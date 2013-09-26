/**
  unfortunate sideffect is to make total execution time nearly a multiple of nanosleep time becasue the
nanosleep must complete before energy counting terminates.

 * RCR logger implementation. 
 * By David O'Brien based on code by Min Yeol Lim 
 *
 * Includes some code based on the UDP server example from the Linux Gazette
 *   (http://www.linuxgazette.com/node/8758)
 * Uses some code based on the libpfm4 examples from Perfmon2
 *   (http://perfmon2.sourceforge.net/)
 */


#include <inttypes.h>  // for int64_t
#include <stdio.h>     // for printf
#include <time.h>      // for clock_gettime

#include "RCR.bb.h"    // for RCRblackboard types

static struct timespec saveTime, finTime;
static struct timespec baseTime;      // the NEVER RESET base Time
double initEnergy;
double baseEnergy;   // the NEVER RESET init Energy 

static struct timespec timerStart;  // designed for polybench: startComputing
static struct timespec timerStop;    // designed for end of computing
double meterStart;                  // energy meter when compute starts
double meterStop;                    // energy meter when compute ends
 
typedef struct energyStatForLoops{
  uint32_t loopID;                  // loop id: obtained from application
  uint64_t numI;                    // number of iterations
				    // (may not know at compiler time, even at
				    // the begining of runtime)
  uint64_t sampleCounter;	    // counts how many times energy was sampled
  double timeSum;                   // total time executing the loop
  double energySum;		    // energy consumed
  double *energyPerIter;            // energy per iteration
} LOOPEnergy;

LOOPEnergy *loopEnergies;           // using array (assume # of loops < N)
double appTotalEnergy;
double appTotalTime;
double appTotalPowerLevel;

#define MAX_NUM_LOOP_ENTER 30
int loop_mapping[MAX_NUM_LOOP_ENTER]; // store the mapping of loop 0
				      // 1, 2, at index [0][1]...
				      // src_line, loop numbe can be 
				      // the contents
int cur_total_num_loops = 0;               // 
int nested_flag = 0;                 // whether nested Enter/Exit has encountered
static void lapPowerEnter(void);
static void lapPowerExit(uint32_t loop);

/******************
 * internal function to print energy consumption 
 *     resets counters -- so multiple calls give increamental usage
 *****************/

static void printPower()
{
  //if ((energy == NULL) || (initEnergy == NULL)) {
  //  printf("Trying to print energy usage before initialization or after termination\n");
  //  return;
  //}
  clock_gettime(CLOCK_MONOTONIC, &finTime);
  double diffTime = (finTime.tv_sec - saveTime.tv_sec)+((finTime.tv_nsec - saveTime.tv_nsec)*10e-10);
  
  double e; 
  readBlackboard(0, &e); 
  printf("e is : %f\n", e);
  double totalEnergy = e - initEnergy;
  initEnergy = e; // reset value

  // update program total
  appTotalEnergy += totalEnergy;
  appTotalTime += diffTime;
  printf("appTotalEnergy:%f\n", appTotalEnergy);  
  // reset values
  saveTime = finTime;
  return;
}

/************************
 * external interface to start measuring energy
 ***********************/

int energyDaemonInit(uint64_t waitTime)
{
   energyDaemon_initBlackboard();
 
  // dealing with loops -- record the energy consumption of different loops
  uint32_t num_loops = MAX_NUM_LOOP_ENTER;
  int l = 0;
  loopEnergies = (LOOPEnergy*)malloc(sizeof(LOOPEnergy) * num_loops); 
  for (l = 0; l < num_loops; l++) {
    loop_mapping[l] = -1;               // loop enter not used yet
    loopEnergies[l].sampleCounter = 0;
    loopEnergies[l].timeSum = 0.0;
    loopEnergies[l].energySum = 0;
  }
  
  appTotalEnergy = 0;   // why initialize to 0? Energy used setting up daemon? 
  appTotalTime = 0;
  appTotalPowerLevel = 0;

  clock_gettime(CLOCK_MONOTONIC, &saveTime);     // Record time 
  baseTime = saveTime;   // permanently record start time
  readBlackboard(0, &initEnergy);
  printf("initial energy is: %f\n", initEnergy);
  baseEnergy = initEnergy;
  return 0;
}

/************************
 * external call to print energy usage since last call
 ***********************/

void energyDaemonPrint()
{
  printPower();
  return;
}

/************************
 * external call to print energy usage since last print call (and cleans up memory)
 ***********************/

void energyDaemonTerm()
{
  int l;
  printPower(); 
  printf("\n");   // some application does not end with \n
  uint32_t num_loops = cur_total_num_loops;
  for (l = 0; l < num_loops; l++) {
    appTotalEnergy += loopEnergies[l].energySum;
    appTotalTime += loopEnergies[l].timeSum;
    printf("Loop %d <line-%d> - Time %f Total energy consumed %f Ave. Power Level %f\n",
	 l, loop_mapping[l], loopEnergies[l].timeSum, 
	 loopEnergies[l].energySum,
	 loopEnergies[l].energySum/loopEnergies[l].timeSum);
  }
  
  double joule =  appTotalEnergy;
  printf("Application(EnergyStat) - Time %f Total energy consumed %f Ave. Power Level %f "
	 "Final Temperature",
	 appTotalTime, joule, joule/appTotalTime);
  printf("\n");
  
  free (loopEnergies);
  initEnergy = 0;
  baseEnergy = 0;
  return;
}

/*
 * Wrapper around lapPowerEnter()
 * External call by user
 */
void energyDaemonEnter(void) {
  if (nested_flag > 0) {
    printf("Nested Enter/Exit not handled yet!\n");
    exit(-1);
  }
  nested_flag++;
  lapPowerEnter();
  return;
}

/* 
 * called before Entering loop to record energy usage of "OpenMP for" loop
 * must be paired with lapPowerExit() to measure energy consumption of a loop
 */ 
static void lapPowerEnter() {
  //if ((energy == NULL) || (initEnergy == NULL)) {
  //  printf("Trying to lap energy usage before initialization or after termination\n");
  //  return;
  //}
  clock_gettime(CLOCK_MONOTONIC, &finTime);
  // Print values
  double diffTime = (finTime.tv_sec - saveTime.tv_sec)+((finTime.tv_nsec - saveTime.tv_nsec)*10e-10);
  
  double e; 
  readBlackboard(0, &e);
  appTotalEnergy += e - initEnergy;
  initEnergy = e; // reset value

  appTotalTime += diffTime; 
  // reset values
  saveTime = finTime;
  return;
}

// wrapper around lapPowerExit();
//Todo: add <file_name, src_line_no> mapping to loop number
void energyDaemonExit(char *file_name, uint32_t src_line_no) {
  int i;

  if (nested_flag != 1) {
    printf("nested enter/exit encountered!\n");
    exit(-1);
  }
  nested_flag--;
 
  for (i=0; i<cur_total_num_loops; i++) {
    // found existing loop already encountered 
    if (loop_mapping[i] == src_line_no) {
	lapPowerExit(i);
	return;
    }
  }

  // new loops to measure power
  loop_mapping[cur_total_num_loops] = src_line_no; 
  lapPowerExit(cur_total_num_loops);
  cur_total_num_loops++; 
  if (cur_total_num_loops >  MAX_NUM_LOOP_ENTER) {
    printf("increase MAX_NUM_LOOP_ENTER in energyStatDaemon please\n");
    exit(-1);
  }

}

/*
 * called immediately after Exiting loop to "lap" the energy 
 * lapPowerEnter() must have been called
 */ 
static void lapPowerExit(uint32_t loop) {
  //if ((energy == NULL) || (initEnergy == NULL)) {
  //  printf("Trying to lap energy usage before initialization or after termination\n");
  //  return;
  //}
  clock_gettime(CLOCK_MONOTONIC, &finTime);
  double diffTime = (finTime.tv_sec - saveTime.tv_sec)+((finTime.tv_nsec - saveTime.tv_nsec)*10e-10);
 
  loopEnergies[loop].sampleCounter++;
 
  double e ;
  readBlackboard(0, &e);							
  loopEnergies[loop].energySum += e - initEnergy;
  initEnergy = e; // reset value

  loopEnergies[loop].timeSum += diffTime; 
  // reset values
  saveTime = finTime;
  return;

}

// start measure time and energy
void energyDaemonTEStart() {
    //if (baseEnergy == NULL) {
    //    printf("Trying to lap energy usage before initialization or after termination\n");
    //    return;
    //}
    clock_gettime(CLOCK_MONOTONIC, &timerStart);
    double e ; readBlackboard(0, &e);							
    meterStart = e;
}

// stop the watch and energy meter, compute difference to energyDaemonStart
void energyDaemonTEStop() {
    //if (baseEnergy == NULL) {
    //    printf("Trying to lap energy usage before initialization or after termination\n");
    //    return;
    //}
    clock_gettime(CLOCK_MONOTONIC, &timerStop);
    double e ; readBlackboard(0, &e);							
    meterStop = e;

    double diffTime = (timerStop.tv_sec - timerStart.tv_sec)+((timerStop.tv_nsec - timerStart.tv_nsec)*10e-10);
    double joule =  (meterStop - meterStart) ; 
    
    printf("\nRegional(EnergyStat) - Time %f Total energy consumed %f Ave. Power Level %f \n",
	 diffTime, joule, joule/diffTime);
}
