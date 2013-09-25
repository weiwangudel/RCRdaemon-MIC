#include <stdint.h> // for int64_t
#include <stdio.h>  // for printf
#include <time.h>   // for time functions
#include "host_micpower_basic.h"
#include "RCR.bb.h"

volatile int64_t ** tempMeter; // used to store location of current socket tempature
volatile int64_t ** energyMeter; // used to store location of energy used to date

uint64_t *energyWrap;
uint64_t *energySave;

struct timespec adaptTimeStart;  // used for measured time
struct timespec adaptTimeStop;  // used for measured time


int main(int argc, char *argv[]) {
	int logPrint = 0;

	if (argc > 1) {
		printf("enable log printing\n");
		logPrint = 1;
	}
	else {
		printf("disable log printing\n");
		logPrint = 0;
	}

	//building shared memory region
    buildBlackboard();
    
    //initializing shared memory region
    initBlackboard();

	struct timespec interval, freq, remainder;
	interval.tv_sec = 0;
	interval.tv_nsec = 1000000;  // currently aim for ~1/1000 second

	//  Overhead of this loop -- ~1% if only the nanosleep
	//                           ~6.7-7.3% if Power Check called  
	//  overhead in PowerCheck seems to be in the pread calls to read the msr's
	//       read only one at a time and that seems to be enforced by the OS
	//       if moved into OS this could be made MUCH lower but I'm not sure I want to 
	//       do that with LINUX(and I'll never get it approved by the powers that run LINUX --
	//       maybe Kitten/LXK

//	while (1) {
		clock_gettime(CLOCK_MONOTONIC, &adaptTimeStart);
		
        RCRMICPowerCheck();

		clock_gettime(CLOCK_MONOTONIC, &adaptTimeStop);

		if ((interval.tv_sec > 0) || (interval.tv_nsec > 0)) {
			nanosleep(&interval, &remainder);
		}

		if (logPrint) {
//			int64_t temp0 = *tempMeter[0];
//			int64_t temp1 = *tempMeter[1];
//			uint64_t energy0 = *energyMeter[0];
//			uint64_t energy1 = *energyMeter[1];
		}
//	}

	return 0;
}
