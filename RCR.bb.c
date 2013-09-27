#include <stdint.h> // for int64_t
#include <stdio.h>  // for printf
#include <stdlib.h>  // for malloc - temp until share memory region allocated
#include <sys/mman.h> // for mmap
#include <fcntl.h>   // for mmap #defs
#include <unistd.h>   // for ftruncate
#include <sys/ioctl.h>
#include <string.h>
#include "RCR.bb.h"

#define RCRFILE_NAME  "RCRMICFile"
#define MAX_RCRFILE_SIZE 1024

void *bbMem;
int allocationHighWater;


char * base() {
	return (char*) bbMem;
}

// allocate next block of shared memory -- track highwater mark

int64_t allocateSharedMemory(int64_t size) {
    int fd = shm_open(RCRFILE_NAME, O_RDWR | O_CREAT, 0);  // TODO should be #def or execution time definable
    if (fd == -1) {
    	perror("shm_open");
    	exit(1);
    }
    
    ftruncate(fd, MAX_RCRFILE_SIZE);
    bbMem = mmap(NULL, MAX_RCRFILE_SIZE, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
    
    printf("in the allocateSharedMemory:%#x\n", bbMem);
	allocationHighWater += size;
    printf("allocationHighWater:%d\n", allocationHighWater);
//    initBlackboard();
    return allocationHighWater;
}

// get current allocation offset -- next allocation starting offset

int64_t getCurOffset() {
	return allocationHighWater;
}

int64_t buildBlackboard() {
	int64_t offset = allocateSharedMemory(sizeof(struct PAPI_MIC_COUNTERS));
	return offset;
}


void initBlackboard() {
    struct PAPI_MIC_COUNTERS *node = (struct PAPI_MIC_COUNTERS *) (base());
    node->total0 = 0.0;
    printf("node->tot0:%f\n", node->total0);
}

void update_blackboard(char* e_name, long long value, long long interval) {
 //*(*)e.() += value*(1e-6)*interval*(1e-6);
  //printf("in reading: %#x:\n", base());
  struct PAPI_MIC_COUNTERS* node = (struct PAPI_MIC_COUNTERS*)(base() );
  //printf("in reading: %#x:\n", node);
  ////printf("energy increment: %f %f\n", interval*(1.0e-6), value*(1.0e-6)*interval*(1.0e-6));
  node->total0 += value*(1.0e-6)*interval*(1.0e-6); 
  printf("accumulative: %f\n", node->total0);
}
double readBlackboard(unsigned int counter, double *value) { 
  struct PAPI_MIC_COUNTERS* node = (struct PAPI_MIC_COUNTERS*)(base() );
  //printf("in readBlackboard:%f\n", node->total0);
  *value = node->total0;
  return *value;
}


void energyDaemon_initBlackboard() {
  int fd = shm_open(RCRFILE_NAME, O_RDONLY, 0);                                  
                                                                                 
  if(fd == -1) {                                                                 
    perror("shm_open "RCRFILE_NAME);                                             
    printf("Please create it at :/dev/shm/RCRMICFile\n");                           
    exit(-1);
  }                                                                              
  ftruncate(fd, MAX_RCRFILE_SIZE);                                               
  bbMem = mmap(NULL, MAX_RCRFILE_SIZE, PROT_READ, MAP_SHARED, fd, 0);            
}
