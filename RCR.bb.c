#include <stdint.h> // for int64_t
#include <stdio.h>  // for printf
#include <stdlib.h>  // for malloc - temp until share memory region allocated
#include <sys/mman.h> // for mmap
#include <fcntl.h>   // for mmap #defs
#include <unistd.h>   // for ftruncate
#include <sys/ioctl.h>
#include <string.h>

#define RCRFILE_NAME  "RCRMICFile"
#define MAX_RCRFILE_SIZE 1024

void *bbMem;
int allocationHighWater;

struct PAPI_MIC_COUNTERS {
        int total0;
        int total1;
        int inst;
        int imax;
        int pcie;
        int c2x3;
        int c2x4;
        int vccp;
        int vddg;
        int vddq;
};

char * base() {
	return (char*) bbMem;
}

// allocate next block of shared memory -- track highwater mark

int64_t allocateSharedMemory(int64_t size) {
	int64_t ret = allocationHighWater;
	if (size == 0) size = 8; // put blank word in if no space required to prevent field overlaps
	if (size + allocationHighWater > MAX_RCRFILE_SIZE) {
		fprintf(stderr, " Out of shared memory\n");
		exit(1);
	}
	if (ret == 0) { // no allocation invent some memory
		int fd = shm_open(RCRFILE_NAME, O_RDWR | O_CREAT, 0);  // TODO should be #def or execution time definable
		if (fd == -1) {
			perror("shm_open");
			exit(1);
		}

		ftruncate(fd, MAX_RCRFILE_SIZE);
		bbMem = mmap(NULL, MAX_RCRFILE_SIZE, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);

	}
	allocationHighWater += size;
	return ret;
}

// get current allocation offset -- next allocation starting offset

int64_t getCurOffset() {
	return allocationHighWater;
}

/******************************
 *  reads current blackboard file and sets it up for access 
 *****************************/
int initBlackBoard() {
	int ret = -1; 

	int fd = shm_open(RCRFILE_NAME, O_RDONLY, 0);

	if (fd == -1) {
		perror("shm_open "RCRFILE_NAME);
		return -1;
	}

	ftruncate(fd, MAX_RCRFILE_SIZE);
	bbMem = mmap(NULL, MAX_RCRFILE_SIZE, PROT_READ, MAP_SHARED, fd, 0);

	return 1;
}

int64_t buildSystem(int64_t numBBMeters, int64_t numNodes) {
	int64_t i;
	int64_t offset = allocateSharedMemory(sizeof(struct PAPI_MIC_COUNTERS));
	return offset;
}

