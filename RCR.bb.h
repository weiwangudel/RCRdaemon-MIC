#include <stdio.h>
int64_t buildBlackboard();
void updateBlackboard();
void initBlackboard() ;
void energyDaemon_initBlackboard() ;
struct PAPI_MIC_COUNTERS {
        double total0;
        double total1;
        double inst;
        double imax;
        double pcie;
        double c2x3;
        double c2x4;
        double vccp;
        double vddg;
        double vddq;
};
