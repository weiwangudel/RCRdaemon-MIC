#!/bin/bash
papi_src=/eecis/cavazos/wwang/papi-5.2.0-micpower-works/src

x86_64-k1om-linux-gcc -g -c MICEnergyStatDaemon.c
x86_64-k1om-linux-gcc -g -c RCR.bb.c
x86_64-k1om-linux-gcc -g -I$papi_src/testlib -I$papi_src -I. -c -o micpower_basic.o micpower_basic.c
x86_64-k1om-linux-gcc -g -I$papi_src/testlib -I$papi_src -I. -c -o RCRdaemonSB.o RCRdaemonSB.c 
x86_64-k1om-linux-gcc -g -I$papi_src/testlib -I$papi_src -I. -o RCRMICPowerDaemon RCR.bb.o micpower_basic.o RCRdaemonSB.o $papi_src/testlib/do_loops.o $papi_src/testlib/test_utils.o $papi_src/testlib/dummy.o $papi_src/libpapi.a -ldl -lpthread -lrt 
