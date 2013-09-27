#!/bin/bash
papi_src=/eecis/cavazos/wwang/papi-5.2.0-micpower-works/src

x86_64-k1om-linux-gcc -g -c MICEnergyStatDaemon.c
x86_64-k1om-linux-gcc -g -c RCR.bb.c
x86_64-k1om-linux-gcc -g -DSTATIC_PAPI_EVENTS_TABLE -DPEINCLUDE="libpfm4/include/perfmon/perf_event.h" -D_REENTRANT -D_GNU_SOURCE -DUSE_COMPILER_TLS  -Wall -Ilibpfm4/include -Wextra -D MICACCESSAPI -D LINUX -I/opt/intel/mic/sysmgmt/sdk/include -I$papi_src/testlib -I$papi_src -I. -c -o host_micpower_basic.o host_micpower_basic.c
x86_64-k1om-linux-gcc -g -DSTATIC_PAPI_EVENTS_TABLE -DPEINCLUDE="libpfm4/include/perfmon/perf_event.h" -D_REENTRANT -D_GNU_SOURCE -DUSE_COMPILER_TLS  -Wall -Ilibpfm4/include -Wextra -D MICACCESSAPI -D LINUX -I/opt/intel/mic/sysmgmt/sdk/include -I$papi_src/testlib -I$papi_src -I. -c -o RCRdaemonSB.o RCRdaemonSB.c 
x86_64-k1om-linux-gcc -g -DSTATIC_PAPI_EVENTS_TABLE -DPEINCLUDE="libpfm4/include/perfmon/perf_event.h" -D_REENTRANT -D_GNU_SOURCE -DUSE_COMPILER_TLS  -Wall -Ilibpfm4/include -Wextra -D MICACCESSAPI -D LINUX -I/opt/intel/mic/sysmgmt/sdk/include -I$papi_src/testlib -I$papi_src -I. -o RCRMICPowerDaemon RCR.bb.o host_micpower_basic.o RCRdaemonSB.o $papi_src/testlib/do_loops.o $papi_src/testlib/test_utils.o $papi_src/testlib/dummy.o $papi_src/libpapi.a -ldl -lpthread -lrt 