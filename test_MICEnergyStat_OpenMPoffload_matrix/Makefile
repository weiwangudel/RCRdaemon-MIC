all:
	/opt/intel/bin/icc -mkl -O3 -openmp -openmp-offload -Wno-unknown-pragmas -std=c99 -vec-report3 matrix.c ../MICEnergyStatDaemon.o ../RCR.bb.o -o matrix.off  -lrt
