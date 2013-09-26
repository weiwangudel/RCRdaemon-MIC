GCC = gcc

#OpenCL_SDK=
#INCLUDE=-I/usr/local/stow/cuda_5.0.35/include/
#LIBPATH=-L${OpenCL_SDK}/OpenCL/common/lib -L${OpenCL_SDK}/shared/lib
LIB= -g -O3 -lOpenCL -lm -lrt

FILES = brdr2d.c  ../RCR.bb.o ../MICEnergyStatDaemon.o

all: brdr2d_gpu

brdr2d_gpu:
	${GCC}  ${LIBPATH} ${INCLUDE} ${LIB} -o brdr2d.exe ${FILES}

clean:
	rm -f *~ *.exe

