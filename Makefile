#
# Copyright (C) Intel Corporation (2012)
#
# This file is subject to the Intel Sample Source Code License.
# A copy of the Intel Sample Source Code License is included. 
#
# Linux OS:	/opt/intel/mic/LICENSE/
# Windows:	C:\Program Files\Intel\MPSS\
#
#


DEFINES += -D MICACCESSAPI -D LINUX
CFLAGS += -m64 -g -O0 $(DEFINES)
LIBPATH = -L/opt/intel/mic/sysmgmt/sdk/lib/Linux

LDFLAGS = $(LIBPATH) $(SCIF_LIBPATH) -lMicAccessSDK -lscif -lpthread

INCDIR = -I/opt/intel/mic/sysmgmt/sdk/include

SOURCES= example-powerusage.c example-gettemperature.c

OBJS = $(SOURCES:.c=.o)
BINARIES = $(SOURCES:.c=)

all: $(SOURCES) $(BINARIES)

$(BINARIES): $(OBJS)
	$(CC) $@.o -o $@ $(LDFLAGS)

.cpp.o:
	$(CC) $(CFLAGS) $(INCDIR) -c $< -o $@

.c.o:
	$(CC) $(CFLAGS) $(INCDIR) -c $< -o $@

clean:
	$(RM) $(OBJS) $(BINARIES)
