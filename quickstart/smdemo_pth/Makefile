# Makefile for Pthreads test programs

# Makefile gfortran compiler with MacOS X

#FC90 = gfortran
#CC = gcc

#OPTS90 = -O3
#CCOPTS = -O3

#PTHREADS
#MPOBJS = LnxMP.o
#PTLIBS = $(CARBON) -lSystemStubs -lpthread

# Makefile gfortran compiler with Linux

FC90 = gfortran
CC = gcc

OPTS90 = -O3
CCOPTS = -O3

#PTHREADS
MPOBJS = LnxMP.o LPProcessors.o
PTLIBS = -lpthread

#

PTOBJS = cpthlib.o $(MPOBJS) dtimer.o

FPTOBJS = fpthlib.o $(MPOBJS) dtimer.o

# Linkage rules

all: pthreads

pthreads: cpthtest fpthtest

#PTHREADS
cpthtest : cpthtest.o $(PTOBJS)
	$(CC) $(CCOPTS) $(LOPTS) -o cpthtest cpthtest.o \
      $(PTOBJS) $(PTLIBS)

fpthtest : fpthtest.o $(FPTOBJS)
	$(FC90) $(OPTS90) $(LOPTS) -o fpthtest fpthtest.o \
      $(FPTOBJS) $(PTLIBS)

# Compilation rules

dtimer.o : dtimer.c
	$(CC) $(CCOPTS) -c dtimer.c

#PTHREADS
LnxMP.o : LnxMP.c
	$(CC) $(CCOPTS) -c LnxMP.c

LPProcessors.o : LPProcessors.c
	$(CC) $(CCOPTS) -c LPProcessors.c

cpthlib.o : pthlib.c
	$(CC) $(CCOPTS) -o cpthlib.o -c pthlib.c

fpthlib.o : pthlib.f
	$(FC90) $(OPTS90) -o fpthlib.o -c pthlib.f

cpthtest.o : cpthtest.c
	$(CC) $(CCOPTS) -c cpthtest.c

fpthtest.o : fpthtest.f90
	$(FC90) $(OPTS90) -c fpthtest.f90

clean:
	rm -f *.o *.mod

clobber: clean
	rm -f *pthtest
