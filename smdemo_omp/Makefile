# Makefile for OpenMP test program

# Makefile gfortran compiler with MacOS X

#FC90 = gfortran
#CC = gcc

#OPTS90 = -O3 -fopenmp
#CCOPTS = -O3 -fopenmp

#LOPTS = -lgomp

# Makefile gfortran compiler with Linux

FC90 = gfortran
CC = gcc

OPTS90 = -O3 -fopenmp
CCOPTS = -O3 -fopenmp

#LOPTS = -lgomp

#

OMPOBJS = complib.o dtimer.o

FOMPOBJS = fomplib.o dtimer.o

# Linkage rules

all: openmp

openmp: comptest fomptest

#OPENMP
comptest : comptest.o $(OMPOBJS)
	$(CC) $(CCOPTS) -o comptest comptest.o $(OMPOBJS)

fomptest : fomptest.o $(FOMPOBJS)
	$(FC90) $(OPTS90) -o fomptest fomptest.o $(FOMPOBJS)

# Compilation rules

dtimer.o : dtimer.c
	$(CC) $(CCOPTS) -c dtimer.c

#OPENMP
complib.o : omplib.c
	$(CC) $(CCOPTS) -o complib.o -c omplib.c

fomplib.o : omplib.f
	$(FC90) $(OPTS90) -o fomplib.o -c omplib.f

comptest.o : comptest.c
	$(CC) $(CCOPTS) -c comptest.c

fomptest.o : fomptest.f90
	$(FC90) $(OPTS90) -c fomptest.f90

clean:
	rm -f *.o *.mod

clobber: clean
	rm -f *omptest
