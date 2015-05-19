# Makefile for MPI test programs

# Makefile gfortran compiler with MacOS X

CARBON = /System/Library/Frameworks/Carbon.framework/Carbon

#FC90 = gfortran
#CC = gcc

#OPTS90 = -O3
#CCOPTS = -O3

#MPI
#LOPTS =

# Makefile gfortran compiler with Linux

FC90 = gfortran
CC = gcc

OPTS90 = -O3
CCOPTS = -O3

#MPI
LOPTS =

#

MPIOBJS = cmpimain.o dtimer.o

FMPIOBJS = fmpimain.o dtimer.o

# Linkage rules

all: mpi

mpi: cmpitest fmpitest

#MPI
cmpitest : cmpitest.o $(MPIOBJS)
	mpicc $(CCOPTS) $(LOPTS) -o cmpitest cmpitest.o \
      $(MPIOBJS)

fmpitest : fmpitest.o $(FMPIOBJS)
	mpif90 $(OPTS90) $(LOPTS) -o fmpitest fmpitest.o \
      $(FMPIOBJS)

# Compilation rules

dtimer.o : dtimer.c
	$(CC) $(CCOPTS) -c dtimer.c

#MPI
cmpimain.o : mpimain.c
	mpicc $(CCOPTS) -o cmpimain.o -c mpimain.c

fmpimain.o : mpimain.f
	mpif90 $(OPTS90) -o fmpimain.o -c mpimain.f

cmpitest.o : cmpitest.c
	mpicc $(CCOPTS) -c cmpitest.c

fmpitest.o : fmpitest.f90
	mpif90 $(OPTS90) -c fmpitest.f90

clean:
	rm -f *.o *.mod

clobber: clean
	rm -f *mpitest
