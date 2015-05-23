# Makefile for MPI test programs

# Makefile gfortran compiler with MacOS X

#MPIFC = mpif90
#MPICC = mpicc

#OPTS90 = -O3
#CCOPTS = -O3

# Makefile gfortran compiler with Linux

MPIFC = mpif90
MPICC = mpicc

OPTS90 = -O3
CCOPTS = -O3

#

MPIOBJS = cmpilib.o dtimer.o

FMPIOBJS = fmpilib.o dtimer.o

# Linkage rules

all: mpi

mpi: cmpitest fmpitest

#MPI
cmpitest : cmpitest.o $(MPIOBJS)
	$(MPICC) $(CCOPTS) $(LOPTS) -o cmpitest cmpitest.o \
      $(MPIOBJS)

fmpitest : fmpitest.o $(FMPIOBJS)
	$(MPIFC) $(OPTS90) $(LOPTS) -o fmpitest fmpitest.o \
      $(FMPIOBJS)

# Compilation rules

dtimer.o : dtimer.c
	$(CC) $(CCOPTS) -c dtimer.c

#MPI
cmpilib.o : mpilib.c
	$(MPICC) $(CCOPTS) -o cmpilib.o -c mpilib.c

fmpilib.o : mpilib.f
	$(MPIFC) $(OPTS90) -o fmpilib.o -c mpilib.f

cmpitest.o : cmpitest.c
	$(MPICC) $(CCOPTS) -c cmpitest.c

fmpitest.o : fmpitest.f90
	$(MPIFC) $(OPTS90) -c fmpitest.f90

clean:
	rm -f *.o *.mod

clobber: clean
	rm -f *mpitest
