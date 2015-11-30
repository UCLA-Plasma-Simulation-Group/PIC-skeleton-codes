#Makefile for 3D Electromagnetic MPI/OpenMP PIC codes

# Makefile gfortran compiler with MacOS X

#OpenMPI
#MPIFC = mpif90 -fopenmp
#MPICC = mpicc -fopenmp

#CC = gcc

#OPTS90 = -O3
#OPTS90 = -O3 -fdefault-real-8 -fdefault-double-8
#OPTS90 = -O3 -fcheck=bounds -fdefault-real-8 -fdefault-double-8 -Wall -std=f95

#OPTS77 = -O3
#OPTS77 = -O3 -fdefault-real-8 -fdefault-double-8
#OPTS77 = -O3 -fcheck=bounds -fdefault-real-8 -fdefault-double-8 -Wall

#CCOPTS = -O3 -std=c99
#CCOPTS = -O3 -Wall -std=c99

#LOPTS = -lgomp
#LEGACY =

# Makefile Intel compiler with Mac OS X

# OpenMPI
#MPIFC = mpif90 -openmp
#MPICC = mpicc -fopenmp

#CC = gcc

#OPTS90 = -O3
#OPTS90 = -O3 -r8
#OPTS90 = -O3 -CB -r8 -warn all -std90

#OPTS77 = -O3
#OPTS77 = -O3 -r8
#OPTS77 = -O3 -CB -r8 -warn all -std77

#CCOPTS = -O3 -std=c99
#CCOPTS = -O3 -Wall -std=c99
#LEGACY = -nofor_main

# Makefile Intel compiler with Linux

# MPI
#MPIFC = mpiifort -openmp
#MPICC = mpiicc -openmp

#CC = icc

#OPTS90 = -O3
#OPTS90 = -O3 -r8
#OPTS90 = -O3 -CB -r8 -warn all -std90

#OPTS77 = -O3
#OPTS77 = -O3 -r8
#OPTS77 = -O3 -CB -r8 -warn all -std77

#CCOPTS = -O3 -std=c99
#CCOPTS = -O3 -no-vec -Wall -std=c99
#LEGACY = -nofor_main

# Makefile gfortran compiler with Linux

# OpenMPI
MPIFC = mpif90 -fopenmp
MPICC = mpicc -fopenmp

CC = gcc

OPTS90 = -O3
#OPTS90 = -O3 -fdefault-real-8 -fdefault-double-8
#OPTS90 = -O3 -fbounds-check -fdefault-real-8 -fdefault-double-8 -Wall -std=f95

OPTS77 = -O3
#OPTS77 = -O3 -fdefault-real-8 -fdefault-double-8
#OPTS77 = -O3 -fbounds-check -fdefault-real-8 -fdefault-double-8 -Wall

CCOPTS = -O3 -std=c99
#CCOPTS = -O3 -Wall -std=c99

#LOPTS = -lgomp
#LEGACY =

# Makefile PGI compiler with Linux

# OpenMPI
#MPIFC = mpipgf90
#MPICC = mpicc -fopenmp

#CC = gcc

#OPTS90 = -O3
#OPTS90 = -O3 -r8
#OPTS90 = -O3 -Mbounds -r8 -Mstandard

#OPTS77 = -O3
#OPTS77 = -O3 -r8
#OPTS77 = -O3 -Mbounds -r8 -Mstandard

#CCOPTS = -O3 -std=c99
#CCOPTS = -O3 -Wall -std=c99
#LEGACY = -Mnomain

# Makefile Cray compiler with Linux

# MPI
#MPIFC = ftn -h omp
#MPICC = cc -h omp

#FC90 = ftn
#CC = cc

#OPTS90 = -O 3
#OPTS90 = -O 3 -s real64
#OPTS90 = -O 3 -R b -s real64 -en

#OPTS77 = -O 3
#OPTS77 = -O 3 -s real64
#OPTS77 = -O 3 -R b -s real64 -en

#CCOPTS = -O 3 -h c99 -h conform
#LEGACY =

#

# Linkage rules

all : fmpbpic3

special: cmpbpic3_f

# Version using Fortran77 mpplib3.f
#fmpbpic3 : fmpbpic3.o fmpbpush3.o fmpplib3.o fomplib.o dtimer.o
#	$(MPIFC) $(OPTS90) $(LOPTS) -o fmpbpic3 \
#   fmpbpic3.o fmpbpush3.o fmpplib3.o fomplib.o mpbpush3_h.o mpplib3_h.o \
#   omplib_h.o dtimer.o

# Version using Fortran90 mpplib3.f90
fmpbpic3 : fmpbpic3.o fmpbpush3.o fmpplib3.o fomplib.o dtimer.o
	$(MPIFC) $(OPTS90) $(LOPTS) -o fmpbpic3 \
    fmpbpic3.o fmpbpush3.o f90mpplib3.o fomplib.o mpbpush3_h.o omplib_h.o \
    dtimer.o

cmpbpic3_f : cmpbpic3.o cmpbpush3_f.o cmpplib3_f.o complib_f.o fmpbpush3.o \
             fmpplib3.o fomplib.o dtimer.o
	$(MPIFC) $(OPTS90) $(LOPTS) $(LEGACY) -o cmpbpic3_f \
    cmpbpic3.o cmpbpush3_f.o cmpplib3_f.o complib_f.o fmpbpush3.o \
    fmpplib3.o fomplib.o dtimer.o

# Compilation rules

dtimer.o : dtimer.c
	$(CC) $(CCOPTS) -c dtimer.c

#OPENMP
fomplib.o : omplib.f
	$(MPIFC) $(OPTS90) -o fomplib.o -c omplib.f

omplib_h.o : omplib_h.f90
	$(MPIFC) $(OPTS90) -o omplib_h.o -c omplib_h.f90

complib.o : omplib.c
	$(MPICC) $(CCOPTS) -o complib.o -c omplib.c

complib_f.o : omplib_f.c
	$(MPICC) $(CCOPTS) -o complib_f.o -c omplib_f.c

fmpplib3.o : mpplib3.f
	$(MPIFC) $(OPTS77) -o fmpplib3.o -c mpplib3.f

mpplib3_h.o : mpplib3_h.f90
	$(MPIFC) $(OPTS90) -o mpplib3_h.o -c mpplib3_h.f90

f90mpplib3.o : mpplib3.f90
	$(MPIFC) $(OPTS90) -o f90mpplib3.o -c mpplib3.f90

cmpplib3_f.o : mpplib3_f.c
	$(MPICC) $(CCOPTS) -o cmpplib3_f.o -c mpplib3_f.c

fmpbpush3.o : mpbpush3.f
	$(MPIFC) $(OPTS90) -o fmpbpush3.o -c mpbpush3.f

mpbpush3_h.o : mpbpush3_h.f90
	$(MPIFC) $(OPTS90) -o mpbpush3_h.o -c mpbpush3_h.f90

# Version using Fortran77 mpplib3.f
#fmpbpic3.o : mpbpic3.f90 mpbpush3_h.o mpplib3_h.o omplib_h.o
#	$(MPIFC) $(OPTS90) -o fmpbpic3.o -c mpbpic3.f90

# Version using Fortran90 mpplib3.f90
fmpbpic3.o : mpbpic3.f90 f90mpplib3.o mpbpush3_h.o omplib_h.o
	$(MPIFC) $(OPTS90) -o fmpbpic3.o -c mpbpic3.f90

cmpbpush3_f.o : mpbpush3_f.c
	$(MPICC) $(CCOPTS) -o cmpbpush3_f.o -c mpbpush3_f.c

cmpbpic3.o : mpbpic3.c
	$(MPICC) $(CCOPTS) -o cmpbpic3.o -c mpbpic3.c

clean :
	rm -f *.o *.mod

clobber: clean
	rm -f fmpbpic3 cmpbpic3_f
