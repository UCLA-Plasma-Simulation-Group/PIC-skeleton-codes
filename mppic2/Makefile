#Makefile for 2D Electrostatic MPI/OpenMP PIC codes

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
#OPTS77 = -O3 -fcheck=bounds -fdefault-real-8 -Wall

#CCOPTS = -O3 -std=c99
#CCOPTS = -O3 -Wall -std=c99
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
MPIFC = mpiifort -openmp
MPICC = mpiicc -openmp

CC = icc

OPTS90 = -O3
#OPTS90 = -O3 -r8
#OPTS90 = -O3 -CB -r8 -warn all -std90

OPTS77 = -O3
#OPTS77 = -O3 -r8
#OPTS77 = -O3 -CB -r8 -warn all -std77

CCOPTS = -O3 -std=c99
#CCOPTS = -O3 -no-vec -Wall -std=c99
LEGACY = -nofor_main

# Makefile gfortran compiler with Linux

# OpenMPI
#MPIFC = mpif90 -fopenmp
#MPICC = mpicc -fopenmp

#CC = gcc

#OPTS90 = -O3
#OPTS90 = -O3 -fdefault-real-8 -fdefault-double-8
#OPTS90 = -O3 -fbounds-check -fdefault-real-8 -fdefault-double-8 -Wall -std=f95

#OPTS77 = -O3
#OPTS77 = -O3 -fdefault-real-8 -fdefault-double-8
#OPTS77 = -O3 -fbounds-check -fdefault-real-8 -fdefault-double-8 -Wall

#CCOPTS = -O3 -std=c99
#CCOPTS = -O3 -Wall -std=c99
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

all : fmppic2 cmppic2

special: fmppic2_c cmppic2_f

# Version using Fortran77 mpplib2.f
#fmppic2 : fmppic2.o fmppush2.o fmpplib2.o fomplib.o dtimer.o
#	$(MPIFC) $(OPTS90) $(LOPTS) -o fmppic2 \
#   fmppic2.o fmppush2.o fmpplib2.o fomplib.o mppush2_h.o mpplib2_h.o \
#   omplib_h.o dtimer.o

# Version using Fortran90 mpplib2.f90
fmppic2 : fmppic2.o fmppush2.o fmpplib2.o fomplib.o dtimer.o
	$(MPIFC) $(OPTS90) $(LOPTS) -o fmppic2 \
    fmppic2.o fmppush2.o f90mpplib2.o fomplib.o mppush2_h.o omplib_h.o \
    dtimer.o

cmppic2 : cmppic2.o cmppush2.o cmpplib2.o complib.o dtimer.o
	$(MPICC) $(CCOPTS) $(LOPTS) -o cmppic2 \
    cmppic2.o cmppush2.o cmpplib2.o complib.o dtimer.o

fmppic2_c : fmppic2_c.o cmppush2.o cmpplib2.o complib.o dtimer.o
	$(MPIFC) $(OPTS90) $(LOPTS) -o fmppic2_c \
    fmppic2_c.o cmppush2.o cmpplib2.o complib.o dtimer.o

cmppic2_f : cmppic2.o cmppush2_f.o cmpplib2_f.o complib_f.o fmppush2.o \
            fmpplib2.o fomplib.o dtimer.o
	$(MPIFC) $(OPTS90) $(LOPTS) $(LEGACY) -o cmppic2_f \
    cmppic2.o cmppush2_f.o cmpplib2_f.o complib_f.o fmppush2.o \
    fmpplib2.o fomplib.o dtimer.o

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

fmpplib2.o : mpplib2.f
	$(MPIFC) $(OPTS77) -o fmpplib2.o -c mpplib2.f

mpplib2_h.o : mpplib2_h.f90
	$(MPIFC) $(OPTS90) -o mpplib2_h.o -c mpplib2_h.f90

f90mpplib2.o : mpplib2.f90
	$(MPIFC) $(OPTS90) -o f90mpplib2.o -c mpplib2.f90

cmpplib2.o : mpplib2.c
	$(MPICC) $(CCOPTS) -o cmpplib2.o -c mpplib2.c

cmpplib2_f.o : mpplib2_f.c
	$(MPICC) $(CCOPTS) -o cmpplib2_f.o -c mpplib2_f.c

fmppush2.o : mppush2.f
	$(MPIFC) $(OPTS90) -o fmppush2.o -c mppush2.f

mppush2_h.o : mppush2_h.f90
	$(MPIFC) $(OPTS90) -o mppush2_h.o -c mppush2_h.f90

cmppush2.o : mppush2.c
	$(MPICC) $(CCOPTS) -o cmppush2.o -c mppush2.c

# Version using Fortran77 mpplib2.f
#fmppic2.o : mppic2.f90 mppush2_h.o mpplib2_h.o omplib_h.o
#	$(MPIFC) $(OPTS90) -o fmppic2.o -c mppic2.f90

# Version using Fortran90 mpplib2.f90
fmppic2.o : mppic2.f90 f90mpplib2.o mppush2_h.o omplib_h.o
	$(MPIFC) $(OPTS90) -o fmppic2.o -c mppic2.f90

cmppush2_f.o : mppush2_f.c
	$(MPICC) $(CCOPTS) -o cmppush2_f.o -c mppush2_f.c

cmppic2.o : mppic2.c
	$(MPICC) $(CCOPTS) -o cmppic2.o -c mppic2.c

fmppic2_c.o : mppic2_c.f90
	$(MPIFC) $(OPTS90) -o fmppic2_c.o -c mppic2_c.f90

clean :
	rm -f *.o *.mod

clobber: clean
	rm -f fmppic2 cmppic2 fmppic2_c cmppic2_f
