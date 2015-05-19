#Makefile for 2-1/2D Darwin MPI/OpenMP PIC codes

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
#OPTS77 = -O3 -fcheck=bounds -fdefault-real-8 -fdefault-double-8 -Wall -std=f95

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
#OPTS77 = -O3 -fbounds-check -fdefault-real-8 -fdefault-double-8 -Wall -std=f95

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

all : fmpdpic2 cmpdpic2

special: fmpdpic2_c cmpdpic2_f

# Version using Fortran77 mpplib2.f
#fmpdpic2 : fmpdpic2.o fmpdpush2.o fmpplib2.o fomplib.o dtimer.o
#	$(MPIFC) $(OPTS90) $(LOPTS) -o fmpdpic2 \
#   fmpdpic2.o fmpdpush2.o fmpplib2.o fomplib.o mpdpush2_h.o mpplib2_h.o \
#   omplib_h.o dtimer.o

# Version using Fortran90 mpplib2.f90
fmpdpic2 : fmpdpic2.o fmpdpush2.o fmpplib2.o fomplib.o dtimer.o
	$(MPIFC) $(OPTS90) $(LOPTS) -o fmpdpic2 \
    fmpdpic2.o fmpdpush2.o f90mpplib2.o fomplib.o mpdpush2_h.o omplib_h.o \
    dtimer.o

cmpdpic2 : cmpdpic2.o cmpdpush2.o cmpplib2.o complib.o dtimer.o
	$(MPICC) $(CCOPTS) $(LOPTS) -o cmpdpic2 \
    cmpdpic2.o cmpdpush2.o cmpplib2.o complib.o dtimer.o

fmpdpic2_c : fmpdpic2_c.o cmpdpush2.o cmpplib2.o complib.o dtimer.o
	$(MPIFC) $(OPTS90) $(LOPTS) -o fmpdpic2_c \
    fmpdpic2_c.o cmpdpush2.o cmpplib2.o complib.o dtimer.o

cmpdpic2_f : cmpdpic2.o cmpdpush2_f.o cmpplib2_f.o complib_f.o fmpdpush2.o \
             fmpplib2.o fomplib.o dtimer.o
	$(MPIFC) $(OPTS90) $(LOPTS) $(LEGACY) -o cmpdpic2_f \
    cmpdpic2.o cmpdpush2_f.o cmpplib2_f.o complib_f.o fmpdpush2.o fmpplib2.o \
    fomplib.o dtimer.o

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

fmpdpush2.o : mpdpush2.f
	$(MPIFC) $(OPTS90) -o fmpdpush2.o -c mpdpush2.f

mpdpush2_h.o : mpdpush2_h.f90
	$(MPIFC) $(OPTS90) -o mpdpush2_h.o -c mpdpush2_h.f90

cmpdpush2.o : mpdpush2.c
	$(MPICC) $(CCOPTS) -o cmpdpush2.o -c mpdpush2.c

cmpdpush2_f.o : mpdpush2_f.c
	$(MPICC) $(CCOPTS) -o cmpdpush2_f.o -c mpdpush2_f.c

# Version using Fortran77 mpplib2.f
#fmpdpic2.o : mpdpic2.f90 mpdpush2_h.o mpplib2_h.o omplib_h.o
#	$(MPIFC) $(OPTS90) -o fmpdpic2.o -c mpdpic2.f90

# Version using Fortran90 mpplib2.f90
fmpdpic2.o : mpdpic2.f90 f90mpplib2.o mpdpush2_h.o omplib_h.o
	$(MPIFC) $(OPTS90) -o fmpdpic2.o -c mpdpic2.f90

cmpdpic2.o : mpdpic2.c
	$(MPICC) $(CCOPTS) -o cmpdpic2.o -c mpdpic2.c

fmpdpic2_c.o : mpdpic2_c.f90
	$(MPIFC) $(OPTS90) -o fmpdpic2_c.o -c mpdpic2_c.f90

clean :
	rm -f *.o *.mod

clobber: clean
	rm -f fmpdpic2 cmpdpic2 fmpdpic2_c cmpdpic2_f
