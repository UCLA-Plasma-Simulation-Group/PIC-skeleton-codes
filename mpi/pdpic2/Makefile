#Makefile for 2-1/2D Darwin MPI PIC codes

# Makefile gfortran compiler with MacOS X

#OpenMPI
#MPIFC = mpif90
#MPICC = mpicc

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
#LOPTS = -lm

# Makefile Intel compiler with Mac OS X

# OpenMPI
#MPIFC = mpif90
#MPICC = mpicc

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
MPIFC = mpiifort
MPICC = mpiicc

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
#MPIFC = mpif90
#MPICC = mpicc

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
#LOPTS = -lm

# Makefile PGI compiler with Linux

# OpenMPI
#MPIFC = mpipgf90
#MPICC = mpicc

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

#MPIFC = ftn
#MPICC = cc

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

all : fpdpic2 cpdpic2

special: fpdpic2_c cpdpic2_f

# Version using Fortran77 pplib2.f
#fpdpic2 : fpdpic2.o fpdpush2.o fpplib2.o dtimer.o
#	$(MPIFC) $(OPTS90) $(LOPTS) -o fpdpic2 \
#       fpdpic2.o fpdpush2.o fpplib2.o pdpush2_h.o pplib2_h.o dtimer.o

# Version using Fortran90 pplib2.f90
fpdpic2 : fpdpic2.o fpdpush2.o dtimer.o
	$(MPIFC) $(OPTS90) $(LOPTS) -o fpdpic2 \
        fpdpic2.o fpdpush2.o f90pplib2.o pdpush2_h.o dtimer.o

cpdpic2 : cpdpic2.o cpdpush2.o cpplib2.o dtimer.o
	$(MPICC) $(CCOPTS) $(LOPTS) -o cpdpic2 \
        cpdpic2.o cpdpush2.o cpplib2.o dtimer.o

fpdpic2_c : fpdpic2_c.o cpdpush2.o cpplib2.o dtimer.o
	$(MPIFC) $(OPTS90) $(LOPTS) -o fpdpic2_c \
        fpdpic2_c.o cpdpush2.o cpplib2.o dtimer.o

cpdpic2_f : cpdpic2.o cpdpush2_f.o cpplib2_f.o fpdpush2.o fpplib2.o dtimer.o
	$(MPIFC) $(OPTS90) $(LOPTS) $(LEGACY) -o cpdpic2_f \
        cpdpic2.o cpdpush2_f.o cpplib2_f.o fpdpush2.o fpplib2.o dtimer.o

# Compilation rules

dtimer.o : dtimer.c
	$(CC) $(CCOPTS) -c dtimer.c

fpplib2.o : pplib2.f
	$(MPIFC) $(OPTS77) -o fpplib2.o -c pplib2.f

pplib2_h.o : pplib2_h.f90
	$(MPIFC) $(OPTS90) -o pplib2_h.o -c pplib2_h.f90

f90pplib2.o : pplib2.f90
	$(MPIFC) $(OPTS90) -o f90pplib2.o -c pplib2.f90

cpplib2.o : pplib2.c
	$(MPICC) $(CCOPTS) -o cpplib2.o -c pplib2.c

cpplib2_f.o : pplib2_f.c
	$(MPICC) $(CCOPTS) -o cpplib2_f.o -c pplib2_f.c

fpdpush2.o : pdpush2.f
	$(MPIFC) $(OPTS90) -o fpdpush2.o -c pdpush2.f

pdpush2_h.o : pdpush2_h.f90
	$(MPIFC) $(OPTS90) -o pdpush2_h.o -c pdpush2_h.f90

cpdpush2.o : pdpush2.c
	$(MPICC) $(CCOPTS) -o cpdpush2.o -c pdpush2.c

cpdpush2_f.o : pdpush2_f.c
	$(MPICC) $(CCOPTS) -o cpdpush2_f.o -c pdpush2_f.c

# Version using Fortran77 pplib2.f
#fpdpic2.o : pdpic2.f90 pdpush2_h.o pplib2_h.o
#	$(MPIFC) $(OPTS90) -o fpdpic2.o -c pdpic2.f90

# Version using Fortran90 pplib2.f90
fpdpic2.o : pdpic2.f90 f90pplib2.o pdpush2_h.o
	$(MPIFC) $(OPTS90) -o fpdpic2.o -c pdpic2.f90

cpdpic2.o : pdpic2.c
	$(MPICC) $(CCOPTS) -o cpdpic2.o -c pdpic2.c

fpdpic2_c.o : pdpic2_c.f90
	$(MPIFC) $(OPTS90) -o fpdpic2_c.o -c pdpic2_c.f90

clean :
	rm -f *.o *.mod

clobber: clean
	rm -f fpdpic2 cpdpic2 fpdpic2_c cpdpic2_f
