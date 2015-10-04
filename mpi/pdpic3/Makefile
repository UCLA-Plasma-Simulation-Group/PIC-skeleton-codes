#Makefile for 3D Darwin MPI PIC codes

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
#LOPTS =

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
#MPIFC = mpiifort
#MPICC = mpiicc

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
MPIFC = mpif90
MPICC = mpicc

CC = gcc

OPTS90 = -O3
#OPTS90 = -O3 -fdefault-real-8 -fdefault-double-8
#OPTS90 = -O3 -fbounds-check -fdefault-real-8 -fdefault-double-8 -Wall -std=f95

OPTS77 = -O3
#OPTS77 = -O3 -fdefault-real-8 -fdefault-double-8
#OPTS77 = -O3 -fbounds-check -fdefault-real-8 -fdefault-double-8 -Wall -std=f95

CCOPTS = -O3 -std=c99
#CCOPTS = -O3 -Wall -std=c99
#LEGACY =
LOPTS =

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

all : fpdpic3

special: cpdpic3_f

# Version using Fortran77 pplib3.f
#fpdpic3 : fpdpic3.o fpdpush3.o fpplib3.o dtimer.o
#	$(MPIFC) $(OPTS90) $(LOPTS) -o fpdpic3 \
#       fpdpic3.o fpdpush3.o fpplib3.o pdpush3_h.o pplib3_h.o dtimer.o -lm

# Version using Fortran90 pplib3.f90
fpdpic3 : fpdpic3.o fpdpush3.o dtimer.o
	$(MPIFC) $(OPTS90) $(LOPTS) -o fpdpic3 \
        fpdpic3.o fpdpush3.o f90pplib3.o pdpush3_h.o dtimer.o

cpdpic3_f : cpdpic3.o cpdpush3_f.o cpplib3_f.o fpdpush3.o fpplib3.o dtimer.o
	$(MPIFC) $(OPTS90) $(LOPTS) $(LEGACY) -o cpdpic3_f \
        cpdpic3.o cpdpush3_f.o cpplib3_f.o fpdpush3.o fpplib3.o dtimer.o -lm

# Compilation rules

dtimer.o : dtimer.c
	$(CC) $(CCOPTS) -c dtimer.c

fpplib3.o : pplib3.f
	$(MPIFC) $(OPTS77) -o fpplib3.o -c pplib3.f

pplib3_h.o : pplib3_h.f90
	$(MPIFC) $(OPTS90) -o pplib3_h.o -c pplib3_h.f90

f90pplib3.o : pplib3.f90
	$(MPIFC) $(OPTS90) -o f90pplib3.o -c pplib3.f90

cpplib3_f.o : pplib3_f.c
	$(MPICC) $(CCOPTS) -o cpplib3_f.o -c pplib3_f.c

fpdpush3.o : pdpush3.f
	$(MPIFC) $(OPTS90) -o fpdpush3.o -c pdpush3.f

pdpush3_h.o : pdpush3_h.f90
	$(MPIFC) $(OPTS90) -o pdpush3_h.o -c pdpush3_h.f90

cpdpush3_f.o : pdpush3_f.c
	$(MPICC) $(CCOPTS) -o cpdpush3_f.o -c pdpush3_f.c

# Version using Fortran77 pplib3.f
#fpdpic3.o : pdpic3.f90 pdpush3_h.o pplib3_h.o
#	$(MPIFC) $(OPTS90) -o fpdpic3.o -c pdpic3.f90

# Version using Fortran90 pplib3.f90
fpdpic3.o : pdpic3.f90 f90pplib3.o pdpush3_h.o
	$(MPIFC) $(OPTS90) -o fpdpic3.o -c pdpic3.f90

cpdpic3.o : pdpic3.c
	$(MPICC) $(CCOPTS) -o cpdpic3.o -c pdpic3.c

clean :
	rm -f *.o *.mod

clobber: clean
	rm -f fpdpic3 cpdpic3_f
