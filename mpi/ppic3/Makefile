#Makefile for 3D Electrostatic MPI PIC codes

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
#OPTS77 = -O3 -fcheck=bounds -fdefault-real-8 -fdefault-double-8 -Wall

#CCOPTS = -O3 -std=c99
#CCOPTS = -O3 -Wall -std=c99
#LEGACY =

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
#OPTS77 = -O3 -fbounds-check -fdefault-real-8 -fdefault-double-8 -Wall

CCOPTS = -O3 -std=c99
#CCOPTS = -O3 -Wall -std=c99
#LEGACY =

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

# MPI
#MPIFC = ftn
#MPICC = cc

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

all : fppic3

special: cppic3_f

# Version using Fortran77 pplib3.f
#fppic3 : fppic3.o fppush3.o fpplib3.o dtimer.o
#	$(MPIFC) $(OPTS90) $(LOPTS) -o fppic3 \
#       fppic3.o fppush3.o fpplib3.o ppush3_h.o pplib3_h.o dtimer.o -lm

# Version using Fortran90 pplib3.f90
fppic3 : fppic3.o fppush3.o f90pplib3.o dtimer.o
	$(MPIFC) $(OPTS90) $(LOPTS) -o fppic3 \
        fppic3.o fppush3.o f90pplib3.o ppush3_h.o dtimer.o

cppic3_f : cppic3.o cppush3_f.o cpplib3_f.o fppush3.o fpplib3.o dtimer.o
	$(MPIFC) $(OPTS90) $(LOPTS) $(LEGACY) -o cppic3_f \
        cppic3.o cppush3_f.o cpplib3_f.o fppush3.o fpplib3.o dtimer.o -lm

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

fppush3.o : ppush3.f
	$(MPIFC) $(OPTS90) -o fppush3.o -c ppush3.f

ppush3_h.o : ppush3_h.f90
	$(MPIFC) $(OPTS90) -o ppush3_h.o -c ppush3_h.f90

cppush3_f.o : ppush3_f.c
	$(MPICC) $(CCOPTS) -o cppush3_f.o -c ppush3_f.c

# Version using Fortran77 pplib3.f
#fppic3.o : ppic3.f90 ppush3_h.o pplib3_h.o
#	$(MPIFC) $(OPTS90) -o fppic3.o -c ppic3.f90

# Version using Fortran90 pplib3.f90
fppic3.o : ppic3.f90 f90pplib3.o ppush3_h.o
	$(MPIFC) $(OPTS90) -o fppic3.o -c ppic3.f90

cppic3.o : ppic3.c
	$(MPICC) $(CCOPTS) -o cppic3.o -c ppic3.c

clean :
	rm -f *.o *.mod

clobber: clean
	rm -f fppic3 cppic3_f
