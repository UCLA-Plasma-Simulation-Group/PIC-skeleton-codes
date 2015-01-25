#Makefile for 2D Electrostatic MPI PIC codes

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
#OPTS77 = -O3 -fcheck=bounds -fdefault-real-8 -Wall

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
#LOPTS = -lm

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
LOPTS = -lm

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

all : fppic2 cppic2

special: fppic2_c cppic2_f

# Version using Fortran77 pplib2.f
#fppic2 : fppic2.o fppush2.o fpplib2.o dtimer.o
#	$(MPIFC) $(OPTS90) $(LOPTS) -o fppic2 \
#       fppic2.o fppush2.o fpplib2.o ppush2_h.o pplib2_h.o dtimer.o

# Version using Fortran90 pplib2.f90
fppic2 : fppic2.o fppush2.o fpplib2.o dtimer.o
	$(MPIFC) $(OPTS90) $(LOPTS) -o fppic2 \
        fppic2.o fppush2.o f90pplib2.o ppush2_h.o dtimer.o

cppic2 : cppic2.o cppush2.o cpplib2.o dtimer.o
	$(MPICC) $(CCOPTS) $(LOPTS) -o cppic2 \
        cppic2.o cppush2.o cpplib2.o dtimer.o

fppic2_c : fppic2_c.o cppush2.o cpplib2.o dtimer.o
	$(MPIFC) $(OPTS90) $(LOPTS) -o fppic2_c \
        fppic2_c.o cppush2.o cpplib2.o dtimer.o

cppic2_f : cppic2.o cppush2_f.o cpplib2_f.o fppush2.o fpplib2.o dtimer.o
	$(MPIFC) $(OPTS90) $(LOPTS) $(LEGACY) -o cppic2_f \
        cppic2.o cppush2_f.o cpplib2_f.o fppush2.o fpplib2.o dtimer.o

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

fppush2.o : ppush2.f
	$(MPIFC) $(OPTS90) -o fppush2.o -c ppush2.f

ppush2_h.o : ppush2_h.f90
	$(MPIFC) $(OPTS90) -o ppush2_h.o -c ppush2_h.f90

cppush2.o : ppush2.c
	$(MPICC) $(CCOPTS) -o cppush2.o -c ppush2.c

# Version using Fortran77 pplib2.f
#fppic2.o : ppic2.f90 ppush2_h.o pplib2_h.o
#	$(MPIFC) $(OPTS90) -o fppic2.o -c ppic2.f90

# Version using Fortran90 pplib2.f90
fppic2.o : ppic2.f90 f90pplib2.o ppush2_h.o
	$(MPIFC) $(OPTS90) -o fppic2.o -c ppic2.f90

cppush2_f.o : ppush2_f.c
	$(MPICC) $(CCOPTS) -o cppush2_f.o -c ppush2_f.c

cppic2.o : ppic2.c
	$(MPICC) $(CCOPTS) -o cppic2.o -c ppic2.c

fppic2_c.o : ppic2_c.f90
	$(MPIFC) $(OPTS90) -o fppic2_c.o -c ppic2_c.f90

clean :
	rm -f *.o *.mod

clobber: clean
	rm -f fppic2 cppic2 fppic2_c cppic2_f
