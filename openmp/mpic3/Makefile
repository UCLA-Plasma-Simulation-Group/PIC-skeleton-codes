#Makefile for 3D Electrostatic OpenMP PIC codes

# Makefile gfortran compiler with MacOS X

#OpenMP
#MPFC = gfortran -fopenmp
#MPCC = gcc -fopenmp
#NoOpenMP
#MPFC = gfortran
#MPCC = gcc

#FC90 = gfortran
#CC = gcc

#OPTS90 = -O3
#OPTS90 = -O3 -fdefault-real-8 -fdefault-double-8
#OPTS90 = -O3 -fcheck=bounds -fdefault-real-8 -Wall -std=f95

#CCOPTS = -O3 -Wall -std=c99

#LOPTS = -lgomp
#LEGACY =

# Makefile Intel compiler with Mac OS X

# OpenMP
#MPFC = ifort -openmp
#MPCC = gcc -fopenmp
#NoOpenMP
#MPFC = ifort
#MPCC = gcc

#FC90 = ifort
#CC = icc

#OPTS90 = -O3
#OPTS90 = -O3 -r8
#OPTS90 = -O3 -CB -r8 -warn all -std90

#CCOPTS = -O3 -std=c99
#LEGACY = -nofor_main

# Makefile Intel compiler with Linux

# OpenMP
#MPFC = ifort -openmp
#MPCC = icc -openmp
#NoOpenMP
#MPFC = ifort
#MPCC = icc

#FC90 = ifort
#CC = icc

#OPTS90 = -O3
#OPTS90 = -O3 -r8
#OPTS90 = -O3 -CB -r8 -warn all -std90

#CCOPTS = -O3 -std=c99
#CCOPTS = -O3 -no-vec -Wall -std=c99
#LEGACY = -nofor_main

# Makefile gfortran compiler with Linux

# OpenMP
MPFC = gfortran -fopenmp
MPCC = gcc -fopenmp
#NoOpenMP
#MPFC = gfortran
#MPCC = gcc

FC90 = gfortran
CC = gcc

OPTS90 = -O3
#OPTS90 = -O3 -fdefault-real-8 -fdefault-double-8
#OPTS90 = -O3 -fbounds-check -fdefault-real-8 -Wall -std=f95

CCOPTS = -O3 -Wall -std=c99
#CCOPTS = -O3 -fno-tree-vectorize -Wall -std=c99
LEGACY =

# Makefile PGI compiler with Linux

# OpenMP
#MPFC = pgf90
#MPCC = gcc -fopenmp

#FC90 = pgf90
#CC = gcc

#OPTS90 = -O3
#OPTS90 = -O3 -r8
#OPTS90 = -O3 -Mbounds -r8 -Mstandard
#LEGACY = -Mnomain

# Makefile Cray compiler with Linux

# OpenMP
#MPFC = ftn -h omp
#MPCC = cc -h omp
#NoOpenMP
#MPFC = ftn -h noomp
#MPCC = cc -h noomp

#FC90 = ftn
#CC = cc

#OPTS90 = -O 3
#OPTS90 = -O 3 -s real64
#OPTS90 = -O 3 -R b -s real64 -en

#CCOPTS = -O 3 -h c99 -h conform
#LEGACY =

#

# Linkage rules

all : fmpic3 cmpic3

special: fmpic3_c cmpic3_f

fmpic3 : fmpic3.o fmpush3.o fomplib.o dtimer.o
	$(MPFC) $(OPTS90) -o fmpic3 fmpic3.o fmpush3.o fomplib.o mpush3_h.o \
        omplib_h.o dtimer.o

cmpic3 : cmpic3.o cmpush3.o complib.o dtimer.o
	$(MPCC) $(CCOPTS) -o cmpic3 cmpic3.o cmpush3.o complib.o dtimer.o -lm

fmpic3_c : fmpic3_c.o cmpush3.o complib.o dtimer.o
	$(MPFC) $(OPTS90) -o fmpic3_c fmpic3_c.o cmpush3.o complib.o dtimer.o 

cmpic3_f : cmpic3.o cmpush3_f.o complib_f.o fmpush3.o fomplib.o dtimer.o
	$(MPFC) $(CCOPTS) $(LEGACY) -o cmpic3_f cmpic3.o cmpush3_f.o complib_f.o \
	fmpush3.o fomplib.o dtimer.o -lm

# Compilation rules

dtimer.o : dtimer.c
	$(CC) $(CCOPTS) -c dtimer.c

#OPENMP
fomplib.o : omplib.f
	$(MPFC) $(OPTS90) -o fomplib.o -c omplib.f

omplib_h.o : omplib_h.f90
	$(FC90) $(OPTS90) -o omplib_h.o -c omplib_h.f90

complib.o : omplib.c
	$(MPCC) $(CCOPTS) -o complib.o -c omplib.c

complib_f.o : omplib_f.c
	$(MPCC) $(CCOPTS) -o complib_f.o -c omplib_f.c

fmpush3.o : mpush3.f
	$(MPFC) $(OPTS90) -o fmpush3.o -c mpush3.f

mpush3_h.o : mpush3_h.f90
	$(FC90) $(OPTS90) -o mpush3_h.o -c mpush3_h.f90

cmpush3.o : mpush3.c
	$(MPCC) $(CCOPTS) -o cmpush3.o -c mpush3.c

fmpic3.o : mpic3.f90 mpush3_h.o omplib_h.o
	$(FC90) $(OPTS90) -o fmpic3.o -c mpic3.f90

cmpush3_f.o : mpush3_f.c
	$(MPCC) $(CCOPTS) -o cmpush3_f.o -c mpush3_f.c

cmpic3.o : mpic3.c
	$(CC) $(CCOPTS) -o cmpic3.o -c mpic3.c

fmpic3_c.o : mpic3_c.f90
	$(FC90) $(OPTS90) -o fmpic3_c.o -c mpic3_c.f90

clean :
	rm -f *.o *.mod

clobber: clean
	rm -f fmpic3 cmpic3 fmpic3_c cmpic3_f
