#Makefile for 2D Electrostatic GPU Emulation OpenMP/AVX PIC codes

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
#CC = gcc

#OPTS90 = -O3
#OPTS90 = -O3 -r8
#OPTS90 = -O3 -CB -r8 -warn all -std90

#CCOPTS = -O3 -Wall -std=c99
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
#OPTS90 = -O3 -no-vec
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

LOPTS = -lgomp
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

#CCOPTS = -O3 -Wall -std=c99
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

all : fmpic2 cmpic2

special: fmpic2_c cmpic2_f

fmpic2 : fmpic2.o fmpush2.o fomplib.o dtimer.o
	$(MPFC) $(OPTS90) -o fmpic2 fmpic2.o fmpush2.o fomplib.o mpush2_h.o \
    omplib_h.o dtimer.o

cmpic2 : cmpic2.o cmpush2.o complib.o dtimer.o
	$(MPCC) $(CCOPTS) -o cmpic2 cmpic2.o cmpush2.o complib.o \
    dtimer.o -lm

fmpic2_c : fmpic2_c.o cmpush2.o complib.o dtimer.o
	$(MPFC) $(OPTS90) -o fmpic2_c fmpic2_c.o cmpush2.o complib.o \
    dtimer.o

cmpic2_f : cmpic2.o cmpush2_f.o complib_f.o fmpush2.o fomplib.o dtimer.o 
	$(MPFC) $(OPTS90) $(LEGACY) -o cmpic2_f cmpic2.o cmpush2_f.o complib_f.o \
    fmpush2.o fomplib.o dtimer.o -lm

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

fmpush2.o : mpush2.f
	$(MPFC) $(OPTS90) -o fmpush2.o -c mpush2.f

mpush2_h.o : mpush2_h.f90
	$(FC90) $(OPTS90) -o mpush2_h.o -c mpush2_h.f90

cmpush2.o : mpush2.c
	$(MPCC) $(CCOPTS) -o cmpush2.o -c mpush2.c

fmpic2.o : mpic2.f90 mpush2_h.o omplib_h.o
	$(FC90) $(OPTS90) -o fmpic2.o -c mpic2.f90

cmpush2_f.o : mpush2_f.c
	$(MPCC) $(CCOPTS) -o cmpush2_f.o -c mpush2_f.c

cmpic2.o : mpic2.c
	$(CC) $(CCOPTS) -o cmpic2.o -c mpic2.c

fmpic2_c.o : mpic2_c.f90
	$(FC90) $(OPTS90) -o fmpic2_c.o -c mpic2_c.f90

clean :
	rm -f *.o *.mod

clobber: clean
	rm -f fmpic2 cmpic2 fmpic2_c cmpic2_f
