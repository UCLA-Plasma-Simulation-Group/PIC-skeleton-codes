#Makefile for 1-2/2D Electromagnetic OpenMP PIC codes

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

#LOPTS = -lgomp#Makefile for 1D Electrostatic OpenMP PIC codes

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

all : fmbpic1 cmbpic1_f

special: cmbpic1_f

fmbpic1 : fmbpic1.o fmbpush1.o fomplib.o dtimer.o
	$(MPFC) $(OPTS90) -o fmbpic1 fmbpic1.o fmbpush1.o fomplib.o mbpush1_h.o \
    omplib_h.o dtimer.o

cmbpic1_f : cmbpic1.o cmbpush1_f.o complib_f.o fmbpush1.o fomplib.o dtimer.o 
	$(MPFC) $(OPTS90) $(LEGACY) -o cmbpic1_f cmbpic1.o cmbpush1_f.o complib_f.o \
    fmbpush1.o fomplib.o dtimer.o -lm

# Compilation rules

dtimer.o : dtimer.c
	$(CC) $(CCOPTS) -c dtimer.c

#OPENMP
fomplib.o : omplib.f
	$(MPFC) $(OPTS90) -o fomplib.o -c omplib.f

omplib_h.o : omplib_h.f90
	$(FC90) $(OPTS90) -o omplib_h.o -c omplib_h.f90

complib_f.o : omplib_f.c
	$(MPCC) $(CCOPTS) -o complib_f.o -c omplib_f.c

fmbpush1.o : mbpush1.f
	$(MPFC) $(OPTS90) -o fmbpush1.o -c mbpush1.f

mbpush1_h.o : mbpush1_h.f90
	$(FC90) $(OPTS90) -o mbpush1_h.o -c mbpush1_h.f90

cmbpush1_f.o : mbpush1_f.c
	$(MPCC) $(CCOPTS) -o cmbpush1_f.o -c mbpush1_f.c

fmbpic1.o : mbpic1.f90 mbpush1_h.o  omplib_h.o
	$(FC90) $(OPTS90) -o fmbpic1.o -c mbpic1.f90

cmbpic1.o : mbpic1.c
	$(CC) $(CCOPTS) -o cmbpic1.o -c mbpic1.c

clean :
	rm -f *.o *.mod

clobber: clean
	rm -f fmbpic1 cmbpic1_f
