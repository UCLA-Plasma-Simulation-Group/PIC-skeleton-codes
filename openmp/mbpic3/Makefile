#Makefile for 3D Electromagnetic OpenMP PIC codes

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

all : fmbpic3 cmbpic3

special : fmbpic3_c cmbpic3_f

fmbpic3 : fmbpic3.o fmbpush3.o fomplib.o dtimer.o
	$(MPFC) $(OPTS90) -o fmbpic3 fmbpic3.o fmbpush3.o fomplib.o mbpush3_h.o \
        dtimer.o

cmbpic3 : cmbpic3.o cmbpush3.o complib.o dtimer.o
	$(MPCC) $(CCOPTS) -o cmbpic3 cmbpic3.o cmbpush3.o complib.o dtimer.o -lm

fmbpic3_c : fmbpic3_c.o cmbpush3.o complib.o dtimer.o
	$(MPFC) $(OPTS90) -o fmbpic3_c fmbpic3_c.o cmbpush3.o complib.o dtimer.o

cmbpic3_f : cmbpic3.o cmbpush3_f.o complib_f.o fmbpush3.o fomplib.o dtimer.o
	$(MPFC) $(CCOPTS) $(LEGACY) -o cmbpic3_f cmbpic3.o cmbpush3_f.o \
	    complib_f.o fmbpush3.o fomplib.o dtimer.o -lm

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

fmbpush3.o : mbpush3.f
	$(MPFC) $(OPTS90) -o fmbpush3.o -c mbpush3.f

mbpush3_h.o : mbpush3_h.f90
	$(FC90) $(OPTS90) -o mbpush3_h.o -c mbpush3_h.f90

cmbpush3.o : mbpush3.c
	$(MPCC) $(CCOPTS) -o cmbpush3.o -c mbpush3.c

fmbpic3.o : mbpic3.f90 mbpush3_h.o omplib_h.o
	$(FC90) $(OPTS90) -o fmbpic3.o -c mbpic3.f90

cmbpush3_f.o : mbpush3_f.c
	$(MPCC) $(CCOPTS) -o cmbpush3_f.o -c mbpush3_f.c

cmbpic3.o : mbpic3.c
	$(CC) $(CCOPTS) -o cmbpic3.o -c mbpic3.c

fmbpic3_c.o : mbpic3_c.f90
	$(FC90) $(OPTS90) -o fmbpic3_c.o -c mbpic3_c.f90

clean :
	rm -f *.o *.mod

clobber: clean
	rm -f fmbpic3 cmbpic3 fmbpic3_c cmbpic3_f
