#Makefile for 3D Darwin OpenMP PIC codes

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

all : fmdpic3

special : cmdpic3_f

fmdpic3 : fmdpic3.o fmdpush3.o fomplib.o dtimer.o
	$(MPFC) $(OPTS90) -o fmdpic3 fmdpic3.o fmdpush3.o fomplib.o mdpush3_h.o \
        dtimer.o

cmdpic3 : cmdpic3.o cmdpush3.o complib.o dtimer.o
	$(MPCC) $(CCOPTS) -o cmdpic3 cmbdpic3.o mdpush3.o complib.o dtimer.o -lm

fmdpic3_c : fmdpic3_c.o cmdpush3.o complib.o dtimer.o
	$(MPFC) $(OPTS90) -o fmdpic3_c fmdpic3_c.o cmdpush3.o complib.o dtimer.o

cmdpic3_f : cmdpic3.o cmdpush3_f.o complib_f.o fmdpush3.o fomplib.o dtimer.o
	$(MPFC) $(CCOPTS) $(LEGACY) -o cmdpic3_f cmdpic3.o cmdpush3_f.o \
	    complib_f.o fmdpush3.o fomplib.o dtimer.o -lm

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

fmdpush3.o : mdpush3.f
	$(MPFC) $(OPTS90) -o fmdpush3.o -c mdpush3.f

mdpush3_h.o : mdpush3_h.f90
	$(FC90) $(OPTS90) -o mdpush3_h.o -c mdpush3_h.f90

cmdpush3.o : mdpush3.c
	$(MPCC) $(CCOPTS) -o cmdpush3.o -c mdpush3.c

fmdpic3.o : mdpic3.f90 mdpush3_h.o omplib_h.o
	$(FC90) $(OPTS90) -o fmdpic3.o -c mdpic3.f90

cmdpush3_f.o : mdpush3_f.c
	$(MPCC) $(CCOPTS) -o cmdpush3_f.o -c mdpush3_f.c

cmdpic3.o : mdpic3.c
	$(CC) $(CCOPTS) -o cmdpic3.o -c mdpic3.c

fmbpic3_c.o : mbpic3_c.f90
	$(FC90) $(OPTS90) -o fmdpic3_c.o -c mdpic3_c.f90

clean :
	rm -f *.o *.mod

clobber: clean
	rm -f fmdpic3 cmdpic3 fmdpic3_c cmdpic3_f
