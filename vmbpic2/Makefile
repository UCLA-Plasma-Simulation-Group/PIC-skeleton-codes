#Makefile for 2-1/2D Electromagnetic OpenMP/Vector PIC codes

# Makefile gfortran compiler with MacOS X

#OpenMP
#MPFC = gfortran -fopenmp
#MPCC = gcc -fopenmp
#NoOpenMP
#MPFC = gfortran
#MPCC = gcc

#FC90 = gfortran
#FC03 = gfortran
#CC = gcc

#OPTS90 = -O3 -ftree-vectorize
#OPTS90 = -O3 -ftree-vectorizer-verbose=6
#OPTS90 = -O3 -fno-tree-vectorize
#OPTS90 = -O3 -fdefault-real-8 -fdefault-double-8
#OPTS90 = -O3 -fcheck=bounds -fdefault-real-8 -fdefault-double-8 -Wall -std=f95

#OPTS03 = -O3 -ftree-vectorize
#OPTS03 = -O3 -ftree-vectorizer-verbose=6
#OPTS03 = -O3 -fno-tree-vectorize
#OPTS03 = -O3 -fdefault-real-8 -fdefault-double-8
#OPTS03 = -O3 -fcheck=bounds -fdefault-real-8 -fdefault-double-8 -Wall -std=f2003

#CCOPTS = -O3 -Wall -std=c99
#CCOPTS = -O3 -fno-tree-vectorize -Wall -std=c99

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
#FC03 = ifort
#CC = gcc

#OPTS90 = -O3
#OPTS90 = -O3 -vec-report6
#OPTS90 = -O3 -no-vec
#OPTS90 = -O3 -r8
#OPTS90 = -O3 -CB -r8 -warn all -std90

#FF03 = -Tf
#OPTS03 = -O3 -free -std03
#OPTS03 = -O3 -free -vec-report6 -std03
#OPTS03 = -O3 -r8 -free -std03
#OPTS03 = -O3 -CB -r8 -warn all -free -std03

#CCOPTS = -O3 -Wall -std=c99
#CCOPTS = -O3 -ftree-vectorize -Wall -std=c99
#CCOPTS = -O3 -fno-tree-vectorize -Wall -std=c99
#LEGACY = -nofor_main

# Makefile Intel compiler with Linux

# OpenMP
MPFC = ifort -openmp
MPCC = icc -openmp
#NoOpenMP
#MPFC = ifort
#MPCC = icc

FC90 = ifort
FC03 = ifort
CC = icc

#OPTS90 = -O3
OPTS90 = -O3 -axSSE2
#OPTS90 = -O3 -axAVX
#OPTS90 = -O3 -vec-report6
#OPTS90 = -O3 -no-vec
#OPTS90 = -O3 -r8
#OPTS90 = -O3 -CB -r8 -warn all -std90

FF03 = -Tf
#OPTS03 = -O3 -free -std03
OPTS03 = -O3 -axSSE2 -free -std03
#OPTS03 = -O3 -axAVX -free -std03
#OPTS03 = -O3 -free -vec-report6 -std03
#OPTS03 = -O3 -free -no-vec -std03
#OPTS03 = -O3 -r8 -free -std03
#OPTS03 = -O3 -CB -r8 -warn all -free -std03

CCOPTS = -O3 -std=c99
#CCOPTS = -O3 -axSSE2 -std=c99
#CCOPTS = -O3 -axAVX -std=c99
#CCOPTS = -O3 -no-vec -Wall -std=c99
LEGACY = -nofor_main

# Makefile gfortran compiler with Linux

# OpenMP
#MPFC = gfortran -fopenmp
#MPCC = gcc -fopenmp
#NoOpenMP
#MPFC = gfortran
#MPCC = gcc

#FC90 = gfortran
#FC03 = gfortran
#CC = gcc

#OPTS90 = -O3 -ftree-vectorize
#OPTS90 = -O3 -ftree-vectorizer-verbose=6
#OPTS90 = -O3 -fno-tree-vectorize
#OPTS90 = -O3 -fdefault-real-8 -fdefault-double-8
#OPTS90 = -O3 -fbounds-check -fdefault-real-8 -fdefault-double-8 -Wall -std=f95

#OPTS03 = -O3 -ftree-vectorize
#OPTS03 = -O3 -ftree-vectorizer-verbose=6
#OPTS03 = -O3 -fno-tree-vectorize
#OPTS03 = -O3 -fdefault-real-8 -fdefault-double-8
#OPTS03 = -O3 -fbounds-check -fdefault-real-8 -fdefault-double-8 -Wall -std=f2003

#CCOPTS = -O3 -Wall -std=c99
#CCOPTS = -O3 -ftree-vectorize -Wall -std=c99
#CCOPTS = -O3 -fno-tree-vectorize -Wall -std=c99

#LOPTS = -lgomp
#LEGACY =

# Makefile PGI compiler with Linux

# OpenMP
#MPFC = pgf90
#MPCC = gcc -fopenmp

#FC90 = pgf90
#FC03 = pgf90
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

all : fvmbpic2 cvmbpic2 f03vmbpic2

fvmbpic2 : fvmbpic2.o fvmbpush2.o fomplib.o csselib2.o csseflib2.o cssembpush2.o \
           dtimer.o
	$(MPFC) $(OPTS90) -o fvmbpic2 fvmbpic2.o fvmbpush2.o fomplib.o csselib2.o \
    csseflib2.o cssembpush2.o sselib2_h.o sseflib2_h.o ssembpush2_h.o vmbpush2_h.o \
    omplib_h.o dtimer.o

cvmbpic2 : cvmbpic2.o cvmbpush2.o complib.o csselib2.o cssembpush2.o dtimer.o
	$(MPCC) $(CCOPTS) -o cvmbpic2 cvmbpic2.o cvmbpush2.o complib.o csselib2.o \
    cssembpush2.o dtimer.o -lm

f03vmbpic2 : f03vmbpic2.o fvmbpush2.o fomplib.o csselib2.o cssembpush2.o dtimer.o
	$(MPFC) $(OPTS03) -o f03vmbpic2 f03vmbpic2.o fvmbpush2.o fomplib.o csselib2.o \
    cssembpush2.o sselib2_c.o ssembpush2_c.o vmbpush2_h.o omplib_h.o dtimer.o

# Compilation rules

dtimer.o : dtimer.c
	$(CC) $(CCOPTS) -c dtimer.c

fomplib.o : omplib.f
	$(MPFC) $(OPTS90) -o fomplib.o -c omplib.f

omplib_h.o : omplib_h.f90
	$(FC90) $(OPTS90) -o omplib_h.o -c omplib_h.f90

complib.o : omplib.c
	$(MPCC) $(CCOPTS) -o complib.o -c omplib.c

complib_f.o : omplib_f.c
	$(MPCC) $(CCOPTS) -o complib_f.o -c omplib_f.c

fvmbpush2.o : vmbpush2.f
	$(MPFC) $(OPTS90) -o fvmbpush2.o -c vmbpush2.f

vmbpush2_h.o : vmbpush2_h.f90
	$(FC90) $(OPTS90) -o vmbpush2_h.o -c vmbpush2_h.f90

cvmbpush2.o : vmbpush2.c
	$(MPCC) $(CCOPTS) -o cvmbpush2.o -c vmbpush2.c

csselib2.o : sselib2.c
	$(CC) $(CCOPTS) -o csselib2.o -c sselib2.c

csseflib2.o : sseflib2.c
	$(CC) $(CCOPTS) -o csseflib2.o -c sseflib2.c

cssembpush2.o : ssembpush2.c
	$(MPCC) $(CCOPTS) -o cssembpush2.o -c ssembpush2.c

sselib2_h.o : sselib2_h.f90
	$(FC90) $(OPTS90) -o sselib2_h.o -c sselib2_h.f90

sseflib2_h.o : sseflib2_h.f90 sselib2_h.o
	$(FC90) $(OPTS90) -o sseflib2_h.o -c sseflib2_h.f90

ssembpush2_h.o : ssembpush2_h.f90
	$(FC90) $(OPTS90) -o ssembpush2_h.o -c ssembpush2_h.f90

sselib2_c.o : sselib2_c.f03
	$(FC03) $(OPTS03) -o sselib2_c.o -c $(FF03) sselib2_c.f03

ssembpush2_c.o : ssembpush2_c.f03
	$(FC03) $(OPTS03) -o ssembpush2_c.o -c $(FF03) ssembpush2_c.f03

fvmbpic2.o : vmbpic2.f90 sseflib2_h.o ssembpush2_h.o vmbpush2_h.o omplib_h.o
	$(FC90) $(OPTS90) -o fvmbpic2.o -c vmbpic2.f90

cvmbpic2.o : vmbpic2.c
	$(CC) $(CCOPTS) -o cvmbpic2.o -c vmbpic2.c

f03vmbpic2.o : vmbpic2.f03 sselib2_c.o ssembpush2_c.o vmbpush2_h.o omplib_h.o
	$(FC03) $(OPTS03) -o f03vmbpic2.o -c $(FF03) vmbpic2.f03

clean :
	rm -f *.o *.mod

clobber: clean
	rm -f fvmbpic2 cvmbpic2 f03vmbpic2
