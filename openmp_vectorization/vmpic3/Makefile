#Makefile for 3D Electrostatic OpenMP/Vector PIC codes

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

all : fvmpic3 f03vmpic3

special: cvmpic3_f

fvmpic3 : fvmpic3.o fvmpush3.o fomplib.o csselib3.o csseflib3.o dtimer.o
	$(MPFC) $(OPTS90) -o fvmpic3 fvmpic3.o fvmpush3.o fomplib.o csselib3.o \
	csseflib3.o sselib3_h.o sseflib3_h.o vmpush3_h.o omplib_h.o dtimer.o

cvmpic3 : cvmpic3.o cvmpush3.o complib.o csselib3.o dtimer.o
	$(MPCC) $(CCOPTS) -o cvmpic3 cvmpic3.o cvmpush3.o complib.o csselib3.o \
    dtimer.o -lm

f03vmpic3 : f03vmpic3.o fvmpush3.o fomplib.o csselib3.o dtimer.o
	$(MPFC) $(OPTS03) -o f03vmpic3 f03vmpic3.o fvmpush3.o fomplib.o csselib3.o \
	sselib3_c.o vmpush3_h.o omplib_h.o dtimer.o

fvmpic3_c : fvmpic3_c.o cvmpush3.o complib.o dtimer.o
	$(MPFC) $(OPTS90) -o fvmpic3_c fvmpic3_c.o cvmpush3.o complib.o dtimer.o

cvmpic3_f : cvmpic3.o cvmpush3_f.o complib_f.o fvmpush3.o fomplib.o csselib3.o \
            dtimer.o
	$(MPFC) $(OPTS90) $(LEGACY) -o cvmpic3_f cvmpic3.o cvmpush3_f.o complib_f.o \
	fvmpush3.o fomplib.o csselib3.o dtimer.o -lm

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

fvmpush3.o : vmpush3.f
	$(MPFC) $(OPTS90) -o fvmpush3.o -c vmpush3.f

vmpush3_h.o : vmpush3_h.f90
	$(FC90) $(OPTS90) -o vmpush3_h.o -c vmpush3_h.f90

cvmpush3.o : vmpush3.c
	$(MPCC) $(CCOPTS) -o cvmpush3.o -c vmpush3.c

csselib3.o : sselib3.c
	$(CC) $(CCOPTS) -o csselib3.o -c sselib3.c

csseflib3.o : sseflib3.c
	$(CC) $(CCOPTS) -o csseflib3.o -c sseflib3.c

sselib3_h.o : sselib3_h.f90
	$(FC90) $(OPTS90) -o sselib3_h.o -c sselib3_h.f90

sseflib3_h.o : sseflib3_h.f90 sselib3_h.o
	$(FC90) $(OPTS90) -o sseflib3_h.o -c sseflib3_h.f90

sselib3_c.o : sselib3_c.f03
	$(FC03) $(OPTS03) -o sselib3_c.o -c $(FF03) sselib3_c.f03

fvmpic3.o : vmpic3.f90 sseflib3_h.o vmpush3_h.o omplib_h.o
	$(FC90) $(OPTS90) -o fvmpic3.o -c vmpic3.f90

cvmpush3_f.o : vmpush3_f.c
	$(MPCC) $(CCOPTS) -o cvmpush3_f.o -c vmpush3_f.c

cvmpic3.o : vmpic3.c
	$(CC) $(CCOPTS) -o cvmpic3.o -c vmpic3.c

f03vmpic3.o : vmpic3.f03 sselib3_c.o vmpush3_h.o omplib_h.o
	$(FC03) $(OPTS03) -o f03vmpic3.o -c $(FF03) vmpic3.f03

clean :
	rm -f *.o *.mod

clobber: clean
	rm -f fvmpic3 cvmpic3 f03vmpic3 fvmpic3_c cvmpic3_f
