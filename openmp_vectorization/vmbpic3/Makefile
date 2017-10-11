#Makefile for 3D Electromagnetic OpenMP/Vector PIC codes

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
#OPTS90 = -O3 -axAVX -vec-report6
#OPTS90 = -O3 -no-vec
#OPTS90 = -O3 -r8
#OPTS90 = -O3 -CB -r8 -warn all -std90

#FF03 = -Tf
#OPTS03 = -O3 -free -std03
#OPTS03 = -O3 -free -vec-report6 -std03
#OPTS03 = -O3 -free -axAVX -vec-report6 -std03
#OPTS03 = -O3 -r8 -free -std03
#OPTS03 = -O3 -CB -r8 -warn all -free -std03

#CCOPTS = -O3 -Wall -std=c99
#CCOPTS = -O3 -ftree-vectorize -Wall -std=c99
#CCOPTS = -O3 -fno-tree-vectorize -Wall -std=c99
#LEGACY = -nofor_main

# Makefile Intel compiler with Linux

# OpenMP
MPFC = ifort -qopenmp
MPCC = icc -qopenmp
#NoOpenMP
#MPFC = ifort
#MPCC = icc

FC90 = ifort
FC03 = ifort
CC = icc

#OPTS90 = -O3
OPTS90 = -O3 -qopt-report=5 -qopt-report-phase=vec
#OPTS90 = -O3 -axAVX -qopt-report=5 -qopt-report-phase=vec
#OPTS90 = -O3 -mmic -qopt-report=5 -qopt-report-phase=vec
#OPTS90 = -O3 -no-vec
#OPTS90 = -O3 -no-vec -no-simd
#OPTS90 = -O3 -r8
#OPTS90 = -O3 -CB -r8 -warn all -std90

FF03 = -Tf
#OPTS03 = -O3 -free -std03
OPTS03 = -O3 -free -qopt-report=5 -qopt-report-phase=vec -std03
#OPTS03 = -O3 -axAVX -qopt-report=5 -qopt-report-phase=vec -free -std03
#OPTS03 = -O3 -mmic -qopt-report=5 -qopt-report-phase=vec -free -std03
#OPTS03 = -O3 -free -no-vec -std03
#OPTS03 = -O3 -r8 -free -std03
#OPTS03 = -O3 -CB -r8 -warn all -free -std03

CCOPTS = -O3 -std=c99
#CCOPTS = -O3 -axAVX -std=c99
#CCOPTS = -O3 -mmic -std=c99
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

all : fvmbpic3 cvmbpic3 f03vmbpic3

special: fvmbpic3_c cvmbpic3_f

fvmbpic3 : fvmbpic3.o fvmbpush3.o fomplib.o cavx512lib3.o cavx512flib3.o \
           ckncmbpush3.o dtimer.o
	$(MPFC) $(OPTS90) -o fvmbpic3 fvmbpic3.o fvmbpush3.o fomplib.o cavx512flib3.o \
 	cavx512lib3.o ckncmbpush3.o avx512lib3_h.o avx512flib3_h.o kncmbpush3_h.o \
 	vmbpush3_h.o omplib_h.o dtimer.o

cvmbpic3 : cvmbpic3.o cvmbpush3.o complib.o cavx512lib3.o ckncmbpush3.o dtimer.o
	$(MPCC) $(CCOPTS) -o cvmbpic3 cvmbpic3.o cvmbpush3.o complib.o cavx512lib3.o \
    ckncmbpush3.o dtimer.o -lm

f03vmbpic3 : f03vmbpic3.o fvmbpush3.o fomplib.o cavx512lib3.o ckncmbpush3.o dtimer.o
	$(MPFC) $(OPTS03) -o f03vmbpic3 f03vmbpic3.o fvmbpush3.o fomplib.o cavx512lib3.o \
	avx512lib3_c.o kncmbpush3_c.o vmbpush3_h.o omplib_h.o ckncmbpush3.o dtimer.o

fvmbpic3_c : fvmbpic3_c.o cvmbpush3.o complib.o cavx512flib3.o cavx512lib3.o dtimer.o
	$(MPFC) $(OPTS90) -o fvmbpic3_c fvmbpic3_c.o cvmbpush3.o complib.o \
 	cavx512flib3.o cavx512lib3.o avx512flib3_h.o dtimer.o

cvmbpic3_f : cvmbpic3.o cvmbpush3_f.o complib_f.o fvmbpush3.o fomplib.o \
             cavx512lib3.o dtimer.o
	$(MPFC) $(OPTS90) $(LEGACY) -o cvmbpic3_f cvmbpic3.o cvmbpush3_f.o complib_f.o \
	fvmbpush3.o fomplib.o cavx512lib3.o dtimer.o -lm

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

fvmbpush3.o : vmbpush3.f
	$(MPFC) $(OPTS90) -o fvmbpush3.o -c vmbpush3.f

vmbpush3_h.o : vmbpush3_h.f90
	$(FC90) $(OPTS90) -o vmbpush3_h.o -c vmbpush3_h.f90

cvmbpush3.o : vmbpush3.c
	$(MPCC) $(CCOPTS) -o cvmbpush3.o -c vmbpush3.c

cavx512lib3.o : avx512lib3.c
	$(CC) $(CCOPTS) -o cavx512lib3.o -c avx512lib3.c

cavx512flib3.o : avx512flib3.c
	$(CC) $(CCOPTS) -o cavx512flib3.o -c avx512flib3.c

ckncmbpush3.o : kncmbpush3.c
	$(MPCC) $(CCOPTS) -o ckncmbpush3.o -c kncmbpush3.c

avx512lib3_h.o : avx512lib3_h.f90
	$(FC90) $(OPTS90) -o avx512lib3_h.o -c avx512lib3_h.f90

avx512flib3_h.o : avx512flib3_h.f90 avx512lib3_h.o
	$(FC90) $(OPTS90) -o avx512flib3_h.o -c avx512flib3_h.f90

kncmbpush3_h.o : kncmbpush3_h.f90
	$(FC90) $(OPTS90) -o kncmbpush3_h.o -c kncmbpush3_h.f90

avx512lib3_c.o : avx512lib3_c.f03
	$(FC03) $(OPTS03) -o avx512lib3_c.o -c $(FF03) avx512lib3_c.f03

kncmbpush3_c.o : kncmbpush3_c.f03
	$(FC03) $(OPTS03) -o kncmbpush3_c.o -c $(FF03) kncmbpush3_c.f03

fvmbpic3.o : vmbpic3.f90 avx512flib3_h.o kncmbpush3_h.o vmbpush3_h.o omplib_h.o
	$(FC90) $(OPTS90) -o fvmbpic3.o -c vmbpic3.f90

cvmbpush3_f.o : vmbpush3_f.c
	$(MPCC) $(CCOPTS) -o cvmbpush3_f.o -c vmbpush3_f.c

cvmbpic3.o : vmbpic3.c
	$(CC) $(CCOPTS) -o cvmbpic3.o -c vmbpic3.c

f03vmbpic3.o : vmbpic3.f03 avx512lib3_c.o kncmbpush3_c.o vmbpush3_h.o omplib_h.o
	$(FC03) $(OPTS03) -o f03vmbpic3.o -c $(FF03) vmbpic3.f03

fvmbpic3_c.o : vmbpic3_c.f90 avx512flib3_h.o
	$(FC90) $(OPTS90) -o fvmbpic3_c.o -c vmbpic3_c.f90

clean :
	rm -f *.o *.mod

clobber: clean
	rm -f fvmbpic3 cvmbpic3 f03vmbpic3 fvmbpic3_c cvmbpic3_f
