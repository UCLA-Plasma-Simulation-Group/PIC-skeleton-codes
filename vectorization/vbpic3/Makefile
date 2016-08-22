#Makefile for 3D Electromagnetic Vector PIC codes

# Makefile gfortran compiler with MacOS X

FC90 = gfortran
FC03 = gfortran
CC = gcc

OPTS90 = -O3 -ftree-vectorize
#OPTS90 = -O3 -ftree-vectorizer-verbose=6
#OPTS90 = -O3 -fno-tree-vectorize
#OPTS90 = -O3 -fdefault-real-8 -fdefault-double-8
#OPTS90 = -O3 -fcheck=bounds -fdefault-real-8 -Wall -std=f95

OPTS03 = -O3 -ftree-vectorize
#OPTS03 = -O3 -ftree-vectorizer-verbose=6
#OPTS03 = -O3 -fno-tree-vectorize
#OPTS03 = -O3 -fdefault-real-8 -fdefault-double-8
#OPTS03 = -O3 -fcheck=bounds -fdefault-real-8 -fdefault-double-8 -Wall -std=f2003

CCOPTS = -O3 -Wall -std=c99
#CCOPTS = -O3 -fno-tree-vectorize -Wall -std=c99
LEGACY =

# Makefile Intel compiler with Mac OS X

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

#CCOPTS = -O3 -std=c99
#CCOPTS = -O3 -fno-tree-vectorize -Wall -std=c99
#LEGACY = -nofor_main

# Makefile Intel compiler with Linux

#FC90 = ifort
#FC03 = ifort
#CC = icc

#OPTS90 = -O3
#OPTS90 = -O3 -vec-report6
#OPTS90 = -O3 -axAVX -vec-report6
#OPTS90 = -O3 -mmic -vec-report6
#OPTS90 = -O3 -no-vec
#OPTS90 = -O3 -no-vec -no-simd
#OPTS90 = -O3 -r8
#OPTS90 = -O3 -CB -r8 -warn all -std90

#FF03 = -Tf
#OPTS03 = -O3 -free -std03
#OPTS03 = -O3 -free -vec-report6 -std03
#OPTS03 = -O3 -free -axAVX -vec-report6 -std03
#OPTS03 = -O3 -free -mmic -vec-report6 -std03
#OPTS03 = -O3 -free -no-vec -std03
#OPTS03 = -O3 -r8 -free -std03
#OPTS03 = -O3 -CB -r8 -warn all -free -std03

#CCOPTS = -O3 -std=c99
#CCOPTS = -O3 -mmic -std=c99
#CCOPTS = -O3 -no-vec -std=c99
#LEGACY = -nofor_main

# Makefile gfortran compiler with Linux

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
#CCOPTS = -O3 -fno-tree-vectorize -Wall -std=c99
#LEGACY =

# Makefile PGI compiler with Linux

#FC90 = pgf90
#FC03 = pgf90
#CC = gcc

#OPTS90 = -O3
#OPTS90 = -O3 -r8
#OPTS90 = -O3 -Mbounds -r8 -Mstandard

#OPTS03 = -O3
#OPTS03 = -O3 -r8
#OPTS03 = -O3 -Mbounds -r8 -Mstandard

#CCOPTS = -O3 -Wall -std=c99
#LEGACY = -Mnomain

# Makefile Cray compiler with Linux

#FC90 = ftn
#FC03 = ftn
#CC = cc

#OPTS90 = -O 3
#OPTS90 = -O 3 -s real64
#OPTS90 = -O 3 -R b -s real64 -en

#OPTS03 = -O 3
#OPTS03 = -O 3 -s real64
#OPTS03 = -O 3 -R b -s real64 -en

#CCOPTS = -O 3 -h c99 -h conform
#LEGACY =

#

# Linkage rules

#all : fvbpic3 cvbpic3 f03vbpic3
all : fvbpic3 cvbpic3

special: fvbpic3_c cvbpic3_f

fvbpic3 : fvbpic3.o fvbpush3.o cavx512lib3.o cavx512flib3.o ckncbpush3.o dtimer.o
	$(FC90) $(OPTS90) -o fvbpic3 fvbpic3.o fvbpush3.o cavx512lib3.o \
	cavx512flib3.o ckncbpush3.o avx512lib3_h.o avx512flib3_h.o vbpush3_h.o \
	dtimer.o

cvbpic3 : cvbpic3.o cvbpush3.o cavx512lib3.o ckncbpush3.o dtimer.o
	$(CC) $(CCOPTS) -o cvbpic3 cvbpic3.o cvbpush3.o cavx512lib3.o ckncbpush3.o \
	dtimer.o -lm

f03vbpic3 : f03vbpic3.o fvbpush3.o cavx512lib3.o ckncbpush3.o dtimer.o
	$(FC03) $(OPTS03) -o f03vbpic3 f03vbpic3.o fvbpush3.o cavx512lib3.o \
	ckncbpush3.o avx512lib3_c.o kncbpush3_c.o vbpush3_h.o dtimer.o

fvbpic3_c : fvbpic3_c.o cvbpush3.o cavx512lib3.o cavx512flib3.o dtimer.o
	$(FC90) $(OPTS90) -o fvbpic3_c fvbpic3_c.o cvbpush3.o cavx512lib3.o \
	cavx512flib3.o avx512flib3_h.o dtimer.o

cvbpic3_f : cvbpic3.o cvbpush3_f.o fvbpush3.o cavx512lib3.o dtimer.o
	$(FC90) $(OPTS90) $(LEGACY) -o cvbpic3_f cvbpic3.o cvbpush3_f.o fvbpush3.o \
    cavx512lib3.o dtimer.o -lm

# Compilation rules

dtimer.o : dtimer.c
	$(CC) $(CCOPTS) -c dtimer.c

fvbpush3.o : vbpush3.f
	$(FC90) $(OPTS90) -o fvbpush3.o -c vbpush3.f

vbpush3_h.o : vbpush3_h.f90
	$(FC90) $(OPTS90) -o vbpush3_h.o -c vbpush3_h.f90

cvbpush3.o : vbpush3.c
	$(CC) $(CCOPTS) -o cvbpush3.o -c vbpush3.c

cavx512lib3.o : avx512lib3.c
	$(CC) $(CCOPTS) -o cavx512lib3.o -c avx512lib3.c

cavx512flib3.o : avx512flib3.c
	$(CC) $(CCOPTS) -o cavx512flib3.o -c avx512flib3.c

ckncbpush3.o : kncbpush3.c
	$(CC) $(CCOPTS) -o ckncbpush3.o -c kncbpush3.c

avx512lib3_h.o : avx512lib3_h.f90
	$(FC90) $(OPTS90) -o avx512lib3_h.o -c avx512lib3_h.f90

avx512flib3_h.o : avx512flib3_h.f90 avx512lib3_h.o
	$(FC90) $(OPTS90) -o avx512flib3_h.o -c avx512flib3_h.f90

kncbpush3_h.o : kncbpush3_h.f90
	$(FC90) $(OPTS90) -o kncbpush3_h.o -c kncbpush3_h.f90

avx512lib3_c.o : avx512lib3_c.f03
	$(FC03) $(OPTS03) -o avx512lib3_c.o -c $(FF03) avx512lib3_c.f03

kncbpush3_c.o : kncbpush3_c.f03
	$(FC03) $(OPTS03) -o kncbpush3_c.o -c $(FF03) kncbpush3_c.f03

cvbpush3_f.o : vbpush3_f.c
	$(CC) $(CCOPTS) -o cvbpush3_f.o -c vbpush3_f.c

fvbpic3.o : vbpic3.f90 avx512flib3_h.o kncbpush3_h.o vbpush3_h.o
	$(FC90) $(OPTS90) -o fvbpic3.o -c vbpic3.f90

cvbpic3.o : vbpic3.c
	$(CC) $(CCOPTS) -o cvbpic3.o -c vbpic3.c

f03vbpic3.o : vbpic3.f03 avx512lib3_c.o kncbpush3_c.o vbpush3_h.o
	$(FC03) $(OPTS03) -o f03vbpic3.o -c $(FF03) vbpic3.f03

fvbpic3_c.o : vbpic3_c.f90 avx512flib3_h.o
	$(FC90) $(OPTS90) -o fvbpic3_c.o -c vbpic3_c.f90

clean :
	rm -f *.o *.mod

clobber: clean
	rm -f fvbpic3 cvbpic3 fvbpic3_c cvbpic3_f
