#Makefile for 3D Electrostatic Vector PIC codes

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

all : fvpic3 cvpic3 f03vpic3

special: fvpic3_c cvpic3_f

fvpic3 : fvpic3.o fvpush3.o cavx512lib3.o cavx512flib3.o ckncpush3.o dtimer.o
	$(FC90) $(OPTS90) -o fvpic3 fvpic3.o fvpush3.o cavx512lib3.o \
	cavx512flib3.o ckncpush3.o avx512lib3_h.o avx512flib3_h.o vpush3_h.o \
	dtimer.o

cvpic3 : cvpic3.o cvpush3.o cavx512lib3.o ckncpush3.o dtimer.o
	$(CC) $(CCOPTS) -o cvpic3 cvpic3.o cvpush3.o cavx512lib3.o ckncpush3.o \
	dtimer.o -lm

f03vpic3 : f03vpic3.o fvpush3.o cavx512lib3.o ckncpush3.o dtimer.o
	$(FC03) $(OPTS03) -o f03vpic3 f03vpic3.o fvpush3.o cavx512lib3.o \
	ckncpush3.o avx512lib3_c.o kncpush3_c.o vpush3_h.o dtimer.o

fvpic3_c : fvpic3_c.o cvpush3.o cavx512lib3.o cavx512flib3.o dtimer.o
	$(FC90) $(OPTS90) -o fvpic3_c fvpic3_c.o cvpush3.o cavx512lib3.o \
	cavx512flib3.o avx512flib3_h.o dtimer.o

cvpic3_f : cvpic3.o cvpush3_f.o fvpush3.o cavx512lib3.o cdtimer.o
	$(FC90) $(OPTS90) $(LEGACY) -o cvpic3_f cvpic3.o cvpush3_f.o fvpush3.o \
    cavx512lib3.o dtimer.o -lm

# Compilation rules

dtimer.o : dtimer.c
	$(CC) $(CCOPTS) -c dtimer.c

fvpush3.o : vpush3.f
	$(FC90) $(OPTS90) -o fvpush3.o -c vpush3.f

vpush3_h.o : vpush3_h.f90
	$(FC90) $(OPTS90) -o vpush3_h.o -c vpush3_h.f90

cvpush3.o : vpush3.c
	$(CC) $(CCOPTS) -o cvpush3.o -c vpush3.c

cavx512lib3.o : avx512lib3.c
	$(CC) $(CCOPTS) -o cavx512lib3.o -c avx512lib3.c

cavx512flib3.o : avx512flib3.c
	$(CC) $(CCOPTS) -o cavx512flib3.o -c avx512flib3.c

ckncpush3.o : kncpush3.c
	$(CC) $(CCOPTS) -o ckncpush3.o -c kncpush3.c

avx512lib3_h.o : avx512lib3_h.f90
	$(FC90) $(OPTS90) -o avx512lib3_h.o -c avx512lib3_h.f90

avx512flib3_h.o : avx512flib3_h.f90 avx512lib3_h.o
	$(FC90) $(OPTS90) -o avx512flib3_h.o -c avx512flib3_h.f90

kncpush3_h.o : kncpush3_h.f90
	$(FC90) $(OPTS90) -o kncpush3_h.o -c kncpush3_h.f90

avx512lib3_c.o : avx512lib3_c.f03
	$(FC03) $(OPTS03) -o avx512lib3_c.o -c $(FF03) avx512lib3_c.f03

kncpush3_c.o : kncpush3_c.f03
	$(FC03) $(OPTS03) -o kncpush3_c.o -c $(FF03) kncpush3_c.f03

cvpush3_f.o : vpush3_f.c
	$(CC) $(CCOPTS) -o cvpush3_f.o -c vpush3_f.c

fvpic3.o : vpic3.f90 avx512flib3_h.o kncpush3_h.o vpush3_h.o
	$(FC90) $(OPTS90) -o fvpic3.o -c vpic3.f90

cvpic3.o : vpic3.c
	$(CC) $(CCOPTS) -o cvpic3.o -c vpic3.c

f03vpic3.o : vpic3.f03 avx512lib3_c.o kncpush3_c.o vpush3_h.o
	$(FC03) $(OPTS03) -o f03vpic3.o -c $(FF03) vpic3.f03

fvpic3_c.o : vpic3_c.f90 avx512flib3_h.o
	$(FC90) $(OPTS90) -o fvpic3_c.o -c vpic3_c.f90

clean :
	rm -f *.o *.mod

clobber: clean
	rm -f fvpic3 cvpic3 fvpic3_c cvpic3_f
