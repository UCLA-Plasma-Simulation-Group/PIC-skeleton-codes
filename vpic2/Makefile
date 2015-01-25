#Makefile for 2D Electrostatic Vector PIC codes

# Makefile gfortran compiler with MacOS X

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

# Makefile Intel compiler with Mac OS X

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
#CCOPTS = -O3 -fno-tree-vectorize -Wall -std=c99

# Makefile Intel compiler with Linux

FC90 = ifort
FC03 = ifort
CC = icc

OPTS90 = -O3
#OPTS90 = -O3 -vec-report6
#OPTS90 = -O3 -no-vec
#OPTS90 = -O3 -r8
#OPTS90 = -O3 -CB -r8 -warn all -std90

FF03 = -Tf
OPTS03 = -O3 -free -std03
#OPTS03 = -O3 -free -vec-report6 -std03
#OPTS03 = -O3 -free -no-vec -std03
#OPTS03 = -O3 -r8 -free -std03
#OPTS03 = -O3 -CB -r8 -warn all -free -std03

CCOPTS = -O3 -std=c99
#CCOPTS = -O3 -no-vec -std=c99

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

#CCOPTS = -O3 -ftree-vectorize -Wall -std=c99
#CCOPTS = -O3 -fno-tree-vectorize -Wall -std=c99

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

#

# Linkage rules

all : fvpic2 cvpic2 f03vpic2

fvpic2 : fvpic2.o fvpush2.o csselib2.o csseflib2.o cssepush2.o dtimer.o
	$(FC90) $(OPTS90) -o fvpic2 fvpic2.o fvpush2.o csselib2.o csseflib2.o \
    cssepush2.o sselib2_h.o sseflib2_h.o ssepush2_h.o vpush2_h.o dtimer.o

cvpic2 : cvpic2.o cvpush2.o csselib2.o cssepush2.o dtimer.o
	$(CC) $(CCOPTS) -o cvpic2 cvpic2.o cvpush2.o csselib2.o cssepush2.o \
	dtimer.o -lm

f03vpic2 : f03vpic2.o fvpush2.o csselib2.o cssepush2.o dtimer.o
	$(FC03) $(OPTS03) -o f03vpic2 f03vpic2.o fvpush2.o csselib2.o \
    cssepush2.o sselib2_c.o ssepush2_c.o vpush2_h.o dtimer.o

# Compilation rules

dtimer.o : dtimer.c
	$(CC) $(CCOPTS) -c dtimer.c

fvpush2.o : vpush2.f
	$(FC90) $(OPTS90) -o fvpush2.o -c vpush2.f

vpush2_h.o : vpush2_h.f90
	$(FC90) $(OPTS90) -o vpush2_h.o -c vpush2_h.f90

cvpush2.o : vpush2.c
	$(CC) $(CCOPTS) -o cvpush2.o -c vpush2.c

csselib2.o : sselib2.c
	$(CC) $(CCOPTS) -o csselib2.o -c sselib2.c

csseflib2.o : sseflib2.c
	$(CC) $(CCOPTS) -o csseflib2.o -c sseflib2.c

cssepush2.o : ssepush2.c
	$(CC) $(CCOPTS) -o cssepush2.o -c ssepush2.c

sselib2_h.o : sselib2_h.f90
	$(FC90) $(OPTS90) -o sselib2_h.o -c sselib2_h.f90

sseflib2_h.o : sseflib2_h.f90 sselib2_h.o
	$(FC90) $(OPTS90) -o sseflib2_h.o -c sseflib2_h.f90

ssepush2_h.o : ssepush2_h.f90
	$(FC90) $(OPTS90) -o ssepush2_h.o -c ssepush2_h.f90

sselib2_c.o : sselib2_c.f03
	$(FC03) $(OPTS03) -o sselib2_c.o -c $(FF03) sselib2_c.f03

ssepush2_c.o : ssepush2_c.f03
	$(FC03) $(OPTS03) -o ssepush2_c.o -c $(FF03) ssepush2_c.f03

fvpic2.o : vpic2.f90 sseflib2_h.o ssepush2_h.o vpush2_h.o
	$(FC90) $(OPTS90) -o fvpic2.o -c vpic2.f90

cvpic2.o : vpic2.c
	$(CC) $(CCOPTS) -o cvpic2.o -c vpic2.c

f03vpic2.o : vpic2.f03 sselib2_c.o ssepush2_c.o vpush2_h.o 
	$(FC03) $(OPTS03) -o f03vpic2.o -c $(FF03) vpic2.f03

clean :
	rm -f *.o *.mod

clobber: clean
	rm -f fvpic2 cvpic2 f03vpic2
