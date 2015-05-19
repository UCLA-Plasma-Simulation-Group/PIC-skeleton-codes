#Makefile for 2-1/2D Electromagnetic Vector PIC codes

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

all : fvbpic2 cvbpic2  f03vbpic2

fvbpic2 : fvbpic2.o fvbpush2.o csselib2.o csseflib2.o cssebpush2.o \
          dtimer.o
	$(FC90) $(OPTS90) -o fvbpic2 fvbpic2.o fvbpush2.o csselib2.o \
	csseflib2.o cssebpush2.o sselib2_h.o sseflib2_h.o ssebpush2_h.o \
	vbpush2_h.o dtimer.o

cvbpic2 : cvbpic2.o cvbpush2.o csselib2.o cssebpush2.o dtimer.o
	$(CC) $(CCOPTS) -o cvbpic2 cvbpic2.o cvbpush2.o csselib2.o \
	cssebpush2.o dtimer.o -lm

f03vbpic2 : f03vbpic2.o fvbpush2.o csselib2.o cssebpush2.o dtimer.o
	$(FC03) $(OPTS03) -o f03vbpic2 f03vbpic2.o fvbpush2.o csselib2.o \
	cssebpush2.o sselib2_c.o ssebpush2_c.o vbpush2_h.o dtimer.o

# Compilation rules

dtimer.o : dtimer.c
	$(CC) $(CCOPTS) -c dtimer.c

fvbpush2.o : vbpush2.f
	$(FC90) $(OPTS90) -o fvbpush2.o -c vbpush2.f

vbpush2_h.o : vbpush2_h.f90
	$(FC90) $(OPTS90) -o vbpush2_h.o -c vbpush2_h.f90

cvbpush2.o : vbpush2.c
	$(CC) $(CCOPTS) -o cvbpush2.o -c vbpush2.c

csselib2.o : sselib2.c
	$(CC) $(CCOPTS) -o csselib2.o -c sselib2.c

csseflib2.o : sseflib2.c
	$(CC) $(CCOPTS) -o csseflib2.o -c sseflib2.c

cssebpush2.o : ssebpush2.c
	$(CC) $(CCOPTS) -o cssebpush2.o -c ssebpush2.c

sselib2_h.o : sselib2_h.f90
	$(FC90) $(OPTS90) -o sselib2_h.o -c sselib2_h.f90

sseflib2_h.o : sseflib2_h.f90 sselib2_h.o
	$(FC90) $(OPTS90) -o sseflib2_h.o -c sseflib2_h.f90

ssebpush2_h.o : ssebpush2_h.f90
	$(FC90) $(OPTS90) -o ssebpush2_h.o -c ssebpush2_h.f90

sselib2_c.o : sselib2_c.f03
	$(FC03) $(OPTS03) -o sselib2_c.o -c $(FF03) sselib2_c.f03

ssebpush2_c.o : ssebpush2_c.f03
	$(FC03) $(OPTS03) -o ssebpush2_c.o -c $(FF03) ssebpush2_c.f03

fvbpic2.o : vbpic2.f90 sseflib2_h.o ssebpush2_h.o vbpush2_h.o
	$(FC90) $(OPTS90) -o fvbpic2.o -c vbpic2.f90

cvbpic2.o : vbpic2.c
	$(CC) $(CCOPTS) -o cvbpic2.o -c vbpic2.c

f03vbpic2.o : vbpic2.f03 sselib2_c.o ssebpush2_c.o vbpush2_h.o
	$(FC03) $(OPTS03) -o f03vbpic2.o -c $(FF03) vbpic2.f03

clean :
	rm -f *.o *.mod

clobber: clean
	rm -f fvbpic2 cvbpic2 f03vbpic2
