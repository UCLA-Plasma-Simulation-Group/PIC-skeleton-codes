#Makefile for 3D Electromagnetic PIC codes

# Makefile gfortran compiler with MacOS X

#FC90 = gfortran
#CC = gcc

#OPTS90 = -O3
#OPTS90 = -O3 -fdefault-real-8 -fdefault-double-8
#OPTS90 = -O3 -fcheck=bounds -fdefault-real-8 -Wall -std=f95

#CCOPTS = -O3 -Wall -std=c99

# Makefile Intel compiler with Mac OS X

#FC90 = ifort
#CC = icc

#OPTS90 = -O3
#OPTS90 = -O3 -r8
#OPTS90 = -O3 -CB -r8 -warn all -std90

#CCOPTS = -O3 -std=c99

# Makefile Intel compiler with Linux

#FC90 = ifort
#CC = icc

#OPTS90 = -O3
#OPTS90 = -O3 -r8
#OPTS90 = -O3 -CB -r8 -warn all -std90

#CCOPTS = -O3 -std=c99

# Makefile gfortran compiler with Linux

FC90 = gfortran
CC = gcc

OPTS90 = -O3
#OPTS90 = -O3 -fdefault-real-8 -fdefault-double-8
#OPTS90 = -O3 -fbounds-check -fdefault-real-8 -Wall -std=f95

CCOPTS = -O3 -Wall -std=c99

# Makefile PGI compiler with Linux

#FC90 = pgf90
#CC = gcc

#OPTS90 = -O3
#OPTS90 = -O3 -r8
#OPTS90 = -O3 -Mbounds -r8 -Mstandard

#CCOPTS = -O3 -Wall -std=c99
#LEGACY =

# Makefile Cray compiler with Linux

#FC90 = ftn
#CC = cc

#OPTS90 = -O 3
#OPTS90 = -O 3 -s real64
#OPTS90 = -O 3 -R b -s real64 -en

#CCOPTS = -O 3 -h c99 -h conform
#LEGACY =

#

# Linkage rules

all : fbpic3 cbpic3

special: cbpic3_f fbpic3_c

fbpic3 : fbpic3.o fbpush3.o dtimer.o
	$(FC90) $(OPTS90) -o fbpic3 fbpic3.o fbpush3.o bpush3_h.o \
        dtimer.o

cbpic3 : cbpic3.o cbpush3.o dtimer.o
	$(CC) $(CCOPTS) -o cbpic3 cbpic3.o cbpush3.o dtimer.o -lm

fbpic3_c : fbpic3_c.o cbpush3.o dtimer.o 
	$(FC90) $(OPTS90) -o fbpic3_c fbpic3_c.o cbpush3.o dtimer.o 

cbpic3_f : cbpic3.o cbpush3_f.o fbpush3.o dtimer.o
	$(FC90) $(CCOPTS) -o cbpic3_f cbpic3.o cbpush3_f.o fbpush3.o \
        dtimer.o -lm

# Compilation rules

dtimer.o : dtimer.c
	$(CC) $(CCOPTS) -c dtimer.c

fbpush3.o : bpush3.f
	$(FC90) $(OPTS90) -o fbpush3.o -c bpush3.f

bpush3_h.o : bpush3_h.f90
	$(FC90) $(OPTS90) -o bpush3_h.o -c bpush3_h.f90

cbpush3.o : bpush3.c
	$(CC) $(CCOPTS) -o cbpush3.o -c bpush3.c

cbpush3_f.o : bpush3_f.c
	$(CC) $(CCOPTS) -o cbpush3_f.o -c bpush3_f.c

fbpic3.o : bpic3.f90 bpush3_h.o
	$(FC90) $(OPTS90) -o fbpic3.o -c bpic3.f90

cbpic3.o : bpic3.c
	$(CC) $(CCOPTS) -o cbpic3.o -c bpic3.c

fbpic3_c.o : bpic3_c.f90
	$(FC90) $(OPTS90) -o fbpic3_c.o -c bpic3_c.f90

clean :
	rm -f *.o *.mod

clobber: clean
	rm -f fbpic3 cbpic3 fbpic3_c cbpic3_f
