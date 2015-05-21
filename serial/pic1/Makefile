#Makefile for 1D Electrostatic PIC codes

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

all : fpic1

special: cpic1_f

fpic1 : fpic1.o fpush1.o dtimer.o
	$(FC90) $(OPTS90) -o fpic1 fpic1.o fpush1.o push1_h.o \
        dtimer.o

cpic1 : cpic1.o cpush1.o dtimer.o
	$(CC) $(CCOPTS) -o cpic1 cpic1.o cpush1.o dtimer.o -lm

fpic1_c : fpic1_c.o cpush1.o dtimer.o
	$(FC90) $(OPTS90) -o fpic1_c fpic1_c.o cpush1.o dtimer.o

cpic1_f : cpic1.o cpush1_f.o fpush1.o dtimer.o
	$(FC90) $(CCOPTS) -o cpic1_f cpic1.o cpush1_f.o fpush1.o \
        dtimer.o -lm

# Compilation rules

dtimer.o : dtimer.c
	$(CC) $(CCOPTS) -c dtimer.c

fpush1.o : push1.f
	$(FC90) $(OPTS90) -o fpush1.o -c push1.f

push1_h.o : push1_h.f90
	$(FC90) $(OPTS90) -o push1_h.o -c push1_h.f90

cpush1.o : push1.c
	$(CC) $(CCOPTS) -o cpush1.o -c push1.c

cpush1_f.o : push1_f.c
	$(CC) $(CCOPTS) -o cpush1_f.o -c push1_f.c

fpic1.o : pic1.f90 push1_h.o
	$(FC90) $(OPTS90) -o fpic1.o -c pic1.f90

cpic1.o : pic1.c
	$(CC) $(CCOPTS) -o cpic1.o -c pic1.c

fpic1_c.o : pic1_c.f90
	$(FC90) $(OPTS90) -o fpic1_c.o -c pic1_c.f90

clean :
	rm -f *.o *.mod

clobber: clean
	rm -f fpic1 cpic1 fpic1_c cpic1_f
