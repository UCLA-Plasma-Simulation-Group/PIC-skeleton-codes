#Makefile for 3D Electrostatic PIC codes

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

all : fpic3

special: cpic3_f

fpic3 : fpic3.o fpush3.o dtimer.o
	$(FC90) $(OPTS90) -o fpic3 fpic3.o fpush3.o push3_h.o \
        dtimer.o

cpic3 : cpic3.o cpush3.o dtimer.o
	$(CC) $(CCOPTS) -o cpic3 cpic3.o cpush3.o dtimer.o -lm

fpic3_c : fpic3_c.o cpush3.o dtimer.o
	$(FC90) $(OPTS90) -o fpic3_c fpic3_c.o cpush3.o dtimer.o

cpic3_f : cpic3.o cpush3_f.o fpush3.o dtimer.o
	$(FC90) $(CCOPTS) -o cpic3_f cpic3.o cpush3_f.o fpush3.o \
        dtimer.o -lm

# Compilation rules

dtimer.o : dtimer.c
	$(CC) $(CCOPTS) -c dtimer.c

fpush3.o : push3.f
	$(FC90) $(OPTS90) -o fpush3.o -c push3.f

push3_h.o : push3_h.f90
	$(FC90) $(OPTS90) -o push3_h.o -c push3_h.f90

cpush3.o : push3.c
	$(CC) $(CCOPTS) -o cpush3.o -c push3.c

cpush3_f.o : push3_f.c
	$(CC) $(CCOPTS) -o cpush3_f.o -c push3_f.c

fpic3.o : pic3.f90 push3_h.o
	$(FC90) $(OPTS90) -o fpic3.o -c pic3.f90

cpic3.o : pic3.c
	$(CC) $(CCOPTS) -o cpic3.o -c pic3.c

fpic3_c.o : pic3_c.f90
	$(FC90) $(OPTS90) -o fpic3_c.o -c pic3_c.f90

clean :
	rm -f *.o *.mod

clobber: clean
	rm -f fpic3 cpic3 fpic3_c cpic3_f
