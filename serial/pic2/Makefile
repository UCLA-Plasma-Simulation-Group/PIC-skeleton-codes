#Makefile for 2D Electrostatic PIC codes

# Makefile gfortran compiler with MacOS X

#FC90 = gfortran
#CC = gcc

#OPTS90 = -O3
#OPTS90 = -O3 -fdefault-real-8 -fdefault-double-8
#OPTS90 = -O3 -fcheck=bounds -fdefault-real-8 -Wall -std=f95

#CCOPTS = -O3 -Wall -std=c99

#LEGACY =

# Makefile Intel compiler with Mac OS X

#FC90 = ifort
#CC = gcc

#OPTS90 = -O3
#OPTS90 = -O3 -r8
#OPTS90 = -O3 -CB -r8 -warn all -std90

#CCOPTS = -O3 -std=c99

#LEGACY = -nofor_main

# Makefile Intel compiler with Linux

#FC90 = ifort
#CC = icc

#OPTS90 = -O3
#OPTS90 = -O3 -r8
#OPTS90 = -O3 -CB -r8 -warn all -std90

#CCOPTS = -O3 -std=c99

#LEGACY = -nofor_main

# Makefile gfortran compiler with Linux

FC90 = gfortran
CC = gcc

OPTS90 = -O3
#OPTS90 = -O3 -fdefault-real-8 -fdefault-double-8
#OPTS90 = -O3 -fbounds-check -fdefault-real-8 -Wall -std=f95

CCOPTS = -O3 -Wall -std=c99

LEGACY =

# Makefile PGI compiler with Linux

#FC90 = pgf90
#CC = gcc

#OPTS90 = -O3
#OPTS90 = -O3 -r8
#OPTS90 = -O3 -Mbounds -r8 -Mstandard

#CCOPTS = -O3 -Wall -std=c99

#LEGACY = -Mnomain

# Makefile Nag compiler with Linux

#FC90 = nagfor
#CC = gcc

#OPTS90 = -O3
#OPTS90 = -O3 -default_kinds:r=64
#OPTS90 = -O3 -default_kinds:r=64 -C=array -nan -w=all -f95

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

all : fpic2 cpic2

special: fpic2_c cpic2_f

fpic2 : fpic2.o fpush2.o dtimer.o
	$(FC90) $(OPTS90) -o fpic2 fpic2.o fpush2.o push2_h.o \
        dtimer.o

cpic2 : cpic2.o cpush2.o dtimer.o
	$(CC) $(CCOPTS) -o cpic2 cpic2.o cpush2.o dtimer.o -lm

fpic2_c : fpic2_c.o cpush2.o dtimer.o
	$(FC90) $(OPTS90) -o fpic2_c fpic2_c.o cpush2.o dtimer.o

cpic2_f : cpic2.o cpush2_f.o fpush2.o dtimer.o
	$(FC90) $(OPTS90) $(LEGACY) -o cpic2_f cpic2.o cpush2_f.o fpush2.o \
        dtimer.o -lm

# Compilation rules

dtimer.o : dtimer.c
	$(CC) $(CCOPTS) -c dtimer.c

fpush2.o : push2.f
	$(FC90) $(OPTS90) -o fpush2.o -c push2.f

push2_h.o : push2_h.f90
	$(FC90) $(OPTS90) -o push2_h.o -c push2_h.f90

cpush2.o : push2.c
	$(CC) $(CCOPTS) -o cpush2.o -c push2.c

cpush2_f.o : push2_f.c
	$(CC) $(CCOPTS) -o cpush2_f.o -c push2_f.c

fpic2.o : pic2.f90 push2_h.o
	$(FC90) $(OPTS90) -o fpic2.o -c pic2.f90

cpic2.o : pic2.c
	$(CC) $(CCOPTS) -o cpic2.o -c pic2.c

fpic2_c.o : pic2_c.f90
	$(FC90) $(OPTS90) -o fpic2_c.o -c pic2_c.f90

clean :
	rm -f *.o *.mod

clobber: clean
	rm -f fpic2 cpic2 fpic2_c cpic2_f
