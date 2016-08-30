#Makefile for 1-2/2D Electromagnetic PIC codes

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

all : fbpic1

special: cbpic1_f

fbpic1 : fbpic1.o fbpush1.o dtimer.o
	$(FC90) $(OPTS90) -o fbpic1 fbpic1.o fbpush1.o bpush1_h.o \
        dtimer.o

cbpic1 : cbpic1.o cbpush1.o dtimer.o
	$(CC) $(CCOPTS) -o cbpic1 cbpic1.o cbpush1.o dtimer.o -lm

fbpic1_c : fbpic1_c.o cbpush1.o dtimer.o
	$(FC90) $(OPTS90) -o fbpic1_c fbpic1_c.o cbpush1.o dtimer.o

cbpic1_f : cbpic1.o cbpush1_f.o fbpush1.o dtimer.o
	$(FC90) $(OPTS90) $(LEGACY) -o cbpic1_f cbpic1.o cbpush1_f.o fbpush1.o \
        dtimer.o -lm

# Compilation rules

dtimer.o : dtimer.c
	$(CC) $(CCOPTS) -c dtimer.c

fbpush1.o : bpush1.f
	$(FC90) $(OPTS90) -o fbpush1.o -c bpush1.f

bpush1_h.o : bpush1_h.f90
	$(FC90) $(OPTS90) -o bpush1_h.o -c bpush1_h.f90

cbpush1.o : bpush1.c
	$(CC) $(CCOPTS) -o cbpush1.o -c bpush1.c

cbpush1_f.o : bpush1_f.c
	$(CC) $(CCOPTS) -o cbpush1_f.o -c bpush1_f.c

fbpic1.o : bpic1.f90 bpush1_h.o
	$(FC90) $(OPTS90) -o fbpic1.o -c bpic1.f90

cbpic1.o : bpic1.c
	$(CC) $(CCOPTS) -o cbpic1.o -c bpic1.c

fbpic1_c.o : bpic1_c.f90
	$(FC90) $(OPTS90) -o fbpic1_c.o -c bpic1_c.f90

clean :
	rm -f *.o *.mod

clobber: clean
	rm -f fbpic1 cbpic1 fbpic1_c cbpic1_f
