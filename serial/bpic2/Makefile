#Makefile for 2-1/2D Electromagnetic PIC codes

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

all : fbpic2 cbpic2

special: fbpic2_c cbpic2_f

fbpic2 : fbpic2.o fbpush2.o dtimer.o
	$(FC90) $(OPTS90) -o fbpic2 fbpic2.o fbpush2.o bpush2_h.o \
        dtimer.o

cbpic2 : cbpic2.o cbpush2.o dtimer.o
	$(CC) $(CCOPTS) -o cbpic2 cbpic2.o cbpush2.o dtimer.o -lm

fbpic2_c : fbpic2_c.o cbpush2.o dtimer.o
	$(FC90) $(OPTS90) -o fbpic2_c fbpic2_c.o cbpush2.o dtimer.o

cbpic2_f : cbpic2.o cbpush2_f.o fbpush2.o dtimer.o
	$(FC90) $(CCOPTS) -o cbpic2_f cbpic2.o cbpush2_f.o fbpush2.o \
        dtimer.o -lm

# Compilation rules

dtimer.o : dtimer.c
	$(CC) $(CCOPTS) -c dtimer.c

fbpush2.o : bpush2.f
	$(FC90) $(OPTS90) -o fbpush2.o -c bpush2.f

bpush2_h.o : bpush2_h.f90
	$(FC90) $(OPTS90) -o bpush2_h.o -c bpush2_h.f90

cbpush2.o : bpush2.c
	$(CC) $(CCOPTS) -o cbpush2.o -c bpush2.c

cbpush2_f.o : bpush2_f.c
	$(CC) $(CCOPTS) -o cbpush2_f.o -c bpush2_f.c

fbpic2.o : bpic2.f90 bpush2_h.o
	$(FC90) $(OPTS90) -o fbpic2.o -c bpic2.f90

cbpic2.o : bpic2.c
	$(CC) $(CCOPTS) -o cbpic2.o -c bpic2.c

fbpic2_c.o : bpic2_c.f90
	$(FC90) $(OPTS90) -o fbpic2_c.o -c bpic2_c.f90

clean :
	rm -f *.o *.mod

clobber: clean
	rm -f fbpic2 cbpic2 fbpic2_c cbpic2_f
