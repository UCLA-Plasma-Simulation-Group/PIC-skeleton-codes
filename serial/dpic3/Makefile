#Makefile for 3D Darwin PIC codes

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

#OPTS90 = -O3 -dusty
#OPTS90 = -O3 -default_kinds:r=64 -dusty
#OPTS90 = -O3 -default_kinds:r=64 -C=array -nan -w=all -f95 -dusty

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

all : fdpic3 cdpic3

special: cdpic3_f fdpic3_c

fdpic3 : fdpic3.o fdpush3.o dtimer.o
	$(FC90) $(OPTS90) -o fdpic3 fdpic3.o fdpush3.o dpush3_h.o \
        dtimer.o

cdpic3 : cdpic3.o cdpush3.o dtimer.o
	$(CC) $(CCOPTS) -o cdpic3 cdpic3.o cdpush3.o dtimer.o -lm

fdpic3_c : fdpic3_c.o cdpush3.o dtimer.o
	$(FC90) $(OPTS90) -o fdpic3_c fdpic3_c.o cdpush3.o dtimer.o

cdpic3_f : cdpic3.o cdpush3_f.o fdpush3.o dtimer.o
	$(FC90) $(OPTS90) $(LEGACY) -o cdpic3_f cdpic3.o cdpush3_f.o fdpush3.o \
        dtimer.o -lm

# Compilation rules

dtimer.o : dtimer.c
	$(CC) $(CCOPTS) -c dtimer.c

fdpush3.o : dpush3.f
	$(FC90) $(OPTS90) -o fdpush3.o -c dpush3.f

dpush3_h.o : dpush3_h.f90
	$(FC90) $(OPTS90) -o dpush3_h.o -c dpush3_h.f90

cdpush3.o : dpush3.c
	$(CC) $(CCOPTS) -o cdpush3.o -c dpush3.c

cdpush3_f.o : dpush3_f.c
	$(CC) $(CCOPTS) -o cdpush3_f.o -c dpush3_f.c

fdpic3.o : dpic3.f90 dpush3_h.o
	$(FC90) $(OPTS90) -o fdpic3.o -c dpic3.f90

cdpic3.o : dpic3.c
	$(CC) $(CCOPTS) -o cdpic3.o -c dpic3.c

fdpic3_c.o : dpic3_c.f90
	$(FC90) $(OPTS90) -o fdpic3_c.o -c dpic3_c.f90

clean :
	rm -f *.o *.mod

clobber: clean
	rm -f fdpic3 cdpic3 fdpic3_c cdpic3_f
