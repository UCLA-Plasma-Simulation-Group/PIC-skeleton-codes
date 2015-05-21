#Makefile for 1-2/2D Darwin PIC codes

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

all : fdpic1

special: cdpic1_f

fdpic1 : fdpic1.o fdpush1.o dtimer.o
	$(FC90) $(OPTS90) -o fdpic1 fdpic1.o fdpush1.o dpush1_h.o \
        dtimer.o

cdpic1 : cdpic1.o cdpush1.o dtimer.o
	$(CC) $(CCOPTS) -o cdpic1 cdpic1.o cdpush1.o dtimer.o -lm

fdpic1_c : fdpic1_c.o cdpush1.o dtimer.o
	$(FC90) $(OPTS90) -o fdpic1_c fdpic1_c.o cdpush1.o dtimer.o

cdpic1_f : cdpic1.o cdpush1_f.o fdpush1.o dtimer.o
	$(FC90) $(CCOPTS) -o cdpic1_f cdpic1.o cdpush1_f.o fdpush1.o \
        dtimer.o -lm

# Compilation rules

dtimer.o : dtimer.c
	$(CC) $(CCOPTS) -c dtimer.c

fdpush1.o : dpush1.f
	$(FC90) $(OPTS90) -o fdpush1.o -c dpush1.f

dpush1_h.o : dpush1_h.f90
	$(FC90) $(OPTS90) -o dpush1_h.o -c dpush1_h.f90

cdpush1.o : dpush1.c
	$(CC) $(CCOPTS) -o cdpush1.o -c dpush1.c

cdpush1_f.o : dpush1_f.c
	$(CC) $(CCOPTS) -o cdpush1_f.o -c dpush1_f.c

fdpic1.o : dpic1.f90 dpush1_h.o
	$(FC90) $(OPTS90) -o fdpic1.o -c dpic1.f90

cdpic1.o : dpic1.c
	$(CC) $(CCOPTS) -o cdpic1.o -c dpic1.c

fdpic1_c.o : dpic1_c.f90
	$(FC90) $(OPTS90) -o fdpic1_c.o -c dpic1_c.f90

clean :
	rm -f *.o *.mod

clobber: clean
	rm -f fdpic1 cdpic1 fdpic1_c cdpic1_f
