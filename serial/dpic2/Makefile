#Makefile for 2-1/2D Darwin PIC codes

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

all : fdpic2 cdpic2

special: cdpic2_f fdpic2_c

fdpic2 : fdpic2.o fdpush2.o dtimer.o
	$(FC90) $(OPTS90) -o fdpic2 fdpic2.o fdpush2.o dpush2_h.o \
        dtimer.o

cdpic2 : cdpic2.o cdpush2.o dtimer.o
	$(CC) $(CCOPTS) -o cdpic2 cdpic2.o cdpush2.o dtimer.o -lm

fdpic2_c : fdpic2_c.o cdpush2.o dtimer.o
	$(FC90) $(OPTS90) -o fdpic2_c fdpic2_c.o cdpush2.o dtimer.o

cdpic2_f : cdpic2.o cdpush2_f.o fdpush2.o dtimer.o
	$(FC90) $(CCOPTS) -o cdpic2_f cdpic2.o cdpush2_f.o fdpush2.o \
        dtimer.o -lm

# Compilation rules

dtimer.o : dtimer.c
	$(CC) $(CCOPTS) -c dtimer.c

fdpush2.o : dpush2.f
	$(FC90) $(OPTS90) -o fdpush2.o -c dpush2.f

dpush2_h.o : dpush2_h.f90
	$(FC90) $(OPTS90) -o dpush2_h.o -c dpush2_h.f90

cdpush2.o : dpush2.c
	$(CC) $(CCOPTS) -o cdpush2.o -c dpush2.c

cdpush2_f.o : dpush2_f.c
	$(CC) $(CCOPTS) -o cdpush2_f.o -c dpush2_f.c

fdpic2.o : dpic2.f90 dpush2_h.o
	$(FC90) $(OPTS90) -o fdpic2.o -c dpic2.f90

cdpic2.o : dpic2.c
	$(CC) $(CCOPTS) -o cdpic2.o -c dpic2.c

fdpic2_c.o : dpic2_c.f90
	$(FC90) $(OPTS90) -o fdpic2_c.o -c dpic2_c.f90

clean :
	rm -f *.o *.mod

clobber: clean
	rm -f fdpic2 cdpic2 fdpic2_c cdpic2_f
