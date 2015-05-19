# Makefile for CUDA Fortran test program

# Makefile gfortran compiler with Linux

CC = gcc

OPTS90 = -O3
CCOPTS = -O3

#CUDA
LOPTS =
CULIBS = -L/usr/local/cuda/lib64 -lcuda -lcudart
#CULIBS = -L/usr/local/cuda/lib -lcuda -lcudart
#normal CUDA
NVOPTS = -O3
#debug CUDA
#NVOPTS = -O -deviceemu
#

FCUOBJS = gpumain_cuf.o dtimer.o

# Linkage rules

all: fgputest_cuf

#CUDA FORTRAN
fgputest_cuf : fgputest_cuf.o $(FCUOBJS)
	pgf90 -Mcuda $(OPTS90) $(LOPTS) -o fgputest_cuf fgputest_cuf.o \
      $(FCUOBJS) $(CULIBS)

# Compilation rules

dtimer.o : dtimer.c
	$(CC) $(CCOPTS) -c dtimer.c

#CUDA FORTRAN
gpumain_cuf.o : gpumain_cuf.cuf
	pgf90 $(OPTS90) -c gpumain_cuf.cuf

fgputest_cuf.o : fgputest_cuf.f90 gpumain_cuf.o
	pgf90 -Mcuda $(OPTS90) -c fgputest_cuf.f90

clean:
	rm -f *.o *.mod

clobber: clean
	rm -f fgputest_cuf
