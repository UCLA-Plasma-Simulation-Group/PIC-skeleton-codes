# Makefile for CUDA test programs

# Makefile gfortran compiler with MacOS X

CARBON = /System/Library/Frameworks/Carbon.framework/Carbon

#FC90 = gfortran
#CC = gcc

#OPTS90 = -O3
#CCOPTS = -O3

#CUDA
#LOPTS = -L/usr/lib/gcc/i686-apple-darwin9/4.0.1 -lstdc++
#LOPTS = -L/usr/lib/gcc/i686-apple-darwin9/4.2.1 -lstdc++
#LOPTS = -L/usr/lib/gcc/i686-apple-darwin10/4.2.1 -lstdc++
#LOPTS = -L/usr/lib/gcc/i686-apple-darwin10/4.2.1/x86_64 -lstdc++
#CULIBS = $(CARBON) -lSystemStubs -L/usr/local/cuda/lib -lcuda -lcudart
#normal CUDA
#NVOPTS = -O3
#NVOPTS = -O3 -m64
#debug CUDA
#NVOPTS = -O -deviceemu

# Makefile gfortran compiler with Linux

FC90 = gfortran
CC = gcc

OPTS90 = -O3
CCOPTS = -O3

#CUDA
LOPTS =
CULIBS = -L/usr/local/cuda/lib64 -lcuda -lcudart
#CULIBS = -L/usr/local/cuda/lib -lcuda -lcudart
#CULIBS = -L/u/local/cuda/4.0.17/lib64 -lcuda -lcudart
#normal CUDA
NVOPTS = -O3
#debug CUDA
#NVOPTS = -O -deviceemu

#

CUOBJS = gpumain_cu.o dtimer.o

FCUOBJS = gpumain_cuf.o dtimer.o


# Linkage rules

all: cuda

cuda: cgputest_cu fgputest_cu

#CUDA
cgputest_cu : cgputest_cu.o $(CUOBJS)
	$(CC) $(CCOPTS) $(LOPTS) -o cgputest_cu cgputest_cu.o \
      $(CUOBJS) $(CULIBS)

fgputest_cu : fgputest_cu.o $(CUOBJS)
	$(FC90) $(OPTS90) $(LOPTS) -o fgputest_cu fgputest_cu.o \
      $(CUOBJS) $(CULIBS)

# Compilation rules

dtimer.o : dtimer.c
	$(CC) $(CCOPTS) -c dtimer.c

#CUDA
gpumain_cu.o : gpumain_cu.cu
	nvcc $(NVOPTS) -c gpumain_cu.cu -I/usr/local/cuda/include

cgputest_cu.o : cgputest_cu.c
	$(CC) $(CCOPTS) -c cgputest_cu.c

fgputest_cu.o : fgputest_cu.f90
	$(FC90) $(OPTS90) -c fgputest_cu.f90

clean:
	rm -f *.o *.mod

clobber: clean
	rm -f *gputest_cu
