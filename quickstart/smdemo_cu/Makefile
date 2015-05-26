# Makefile for CUDA test programs

# Makefile gfortran compiler with MacOS X

#FC90 = gfortran
#CC = gcc

#OPTS90 = -O3
#CCOPTS = -O3

#CUDA
#CUDAHOME = /usr/local/cuda
#LOPTS = -L/usr/lib/gcc/i686-apple-darwin9/4.0.1 -lstdc++
#LOPTS = -L/usr/lib/gcc/i686-apple-darwin9/4.2.1 -lstdc++
#LOPTS = -L/usr/lib/gcc/i686-apple-darwin10/4.2.1 -lstdc++
#LOPTS = -L/usr/lib/gcc/i686-apple-darwin10/4.2.1/x86_64 -lstdc++
#CULIBS = -lSystemStubs -L$(CUDAHOME)/lib -lcuda -lcudart
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
CUDAHOME = /u/local/cuda/5.5
LOPTS =
CULIBS = -L$(CUDAHOME)/lib64 -lcuda -lcudart
#CULIBS = -L$(CUDAHOME)/lib -lcuda -lcudart
#normal CUDA
NVOPTS = -O3
#debug CUDA
#NVOPTS = -O -deviceemu

#

CUOBJS = gpulib_cu.o dtimer.o

# Linkage rules

all: cuda

cuda: cgputest_cu fgputest_cu

#CUDA
cgputest_cu : cgputest_cu.o $(CUOBJS)
	$(CC) $(CCOPTS) $(LOPTS) -o cgputest_cu cgputest_cu.o \
      $(CUOBJS) $(CULIBS)

fgputest_cu : fgputest_cu.o $(CUOBJS)
	$(FC90) $(OPTS90) $(LOPTS) -o fgputest_cu fgputest_cu.o gpulib_c.o \
      $(CUOBJS) $(CULIBS)

# Compilation rules

dtimer.o : dtimer.c
	$(CC) $(CCOPTS) -c dtimer.c

#CUDA
gpulib_cu.o : gpulib_cu.cu
	nvcc $(NVOPTS) -c gpulib_cu.cu -I$(CUDAHOME)/include

gpulib_c.o : gpulib_c.f03
	$(FC90) $(OPTS90) -o gpulib_c.o -c gpulib_c.f03

cgputest_cu.o : cgputest_cu.c
	$(CC) $(CCOPTS) -c cgputest_cu.c

fgputest_cu.o : fgputest_cu.f03 gpulib_c.o
	$(FC90) $(OPTS90) -c fgputest_cu.f03

clean:
	rm -f *.o *.mod

clobber: clean
	rm -f *gputest_cu
