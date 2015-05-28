#Makefile for GPU Tutorial

# Makefile gfortran compiler with Linux

FC90 = gfortran
CC = gcc

OPTS90 = -O3
#OPTS90 = -O3 -fcheck=bounds -Wall -std=f95
#OPTS90 = -O3 -fdefault-real-8
#OPTS90 = -O3 -fcheck=bounds -fdefault-real-8 -Wall -std=f95

CCOPTS = -O3 -Wall

# Makefile Cuda Fortran PGI compiler with Linux

FCU90 = pgf90 -Mcuda 

OPTS90 = -O3
#OPTS90 = -O3 -Mbounds -Mstandard
#OPTS90 = -O3 -r8
#OPTS90 = -O3 -Mbounds -r8 -Mstandard

#CUDA for MacOS
#CUDAHOME = /usr/local/cuda
#FCUDAHOME = /usr/local/cuda

#NVOPTS = -O3 -arch sm_20
#NVOPTS = -O3 -m64 -arch sm_20

#LOPTS = -L/usr/lib/gcc/i686-apple-darwin9/4.0.1 -lstdc++
#LOPTS = -L/usr/lib/gcc/i686-apple-darwin9/4.2.1 -lstdc++
#LOPTS = -L/usr/lib/gcc/i686-apple-darwin10/4.2.1 -lstdc++
#LOPTS = -L/usr/lib/gcc/i686-apple-darwin10/4.2.1/x86_64 -lstdc++
#CULIBS = -lSystemStubs -L$(CUDAHOME)/lib -lcuda -lcudart

#CUDA for Linux
CUDAHOME = /u/local/cuda/5.5

FNVOPTS = -O3 -ta=nvidia:cc20
NVOPTS = -O3 -arch sm_20

CULIBS = -L$(CUDAHOME)/lib64 -lcuda -lcudart
#CULIBS = -L$(CUDAHOME)/lib -lcuda -lcudart

#

# Linkage rules

all : cuda fccuda

cuda : cexample1 cexample2  cexample3

fccuda : fcexample1 fcexample2 fcexample3

cudaf : fexample1 fexample2 fexample3

#CUDA C
cexample1 : cexample1.o ccopy.o gpulib2.o dtimer.o
	$(CC) $(LOPTS) -o cexample1 cexample1.o ccopy.o gpulib2.o \
	dtimer.o $(CULIBS)

cexample2 : cexample2.o ctranspose.o gpulib2.o dtimer.o
	$(CC) $(LOPTS) -o cexample2 cexample2.o ctranspose.o gpulib2.o \
	dtimer.o $(CULIBS)

cexample3 : cexample3.o credux.o gpulib2.o dtimer.o
	$(CC) $(LOPTS) -o cexample3 cexample3.o credux.o gpulib2.o \
	dtimer.o $(CULIBS)

#Fortran2003 with CUDA C
fcexample1 : fcexample1.o gpulib2.o dtimer.o
	$(FC90) $(LOPTS) -o fcexample1 fcexample1.o fccopy.o gpulib2_c.o \
	gpulib2.o dtimer.o $(CULIBS)

fcexample2 : fcexample2.o gpulib2.o dtimer.o
	$(FC90) $(LOPTS) -o fcexample2 fcexample2.o fctranspose.o gpulib2_c.o \
	gpulib2.o dtimer.o $(CULIBS)

fcexample3 : fcexample3.o gpulib2.o dtimer.o
	$(FC90) $(LOPTS) -o fcexample3 fcexample3.o fcredux.o gpulib2_c.o \
	gpulib2.o dtimer.o $(CULIBS)

#CUDA Fortran
fexample1 : fexample1.o dtimer.o
	$(FCU90) -o fexample1 fexample1.o fcopy.o gpuflib2.o \
	dtimer.o

fexample2 : fexample2.o dtimer.o
	$(FCU90) -o fexample2 fexample2.o ftranspose.o gpuflib2.o \
	dtimer.o

fexample3 : fexample3.o dtimer.o
	$(FCU90) -o fexample3 fexample3.o fredux.o gpuflib2.o \
	dtimer.o

# Compilation rules

dtimer.o : dtimer.c
	$(CC) $(CCOPTS) -c dtimer.c

#CUDA C
gpulib2.o : gpulib2.cu
	nvcc $(NVOPTS) -c gpulib2.cu -I$(CUDAHOME)/include

ccopy.o : copy.c
	$(CC) $(CCOPTS) -o ccopy.o -c copy.c

ctranspose.o : transpose.c
	$(CC) $(CCOPTS) -o ctranspose.o -c transpose.c

credux.o : redux.c
	$(CC) $(CCOPTS) -o credux.o -c redux.c

cexample1.o : example1.c
	$(CC) $(CCOPTS) -o cexample1.o -c example1.c

cexample2.o : example2.c
	$(CC) $(CCOPTS) -o cexample2.o -c example2.c

cexample3.o : example3.c
	$(CC) $(CCOPTS) -o cexample3.o -c example3.c

#Fortran with CUDA C
gpulib2_h.o : gpulib2_h.f90
	$(FC90) $(OPTS90) -c gpulib2_h.f90

gpulib2_c.o : gpulib2_c.f03
	$(FC90) $(OPTS90) -c gpulib2_c.f03

fccopy.o : copy.f90
	$(FC90) $(OPTS90) -o fccopy.o -c copy.f90

fctranspose.o : transpose.f90
	$(FC90) $(OPTS90) -o fctranspose.o -c transpose.f90

fcredux.o : redux.f90
	$(FC90) $(OPTS90) -o fcredux.o -c redux.f90

fcexample1.o : fcexample1.f03 fccopy.o gpulib2_c.o
	$(FC90) $(OPTS90) -c fcexample1.f03

fcexample2.o : fcexample2.f03 fctranspose.o gpulib2_c.o
	$(FC90) $(OPTS90) -c fcexample2.f03
	
fcexample3.o : fcexample3.f03 fcredux.o gpulib2_c.o
	$(FC90) $(OPTS90) -c fcexample3.f03

#CUDA Fortran
gpuflib2.o : gpuflib2.cuf
	$(FCU90) $(FNVOPTS) -c gpuflib2.cuf -I$(CUDAHOME)/include

fcopy.o : copy.f90
	$(FCU90) $(OPTS90) -o fcopy.o -c copy.f90

ftranspose.o : transpose.f90
	$(FCU90) $(OPTS90) -o ftranspose.o -c transpose.f90

fredux.o : redux.f90
	$(FCU90) $(OPTS90) -o fredux.o -c redux.f90

fexample1.o : example1.f90 fcopy.o gpuflib2.o
	$(FCU90) $(OPTS90) -o fexample1.o -c example1.f90

fexample2.o : example2.f90 ftranspose.o gpuflib2.o
	$(FCU90) $(OPTS90) -o fexample2.o -c example2.f90

fexample3.o : example3.f90 fredux.o gpuflib2.o
	$(FCU90) $(OPTS90) -o fexample3.o -c example3.f90

clean :
	rm -f *.o *.mod

clobber: clean
	rm -f cexample1 cexample2 cexample3 \
	      fcexample1 fcexample2 fcexample3 \
	      fexample1 fexample2 fexample3
