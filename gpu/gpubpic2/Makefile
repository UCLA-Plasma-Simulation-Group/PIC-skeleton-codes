#Makefile for 2-1/2D Electromagnetic GPU PIC codes

# Makefile gfortran compiler with MacOS X

#CARBON = /System/Library/Frameworks/Carbon.framework/Carbon

#FC90 = gfortran
#FC03 = gfortran
#CC = gcc

#OPTS90 = -O3
#OPTS90 = -O3 -fdefault-real-8
#OPTS90 = -O3 -fcheck=bounds -fdefault-real-8 -Wall -std=f95

#OPTS03 = -O3
#OPTS03 = -O3 -fdefault-real-8 -fdefault-double-8
#OPTS03 = -O3 -fcheck=bounds -fdefault-real-8 -fdefault-double-8 -Wall -std=f2003

#CCOPTS = -O3 -Wall -std=c99

#CUDA
#LOPTS = -L/usr/lib/gcc/i686-apple-darwin9/4.0.1 -lstdc++
#LOPTS = -L/usr/lib/gcc/i686-apple-darwin9/4.2.1 -lstdc++
#LOPTS = -L/usr/lib/gcc/i686-apple-darwin10/4.2.1 -lstdc++
#LOPTS = -L/usr/lib/gcc/i686-apple-darwin10/4.2.1/x86_64 -lstdc++
#CULIBS = $(CARBON) -lSystemStubs -L/usr/local/cuda/lib -lcuda -lcudart \
#-lcufft

#NVOPTS = -O3 -Xptxas=-v
#NVOPTS = -O3 -arch sm_12
#NVOPTS = -O3 -arch sm_20
#NVOPTS = -O3 -m64 -arch sm_12
#NVOPTS = -O3 -m64 -arch sm_20

# Makefile Intel compiler with Mac OS X

#CARBON = /System/Library/Frameworks/Carbon.framework/Carbon

#FC90 = ifort
#FC03 = ifort
#CC = gcc

#OPTS90 = -O3
#OPTS90 = -O3 -r8
#OPTS90 = -O3 -CB -r8 -warn all -std90

#FF03 = -Tf
#OPTS03 = -O3 -free -std03
#OPTS03 = -O3 -r8 -free -std03
#OPTS03 = -O3 -CB -r8 -warn all -free -std03

#CCOPTS = -O3 -Wall -std=c99

#CUDA
#LOPTS = -L/usr/lib/gcc/i686-apple-darwin9/4.0.1 -lstdc++
#LOPTS = -L/usr/lib/gcc/i686-apple-darwin9/4.2.1 -lstdc++
#LOPTS = -L/usr/lib/gcc/i686-apple-darwin10/4.2.1 -lstdc++
#LOPTS = -L/usr/lib/gcc/i686-apple-darwin10/4.2.1/x86_64 -lstdc++
#CULIBS = $(CARBON) -lSystemStubs -L/usr/local/cuda/lib -lcuda -lcudart \
#-lcufft

#NVOPTS = -O3 -arch sm_12
#NVOPTS = -O3 -arch sm_20
#NVOPTS = -O3 -m64 -arch sm_12
#NVOPTS = -O3 -m64 -arch sm_20

# Makefile Intel compiler with Linux

#FC90 = ifort
#FC03 = ifort
#CC = gcc

#OPTS90 = -O3
#OPTS90 = -O3 -r8
#OPTS90 = -O3 -CB -r8 -warn all -std90

#FF03 = -Tf
#OPTS03 = -O3 -free -std03
#OPTS03 = -O3 -r8 -free -std03
#OPTS03 = -O3 -CB -r8 -warn all -free -std03

#CCOPTS = -O3 -Wall -std=c99

#LOPTS =
#CULIBS = -L/u/local/cuda/5.0/lib64 -lcudart -lcufft
#NVOPTS = -O3 -arch sm_20

# Makefile gfortran compiler with Linux

FC90 = gfortran
FC03 = gfortran
CC = gcc

OPTS90 = -O3
#OPTS90 = -O3 -fdefault-real-8
#OPTS90 = -O3 -fbounds-check -fdefault-real-8 -Wall -std=f95

OPTS03 = -O3
#OPTS03 = -O3 -fdefault-real-8
#OPTS03 = -O3 -fbounds-check -fdefault-real-8 -fdefault-double-8 -Wall -std=f2003

CCOPTS = -O3 -Wall -std=c99

#CUDA
LOPTS =
#CULIBS = -L/usr/local/cuda/lib64 -lcuda -lcudart -lcufft
#CULIBS = -L/usr/local/cuda/lib -lcuda -lcudart -lcufft
#CULIBS = -L/u/local/cuda/4.0.17/lib64 -lcuda -lcudart -lcufft
CULIBS = -L/u/local/cuda/5.0/lib64 -lcudart -lcufft
#CULIBS = -lcudart -lcufft
#NVOPTS = -O3 -Xptxas=-v
#NVOPTS = -O3 -arch sm_12
NVOPTS = -O3 -arch sm_20

# Makefile PGI compiler with Linux

#FC90 = pgf90
#FC03 = pgf90
#CC = gcc

#OPTS90 = -O3
#OPTS90 = -O3 -r8
#OPTS90 = -O3 -Mbounds -r8 -Mstandard

#OPTS03 = -O3
#OPTS03 = -O3 -r8
#OPTS03 = -O3 -Mbounds -r8 -Mstandard

#CCOPTS = -O3 -std=c99

#CUDA
#LOPTS =
#CULIBS = -L/u/local/cuda/5.0/lib64 -lcudart -lcufft
#NVOPTS = -O3 -arch sm_20

# Makefile Cray compiler with Linux
 
#FC90 = ftn
#FC03 = ftn
#CC = nvcc

#OPTS90 = -O 3
#OPTS90 = -O 3 -s real64
#OPTS90 = -O 3 -R b -s real64 -en

#OPTS03 = -O 3
#OPTS03 = -O 3 -s real64
#OPTS03 = -O 3 -R b -s real64 -en

#CCOPTS = -O3
#CULIBS = -lcufft
#NVOPTS = -O3 -arch sm_20

#

# Linkage rules

all: cuda

cuda : fgpubpic2 cgpubpic2 f03gpubpic2

fgpubpic2 : fgpubpic2.o gpulib2.o gpubpush2.o gpufft2.o fbpush2.o dtimer.o
	$(FC90) $(OPTS90) $(LOPTS) -o fgpubpic2 fgpubpic2.o \
	gpulib2.o gpubpush2.o gpufft2.o fbpush2.o gpulib2_h.o gpubpush2_h.o \
	gpufft2_h.o bpush2_h.o dtimer.o $(CULIBS)

cgpubpic2 : cgpubpic2.o gpulib2.o gpubpush2.o gpufft2.o cbpush2.o dtimer.o
	$(CC) $(CCOPTS) $(LOPTS) -o cgpubpic2 cgpubpic2.o \
	gpulib2.o gpubpush2.o gpufft2.o cbpush2.o dtimer.o $(CULIBS) -lm

f03gpubpic2 : f03gpubpic2.o gpulib2.o gpubpush2.o gpufft2.o fbpush2.o \
              dtimer.o
	$(FC03) $(OPTS03) $(LOPTS) -o f03gpubpic2 f03gpubpic2.o \
	gpulib2.o gpubpush2.o gpufft2.o fbpush2.o gpulib2_c.o gpubpush2_c.o \
	gpufft2_c.o bpush2_h.o dtimer.o $(CULIBS)

# Compilation rules

dtimer.o : dtimer.c
	$(CC) $(CCOPTS) -c dtimer.c

gpulib2.o : gpulib2.cu
	nvcc $(NVOPTS) -c gpulib2.cu -I/usr/local/cuda/include

gpubpush2.o : gpubpush2.cu
	nvcc $(NVOPTS) -c gpubpush2.cu -I/usr/local/cuda/include

gpufft2.o : gpufft2.cu
	nvcc $(NVOPTS) -c gpufft2.cu -I/usr/local/cuda/include

gpulib2_h.o : gpulib2_h.f90
	$(FC90) $(OPTS90) -c gpulib2_h.f90

gpubpush2_h.o : gpubpush2_h.f90
	$(FC90) $(OPTS90) -c gpubpush2_h.f90

gpufft2_h.o : gpufft2_h.f90
	$(FC90) $(OPTS90) -c gpufft2_h.f90

gpulib2_c.o : gpulib2_c.f03
	$(FC03) $(OPTS03) -c $(FF03) gpulib2_c.f03

gpubpush2_c.o : gpubpush2_c.f03
	$(FC03) $(OPTS03) -c $(FF03) gpubpush2_c.f03

gpufft2_c.o : gpufft2_c.f03
	$(FC03) $(OPTS03) -c $(FF03) gpufft2_c.f03

fbpush2.o : bpush2.f
	$(FC90) $(OPTS90) -o fbpush2.o -c bpush2.f

bpush2_h.o : bpush2_h.f90
	$(FC90) $(OPTS90) -o bpush2_h.o -c bpush2_h.f90

cbpush2.o : bpush2.c
	$(CC) $(CCOPTS) -o cbpush2.o -c bpush2.c

fgpubpic2.o : gpubpic2.f90 gpulib2_h.o gpubpush2_h.o gpufft2_h.o bpush2_h.o
	$(FC90) $(OPTS90) -o fgpubpic2.o -c gpubpic2.f90

cgpubpic2.o : gpubpic2.c
	$(CC) $(CCOPTS) -o cgpubpic2.o -c gpubpic2.c

f03gpubpic2.o : gpubpic2.f03 gpulib2_c.o gpubpush2_c.o gpufft2_c.o bpush2_h.o
	$(FC03) $(OPTS03) -o f03gpubpic2.o -c $(FF03) gpubpic2.f03

clean :
	rm -f *.o *.mod

clobber: clean
	rm -f fgpubpic2 cgpubpic2 f03gpubpic2
