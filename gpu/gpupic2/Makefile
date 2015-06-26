#Makefile for 2D Electrostatic GPU PIC codes

# Makefile gfortran compiler with MacOS X

#CARBON = /System/Library/Frameworks/Carbon.framework/Carbon

#FC90 = gfortran
#FC03 = gfortran
#CC = gcc

#OPTS90 = -O3
#OPTS90 = -O3 -fdefault-real-8 -fdefault-double-8
#OPTS90 = -O3 -fcheck=bounds -fdefault-real-8 -fdefault-double-8 -Wall -std=f95

#OPTS03 = -O3
#OPTS03 = -O3 -fdefault-real-8 -fdefault-double-8
#OPTS03 = -O3 -fcheck=bounds -fdefault-real-8 -fdefault-double-8 -Wall -std=f2003

#CCOPTS = -O3 -Wall -std=c99

#CUDA
#CUDAHOME = /usr/local/cuda
#LOPTS = -L/usr/lib/gcc/i686-apple-darwin9/4.0.1 -lstdc++
#LOPTS = -L/usr/lib/gcc/i686-apple-darwin9/4.2.1 -lstdc++
#LOPTS = -L/usr/lib/gcc/i686-apple-darwin10/4.2.1 -lstdc++
#LOPTS = -L/usr/lib/gcc/i686-apple-darwin10/4.2.1/x86_64 -lstdc++
#CULIBS = $(CARBON) -lSystemStubs -L$(CUDAHOME)/lib -lcuda -lcudart \
#-lcufft

#NVOPTS = -O3 -Xptxas=-v
#NVOPTS = -O3 -arch sm_20
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
#CUDAHOME = /usr/local/cuda
#LOPTS = -L/usr/lib/gcc/i686-apple-darwin9/4.0.1 -lstdc++
#LOPTS = -L/usr/lib/gcc/i686-apple-darwin9/4.2.1 -lstdc++
#LOPTS = -L/usr/lib/gcc/i686-apple-darwin10/4.2.1 -lstdc++
#LOPTS = -L/usr/lib/gcc/i686-apple-darwin10/4.2.1/x86_64 -lstdc++
#CULIBS = $(CARBON) -lSystemStubs -L$(CUDAHOME)/lib -lcuda -lcudart \
#-lcufft

#NVOPTS = -O3 -arch sm_20
#NVOPTS = -O3 -m64 -arch sm_20

# Makefile Intel compiler with Linux

#FC90 = ifort
#FC03 = ifort
#CC = icc

#OPTS90 = -O3
#OPTS90 = -O3 -r8
#OPTS90 = -O3 -CB -r8 -warn all -std90

#FF03 = -Tf
#OPTS03 = -O3 -free -std03
#OPTS03 = -O3 -r8 -free -std03
#OPTS03 = -O3 -CB -r8 -warn all -free -std03

#CCOPTS = -O3 -Wall -std=c99

#CUDA
#CUDAHOME = /u/local/cuda/5.5
#LOPTS =
#CULIBS = -L$(CUDAHOME)/lib64 -lcudart -lcufft
#NVOPTS = -O3 -arch sm_20

# Makefile gfortran compiler with Linux

FC90 = gfortran
FC03 = gfortran
CC = gcc

OPTS90 = -O3
#OPTS90 = -O3 -fdefault-real-8 -fdefault-double-8
#OPTS90 = -O3 -fbounds-check -fdefault-real-8 -Wall -std=f95

OPTS03 = -O3
#OPTS03 = -O3 -fdefault-real-8
#OPTS03 = -O3 -fbounds-check -fdefault-real-8 -fdefault-double-8 -Wall -std=f2003

CCOPTS = -O3 -Wall -std=c99

#CUDA
CUDAHOME = /u/local/cuda/5.5
LOPTS =
CULIBS = -L$(CUDAHOME)/lib64 -lcuda -lcudart -lcufft
#CULIBS = -L$(CUDAHOME)/lib -lcuda -lcudart -lcufft
#NVOPTS = -O3 -Xptxas=-v
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

#CCOPTS = -O3 -Wall -std=c99

#CUDA
#CUDAHOME = /u/local/cuda/5.5
#LOPTS =
#CULIBS = -L$(CUDAHOME)/lib64 -lcudart -lcufft
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

all : fgpupic2 cgpupic2 f03gpupic2

fgpupic2 : fgpupic2.o gpulib2.o gpupush2.o gpufft2.o fpush2.o \
           dtimer.o
	$(FC90) $(OPTS90) $(LOPTS) -o fgpupic2 fgpupic2.o \
	gpulib2.o gpupush2.o gpufft2.o fpush2.o gpulib2_h.o \
	gpupush2_h.o gpufft2_h.o push2_h.o dtimer.o $(CULIBS)

cgpupic2 : cgpupic2.o gpulib2.o gpupush2.o gpufft2.o cpush2.o \
           dtimer.o
	$(CC) $(CCOPTS) $(LOPTS) -o cgpupic2 cgpupic2.o \
	gpulib2.o gpupush2.o gpufft2.o cpush2.o dtimer.o $(CULIBS) -lm

f03gpupic2 : f03gpupic2.o gpulib2.o gpupush2.o gpufft2.o fpush2.o \
             dtimer.o
	$(FC03) $(OPTS03) $(LOPTS) -o f03gpupic2 f03gpupic2.o \
	gpulib2.o gpupush2.o gpufft2.o fpush2.o gpulib2_c.o \
	gpupush2_c.o gpufft2_c.o push2_h.o dtimer.o $(CULIBS)

# Compilation rules

dtimer.o : dtimer.c
	$(CC) $(CCOPTS) -c dtimer.c

gpulib2.o : gpulib2.cu
	nvcc $(NVOPTS) -c gpulib2.cu -I$(CUDAHOME)/include

gpupush2.o : gpupush2.cu
	nvcc $(NVOPTS) -c gpupush2.cu -I$(CUDAHOME)/include

gpufft2.o : gpufft2.cu
	nvcc $(NVOPTS) -c gpufft2.cu -I$(CUDAHOME)/include

gpulib2_h.o : gpulib2_h.f90
	$(FC90) $(OPTS90) -c gpulib2_h.f90

gpupush2_h.o : gpupush2_h.f90
	$(FC90) $(OPTS90) -c gpupush2_h.f90

gpufft2_h.o : gpufft2_h.f90
	$(FC90) $(OPTS90) -c gpufft2_h.f90

gpulib2_c.o : gpulib2_c.f03
	$(FC03) $(OPTS03) -c $(FF03) gpulib2_c.f03

gpupush2_c.o : gpupush2_c.f03
	$(FC03) $(OPTS03) -c $(FF03) gpupush2_c.f03

gpufft2_c.o : gpufft2_c.f03
	$(FC03) $(OPTS03) -c $(FF03) gpufft2_c.f03

fpush2.o : push2.f
	$(FC90) $(OPTS90) -o fpush2.o -c push2.f

push2_h.o : push2_h.f90
	$(FC90) $(OPTS90) -o push2_h.o -c push2_h.f90

cpush2.o : push2.c
	$(CC) $(CCOPTS) -o cpush2.o -c push2.c

fgpupic2.o : gpupic2.f90 gpulib2_h.o gpupush2_h.o gpufft2_h.o push2_h.o
	$(FC90) $(OPTS90) -o fgpupic2.o -c gpupic2.f90

cgpupic2.o : gpupic2.c
	$(CC) $(CCOPTS) -o cgpupic2.o -c gpupic2.c

f03gpupic2.o : gpupic2.f03 gpulib2_c.o gpupush2_c.o gpufft2_c.o push2_h.o
	$(FC03) $(OPTS03) -o f03gpupic2.o -c $(FF03) gpupic2.f03

clean :
	rm -f *.o *.mod

clobber: clean
	rm -f fgpupic2 cgpupic2 f03gpupic2
