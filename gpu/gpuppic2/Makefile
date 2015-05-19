#Makefile for 2D Electrostatic GPU-MPI PIC codes

# Makefile gfortran compiler with MacOS X

#CARBON = /System/Library/Frameworks/Carbon.framework/Carbon

#OpenMPI
#MPIFC = mpif90
#MPICC = mpicc

#CC = gcc

#OPTS90 = -O3
#OPTS90 = -O3 -fdefault-real-8 -fdefault-double-8
#OPTS90 = -O3 -fcheck=bounds -fdefault-real-8 -fdefault-double-8 -Wall -std=f95

#OPTS03 = -O3
#OPTS03 = -O3 -fdefault-real-8 -fdefault-double-8
#OPTS03 = -O3 -fcheck=bounds -fdefault-real-8 -fdefault-double-8 -Wall -std=f2003

#CCOPTS = -O3 -std=c99
#CCOPTS = -O3 -Wall -std=c99

#CUDA
#LOPTS = -L/usr/lib/gcc/i686-apple-darwin9/4.0.1 -lstdc++
#LOPTS = -L/usr/lib/gcc/i686-apple-darwin9/4.2.1 -lstdc++
#LOPTS = -L/usr/lib/gcc/i686-apple-darwin10/4.2.1 -lstdc++
#LOPTS = -L/usr/lib/gcc/i686-apple-darwin10/4.2.1/x86_64 -lstdc++
#CULIBS = $(CARBON) -lSystemStubs -L/usr/local/cuda/lib -lcuda -lcudart \
#-lcufft

#NVOPTS = -O3 -arch sm_20
#NVOPTS = -O3 -m64 -arch sm_20

# Makefile Intel compiler with Mac OS X

#CARBON = /System/Library/Frameworks/Carbon.framework/Carbon

# OpenMPI
#MPIFC = mpif90
#MPICC = mpicc

#CC = gcc

#OPTS90 = -O3
#OPTS90 = -O3 -r8
#OPTS90 = -O3 -CB -r8 -warn all -std90

#FF03 = -Tf
#OPTS03 = -O3
#OPTS03 = -O3 -r8
#OPTS03 = -O3 -CB -r8 -warn all -std03

#CCOPTS = -O3 -std=c99
#CCOPTS = -O3 -Wall -std=c99

#CUDA
#LOPTS = -L/usr/lib/gcc/i686-apple-darwin9/4.0.1 -lstdc++
#LOPTS = -L/usr/lib/gcc/i686-apple-darwin9/4.2.1 -lstdc++
#LOPTS = -L/usr/lib/gcc/i686-apple-darwin10/4.2.1 -lstdc++
#LOPTS = -L/usr/lib/gcc/i686-apple-darwin10/4.2.1/x86_64 -lstdc++
#CULIBS = $(CARBON) -lSystemStubs -L/usr/local/cuda/lib -lcuda -lcudart \
#-lcufft

#NVOPTS = -O3 -arch sm_20
#NVOPTS = -O3 -m64 -arch sm_20

# Makefile Intel compiler with Linux

# MPI
MPIFC = mpiifort
MPICC = mpiicc

CC = icc

OPTS90 = -O3
#OPTS90 = -O3 -r8
#OPTS90 = -O3 -CB -r8 -warn all -std90

FF03 = -Tf
OPTS03 = -O3 -free -std03
#OPTS03 = -O3 -r8 -free -std03
#OPTS03 = -O3 -CB -r8 -warn all -free -std03

CCOPTS = -O3 -Wall -std=c99

#CUDA
LOPTS =
#CULIBS = -L/usr/local/cuda/lib64 -lcuda -lcudart -lcufft
#CULIBS = -L/usr/local/cuda/lib -lcuda -lcudart -lcufft
#CULIBS = -L/u/local/cuda/4.0.17/lib64 -lcuda -lcudart -lcufft
#CULIBS = -L/u/local/cuda/4.0.17/lib64 -lcudart -lcufft
CULIBS = -L/u/local/cuda/5.0/lib64 -lcudart -lcufft
#CULIBS = -lcudart -lcufft

NVOPTS = -O3 -arch sm_20

# Makefile gfortran compiler with Linux

# OpenMPI
#MPIFC = mpif90
#MPICC = mpicc

#CC = gcc

#OPTS90 = -O3
#OPTS90 = -O3 -fdefault-real-8 -fdefault-double-8
#OPTS90 = -O3 -fbounds-check -fdefault-real-8 -fdefault-double-8 -Wall -std=f95

#OPTS03 = -O3
#OPTS03 = -O3 -fdefault-real-8 -fdefault-double-8
#OPTS03 = -O3 -fbounds-check -fdefault-real-8 -fdefault-double-8 -Wall -std=f2003

#CCOPTS = -O3 -Wall -std=c99

#CUDA
#LOPTS =
#CULIBS = -L/usr/local/cuda/lib64 -lcuda -lcudart -lcufft
#CULIBS = -L/usr/local/cuda/lib -lcuda -lcudart -lcufft
#CULIBS = -L/u/local/cuda/4.0.17/lib64 -lcuda -lcudart -lcufft
#CULIBS = -L/u/local/cuda/5.0/lib64 -lcudart -lcufft
#CULIBS = -lcudart -lcufft

#NVOPTS = -O3 -arch sm_20

# Makefile PGI compiler with Linux

# OpenMPI
#MPIFC = mpipgf90
#MPICC = mpicc

#CC = gcc

#OPTS90 = -O3
#OPTS90 = -O3 -r8
#OPTS90 = -O3 -Mbounds -r8 -Mstandard

#OPTS03 = -O3
#OPTS03 = -O3 -r8
#OPTS03 = -O3 -Mbounds -r8 -Mstandard

#CCOPTS = -O3 -Wall -std=c99

#CUDA
#LOPTS =
#CULIBS = -L/usr/local/cuda/lib64 -lcuda -lcudart -lcufft
#CULIBS = -L/usr/local/cuda/lib -lcuda -lcudart -lcufft
#CULIBS = -L/u/local/cuda/4.0.17/lib64 -lcuda -lcudart -lcufft
#CULIBS = -L/u/local/cuda/5.0/lib64 -lcudart -lcufft
#CULIBS = -lcudart -lcufft

#NVOPTS = -O3 -arch sm_20

# Makefile Cray compiler with Linux
 
#MPICC = nvcc

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

all : fgpuppic2 cgpuppic2 f03gpuppic2

fgpuppic2 : fgpuppic2.o gpulib2.o gpuppush2.o gpupfft2.o fppush2.o \
            fpplib2.o fgpplib2.o dtimer.o
	$(MPIFC) $(OPTS90) $(LOPTS) -o fgpuppic2 fgpuppic2.o \
    gpulib2.o gpuppush2.o gpupfft2.o fppush2.o fpplib2.o fgpplib2.o \
    gpulib2_h.o gpuppush2_h.o gpupfft2_h.o ppush2_h.o dtimer.o \
    $(CULIBS)

cgpuppic2 : cgpuppic2.o gpulib2.o gpuppush2.o gpupfft2.o cppush2.o \
            cpplib2.o cgpplib2.o dtimer.o
	$(MPICC) $(CCOPTS) $(LOPTS) -o cgpuppic2 cgpuppic2.o \
	gpulib2.o gpuppush2.o gpupfft2.o cppush2.o cpplib2.o cgpplib2.o \
    dtimer.o $(CULIBS)

f03gpuppic2 : f03gpuppic2.o gpulib2.o gpuppush2.o gpupfft2.o fppush2.o \
              fpplib2.o f03gpplib2.o dtimer.o
	$(MPIFC) $(OPTS03) $(LOPTS) -o f03gpuppic2 f03gpuppic2.o \
    gpulib2.o gpuppush2.o gpupfft2.o fppush2.o fpplib2.o f03gpplib2.o \
    gpulib2_c.o gpuppush2_c.o gpupfft2_c.o ppush2_h.o dtimer_c.o dtimer.o \
    $(CULIBS)

# Compilation rules

dtimer.o : dtimer.c
	$(CC) $(CCOPTS) -c dtimer.c

dtimer_c.o : dtimer_c.f03
	$(MPIFC) $(OPTS03) -c $(FF03) dtimer_c.f03

fgpplib2.o : gpplib2.f90 gpulib2_h.o gpuppush2_h.o gpupfft2_h.o fpplib2.o
	$(MPIFC) $(OPTS90) -o fgpplib2.o -c gpplib2.f90

cgpplib2.o : gpplib2.c
	$(MPICC) $(CCOPTS) -o cgpplib2.o -c gpplib2.c

f03gpplib2.o : gpplib2.f03 gpulib2_c.o gpuppush2_c.o gpupfft2_c.o fpplib2.o \
               dtimer_c.o
	$(MPIFC) $(OPTS03) -o f03gpplib2.o -c $(FF03) gpplib2.f03

gpulib2.o : gpulib2.cu
	nvcc $(NVOPTS) -c gpulib2.cu -I/usr/local/cuda/include

gpuppush2.o : gpuppush2.cu
	nvcc $(NVOPTS) -c gpuppush2.cu -I/usr/local/cuda/include

gpupfft2.o : gpupfft2.cu
	nvcc $(NVOPTS) -c gpupfft2.cu -I/usr/local/cuda/include

fpplib2.o : pplib2.f90
	$(MPIFC) $(OPTS90) -o fpplib2.o -c pplib2.f90

cpplib2.o : pplib2.c
	$(MPICC) $(CCOPTS) -o cpplib2.o -c pplib2.c

gpulib2_h.o : gpulib2_h.f90
	$(MPIFC) $(OPTS90) -c gpulib2_h.f90

gpuppush2_h.o : gpuppush2_h.f90
	$(MPIFC) $(OPTS90) -c gpuppush2_h.f90

gpupfft2_h.o : gpupfft2_h.f90
	$(MPIFC) $(OPTS90) -c gpupfft2_h.f90

gpulib2_c.o : gpulib2_c.f03
	$(MPIFC) $(OPTS03) -c $(FF03) gpulib2_c.f03

gpuppush2_c.o : gpuppush2_c.f03
	$(MPIFC) $(OPTS03) -c $(FF03) gpuppush2_c.f03

gpupfft2_c.o : gpupfft2_c.f03
	$(MPIFC) $(OPTS03) -c $(FF03) gpupfft2_c.f03

fppush2.o : ppush2.f
	$(MPIFC) $(OPTS90) -o fppush2.o -c ppush2.f

ppush2_h.o : ppush2_h.f90
	$(MPIFC) $(OPTS90) -o ppush2_h.o -c ppush2_h.f90

cppush2.o : ppush2.c
	$(MPICC) $(CCOPTS) -o cppush2.o -c ppush2.c

fgpuppic2.o : gpuppic2.f90 fgpplib2.o gpulib2_h.o gpuppush2_h.o gpupfft2_h.o \
              ppush2_h.o fpplib2.o 
	$(MPIFC) $(OPTS90) -o fgpuppic2.o -c gpuppic2.f90

cgpuppic2.o : gpuppic2.c
	$(MPICC) $(CCOPTS) -o cgpuppic2.o -c gpuppic2.c

f03gpuppic2.o : gpuppic2.f03 f03gpplib2.o gpulib2_c.o gpuppush2_c.o \
                gpupfft2_c.o ppush2_h.o fpplib2.o
	$(MPIFC) $(OPTS03) -o f03gpuppic2.o -c $(FF03) gpuppic2.f03

clean :
	rm -f *.o *.mod

clobber: clean
	rm -f fgpuppic2 cgpuppic2 f03gpuppic2
