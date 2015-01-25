#Makefile for 2-1/2D Electromagnetic GPU-MPI PIC codes

# Makefile CUDA Fortran PGI compiler with Linux

# OpenMPI
MPIFC = pgf90 -I/u/local/compilers/intel-cs/2013.0.028/mpi/intel64/include

FCU90 = pgf90
#FCU90 = pgf90 -Mcuda=emu
CC = gcc

OPTS90 = -O3
#OPTS90 = -O3 -r8
#OPTS90 = -O3 -Mbounds -r8 -Mstandard

CCOPTS = -O3 -Wall -std=c99
LEGACY = -Mnomain

#CUDA
LOPTS =
#CULIBS = -L/usr/local/cuda/lib64 -lcuda -lcudart -lcufft
#CULIBS = -L/usr/local/cuda/lib -lcuda -lcudart -lcufft
#CULIBS = -L/u/local/cuda/4.0.17/lib64 -lcuda -lcudart -lcufft
#CULIBS = -L/u/local/cuda/4.0.17/lib64 -lcudart -lcufft
CULIBS = -lcudart -lcufft
LMPILIBS= -L/u/local/compilers/intel-cs/2013.0.028/mpi/intel64/lib \
-lmpigf -lmpi -lmpigi -ldl -lrt -lpthread

FNVOPTS = -O3 -ta=nvidia:cc20
#FNVOPTS = -O3 -r8 -ta=nvidia:cc20
#FNVOPTS = -O3 -ta=nvidia:cc20 -Mbounds -Mstandard
#FNVOPTS = -O3 -r8 -ta=nvidia:cc20 -Mbounds -Mstandard
NVOPTS = -O3 -arch sm_20

#

# Linkage rules

all: fgpufpbpic2

special : cgpupbpic2_f

fgpufpbpic2 : fgpufpbpic2.o fpbpush2.o dtimer.o
	$(MPIFC) -Mcuda $(OPTS90) $(LOPTS) -o fgpufpbpic2 fgpufpbpic2.o \
        fgpulib2.o fgpupbpush2.o fgpupfft2.o fgpplib2.o fpplib2.o \
        fpbpush2.o pbpush2_h.o dtimer_c.o dtimer.o $(CULIBS) $(LMPILIBS)

cgpupbpic2_f : cgpupbpic2_f.o fgpulib2_c.o fgpupbpush2_c.o fgpupfft2_c.o \
               fgpplib2_c.o fpplib2_c.o cpbpush2.o gpulib2s.o dtimer.o
	$(MPIFC) -Mcuda $(OPTS90) $(LOPTS) $(LEGACY) -o cgpupbpic2_f cgpupbpic2_f.o \
    fgpulib2_c.o fgpupbpush2_c.o fgpupfft2_c.o fgpplib2_c.o fpplib2_c.o \
    fgpulib2.o fgpupbpush2.o fgpupfft2.o fgpplib2.o fpplib2.o cpbpush2.o \
    gpulib2s.o dtimer_c.o dtimer.o $(CULIBS) $(LMPILIBS) -lm

# Compilation rules

dtimer.o : dtimer.c
	$(CC) $(CCOPTS) -c dtimer.c

dtimer_c.o : dtimer_c.f03
	$(MPIFC) $(OPTS03) -c $(FF03) dtimer_c.f03

fgpplib2.o : gpplib2.cuf fgpupfft2.o fpplib2.o dtimer_c.o
	$(FCU90) $(FNVOPTS) -o fgpplib2.o -c gpplib2.cuf

fgpulib2.o : gpulib2.cuf
	$(FCU90) $(FNVOPTS) -o fgpulib2.o -c gpulib2.cuf -I/usr/local/cuda/include

fgpupbpush2.o : gpupbpush2.cuf fgpulib2.o
	$(FCU90) $(FNVOPTS) -o fgpupbpush2.o -c gpupbpush2.cuf -I/usr/local/cuda/include

fgpupfft2.o : gpupfft2.cuf fgpupbpush2.o
	$(FCU90) $(FNVOPTS) -o fgpupfft2.o -c gpupfft2.cuf -I/usr/local/cuda/include

# Version using Fortran77 pplib2.f
#fpplib2.o : pplib2.f
#	$(MPIFC) $(OPTS90) -o fpplib2.o -c pplib2.f
#
#pplib2_h.o : pplib2_h.f90
#	$(FCU90) $(OPTS90) -o pplib2_h.o -c pplib2_h.f90

# Version using Fortran90 pplib2.f90
fpplib2.o : pplib2.f90
	$(MPIFC) $(OPTS90) -o fpplib2.o -c pplib2.f90

fpplib2_c.o : pplib2_c.f03 fpplib2.o
	$(FCU90) $(OPTS90) -o fpplib2_c.o -c pplib2_c.f03

fpbpush2.o : pbpush2.f
	$(FCU90) $(OPTS90) -o fpbpush2.o -c pbpush2.f

pbpush2_h.o : pbpush2_h.f90
	$(FCU90) $(OPTS90) -o pbpush2_h.o -c pbpush2_h.f90

cpbpush2.o : pbpush2.c
	$(CC) $(CCOPTS) -o cpbpush2.o -c pbpush2.c

fgpufpbpic2.o : gpupbpic2.cuf fgpplib2.o pbpush2_h.o
	$(FCU90) $(FNVOPTS) -o fgpufpbpic2.o -c gpupbpic2.cuf

fgpulib2_c.o : gpulib2_c.cuf fgpulib2.o
	$(FCU90) $(FNVOPTS) -o fgpulib2_c.o -c gpulib2_c.cuf

fgpupbpush2_c.o : gpupbpush2_c.cuf fgpupbpush2.o
	$(FCU90) $(FNVOPTS) -o fgpupbpush2_c.o -c gpupbpush2_c.cuf

fgpupfft2_c.o : gpupfft2_c.cuf fgpupfft2.o
	$(FCU90) $(FNVOPTS) -o fgpupfft2_c.o -c gpupfft2_c.cuf

fgpplib2_c.o : gpplib2_c.cuf fgpplib2.o
	$(FCU90) $(FNVOPTS) -o fgpplib2_c.o -c gpplib2_c.cuf

gpulib2s.o : gpulib2s.cu
	nvcc $(NVOPTS) -c gpulib2s.cu -I/usr/local/cuda/include

cgpupbpic2_f.o : gpupbpic2_f.c
	$(CC) $(CCOPTS) -o cgpupbpic2_f.o -c gpupbpic2_f.c

clean :
	rm -f *.o *.mod

clobber: clean
	rm -f fgpufpbpic2 cgpupbpic2_f
