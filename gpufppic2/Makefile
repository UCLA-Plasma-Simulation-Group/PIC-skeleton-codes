#Makefile for 2D Electrostatic GPU-MPI PIC codes

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

all : fgpufppic2

special : cgpuppic2_f

fgpufppic2 : fgpufppic2.o fppush2.o dtimer.o
	$(MPIFC) -Mcuda $(OPTS90) $(LOPTS) -o fgpufppic2 fgpufppic2.o \
        fgpulib2.o fgpuppush2.o fgpupfft2.o fgpplib2.o fpplib2.o \
        fppush2.o ppush2_h.o dtimer_c.o dtimer.o $(CULIBS) $(LMPILIBS)

cgpuppic2_f : cgpuppic2_f.o fgpulib2_c.o fgpuppush2_c.o fgpupfft2_c.o \
              fgpplib2_c.o fpplib2_c.o cppush2.o gpulib2s.o dtimer.o
	$(MPIFC) -Mcuda $(OPTS90) $(LOPTS) $(LEGACY) -o cgpuppic2_f cgpuppic2_f.o \
    fgpulib2_c.o fgpuppush2_c.o fgpupfft2_c.o fgpplib2_c.o fpplib2_c.o \
    fgpulib2.o fgpuppush2.o fgpupfft2.o fgpplib2.o fpplib2.o cppush2.o \
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

fgpuppush2.o : gpuppush2.cuf fgpulib2.o
	$(FCU90) $(FNVOPTS) -o fgpuppush2.o -c gpuppush2.cuf -I/usr/local/cuda/include

fgpupfft2.o : gpupfft2.cuf fgpuppush2.o
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

fppush2.o : ppush2.f
	$(FCU90) $(OPTS90) -o fppush2.o -c ppush2.f

ppush2_h.o : ppush2_h.f90
	$(FCU90) $(OPTS90) -o ppush2_h.o -c ppush2_h.f90

cppush2.o : ppush2.c
	$(CC) $(CCOPTS) -o cppush2.o -c ppush2.c

fgpufppic2.o : gpuppic2.cuf fgpplib2.o ppush2_h.o
	$(FCU90) $(FNVOPTS) -o fgpufppic2.o -c gpuppic2.cuf

fgpulib2_c.o : gpulib2_c.cuf fgpulib2.o
	$(FCU90) $(FNVOPTS) -o fgpulib2_c.o -c gpulib2_c.cuf

fgpuppush2_c.o : gpuppush2_c.cuf fgpuppush2.o
	$(FCU90) $(FNVOPTS) -o fgpuppush2_c.o -c gpuppush2_c.cuf

fgpupfft2_c.o : gpupfft2_c.cuf fgpupfft2.o
	$(FCU90) $(FNVOPTS) -o fgpupfft2_c.o -c gpupfft2_c.cuf

fgpplib2_c.o : gpplib2_c.cuf fgpplib2.o
	$(FCU90) $(FNVOPTS) -o fgpplib2_c.o -c gpplib2_c.cuf

gpulib2s.o : gpulib2s.cu
	nvcc $(NVOPTS) -c gpulib2s.cu -I/usr/local/cuda/include

cgpuppic2_f.o : gpuppic2_f.c
	$(CC) $(CCOPTS) -o cgpuppic2_f.o -c gpuppic2_f.c

clean :
	rm -f *.o *.mod

clobber: clean
	rm -f fgpufppic2 cgpuppic2_f
