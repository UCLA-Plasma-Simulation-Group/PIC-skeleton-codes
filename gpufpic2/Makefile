#Makefile for 2D Electrostatic GPU PIC codes

# Makefile Cuda Fortran PGI compiler with Linux

FCU90 = pgf90 -Mcuda
#FCU90 = pgf90 -Mcuda=emu
FC90 = pgf90
CC = gcc

FCUOPTS90 = -O3
#FCUOPTS90 = -O3 -r8
#FCUOPTS90 = -O3 -Mbounds -Mstandard
#FCUOPTS90 = -O3 -r8 -Mbounds -Mstandard

OPTS90 = -O3
#OPTS90 = -O3 -r8

CCOPTS = -O3 -Wall -std=c99

FNVOPTS = -O3 -ta=nvidia:cc20
#FNVOPTS = -O3 -r8 -ta=nvidia:cc20
#FNVOPTS = -O3 -ta=nvidia:cc20 -Mbounds -Mstandard
#FNVOPTS = -O3 -r8 -ta=nvidia:cc20 -Mbounds -Mstandard
NVOPTS = -O3 -arch sm_20

LEGACY = -Mnomain

CULIBS = -L/u/local/cuda/5.0/lib64 -lcudart -lcufft

#

# Linkage rules

all : fgpufpic2

special : cgpupic2_f

fgpufpic2 : fgpufpic2.o fpush2.o dtimer.o
	$(FCU90) $(FCUOPTS90) $(LOPTS) -o fgpufpic2 fgpufpic2.o \
	fgpupush2.o fgpulib2.o fgpufft2.o fpush2.o push2_h.o dtimer.o $(CULIBS)

cgpupic2_f : cgpupic2_f.o fgpulib2_c.o fgpupush2_c.o fgpufft2_c.o \
             fgpulib2.o fgpupush2.o fgpufft2.o cpush2.o gpulib2s.o dtimer.o
	$(FCU90) $(FCUOPTS90) $(LEGACY) $(LOPTS) -o cgpupic2_f cgpupic2_f.o \
    fgpulib2_c.o fgpupush2_c.o fgpufft2_c.o fgpulib2.o fgpupush2.o \
    fgpufft2.o cpush2.o gpulib2s.o dtimer.o $(CULIBS) -lm

# Compilation rules

dtimer.o : dtimer.c
	$(CC) $(CCOPTS) -c dtimer.c

fgpulib2.o : gpulib2.cuf
	$(FCU90) $(FNVOPTS) -o fgpulib2.o -c gpulib2.cuf -I/usr/local/cuda/include

fgpupush2.o : gpupush2.cuf fgpulib2.o
	$(FCU90) $(FNVOPTS) -o fgpupush2.o -c gpupush2.cuf -I/usr/local/cuda/include

fgpufft2.o : gpufft2.cuf fgpupush2.o
	$(FCU90) $(FNVOPTS) -o fgpufft2.o -c gpufft2.cuf -I/usr/local/cuda/include

fpush2.o : push2.f
	$(FC90) $(OPTS90) -o fpush2.o -c push2.f

push2_h.o : push2_h.f90
	$(FC90) $(OPTS90) -o push2_h.o -c push2_h.f90

cpush2.o : push2.c
	$(CC) $(CCOPTS) -o cpush2.o -c push2.c

fgpufpic2.o : gpupic2.cuf push2_h.o fgpupush2.o fgpulib2.o fgpufft2.o
	$(FCU90) $(FCUOPTS90) -o fgpufpic2.o -c gpupic2.cuf

fgpulib2_c.o : gpulib2_c.cuf fgpulib2.o
	$(FCU90) $(OPTS90) -o fgpulib2_c.o -c gpulib2_c.cuf

fgpupush2_c.o : gpupush2_c.cuf fgpupush2.o
	$(FCU90) $(OPTS90) -o fgpupush2_c.o -c gpupush2_c.cuf

fgpufft2_c.o : gpufft2_c.cuf fgpufft2.o
	$(FCU90) $(OPTS90) -o fgpufft2_c.o -c gpufft2_c.cuf

gpulib2s.o : gpulib2s.cu
	nvcc $(NVOPTS) -c gpulib2s.cu -I/usr/local/cuda/include

cgpupic2_f.o : gpupic2_f.c
	$(CC) $(CCOPTS) -o cgpupic2_f.o -c gpupic2_f.c

clean :
	rm -f *.o *.mod

clobber: clean
	rm -f fgpufpic2 cgpupic2_f
