#Makefile for 2-1/2D Electromagnetic GPU PIC codes

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

all : fgpufbpic2

special : cgpubpic2_f

fgpufbpic2 : fgpufbpic2.o fbpush2.o dtimer.o
	$(FCU90) $(FCUOPTS90) $(LOPTS) -o fgpufbpic2 fgpufbpic2.o \
	fgpubpush2.o fgpulib2.o fgpufft2.o fbpush2.o bpush2_h.o dtimer.o $(CULIBS)

cgpubpic2_f : cgpubpic2_f.o fgpulib2_c.o fgpubpush2_c.o fgpufft2_c.o \
              fgpulib2.o fgpubpush2.o fgpufft2.o cbpush2.o gpulib2s.o dtimer.o
	$(FCU90) $(FCUOPTS90) $(LEGACY) $(LOPTS) -o cgpubpic2_f cgpubpic2_f.o \
    fgpulib2_c.o fgpubpush2_c.o fgpufft2_c.o fgpulib2.o fgpubpush2.o \
    fgpufft2.o cbpush2.o gpulib2s.o dtimer.o $(CULIBS) -lm

# Compilation rules

dtimer.o : dtimer.c
	$(CC) $(CCOPTS) -c dtimer.c

fgpulib2.o : gpulib2.cuf
	$(FCU90) $(FNVOPTS) -o fgpulib2.o -c gpulib2.cuf -I/usr/local/cuda/include

fgpubpush2.o : gpubpush2.cuf fgpulib2.o
	$(FCU90) $(FNVOPTS) -o fgpubpush2.o -c gpubpush2.cuf -I/usr/local/cuda/include

fgpufft2.o : gpufft2.cuf fgpubpush2.o
	$(FCU90) $(FNVOPTS) -o fgpufft2.o -c gpufft2.cuf -I/usr/local/cuda/include

fbpush2.o : bpush2.f
	$(FC90) $(OPTS90) -o fbpush2.o -c bpush2.f

bpush2_h.o : bpush2_h.f90
	$(FC90) $(OPTS90) -o bpush2_h.o -c bpush2_h.f90

cbpush2.o : bpush2.c
	$(CC) $(CCOPTS) -o cbpush2.o -c bpush2.c

fgpufbpic2.o : gpubpic2.cuf bpush2_h.o fgpubpush2.o fgpulib2.o fgpufft2.o
	$(FCU90) $(FCUOPTS90) -o fgpufbpic2.o -c gpubpic2.cuf

fgpulib2_c.o : gpulib2_c.cuf fgpulib2.o
	$(FCU90) $(OPTS90) -o fgpulib2_c.o -c gpulib2_c.cuf

fgpubpush2_c.o : gpubpush2_c.cuf fgpubpush2.o
	$(FCU90) $(OPTS90) -o fgpubpush2_c.o -c gpubpush2_c.cuf

fgpufft2_c.o : gpufft2_c.cuf fgpufft2.o
	$(FCU90) $(OPTS90) -o fgpufft2_c.o -c gpufft2_c.cuf

gpulib2s.o : gpulib2s.cu
	nvcc $(NVOPTS) -c gpulib2s.cu -I/usr/local/cuda/include

cgpubpic2_f.o : gpubpic2_f.c
	$(CC) $(CCOPTS) -o cgpubpic2_f.o -c gpubpic2_f.c

clean :
	rm -f *.o *.mod

clobber: clean
	rm -f fgpufbpic2 cgpubpic2_f
