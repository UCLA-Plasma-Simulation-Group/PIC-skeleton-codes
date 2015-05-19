# Makefile for OpenCL test program

# Makefile gfortran compiler with MacOS X

CARBON = /System/Library/Frameworks/Carbon.framework/Carbon

#FC90 = gfortran
#CC = gcc

#OPTS90 = -O3
#CCOPTS = -O3

#OpenCL
#CLI =
#CLLIBS = -framework OpenCL

# Makefile gfortran compiler with Linux

FC90 = gfortran
CC = gcc

OPTS90 = -O3
CCOPTS = -O3

#OpenCL
CLI = -I/u/local/cuda/current/include
#CLI = -I/usr/local/atisdk/include
CLLIBS = -lOpenCL
#CLLIBS = -L/usr/local/atisdk/lib/x86_64 -lOpenCL

#

CLOBJS = gpumain_cl.o dtimer.o

# Linkage rules

all: opencl

opencl: cgputest_cl fgputest_cl

#OPENCL
cgputest_cl : cgputest_cl.o $(CLOBJS)
	$(CC) $(CCOPTS) $(LOPTS) -o cgputest_cl cgputest_cl.o \
      $(CLOBJS) $(CLLIBS)

fgputest_cl : fgputest_cl.o $(CLOBJS)
	$(FC90) $(OPTS90) $(LOPTS) -o fgputest_cl fgputest_cl.o \
      $(CLOBJS) $(CLLIBS)

# Compilation rules

dtimer.o : dtimer.c
	$(CC) $(CCOPTS) -c dtimer.c

#OPENCL
gpumain_cl.o : gpumain_cl.c
	$(CC) $(CCOPTS) -c gpumain_cl.c $(CLI)

cgputest_cl.o : cgputest_cl.c
	$(CC) $(CCOPTS) -c cgputest_cl.c $(CLI)

fgputest_cl.o : fgputest_cl.f90
	$(FC90) $(OPTS90) -c fgputest_cl.f90

clean:
	rm -f *.o *.mod

clobber: clean
	rm -f *gputest_cl
