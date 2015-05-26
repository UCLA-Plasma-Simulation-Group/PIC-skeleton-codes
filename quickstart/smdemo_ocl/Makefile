# Makefile for OpenCL test program

# Makefile gfortran compiler with MacOS X

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
OCLHOME = /u/local/cuda/current
CLI = -I$(OCLHOME)/include
#CLI = -I/usr/local/atisdk/include
CLLIBS = -lOpenCL
#OCLHOME = /usr/local/atisdk/lib/x86_64
#CLLIBS = -L$(OCLHOME) -lOpenCL

#

CLOBJS = gpulib_cl.o dtimer.o

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
gpulib_cl.o : gpulib_cl.c
	$(CC) $(CCOPTS) -c gpulib_cl.c $(CLI)

ocllib_c.o : ocllib_c.f03
	$(FC90) $(OPTS90) -o ocllib_c.o -c ocllib_c.f03

cgputest_cl.o : cgputest_cl.c
	$(CC) $(CCOPTS) -c cgputest_cl.c $(CLI)

fgputest_cl.o : fgputest_cl.f03 ocllib_c.o
	$(FC90) $(OPTS90) -c fgputest_cl.f03

clean:
	rm -f *.o *.mod

clobber: clean
	rm -f *gputest_cl
