# Makefile for SSE test program

# Makefile gfortran compiler with Linux

FC90 = gfortran
CC = gcc

OPTS90 = -O3 -fno-tree-vectorize
CCOPTS = -O3 -fno-tree-vectorize

# Makefile Intel compiler with Linux

#FC90 = ifort
#CC = icc

#FF03 = -Tf
#OPTS90 = -O3 -no-vec
#CCOPTS = -O3 -no-vec

#

SSEOBJS = sselib.o dtimer.o

# Linkage rules

all: sse

sse: cssetest fssetest

#SSE
cssetest : cssetest.o $(SSEOBJS)
	$(CC) $(CCOPTS) $(LOPTS) -o cssetest cssetest.o \
      $(SSEOBJS)

fssetest : fssetest.o sselib_c.o $(SSEOBJS)
	$(FC90) $(OPTS90) $(LOPTS) -o fssetest fssetest.o \
      sselib_c.o $(SSEOBJS)

# Compilation rules

dtimer.o : dtimer.c
	$(CC) $(CCOPTS) -c dtimer.c

#SSE
sselib.o : sselib.c
	$(CC) $(CCOPTS) -o sselib.o -c sselib.c

sselib_c.o : sselib_c.f03
	$(FC90) $(OPTS90) -o sselib_c.o -c $(FF03) sselib_c.f03

cssetest.o : cssetest.c
	$(CC) $(CCOPTS) -c cssetest.c

fssetest.o : fssetest.f90 sselib_c.o
	$(FC90) $(OPTS90) -c fssetest.f90

clean:
	rm -f *.o *.mod

clobber: clean
	rm -f *ssetest
