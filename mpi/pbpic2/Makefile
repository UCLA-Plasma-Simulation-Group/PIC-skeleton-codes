#Makefile for 2-1/2D Electromagnetic MPI PIC codes

# Makefile gfortran compiler with MacOS X

#OpenMPI
#MPIFC = mpif90
#MPICC = mpicc

#CC = gcc

#OPTS90 = -O3
#OPTS90 = -O3 -fdefault-real-8 -fdefault-double-8
#OPTS90 = -O3 -fcheck=bounds -fdefault-real-8 -fdefault-double-8 -Wall -std=f95

#OPTS77 = -O3
#OPTS77 = -O3 -fdefault-real-8 -fdefault-double-8
#OPTS77 = -O3 -fcheck=bounds -fdefault-real-8 -fdefault-double-8 -Wall -std=f95

#CCOPTS = -O3 -std=c99
#CCOPTS = -O3 -Wall -std=c99
#LEGACY =
#LOPTS = -lm

# Makefile Intel compiler with Mac OS X

# OpenMPI
#MPIFC = mpif90
#MPICC = mpicc

#CC = gcc

#OPTS90 = -O3
#OPTS90 = -O3 -r8
#OPTS90 = -O3 -CB -r8 -warn all -std90

#OPTS77 = -O3
#OPTS77 = -O3 -r8
#OPTS77 = -O3 -CB -r8 -warn all -std77

#CCOPTS = -O3 -std=c99
#CCOPTS = -O3 -Wall -std=c99
#LEGACY = -nofor_main

# Makefile Intel compiler with Linux

# MPI
#MPIFC = mpiifort
#MPICC = mpiicc

#CC = icc

#OPTS90 = -O3
#OPTS90 = -O3 -r8
#OPTS90 = -O3 -CB -r8 -warn all -std90

#OPTS77 = -O3
#OPTS77 = -O3 -r8
#OPTS77 = -O3 -CB -r8 -warn all -std77

#CCOPTS = -O3 -std=c99
#CCOPTS = -O3 -no-vec -Wall -std=c99
#LEGACY = -nofor_main

# Makefile gfortran compiler with Linux

# OpenMPI
MPIFC = mpif90
MPICC = mpicc

CC = gcc

OPTS90 = -O3
#OPTS90 = -O3 -fdefault-real-8 -fdefault-double-8
#OPTS90 = -O3 -fbounds-check -fdefault-real-8 -fdefault-double-8 -Wall -std=f95

OPTS77 = -O3
#OPTS77 = -O3 -fdefault-real-8 -fdefault-double-8
#OPTS77 = -O3 -fbounds-check -fdefault-real-8 -fdefault-double-8 -Wall -std=f95

CCOPTS = -O3 -std=c99
#CCOPTS = -O3 -Wall -std=c99
#LEGACY =
LOPTS = -lm

# Makefile PGI compiler with Linux

# OpenMPI
#MPIFC = mpipgf90
#MPICC = mpicc

#CC = gcc

#OPTS90 = -O3
#OPTS90 = -O3 -r8
#OPTS90 = -O3 -Mbounds -r8 -Mstandard

#OPTS77 = -O3
#OPTS77 = -O3 -r8
#OPTS77 = -O3 -Mbounds -r8 -Mstandard

#CCOPTS = -O3 -std=c99
#CCOPTS = -O3 -Wall -std=c99
#LEGACY = -Mnomain

# Makefile Cray compiler with Linux

#MPIFC = ftn
#MPICC = cc

#OPTS90 = -O 3
#OPTS90 = -O 3 -s real64
#OPTS90 = -O 3 -R b -s real64 -en

#OPTS77 = -O 3
#OPTS77 = -O 3 -s real64
#OPTS77 = -O 3 -R b -s real64 -en

#CCOPTS = -O 3 -h c99 -h conform
#LEGACY =

#

# Linkage rules

all : fpbpic2 cpbpic2

special: fpbpic2_c cpbpic2_f

# Version using Fortran77 pplib2.f
#fpbpic2 : fpbpic2.o fpbpush2.o fpplib2.o dtimer.o
#	$(MPIFC) $(OPTS90) $(LOPTS) -o fpbpic2 \
#       fpbpic2.o fpbpush2.o fpplib2.o pbpush2_h.o pplib2_h.o dtimer.o

# Version using Fortran90 pplib2.f90
fpbpic2 : fpbpic2.o fpbpush2.o dtimer.o
	$(MPIFC) $(OPTS90) $(LOPTS) -o fpbpic2 \
        fpbpic2.o fpbpush2.o f90pplib2.o pbpush2_h.o dtimer.o

cpbpic2 : cpbpic2.o cpbpush2.o cpplib2.o dtimer.o
	$(MPICC) $(CCOPTS) $(LOPTS) -o cpbpic2 \
        cpbpic2.o cpbpush2.o cpplib2.o dtimer.o

fpbpic2_c : fpbpic2_c.o cpbpush2.o cpplib2.o dtimer.o
	$(MPIFC) $(OPTS90) $(LOPTS) -o fpbpic2_c \
        fpbpic2_c.o cpbpush2.o cpplib2.o dtimer.o

cpbpic2_f : cpbpic2.o cpbpush2_f.o cpplib2_f.o fpbpush2.o fpplib2.o dtimer.o
	$(MPIFC) $(OPTS90) $(LOPTS) $(LEGACY) -o cpbpic2_f \
        cpbpic2.o cpbpush2_f.o cpplib2_f.o fpbpush2.o fpplib2.o dtimer.o

# Compilation rules

dtimer.o : dtimer.c
	$(CC) $(CCOPTS) -c dtimer.c

fpplib2.o : pplib2.f
	$(MPIFC) $(OPTS77) -o fpplib2.o -c pplib2.f

pplib2_h.o : pplib2_h.f90
	$(MPIFC) $(OPTS90) -o pplib2_h.o -c pplib2_h.f90

f90pplib2.o : pplib2.f90
	$(MPIFC) $(OPTS90) -o f90pplib2.o -c pplib2.f90

cpplib2.o : pplib2.c
	$(MPICC) $(CCOPTS) -o cpplib2.o -c pplib2.c

cpplib2_f.o : pplib2_f.c
	$(MPICC) $(CCOPTS) -o cpplib2_f.o -c pplib2_f.c

fpbpush2.o : pbpush2.f
	$(MPIFC) $(OPTS90) -o fpbpush2.o -c pbpush2.f

pbpush2_h.o : pbpush2_h.f90
	$(MPIFC) $(OPTS90) -o pbpush2_h.o -c pbpush2_h.f90

cpbpush2.o : pbpush2.c
	$(MPICC) $(CCOPTS) -o cpbpush2.o -c pbpush2.c

cpbpush2_f.o : pbpush2_f.c
	$(MPICC) $(CCOPTS) -o cpbpush2_f.o -c pbpush2_f.c

# Version using Fortran77 pplib2.f
#fpbpic2.o : pbpic2.f90 pbpush2_h.o pplib2_h.o
#	$(MPIFC) $(OPTS90) -o fpbpic2.o -c pbpic2.f90

# Version using Fortran90 pplib2.f90
fpbpic2.o : pbpic2.f90 f90pplib2.o pbpush2_h.o
	$(MPIFC) $(OPTS90) -o fpbpic2.o -c pbpic2.f90

cpbpic2.o : pbpic2.c
	$(MPICC) $(CCOPTS) -o cpbpic2.o -c pbpic2.c

fpbpic2_c.o : pbpic2_c.f90
	$(MPIFC) $(OPTS90) -o fpbpic2_c.o -c pbpic2_c.f90

clean :
	rm -f *.o *.mod

clobber: clean
	rm -f fpbpic2 cpbpic2 fpbpic2_c cpbpic2_f
