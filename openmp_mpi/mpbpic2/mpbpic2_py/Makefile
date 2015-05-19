#Makefile for 2-1/2D Electromagnetic Python MPI/OpenMP PIC codes

# Makefile gfortran compiler with MacOS X

#PFC = f2py --fcompiler=gnu95
#PCC = f2py --compiler=unix

#OPTS90 = -O3 -fopenmp
#OPTS90 = -O3 -fopenmp -fdefault-real-8 -fdefault-double-8
#OPTS90 = -O3 -fopenmp -fcheck=bounds -fdefault-real-8 -Wall -std=f95

#export CFLAGS = -O3 -fopenmp
#export CFLAGS = -O3 -fopenmp -Wall

#LOPTS = -lgomp
#LEGACY =

# Makefile Intel compiler with Mac OS X

#PFC = f2py --fcompiler=intelem
#PCC = f2py --compiler=unix

#OPTS90 = -O3 -openmp
#OPTS90 = -O3 -openmp -r8
#OPTS90 = -O3 -openmp -CB -r8 -warn all -std90

#export CFLAGS = -O3 -fopenmp
#export CFLAGS = -O3 -fopenmp -Wall

#LOPTS = -liomp5 -lirc
#LEGACY = -nofor_main

# Makefile Intel compiler with Linux

PFC = f2py --fcompiler=intelem
PCC = f2py --compiler=unix

OPTS90 = -O3 -openmp
#OPTS90 = -O3 -openmp -r8
#OPTS90 = -O3 -openmp -CB -r8 -warn all -std90

export CFLAGS = -O3 -fopenmp
#export CFLAGS = -O3 -fopenmp -Wall

# to find where MPI is located, execute: mpiifort -show
MPIHOME = /u/local/compilers/intel-cs/2013.0.028/mpi
MPIINC = -I$(MPIHOME)/intel64/include
MLOPTS = -L$(MPIHOME)/intel64/lib -lmpigf -lmpi

LOPTS = -liomp5 -lirc
LEGACY = -nofor_main

# Makefile gfortran compiler with Linux

#PFC = f2py --fcompiler=gnu95
#PCC = f2py --compiler=unix

#OPTS90 = -O3 -fopenmp
#OPTS90 = -O3 -fopenmp -fdefault-real-8 -fdefault-double-8
#OPTS90 = -O3 -fopenmp -fbounds-check -fdefault-real-8 -Wall -std=f95

#export CFLAGS = -O3 -fopenmp
#export CFLAGS = -O3 -fopenmp -Wall

# to find where MPI is located, execute: mpif90 -show
#MPIHOME = /u/local/compilers/intel-cs/2013.0.028/mpi
#MPIINC = -I$(MPIHOME)/intel64/include
#MLOPTS = -L$(MPIHOME)/intel64/lib -lmpigf -lmpi

#LOPTS = -lgomp
#LEGACY =

# Makefile PGI compiler with Linux

#PFC = f2py --fcompiler=pg
#PCC = f2py --compiler=unix

#OPTS90 = -O3
#OPTS90 = -O3 -r8
#OPTS90 = -O3 -Mbounds -r8 -Mstandard

#export CFLAGS = -O3 -fopenmp
#export CFLAGS = -O3 -fopenmp -Wall

#LOPTS =
#LEGACY =

#

# Linkage rules

all : fmpbpic2.py cmpbpic2.py

fmpbpic2.py : fmpbpush2.so dtimer.so fomplib.so

cmpbpic2.py : cmpbpush2.so dtimer.so complib.so

# Compilation rules

dtimer.so : dtimer_h.pyf dtimer.c
	$(PCC) -m dtimer -c dtimer_h.pyf dtimer.c

dtimer_h.pyf : dtimer_hpy.f90
	f2py --overwrite-signature -m dtimer -h dtimer_h.pyf dtimer_hpy.f90

fomplib.so : omplib_h.pyf omplib.f
	$(PFC) --opt="$(OPTS90)" $(LOPTS) -m fomplib -c omplib_h.pyf omplib.f

omplib_h.pyf : omplib_hpy.f90
	f2py --overwrite-signature -m fomplib -h omplib_h.pyf omplib_hpy.f90

complib.so : complib_h.pyf omplib.c
	$(PCC) $(LOPTS) -m complib -c complib_h.pyf omplib.c

complib_h.pyf : complib_hpy.f90
	f2py --overwrite-signature -m complib -h complib_h.pyf complib_hpy.f90

#fprecision.so : fprecision.f90
#	$(PFC) --opt="$(OPTS90)" -m fprecision -c fprecision.f90

# Version using Fortran77 mpplib2.f
#fmpbpush2.so : mpbpush2_h.pyf mpbpush2.f mpplib2.f
#	$(PFC) --opt="$(OPTS90)" $(MLOPTS) $(LOPTS) -m fmpbpush2 -c mpbpush2_h.pyf \
#   mpbpush2.f $(MPIINC) mpplib2.f

# Version using Fortran90 mpplib2.f90
fmpbpush2.so : mpbpush2_h.pyf mpbpush2.f mpplib2.f90
	$(PFC) --opt="$(OPTS90)" $(MLOPTS) $(LOPTS) -m fmpbpush2 -c mpbpush2_h.pyf \
	mpbpush2.f $(MPIINC) mpplib2.f90

mpbpush2_h.pyf : mpbpush2pp_hpy.f90
	f2py --overwrite-signature -m fmpbpush2 -h mpbpush2_h.pyf mpbpush2pp_hpy.f90

cmpbpush2.so : cmpbpush2_h.pyf mpbpush2.c mpplib2.c
	$(PCC) $(MLOPTS) $(LOPTS) -m cmpbpush2 -c cmpbpush2_h.pyf mpbpush2.c \
	$(MPIINC) mpplib2.c

cmpbpush2_h.pyf : cmpbpush2pp_hpy.f90
	f2py --overwrite-signature -m cmpbpush2 -h cmpbpush2_h.pyf cmpbpush2pp_hpy.f90

clean :
	rm -f *.so *.pyf


