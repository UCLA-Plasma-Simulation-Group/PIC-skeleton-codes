#Makefile for 3D Electromagnetic Python MPI/OpenMP PIC codes

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

# Makefile PGI compiler with Linux

#PFC = f2py --fcompiler=pg
#PCC = f2py --compiler=unix

#OPTS90 = -O3
#OPTS90 = -O3 -r8
#OPTS90 = -O3 -Mbounds -r8 -Mstandard

#export CFLAGS = -O3 -fopenmp
#export CFLAGS = -O3 -fopenmp -Wall

#LOPTS =

#

# Linkage rules

all : fmpbpic3.py

fmpbpic3.py : fmpbpush3.so dtimer.so fomplib.so

# Compilation rules

dtimer.so : dtimer_h.pyf dtimer.c
	$(PCC) -m dtimer -c dtimer_h.pyf dtimer.c

dtimer_h.pyf : dtimer_hpy.f90
	f2py --overwrite-signature -m dtimer -h dtimer_h.pyf dtimer_hpy.f90

fomplib.so : omplib_h.pyf omplib.f
	$(PFC) --opt="$(OPTS90)" $(LOPTS) -m fomplib -c omplib_h.pyf omplib.f

omplib_h.pyf : omplib_hpy.f90
	f2py --overwrite-signature -m fomplib -h omplib_h.pyf omplib_hpy.f90

#fprecision.so : fprecision.f90
#	$(PFC) --opt="$(OPTS90)" -m fprecision -c fprecision.f90

# Version using Fortran77 mpplib3.f
#fmpbpush3.so : mpbpush3_h.pyf mpbpush3.f mpplib3.f
#	$(PFC) --opt="$(OPTS90)" $(MLOPTS) $(LOPTS) -m fmpbpush3 -c mpbpush3_h.pyf \
#   mpbpush3.f $(MPIINC) mpplib3.f

# Version using Fortran90 mpplib3.f90
fmpbpush3.so : mpbpush3_h.pyf mpbpush3.f mpplib3.f90
	$(PFC) --opt="$(OPTS90)" $(MLOPTS) $(LOPTS) -m fmpbpush3 -c mpbpush3_h.pyf \
	mpbpush3.f $(MPIINC) mpplib3.f90

mpbpush3_h.pyf : mpbpush3pp_hpy.f90
	f2py --overwrite-signature -m fmpbpush3 -h mpbpush3_h.pyf mpbpush3pp_hpy.f90

clean :
	rm -f *.so *.pyf


