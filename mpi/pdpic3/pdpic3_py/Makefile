#Makefile for 3D Darwin Python MPI PIC codes

# Makefile gfortran compiler with MacOS X

#PFC = f2py --fcompiler=gnu95
#PCC = f2py --compiler=unix

#OPTS90 = -O3
#OPTS90 = -O3 -fdefault-real-8 -fdefault-double-8
#OPTS90 = -O3 -fcheck=bounds -fdefault-real-8 -Wall -std=f95

#export CFLAGS = -O3
#export CFLAGS = -O3 -Wall

# Makefile Intel compiler with Mac OS X

#PFC = f2py --fcompiler=intelem
#PCC = f2py --compiler=unix

#OPTS90 = -O3
#OPTS90 = -O3 -r8
#OPTS90 = -O3 -CB -r8 -warn all -std90

#export CFLAGS = -O3
#export CFLAGS = -O3 -Wall

# Makefile Intel compiler with Linux

PFC = f2py --fcompiler=intelem
PCC = f2py --compiler=intelem

OPTS90 = -O3
#OPTS90 = -O3 -r8
#OPTS90 = -O3 -CB -r8 -warn all -std90

export CFLAGS = -O3
#export CFLAGS = -O3 -Wall

# to find where MPI is located, execute: mpiifort -show
MPIHOME = /u/local/compilers/intel-cs/2013.0.028/mpi
MPIINC = -I$(MPIHOME)/intel64/include
LOPTS = -L$(MPIHOME)/intel64/lib -lmpigf -lmpi

# Makefile gfortran compiler with Linux

#PFC = f2py --fcompiler=gnu95
#PCC = f2py --compiler=unix

#OPTS90 = -O3
#OPTS90 = -O3 -fdefault-real-8 -fdefault-double-8
#OPTS90 = -O3 -fbounds-check -fdefault-real-8 -Wall -std=f95

#export CFLAGS = -O3
#export CFLAGS = -O3 -Wall

# to find where MPI is located, execute: mpif90 -show
#MPIHOME = /u/local/compilers/intel-cs/2013.0.028/mpi
#MPIINC = -I$(MPIHOME)/intel64/include
#LOPTS = -L$(MPIHOME)/intel64/lib -lmpigf -lmpi

# Makefile PGI compiler with Linux

#PFC = f2py --fcompiler=pg
#PCC = f2py --compiler=unix

#OPTS90 = -O3
#OPTS90 = -O3 -r8
#OPTS90 = -O3 -Mbounds -r8 -Mstandard

#export CFLAGS = -O3
#export CFLAGS = -O3 -Wall

#

# Linkage rules

all : fpdpic3.py

fpdpic3.py : fpdpush3.so dtimer.so

# Compilation rules

dtimer.so : dtimer_h.pyf dtimer.c
	$(PCC) -m dtimer -c dtimer_h.pyf dtimer.c

dtimer_h.pyf : dtimer_hpy.f90
	f2py --overwrite-signature -m dtimer -h dtimer_h.pyf dtimer_hpy.f90

#fprecision.so : fprecision.f90
#	$(PFC) --opt="$(OPTS90)" -m fprecision -c fprecision.f90

# Version using Fortran77 pplib3.f
#fpdpush3.so : pdpush3_h.pyf pdpush3.f pplib3.f
#	$(PFC) --opt="$(OPTS90)" $(LOPTS) -m fpdpush3 -c pdpush3_h.pyf pdpush3.f \
# 	$(MPIINC) pplib3.f

# Version using Fortran90 pplib3.f90
fpdpush3.so : pdpush3_h.pyf pdpush3.f pplib3.f90
	$(PFC) --opt="$(OPTS90)" $(LOPTS) -m fpdpush3 -c pdpush3_h.pyf pdpush3.f \
	$(MPIINC) pplib3.f90

pdpush3_h.pyf : pdpush3pp_hpy.f90
	f2py --overwrite-signature -m fpdpush3 -h pdpush3_h.pyf pdpush3pp_hpy.f90

clean :
	rm -f *.so *.pyf
