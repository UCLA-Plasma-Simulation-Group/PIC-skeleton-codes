#Makefile for 2D Electrostatic Python OpenMP PIC codes

# Makefile gfortran compiler with MacOS X

#PFC = f2py --fcompiler=gnu95
#PCC = f2py --compiler=unix

#OPTS90 = -O3 -fopenmp
#OPTS90 = -O3 -fopenmp -fdefault-real-8 -fdefault-double-8
#OPTS90 = -O3 -fopenmp -fcheck=bounds -fdefault-real-8 -Wall -std=f95

#export CFLAGS = -O3 -fopenmp
#export CFLAGS = -O3 -fopenmp -Wall

#LOPTS = -lgomp

# Makefile Intel compiler with Mac OS X

#PFC = f2py --fcompiler=intelem
#PCC = f2py --compiler=unix

#OPTS90 = -O3 -openmp
#OPTS90 = -O3 -openmp -r8
#OPTS90 = -O3 -openmp -CB -r8 -warn all -std90

#export CFLAGS = -O3 -fopenmp
#export CFLAGS = -O3 -fopenmp -Wall

# Makefile Intel compiler with Linux

#PFC = f2py --fcompiler=intelem
#PCC = f2py --compiler=intelem

#OPTS90 = -O3 -openmp
#OPTS90 = -O3 -openmp -r8
#OPTS90 = -O3 -openmp -CB -r8 -warn all -std90

#export CFLAGS = -O3 -openmp
#export CFLAGS = -O3 -openmp -Wall

# Makefile gfortran compiler with Linux

PFC = f2py --fcompiler=gnu95
PCC = f2py --compiler=unix

OPTS90 = -O3 -fopenmp
#OPTS90 = -O3 -fopenmp -fdefault-real-8 -fdefault-double-8
#OPTS90 = -O3 -fopenmp -fbounds-check -fdefault-real-8 -Wall -std=f95

export CFLAGS = -O3 -fopenmp
#export CFLAGS = -O3 -fopenmp -Wall

LOPTS = -lgomp

# Makefile PGI compiler with Linux

#PFC = f2py --fcompiler=pg
#PCC = f2py --compiler=unix

#OPTS90 = -O3
#OPTS90 = -O3 -r8
#OPTS90 = -O3 -Mbounds -r8 -Mstandard

#export CFLAGS = -O3 -fopenmp
#export CFLAGS = -O3 -fopenmp -Wall

#

# Linkage rules

all : fmpic2.py cmpic2.py

fmpic2.py : fmpush2.so dtimer.so fomplib.so

cmpic2.py : cmpush2.so dtimer.so complib.so

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

fmpush2.so : mpush2_h.pyf mpush2.f
	$(PFC) --opt="$(OPTS90)" $(LOPTS) -m fmpush2 -c mpush2_h.pyf mpush2.f

mpush2_h.pyf : mpush2_hpy.f90
	f2py --overwrite-signature -m fmpush2 -h mpush2_h.pyf mpush2_hpy.f90

cmpush2.so : cmpush2_h.pyf mpush2.c
	$(PCC) $(LOPTS) -m cmpush2 -c cmpush2_h.pyf mpush2.c

cmpush2_h.pyf : cmpush2_hpy.f90
	f2py --overwrite-signature -m cmpush2 -h cmpush2_h.pyf cmpush2_hpy.f90

clean :
	rm -f *.so *.pyf


