#Makefile for 2-1/2D Electromagnetic Python OpenMP PIC codes

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

all : fmbpic2.py cmbpic2.py

fmbpic2.py : fmbpush2.so dtimer.so fomplib.so

cmbpic2.py : cmbpush2.so dtimer.so complib.so

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

fmbpush2.so : mbpush2_h.pyf mbpush2.f
	$(PFC) --opt="$(OPTS90)" $(LOPTS) -m fmbpush2 -c mbpush2_h.pyf mbpush2.f

mbpush2_h.pyf : mbpush2_hpy.f90
	f2py --overwrite-signature -m fmbpush2 -h mbpush2_h.pyf mbpush2_hpy.f90

cmbpush2.so : cmbpush2_h.pyf mbpush2.c
	$(PCC) $(LOPTS) -m cmbpush2 -c cmbpush2_h.pyf mbpush2.c

cmbpush2_h.pyf : cmbpush2_hpy.f90
	f2py --overwrite-signature -m cmbpush2 -h cmbpush2_h.pyf cmbpush2_hpy.f90

clean :
	rm -f *.so *.pyf


