#Makefile for 2-1/2D Darwin Python PIC codes

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

#PFC = f2py --fcompiler=intelem
#PCC = f2py --compiler=intelem

#OPTS90 = -O3
#OPTS90 = -O3 -r8
#OPTS90 = -O3 -CB -r8 -warn all -std90

#export CFLAGS = -O3
#export CFLAGS = -O3 -Wall

# Makefile gfortran compiler with Linux

PFC = f2py --fcompiler=gnu95
PCC = f2py --compiler=unix

OPTS90 = -O3
#OPTS90 = -O3 -fdefault-real-8 -fdefault-double-8
#OPTS90 = -O3 -fbounds-check -fdefault-real-8 -Wall -std=f95

#export CFLAGS = -O3
#export CFLAGS = -O3 -Wall

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

all : fdpic2.py cdpic2.py

fdpic2.py : fdpush2.so dtimer.so

cdpic2.py : cdpush2.so dtimer.so

# Compilation rules

dtimer.so : dtimer_h.pyf dtimer.c
	$(PCC) -m dtimer -c dtimer_h.pyf dtimer.c

dtimer_h.pyf : dtimer_hpy.f90
	f2py --overwrite-signature -m dtimer -h dtimer_h.pyf dtimer_hpy.f90

#fprecision.so : fprecision.f90
#	$(PFC) --opt="$(OPTS90)" -m fprecision -c fprecision.f90

fdpush2.so : dpush2_h.pyf dpush2.f
	$(PFC) --opt="$(OPTS90)" -m fdpush2 -c dpush2_h.pyf dpush2.f

dpush2_h.pyf : dpush2_hpy.f90
	f2py --overwrite-signature -m fdpush2 -h dpush2_h.pyf dpush2_hpy.f90

cdpush2.so : cdpush2_h.pyf dpush2.c
	$(PCC) -m cdpush2 -c cdpush2_h.pyf dpush2.c

cdpush2_h.pyf : cdpush2_hpy.f90
	f2py --overwrite-signature -m cdpush2 -h cdpush2_h.pyf cdpush2_hpy.f90

clean :
	rm -f *.so *.pyf


