#Makefile for 2D Electrostatic Python PIC codes

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

all : fpic2.py cpic2.py

fpic2.py : fpush2.so dtimer.so

cpic2.py : cpush2.so dtimer.so

# Compilation rules

dtimer.so : dtimer_h.pyf dtimer.c
	$(PCC) -m dtimer -c dtimer_h.pyf dtimer.c

dtimer_h.pyf : dtimer_hpy.f90
	f2py --overwrite-signature -m dtimer -h dtimer_h.pyf dtimer_hpy.f90

#fprecision.so : fprecision.f90
#	$(PFC) --opt="$(OPTS90)" -m fprecision -c fprecision.f90

fpush2.so : push2_h.pyf push2.f
	$(PFC) --opt="$(OPTS90)" -m fpush2 -c push2_h.pyf push2.f

push2_h.pyf : push2_hpy.f90
	f2py --overwrite-signature -m fpush2 -h push2_h.pyf push2_hpy.f90

cpush2.so : cpush2_h.pyf push2.c
	$(PCC) -m cpush2 -c cpush2_h.pyf push2.c

cpush2_h.pyf : cpush2_hpy.f90
	f2py --overwrite-signature -m cpush2 -h cpush2_h.pyf cpush2_hpy.f90

clean :
	rm -f *.so *.pyf


