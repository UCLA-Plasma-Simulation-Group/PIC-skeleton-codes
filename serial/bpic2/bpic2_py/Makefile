#Makefile for 2-1/2D Electromagnetic Python PIC codes

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

all : fbpic2.py cbpic2.py

fbpic2.py : fbpush2.so dtimer.so

cbpic2.py : cbpush2.so dtimer.so

# Compilation rules

dtimer.so : dtimer_h.pyf dtimer.c
	$(PCC) -m dtimer -c dtimer_h.pyf dtimer.c

dtimer_h.pyf : dtimer_hpy.f90
	f2py --overwrite-signature -m dtimer -h dtimer_h.pyf dtimer_hpy.f90

#fprecision.so : fprecision.f90
#	$(PFC) --opt="$(OPTS90)" -m fprecision -c fprecision.f90

fbpush2.so : bpush2_h.pyf bpush2.f
	$(PFC) --opt="$(OPTS90)" -m fbpush2 -c bpush2_h.pyf bpush2.f

bpush2_h.pyf : bpush2_hpy.f90
	f2py --overwrite-signature -m fbpush2 -h bpush2_h.pyf bpush2_hpy.f90

cbpush2.so : cbpush2_h.pyf bpush2.c
	$(PCC) -m cbpush2 -c cbpush2_h.pyf bpush2.c

cbpush2_h.pyf : cbpush2_hpy.f90
	f2py --overwrite-signature -m cbpush2 -h cbpush2_h.pyf cbpush2_hpy.f90

clean :
	rm -f *.so *.pyf


