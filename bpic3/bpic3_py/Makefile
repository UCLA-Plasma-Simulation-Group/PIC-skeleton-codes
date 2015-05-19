#Makefile for 3D Electromagnetic Python PIC codes

# Makefile gfortran compiler with MacOS X

#PFC = --fcompiler=gnu95

#OPTS90 = -O3
#OPTS90 = -O3 -fdefault-real-8 -fdefault-double-8
#OPTS90 = -O3 -fcheck=bounds -fdefault-real-8 -Wall -std=f95

#CCOPTS = -O3
#CCOPTS = -O3 -Wall

# Makefile Intel compiler with Mac OS X

#PFC = --fcompiler=intel

#OPTS90 = -O3
#OPTS90 = -O3 -r8
#OPTS90 = -O3 -CB -r8 -warn all -std90

#CCOPTS = -O3
#CCOPTS = -O3 -Wall

# Makefile Intel compiler with Linux

#PFC = --fcompiler=intel

#OPTS90 = -O3
#OPTS90 = -O3 -r8
#OPTS90 = -O3 -CB -r8 -warn all -std90

#CCOPTS = -O3
#CCOPTS = -O3 -Wall

# Makefile gfortran compiler with Linux

PFC = --fcompiler=gnu95

OPTS90 = -O3
#OPTS90 = -O3 -fdefault-real-8 -fdefault-double-8
#OPTS90 = -O3 -fbounds-check -fdefault-real-8 -Wall -std=f95

CCOPTS = -O3
#CCOPTS = -O3 -Wall

# Makefile PGI compiler with Linux

#PFC = --fcompiler=pg

#OPTS90 = -O3
#OPTS90 = -O3 -r8
#OPTS90 = -O3 -Mbounds -r8 -Mstandard

#CCOPTS = -O3
#CCOPTS = -O3 -Wall

#

# Linkage rules

all : fbpic3.py

fbpic3.py : fbpush3.so dtimer.so

# Compilation rules

dtimer.so : dtimer_h.pyf dtimer.c
	f2py --opt="$(CCOPTS)" -m dtimer -c dtimer_h.pyf dtimer.c

dtimer_h.pyf : dtimer_hpy.f90
	f2py --overwrite-signature -m dtimer -h dtimer_h.pyf dtimer_hpy.f90

#fprecision.so : fprecision.f90
#	f2py --opt="$(OPTS90)" -m fprecision -c fprecision.f90

fbpush3.so : bpush3_h.pyf bpush3.f
	f2py $(PFC) --opt="$(OPTS90)" -m fbpush3 -c bpush3_h.pyf bpush3.f

bpush3_h.pyf : bpush3_hpy.f90
	f2py --overwrite-signature -m fbpush3 -h bpush3_h.pyf bpush3_hpy.f90

clean :
	rm -f *.so *.pyf


