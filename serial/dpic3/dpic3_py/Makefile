#Makefile for 3D Darwin Python PIC codes

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

all : fdpic3.py

fdpic3.py : fdpush3.so dtimer.so

# Compilation rules

dtimer.so : dtimer_h.pyf dtimer.c
	f2py --opt="$(CCOPTS)" -m dtimer -c dtimer_h.pyf dtimer.c

dtimer_h.pyf : dtimer_hpy.f90
	f2py --overwrite-signature -m dtimer -h dtimer_h.pyf dtimer_hpy.f90

#fprecision.so : fprecision.f90
#	f2py --opt="$(OPTS90)" -m fprecision -c fprecision.f90

fdpush3.so : dpush3_h.pyf dpush3.f
	f2py $(PFC) --opt="$(OPTS90)" -m fdpush3 -c dpush3_h.pyf dpush3.f

dpush3_h.pyf : dpush3_hpy.f90
	f2py --overwrite-signature -m fdpush3 -h dpush3_h.pyf dpush3_hpy.f90

clean :
	rm -f *.so *.pyf


