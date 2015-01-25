Skeleton 2D Electrostatic MPI Particle-in-Cell (PIC) codes
by Viktor K. Decyk
copyright 2000-2013, regents of the university of california

This program contains sample codes for illustrating the basic structure
of a 2D Electrostatic MPI Particle-in-Cell (PIC) code, in both Fortran
and C.  The codes have no diagnosics except for initial and final
energies.  Their primary purpose is to provide example codes for
physical science students learning about MPI PIC codes.  They are also
intended as benchmark reference codes to aid in developing new codes and
in evaluating new computer architectures.  A serial version of this code
with the same structure (pic2) also exists, and can be compared to this
code in order to understand the parallel algorithms.

PIC codes are widely used in plasma physics.  They model plasmas as
particles which interact self-consistently via the electromagnetic
fields they themselves produce.  PIC codes generally have three
important procedures in the main iteration loop.  The first is the
deposit, where some particle quantity, such as a charge, is accumulated
on a grid via interpolation to produce a source density.  The second
important procedure is the field solver, which solves Maxwell’s equation
or a subset to obtain the electric and/or magnetic fields from the
source densities.  Finally, once the fields are obtained, the particle
forces are found by interpolation from the grid, and the particle
co-ordinates are updated, using Newton’s second law and the Lorentz
force.  The particle processing parts dominate over the field solving
parts in a typical PIC application. 

More details about PIC codes can be found in the texts by C. K. Birdsall
and A. B. Langdon, Plasma Physics via Computer Simulation, 1985,
R. W. Hockney and J. W. Eastwood, Computer Simulation Using Particles,
1981, and John M. Dawson, "Particle simulation of plasmas", Rev. Mod.
Phys. 55, 403 (1983).  Details about the mathematical equations and
units used in this code is given in the companion article,
"Description of Electrostatic Spectral Code from the UPIC Framework" by
Viktor K. Decyk, UCLA, in the file ESModels.pdf.

Details abut MPI can be found in the book by William Gropp, Ewing Lusk,
and Anthony Skjellum, Using MPI: Portable Parallel Programming with the
Message-Passing Interface, The MIT Press, 1994.

No warranty for proper operation of this software is given or implied.
Software or information may be copied, distributed, and used at own
risk; it may not be distributed without this notice included verbatim
with each file.  If use of these codes results in a publication, an
acknowledgement is requested.

The code here uses the simplest force, the electrostatic Coulomb
interaction, obtained by solving a Poisson equation.  A spectral method
using Fast Fourier Transforms (FFTs) is used to solve the Poisson
equation.  A real to complex FFT is used, and the data in Fourier space
is stored in a packed format, where the input and output sizes are the
same.  The boundary conditions are periodic, only electron species are
included, and linear interpolation is used.

For parallelization, the code uses a simple domain decomposition scheme,
where the field quantities (electric field, charge density) are divided
among the computational nodes.  The primary decomposition divides the y
values evenly, that is, each node has all the x values for some y.  The
particles are distributed so that the y co-ordinates of the particles
have a value within the domain.  This simple decomposition works if the
particles are uniformly distributed in space.  Particles at the edge of
the domain may need information from the next domain in order to
interpolate the fields.  To avoid unnecessary communication, one extra
guard cell in y is added at the end of each domain that replicates the
first y value in the next domain.  After particles are updated, some
particles may move to a neighboring domain.  A particle manager
(PPMOVE2) is responsible for moving such particles to the appropriate
domain.  The FFT is performed in 3 steps.  In going from real space to
Fourier space, the FFT is first performed in x for the y values in the
primary domain.  The data is then transposed to a secondary domain
decomposition, where each node has all the y values for some x.  The FFT
is then performed in the y direction for the x values in the secondary
domain.  Poisson's equation is solved using this secondary decomposition
There are four main communication procedures which use MPI.  The first
adds the guard cells for the charge density, the second copies the guard
cells for the electric field.  The third is the particle manager, and
the fourth transposes the data between primary and secondary
decompositions using an all to all communication pattern.  Further
information about the domain decomposition parallel algorithms used can
be found in the companion presentation Dcomp.pdf and in the article:
p. c. liewer and v. k. decyk, j. computational phys. 85, 302 (1989).

Important differences between the push and deposit procedures (in
ppush2.f and ppush2.c) and the serial versions (in push2.f and push2.c
in the pic2 directory) are highlighted in the files dppush2_f.pdf and
dppush2_c.pdf, respectively.

Differences between the main codes (ppic2.f90 and ppic2.c) and the main
codes in the serial versions in the pic2 directory (pic2.f90 and pic2.c)
are highlighted in the files dppic2_f90.pdf and dppic2_c.pdf,
respectively. 

Particles are initialized with a uniform distribution in space and a 
gaussian distribution in velocity space.  This describes a plasma in
thermal equilibrium.  The inner loop contains a charge deposit, add
guard cell procedures, a scalar FFT, a Poisson solver, a vector FFT,
copy guard cell procedures, a particle push, a particle manager, and a
particle sorting procedure.  The final energy and timings are printed.
A sample output file for the default input parameters is included in the
file output.

In more detail, the inner loop of the code contains the following
procedures in Fortran (C):

Deposit section:
   PPGPOST2L (cppgpost2l): deposit charge density
   PPAGUARD2XL (cppaguard2xl): add charge density guard cells in x on
                               local processor
   PPNAGUARD2L (cppnaguard2l): add charge density guard cells in y from
                               remote processor 

Field solve section:
   WPPFFT2R (cwppfft2r): FFT charge density to fourier space
   PPPOIS22 (cppois22): calculate smoothed longitudinal electric field
                        in fourier space.
   WPPFFT2R2 (cwppfft2r2): FFT smoothed electric field to real space

Particle Push section:
   PPNCGUARD2L (cppncguard2l): fill in guard cells for smoothed electric
                               field in y from remote processor
   PPCGUARD2XL (cppcguard2xl): fill in guard cells for smoothed electric
                               field in x field on local processor
   PPGPUSH2L (cppgpush2l): update particle co-ordinates with smoothed
                           electric field. also calculate locations of
                           particles leaving processor for PPMOVE2.
                           x(t)->x(t+dt); v(t-dt/2)->v(t+dt/2)
   PPMOVE2 (cppmove2): moves particles to appropriate processor from
                       from list supplied by PPGPUSH2L
   PPDSORTP2YL (cppdsortp2yl) : sort particles by cell

The inputs to the code are the grid parameters indx, indy, the particle
number parameters npx, npy, the time parameters tend, dt, and the
velocity paramters vtx, vty, vx0, vy0, and the sorting parameter sortime.

In more detail:
indx = exponent which determines length in x direction, nx=2**indx.
indy = exponent which determines length in y direction, ny=2**indy.
   These ensure the system lengths are a power of 2.
npx = number of electrons distributed in x direction.
npy = number of electrons distributed in y direction.
   The total number of particles in the simulation is npx*npy.
tend = time at end of simulation, in units of plasma frequency.
dt = time interval between successive calculations.
   The number of time steps the code runs is given by tend/dt.
   dt should be less than .2 for the electrostatic code.
vtx/vty = thermal velocity of electrons in x/y direction
   a typical value is 1.0.
vx0/vy0 = drift velocity of electrons in x/y direction.
sortime = number of time steps between electron sorting.
   This is used to improve cache performance.  sortime=0 to suppress.

The major program files contained here include:
ppic2.f90    Fortran90 main program 
ppic2.c      C main program
pplib2.f     Fortran77 MPI communications library
pplib2_h.f90 Fortran90 MPI communications interface (header) library
pplib2.f90   Fortran90 MPI communications library
pplib2.c     C MPI communications library
pplib2.h     C MPI communications header library
ppush2.f     Fortran77 procedure library
ppush2_h.f90 Fortran90 procedure interface (header) library
ppush2.c     C procedure library
ppush2.h     C procedure header library
dtimer.c     C timer function, used by both C and Fortran

Files with the suffix .f90 adhere to the Fortran 90 standard, files with
the suffix .f adhere to the Fortran77 standard, files with the suffix .c
and .h adhere to the C99 standard.

The makefile is setup to use gcc and gfortran with Linux.  A version for
Mac OS X is also in the Makefile, but is commented out.  

Two executables can be created, fppic2 for Fortran, and cppic2 for C.

To compile program, execute:

Make program_name

where program_name is either: fppic2 or cppic2, or execute:

make

to create both programs.

To execute, type the name of the executable:

mpirun -np nproc ./program_name

where program_name is either fppic2 or cppic2, and
where nproc is the number of processors to be used.

There is one restriction on the number of processors which can be used:
this simple skeleton code does not support the case where MPI nodes have
zero grid points.  This special case can happen for certain combinations
of the grid size in y (set by the parameter indy) and the number of
processors chosen.  If this happens the code will exit with an error
message.  This special case will never occur if the grid size in y is an
exact multiple of the number of processors.

The file output contains the results produced for the default parameters.
Typical timing results are shown in the file fppic2_bench.pdf.

The Fortran version can be compiled to run with double precision by
changing the Makefile (typically by setting the compiler options flags
-r8).

The libraries pplib2.c and ppush2.c contain wrapper functions to allow
the C libraries to be called from Fortran.  The libraries pplib2_f.c and
ppush2_f.c contain wrapper functions to allow the Fortran libraries to
be called from C.

