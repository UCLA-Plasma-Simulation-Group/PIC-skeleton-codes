(This is the Markdown annotated version of the plaintext README and contains the same information)

Skeleton 2D Electrostatic Particle-in-Cell (PIC) codes
by Viktor K. Decyk
copyright 1994-2013, regents of the university of california

This program contains sample codes for illustrating the basic structure
of a 2D Electrostatic Particle-in-Cell (PIC) code, in both Fortran
and C.  The codes have no diagnosics except for initial and final
energies.  Their primary purpose is to provide example codes for
physical science students learning about PIC codes.  They are also
intended as benchmark reference codes to aid in developing new codes and
in evaluating new computer architectures.

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

Particles are initialized with a uniform distribution in space and a
gaussian distribution in velocity space.  This describes a plasma in
thermal equilibrium.  The inner loop contains a charge deposit, an add
guard cell procedure, a scalar FFT, a Poisson solver, a vector FFT, a
copy guard cell procedure, a particle push, and a particle sorting
procedure.  The final energy and timings are printed.  A sample output
file for the default input parameters is included in the file output.

#Code Description
In more detail, the inner loop of the code contains the following
procedures in Fortran (C):

__Deposit section__:

|Fortran Name|C Name| description|
| ---------- | ---- | ---------- |
|GPOST2L|cgpost2l|deposit charge density|
|AGUARD2L|caguard2)|add charge density guard cells|

__Field solve section__:

|Fortran Name|C Name| description|
| ---------- | ---- | ---------- |
|WFFT2RX|cwfft2rx|FFT charge density to fourier space
|POIS22|cpois22|calculate smoothed longitudinal electric field in fourier space.|
|WFFT2R2|cwfft2r2|FFT smoothed electric field to real space|

#Particle Push section:
|Fortran Name|C Name| description|
| ---------- | ---- | ---------- |
|CGUARD2L|ccguard2l|fill in guard cells for smoothed electric field|
|GPUSH2L|cgpush2l|update particle co-ordinates with smoothed electric field: <br> x(t)->x(t+dt); v(t-dt/2)->v(t+dt/2)|
|DSORTP2YL|cdsortp2yl|sort particles by cell|

#Input Paramters
There is no input file. For simplicity, input parameters are set at within the source code near at the start of the _main_ sections of each programs (i.e at the top of the the file `pic2.c` for the C version and `pic2.f90` for the Fortran version ). Repository versions come with resonable defaults and do not require editing if you wish to jsut run a simulation out-of-the-box.

The inputs to the code are the grid parameters (_indx_, _indy_) the particle
number parameters (_npx_, _npy_), the time parameters (_tend_, _dt_), the velocity
parameters (_vtx_, _vty_, _vx0_, _vy0_), and the sorting parameter _sortime_.

In more detail:

|Parameter|description|
| ------- | --------- |
|indx|exponent which determines length in x direction, nx=2**indx|
|indy|exponent which determines length in y direction, ny=2**indy|

   These ensure the system lengths are a power of 2.

|Parameter|description|
| ------- | --------- |
|npx|exponent which determines length in x direction, nx=2**indx|
|npy|exponent which determines length in y direction, ny=2**indy|

   The total number of particles in the simulation is npx*npy.

|Parameter|description|
| ------- | --------- |
|tend| time at end of simulation, in units of plasma frequency.|
|dt|time interval between successive calculations.|

   The number of time steps the code runs is given by tend/dt.
   dt should be less than .2 for the electrostatic code.

|Parameter|description|
| ------- | --------- |
|vtx/vty|thermal velocity of electrons in x/y direction. A typical value is 1.0. |
|vx0/vy0|drift velocity of electrons in x/y direction.|
|sortime|number of time steps between electron sorting. This is used to improve cache performance.  sortime=0 to suppress.|

# File listing
The major program files contained here include:


|File Name|description|
| ------- | --------- |
|pic2.f90|Fortran90 main program|
|pic2.c| C main program|
|push2.f|Fortran77 procedure library|
|push2_h.f90|Fortran90 procedure interface (header) library|
|push2.c| C procedure library|
|push2.h| C procedure header library|
|dtimer.c|C timer function, used by both C and Fortran|


Files with the suffix .f90 adhere to the Fortran 90 standard, files with
the suffix .f adhere to the Fortran77 standard, files with the suffix .c
and .h adhere to the C99 standard.

The makefile is setup to use gcc and gfortran with Linux.  A version for
Mac OS X is also in the Makefile, but is commented out.  

Two executables can be created, _fpic2_ for Fortran and _cpic2_ for C.

To compile program, execute:

```
make program_name
```

where program_name is either: _fpic2_ or _cpic2_, or execute:

```
make
```

to create both programs.

To execute, type the name of the executable:

```
./program_name
```

where program_name is either _fpic2_ or _cpic2_.

The file output contains the results produced for the default parameters.

The Fortran version can be compiled to run with double precision by
changing the Makefile (typically by setting the compiler options flags
-r8).

The library `push2.c` contains wrapper functions to allow the C library to
be called from Fortran. The library `push2_f.c` contains wrapper functions
to allow the Fortran library to be called from C.
