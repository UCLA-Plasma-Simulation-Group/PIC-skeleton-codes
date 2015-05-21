#-----------------------------------------------------------------------
# Skeleton 1-2/2D Darwin PIC code
# written by Viktor K. Decyk, Adam Tableman, and Qiyang Hu, UCLA
import math
import numpy
from fdpush1 import *
from dtimer import *

int_type = numpy.int32
double_type = numpy.float64
float_type = numpy.float32
complex_type = numpy.complex64

# indx = exponent which determines grid points in x direction:
# nx = 2**indx.
indx =   9
# npx = number of electrons distributed in x direction.
npx =  18432
# tend = time at end of simulation, in units of plasma frequency.
# dt = time interval between successive calculations.
# qme = charge on electron, in units of e.
tend = 10.0; dt = 0.1; qme = -1.0
# vtx/vty = thermal velocity of electrons in x/y direction
# vx0/vy0 = drift velocity of electrons in x/y direction.
vtx = 1.0; vty = 1.0; vx0 = 0.0; vy0 = 0.0
# vtx/vz0 = thermal/drift velocity of electrons in z direction
vtz = 1.0; vz0 = 0.0
# ax = smoothed particle size in x direction
# ci = reciprocal of velocity of light.
ax = .912871; ci = 0.1
# idimp = number of particle coordinates = 4
# ipbc = particle boundary condition: 1 = periodic
# sortime = number of time steps between standard electron sorting
idimp = 4; ipbc = 1; sortime = 50
# omx/omy/omz = magnetic field electron cyclotron frequency in x/y/z 
omx = 0.4; omy = 0.0; omz = 0.0
# ndc = number of corrections in darwin iteration
ndc = 1
# wke/we = particle kinetic/electrostatic field energy
# wf/wm/wt = magnetic field/transverse electric field/total energy
wke = numpy.zeros((1),float_type)
we = numpy.zeros((1),float_type)
wf = numpy.zeros((1),float_type)
wm = numpy.zeros((1),float_type)
wt = numpy.zeros((1),float_type)
zero = 0.0
# declare scalars for standard code
wpmax = numpy.empty((1),float_type)
wpmin = numpy.empty((1),float_type)

# declare and initialize timing data
itime = numpy.empty((4),numpy.int32)
tdpost = 0.0; tguard = 0.0; tfft = 0.0; tfield = 0.0
tdjpost = 0.0; tdcjpost = 0.0; tpush = 0.0; tsort = 0.0
dtime = numpy.empty((1),double_type)

# initialize scalars for standard code
# np = total number of particles in simulation
# nx = number of grid points in x direction
np = npx; nx = int(math.pow(2,indx)); nxh = int(nx/2)
nxe = nx + 2; nxeh = nxe/2; nx1 = nx + 1
# nloop = number of time steps in simulation
# ntime = current time step
nloop = int(tend/dt + .0001); ntime = 0
qbme = qme
affp = float(nx)/float(np)

# allocate data for standard code
# part, part2 = particle arrays
part = numpy.empty((idimp,np),float_type,'F')
if (sortime > 0):
   part2 = numpy.empty((idimp,np),float_type,'F')
# qe = electron charge density with guard cells
qe = numpy.empty((nxe),float_type,'F')
# fxe = smoothed longitudinal electric field with guard cells
fxe = numpy.empty((nxe),float_type,'F')
# cue = electron current density with guard cells
cue = numpy.empty((2,nxe),float_type,'F')
# dcu = acceleration density with guard cells
dcu = numpy.empty((2,nxe),float_type,'F')
# cus = transverse electric field with guard cells
cus = numpy.empty((2,nxe),float_type,'F')
# amu = momentum flux with guard cells
amu = numpy.empty((2,nxe),float_type,'F')
# exyze = smoothed total electric field with guard cells
exyze = numpy.empty((3,nxe),float_type,'F')
# byze = smoothed magnetic field with guard cells
byze = numpy.empty((2,nxe),float_type,'F')
# ffc, ffe = form factor arrays for poisson solvers
ffc = numpy.empty((nxh),complex_type,'F')
ffe = numpy.empty((nxh),complex_type,'F')
# mixup = bit reverse table for FFT
mixup = numpy.empty((nxh),int_type,'F')
# sct = sine/cosine table for FFT
sct = numpy.empty((nxh),complex_type,'F')
# npic = scratch array for reordering particles
npic = numpy.empty((nx1),int_type,'F')
# gxe, gyze = scratch arrays for fft
gxe = numpy.empty((nxe),float_type,'F')
gyze = numpy.empty((2,nxe),float_type,'F')

# prepare fft tables
wfft1rinit(mixup,sct,indx,nxh)
# calculate form factor: ffc
isign = 0
pois1(qe,fxe,isign,ffc,ax,affp,we,nx)
# initialize electrons
distr1h(part,vtx,vty,vtz,vx0,vy0,vz0,npx,idimp,np,nx,ipbc)

# find maximum and minimum initial electron density
qe.fill(0.0)
gpost1l(part,qe,qme,np,idimp,nxe)
aguard1l(qe,nx,nxe)
fwpminmx1(qe,qbme,wpmax,wpmin,nx,nxe)
wpm = 0.5*(wpmax[0] + wpmin[0])*affp
# accelerate convergence: update wpm
if (wpm <= 10.0):
   wpm = 0.75*wpm
print "wpm=",wpm
q2m0 = wpm/affp
# calculate form factor: ffe
isign = 0
epois13(dcu,cus,isign,ffe,ax,affp,wpm,ci,wf,nx,nxeh,nxh)

# initialize transverse electric field
cus.fill(0.0)

# * * * start main iteration loop * * *

for ntime in xrange(0,nloop):
#  print "ntime = ", ntime

# deposit current with standard procedure: updates cue
   dtimer(dtime,itime,-1)
   cue.fill(0.0)
   gjpost1l(part,cue,qme,zero,np,idimp,nx,nxe,ipbc)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tdjpost = tdjpost + time

# deposit charge with standard procedure: updates qe
   dtimer(dtime,itime,-1)
   qe.fill(0.0)
   gpost1l(part,qe,qme,np,idimp,nxe)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tdpost = tdpost + time

# add guard cells with standard procedure: updates qe, cue
   dtimer(dtime,itime,-1)
   aguard1l(qe,nx,nxe)
   acguard1l(cue,nx,nxe)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tguard = tguard + time

# transform charge to fourier space with standard procedure:
# updates qe, gxe
   dtimer(dtime,itime,-1)
   isign = -1
   fft1rxx(qe,gxe,isign,mixup,sct,indx,nxe,nxh)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tfft = tfft + time

# calculate longitudinal force/charge in fourier space with standard
# procedure: updates fxe, we
   dtimer(dtime,itime,-1)
   isign = -1
   pois1(qe,fxe,isign,ffc,ax,affp,we,nx)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tfield = tfield + time

# transform longitudinal electric force to real space with standard
# procedure: updates fxe, gxe
   dtimer(dtime,itime,-1)
   isign = 1
   fft1rxx(fxe,gxe,isign,mixup,sct,indx,nxe,nxh)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tfft = tfft + time

# transform current to fourier space with standard procedure:
# updates cue, gyze
   dtimer(dtime,itime,-1)
   isign = -1
   fft1r2x(cue,gyze,isign,mixup,sct,indx,nxe,nxh)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tfft = tfft + time

# calculate magnetic field in fourier space with standard procedure:
# updates byze, wm
   dtimer(dtime,itime,-1)
   bbpois13(cue,byze,ffc,ci,wm,nx,nxeh,nxh)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tfield = tfield + time

# transform magnetic force to real space with standard procedure:
# updates byze, gyze
   dtimer(dtime,itime,-1)
   isign = 1
   fft1r2x(byze,gyze,isign,mixup,sct,indx,nxe,nxh)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tfft = tfft + time

# add constant to magnetic field with standard procedure: updates byze
   dtimer(dtime,itime,-1)
   baddext1(byze,omy,omz,nx,nxe)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tfield = tfield + time

# copy guard cells with standard procedure: updates fxe, byze
   dtimer(dtime,itime,-1)
   dguard1l(fxe,nx,nxe)
   cguard1l(byze,nx,nxe)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tguard = tguard + time

# add longitudinal and old transverse electric fields with standard
# procedure: updates exyze
   dtimer(dtime,itime,-1)
   addvrfield13(exyze,cus,fxe,nxe)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tfield = tfield + time

# deposit electron acceleration density and momentum flux with standard
# procedure: updates dcu, amu
   dtimer(dtime,itime,-1)
   dcu.fill(0.0); amu.fill(0.0)
   gdjpost1l(part,exyze,byze,dcu,amu,omx,qme,qbme,dt,idimp,np,nxe)
# add old scaled electric field with standard procedure: updates dcu
   ascfguard1l(dcu,cus,q2m0,nx,nxe)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tdcjpost = tdcjpost + time

# add guard cells with standard procedure: updates dcu, amu
   dtimer(dtime,itime,-1)
   acguard1l(dcu,nx,nxe)
   acguard1l(amu,nx,nxe)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tguard = tguard + time

# transform acceleration density and momentum flux to fourier space
# with standard procedure: updates dcu, amu, gyze
   dtimer(dtime,itime,-1)
   isign = -1
   fft1r2x(dcu,gyze,isign,mixup,sct,indx,nxe,nxh)
   fft1r2x(amu,gyze,isign,mixup,sct,indx,nxe,nxh)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tfft = tfft + time

# take transverse part of time derivative of current with standard
# procedure: updates dcu
   dtimer(dtime,itime,-1)
   adcuperp13(dcu,amu,nx,nxeh)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tfield = tfield + time

# calculate transverse electric field with standard procedure:
# updates cus, wf
   dtimer(dtime,itime,-1)
   isign = -1
   epois13(dcu,cus,isign,ffe,ax,affp,wpm,ci,wf,nx,nxeh,nxh)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tfield = tfield + time

# transform transverse electric field to real space with standard
# procedure: updates cus, gyze
   dtimer(dtime,itime,-1)
   isign = 1
   fft1r2x(cus,gyze,isign,mixup,sct,indx,nxe,nxh)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tfft = tfft + time

# copy guard cells with standard procedure: updates cus
   dtimer(dtime,itime,-1)
   cguard1l(cus,nx,nxe)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tguard = tguard + time

# add longitudinal and transverse electric fields with standard
# procedure: exyze = cus + fxe, updates exyze
# cus needs to be retained for next time step
   dtimer(dtime,itime,-1)
   addvrfield13(exyze,cus,fxe,nxe)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tfield = tfield + time

# inner iteration loop
   for k in xrange(0,ndc):

# deposit electron current and acceleration density and momentum flux
# with standard procedure: updates cue, dcu, amu
      dtimer(dtime,itime,-1)
      cue.fill(0.0); dcu.fill(0.0); amu.fill(0.0)
      gdcjpost1l(part,exyze,byze,cue,dcu,amu,omx,qme,qbme,dt,idimp,np,
                 nxe)
# add scaled electric field with standard procedure: updates dcu
      ascfguard1l(dcu,cus,q2m0,nx,nxe)
      dtimer(dtime,itime,1)
      time = float(dtime)
      tdcjpost = tdcjpost + time

# add guard cells for current, acceleration density, and momentum flux
# with standard procedure: updates cue, dcu, amu
      dtimer(dtime,itime,-1)
      acguard1l(cue,nx,nxe)
      acguard1l(dcu,nx,nxe)
      acguard1l(amu,nx,nxe)
      dtimer(dtime,itime,1)
      time = float(dtime)
      tguard = tguard + time

# transform current to fourier space with standard procedure:
# update cue, gyze
      dtimer(dtime,itime,-1)
      isign = -1
      fft1r2x(cue,gyze,isign,mixup,sct,indx,nxe,nxh)
      dtimer(dtime,itime,1)
      time = float(dtime)
      tfft = tfft + time

# calculate magnetic field in fourier space with standard procedure:
# updates byze, wm
      dtimer(dtime,itime,-1)
      bbpois13(cue,byze,ffc,ci,wm,nx,nxeh,nxh)
      dtimer(dtime,itime,1)
      time = float(dtime)
      tfield = tfield + time

# transform magnetic force to real space with standard procedure:
# updates byze, gyze
      dtimer(dtime,itime,-1)
      isign = 1
      fft1r2x(byze,gyze,isign,mixup,sct,indx,nxe,nxh)
      dtimer(dtime,itime,1)
      time = float(dtime)
      tfft = tfft + time

# add constant to magnetic field with standard procedure: updates bzye
      dtimer(dtime,itime,-1)
      baddext1(byze,omy,omz,nx,nxe)
      dtimer(dtime,itime,1)
      time = float(dtime)
      tfield = tfield + time

# transform acceleration density and momentum flux to fourier space
# with standard procedure: updates dcu, amu, gyze
      dtimer(dtime,itime,-1)
      isign = -1
      fft1r2x(dcu,gyze,isign,mixup,sct,indx,nxe,nxh)
      fft1r2x(amu,gyze,isign,mixup,sct,indx,nxe,nxh)
      dtimer(dtime,itime,1)
      time = float(dtime)
      tfft = tfft + time

# take transverse part of time derivative of current with standard
# procedure: updates dcu
      dtimer(dtime,itime,-1)
      adcuperp13(dcu,amu,nx,nxeh)
      dtimer(dtime,itime,1)
      time = float(dtime)
      tfield = tfield + time

# calculate transverse electric field with standard procedure:
# updates cus, wf
      dtimer(dtime,itime,-1)
      isign = -1
      epois13(dcu,cus,isign,ffe,ax,affp,wpm,ci,wf,nx,nxeh,nxh)
      dtimer(dtime,itime,1)
      time = float(dtime)
      tfield = tfield + time

# transform transverse electric field to real space with standard
      dtimer(dtime,itime,-1)
      isign = 1
      fft1r2x(cus,gyze,isign,mixup,sct,indx,nxe,nxh)
      dtimer(dtime,itime,1)
      time = float(dtime)
      tfft = tfft + time

# copy guard cells with standard procedure: updates byze, cus
      dtimer(dtime,itime,-1)
      cguard1l(byze,nx,nxe)
      cguard1l(cus,nx,nxe)
      dtimer(dtime,itime,1)
      time = float(dtime)
      tguard = tguard + time

# add longitudinal and transverse electric fields with standard
# procedure: exyze = cus + fxyze, updates exyze
# cus needs to be retained for next time step
      dtimer(dtime,itime,-1)
      addvrfield13(exyze,cus,fxe,nxe)
      dtimer(dtime,itime,1)
      time = float(dtime)
      tfield = tfield + time
   pass

# push particles with standard procedure: updates part, wke
   wke[0] = 0.0
   dtimer(dtime,itime,-1)
   gbpush13l(part,exyze,byze,omx,qbme,dt,dt,wke,idimp,np,nx,nxe,ipbc)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tpush = tpush + time

# sort particles by cell for standard procedure
   if (sortime > 0):
      if (ntime%sortime==0):
         dtimer(dtime,itime,-1)
         dsortp1xl(part,part2,npic,idimp,np,nx1)
         tpart = part
         part = part2
         part2 = tpart
         dtimer(dtime,itime,1)
         time = float(dtime)
         tsort = tsort + time
      pass
   pass

   if (ntime==0):
      wt = we + wm
      print "Initial Total Field, Kinetic and Total Energies:"
      print "%14.7e %14.7e %14.7e" % (wt, wke, wke + wt)
      print "Initial Electrostatic, Transverse Electric and Magnetic " \
      "Field Energies:"
      print "%14.7e %14.7e %14.7e" % (we, wf, wm)
ntime = ntime + 1

# * * * end main iteration loop * * *

print "ntime, ndc = ", ntime, ndc
wt = we + wm
print "Final Total Field, Kinetic and Total Energies:"
print "%14.7e %14.7e %14.7e" % (wt, wke, wke + wt)
print "Final Electrostatic, Transverse Electric and Magnetic Field " \
"Energies:"
print "%14.7e %14.7e %14.7e" % (we, wf, wm)

print ""
print "deposit time = ", tdpost
print "current deposit time = ", tdjpost
print "current derivative deposit time = ", tdcjpost
tdpost = tdpost + tdjpost + tdcjpost
print "total deposit time = ", tdpost
print "guard time = ", tguard
print "solver time = ", tfield
print "fft time = ", tfft
print "push time = ", tpush
print "sort time = ", tsort
tfield = tfield + tguard + tfft
print "total solver time = ", tfield
time = tdpost + tpush + tsort
print "total particle time = ", time
wt = time + tfield
print "total time = ", wt
print ""

wt = 1.0e+09/(float(nloop)*float(np))
print "Push Time (nsec) = ", tpush*wt
print "Deposit Time (nsec) = ", tdpost*wt
print "Sort Time (nsec) = ", tsort*wt
print "Total Particle Time (nsec) = ", time*wt

