#-----------------------------------------------------------------------
# Skeleton 1-2/2D Electromagnetic PIC code
# written by Viktor K. Decyk, Adam Tableman, and Qiyang Hu, UCLA
import math
import numpy
from fbpush1 import *
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
tend = 10.0; dt = 0.05; qme = -1.0
# vtx/vty = thermal velocity of electrons in x/y direction
# vx0/vy0 = drift velocity of electrons in x/y direction.
vtx = 1.0; vty = 1.0; vx0 = 0.0; vy0 = 0.0
# vtx/vz0 = thermal/drift velocity of electrons in z direction
vtz = 1.0; vz0 = 0.0
# omx = magnetic field electron cyclotron frequency in x
omx = 0.0
# ax = smoothed particle size in x direction
# ci = reciprocal of velocity of light.
ax = .912871; ci = 0.1
# idimp = number of particle coordinates = 4
# ipbc = particle boundary condition: 1 = periodic
# sortime = number of time steps between standard electron sorting
# relativity = (no,yes) = (0,1) = relativity is used
idimp = 4; ipbc = 1; sortime = 50; relativity = 1
# wke/we = particle kinetic/electrostatic field energy
# wf/wm/wt = magnetic field/transverse electric field/total energy
wke = numpy.zeros((1),float_type)
we = numpy.zeros((1),float_type)
wf = numpy.zeros((1),float_type)
wm = numpy.zeros((1),float_type)
wt = numpy.zeros((1),float_type)

# declare and initialize timing data
itime = numpy.empty((4),numpy.int32)
tdpost = 0.0; tguard = 0.0; tfft = 0.0; tfield = 0.0
tdjpost = 0.0; tpush = 0.0; tsort = 0.0
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
dth = 0.0

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
# fxyze = smoothed electric field with guard cells
fxyze = numpy.empty((3,nxe),float_type,'F')
# byze = smoothed magnetic field with guard cells
byze = numpy.empty((2,nxe),float_type,'F')
# eyz = transverse electric field in fourier space
eyz = numpy.empty((2,nxeh),complex_type,'F')
# byz = magnetic field in fourier space
byz = numpy.empty((2,nxeh),complex_type,'F')
# ffc = form factor array for poisson solver
ffc = numpy.empty((nxh),complex_type,'F')
# mixup = bit reverse table for FFT
mixup = numpy.empty((nxh),int_type,'F')
# sct = sine/cosine table for FFT
sct = numpy.empty((nxh),complex_type,'F')
# npic = scratch array for reordering particles
npic = numpy.empty((nx1),int_type,'F')
# gxyze = scratch array for fft
gxyze = numpy.empty((3,nxe),float_type,'F')

# prepare fft tables
wfft1rinit(mixup,sct,indx,nxh)
# calculate form factors
isign = 0
pois1(qe,fxe,isign,ffc,ax,affp,we,nx)
# initialize electrons
distr1h(part,vtx,vty,vtz,vx0,vy0,vz0,npx,idimp,np,nx,ipbc)

# initialize transverse electromagnetic fields
eyz.fill(numpy.complex(0.0,0.0))
byz.fill(numpy.complex(0.0,0.0))

if (dt > 0.64*ci):
   print "Warning: Courant condition may be exceeded!"

# * * * start main iteration loop * * *

for ntime in xrange(0,nloop):
#  print "ntime = ", ntime

# deposit current with standard procedure: updates part, cue
   dtimer(dtime,itime,-1)
   cue.fill(0.0)
   if (relativity==1):
      grjpost1l(part,cue,qme,dth,ci,np,idimp,nx,nxe,ipbc)
   else:
      gjpost1l(part,cue,qme,dth,np,idimp,nx,nxe,ipbc)
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

# add guard cells with standard procedure: updates cue, qe
   dtimer(dtime,itime,-1)
   acguard1l(cue,nx,nxe)
   aguard1l(qe,nx,nxe)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tguard = tguard + time

# transform charge to fourier space with standard procedure:
# updates qe, fxe
   dtimer(dtime,itime,-1)
   isign = -1
   fft1rxx(qe,fxe,isign,mixup,sct,indx,nxe,nxh)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tfft = tfft + time

# transform current to fourier space with standard procedure:
# updates cue, byze
   dtimer(dtime,itime,-1)
   isign = -1
   fft1r2x(cue,byze,isign,mixup,sct,indx,nxe,nxh)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tfft = tfft + time

# calculate electromagnetic fields in fourier space with standard
# procedure: updates eyz, byz
   dtimer(dtime,itime,-1)
   if (ntime==0):
      ibpois13(cue,byz,ffc,ci,wm,nx,nxeh,nxh)
      wf[0] = 0.0
      dth = 0.5*dt
   else:
      maxwel1(eyz,byz,cue,ffc,ci,dt,wf,wm,nx,nxeh,nxh)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tfield = tfield + time

# calculate force/charge in fourier space with standard procedure:
# updates fxe, we
   dtimer(dtime,itime,-1)
   isign = -1
   pois1(qe,fxe,isign,ffc,ax,affp,we,nx)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tfield = tfield + time

# add longitudinal and transverse electric fields with standard
# procedure: updates fxyze
   dtimer(dtime,itime,-1)
   isign = 1
   emfield1(fxyze,fxe,eyz,ffc,nx,nxeh,nxh)
# copy magnetic field with standard procedure: updates byze
   isign = -1
   bmfield1(byze,byz,ffc,nx,nxeh,nxh)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tfield = tfield + time

# transform electric force to real space with standard procedure:
# updates fxyze, gxyze
   dtimer(dtime,itime,-1)
   isign = 1
   fft1r3x(fxyze,gxyze,isign,mixup,sct,indx,nxe,nxh)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tfft = tfft + time

# transform magnetic force to real space with standard procedure:
# updates byze, cue
   dtimer(dtime,itime,-1)
   isign = 1
   fft1r2x(byze,cue,isign,mixup,sct,indx,nxe,nxh)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tfft = tfft + time

# copy guard cells with standard procedure: updates fxyze, byze
   dtimer(dtime,itime,-1)
   bguard1l(fxyze,nx,nxe)
   cguard1l(byze,nx,nxe)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tguard = tguard + time

# push particles with standard procedure: updates part, wke
   wke[0] = 0.0
   dtimer(dtime,itime,-1)
   if (relativity==1):
      grbpush13l(part,fxyze,byze,omx,qbme,dt,dth,ci,wke,idimp,np,nx,nxe,
                 ipbc)
   else:
      gbpush13l(part,fxyze,byze,omx,qbme,dt,dth,wke,idimp,np,nx,nxe,
                ipbc)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tpush = tpush + time

# sort particles by cell for standard procedure
   if (sortime > 0):
      if (ntime%sortime==0):
         dtimer(dtime,itime,-1)
         dsortp1xl(part,part2,npic,idimp,np,nx1)
# exchange pointers
         tpart = part
         part = part2
         part2 = tpart
         dtimer(dtime,itime,1)
         time = float(dtime)
         tsort = tsort + time
      pass
   pass

   if (ntime==0):
      wt = we + wf + wm
      print "Initial Total Field, Kinetic and Total Energies:"
      print "%14.7e %14.7e %14.7e" % (wt, wke, wke + wt)
      print "Initial Electrostatic, Transverse Electric and Magnetic " \
      "Field Energies:"
      print "%14.7e %14.7e %14.7e" % (we, wf, wm)
ntime = ntime + 1

# * * * end main iteration loop * * *

print "ntime, relativity = ", ntime, relativity
wt = we + wf + wm
print "Final Total Field, Kinetic and Total Energies:"
print "%14.7e %14.7e %14.7e" % (wt, wke, wke + wt)
print "Final Electrostatic, Transverse Electric and Magnetic Field " \
"Energies:"
print "%14.7e %14.7e %14.7e" % (we, wf, wm)

print ""
print "deposit time = ", tdpost
print "current deposit time = ", tdjpost
tdpost = tdpost + tdjpost
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

