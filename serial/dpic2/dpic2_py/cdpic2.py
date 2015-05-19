#-----------------------------------------------------------------------
# Skeleton 2-1/2D Darwin PIC code
# written by Viktor K. Decyk, Adam Tableman, and Qiyang Hu, UCLA
import math
import numpy
from cdpush2 import *
from dtimer import *

int_type = numpy.int32
double_type = numpy.float64
float_type = numpy.float32
complex_type = numpy.complex64

# indx/indy = exponent which determines grid points in x/y direction:
# nx = 2**indx, ny = 2**indy.
indx =   9; indy =   9
# npx/npy = number of electrons distributed in x/y direction.
npx =  3072; npy =   3072
# ndim = number of velocity coordinates = 3
ndim = 3
# tend = time at end of simulation, in units of plasma frequency.
# dt = time interval between successive calculations.
# qme = charge on electron, in units of e.
tend = 10.0; dt = 0.1; qme = -1.0
# vtx/vty = thermal velocity of electrons in x/y direction
# vx0/vy0 = drift velocity of electrons in x/y direction.
vtx = 1.0; vty = 1.0; vx0 = 0.0; vy0 = 0.0
# vtx/vz0 = thermal/drift velocity of electrons in z direction
vtz = 1.0; vz0 = 0.0
# ax/ay = smoothed particle size in x/y direction
# ci = reciprocal of velocity of light.
ax = .912871; ay = .912871; ci = 0.1
# idimp = number of particle coordinates = 5
# ipbc = particle boundary condition: 1 = periodic
# sortime = number of time steps between standard electron sorting
idimp = 5; ipbc = 1; sortime = 50
# omx/omy/omz = magnetic field electron cyclotron frequency in x/y/z 
omx = 0.4; omy = 0.0; omz = 0.0
# ndc = number of corrections in darwin iteration
ndc = 1
# wke/we/wt = particle kinetic/electric field/total energy
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
# nx/ny = number of grid points in x/y direction
np = npx*npy; nx = int(math.pow(2,indx)); ny = int(math.pow(2,indy))
nxh = int(nx/2); nyh = max(1,int(ny/2))
nxe = nx + 2; nye = ny + 1; nxeh = int(nxe/2)
nxyh = int(max(nx,ny)/2); nxhy = max(nxh,ny); ny1 = ny + 1
# nloop = number of time steps in simulation
# ntime = current time step
nloop = int(tend/dt + .0001); ntime = 0
# mdim = dimension of amu array
mdim = 2*ndim - 2
qbme = qme
affp = float(nx*ny)/float(np)

# allocate data for standard code
# part, part2 = particle arrays
part = numpy.empty((idimp,np),float_type,'F')
if (sortime > 0):
   part2 = numpy.empty((idimp,np),float_type,'F')
# qe = electron charge density with guard cells
qe = numpy.empty((nxe,nye),float_type,'F')
# cue = electron current density with guard cells
cue = numpy.empty((ndim,nxe,nye),float_type,'F')
# dcu = acceleration density with guard cells
dcu = numpy.empty((ndim,nxe,nye),float_type,'F')
# cus = smoothed transverse electric field with guard cells
cus = numpy.empty((ndim,nxe,nye),float_type,'F')
# amu = momentum flux with guard cells
amu = numpy.empty((mdim,nxe,nye),float_type,'F')
# exyze = smoothed total electric field with guard cells
exyze = numpy.empty((ndim,nxe,nye),float_type,'F')
# fxyze = smoothed longitudinal electric field with guard cells
fxyze = numpy.empty((ndim,nxe,nye),float_type,'F')
# bxyze = smoothed magnetic field with guard cells
bxyze = numpy.empty((ndim,nxe,nye),float_type,'F')
# ffc, ffe = form factor arrays for poisson solvers
ffc = numpy.empty((nxh,nyh),complex_type,'F')
ffe = numpy.empty((nxh,nyh),complex_type,'F')
# mixup = bit reverse table for FFT
mixup = numpy.empty((nxhy),int_type,'F')
# sct = sine/cosine table for FFT
sct = numpy.empty((nxyh),complex_type,'F')
# npicy = scratch array for reordering particles
npicy = numpy.empty((ny1),int_type,'F')
# ss = scratch array for WFFT2RN
ss = numpy.empty((mdim,nxeh),complex_type,'F')

# prepare fft tables
cwfft2rinit(mixup,sct,indx,indy,nxhy,nxyh)
# calculate form factors: ffc
isign = 0
cpois23(qe,fxyze,isign,ffc,ax,ay,affp,we,nx,ny,nxeh,nye,nxh,nyh)
# initialize electrons
cdistr2h(part,vtx,vty,vtz,vx0,vy0,vz0,npx,npy,idimp,np,nx,ny,ipbc)

# find maximum and minimum initial electron density
qe.fill(0.0)
cgpost2l(part,qe,qme,np,idimp,nxe,nye)
caguard2l(qe,nx,ny,nxe,nye)
cfwpminmx2(qe,qbme,wpmax,wpmin,nx,ny,nxe,nye)
wpm = 0.5*(wpmax[0] + wpmin[0])*affp
# accelerate convergence: update wpm
if (wpm <= 10.0):
   wpm = 0.75*wpm
print "wpm=",wpm
q2m0 = wpm/affp
# calculate form factor: ffe
isign = 0
cepois23(dcu,cus,isign,ffe,ax,ay,affp,wpm,ci,wf,nx,ny,nxeh,nye,
         nxh,nyh)

# initialize transverse electric field
cus.fill(0.0)

# * * * start main iteration loop * * *

for ntime in xrange(0,nloop):
#  print "ntime = ", ntime

# deposit current with standard procedure: updates cue
   dtimer(dtime,itime,-1)
   cue.fill(0.0)
   cgjpost2l(part,cue,qme,zero,np,idimp,nx,ny,nxe,nye,ipbc)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tdjpost = tdjpost + time

# deposit charge with standard procedure: updates qe
   dtimer(dtime,itime,-1)
   qe.fill(0.0)
   cgpost2l(part,qe,qme,np,idimp,nxe,nye)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tdpost = tdpost + time

# add guard cells with standard procedure: updates qe, cue
   dtimer(dtime,itime,-1)
   caguard2l(qe,nx,ny,nxe,nye)
   cacguard2l(cue,nx,ny,nxe,nye)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tguard = tguard + time

# transform charge to fourier space with standard procedure: updates qe
   dtimer(dtime,itime,-1)
   isign = -1
   cwfft2rx(qe,isign,mixup,sct,indx,indy,nxeh,nye,nxhy,nxyh)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tfft = tfft + time

# calculate longitudinal force/charge in fourier space with standard
# procedure: updates fxyze, we
   dtimer(dtime,itime,-1)
   isign = -1
   cpois23(qe,fxyze,isign,ffc,ax,ay,affp,we,nx,ny,nxeh,nye,nxh,nyh)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tfield = tfield + time

# transform longitudinal electric force to real space with standard
# procedure: updates fxyze
   dtimer(dtime,itime,-1)
   isign = 1
   cwfft2r3(fxyze,isign,mixup,sct,indx,indy,nxeh,nye,nxhy,nxyh)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tfft = tfft + time

# transform current to fourier space with standard procedure: update cue
   dtimer(dtime,itime,-1)
   isign = -1
   cwfft2r3(cue,isign,mixup,sct,indx,indy,nxeh,nye,nxhy,nxyh)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tfft = tfft + time

# take transverse part of current with standard procedure: updates cue
   dtimer(dtime,itime,-1)
   ccuperp2(cue,nx,ny,nxeh,nye)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tfield = tfield + time

# calculate magnetic field in fourier space with standard procedure:
# updates bxyze, wm
   dtimer(dtime,itime,-1)
   cbbpois23(cue,bxyze,ffc,ci,wm,nx,ny,nxeh,nye,nxh,nyh)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tfield = tfield + time

# transform magnetic force to real space with standard procedure:
# updates bxyze
   dtimer(dtime,itime,-1)
   isign = 1
   cwfft2r3(bxyze,isign,mixup,sct,indx,indy,nxeh,nye,nxhy,nxyh)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tfft = tfft + time

# add constant to magnetic field with standard procedure: updates bxyze
   dtimer(dtime,itime,-1)
   cbaddext2(bxyze,omx,omy,omz,nx,ny,nxe,nye)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tfield = tfield + time

# copy guard cells with standard procedure: updates fxyze, bxyze
   dtimer(dtime,itime,-1)
   cbguard2l(fxyze,nx,ny,nxe,nye)
   cbguard2l(bxyze,nx,ny,nxe,nye)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tguard = tguard + time

# add longitudinal and old transverse electric fields with standard
# procedure: updates exyze
   dtimer(dtime,itime,-1)
   caddvrfield2(exyze,cus,fxyze,ndim,nxe,nye)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tfield = tfield + time

# deposit electron acceleration density and momentum flux with standard
# procedure: updates dcu, amu
   dtimer(dtime,itime,-1)
   dcu.fill(0.0); amu.fill(0.0)
   cgdjpost2l(part,exyze,bxyze,dcu,amu,qme,qbme,dt,idimp,np,nxe,nye)
# add old scaled electric field with standard procedure: updates dcu
   cascfguard2l(dcu,cus,q2m0,nx,ny,nxe,nye)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tdcjpost = tdcjpost + time

# add guard cells with standard procedure: updates dcu, amu
   dtimer(dtime,itime,-1)
   cacguard2l(dcu,nx,ny,nxe,nye)
   camcguard2l(amu,nx,ny,nxe,nye,mdim)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tguard = tguard + time

# transform acceleration density and momentum flux to fourier space
# with standard procedure: updates dcu, amu
   dtimer(dtime,itime,-1)
   isign = -1
   cwfft2r3(dcu,isign,mixup,sct,indx,indy,nxeh,nye,nxhy,nxyh)
   cwfft2rn(amu,ss,isign,mixup,sct,indx,indy,nxeh,nye,mdim,nxhy,nxyh)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tfft = tfft + time

# take transverse part of time derivative of current with standard
# procedure: updates dcu
   dtimer(dtime,itime,-1)
   cadcuperp23(dcu,amu,nx,ny,nxeh,nye)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tfield = tfield + time

# calculate transverse electric field with standard procedure:
# updates cus, wf
   dtimer(dtime,itime,-1)
   isign = -1
   cepois23(dcu,cus,isign,ffe,ax,ay,affp,wpm,ci,wf,nx,ny,nxeh,nye,nxh,
            nyh)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tfield = tfield + time

# transform transverse electric field to real space with standard
# procedure: updates cus
   dtimer(dtime,itime,-1)
   isign = 1
   cwfft2r3(cus,isign,mixup,sct,indx,indy,nxeh,nye,nxhy,nxyh)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tfft = tfft + time

# copy guard cells with standard procedure: updates cus
   dtimer(dtime,itime,-1)
   cbguard2l(cus,nx,ny,nxe,nye)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tguard = tguard + time

# add longitudinal and transverse electric fields with standard
# procedure: exyze = cus + fxyze, updates exyze
# cus needs to be retained for next time step
   dtimer(dtime,itime,-1)
   caddvrfield2(exyze,cus,fxyze,ndim,nxe,nye)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tfield = tfield + time

# inner iteration loop
   for k in xrange(0,ndc):

# deposit electron current and acceleration density and momentum flux
# with standard procedure: updates cue, dcu, amu
      dtimer(dtime,itime,-1)
      cue.fill(0.0); dcu.fill(0.0); amu.fill(0.0)
      cgdcjpost2l(part,exyze,bxyze,cue,dcu,amu,qme,qbme,dt,idimp,np,nxe,\
                 nye)
# add scaled electric field with standard procedure: updates dcu
      cascfguard2l(dcu,cus,q2m0,nx,ny,nxe,nye)
      dtimer(dtime,itime,1)
      time = float(dtime)
      tdcjpost = tdcjpost + time

# add guard cells for current, acceleration density, and momentum flux
# with standard procedure: updates cue, dcu, amu
      dtimer(dtime,itime,-1)
      cacguard2l(cue,nx,ny,nxe,nye)
      cacguard2l(dcu,nx,ny,nxe,nye)
      camcguard2l(amu,nx,ny,nxe,nye,mdim)
      dtimer(dtime,itime,1)
      time = float(dtime)
      tguard = tguard + time

# transform current to fourier space with standard procedure: update cue
      dtimer(dtime,itime,-1)
      isign = -1
      cwfft2r3(cue,isign,mixup,sct,indx,indy,nxeh,nye,nxhy,nxyh)
      dtimer(dtime,itime,1)
      time = float(dtime)
      tfft = tfft + time

# take transverse part of current with standard procedure: updates cue
      dtimer(dtime,itime,-1)
      ccuperp2(cue,nx,ny,nxeh,nye)
      dtimer(dtime,itime,1)
      time = float(dtime)
      tfield = tfield + time

# calculate magnetic field in fourier space with standard procedure:
# updates bxyze, wm
      dtimer(dtime,itime,-1)
      cbbpois23(cue,bxyze,ffc,ci,wm,nx,ny,nxeh,nye,nxh,nyh)
      dtimer(dtime,itime,1)
      time = float(dtime)
      tfield = tfield + time

# transform magnetic force to real space with standard procedure:
# updates bxyze
      dtimer(dtime,itime,-1)
      isign = 1
      cwfft2r3(bxyze,isign,mixup,sct,indx,indy,nxeh,nye,nxhy,nxyh)
      dtimer(dtime,itime,1)
      time = float(dtime)
      tfft = tfft + time

# add constant to magnetic field with standard procedure: updates bxzye
      dtimer(dtime,itime,-1)
      cbaddext2(bxyze,omx,omy,omz,nx,ny,nxe,nye)
      dtimer(dtime,itime,1)
      time = float(dtime)
      tfield = tfield + time

# transform acceleration density and momentum flux to fourier space
# with standard procedure: updates dcu and amu
      dtimer(dtime,itime,-1)
      isign = -1
      cwfft2r3(dcu,isign,mixup,sct,indx,indy,nxeh,nye,nxhy,nxyh)
      cwfft2rn(amu,ss,isign,mixup,sct,indx,indy,nxeh,nye,mdim,nxhy,nxyh)
      dtimer(dtime,itime,1)
      time = float(dtime)
      tfft = tfft + time

# take transverse part of time derivative of current with standard
# procedure: updates dcu
      dtimer(dtime,itime,-1)
      cadcuperp23(dcu,amu,nx,ny,nxeh,nye)
      dtimer(dtime,itime,1)
      time = float(dtime)
      tfield = tfield + time

# calculate transverse electric field with standard procedure:
# updates cus, wf
      dtimer(dtime,itime,-1)
      isign = -1
      cepois23(dcu,cus,isign,ffe,ax,ay,affp,wpm,ci,wf,nx,ny,nxeh,nye,
               nxh,nyh)
      dtimer(dtime,itime,1)
      time = float(dtime)
      tfield = tfield + time

# transform transverse electric field to real space with standard
# procedure: updates cus
      dtimer(dtime,itime,-1)
      isign = 1
      cwfft2r3(cus,isign,mixup,sct,indx,indy,nxeh,nye,nxhy,nxyh)
      dtimer(dtime,itime,1)
      time = float(dtime)
      tfft = tfft + time

# copy guard cells with standard procedure: updates bxyze, cus
      dtimer(dtime,itime,-1)
      cbguard2l(bxyze,nx,ny,nxe,nye)
      cbguard2l(cus,nx,ny,nxe,nye)
      dtimer(dtime,itime,1)
      time = float(dtime)
      tguard = tguard + time

# add longitudinal and transverse electric fields with standard
# procedure: exyze = cus + fxyze, updates exyze
# cus needs to be retained for next time step
      dtimer(dtime,itime,-1)
      caddvrfield2(exyze,cus,fxyze,ndim,nxe,nye)
      dtimer(dtime,itime,1)
      time = float(dtime)
      tfield = tfield + time
   pass

# push particles with standard procedure: updates part, wke
   wke[0] = 0.0
   dtimer(dtime,itime,-1)
   cgbpush23l(part,exyze,bxyze,qbme,dt,dt,wke,idimp,np,nx,ny,nxe,nye,
              ipbc)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tpush = tpush + time

# sort particles by cell for standard procedure
   if (sortime > 0):
      if (ntime%sortime==0):
         dtimer(dtime,itime,-1)
         cdsortp2yl(part,part2,npicy,idimp,np,ny1)
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
