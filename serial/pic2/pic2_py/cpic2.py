#-----------------------------------------------------------------------
# Skeleton 2D Electrostatic PIC code
# written by Viktor K. Decyk, Adam Tableman, and Qiyang Hu, UCLA
import math
import numpy
from cpush2 import *
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
# ndim = number of velocity coordinates = 2
ndim = 2
# tend = time at end of simulation, in units of plasma frequency.
# dt = time interval between successive calculations.
# qme = charge on electron, in units of e.
tend = 10.0; dt = 0.1; qme = -1.0
# vtx/vty = thermal velocity of electrons in x/y direction
# vx0/vy0 = drift velocity of electrons in x/y direction.
vtx = 1.0; vty = 1.0; vx0 = 0.0; vy0 = 0.0
# ax/ay = smoothed particle size in x/y direction
ax = .912871; ay = .912871
# idimp = number of particle coordinates = 4
# ipbc = particle boundary condition: 1 = periodic
# sortime = number of time steps between standard electron sorting
idimp = 4; ipbc = 1; sortime = 50
# wke/we/wt = particle kinetic/electric field/total energy
wke = numpy.zeros((1),float_type)
we = numpy.zeros((1),float_type)
wt = numpy.zeros((1),float_type)

# declare and initialize timing data
itime = numpy.empty((4),numpy.int32)
tdpost = 0.0; tguard = 0.0; tfft = 0.0; tfield = 0.0
tpush = 0.0; tsort = 0.0
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
qbme = qme
affp = float(nx*ny)/float(np)

# allocate data for standard code
# part, part2 = particle arrays
part = numpy.empty((idimp,np),float_type,'F')
if (sortime > 0):
   part2 = numpy.empty((idimp,np),float_type,'F')
# qe = electron charge density with guard cells
qe = numpy.empty((nxe,nye),float_type,'F')
# fxye = smoothed electric field with guard cells
fxye = numpy.empty((ndim,nxe,nye),float_type,'F')
# ffc = form factor array for poisson solver
ffc = numpy.empty((nxh,nyh),complex_type,'F')
# mixup = bit reverse table for FFT
mixup = numpy.empty((nxhy),int_type,'F')
# sct = sine/cosine table for FFT
sct = numpy.empty((nxyh),complex_type,'F')
# npicy = scratch array for reordering particles
npicy = numpy.empty((ny1),int_type,'F')

# prepare fft tables
cwfft2rinit(mixup,sct,indx,indy,nxhy,nxyh)
# calculate form factors
isign = 0
cpois22(qe,fxye,isign,ffc,ax,ay,affp,we,nx,ny,nxeh,nye,nxh,nyh)
# initialize electrons
cdistr2(part,vtx,vty,vx0,vy0,npx,npy,idimp,np,nx,ny,ipbc)

# * * * start main iteration loop * * *

for ntime in xrange(0,nloop):
#  print "ntime = ", ntime

# deposit charge with standard procedure: updates qe
   dtimer(dtime,itime,-1)
   qe.fill(0.0)
   cgpost2l(part,qe,qme,np,idimp,nxe,nye)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tdpost = tdpost + time

# add guard cells with standard procedure: updates qe
   dtimer(dtime,itime,-1)
   caguard2l(qe,nx,ny,nxe,nye)
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

# calculate force/charge in fourier space with standard procedure:
# updates fxye, we
   dtimer(dtime,itime,-1)
   isign = -1
   cpois22(qe,fxye,isign,ffc,ax,ay,affp,we,nx,ny,nxeh,nye,nxh,nyh)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tfield = tfield + time

# transform force to real space with standard procedure: updates fxye
   dtimer(dtime,itime,-1)
   isign = 1
   cwfft2r2(fxye,isign,mixup,sct,indx,indy,nxeh,nye,nxhy,nxyh)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tfft = tfft + time

# copy guard cells with standard procedure: updates fxye
   dtimer(dtime,itime,-1)
   ccguard2l(fxye,nx,ny,nxe,nye)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tguard = tguard + time

# push particles with standard procedure: updates part, wke
   wke[0] = 0.0
   dtimer(dtime,itime,-1)
   cgpush2l(part,fxye,qbme,dt,wke,idimp,np,nx,ny,nxe,nye,ipbc)
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
      print "Initial Field, Kinetic and Total Energies:"
      print "%14.7e %14.7e %14.7e" % (we, wke, wke + we)
ntime = ntime + 1

# * * * end main iteration loop * * *

print "ntime = ", ntime
print "Final Field, Kinetic and Total Energies:"
print "%14.7e %14.7e %14.7e" % (we, wke, wke + we)

print ""
print "deposit time = ", tdpost
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
