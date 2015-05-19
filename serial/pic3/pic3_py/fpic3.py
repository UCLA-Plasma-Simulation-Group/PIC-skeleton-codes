#-----------------------------------------------------------------------
# Skeleton 3D Electrostatic PIC code
# written by Viktor K. Decyk, Adam Tableman, and Qiyang Hu, UCLA
import math
import numpy
from fpush3 import *
from dtimer import *

int_type = numpy.int32
double_type = numpy.float64
float_type = numpy.float32
complex_type = numpy.complex64

# indx/indy/indz = exponent which determines grid points in x/y/z
# direction: nx = 2**indx, ny = 2**indy, nz = 2**indz.
indx =   7; indy =   7; indz =   7
# npx/npy/npz = number of electrons distributed in x/y/z direction.
npx =  384; npy =   384; npz =   384
# ndim = number of velocity coordinates = 3
ndim = 3
# tend = time at end of simulation, in units of plasma frequency.
# dt = time interval between successive calculations.
# qme = charge on electron, in units of e.
tend = 10.0; dt = 0.1; qme = -1.0
# vtx/vty/vtz = thermal velocity of electrons in x/y/z direction
vtx = 1.0; vty = 1.0; vtz = 1.0
# vx0/vy0/vz0 = drift velocity of electrons in x/y/z direction
vx0 = 0.0; vy0 = 0.0; vz0 = 0.0
# ax/ay/az = smoothed particle size in x/y/z direction
ax = .912871; ay = .912871; az = .912871
# idimp = number of particle coordinates = 6
# ipbc = particle boundary condition: 1 = periodic
# sortime = number of time steps between standard electron sorting
idimp = 6; ipbc = 1; sortime = 20
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
# nx/ny/nz = number of grid points in x/y/z direction
np = npx*npy*npz; nx = int(math.pow(2,indx));
ny = int(math.pow(2,indy)); nz = int(math.pow(2,indz))
nxh = int(nx/2); nyh = max(1,int(ny/2)); nzh = max(1,int(nz/2))
nxe = nx + 2; nye = ny + 1; nze = nz + 1; nxeh = int(nxe/2)
nxyzh = int(max(nx,ny,nz)/2); nxhyz = max(nxh,ny,nz)
ny1 = ny + 1; nyz1 = ny1*(nz + 1)
# nloop = number of time steps in simulation
# ntime = current time step
nloop = int(tend/dt + .0001); ntime = 0
qbme = qme
affp = float(nx)*float(ny)*float(nz)/float(np)

# allocate data for standard code
# part, part2 = particle arrays
part = numpy.empty((idimp,np),float_type,'F')
if (sortime > 0):
   part2 = numpy.empty((idimp,np),float_type,'F')
# qe = electron charge density with guard cells
qe = numpy.empty((nxe,nye,nze),float_type,'F')
# fxyze = smoothed electric field with guard cells
fxyze = numpy.empty((ndim,nxe,nye,nze),float_type,'F')
# ffc = form factor array for poisson solver
ffc = numpy.empty((nxh,nyh,nzh),complex_type,'F')
# mixup = bit reverse table for FFT
mixup = numpy.empty((nxhyz),int_type,'F')
# sct = sine/cosine table for FFT
sct = numpy.empty((nxyzh),complex_type,'F')
# npic = scratch array for reordering particles
npic = numpy.empty((nyz1),int_type,'F')

# prepare fft tables
wfft3rinit(mixup,sct,indx,indy,indz,nxhyz,nxyzh)
# calculate form factors
isign = 0
pois33(qe,fxyze,isign,ffc,ax,ay,az,affp,we,nx,ny,nz,nxeh,nye,nze,nxh,
       nyh,nzh)
# initialize electrons
distr3(part,vtx,vty,vtz,vx0,vy0,vz0,npx,npy,npz,idimp,np,nx,ny,nz,ipbc)

# * * * start main iteration loop * * *

for ntime in xrange(0,nloop):
#  print "ntime = ", ntime

# deposit charge with standard procedure: updates qe
   dtimer(dtime,itime,-1)
   qe.fill(0.0)
   gpost3l(part,qe,qme,np,idimp,nxe,nye,nze)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tdpost = tdpost + time

# add guard cells with standard procedure: updates qe
   dtimer(dtime,itime,-1)
   aguard3l(qe,nx,ny,nz,nxe,nye,nze)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tguard = tguard + time

# transform charge to fourier space with standard procedure: updates qe
   dtimer(dtime,itime,-1)
   isign = -1
   wfft3rx(qe,isign,mixup,sct,indx,indy,indz,nxeh,nye,nze,nxhyz,nxyzh)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tfft = tfft + time

# calculate force/charge in fourier space with standard procedure:
# updates fxyze
   dtimer(dtime,itime,-1)
   isign = -1
   pois33(qe,fxyze,isign,ffc,ax,ay,az,affp,we,nx,ny,nz,nxeh,nye,nze,nxh,
          nyh,nzh)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tfield = tfield + time

# transform force to real space with standard procedure: updates fxyze
   dtimer(dtime,itime,-1)
   isign = 1
   wfft3r3(fxyze,isign,mixup,sct,indx,indy,indz,nxeh,nye,nze,nxhyz,
           nxyzh)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tfft = tfft + time

# copy guard cells with standard procedure: updates fxyze
   dtimer(dtime,itime,-1)
   cguard3l(fxyze,nx,ny,nz,nxe,nye,nze)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tguard = tguard + time

# push particles with standard procedure: updates part, wke
   wke[0] = 0.0
   dtimer(dtime,itime,-1)
   gpush3l(part,fxyze,qbme,dt,wke,idimp,np,nx,ny,nz,nxe,nye,nze,ipbc)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tpush = tpush + time

# sort particles by cell for standard procedure
   if (sortime > 0):
      if (ntime%sortime==0):
         dtimer(dtime,itime,-1)
         dsortp3yzl(part,part2,npic,idimp,np,ny1,nyz1)
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

