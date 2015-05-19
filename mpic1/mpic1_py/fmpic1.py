#-----------------------------------------------------------------------
# Skeleton 1D Electrostatic OpenMP PIC code
# written by Viktor K. Decyk, Adam Tableman, and Qiyang Hu, UCLA
import math
import numpy
from fmpush1 import *
from dtimer import *
from fomplib import *

int_type = numpy.int32
double_type = numpy.float64
float_type = numpy.float32
complex_type = numpy.complex64

# indx = exponent which determines grid points in x direction:
# nx = 2**indx.
indx =   9
# npx = number of electrons distributed in x direction.
npx = 18432
# tend = time at end of simulation, in units of plasma frequency.
# dt = time interval between successive calculations.
# qme = charge on electron, in units of e.
tend = 10.0; dt = 0.1; qme = -1.0
# vtx = thermal velocity of electrons in x direction
# vx0 = drift velocity of electrons in x direction.
vtx = 1.0; vx0 = 0.0
# ax = smoothed particle size in x direction
ax = .912871
# idimp = number of particle coordinates = 2
# ipbc = particle boundary condition: 1 = periodic
idimp = 2; ipbc = 1
# wke/we/wt = particle kinetic/electric field/total energy
wke = numpy.zeros((1),float_type)
we = numpy.zeros((1),float_type)
wt = numpy.zeros((1),float_type)
# mx = number of grids in x in sorting tiles
mx = 32
# xtras = fraction of extra particles needed for particle management
xtras = 0.2

# declare scalars for OpenMP code
nppmx = numpy.empty((1),int_type)
irc = numpy.zeros((1),int_type)

# declare and initialize timing data
itime = numpy.empty((4),numpy.int32)
tdpost = 0.0; tguard = 0.0; tfft = 0.0; tfield = 0.0
tpush = 0.0; tsort = 0.0
dtime = numpy.empty((1),double_type)

# nvp = number of shared memory nodes (0=default)
nvp = 0
#nvp = int(input("enter number of nodes: "))
# initialize for shared memory parallel processing
init_omp(nvp)

# initialize scalars for standard code
# np = total number of particles in simulation
# nx = number of grid points in x direction
np = npx; nx = int(math.pow(2,indx)); nxh = int(nx/2)
nxe = nx + 2
# mx1 = number of tiles in x direction
mx1 = int((nx - 1)/mx + 1)
# nloop = number of time steps in simulation
# ntime = current time step
nloop = int(tend/dt + .0001); ntime = 0
qbme = qme
affp = float(nx)/float(np)

# allocate data for standard code
# part = particle array
part = numpy.empty((idimp,np),float_type,'F')
# qe = electron charge density with guard cells
qe = numpy.empty((nxe),float_type,'F')
# fxe = smoothed electric field with guard cells
fxe = numpy.empty((nxe),float_type,'F')
# ffc = form factor array for poisson solver
ffc = numpy.empty((nxh),complex_type,'F')
# mixup = bit reverse table for FFT
mixup = numpy.empty((nxh),int_type,'F')
# sct = sine/cosine table for FFT
sct = numpy.empty((nxh),complex_type,'F')
# kpic = number of particles in each tile
kpic = numpy.empty((mx1),int_type,'F')

# prepare fft tables
wfft1rinit(mixup,sct,indx,nxh)
# calculate form factors
isign = 0
pois1(qe,fxe,isign,ffc,ax,affp,we,nx)
# initialize electrons
distr1(part,vtx,vx0,npx,idimp,np,nx,ipbc)

# find number of particles in each of mx, tiles: updates kpic, nppmx
dblkp1l(part,kpic,nppmx,idimp,np,mx,mx1,irc)
if (irc[0] != 0):
   print "dblkp1l error, irc=", irc[0]
   exit(0)
# allocate vector particle data
nppmx0 = int((1.0 + xtras)*nppmx)
ntmax = int(xtras*nppmx)
npbmx = int(xtras*nppmx)
# ppart = tiled particle array
ppart = numpy.empty((idimp,nppmx0,mx1),float_type,'F')
# ppbuff = buffer array for reordering tiled particle array
ppbuff = numpy.empty((idimp,npbmx,mx1),float_type,'F')
# ncl = number of particles departing tile in each direction
ncl = numpy.empty((2,mx1),int_type,'F')
# ihole = location/destination of each particle departing tile
ihole = numpy.empty((2,ntmax+1,mx1),int_type,'F')
# copy ordered particle data for OpenMP: updates ppart and kpic
ppmovin1l(part,ppart,kpic,nppmx0,idimp,np,mx,mx1,irc)
if (irc[0] != 0):
   print "ppmovin1l overflow error, irc=", irc[0]
   exit(0)
# sanity check
ppcheck1l(ppart,kpic,idimp,nppmx0,nx,mx,mx1,irc)
if (irc[0] != 0):
   print "ppcheck1l error, irc=", irc[0]
   exit(0)

# * * * start main iteration loop * * *

for ntime in xrange(0,nloop):
#  print "ntime = ", ntime

# deposit charge with OpenMP: updates qe
   dtimer(dtime,itime,-1)
   qe.fill(0.0)
   gppost1l(ppart,qe,kpic,qme,nppmx0,idimp,mx,nxe,mx1)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tdpost = tdpost + time

# add guard cells with standard procedure: updates qe
   dtimer(dtime,itime,-1)
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

# calculate force/charge in fourier space with standard procedure:
# updates fxe, we
   dtimer(dtime,itime,-1)
   isign = -1
   pois1(qe,fxe,isign,ffc,ax,affp,we,nx)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tfield = tfield + time

# transform force to real space with standard procedure: updates fxe, qe
   dtimer(dtime,itime,-1)
   isign = 1
   fft1rxx(fxe,qe,isign,mixup,sct,indx,nxe,nxh)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tfft = tfft + time

# copy guard cells with standard procedure: updates fxe
   dtimer(dtime,itime,-1)
   cguard1l(fxe,nx,nxe)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tguard = tguard + time

# push particles with OpenMP:
   wke[0] = 0.0
   dtimer(dtime,itime,-1)
# updates part, wke
#  gppush1l(ppart,fxe,kpic,qbme,dt,wke,idimp,nppmx0,nx,mx,nxe,mx1,ipbc)
# updates ppart, ncl, ihole, wke, irc
   gppushf1l(ppart,fxe,kpic,ncl,ihole,qbme,dt,wke,idimp,nppmx0,nx,mx,
             nxe,mx1,ntmax,irc)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tpush = tpush + time
   if (irc[0] != 0):
      print "gppushf1l error, irc=", irc[0]
      exit(0)

# reorder particles by tile with OpenMP:
   dtimer(dtime,itime,-1)
# updates ppart, ppbuff, kpic, ncl, ihole, and irc
#  pporder1l(ppart,ppbuff,kpic,ncl,ihole,idimp,nppmx0,nx,mx,mx1,npbmx,
#            ntmax,irc)
# updates ppart, ppbuff, kpic, ncl, and irc
   pporderf1l(ppart,ppbuff,kpic,ncl,ihole,idimp,nppmx0,mx1,npbmx,ntmax,
              irc)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tsort = tsort + time
   if (irc[0] != 0):
      print "pporderf1l error, ntmax, irc=", ntmax, irc[0]
      exit(0)

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
