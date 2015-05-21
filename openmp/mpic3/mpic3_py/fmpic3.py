#-----------------------------------------------------------------------
# Skeleton 3D Electrostatic OpenMP PIC code
# written by Viktor K. Decyk, Adam Tableman, and Qiyang Hu, UCLA
import math
import numpy
from fmpush3 import *
from dtimer import *
from fomplib import *

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
idimp = 6; ipbc = 1
# wke/we/wt = particle kinetic/electric field/total energy
wke = numpy.zeros((1),float_type)
we = numpy.zeros((1),float_type)
wt = numpy.zeros((1),float_type)
# mx/my/mz = number of grids in x/y/z in sorting tiles
mx = 8; my = 8; mz = 8
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
# nx/ny/nz = number of grid points in x/y/z direction
np = npx*npy*npz; nx = int(math.pow(2,indx));
ny = int(math.pow(2,indy)); nz = int(math.pow(2,indz))
nxh = int(nx/2); nyh = max(1,int(ny/2)); nzh = max(1,int(nz/2))
nxe = nx + 2; nye = ny + 1; nze = nz + 1; nxeh = int(nxe/2)
nxyzh = int(max(nx,ny,nz)/2); nxhyz = max(nxh,ny,nz)
# mx1/my1/mz1 = number of tiles in x/y/z direction
mx1 = int((nx - 1)/mx + 1); my1 = int((ny - 1)/my + 1);
mz1 = int((nz - 1)/mz + 1); mxyz1 = mx1*my1*mz1
# nloop = number of time steps in simulation
# ntime = current time step
nloop = int(tend/dt + .0001); ntime = 0
qbme = qme
affp = float(nx)*float(ny)*float(nz)/float(np)

# allocate data for standard code
# part, part2 = particle arrays
part = numpy.empty((idimp,np),float_type,'F')
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
# kpic = number of particles in each tile
kpic = numpy.empty((mxyz1),int_type,'F')

# prepare fft tables
wfft3rinit(mixup,sct,indx,indy,indz,nxhyz,nxyzh)
# calculate form factors
isign = 0
mpois33(qe,fxyze,isign,ffc,ax,ay,az,affp,we,nx,ny,nz,nxeh,nye,nze,nxh,
        nyh,nzh)
# initialize electrons
distr3(part,vtx,vty,vtz,vx0,vy0,vz0,npx,npy,npz,idimp,np,nx,ny,nz,ipbc)

# find number of particles in each of mx, my, mz, tiles:
# updates kpic, nppmx
dblkp3l(part,kpic,nppmx,idimp,np,mx,my,mz,mx1,my1,mxyz1,irc)
if (irc[0] != 0):
   print "dblkp3l error, irc=", irc[0]
   exit(0)
# allocate vector particle data
nppmx0 = int((1.0 + xtras)*nppmx)
ntmax = int(xtras*nppmx)
npbmx = int(xtras*nppmx)
# ppart = tiled particle array
ppart = numpy.empty((idimp,nppmx0,mxyz1),float_type,'F')
# ppbuff = buffer array for reordering tiled particle array
ppbuff = numpy.empty((idimp,npbmx,mxyz1),float_type,'F')
# ncl = number of particles departing tile in each direction
ncl = numpy.empty((26,mxyz1),int_type,'F')
# ihole = location/destination of each particle departing tile
ihole = numpy.empty((2,ntmax+1,mxyz1),int_type,'F')
# copy ordered particle data for OpenMP: updates ppart and kpic
ppmovin3l(part,ppart,kpic,nppmx0,idimp,np,mx,my,mz,mx1,my1,mxyz1,irc)
if (irc[0] != 0):
   print "ppmovin3l overflow error, irc=", irc[0]
   exit(0)
# sanity check
ppcheck3l(ppart,kpic,idimp,nppmx0,nx,ny,nz,mx,my,mz,mx1,my1,mz1,irc)
if (irc[0] != 0):
   print "ppcheck3l error, irc=", irc[0]
   exit(0)

# * * * start main iteration loop * * *

for ntime in xrange(0,nloop):
#  print "ntime = ", ntime

# deposit charge with OpenMP: updates qe
   dtimer(dtime,itime,-1)
   qe.fill(0.0)
   gppost3l(ppart,qe,kpic,qme,nppmx0,idimp,mx,my,mz,nxe,nye,nze,mx1,my1,
            mxyz1)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tdpost = tdpost + time

# add guard cells with OpenMP: updates qe
   dtimer(dtime,itime,-1)
   aguard3l(qe,nx,ny,nz,nxe,nye,nze)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tguard = tguard + time

# transform charge to fourier space with OpenMP: updates qe
   dtimer(dtime,itime,-1)
   isign = -1
   wfft3rmx(qe,isign,mixup,sct,indx,indy,indz,nxeh,nye,nze,nxhyz,nxyzh)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tfft = tfft + time

# calculate force/charge in fourier space with OpenMP: updates fxyze, we
   dtimer(dtime,itime,-1)
   isign = -1
   mpois33(qe,fxyze,isign,ffc,ax,ay,az,affp,we,nx,ny,nz,nxeh,nye,nze,nxh,
           nyh,nzh)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tfield = tfield + time

# transform force to real space with OpenMP: updates fxyze
   dtimer(dtime,itime,-1)
   isign = 1
   wfft3rm3(fxyze,isign,mixup,sct,indx,indy,indz,nxeh,nye,nze,nxhyz,
            nxyzh)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tfft = tfft + time

# copy guard cells with OpenMP: updates fxyze
   dtimer(dtime,itime,-1)
   cguard3l(fxyze,nx,ny,nz,nxe,nye,nze)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tguard = tguard + time

# push particles with OpenMP:
   wke[0] = 0.0
   dtimer(dtime,itime,-1)
# updates ppart, wke
#  gppush3l(ppart,fxyze,kpic,qbme,dt,wke,idimp,nppmx0,nx,ny,nz,mx,my,
#           mz,nxe,nye,nze,mx1,my1,mxyz1,ipbc)
# updates ppart, ncl, ihole, wke, irc
   gppushf3l(ppart,fxyze,kpic,ncl,ihole,qbme,dt,wke,idimp,nppmx0,nx,ny,
             nz,mx,my,mz,nxe,nye,nze,mx1,my1,mxyz1,ntmax,irc)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tpush = tpush + time
   if (irc[0] != 0):
      print "gppushf3l error, irc=", irc[0]
      exit(0)

# reorder particles by cell with OpenMP:
   dtimer(dtime,itime,-1)
# updates ppart, ppbuff, kpic, ncl, ihole, and irc
#  pporder3l(ppart,ppbuff,kpic,ncl,ihole,idimp,nppmx0,nx,ny,nz,mx,my,mz,
#            mx1,my1,mz1,npbmx,ntmax,irc)
# updates ppart, ppbuff, kpic, ncl, and irc
   pporderf3l(ppart,ppbuff,kpic,ncl,ihole,idimp,nppmx0,mx1,my1,mz1,
              npbmx,ntmax,irc)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tsort = tsort + time
   if (irc[0] != 0):
      print "pporderf3l error, ntmax, irc=", ntmax, irc[0]
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
print ""

