#-----------------------------------------------------------------------
# Skeleton 3D Electromagnetic MPI PIC code
# written by Viktor K. Decyk, Adam Tableman, and Qiyang Hu, UCLA
import math
import numpy
from fppush3 import *
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
# idps = number of partition boundaries = 4
# idds = dimensionality of domain decomposition = 2
idps = 4; idds =    2
# wke/we/wt = particle kinetic/electric field/total energy
wke = numpy.zeros((1),float_type)
we = numpy.zeros((1),float_type)
wt = numpy.zeros((1),float_type)
# declare scalars for standard code
ierr = numpy.empty((1),int_type)

# declare scalars for MPI code
ntpose = 1
nnvp = numpy.empty((1),int_type)
nidproc = numpy.empty((1),int_type)
nnvpy = numpy.empty((1),int_type)
nnvpz = numpy.empty((1),int_type)
nnpp = numpy.empty((1),int_type)
nnypmx = numpy.empty((1),int_type)
nnzpmx = numpy.empty((1),int_type)
nnypmn = numpy.empty((1),int_type)
nnzpmn = numpy.empty((1),int_type)

# declare and initialize timing data
itime = numpy.empty((4),numpy.int32)
tdpost = 0.0; tguard = 0.0; tfield = 0.0
tpush = 0.0; tsort = 0.0; tmov = 0.0
ttp = numpy.zeros((1),float_type)
tfft = numpy.zeros((2),float_type)
dtime = numpy.empty((1),double_type)

# initialize scalars for standard code
# np = total number of particles in simulation
np =  float(npx)*float(npy)*float(npz)
# nx/ny/nz = number of grid points in x/y/z direction
nx = int(math.pow(2,indx)); ny = int(math.pow(2,indy))
nz = int(math.pow(2,indz))
nxh = int(nx/2); nyh = max(1,int(ny/2)); nzh = max(1,int(nz/2))
nxe = nx + 2; nye = ny + 2; nze = nz + 2
nxeh = int(nxe/2); nnxe = ndim*nxe
nxyzh = int(max(nx,ny,nz)/2); nxhyz = max(nxh,ny,nz)
# nloop = number of time steps in simulation
# ntime = current time step
nloop = int(tend/dt + .0001); ntime = 0
qbme = qme
affp = float(nx)*float(ny)*float(nz)/np

# nvp = number of MPI ranks
# initialize for distributed memory parallel processing
ppinit2(nidproc,nnvp)
idproc = nidproc[0]; nvp = nnvp[0]
kstrt = idproc + 1
# obtain 2D partition (nvpy,nvpz) from nvp:
# nvpy/nvpz = number of processors in y/z
fcomp32(nvp,nx,ny,nz,nnvpy,nnvpz,ierr)
nvpy = nnvpy[0]; nvpz = nnvpz[0]
if (ierr[0] != 0):
   if (kstrt==1):
      print "fcomp32 error: nvp,nvpy,nvpz=", nvp, nvpy, nvpz
   ppexit()
   exit(0)

# initialize data for MPI code
edges = numpy.empty((idps),float_type,'F')
nyzp = numpy.empty((idds),int_type)
noff = numpy.empty((idds),int_type)
# calculate partition variables:
# edges, nyzp, noff, nypmx, nzpmx, nypmn, nzpmn
# edges[0:1] = lower:upper boundary of particle partition in y
# edges[2:3] = back:front boundary of particle partition in z
# nyzp[0:1] = number of primary (complete) gridpoints in y/z
# noff[0:1] = lowermost global gridpoint in y/z in particle partition
# nypmx = maximum size of particle partition in y, including guard cells
# nzpmx = maximum size of particle partition in z, including guard cells
# nypmn = minimum value of nyzp[0]
# nzpmn = minimum value of nyzp[1]
pdicomp32l(edges,nyzp,noff,nnypmx,nnzpmx,nnypmn,nnzpmn,ny,nz,kstrt,nvpy,
           nvpz,idps,idds)
nypmx = nnypmx[0]; nzpmx = nnzpmx[0]
nypmn = nnypmn[0]; nzpmn = nnzpmn[0]
if (kstrt==1):
   if (nypmn < 1):
      print "combination not supported nvpy, ny =",nvpy,ny
   if (nzpmn < 1):
      print "combination not supported nvpz, nz =",nvpz,nz
if ((nypmn < 1) or (nzpmn < 1)):
   ppexit()
   exit(0)
# initialize additional scalars for MPI code
# kyp = number of complex grids in each field partition in y direction
kyp = int((ny - 1)/nvpy) + 1
# kzp = number of complex grids in each field partition in z direction
kzp = int((nz - 1)/nvpz) + 1
# kxyp = number of complex grids in each field partition in x direction
# in transposed data
kxyp = int((nxh - 1)/nvpy) + 1
# kyzp = number of complex grids in each field partition in y direction,
# in transposed data
kyzp = int((ny - 1)/nvpz) + 1; kzyp = max(kyzp,kyp)
# dimension for scratch array for reordering particles
nyzpm1 = (kyp + 1)*(kzp + 1)
# npmax = maximum number of electrons in each partition
npmax = int((np/float(nvp))*1.25)
# nbmax = size of buffer for passing particles between processors
nbmax = int(0.1*float(npmax))
# ntmax = size of ihole buffer for particles leaving processor
ntmax = 2*nbmax

# allocate data for standard code
# part, part2 = particle arrays
part = numpy.empty((idimp,npmax),float_type,'F')
if (sortime > 0):
   part2 = numpy.empty((idimp,npmax),float_type,'F')
# qe = electron charge density with guard cells
qe = numpy.empty((nxe,nypmx,nzpmx),float_type,'F')
# fxyze = smoothed longitudinal electric field with guard cells
fxyze = numpy.empty((ndim,nxe,nypmx,nzpmx),float_type,'F')
# qt, qs = scalar charge density field arrays in fourier space
qt = numpy.empty((nze,kxyp,kyzp),complex_type,'F')
qs = numpy.empty((nye,kxyp,nzpmx),complex_type,'F')
# fxyzt, fxyzs = vector longitudinal electric field in fourier space
fxyzt  = numpy.empty((ndim,nze,kxyp,kyzp),complex_type,'F')
fxyzs = numpy.empty((ndim,nye,kxyp,nzpmx),complex_type,'F')
# ffc = form factor array for poisson solvers
ffc = numpy.empty((nzh,kxyp,kyzp),complex_type,'F')
# mixup = bit reverse table for FFT
mixup = numpy.empty((nxhyz),int_type,'F')
# sct = sine/cosine table for FFT
sct = numpy.empty((nxyzh),complex_type,'F')
# ihole = location of hole left in particle arrays
ihole = numpy.empty((ntmax+1,2),int_type,'F')
# npic = scratch array for reordering particles
npic = numpy.empty((nyzpm1),int_type,'F')
wtot = numpy.empty((4),double_type)
work = numpy.empty((4),double_type)
info = numpy.empty((7),int_type)

# allocate data for MPI code
# bs/br = complex send/receive buffers for data transpose
bs = numpy.empty((ndim,kxyp*kzyp,kzp),complex_type,'F')
br = numpy.empty((ndim,kxyp*kzyp,kzp),complex_type,'F')
# sbufl/sbufr = particle buffers sent to nearby processors
sbufl = numpy.empty((idimp,nbmax),float_type,'F')
sbufr = numpy.empty((idimp,nbmax),float_type,'F')
# rbufl/rbufr = particle buffers received from nearby processors
rbufl = numpy.empty((idimp,nbmax),float_type,'F')
rbufr = numpy.empty((idimp,nbmax),float_type,'F')
# scr/scs = guard cell buffers received/sent from nearby processors
scr = numpy.empty((nxe,nypmx),float_type,'F')
scs = numpy.empty((nnxe,2*nzpmx),float_type,'F')

# prepare fft tables
wpfft32rinit(mixup,sct,indx,indy,indz,nxhyz,nxyzh)
# calculate form factors
isign = 0
ppois332(qt,fxyzt,isign,ffc,ax,ay,az,affp,we,nx,ny,nz,kstrt,nvpy,nvpz,
         nze,kxyp,kyzp,nzh)

# initialize electrons
nps = 1
nnpp[0] = 0
pdistr32(part,edges,nnpp,nps,vtx,vty,vtz,vx0,vy0,vz0,npx,npy,npz,nx,ny,
         nz,idimp,npmax,idps,ipbc,ierr)
npp = nnpp[0]
# check for particle initialization error
if (ierr[0] != 0):
   if (kstrt==1):
      print "particle initialization error: ierr=", ierr[0]
   ppexit()
   exit(0)

# * * * start main iteration loop * * *

for ntime in xrange(0,nloop):
#  if (kstrt==1):
#     print "ntime = ", ntime

# deposit charge with standard procedure: updates qe
   dtimer(dtime,itime,-1)
   qe.fill(0.0)
   ppgpost32l(part,qe,npp,noff,qme,idimp,npmax,nxe,nypmx,nzpmx,idds)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tdpost = tdpost + time

# add guard cells with standard procedure: updates qe
   dtimer(dtime,itime,-1)
   ppaguard32xl(qe,nyzp,nx,nxe,nypmx,nzpmx,idds)
   ppnaguard32l(qe,scs,scr,nyzp,kstrt,nvpy,nvpz,nx,nxe,nypmx,nzpmx,idds)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tguard = tguard + time

# transform charge to fourier space with standard procedure: updates qt
# modifies qe
   dtimer(dtime,itime,-1)
   isign = -1
   wppfft32r(qe,qs,qt,bs,br,isign,ntpose,mixup,sct,ttp,indx,indy,indz,
             kstrt,nvpy,nvpz,nxeh,nye,nze,kxyp,kyp,kyzp,kzp,kxyp,nypmx,
             kyzp,nzpmx,kzyp,nxhyz,nxyzh)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tfft[0] = tfft[0] + time
   tfft[1] = tfft[1] + ttp[0]

# calculate force/charge in fourier space with standard procedure:
# updates fxyzt, we
   dtimer(dtime,itime,-1)
   isign = -1
   ppois332(qt,fxyzt,isign,ffc,ax,ay,az,affp,we,nx,ny,nz,kstrt,nvpy,
            nvpz,nze,kxyp,kyzp,nzh)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tfield = tfield + time

# transform force to real space with standard procedure: updates fxyze
# modifies fxyzt
   dtimer(dtime,itime,-1)
   isign = 1
   wppfft32r3(fxyze,fxyzs,fxyzt,bs,br,isign,ntpose,mixup,sct,ttp,indx,
              indy,indz,kstrt,nvpy,nvpz,nxeh,nye,nze,kxyp,kyp,kyzp,kzp,
              kxyp,nypmx,kyzp,nzpmx,kzyp,nxhyz,nxyzh)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tfft[0] = tfft[0] + time
   tfft[1] = tfft[1] + ttp[0]

# copy guard cells with standard procedure: updates fxyze
   dtimer(dtime,itime,-1)
   ppncguard32l(fxyze,scs,nyzp,kstrt,nvpy,nvpz,nnxe,nypmx,nzpmx,idds)
   ppcguard32xl(fxyze,nyzp,nx,ndim,nxe,nypmx,nzpmx,idds)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tguard = tguard + time

# push particles with standard procedure: updates part, wke, and ihole
   wke[0] = 0.0
   dtimer(dtime,itime,-1)
   ppgpush32l(part,fxyze,edges,npp,noff,ihole,qbme,dt,wke,nx,ny,nz,
              idimp,npmax,nxe,nypmx,nzpmx,idps,idds,ntmax,ipbc)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tpush = tpush + time
# check for ihole overflow error
   if (ihole[0,0] < 0):
      ierr[0] = -ihole[0,0]
      print kstrt,"ihole overflow error: ntmax, ih=", ntmax, ierr[0]
      ppabort()
      ppexit()
      exit(0)

# move electrons into appropriate spatial regions: updates part, npp
   dtimer(dtime,itime,-1)
   ppmove32(part,edges,nnpp,sbufr,sbufl,rbufr,rbufl,ihole,ny,nz,kstrt,
            nvpy,nvpz,idimp,npmax,idps,nbmax,ntmax,info)
   npp = nnpp[0]
   dtimer(dtime,itime,1)
   time = float(dtime)
   tmov = tmov + time
# check for particle manager error
   if (info[0] != 0):
      ierr[0] = info[0]
      if (kstrt==1):
         print "push particle manager error: ierr=", ierr[0]
      ppexit()
      exit(0)

# sort particles by cell for standard procedure
   if (sortime > 0):
      if (ntime%sortime==0):
         dtimer(dtime,itime,-1)
         ppdsortp32yzl(part,part2,npic,npp,noff,nyzp,idimp,npmax,nyzpm1,
                       idds)
# exchange pointers
         tpart = part
         part = part2
         part2 = tpart
         dtimer(dtime,itime,1)
         time = float(dtime)
         tsort = tsort + time
      pass
   pass

# energy diagnostic
   wtot[0] = we
   wtot[1] = wke[0]
   wtot[2] = 0.0
   wtot[3] = we + wke[0]
   ppdsum(wtot,work,4)
   we[0] = wtot[0]
   wke[0] = wtot[1]
   if (ntime==0):
      if (kstrt==1):
         print "Initial Field, Kinetic and Total Energies:"
         print "%14.7e %14.7e %14.7e" % (we, wke, wke + we)
ntime = ntime + 1

# * * * end main iteration loop * * *

if (kstrt==1):
   print "ntime = ", ntime
   print "MPI nodes nvpy, nvpz = ", nvpy, nvpz
   print "Final Field, Kinetic and Total Energies:"
   print "%14.7e %14.7e %14.7e" % (we, wke, wke + we)

   print ""
   print "deposit time = ", tdpost
   print "guard time = ", tguard
   print "solver time = ", tfield
   print "fft and transpose time = ", tfft[0], tfft[1]
   print "push time = ", tpush
   print "particle move time = ", tmov
   print "sort time = ", tsort
   tfield = tfield + tguard + tfft[0]
   print "total solver time = ", tfield
   tsort = tsort + tmov
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

ppexit()
