#-----------------------------------------------------------------------
# Skeleton 2D Electrostatic MPI PIC code
# written by Viktor K. Decyk, Adam Tableman, and Qiyang Hu, UCLA
import math
import numpy
from cppush2 import *
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
# idps = number of partition boundaries
idps = 2
# wke/we/wt = particle kinetic/electric field/total energy
wke = numpy.zeros((1),float_type)
we = numpy.zeros((1),float_type)
wt = numpy.zeros((1),float_type)

# declare scalars for MPI code
ntpose = 1
argv = numpy.empty((1),numpy.str)
nnvp = numpy.empty((1),int_type)
nidproc = numpy.empty((1),int_type)
nnyp = numpy.empty((1),int_type)
nnoff = numpy.empty((1),int_type)
nnpp = numpy.empty((1),int_type)
nnypmx = numpy.empty((1),int_type)
nnypmn = numpy.empty((1),int_type)
ierr = numpy.empty((1),int_type)

# declare and initialize timing data
itime = numpy.empty((4),numpy.int32)
tdpost = 0.0; tguard = 0.0; tfield = 0.0
tpush = 0.0; tsort = 0.0; tmov = 0.0
ttp = numpy.zeros((1),float_type)
tfft = numpy.zeros((2),float_type)
dtime = numpy.empty((1),double_type)

# initialize scalars for standard code
# np = total number of particles in simulation
np =  float(npx)*float(npy)
# nx/ny = number of grid points in x/y direction
nx = int(math.pow(2,indx)); ny = int(math.pow(2,indy))
nxh = int(nx/2); nyh = max(1,int(ny/2))
nxe = nx + 2; nye = ny + 2; nxeh = int(nxe/2); nnxe = ndim*nxe
nxyh = int(max(nx,ny)/2); nxhy = max(nxh,ny); ny1 = ny + 1
# nloop = number of time steps in simulation
# ntime = current time step
nloop = int(tend/dt + .0001); ntime = 0
qbme = qme
affp = float(nx)*float(ny)/np

# nvp = number of MPI ranks
# initialize for distributed memory parallel processing
cppinit2(nidproc,nnvp,0,argv)
idproc = nidproc[0]; nvp = nnvp[0]
kstrt = idproc + 1
# check if too many processors
if (nvp > ny):
   if (kstrt==1):
      print "Too many processors requested: ny, nvp=", ny, nvp
   ppexit()
   exit(0)

# initialize data for MPI code
edges = numpy.empty((idps),float_type,'F')
# calculate partition variables: edges, nyp, noff, nypmx
# edges[0:1] = lower:upper boundary of particle partition
# nyp = number of primary (complete) gridpoints in particle partition
# noff = lowermost global gridpoint in particle partition
# nypmx = maximum size of particle partition, including guard cells
# nypmn = minimum value of nyp
cpdicomp2l(edges,nnyp,nnoff,nnypmx,nnypmn,ny,kstrt,nvp,idps)
nyp = nnyp[0]; noff = nnoff[0]; nypmx = nnypmx[0]; nypmn = nnypmn[0]
if (nypmn < 1):
   if (kstrt==1):
      print "combination not supported nvp, ny =",nvp,ny
   cppexit()
   exit(0)
# initialize additional scalars for MPI code
# kxp = number of complex grids in each field partition in x direction
kxp = int((nxh - 1)/nvp) + 1
# kyp = number of complex grids in each field partition in y direction
kyp = int((ny - 1)/nvp) + 1
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
qe = numpy.empty((nxe,nypmx),float_type,'F')
# fxye = smoothed longitudinal electric field with guard cells
fxye = numpy.empty((ndim,nxe,nypmx),float_type,'F')
# qt = scalar charge density field array in fourier space
qt = numpy.empty((nye,kxp),complex_type,'F')
# fxyt = vector longitudinal electric field in fourier space
fxyt = numpy.empty((ndim,nye,kxp),complex_type,'F')
# ffc = form factor array for poisson solver
ffc = numpy.empty((nyh,kxp),complex_type,'F')
# mixup = bit reverse table for FFT
mixup = numpy.empty((nxhy),int_type,'F')
# sct = sine/cosine table for FFT
sct = numpy.empty((nxyh),complex_type,'F')
# ihole = location of hole left in particle arrays
ihole = numpy.empty((ntmax+1),int_type,'F')
# npic = scratch array for reordering particles
npic = numpy.empty((nypmx),int_type,'F')
wtot = numpy.empty((4),double_type)
work = numpy.empty((4),double_type)
info = numpy.empty((7),int_type)

# allocate data for MPI code
# bs/br = complex send/receive buffers for data transpose
bs = numpy.empty((ndim,kxp,kyp),complex_type,'F')
br = numpy.empty((ndim,kxp,kyp),complex_type,'F')
# sbufl/sbufr = particle buffers sent to nearby processors
sbufl = numpy.empty((idimp,nbmax),float_type,'F')
sbufr = numpy.empty((idimp,nbmax),float_type,'F')
# rbufl/rbufr = particle buffers received from nearby processors
rbufl = numpy.empty((idimp,nbmax),float_type,'F')
rbufr = numpy.empty((idimp,nbmax),float_type,'F')
# scr = guard cell buffer received from nearby processors
scr = numpy.empty((nxe),float_type,'F')

# prepare fft tables
cwpfft2rinit(mixup,sct,indx,indy,nxhy,nxyh)
# calculate form factors
isign = 0
cppois22(qt,fxyt,isign,ffc,ax,ay,affp,we,nx,ny,kstrt,nye,kxp,nyh)
# initialize electrons
nps = 1
nnpp[0] = 0
cpdistr2(part,edges,nnpp,nps,vtx,vty,vx0,vy0,npx,npy,nx,ny,idimp,npmax,
         idps,ipbc,ierr)
npp = nnpp[0]
# check for particle initialization error
if (ierr[0] != 0):
   if (kstrt==1):
      print "particle initialization error: ierr=", ierr[0]
   cppexit()
   exit(0)

# * * * start main iteration loop * * *

for ntime in xrange(0,nloop):
#  if (kstrt==1):
#     print "ntime = ", ntime

# deposit charge with standard procedure: updates qe
   dtimer(dtime,itime,-1)
   qe.fill(0.0)
   cppgpost2l(part,qe,npp,noff,qme,idimp,npmax,nxe,nypmx)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tdpost = tdpost + time

# add guard cells with standard procedure: updates qe
   dtimer(dtime,itime,-1)
   cppaguard2xl(qe,nyp,nx,nxe,nypmx)
   cppnaguard2l(qe,scr,nyp,nx,kstrt,nvp,nxe,nypmx)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tguard = tguard + time

# transform charge to fourier space with standard procedure: updates qt
# modifies qe
   dtimer(dtime,itime,-1)
   isign = -1
   cwppfft2r(qe,qt,bs,br,isign,ntpose,mixup,sct,ttp,indx,indy,kstrt,nvp,
             nxeh,nye,kxp,kyp,nypmx,nxhy,nxyh)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tfft[0] = tfft[0] + time
   tfft[1] = tfft[1] + ttp[0]

# calculate force/charge in fourier space with standard procedure:
# updates fxyt, we
   dtimer(dtime,itime,-1)
   isign = -1
   cppois22(qt,fxyt,isign,ffc,ax,ay,affp,we,nx,ny,kstrt,nye,kxp,nyh)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tfield = tfield + time

# transform force to real space with standard procedure: updates fxye
# modifies fxyt
   dtimer(dtime,itime,-1)
   isign = 1
   cwppfft2r2(fxye,fxyt,bs,br,isign,ntpose,mixup,sct,ttp,indx,indy,kstrt,
              nvp,nxeh,nye,kxp,kyp,nypmx,nxhy,nxyh)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tfft[0] = tfft[0] + time
   tfft[1] = tfft[1] + ttp[0]

# copy guard cells with standard procedure: updates fxye
   dtimer(dtime,itime,-1)
   cppncguard2l(fxye,nyp,kstrt,nvp,nnxe,nypmx)
   cppcguard2xl(fxye,nyp,nx,ndim,nxe,nypmx)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tguard = tguard + time

# push particles: updates part, wke, and ihole
   wke[0] = 0.0
   dtimer(dtime,itime,-1)
   cppgpush2l(part,fxye,edges,npp,noff,ihole,qbme,dt,wke,nx,ny,idimp,
              npmax,nxe,nypmx,idps,ntmax,ipbc)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tpush = tpush + time
# check for ihole overflow error
   if (ihole[0] < 0):
      ierr[0] = -ihole[0]
      print kstrt,"ihole overflow error: ntmax, ih=", ntmax, ierr[0]
      cppabort()
      cppexit()
      exit(0)

# move electrons into appropriate spatial regions: updates part, npp
   dtimer(dtime,itime,-1)
   cppmove2(part,edges,nnpp,sbufr,sbufl,rbufr,rbufl,ihole,ny,kstrt,nvp,
            idimp,npmax,idps,nbmax,ntmax,info)
   npp = nnpp[0]
   dtimer(dtime,itime,1)
   time = float(dtime)
   tmov = tmov + time
# check for particle manager error
   if (info[0] != 0):
      ierr[0] = info[0]
      if (kstrt==1):
         print "push particle manager error: ierr=", ierr[0]
      cppexit()
      exit(0)

# sort particles for standard code: updates part
   if (sortime > 0):
      if (ntime%sortime==0):
         dtimer(dtime,itime,-1)
         cppdsortp2yl(part,part2,npic,npp,noff,nyp,idimp,npmax,nypmx)
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
   cppdsum(wtot,work,4)
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
   print "MPI nodes nvp = ", nvp
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

cppexit()
