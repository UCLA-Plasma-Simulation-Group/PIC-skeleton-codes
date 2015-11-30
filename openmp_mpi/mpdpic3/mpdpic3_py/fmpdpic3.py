#-----------------------------------------------------------------------
# Skeleton 3D Darwin MPI/OpenMP PIC code
# written by Viktor K. Decyk, Adam Tableman, and Qiyang Hu, UCLA
import math
import numpy
from fmpdpush3 import *
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
# ci = reciprocal of velocity of light.
ax = .912871; ay = .912871; az = .912871; ci = 0.1
# idimp = number of particle coordinates = 6
# ipbc = particle boundary condition: 1 = periodic
idimp = 6; ipbc = 1
# omx/omy/omz = magnetic field electron cyclotron frequency in x/y/z 
omx = 0.4; omy = 0.0; omz = 0.0
# ndc = number of corrections in darwin iteration
ndc = 1
# idps = number of partition boundaries = 4
# idds = dimensionality of domain decomposition = 2
idps = 4; idds =    2
# wke/we = particle kinetic/electrostatic field energy
# wf/wm/wt = magnetic field/transverse electric field/total energy
wke = numpy.zeros((1),float_type)
we = numpy.zeros((1),float_type)
wf = numpy.zeros((1),float_type)
wm = numpy.zeros((1),float_type)
wt = numpy.zeros((1),float_type)
zero = 0.0
# mx/my/mz = number of grids in x/y/z in sorting tiles
# sorting tiles, should be less than or equal to 16
mx = 8; my = 8; mz = 8
# fraction of extra particles needed for particle management
xtras = 0.2
# declare scalars for standard code
wpmax = numpy.empty((1),float_type)
wpmin = numpy.empty((1),float_type)

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
ierr = numpy.empty((1),int_type)

# declare scalars for OpenMP code
nppmx = numpy.empty((1),int_type)
irc = numpy.zeros((1),int_type)
nvpp = numpy.empty((1),int_type)

# declare and initialize timing data
itime = numpy.empty((4),numpy.int32)
tdpost = 0.0; tguard = 0.0; tfield = 0.0; tdjpost = 0.0
tdcjpost = 0.0; tpush = 0.0; tsort = 0.0; tmov = 0.0
ttp = numpy.zeros((1),float_type)
tfft = numpy.zeros((2),float_type)
dtime = numpy.empty((1),double_type)

# nvpp = number of shared memory nodes (0=default)
nvpp = 0
#nvpp = int(input("enter number of nodes: "))
# initialize for shared memory parallel processing
init_omp(nvpp)

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
# mx1 = number of tiles in x direction
mx1 = int((nx - 1)/mx) + 1
# nloop = number of time steps in simulation
# ntime = current time step
nloop = int(tend/dt + .0001); ntime = 0
# mdim = dimension of amu array
mdim = 2*ndim
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
# npmax = maximum number of electrons in each partition
npmax = int((np/float(nvp))*1.25)
# myp1/mzp1 = number of tiles in y/z direction
myp1 = int((nyzp[0] - 1)/my) + 1; mzp1 = int((nyzp[1] - 1)/mz) + 1
# mxzyp1 = mx1*max(max(mzp1),max(myp1))
mxzyp1 = mx1*max((nzpmx-2)/mz+1,(nypmx-2)/my+1)
mxyzp1 = mx1*myp1*mzp1

# allocate data for standard code
# part = particle array
part = numpy.empty((idimp,npmax),float_type,'F')
# qe = electron charge density with guard cells
qe = numpy.empty((nxe,nypmx,nzpmx),float_type,'F')
# cue = electron current density with guard cells
cue = numpy.empty((ndim,nxe,nypmx,nzpmx),float_type,'F')
# dcu = acceleration density with guard cells
dcu = numpy.empty((ndim,nxe,nypmx,nzpmx),float_type,'F')
# cus = smoothed transverse electric field with guard cells
cus = numpy.empty((ndim,nxe,nypmx,nzpmx),float_type,'F')
# amu = momentum flux with guard cells
amu = numpy.empty((mdim,nxe,nypmx,nzpmx),float_type,'F')
# exyze = smoothed total electric field with guard cells
exyze = numpy.empty((ndim,nxe,nypmx,nzpmx),float_type,'F')
# fxyze = smoothed longitudinal electric field with guard cells
fxyze = numpy.empty((ndim,nxe,nypmx,nzpmx),float_type,'F')
# bxyze = smoothed magnetic field with guard cells
bxyze = numpy.empty((ndim,nxe,nypmx,nzpmx),float_type,'F')
# ss = scratch array for WPPFFT32RMN
ss = numpy.empty((mdim*nxeh,nzpmx),complex_type,'F')
# qt, qs = scalar charge density field arrays in fourier space
qt = numpy.empty((nze,kxyp,kyzp),complex_type,'F')
qs = numpy.empty((nye,kxyp,nzpmx),complex_type,'F')
# cut = vector current density in fourier space
cut = numpy.empty((ndim,nze,kxyp,kyzp),complex_type,'F')
# dcut = vector acceleration density in fourier space
dcut = numpy.empty((ndim,nze,kxyp,kyzp),complex_type,'F')
# exyzt = vector transverse electric field in fourier space
exyzt = numpy.empty((ndim,nze,kxyp,kyzp),complex_type,'F')
# amut = tensor momentum flux in fourier space
amut = numpy.empty((mdim,nze,kxyp,kyzp),complex_type,'F')
# fxyzs = vector field array in fourier space
fxyzs = numpy.empty((ndim,nye,kxyp,nzpmx),complex_type,'F')
# amus = tensor field array in fourier space
amus = numpy.empty((mdim,nye,kxyp,nzpmx),complex_type,'F')
# fxyzt = vector longitudinal electric field in fourier space
fxyzt = numpy.empty((ndim,nze,kxyp,kyzp),complex_type,'F')
# bxyzt = vector magnetic field array in fourier space
bxyzt = numpy.empty((ndim,nze,kxyp,kyzp),complex_type,'F')
# ffc, ffe = form factor arrays for poisson solvers
ffc = numpy.empty((nzh,kxyp,kyzp),complex_type,'F')
ffe = numpy.empty((nzh,kxyp,kyzp),complex_type,'F')
# mixup = bit reverse table for FFT
mixup = numpy.empty((nxhyz),int_type,'F')
# sct = sine/cosine table for FFT
sct = numpy.empty((nxyzh),complex_type,'F')
# kpic = number of particles in each tile
kpic = numpy.empty((mxyzp1),int_type,'F')
wtot = numpy.empty((7),double_type)
work = numpy.empty((7),double_type)

# allocate data for MPI code
# bs/br = complex send/receive buffers for data transpose
bs = numpy.empty((mdim,kxyp*kzyp,kzp),complex_type,'F')
br = numpy.empty((mdim,kxyp*kzyp,kzp),complex_type,'F')
# scr/scs = guard cell buffers received/sent from nearby processors
scr = numpy.empty((mdim*nxe,nypmx),float_type,'F')
scs = numpy.empty((mdim*nxe,2*nzpmx),float_type,'F')

# prepare fft tables
wpfft32rinit(mixup,sct,indx,indy,indz,nxhyz,nxyzh)
# calculate form factor: ffc
isign = 0
mppois332(qt,fxyzt,isign,ffc,ax,ay,az,affp,we,nx,ny,nz,kstrt,nvpy,nvpz,
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

# find number of particles in each of mx, my, mz tiles:
# updates kpic, nppmx
ppdblkp3l(part,kpic,npp,noff,nppmx,idimp,npmax,mx,my,mz,mx1,myp1,mxyzp1,
          idds,irc)
if (irc[0] != 0):
   print "ppdblkp3l error, irc=", irc[0]
   ppabort()
   exit(0)

# allocate vector particle data
nppmx0 = int((1.0 + xtras)*nppmx)
ntmaxp = int(xtras*nppmx)
npbmx = int(xtras*nppmx)
nbmaxp = int(0.125*mxzyp1*npbmx)
# sbufl/sbufr = particle buffers sent to nearby processors
sbufl = numpy.empty((idimp,nbmaxp,2),float_type,'F')
sbufr = numpy.empty((idimp,nbmaxp,2),float_type,'F')
# rbufl/rbufr = particle buffers received from nearby processors
rbufl = numpy.empty((idimp,nbmaxp,2),float_type,'F')
rbufr = numpy.empty((idimp,nbmaxp,2),float_type,'F')
# ppart = tiled particle array
ppart = numpy.empty((idimp,nppmx0,mxyzp1),float_type,'F')
# ppbuff = buffer array for reordering tiled particle array
ppbuff = numpy.empty((idimp,npbmx,mxyzp1),float_type,'F')
# ncl = number of particles departing tile in each direction
ncl = numpy.empty((26,mxyzp1),int_type,'F')
# iholep = location/destination of each particle departing tile
iholep = numpy.empty((2,ntmaxp+1,mxyzp1),int_type,'F')
# ncll/nclr/mcll/mclr = number offsets send/received from processors
ncll = numpy.empty((3,mxzyp1,3,2),int_type,'F')
nclr = numpy.empty((3,mxzyp1,3,2),int_type,'F')
mcll = numpy.empty((3,mxzyp1,3,2),int_type,'F')
mclr = numpy.empty((3,mxzyp1,3,2),int_type,'F')
mcls = numpy.empty((3,mx1+1,4),int_type,'F')

# copy ordered particle data for OpenMP: updates ppart and kpic
pppmovin3l(part,ppart,kpic,npp,noff,nppmx0,idimp,npmax,mx,my,mz,mx1,
           myp1,mxyzp1,idds,irc)
if (irc[0] != 0):
   print "pppmovin3l overflow error, irc=", irc[0]
   ppabort()
   exit(0)
# sanity check
pppcheck3l(ppart,kpic,noff,nyzp,idimp,nppmx0,nx,mx,my,mz,mx1,myp1,mzp1,
           idds,irc)
if (irc[0] != 0):
   print "pppcheck2l error, irc=", irc[0]
   ppabort()
   exit(0)

# find maximum and minimum initial electron density
qe.fill(0.0)
ppgppost32l(ppart,qe,kpic,noff,qme,nppmx0,idimp,mx,my,mz,nxe,nypmx,
            nzpmx,mx1,myp1,mxyzp1,idds)
ppaguard32xl(qe,nyzp,nx,nxe,nypmx,nzpmx,idds)
ppnaguard32l(qe,scs,scr,nyzp,kstrt,nvpy,nvpz,nx,nxe,nypmx,nzpmx,idds)
ppfwpminmx32(qe,nyzp,qbme,wpmax,wpmin,nx,nxe,nypmx,nzpmx,idds)
wtot[0] = wpmax[0]
wtot[1] = -wpmin[0]
ppdmax(wtot,work,2)
wpmax[0] = wtot[0]
wpmin[0] = -wtot[1]
wpm = 0.5*(wpmax[0] + wpmin[0])*affp
# accelerate convergence: update wpm
if (wpm <= 10.0):
   wpm = 0.75*wpm
if (kstrt==1):
   print "wpm=",wpm
q2m0 = wpm/affp
# calculate form factor: ffe
isign = 0
mppepoisp332(dcut,exyzt,isign,ffe,ax,ay,az,affp,wpm,ci,wf,nx,ny,nz,
             kstrt,nvpy,nvpz,nze,kxyp,kyzp,nzh)

# initialize transverse electric field
cus.fill(0.0)

# * * * start main iteration loop * * *

for ntime in xrange(0,nloop):
#  if (kstrt==1):
#     print "ntime = ", ntime

# deposit current with OpenMP: updates cue
   dtimer(dtime,itime,-1)
   cue.fill(0.0)
   ppgjppost32l(ppart,cue,kpic,noff,qme,zero,nppmx0,idimp,nx,ny,nz,mx,
                my,mz,nxe,nypmx,nzpmx,mx1,myp1,mxyzp1,idds,ipbc)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tdjpost = tdjpost + time

# deposit charge with OpenMP: updates qe
   dtimer(dtime,itime,-1)
   qe.fill(0.0)
   ppgppost32l(ppart,qe,kpic,noff,qme,nppmx0,idimp,mx,my,mz,nxe,nypmx,
               nzpmx,mx1,myp1,mxyzp1,idds)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tdpost = tdpost + time

# add guard cells with OpenMP: updates qe, cue
   dtimer(dtime,itime,-1)
   ppaguard32xl(qe,nyzp,nx,nxe,nypmx,nzpmx,idds)
   ppnaguard32l(qe,scs,scr,nyzp,kstrt,nvpy,nvpz,nx,nxe,nypmx,nzpmx,idds)
   ppacguard32xl(cue,nyzp,nx,ndim,nxe,nypmx,nzpmx,idds)
   ppnacguard32l(cue,scs,scr,nyzp,ndim,kstrt,nvpy,nvpz,nx,nxe,nypmx,
                 nzpmx,idds)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tguard = tguard + time

# transform charge to fourier space with OpenMP: updates qt, modifies qe
   dtimer(dtime,itime,-1)
   isign = -1
   wppfft32rm(qe,qs,qt,bs,br,isign,ntpose,mixup,sct,ttp,indx,indy,indz,
              kstrt,nvpy,nvpz,nxeh,nye,nze,kxyp,kyp,kyzp,kzp,kxyp,nypmx,
              kyzp,nzpmx,kzyp,nxhyz,nxyzh)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tfft[0] = tfft[0] + time
   tfft[1] = tfft[1] + ttp[0]

# calculate longitudinal force/charge in fourier space with OpenMP:
# updates fxyzt, we
   dtimer(dtime,itime,-1)
   isign = -1
   mppois332(qt,fxyzt,isign,ffc,ax,ay,az,affp,we,nx,ny,nz,kstrt,nvpy,
             nvpz,nze,kxyp,kyzp,nzh)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tfield = tfield + time

# transform longitudinal electric force to real space with OpenMP:
# updates fxyze, modifies fxyzt
   dtimer(dtime,itime,-1)
   isign = 1
   wppfft32rm3(fxyze,fxyzs,fxyzt,bs,br,isign,ntpose,mixup,sct,ttp,indx,
               indy,indz,kstrt,nvpy,nvpz,nxeh,nye,nze,kxyp,kyp,kyzp,kzp,
               kxyp,nypmx,kyzp,nzpmx,kzyp,nxhyz,nxyzh)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tfft[0] = tfft[0] + time
   tfft[1] = tfft[1] + ttp[0]

# transform current to fourier space with OpenMP: update cut
# modifies cue
   dtimer(dtime,itime,-1)
   isign = -1
   wppfft32rm3(cue,fxyzs,cut,bs,br,isign,ntpose,mixup,sct,ttp,indx,indy,
               indz,kstrt,nvpy,nvpz,nxeh,nye,nze,kxyp,kyp,kyzp,kzp,kxyp,
               nypmx,kyzp,nzpmx,kzyp,nxhyz,nxyzh)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tfft[0] = tfft[0] + time
   tfft[1] = tfft[1] + ttp[0]

# take transverse part of current with OpenMP: updates cut
   dtimer(dtime,itime,-1)
   mppcuperp32(cut,nx,ny,nz,kstrt,nvpy,nvpz,nze,kxyp,kyzp)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tfield = tfield + time

# calculate magnetic field in fourier space with OpenMP:
# updates bxyzt, wm
   dtimer(dtime,itime,-1)
   mppbbpoisp332(cut,bxyzt,ffc,ci,wm,nx,ny,nz,kstrt,nvpy,nvpz,nze,kxyp,
                 kyzp,nzh)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tfield = tfield + time

# transform magnetic force to real space with OpenMP: updates bxyze
# modifies bxyzt
   dtimer(dtime,itime,-1)
   isign = 1
   wppfft32rm3(bxyze,fxyzs,bxyzt,bs,br,isign,ntpose,mixup,sct,ttp,indx,
               indy,indz,kstrt,nvpy,nvpz,nxeh,nye,nze,kxyp,kyp,kyzp,kzp,
               kxyp,nypmx,kyzp,nzpmx,kzyp,nxhyz,nxyzh)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tfft[0] = tfft[0] + time
   tfft[1] = tfft[1] + ttp[0]

# add constant to magnetic field with OpenMP: updates bxyze
   dtimer(dtime,itime,-1)
   mppbaddext32(bxyze,nyzp,omx,omy,omz,nx,nxe,nypmx,nzpmx,idds)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tfield = tfield + time

# copy guard cells with OpenMP: updates fxyze, bxyze
   dtimer(dtime,itime,-1)
   ppncguard32l(fxyze,scs,nyzp,kstrt,nvpy,nvpz,nnxe,nypmx,nzpmx,idds)
   ppcguard32xl(fxyze,nyzp,nx,ndim,nxe,nypmx,nzpmx,idds)
   ppncguard32l(bxyze,scs,nyzp,kstrt,nvpy,nvpz,nnxe,nypmx,nzpmx,idds)
   ppcguard32xl(bxyze,nyzp,nx,ndim,nxe,nypmx,nzpmx,idds)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tguard = tguard + time

# add longitudinal and old transverse electric fields with OpenMP:
# updates exyze
   dtimer(dtime,itime,-1)
   mppaddvrfield32(exyze,cus,fxyze,ndim,nxe,nypmx,nzpmx)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tfield = tfield + time

# deposit electron acceleration density and momentum flux with OpenMP:
# updates dcu, amu
   dtimer(dtime,itime,-1)
   dcu.fill(0.0); amu.fill(0.0)
   ppgdjppost32l(ppart,exyze,bxyze,kpic,noff,nyzp,dcu,amu,qme,qbme,dt,
                 idimp,nppmx0,nx,mx,my,mz,nxe,nypmx,nzpmx,mx1,myp1,
                 mxyzp1,idds)
# add old scaled electric field with standard procedure: updates dcu
   mppascfguard32l(dcu,cus,nyzp,q2m0,nx,nxe,nypmx,nzpmx,idds)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tdcjpost = tdcjpost + time

# add guard cells with OpenMP: updates dcu, amu
   dtimer(dtime,itime,-1)
   ppacguard32xl(dcu,nyzp,nx,ndim,nxe,nypmx,nzpmx,idds)
   ppnacguard32l(dcu,scs,scr,nyzp,ndim,kstrt,nvpy,nvpz,nx,nxe,nypmx,
                 nzpmx,idds)
   ppacguard32xl(amu,nyzp,nx,mdim,nxe,nypmx,nzpmx,idds)
   ppnacguard32l(amu,scs,scr,nyzp,mdim,kstrt,nvpy,nvpz,nx,nxe,nypmx,
                 nzpmx,idds)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tguard = tguard + time

# transform acceleration density and momentum flux to fourier space
# with OpenMP: updates dcut, amut, modifies dcu, amu
   dtimer(dtime,itime,-1)
   isign = -1
   wppfft32rm3(dcu,fxyzs,dcut,bs,br,isign,ntpose,mixup,sct,ttp,indx,
               indy,indz,kstrt,nvpy,nvpz,nxeh,nye,nze,kxyp,kyp,kyzp,kzp,
               kxyp,nypmx,kyzp,nzpmx,kzyp,nxhyz,nxyzh)
   tfft[1] = tfft[1] + ttp[0]
   wppfft32rmn(amu,amus,amut,bs,br,ss,isign,ntpose,mixup,sct,ttp,indx,
               indy,indz,kstrt,nvpy,nvpz,nxeh,nye,nze,kxyp,kyp,kyzp,kzp,
               kxyp,nypmx,kyzp,nzpmx,kzyp,mdim,nxhyz,nxyzh)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tfft[0] = tfft[0] + time
   tfft[1] = tfft[1] + ttp[0]

# take transverse part of time derivative of current with OpenMP:
# updates dcut
   dtimer(dtime,itime,-1)
   mppadcuperp32(dcut,amut,nx,ny,nz,kstrt,nvpy,nvpz,nze,kxyp,kyzp)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tfield = tfield + time

# calculate transverse electric field with OpenMP:
# updates exyzt, wf
   dtimer(dtime,itime,-1)
   isign = -1
   mppepoisp332(dcut,exyzt,isign,ffe,ax,ay,az,affp,wpm,ci,wf,nx,ny,nz,
                kstrt,nvpy,nvpz,nze,kxyp,kyzp,nzh)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tfield = tfield + time

# transform transverse electric field to real space with OpenMP:
# updates cus, modifies exyzt
   dtimer(dtime,itime,-1)
   isign = 1
   wppfft32rm3(cus,fxyzs,exyzt,bs,br,isign,ntpose,mixup,sct,ttp,indx,
               indy,indz,kstrt,nvpy,nvpz,nxeh,nye,nze,kxyp,kyp,kyzp,kzp,
               kxyp,nypmx,kyzp,nzpmx,kzyp,nxhyz,nxyzh)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tfft[0] = tfft[0] + time
   tfft[1] = tfft[1] + ttp[0]

# copy guard cells with OpenMP: updates cus
   dtimer(dtime,itime,-1)
   ppncguard32l(cus,scs,nyzp,kstrt,nvpy,nvpz,nnxe,nypmx,nzpmx,idds)
   ppcguard32xl(cus,nyzp,nx,ndim,nxe,nypmx,nzpmx,idds)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tguard = tguard + time

# add longitudinal and transverse electric fields with OpenMP:
# exyze = cus + fxyze, updates exyze
# cus needs to be retained for next time step
   dtimer(dtime,itime,-1)
   mppaddvrfield32(exyze,cus,fxyze,ndim,nxe,nypmx,nzpmx)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tfield = tfield + time

# inner iteration loop
   for k in xrange(0,ndc):

# deposit electron current and acceleration density and momentum flux
# with OpenMP: updates cue, dcu, amu
      dtimer(dtime,itime,-1)
      cue.fill(0.0); dcu.fill(0.0); amu.fill(0.0)
      ppgdcjppost32l(ppart,exyze,bxyze,kpic,noff,nyzp,cue,dcu,amu,qme,
                     qbme,dt,idimp,nppmx0,nx,mx,my,mz,nxe,nypmx,nzpmx,
                     mx1,myp1,mxyzp1,idds)
# add scaled electric field with standard procedure: updates dcu
      mppascfguard32l(dcu,cus,nyzp,q2m0,nx,nxe,nypmx,nzpmx,idds)
      dtimer(dtime,itime,1)
      time = float(dtime)
      tdcjpost = tdcjpost + time

# add guard cells for current, acceleration density, and momentum flux
# with OpenMP: updates cue, dcu, amu
      dtimer(dtime,itime,-1)
      ppacguard32xl(cue,nyzp,nx,ndim,nxe,nypmx,nzpmx,idds)
      ppnacguard32l(cue,scs,scr,nyzp,ndim,kstrt,nvpy,nvpz,nx,nxe,nypmx,
                    nzpmx,idds)
      ppacguard32xl(dcu,nyzp,nx,ndim,nxe,nypmx,nzpmx,idds)
      ppnacguard32l(dcu,scs,scr,nyzp,ndim,kstrt,nvpy,nvpz,nx,nxe,nypmx,
                    nzpmx,idds)
      ppacguard32xl(amu,nyzp,nx,mdim,nxe,nypmx,nzpmx,idds)
      ppnacguard32l(amu,scs,scr,nyzp,mdim,kstrt,nvpy,nvpz,nx,nxe,nypmx,
                    nzpmx,idds)
      dtimer(dtime,itime,1)
      time = float(dtime)
      tguard = tguard + time

# transform current to fourier space with OpenMP: update cut
# modifies cue
      dtimer(dtime,itime,-1)
      isign = -1
      wppfft32rm3(cue,fxyzs,cut,bs,br,isign,ntpose,mixup,sct,ttp,indx,
                  indy,indz,kstrt,nvpy,nvpz,nxeh,nye,nze,kxyp,kyp,kyzp,
                  kzp,kxyp,nypmx,kyzp,nzpmx,kzyp,nxhyz,nxyzh)
      dtimer(dtime,itime,1)
      time = float(dtime)
      tfft[0] = tfft[0] + time
      tfft[1] = tfft[1] + ttp[0]

# take transverse part of current with OpenMP: updates cut
      dtimer(dtime,itime,-1)
      mppcuperp32(cut,nx,ny,nz,kstrt,nvpy,nvpz,nze,kxyp,kyzp)
      dtimer(dtime,itime,1)
      time = float(dtime)
      tfield = tfield + time

# calculate magnetic field in fourier space with OpenMP:
# updates bxyzt, wm
      dtimer(dtime,itime,-1)
      mppbbpoisp332(cut,bxyzt,ffc,ci,wm,nx,ny,nz,kstrt,nvpy,nvpz,nze,
                    kxyp,kyzp,nzh)
      dtimer(dtime,itime,1)
      time = float(dtime)
      tfield = tfield + time

# transform magnetic force to real space with OpenMP: updates bxyze
# modifies bxyzt
      dtimer(dtime,itime,-1)
      isign = 1
      wppfft32rm3(bxyze,fxyzs,bxyzt,bs,br,isign,ntpose,mixup,sct,ttp,
                  indx,indy,indz,kstrt,nvpy,nvpz,nxeh,nye,nze,kxyp,kyp,
                  kyzp,kzp,kxyp,nypmx,kyzp,nzpmx,kzyp,nxhyz,nxyzh)
      dtimer(dtime,itime,1)
      time = float(dtime)
      tfft[0] = tfft[0] + time
      tfft[1] = tfft[1] + ttp[0]

# add constant to magnetic field with OpenMP: updates bxyze
      dtimer(dtime,itime,-1)
      mppbaddext32(bxyze,nyzp,omx,omy,omz,nx,nxe,nypmx,nzpmx,idds)
      dtimer(dtime,itime,1)
      time = float(dtime)
      tfield = tfield + time

# transform acceleration density and momentum flux to fourier space
# with OpenMP: updates dcut, amut, modifies dcu, amu
      dtimer(dtime,itime,-1)
      isign = -1
      wppfft32rm3(dcu,fxyzs,dcut,bs,br,isign,ntpose,mixup,sct,ttp,indx,
                  indy,indz,kstrt,nvpy,nvpz,nxeh,nye,nze,kxyp,kyp,kyzp,
                  kzp,kxyp,nypmx,kyzp,nzpmx,kzyp,nxhyz,nxyzh)
      tfft[1] = tfft[1] + ttp[0]
      wppfft32rmn(amu,amus,amut,bs,br,ss,isign,ntpose,mixup,sct,ttp,
                  indx,indy,indz,kstrt,nvpy,nvpz,nxeh,nye,nze,kxyp,kyp,
                  kyzp,kzp,kxyp,nypmx,kyzp,nzpmx,kzyp,mdim,nxhyz,nxyzh)
      dtimer(dtime,itime,1)
      time = float(dtime)
      tfft[0] = tfft[0] + time
      tfft[1] = tfft[1] + ttp[0]

# take transverse part of time derivative of current with OpenMP:
# updates dcut
      dtimer(dtime,itime,-1)
      mppadcuperp32(dcut,amut,nx,ny,nz,kstrt,nvpy,nvpz,nze,kxyp,kyzp)
      dtimer(dtime,itime,1)
      time = float(dtime)
      tfield = tfield + time

# calculate transverse electric field with OpenMP:
# updates exyzt, wf
      dtimer(dtime,itime,-1)
      isign = -1
      mppepoisp332(dcut,exyzt,isign,ffe,ax,ay,az,affp,wpm,ci,wf,nx,ny,
                   nz,kstrt,nvpy,nvpz,nze,kxyp,kyzp,nzh)
      dtimer(dtime,itime,1)
      time = float(dtime)
      tfield = tfield + time

# transform transverse electric field to real space with OpenMP:
# updates cus, modifies exyzt
      dtimer(dtime,itime,-1)
      isign = 1
      wppfft32rm3(cus,fxyzs,exyzt,bs,br,isign,ntpose,mixup,sct,ttp,indx,
                  indy,indz,kstrt,nvpy,nvpz,nxeh,nye,nze,kxyp,kyp,kyzp,
                  kzp,kxyp,nypmx,kyzp,nzpmx,kzyp,nxhyz,nxyzh)
      dtimer(dtime,itime,1)
      time = float(dtime)
      tfft[0] = tfft[0] + time
      tfft[1] = tfft[1] + ttp[0]

# copy guard cells with OpenMP: updates bxyze, cus
      dtimer(dtime,itime,-1)
      ppncguard32l(bxyze,scs,nyzp,kstrt,nvpy,nvpz,nnxe,nypmx,nzpmx,idds)
      ppcguard32xl(bxyze,nyzp,nx,ndim,nxe,nypmx,nzpmx,idds)
      ppncguard32l(cus,scs,nyzp,kstrt,nvpy,nvpz,nnxe,nypmx,nzpmx,idds)
      ppcguard32xl(cus,nyzp,nx,ndim,nxe,nypmx,nzpmx,idds)
      dtimer(dtime,itime,1)
      time = float(dtime)
      tguard = tguard + time

# add longitudinal and transverse electric fields with OpenMP:
# exyze = cus + fxyze, updates exyze
# cus needs to be retained for next time step
      dtimer(dtime,itime,-1)
      mppaddvrfield32(exyze,cus,fxyze,ndim,nxe,nypmx,nzpmx)
      dtimer(dtime,itime,1)
      time = float(dtime)
      tfield = tfield + time
   pass

# push particles with OpenMP: updates part, wke, and ihole
   wke[0] = 0.0
   dtimer(dtime,itime,-1)
# updates ppart, wke
#  ppgbppush32l(ppart,exyze,bxyze,kpic,noff,nyzp,qbme,dt,dt,wke,idimp,
#               nppmx0,nx,ny,nz,mx,my,mz,nxe,nypmx,nzpmx,mx1,myp1,
#               mxyzp1,idds,ipbc)
# updates ppart, ncl, iholep, wke, irc
   ppgbppushf32l(ppart,exyze,bxyze,kpic,ncl,iholep,noff,nyzp,qbme,dt,dt,
                 wke,idimp,nppmx0,nx,ny,nz,mx,my,mz,nxe,nypmx,nzpmx,mx1,
                 myp1,mxyzp1,ntmaxp,idds,irc)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tpush = tpush + time
   if (irc[0] != 0):
      print kstrt, "ppgbppushf32l error: irc=", irc[0]
      ppabort()
      exit(0)

# reorder particles by tile with OpenMP
# first part of particle reorder on x, y and z cell
# with mx, my, mz tiles:
   dtimer(dtime,itime,-1)
# updates ppart, ppbuff, sbufl, sbufr, ncl, iholep, ncll, nclr, irc
#  ppporder32la(ppart,ppbuff,sbufl,sbufr,kpic,ncl,iholep,ncll,nclr,noff,
#               nyzp,idimp,nppmx0,nx,ny,nz,mx,my,mz,mx1,myp1,mzp1,
#               mxzyp1,npbmx,ntmaxp,nbmaxp,idds,irc)
# updates: ppart, ppbuff, sbufl, sbufr, ncl, ncll, nclr, irc
   ppporderf32la(ppart,ppbuff,sbufl,sbufr,ncl,iholep,ncll,nclr,idimp,
                 nppmx0,mx1,myp1,mzp1,mxzyp1,npbmx,ntmaxp,nbmaxp,irc)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tsort = tsort + time
   if (irc[0] != 0):
      print kstrt, "ppporderf32la error: ntmaxp, irc=", ntmaxp, irc[0]
      ppabort()
      exit(0)
# move particles into appropriate spatial regions:
# updates rbufr, rbufl, mcll, mclr, mcls
   dtimer(dtime,itime,-1)
   pppmove32(sbufr,sbufl,rbufr,rbufl,ncll,nclr,mcll,mclr,mcls,kstrt,
             nvpy,nvpz,idimp,nbmaxp,mx1,myp1,mzp1,mxzyp1,irc)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tmov = tmov + time
   if (irc[0] != 0):
      if (kstrt==1):
         print kstrt, "pppmove32 error: nbmaxp, irc=", nbmaxp, irc[0]
         ppabort()
         exit(0)

# second part of particle reorder on x and y cell with mx, my, mz tiles:
# updates ppart, kpic
   dtimer(dtime,itime,-1)
   ppporder32lb(ppart,ppbuff,rbufl,rbufr,kpic,ncl,iholep,mcll,mclr,mcls,
                idimp,nppmx0,mx1,myp1,mzp1,mxzyp1,npbmx,ntmaxp,nbmaxp,
                irc)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tsort = tsort + time
   if (irc[0] != 0):
      print kstrt, "ppporder32lb error: nppmx0, irc=", nppmx0, irc[0]
      ppabort()
      exit(0)

# energy diagnostic
   wt = we + wm
   wtot[0] = wt
   wtot[1] = wke[0]
   wtot[2] = 0.0
   wtot[3] = wke[0] + wt
   wtot[4] = we[0]
   wtot[5] = wf[0]
   wtot[6] = wm[0]
   ppdsum(wtot,work,7)
   wke[0] = wtot[1]
   we[0] = wtot[4]
   wf[0] = wtot[5]
   wm[0] = wtot[6]
   if (ntime==0):
      if (kstrt==1):
         wt = we + wm
         print "Initial Total Field, Kinetic and Total Energies:"
         print "%14.7e %14.7e %14.7e" % (wt, wke, wke + wt)
         print "Initial Electrostatic, Transverse Electric and Magnetic " \
         "Field Energies:"
         print "%14.7e %14.7e %14.7e" % (we, wf, wm)
ntime = ntime + 1

# * * * end main iteration loop * * *

if (kstrt==1):
   print "ntime, ndc = ", ntime, ndc
   print "MPI nodes nvpy, nvpz = ", nvpy, nvpz
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

