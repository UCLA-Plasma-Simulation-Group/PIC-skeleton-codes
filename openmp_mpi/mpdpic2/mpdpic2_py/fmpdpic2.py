#-----------------------------------------------------------------------
# Skeleton 2-1/2D Darwin MPI/OpenMP PIC code
# written by Viktor K. Decyk, Adam Tableman, and Qiyang Hu, UCLA
import math
import numpy
from fmpdpush2 import *
from dtimer import *
from fomplib import *

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
# idimp = dimension of phase space = 5
# ipbc = particle boundary condition: 1 = periodic
idimp = 5; ipbc = 1
# omx/omy/omz = magnetic field electron cyclotron frequency in x/y/z 
omx = 0.4; omy = 0.0; omz = 0.0
# ndc = number of corrections in darwin iteration
ndc = 1
# idps = number of partition boundaries
idps = 2
# wke/we = particle kinetic/electrostatic field energy
# wf/wm/wt = magnetic field/transverse electric field/total energy
wke = numpy.zeros((1),float_type)
we = numpy.zeros((1),float_type)
wf = numpy.zeros((1),float_type)
wm = numpy.zeros((1),float_type)
wt = numpy.zeros((1),float_type)
zero = 0.0
# sorting tiles, should be less than or equal to 32
mx = 16; my = 16
# fraction of extra particles needed for particle management
xtras = 0.2
# declare scalars for standard code
wpmax = numpy.empty((1),float_type)
wpmin = numpy.empty((1),float_type)

# declare scalars for MPI code
ntpose = 1
nnvp = numpy.empty((1),int_type)
nidproc = numpy.empty((1),int_type)
nnyp = numpy.empty((1),int_type)
nnoff = numpy.empty((1),int_type)
nnpp = numpy.empty((1),int_type)
nnypmx = numpy.empty((1),int_type)
nnypmn = numpy.empty((1),int_type)
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
np =  float(npx)*float(npy)
# nx/ny = number of grid points in x/y direction
nx = int(math.pow(2,indx)); ny = int(math.pow(2,indy))
nxh = int(nx/2); nyh = max(1,int(ny/2))
nxe = nx + 2; nye = ny + 2; nxeh = int(nxe/2); nnxe = ndim*nxe
nxyh = int(max(nx,ny)/2); nxhy = max(nxh,ny)
# mx1 = number of tiles in x direction
mx1 = int((nx - 1)/mx) + 1
# nloop = number of time steps in simulation
# ntime = current time step
nloop = int(tend/dt + .0001); ntime = 0
# mdim = dimension of amu array
mdim = 2*ndim - 2
qbme = qme
affp = float(nx)*float(ny)/np

# nvp = number of distributed memory nodes
# initialize for distributed memory parallel processing
ppinit2(nidproc,nnvp)
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
pdicomp2l(edges,nnyp,nnoff,nnypmx,nnypmn,ny,kstrt,nvp,idps)
nyp = nnyp[0]; noff = nnoff[0]; nypmx = nnypmx[0]; nypmn = nnypmn[0]
if (nypmn < 1):
   if (kstrt==1):
      print "combination not supported nvp, ny =",nvp,ny
   ppexit()
   exit(0)

# initialize additional scalars for MPI code
# kxp = number of complex grids in each field partition in x direction
kxp = int((nxh - 1)/nvp) + 1
# kyp = number of complex grids in each field partition in y direction
kyp = int((ny - 1)/nvp) + 1
# npmax = maximum number of electrons in each partition
npmax = int((np/float(nvp))*1.25)
# myp1 = number of tiles in y direction
myp1 = int((nyp - 1)/my) + 1; mxyp1 = mx1*myp1

# allocate data for standard code
# part = particle array
part = numpy.empty((idimp,npmax),float_type,'F')
# qe = electron charge density with guard cells
qe = numpy.empty((nxe,nypmx),float_type,'F')
# cue = electron current density with guard cells
cue = numpy.empty((ndim,nxe,nypmx),float_type,'F')
# dcu = acceleration density with guard cells
dcu = numpy.empty((ndim,nxe,nypmx),float_type,'F')
# cus = smoothed transverse electric field with guard cells
cus = numpy.empty((ndim,nxe,nypmx),float_type,'F')
# amu = momentum flux with guard cells
amu = numpy.empty((mdim,nxe,nypmx),float_type,'F')
# exyze = smoothed total electric field with guard cells
exyze = numpy.empty((ndim,nxe,nypmx),float_type,'F')
# fxyze = smoothed longitudinal electric field with guard cells
fxyze = numpy.empty((ndim,nxe,nypmx),float_type,'F')
# bxyze = smoothed magnetic field with guard cells
bxyze = numpy.empty((ndim,nxe,nypmx),float_type,'F')
# ss = scratch array for WPPFFT2RN
ss = numpy.empty((mdim*nxeh,nypmx),complex_type,'F')
# qt = scalar charge density field array in fourier space
qt = numpy.empty((nye,kxp),complex_type,'F')
# cut = vector current density in fourier space
cut = numpy.empty((ndim,nye,kxp),complex_type,'F')
# dcut = vector acceleration density in fourier space
dcut = numpy.empty((ndim,nye,kxp),complex_type,'F')
# cur = vector transverse electric field in fourier space
cur = numpy.empty((ndim,nye,kxp),complex_type,'F')
# amut = tensor momentum flux in fourier space
amut = numpy.empty((mdim,nye,kxp),complex_type,'F')
# exyt = vector total electric field in fourier space
exyt = numpy.empty((ndim,nye,kxp),complex_type,'F')
# fxyt = vector longitudinal electric field in fourier space
fxyt = numpy.empty((ndim,nye,kxp),complex_type,'F')
# bxyt = vector magnetic field in fourier space
bxyt = numpy.empty((ndim,nye,kxp),complex_type,'F')
# ffc, ffe = form factor arrays for poisson solvers
ffc = numpy.empty((nyh,kxp),complex_type,'F')
ffe = numpy.empty((nyh,kxp),complex_type,'F')
# mixup = bit reverse table for FFT
mixup = numpy.empty((nxhy),int_type,'F')
# sct = sine/cosine table for FFT
sct = numpy.empty((nxyh),complex_type,'F')
# kpic = number of particles in each tile
kpic = numpy.empty((mxyp1),int_type,'F')
wtot = numpy.empty((7),double_type)
work = numpy.empty((7),double_type)

# allocate data for MPI code
# bs/br = complex send/receive buffers for data transpose
bs = numpy.empty((mdim,kxp,kyp),complex_type,'F')
br = numpy.empty((mdim,kxp,kyp),complex_type,'F')
# scs/scr = guard cell buffers received from nearby processors
scs = numpy.empty((mdim*nxe),float_type,'F')
scr = numpy.empty((mdim*nxe),float_type,'F')

# prepare fft tables
wpfft2rinit(mixup,sct,indx,indy,nxhy,nxyh)
# calculate form factor: ffc
isign = 0
mppois23(qt,fxyt,isign,ffc,ax,ay,affp,we,nx,ny,kstrt,nye,kxp,nyh)
# initialize electrons
nps = 1
nnpp[0] = 0
pdistr2h(part,edges,nnpp,nps,vtx,vty,vtz,vx0,vy0,vz0,npx,npy,nx,ny,
         idimp,npmax,idps,ipbc,ierr)
npp = nnpp[0]
# check for particle initialization error
if (ierr[0] != 0):
   if (kstrt==1):
      print "particle initialization error: ierr=", ierr[0]
   ppexit()
   exit(0)

# find number of particles in each of mx, my tiles: updates kpic, nppmx
ppdblkp2l(part,kpic,npp,noff,nppmx,idimp,npmax,mx,my,mx1,mxyp1,irc)
if (irc[0] != 0):
   print "ppdblkp2l error, irc=", irc[0]
   ppabort()
   exit(0)

# allocate vector particle data
nppmx0 = int((1.0 + xtras)*nppmx)
ntmaxp = int(xtras*nppmx)
npbmx = int(xtras*nppmx)
nbmaxp = int(0.25*mx1*npbmx)
# sbufl/sbufr = particle buffers sent to nearby processors
sbufl = numpy.empty((idimp,nbmaxp),float_type,'F')
sbufr = numpy.empty((idimp,nbmaxp),float_type,'F')
# rbufl/rbufr = particle buffers received from nearby processors
rbufl = numpy.empty((idimp,nbmaxp),float_type,'F')
rbufr = numpy.empty((idimp,nbmaxp),float_type,'F')
# ppart = tiled particle array
ppart = numpy.empty((idimp,nppmx0,mxyp1),float_type,'F')
# ppbuff = buffer array for reordering tiled particle array
ppbuff = numpy.empty((idimp,npbmx,mxyp1),float_type,'F')
# ncl = number of particles departing tile in each direction
ncl = numpy.empty((8,mxyp1),int_type,'F')
# iholep = location/destination of each particle departing tile
iholep = numpy.empty((2,ntmaxp+1,mxyp1),int_type,'F')
# ncll/nclr/mcll/mclr = number offsets send/received from processors
ncll = numpy.empty((3,mx1),int_type,'F')
nclr = numpy.empty((3,mx1),int_type,'F')
mcll = numpy.empty((3,mx1),int_type,'F')
mclr = numpy.empty((3,mx1),int_type,'F')

# copy ordered particle data for OpenMP
pppmovin2l(part,ppart,kpic,npp,noff,nppmx0,idimp,npmax,mx,my,mx1,mxyp1,
           irc)
if (irc[0] != 0):
   print "pppmovin2l overflow error, irc=", irc[0]
   ppabort()
   exit(0)

# sanity check
pppcheck2l(ppart,kpic,noff,nyp,idimp,nppmx0,nx,mx,my,mx1,myp1,irc)
if (irc[0] != 0):
   print "pppcheck2l error, irc=", irc[0]
   ppabort()
   exit(0)

# find maximum and minimum initial electron density
qe.fill(0.0)
ppgppost2l(ppart,qe,kpic,noff,qme,idimp,nppmx0,mx,my,nxe,nypmx,mx1,
           mxyp1)
ppaguard2xl(qe,nyp,nx,nxe,nypmx)
ppnaguard2l(qe,scr,nyp,nx,kstrt,nvp,nxe,nypmx)
ppfwpminmx2(qe,nyp,qbme,wpmax,wpmin,nx,nxe,nypmx)
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
mppepoisp23(dcut,cur,isign,ffe,ax,ay,affp,wpm,ci,wf,nx,ny,kstrt,nye,kxp,
            nyh)

# initialize transverse electric field
cus.fill(0.0)

# * * * start main iteration loop * * *

for ntime in xrange(0,nloop):
#  if (kstrt==1):
#     print "ntime = ", ntime

# deposit current with OpenMP: updates cue
   dtimer(dtime,itime,-1)
   cue.fill(0.0)
   ppgjppost2l(ppart,cue,kpic,noff,qme,zero,nppmx0,idimp,nx,ny,mx,my,
               nxe,nypmx,mx1,mxyp1,ipbc)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tdjpost = tdjpost + time

# deposit charge with OpenMP: updates qe
   dtimer(dtime,itime,-1)
   qe.fill(0.0)
   ppgppost2l(ppart,qe,kpic,noff,qme,idimp,nppmx0,mx,my,nxe,nypmx,mx1,
              mxyp1)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tdpost = tdpost + time

# add guard cells with OpenMP: updates cue, qe
   dtimer(dtime,itime,-1)
   ppaguard2xl(qe,nyp,nx,nxe,nypmx)
   ppnaguard2l(qe,scr,nyp,nx,kstrt,nvp,nxe,nypmx)
   ppacguard2xl(cue,nyp,nx,ndim,nxe,nypmx)
   ppnacguard2l(cue,scr,nyp,nx,ndim,kstrt,nvp,nxe,nypmx)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tguard = tguard + time

# transform charge to fourier space with OpenMP: updates qt, modifies qe
   dtimer(dtime,itime,-1)
   isign = -1
   wppfft2rm(qe,qt,bs,br,isign,ntpose,mixup,sct,ttp,indx,indy,kstrt,nvp,
             nxeh,nye,kxp,kyp,nypmx,nxhy,nxyh)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tfft[0] = tfft[0] + time
   tfft[1] = tfft[1] + ttp[0]

# calculate longitudinal force/charge in fourier space with OpenMP:
# updates fxyt, we
   dtimer(dtime,itime,-1)
   isign = -1
   mppois23(qt,fxyt,isign,ffc,ax,ay,affp,we,nx,ny,kstrt,nye,kxp,nyh)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tfield = tfield + time

# transform longitudinal electric force to real space with OpenMP:
# updates fxyze, modifies fxyt
   dtimer(dtime,itime,-1)
   isign = 1
   wppfft2rm3(fxyze,fxyt,bs,br,isign,ntpose,mixup,sct,ttp,indx,indy,
              kstrt,nvp,nxeh,nye,kxp,kyp,nypmx,nxhy,nxyh)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tfft[0] = tfft[0] + time
   tfft[1] = tfft[1] + ttp[0]

# transform current to fourier space with OpenMP: updates cut
# modifies cue
   dtimer(dtime,itime,-1)
   isign = -1
   wppfft2rm3(cue,cut,bs,br,isign,ntpose,mixup,sct,ttp,indx,indy,kstrt,
              nvp,nxeh,nye,kxp,kyp,nypmx,nxhy,nxyh)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tfft[0] = tfft[0] + time
   tfft[1] = tfft[1] + ttp[0]

# take transverse part of current with OpenMP: updates cut
   dtimer(dtime,itime,-1)
   mppcuperp2(cut,nx,ny,kstrt,nye,kxp)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tfield = tfield + time

# calculate magnetic field in fourier space with standard procedure:
# updates bxyt, wm
   dtimer(dtime,itime,-1)
   mppbbpoisp23(cut,bxyt,ffc,ci,wm,nx,ny,kstrt,nye,kxp,nyh)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tfield = tfield + time

# transform magnetic field to real space with OpenMP: updates bxyze
# modifies bxyt
   dtimer(dtime,itime,-1)
   isign = 1
   wppfft2rm3(bxyze,bxyt,bs,br,isign,ntpose,mixup,sct,ttp,indx,indy,
              kstrt,nvp,nxeh,nye,kxp,kyp,nypmx,nxhy,nxyh)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tfft[0] = tfft[0] + time
   tfft[1] = tfft[1] + ttp[0]

# add constant to magnetic field with standard procedure: updates bxyze
   dtimer(dtime,itime,-1)
   ppbaddext2(bxyze,nyp,omx,omy,omz,nx,nxe,nypmx)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tfield = tfield + time

# copy guard cells with OpenMP: updates fxyze, bxyze
   dtimer(dtime,itime,-1)
   ppncguard2l(fxyze,nyp,kstrt,nvp,nnxe,nypmx)
   ppcguard2xl(fxyze,nyp,nx,ndim,nxe,nypmx)
   ppncguard2l(bxyze,nyp,kstrt,nvp,nnxe,nypmx)
   ppcguard2xl(bxyze,nyp,nx,ndim,nxe,nypmx)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tguard = tguard + time

# add longitudinal and old transverse electric fields with standard
# procedure: updates exyze
   dtimer(dtime,itime,-1)
   ppaddvrfield2(exyze,cus,fxyze,ndim,nxe,nypmx)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tfield = tfield + time

# deposit electron acceleration density and momentum flux with standard
# procedure: updates dcu, amu
   dtimer(dtime,itime,-1)
   dcu.fill(0.0); amu.fill(0.0)
   ppgdjppost2l(ppart,exyze,bxyze,dcu,amu,kpic,noff,nyp,qme,qbme,dt,
                idimp,nppmx0,nx,mx,my,nxe,nypmx,mx1,mxyp1)
# add old scaled electric field with standard procedure: updates dcu
   ppascfguard2l(dcu,cus,nyp,q2m0,nx,nxe,nypmx)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tdcjpost = tdcjpost + time

# add guard cells with standard procedure: updates dcu, amu
   dtimer(dtime,itime,-1)
   ppacguard2xl(dcu,nyp,nx,ndim,nxe,nypmx)
   ppnacguard2l(dcu,scr,nyp,nx,ndim,kstrt,nvp,nxe,nypmx)
   ppacguard2xl(amu,nyp,nx,mdim,nxe,nypmx)
   ppnacguard2l(amu,scr,nyp,nx,mdim,kstrt,nvp,nxe,nypmx)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tguard = tguard + time

# transform acceleration density and momentum flux to fourier space
# with OpenMP: updates dcut, amut, modifies dcu, amu
   dtimer(dtime,itime,-1)
   isign = -1
   wppfft2rm3(dcu,dcut,bs,br,isign,ntpose,mixup,sct,ttp,indx,indy,kstrt,
              nvp,nxeh,nye,kxp,kyp,nypmx,nxhy,nxyh)
   tfft[1] = tfft[1] + ttp[0]
   wppfft2rmn(amu,amut,bs,br,ss,isign,ntpose,mixup,sct,ttp,indx,indy,
              kstrt,nvp,nxeh,nye,kxp,kyp,nypmx,mdim,nxhy,nxyh)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tfft[0] = tfft[0] + time
   tfft[1] = tfft[1] + ttp[0]

# take transverse part of time derivative of current with standard
# procedure: updates dcut
   dtimer(dtime,itime,-1)
   mppadcuperp23(dcut,amut,nx,ny,kstrt,nye,kxp)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tfield = tfield + time

# calculate transverse electric field with standard procedure:
# updates cur, wf
   dtimer(dtime,itime,-1)
   isign = -1
   mppepoisp23(dcut,cur,isign,ffe,ax,ay,affp,wpm,ci,wf,nx,ny,kstrt,nye,
               kxp,nyh)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tfield = tfield + time

# transform transverse electric field to real space with OpenMP:
# updates cus, modifies cur
   dtimer(dtime,itime,-1)
   isign = 1
   wppfft2rm3(cus,cur,bs,br,isign,ntpose,mixup,sct,ttp,indx,indy,kstrt,
              nvp,nxeh,nye,kxp,kyp,nypmx,nxhy,nxyh)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tfft[0] = tfft[0] + time
   tfft[1] = tfft[1] + ttp[0]

# copy guard cells with OpenMP: updates cus
   dtimer(dtime,itime,-1)
   ppncguard2l(cus,nyp,kstrt,nvp,nnxe,nypmx)
   ppcguard2xl(cus,nyp,nx,ndim,nxe,nypmx)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tguard = tguard + time

# add longitudinal and transverse electric fields with standard
# procedure: exyze = cus + fxyze, updates exyze
# cus needs to be retained for next time step
   dtimer(dtime,itime,-1)
   ppaddvrfield2(exyze,cus,fxyze,ndim,nxe,nypmx)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tfield = tfield + time

# inner iteration loop
   for k in xrange(0,ndc):

# deposit electron current and acceleration density and momentum flux
# with standard procedure: updates cue, dcu, amu
      dtimer(dtime,itime,-1)
      cue.fill(0.0); dcu.fill(0.0); amu.fill(0.0)
      ppgdcjppost2l(ppart,exyze,bxyze,cue,dcu,amu,kpic,noff,nyp,qme,
                    qbme,dt,idimp,nppmx0,nx,mx,my,nxe,nypmx,mx1,mxyp1)
# add scaled electric field with standard procedure: updates dcu
      ppascfguard2l(dcu,cus,nyp,q2m0,nx,nxe,nypmx)
      dtimer(dtime,itime,1)
      time = float(dtime)
      tdcjpost = tdcjpost + time

# add guard cells for current, acceleration density, and momentum flux
# with OpenMP: updates cue, dcu, amu
      dtimer(dtime,itime,-1)
      ppacguard2xl(cue,nyp,nx,ndim,nxe,nypmx)
      ppnacguard2l(cue,scr,nyp,nx,ndim,kstrt,nvp,nxe,nypmx)
      ppacguard2xl(dcu,nyp,nx,ndim,nxe,nypmx)
      ppnacguard2l(dcu,scr,nyp,nx,ndim,kstrt,nvp,nxe,nypmx)
      ppacguard2xl(amu,nyp,nx,mdim,nxe,nypmx)
      ppnacguard2l(amu,scr,nyp,nx,mdim,kstrt,nvp,nxe,nypmx)
      dtimer(dtime,itime,1)
      time = float(dtime)
      tguard = tguard + time

# transform current to fourier space with OpenMP: updates cut
# modifies cue
      dtimer(dtime,itime,-1)
      isign = -1
      wppfft2rm3(cue,cut,bs,br,isign,ntpose,mixup,sct,ttp,indx,indy,
                 kstrt,nvp,nxeh,nye,kxp,kyp,nypmx,nxhy,nxyh)
      dtimer(dtime,itime,1)
      time = float(dtime)
      tfft[0] = tfft[0] + time
      tfft[1] = tfft[1] + ttp[0]

# take transverse part of current with OpenMP: updates cut
      dtimer(dtime,itime,-1)
      mppcuperp2(cut,nx,ny,kstrt,nye,kxp)
      dtimer(dtime,itime,1)
      time = float(dtime)
      tfield = tfield + time

# calculate magnetic field in fourier space with standard procedure:
# updates bxyt, wm
      dtimer(dtime,itime,-1)
      mppbbpoisp23(cut,bxyt,ffc,ci,wm,nx,ny,kstrt,nye,kxp,nyh)
      dtimer(dtime,itime,1)
      time = float(dtime)
      tfield = tfield + time

# transform magnetic field to real space with OpenMP: updates bxyze
# modifies bxyt
      dtimer(dtime,itime,-1)
      isign = 1
      wppfft2rm3(bxyze,bxyt,bs,br,isign,ntpose,mixup,sct,ttp,indx,indy,
                 kstrt,nvp,nxeh,nye,kxp,kyp,nypmx,nxhy,nxyh)
      dtimer(dtime,itime,1)
      time = float(dtime)
      tfft[0] = tfft[0] + time
      tfft[1] = tfft[1] + ttp[0]

# add constant to magnetic field with standard procedure: updates bxyze
      dtimer(dtime,itime,-1)
      ppbaddext2(bxyze,nyp,omx,omy,omz,nx,nxe,nypmx)
      dtimer(dtime,itime,1)
      time = float(dtime)
      tfield = tfield + time

# transform acceleration density and momentum flux to fourier space
# with OpenMP: updates dcut, amut, modifies dcu, amu
      dtimer(dtime,itime,-1)
      isign = -1
      wppfft2rm3(dcu,dcut,bs,br,isign,ntpose,mixup,sct,ttp,indx,indy,
                 kstrt,nvp,nxeh,nye,kxp,kyp,nypmx,nxhy,nxyh)
      tfft[1] = tfft[1] + ttp[0]
      wppfft2rmn(amu,amut,bs,br,ss,isign,ntpose,mixup,sct,ttp,indx,indy,
                 kstrt,nvp,nxeh,nye,kxp,kyp,nypmx,mdim,nxhy,nxyh)
      dtimer(dtime,itime,1)
      time = float(dtime)
      tfft[0] = tfft[0] + time
      tfft[1] = tfft[1] + ttp[0]

# take transverse part of time derivative of current with standard
# procedure: updates dcut
      dtimer(dtime,itime,-1)
      mppadcuperp23(dcut,amut,nx,ny,kstrt,nye,kxp)
      dtimer(dtime,itime,1)
      time = float(dtime)
      tfield = tfield + time

# calculate transverse electric field with standard procedure:
# updates cur, wf
      dtimer(dtime,itime,-1)
      isign = -1
      mppepoisp23(dcut,cur,isign,ffe,ax,ay,affp,wpm,ci,wf,nx,ny,kstrt,
                  nye,kxp,nyh)
      dtimer(dtime,itime,1)
      time = float(dtime)
      tfield = tfield + time

# transform transverse electric field to real space with OpenMP:
# updates cus, modifies cur
      dtimer(dtime,itime,-1)
      isign = 1
      wppfft2rm3(cus,cur,bs,br,isign,ntpose,mixup,sct,ttp,indx,indy,
                 kstrt,nvp,nxeh,nye,kxp,kyp,nypmx,nxhy,nxyh)
      dtimer(dtime,itime,1)
      time = float(dtime)
      tfft[0] = tfft[0] + time
      tfft[1] = tfft[1] + ttp[0]

# copy guard cells with OpenMP: updates bxyze, cus
      dtimer(dtime,itime,-1)
      ppncguard2l(bxyze,nyp,kstrt,nvp,nnxe,nypmx)
      ppcguard2xl(bxyze,nyp,nx,ndim,nxe,nypmx)
      ppncguard2l(cus,nyp,kstrt,nvp,nnxe,nypmx)
      ppcguard2xl(cus,nyp,nx,ndim,nxe,nypmx)
      dtimer(dtime,itime,1)
      time = float(dtime)
      tguard = tguard + time

# add longitudinal and transverse electric fields with standard
# procedure: exyze = cus + fxyze, updates exyze
# cus needs to be retained for next time step
      dtimer(dtime,itime,-1)
      ppaddvrfield2(exyze,cus,fxyze,ndim,nxe,nypmx)
      dtimer(dtime,itime,1)
      time = float(dtime)
      tfield = tfield + time
   pass

# push particles with OpenMP:
   dtimer(dtime,itime,-1)
   wke[0] = 0.0
# updates ppart and wke
#  ppgbppush23l(ppart,exyze,bxyze,kpic,noff,nyp,qbme,dt,dt,wke,idimp,
#               nppmx0,nx,ny,mx,my,nxe,nypmx,mx1,mxyp1,ipbc)
# updates ppart, ncl, iholep, ek, irc
   ppgbppushf23l(ppart,exyze,bxyze,kpic,ncl,iholep,noff,nyp,qbme,dt,dt,
                 wke,idimp,nppmx0,nx,ny,mx,my,nxe,nypmx,mx1,mxyp1,
                 ntmaxp,irc)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tpush = tpush + time
   if (irc[0] != 0):
      print kstrt, "ppgbppushf23l error: irc=", irc[0]
      ppabort()
      exit(0)

# reorder particles by tile with OpenMP
# first part of particle reorder on x and y cell with mx, my tiles:
   dtimer(dtime,itime,-1)
# updates ppart, ppbuff, sbufl, sbufr, ncl, iholep, ncll, nclr, irc
#  ppporder2la(ppart,ppbuff,sbufl,sbufr,kpic,ncl,iholep,ncll,nclr,noff,
#              nyp,idimp,nppmx0,nx,ny,mx,my,mx1,myp1,npbmx,ntmaxp,
#              nbmaxp,irc)
# updates ppart, ppbuff, sbufl, sbufr, ncl, ncll, nclr, irc
   ppporderf2la(ppart,ppbuff,sbufl,sbufr,ncl,iholep,ncll,nclr,idimp,
                nppmx0,mx1,myp1,npbmx,ntmaxp,nbmaxp,irc)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tsort = tsort + time
   if (irc[0] != 0):
      print kstrt, "ppporderf2la error: ntmaxp, irc=", ntmaxp, irc[0]
      ppabort()
      exit(0)
# move particles into appropriate spatial regions:
# updates rbufr, rbufl, mcll, mclr
   dtimer(dtime,itime,-1)
   pppmove2(sbufr,sbufl,rbufr,rbufl,ncll,nclr,mcll,mclr,kstrt,nvp,idimp,
            nbmaxp,mx1)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tmov = tmov + time
# second part of particle reorder on x and y cell with mx, my tiles:
# updates ppart, kpic
   dtimer(dtime,itime,-1)
   ppporder2lb(ppart,ppbuff,rbufl,rbufr,kpic,ncl,iholep,mcll,mclr,idimp,
               nppmx0,mx1,myp1,npbmx,ntmaxp,nbmaxp,irc)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tsort = tsort + time
   if (irc[0] != 0):
      print kstrt, "ppporder2lb error: nppmx0, irc=", nppmx0, irc[0]
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
   print "MPI nodes nvp = ", nvp
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
