#-----------------------------------------------------------------------
# Skeleton 2-1/2D Darwin OpenMP PIC code
# written by Viktor K. Decyk, Adam Tableman, and Qiyang Hu, UCLA
import math
import numpy
from fmdpush2 import *
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
# idimp = number of particle coordinates = 5
# ipbc = particle boundary condition: 1 = periodic
idimp = 5; ipbc = 1
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
# mx/my = number of grids in x/y in sorting tiles
mx = 16; my = 16
# xtras = fraction of extra particles needed for particle management
xtras = 0.2
# declare scalars for standard code
wpmax = numpy.empty((1),float_type)
wpmin = numpy.empty((1),float_type)

# declare scalars for OpenMP code
nppmx = numpy.empty((1),int_type)
irc = numpy.zeros((1),int_type)

# declare and initialize timing data
itime = numpy.empty((4),numpy.int32)
tdpost = 0.0; tguard = 0.0; tfft = 0.0; tfield = 0.0
tdjpost = 0.0; tdcjpost = 0.0; tpush = 0.0; tsort = 0.0
dtime = numpy.empty((1),double_type)

# nvp = number of shared memory nodes (0=default)
nvp = 0
#nvp = int(input("enter number of nodes: "))
# initialize for shared memory parallel processing
init_omp(nvp)

# initialize scalars for standard code
# np = total number of particles in simulation
# nx/ny = number of grid points in x/y direction
np = npx*npy; nx = int(math.pow(2,indx)); ny = int(math.pow(2,indy))
nxh = int(nx/2); nyh = max(1,int(ny/2))
nxe = nx + 2; nye = ny + 1; nxeh = int(nxe/2)
nxyh = int(max(nx,ny)/2); nxhy = max(nxh,ny)
# mx1/my1 = number of tiles in x/y direction
mx1 = int((nx - 1)/mx + 1); my1 = int((ny - 1)/my + 1); mxy1 = mx1*my1
# nloop = number of time steps in simulation
# ntime = current time step
nloop = int(tend/dt + .0001); ntime = 0
# mdim = dimension of amu array
mdim = 2*ndim - 2
qbme = qme
affp = float(nx*ny)/float(np)

# allocate data for standard code
# part = particle array
part = numpy.empty((idimp,np),float_type,'F')
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
# kpic = number of particles in each tile
kpic = numpy.empty((mxy1),int_type,'F')
# ss = scratch array for WFFT2RN
ss = numpy.empty((mdim*nxeh,nye),complex_type,'F')

# prepare fft tables
wfft2rinit(mixup,sct,indx,indy,nxhy,nxyh)
# calculate form factors: ffc
isign = 0
mpois23(qe,fxyze,isign,ffc,ax,ay,affp,we,nx,ny,nxeh,nye,nxh,nyh)
# initialize electrons
distr2h(part,vtx,vty,vtz,vx0,vy0,vz0,npx,npy,idimp,np,nx,ny,ipbc)

# find number of particles in each of mx, my tiles: updates kpic, nppmx
dblkp2l(part,kpic,nppmx,idimp,np,mx,my,mx1,mxy1,irc)
if (irc[0] != 0):
   print "dblkp2l error, irc=", irc[0]
   exit(0)
# allocate vector particle data
nppmx0 = int((1.0 + xtras)*nppmx)
ntmax = int(xtras*nppmx)
npbmx = int(xtras*nppmx)
# ppart = tiled particle array
ppart = numpy.empty((idimp,nppmx0,mxy1),float_type,'F')
# ppbuff = buffer array for reordering tiled particle array
ppbuff = numpy.empty((idimp,npbmx,mxy1),float_type,'F')
# ncl = number of particles departing tile in each direction
ncl = numpy.empty((8,mxy1),int_type,'F')
# ihole = location/destination of each particle departing tile
ihole = numpy.empty((2,ntmax+1,mxy1),int_type,'F')
# copy ordered particle data for OpenMP: updates ppart and kpic
ppmovin2l(part,ppart,kpic,nppmx0,idimp,np,mx,my,mx1,mxy1,irc)
if (irc[0] != 0):
   print "ppmovin2l overflow error, irc=", irc[0]
   exit(0)
# sanity check
ppcheck2l(ppart,kpic,idimp,nppmx0,nx,ny,mx,my,mx1,my1,irc)
if (irc[0] != 0):
   print "ppcheck2l error, irc=", irc[0]
   exit(0)

# find maximum and minimum initial electron density
qe.fill(0.0)
gppost2l(ppart,qe,kpic,qme,nppmx0,idimp,mx,my,nxe,nye,mx1,mxy1)
aguard2l(qe,nx,ny,nxe,nye)
fwpminmx2(qe,qbme,wpmax,wpmin,nx,ny,nxe,nye)
wpm = 0.5*(wpmax[0] + wpmin[0])*affp
# accelerate convergence: update wpm
if (wpm <= 10.0):
   wpm = 0.75*wpm
print "wpm=",wpm
q2m0 = wpm/affp
# calculate form factor: ffe
isign = 0
mepois23(dcu,cus,isign,ffe,ax,ay,affp,wpm,ci,wf,nx,ny,nxeh,nye,
         nxh,nyh)

# initialize transverse electric field
cus.fill(0.0)

# * * * start main iteration loop * * *

for ntime in xrange(0,nloop):
#  print "ntime = ", ntime

# deposit current with OpenMP: updates cue
   dtimer(dtime,itime,-1)
   cue.fill(0.0)
   gjppost2l(ppart,cue,kpic,qme,zero,nppmx0,idimp,nx,ny,mx,my,nxe,nye,
             mx1,mxy1,ipbc)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tdjpost = tdjpost + time

# deposit charge with OpenMP: updates qe
   dtimer(dtime,itime,-1)
   qe.fill(0.0)
   gppost2l(ppart,qe,kpic,qme,nppmx0,idimp,mx,my,nxe,nye,mx1,mxy1)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tdpost = tdpost + time

# add guard cells with OpenMP: updates qe, cue
   dtimer(dtime,itime,-1)
   aguard2l(qe,nx,ny,nxe,nye)
   acguard2l(cue,nx,ny,nxe,nye)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tguard = tguard + time

# transform charge to fourier space with OpenMP: updates qe
   dtimer(dtime,itime,-1)
   isign = -1
   wfft2rmx(qe,isign,mixup,sct,indx,indy,nxeh,nye,nxhy,nxyh)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tfft = tfft + time

# calculate longitudinal force/charge in fourier space with OpenMP:
# updates fxyze, we
   dtimer(dtime,itime,-1)
   isign = -1
   mpois23(qe,fxyze,isign,ffc,ax,ay,affp,we,nx,ny,nxeh,nye,nxh,nyh)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tfield = tfield + time

# transform longitudinal electric force to real space with OpenMP:
# updates fxyze
   dtimer(dtime,itime,-1)
   isign = 1
   wfft2rm3(fxyze,isign,mixup,sct,indx,indy,nxeh,nye,nxhy,nxyh)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tfft = tfft + time

# transform current to fourier space with OpenMP: update cue
   dtimer(dtime,itime,-1)
   isign = -1
   wfft2rm3(cue,isign,mixup,sct,indx,indy,nxeh,nye,nxhy,nxyh)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tfft = tfft + time

# take transverse part of current with OpenMP: updates cue
   dtimer(dtime,itime,-1)
   mcuperp2(cue,nx,ny,nxeh,nye)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tfield = tfield + time

# calculate magnetic field in fourier space with OpenMP:
# updates bxyze, wm
   dtimer(dtime,itime,-1)
   mbbpois23(cue,bxyze,ffc,ci,wm,nx,ny,nxeh,nye,nxh,nyh)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tfield = tfield + time

# transform magnetic force to real space with OpenMP: updates bxyze
   dtimer(dtime,itime,-1)
   isign = 1
   wfft2rm3(bxyze,isign,mixup,sct,indx,indy,nxeh,nye,nxhy,nxyh)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tfft = tfft + time

# add constant to magnetic field with OpenMP: updates bxyze
   dtimer(dtime,itime,-1)
   baddext2(bxyze,omx,omy,omz,nx,ny,nxe,nye)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tfield = tfield + time

# copy guard cells with OpenMP: updates fxyze, bxyze
   dtimer(dtime,itime,-1)
   bguard2l(fxyze,nx,ny,nxe,nye)
   bguard2l(bxyze,nx,ny,nxe,nye)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tguard = tguard + time

# add longitudinal and old transverse electric fields with OpenMP:
# updates exyze
   dtimer(dtime,itime,-1)
   addvrfield2(exyze,cus,fxyze,ndim,nxe,nye)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tfield = tfield + time

# deposit electron acceleration density and momentum flux with OpenMP:
# updates dcu, amu
   dtimer(dtime,itime,-1)
   dcu.fill(0.0); amu.fill(0.0)
   gdjppost2l(ppart,exyze,bxyze,dcu,amu,kpic,qme,qbme,dt,idimp,nppmx0,
              nx,ny,mx,my,nxe,nye,mx1,mxy1)
# add old scaled electric field with OpenMP: updates dcu
   ascfguard2l(dcu,cus,q2m0,nx,ny,nxe,nye)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tdcjpost = tdcjpost + time

# add guard cells with OpenMP: updates dcu, amu
   dtimer(dtime,itime,-1)
   acguard2l(dcu,nx,ny,nxe,nye)
   amcguard2l(amu,nx,ny,nxe,nye,mdim)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tguard = tguard + time

# transform acceleration density and momentum flux to fourier space
# with OpenMP: updates dcu, amu
   dtimer(dtime,itime,-1)
   isign = -1
   wfft2rm3(dcu,isign,mixup,sct,indx,indy,nxeh,nye,nxhy,nxyh)
   wfft2rmn(amu,ss,isign,mixup,sct,indx,indy,nxeh,nye,mdim,nxhy,nxyh)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tfft = tfft + time

# take transverse part of time derivative of current with OpenMP:
# updates dcu
   dtimer(dtime,itime,-1)
   madcuperp23(dcu,amu,nx,ny,nxeh,nye)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tfield = tfield + time

# calculate transverse electric field with OpenMP: updates cus, wf
   dtimer(dtime,itime,-1)
   isign = -1
   mepois23(dcu,cus,isign,ffe,ax,ay,affp,wpm,ci,wf,nx,ny,nxeh,nye,nxh,
            nyh)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tfield = tfield + time

# transform transverse electric field to real space with OpenMP:
# updates cus
   dtimer(dtime,itime,-1)
   isign = 1
   wfft2rm3(cus,isign,mixup,sct,indx,indy,nxeh,nye,nxhy,nxyh)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tfft = tfft + time

# copy guard cells with OpenMP: updates cus
   dtimer(dtime,itime,-1)
   bguard2l(cus,nx,ny,nxe,nye)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tguard = tguard + time

# add longitudinal and transverse electric fields with OpenMP:
# exyze = cus + fxyze, updates exyze
# cus needs to be retained for next time step
   dtimer(dtime,itime,-1)
   addvrfield2(exyze,cus,fxyze,ndim,nxe,nye)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tfield = tfield + time

# inner iteration loop
   for k in xrange(0,ndc):

# deposit electron current and acceleration density and momentum flux
# with OpenMP: updates cue, dcu, amu
      dtimer(dtime,itime,-1)
      cue.fill(0.0); dcu.fill(0.0); amu.fill(0.0)
      gdcjppost2l(ppart,exyze,bxyze,cue,dcu,amu,kpic,qme,qbme,dt,idimp,
                  nppmx0,nx,ny,mx,my,nxe,nye,mx1,mxy1)
# add scaled electric field with OpenMP: updates dcu
      ascfguard2l(dcu,cus,q2m0,nx,ny,nxe,nye)
      dtimer(dtime,itime,1)
      time = float(dtime)
      tdcjpost = tdcjpost + time

# add guard cells for current, acceleration density, and momentum flux
# with OpenMP: updates cue, dcu, amu
      dtimer(dtime,itime,-1)
      acguard2l(cue,nx,ny,nxe,nye)
      acguard2l(dcu,nx,ny,nxe,nye)
      amcguard2l(amu,nx,ny,nxe,nye,mdim)
      dtimer(dtime,itime,1)
      time = float(dtime)
      tguard = tguard + time

# transform current to fourier space with OpenMP: update cue
      dtimer(dtime,itime,-1)
      isign = -1
      wfft2rm3(cue,isign,mixup,sct,indx,indy,nxeh,nye,nxhy,nxyh)
      dtimer(dtime,itime,1)
      time = float(dtime)
      tfft = tfft + time

# take transverse part of current with OpenMP: updates cue
      dtimer(dtime,itime,-1)
      mcuperp2(cue,nx,ny,nxeh,nye)
      dtimer(dtime,itime,1)
      time = float(dtime)
      tfield = tfield + time

# calculate magnetic field in fourier space with OpenMP:
# updates bxyze, wm
      dtimer(dtime,itime,-1)
      mbbpois23(cue,bxyze,ffc,ci,wm,nx,ny,nxeh,nye,nxh,nyh)
      dtimer(dtime,itime,1)
      time = float(dtime)
      tfield = tfield + time

# transform magnetic force to real space with OpenMP: updates bxyze
      dtimer(dtime,itime,-1)
      isign = 1
      wfft2rm3(bxyze,isign,mixup,sct,indx,indy,nxeh,nye,nxhy,nxyh)
      dtimer(dtime,itime,1)
      time = float(dtime)
      tfft = tfft + time

# add constant to magnetic field with OpenMP: updates bxzye
      dtimer(dtime,itime,-1)
      baddext2(bxyze,omx,omy,omz,nx,ny,nxe,nye)
      dtimer(dtime,itime,1)
      time = float(dtime)
      tfield = tfield + time

# transform acceleration density and momentum flux to fourier space
# with OpenMP: updates dcu and amu
      dtimer(dtime,itime,-1)
      isign = -1
      wfft2rm3(dcu,isign,mixup,sct,indx,indy,nxeh,nye,nxhy,nxyh)
      wfft2rmn(amu,ss,isign,mixup,sct,indx,indy,nxeh,nye,mdim,nxhy,nxyh)
      dtimer(dtime,itime,1)
      time = float(dtime)
      tfft = tfft + time

# take transverse part of time derivative of current with OpenMP:
# updates dcu
      dtimer(dtime,itime,-1)
      madcuperp23(dcu,amu,nx,ny,nxeh,nye)
      dtimer(dtime,itime,1)
      time = float(dtime)
      tfield = tfield + time

# calculate transverse electric field with OpenMP: updates cus, wf
      dtimer(dtime,itime,-1)
      isign = -1
      mepois23(dcu,cus,isign,ffe,ax,ay,affp,wpm,ci,wf,nx,ny,nxeh,nye,
               nxh,nyh)
      dtimer(dtime,itime,1)
      time = float(dtime)
      tfield = tfield + time

# transform transverse electric field to real space with OpenMP:
# updates cus
      dtimer(dtime,itime,-1)
      isign = 1
      wfft2rm3(cus,isign,mixup,sct,indx,indy,nxeh,nye,nxhy,nxyh)
      dtimer(dtime,itime,1)
      time = float(dtime)
      tfft = tfft + time

# copy guard cells with OpenMP: updates bxyze, cus
      dtimer(dtime,itime,-1)
      bguard2l(bxyze,nx,ny,nxe,nye)
      bguard2l(cus,nx,ny,nxe,nye)
      dtimer(dtime,itime,1)
      time = float(dtime)
      tguard = tguard + time

# add longitudinal and transverse electric fields with OpenMP:
# exyze = cus + fxyze, updates exyze
# cus needs to be retained for next time step
      dtimer(dtime,itime,-1)
      addvrfield2(exyze,cus,fxyze,ndim,nxe,nye)
      dtimer(dtime,itime,1)
      time = float(dtime)
      tfield = tfield + time
   pass

# push particles with OpenMP:
   wke[0] = 0.0
   dtimer(dtime,itime,-1)
# updates ppart, wke
#  gbppush23l(ppart,exyze,bxyze,kpic,qbme,dt,dt,wke,idimp,nppmx0,nx,ny,
#             mx,my,nxe,nye,mx1,mxy1,ipbc)
# updates ppart, ncl, ihole, wke, irc
   gbppushf23l(ppart,exyze,bxyze,kpic,ncl,ihole,qbme,dt,dt,wke,idimp,
               nppmx0,nx,ny,mx,my,nxe,nye,mx1,mxy1,ntmax,irc)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tpush = tpush + time
   if (irc[0] != 0):
      print "gbppushf23l error, irc=", irc[0]
      exit(0)

# reorder particles by cell with OpenMP:
   dtimer(dtime,itime,-1)
# updates ppart, ppbuff, kpic, ncl, ihole, and irc
#  pporder2l(ppart,ppbuff,kpic,ncl,ihole,idimp,nppmx0,nx,ny,mx,my,mx1,
#            my1,npbmx,ntmax,irc)
# updates ppart, ppbuff, kpic, ncl, and irc
   pporderf2l(ppart,ppbuff,kpic,ncl,ihole,idimp,nppmx0,mx1,my1,npbmx,
              ntmax,irc)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tsort = tsort + time
   if (irc[0] != 0):
      print "pporderf2l error, ntmax, irc=", ntmax, irc[0]
      exit(0)

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
