#-----------------------------------------------------------------------
# Skeleton 3D Darwin OpenMP PIC code
# written by Viktor K. Decyk, Adam Tableman, and Qiyang Hu, UCLA
import math
import numpy
from fmdpush3 import *
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
# wke/we = particle kinetic/electrostatic field energy
# wf/wm/wt = magnetic field/transverse electric field/total energy
wke = numpy.zeros((1),float_type)
we = numpy.zeros((1),float_type)
wf = numpy.zeros((1),float_type)
wm = numpy.zeros((1),float_type)
wt = numpy.zeros((1),float_type)
zero = 0.0
# mx/my/mz = number of grids in x/y/z in sorting tiles
mx = 8; my = 8; mz = 8
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
# nx/ny/nz = number of grid points in x/y/z direction
np = npx*npy*npz; nx = int(math.pow(2,indx))
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
# mdim = dimension of amu array
mdim = 2*ndim
qbme = qme
affp = float(nx)*float(ny)*float(nz)/float(np)

# allocate data for standard code
# part = particle array
part = numpy.empty((idimp,np),float_type,'F')
# qe = electron charge density with guard cells
qe = numpy.empty((nxe,nye,nze),float_type,'F')
# cue = electron current density with guard cells
cue = numpy.empty((ndim,nxe,nye,nze),float_type,'F')
# dcu = acceleration density with guard cells
dcu = numpy.empty((ndim,nxe,nye,nze),float_type,'F')
# cus = smoothed transverse electric field with guard cells
cus = numpy.empty((ndim,nxe,nye,nze),float_type,'F')
# amu = momentum flux with guard cells
amu = numpy.empty((mdim,nxe,nye,nze),float_type,'F')
# exyze = smoothed total electric field with guard cells
exyze = numpy.empty((ndim,nxe,nye,nze),float_type,'F')
# fxyze = smoothed longitudinal electric field with guard cells
fxyze = numpy.empty((ndim,nxe,nye,nze),float_type,'F')
# bxyze = smoothed magnetic field with guard cells
bxyze = numpy.empty((ndim,nxe,nye,nze),float_type,'F')
# ffc, ffe = form factor arrays for poisson solvers
ffc = numpy.empty((nxh,nyh,nzh),complex_type,'F')
ffe = numpy.empty((nxh,nyh,nzh),complex_type,'F')
# mixup = bit reverse table for FFT
mixup = numpy.empty((nxhyz),int_type,'F')
# sct = sine/cosine table for FFT
sct = numpy.empty((nxyzh),complex_type,'F')
# kpic = number of particles in each tile
kpic = numpy.empty((mxyz1),int_type,'F')
# ss = scratch array for WFFT3RMN
ss = numpy.empty((mdim*nxeh,nze),complex_type,'F')

# prepare fft tables
wfft3rinit(mixup,sct,indx,indy,indz,nxhyz,nxyzh)
# calculate form factor: ffc
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

# find maximum and minimum initial electron density
qe.fill(0.0)
gppost3l(ppart,qe,kpic,qme,nppmx0,idimp,mx,my,mz,nxe,nye,nze,mx1,my1,
         mxyz1)
aguard3l(qe,nx,ny,nz,nxe,nye,nze)
fwpminmx3(qe,qbme,wpmax,wpmin,nx,ny,nz,nxe,nye,nze)
wpm = 0.5*(wpmax[0] + wpmin[0])*affp
# accelerate convergence: update wpm
if (wpm <= 10.0):
   wpm = 0.75*wpm
print "wpm=",wpm
q2m0 = wpm/affp
# calculate form factor: ffe
isign = 0
mepois33(dcu,cus,isign,ffe,ax,ay,az,affp,wpm,ci,wf,nx,ny,nz,nxeh,nye,nze,
         nxh,nyh,nzh)

# initialize transverse electric field
cus.fill(0.0)

# * * * start main iteration loop * * *

for ntime in xrange(0,nloop):
#  print "ntime = ", ntime

# deposit current with OpenMP: updates cue
   dtimer(dtime,itime,-1)
   cue.fill(0.0)
   gjppost3l(ppart,cue,kpic,qme,zero,nppmx0,idimp,nx,ny,nz,mx,my,mz,nxe,
             nye,nze,mx1,my1,mxyz1,ipbc)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tdjpost = tdjpost + time

# deposit charge with OpenMP: updates qe
   dtimer(dtime,itime,-1)
   qe.fill(0.0)
   gppost3l(ppart,qe,kpic,qme,nppmx0,idimp,mx,my,mz,nxe,nye,nze,mx1,my1,
            mxyz1)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tdpost = tdpost + time

# add guard cells with OpenMP: updates qe, cue
   dtimer(dtime,itime,-1)
   aguard3l(qe,nx,ny,nz,nxe,nye,nze)
   acguard3l(cue,nx,ny,nz,nxe,nye,nze)
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

# calculate longitudinal force/charge in fourier space with OpenMP:
# updates fxyze, we
   dtimer(dtime,itime,-1)
   isign = -1
   mpois33(qe,fxyze,isign,ffc,ax,ay,az,affp,we,nx,ny,nz,nxeh,nye,nze,nxh,
           nyh,nzh)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tfield = tfield + time

# transform longitudinal electric force to real space with OpenMP:
# updates fxyze
   dtimer(dtime,itime,-1)
   isign = 1
   wfft3rm3(fxyze,isign,mixup,sct,indx,indy,indz,nxeh,nye,nze,nxhyz,
            nxyzh)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tfft = tfft + time

# transform current to fourier space with OpenMP: update cue
   dtimer(dtime,itime,-1)
   isign = -1
   wfft3rm3(cue,isign,mixup,sct,indx,indy,indz,nxeh,nye,nze,nxhyz,nxyzh)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tfft = tfft + time

# take transverse part of current with OpenMP: updates cue
   dtimer(dtime,itime,-1)
   mcuperp3(cue,nx,ny,nz,nxeh,nye,nze)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tfield = tfield + time

# calculate magnetic field in fourier space with OpenMP:
# updates bxyze, wm
   dtimer(dtime,itime,-1)
   mbbpois33(cue,bxyze,ffc,ci,wm,nx,ny,nz,nxeh,nye,nze,nxh,nyh,nzh)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tfield = tfield + time

# transform magnetic force to real space with OpenMP: updates bxyze
   dtimer(dtime,itime,-1)
   isign = 1
   wfft3rm3(bxyze,isign,mixup,sct,indx,indy,indz,nxeh,nye,nze,nxhyz,
            nxyzh)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tfft = tfft + time

# add constant to magnetic field with OpenMP: updates bxyze
   dtimer(dtime,itime,-1)
   baddext3(bxyze,omx,omy,omz,nx,ny,nz,nxe,nye,nze)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tfield = tfield + time

# copy guard cells with OpenMP: updates fxyze, bxyze
   dtimer(dtime,itime,-1)
   cguard3l(fxyze,nx,ny,nz,nxe,nye,nze)
   cguard3l(bxyze,nx,ny,nz,nxe,nye,nze)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tguard = tguard + time

# add longitudinal and old transverse electric fields with OpenMP:
# updates exyze
   dtimer(dtime,itime,-1)
   addvrfield3(exyze,cus,fxyze,ndim,nxe,nye,nze)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tfield = tfield + time

# deposit electron acceleration density and momentum flux with OpenMP:
# updates dcu, amu
   dtimer(dtime,itime,-1)
   dcu.fill(0.0); amu.fill(0.0)
   gdjppost3l(ppart,exyze,bxyze,kpic,dcu,amu,qme,qbme,dt,idimp,nppmx0,
              nx,ny,nz,mx,my,mz,nxe,nye,nze,mx1,my1,mxyz1)
# add old scaled electric field with OpenMP: updates dcu
   ascfguard3l(dcu,cus,q2m0,nx,ny,nz,nxe,nye,nze)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tdcjpost = tdcjpost + time

# add guard cells with OpenMP: updates dcu, amu
   dtimer(dtime,itime,-1)
   acguard3l(dcu,nx,ny,nz,nxe,nye,nze)
   amcguard3l(amu,nx,ny,nz,nxe,nye,nze,mdim)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tguard = tguard + time

# transform acceleration density and momentum flux to fourier space
# with OpenMP: updates dcu, amu
   dtimer(dtime,itime,-1)
   isign = -1
   wfft3rm3(dcu,isign,mixup,sct,indx,indy,indz,nxeh,nye,nze,nxhyz,nxyzh)
   wfft3rmn(amu,ss,isign,mixup,sct,indx,indy,indz,nxeh,nye,nze,mdim,
            nxhyz,nxyzh)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tfft = tfft + time

# take transverse part of time derivative of current with OpenMP:
# updates dcu
   dtimer(dtime,itime,-1)
   madcuperp3(dcu,amu,nx,ny,nz,nxeh,nye,nze)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tfield = tfield + time

# calculate transverse electric field with OpenMP: updates cus, wf
   dtimer(dtime,itime,-1)
   isign = -1
   mepois33(dcu,cus,isign,ffe,ax,ay,az,affp,wpm,ci,wf,nx,ny,nz,nxeh,nye,
            nze,nxh,nyh,nzh)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tfield = tfield + time

# transform transverse electric field to real space with OpenMP:
# updates cus
   dtimer(dtime,itime,-1)
   isign = 1
   wfft3rm3(cus,isign,mixup,sct,indx,indy,indz,nxeh,nye,nze,nxhyz,nxyzh)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tfft = tfft + time

# copy guard cells with OpenMP: updates cus
   dtimer(dtime,itime,-1)
   cguard3l(cus,nx,ny,nz,nxe,nye,nze)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tguard = tguard + time

# add longitudinal and transverse electric fields with OpenMP:
# exyze = cus + fxyze, updates exyze
# cus needs to be retained for next time step
   dtimer(dtime,itime,-1)
   addvrfield3(exyze,cus,fxyze,ndim,nxe,nye,nze)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tfield = tfield + time

# inner iteration loop
   for k in xrange(0,ndc):

# deposit electron current and acceleration density and momentum flux
# with OpenMP: updates cue, dcu, amu
      dtimer(dtime,itime,-1)
      cue.fill(0.0); dcu.fill(0.0); amu.fill(0.0)
      gdcjppost3l(ppart,exyze,bxyze,kpic,cue,dcu,amu,qme,qbme,dt,idimp,
                  nppmx0,nx,ny,nz,mx,my,mz,nxe,nye,nze,mx1,my1,mxyz1)
# add caled electric field with OpenMP: updates dcu
      ascfguard3l(dcu,cus,q2m0,nx,ny,nz,nxe,nye,nze)
      dtimer(dtime,itime,1)
      time = float(dtime)
      tdcjpost = tdcjpost + time

# add guard cells for current, acceleration density, and momentum flux
# with OpenMP: updates cue, dcu, amu
      dtimer(dtime,itime,-1)
      acguard3l(cue,nx,ny,nz,nxe,nye,nze)
      acguard3l(dcu,nx,ny,nz,nxe,nye,nze)
      amcguard3l(amu,nx,ny,nz,nxe,nye,nze,mdim)
      dtimer(dtime,itime,1)
      time = float(dtime)
      tguard = tguard + time

# transform current to fourier space with OpenMP: update cue
      dtimer(dtime,itime,-1)
      isign = -1
      wfft3rm3(cue,isign,mixup,sct,indx,indy,indz,nxeh,nye,nze,nxhyz,
               nxyzh)
      dtimer(dtime,itime,1)
      time = float(dtime)
      tfft = tfft + time

# take transverse part of current with OpenMP: updates cue
      dtimer(dtime,itime,-1)
      mcuperp3(cue,nx,ny,nz,nxeh,nye,nze)
      dtimer(dtime,itime,1)
      time = float(dtime)
      tfield = tfield + time

# calculate magnetic field in fourier space with OpenMP:
# updates bxyze, wm
      dtimer(dtime,itime,-1)
      mbbpois33(cue,bxyze,ffc,ci,wm,nx,ny,nz,nxeh,nye,nze,nxh,nyh,nzh)
      dtimer(dtime,itime,1)
      time = float(dtime)
      tfield = tfield + time

# transform magnetic force to real space with OpenMP: updates bxyze
      dtimer(dtime,itime,-1)
      isign = 1
      wfft3rm3(bxyze,isign,mixup,sct,indx,indy,indz,nxeh,nye,nze,nxhyz,
               nxyzh)
      dtimer(dtime,itime,1)
      time = float(dtime)
      tfft = tfft + time

# add constant to magnetic field with OpenMP: updates bxyze
      dtimer(dtime,itime,-1)
      baddext3(bxyze,omx,omy,omz,nx,ny,nz,nxe,nye,nze)
      dtimer(dtime,itime,1)
      time = float(dtime)
      tfield = tfield + time

# transform acceleration density and momentum flux to fourier space
# with OpenMP: updates dcu, amu
      dtimer(dtime,itime,-1)
      isign = -1
      wfft3rm3(dcu,isign,mixup,sct,indx,indy,indz,nxeh,nye,nze,nxhyz,
               nxyzh)
      wfft3rmn(amu,ss,isign,mixup,sct,indx,indy,indz,nxeh,nye,nze,mdim,
               nxhyz,nxyzh)
      dtimer(dtime,itime,1)
      time = float(dtime)
      tfft = tfft + time

# take transverse part of time derivative of current with OpenMP:
# updates dcu
      dtimer(dtime,itime,-1)
      madcuperp3(dcu,amu,nx,ny,nz,nxeh,nye,nze)
      dtimer(dtime,itime,1)
      time = float(dtime)
      tfield = tfield + time

# calculate transverse electric field with OpenMP: updates cus, wf
      dtimer(dtime,itime,-1)
      isign = -1
      mepois33(dcu,cus,isign,ffe,ax,ay,az,affp,wpm,ci,wf,nx,ny,nz,nxeh,
               nye,nze,nxh,nyh,nzh)
      dtimer(dtime,itime,1)
      time = float(dtime)
      tfield = tfield + time

# transform transverse electric field to real space with OpenMP:
# updates cus
      dtimer(dtime,itime,-1)
      isign = 1
      wfft3rm3(cus,isign,mixup,sct,indx,indy,indz,nxeh,nye,nze,nxhyz,
               nxyzh)
      dtimer(dtime,itime,1)
      time = float(dtime)
      tfft = tfft + time

# copy guard cells with OpenMP: updates bxyze, cus
      dtimer(dtime,itime,-1)
      cguard3l(bxyze,nx,ny,nz,nxe,nye,nze)
      cguard3l(cus,nx,ny,nz,nxe,nye,nze)
      dtimer(dtime,itime,1)
      time = float(dtime)
      tguard = tguard + time

# add longitudinal and transverse electric fields with OpenMP:
# exyze = cus + fxyze, updates exyze
# cus needs to be retained for next time step
      dtimer(dtime,itime,-1)
      addvrfield3(exyze,cus,fxyze,ndim,nxe,nye,nze)
      dtimer(dtime,itime,1)
      time = float(dtime)
      tfield = tfield + time
   pass

# push particles with OpenMP:
   wke[0] = 0.0
   dtimer(dtime,itime,-1)
# updates ppart, wke
#  gbppush3l(ppart,exyze,bxyze,kpic,qbme,dt,dt,wke,idimp,nppmx0,nx,ny,
#            nz,mx,my,mz,nxe,nye,nze,mx1,my1,mxyz1,ipbc)
# updates ppart, ncl, ihole, wke, irc
   gbppushf3l(ppart,exyze,bxyze,kpic,ncl,ihole,qbme,dt,dt,wke,idimp,
              nppmx0,nx,ny,nz,mx,my,mz,nxe,nye,nze,mx1,my1,mxyz1,ntmax,
              irc)
   dtimer(dtime,itime,1)
   time = float(dtime)
   tpush = tpush + time
   if (irc[0] != 0):
      print "gbppushf3l error, irc=", irc[0]
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
print ""
