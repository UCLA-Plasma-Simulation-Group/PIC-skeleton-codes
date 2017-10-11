/*---------------------------------------------------------------------*/
/* Skeleton 3D Electromagnetic OpenMP PIC code */
/* written by Viktor K. Decyk, UCLA and Ricardo Fonseca, ISCTE */
#include <stdlib.h>
#include <stdio.h>
#include <complex.h>
#include <sys/time.h>
#include "vmbpush3.h"
#include "omplib.h"
#include "avx512lib3.h"
#include "kncmbpush3.h"

void dtimer(double *time, struct timeval *itime, int icntrl);

int main(int argc, char *argv[]) {
/* indx/indy/indz = exponent which determines grid points in x/y/z: */
/* direction: nx = 2**indx, ny = 2**indy, nz = 2**indz */
   int indx =   7, indy =   7, indz =   7;
/* npx/npy/npz = number of electrons distributed in x/y/z direction */
   int npx =  384, npy =   384, npz =   384;
/* ndim = number of velocity coordinates = 3 */
   int ndim = 4;
/* tend = time at end of simulation, in units of plasma frequency */
/* dt = time interval between successive calculations */
/* qme = charge on electron, in units of e */
   float tend = 10.0, dt = 0.035, qme = -1.0;
/* vtx/vty/vtz = thermal velocity of electrons in x/y/z direction */
   float vtx = 1.0, vty = 1.0, vtz = 1.0;
/* vx0/vy0/vz0 = drift velocity of electrons in x/y/z direction */
   float vx0 = 0.0, vy0 = 0.0, vz0 = 0.0;
/* ax/ay/az = smoothed particle size in x/y/z direction */
/* ci = reciprocal of velocity of light */
   float ax = .912871, ay = .912871, az = .912871, ci = 0.1;
/* idimp = number of particle coordinates = 6 */
/* ipbc = particle boundary condition: 1 = periodic */
/* relativity = (no,yes) = (0,1) = relativity is used */
   int idimp = 6, ipbc = 1, relativity = 1;
/* wke/we = particle kinetic/electrostatic field energy             */
/* wf/wm/wt = magnetic field/transverse electric field/total energy */
   float wke = 0.0, we = 0.0, wf = 0.0, wm = 0.0, wt = 0.0;
/* mx/my/mz = number of grids in x/y/z in sorting tiles */
   int mx = 8, my = 8, mz = 8;
/* xtras = fraction of extra particles needed for particle management */
   float xtras = 0.2;
/* kvec = (1,2) = run (autovector,KNC) version */
   int kvec = 1;

/* declare scalars for standard code */
   int j;
   int np, nx, ny, nz, nxh, nyh, nzh, nxe, nye, nze, nxeh;
   int nxyzh, nxhyz, mx1, my1, mz1, mxyz1;
   int ntime, nloop, isign, lvect;
   int irc = 0;
   float qbme, affp, dth;

/* declare scalars for OpenMP code */
   int nppmx, nppmx0, ntmax, npbmx;
   int nvp;

/* declare arrays for standard code: */
/* part = original particle array */
   float *part = NULL;
/* qe = electron charge density with guard cells */
/* cue = electron current density with guard cells */
/* fxyze/bxyze = smoothed electric/magnetic field with guard cells */
   float *qe = NULL, *cue = NULL, *fxyze = NULL, *bxyze = NULL;
/* exyz/bxyz = transverse electric/magnetic field in fourier space */
   float complex *exyz = NULL, *bxyz = NULL;
/* ffc = form factor array for poisson solver */
/* sct = sine/cosine table for FFT */
   float complex *ffc = NULL, *sct = NULL;
/* mixup = bit reverse table for FFT */
   int *mixup = NULL;

/* declare arrays for OpenMP (tiled) code: */
/* ppartt = tiled particle array */
/* ppbuff = buffer array for reordering tiled particle array */
   float *ppartt = NULL, *ppbuff = NULL;
/* kpic = number of particles in each tile */
   int *kpic = NULL;
/* ncl = number of particles departing tile in each direction */
   int *ncl = NULL;
/* ihole = location/destination of each particle departing tile */
   int *ihole = NULL;
/* kp = original location of reordered particle */
   int *kp = NULL;

/* declare and initialize timing data */
   float time;
   struct timeval itime;
   float tdpost = 0.0, tguard = 0.0, tfft = 0.0, tfield = 0.0;
   float tdjpost = 0.0, tpush = 0.0, tsort = 0.0;
   double dtime;

   irc = 0;
/* nvp = number of shared memory nodes  (0=default) */
   nvp = 0;
/* printf("enter number of nodes:\n"); */
/* scanf("%i",&nvp);                   */
/* initialize for shared memory parallel processing */
   cinit_omp(nvp);

/* initialize scalars for standard code */
/* np = total number of particles in simulation */
/* nx/ny/nz = number of grid points in x/y direction */
   np = npx*npy*npz; nx = 1L<<indx; ny = 1L<<indy; nz = 1L<<indz;
   nxh = nx/2; nyh = 1 > ny/2 ? 1 : ny/2; nzh = 1 > nz/2 ? 1 : nz/2;
   nxe = nx + 2; nye = ny + 1; nze = nz + 1; nxeh = nxe/2;
   nxyzh = (nx > ny ? nx : ny); nxyzh = (nxyzh > nz ? nxyzh : nz)/2;
   nxhyz = nxh > ny ? nxh : ny; nxhyz = nxhyz > nz ? nxhyz : nz;
/* mx1/my1/mz1 = number of tiles in x/y/z direction */
   mx1 = (nx - 1)/mx + 1; my1 = (ny - 1)/my + 1;
   mz1 = (nz - 1)/mz + 1; mxyz1 = mx1*my1*mz1;
/* nloop = number of time steps in simulation */
/* ntime = current time step */
   nloop = tend/dt + .0001; ntime = 0;
   qbme = qme;
   affp = ((float) nx)*((float) ny)*((float) nz)/(float ) np;
   dth = 0.0;

/* allocate data for standard code */
   part = (float *) malloc(idimp*np*sizeof(float));
   mixup = (int *) malloc(nxhyz*sizeof(int));
   sct = (float complex *) malloc(nxyzh*sizeof(float complex));
   kpic = (int *) malloc(mxyz1*sizeof(int));

   lvect = 16;
/* allocate vector field data */
   nxe = lvect*((nxe - 1)/lvect + 1);
   nxeh = nxe/2;
   avx512_fallocate(&qe,nxe*nye*nze,&irc);
   avx512_fallocate(&cue,ndim*nxe*nye*nze,&irc);
   avx512_fallocate(&fxyze,ndim*nxe*nye*nze,&irc);
   avx512_fallocate(&bxyze,ndim*nxe*nye*nze,&irc);
   avx512_callocate(&exyz,ndim*nxeh*nye*nze,&irc);
   avx512_callocate(&bxyz,ndim*nxeh*nye*nze,&irc);
   avx512_callocate(&ffc,nxh*nyh*nzh,&irc);
   if (irc != 0) {
      printf("aligned field allocation error: irc = %d\n",irc);
   }

/* prepare fft tables */
   cwfft3rinit(mixup,sct,indx,indy,indz,nxhyz,nxyzh);
/* calculate form factors */
   isign = 0;
   cvmpois33((float complex *)qe,(float complex *)fxyze,isign,ffc,ax,ay,
             az,affp,&we,nx,ny,nz,nxeh,nye,nze,nxh,nyh,nzh);
/* initialize electrons */
   cdistr3(part,vtx,vty,vtz,vx0,vy0,vz0,npx,npy,npz,idimp,np,nx,ny,nz,
           ipbc);

/* find number of particles in each of mx, my, mz tiles: */
/* updates kpic, nppmx */
   cdblkp3l(part,kpic,&nppmx,idimp,np,mx,my,mz,mx1,my1,mxyz1,&irc);
   if (irc != 0) { 
      printf("cdblkp3l error, irc=%d\n",irc);
      exit(1);
   }
/* allocate vector particle data */
   nppmx0 = (1.0 + xtras)*nppmx;
   ntmax = xtras*nppmx;
   npbmx = xtras*nppmx;
/* align data for Vector Processor */
   nppmx0 = lvect*((nppmx0 - 1)/lvect + 1);
   ntmax = lvect*(ntmax/lvect + 1);
   npbmx = lvect*((npbmx - 1)/lvect + 1);
   avx512_fallocate(&ppartt,nppmx0*idimp*mxyz1,&irc);
   avx512_fallocate(&ppbuff,npbmx*idimp*mxyz1,&irc);
   ncl = (int *) malloc(26*mxyz1*sizeof(int));
   ihole = (int *) malloc(2*(ntmax+1)*mxyz1*sizeof(int));
   kp = (int *) malloc(nppmx0*mxyz1*sizeof(int));
   if (irc != 0) {
      printf("aligned particle allocation error: irc = %d\n",irc);
   }

/* copy ordered particle data for OpenMP: updates ppart and kpic */
   cppmovin3ltp(part,ppartt,kpic,kp,nppmx0,idimp,np,mx,my,mz,mx1,my1,
                mxyz1,&irc);
   if (irc != 0) { 
      printf("cppmovin3ltp overflow error, irc=%d\n",irc);
      exit(1);
   }
/* sanity check */
   cppcheck3lt(ppartt,kpic,idimp,nppmx0,nx,ny,nz,mx,my,mz,mx1,my1,mz1,
               &irc);
   if (irc != 0) {
      printf("%d,cppcheck3lt error: irc=%d\n",ntime,irc);
      exit(1);
   }

/* initialize transverse electromagnetic fields */
   cset_cvzero3(exyz,nx,ny,nz,ndim,nxeh,nye,nze);
   cset_cvzero3(bxyz,nx,ny,nz,ndim,nxeh,nye,nze);

   if (dt > 0.37*ci) {
      printf("Warning: Courant condition may be exceeded!\n");
   }

/* * * * start main iteration loop * * * */
 
L500: if (nloop <= ntime)
         goto L2000;
/*    printf("ntime = %i\n",ntime); */
 
/* deposit current with OpenMP: */
      dtimer(&dtime,&itime,-1);
      cset_vzero3(cue,mx,my,mz,ndim,nxe,nye,nze,mx1,my1,mxyz1);
      if (relativity==1) {
         if (kvec==1)
/* updates ppart, cue */
            cvgrjppost3lt(ppartt,cue,kpic,qme,dth,ci,nppmx0,idimp,nx,
                          ny,nz,mx,my,mz,nxe,nye,nze,mx1,my1,mxyz1,
                          ipbc);
/* updates ppart, cue, ncl, ihole, irc */
/*          cvgrjppostf3lt(ppartt,cue,kpic,ncl,ihole,qme,dth,ci,nppmx0, */
/*                      idimp,nx,ny,nz,mx,my,mz,nxe,nye,nze,mx1,my1,    */
/*                      mxyz1,ntmax,&irc);                              */
/* KNC function */
         else if (kvec==2)
/* updates ppart, cue */
/*          ckncgrjppost3lt(ppartt,cue,kpic,qme,dth,ci,nppmx0,idimp, */
/*                          nx,ny,nz,mx,my,mz,nxe,nye,nze,mx1,my1,   */
/*                          mxyz1,ipbc);                             */
            cknc2grjppost3lt(ppartt,cue,kpic,qme,dth,ci,nppmx0,idimp,
                             nx,ny,nz,mx,my,mz,nxe,nye,nze,mx1,my1,
                             mxyz1,ipbc);
/* updates ppart, cue, ncl, ihole, irc */
/*          ckncgrjppostf3lt(ppartt,cue,kpic,ncl,ihole,qme,dth,ci,   */
/*                           nppmx0,idimp,nx,ny,nz,mx,my,mz,nxe,nye, */
/*                           nze,mx1,my1,mxyz1,ntmax,&irc);          */
      }
      else {
         if (kvec==1)
/* updates ppart, cue */
            cvgjppost3lt(ppartt,cue,kpic,qme,dth,nppmx0,idimp,nx,ny,
                         nz,mx,my,mz,nxe,nye,nze,mx1,my1,mxyz1,ipbc);
/* updates ppart, cue, ncl, ihole, irc */
/*          cvgjppostf3lt(ppartt,cue,kpic,ncl,ihole,qme,dth,nppmx0,    */
/*                        idimp,nx,ny,nz,mx,my,mz,nxe,nye,nze,mx1,my1, */
/*                        mxyz1,ntmax,&irc);                           */
/* KNC function */
         else if (kvec==2)
/* updates ppart, cue */
/*          ckncgjppost3lt(ppartt,cue,kpic,qme,dth,nppmx0,idimp,nx,ny,  */
/*                         nz,mx,my,mz,nxe,nye,nze,mx1,my1,mxyz1,ipbc); */
            cknc2gjppost3lt(ppartt,cue,kpic,qme,dth,nppmx0,idimp,nx,ny,
                            nz,mx,my,mz,nxe,nye,nze,mx1,my1,mxyz1,ipbc);
/* updates ppart, cue, ncl, ihole, irc */
/*          ckncgjppostf3lt(ppartt,cue,kpic,ncl,ihole,qme,dth,nppmx0, */
/*                          idimp,nx,ny,nz,mx,my,mz,nxe,nye,nze,mx1,  */
/*                          my1 mxyz1,ntmax,&irc);                    */
      }
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tdjpost += time;
      if (irc != 0) {
         if (relativity==1) {
            printf("cvgjppostf3lt error: irc=%d\n",irc);
         }
         else {
            printf("cvgrjppostf3lt error: irc=%d\n",irc);
         }
         exit(1);
      }

/* reorder particles by tile with OpenMP: */
      dtimer(&dtime,&itime,-1);
      if (kvec==1)
/* updates ppartt, ppbuff, kpic, ncl, ihole, and irc */
         cvpporder3lt(ppartt,ppbuff,kpic,ncl,ihole,idimp,nppmx0,nx,ny,
                      nz,mx,my,mz,mx1,my1,mz1,npbmx,ntmax,&irc);
/* updates ppartt, ppbuff, kpic, ncl, and irc */
/*       cvpporderf3lt(ppartt,ppbuff,kpic,ncl,ihole,idimp,nppmx0,mx1, */
/*                     my1,mz1,npbmx,ntmax,&irc);                     */
/* KNC function */
      else if (kvec==2)
/* updates ppartt, ppbuff, kpic, ncl, ihole, and irc */
         ckncpporder3lt(ppartt,ppbuff,kpic,ncl,ihole,idimp,nppmx0,nx,
                        ny,nz,mx,my,mz,mx1,my1,mz1,npbmx,ntmax,&irc);
/* updates ppartt, ppbuff, kpic, ncl, and irc */
/*       ckncpporderf3lt(ppartt,ppbuff,kpic,ncl,ihole,idimp,nppmx0, */
/*                       mx1,my1,mz1,npbmx,ntmax,&irc);             */
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tsort += time;
      if (irc != 0) {
         printf("cvpporderf3lt error: ntmax, irc=%d,%d\n",ntmax,irc);
         exit(1);
      }

/* deposit charge with OpenMP: updates qe */
      dtimer(&dtime,&itime,-1);
      cset_szero3(qe,mx,my,mz,nxe,nye,nze,mx1,my1,mxyz1);
      if (kvec==1)
         cvgppost3lt(ppartt,qe,kpic,qme,nppmx0,idimp,mx,my,mz,nxe,nye,
                     nze,mx1,my1,mxyz1);
/* KNC function */
      else if (kvec==2)
         cknc2gppost3lt(ppartt,qe,kpic,qme,nppmx0,idimp,mx,my,mz,nxe,
                        nye,nze,mx1,my1,mxyz1);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tdpost += time;

/* add guard cells with OpenMP: updates cue, qe */
      dtimer(&dtime,&itime,-1);
      if (kvec==1) {
         cacguard3l(cue,nx,ny,nz,nxe,nye,nze);
         caguard3l(qe,nx,ny,nz,nxe,nye,nze);
      }
/* KNC function */
      else if (kvec==2) {
         ckncacguard3l(cue,nx,ny,nz,nxe,nye,nze);
         ckncaguard3l(qe,nx,ny,nz,nxe,nye,nze);
      }
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tguard += time;

/* transform charge to fourier space with OpenMP: updates qe */
      dtimer(&dtime,&itime,-1);
      isign = -1;
      if (kvec==1)
         cwfft3rvmx((float complex *)qe,isign,mixup,sct,indx,indy,indz,
                    nxeh,nye,nze,nxhyz,nxyzh);
/* KNC function */
      else if (kvec==2)
         ckncwfft3rmx((float complex *)qe,isign,mixup,sct,indx,indy,
                      indz,nxeh,nye,nze,nxhyz,nxyzh);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tfft += time;

/* transform current to fourier space with OpenMP: update cue */
      dtimer(&dtime,&itime,-1);
      isign = -1;
      if (kvec==1)
         cwfft3rvm3((float complex *)cue,isign,mixup,sct,indx,indy,indz,
                     nxeh,nye,nze,nxhyz,nxyzh);
/* KNC function */
      else if (kvec==2)
         ckncwfft3rm3((float complex *)cue,isign,mixup,sct,indx,indy,
                      indz,nxeh,nye,nze,nxhyz,nxyzh);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tfft += time;

/* take transverse part of current with OpenMP: updates cue */
      dtimer(&dtime,&itime,-1);
      if (kvec==1)
         cmcuperp3((float complex *)cue,nx,ny,nz,nxeh,nye,nze);
/* KNC function */
      else if (kvec==2)
         ckncmcuperp3((float complex *)cue,nx,ny,nz,nxeh,nye,nze);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tfield += time;

/* calculate electromagnetic fields in fourier space with OpenMP: */
/* updates exyz, bxyz, wf, wm                                     */
      dtimer(&dtime,&itime,-1);
      if (ntime==0) {
         if (kvec==1)
            cvmibpois33((float complex *)cue,bxyz,ffc,ci,&wm,nx,ny,nz,
                        nxeh,nye,nze,nxh,nyh,nzh);
/* KNC function */
         else if (kvec==2)
            ckncmibpois33((float complex *)cue,bxyz,ffc,ci,&wm,nx,ny,nz,
                           nxeh,nye,nze,nxh,nyh,nzh);
         wf = 0.0;
         dth = 0.5*dt;
      }
      else {
         if (kvec==1)
            cvmmaxwel3(exyz,bxyz,(float complex *)cue,ffc,ci,dt,&wf,&wm,
                       nx,ny,nz,nxeh,nye,nze,nxh,nyh,nzh);
/* KNC function */
         else if (kvec==2)
            ckncmmaxwel3(exyz,bxyz,(float complex *)cue,ffc,ci,dt,&wf,
                         &wm,nx,ny,nz,nxeh,nye,nze,nxh,nyh,nzh);
      }
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tfield += time;

/* calculate force/charge in fourier space with OpenMP: */
/* updates fxyze, we */
      dtimer(&dtime,&itime,-1);
      isign = -1;
      if (kvec==1)
         cvmpois33((float complex *)qe,(float complex *)fxyze,isign,ffc,
                   ax,ay,az,affp,&we,nx,ny,nz,nxeh,nye,nze,nxh,nyh,nzh);
/* KNC function */
      else if (kvec==2)
         ckncmpois33((float complex *)qe,(float complex *)fxyze,isign,
                     ffc,ax,ay,az,affp,&we,nx,ny,nz,nxeh,nye,nze,nxh,
                     nyh,nzh);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tfield += time;

/* add longitudinal and transverse electric fields with OpenMP: */
/* updates fxyze */
      dtimer(&dtime,&itime,-1);
      isign = 1;
      if (kvec==1)
         cvmemfield3((float complex *)fxyze,exyz,ffc,isign,nx,ny,nz,
                     nxeh,nye,nze,nxh,nyh,nzh);
/* KNC function */
      else if (kvec==2)
         ckncmemfield3((float complex *)fxyze,exyz,ffc,isign,nx,ny,nz,
                        nxeh,nye,nze,nxh,nyh,nzh);
/* copy magnetic field with OpenMP: updates bxyze */
      isign = -1;
      if (kvec==1)
         cvmemfield3((float complex *)bxyze,bxyz,ffc,isign,nx,ny,nz,
                     nxeh,nye,nze,nxh,nyh,nzh);
/* KNC function */
      else if (kvec==2)
         ckncmemfield3((float complex *)bxyze,bxyz,ffc,isign,nx,ny,nz,
                        nxeh,nye,nze,nxh,nyh,nzh);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tfield += time;

/* transform force to real space with OpenMP: updates fxyze */
      dtimer(&dtime,&itime,-1);
      isign = 1;
      if (kvec==1)
         cwfft3rvm3((float complex *)fxyze,isign,mixup,sct,indx,indy,
                    indz,nxeh,nye,nze,nxhyz,nxyzh);
/* KNC function */
      else if (kvec==2)
         ckncwfft3rm3((float complex *)fxyze,isign,mixup,sct,indx,indy,
                       indz,nxeh,nye,nze,nxhyz,nxyzh);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tfft += time;

/* transform magnetic force to real space with OpenMP: updates bxyze */
      dtimer(&dtime,&itime,-1);
      isign = 1;
      if (kvec==1)
         cwfft3rvm3((float complex *)bxyze,isign,mixup,sct,indx,indy,
                    indz,nxeh,nye,nze,nxhyz,nxyzh);
/* KNC function */
      else if (kvec==2)
         ckncwfft3rm3((float complex *)bxyze,isign,mixup,sct,indx,indy,
                      indz,nxeh,nye,nze,nxhyz,nxyzh);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tfft += time;

/* copy guard cells with OpenMP: updates fxyze, bxyze */
      dtimer(&dtime,&itime,-1);
      if (kvec==1) {
         ccguard3l(fxyze,nx,ny,nz,nxe,nye,nze);
         ccguard3l(bxyze,nx,ny,nz,nxe,nye,nze);
      }
/* KNC function */
      else if (kvec==2) {
         cknccguard3l(fxyze,nx,ny,nz,nxe,nye,nze);
         cknccguard3l(bxyze,nx,ny,nz,nxe,nye,nze);
      }
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tguard += time;

/* push particles with OpenMP: */
      wke = 0.0;
      dtimer(&dtime,&itime,-1);
      if (relativity==1) {
         if (kvec==1)
/* updates ppart, wke */
/*          cvgrbppush3lt(ppartt,fxyze,bxyze,kpic,qbme,dt,dth,ci,&wke, */
/*                     idimp,nppmx0,nx,ny,nz,mx,my,mz,nxe,nye,nze,mx1, */
/*                     my1, mxyz1,ipbc);                               */
            cv2grbppush3lt(ppartt,fxyze,bxyze,kpic,qbme,dt,dth,ci,&wke,
                        idimp,nppmx0,nx,ny,nz,mx,my,mz,nxe,nye,nze,mx1,
                        my1, mxyz1,ipbc);
/* updates ppart, ncl, ihole, wke, irc */
/*          cvgrbppushf3lt(ppartt,fxyze,bxyze,kpic,ncl,ihole,qbme,dt, */
/*                         dth, ci,&wke,idimp,nppmx0,nx,ny,nz,mx,my,  */
/*                         mz,nxe,nye,nze,mx1,my1,mxyz1,ntmax,&irc);  */
/*          cv2grbppushf3lt(ppartt,fxyze,bxyze,kpic,ncl,ihole,qbme,dt, */
/*                          dth, ci,&wke,idimp,nppmx0,nx,ny,nz,mx,my,  */
/*                          mz,nxe,nye,nze,mx1,my1,mxyz1,ntmax,&irc);  */
/* KNC function */
         else if (kvec==2)
/* updates ppart, wke */
            ckncgrbppush3lt(ppartt,fxyze,bxyze,kpic,qbme,dt,dth,ci,
                            &wke,idimp,nppmx0,nx,ny,nz,mx,my,mz,nxe,
                            nye,nze,mx1,my1,mxyz1,ipbc);
/* updates ppart, ncl, ihole, wke, irc */
/*          ckncgrbppushf3lt(ppartt,fxyze,bxyze,kpic,ncl,ihole,qbme,dt, */
/*                           dth,ci,&wke,idimp,nppmx0,nx,ny,nz,mx,my,   */
/*                           mz,nxe,nye,nze,mx1,my1,mxyz1,ntmax,&irc);  */
      }
      else {
         if (kvec==1)
/* updates ppart, wke */
            cvgbppush3lt(ppartt,fxyze,bxyze,kpic,qbme,dt,dth,&wke,
                         idimp,nppmx0,nx,ny,nz,mx,my,mz,nxe,nye,nze,
                         mx1,my1,mxyz1,ipbc);
/*          cv2gbppush3lt(ppartt,fxyze,bxyze,kpic,qbme,dt,dth,&wke,   */
/*                        idimp,nppmx0,nx,ny,nz,mx,my,mz,nxe,nye,nze, */
/*                        mx1,my1,mxyz1,ipbc);                        */
/* updates ppart, ncl, ihole, wke, irc */
/*          cvgbppushf3lt(ppartt,fxyze,bxyze,kpic,ncl,ihole,qbme,dt,   */
/*                        dth,&wke,idimp,nppmx0,nx,ny,nz,mx,my,mz,nxe, */
/*                        nye,nze,mx1,my1,mxyz1,ntmax,&irc);           */
/*          cv2gbppushf3lt(ppartt,fxyze,bxyze,kpic,ncl,ihole,qbme,dt,   */
/*                         dth,&wke,idimp,nppmx0,nx,ny,nz,mx,my,mz,nxe, */
/*                         nye,nze,mx1,my1,mxyz1,ntmax,&irc);           */
/* KNC function */
         else if (kvec==2)
/* updates ppart, wke */
            ckncgbppush3lt(ppartt,fxyze,bxyze,kpic,qbme,dt,dth,&wke,
                           idimp,nppmx0,nx,ny,nz,mx,my,mz,nxe,nye,nze,
                           mx1,my1,mxyz1,ipbc);
/* updates ppart, ncl, ihole, wke, irc */
/*          ckncgbppushf3lt(ppartt,fxyze,bxyze,kpic,ncl,ihole,qbme,dt,  */
/*                          dt,dth,&wke,idimp,nppmx0,nx,ny,nz,mx,my,mz, */
/*                          nxe,nye,nze,mx1,my1,mxyz1,ntmax,&irc);      */
      }
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tpush += time;
      if (irc != 0) {
         if (relativity==1) {
            printf("cgrbppushf3l error: irc=%d\n",irc);
         }
         else {
            printf("cgbppushf3l error: irc=%d\n",irc);
         }
         exit(1);
      }

/* reorder particles by tile with OpenMP: */
      dtimer(&dtime,&itime,-1);
      if (kvec==1)
/* updates ppartt, ppbuff, kpic, ncl, ihole, and irc */
         cvpporder3lt(ppartt,ppbuff,kpic,ncl,ihole,idimp,nppmx0,nx,ny,
                      nz,mx,my,mz,mx1,my1,mz1,npbmx,ntmax,&irc);
/* updates ppartt, ppbuff, kpic, ncl, and irc */
/*       cvpporderf3lt(ppartt,ppbuff,kpic,ncl,ihole,idimp,nppmx0,mx1, */
/*                     my1,mz1,npbmx,ntmax,&irc);                     */
/* KNC function */
      else if (kvec==2)
/* updates ppartt, ppbuff, kpic, ncl, ihole, and irc */
         ckncpporder3lt(ppartt,ppbuff,kpic,ncl,ihole,idimp,nppmx0,nx,
                        ny,nz,mx,my,mz,mx1,my1,mz1,npbmx,ntmax,&irc);
/* updates ppartt, ppbuff, kpic, ncl, and irc */
/*       ckncpporderf3lt(ppartt,ppbuff,kpic,ncl,ihole,idimp,nppmx0, */
/*                       mx1,my1,mz1,npbmx,ntmax,&irc);             */
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tsort += time;
      if (irc != 0) {
         printf("cvpporderf3lt error: ntmax, irc=%d,%d\n",ntmax,irc);
         exit(1);
      }

      if (ntime==0) {
         wt = we + wf + wm;
         printf("Initial Total Field, Kinetic and Total Energies:\n");
         printf("%e %e %e\n",wt,wke,wke+wt);
         printf("Initial Electrostatic, Transverse Electric and Magnetic \
Field Energies:\n");
         printf("%e %e %e\n",we,wf,wm);
      }
      ntime += 1;
      goto L500;
L2000:

/* * * * end main iteration loop * * * */

   printf("ntime, relativity = %i,%i\n",ntime,relativity);
   printf("kvec = %i\n",kvec);
   wt = we + wf + wm;
   printf("Final Total Field, Kinetic and Total Energies:\n");
   printf("%e %e %e\n",wt,wke,wke+wt);
   printf("Final Electrostatic, Transverse Electric and Magnetic Field \
Energies:\n");
   printf("%e %e %e\n",we,wf,wm);

   printf("\n");
   printf("deposit time = %f\n",tdpost);
   printf("current deposit time = %f\n",tdjpost);
   tdpost += tdjpost;
   printf("total deposit time = %f\n",tdpost);
   printf("guard time = %f\n",tguard);
   printf("solver time = %f\n",tfield);
   printf("fft time = %f\n",tfft);
   printf("push time = %f\n",tpush);
   printf("sort time = %f\n",tsort);
   tfield += tguard + tfft;
   printf("total solver time = %f\n",tfield);
   time = tdpost + tpush + tsort;
   printf("total particle time = %f\n",time);
   wt = time + tfield;
   printf("total time = %f\n",wt);
   printf("\n");

   wt = 1.0e+09/(((float) nloop)*((float) np));
   printf("Push Time (nsec) = %f\n",tpush*wt);
   printf("Deposit Time (nsec) = %f\n",tdpost*wt);
   printf("Sort Time (nsec) = %f\n",tsort*wt);
   printf("Total Particle Time (nsec) = %f\n",time*wt);
   printf("\n");

   avx512_deallocate(ppartt);
   avx512_deallocate(ppbuff);
   avx512_deallocate(ffc);
   avx512_deallocate(bxyz);
   avx512_deallocate(exyz);
   avx512_deallocate(bxyze);
   avx512_deallocate(fxyze);
   avx512_deallocate(cue);
   avx512_deallocate(qe);

   return 0;
}
