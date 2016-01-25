/*--------------------------------------------------------------------*/
/* Skeleton 2-1/2D Electromagnetic MPI/OpenMP PIC code */
/* written by Viktor K. Decyk, UCLA */
#include <stdlib.h>
#include <stdio.h>
#include <complex.h>
#include <sys/time.h>
#include "mpbpush2.h"
#include "mpplib2.h"
#include "omplib.h"

void dtimer(double *time, struct timeval *itime, int icntrl);

int main(int argc, char *argv[]) {
/* indx/indy = exponent which determines grid points in x/y direction: */
/* nx = 2**indx, ny = 2**indy */
   int indx =   9, indy =   9;
/* npx/npy = number of electrons distributed in x/y direction */
   int npx =  3072, npy =   3072;
/* ndim = number of velocity coordinates = 3 */
   int ndim = 3;
/* tend = time at end of simulation, in units of plasma frequency */
/* dt = time interval between successive calculations */
/* qme = charge on electron, in units of e */
   float tend = 10.0, dt = 0.04, qme = -1.0;
/* float tend = 10.0, dt = 0.025, qme = -1.0;        */
/* vtx/vty = thermal velocity of electrons in x/y direction */
/* vx0/vy0 = drift velocity of electrons in x/y direction */
   float vtx = 1.0, vty = 1.0, vx0 = 0.0, vy0 = 0.0;
/* vtx/vz0 = thermal/drift velocity of electrons in z direction */
   float vtz = 1.0, vz0 = 0.0;
/* ax/ay = smoothed particle size in x/y direction */
/* ci = reciprocal of velocity of light */
   float ax = .912871, ay = .912871, ci = 0.1;
/* idimp = dimension of phase space = 5 */
/* ipbc = particle boundary condition: 1 = periodic */
/* relativity = (no,yes) = (0,1) = relativity is used */
   int idimp = 5, ipbc = 1, relativity = 1;
/* idps = number of partition boundaries */
   int idps = 2;
/* wke/we = particle kinetic/electrostatic field energy */
/* wf/wm/wt = magnetic field/transverse electric field/total energy */
   float wke = 0.0, we = 0.0, wf = 0.0, wm = 0.0, wt = 0.0;
/* sorting tiles, should be less than or equal to 32 */
   int mx = 16, my = 16;
/* fraction of extra particles needed for particle management */
   float xtras = 0.2;
/* declare scalars for standard code */
   int j;
   int nx, ny, nxh, nyh, nxe, nye, nxeh, nnxe, nxyh, nxhy;
   int mx1, ntime, nloop, isign, ierr;
   float qbme, affp, dth;
   double np;

/* declare scalars for MPI code */
   int ntpose = 1;
   int nvp, idproc, kstrt, npmax, kxp, kyp, nypmx, nypmn;
   int nyp, noff, npp, nps, myp1, mxyp1;

/* declare scalars for OpenMP code */
   int nppmx, nppmx0, nbmaxp, ntmaxp, npbmx, irc;
   int nvpp;

/* declare arrays for standard code */
/* part = particle array */
   float *part = NULL;
/* qe = electron charge density with guard cells */
   float *qe = NULL;
/* cue = electron current density with guard cells */
/* fxyze/bxyze = smoothed electric/magnetic field with guard cells */
   float *cue = NULL, *fxyze = NULL, *bxyze = NULL;
/* exyz/bxyz = transverse electric/magnetic field in fourier space */
   float complex *exyz = NULL, *bxyz = NULL;
/* qt = scalar charge density field array in fourier space */
   float complex *qt = NULL;
/* cut = vector current density field array in fourier space */
/* fxyt/bxyt = vector electric/magnetic field in fourier space */
   float complex *cut = NULL, *fxyt = NULL, *bxyt = NULL;
/* ffc = form factor array for poisson solver */
   float complex *ffc = NULL;
/* mixup = bit reverse table for FFT */
   int *mixup = NULL;
/* sct = sine/cosine table for FFT */
   float complex *sct = NULL;
   double wtot[7], work[7];

/* declare arrays for MPI code */
/* bs/br = complex send/receive buffers for data transpose */
   float complex *bs = NULL, *br = NULL;
/* sbufl/sbufr = particle buffers sent to nearby processors */
/* rbufl/rbufr = particle buffers received from nearby processors */
   float *sbufl = NULL, *sbufr = NULL, *rbufl = NULL, *rbufr = NULL;
/* edges[0:1] = lower:upper y boundaries of particle partition */
   float *edges = NULL;
/* scr = guard cell buffer received from nearby processors */
   float *scr = NULL;

/* declare arrays for OpenMP code */
/* ppart = tiled particle array */
/* ppbuff = buffer array for reordering tiled particle array */
   float *ppart = NULL, *ppbuff = NULL;
/* kpic = number of particles in each tile */
   int *kpic = NULL;
/* ncl = number of particles departing tile in each direction */
/* iholep = location/destination of each particle departing tile */
   int *ncl = NULL, *iholep = NULL;
/* ncll/nclr/mcll/mclr = number offsets send/received from processors */
   int *ncll = NULL, *nclr = NULL, *mcll = NULL, *mclr = NULL;

/* declare and initialize timing data */
   float time;
   struct timeval itime;
   float tdpost = 0.0, tguard = 0.0, ttp = 0.0, tfield = 0.0;
   float tdjpost = 0.0, tpush = 0.0, tsort = 0.0, tmov = 0.0;
   float tfft[2] = {0.0,0.0};
   double dtime;

   irc = 0;
/* nvpp = number of shared memory nodes  (0=default) */
   nvpp = 0;
/* printf("enter number of nodes:\n"); */
/* scanf("%i",&nvpp);                   */
/* initialize for shared memory parallel processing */
   cinit_omp(nvpp);

/* initialize scalars for standard code */
/* np = total number of particles in simulation */
   np =  (double) npx*(double) npy;
/* nx/ny = number of grid points in x/y direction */
   nx = 1L<<indx; ny = 1L<<indy;
   nxh = nx/2; nyh = 1 > ny/2 ? 1 : ny/2;
   nxe = nx + 2; nye = ny + 2; nxeh = nxe/2; nnxe = ndim*nxe;
   nxyh = (nx > ny ? nx : ny)/2; nxhy = nxh > ny ? nxh : ny;
/* mx1 = number of tiles in x direction */
   mx1 = (nx - 1)/mx + 1;
/* nloop = number of time steps in simulation */
/* ntime = current time step */
   nloop = tend/dt + .0001; ntime = 0;
   qbme = qme;
   affp = (double) nx*(double) ny/np;
   dth = 0.0;

/* nvp = number of distributed memory nodes */
/* initialize for distributed memory parallel processing */
   cppinit2(&idproc,&nvp,argc,argv);
   kstrt = idproc + 1;
/* check if too many processors */
   if (nvp > ny) {
      if (kstrt==1) {
         printf("Too many processors requested: ny, nvp=%d,%d\n",ny,nvp);
      }
      goto L3000;
   }

/* initialize data for MPI code */
   edges = (float *) malloc(idps*sizeof(float));
/* calculate partition variables: edges, nyp, noff, nypmx              */
/* edges[0:1] = lower:upper boundary of particle partition             */
/* nyp = number of primary (complete) gridpoints in particle partition */
/* noff = lowermost global gridpoint in particle partition             */
/* nypmx = maximum size of particle partition, including guard cells   */
/* nypmn = minimum value of nyp                                        */
   cpdicomp2l(edges,&nyp,&noff,&nypmx,&nypmn,ny,kstrt,nvp,idps);
   if (nypmn < 1) {
      if (kstrt==1) {
         printf("combination not supported nvp, ny = %d,%d\n",nvp,ny);
      }
      goto L3000;
   }

/* initialize additional scalars for MPI code */
/* kxp = number of complex grids in each field partition in x direction */
   kxp = (nxh - 1)/nvp + 1;
/* kyp = number of complex grids in each field partition in y direction */
   kyp = (ny - 1)/nvp + 1;
/* npmax = maximum number of electrons in each partition */
   npmax = (np/nvp)*1.25;
/* myp1 = number of tiles in y direction */
   myp1 = (nyp - 1)/my + 1; mxyp1 = mx1*myp1;

/* allocate and initialize data for standard code */
   part = (float *) malloc(idimp*npmax*sizeof(float));
   qe = (float *) malloc(nxe*nypmx*sizeof(float));
   fxyze = (float *) malloc(ndim*nxe*nypmx*sizeof(float));
   cue = (float *) malloc(ndim*nxe*nypmx*sizeof(float));
   bxyze = (float *) malloc(ndim*nxe*nypmx*sizeof(float));
   exyz = (float complex *) malloc(ndim*nye*kxp*sizeof(float complex));
   bxyz = (float complex *) malloc(ndim*nye*kxp*sizeof(float complex));
   qt = (float complex *) malloc(nye*kxp*sizeof(float complex));
   fxyt = (float complex *) malloc(ndim*nye*kxp*sizeof(float complex));
   cut = (float complex *) malloc(ndim*nye*kxp*sizeof(float complex));
   bxyt = (float complex *) malloc(ndim*nye*kxp*sizeof(float complex));
   ffc = (float complex *) malloc(nyh*kxp*sizeof(float complex));
   mixup = (int *) malloc(nxhy*sizeof(int));
   sct = (float complex *) malloc(nxyh*sizeof(float complex));
   kpic = (int *) malloc(mxyp1*sizeof(int));

/* allocate and initialize data for MPI code */
   bs = (float complex *) malloc(ndim*kxp*kyp*sizeof(float complex));
   br = (float complex *) malloc(ndim*kxp*kyp*sizeof(float complex));
   scr = (float *) malloc(ndim*nxe*sizeof(float));

/* prepare fft tables */
   cwpfft2rinit(mixup,sct,indx,indy,nxhy,nxyh);
/* calculate form factors */
   isign = 0;
   cmppois23(qt,fxyt,isign,ffc,ax,ay,affp,&we,nx,ny,kstrt,nye,kxp,nyh);
/* initialize electrons */
   nps = 1;
   npp = 0;
   cpdistr2h(part,edges,&npp,nps,vtx,vty,vtz,vx0,vy0,vz0,npx,npy,nx,ny,
             idimp,npmax,idps,ipbc,&ierr);
/* check for particle initialization error */
   if (ierr != 0) {
      if (kstrt==1) {
         printf("particle initialization error: ierr=%d\n",ierr);
      }
      goto L3000;
   }

/* initialize transverse electromagnetic fields */
   for (j = 0; j < ndim*nye*kxp; j++) {
      exyz[j] = 0.0 + 0.0*_Complex_I;
      bxyz[j] = 0.0 + 0.0*_Complex_I;
   }

/* find number of particles in each of mx, my tiles: updates kpic, nppmx */
   cppdblkp2l(part,kpic,npp,noff,&nppmx,idimp,npmax,mx,my,mx1,mxyp1,
              &irc);
   if (irc != 0) {
      printf("%d,cppdblkp2l error, irc=%d\n",kstrt,irc);
      cppabort();
      exit(1);
   }
/* allocate vector particle data */
   nppmx0 = (1.0 + xtras)*nppmx;
   ntmaxp = xtras*nppmx;
   npbmx = xtras*nppmx;
   nbmaxp = 0.25*mx1*npbmx;
   sbufl = (float *) malloc(idimp*nbmaxp*sizeof(float));
   sbufr = (float *) malloc(idimp*nbmaxp*sizeof(float));
   rbufl = (float *) malloc(idimp*nbmaxp*sizeof(float));
   rbufr = (float *) malloc(idimp*nbmaxp*sizeof(float));
   ppart = (float *) malloc(idimp*nppmx0*mxyp1*sizeof(float));
   ppbuff = (float *) malloc(idimp*npbmx*mxyp1*sizeof(float));
   ncl = (int *) malloc(8*mxyp1*sizeof(int));
   iholep = (int *) malloc(2*(ntmaxp+1)*mxyp1*sizeof(int));
   ncll = (int *) malloc(3*mxyp1*sizeof(int));
   nclr = (int *) malloc(3*mxyp1*sizeof(int));
   mcll = (int *) malloc(3*mxyp1*sizeof(int));
   mclr = (int *) malloc(3*mxyp1*sizeof(int));

/* copy ordered particle data for OpenMP */
   cpppmovin2l(part,ppart,kpic,npp,noff,nppmx0,idimp,npmax,mx,my,mx1,
                  mxyp1,&irc);
   if (irc != 0) { 
      printf("%d,cpppmovin2l overflow error, irc=%d\n",kstrt,irc);
      cppabort();
      exit(1);
   }
/* sanity check */
   cpppcheck2l(ppart,kpic,noff,nyp,idimp,nppmx0,nx,mx,my,mx1,myp1,&irc);
   if (irc != 0) { 
      printf("%d,cpppcheck2l error: irc=%d\n",kstrt,irc);
      cppabort();
      exit(1);
   }

   if (dt > 0.45*ci) {
      if (kstrt==1) {
         printf("Warning: Courant condition may be exceeded!\n");
      }
   }

/* * * * start main iteration loop * * * */

L500: if (nloop <= ntime)
         goto L2000;
/*    if (kstrt==1) printf("ntime = %i\n",ntime); */

/* deposit current with OpenMP: */
      dtimer(&dtime,&itime,-1);
      for (j = 0; j < ndim*nxe*nypmx; j++) {
         cue[j] = 0.0;
      }
      if (relativity==1) {
/* updates ppart, cue */
/*       cppgrjppost2l(ppart,cue,kpic,noff,qme,dth,ci,nppmx0,idimp,nx, */
/*                     ny,mx,my,nxe,nypmx,mx1,mxyp1,ipbc);             */
/* updates ppart, cue, ncl, iholep, irc */
         cppgrjppostf2l(ppart,cue,kpic,ncl,iholep,noff,nyp,qme,dth,
                        ci,nppmx0,idimp,nx,ny,mx,my,nxe,nypmx,mx1,
                        mxyp1,ntmaxp,&irc);
      }
      else {
/* updates ppart, cue */
/*       cppgjppost2l(ppart,cue,kpic,noff,qme,dth,nppmx0,idimp,nx,ny, */
/*                    mx,my,nxe,nypmx,mx1,mxyp1,ipbc);                */
/* updates ppart, cue, ncl, iholep, irc */
         cppgjppostf2l(ppart,cue,kpic,ncl,iholep,noff,nyp,qme,dth,
                       nppmx0,idimp,nx,ny,mx,my,nxe,nypmx,mx1,mxyp1,
                       ntmaxp,&irc);
      }
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tdjpost += time;
      if (irc != 0) {
         if (relativity==1) {
            printf("%d,cppgrjppostf2l error: irc=%d\n",kstrt,irc);
         }
         else {
            printf("%d,cppgjppostf2l error: irc=%d\n",kstrt,irc);
         }
         cppabort();
         exit(1);
      }

/* reorder particles by tile with OpenMP */
/* first part of particle reorder on x and y cell with mx, my tiles: */
      dtimer(&dtime,&itime,-1);
/* updates ppart, ppbuff, sbufl, sbufr, ncl, iholep, ncll, nclr, irc */
/*    cppporder2la(ppart,ppbuff,sbufl,sbufr,kpic,ncl,iholep,ncll,nclr, */
/*                 noff,nyp,idimp,nppmx0,nx,ny,mx,my,mx1,myp1,npbmx,   */
/*                 ntmaxp,nbmaxp,&irc);                                */
/* updates ppart, ppbuff, sbufl, sbufr, ncl, ncll, nclr, irc */
      cppporderf2la(ppart,ppbuff,sbufl,sbufr,ncl,iholep,ncll,nclr,
                    idimp,nppmx0,mx1,myp1,npbmx,ntmaxp,nbmaxp,&irc);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tsort += time;
      if (irc != 0) {
         printf("%d,cppporderf2la error:ntmaxp,irc=%d,%d\n",kstrt,
                ntmaxp,irc);
         cppabort();
         exit(1);
      }
/* move particles into appropriate spatial regions: */
/* updates rbufr, rbufl, mcll, mclr */
      dtimer(&dtime,&itime,-1);
      cpppmove2(sbufr,sbufl,rbufr,rbufl,ncll,nclr,mcll,mclr,kstrt,nvp,
                idimp,nbmaxp,mx1);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tmov += time;
/* second part of particle reorder on x and y cell with mx, my tiles: */
/* updates ppart, kpic */
      dtimer(&dtime,&itime,-1);
      cppporder2lb(ppart,ppbuff,rbufl,rbufr,kpic,ncl,iholep,mcll,mclr,
                   idimp,nppmx0,mx1,myp1,npbmx,ntmaxp,nbmaxp,&irc);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tsort += time;
      if (irc != 0) {
         printf("%d,cppporder2lb error:nppmx0,irc=%d,%d\n",kstrt,nppmx0,
                irc);
         cppabort();
         exit(1);
      }

/* deposit charge with OpenMP: updates qe */
      dtimer(&dtime,&itime,-1);
      for (j = 0; j < nxe*nypmx; j++) {
         qe[j] = 0.0;
      }
      cppgppost2l(ppart,qe,kpic,noff,qme,idimp,nppmx0,mx,my,nxe,nypmx,
                  mx1,mxyp1);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tdpost += time;

/* add guard cells with OpenMP: updates cue, qe */
      dtimer(&dtime,&itime,-1);
      cppacguard2xl(cue,nyp,nx,ndim,nxe,nypmx);
      cppnacguard2l(cue,scr,nyp,nx,ndim,kstrt,nvp,nxe,nypmx);
      cppaguard2xl(qe,nyp,nx,nxe,nypmx);
      cppnaguard2l(qe,scr,nyp,nx,kstrt,nvp,nxe,nypmx);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tguard += time;

/* transform charge to fourier space with OpenMP: updates qt */
/* modifies qe */
      dtimer(&dtime,&itime,-1);
      isign = -1;
      cwppfft2rm((float complex *)qe,qt,bs,br,isign,ntpose,mixup,sct,&ttp,
                 indx,indy,kstrt,nvp,nxeh,nye,kxp,kyp,nypmx,nxhy,nxyh);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tfft[0] += time;
      tfft[1] += ttp;

/* transform current to fourier space with OpenMP: updates cut */
/* modifies cue */
      dtimer(&dtime,&itime,-1);
      isign = -1;
      cwppfft2rm3((float complex *)cue,cut,bs,br,isign,ntpose,mixup,sct,
                  &ttp,indx,indy,kstrt,nvp,nxeh,nye,kxp,kyp,nypmx,nxhy,
                  nxyh);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tfft[0] += time;
      tfft[1] += ttp;

/* take transverse part of current with OpenMP: updates cut */
/* modifies cue */
      dtimer(&dtime,&itime,-1);
      cmppcuperp2(cut,nx,ny,kstrt,nye,kxp);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tfield += time;

/* calculate electromagnetic fields in fourier space with OpenMP: */
/* updates exyz, bxyz, wf, wm */
      dtimer(&dtime,&itime,-1);
      if (ntime==0) {
         cmippbpoisp23(cut,bxyz,ffc,ci,&wm,nx,ny,kstrt,nye,kxp,nyh);
         wf = 0.0;
         dth = 0.5*dt;
      }
      else {
         cmppmaxwel2(exyz,bxyz,cut,ffc,affp,ci,dt,&wf,&wm,nx,ny,kstrt,
                     nye,kxp,nyh);
      }
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tfield += time;

/* calculate force/charge in fourier space with OpenMP: */
/* updates fxyt, we */
      dtimer(&dtime,&itime,-1);
      isign = -1;
      cmppois23(qt,fxyt,isign,ffc,ax,ay,affp,&we,nx,ny,kstrt,nye,kxp,
                nyh);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tfield += time;

/* add longitudinal and transverse electric fields with OpenMP: */
/* updates fxyt */
      dtimer(&dtime,&itime,-1);
      isign = 1;
      cmppemfield2(fxyt,exyz,ffc,isign,nx,ny,kstrt,nye,kxp,nyh);
/* copy magnetic field with OpenMP: updates bxyt */
      isign = -1;
      cmppemfield2(bxyt,bxyz,ffc,isign,nx,ny,kstrt,nye,kxp,nyh);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tfield += time;

/* transform force to real space with OpenMP: updates fxyze */
/* modifies fxyt */
      dtimer(&dtime,&itime,-1);
      isign = 1;
      cwppfft2rm3((float complex *)fxyze,fxyt,bs,br,isign,ntpose,mixup,
                  sct,&ttp,indx,indy,kstrt,nvp,nxeh,nye,kxp,kyp,nypmx,
                  nxhy,nxyh);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tfft[0] += time;
      tfft[1] += ttp;

/* transform magnetic field to real space with OpenMP: */
/* updates bxyze, modifies bxyt */
      dtimer(&dtime,&itime,-1);
      isign = 1;
      cwppfft2rm3((float complex *)bxyze,bxyt,bs,br,isign,ntpose,mixup,
                  sct,&ttp,indx,indy,kstrt,nvp,nxeh,nye,kxp,kyp,nypmx,
                  nxhy,nxyh);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tfft[0] += time;
      tfft[1] += ttp;

/* copy guard cells with OpenMP: updates fxyze, bxyze */
      dtimer(&dtime,&itime,-1);
      cppncguard2l(fxyze,nyp,kstrt,nvp,nnxe,nypmx);
      cppcguard2xl(fxyze,nyp,nx,ndim,nxe,nypmx);
      cppncguard2l(bxyze,nyp,kstrt,nvp,nnxe,nypmx);
      cppcguard2xl(bxyze,nyp,nx,ndim,nxe,nypmx);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tguard += time;

/* push particles with OpenMP: */
      dtimer(&dtime,&itime,-1);
      wke = 0.0;
      if (relativity==1) {
/* updates ppart and wke */
/*       cppgrbppush23l(ppart,fxyze,bxyze,kpic,noff,nyp,qbme,dt,dth,ci, */
/*                      &wke,idimp,nppmx0,nx,ny,mx,my,nxe,nypmx,mx1,    */
/*                      mxyp1,ipbc);                                    */
/* updates ppart, ncl, iholep, wke, irc */
         cppgrbppushf23l(ppart,fxyze,bxyze,kpic,ncl,iholep,noff,nyp,
                         qbme,dt,dth,ci,&wke,idimp,nppmx0,nx,ny,mx,my,
                         nxe,nypmx,mx1,mxyp1,ntmaxp,&irc);
      }
      else {
/* updates ppart and wke */
/*       cppgbppush23l(ppart,fxyze,bxyze,kpic,noff,nyp,qbme,dt,dth,&wke, */
/*                     idimp,nppmx0,nx,ny,mx,my,nxe,nypmx,mx1,mxyp1,     */
/*                     ipbc);                                            */
/* updates ppart, ncl, iholep, wke, irc */
         cppgbppushf23l(ppart,fxyze,bxyze,kpic,ncl,iholep,noff,nyp,qbme,
                        dt,dth,&wke,idimp,nppmx0,nx,ny,mx,my,nxe,nypmx,
                        mx1,mxyp1,ntmaxp,&irc);
      }
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tpush += time;
      if (irc != 0) {
         if (relativity==1) {
            printf("%d,cppgrbppushf23l error: irc=%d\n",kstrt,irc);
         }
         else {
            printf("%d,cppgbppushf23l error: irc=%d\n",kstrt,irc);
         }
         cppabort();
         exit(1);
      }

/* reorder particles by tile with OpenMP */
/* first part of particle reorder on x and y cell with mx, my tiles: */
      dtimer(&dtime,&itime,-1);
/* updates ppart, ppbuff, sbufl, sbufr, ncl, iholep, ncll, nclr, irc */
/*    cppporder2la(ppart,ppbuff,sbufl,sbufr,kpic,ncl,iholep,ncll,nclr, */
/*                 noff,nyp,idimp,nppmx0,nx,ny,mx,my,mx1,myp1,npbmx,   */
/*                 ntmaxp,nbmaxp,&irc);                                */
/* updates ppart, ppbuff, sbufl, sbufr, ncl, ncll, nclr, irc */
      cppporderf2la(ppart,ppbuff,sbufl,sbufr,ncl,iholep,ncll,nclr,
                    idimp,nppmx0,mx1,myp1,npbmx,ntmaxp,nbmaxp,&irc);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tsort += time;
      if (irc != 0) {
         printf("%d,cppporderf2la error:ntmaxp,irc=%d,%d\n",kstrt,
                ntmaxp,irc);
         cppabort();
         exit(1);
      }
/* move particles into appropriate spatial regions: */
/* updates rbufr, rbufl, mcll, mclr */
      dtimer(&dtime,&itime,-1);
      cpppmove2(sbufr,sbufl,rbufr,rbufl,ncll,nclr,mcll,mclr,kstrt,nvp,
                idimp,nbmaxp,mx1);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tmov += time;
/* second part of particle reorder on x and y cell with mx, my tiles: */
/* updates ppart, kpic */
      dtimer(&dtime,&itime,-1);
      cppporder2lb(ppart,ppbuff,rbufl,rbufr,kpic,ncl,iholep,mcll,mclr,
                   idimp,nppmx0,mx1,myp1,npbmx,ntmaxp,nbmaxp,&irc);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tsort += time;
      if (irc != 0) {
         printf("%d,cppporder2lb error:nppmx0,irc=%d,%d\n",kstrt,nppmx0,
                irc);
         cppabort();
         exit(1);
      }

/* energy diagnostic */
      wt = we + wf + wm;
      wtot[0] = wt;
      wtot[1] = wke;
      wtot[2] = 0.0;
      wtot[3] = wke + wt;
      wtot[4] = we;
      wtot[5] = wf;
      wtot[6] = wm;
      cppdsum(wtot,work,7);
      wke = wtot[1];
      we = wtot[4];
      wf = wtot[5];
      wm = wtot[6];
      if (ntime==0) {
         if (kstrt==1) {
            wt = we + wf + wm;
            printf("Initial Total Field, Kinetic and Total Energies:\n");
            printf("%e %e %e\n",wt,wke,wke+wt);
            printf("Initial Electrostatic, Transverse Electric and \
Magnetic Field Energies:\n");
            printf("%e %e %e\n",we,wf,wm);
         }
      }
      ntime += 1;
      goto L500;
L2000:

/* * * * end main iteration loop * * * */

   if (kstrt==1) {
      printf("ntime, relativity = %i,%i\n",ntime,relativity);
      printf("MPI nodes nvp = %i\n",nvp);
      wt = we + wf + wm;
      printf("Final Total Field, Kinetic and Total Energies:\n");
      printf("%e %e %e\n",wt,wke,wke+wt);
      printf("Final Electrostatic, Transverse Electric and Magnetic \
Field Energies:\n");
      printf("%e %e %e\n",we,wf,wm);

      printf("\n");
      printf("deposit time = %f\n",tdpost);
      printf("current deposit time = %f\n",tdjpost);
      tdpost += tdjpost;
      printf("total deposit time = %f\n",tdpost);
      printf("guard time = %f\n",tguard);
      printf("solver time = %f\n",tfield);
      printf("fft and transpose time = %f,%f\n",tfft[0],tfft[1]);
      printf("push time = %f\n",tpush);
      printf("particle move time = %f\n",tmov);
      printf("sort time = %f\n",tsort);
      tfield += tguard + tfft[0];
      printf("total solver time = %f\n",tfield);
      tsort += tmov;
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
   }

L3000:
   cppexit();
   return 0;
}
