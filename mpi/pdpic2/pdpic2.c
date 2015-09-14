/*--------------------------------------------------------------------*/
/* Skeleton 2-1/2D Darwin MPI PIC code */
/* written by Viktor K. Decyk, UCLA */
#include <stdlib.h>
#include <stdio.h>
#include <complex.h>
#include <sys/time.h>
#include "pdpush2.h"
#include "pplib2.h"

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
   float tend = 10.0, dt = 0.1, qme = -1.0;
/* vtx/vty = thermal velocity of electrons in x/y direction */
/* vx0/vy0 = drift velocity of electrons in x/y direction */
   float vtx = 1.0, vty = 1.0, vx0 = 0.0, vy0 = 0.0;
/* vtx/vz0 = thermal/drift velocity of electrons in z direction */
   float vtz = 1.0, vz0 = 0.0;
/* ax/ay = smoothed particle size in x/y direction */
/* ci = reciprocal of velocity of light */
   float ax = .912871, ay = .912871, ci = 0.1;
/* idimp = number of particle coordinates = 5 */
/* ipbc = particle boundary condition: 1 = periodic */
/* sortime = number of time steps between standard electron sorting */
   int idimp = 5, ipbc = 1, sortime = 50;
/* omx/omy/omz = magnetic field electron cyclotron frequency in x/y/z */
   float omx = 0.4, omy = 0.0, omz = 0.0;
/* ndc = number of corrections in darwin iteration */
   int ndc = 1;
/* idps = number of partition boundaries */
   int idps = 2;
/* wke/we = particle kinetic/electrostatic field energy */
/* wf/wm/wt = magnetic field/transverse electric field/total energy */
   float wke = 0.0, we = 0.0, wf = 0.0, wm = 0.0, wt = 0.0;
   float zero = 0.0;
/* declare scalars for standard code */
   int j, k;
   int nx, ny, nxh, nyh, nxe, nye, nxeh, nnxe, nxyh, nxhy;
   int mdim, ny1, ntime, nloop, isign, ierr;
   float qbme, affp, q2m0, wpm, wpmax, wpmin;
   double np;

/* declare scalars for MPI code */
   int ntpose = 1;
   int nvp, idproc, kstrt, npmax, kxp, kyp, nypmx, nypmn;
   int nyp, noff, npp, nps, nbmax, ntmax;

/* declare arrays for standard code: */
/* part, part2 = particle arrays */
   float *part = NULL, *part2 = NULL, *tpart = NULL;
/* qe = electron charge density with guard cells */
   float *qe = NULL;
/* cue = electron current density with guard cells */
/* dcu = acceleration density with guard cells */
/* cus = transverse electric field with guard cells */
/* amu = momentum flux with guard cells */
   float *cue = NULL, *dcu = NULL, *cus = NULL, *amu = NULL;
/* exyze = smoothed total electric field with guard cells */
/* fxyze = smoothed longitudinal electric field with guard cells */
/* bxyze = smoothed magnetic field with guard cells */
   float *fxyze = NULL, *exyze = NULL, *bxyze = NULL;
/* ss = scratch array for cwppfft2rn */
   float complex *ss = NULL;
/* qt = scalar charge density field array in fourier space */
   float complex *qt = NULL;
/* cut = vector current density in fourier space */
/* dcut = vector acceleration density in fourier space */
/* exyt = vector transverse electric field in fourier space */
/* amut = tensor momentum flux in fourier space */
   float complex *cut = NULL, *dcut = NULL, *exyt = NULL, *amut = NULL;
/* fxyt = vector longitudinal electric field in fourier space */
/* bxyt = vector magnetic field in fourier space */
   float complex *fxyt = NULL, *bxyt = NULL;
/* ffc, ffe = form factor arrays for poisson solvers */
/* sct = sine/cosine table for FFT */
   float complex *ffc = NULL, *ffe = NULL;
/* mixup = bit reverse table for FFT */
   int *mixup = NULL;
/* sct = sine/cosine table for FFT */
   float complex *sct = NULL;
/* ihole = location of hole left in particle arrays */
   int *ihole = NULL;
/* npic = scratch array for reordering particles */
   int *npic = NULL;
   double wtot[7], work[7];
   int info[7];

/* declare arrays for MPI code: */
/* bs/br = complex send/receive buffers for data transpose */
   float complex *bs = NULL, *br = NULL;
/* sbufl/sbufr = particle buffers sent to nearby processors */
/* rbufl/rbufr = particle buffers received from nearby processors */
   float *sbufl = NULL, *sbufr = NULL, *rbufl = NULL, *rbufr = NULL;
/* edges[0:1] = lower:upper y boundaries of particle partition */
   float *edges = NULL;
/* scr = guard cell buffer received from nearby processors */
   float *scr = NULL;

/* declare and initialize timing data */
   float time;
   struct timeval itime;
   float tdpost = 0.0, tguard = 0.0, ttp = 0.0, tfield = 0.0;
   float tdjpost = 0.0, tdcjpost = 0.0, tpush = 0.0, tsort = 0.0;
   float tmov = 0.0;
   float tfft[2] = {0.0,0.0};
   double dtime;

/* initialize scalars for standard code */
/* np = total number of particles in simulation */
   np =  (double) npx*(double) npy;
/* nx/ny = number of grid points in x/y direction */
   nx = 1L<<indx; ny = 1L<<indy;
   nxh = nx/2; nyh = 1 > ny/2 ? 1 : ny/2;
   nxe = nx + 2; nye = ny + 2; nxeh = nxe/2; nnxe = ndim*nxe;
   nxyh = (nx > ny ? nx : ny)/2; nxhy = nxh > ny ? nxh : ny;
   ny1 = ny + 1;
/* nloop = number of time steps in simulation */
/* ntime = current time step */
   nloop = tend/dt + .0001; ntime = 0;
/* mdim = dimension of amu array */
   mdim = 2*ndim - 2;
   qbme = qme;
   affp = (double) nx*(double) ny/np;

/* nvp = number of MPI ranks */
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
/* nbmax = size of buffer for passing particles between processors */
   nbmax = 0.1*npmax;
/* ntmax = size of ihole buffer for particles leaving processor */
   ntmax = 2*nbmax;

/* allocate data for standard code */
   part = (float *) malloc(idimp*npmax*sizeof(float));
   if (sortime > 0)
      part2 = (float *) malloc(idimp*npmax*sizeof(float));
   qe = (float *) malloc(nxe*nypmx*sizeof(float));
   fxyze = (float *) malloc(ndim*nxe*nypmx*sizeof(float));
   exyze = (float *) malloc(ndim*nxe*nypmx*sizeof(float));
   cue = (float *) malloc(ndim*nxe*nypmx*sizeof(float));
   dcu = (float *) malloc(ndim*nxe*nypmx*sizeof(float));
   cus = (float *) malloc(ndim*nxe*nypmx*sizeof(float));
   amu = (float *) malloc(mdim*nxe*nypmx*sizeof(float));
   bxyze = (float *) malloc(ndim*nxe*nypmx*sizeof(float));
   qt = (float complex *) malloc(nye*kxp*sizeof(float complex));
   fxyt = (float complex *) malloc(ndim*nye*kxp*sizeof(float complex));
   cut = (float complex *) malloc(ndim*nye*kxp*sizeof(float complex));
   dcut = (float complex *) malloc(ndim*nye*kxp*sizeof(float complex));
   exyt = (float complex *) malloc(ndim*nye*kxp*sizeof(float complex));
   amut = (float complex *) malloc(mdim*nye*kxp*sizeof(float complex));
   bxyt = (float complex *) malloc(ndim*nye*kxp*sizeof(float complex));
   ffc = (float complex *) malloc(nyh*kxp*sizeof(float complex));
   ffe = (float complex *) malloc(nyh*kxp*sizeof(float complex));
   mixup = (int *) malloc(nxhy*sizeof(int));
   sct = (float complex *) malloc(nxyh*sizeof(float complex));
   ihole = (int *) malloc((ntmax+1)*sizeof(int));
   npic = (int *) malloc(nypmx*sizeof(int));
   ss = (float complex *) malloc(mdim*nxeh*sizeof(float complex));

/* allocate data for MPI code */
   bs = (float complex *) malloc(mdim*kxp*kyp*sizeof(float complex));
   br = (float complex *) malloc(mdim*kxp*kyp*sizeof(float complex));
   sbufl = (float *) malloc(idimp*nbmax*sizeof(float));
   sbufr = (float *) malloc(idimp*nbmax*sizeof(float));
   rbufl = (float *) malloc(idimp*nbmax*sizeof(float));
   rbufr = (float *) malloc(idimp*nbmax*sizeof(float));
   scr = (float *) malloc(mdim*nxe*sizeof(float));

/* prepare fft tables */
   cwpfft2rinit(mixup,sct,indx,indy,nxhy,nxyh);
/* calculate form factor: ffc */
   isign = 0;
   cppois23(qt,fxyt,isign,ffc,ax,ay,affp,&we,nx,ny,kstrt,nye,kxp,nyh);
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

/* find maximum and minimum initial electron density */
   for (j = 0; j < nxe*nypmx; j++) {
      qe[j] = 0.0;
   }
   cppgpost2l(part,qe,npp,noff,qme,idimp,npmax,nxe,nypmx);
   cppaguard2xl(qe,nyp,nx,nxe,nypmx);
   cppnaguard2l(qe,scr,nyp,nx,kstrt,nvp,nxe,nypmx);
   cppfwpminmx2(qe,nyp,qbme,&wpmax,&wpmin,nx,nxe,nypmx);
   wtot[0] = wpmax;
   wtot[1] = -wpmin;
   cppdmax(wtot,work,2);
   wpmax = wtot[0];
   wpmin = -wtot[1];
   wpm = 0.5*(wpmax + wpmin)*affp;
/* accelerate convergence: update wpm */
   if (wpm <= 10.0)
      wpm = 0.75*wpm;
   if (kstrt==1)
      printf("wpm=%f\n",wpm);
   q2m0 = wpm/affp;
/* calculate form factor: ffe */
   isign = 0;
   cppepoisp23(dcut,exyt,isign,ffe,ax,ay,affp,wpm,ci,&wf,nx,ny,kstrt,
               nye,kxp,nyh);

/* initialize transverse electric field */
   for (j = 0; j < ndim*nxe*nypmx; j++) {
      cus[j] = 0.0;
   }

/* * * * start main iteration loop * * * */

L500: if (nloop <= ntime)
         goto L2000;
/*    if (kstrt==1) printf("ntime = %i\n",ntime); */

/* deposit current with standard procedure: updates cue, ihole */
      dtimer(&dtime,&itime,-1);
      for (j = 0; j < ndim*nxe*nypmx; j++) {
         cue[j] = 0.0;
      }
      cppgjpost2l(part,cue,edges,npp,noff,ihole,qme,zero,nx,ny,idimp,
                  npmax,nxe,nypmx,idps,ntmax,ipbc);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tdjpost += time;

/* deposit charge with standard procedure: updates qe */
      dtimer(&dtime,&itime,-1);
      for (j = 0; j < nxe*nypmx; j++) {
         qe[j] = 0.0;
      }
      cppgpost2l(part,qe,npp,noff,qme,idimp,npmax,nxe,nypmx);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tdpost += time;

/* add guard cells with standard procedure: updates cue, qe */
      dtimer(&dtime,&itime,-1);
      cppaguard2xl(qe,nyp,nx,nxe,nypmx);
      cppnaguard2l(qe,scr,nyp,nx,kstrt,nvp,nxe,nypmx);
      cppacguard2xl(cue,nyp,nx,ndim,nxe,nypmx);
      cppnacguard2l(cue,scr,nyp,nx,ndim,kstrt,nvp,nxe,nypmx);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tguard += time;

/* transform charge to fourier space with standard procedure: updates qt */
/* modifies qe */
      dtimer(&dtime,&itime,-1);
      isign = -1;
      cwppfft2r((float complex *)qe,qt,bs,br,isign,ntpose,mixup,sct,&ttp,
                indx,indy,kstrt,nvp,nxeh,nye,kxp,kyp,nypmx,nxhy,nxyh);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tfft[0] += time;
      tfft[1] += ttp;

/* calculate longitudinal force/charge in fourier space with standard */
/* procedure: updates fxyt, we */
      dtimer(&dtime,&itime,-1);
      isign = -1;
      cppois23(qt,fxyt,isign,ffc,ax,ay,affp,&we,nx,ny,kstrt,nye,kxp,
               nyh);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tfield += time;

/* transform longitudinal electric force to real space with standard */
/* procedure: updates fxyze, modifies fxyt */
      dtimer(&dtime,&itime,-1);
      isign = 1;
      cwppfft2r3((float complex *)fxyze,fxyt,bs,br,isign,ntpose,mixup,
                 sct,&ttp,indx,indy,kstrt,nvp,nxeh,nye,kxp,kyp,nypmx,
                 nxhy,nxyh);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tfft[0] += time;
      tfft[1] += ttp;

/* transform current to fourier space with standard procedure: updates cut */
/* modifies cue */
      dtimer(&dtime,&itime,-1);
      isign = -1;
      cwppfft2r3((float complex *)cue,cut,bs,br,isign,ntpose,mixup,sct,
                 &ttp,indx,indy,kstrt,nvp,nxeh,nye,kxp,kyp,nypmx,nxhy,
                 nxyh);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tfft[0] += time;
      tfft[1] += ttp;

/* take transverse part of current with standard procedure: updates cut */
      dtimer(&dtime,&itime,-1);
      cppcuperp2(cut,nx,ny,kstrt,nye,kxp);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tfield += time;

/* calculate magnetic field in fourier space with standard procedure: */
/* updates bxyt, wm */
      dtimer(&dtime,&itime,-1);
      cppbbpoisp23(cut,bxyt,ffc,ci,&wm,nx,ny,kstrt,nye,kxp,nyh);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tfield += time;

/* transform magnetic field to real space with standard procedure: */
/* updates bxyze, modifies bxyt */
      dtimer(&dtime,&itime,-1);
      isign = 1;
      cwppfft2r3((float complex *)bxyze,bxyt,bs,br,isign,ntpose,mixup,
                 sct,&ttp,indx,indy,kstrt,nvp,nxeh,nye,kxp,kyp,nypmx,
                 nxhy,nxyh);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tfft[0] += time;
      tfft[1] += ttp;

/* add constant to magnetic field with standard procedure: updates bxyze */
      dtimer(&dtime,&itime,-1);
      cppbaddext2(bxyze,nyp,omx,omy,omz,nx,nxe,nypmx);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tfield += time;

/* copy guard cells with standard procedure: updates fxyze, bxyze */
      dtimer(&dtime,&itime,-1);
      cppncguard2l(fxyze,nyp,kstrt,nvp,nnxe,nypmx);
      cppcguard2xl(fxyze,nyp,nx,ndim,nxe,nypmx);
      cppncguard2l(bxyze,nyp,kstrt,nvp,nnxe,nypmx);
      cppcguard2xl(bxyze,nyp,nx,ndim,nxe,nypmx);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tguard += time;

/* add longitudinal and old transverse electric fields with standard */
/* procedure: updates exyze */
      dtimer(&dtime,&itime,-1);
      cppaddvrfield2(exyze,cus,fxyze,ndim,nxe,nypmx);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tfield += time;

/* deposit electron acceleration density and momentum flux with */
/* standard procedure: updates dcu, amu                         */
      dtimer(&dtime,&itime,-1);
      for (j = 0; j < ndim*nxe*nypmx; j++) {
         dcu[j] = 0.0;
      }
      for (j = 0; j < mdim*nxe*nypmx; j++) {
         amu[j] = 0.0;
      }
      cppgdjpost2l(part,exyze,bxyze,npp,noff,dcu,amu,qme,qbme,dt,idimp,
                   npmax,nxe,nypmx);
/* add old scaled electric field with standard procedure: updates dcu */
      cppascfguard2l(dcu,cus,nyp,q2m0,nx,nxe,nypmx);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tdcjpost += time;

/* add guard cells with standard procedure: updates dcu, amu */
      dtimer(&dtime,&itime,-1);
      cppacguard2xl(dcu,nyp,nx,ndim,nxe,nypmx);
      cppnacguard2l(dcu,scr,nyp,nx,ndim,kstrt,nvp,nxe,nypmx);
      cppacguard2xl(amu,nyp,nx,mdim,nxe,nypmx);
      cppnacguard2l(amu,scr,nyp,nx,mdim,kstrt,nvp,nxe,nypmx);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tguard += time;

/* transform acceleration density and momentum flux to fourier space */
/* with standard procedure: updates dcut, amut, modifies dcu, amu */
      dtimer(&dtime,&itime,-1);
      isign = -1;
      cwppfft2r3((float complex *)dcu,dcut,bs,br,isign,ntpose,mixup,sct,
                 &ttp,indx,indy,kstrt,nvp,nxeh,nye,kxp,kyp,nypmx,nxhy,
                 nxyh);
      tfft[1] += ttp;
      cwppfft2rn((float complex *)amu,amut,bs,br,ss,isign,ntpose,mixup,
                 sct,&ttp,indx,indy,kstrt,nvp,nxeh,nye,kxp,kyp,nypmx,
                 mdim,nxhy,nxyh);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tfft[0] += time;
      tfft[1] += ttp;

/* take transverse part of time derivative of current with standard */
/* procedure: updates dcut                                          */
      dtimer(&dtime,&itime,-1);
      cppadcuperp23(dcut,amut,nx,ny,kstrt,nye,kxp);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tfield += time;

/* calculate transverse electric field with standard procedure: */
/* updates exyt, wf                                             */
      dtimer(&dtime,&itime,-1);
      isign = -1;
      cppepoisp23(dcut,exyt,isign,ffe,ax,ay,affp,wpm,ci,&wf,nx,ny,kstrt,
                  nye,kxp,nyh);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tfield += time;

/* transform transverse electric field to real space with standard */
/* procedure: updates cus, modifies exyt */
      dtimer(&dtime,&itime,-1);
      isign = 1;
      cwppfft2r3((float complex *)cus,exyt,bs,br,isign,ntpose,mixup,sct,
                 &ttp,indx,indy,kstrt,nvp,nxeh,nye,kxp,kyp,nypmx,nxhy,
                 nxyh);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tfft[0] += time;
      tfft[1] += ttp;

/* copy guard cells with standard procedure: updates cus */
      dtimer(&dtime,&itime,-1);
      cppncguard2l(cus,nyp,kstrt,nvp,nnxe,nypmx);
      cppcguard2xl(cus,nyp,nx,ndim,nxe,nypmx);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tguard += time;

/* add longitudinal and transverse electric fields with standard */
/* procedure: exyze = cus + fxyze, updates exyze                 */
/* cus needs to be retained for next time step */
      dtimer(&dtime,&itime,-1);
      cppaddvrfield2(exyze,cus,fxyze,ndim,nxe,nypmx);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tfield += time;

/* inner iteration loop */
      for (k = 0; k < ndc; k++) {

/* deposit electron current and acceleration density and momentum flux */
/* with standard procedure: updates cue, dcu, amu                      */
         dtimer(&dtime,&itime,-1);
         for (j = 0; j < ndim*nxe*nypmx; j++) {
            cue[j] = 0.0;
            dcu[j] = 0.0;
         }
         for (j = 0; j < mdim*nxe*nypmx; j++) {
            amu[j] = 0.0;
         }
         cppgdcjpost2l(part,exyze,bxyze,npp,noff,cue,dcu,amu,qme,qbme,
                       dt,idimp,npmax,nxe,nypmx);
/* add scaled electric field with standard procedure: updates dcu */
         cppascfguard2l(dcu,cus,nyp,q2m0,nx,nxe,nypmx);
         dtimer(&dtime,&itime,1);
         time = (float) dtime;
         tdcjpost += time;

/* add guard cells for current, acceleration density, and momentum flux */
/* with standard procedure: updates cue, dcu, amu                       */
         dtimer(&dtime,&itime,-1);
         cppacguard2xl(cue,nyp,nx,ndim,nxe,nypmx);
         cppnacguard2l(cue,scr,nyp,nx,ndim,kstrt,nvp,nxe,nypmx);
         cppacguard2xl(dcu,nyp,nx,ndim,nxe,nypmx);
         cppnacguard2l(dcu,scr,nyp,nx,ndim,kstrt,nvp,nxe,nypmx);
         cppacguard2xl(amu,nyp,nx,mdim,nxe,nypmx);
         cppnacguard2l(amu,scr,nyp,nx,mdim,kstrt,nvp,nxe,nypmx);
         dtimer(&dtime,&itime,1);
         time = (float) dtime;
         tguard += time;

/* transform current to fourier space with standard procedure: update cut */
/* modifies cue */
         dtimer(&dtime,&itime,-1);
         isign = -1;
         cwppfft2r3((float complex *)cue,cut,bs,br,isign,ntpose,mixup,
                    sct,&ttp,indx,indy,kstrt,nvp,nxeh,nye,kxp,kyp,nypmx,
                    nxhy,nxyh);
         dtimer(&dtime,&itime,1);
         time = (float) dtime;
         tfft[0] += time;
         tfft[1] += ttp;

/* take transverse part of current with standard procedure: updates cut */
         dtimer(&dtime,&itime,-1);
         cppcuperp2(cut,nx,ny,kstrt,nye,kxp);
         dtimer(&dtime,&itime,1);
         time = (float) dtime;
         tfield += time;

/* calculate magnetic field in fourier space with standard procedure: */
/* updates bxyt, wm                                                   */
         dtimer(&dtime,&itime,-1);
         cppbbpoisp23(cut,bxyt,ffc,ci,&wm,nx,ny,kstrt,nye,kxp,nyh);
         dtimer(&dtime,&itime,1);
         time = (float) dtime;
         tfield += time;

/* transform magnetic force to real space with standard procedure: */
/* updates bxyze, modifies bxyt                                    */
         dtimer(&dtime,&itime,-1);
         isign = 1;
         cwppfft2r3((float complex *)bxyze,bxyt,bs,br,isign,ntpose,mixup,
                    sct,&ttp,indx,indy,kstrt,nvp,nxeh,nye,kxp,kyp,nypmx,
                    nxhy,nxyh);
         dtimer(&dtime,&itime,1);
         time = (float) dtime;
         tfft[0] += time;
         tfft[1] += ttp;

/* add constant to magnetic field with standard procedure: updates bxyze */
         dtimer(&dtime,&itime,-1);
         cppbaddext2(bxyze,nyp,omx,omy,omz,nx,nxe,nypmx);
         dtimer(&dtime,&itime,1);
         time = (float) dtime;
         tfield += time;

/* transform acceleration density and momentum flux to fourier space */
/* with standard procedure: updates dcut, amut, modifies dcu, amu */
         dtimer(&dtime,&itime,-1);
         isign = -1;
         cwppfft2r3((float complex *)dcu,dcut,bs,br,isign,ntpose,mixup,
                    sct,&ttp,indx,indy,kstrt,nvp,nxeh,nye,kxp,kyp,nypmx,
                    nxhy,nxyh);
         tfft[1] += ttp;
         cwppfft2rn((float complex *)amu,amut,bs,br,ss,isign,ntpose,
                    mixup,sct,&ttp,indx,indy,kstrt,nvp,nxeh,nye,kxp,kyp,
                    nypmx,mdim,nxhy,nxyh);
         dtimer(&dtime,&itime,1);
         time = (float) dtime;
         tfft[0] += time;
         tfft[1] += ttp;
 
/* take transverse part of time derivative of current with standard */
/* procedure: updates dcut                                          */
         dtimer(&dtime,&itime,-1);
         cppadcuperp23(dcut,amut,nx,ny,kstrt,nye,kxp);
         dtimer(&dtime,&itime,1);
         time = (float) dtime;
         tfield += time;
 
/* calculate convective part of transverse electric field with standard */
/* procedure: updates exyt, wf                                          */
         dtimer(&dtime,&itime,-1);
         isign = -1;
         cppepoisp23(dcut,exyt,isign,ffe,ax,ay,affp,wpm,ci,&wf,nx,ny,
                     kstrt,nye,kxp,nyh);
         dtimer(&dtime,&itime,1);
         time = (float) dtime;
         tfield += time;
 
/* transform transverse electric field to real space with standard */
/* procedure: updates cus, modifies exyt */
         dtimer(&dtime,&itime,-1);
         isign = 1;
         cwppfft2r3((float complex *)cus,exyt,bs,br,isign,ntpose,mixup,
                    sct,&ttp,indx,indy,kstrt,nvp,nxeh,nye,kxp,kyp,nypmx,
                    nxhy,nxyh);
         dtimer(&dtime,&itime,1);
         time = (float) dtime;
         tfft[0] += time;
         tfft[1] += ttp;
 
/* copy guard cells with standard procedure: updates bxyze, cus */
         dtimer(&dtime,&itime,-1);
         cppncguard2l(bxyze,nyp,kstrt,nvp,nnxe,nypmx);
         cppcguard2xl(bxyze,nyp,nx,ndim,nxe,nypmx);
         cppncguard2l(cus,nyp,kstrt,nvp,nnxe,nypmx);
         cppcguard2xl(cus,nyp,nx,ndim,nxe,nypmx);
         dtimer(&dtime,&itime,1);
         time = (float) dtime;
         tguard += time;
 
/* add longitudinal and transverse electric fields with standard */
/* procedure: exyze = cus + fxyze, updates exyze                 */
/* cus needs to be retained for next time step                   */
         dtimer(&dtime,&itime,-1);
         cppaddvrfield2(exyze,cus,fxyze,ndim,nxe,nypmx);
         dtimer(&dtime,&itime,1);
         time = (float) dtime;
         tfield += time;

      }

/* push particles: updates part, wke, and ihole */
      dtimer(&dtime,&itime,-1);
      wke = 0.0;
      cppgbpush23l(part,exyze,bxyze,edges,npp,noff,ihole,qbme,dt,dt,
                   &wke,nx,ny,idimp,npmax,nxe,nypmx,idps,ntmax,ipbc);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tpush += time;
/* check for ihole overflow error */
      if (ihole[0] < 0) {
         ierr = -ihole[0];
         printf("ihole overflow error: ntmax,ih=%d,%d\n",ntmax,ierr);
         cppabort();
         goto L3000;
      }

/* move electrons into appropriate spatial regions: updates part, npp */
      dtimer(&dtime,&itime,-1);
      cppmove2(part,edges,&npp,sbufr,sbufl,rbufr,rbufl,ihole,ny,kstrt,nvp,
               idimp,npmax,idps,nbmax,ntmax,info);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tmov += time;
/* check for particle manager error */
      if (info[0] != 0) {
         ierr = info[0];
         if (kstrt==1) {
            printf("push particle manager error: ierr=%d\n",ierr);
         }
         goto L3000;
      }

/* sort particles for standard code: updates part */
      if (sortime > 0) {
         if (ntime%sortime==0) {
            dtimer(&dtime,&itime,-1);
            cppdsortp2yl(part,part2,npic,npp,noff,nyp,idimp,npmax,nypmx);
/* exchange pointers */
            tpart = part;
            part = part2;
            part2 = tpart;
            dtimer(&dtime,&itime,1);
            time = (float) dtime;
            tsort += time;
         }
      }

/* energy diagnostic */
      wt = we + wm;
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
            wt = we + wm;
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
      printf("ntime, ndc = %i,%i\n",ntime,ndc);
      printf("MPI nodes nvp = %i\n",nvp);
      wt = we + wm;
      printf("Final Total Field, Kinetic and Total Energies:\n");
      printf("%e %e %e\n",wt,wke,wke+wt);
      printf("Final Electrostatic, Transverse Electric and Magnetic \
Field Energies:\n");
      printf("%e %e %e\n",we,wf,wm);

      printf("\n");
      printf("deposit time = %f\n",tdpost);
      printf("current deposit time = %f\n",tdjpost);
      printf("current derivative deposit time = %f\n",tdcjpost);
      tdpost += tdjpost + tdcjpost;
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
