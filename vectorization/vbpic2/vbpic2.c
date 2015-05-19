/*---------------------------------------------------------------------*/
/* Skeleton 2-1/2D Electromagnetic Vector PIC code */
/* written by Viktor K. Decyk, UCLA and Ricardo Fonseca, ISCTE */
#include <stdlib.h>
#include <stdio.h>
#include <complex.h>
#include <sys/time.h>
#include "vbpush2.h"
#include "sselib2.h"
#include "ssebpush2.h"

void dtimer(double *time, struct timeval *itime, int icntrl);

int main(int argc, char *argv[]) {
/* indx/indy = exponent which determines grid points in x/y direction: */
/* nx = 2**indx, ny = 2**indy */
   int indx =   9, indy =   9;
/* npx/npy = number of electrons distributed in x/y direction */
   int npx =  3072, npy =   3072;
/* ndim = number of velocity coordinates = 3 */
   int ndim = 4;
/* tend = time at end of simulation, in units of plasma frequency */
/* dt = time interval between successive calculations */
/* qme = charge on electron, in units of e */
   float tend = 10.0, dt = 0.04, qme = -1.0;
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
/* relativity = (no,yes) = (0,1) = relativity is used */
   int idimp = 5, ipbc = 1, sortime = 50, relativity = 1;
   float wke = 0.0, we = 0.0, wf = 0.0, wm = 0.0, wt = 0.0;
/* kvec = (1,2) = run (autovector,SSE2) version */
   int kvec = 1;

/* declare scalars for standard code */
   int j;
   int np, nx, ny, nxh, nyh, nxe, nye, nxeh, nxyh, nxhy;
   int npe, ny1, ntime, nloop, isign;
   int irc = 0;
   float qbme, affp, dth;

/* declare arrays for standard code: */
/* partt, partt2 = transposed particle arrays */
   float *partt = NULL, *partt2 = NULL, *tpartt = NULL;
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
/* npicy = scratch array for reordering particles */
   int *mixup = NULL, *npicy = NULL;

/* declare and initialize timing data */
   float time;
   struct timeval itime;
   float tdpost = 0.0, tguard = 0.0, tfft = 0.0, tfield = 0.0;
   float tdjpost = 0.0, tpush = 0.0, tsort = 0.0;
   double dtime;

/* initialize scalars for standard code */
/* np = total number of particles in simulation */
/* nx/ny = number of grid points in x/y direction */
   np = npx*npy; nx = 1L<<indx; ny = 1L<<indy; nxh = nx/2; nyh = ny/2;
   nxe = nx + 2; nye = ny + 1; nxeh = nxe/2;
   nxyh = (nx > ny ? nx : ny)/2; nxhy = nxh > ny ? nxh : ny;
   ny1 = ny + 1;
/* nloop = number of time steps in simulation */
/* ntime = current time step */
   nloop = tend/dt + .0001; ntime = 0;
   qbme = qme;
   affp = (float) (nx*ny)/(float ) np;
   dth = 0.0;

/* allocate data for standard code */
   mixup = (int *) malloc(nxhy*sizeof(int));
   sct = (float complex *) malloc(nxyh*sizeof(float complex));

/* align memory for SSE */
   npe = 4*((np - 1)/4 + 1);
   nxe = 4*((nxe - 1)/4 + 1);
   nxeh = nxe/2;
   sse_fallocate(&partt,npe*idimp,&irc);
   if (sortime > 0)
      sse_fallocate(&partt2,npe*idimp,&irc);
   sse_fallocate(&qe,nxe*nye,&irc);
   sse_fallocate(&cue,ndim*nxe*nye,&irc);
   sse_fallocate(&fxyze,ndim*nxe*nye,&irc);
   sse_fallocate(&bxyze,ndim*nxe*nye,&irc);
   sse_callocate(&exyz,ndim*nxeh*nye,&irc);
   sse_callocate(&bxyz,ndim*nxeh*nye,&irc);
   sse_callocate(&ffc,nxh*nyh,&irc);
   sse_iallocate(&npicy,ny1,&irc);
   if (irc != 0) {
      printf("aligned allocation error: irc = %d\n",irc);
   }

/* prepare fft tables */
   cwfft2rinit(mixup,sct,indx,indy,nxhy,nxyh);
/* calculate form factors */
   isign = 0;
   cvpois23((float complex *)qe,(float complex *)fxyze,isign,ffc,ax,ay,
             affp,&we,nx,ny,nxeh,nye,nxh,nyh);
/* initialize electrons */
   cdistr2ht(partt,vtx,vty,vtz,vx0,vy0,vz0,npx,npy,idimp,npe,nx,ny,ipbc);

/* initialize transverse electromagnetic fields */
   for (j = 0; j < ndim*nxeh*nye; j++) {
      exyz[j] = 0.0 + 0.0*_Complex_I;
      bxyz[j] = 0.0 + 0.0*_Complex_I;
   }

   if (dt > 0.45*ci) {
      printf("Warning: Courant condition may be exceeded!\n");
   }

/* * * * start main iteration loop * * * */
 
L500: if (nloop <= ntime)
         goto L2000;
/*    printf("ntime = %i\n",ntime); */
 
/* deposit current with standard procedure: updates part, cue */
      dtimer(&dtime,&itime,-1);
      for (j = 0; j < ndim*nxe*nye; j++) {
         cue[j] = 0.0;
      }
      if (relativity==1) {
         if (kvec==1)
            cvgrjpost2lt(partt,cue,qme,dth,ci,np,npe,idimp,nx,ny,nxe,
                         nye,ipbc);
/* SSE2 function */
         else if (kvec==2)
            csse2grjpost2lt(partt,cue,qme,dth,ci,np,npe,idimp,nx,ny,nxe,
                            nye,ipbc);
      }
      else {
         if (kvec==1)
            cvgjpost2lt(partt,cue,qme,dth,np,npe,idimp,nx,ny,nxe,nye,
                        ipbc);
/* SSE2 function */
         else if (kvec==2)
           csse2gjpost2lt(partt,cue,qme,dth,np,npe,idimp,nx,ny,nxe,nye,
                          ipbc);
      }
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tdjpost += time;

/* deposit charge with standard procedure: updates qe */
      dtimer(&dtime,&itime,-1);
      for (j = 0; j < nxe*nye; j++) {
         qe[j] = 0.0;
      }
      if (kvec==1)
         cvgpost2lt(partt,qe,qme,np,npe,idimp,nxe,nye);
/* SSE2 function */
      else if (kvec==2)
         csse2gpost2lt(partt,qe,qme,np,npe,idimp,nxe,nye);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tdpost += time;

/* add guard cells with standard procedure: updates cue, qe */
      dtimer(&dtime,&itime,-1);
      if (kvec==1) {
         cacguard2l(cue,nx,ny,nxe,nye);
         caguard2l(qe,nx,ny,nxe,nye);
      }
/* SSE2 function */
      else if (kvec==2) {
         csse2acguard2l(cue,nx,ny,nxe,nye);
         csse2aguard2l(qe,nx,ny,nxe,nye);
      }
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tguard += time;

/* transform charge to fourier space with standard procedure: updates qe */
      dtimer(&dtime,&itime,-1);
      isign = -1;
      if (kvec==1)
         cwfft2rvx((float complex *)qe,isign,mixup,sct,indx,indy,nxeh,
                   nye,nxhy,nxyh);
/* SSE2 function */
      else if (kvec==2)
         csse2wfft2rx((float complex *)qe,isign,mixup,sct,indx,indy,
                      nxeh,nye,nxhy,nxyh);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tfft += time;

/* transform current to fourier space with standard procedure: update cue */
      dtimer(&dtime,&itime,-1);
      isign = -1;
      if (kvec==1)
         cwfft2rv3((float complex *)cue,isign,mixup,sct,indx,indy,nxeh,
                   nye,nxhy,nxyh);
/* SSE2 function */
      else if (kvec==2)
         csse2wfft2r3((float complex *)cue,isign,mixup,sct,indx,indy,
                      nxeh,nye,nxhy,nxyh);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tfft += time;

/* take transverse part of current with standard procedure: updates cue */
      dtimer(&dtime,&itime,-1);
      if (kvec==1)
         ccuperp2((float complex *)cue,nx,ny,nxeh,nye);
/* SSE2 function */
      else if (kvec==2)
         csse2cuperp2((float complex *)cue,nx,ny,nxeh,nye);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tfield += time;

/* calculate electromagnetic fields in fourier space with standard */
/* procedure: updates exyz, bxyz                                   */
      dtimer(&dtime,&itime,-1);
      if (ntime==0) {
         if (kvec==1)
            cvibpois23((float complex *)cue,bxyz,ffc,ci,&wm,nx,ny,nxeh,
                       nye,nxh,nyh);
/* SSE2 function */
         else if (kvec==2)
            csse2ibpois23((float complex *)cue,bxyz,ffc,ci,&wm,nx,ny,
                          nxeh,nye,nxh,nyh);
         wf = 0.0;
         dth = 0.5*dt;
      }
      else {
         if (kvec==1)
            cvmaxwel2(exyz,bxyz,(float complex *)cue,ffc,ci,dt,&wf,&wm,
                      nx,ny,nxeh,nye,nxh,nyh);
/* SSE2 function */
         else if (kvec==2)
            csse2maxwel2(exyz,bxyz,(float complex *)cue,ffc,ci,dt,&wf,
                         &wm,nx,ny,nxeh,nye,nxh,nyh);
      }
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tfield += time;

/* calculate force/charge in fourier space with standard procedure: */
/* updates fxyze                                                    */
      dtimer(&dtime,&itime,-1);
      isign = -1;
      if (kvec==1)
         cvpois23((float complex *)qe,(float complex *)fxyze,isign,ffc,
                  ax,ay,affp,&we,nx,ny,nxeh,nye,nxh,nyh);
/* SSE2 function */
      else if (kvec==2)
         csse2pois23((float complex *)qe,(float complex *)fxyze,isign,
                     ffc,ax,ay,affp,&we,nx,ny,nxeh,nye,nxh,nyh);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tfield += time;

/* add longitudinal and transverse electric fields with standard */
/* procedure: updates fxyze */
      dtimer(&dtime,&itime,-1);
      isign = 1;
      if (kvec==1)
         cvemfield2((float complex *)fxyze,exyz,ffc,isign,nx,ny,nxeh,
                    nye,nxh,nyh);
/* SSE2 function */
      else if (kvec==2)
         csse2emfield2((float complex *)fxyze,exyz,ffc,isign,nx,ny,nxeh,
                       nye,nxh,nyh);
/* copy magnetic field with standard procedure: updates bxyze */
      isign = -1;
      if (kvec==1)
         cvemfield2((float complex *)bxyze,bxyz,ffc,isign,nx,ny,nxeh,
                    nye,nxh,nyh);
/* SSE2 function */
      else if (kvec==2)
         csse2emfield2((float complex *)bxyze,bxyz,ffc,isign,nx,ny,nxeh,
                       nye,nxh,nyh);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tfield += time;

/* transform electric force to real space with standard procedure: */
/* updates fxyze                                                   */
      dtimer(&dtime,&itime,-1);
      isign = 1;
      if (kvec==1)
         cwfft2rv3((float complex *)fxyze,isign,mixup,sct,indx,indy,
                   nxeh,nye,nxhy,nxyh);
/* SSE2 function */
      else if (kvec==2)
         csse2wfft2r3((float complex *)fxyze,isign,mixup,sct,indx,indy,
                      nxeh,nye,nxhy,nxyh);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tfft += time;

/* transform magnetic force to real space with standard procedure: */
/* updates bxyze                                                   */
      dtimer(&dtime,&itime,-1);
      isign = 1;
      if (kvec==1)
         cwfft2rv3((float complex *)bxyze,isign,mixup,sct,indx,indy,
                   nxeh,nye,nxhy,nxyh);
/* SSE2 function */
      else if (kvec==2)
         csse2wfft2r3((float complex *)bxyze,isign,mixup,sct,indx,indy,
                      nxeh,nye,nxhy,nxyh);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tfft += time;

/* copy guard cells with standard procedure: updates fxyze, bxyze */
      dtimer(&dtime,&itime,-1);
      if (kvec==1) {
         cbguard2l(fxyze,nx,ny,nxe,nye);
         cbguard2l(bxyze,nx,ny,nxe,nye);
      }
/* SSE2 function */
      else if (kvec==2) {
         csse2bguard2l(fxyze,nx,ny,nxe,nye);
         csse2bguard2l(bxyze,nx,ny,nxe,nye);
      }
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tguard += time;

/* push particles with standard procedure: updates part, wke */
      wke = 0.0;
      dtimer(&dtime,&itime,-1);
      if (relativity==1) {
         if (kvec==1)
            cvgrbpush23lt(partt,fxyze,bxyze,qbme,dt,dth,ci,&wke,idimp,
                          np,npe,nx,ny,nxe,nye,ipbc);
/* SSE2 function */
        else if (kvec==2)
            csse2grbpush23lt(partt,fxyze,bxyze,qbme,dt,dth,ci,&wke,
                            idimp,np,npe,nx,ny,nxe,nye,ipbc);
      }
      else {
         if (kvec==1)
            cvgbpush23lt(partt,fxyze,bxyze,qbme,dt,dth,&wke,idimp,np,
                         npe,nx,ny,nxe,nye,ipbc);
/* SSE2 function */
         else if (kvec==2)
            csse2gbpush23lt(partt,fxyze,bxyze,qbme,dt,dth,&wke,idimp,np,
                            npe,nx,ny,nxe,nye,ipbc);
      }
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tpush += time;

/* sort particles by cell for standard procedure */
      if (sortime > 0) {
         if (ntime%sortime==0) {
            dtimer(&dtime,&itime,-1);
            if (kvec==1)
               cdsortp2ylt(partt,partt2,npicy,idimp,np,npe,ny1);
/* SSE2 function */
            else if (kvec==2)
               csse2dsortp2ylt(partt,partt2,npicy,idimp,np,npe,ny1);
/* exchange pointers */
            tpartt = partt;
            partt = partt2;
            partt2 = tpartt;
            dtimer(&dtime,&itime,1);
            time = (float) dtime;
            tsort += time;
         }
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

   sse_deallocate(npicy);
   sse_deallocate(ffc);
   sse_deallocate(bxyz);
   sse_deallocate(exyz);
   sse_deallocate(bxyze);
   sse_deallocate(fxyze);
   sse_deallocate(cue);
   sse_deallocate(qe);
   if (sortime > 0)
      sse_deallocate(partt2);
   sse_deallocate(partt);

   return 0;
}
