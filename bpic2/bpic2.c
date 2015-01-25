/*---------------------------------------------------------------------*/
/* Skeleton 2-1/2D Electromagnetic PIC code */
/* written by Viktor K. Decyk, UCLA */
#include <stdlib.h>
#include <stdio.h>
#include <complex.h>
#include <sys/time.h>
#include "bpush2.h"

void dtimer(double *time, struct timeval *itime, int icntrl);

int main(int argc, char *argv[]) {
   int indx =   9, indy =   9;
   int npx =  3072, npy =   3072;
   int ndim = 3;
   float tend = 10.0, dt = 0.04, qme = -1.0;
   float vtx = 1.0, vty = 1.0, vx0 = 0.0, vy0 = 0.0;
   float vtz = 1.0, vz0 = 0.0;
   float ax = .912871, ay = .912871, ci = 0.1;
/* idimp = dimension of phase space = 5 */
/* sortime = number of time steps between standard electron sorting */
/* relativity = (no,yes) = (0,1) = relativity is used */
   int idimp = 5, ipbc = 1, sortime = 50, relativity = 1;
   float wke = 0.0, we = 0.0, wf = 0.0, wm = 0.0, wt = 0.0;

/* declare scalars for standard code */
   int j;
   int np, nx, ny, nxh, nyh, nxe, nye, nxeh, nxyh, nxhy;
   int ny1, ntime, nloop, isign;
   float qbme, affp, dth;

/* declare arrays for standard code */
   float *part = NULL, *part2 = NULL, *tpart = NULL;
   float *qe = NULL, *cue = NULL, *fxyze = NULL, *bxyze = NULL;
   float complex *exyz = NULL, *bxyz = NULL;
   float complex *ffc = NULL, *sct = NULL;
   int *mixup = NULL, *npicy = NULL;

/* declare and initialize timing data */
   float time;
   struct timeval itime;
   float tdpost = 0.0, tguard = 0.0, tfft = 0.0, tfield = 0.0;
   float tdjpost = 0.0, tpush = 0.0, tsort = 0.0;
   double dtime;

/* initialize scalars for standard code */
   np = npx*npy; nx = 1L<<indx; ny = 1L<<indy; nxh = nx/2; nyh = ny/2;
   nxe = nx + 2; nye = ny + 1; nxeh = nxe/2;
   nxyh = (nx > ny ? nx : ny)/2; nxhy = nxh > ny ? nxh : ny;
   ny1 = ny + 1;
   nloop = tend/dt + .0001; ntime = 0;
   qbme = qme;
   affp = (float) (nx*ny)/(float ) np;
   dth = 0.0;

/* allocate and initialize data for standard code */
   part = (float *) malloc(idimp*np*sizeof(float));
   part2 = (float *) malloc(idimp*np*sizeof(float));
   qe = (float *) malloc(nxe*nye*sizeof(float));
   fxyze = (float *) malloc(ndim*nxe*nye*sizeof(float));
   cue = (float *) malloc(ndim*nxe*nye*sizeof(float));
   bxyze = (float *) malloc(ndim*nxe*nye*sizeof(float));
   exyz = (float complex *) malloc(ndim*nxeh*nye*sizeof(float complex));
   bxyz = (float complex *) malloc(ndim*nxeh*nye*sizeof(float complex));
   ffc = (float complex *) malloc(nxh*nyh*sizeof(float complex));
   mixup = (int *) malloc(nxhy*sizeof(int));
   sct = (float complex *) malloc(nxyh*sizeof(float complex));
   npicy = (int *) malloc(ny1*sizeof(int));

/* prepare fft tables */
   cwfft2rinit(mixup,sct,indx,indy,nxhy,nxyh);
/* calculate form factors */
   isign = 0;
   cpois23((float complex *)qe,(float complex *)fxyze,isign,ffc,ax,ay,
            affp,&we,nx,ny,nxeh,nye,nxh,nyh);
/* initialize electrons */
   cdistr2h(part,vtx,vty,vtz,vx0,vy0,vz0,npx,npy,idimp,np,nx,ny,ipbc);

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
      if (relativity==1)
         cgrjpost2l(part,cue,qme,dth,ci,np,idimp,nx,ny,nxe,nye,ipbc);
      else
         cgjpost2l(part,cue,qme,dth,np,idimp,nx,ny,nxe,nye,ipbc);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tdjpost += time;

/* deposit charge with standard procedure: updates qe */
      dtimer(&dtime,&itime,-1);
      for (j = 0; j < nxe*nye; j++) {
         qe[j] = 0.0;
      }
      cgpost2l(part,qe,qme,np,idimp,nxe,nye);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tdpost += time;

/* add guard cells with standard procedure: updates cue, qe */
      dtimer(&dtime,&itime,-1);
      cacguard2l(cue,nx,ny,nxe,nye);
      caguard2l(qe,nx,ny,nxe,nye);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tguard += time;

/* transform charge to fourier space with standard procedure: updates qe */
      dtimer(&dtime,&itime,-1);
      isign = -1;
      cwfft2rx((float complex *)qe,isign,mixup,sct,indx,indy,nxeh,nye,
               nxhy,nxyh);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tfft += time;

/* transform current to fourier space with standard procedure: update cue */
      dtimer(&dtime,&itime,-1);
      isign = -1;
      cwfft2r3((float complex *)cue,isign,mixup,sct,indx,indy,nxeh,nye,
               nxhy,nxyh);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tfft += time;

/* take transverse part of current with standard procedure: updates cue */
      dtimer(&dtime,&itime,-1);
      ccuperp2((float complex *)cue,nx,ny,nxeh,nye);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tfield += time;

/* calculate electromagnetic fields in fourier space with standard */
/* procedure: updates exyz, bxyz                                   */
      dtimer(&dtime,&itime,-1);
      if (ntime==0) {
         cibpois23((float complex *)cue,bxyz,ffc,ci,&wm,nx,ny,nxeh,nye,
                   nxh,nyh);
         wf = 0.0;
         dth = 0.5*dt;
      }
      else {
         cmaxwel2(exyz,bxyz,(float complex *)cue,ffc,ci,dt,&wf,&wm,nx,ny,
                  nxeh,nye,nxh,nyh);
      }
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tfield += time;

/* calculate force/charge in fourier space with standard procedure: */
/* updates fxyze                                                    */
      dtimer(&dtime,&itime,-1);
      isign = -1;
      cpois23((float complex *)qe,(float complex *)fxyze,isign,ffc,ax,ay,
              affp,&we,nx,ny,nxeh,nye,nxh,nyh);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tfield += time;

/* add longitudinal and transverse electric fields with standard */
/* procedure: updates fxyze */
      dtimer(&dtime,&itime,-1);
      isign = 1;
      cemfield2((float complex *)fxyze,exyz,ffc,isign,nx,ny,nxeh,nye,nxh,
                nyh);

/* copy magnetic field with standard procedure: updates bxyze */
      isign = -1;
      cemfield2((float complex *)bxyze,bxyz,ffc,isign,nx,ny,nxeh,nye,nxh,
                nyh);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tfield += time;

/* transform electric force to real space with standard procedure: */
/* updates fxyze                                                   */
      dtimer(&dtime,&itime,-1);
      isign = 1;
      cwfft2r3((float complex *)fxyze,isign,mixup,sct,indx,indy,nxeh,nye,
               nxhy,nxyh);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tfft += time;

/* transform magnetic force to real space with standard procedure: */
/* updates bxyze                                                   */
      dtimer(&dtime,&itime,-1);
      isign = 1;
      cwfft2r3((float complex *)bxyze,isign,mixup,sct,indx,indy,nxeh,nye,
               nxhy,nxyh);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tfft += time;

/* copy guard cells with standard procedure: updates fxyze, bxyze */
      dtimer(&dtime,&itime,-1);
      cbguard2l(fxyze,nx,ny,nxe,nye);
      cbguard2l(bxyze,nx,ny,nxe,nye);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tguard += time;

/* push particles with standard procedure: updates part, wke */
      wke = 0.0;
      dtimer(&dtime,&itime,-1);
      if (relativity==1)
         cgrbpush23l(part,fxyze,bxyze,qbme,dt,dth,ci,&wke,idimp,np,nx,ny,
                     nxe,nye,ipbc);
      else
         cgbpush23l(part,fxyze,bxyze,qbme,dt,dth,&wke,idimp,np,nx,ny,nxe,
                    nye,ipbc);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tpush += time;

/* sort particles by cell for standard procedure */
      if (sortime > 0) {
         if (ntime%sortime==0) {
            dtimer(&dtime,&itime,-1);
            cdsortp2yl(part,part2,npicy,idimp,np,ny1);
/* exchange pointers */
            tpart = part;
            part = part2;
            part2 = tpart;
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

   return 0;
}
