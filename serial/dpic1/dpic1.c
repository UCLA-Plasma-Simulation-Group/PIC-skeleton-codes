/*---------------------------------------------------------------------*/
/* Skeleton 1-2/2D Darwin PIC code */
/* written by Viktor K. Decyk, UCLA */
#include <stdlib.h>
#include <stdio.h>
#include <complex.h>
#include <sys/time.h>
#include "dpush1.h"

void dtimer(double *time, struct timeval *itime, int icntrl);

int main(int argc, char *argv[]) {
/* indx = exponent which determines grid points in x direction: */
/* nx = 2**indx */
   int indx =   9;
/* npx = number of electrons distributed in x direction */
   int npx =  18432;
/* tend = time at end of simulation, in units of plasma frequency */
/* dt = time interval between successive calculations */
/* qme = charge on electron, in units of e */
   float tend = 10.0, dt = 0.1, qme = -1.0;
/* vtx/vty = thermal velocity of electrons in x/y direction */
/* vx0/vy0 = drift velocity of electrons in x/y direction */
   float vtx = 1.0, vty = 1.0, vx0 = 0.0, vy0 = 0.0;
/* vtx/vz0 = thermal/drift velocity of electrons in z direction */
   float vtz = 1.0, vz0 = 0.0;
/* ax = smoothed particle size in x direction */
/* ci = reciprocal of velocity of light */
   float ax = .912871, ci = 0.1;
/* idimp = number of particle coordinates = 4 */
/* ipbc = particle boundary condition: 1 = periodic */
/* sortime = number of time steps between standard electron sorting */
   int idimp = 4, ipbc = 1, sortime = 50;
/* omx/omy/omz = magnetic field electron cyclotron frequency in x/y/z */
   float omx = 0.4, omy = 0.0, omz = 0.0;
/* ndc = number of corrections in darwin iteration */
   int ndc = 1;
/* wke/we = particle kinetic/electrostatic field energy */
/* wf/wm/wt = magnetic field/transverse electric field/total energy */
   float wke = 0.0, we = 0.0, wf = 0.0, wm = 0.0, wt = 0.0;
   float zero = 0.0;
/* declare scalars for standard code */
   int j, k;
   int np, nx, nxh, nxe, nxeh;
   int nx1, ntime, nloop, isign;
   float qbme, affp, q2m0, wpm, wpmax, wpmin;

/* declare arrays for standard code: */
/* part, part2 = particle arrays */
   float *part = NULL, *part2 = NULL, *tpart = NULL;
/* qe = electron charge density with guard cells */
/* fxe = smoothed longitudinal electric field with guard cells */
   float *qe = NULL, *fxe = NULL;
/* cue = electron current density with guard cells */
/* dcu = acceleration density with guard cells */
/* cus = transverse electric field */
/* amu = momentum flux with guard cells */
   float *cue = NULL, *dcu = NULL, *cus = NULL, *amu = NULL;
/* exyze/byze = smoothed electric/magnetic field with guard cells */
   float *exyze = NULL, *byze = NULL;
/* ffc, ffe = form factor arrays for poisson solvers */
/* sct = sine/cosine table for FFT */
   float complex *ffc = NULL, *ffe = NULL, *sct = NULL;
/* mixup = bit reverse table for FFT */
/* npic = scratch array for reordering particles */
   int *mixup = NULL, *npic = NULL;
/* gxe, gyze, gxyze = scratch arrays for fft */
   float *gxe = NULL, *gyze = NULL;

/* declare and initialize timing data */
   float time;
   struct timeval itime;
   float tdpost = 0.0, tguard = 0.0, tfft = 0.0, tfield = 0.0;
   float tdjpost = 0.0, tdcjpost = 0.0, tpush = 0.0, tsort = 0.0;
   double dtime;

/* initialize scalars for standard code */
/* np = total number of particles in simulation */
/* nx = number of grid points in x direction */
   np = npx; nx = 1L<<indx; nxh = nx/2;
   nxe = nx + 2; nxeh = nxe/2; nx1 = nx + 1;
/* nloop = number of time steps in simulation */
/* ntime = current time step */
   nloop = tend/dt + .0001; ntime = 0;
   qbme = qme;
   affp = (float) nx/(float ) np;

/* allocate data for standard code */
   part = (float *) malloc(idimp*np*sizeof(float));
   if (sortime > 0)
      part2 = (float *) malloc(idimp*np*sizeof(float));
   qe = (float *) malloc(nxe*sizeof(float));
   fxe = (float *) malloc(nxe*sizeof(float));
   cue = (float *) malloc(2*nxe*sizeof(float));
   dcu = (float *) malloc(2*nxe*sizeof(float));
   cus = (float *) malloc(2*nxe*sizeof(float));
   amu = (float *) malloc(2*nxe*sizeof(float));
   exyze = (float *) malloc(3*nxe*sizeof(float));
   byze = (float *) malloc(2*nxe*sizeof(float));
   ffc = (float complex *) malloc(nxh*sizeof(float complex));
   ffe = (float complex *) malloc(nxh*sizeof(float complex));
   mixup = (int *) malloc(nxh*sizeof(int));
   sct = (float complex *) malloc(nxh*sizeof(float complex));
   npic = (int *) malloc(nx1*sizeof(int));
   gxe = (float *) malloc(nxe*sizeof(float));
   gyze = (float *) malloc(2*nxe*sizeof(float));

/* prepare fft tables */
   cwfft1rinit(mixup,sct,indx,nxh);
/* calculate form factor: ffc */
   isign = 0;
   cpois1((float complex *)qe,(float complex *)fxe,isign,ffc,ax,affp,
          &we,nx);
/* initialize electrons */
   cdistr1h(part,vtx,vty,vtz,vx0,vy0,vz0,npx,idimp,np,nx,ipbc);

/* find maximum and minimum initial electron density */
   for (j = 0; j < nxe; j++) {
      qe[j] = 0.0;
   }
   cgpost1l(part,qe,qme,np,idimp,nxe);
   caguard1l(qe,nx,nxe);
   cfwpminmx1(qe,qbme,&wpmax,&wpmin,nx,nxe);
   wpm = 0.5*(wpmax + wpmin)*affp;
/* accelerate convergence: update wpm */
   if (wpm <= 10.0)
      wpm = 0.75*wpm;
   printf("wpm=%f\n",wpm);
   q2m0 = wpm/affp;
/* calculate form factor: ffe */
   isign = 0;
   cepois13((float complex *)dcu,(float complex *)cus,isign,ffe,ax,
            affp,wpm,ci,&wf,nx,nxeh,nxh);

/* initialize transverse electric field */
   for (j = 0; j < 2*nxe; j++) {
      cus[j] = 0.0;
   }

/* * * * start main iteration loop * * * */
 
L500: if (nloop <= ntime)
         goto L2000;
/*    printf("ntime = %i\n",ntime); */
 
/* deposit current with standard procedure: updates cue */
      dtimer(&dtime,&itime,-1);
      for (j = 0; j < 2*nxe; j++) {
         cue[j] = 0.0;
      }
      cgjpost1l(part,cue,qme,zero,np,idimp,nx,nxe,ipbc);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tdjpost += time;

/* deposit charge with standard procedure: updates qe */
      dtimer(&dtime,&itime,-1);
      for (j = 0; j < nxe; j++) {
         qe[j] = 0.0;
      }
      cgpost1l(part,qe,qme,np,idimp,nxe);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tdpost += time;

/* add guard cells with standard procedure: updates qe, cue */
      dtimer(&dtime,&itime,-1);
      caguard1l(qe,nx,nxe);
      cacguard1l(cue,nx,nxe);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tguard += time;

/* transform charge to fourier space with standard procedure: */
/* updates qe, gxe */
      dtimer(&dtime,&itime,-1);
      isign = -1;
      cfft1rxx((float complex *)qe,(float complex *)gxe,isign,mixup,sct,
               indx,nxe,nxh);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tfft += time;

/* calculate longitudinal force/charge in fourier space with standard */
/* procedure: updates fxe, we                                         */
      dtimer(&dtime,&itime,-1);
      isign = -1;
      cpois1((float complex *)qe,(float complex *)fxe,isign,ffc,ax,affp,
             &we,nx);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tfield += time;

/* transform longitudinal electric force to real space with standard */
/* procedure: updates fxe, gxe                                       */
      dtimer(&dtime,&itime,-1);
      isign = 1;
      cfft1rxx((float complex *)fxe,(float complex *)gxe,isign,mixup,
               sct,indx,nxe,nxh);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tfft += time;

/* transform current to fourier space with standard procedure: */
/* updates cue, gyze                                           */
      dtimer(&dtime,&itime,-1);
      isign = -1;
      cfft1r2x((float complex *)cue,(float complex *)gyze,isign,mixup,
               sct,indx,nxe,nxh);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tfft += time;

/* calculate magnetic field in fourier space with standard procedure: */
/* updates byze, wm                                                   */
      dtimer(&dtime,&itime,-1);
      cbbpois13((float complex *)cue,(float complex *)byze,ffc,ci,&wm,
                nx,nxeh,nxh);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tfield += time;

/* transform magnetic force to real space with standard procedure: */
/* updates byze, gyze */
      dtimer(&dtime,&itime,-1);
      isign = 1;
      cfft1r2x((float complex *)byze,(float complex *)gyze,isign,mixup,
               sct,indx,nxe,nxh);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tfft += time;

/* add constant to magnetic field with standard procedure: updates byze */
      dtimer(&dtime,&itime,-1);
      cbaddext1(byze,omy,omz,nx,nxe);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tfield += time;

/* copy guard cells with standard procedure: updates fxe, byze */
      dtimer(&dtime,&itime,-1);
      cdguard1l(fxe,nx,nxe);
      ccguard1l(byze,nx,nxe);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tguard += time;

/* add longitudinal and old transverse electric fields with standard */
/* procedure: updates exyze                                          */
      dtimer(&dtime,&itime,-1);
      caddvrfield13(exyze,cus,fxe,nxe);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tfield += time;

/* deposit electron acceleration density and momentum flux with */
/* standard procedure: updates dcu, amu                         */
      dtimer(&dtime,&itime,-1);
      for (j = 0; j < 2*nxe; j++) {
         dcu[j] = 0.0;
         amu[j] = 0.0;
      }
      cgdjpost1l(part,exyze,byze,dcu,amu,omx,qme,qbme,dt,idimp,np,nxe);
/* add old scaled electric field with standard procedure: updates dcu */
      cascfguard1l(dcu,cus,q2m0,nx,nxe);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tdcjpost += time;

/* add guard cells with standard procedure: updates dcu, amu */
      dtimer(&dtime,&itime,-1);
      cacguard1l(dcu,nx,nxe);
      cacguard1l(amu,nx,nxe);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tguard += time;

/* transform acceleration density and momentum flux to fourier space */
/* with standard procedure: updates dcu, amu, gyze                   */
      dtimer(&dtime,&itime,-1);
      isign = -1;
      cfft1r2x((float complex *)dcu,(float complex *)gyze,isign,mixup,
               sct,indx,nxe,nxh);
      cfft1r2x((float complex *)amu,(float complex *)gyze,isign,mixup,
               sct,indx,nxe,nxh);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tfft += time;

/* take transverse part of time derivative of current with standard */
/* procedure: updates dcu                                           */
      dtimer(&dtime,&itime,-1);
      cadcuperp13((float complex *)dcu,(float complex *)amu,nx,nxeh);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tfield += time;

/* calculate transverse electric field with standard procedure: */
/* updates cus, wf                                              */
      dtimer(&dtime,&itime,-1);
      isign = -1;
      cepois13((float complex *)dcu,(float complex *)cus,isign,ffe,ax,
               affp,wpm,ci,&wf,nx,nxeh,nxh);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tfield += time;

/* transform transverse electric field to real space with standard */
/* procedure: updates cus, gyze                                    */
      dtimer(&dtime,&itime,-1);
      isign = 1;
      cfft1r2x((float complex *)cus,(float complex *)gyze,isign,mixup,
               sct,indx,nxe,nxh);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tfft += time;

/* copy guard cells with standard procedure: updates cus */
      dtimer(&dtime,&itime,-1);
      ccguard1l(cus,nx,nxe);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tguard += time;

/* add longitudinal and transverse electric fields with standard */
/* procedure: exyze = cus + fxe, updates exyze                   */
/* cus needs to be retained for next time step */
      dtimer(&dtime,&itime,-1);
      caddvrfield13(exyze,cus,fxe,nxe);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tfield += time;

/* inner iteration loop */
      for (k = 0; k < ndc; k++) {

/* deposit electron current and acceleration density and momentum flux */
/* with standard procedure: updates cue, dcu, amu                      */
         dtimer(&dtime,&itime,-1);
         for (j = 0; j < 2*nxe; j++) {
            cue[j] = 0.0;
            dcu[j] = 0.0;
            amu[j] = 0.0;
         }
         cgdcjpost1l(part,exyze,byze,cue,dcu,amu,omx,qme,qbme,dt,idimp,
                    np,nxe);
/* add scaled electric field with standard procedure: updates dcu */
         cascfguard1l(dcu,cus,q2m0,nx,nxe);
         dtimer(&dtime,&itime,1);
         time = (float) dtime;
         tdcjpost += time;

/* add guard cells for current, acceleration density, and momentum flux */
/* with standard procedure: updates cue, dcu, amu                       */
         dtimer(&dtime,&itime,-1);
         cacguard1l(cue,nx,nxe);
         cacguard1l(dcu,nx,nxe);
         cacguard1l(amu,nx,nxe);
         dtimer(&dtime,&itime,1);
         time = (float) dtime;
         tguard += time;

/* transform current to fourier space with standard procedure: */
/* update cue, gyze                                            */
         dtimer(&dtime,&itime,-1);
         isign = -1;
         cfft1r2x((float complex *)cue,(float complex *)gyze,isign,mixup,
                  sct,indx,nxe,nxh);
         dtimer(&dtime,&itime,1);
         time = (float) dtime;
         tfft += time;

/* calculate magnetic field in fourier space with standard procedure: */
/* updates byze, wm                                                  */
         dtimer(&dtime,&itime,-1);
         cbbpois13((float complex *)cue,(float complex *)byze,ffc,ci,&wm,
                   nx,nxeh,nxh);
         dtimer(&dtime,&itime,1);
         time = (float) dtime;
         tfield += time;

/* transform magnetic force to real space with standard procedure: */
/* updates byze, gyze                                             */
         dtimer(&dtime,&itime,-1);
         isign = 1;
         cfft1r2x((float complex *)byze,(float complex *)gyze,isign,
                  mixup,sct,indx,nxe,nxh);
         dtimer(&dtime,&itime,1);
         time = (float) dtime;
         tfft += time;

/* add constant to magnetic field with standard procedure: updates byze */
         dtimer(&dtime,&itime,-1);
         cbaddext1(byze,omy,omz,nx,nxe);
         dtimer(&dtime,&itime,1);
         time = (float) dtime;
         tfield += time;

/* transform acceleration density and momentum flux to fourier space */
/* with standard procedure: updates dcu amu, gyze                    */
         dtimer(&dtime,&itime,-1);
         isign = -1;
         cfft1r2x((float complex *)dcu,(float complex *)gyze,isign,mixup,
                  sct,indx,nxe,nxh);
         cfft1r2x((float complex *)amu,(float complex *)gyze,isign,mixup,
                  sct,indx,nxe,nxh);
         dtimer(&dtime,&itime,1);
         time = (float) dtime;
         tfft += time;
 
/* take transverse part of time derivative of current with standard */
/* procedure: updates dcu                                           */
         dtimer(&dtime,&itime,-1);
         cadcuperp13((float complex *)dcu,(float complex *)amu,nx,nxeh);
         dtimer(&dtime,&itime,1);
         time = (float) dtime;
         tfield += time;
 
/* calculate convective part of transverse electric field with standard */
/* procedure: updates cus, wf                                           */
         dtimer(&dtime,&itime,-1);
         isign = -1;
         cepois13((float complex *)dcu,(float complex *)cus,isign,ffe,ax,
                   affp,wpm,ci,&wf,nx,nxeh,nxh);
         dtimer(&dtime,&itime,1);
         time = (float) dtime;
         tfield += time;
 
/* transform transverse electric field to real space with standard */
/* procedure: updates cus, gyze                                    */
         dtimer(&dtime,&itime,-1);
         isign = 1;
         cfft1r2x((float complex *)cus,(float complex *)gyze,isign,mixup,
                  sct,indx,nxe,nxh);
         dtimer(&dtime,&itime,1);
         time = (float) dtime;
         tfft += time;
 
/* copy guard cells with standard procedure: updates byze, cus */
         dtimer(&dtime,&itime,-1);
         ccguard1l(byze,nx,nxe);
         ccguard1l(cus,nx,nxe);
         dtimer(&dtime,&itime,1);
         time = (float) dtime;
         tguard += time;
 
/* add longitudinal and transverse electric fields with standard */
/* procedure: exyze = cus + fxe, updates exyze                   */
/* cus needs to be retained for next time step                   */
         dtimer(&dtime,&itime,-1);
         caddvrfield13(exyze,cus,fxe,nxe);
         dtimer(&dtime,&itime,1);
         time = (float) dtime;
         tfield += time;

      }

/* push particles with standard procedure: updates part, wke */
      wke = 0.0;
      dtimer(&dtime,&itime,-1);
      cgbpush13l(part,exyze,byze,omx,qbme,dt,dt,&wke,idimp,np,nx,nxe,
                 ipbc);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tpush += time;

/* sort particles by cell for standard procedure */
      if (sortime > 0) {
         if (ntime%sortime==0) {
            dtimer(&dtime,&itime,-1);
            cdsortp1xl(part,part2,npic,idimp,np,nx1);
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
         wt = we + wm;
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

   printf("ntime, ndc = %i,%i\n",ntime,ndc);
   wt = we + wm;
   printf("Final Total Field, Kinetic and Total Energies:\n");
   printf("%e %e %e\n",wt,wke,wke+wt);
   printf("Final Electrostatic, Transverse Electric and Magnetic Field \
Energies:\n");
   printf("%e %e %e\n",we,wf,wm);

   printf("\n");
   printf("deposit time = %f\n",tdpost);
   printf("current deposit time = %f\n",tdjpost);
   printf("current derivative deposit time = %f\n",tdcjpost);
   tdpost += tdjpost + tdcjpost;
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
