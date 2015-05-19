/*---------------------------------------------------------------------*/
/* Skeleton 1D Electrostatic PIC code */
/* written by Viktor K. Decyk, UCLA */
#include <stdlib.h>
#include <stdio.h>
#include <complex.h>
#include <sys/time.h>
#include "push1.h"

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
/* vtx = thermal velocity of electrons in x direction */
/* vx0 = drift velocity of electrons in x direction */
   float vtx = 1.0, vx0 = 0.0;
/* ax = smoothed particle size in x direction */
   float ax = .912871;
/* idimp = number of particle coordinates = 2 */
/* ipbc = particle boundary condition: 1 = periodic */
/* sortime = number of time steps between standard electron sorting */
   int idimp = 2, ipbc = 1, sortime = 50;
/* wke/we/wt = particle kinetic/electric field/total energy */
   float wke = 0.0, we = 0.0, wt = 0.0;
/* declare scalars for standard code */
   int j;
   int np, nx, nxh, nxe;
   int nx1, ntime, nloop, isign;
   float qbme, affp;

/* declare arrays for standard code: */
/* part, part2 = particle arrays */
   float *part = NULL, *part2 = NULL, *tpart = NULL;
/* qe = electron charge density with guard cells */
   float *qe = NULL;
/* fxe = smoothed electric field with guard cells */
   float *fxe = NULL;
/* ffc = form factor array for poisson solver */
   float complex *ffc = NULL;
/* mixup = bit reverse table for FFT */
   int *mixup = NULL;
/* sct = sine/cosine table for FFT */
   float complex *sct = NULL;
/* npic = scratch array for reordering particles */
   int *npic = NULL;
  
/* declare and initialize timing data */
   float time;
   struct timeval itime;
   float tdpost = 0.0, tguard = 0.0, tfft = 0.0, tfield = 0.0;
   float tpush = 0.0, tsort = 0.0;
   double dtime;

/* initialize scalars for standard code */
/* np = total number of particles in simulation */
/* nx = number of grid points in x direction */
   np = npx; nx = 1L<<indx; nxh = nx/2;
   nxe = nx + 2; nx1 = nx + 1;
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
   ffc = (float complex *) malloc(nxh*sizeof(float complex));
   mixup = (int *) malloc(nxh*sizeof(int));
   sct = (float complex *) malloc(nxh*sizeof(float complex));
   npic = (int *) malloc(nx1*sizeof(int));

/* prepare fft tables */
   cwfft1rinit(mixup,sct,indx,nxh);
/* calculate form factors */
   isign = 0;
   cpois1((float complex *)qe,(float complex *)fxe,isign,ffc,ax,affp,
          &we,nx);
/* initialize electrons */
   cdistr1(part,vtx,vx0,npx,idimp,np,nx,ipbc);

/* * * * start main iteration loop * * * */
 
L500: if (nloop <= ntime)
         goto L2000;
/*    printf("ntime = %i\n",ntime); */
 
/* deposit charge with standard procedure: updates qe */
      dtimer(&dtime,&itime,-1);
      for (j = 0; j < nxe; j++) {
         qe[j] = 0.0;
      }
      cgpost1l(part,qe,qme,np,idimp,nxe);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tdpost += time;

/* add guard cells with standard procedure: updates qe */
      dtimer(&dtime,&itime,-1);
      caguard1l(qe,nx,nxe);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tguard += time;

/* transform charge to fourier space with standard procedure: */
/* updates qe, fxe */
      dtimer(&dtime,&itime,-1);
      isign = -1;
      cfft1rxx((float complex *)qe,(float complex *)fxe,isign,mixup,sct,
               indx,nxe,nxh);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tfft += time;

/* calculate force/charge in fourier space with standard procedure: */
/* updates fxe, we */
      dtimer(&dtime,&itime,-1);
      isign = -1;
      cpois1((float complex *)qe,(float complex *)fxe,isign,ffc,ax,affp,
             &we,nx);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tfield += time;

/* transform force to real space with standard procedure: */
/* updates fxe, qe */
      dtimer(&dtime,&itime,-1);
      isign = 1;
      cfft1rxx((float complex *)fxe,(float complex *)qe,isign,mixup,sct,
               indx,nxe,nxh);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tfft += time;

/* copy guard cells with standard procedure: updates fxe */
      dtimer(&dtime,&itime,-1);
      ccguard1l(fxe,nx,nxe);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tguard += time;

/* push particles with standard procedure: updates part, wke */
      wke = 0.0;
      dtimer(&dtime,&itime,-1);
      cgpush1l(part,fxe,qbme,dt,&wke,idimp,np,nx,nxe,ipbc);
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
         printf("Initial Field, Kinetic and Total Energies:\n");
         printf("%e %e %e\n",we,wke,wke+we);
      }
      ntime += 1;
      goto L500;
L2000:

/* * * * end main iteration loop * * * */

   printf("ntime = %i\n",ntime);
   printf("Final Field, Kinetic and Total Energies:\n");
   printf("%e %e %e\n",we,wke,wke+we);

   printf("\n");
   printf("deposit time = %f\n",tdpost);
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
