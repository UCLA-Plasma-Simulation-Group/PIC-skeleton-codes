/*---------------------------------------------------------------------*/
/* Skeleton 1D Electrostatic OpenMP PIC code */
/* written by Viktor K. Decyk, UCLA */
#include <stdlib.h>
#include <stdio.h>
#include <complex.h>
#include <sys/time.h>
#include "mpush1.h"
#include "omplib.h"

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
   int idimp = 2, ipbc = 1;
/* wke/we/wt = particle kinetic/electric field/total energy */
   float wke = 0.0, we = 0.0, wt = 0.0;
/* mx = number of grids in x in sorting tiles */
   int mx = 32;
/* xtras = fraction of extra particles needed for particle management */
   float xtras = 0.2;
/* declare scalars for standard code */
   int j;
   int np, nx, nxh, nxe;
   int mx1, ntime, nloop, isign;
   float qbme, affp;

/* declare scalars for OpenMP code */
   int nppmx, nppmx0, ntmax, npbmx, irc;
   int nvp;

/* declare arrays for standard code: */
/* part = particle array */
   float *part = NULL;
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

/* declare arrays for OpenMP (tiled) code: */
/* ppart = tiled particle array */
/* ppbuff = buffer array for reordering tiled particle array */
   float *ppart = NULL, *ppbuff = NULL;
/* kpic = number of particles in each tile */
   int *kpic = NULL;
/* ncl = number of particles departing tile in each direction */
   int *ncl = NULL;
/* ihole = location/destination of each particle departing tile */
   int *ihole = NULL;

/* declare and initialize timing data */
   float time;
   struct timeval itime;
   float tdpost = 0.0, tguard = 0.0, tfft = 0.0, tfield = 0.0;
   float tpush = 0.0, tsort = 0.0;
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
/* nx = number of grid points in x direction */
   np = npx; nx = 1L<<indx; nxh = nx/2;
   nxe = nx + 2;
/* mx1 = number of tiles in x direction */
   mx1 = (nx - 1)/mx + 1;
/* nloop = number of time steps in simulation */
/* ntime = current time step */
   nloop = tend/dt + .0001; ntime = 0;
   qbme = qme;
   affp = (float) nx/(float ) np;

/* allocate data for standard code */
   part = (float *) malloc(idimp*np*sizeof(float));
   qe = (float *) malloc(nxe*sizeof(float));
   fxe = (float *) malloc(nxe*sizeof(float));
   ffc = (float complex *) malloc(nxh*sizeof(float complex));
   mixup = (int *) malloc(nxh*sizeof(int));
   sct = (float complex *) malloc(nxh*sizeof(float complex));
   kpic = (int *) malloc(mx1*sizeof(int));

/* prepare fft tables */
   cwfft1rinit(mixup,sct,indx,nxh);
/* calculate form factors */
   isign = 0;
   cpois1((float complex *)qe,(float complex *)fxe,isign,ffc,ax,affp,
          &we,nx);
/* initialize electrons */
   cdistr1(part,vtx,vx0,npx,idimp,np,nx,ipbc);

/* find number of particles in each of mx, tiles: updates kpic, nppmx */
   cdblkp1l(part,kpic,&nppmx,idimp,np,mx,mx1,&irc);
   if (irc != 0) { 
      printf("cdblkp1l error, irc=%d\n",irc);
      exit(1);
   }
/* allocate vector particle data */
   nppmx0 = (1.0 + xtras)*nppmx;
   ntmax = xtras*nppmx;
   npbmx = xtras*nppmx;
   ppart = (float *) malloc(idimp*nppmx0*mx1*sizeof(float));
   ppbuff = (float *) malloc(idimp*npbmx*mx1*sizeof(float));
   ncl = (int *) malloc(2*mx1*sizeof(int));
   ihole = (int *) malloc(2*(ntmax+1)*mx1*sizeof(int));
/* copy ordered particle data for OpenMP: updates ppart and kpic */
   cppmovin1l(part,ppart,kpic,nppmx0,idimp,np,mx,mx1,&irc);
   if (irc != 0) { 
      printf("cppmovin1l overflow error, irc=%d\n",irc);
      exit(1);
   }
/* sanity check */
   cppcheck1l(ppart,kpic,idimp,nppmx0,nx,mx,mx1,&irc);
   if (irc != 0) {
      printf("%d,cppcheck1l error: irc=%d\n",ntime,irc);
      exit(1);
   }

/* * * * start main iteration loop * * * */
 
L500: if (nloop <= ntime)
         goto L2000;
/*    printf("ntime = %i\n",ntime); */
 
/* deposit charge with OpenMP: updates qe */
      dtimer(&dtime,&itime,-1);
      for (j = 0; j < nxe; j++) {
         qe[j] = 0.0;
      }
      cgppost1l(ppart,qe,kpic,qme,nppmx0,idimp,mx,nxe,mx1);
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

/* push particles with OpenMP: updates part, wke */
      wke = 0.0;
      dtimer(&dtime,&itime,-1);
/* updates ppart, wke */
/*    cgppush1l(ppart,fxe,kpic,qbme,dt,&wke,idimp,nppmx0,nx,mx,nxe, */
/*              mx1,ipbc);                                          */
/* updates ppart, ncl, ihole, wke, irc */
      cgppushf1l(ppart,fxe,kpic,ncl,ihole,qbme,dt,&wke,idimp,nppmx0,
                 nx,mx,nxe,mx1,ntmax,&irc);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tpush += time;
      if (irc != 0) {
         printf("cgppushf1l error: irc=%d\n",irc);
         exit(1);
      }

/* reorder particles by tile with OpenMP: */
      dtimer(&dtime,&itime,-1);
/* updates ppart, ppbuff, kpic, ncl, ihole, and irc */
/*    cpporder1l(ppart,ppbuff,kpic,ncl,ihole,idimp,nppmx0,nx,mx,mx1, */
/*               npbmx,ntmax,&irc);                                  */
/* updates ppart, ppbuff, kpic, ncl, and irc */
      cpporderf1l(ppart,ppbuff,kpic,ncl,ihole,idimp,nppmx0,mx1,npbmx,
                  ntmax,&irc);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tsort += time;
      if (irc != 0) {
         printf("cpporderf1l error: ntmax, irc=%d,%d\n",ntmax,irc);
         exit(1);
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
   printf("\n");

   return 0;
}
