/*---------------------------------------------------------------------*/
/* Skeleton 1-2/2D Electromagnetic OpenMP PIC code */
/* written by Viktor K. Decyk, UCLA */
#include <stdlib.h>
#include <stdio.h>
#include <complex.h>
#include <sys/time.h>
#include "mbpush1.h"
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
   float tend = 10.0, dt = 0.05, qme = -1.0;
/* vtx/vty = thermal velocity of electrons in x/y direction */
/* vx0/vy0 = drift velocity of electrons in x/y direction */
   float vtx = 1.0, vty = 1.0, vx0 = 0.0, vy0 = 0.0;
/* vtx/vz0 = thermal/drift velocity of electrons in z direction */
   float vtz = 1.0, vz0 = 0.0;
/* omx = magnetic field electron cyclotron frequency in x */
   float omx = 0.0;
/* ax = smoothed particle size in x direction */
/* ci = reciprocal of velocity of light */
   float ax = .912871, ci = 0.1;
/* idimp = number of particle coordinates = 4 */
/* ipbc = particle boundary condition: 1 = periodic */
/* relativity = (no,yes) = (0,1) = relativity is used */
   int idimp = 4, ipbc = 1, relativity = 1;
   float wke = 0.0, we = 0.0, wf = 0.0, wm = 0.0, wt = 0.0;
/* mx = number of grids in x in sorting tiles */
   int mx = 32;
/* xtras = fraction of extra particles needed for particle management */
   float xtras = 0.2;
/* declare scalars for standard code */
   int j;
   int np, nx, nxh, nxe, nxeh;
   int mx1, ntime, nloop, isign;
   float qbme, affp, dth;

/* declare scalars for OpenMP code */
   int nppmx, nppmx0, ntmax, npbmx, irc;
   int nvp;

/* declare arrays for standard code: */
/* part, part2 = particle arrays */
   float *part = NULL;
/* qe = electron charge density with guard cells */
/* fxe = smoothed longitudinal electric field with guard cells */
   float *qe = NULL, *fxe = NULL;
/* cue = electron current density with guard cells */
/* fxyze/byze = smoothed electric/magnetic field with guard cells */
   float *cue = NULL, *fxyze = NULL, *byze = NULL;
/* eyz/byz = transverse electric/magnetic field in fourier space */
   float complex *eyz = NULL, *byz = NULL;
/* ffc = form factor array for poisson solver */
/* sct = sine/cosine table for FFT */
   float complex *ffc = NULL, *sct = NULL;
/* mixup = bit reverse table for FFT */
   int *mixup = NULL;
/* gxyze = scratch array for fft */
   float *gxyze = NULL;

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
/* nx = number of grid points in x direction */
   np = npx; nx = 1L<<indx; nxh = nx/2;
   nxe = nx + 2; nxeh = nxe/2;
/* mx1 = number of tiles in x direction */
   mx1 = (nx - 1)/mx + 1;
/* nloop = number of time steps in simulation */
/* ntime = current time step */
   nloop = tend/dt + .0001; ntime = 0;
   qbme = qme;
   affp = (float) nx/(float ) np;
   dth = 0.0;

/* allocate data for standard code */
   part = (float *) malloc(idimp*np*sizeof(float));
   qe = (float *) malloc(nxe*sizeof(float));
   fxe = (float *) malloc(nxe*sizeof(float));
   fxyze = (float *) malloc(3*nxe*sizeof(float));
   cue = (float *) malloc(2*nxe*sizeof(float));
   byze = (float *) malloc(2*nxe*sizeof(float));
   eyz = (float complex *) malloc(2*nxeh*sizeof(float complex));
   byz = (float complex *) malloc(2*nxeh*sizeof(float complex));
   ffc = (float complex *) malloc(nxh*sizeof(float complex));
   mixup = (int *) malloc(nxh*sizeof(int));
   sct = (float complex *) malloc(nxh*sizeof(float complex));
   kpic = (int *) malloc(mx1*sizeof(int));
   gxyze = (float *) malloc(3*nxe*sizeof(float));

/* prepare fft tables */
   cwfft1rinit(mixup,sct,indx,nxh);
/* calculate form factors */
   isign = 0;
   cpois1((float complex *)qe,(float complex *)fxe,isign,ffc,ax,affp,
          &we,nx);
/* initialize electrons */
   cdistr1h(part,vtx,vty,vtz,vx0,vy0,vz0,npx,idimp,np,nx,ipbc);

/* initialize transverse electromagnetic fields */
   for (j = 0; j < 2*nxeh; j++) {
      eyz[j] = 0.0 + 0.0*_Complex_I;
      byz[j] = 0.0 + 0.0*_Complex_I;
   }

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

   if (dt > 0.64*ci) {
      printf("Warning: Courant condition may be exceeded!\n");
   }

/* * * * start main iteration loop * * * */
 
L500: if (nloop <= ntime)
         goto L2000;
/*    printf("ntime = %i\n",ntime); */

/* deposit current with OpenMP: */
      dtimer(&dtime,&itime,-1);
      for (j = 0; j < 2*nxe; j++) {
         cue[j] = 0.0;
      }
      if (relativity==1) {
/* updates ppart, cue */
/*       cgrjppost1l(ppart,cue,kpic,qme,dth,ci,nppmx0,idimp,nx,mx,nxe, */
/*                   mx1,ipbc);                                        */
/* updates ppart, cue, ncl, ihole, irc */
         cgrjppostf1l(ppart,cue,kpic,ncl,ihole,qme,dth,ci,nppmx0,idimp,
                      nx,mx,nxe,mx1,ntmax,&irc);
      }
      else {
/* updates ppart, cue */
/*       cgjppost1l(ppart,cue,kpic,qme,dth,nppmx0,idimp,nx,mx,nxe,mx1, */
/*                  ipbc);                                             */
/* updates ppart, cue, ncl, ihole, irc */
         cgjppostf1l(ppart,cue,kpic,ncl,ihole,qme,dth,nppmx0,idimp,nx,
                      mx,nxe,mx1,ntmax,&irc);
      }
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tdjpost += time;
      if (irc != 0) {
         if (relativity==1) {
            printf("cgrjppostf1l error: irc=%d\n",irc);
         }
         else {
            printf("cgjppostf1l error: irc=%d\n",irc);
         }
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

/* deposit charge with OpenMP: updates qe */
      dtimer(&dtime,&itime,-1);
      for (j = 0; j < nxe; j++) {
         qe[j] = 0.0;
      }
      cgppost1l(ppart,qe,kpic,qme,nppmx0,idimp,mx,nxe,mx1);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tdpost += time;

/* add guard cells with standard procedure: updates cue, qe */
      dtimer(&dtime,&itime,-1);
      cacguard1l(cue,nx,nxe);
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

/* transform current to fourier space with standard procedure: */
/*updates cue, byze */
      dtimer(&dtime,&itime,-1);
      isign = -1;
      cfft1r2x((float complex *)cue,(float complex *)byze,isign,mixup,
               sct,indx,nxe,nxh);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tfft += time;

/* calculate electromagnetic fields in fourier space with standard */
/* procedure: updates eyz, byz                                     */
      dtimer(&dtime,&itime,-1);
      if (ntime==0) {
         cibpois13((float complex *)cue,byz,ffc,ci,&wm,nx,nxeh,nxh);
         wf = 0.0;
         dth = 0.5*dt;
      }
      else {
         cmaxwel1(eyz,byz,(float complex *)cue,ffc,ci,dt,&wf,&wm,nx,
                  nxeh,nxh);
      }
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tfield += time;

/* calculate force/charge in fourier space with standard procedure: */
/* updates fxe, we */
      dtimer(&dtime,&itime,-1);
      isign = -1;
      cpois1((float complex *)qe,(float complex *)fxe,isign,ffc,ax,affp,
             &we,nx);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tfield += time;

/* add longitudinal and transverse electric fields with standard */
/* procedure: updates fxyze */
      dtimer(&dtime,&itime,-1);
      isign = 1;
      cemfield1((float complex *)fxyze,(float complex *)fxe,eyz,ffc,nx,
                nxeh,nxh);
/* copy magnetic field with standard procedure: updates byze */
      isign = -1;
      cbmfield1((float complex *)byze,byz,ffc,nx,nxeh,nxh);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tfield += time;

/* transform electric force to real space with standard procedure: */
/* updates fxyze, gxyze */
      dtimer(&dtime,&itime,-1);
      isign = 1;
      cfft1r3x((float complex *)fxyze,(float complex *)gxyze,isign,
                mixup,sct,indx,nxe,nxh);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tfft += time;

/* transform magnetic force to real space with standard procedure: */
/* updates byze, cue */
      dtimer(&dtime,&itime,-1);
      isign = 1;
      cfft1r2x((float complex *)byze,(float complex *)cue,isign,mixup,
               sct,indx,nxe,nxh);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tfft += time;

/* copy guard cells with standard procedure: updates fxyze, byze */
      dtimer(&dtime,&itime,-1);
      cbguard1l(fxyze,nx,nxe);
      ccguard1l(byze,nx,nxe);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tguard += time;

/* push particles with OpenMP: */
      wke = 0.0;
      dtimer(&dtime,&itime,-1);
      if (relativity==1) {
/* updates ppart, wke */
/*       cgrbppush13l(ppart,fxyze,byze,kpic,omx,qbme,dt,dth,ci,&wke, */
/*                    idimp,nppmx0,nx,mx,nxe,mx1,ipbc);              */
/* updates ppart, ncl, ihole, wke, irc */
         cgrbppushf13l(ppart,fxyze,byze,kpic,ncl,ihole,omx,qbme,dt,dth,
                       ci,&wke,idimp,nppmx0,nx,mx,nxe,mx1,ntmax,&irc);
      }
      else {
/* updates ppart, wke */
/*       cgbppush13l(ppart,fxyze,byze,kpic,omx,qbme,dt,dth,&wke,idimp, */
/*                   nppmx0,nx,mx,nxe,mx1,ipbc);                       */
/* updates ppart, ncl, ihole, wke, irc */
         cgbppushf13l(ppart,fxyze,byze,kpic,ncl,ihole,omx,qbme,dt,dth,
                      &wke,idimp,nppmx0,nx,mx,nxe,mx1,ntmax,&irc);
      }
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tpush += time;
      if (irc != 0) {
         if (relativity==1) {
            printf("cgrbppushf13l error: irc=%d\n",irc);
         }
         else {
            printf("cgbppushf13l error: irc=%d\n",irc);
         }
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
