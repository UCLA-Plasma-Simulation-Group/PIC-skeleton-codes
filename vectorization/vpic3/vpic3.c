/*---------------------------------------------------------------------*/
/* Skeleton 3D Electrostatic Vector PIC code */
/* written by Viktor K. Decyk, UCLA and Ricardo Fonseca, ISCTE */
#include <stdlib.h>
#include <stdio.h>
#include <complex.h>
#include <sys/time.h>
#include "vpush3.h"
#include "avx512lib3.h"
#include "kncpush3.h"

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
   float tend = 10.0, dt = 0.1, qme = -1.0;
/* vtx/vty/vtz = thermal velocity of electrons in x/y/z direction */
   float vtx = 1.0, vty = 1.0, vtz = 1.0;
/* vx0/vy0/vz0 = drift velocity of electrons in x/y/z direction */
   float vx0 = 0.0, vy0 = 0.0, vz0 = 0.0;
/* ax/ay/az = smoothed particle size in x/y/z direction */
   float ax = .912871, ay = .912871, az = .912871;
/* idimp = number of particle coordinates = 6 */
/* ipbc = particle boundary condition: 1 = periodic */
/* sortime = number of time steps between standard electron sorting */
   int idimp = 6, ipbc = 1, sortime = 20;
/* wke/we/wt = particle kinetic/electric field/total energy */
   float wke = 0.0, we = 0.0, wt = 0.0;
/* kvec = (1,2) = run (autovector,KNC) version */
   int kvec = 1;
   
/* declare scalars for standard code */
   int j;
   int np, nx, ny, nz, nxh, nyh, nzh, nxe, nye, nze, nxeh;
   int nxyzh, nxhyz, npe, ny1, nyz1, ntime, nloop, isign;
   int irc = 0;
   float qbme, affp;

/* declare arrays for standard code: */
/* partt, partt2 = transposed particle arrays */
   float *partt = NULL, *partt2 = NULL, *tpartt = NULL;
/* qe = electron charge density with guard cells */
   float *qe = NULL;
/* fxyze = smoothed electric field with guard cells */
   float *fxyze = NULL;
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
/* nx/ny/nz = number of grid points in x/y direction */
   np = npx*npy*npz; nx = 1L<<indx; ny = 1L<<indy; nz = 1L<<indz;
   nxh = nx/2; nyh = 1 > ny/2 ? 1 : ny/2; nzh = 1 > nz/2 ? 1 : nz/2;
   nxe = nx + 2; nye = ny + 1; nze = nz + 1; nxeh = nxe/2;
   nxyzh = (nx > ny ? nx : ny); nxyzh = (nxyzh > nz ? nxyzh : nz)/2;
   nxhyz = nxh > ny ? nxh : ny; nxhyz = nxhyz > nz ? nxhyz : nz;
   ny1 = ny + 1; nyz1 = ny1*(nz + 1);
/* nloop = number of time steps in simulation */
/* ntime = current time step */
   nloop = tend/dt + .0001; ntime = 0;
   qbme = qme;
   affp = ((float) nx)*((float) ny)*((float) nz)/(float ) np;

/* allocate data for standard code */
   mixup = (int *) malloc(nxhyz*sizeof(int));
   sct = (float complex *) malloc(nxyzh*sizeof(float complex));

/* align memory for avx512 */
   npe = 16*((np - 1)/16 + 1);
   nxe = 16*((nxe - 1)/16 + 1);
   nxeh = nxe/2;
   avx512_fallocate(&partt,npe*idimp,&irc);
   if (sortime > 0)
      avx512_fallocate(&partt2,npe*idimp,&irc);
   avx512_fallocate(&qe,nxe*nye*nze,&irc);
   avx512_fallocate(&fxyze,ndim*nxe*nye*nze,&irc);
   avx512_callocate(&ffc,nxh*nyh*nzh,&irc);
   avx512_iallocate(&npic,nyz1,&irc);
   if (irc != 0) {
      printf("aligned allocation error: irc = %d\n",irc);
   }

/* prepare fft tables */
   cwfft3rinit(mixup,sct,indx,indy,indz,nxhyz,nxyzh);
/* calculate form factors */
   isign = 0;
   cvpois33((float complex *)qe,(float complex *)fxyze,isign,ffc,ax,ay,
            az,affp,&we,nx,ny,nz,nxeh,nye,nze,nxh,nyh,nzh);
/* initialize electrons */
   cdistr3t(partt,vtx,vty,vtz,vx0,vy0,vz0,npx,npy,npz,idimp,np,nx,ny,nz,
            ipbc);

/* * * * start main iteration loop * * * */
 
L500: if (nloop <= ntime)
         goto L2000;
/*    printf("ntime = %i\n",ntime); */
 
/* deposit charge with standard procedure: updates qe */
      dtimer(&dtime,&itime,-1);
      for (j = 0; j < nxe*nye; j++) {
         qe[j] = 0.0;
      }
      if (kvec==1)
         cgpost3lt(partt,qe,qme,np,npe,idimp,nxe,nye,nze);
/*       cvgpost3lt(partt,qe,qme,np,npe,idimp,nxe,nye,nze); */
/* KNC function */
      else if (kvec==2)
         cknc2gpost3lt(partt,qe,qme,np,npe,idimp,nxe,nye,nze);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tdpost += time;

/* add guard cells with standard procedure: updates qe */
      dtimer(&dtime,&itime,-1);
      if (kvec==1)
         caguard3l(qe,nx,ny,nz,nxe,nye,nze);
/* KNC function */
      else if (kvec==2)
         ckncaguard3l(qe,nx,ny,nz,nxe,nye,nze);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tguard += time;

/* transform charge to fourier space with standard procedure: updates qe */
      dtimer(&dtime,&itime,-1);
      isign = -1;
      if (kvec==1)
         cwfft3rvx((float complex *)qe,isign,mixup,sct,indx,indy,indz,
                   nxeh,nye,nze,nxhyz,nxyzh);
/* KNC function */
      else if (kvec==2)
         ckncwfft3rvx((float complex *)qe,isign,mixup,sct,indx,indy,
                      indz,nxeh,nye,nze,nxhyz,nxyzh);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tfft += time;

/* calculate force/charge in fourier space with standard procedure: */
/* updates fxyze, we                                                */
      dtimer(&dtime,&itime,-1);
      isign = -1;
      if (kvec==1)
         cvpois33((float complex *)qe,(float complex *)fxyze,isign,ffc,
                   ax,ay,az,affp,&we,nx,ny,nz,nxeh,nye,nze,nxh,nyh,nzh);
/* KNC function */
      else if (kvec==2)
         ckncpois33((float complex *)qe,(float complex *)fxyze,isign,
                    ffc,ax,ay,az,affp,&we,nx,ny,nz,nxeh,nye,nze,nxh,nyh,
                    nzh);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tfield += time;

/* transform force to real space with standard procedure: updates fxyze */
      dtimer(&dtime,&itime,-1);
      isign = 1;
      if (kvec==1)
         cwfft3rv3((float complex *)fxyze,isign,mixup,sct,indx,indy,
                   indz,nxeh,nye,nze,nxhyz,nxyzh);
/* KNC function */
      else if (kvec==2)
         ckncwfft3rv3((float complex *)fxyze,isign,mixup,sct,indx,indy,
                      indz,nxeh,nye,nze,nxhyz,nxyzh);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tfft += time;

/* copy guard cells with standard procedure: updates fxyze */
      dtimer(&dtime,&itime,-1);
      if (kvec==1)
         ccguard3l(fxyze,nx,ny,nz,nxe,nye,nze);
/* KNC function */
      else if (kvec==2)
         cknccguard3l(fxyze,nx,ny,nz,nxe,nye,nze);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tguard += time;

/* push particles with standard procedure: updates part, wke */
      wke = 0.0;
      dtimer(&dtime,&itime,-1);
      if (kvec==1)
/*       cgpush3lt(partt,fxyze,qbme,dt,&wke,idimp,np,npe,nx,ny,nz,nxe, */
/*                 nye,nze,ipbc);                                      */
         cvgpush3lt(partt,fxyze,qbme,dt,&wke,idimp,np,npe,nx,ny,nz,nxe,
                    nye,nze,ipbc);
/* KNC function */
      else if (kvec==2)
         ckncgpush3lt(partt,fxyze,qbme,dt,&wke,idimp,np,npe,nx,ny,nz,
                      nxe,nye,nze,ipbc);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tpush += time;

/* sort particles by cell for standard procedure */
      if (sortime > 0) {
         if (ntime%sortime==0) {
            dtimer(&dtime,&itime,-1);
            if (kvec==1)
               cdsortp3yzlt(partt,partt2,npic,idimp,np,npe,ny1,nyz1);
/* KNC function */
            else if (kvec==2)
               ckncdsortp3yzlt(partt,partt2,npic,idimp,np,npe,ny1,nyz1);
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
         printf("Initial Field, Kinetic and Total Energies:\n");
         printf("%e %e %e\n",we,wke,wke+we);
      }
      ntime += 1;
      goto L500;
L2000:

/* * * * end main iteration loop * * * */

   printf("ntime = %i, kvec = %i\n",ntime,kvec);
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

   avx512_deallocate(npic);
   avx512_deallocate(ffc);
   avx512_deallocate(fxyze);
   avx512_deallocate(qe);
   if (sortime > 0)
      avx512_deallocate(partt2);
   avx512_deallocate(partt);

   return 0;
}
