/*---------------------------------------------------------------------*/
/* Skeleton 3D Electromagnetic Vector PIC code */
/* written by Viktor K. Decyk, UCLA and Ricardo Fonseca, ISCTE */
#include <stdlib.h>
#include <stdio.h>
#include <complex.h>
#include <sys/time.h>
#include "vbpush3.h"
#include "avx512lib3.h"
#include "kncbpush3.h"

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
/* sortime = number of time steps between standard electron sorting */
/* relativity = (no,yes) = (0,1) = relativity is used */
   int idimp = 6, ipbc = 1, sortime = 20, relativity = 1;
/* wke/we = particle kinetic/electrostatic field energy             */
/* wf/wm/wt = magnetic field/transverse electric field/total energy */
   float wke = 0.0, we = 0.0, wf = 0.0, wm = 0.0, wt = 0.0;
/* kvec = (1,2) = run (autovector,KNC) version */
   int kvec = 1;

/* declare scalars for standard code */
   int j;
   int np, nx, ny, nz, nxh, nyh, nzh, nxe, nye, nze, nxeh;
   int nxyzh, nxhyz, npe, ny1, nyz1, ntime, nloop, isign;
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
/* npic = scratch array for reordering particles */
   int *mixup = NULL, *npic = NULL;

/* declare and initialize timing data */
   float time;
   struct timeval itime;
   float tdpost = 0.0, tguard = 0.0, tfft = 0.0, tfield = 0.0;
   float tdjpost = 0.0, tpush = 0.0, tsort = 0.0;
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
   dth = 0.0;

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
   avx512_fallocate(&cue,ndim*nxe*nye*nze,&irc);
   avx512_fallocate(&fxyze,ndim*nxe*nye*nze,&irc);
   avx512_fallocate(&bxyze,ndim*nxe*nye*nze,&irc);
   avx512_callocate(&exyz,ndim*nxeh*nye*nze,&irc);
   avx512_callocate(&bxyz,ndim*nxeh*nye*nze,&irc);
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
   cdistr3t(partt,vtx,vty,vtz,vx0,vy0,vz0,npx,npy,npz,idimp,npe,nx,ny,
            nz,ipbc);

/* initialize transverse electromagnetic fields */
   for (j = 0; j < ndim*nxeh*nye*nze; j++) {
      exyz[j] = 0.0 + 0.0*_Complex_I;
      bxyz[j] = 0.0 + 0.0*_Complex_I;
   }

   if (dt > 0.37*ci) {
      printf("Warning: Courant condition may be exceeded!\n");
   }

/* * * * start main iteration loop * * * */
 
L500: if (nloop <= ntime)
         goto L2000;
/*    printf("ntime = %i\n",ntime); */
 
/* deposit current with standard procedure: updates part, cue */
      dtimer(&dtime,&itime,-1);
      for (j = 0; j < ndim*nxe*nye*nze; j++) {
         cue[j] = 0.0;
      }
      if (relativity==1) {
         if (kvec==1)
/*          cgrjpost3lt(partt,cue,qme,dth,ci,np,npe,idimp,nx,ny,nz, */
/*                      nxe,nye,nze,ipbc);                          */
            cvgrjpost3lt(partt,cue,qme,dth,ci,np,npe,idimp,nx,ny,nz,
                         nxe,nye,nze,ipbc);
/* KNC function */
         else if (kvec==2)
             ckncgrjpost3lt(partt,cue,qme,dth,ci,np,npe,idimp,nx,ny,nz,
                            nxe,nye,nze,ipbc);
      }
      else {
         if (kvec==1)
/*          cgjpost3lt(partt,cue,qme,dth,np,npe,idimp,nx,ny,nz,nxe, */
/*                     nye,nze,ipbc);                               */
            cvgjpost3lt(partt,cue,qme,dth,np,npe,idimp,nx,ny,nz,nxe,
                        nye,nze,ipbc);
/* KNC function */
         else if (kvec==2)
            ckncgjpost3lt(partt,cue,qme,dth,np,npe,idimp,nx,ny,nz,nxe,
                          nye,nze,ipbc);
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
         cgpost3lt(partt,qe,qme,np,npe,idimp,nxe,nye,nze);
/*       cvgpost3lt(partt,qe,qme,np,npe,idimp,nxe,nye,nze); */
/* KNC function */
      else if (kvec==2)
         cknc2gpost3lt(partt,qe,qme,np,npe,idimp,nxe,nye,nze);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tdpost += time;

/* add guard cells with standard procedure: updates cue, qe */
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

/* transform current to fourier space with standard procedure: update cue */
      dtimer(&dtime,&itime,-1);
      isign = -1;
      if (kvec==1)
         cwfft3rv3((float complex *)cue,isign,mixup,sct,indx,indy,indz,
                nxeh,nye,nze,nxhyz,nxyzh);
/* KNC function */
      else if (kvec==2)
         ckncwfft3rv3((float complex *)cue,isign,mixup,sct,indx,indy,
                      indz,nxeh,nye,nze,nxhyz,nxyzh);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tfft += time;

/* take transverse part of current with standard procedure: updates cue */
      dtimer(&dtime,&itime,-1);
      if (kvec==1)
         ccuperp3((float complex *)cue,nx,ny,nz,nxeh,nye,nze);
/* KNC function */
      else if (kvec==2)
         cknccuperp3((float complex *)cue,nx,ny,nz,nxeh,nye,nze);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tfield += time;

/* calculate electromagnetic fields in fourier space with standard */
/* procedure: updates exyz, bxyz                                   */
      dtimer(&dtime,&itime,-1);
      if (ntime==0) {
         if (kvec==1)
            cvibpois33((float complex *)cue,bxyz,ffc,ci,&wm,nx,ny,nz,
                       nxeh,nye,nze,nxh,nyh,nzh);
/* KNC function */
         else if (kvec==2)
            ckncibpois33((float complex *)cue,bxyz,ffc,ci,&wm,nx,ny,nz,
                         nxeh,nye,nze,nxh,nyh,nzh);
         wf = 0.0;
         dth = 0.5*dt;
      }
      else {
         if (kvec==1)
            cvmaxwel3(exyz,bxyz,(float complex *)cue,ffc,ci,dt,&wf,&wm,
                      nx,ny,nz,nxeh,nye,nze,nxh,nyh,nzh);
/* KNC function */
         else if (kvec==2)
            ckncmaxwel3(exyz,bxyz,(float complex *)cue,ffc,ci,dt,&wf,&wm,
                        nx,ny,nz,nxeh,nye,nze,nxh,nyh,nzh);
      }
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tfield += time;

/* calculate force/charge in fourier space with standard procedure: */
/* updates fxyze                                                    */
      dtimer(&dtime,&itime,-1);
      isign = -1;
      if (kvec==1)
         cvpois33((float complex *)qe,(float complex *)fxyze,isign,ffc,
                  ax,ay,az,affp,&we,nx,ny,nz,nxeh,nye,nze,nxh,nyh,nzh);
/* KNC function */
      else if (kvec==2)
         ckncpois33((float complex *)qe,(float complex *)fxyze,isign,ffc,
                     ax,ay,az,affp,&we,nx,ny,nz,nxeh,nye,nze,nxh,nyh,
                     nzh);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tfield += time;

/* add longitudinal and transverse electric fields with standard */
/* procedure: updates fxyze */
      dtimer(&dtime,&itime,-1);
      isign = 1;
      if (kvec==1)
         cvemfield3((float complex *)fxyze,exyz,ffc,isign,nx,ny,nz,nxeh,
                    nye,nze,nxh,nyh,nzh);
/* KNC function */
      else if (kvec==2)
         ckncemfield3((float complex *)fxyze,exyz,ffc,isign,nx,ny,nz,
                      nxeh,nye,nze,nxh,nyh,nzh);
/* copy magnetic field with standard procedure: updates bxyze */
      isign = -1;
      if (kvec==1)
         cvemfield3((float complex *)bxyze,bxyz,ffc,isign,nx,ny,nz,nxeh,
                    nye,nze,nxh,nyh,nzh);
/* KNC function */
      else if (kvec==2)
         ckncemfield3((float complex *)bxyze,bxyz,ffc,isign,nx,ny,nz,
                      nxeh,nye,nze,nxh,nyh,nzh);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tfield += time;

/* transform electric force to real space with standard procedure: */
/* updates fxyze                                                   */
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

/* transform magnetic force to real space with standard procedure: */
/* updates bxyze                                                   */
      dtimer(&dtime,&itime,-1);
      isign = 1;
      if (kvec==1)
         cwfft3rv3((float complex *)bxyze,isign,mixup,sct,indx,indy,
                   indz,nxeh,nye,nze,nxhyz,nxyzh);
/* KNC function */
      else if (kvec==2)
         ckncwfft3rv3((float complex *)bxyze,isign,mixup,sct,indx,indy,
                      indz,nxeh,nye,nze,nxhyz,nxyzh);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tfft += time;

/* copy guard cells with standard procedure: updates fxyze, bxyze */
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

/* push particles with standard procedure: updates part, wke */
      wke = 0.0;
      dtimer(&dtime,&itime,-1);
      if (relativity==1) {
         if (kvec==1)
/*          cgrbpush3lt(partt,fxyze,bxyze,qbme,dt,dth,ci,&wke,idimp, */
/*                      np,npe,nx,ny,nz,nxe,nye,nze,ipbc);           */
            cvgrbpush3lt(partt,fxyze,bxyze,qbme,dt,dth,ci,&wke,idimp,
                         np,npe,nx,ny,nz,nxe,nye,nze,ipbc);
/* KNC function */
         else if (kvec==2)
            ckncgrbpush3lt(partt,fxyze,bxyze,qbme,dt,dth,ci,&wke,idimp,
                           np,npe,nx,ny,nz,nxe,nye,nze,ipbc);
      }
      else {
         if (kvec==1)
/*          cgbpush3lt(partt,fxyze,bxyze,qbme,dt,dth,&wke,idimp,np, */
/*                     npe,nx,ny,nz,nxe,nye,nze,ipbc);              */
            cvgbpush3lt(partt,fxyze,bxyze,qbme,dt,dth,&wke,idimp,np,
                       npe,nx,ny,nz,nxe,nye,nze,ipbc);
/* KNC function */
         else if (kvec==2)
            ckncgbpush3lt(partt,fxyze,bxyze,qbme,dt,dth,&wke,idimp,np,
                          npe,nx,ny,nz,nxe,nye,nze,ipbc);
      }
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

   avx512_deallocate(npic);
   avx512_deallocate(ffc);
   avx512_deallocate(bxyz);
   avx512_deallocate(exyz);
   avx512_deallocate(bxyze);
   avx512_deallocate(fxyze);
   avx512_deallocate(cue);
   avx512_deallocate(qe);
   if (sortime > 0)
      avx512_deallocate(partt2);
   avx512_deallocate(partt);

   return 0;
}
