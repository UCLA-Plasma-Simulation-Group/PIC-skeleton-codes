/*---------------------------------------------------------------------*/
/* Skeleton 3D Electrostatic MPI/OpenMP PIC code */
/* written by Viktor K. Decyk, UCLA */
#include <stdlib.h>
#include <stdio.h>
#include <complex.h>
#include <sys/time.h>
#include "mppush3.h"
#include "mpplib3.h"
#include "omplib.h"

void dtimer(double *time, struct timeval *itime, int icntrl);

int main(int argc, char *argv[]) {
/* indx/indy/indz = exponent which determines grid points in x/y/z: */
/* direction: nx = 2**indx, ny = 2**indy, nz = 2**indz */
   int indx =   7, indy =   7, indz =   7;
/* npx/npy/npz = number of electrons distributed in x/y/z direction */
   int npx =  384, npy =   384, npz =   384;
/* ndim = number of velocity coordinates = 3 */
   int ndim = 3;
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
   int idimp = 6, ipbc = 1;
/* idps = number of partition boundaries = 4 */
/* idds = dimensionality of domain decomposition = 2 */
   int idps = 4, idds =    2;
/* wke/we/wt = particle kinetic/electric field/total energy */
   float wke = 0.0, we = 0.0, wt = 0.0;
/* mx/my/mz = number of grids in x/y/z in sorting tiles */
/* sorting tiles, should be less than or equal to 16    */
   int mx = 8, my = 8, mz = 8;
/* fraction of extra particles needed for particle management */
   float xtras = 0.2;
/* declare scalars for standard code */
   int j;
   int nx, ny, nz, nxh, nyh, nzh, nxe, nye, nze, nxeh, nnxe;
   int nxyzh, nxhyz, mx1, ntime, nloop, isign, ierr;
   float qbme, affp;
   double np;

/* declare scalars for MPI code */
   int ntpose = 1;
   int nvpy, nvpz, nvp, idproc, kstrt, npmax, kyp, kzp;
   int kxyp, kyzp, kzyp, nypmx, nzpmx, nypmn, nzpmn;
   int npp, nps, myp1, mzp1, mxyzp1, mxzyp1;

/* declare scalars for OpenMP code */
   int nppmx, nppmx0, nbmaxp, ntmaxp, npbmx, irc;
   int nvpp;

/* declare arrays for standard code: */
/* part = particle array */
   float *part = NULL;
/* qe = electron charge density with guard cells */
   float *qe = NULL;
/* fxyze = smoothed electric field with guard cells */
   float *fxyze = NULL;
/* qt, qs = scalar charge density field arrays in fourier space */
   float complex *qt = NULL, *qs = NULL;
/* fxyzt, fxyzs = vector electric field arrays in fourier space */
   float complex *fxyzt = NULL, *fxyzs = NULL;
/* ffc = form factor array for poisson solver */
   float complex *ffc = NULL;
/* mixup = bit reverse table for FFT */
   int *mixup = NULL;
/* sct = sine/cosine table for FFT */
   float complex *sct = NULL;
   double wtot[4], work[4];

/* declare arrays for MPI code: */
/* bs/br = complex send/receive buffers for data transpose */
   float complex *bs = NULL, *br = NULL;
/* sbufl/sbufr = particle buffers sent to nearby processors */
/* rbufl/rbufr = particle buffers received from nearby processors */
   float *sbufl = NULL, *sbufr = NULL, *rbufl = NULL, *rbufr = NULL;
/* edges[0:1] = lower:upper y boundaries of particle partition */
/* edges[2:3] = back:front z boundaries of particle partition */
   float *edges = NULL;
/* nyzp[0:1] = number of primary (complete) gridpoints in y/z */
/* noff[0:1] = lowermost global gridpoint in y/z */
   int *nyzp = NULL, *noff = NULL;
/* scr/scs = guard cell buffers received/sent from nearby processors */
   float *scr = NULL, *scs = NULL;

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
/* mcls = number offsets received from corner processors */
   int *mcls = NULL;

/* declare and initialize timing data */
   float time;
   struct timeval itime;
   float tdpost = 0.0, tguard = 0.0, ttp = 0.0, tfield = 0.0;
   float tpush = 0.0, tsort = 0.0, tmov = 0.0;
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
   np = ((double) npx)*((double) npy)*((double) npz);
/* nx/ny/nz = number of grid points in x/y direction */
   nx = 1L<<indx; ny = 1L<<indy; nz = 1L<<indz;
   nxh = nx/2; nyh = 1 > ny/2 ? 1 : ny/2; nzh = 1 > nz/2 ? 1 : nz/2;
   nxe = nx + 2; nye = ny + 2; nze = nz + 2;
   nxeh = nxe/2;  nnxe = ndim*nxe;
   nxyzh = (nx > ny ? nx : ny); nxyzh = (nxyzh > nz ? nxyzh : nz)/2;
   nxhyz = nxh > ny ? nxh : ny; nxhyz = nxhyz > nz ? nxhyz : nz;
/* mx1 = number of tiles in x direction */
   mx1 = (nx - 1)/mx + 1;
/* nloop = number of time steps in simulation */
/* ntime = current time step */
   nloop = tend/dt + .0001; ntime = 0;
   qbme = qme;
   affp = ((double) nx)*((double) ny)*((double) nz)/np;

/* nvp = number of MPI ranks */
/* initialize for distributed memory parallel processing */
   cppinit2(&idproc,&nvp,argc,argv);
   kstrt = idproc + 1;
/* obtain 2D partition (nvpy,nvpz) from nvp: */
/* nvpy/nvpz = number of processors in y/z */
   cfcomp32(nvp,nx,ny,nz,&nvpy,&nvpz,&ierr);
   if (ierr != 0) {
      if (kstrt==1) {
         printf("cfcomp32 error: nvp,nvpy,nvpz=%d,%d,%d\n",nvp,nvpy,nvpz);
      }
      goto L3000;
   }

/* initialize data for MPI code */
   edges = (float *) malloc(idps*sizeof(float));
   nyzp = (int *) malloc(idds*sizeof(int));
   noff = (int *) malloc(idds*sizeof(int));
/* calculate partition variables: */
/* edges, nyzp, noff, nypmx, nzpmx, nypmn, nzpmn                          */
/* edges[0:1] = lower:upper boundary of particle partition in y           */
/* edges[2:3] = back:front boundary of particle partition in z            */
/* nyzp[0:1] = number of primary (complete) gridpoints in y/z             */
/* noff[0:1] = lowermost global gridpoint in y/z in particle partition    */
/* nypmx = maximum size of particle partition in y, including guard cells */
/* nzpmx = maximum size of particle partition in z, including guard cells */
/* nypmn = minimum value of nyzp[0]                                       */
/* nzpmn = minimum value of nyzp[1]                                       */
   cpdicomp32l(edges,nyzp,noff,&nypmx,&nzpmx,&nypmn,&nzpmn,ny,nz,kstrt,
               nvpy,nvpz,idps,idds);
   if (kstrt==1) {
      if (nypmn < 1) {
                     printf("combination not supported nvpy,ny= %d,%d\n",
                            nvpy,ny);
      }
      if (nzpmn < 1) {
                     printf("combination not supported nvpz,nz= %d,%d\n",
                            nvpz,nz);
      }
   }
   if ((nypmn < 1) || (nzpmn < 1)) goto L3000;
/* initialize additional scalars for MPI code */
/* kyp = number of complex grids in each field partition in y direction */
   kyp = (ny - 1)/nvpy + 1;
/* kzp = number of complex grids in each field partition in z direction */
   kzp = (nz - 1)/nvpz + 1;
/* kxyp = number of complex grids in each field partition in x direction */
/* in transposed data */
   kxyp = (nxh - 1)/nvpy + 1;
/* kyzp = number of complex grids in each field partition in y direction, */
/* in transposed data */
   kyzp = (ny - 1)/nvpz + 1; kzyp = kyzp > kyp ? kyzp : kyp;
/* npmax = maximum number of electrons in each partition */
   npmax = (np/nvp)*1.25;
/* myp1/mzp1 = number of tiles in y/z direction */
   myp1 = (nyzp[0] - 1)/my + 1; mzp1 = (nyzp[1] - 1)/mz + 1;
/* mxzyp1 = mx1*max(max(mzp1),max(myp1)) */
   mxzyp1 = (nzpmx-2)/mz + 1; mxyzp1 = (nypmx-2)/my + 1;
   mxzyp1 = mx1*(mxzyp1 > mxyzp1 ? mxzyp1 : mxyzp1);
   mxyzp1 = mx1*myp1*mzp1;

/* allocate data for standard code */
   part = (float *) malloc(idimp*npmax*sizeof(float));
   qe = (float *) malloc(nxe*nypmx*nzpmx*sizeof(float));
   fxyze = (float *) malloc(ndim*nxe*nypmx*nzpmx*sizeof(float));
   qt = (float complex *) malloc(nze*kxyp*kyzp*sizeof(float complex));
   qs = (float complex *) malloc(nye*kxyp*nzpmx*sizeof(float complex));
   fxyzt = (float complex *) malloc(ndim*nze*kxyp*kyzp
                                    *sizeof(float complex));
   fxyzs = (float complex *) malloc(ndim*nye*kxyp*nzpmx
                                    *sizeof(float complex));
   ffc = (float complex *) malloc(nzh*kxyp*kyzp*sizeof(float complex));
   mixup = (int *) malloc(nxhyz*sizeof(int));
   sct = (float complex *) malloc(nxyzh*sizeof(float complex));
   kpic = (int *) malloc(mxyzp1*sizeof(int));

/* allocate data for MPI code */
   bs = (float complex *) malloc(ndim*kxyp*kzyp*kzp*sizeof(float complex));
   br = (float complex *) malloc(ndim*kxyp*kzyp*kzp*sizeof(float complex));
   scr = (float *) malloc(nxe*nypmx*sizeof(float));
   scs = (float *) malloc(nnxe*2*nzpmx*sizeof(float));

/* prepare fft tables */
   cwpfft32rinit(mixup,sct,indx,indy,indz,nxhyz,nxyzh);
/* calculate form factors */
   isign = 0;
   cmppois332(qt,fxyzt,isign,ffc,ax,ay,az,affp,&we,nx,ny,nz,kstrt,nvpy,
              nvpz,nze,kxyp,kyzp,nzh);
/* initialize electrons */
   nps = 1;
   npp = 0;
   cpdistr32(part,edges,&npp,nps,vtx,vty,vtz,vx0,vy0,vz0,npx,npy,npz,nx,
             ny,nz,idimp,npmax,idps,ipbc,&ierr);
/* check for particle initialization error */
   if (ierr != 0) {
      if (kstrt==1) {
         printf("particle initialization error: ierr=%d\n",ierr);
      }
      goto L3000;
   }

/* find number of particles in each of mx, my, mz tiles: */
/* updates kpic, nppmx                                   */
   cppdblkp3l(part,kpic,npp,noff,&nppmx,idimp,npmax,mx,my,mz,mx1,myp1,
              mxyzp1,idds,&irc);
   if (irc != 0) {
      printf("%d,cppdblkp3l error, irc=%d\n",kstrt,irc);
      cppabort();
      exit(1);
   }
/* allocate vector particle data */
   nppmx0 = (1.0 + xtras)*nppmx;
   ntmaxp = xtras*nppmx;
   npbmx = xtras*nppmx;
   nbmaxp = 0.125*mxzyp1*npbmx;
   sbufl = (float *) malloc(idimp*nbmaxp*2*sizeof(float));
   sbufr = (float *) malloc(idimp*nbmaxp*2*sizeof(float));
   rbufl = (float *) malloc(idimp*nbmaxp*2*sizeof(float));
   rbufr = (float *) malloc(idimp*nbmaxp*2*sizeof(float));
   ppart = (float *) malloc(idimp*nppmx0*mxyzp1*sizeof(float));
   ppbuff = (float *) malloc(idimp*npbmx*mxyzp1*sizeof(float));
   ncl = (int *) malloc(26*mxyzp1*sizeof(int));
   iholep = (int *) malloc(2*(ntmaxp+1)*mxyzp1*sizeof(int));
   ncll = (int *) malloc(3*mxzyp1*6*sizeof(int));
   nclr = (int *) malloc(3*mxzyp1*6*sizeof(int));
   mcll = (int *) malloc(3*mxzyp1*6*sizeof(int));
   mclr = (int *) malloc(3*mxzyp1*6*sizeof(int));
   mcls = (int *) malloc(3*(mx1+1)*4*sizeof(int));

/* copy ordered particle data for OpenMP */
   cpppmovin3l(part,ppart,kpic,npp,noff,nppmx0,idimp,npmax,mx,my,mz,mx1,
               myp1,mxyzp1,idds,&irc);
   if (irc != 0) { 
      printf("%d,cpppmovin3l overflow error, irc=%d\n",kstrt,irc);
      cppabort();
      exit(1);
   }
/* sanity check */
   cpppcheck3l(ppart,kpic,noff,nyzp,idimp,nppmx0,nx,mx,my,mz,mx1,myp1,
               mzp1,idds,&irc);
   if (irc != 0) { 
      printf("%d,cpppcheck3l error: irc=%d\n",kstrt,irc);
      cppabort();
      exit(1);
   }

/* * * * start main iteration loop * * * */
 
L500: if (nloop <= ntime)
         goto L2000;
/*    if (kstrt==1) printf("ntime = %i\n",ntime); */
 
/* deposit charge with OpenMP: updates qe */
      dtimer(&dtime,&itime,-1);
      for (j = 0; j < nxe*nypmx*nzpmx; j++) {
         qe[j] = 0.0;
      }
      cppgppost32l(ppart,qe,kpic,noff,qme,nppmx0,idimp,mx,my,mz,nxe,
                   nypmx,nzpmx,mx1,myp1,mxyzp1,idds);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tdpost += time;

/* add guard cells with OpenMP: updates qe */
      dtimer(&dtime,&itime,-1);
      cppaguard32xl(qe,nyzp,nx,nxe,nypmx,nzpmx,idds);
      cppnaguard32l(qe,scs,scr,nyzp,kstrt,nvpy,nvpz,nx,nxe,nypmx,nzpmx,
                    idds);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tguard += time;

/* transform charge to fourier space with OpenMP: updates qt, modifies qe */
      dtimer(&dtime,&itime,-1);
      isign = -1;
      cwppfft32rm((float complex *)qe,qs,qt,bs,br,isign,ntpose,mixup,
                  sct,&ttp,indx,indy,indz,kstrt,nvpy,nvpz,nxeh,nye,nze,
                  kxyp,kyp,kyzp,kzp,kxyp,nypmx,kyzp,nzpmx,kzyp,nxhyz,
                  nxyzh);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tfft[0] += time;
      tfft[1] += ttp;

/* calculate force/charge in fourier space with OpenMP: */
/* updates fxyzt, we                                    */
      dtimer(&dtime,&itime,-1);
      isign = -1;
      cmppois332(qt,fxyzt,isign,ffc,ax,ay,az,affp,&we,nx,ny,nz,kstrt,
                 nvpy,nvpz,nze,kxyp,kyzp,nzh);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tfield += time;

/* transform force to real space with OpenMP: updates fxyze */
/* modifies fxyzt */
      dtimer(&dtime,&itime,-1);
      isign = 1;
      cwppfft32rm3((float complex *)fxyze,fxyzs,fxyzt,bs,br,isign,
                   ntpose,mixup,sct,&ttp,indx,indy,indz,kstrt,nvpy,nvpz,
                   nxeh,nye,nze,kxyp,kyp,kyzp,kzp,kxyp,nypmx,kyzp,nzpmx,
                   kzyp,nxhyz,nxyzh);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tfft[0] += time;
      tfft[1] += ttp;

/* copy guard cells with OpenMP: updates fxyze */
      dtimer(&dtime,&itime,-1);
      cppncguard32l(fxyze,scs,nyzp,kstrt,nvpy,nvpz,nnxe,nypmx,nzpmx,
                    idds);
      cppcguard32xl(fxyze,nyzp,nx,ndim,nxe,nypmx,nzpmx,idds);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tguard += time;

/* push particles with OpenMP: updates part, wke */
      wke = 0.0;
      dtimer(&dtime,&itime,-1);
/* updates part, wke */
/*    cppgppush32l(ppart,fxyze,kpic,noff,nyzp,qbme,dt,&wke,idimp,nppmx0, */
/*                 nx,ny,nz,mx,my,mz,nxe,nypmx,nzpmx,mx1,myp1,mxyzp1,    */
/*                 idds,ipbc);                                           */
/* updates ppart, wke, ncl, iholep, irc */
      cppgppushf32l(ppart,fxyze,kpic,ncl,iholep,noff,nyzp,qbme,dt,&wke,
                    idimp,nppmx0,nx,ny,nz,mx,my,mz,nxe,nypmx,nzpmx,mx1,
                    myp1,mxyzp1,ntmaxp,idds,&irc);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tpush += time;
      if (irc != 0) { 
         printf("%d,cppgppushf32l error: irc=%d\n",kstrt,irc);
         cppabort();
         exit(1);
      }

/* reorder particles by tile with OpenMP */
/* first part of particle reorder on x, y and z cell */
/* with mx, my, mz tiles: */
      dtimer(&dtime,&itime,-1);
/* updates ppart, ppbuff, sbufl, sbufr, ncl, iholep, ncll, nclr, irc */
/*    cppporder32la(ppart,ppbuff,sbufl,sbufr,kpic,ncl,iholep,ncll,     */
/*                  nclr,noff,nyzp,idimp,nppmx0,nx,ny,nz,mx,my,mz,mx1, */
/*                  myp1,mzp1,mxzyp1,npbmx,ntmaxp,nbmaxp,idds,&irc);   */
/* updates: ppart, ppbuff, sbufl, sbufr, ncl, ncll, nclr, irc */
      cppporderf32la(ppart,ppbuff,sbufl,sbufr,ncl,iholep,ncll,nclr,
                     idimp,nppmx0,mx1,myp1,mzp1,mxzyp1,npbmx,ntmaxp,
                     nbmaxp,&irc);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tsort += time;
      if (irc != 0) {
         printf("%d,cppporderf32la error:ntmaxp,irc=%d,%d\n",kstrt,
                ntmaxp,irc);
         cppabort();
         exit(1);
      }
/* move particles into appropriate spatial regions: */
/* updates rbufr, rbufl, mcll, mclr, mcls */
      dtimer(&dtime,&itime,-1);
      cpppmove32(sbufr,sbufl,rbufr,rbufl,ncll,nclr,mcll,mclr,mcls,kstrt,
                 nvpy,nvpz,idimp,nbmaxp,mx1,myp1,mzp1,mxzyp1,&irc);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tmov += time;
      if (irc != 0) {
         if (kstrt==1) {
            printf("cpppmove32 error: nbmaxp, irc=%d,%d\n",nbmaxp,irc);
            goto L3000;
         }
      }
/* second part of particle reorder on x and y cell */
/* with mx, my, mz tiles: updates ppart, kpic */
      dtimer(&dtime,&itime,-1);
      cppporder32lb(ppart,ppbuff,rbufl,rbufr,kpic,ncl,iholep,mcll,mclr,
                    mcls,idimp,nppmx0,mx1,myp1,mzp1,mxzyp1,npbmx,ntmaxp,
                    nbmaxp,&irc);
      dtimer(&dtime,&itime,1);
      time = (float) dtime;
      tsort += time;
      if (irc != 0) {
         printf("%d,cppporder32lb error:nppmx0,irc=%d,%d\n",kstrt,nppmx0,
                irc);
         cppabort();
         exit(1);
      }

/* energy diagnostic */
      wtot[0] = we;
      wtot[1] = wke;
      wtot[2] = 0.0;
      wtot[3] = we + wke;
      cppdsum(wtot,work,4);
      we = wtot[0];
      wke = wtot[1];
      if (ntime==0) {
         if (kstrt==1) {
            printf("Initial Field, Kinetic and Total Energies:\n");
            printf("%e %e %e\n",we,wke,wke+we);
         }
      }
      ntime += 1;
      goto L500;
L2000:

/* * * * end main iteration loop * * * */

   if (kstrt==1) {
      printf("ntime = %i\n",ntime);
      printf("MPI nodes nvpy, nvpz = %i,%i\n",nvpy,nvpz);
      printf("Final Field, Kinetic and Total Energies:\n");
      printf("%e %e %e\n",we,wke,wke+we);

      printf("\n");
      printf("deposit time = %f\n",tdpost);
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
