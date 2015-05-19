/*--------------------------------------------------------------------*/
/* Basic parallel PIC library for MPI communications
   pplib2.c contains basic communications procedures for 1d partitions
   cppinit2 initializes parallel processing for C, returns
            number of processors and processor id.
   cppexit terminates parallel processing.
   cppabort aborts parallel processing.
   cpwtimera performs parallel local wall clock timing.
   cppsum performs parallel sum of a real vector.
   cppdsum performs parallel sum of a double precision vector.
   cppimax performs parallel maximum of an integer vector.
   cppdmax performs parallel maximum of a double precision vector.
   cppncguard2l copies data to guard cells in y for scalar data, linear
                interpolation, and distributed data with non-uniform
                partition.
   cppnaguard2l adds guard cells in y for scalar array, linear
                interpolation, and distributed data with non-uniform
                partition.
   cppnacguard2lL adds guard cells in y for vector array, linear
                  interpolation, and distributed data with non-uniform
                  partition.
   cpptpose performs a transpose of a complex scalar array, distributed
            in y, to a complex scalar array, distributed in x.
   cppntpose performs a transpose of an n component complex vector array,
            distributed in y, to an n component complex vector array,
            distributed in x.
   cppmove2 moves particles into appropriate spatial regions with periodic
            boundary conditions.  Assumes ihole list has been found.
   written by viktor k. decyk, ucla
   copyright 1995, regents of the university of california
   update: april 21, 2013                                         */

#include <stdlib.h>
#include <stdio.h>
#include <complex.h>
#include "mpi.h"
#include "pplib2.h"
      
/* common block for parallel processing
   nproc = number of real or virtual processors obtained
   lgrp = current communicator
   mreal = default datatype for reals
   mint = default datatype for integers
   mcplx = default datatype for complex type
   mdouble = default double precision type
   lworld = MPI_COMM_WORLD communicator
   msum = MPI_SUM
   mmax = MPI_MAX */

static int nproc;
static MPI_Comm lgrp, lworld;
static MPI_Datatype mreal, mint, mcplx, mdouble;
static MPI_Op msum, mmax;

static FILE *unit2 = NULL;

float vresult(float prec) {
   float vresult;
   vresult = prec;
   return vresult;
}

int iresult(int iprec) {
   int iresult;
   iresult = iprec;
   return iresult;
}

/*--------------------------------------------------------------------*/
void cppinit2(int *idproc, int *nvp, int argc, char *argv[]) {
/* this subroutine initializes parallel processing
   lgrp communicator = MPI_COMM_WORLD
   output: idproc, nvp
   idproc = processor id in lgrp communicator
   nvp = number of real or virtual processors obtained
local data */
   static int ibig = 2147483647;
   static float small = 1.0e-12;
   int ierror, flag, ndprec, idprec, iprec;
   float prec;
   prec = 1.0 + small;
   iprec = ibig + 1;
/* ndprec = (0,1) = (no,yes) use (normal,autodouble) precision */
   if (vresult(prec) > 1.0)
      ndprec = 1;
   else
      ndprec = 0;
/* idprec = (0,1) = (no,yes) use (normal,autodouble) integer precision */
   if (iresult(iprec) > 0)
      idprec = 1;
   else
      idprec = 0;
/* Open error file */
   unit2 = fopen("C.2","w");
/* indicate whether MPI_INIT has been called */
   ierror = MPI_Initialized(&flag);
   if (!flag) {
/* initialize the MPI execution environment */
      ierror = MPI_Init(&argc,&argv);
      if (ierror) exit(1);
   }
   lworld = MPI_COMM_WORLD;
   lgrp = lworld;
/* determine the rank of the calling process in the communicator */
   ierror = MPI_Comm_rank(lgrp,idproc);
/* determine the size of the group associated with a communicator */
   ierror = MPI_Comm_size(lgrp,&nproc);
/* set default datatypes */
   mint = MPI_INT;
   mdouble = MPI_DOUBLE;
/* single precision real */
   if (ndprec==0) {
      mreal = MPI_FLOAT;
      mcplx = MPI_COMPLEX;
   }
/* double precision real */
   else {
      mreal = MPI_DOUBLE;
      mcplx = MPI_DOUBLE_COMPLEX;
   }
/* single precision integer */
/* if (idprec==0)           */
/*    mint = MPI_INT;       */
/* double precision integer */
/* else                     */
/*    mint = MPI_LONG;      */
/* operators */
   msum = MPI_SUM;
   mmax = MPI_MAX;
   *nvp = nproc;
   return;
}

/*--------------------------------------------------------------------*/
void cppexit() {
/* this subroutine terminates parallel processing
local data */
   int ierror, flag;
/* close error file */
   fclose(unit2);
/* indicate whether MPI_INIT has been called */
   ierror = MPI_Initialized(&flag);
   if (flag) {
/* synchronize processes */
      ierror = MPI_Barrier(lworld);
/* terminate MPI execution environment */
      ierror = MPI_Finalize();
   }
   return;
}

/*--------------------------------------------------------------------*/
void cppabort() {
/* this subroutine aborts parallel processing
local data */
   int errorcode, ierror, flag;
/* close error file */
   fclose(unit2);
/* indicate whether MPI_INIT has been called */
   ierror = MPI_Initialized(&flag);
   if (flag) {
      errorcode = 1;
/* terminate MPI execution environment */
      ierror = MPI_Abort(lworld,errorcode);
   }
   return;
}

/*--------------------------------------------------------------------*/
void cpwtimera(int icntrl, float *time, double *dtime) {
/* this subroutine performs local wall clock timing
   input: icntrl, dtime
   icntrl = (-1,0,1) = (initialize,ignore,read) clock
   clock should be initialized before it is read!
   time = elapsed time in seconds
   dtime = current time
   written for mpi
local data */
   double jclock;
/* initialize clock */
   if (icntrl==(-1))
      *dtime = MPI_Wtime();
/* read clock and write time difference from last clock initialization */
   else if (icntrl==1) {
      jclock = *dtime;
      *dtime = MPI_Wtime();
      *time = (float ) (*dtime - jclock);
   }
   return;
}

/*--------------------------------------------------------------------*/
void cppsum(float f[], float g[], int nxp) {
/* this subroutine performs a parallel sum of a vector, that is:
   f[k][j] = sum over k of f[k][j]
   at the end, all processors contain the same summation.
   f = input and output real data
   g = scratch real array
   nxp = number of data values in vector
local data */
   int j, ierr;
/* perform sum */
   ierr = MPI_Allreduce(f,g,nxp,mreal,msum,lgrp);
/* copy output from scratch array */
   for (j = 0; j < nxp; j++) {
      f[j] = g[j];
   }
   return;
}

/*--------------------------------------------------------------------*/
void cppdsum(double f[], double g[], int nxp) {
/* this subroutine performs a parallel sum of a vector, that is:
   f[k][j] = sum over k of f[k][j]
   at the end, all processors contain the same summation.
   f = input and output double precision data
   g = scratch double precision array
   nxp = number of data values in vector
local data */
   int j, ierr;
/* perform sum */
   ierr = MPI_Allreduce(f,g,nxp,mdouble,msum,lgrp);
/* copy output from scratch array */
   for (j = 0; j < nxp; j++) {
      f[j] = g[j];
   }
   return;
}

/*--------------------------------------------------------------------*/
void cppimax(int f[], int g[], int nxp) {
/* this subroutine finds parallel maximum for each element of a vector
   that is, f[k][j] = maximum as a function of k of f[k][j]
   at the end, all processors contain the same maximum.
   f = input and output integer data
   g = scratch integer array
   nxp = number of data values in vector
local data */
   int j, ierr;
/* find maximum */
   ierr = MPI_Allreduce(f,g,nxp,mint,mmax,lgrp);
/* copy output from scratch array */
   for (j = 0; j < nxp; j++) {
      f[j] = g[j];
   }
   return;
}

/*--------------------------------------------------------------------*/
void cppdmax(double f[], double g[], int nxp) {
/* this subroutine finds parallel maximum for each element of a vector
   that is, f[k][j] = maximum as a function of k of f[k][j]
   at the end, all processors contain the same maximum.
   f = input and output double precision data
   g = scratch double precision array
   nxp = number of data values in vector
local data */
   int j, ierr;
/* find maximum */
   ierr = MPI_Allreduce(f,g,nxp,mdouble,mmax,lgrp);
/* copy output from scratch array */
   for (j = 0; j < nxp; j++) {
      f[j] = g[j];
   }
   return;
}

/*--------------------------------------------------------------------*/
void cppncguard2l(float f[], int nyp, int kstrt, int nvp, int nxv,
                  int nypmx) {
/* this subroutine copies data to guard cells in non-uniform partitions
   f[k][j] = real data for grid j,k in particle partition.
   the grid is non-uniform and includes one extra guard cell.
   output: f
   nyp = number of primary gridpoints in field partition
   it is assumed the nyp > 0.
   kstrt = starting data block number
   nvp = number of real or virtual processors
   nxv = first dimension of f, must be >= nx
   nypmx = maximum size of field partition, including guard cell.
   linear interpolation, for distributed data
local data */
   int j, ks, moff, kl, kr, ierr;
   MPI_Request msid;
   MPI_Status istatus;
/* special case for one processor */
   if (nvp==1) {
      for (j = 0; j < nxv; j++) {
        f[j+nxv*nyp] = f[j];
      }
      return;
   }
   ks = kstrt - 1;
   moff = nypmx*nvp + 2;
/* copy guard cells */
   kr = ks + 1;
   if (kr >= nvp)
      kr = kr - nvp;
   kl = ks - 1;
   if (kl < 0)
      kl = kl + nvp;
/* this segment is used for mpi computers */
   ierr = MPI_Irecv(&f[nxv*nyp],nxv,mreal,kr,moff,lgrp,&msid);
   ierr = MPI_Send(f,nxv,mreal,kl,moff,lgrp);
   ierr = MPI_Wait(&msid,&istatus);
   return;
}

/*--------------------------------------------------------------------*/
void cppnaguard2l(float f[], float scr[], int nyp, int nx, int kstrt,
                  int nvp, int nxv, int nypmx) {
/* this subroutine adds data from guard cells in non-uniform partitions
   f[k][j] = real data for grid j,k in particle partition.
   the grid is non-uniform and includes one extra guard cell.
   output: f, scr
   scr[j] = scratch array for particle partition
   nyp = number of primary gridpoints in particle partition
   it is assumed the nyp > 0.
   kstrt = starting data block number
   nvp = number of real or virtual processors
   nx = system length in x direction
   nxv = first dimension of f, must be >= nx
   nypmx = maximum size of field partition, including guard cells.
   linear interpolation, for distributed data
local data */
   int j, nx1, ks, moff, kl, kr, ierr;
   MPI_Request msid;
   MPI_Status istatus;
   nx1 = nx + 1;
/* special case for one processor */
   if (nvp==1) {
      for (j = 0; j < nx1; j++) {
         f[j] += f[j+nxv*nyp];
         f[j+nxv*nyp] = 0.0;
      }
      return;
   }
   ks = kstrt - 1;
   moff = nypmx*nvp + 1;
/* add guard cells */
   kr = ks + 1;
   if (kr >= nvp)
      kr = kr - nvp;
   kl = ks - 1;
   if (kl < 0)
      kl = kl + nvp;
/* this segment is used for mpi computers */
   ierr = MPI_Irecv(scr,nxv,mreal,kl,moff,lgrp,&msid);
   ierr = MPI_Send(&f[nxv*nyp],nxv,mreal,kr,moff,lgrp);
   ierr = MPI_Wait(&msid,&istatus);
/* add up the guard cells */
   for (j = 0; j < nx1; j++) {
      f[j] += scr[j];
      f[j+nxv*nyp] = 0.0;
   }
   return;
}

/*--------------------------------------------------------------------*/
void cppnacguard2l(float f[], float scr[], int nyp, int nx, int ndim,
                   int kstrt, int nvp, int nxv, int nypmx) {
/* this subroutine adds data from guard cells in non-uniform partitions
   f[k][j][ndim] = real data for grid j,k in particle partition.
   the grid is non-uniform and includes one extra guard cell.
   output: f, scr
   scr[j][ndim] = scratch array for particle partition
   nyp = number of primary gridpoints in particle partition
   it is assumed the nyp > 0.
   kstrt = starting data block number
   nvp = number of real or virtual processors
   nx = system length in x direction
   ndim = leading dimension of array f
   nxv = first dimension of f, must be >= nx
   nypmx = maximum size of field partition, including guard cells.
   linear interpolation, for distributed data
local data */
   int j, n, nx1, ks, moff, kl, kr, ierr;
   int nnxv;
   MPI_Request msid;
   MPI_Status istatus;
   nx1 = nx + 1;
/* special case for one processor */
   if (nvp==1) {
      for (j = 0; j < nx1; j++) {
         for (n = 0; n < ndim; n++) {
            f[n+ndim*j] += f[n+ndim*(j+nxv*nyp)];
            f[n+ndim*(j+nxv*nyp)] = 0.0;
         }
      }
      return;
   }
   ks = kstrt - 1;
   moff = nypmx*nvp + 1;
   nnxv = ndim*nxv;
/* add guard cells */
   kr = ks + 1;
   if (kr >= nvp)
      kr = kr - nvp;
   kl = ks - 1;
   if (kl < 0)
      kl = kl + nvp;
/* this segment is used for mpi computers */
   ierr = MPI_Irecv(scr,nnxv,mreal,kl,moff,lgrp,&msid);
   ierr = MPI_Send(&f[nnxv*nyp],nnxv,mreal,kr,moff,lgrp);
   ierr = MPI_Wait(&msid,&istatus);
/* add up the guard cells */
   for (j = 0; j < nx1; j++) {
      for (n = 0; n < ndim; n++) {
         f[n+ndim*j] += scr[n+ndim*j];
         f[n+ndim*(j+nxv*nyp)] = 0.0;
      }
   }
   return;
}

/*--------------------------------------------------------------------*/
void cpptpose(float complex f[], float complex g[], float complex s[],
              float complex t[], int nx, int ny, int kxp, int kyp,
              int kstrt, int nvp, int nxv, int nyv, int kxpd, int kypd) {
/* this subroutine performs a transpose of a matrix f, distributed in y,
   to a matrix g, distributed in x, that is,
   g[l][j][k+kyp*m] = f[m][k][j+kxp*l], where
   0 <= j < kxp, 0 <= k < kyp, 0 <= l < nx/kxp, 0 <= m < ny/kyp
   and where indices l and m can be distributed across processors.
   this subroutine sends and receives one message at a time, either
   synchronously or asynchronously. it uses a minimum of system resources
   f = complex input array
   g = complex output array
   s, t = complex scratch arrays
   nx/ny = number of points in x/y
   kxp/kyp = number of data values per block in x/y
   kstrt = starting data block number
   nvp = number of real or virtual processors
   nxv/nyv = first dimension of f/g
   kypd/kxpd = second dimension of f/g
local data */
   int n, j, k, ks, kxps, kyps, kxyp, id, joff, koff, ld, ierr;
   MPI_Request msid;
   MPI_Status istatus;
   ks = kstrt - 1;
   kxps = nx - kxp*ks;
   kxps = 0 > kxps ? 0 : kxps;
   kxps = kxp < kxps ? kxp : kxps;
   kyps = ny - kyp*ks;
   kyps = 0 > kyps ? 0 : kyps;
   kyps = kyp < kyps ? kyp : kyps;
   kxyp = kxp*kyp;
/* special case for one processor */
   if (nvp==1) {
      for (k = 0; k < kyp; k++) {
         for (j = 0; j < kxp; j++) {
            g[k+nyv*j] = f[j+nxv*k];
         }
      }
      return;
   }
/* this segment is used for shared memory computers */
/* for (m = 0; m < min(ny,nvp); m++) {                          */
/*    koff = kyp*m;                                             */
/*    for (k = 0; k < min(kyp,max(0,ny-koff)); k++) {           */
/*       for (l = 0; l < min(nx,nvp); l++) {                    */
/*          joff = kxp*l;                                       */
/*          for (j = 0; j < min(kxp,max(0,nx-joff)); j++) {     */
/*             g[k+koff+nyv*(j+joff)] = f[j+joff+nxv*(k+koff)]; */
/*          }                                                   */
/*       }                                                      */
/*    }                                                         */
/* }                                                            */
/* this segment is used for mpi computers */
   for (n = 0; n < nvp; n++) {
      id = n - ks;
      if (id < 0)
         id += nvp;
/* extract data to send */
      joff = kxp*id;
      ld = nx - joff;
      ld = 0 > ld ? 0 : ld;
      ld = kxp < ld ? kxp : ld;
      for (k = 0; k < kyps; k++) {
         for (j = 0; j < ld; j++) {
            s[j+ld*k] = f[j+joff+nxv*k];
         }
      }
      ld *= kyps;
/* post receive */
      ierr = MPI_Irecv(t,kxyp,mcplx,id,n,lgrp,&msid);
/* send data */
      ierr = MPI_Send(s,ld,mcplx,id,n,lgrp);
/* receive data */
      ierr = MPI_Wait(&msid,&istatus);
/* insert data received */
      koff = kyp*id;
      ld = ny - koff;
      ld = 0 > ld ? 0 : ld;
      ld = kyp < ld ? kyp : ld;
      for (k = 0; k < ld; k++) {
         for (j = 0; j < kxps; j++) {
            g[k+koff+nyv*j] = t[j+kxps*k];
         }
      }
   }
   return;
}

/*--------------------------------------------------------------------*/
void cppntpose(float complex f[], float complex g[], float complex s[],
               float complex t[], int nx, int ny, int kxp, int kyp,
               int kstrt, int nvp, int ndim, int nxv, int nyv, int kxpd,
               int kypd) {
/* this subroutine performs a transpose of a matrix f, distributed in y,
   to a matrix g, distributed in x, that is,
   g[l][j][k+kyp*m][1:ndim] = f[m][k][j+kxp*l][1:ndim], where
   0 <= j < kxp, 0 <= k < kyp, 0 <= l < nx/kxp, 0 <= m < ny/kyp
   and where indices l and m can be distributed across processors.
   this subroutine sends and receives one message at a time, either
   synchronously or asynchronously. it uses a minimum of system resources
   f = complex input array
   g = complex output array
   s, t = complex scratch arrays
   nx/ny = number of points in x/y
   kxp/kyp = number of data values per block in x/y
   kstrt = starting data block number
   nvp = number of real or virtual processors
   ndim = leading dimension of arrays f and g
   nxv/nyv = first dimension of f/g
   kypd/kxpd = second dimension of f/g
local data */
   int i, n, j, k, ks, kxps, kyps, kxyp, id, joff, koff, ld, ierr;
   int nnxv, nnyv;
   MPI_Request msid;
   MPI_Status istatus;
   ks = kstrt - 1;
   kxps = nx - kxp*ks;
   kxps = 0 > kxps ? 0 : kxps;
   kxps = kxp < kxps ? kxp : kxps;
   kyps = ny - kyp*ks;
   kyps = 0 > kyps ? 0 : kyps;
   kyps = kyp < kyps ? kyp : kyps;
   kxyp = ndim*kxp*kyp;
   nnxv = ndim*nxv;
   nnyv = ndim*nyv;
/* special case for one processor */
   if (nvp==1) {
      for (k = 0; k < kyp; k++) {
         for (j = 0; j < kxp; j++) {
            for (i = 0; i < ndim; i++) {
               g[i+ndim*k+nnyv*j] = f[i+ndim*j+nnxv*k];
            }
         }
      }
      return;
   }
/* this segment is used for shared memory computers */
/* for (m = 0; m < min(ny,nvp); m++) {                      */
/*    koff = kyp*m;                                         */
/*    for (k = 0; k < min(kyp,max(0,ny-koff)); k++) {       */
/*       for (l = 0; l < min(nx,nvp); l++) {                */
/*          joff = kxp*l;                                   */
/*          for (j = 0; j < min(kxp,max(0,nx-joff)); j++) { */
/*             for (i = 0; i < ndim; i++) {                 */
/*                g[i+ndim*(k+koff)+nnyv*(j+joff)] =        */
/*                f[i+ndim*(j+joff)+nnxv*(k+koff)];         */
/*             }                                            */
/*          }                                               */
/*       }                                                  */
/*    }                                                     */
/* }                                                        */
/* this segment is used for mpi computers */
   for (n = 0; n < nvp; n++) {
      id = n - ks;
      if (id < 0)
         id += nvp;
/* extract data to send */
      joff = kxp*id;
      ld = nx - joff;
      ld = 0 > ld ? 0 : ld;
      ld = kxp < ld ? kxp : ld;
      for (k = 0; k < kyps; k++) {
         for (j = 0; j < ld; j++) {
            for (i = 0; i < ndim; i++) {
               s[i+ndim*(j+ld*k)] = f[i+ndim*(j+joff)+nnxv*k];
            }
         }
      }
      ld *= ndim*kyps;
/* post receive */
      ierr = MPI_Irecv(t,kxyp,mcplx,id,n,lgrp,&msid);
/* send data */
      ierr = MPI_Send(s,ld,mcplx,id,n,lgrp);
/* receive data */
      ierr = MPI_Wait(&msid,&istatus);
/* insert data received */
      koff = kyp*id;
      ld = ny - koff;
      ld = 0 > ld ? 0 : ld;
      ld = kyp < ld ? kyp : ld;
      for (k = 0; k < ld; k++) {
         for (j = 0; j < kxps; j++) {
            for (i = 0; i < ndim; i++) {
               g[i+ndim*(k+koff)+nnyv*j] = t[i+ndim*(j+kxps*k)];
            }
         }
      }
   }
   return;
}

/*--------------------------------------------------------------------*/
void cppmove2(float part[], float edges[], int *npp, float sbufr[],
              float sbufl[], float rbufr[], float rbufl[], int ihole[],
              int ny, int kstrt, int nvp, int idimp, int npmax, int idps,
              int nbmax, int ntmax, int info[]) {
/* this subroutine moves particles into appropriate spatial regions
   periodic boundary conditions
   output: part, npp, sbufr, sbufl, rbufr, rbufl, info
   part[n][0] = position x of particle n in partition
   part[n][1] = position y of particle n in partition
   part[n][2] = velocity vx of particle n in partition
   part[n][3] = velocity vy of particle n in partition
   edges[0:1] = lower:lower boundary of particle partition
   npp = number of particles in partition
   sbufl = buffer for particles being sent to lower processor
   sbufr = buffer for particles being sent to upper processor
   rbufl = buffer for particles being received from lower processor
   rbufr = buffer for particles being received from upper processor
   ihole = location of holes left in particle arrays
   ny = system length in y direction
   kstrt = starting data block number
   nvp = number of real or virtual processors
   idimp = size of phase space = 4
   npmax = maximum number of particles in each partition.
   idps = number of partition boundaries
   nbmax =  size of buffers for passing particles between processors
   ntmax =  size of hole array for particles leaving processors
   info = status information
   info[0] = ierr = (0,N) = (no,yes) error condition exists
   info[1] = maximum number of particles per processor
   info[2] = minimum number of particles per processor
   info[3] = maximum number of buffer overflows
   info[4] = maximum number of particle passes required
local data */
/* iy = partitioned co-ordinate */
   int iy = 1;
   int ierr, ks, ih, iter, nps, itg, kl, kr, j, j1, j2, i;
   int joff, jin, nbsize, nter, mter, itermax, mpp;
   float any, yt;
   MPI_Request msid[4];
   MPI_Status istatus;
   int jsl[2], jsr[2], jss[2], ibflg[4], iwork[4];
   any = (float) ny;
   ks = kstrt - 1;
   nbsize = idimp*nbmax;
   iter = 2;
   nter = 0;
   info[0] = 0;
   info[4] = 0;
   ih = ihole[0];
   joff = 1;
   jin = 1;
   itermax = 2000;
   mpp = *npp;
/* ih = number of particles extracted from holes   */
/* joff = next hole location for extraction        */
/* jss[0] = number of holes available to be filled */
/* jin = next hole location to be filled           */
/* start loop */
L10: mter = 0;
   nps = 0;
/* buffer outgoing particles */
   jsl[0] = 0;
   jsr[0] = 0;
/* load particle buffers */
   for (j = 0; j < ih; j++) {
      j1 = ihole[j+joff] - 1;
      yt = part[iy+idimp*j1];
/* particles going down */
      if (yt < edges[0]) {
         if (ks==0)
            yt += any;
         if (jsl[0] < nbmax) {
            for (i = 0; i < idimp; i++) {
               sbufl[i+idimp*jsl[0]] = part[i+idimp*j1];
            }
            sbufl[iy+idimp*jsl[0]] = yt;
            jsl[0] += 1;
         }
         else {
            nps = 1;
            goto L50;
         }
      }
/* particles going up */
      else {
         if (ks==(nvp-1))
            yt -= any;
         if (jsr[0] < nbmax) {
            for (i = 0; i < idimp; i++) {
               sbufr[i+idimp*jsr[0]] = part[i+idimp*j1];
            }
            sbufr[iy+idimp*jsr[0]] = yt;
            jsr[0] += 1;
         }
         else {
            nps = 1;
            goto L50;
         }
      }
   }
L50: jss[0] = jsl[0] + jsr[0];
   joff += jss[0];
   ih -= jss[0];
/* check for full buffer condition */
   ibflg[2] = nps;
/* copy particle buffers */
L60: iter += 2;
   mter += 1;
/* special case for one processor */
   if (nvp==1) {
      jsl[1] = jsr[0];
      for (j = 0; j < jsl[1]; j++) {
         for (i = 0; i < idimp; i++) {
            rbufl[i+idimp*j] = sbufr[i+idimp*j];
         }
      }
      jsr[1] = jsl[0];
      for (j = 0; j < jsr[1]; j++) {
         for (i = 0; i < idimp; i++) {
            rbufr[i+idimp*j] = sbufl[i+idimp*j];
         }
      }
   }
/* this segment is used for mpi computers */
   else {
/* get particles from below and above */
      kr = ks + 1;
      if (kr >= nvp)
         kr -= nvp;
      kl = ks - 1;
      if (kl < 0)
         kl += nvp;
/* post receive */
      itg = iter - 1;
      ierr = MPI_Irecv(rbufl,nbsize,mreal,kl,itg,lgrp,&msid[0]);
      ierr = MPI_Irecv(rbufr,nbsize,mreal,kr,iter,lgrp,&msid[1]);
/* send particles */
      jsr[0] = idimp*jsr[0];
      ierr = MPI_Isend(sbufr,jsr[0],mreal,kr,itg,lgrp,&msid[2]);
      jsl[0] = idimp*jsl[0];
      ierr = MPI_Isend(sbufl,jsl[0],mreal,kl,iter,lgrp,&msid[3]);
/* wait for particles to arrive */
      ierr = MPI_Wait(&msid[0],&istatus);
      ierr = MPI_Get_count(&istatus,mreal,&nps);
      jsl[1] = nps/idimp;
      ierr = MPI_Wait(&msid[1],&istatus);
      ierr = MPI_Get_count(&istatus,mreal,&nps);
      jsr[1] = nps/idimp;
   }
/* check if particles must be passed further */
/* check if any particles coming from above belong here */
   jsl[0] = 0;
   jsr[0] = 0;
   jss[1] = 0;
   for (j = 0; j < jsr[1]; j++) {
      if (rbufr[iy+idimp*j] < edges[0])
         jsl[0] += 1;
      if (rbufr[iy+idimp*j] >= edges[1])
         jsr[0] += 1;
   }
   if (jsr[0] != 0)
      fprintf(unit2,"%d,Info: particles returning up\n",ks+1);
/* check if any particles coming from below belong here */
   for (j = 0; j < jsl[1]; j++) {
      if (rbufl[iy+idimp*j] >= edges[1])
         jsr[0] += 1;
      if (rbufl[iy+idimp*j] < edges[0])
         jss[1] += 1;
   }
   if (jss[1] != 0)
      fprintf(unit2,"%d,Info: particles returning down\n",ks+1);
   nps = jsl[0] + jsr[0] + jss[1];
   ibflg[1] = nps;
/* make sure sbufr and sbufl have been sent */
   if (nvp != 1) {
      ierr = MPI_Wait(&msid[2],&istatus);
      ierr = MPI_Wait(&msid[3],&istatus);
   }
   if (nps==0)
      goto L180;
/* remove particles which do not belong here */
/* first check particles coming from above */
   jsl[0] = 0;
   jsr[0] = 0;
   jss[1] = 0;
   for (j = 0; j < jsr[1]; j++) {
      yt = rbufr[iy+idimp*j];
/* particles going down */
      if (yt < edges[0]) {
         if (ks==0)
            yt += any;
         rbufr[iy+idimp*j] = yt;
         for (i = 0; i < idimp; i++) {
            sbufl[i+idimp*jsl[0]] = rbufr[i+idimp*j];
         }
         jsl[0] += 1;
      }
/* particles going up, should not happen */
      else if (yt >= edges[1]) {
         if (ks==(nvp-1))
            yt -= any;
         rbufr[iy+idimp*j] = yt;
         for (i = 0; i < idimp; i++) {
            sbufr[i+idimp*jsr[0]] = rbufr[i+idimp*j];
         }
         jsr[0] += 1;
      }
/* particles staying here */
      else {
         for (i = 0; i < idimp; i++) {
            rbufr[i+idimp*jss[1]] = rbufr[i+idimp*j];
         }
         jss[1] += 1;
      }
   }
   jsr[1] = jss[1];
/* next check particles coming from below */
   jss[1] = 0;
   for (j = 0; j < jsl[1]; j++) {
      yt = rbufl[iy+idimp*j];
/* particles going up */
      if (yt >= edges[1]) {
         if (jsr[0] < nbmax) {
            if (ks==(nvp-1))
               yt -= any;
            rbufl[iy+idimp*j] = yt;
            for (i = 0; i < idimp; i++) {
               sbufr[i+idimp*jsr[0]] = rbufl[i+idimp*j];
            }
            jsr[0] += 1;
         }
         else {
            jss[1] = 2*npmax;
            goto L170;
         }
      }
/* particles going down, should not happen */
      else if (yt < edges[0]) {
         if (jsl[0] < nbmax) {
            if (ks==0)
               yt += any;
            rbufl[iy+idimp*j] = yt;
            for (i = 0; i < idimp; i++) {
               sbufl[i+idimp*jsl[0]] = rbufl[i+idimp*j];
            }
            jsl[0] += 1;
         }
         else {
            jss[1] = 2*npmax;
            goto L170;
         }
      }
/* particles staying here */
      else {
         for (i = 0; i < idimp; i++) {
            rbufl[i+idimp*jss[1]] = rbufl[i+idimp*j];
         }
         jss[1] += 1;
      }
   }
L170: jsl[1] = jss[1];
/* check if move would overflow particle array */
L180: nps = mpp + jsl[1] + jsr[1] - jss[0];
   ibflg[0] = nps;
   nps = npmax < nps ? npmax : nps;
   ibflg[3] = -nps;
   cppimax(ibflg,iwork,4);
   info[1] = ibflg[0];
   info[2] = -ibflg[3];
   ierr = ibflg[0] - npmax;
   if (ierr > 0) {
      fprintf(unit2,"particle overflow error, ierr = %d\n",ierr);
      info[0] = ierr;
      *npp = mpp;
      return;
   }
/* distribute incoming particles from buffers */
/* distribute particles coming from below into holes */
   jss[1] = jss[0] < jsl[1] ? jss[0] : jsl[1];
   for (j = 0; j < jss[1]; j++) {
      j1 = ihole[j+jin] - 1;
      for (i = 0; i < idimp; i++) {
         part[i+idimp*j1] = rbufl[i+idimp*j];
      }
   }
   jin += jss[1];
   if (jss[0] > jsl[1]) {
      jss[1] = jss[0] - jsl[1];
      jss[1] = jss[1] < jsr[1] ? jss[1] : jsr[1];
   }
   else
      jss[1] = jsl[1] - jss[0];
   for (j = 0; j < jss[1]; j++) {
/* no more particles coming from below */
/* distribute particles coming from above into holes */
      if (jss[0] > jsl[1]) {
         j1 = ihole[j+jin] - 1;
         for (i = 0; i < idimp; i++) {
            part[i+idimp*j1] = rbufr[i+idimp*j];
         }
      }
/* no more holes */
/* distribute remaining particles from below into bottom */
      else {
         for (i = 0; i < idimp; i++) {
            part[i+idimp*(j+mpp)] = rbufl[i+idimp*(j+jss[0])];
         }
      }
   }
   if (jss[0] > jsl[1])
      jin += jss[1];
   nps = jsl[1] + jsr[1];
   if (jss[0] <= jsl[1]) {
      mpp += jsl[1] - jss[0];
      jss[0] = jsl[1];
   }
/* no more holes */
/* distribute remaining particles from above into bottom */
   jsr[1] = nps - jss[0];
   jsr[1] = 0 > jsr[1] ? 0 : jsr[1];
   jss[0] -= jsl[1];
   for (j = 0; j < jsr[1]; j++) {
      for (i = 0; i < idimp; i++) {
         part[i+idimp*(j+mpp)] = rbufr[i+idimp*(j+jss[0])];
     }
   }
   mpp += jsr[1];
/* holes left over */
/* fill up remaining holes in particle array with particles from bottom */
   if (ih==0) {
      jsr[1] = ihole[0] - jin + 1;
      jsr[1] = 0 > jsr[1] ? 0 : jsr[1];
      for (j = 0; j < jsr[1]; j++) {
         j1 = mpp - j - 1;
         j2 = ihole[jsr[1]-j+jin-1] - 1;
         if (j1 > j2) {
/* move particle only if it is below current hole */
            for (i = 0; i < idimp; i++) {
               part[i+idimp*j2] = part[i+idimp*j1];
            }
         }
      }
      jin += jsr[1];
      mpp -= jsr[1];
   }
   jss[0] = 0;
/* check if any particles have to be passed further */
   if (ibflg[2] > 0)
      ibflg[2] = 1;
   info[4] = info[4] > mter ? info[4] : mter;
   if (ibflg[1] > 0) {
      fprintf(unit2,"Info: particles being passed further = %d\n",
             ibflg[1]);
      if (iter < itermax)
         goto L60;
      ierr = -((iter-2)/2);
      fprintf(unit2,"Iteration overflow, iter = %d\n",ierr);
      info[0] = ierr;
      *npp = mpp;
      return;
   }
/* check if buffer overflowed and more particles remain to be checked */
   if (ibflg[2] > 0) {
      nter += 1;
      info[3] = nter;
      goto L10;
   }
   if (nter > 0) {
      fprintf(unit2,"Info: %d buffer overflows, nbmax=%d\n",nter,nbmax);
   }
   *npp = mpp;
   return;
}

/* Interfaces to Fortran */

/*--------------------------------------------------------------------*/
void cppinit2_(int *idproc, int *nvp, int *argc, char *argv[]) {
   cppinit2(idproc,nvp,*argc,argv);
   return;
}

void cppexit_() {
   cppexit();
   return;
}

void cppabort_() {
   cppabort();
   return;
}

/*--------------------------------------------------------------------*/
void cpwtimera_(int *icntrl, float *time, double *dtime) {
   cpwtimera(*icntrl,time,dtime);
   return;
}

/*--------------------------------------------------------------------*/
void cppsum_(float *f, float *g, int *nxp) {
   cppsum(f,g,*nxp);
   return;
}

/*--------------------------------------------------------------------*/
void cppdsum_(double *f, double *g, int *nxp) {
   cppdsum(f,g,*nxp);
   return;
}

/*--------------------------------------------------------------------*/
void cppimax_(int *f, int *g, int *nxp) {
   cppimax(f,g,*nxp);
   return;
}

/*--------------------------------------------------------------------*/
void cppdmax_(double *f, double *g, int *nxp) {
   cppdmax(f,g,*nxp);
   return;
}

/*--------------------------------------------------------------------*/
void cppncguard2l_(float *f, int *nyp, int *kstrt, int *nvp, int *nxv,
                   int *nypmx) {
   cppncguard2l(f,*nyp,*kstrt,*nvp,*nxv,*nypmx);
   return;
}
 
/*--------------------------------------------------------------------*/
void cppnaguard2l_(float *f, float *scr, int *nyp, int *nx, int *kstrt,
                   int *nvp, int *nxv, int *nypmx) {
   cppnaguard2l(f,scr,*nyp,*nx,*kstrt,*nvp,*nxv,*nypmx);
   return;
}

/*--------------------------------------------------------------------*/
void cppnacguard2l_(float *f, float *scr, int *nyp, int *nx, int *ndim,
                    int *kstrt, int *nvp, int *nxv, int *nypmx) {
   cppnacguard2l(f,scr,*nyp,*nx,*ndim,*kstrt,*nvp,*nxv,*nypmx);
   return;
}

/*--------------------------------------------------------------------*/
void cpptpose_(float complex *f, float complex *g, float complex *s,
               float complex *t, int *nx, int *ny, int *kxp, int *kyp,
               int *kstrt, int *nvp, int *nxv, int *nyv, int *kxpd,
               int *kypd) {
   cpptpose(f,g,s,t,*nx,*ny,*kxp,*kyp,*kstrt,*nvp,*nxv,*nyv,*kxpd,*kypd);
   return;
}

/*--------------------------------------------------------------------*/
void cppntpose_(float complex *f, float complex *g, float complex *s,
                float complex *t, int *nx, int *ny, int *kxp, int *kyp,
                int *kstrt, int *nvp, int *ndim, int *nxv, int *nyv,
                int *kxpd, int *kypd) {
   cppntpose(f,g,s,t,*nx,*ny,*kxp,*kyp,*kstrt,*nvp,*ndim,*nxv,*nyv,*kxpd,
             *kypd);
   return;
}

/*--------------------------------------------------------------------*/
void cppmove2_(float *part, float *edges, int *npp, float *sbufr,
               float *sbufl, float *rbufr, float *rbufl, int *ihole,
               int *ny, int *kstrt, int *nvp, int *idimp, int *npmax,
               int *idps, int *nbmax, int *ntmax, int *info) {
   cppmove2(part,edges,npp,sbufr,sbufl,rbufr,rbufl,ihole,*ny,*kstrt,*nvp,
            *idimp,*npmax,*idps,*nbmax,*ntmax,info);
   return;
}
