/*--------------------------------------------------------------------*/
/* Basic parallel PIC library for MPI communications with OpenMP
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
   cpppmove2 moves particles into appropriate spatial regions for tiled
             distributed data.
   written by viktor k. decyk, ucla
   copyright 1995, regents of the university of california
   update: may 10, 2015                                         */

#include <stdlib.h>
#include <stdio.h>
#include <complex.h>
#include "mpi.h"
#include "mpplib2.h"
      
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
#pragma omp parallel for private(j,k)
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
#pragma omp parallel for private(j,k)
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
#pragma omp parallel for private(j,k)
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
#pragma omp parallel for private(i,j,k)
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
#pragma omp parallel for private(i,j,k)
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
#pragma omp parallel for private(i,j,k)
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
void cpppmove2(float sbufr[], float sbufl[], float rbufr[], 
               float rbufl[], int ncll[], int nclr[], int mcll[],
               int mclr[], int kstrt, int nvp, int idimp, int nbmax,
               int mx1) {
/* this subroutine moves particles into appropriate spatial regions
   for distributed data, with 1d domain decomposition in y.
   tiles are assumed to be arranged in 2D linear memory
   output: rbufr, rbufl, mcll, mclr
   sbufl = buffer for particles being sent to lower processor
   sbufr = buffer for particles being sent to upper processor
   rbufl = buffer for particles being received from lower processor
   rbufr = buffer for particles being received from upper processor
   ncll = particle number being sent to lower processor
   nclr = particle number being sent to upper processor
   mcll = particle number being received from lower processor
   mclr = particle number being received from upper processor
   kstrt = starting data block number
   nvp = number of real or virtual processors
   idimp = size of phase space = 4 or 5
   nbmax =  size of buffers for passing particles between processors
   mx1 = (system length in x direction - 1)/mx + 1
local data */
   int ierr, ks, kl, kr, i, j, jsl, jsr;
   int nbsize, ncsize;
   int itg[4] = {3,4,5,6};
   MPI_Request msid[8];
   MPI_Status istatus;
   ks = kstrt - 1;
   nbsize = idimp*nbmax;
   ncsize = 3*mx1;
/* copy particle buffers: update rbufl, rbufr, mcll, mclr */
/* special case for one processor */
   if (nvp==1) {
      for (j = 0; j < mx1; j++) {
         for (i = 0; i < 3; i++) {
            mcll[i+3*j] = nclr[i+3*j];
        }
      }
      for (j = 0; j < mx1; j++) {
         for (i = 0; i < 3; i++) {
            mclr[i+3*j] = ncll[i+3*j];
        }
      }
      for (j = 0; j < nclr[3*mx1-1]; j++) {
         for (i = 0; i < idimp; i++) {
            rbufl[i+idimp*j] = sbufr[i+idimp*j];
         }
      }
      for (j = 0; j < ncll[3*mx1-1]; j++) {
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
/* post receives */
      ierr = MPI_Irecv(mcll,ncsize,mint,kl,itg[0],lgrp,&msid[0]);
      ierr = MPI_Irecv(mclr,ncsize,mint,kr,itg[1],lgrp,&msid[1]);
      ierr = MPI_Irecv(rbufl,nbsize,mreal,kl,itg[2],lgrp,&msid[2]);
      ierr = MPI_Irecv(rbufr,nbsize,mreal,kr,itg[3],lgrp,&msid[3]);
/* send particle number offsets */
      ierr = MPI_Isend(nclr,ncsize,mint,kr,itg[0],lgrp,&msid[4]);
      ierr = MPI_Isend(ncll,ncsize,mint,kl,itg[1],lgrp,&msid[5]);
      ierr = MPI_Wait(&msid[0],&istatus);
      ierr = MPI_Wait(&msid[1],&istatus);
/* send particles */
      jsr = idimp*nclr[3*mx1-1];
      ierr = MPI_Isend(sbufr,jsr,mreal,kr,itg[2],lgrp,&msid[6]);
      jsl = idimp*ncll[3*mx1-1];
      ierr = MPI_Isend(sbufl,jsl,mreal,kl,itg[3],lgrp,&msid[7]);
      ierr = MPI_Wait(&msid[2],&istatus);
      ierr = MPI_Wait(&msid[3],&istatus);
   }
/* make sure sbufr, sbufl, ncll, and nclr have been sent */
   if (nvp != 1) {
      for (i = 0; i < 4; i++) {
         ierr = MPI_Wait(&msid[i+4],&istatus);
      }
   }
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
void cpppmove2_(float *sbufr, float *sbufl, float *rbufr, float *rbufl,
                int *ncll, int *nclr, int *mcll, int *mclr, int *kstrt,
                int *nvp, int *idimp, int *nbmax, int *mx1) {
   cpppmove2(sbufr,sbufl,rbufr, rbufl,ncll,nclr,mcll,mclr,*kstrt,*nvp,
             *idimp,*nbmax,*mx1);
   return;
}
