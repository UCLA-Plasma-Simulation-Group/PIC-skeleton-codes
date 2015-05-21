/*--------------------------------------------------------------------*/
/* Basic parallel PIC library for MPI communications
   pplib2.c contains basic communications procedures for 1d partitions
   cppinit2 initializes parallel processing for C, returns
            number of processors and processor id.
   cppfndgrp finds which MPI nodes have the same host name, creates
             ordered list of such nodes, and returns where current node
             is located in that list 
   cppexit terminates parallel processing.
   cppabort aborts parallel processing.
   cpwtimera performs parallel local wall clock timing.
   cppsum performs parallel sum of a real vector.
   cppdsum performs parallel sum of a double precision vector.
   cppimax performs parallel maximum of an integer vector.
   cpppncguard2l sends/receives guard cells in y for scalar array, linear
                 interpolation, and distributed data with non-uniform
                 partition.  for copying guard cells.
   cpppnaguard2l sends/receives guard cells in y for scalar array, linear
                 interpolation, and distributed data with non-uniform
                 partition.  for adding guard cells.
   cppptpose performs a transpose of a complex scalar array, distributed
             in y, to a complex scalar array, distributed in x.
             optimized for GPU
   cppptposen performs a transpose of an n component complex vector array,
              distributed in y, to an n component complex vector array,
              distributed in x.  optimized for GPU with vector data
   cacsndrec helps perform asynchronous transpose between GPUS
   cpppmove2 moves particles into appropriate spatial regions for tiled
             distributed data.
   written by viktor k. decyk, ucla
   copyright 1995, regents of the university of california
   update: may 7, 2014                                         */

#include <stdlib.h>
#include <stdio.h>
#include <complex.h>
#include <string.h>
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
void cppfndgrp(int locl[], int kstrt, int nvp, int *idev, int *ndev) {
/* this subroutine finds which MPI nodes have the same host name
   creates ordered list of such nodes andreturns where current node is
   located in that list 
   used to ensure that different GPUs on the same host have different ids
   input: all except locl, idev, ndev, output: locl, idev, ndev
   locl = ordered list of ranks on same host
   kstrt = starting data block number
   nvp = number of real or virtual processors
   idev = location in rank list with value kstrt - 1 (-1 if not found)
   ndev = number of nodes found on host
local data */
   int j, l, m, n, nn, ks, id, lname, ierr;
   MPI_Request msid;
   MPI_Status istatus;
   char *hname = NULL, *name = NULL;
   lname = MPI_MAX_PROCESSOR_NAME;
   hname = (char *) malloc(lname*sizeof(char));
   name = (char *) malloc(lname*sizeof(char));
   ks = kstrt - 1;
   ierr =  MPI_Get_processor_name(hname,&l);
   nn = -1;
/* this segment is used for mpi computers  */
/* find and save ranks with same host name */
   for (n = 0; n < nvp; n++) {
      id = n - ks;
      if (id < 0)
         id += nvp;
      if (id != ks) {
/* post receive*/
         ierr = MPI_Irecv(name,lname,MPI_CHAR,id,n,lgrp,&msid);
/* send data */
         ierr = MPI_Send(hname,lname,MPI_CHAR,id,n,lgrp);
/* receive data */
         ierr = MPI_Wait(&msid,&istatus);
      }
      else {
         strcpy(name,hname);
      }
/* save rank if remote name equals local name */
      if (!strcmp(name,hname)) {
         nn += 1;
         locl[nn] = id;
      }
   }
   nn += 1;
/* order rank list */
   for (j = 0; j < (nn-1); j++) {
      l = j;
      m = locl[j];
/* find minimum value and location */
      for (n = j; n < nn; n++) {
         id = locl[n];
         if (id < m) {
            m = id;
            l = n;
         }
      }
/* swap minimum to beginning of array */
      if (l > j) {
         id = locl[j];
         locl[j] = m;
         locl[l] = id;
      }
   }
/* find location in rank list with value kstrt - 1 */
   *idev = -1;
   for (j = 0; j < nn; j++) {
      if (ks==locl[j])
         *idev = j;
   }
/* return number of nodes found on host */
   *ndev = nn;
   free(hname);
   free(name);
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
void cpppncguard2l(float scs[], float scr[], int kstrt, int nvp,
                   int nxv) {
/* this subroutine sends/receives guard cells in y for scalar array
   for copying guard cells, sends to left processor, receives from right
   output: scr
   scs[j] = input data to send
   scr[j] = received data
   kstrt = starting data block number
   nvp = number of real or virtual processors
   nxv = size of array to send
   linear interpolation, for distributed data
local data */
   int ks, moff, kl, kr, ierr;
   MPI_Request msid;
   MPI_Status istatus;
   ks = kstrt - 1;
   moff = nxv*nvp + 2;
/* copy guard cells */
   kr = ks + 1;
   if (kr >= nvp)
      kr = kr - nvp;
   kl = ks - 1;
   if (kl < 0)
      kl = kl + nvp;
/* this segment is used for mpi computers */
   ierr = MPI_Irecv(scr,nxv,mreal,kr,moff,lgrp,&msid);
   ierr = MPI_Send(scs,nxv,mreal,kl,moff,lgrp);
   ierr = MPI_Wait(&msid,&istatus);
   return;
}

/*--------------------------------------------------------------------*/
void cpppnaguard2l(float scs[], float scr[], int kstrt, int nvp,
                   int nxv) {
/* this subroutine sends/receives guard cells in y for scalar array
   for adding guard cells, sends to right processor, receives from left
   output: scr
   scs[j] = input data to send
   scr[j] = received data
   kstrt = starting data block number
   nvp = number of real or virtual processors
   nxv = size of array to send
   linear interpolation, for distributed data
local data */
   int ks, moff, kl, kr, ierr;
   MPI_Request msid;
   MPI_Status istatus;
   ks = kstrt - 1;
   moff = nxv*nvp + 1;
/* add guard cells */
   kr = ks + 1;
   if (kr >= nvp)
      kr = kr - nvp;
   kl = ks - 1;
   if (kl < 0)
      kl = kl + nvp;
/* this segment is used for mpi computers */
   ierr = MPI_Irecv(scr,nxv,mreal,kl,moff,lgrp,&msid);
   ierr = MPI_Send(scs,nxv,mreal,kr,moff,lgrp);
   ierr = MPI_Wait(&msid,&istatus);
   return;
}

/*--------------------------------------------------------------------*/
void cppptpose(float complex sm[], float complex tm[], int nx, int ny,
               int kxp, int kyp, int kstrt, int nvp) {
/* this subroutine sends and receives data between MPI nodes to perform
   a transpose of a matrix distributed in y, to another matrix
   distributed in x.  one message is sent and received at a time.
   optimized for GPU
   ss/tm = complex buffers on host to be sent/received
   nx/ny = number of points in x/y
   kxp/kyp = number of data values per block in x/y
   kstrt = starting data block number
   nvp = number of real or virtual processors
local data */
   int j, n, nn, ks, kyps, kxyp, id, joff, ld, ierr;
   MPI_Request msid;
   MPI_Status istatus;
   ks = kstrt - 1;
   kyps = ny - kyp*ks;
   kyps = 0 > kyps ? 0 : kyps;
   kyps = kyp < kyps ? kyp : kyps;
   kxyp = kxp*kyp;
/* special case for one processor */
   if (nvp==1) {
      for (j = 0; j < kxyp; j++) {
         tm[j] = sm[j];
      }
      return;
   }
   nn = -1;
/* this segment is used for mpi computers */
   for (n = 0; n < nvp; n++) {
      id = n - ks;
      if (id < 0)
         id += nvp;
      if (id != ks) {
/* adjust counter to omit data sent to oneself */
         nn += 1;
/* calculate length of data to send */
         joff = kxp*id;
         ld = nx - joff;
         ld = 0 > ld ? 0 : ld;
         ld = kyps*(kxp < ld ? kxp : ld);
/* post receive */
         ierr = MPI_Irecv(&tm[kxyp*nn],kxyp,mcplx,id,n+1,lgrp,&msid);
/* send data */
         ierr = MPI_Send(&sm[kxyp*nn],ld,mcplx,id,n+1,lgrp);
/* receive data */
         ierr = MPI_Wait(&msid,&istatus);
      }
   }
   return;
}

/*--------------------------------------------------------------------*/
void cppptposen(float complex sm[], float complex tm[], int nx, int ny,
                int kxp, int kyp, int kstrt, int nvp, int ndim) {
/* this subroutine sends and receives data between MPI nodes to perform
   a transpose of an n component matrix distributed in y, to another 
   n component matrix distributed in x.
   one message is sent and received at a time.
   optimized for GPU with vector data
   ss/tm = complex buffers on host to be sent/received
   nx/ny = number of points in x/y
   kxp/kyp = number of data values per block in x/y
   kstrt = starting data block number
   nvp = number of real or virtual processors
local data */
   int j, n, nn, ks, kyps, kxyp, id, joff, ld, ierr;
   MPI_Request msid;
   MPI_Status istatus;
   ks = kstrt - 1;
   kyps = ny - kyp*ks;
   kyps = 0 > kyps ? 0 : kyps;
   kyps = ndim*(kyp < kyps ? kyp : kyps);
   kxyp = kxp*ndim*kyp;
/* special case for one processor */
   if (nvp==1) {
      for (j = 0; j < kxyp; j++) {
         tm[j] = sm[j];
      }
      return;
   }
   nn = -1;
/* this segment is used for mpi computers */
   for (n = 0; n < nvp; n++) {
      id = n - ks;
      if (id < 0)
         id += nvp;
      if (id != ks) {
/* adjust counter to omit data sent to oneself */
         nn += 1;
/* calculate length of data to send */
         joff = kxp*id;
         ld = nx - joff;
         ld = 0 > ld ? 0 : ld;
         ld = kyps*(kxp < ld ? kxp : ld);
/* post receive */
         ierr = MPI_Irecv(&tm[kxyp*nn],kxyp,mcplx,id,n+1,lgrp,&msid);
/* send data */
         ierr = MPI_Send(&sm[kxyp*nn],ld,mcplx,id,n+1,lgrp);
/* receive data */
         ierr = MPI_Wait(&msid,&istatus);
      }
   }
   return;
}

/*--------------------------------------------------------------------*/
void cacsndrec(float complex stm[], int idproc, int nsize, int ntag,
               int mode) {
/* this subroutine is part of a family of procedures for performing
   a transpose by sending and receiving asynchronous data between GPUs
   stm = receive/send buffer
   idproc = processor id for sending/receiving
   nsize = size of data packet in words
   ntag = MPI tag
   mode = (1,2,3) = (post receive, post send, wait for send/receive)
   modes 1 and 2 should be called before mode 3 is called.
   local data */
   int ierr;
   static MPI_Request msid, mrid;
   static MPI_Status istatus;
   if (mode==1) {
      ierr = MPI_Irecv(stm,nsize,mcplx,idproc,ntag,lgrp,&mrid);
   }
   else if (mode==2) {
      ierr = MPI_Isend(stm,nsize,mcplx,idproc,ntag,lgrp,&msid);
   }
   else if (mode==3) {
      ierr = MPI_Wait(&msid,&istatus);
      ierr = MPI_Wait(&mrid,&istatus);
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

/*--------------------------------------------------------------------*/
void cppfndgrp_(int *locl, int *kstrt, int *nvp, int *idev, int *ndev) {
   cppfndgrp(locl,*kstrt,*nvp,idev,ndev);
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
void cpppncguard2l_(float *scs, float *scr, int *kstrt, int *nvp,
                    int *nxv) {
   cpppncguard2l(scs,scr,*kstrt,*nvp,*nxv);
   return;
}

/*--------------------------------------------------------------------*/
void cpppnaguard2l_(float *scs, float *scr, int *kstrt, int *nvp,
                    int *nxv) {
   cpppnaguard2l(scs,scr,*kstrt,*nvp,*nxv);
   return;
}

/*--------------------------------------------------------------------*/
void cppptpose_(float complex *sm, float complex *tm, int *nx, int *ny,
                int *kxp, int *kyp, int *kstrt, int *nvp) {
   cppptpose(sm,tm,*nx,*ny,*kxp,*kyp,*kstrt,*nvp);
   return;
}

/*--------------------------------------------------------------------*/
void cppptposen_(float complex *sm, float complex *tm, int *nx, int *ny,
                 int *kxp, int *kyp, int *kstrt, int *nvp, int *ndim) {
   cppptposen(sm,tm,*nx,*ny,*kxp,*kyp,*kstrt,*nvp,*ndim);
   return;
}

/*--------------------------------------------------------------------*/
void cacsndrec_(float complex *stm, int *idproc, int *nsize, int *ntag,
                int *mode) {
   cacsndrec(stm,*idproc,*nsize,*ntag,*mode); 
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
