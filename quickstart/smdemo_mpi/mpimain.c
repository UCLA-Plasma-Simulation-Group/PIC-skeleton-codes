/* vector add test program for MPI  */
/* written by Viktor K. Decyk, UCLA */

#include <stdlib.h>
#include <stdio.h>
#include "mpi.h"

/* nproc = number of real or virtual processors obtained */
static int nproc;

/* lgrp = current communicator                           */
/* lworld = MPI_COMM_WORLD communicator                  */
static MPI_Comm lgrp, lworld;

/* mreal = default datatype for reals                    */
/* mint = default datatype for integers                  */
/* mcplx = default datatype for complex type             */
/* mdouble = default double precision type               */
static MPI_Datatype mreal, mint, mcplx, mdouble;

/* msum = MPI_SUM */
/* mmax = MPI_MAX */
static MPI_Op msum, mmax;
      
/*--------------------------------------------------------------------*/
void mpadd(float a[], float b[], float c[], int nx) {
   int j;
   for (j = 0; j < nx; j++) {
      a[j] = b[j] + c[j];
   }
   return;
}

/*--------------------------------------------------------------------*/
void init_mpi(int *idproc, int *nvp, int *irc, int argc, char *argv[]) {
/* this subroutine initializes parallel processing     */
/* lgrp communicator = MPI_COMM_WORLD                  */
/* output: idproc, nvp                                 */
/* idproc = processor id in lgrp communicator          */
/* nvp = number of real or virtual processors obtained */
/* error code is modified only if there is an error    */
/* local data */
/* ndprec = (0,1) = (no,yes) use (normal,autodouble) precision         */
/* idprec = (0,1) = (no,yes) use (normal,autodouble) integer precision */
   int ndprec = 0;
/* int idprec = 0; */
   int ierror, flag;
/* indicate whether MPI_INIT has been called */
   ierror = MPI_Initialized(&flag);
   if (!flag) {
/* initialize the MPI execution environment */
      ierror = MPI_Init(&argc,&argv);
      if (ierror) {
         *irc = 1;
         return;
      }
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
void vscatter(float f[],float g[], int nx, int nxp) {
/* this subroutine performs a scatter of a real vector, that is, */
/* successive nxp blocks of source f on processor 0, are sent to */
/* destination g on successive processors:                       */
/* g[0:nxp-1] = f[nxp*idproc:min(nxp*(idproc+1),nx)-1]           */
/* f = input real data                                           */
/* g = output real data                                          */
/* nx = size of source data                                      */
/* nxp = size of destination data                                */
/* local data */
   int is, j, idproc, nvp, nxps, koff, ierr;
   MPI_Status istatus;
/* find processor id */
   ierr = MPI_Comm_rank(lgrp,&idproc);
/* find number of processors */
   ierr = MPI_Comm_size(lgrp,&nvp);
/* copy directly on processor 0 */
   if (idproc==0) {
      for (j = 0; j < nxp; j++) {
         g[j] = f[j];
      }
/* send to remaining processors */
      for (is = 1; is < nvp; is++) {
         nxps = nxp*is;
         koff = nxps < nx-1 ? nxps : nx-1;
         nxps = 0 > nx-nxps ? 0 : nx-nxps;
         nxps = nxp < nxps ? nxp : nxps;
         ierr = MPI_Send(&f[koff],nxps,mreal,is,1,lgrp);
      }
   }
/* receive from processor 0 */
   else {
      nxps = 0 > nx-nxp*idproc ? 0 : nx-nxp*idproc;
      nxps = nxp < nxps ? nxp : nxps;
      ierr = MPI_Recv(g,nxps,mreal,0,1,lgrp,&istatus);
   }
   return;
}

/*--------------------------------------------------------------------*/
void vgather(float f[],float g[], int nx, int nxp) {
/* this subroutine performs a gather of a real vector, that is,        */
/* successive nxp blocks of destination f on processor 0, are received */
/* from source g on successive processors:                             */
/* f[nxp*idproc:min(nxp*(idproc+1),nx)-1] = g[0:nxp-1]                 */
/* f = input real data                                                 */
/* g = output real data                                                */
/* nx = size of destination data                                       */
/* nxp = size of source data                                           */
/* local data */
   int is, j, idproc, nvp, nxps, koff, ierr;
   MPI_Status istatus;
/* find processor id */
   ierr = MPI_Comm_rank(lgrp,&idproc);
/* find number of processors */
   ierr = MPI_Comm_size(lgrp,&nvp);
/* copy directly on processor 0 */
   if (idproc==0) {
      for (j = 0; j < nxp; j++) {
         f[j] = g[j];
      }
/* receive from remaining processors */
      for (is = 1; is < nvp; is++) {
         nxps = nxp*is;
         koff = nxps < nx-1 ? nxps : nx-1;
         nxps = 0 > nx-nxps ? 0 : nx-nxps;
         nxps = nxp < nxps ? nxp : nxps;
         ierr = MPI_Recv(&f[koff],nxps,mreal,is,2,lgrp,&istatus);
      }
   }
/* send to processor 0 */
   else {
      nxps = 0 > nx-nxp*idproc ? 0 : nx-nxp*idproc;
      nxps = nxp < nxps ? nxp : nxps;
      ierr = MPI_Send(g,nxps,mreal,0,2,lgrp);
   }
   return;
}
 
/*--------------------------------------------------------------------*/
void end_mpi() {
/* this subroutine terminates parallel processing */
   int ierror, flag;
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

