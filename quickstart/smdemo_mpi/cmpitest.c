#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>
#include "mpilib.h"

void cadd(float a[], float b[], float c[], int nx) {
   int j;
   for (j = 0; j < nx; j++) {
      a[j] = b[j] + c[j];
   }
   return;
}

void dtimer(double *time, struct timeval *itime, int icntrl);

int main(int argc, char *argv[])
{
/* MPI vector add test program      */
/* written by Viktor K. Decyk, UCLA */
/* nx = size of array */
   int nx = 1048576;
/* idproc = processor id               */
/* nvp = number of processors obtained */
   int idproc, nvp;
   int j, joff, nxp, nxps, irc;
   float eps, epsmax;
   float epsm[1], work[1];
/* timing data */
   double dtime;
   struct timeval itime;
/* global data for C */
   float *a = NULL, *b = NULL, *c = NULL; 
/* local data for MPI */
   float *p_a = NULL, *p_b = NULL, *p_c = NULL; 

/* initialize global Host data */
   a = (float *) malloc(nx*sizeof(float));
   b = (float *) malloc(nx*sizeof(float));
   c = (float *) malloc(nx*sizeof(float));
/* initialize the same global vectors on each MPI node */
   for (j = 0; j < nx; j++) {
      a[j] = 0.0; b[j] = j + 1; c[j] = 2*j + 2;
   }
/* set up MPI */
/* idproc contains local processor id, different for each MPI node */
   irc = 0;
   init_mpi(&idproc,&nvp,&irc,argc,argv);
   if (irc != 0) {
      printf("%i MPI initialization error!\n",idproc);
      exit(1);
   }
   if (idproc==0)
      printf("MPI nodes available = %i\n",nvp);

/* vector data will be partitioned into local arrays among MPI nodes */
/* nxp = maximum size of local array on each MPI node                */
   nxp = (nx - 1)/nvp + 1;
/* allocate local data for MPI */
   p_a = (float *) malloc(nxp*sizeof(float));
   p_b = (float *) malloc(nxp*sizeof(float));
   p_c = (float *) malloc(nxp*sizeof(float));
/* joff = starting location in global array for local array     */
/* different value on each MPI node, since it depends on idproc */
   joff = nxp*idproc;
/* nxps = actual size of local array, if nx not an exact multiple of nvp */
/* possibly different value on each MPI node, since it depends on joff   */
   nxps = 0 > nx-joff ? 0 : joff;
   nxps = nxp < nxps ? nxp : nxps;
/* Copy part of initial vector data to local arrays for MPI */
   for (j = 0; j < nxps; j++) {
      p_b[j] = b[j+joff];
      p_c[j] = c[j+joff];
   }
/* alternative copy data for MPI:                                      */
/* copy from MPI node 0 to other nodes if only node 0 has initial data */
/* vscatter(b,p_b,nx,nxp); */
/* vscatter(c,p_c,nx,nxp); */

/* First execute on Host in C: a = b + c               */
/* Each processor executes the same add of global data */
   dtimer(&dtime,&itime,-1);
   cadd(a,b,c,nx);
   dtimer(&dtime,&itime,1);
   if (idproc==0)
      printf("C add time=%e\n",(float)dtime);

/* Execute with MPI: p_a = p_b + p_c       */
/* Each processor adds only its local data */
   dtimer(&dtime,&itime,-1);
   mpadd(p_a,p_b,p_c,nxps);
   dtimer(&dtime,&itime,1);
   if (idproc==0)
      printf("MPI add time=%e\n",(float)dtime);

/* Check for correctness: compare a and p_a                     */
/* Each MPI node compares global result with local array result */
   epsmax = 0.0;
   for (j = 0; j < nxps; j++) {
      eps = a[j+joff] - p_a[j];
      if (eps < 0.0)
         eps = -eps;
      if (eps > epsmax)
          epsmax = eps;
   }
/* find global maximum error for each of the local maximum errors */
   epsm[0] = epsmax;
   ppmax(epsm,work,1);
   epsmax = epsm[0];
   if (idproc==0)
      printf("maximum difference = %e\n",epsmax);

/* Alternate check for correctness: compare a and p_a                 */
/* copy to c from other MPI nodes to node 0 if only node 0 has answer */
/* vgather(c,p_a,nx,nxp);                         */
/* if (idproc==0) {                               */
/*    epsmax = 0.0;                               */
/*    for (j = 0; j < nx; j++) {                  */
/*       eps = a[j] - c[j];                       */
/*       if (eps < 0.0)                           */
/*          eps = -eps;                           */
/*       if (eps > epsmax)                        */
/*           epsmax = eps;                        */
/*    }                                           */
/*    printf("maximum difference = %e\n",epsmax); */
/* }                                              */

/* deallocate memory for MPI */
   free(p_a);
   free(p_b);
   free(p_c);
/* close down MPI */
   end_mpi();
/* deallocate Host memory */
   free(a);
   free(b);
   free(c);

   return 0;
}

