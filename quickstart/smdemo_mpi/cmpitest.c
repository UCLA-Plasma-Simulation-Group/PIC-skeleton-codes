#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>
#include "mpimain.h"

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
/* mpi vector add test program      */
/* written by Viktor K. Decyk, UCLA */
/* nx = size of array */
   int nx = 1048576;
/* idproc = processor id               */
/* nvp = number of processors obtained */
   int idproc, nvp;
   int j, nxp, nxps, irc;
   float time = 0.0;
   float eps, epsmax;
   double dtime;
   float *a = NULL, *b = NULL, *c = NULL, *d = NULL; 
   struct timeval itime;
   float *p_a = NULL, *p_b = NULL, *p_c = NULL; 

/* set up mpi */
   irc = 0;
   init_mpi(&idproc,&nvp,&irc,argc,argv);
   if (irc != 0) {
      printf("%i MPI initialization error!\n",idproc);
      exit(1);
   }
   if (idproc==0)
      printf("mpi nodes available = %i\n",nvp);
/* initialize Host data */
   if (idproc==0) {
      a = (float *) malloc(nx*sizeof(float));
      b = (float *) malloc(nx*sizeof(float));
      c = (float *) malloc(nx*sizeof(float));
      d = (float *) malloc(nx*sizeof(float));
/* initialize vectors */
      for (j = 0; j < nx; j++) {
         a[j] = 0.0; b[j] = j + 1; c[j] = 2*j + 2; d[j] = -1.0;
      }
   }
/* nxp = maximum size of array on each processor */
   nxp = (nx - 1)/nvp + 1;
/* allocate data for mpi */
   p_a = (float *) malloc(nxp*sizeof(float));
   p_b = (float *) malloc(nxp*sizeof(float));
   p_c = (float *) malloc(nxp*sizeof(float));
/* nxps = size of actual data, if nx is not an exact multiple of nvp */
   nxps = 0 > nx-nxp*idproc ? 0 : nx-nxp*idproc;
   nxps = nxp < nxps ? nxp : nxps;

/* First execute on Host in C: a = b + c */
   if (idproc==0) {
      dtimer(&dtime,&itime,-1);
      cadd(a,b,c,nx);
      dtimer(&dtime,&itime,1);
      printf("C add time=%e\n",(float)dtime);
   }
/* Copy data for mpi */
   dtimer(&dtime,&itime,-1);
   vscatter(b,p_b,nx,nxp);
   vscatter(c,p_c,nx,nxp);
   dtimer(&dtime,&itime,1);
   time += (float)dtime;
   if (idproc==0)
      printf("Copyin time=%e\n",(float)dtime);
/* Execute with mpi: g_a = g_b + g_c */
   dtimer(&dtime,&itime,-1);
   mpadd(p_a,p_b,p_c,nxps);
   dtimer(&dtime,&itime,1);
   time += (float)dtime;
   if (idproc==0)
      printf("mpi add time=%e\n",(float)dtime);
/* Copy data from mpi: d = g_a */
   dtimer(&dtime,&itime,-1);
   vgather(d,p_a,nx,nxp);
   dtimer(&dtime,&itime,1);
   time += (float)dtime;
   if (idproc==0) {
      printf("Copyout time=%e\n",(float)dtime);
      printf("Total mpi time=%e\n",time);
   }

/* Check for correctness: compare a and d */
   if (idproc==0) {
      epsmax = 0.0;
      for (j = 0; j < nx; j++) {
         eps = a[j] - d[j];
         if (eps < 0.0)
         eps = -eps;
         if (eps > epsmax)
            epsmax = eps;
      }
      printf("maximum difference = %e\n",epsmax);
   }

/* deallocate memory for mpi */
   free(p_a);
   free(p_b);
   free(p_c);
/* close down mpi */
   end_mpi();
/* deallocate Host memory */
   if (idproc==0) {
      free(a);
      free(b);
      free(c);
      free(d);
   }

   return 0;
}

