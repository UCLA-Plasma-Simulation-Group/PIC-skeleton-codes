#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>
#include "gpulib_cu.h"

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
/* GPU vector add test program */
/* written by Viktor K. Decyk, UCLA */
/* nx = size of array, nblock = block size on GPU */
   int nx = 1048576, nblock = 64;
   int j, irc;
   float eps, epsmax;
/* timing data */
   double dtime;
   struct timeval itime;
/* data for C */
   float *a = NULL, *b = NULL, *c = NULL; 
/* data for GPU */
   float *g_a = NULL, *g_b = NULL, *g_c = NULL; 

/* initialize C data on Host */
   a = (float *) malloc(nx*sizeof(float));
   b = (float *) malloc(nx*sizeof(float));
   c = (float *) malloc(nx*sizeof(float));
/* initialize vectors */
   for (j = 0; j < nx; j++) {
      a[j] = 0.0; b[j] = j + 1; c[j] = 2*j + 2;
   }
/* set up GPU */
   irc = 0;
   setgbsize(nblock);
   init_cu(0,&irc);
   if (irc != 0) {
      printf("CUDA initialization error!\n");
      exit(1);
   }
/* allocate data on GPU, return address to C */
   gpu_fallocate(&g_a,nx,&irc);
   gpu_fallocate(&g_b,nx,&irc);
   gpu_fallocate(&g_c,nx,&irc);
   if (irc != 0) {
      printf("GPU allocate error!\n");
      exit(1);
   }
/* Copy initial data to GPU */
   gpu_fcopyin(b,g_b,nx);
   gpu_fcopyin(c,g_c,nx);

/* First execute on Host in C: a = b + c */
   dtimer(&dtime,&itime,-1);
   cadd(a,b,c,nx);
   dtimer(&dtime,&itime,1);
   printf("C add time=%e\n",(float)dtime);

/* Execute on GPU: g_a = g_b + g_c */
   dtimer(&dtime,&itime,-1);
   gpadd(g_a,g_b,g_c,nx);
   dtimer(&dtime,&itime,1);
   printf("GPU add time=%e\n",(float)dtime);

/* Check for correctness: compare a and g_a */
/* Copy g_a from GPU to c on Host, then compare a with c */
   gpu_fcopyout(c,g_a,nx);
   epsmax = 0.0;
   for (j = 0; j < nx; j++) {
      eps = a[j] - c[j];
      if (eps < 0.0)
         eps = -eps;
      if (eps > epsmax)
         epsmax = eps;
   }
   printf("maximum difference = %e\n",epsmax);

/* deallocate memory on GPU */
   gpu_deallocate(&g_a,&irc);
   gpu_deallocate(&g_b,&irc);
   gpu_deallocate(&g_c,&irc);
/* close down GPU */
   end_cu();
/* deallocate Host memory */
   free(a);
   free(b);
   free(c);

   return 0;
}

