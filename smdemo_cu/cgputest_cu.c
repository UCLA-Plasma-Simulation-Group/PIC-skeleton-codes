#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>
#include "gpumain_cu.h"

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
   float time = 0.0;
   float eps, epsmax;
   double dtime;
   float *a = NULL, *b = NULL, *c = NULL, *d = NULL; 
   struct timeval itime;
   float *g_a = NULL, *g_b = NULL, *g_c = NULL; 

/* initialize Host data */
   a = (float *) malloc(nx*sizeof(float));
   b = (float *) malloc(nx*sizeof(float));
   c = (float *) malloc(nx*sizeof(float));
   d = (float *) malloc(nx*sizeof(float));
/* initialize vectors */
   for (j = 0; j < nx; j++) {
      a[j] = 0.0; b[j] = j + 1; c[j] = 2*j + 2; d[j] = -1.0;
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
   g_fallocate(&g_a,nx,&irc);
   g_fallocate(&g_b,nx,&irc);
   g_fallocate(&g_c,nx,&irc);
   if (irc != 0) {
      printf("GPU allocate error!\n");
      exit(1);
   }

/* First execute on Host in C: a = b + c */
   dtimer(&dtime,&itime,-1);
   cadd(a,b,c,nx);
   dtimer(&dtime,&itime,1);
   printf("C add time=%e\n",(float)dtime);
/* Copy data to GPU */
   dtimer(&dtime,&itime,-1);
   copyin_gmemptr(b,g_b,nx);
   copyin_gmemptr(c,g_c,nx);
   dtimer(&dtime,&itime,1);
   time += (float)dtime;
   printf("Copyin time=%e\n",(float)dtime);
/* Execute on GPU: g_a = g_b + g_c */
   dtimer(&dtime,&itime,-1);
   gpadd(g_a,g_b,g_c,nx);
   dtimer(&dtime,&itime,1);
   time += (float)dtime;
   printf("GPU add time=%e\n",(float)dtime);
/* Copy data from GPU: d = g_a */
   dtimer(&dtime,&itime,-1);
   copyout_gmemptr(d,g_a,nx);
   dtimer(&dtime,&itime,1);
   time += (float)dtime;
   printf("Copyout time=%e\n",(float)dtime);
   printf("Total GPU time=%e\n",time);

/* Check for correctness: compare a and d */
   epsmax = 0.0;
   for (j = 0; j < nx; j++) {
      eps = a[j] - d[j];
      if (eps < 0.0)
         eps = -eps;
      if (eps > epsmax)
         epsmax = eps;
   }
   printf("maximum difference = %e\n",epsmax);

/* deallocate memory on GPU */
   g_deallocate(&g_a,&irc);
   g_deallocate(&g_b,&irc);
   g_deallocate(&g_c,&irc);
/* close down GPU */
   end_cu();
/* deallocate Host memory */
   free(a);
   free(b);
   free(c);
   free(d);

   return 0;
}

