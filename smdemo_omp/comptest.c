#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>
#include "omplib.h"

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
/* OpenMP vector add test program */
/* written by Viktor K. Decyk, UCLA */
/* nx = size of array, nthreads = number of threads */
   int nx = 1048576, nthreads = 1;
   int j;
   float time = 0.0;
   float eps, epsmax;
/* timing data */
   double dtime;
   struct timeval itime;
/* data for C */
   float *a = NULL, *b = NULL, *c = NULL; 
/* data for OpenMP */
   float *p_a = NULL, *p_b = NULL, *p_c = NULL; 

/* initialize Host data */
   a = (float *) malloc(nx*sizeof(float));
   b = (float *) malloc(nx*sizeof(float));
   c = (float *) malloc(nx*sizeof(float));
/* initialize vectors */
   for (j = 0; j < nx; j++) {
      a[j] = 0.0; b[j] = j + 1; c[j] = 2*j + 2;
   }
/* set up OpenMP */
   init_omp(0);
/* setnthsize(nthreads); */
/* allocate data on OpenMP */
   p_a = (float *) malloc(nx*sizeof(float));
   p_b = (float *) malloc(nx*sizeof(float));
   p_c = (float *) malloc(nx*sizeof(float));
/* Copy initial data for OpenMP */
   for (j = 0; j < nx; j++) {
      p_b[j] = b[j];
      p_c[j] = c[j];
   }

/* First execute on Host in C: a = b + c */
   dtimer(&dtime,&itime,-1);
   cadd(a,b,c,nx);
   dtimer(&dtime,&itime,1);
   printf("C add time=%e\n",(float)dtime);

/* Execute with OpenMP: s_a = s_b + s_c */
   dtimer(&dtime,&itime,-1);
   mpadd(p_a,p_b,p_c,nx);
   dtimer(&dtime,&itime,1);
   time += (float)dtime;
   printf("OpenMP add time=%e\n",(float)dtime);

/* Check for correctness: compare a and p_a */
   epsmax = 0.0;
   for (j = 0; j < nx; j++) {
      eps = a[j] - p_a[j];
      if (eps < 0.0)
         eps = -eps;
      if (eps > epsmax)
         epsmax = eps;
   }
   printf("maximum difference = %e\n",epsmax);

/* deallocate memory for OpenMP */
   free(p_a);
   free(p_b);
   free(p_c);
/* deallocate Host memory */
   free(a);
   free(b);
   free(c);

   return 0;
}

