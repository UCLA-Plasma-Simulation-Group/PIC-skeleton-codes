#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>
#include "pthlib.h"

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
/* Pthreads vector add test program */
/* written by Viktor K. Decyk, UCLA */
/* nx = size of array, nthreads = number of threads */
   int nx = 1048576, nthreads = 1;
   int j, irc;
   float eps, epsmax;
/* timing data */
   double dtime;
   struct timeval itime;
/* data for C */
   float *a = NULL, *b = NULL, *c = NULL; 
/* data for Pthreads */
   float *p_a = NULL, *p_b = NULL, *p_c = NULL; 

/* initialize Host data */
   a = (float *) malloc(nx*sizeof(float));
   b = (float *) malloc(nx*sizeof(float));
   c = (float *) malloc(nx*sizeof(float));
/* initialize vectors */
   for (j = 0; j < nx; j++) {
      a[j] = 0.0; b[j] = j + 1; c[j] = 2*j + 2;
   }
/* set up Pthreads */
   irc = 0;
   init_pt(0,&irc);
   if (irc != 0) {
      printf("Pthreads initialization error!\n");
      exit(1);
   }
/* setnthsize(nthreads); */
/* allocate data on pthreads */
   p_a = (float *) malloc(nx*sizeof(float));
   p_b = (float *) malloc(nx*sizeof(float));
   p_c = (float *) malloc(nx*sizeof(float));
/* Copy initial data for Pthreads */
   for (j = 0; j < nx; j++) {
      p_b[j] = b[j];
      p_c[j] = c[j];
   }

/* First execute on Host in C: a = b + c */
   dtimer(&dtime,&itime,-1);
   cadd(a,b,c,nx);
   dtimer(&dtime,&itime,1);
   printf("C add time=%e\n",(float)dtime);

/* Execute with Pthreads: p_a = p_b + p_c */
   dtimer(&dtime,&itime,-1);
   ptadd(p_a,p_b,p_c,nx,&irc);
   dtimer(&dtime,&itime,1);
   if (irc != 0)
      printf("ptadd error: irc=%i\n",irc);
   printf("pthreads add time=%e\n",(float)dtime);

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

/* deallocate memory for Pthreads */
   free(p_a);
   free(p_b);
   free(p_c);
/* close down Pthreads */
   end_pt();
/* deallocate Host memory */
   free(a);
   free(b);
   free(c);

   return 0;
}

