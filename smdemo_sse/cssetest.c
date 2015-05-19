#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>
#include "sselib.h"

void cadd(float a[], float b[], float c[], int nx) {
/* C add procedure */
   int j;
   for (j = 0; j < nx; j++) {
      a[j] = b[j] + c[j];
   }
   return;
}

void dtimer(double *time, struct timeval *itime, int icntrl);

int main(int argc, char *argv[])
{
/* SSE vector add test program */
/* written by Viktor K. Decyk, UCLA */
/* nx = size of array */
   int nx = 1048576;
   int j, irc = 0;
   float eps, epsmax;
/* timing data */
   double dtime;
   struct timeval itime;
/* data for C */
   float *a = NULL, *b = NULL, *c = NULL; 
/* data for SSE */
   float *s_a = NULL, *s_b = NULL, *s_c = NULL; 

/* initialize data for C */
   a = (float *) malloc(nx*sizeof(float));
   b = (float *) malloc(nx*sizeof(float));
   c = (float *) malloc(nx*sizeof(float));
/* initialize vectors */
   for (j = 0; j < nx; j++) {
      a[j] = 0.0; b[j] = j + 1; c[j] = 2*j + 2;
   }
/* allocate aligned 1d array for SSE */
   sse_fallocate(&s_a,nx,&irc);
   sse_fallocate(&s_b,nx,&irc);
   sse_fallocate(&s_c,nx,&irc);
   if (irc != 0) {
      printf("SSE allocate error!\n");
      exit(1);
   }
/* Copy initial data for SSE */
   for (j = 0; j < nx; j++) {
      s_b[j] = b[j]; s_c[j] = c[j];
   }

/* First execute in C: a = b + c */
   dtimer(&dtime,&itime,-1);
   cadd(a,b,c,nx);
   dtimer(&dtime,&itime,1);
   printf("C add time=%e\n",(float)dtime);

/* Execute on SSE: s_a = s_b + s_c */
   dtimer(&dtime,&itime,-1);
   ssadd(s_a,s_b,s_c,nx);
   dtimer(&dtime,&itime,1);
   printf("SSE add time=%e\n",(float)dtime);

/* Check for correctness: compare a and s_a */
   epsmax = 0.0;
   for (j = 0; j < nx; j++) {
      eps = a[j] - s_a[j];
      if (eps < 0.0)
         eps = -eps;
      if (eps > epsmax)
         epsmax = eps;
   }
   printf("maximum difference = %e\n",epsmax);

/* deallocate memory for SSE */
   sse_deallocate(s_a);
   sse_deallocate(s_b);
   sse_deallocate(s_c);
/* deallocate C memory */
   free(a);
   free(b);
   free(c);

   return 0;
}

