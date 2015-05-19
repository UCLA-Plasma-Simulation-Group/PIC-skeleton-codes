/*--------------------------------------------------------------------*/
/* C GPU Tutorial: Reduction        */
/* written by Viktor K. Decyk, UCLA */

#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>
#include "redux.h"
#include "gpulib2.h"

void dtimer(double *time, struct timeval *itime, int icntrl);

int main(int argc, char *argv[]) {
   int nx = 3000, mx = 128;
   int nblock = 64;
   int j, nbx, nbxs, irc;
   float s, t, eps, epsmax;
   double dtime;
   struct timeval itime;
   float *a1 = NULL, *d1 = NULL;
   float *g_a1 = NULL, *g_d1 = NULL, *g_s = NULL;
   nbx = (nx - 1)/mx + 1;
   nbxs = (nbx - 1)/mx + 1;

/* allocate host data */
   a1 = (float *) malloc(nx*sizeof(float));
   d1 = (float *) malloc(nbx*sizeof(float));

/* set up GPU */
   irc = 0;
   setgbsize(nblock);
   init_cu(0,&irc);
   if (irc != 0) {
      printf("CUDA initialization error!\n");
      exit(1);
   }

/* allocate 1d data on GPU */
   gpu_fallocate(&g_a1,nx,&irc);
   gpu_fallocate(&g_d1,nbx,&irc);
   gpu_fallocate(&g_s,2*nbxs,&irc);

   if (irc != 0) {
      printf("GPU allocate error!\n");
      exit(1);
   }

/* initialize 1d data on host */
   for (j = 0; j < nx; j++) {
      a1[j] = (float) (j + 1);
   }
   for (j = 0; j < nbx; j++) {
      d1[j] = 0.0;
   }
   gpu_fcopyin(d1,g_d1,nbx);
   s = 0.0; t = 0.0;

/* measure overhead time by running empty kernel */
   dtimer(&dtime,&itime,-1);
   emptykernel();
   dtimer(&dtime,&itime,1);
   printf("C empty kernel time=%e\n",(float)dtime);

/* segmented 1d sum on host with block size mx */
   dtimer(&dtime,&itime,-1);
/* sum0(a1,&s,nx);     */
/* sum1(a1,&s,mx,nx);  */
   sum2(a1,d1,mx,nx);
   sum1(d1,&s,mx,nbx);
   dtimer(&dtime,&itime,1);
   printf("C 1d sum time=%e\n",(float)dtime);

/* 1d copy on GPU with block size mx */
   gpu_fcopyin(a1,g_a1,nx);
   dtimer(&dtime,&itime,-1);
/* gpu_sum1(g_a1,g_s,mx,nx);  */
/* gpu_sum2(g_a1,g_d1,mx,nx); */
/* gpu_sum1(g_d1,g_s,mx,nbx); */
   gpu_sum3(g_a1,g_d1,g_s,mx,nx);
   dtimer(&dtime,&itime,1);
   printf("GPU 1d sum time=%e\n",(float)dtime);
   gpu_fcopyout(a1,g_d1,nbx);
   gpu_fcopyout(&t,g_s,1);

/* Check for correctness: compare d1 and g_d1 */
   epsmax = 0.0;
   for (j = 0; j < nbx; j++) {
      eps = d1[j] - a1[j];
      if (eps < 0.0)
         eps = -eps;
      if (eps > epsmax)
         epsmax = eps;
   }
   printf("1d sum maximum difference = %e\n",epsmax);
   printf("s,t = %f,%f\n",s,t);

/* deallocate memory on GPU */
   gpu_deallocate((void *)g_a1,&irc);
   gpu_deallocate((void *)g_d1,&irc);
   gpu_deallocate((void *)g_s,&irc);
   if (irc != 0) {
      printf("GPU deallocate error!\n");
      exit(1);
   }
/* close down GPU */
   end_cu();

   return 0;
}

