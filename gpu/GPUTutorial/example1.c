/*--------------------------------------------------------------------*/
/* CUDA C GPU Tutorial: Copy            */
/* written by Viktor K. Decyk, UCLA */

#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>
#include "copy.h"
#include "gpulib2.h"

void dtimer(double *time, struct timeval *itime, int icntrl);

int main(int argc, char *argv[]) {
/* nx, ny = size of array */
   int nx = 3000, ny = 600;
/* nblock = block size on GPU */
   int nblock = 64;
/* mx, my = data block size */
   int mx, my;
   float s = 0.5;
   int j, k, irc;
   float eps, epsmax;
/* timing data */
   double dtime;
   struct timeval itime;
/* data for C Host */
   float *a1 = NULL, *b1 = NULL;
   float *a2 = NULL, *b2 = NULL;
/* data for GPU */
   float *g_a1 = NULL, *g_b1 = NULL;
   float *g_a2 = NULL, *g_b2 = NULL;
/* allocate Host data */
   a1 = (float *) malloc(nx*sizeof(float));
   b1 = (float *) malloc(nx*sizeof(float));
   a2 = (float *) malloc(nx*ny*sizeof(float));
   b2 = (float *) malloc(nx*ny*sizeof(float));

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
   gpu_fallocate(&g_b1,nx,&irc);

/* allocate 2d data on GPU */
   gpu_fallocate(&g_a2,nx*ny,&irc);
   gpu_fallocate(&g_b2,nx*ny,&irc);

   if (irc != 0) {
      printf("GPU allocate error!\n");
      exit(1);
   }

/* initialize 1d data on Host */
   for (j = 0; j < nx; j++) {
      b1[j] = (float) (j + 1);
      a1[j] = 0.0;
   }

/* initialize 2d data on Host */
   for (k = 0; k < ny; k++) {
      for (j = 0; j < nx; j++) {
         b2[j+nx*k] = (float) (j + nx*k + 1);
         a2[j+nx*k] = 0.0;
      }
   }
/* copy data to GPU */
   gpu_fcopyin(a1,g_a1,nx);
   gpu_fcopyin(b1,g_b1,nx);
   gpu_fcopyin(a2,g_a2,nx*ny);
   gpu_fcopyin(b2,g_b2,nx*ny);

/* measure overhead time by running empty kernel */
   dtimer(&dtime,&itime,-1);
   emptykernel();
   dtimer(&dtime,&itime,1);
   printf("C empty kernel time=%e\n",(float)dtime);

/* 1D copy */
   mx = 128;

/* segmented 1d copy on host with block size mx */
   dtimer(&dtime,&itime,-1);
/* copy0(a1,b1,nx); */
   copy1(a1,b1,mx,nx);
   dtimer(&dtime,&itime,1);
   printf("C 1d copy time=%e\n",(float)dtime);

/* 1d copy on GPU with block size mx */
   dtimer(&dtime,&itime,-1);
   gpu_copy1(g_a1,g_b1,mx,nx);
   dtimer(&dtime,&itime,1);
   printf("GPU 1d copy time=%e\n",(float)dtime);

/* copy data from GPU */
   gpu_fcopyout(b1,g_a1,nx);

/* Check for correctness: compare a1 and g_a1 */
   epsmax = 0.0;
   for (j = 0; j < nx; j++) {
      eps = a1[j] - b1[j];
      if (eps < 0.0)
         eps = -eps;
      if (eps > epsmax)
         epsmax = eps;
   }
   printf("1d copy maximum difference = %e\n",epsmax);

/* 2D copy */
   mx = 16; my = 16;

/* segmented 2d copy on host with block size mx, my */
   dtimer(&dtime,&itime,-1);
/* copy2(a2,b2,mx,nx,ny);    */
/* saxpy2(a2,b2,s,mx,nx,ny); */
   copy3(a2,b2,mx,my,nx,ny);
   dtimer(&dtime,&itime,1);
   printf("C 2d copy time=%e\n",(float)dtime);

/* 2d copy on GPU with block size mx, my */
   dtimer(&dtime,&itime,-1);
/* gpu_copy2a(g_a2,g_b2,mx,nx,ny);   */
/* gpu_copy2b(g_a2,g_b2,mx,nx,ny);   */
/* gpu_saxpy2(g_a2,g_b2,s,mx,nx,ny); */
   gpu_copy3(g_a2,g_b2,mx,my,nx,ny);
   dtimer(&dtime,&itime,1);
   printf("GPU 2d copy time=%e\n",(float)dtime);

/* copy data from GPU */
   gpu_fcopyout(b2,g_a2,nx*ny);

/* Check for correctness: compare a2 and g_a2 */
   epsmax = 0.0;
   for (k = 0; k < ny; k++) {
      for (j = 0; j < nx; j++) {
         eps = a2[j+nx*k] - b2[j+nx*k];
         if (eps < 0.0)
            eps = -eps;
         if (eps > epsmax)
            epsmax = eps;
      }
   }
   printf("2d copy maximum difference = %e\n",epsmax);

/* deallocate memory on GPU */
   gpu_deallocate((void *)g_a1,&irc);
   gpu_deallocate((void *)g_b1,&irc);
   gpu_deallocate((void *)g_a2,&irc);
   gpu_deallocate((void *)g_b2,&irc);
   if (irc != 0) {
      printf("GPU deallocate error!\n");
      exit(1);
   }
/* close down GPU */
   end_cu();

   return 0;
}

