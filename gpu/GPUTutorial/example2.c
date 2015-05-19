/*--------------------------------------------------------------------*/
/* C GPU Tutorial: Transpose        */
/* written by Viktor K. Decyk, UCLA */

#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>
#include "transpose.h"
#include "gpulib2.h"

void dtimer(double *time, struct timeval *itime, int icntrl);

int main(int argc, char *argv[]) {
   int nx = 512, ny = 512, mx = 16, my = 16;
   int nblock = 64;
   int j, k, irc;
   float eps, epsmax;
   double dtime;
   struct timeval itime;
   float *a2 = NULL, *b2 = NULL, *c2 = NULL;
   float *g_a2 = NULL, *g_b2 = NULL;
/* allocate host data */
   a2 = (float *) malloc(ny*nx*sizeof(float));
   b2 = (float *) malloc(nx*ny*sizeof(float));
   c2 = (float *) malloc(ny*nx*sizeof(float));

/* set up GPU */
   irc = 0;
   setgbsize(nblock);
   init_cu(0,&irc);
   if (irc != 0) {
      printf("CUDA initialization error!\n");
      exit(1);
   }

/* allocate 2d data on GPU */
   gpu_fallocate(&g_a2,ny*nx,&irc);
   gpu_fallocate(&g_b2,nx*ny,&irc);

   if (irc != 0) {
      printf("GPU allocate error!\n");
      exit(1);
   }

/* initialize 2d data on host */
   for (k = 0; k < ny; k++) {
      for (j = 0; j < nx; j++) {
         b2[j+nx*k] = (float) (j + nx*k + 1);
         a2[k+ny*j] = 0.0;
      }
   }
   gpu_fcopyin(a2,g_a2,ny*nx);

/* measure overhead time by running empty kernel */
   dtimer(&dtime,&itime,-1);
   emptykernel();
   dtimer(&dtime,&itime,1);
   printf("C empty kernel time=%e\n",(float)dtime);

/* segmented 2d transpose on host with block size mx, my */
   dtimer(&dtime,&itime,-1);
/* transpose0(a2,b2,nx,ny); */
   transpose2(a2,b2,mx,my,nx,ny);
   dtimer(&dtime,&itime,1);
   printf("C 2d transpose time=%e\n",(float)dtime);

/* 2d transpose on GPU with block size mx, mx */
   gpu_fcopyin(b2,g_b2,nx*ny);
   dtimer(&dtime,&itime,-1);
   gpu_transpose2(g_a2,g_b2,mx,nx,ny);
   dtimer(&dtime,&itime,1);
   printf("GPU 2d transpose time=%e\n",(float)dtime);
   gpu_fcopyout(c2,g_a2,nx*ny);

/* Check for correctness: compare a2 and g_a2 */
   epsmax = 0.0;
   for (k = 0; k < ny; k++) {
      for (j = 0; j < nx; j++) {
         eps = a2[j+nx*k] - c2[j+nx*k];
         if (eps < 0.0)
            eps = -eps;
         if (eps > epsmax)
            epsmax = eps;
      }
   }
   printf("2d transpose maximum difference = %e\n",epsmax);

/* deallocate memory on GPU */
   gpu_deallocate((void *)g_a2,&irc);
   gpu_deallocate((void *)g_b2,&irc);
/* close down GPU */
   end_cu();

   return 0;
}

