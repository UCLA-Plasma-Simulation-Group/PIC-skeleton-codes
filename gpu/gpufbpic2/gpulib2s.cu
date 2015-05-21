/*--------------------------------------------------------------------*/
/* CUDA special utility Library */
/* written by Viktor K. Decyk, UCLA */

#include <stdlib.h>
#include <stdio.h>
#include "cuda.h"

static cudaError_t crc;

/*--------------------------------------------------------------------*/
extern "C" void gpu_fallocate(float **g_f, int nsize, int *irc) {
/* allocate global float memory on GPU, return pointer to C */
   void *gptr;
   crc = cudaMalloc(&gptr,sizeof(float)*nsize);
   if (crc) {
      printf("cudaMalloc float Error=%d:%s,l=%d\n",crc,
              cudaGetErrorString(crc),nsize);
      *irc = 1;
   }
   *g_f = (float *)gptr;
   return;
}

/*--------------------------------------------------------------------*/
extern "C" void gpu_iallocate(int **g_i, int nsize, int *irc) {
/* allocate global integer memory on GPU, return pointer to C */
   void *gptr;
   crc = cudaMalloc(&gptr,sizeof(int)*nsize);
   if (crc) {
      printf("cudaMalloc int Error=%d:%s,l=%d\n",crc,
              cudaGetErrorString(crc),nsize);
      *irc = 1;
   }
   *g_i = (int *)gptr;
   return;
}

/*--------------------------------------------------------------------*/
extern "C" void gpu_callocate(float2 **g_c, int nsize, int *irc) {
/* allocate global float2 memory on GPU, return pointer to C */
   void *gptr;
   crc = cudaMalloc(&gptr,sizeof(float2)*nsize);
   if (crc) {
      printf("cudaMalloc float2 Error=%d:%s,l=%d\n",crc,
              cudaGetErrorString(crc),nsize);
      *irc = 1;
   }
   *g_c = (float2 *)gptr;
   return;
}

/*--------------------------------------------------------------------*/
extern "C" void gpu_deallocate(void *g_d, int *irc) {
/* deallocate global memory on GPU */
   crc = cudaFree(g_d);
   if (crc) {
      printf("cudaFree Error=%d:%s\n",crc,cudaGetErrorString(crc));
      *irc = 1;
   }
   return;
}

/* Interfaces to Fortran */

/*--------------------------------------------------------------------*/
extern "C" void gpu_fallocate_(unsigned long *gp_f, int *nsize,
                               int *irc) {
/* allocate global float memory on GPU, return pointer to Fortran */
   float *fptr;
   gpu_fallocate(&fptr,*nsize,irc);
   *gp_f = (long )fptr;
   return;
}

/*--------------------------------------------------------------------*/
extern "C" void gpu_iallocate_(unsigned long *gp_i, int *nsize,
                               int *irc) {
/* allocate global integer memory on GPU, return pointer to Fortran */
   int *iptr;
   gpu_iallocate(&iptr,*nsize,irc);
   *gp_i = (long )iptr;
   return;
}

/*--------------------------------------------------------------------*/
extern "C" void gpu_callocate_(unsigned long *gp_f, int *nsize,
                               int *irc) {
/* allocate global float2 memory on GPU, return pointer */
/* to Fortran */
   float2 *fptr;
   gpu_callocate(&fptr,*nsize,irc);
   *gp_f = (long )fptr;
   return;
}

/*--------------------------------------------------------------------*/
extern "C" void gpu_deallocate_(unsigned long *gp_d, int *irc) {
/* deallocate global memory on GPU, return pointer to Fortran */
   void *d;
   d = (void *)*gp_d;
   gpu_deallocate(d,irc);
   *gp_d = 0;
   return;
}


