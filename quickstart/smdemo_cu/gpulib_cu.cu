/* GPU vector add test program for CUDA */
/* written by Viktor K. Decyk, UCLA */

#include <stdlib.h>
#include <stdio.h>
#include "cuda.h"

static int nblock_size = 64;
static int ngrid_size = 1;
int maxgsx = 65535;
int mmcc = 0;
static int devid;

static cudaError_t crc;

/*--------------------------------------------------------------------*/
__global__ void gadd(float a[], float b[], float c[], int nx) {
   int j;
   j = threadIdx.x+blockDim.x*blockIdx.x;
   if (j < nx)
      a[j] = b[j] + c[j];
   return;
}

/*--------------------------------------------------------------------*/
extern "C" void gpadd(float *a, float *b, float *c, int nx) {
/* Vector Add Interface for C */
   dim3 dimBlock(nblock_size);
   dim3 dimGrid((nx - 1)/nblock_size + 1);
   gadd<<<dimGrid,dimBlock>>>(a,b,c,nx);
   cudaThreadSynchronize();
   return;
}

__global__ void emptyKernel() {}

/*--------------------------------------------------------------------*/
extern "C" void emptykernel() {
   int ngx, ngy;
   ngx  = nblock_size < 32768 ? nblock_size : 32768;
   ngy = (ngrid_size - 1)/ngx + 1;
   dim3 dimBlock(nblock_size,1);
   dim3 dimGrid(ngx,ngy);
   crc = cudaGetLastError();
   emptyKernel<<<dimGrid,dimBlock>>>();
   cudaThreadSynchronize();
   crc = cudaGetLastError();
   if (crc) {
      printf("emptyKernel error=%d:%s\n",crc,cudaGetErrorString(crc));
      exit(1);
   }
   return;
}

/*--------------------------------------------------------------------*/
extern "C" void init_cu(int dev, int *irc) {
/* initialize CUDA with device dev or selects best GPU available       */
/* searches throughs devices, selects the device with the most compute */
/* units, and saves the device id devid                                */
/* if dev is a valid device, it is used, otherwise the GPU with the    */
/* most multi-processors is selected                                   */
/* error code is modified only if there is an error */
   int maxcpus = 0, jm = -1;
   int j, ndevs, maxunits;
   unsigned long msize;
   double z;
   struct cudaDeviceProp prop;
/* returns number of device */
   crc = cudaGetDeviceCount(&ndevs);
   if (crc) {
      printf("cudaGetDeviceCount Error=%i:%s\n",crc,
             cudaGetErrorString(crc));
      *irc = 1;
      return;
   }
/* get information about devices */
   for (j = 0; j < ndevs; j++) {
      crc = cudaGetDeviceProperties(&prop,j);
      if (crc) {
         printf("cudaGetDeviceProperties Error=%i:%s\n",crc,
                cudaGetErrorString(crc));
         prop.name[0] = 0;
      }
      maxunits = prop.multiProcessorCount;
      if (dev <= 0) {
         printf("j=%i:CUDA_DEVICE_NAME=%s,CUDA_MULTIPROCESSOR_COUNT=%i\n",
                j,prop.name,maxunits);
         msize = prop.totalGlobalMem;
         z = ((double) msize)/1073741824.0;
         mmcc = 10*prop.major + prop.minor;
         printf("    CUDA_GLOBAL_MEM_SIZE=%lu(%f GB),Capability=%d\n",
                msize,(float) z,mmcc);
         printf("    Capability=%d\n",mmcc);
         if (maxunits > maxcpus) {
            maxcpus = maxunits;
            jm = j;
         }
      }
   }
   devid = jm;
   if (dev >= 0)
      devid = dev % ndevs;
   printf("using device j=%i\n",devid);
/* get properties for this device */
   crc = cudaGetDeviceProperties(&prop,devid);
   maxgsx = prop.maxGridSize[0];
   mmcc = 10*prop.major + prop.minor;
/* set device */
   crc = cudaSetDevice(devid);
   if (crc) {
      printf("cudaSetDevice Error=%i:%s\n",crc,
             cudaGetErrorString(crc));
      *irc = 1;
      return;
   }
/* run empty kernel */
   emptykernel();
   return;
}

/*--------------------------------------------------------------------*/
extern "C" void end_cu() {
   crc = cudaThreadExit();
   if (crc) {
      printf("cudaThreadExit Error=%d:%s\n",crc,cudaGetErrorString(crc));
   }
   return;
}

/*--------------------------------------------------------------------*/
extern "C" void setgbsize(int nblock) {
   nblock_size = nblock;
   return;
}

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
extern "C" void gpu_deallocate(float **g_f, int *irc) {
/* deallocate global memory on GPU, return pointer to C */
   crc = cudaFree((void *)*g_f);
   if (crc) {
      printf("cudaFree Error=%d:%s\n",crc,cudaGetErrorString(crc));
      *irc = 1;
   }
   *g_f = NULL;
   return;
}

/*--------------------------------------------------------------------*/
extern "C" void gpu_fcopyin(float *f, float *g_f, int nsize) {
/* copy float array from host memory to global GPU memory */
   crc = cudaMemcpy((void *)g_f,f,sizeof(float)*nsize,
                    cudaMemcpyHostToDevice);
   if (crc) {
      printf("cudaMemcpyHostToDevice float Error=%d:%s\n",crc,
              cudaGetErrorString(crc));
      exit(1);
   }
   return;
}

/*--------------------------------------------------------------------*/
extern "C" void gpu_fcopyout(float *f, float *g_f, int nsize) {
/* copy float array from global GPU memory to host memory */
   crc = cudaMemcpy(f,(void *)g_f,sizeof(float)*nsize,
                    cudaMemcpyDeviceToHost);
   if (crc) {
      printf("cudaMemcpyDeviceToHost float Error=%d:%s\n",crc,
              cudaGetErrorString(crc));
      exit(1);
   }
   return;
}

/* Interfaces to Fortran */

/*--------------------------------------------------------------------*/
extern "C" void gpadd_(unsigned long *gp_a, unsigned long *gp_b,
                       unsigned long *gp_c, int *nx) {
/* Vector Add Interface for Fortran */
   float *a, *b, *c;
   a = (float *)*gp_a;
   b = (float *)*gp_b;
   c = (float *)*gp_c;
   gpadd(a,b,c,*nx);
}

/*--------------------------------------------------------------------*/
extern "C" void init_cu_(int *dev, int *irc) {
   init_cu(*dev,irc);
   return;
}

/*--------------------------------------------------------------------*/
extern "C" void end_cu_() {
   end_cu();
   return;
}

/*--------------------------------------------------------------------*/
extern "C" void setgbsize_(int *nblock) {
   setgbsize(*nblock);
   return;
}

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
extern "C" void gpu_deallocate_(unsigned long *gp_f, int *irc) {
/* deallocate global memory on GPU, return pointer to Fortran */
   float *f;
   f = (float *)*gp_f;
   gpu_deallocate(&f,irc);
   *gp_f = 0;
   return;
}

/*--------------------------------------------------------------------*/
extern "C" void gpu_fcopyin_(float *f, unsigned long *gp_f,
                             int *nsize) {
/* copy float array from main memory to global GPU memory */
   float *g_f;
   g_f = (float *)*gp_f;
   gpu_fcopyin(f,g_f,*nsize);
   return;
}

/*--------------------------------------------------------------------*/
extern "C" void gpu_fcopyout_(float *f, unsigned long *gp_f,
                              int *nsize) {
/* copy float array from global GPU memory to main memory */
   float *g_f;
   g_f = (float *)*gp_f;
   gpu_fcopyout(f,g_f,*nsize);
   return;
}
