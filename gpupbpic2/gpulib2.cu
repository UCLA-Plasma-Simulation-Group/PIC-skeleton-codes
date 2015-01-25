/* CUDA utility Library */
/* written by Viktor K. Decyk, UCLA */

#include <stdlib.h>
#include <stdio.h>
#include "cuda.h"

int nblock_size = 64;
int ngrid_size = 1;
int maxgsx = 65535;
int mmcc = 0;
static int devid;

static cudaError_t crc;

#define MAXSTREAMS             4
static cudaStream_t streams[MAXSTREAMS] = {NULL,NULL,NULL,NULL};

/* Prototypes for Fortran function called by C */
extern "C" void getfcptr_(unsigned long *carrayref, float *carray,
                          int *nx);

extern "C" void getf2cptr_(unsigned long *carrayref, float *carray,
                           int *nx, int *ny);

extern "C" void getc2cptr_(unsigned long *carrayref, float2 *carray,
                           int *nx, int *ny);

__global__ void emptyKernel() {}

/*--------------------------------------------------------------------*/
extern "C" void gpu_setgbsize(int nblock) {
/* set blocksize */
   nblock_size = nblock;
   return;
}

/*--------------------------------------------------------------------*/
extern "C" int getmmcc() {
/* get major and minor computer capability */
   return mmcc;
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

/*--------------------------------------------------------------------*/
extern "C" void hpl_fallocate(float **h_f, int nsize, int *irc) {
/* allocate page-locked float memory on host, return pointer to C */
   void *hptr = NULL;
   crc = cudaMallocHost(&hptr,sizeof(float)*nsize);
   if (crc) {
      printf("cudaMallocHost float Error=%d:%s,l=%d\n",crc,
              cudaGetErrorString(crc),nsize);
      *irc = 1;
   }
   *h_f = (float *)hptr;
   return;
}

/*--------------------------------------------------------------------*/
extern "C" void hpl_callocate(float2 **h_c, int nsize, int *irc) {
/* allocate page-locked float2 memory on host, return pointer to C */
   void *hptr = NULL;
   crc = cudaMallocHost(&hptr,sizeof(float2)*nsize);
   if (crc) {
      printf("cudaMallocHost float2 Error=%d:%s,l=%d\n",crc,
              cudaGetErrorString(crc),nsize);
      *irc = 1;
   }
   *h_c = (float2 *)hptr;
   return;
}

/*--------------------------------------------------------------------*/
extern "C" void hpl_deallocate(void *h_d, int *irc) {
/* deallocate page-locked on host */
   crc = cudaFreeHost(h_d);
   if (crc) {
      printf("cudaFreeHost Error=%d:%s\n",crc,cudaGetErrorString(crc));
      *irc = 1;
   }
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

/*--------------------------------------------------------------------*/
extern "C" void gpu_icopyin(int *f, int *g_f, int nsize) {
/* copy int array from host memory to global GPU memory */
   crc = cudaMemcpy((void *)g_f,f,sizeof(int)*nsize,
                    cudaMemcpyHostToDevice);
   if (crc) {
      printf("cudaMemcpyHostToDevice int Error=%d:%s\n",crc,
              cudaGetErrorString(crc));
      exit(1);
   }
   return;
}

/*--------------------------------------------------------------------*/
extern "C" void gpu_icopyout(int *f, int *g_f, int nsize) {
/* copy int array from global GPU memory to host memory */
   crc = cudaMemcpy(f,(void *)g_f,sizeof(int)*nsize,
                    cudaMemcpyDeviceToHost);
   if (crc) {
      printf("cudaMemcpyDeviceToHost int Error=%d:%s\n",crc,
              cudaGetErrorString(crc));
      exit(1);
   }
   return;
}

/*--------------------------------------------------------------------*/
extern "C" void gpu_ccopyin(float2 *f, float2 *g_f, int nsize) {
/* copy float2 array from host memory to global GPU memory */
   crc = cudaMemcpy((void *)g_f,f,sizeof(float2)*nsize,
                    cudaMemcpyHostToDevice);
   if (crc) {
      printf("cudaMemcpyHostToDevice float2 Error=%d:%s\n",crc,
              cudaGetErrorString(crc));
      exit(1);
   }
   return;
}

/*--------------------------------------------------------------------*/
extern "C" void gpu_ccopyout(float2 *f, float2 *g_f, int nsize) {
/* copy float2 array from global GPU memory to host memory */
   crc = cudaMemcpy(f,(void *)g_f,sizeof(float2)*nsize,
                    cudaMemcpyDeviceToHost);
   if (crc) {
      printf("cudaMemcpyDeviceToHost float2 Error=%d:%s\n",crc,
              cudaGetErrorString(crc));
      exit(1);
   }
   return;
}

/*--------------------------------------------------------------------*/
extern "C" void gpu_initstream(int nstream) {
/* Create Stream for requested identifier nstream       */
/* nstream should be between 1 and MAXSTREAMS inclusive */
   if ((nstream < 1) || (nstream > MAXSTREAMS)) {
      printf("gpu_initstream: nstream out of bounds = %d\n",nstream);
      exit(1);
   }
   if (streams[nstream-1] != NULL) {
      printf("gpu_initstream: nstream already used = %d\n",nstream);
      exit(1);
   }
   crc = cudaStreamCreate(&streams[nstream-1]);
   if (crc) {
      printf("cudaStreamCreate Error=%d:%s\n",crc,
             cudaGetErrorString(crc));
      exit(1);
   }
   return;
}

/*--------------------------------------------------------------------*/
extern "C" void gpu_delstream(int nstream) {
/* Destroy Stream for requested identifier nstream      */
/* nstream should be between 1 and MAXSTREAMS inclusive */
   if ((nstream < 1) || (nstream > MAXSTREAMS)) {
      printf("gpu_delstream: nstream out of bounds = %d\n",nstream);
   }
   if (streams[nstream-1] == NULL) {
      printf("gpu_delstream: nstream not allocated = %d\n",nstream);
   }
   crc = cudaStreamDestroy(streams[nstream-1]);
   if (crc) {
      printf("cudaStreamDestroy Error=%d:%s\n",crc,
             cudaGetErrorString(crc));
      exit(1);
   }
   return;
}

/*--------------------------------------------------------------------*/
extern "C" void gpu_waitstream(int nstream) {
/* Synchronize Stream for requested identifier nstream  */
/* nstream should be between 0 and MAXSTREAMS inclusive */
   cudaStream_t stream = NULL;
   if ((nstream >= 0) || (nstream <= MAXSTREAMS)) {
      if (nstream > 0) stream = streams[nstream-1];
   }
   else {
      printf("gpu_waitstream: nstream undefined = %d\n",nstream);
      exit(1);
   }
   crc = cudaStreamSynchronize(stream);
   if (crc) {
      printf("cudaStreamSynchronize Error=%d:%s\n",crc,
             cudaGetErrorString(crc));
      exit(1);
   }
   return;
}

/*--------------------------------------------------------------------*/
extern "C" void gpu_cascopyin(float2 *f, float2 *g_f, int noff, 
                              int nsize, int nstream) {
/* copy float2 array segment from host memory to global GPU memory */
/* asynchronous copy */
   float2 *cptr;
   cudaStream_t stream = NULL;
   cptr = &g_f[noff];
   if ((nstream >= 0) || (nstream <= MAXSTREAMS)) {
      if (nstream > 0) stream = streams[nstream-1];
   }
   else {
      printf("gpu_cascopyin: nstream undefined = %d\n",nstream);
      exit(1);
   }
   crc = cudaMemcpyAsync((void *)cptr,f,sizeof(float2)*nsize,
                         cudaMemcpyHostToDevice,stream);
   if (crc) {
      printf("Async cudaMemcpyHostToDevice float2 Error=%d:%s\n",crc,
              cudaGetErrorString(crc));
      exit(1);
   }
   return;
}

/*--------------------------------------------------------------------*/
extern "C" void gpu_cascopyout(float2 *f, float2 *g_f, int noff,
                               int nsize, int nstream) {
/* copy float2 array segment from global GPU memory to host memory */
/* asynchronous copy */
   float2 *cptr;
   cudaStream_t stream = NULL;
   cptr = &g_f[noff];
   if ((nstream >= 0) || (nstream <= MAXSTREAMS)) {
      if (nstream > 0) stream = streams[nstream-1];
   }
   else {
      printf("gpu_cascopyout: nstream undefined = %d\n",nstream);
      exit(1);
   }
   crc = cudaMemcpyAsync(f,(void *)cptr,sizeof(float2)*nsize,
                         cudaMemcpyDeviceToHost,stream);
   if (crc) {
      printf("Async cudaMemcpyDeviceToHost float2 Error=%d:%s\n",crc,
              cudaGetErrorString(crc));
      exit(1);
   }
   return;
}

/*--------------------------------------------------------------------*/
extern "C" void gpu_zfmem(float *g_f, int nsize) {
/* initialize float array in global GPU memory to zero */
   crc = cudaMemset((void *)g_f,0,sizeof(float)*nsize);
   if (crc) {
      printf("cudaMemset Error=%d:%s\n",crc,cudaGetErrorString(crc));
      exit(1);
   }
   return;
}

/*--------------------------------------------------------------------*/
extern "C" void gpu_zcmem(float2 *g_f, int nsize) {
/* initialize float2 array in global GPU memory to zero */
   crc = cudaMemset((void *)g_f,0,sizeof(float2)*nsize);
   if (crc) {
      printf("cudaMemset Error=%d:%s\n",crc,cudaGetErrorString(crc));
      exit(1);
   }
   return;
}

/*--------------------------------------------------------------------*/
extern "C" void gpu_set_cache_size(int nscache) {
/* request preferred cache size, requires CUDA 3.2 or higher */
/* nscache = (0,1,2) = (no,small,big) cache size */
   cudaFuncCache cpref;
   if ((nscache < 0) || (nscache > 2))
      return;
   if (nscache==0)
      cpref = cudaFuncCachePreferNone;
   else if (nscache==1)
      cpref = cudaFuncCachePreferShared;
   else if (nscache==2)
      cpref = cudaFuncCachePreferL1;
   crc = cudaThreadSetCacheConfig(cpref);
/* crc = cudaDeviceSetCacheConfig(cpref); */
   if (crc) {
      printf("cudaThreadSetCacheConfig error=%d:%s\n",crc,
             cudaGetErrorString(crc));
   }
   return;
}

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

extern "C" void end_cu() {
/* terminate CUDA */
   crc = cudaThreadExit();
   if (crc) {
      printf("cudaThreadExit Error=%d:%s\n",crc,cudaGetErrorString(crc));
   }
   return;
}

/* Interfaces to Fortran */

/*--------------------------------------------------------------------*/
extern "C" void gpu_setgbsize_(int *nblock) {
   gpu_setgbsize(*nblock);
   return;
}

/*--------------------------------------------------------------------*/
extern "C" int getmmcc_() {
/* get major and minor computer capability */
   return getmmcc();
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

/*--------------------------------------------------------------------*/
extern "C" void hpl_f1allocate_(unsigned long *hp_f, int *nx,
                                int *irc) {
/* allocate page-locked 1d real memory on host, assign */
/* data pointer to Fortran pointer object hp_f         */ 
/* This procedure needs an interface in Fortran90      */
/* interface                                 */
/*    subroutine hpl_f1allocate(hp_f,nx,irc) */
/*    implicit none                          */
/*    integer :: nx, irc                     */
/*    real, dimension(:), pointer :: hp_f    */
/*    end subroutine                         */
/* end interface                             */
   int nsize;
   float *fptr;
   nsize = *nx;
/* allocate data on host */
   hpl_fallocate(&fptr,nsize,irc);
/* set reference to C data in real Fortran pointer object */
   getfcptr_(hp_f,fptr,nx);
   return;
}

/*--------------------------------------------------------------------*/
extern "C" void hpl_f2allocate_(unsigned long *hp_f, int *nx, int *ny,
                                int *irc) {
/* allocate page-locked 2d real memory on host, assign */
/* data pointer to Fortran pointer object hp_f         */ 
/* This procedure needs an interface in Fortran90      */
/* interface                                    */
/*    subroutine hpl_f2allocate(hp_f,nx,ny,irc) */
/*    implicit none                             */
/*    integer :: nx, ny, irc                    */
/*    real, dimension(:,:), pointer :: hp_f     */
/*    end subroutine                            */
/* end interface                                */
   int nsize;
   float *fptr;
   nsize = (*nx)*(*ny);
/* allocate data on host */
   hpl_fallocate(&fptr,nsize,irc);
/* set reference to C data in real Fortran pointer object */
   getf2cptr_(hp_f,fptr,nx,ny);
   return;
}

/*--------------------------------------------------------------------*/
extern "C" void hpl_c2allocate_(unsigned long *hp_f, int *nx, int *ny,
                                int *irc) {
/* allocate page-locked 2d complex memory on host, assign */
/* data pointer to Fortran pointer object hp_f            */ 
/* This procedure needs an interface in Fortran90         */
/* interface                                    */
/*    subroutine hpl_c2allocate(hp_f,nx,ny,irc) */
/*    implicit none                             */
/*    integer :: nx, ny, irc                    */
/*    complex, dimension(:,:), pointer :: hp_f  */
/*    end subroutine                            */
/* end interface                                */
   int nsize;
   float2 *cptr;
   nsize = (*nx)*(*ny);
/* allocate data on host */
   hpl_callocate(&cptr,nsize,irc);
/* set reference to C data in complex Fortran pointer object */
   getc2cptr_(hp_f,cptr,nx,ny);
   return;
}

/*--------------------------------------------------------------------*/
extern "C" void hpl_deallocate_(void *h_d, int *irc) {
/* deallocate page-locked memory on host                  */
/* pointer in Fortran should also be nullified            */
   hpl_deallocate(h_d,irc);
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

/*--------------------------------------------------------------------*/
extern "C" void gpu_icopyin_(int *f, unsigned long *gp_f, int *nsize) {
/* copy int array from main memory to global GPU memory */
   int *g_f;
   g_f = (int *)*gp_f;
   gpu_icopyin(f,g_f,*nsize);
   return;
}

/*--------------------------------------------------------------------*/
extern "C" void gpu_icopyout_(int *f, unsigned long *gp_f, int *nsize) {
/* copy int array from global GPU memory to main memory */
   int *g_f;
   g_f = (int *)*gp_f;
   gpu_icopyout(f,g_f,*nsize);
   return;
}

/*--------------------------------------------------------------------*/
extern "C" void gpu_ccopyin_(float2 *f, unsigned long *gp_f,
                             int *nsize) {
/* copy float2 array from main memory to global GPU memory */
   float2 *g_f;
   g_f = (float2 *)*gp_f;
   gpu_ccopyin(f,g_f,*nsize);
   return;
}

/*--------------------------------------------------------------------*/
extern "C" void gpu_ccopyout_(float2 *f, unsigned long *gp_f,
                              int *nsize) {
/* copy float2 array from global GPU memory to main memory */
   float2 *g_f;
   g_f = (float2 *)*gp_f;
   gpu_ccopyout(f,g_f,*nsize);
   return;
}

/*--------------------------------------------------------------------*/
extern "C" void gpu_initstream_(int *nstream) {
   gpu_initstream(*nstream);
   return;
}

/*--------------------------------------------------------------------*/
extern "C" void gpu_delstream_(int *nstream) {
   gpu_delstream(*nstream);
   return;
}

/*--------------------------------------------------------------------*/
extern "C" void gpu_waitstream_(int *nstream) {
   gpu_waitstream(*nstream);
   return;
}

/*--------------------------------------------------------------------*/
extern "C" void gpu_cascopyin_(float2 *f, unsigned long *gp_f,
                               int *noff, int *nsize, int *nstream) {
/* copy float2 array segment from main memory to global GPU memory */
/* asynchronous copy */
   float2 *g_f;
   g_f = (float2 *)*gp_f;
   gpu_cascopyin(f,g_f,*noff,*nsize,*nstream);
   return;
}

/*--------------------------------------------------------------------*/
extern "C" void gpu_cascopyout_(float2 *f, unsigned long *gp_f,
                                int *noff, int *nsize, int *nstream) {
/* copy float2 array segment from global GPU memory to main memory */
/* asynchronous copy */
   float2 *g_f;
   g_f = (float2 *)*gp_f;
   gpu_cascopyout(f,g_f,*noff,*nsize,*nstream);
   return;
}

/*--------------------------------------------------------------------*/
extern "C" void gpu_zfmem_(unsigned long *gp_f, int *nsize) {
   float *g_f;
   g_f = (float *)*gp_f;
   gpu_zfmem(g_f,*nsize);
   return;
}

/*--------------------------------------------------------------------*/
extern "C" void gpu_zcmem_(unsigned long *gp_f, int *nsize) {
   float2 *g_f;
   g_f = (float2 *)*gp_f;
   gpu_zcmem(g_f,*nsize);
   return;
}

/*--------------------------------------------------------------------*/
extern "C" void gpu_set_cache_size_(int *nscache) {
   gpu_set_cache_size(*nscache);
   return;
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
