/* CUDA FFT Library */
/* written by Viktor K. Decyk, UCLA */

#include <stdlib.h>
#include <stdio.h>
#include "cuda.h"
#include <cufft.h>

extern int nblock_size;
extern int maxgsx;

static cudaError_t crc;
static cufftResult cfrc;
static cufftHandle planrx, planxr, planrxn, planxrn;
static cufftHandle plany, planyn;

__global__ void gpuctpose4(float2 f[], float2 g[], int nx, int ny,
                           int nxv, int nyv);

__global__ void gpuctpose4n(float2 fn[], float2 gn[], int nx, int ny,
                            int ndim, int nxv, int nyv);

/*--------------------------------------------------------------------*/
__global__ void gpusctpose4(float2 f[], float2 g[], float ani, int nx,
                            int ny, int nxv, int nyv) {
/* scaled complex transpose using blocking algorithm with gaps */
/* local data */
   int j, k, js, ks, joff, koff, mx, mxv;
   float2 a;
/* The size of the shared memory array is as follows: */
/* float2 shm2[(mx + 1)*mx];                          */
   extern __shared__ float2 shm2[];
   mx = blockDim.x;
   mxv = mx + 1;
   joff = mx*blockIdx.x;
   koff = mx*blockIdx.y;
   js = threadIdx.x;
   ks = threadIdx.y;
/* copy into block */
   j = js + joff;
   k = ks + koff;
   if ((j < nx) && (k < ny)) {
      shm2[js+mxv*ks] = f[j+nxv*k];
   }
   __syncthreads();
/* copy out from block with scaling */
   j = ks + joff;
   k = js + koff;
   if ((j < nx) && (k < ny)) {
      a = shm2[ks+mxv*js];
      a.x = ani*a.x;
      a.y = ani*a.y;
      g[k+nyv*j] = a;
   }
   return;
}

/*--------------------------------------------------------------------*/
__global__ void gpusctpose4n(float2 fn[], float2 gn[], float ani,
                             int nx, int ny, int ndim, int nxv,
                             int nyv) {
/* scaled complex vector transpose using blocking algorithm with gaps */
/* ndim = vector dimension                                            */
/* local data */
   int i, j, k, js, ks, joff, koff, mx, mxv, nmxv, nnxv, nnyv, jj, kk;
   float2 a;
/* The size of the shared memory array is as follows: */
/* float2 shmn2[ndim*(mx + 1)*mx];                    */
   extern __shared__ float2 shmn2[];
   mx = blockDim.x;
   mxv = mx + 1;
   joff = mx*blockIdx.x;
   koff = mx*blockIdx.y;
   js = threadIdx.x;
   ks = threadIdx.y;
   nmxv = ndim*mxv;
   nnxv = ndim*nxv;
   nnyv = ndim*nyv;
/* copy into block */
   j = js + joff;
   k = ks + koff;
   if ((j < nx) && (k < ny)) {
      jj = j + nnxv*k;
      kk = js + nmxv*ks;
      for (i = 0; i < ndim; i++) {
         shmn2[kk+mxv*i] = fn[jj+nxv*i];
      }
   }
   __syncthreads();
/* copy out from block with scaling */
   j = ks + joff;
   k = js + koff;
   if ((j < nx) && (k < ny)) {
      kk = k + nnyv*j;
      jj = ks + nmxv*js;
      for (i = 0; i < ndim; i++) {
         a = shmn2[jj+mxv*i];
         a.x = ani*a.x;
         a.y = ani*a.y;
         gn[kk+nyv*i] = a;
      }
   }
   return;
}

/*--------------------------------------------------------------------*/	
extern "C" void gpufft2rrcuinit(int nx, int ny, int ndim) {
   cfrc = cufftPlan1d(&planrx,nx,CUFFT_R2C,ny);
   if (cfrc) {
      printf("cufftPlan1d planrx error=%d\n",cfrc);
      exit(1);
   }
   cfrc = cufftPlan1d(&planxr,nx,CUFFT_C2R,ny);
   if (cfrc) {
      printf("cufftPlan1d planxr error=%d\n",cfrc);
      exit(1);
   }
   cfrc = cufftPlan1d(&planrxn,nx,CUFFT_R2C,ndim*ny);
   if (cfrc) {
      printf("cufftPlan1d planrxn error=%d\n",cfrc);
      exit(1);
   }
   cfrc = cufftPlan1d(&planxrn,nx,CUFFT_C2R,ndim*ny);
   if (cfrc) {
      printf("cufftPlan1d planxrn error=%d\n",cfrc);
      exit(1);
   }
   return;
}

/*--------------------------------------------------------------------*/	
extern "C" void gpufft2cuinit(int nx, int ny, int ndim) {
   int nxh1;
   nxh1 = nx/2 + 1;
   cfrc = cufftPlan1d(&plany,ny,CUFFT_C2C,nxh1);
   if (cfrc) {
      printf("cufftPlan1d plany error=%d\n",cfrc);
      exit(1);
   }
   cfrc = cufftPlan1d(&planyn,ny,CUFFT_C2C,ndim*nxh1);
   if (cfrc) {
      printf("cufftPlan1d planyn error=%d\n",cfrc);
      exit(1);
   }
   return;
}

/*--------------------------------------------------------------------*/
extern "C" void gpufft2rrcudel() {
   cfrc = cufftDestroy(planrx);
   if (cfrc) {
      printf("cufftDestroy planrx error=%d\n",cfrc);
      exit(1);
   }
   cfrc = cufftDestroy(planxr);
   if (cfrc) {
      printf("cufftDestroy planxr error=%d\n",cfrc);
      exit(1);
   }
   cfrc = cufftDestroy(planrxn);
   if (cfrc) {
      printf("cufftDestroy planrxn error=%d\n",cfrc);
      exit(1);
   }
   cfrc = cufftDestroy(planxrn);
   if (cfrc) {
      printf("cufftDestroy planxrn error=%d\n",cfrc);
      exit(1);
   }
   return;
}

/*--------------------------------------------------------------------*/
extern "C" void gpufft2cudel() {
   cfrc = cufftDestroy(plany);
   if (cfrc) {
      printf("cufftDestroy plany error=%d\n",cfrc);
      exit(1);
   }
   cfrc = cufftDestroy(planyn);
   if (cfrc) {
      printf("cufftDestroy planyn error=%d\n",cfrc);
      exit(1);
   }
   return;
}

/*--------------------------------------------------------------------*/
extern "C" void gpufft2rrcu(float2 f[], float2 g[], int isign,
                            int indx, int indy, int nxh1d, int nyd) {
/* wrapper function for real to complex fft, without packed data */
/* uses 1D real to complex and complex to complex NVIDIA FFTs    */
/* nxh1d must be = nx/2+1                                        */
/* local data */
   int nx, nxh1, ny, ns;
   int mx = 16;
   float ani;
   dim3 dimBlock(nblock_size);
   dim3 dimBlockt(mx,mx);
/* calculate range of indices */
   nx = 1L<<indx;
   nxh1 = nx/2 + 1;
   ny = 1L<<indy;
   dim3 dimGridtx((nxh1-1)/mx+1,(ny-1)/mx+1);
   dim3 dimGridty((ny-1)/mx+1,(nxh1-1)/mx+1);
   ns = (mx+1)*mx*sizeof(float2);
/* inverse fourier transform */
   if (isign < 0) {
/* perform x fft */
      cfrc = cufftExecR2C(planrx,(cufftReal *)f,(cufftComplex *)f);
/*    cudaThreadSynchronize(); */
      if (cfrc) {
         printf("cufftExecR2C(-1) planrx error=%d\n",cfrc);
         exit(1);
      }
/* transpose f to g and normalize */
      ani = 1.0f/(((float) nx)*((float) ny));
      crc = cudaGetLastError();
      gpusctpose4<<<dimGridtx,dimBlockt,ns>>>(f,g,ani,nxh1,ny,nxh1d,
                                              nyd);
/*    cudaThreadSynchronize(); */
      crc = cudaGetLastError();
      if (crc) {
         printf("gpusctpose4 error=%d:%s\n",crc,
                cudaGetErrorString(crc));
         exit(1);
      }
/* perform y fft */
      cfrc = cufftExecC2C(plany,(cufftComplex *)g,(cufftComplex *)g,
                          CUFFT_FORWARD);
      cudaThreadSynchronize();
      if (cfrc) {
         printf("cufftExecC2C(-1) plany error=%d\n",cfrc);
         exit(1);
      }
   }
/* forward fourier transform */
   else if (isign > 0) {
/* perform y fft */
      cfrc = cufftExecC2C(plany,(cufftComplex *)g,(cufftComplex *)g,
                          CUFFT_INVERSE);
/*    cudaThreadSynchronize(); */
      if (cfrc) {
         printf("cufftExecC2C(1) plany error=%d\n",cfrc);
         exit(1);
      }
/* transpose g to f */
      crc = cudaGetLastError();
      gpuctpose4<<<dimGridty,dimBlockt,ns>>>(g,f,ny,nxh1,nyd,nxh1d);
/*    cudaThreadSynchronize(); */
      crc = cudaGetLastError();
      if (crc) {
         printf("gpuctpose4 error=%d:%s\n",crc,cudaGetErrorString(crc));
         exit(1);
      }
/* perform x fft */
      cfrc = cufftExecC2R(planxr,(cufftComplex *)f,(cufftReal *)f);
      cudaThreadSynchronize();
      if (cfrc) {
         printf("cufftExecC2R(1) planxr error=%d\n",cfrc);
         exit(1);
      }
   }
   return;
}

/*--------------------------------------------------------------------*/
extern "C" void gpufft2rrcun(float2 fn[], float2 gn[], int isign,
                             int indx, int indy, int ndim, int nxh1d,
                             int nyd) {
/* wrapper function for real to complex fft, without packed data */
/* for vector data                                               */
/* uses 1D real to complex and complex to complex NVIDIA FFTs    */
/* ndim = vector dimension                                       */
/* nxh1d must be = nx/2+1                                        */
/* local data */
   int nx, nxh1, ny, ns;
   int mx = 16;
   float ani;
   dim3 dimBlock(nblock_size);
   dim3 dimBlockt(mx,mx);
/* calculate range of indices */
   nx = 1L<<indx;
   nxh1 = nx/2 + 1;
   ny = 1L<<indy;
   dim3 dimGridtx((nxh1-1)/mx+1,(ny-1)/mx+1);
   dim3 dimGridty((ny-1)/mx+1,(nxh1-1)/mx+1);
   ns = ndim*(mx+1)*mx*sizeof(float2);
/* inverse fourier transform */
   if (isign < 0) {
/* perform x fft */
      cfrc = cufftExecR2C(planrxn,(cufftReal *)fn,(cufftComplex *)fn);
/*    cudaThreadSynchronize(); */
      if (cfrc) {
         printf("cufftExecR2C(-1) planrxn error=%d\n",cfrc);
         exit(1);
      }
/* transpose f to g and normalize */
      ani = 1.0f/(((float) nx)*((float) ny));
      crc = cudaGetLastError();
      gpusctpose4n<<<dimGridtx,dimBlockt,ns>>>(fn,gn,ani,nxh1,ny,ndim,
                                               nxh1d,nyd);
/*    cudaThreadSynchronize(); */
      crc = cudaGetLastError();
      if (crc) {
         printf("gpusctpose4n error=%d:%s\n",crc,
                cudaGetErrorString(crc));
         exit(1);
      }
/* perform y fft */
      cfrc = cufftExecC2C(planyn,(cufftComplex *)gn,(cufftComplex *)gn,
                          CUFFT_FORWARD);
      cudaThreadSynchronize();
      if (cfrc) {
         printf("cufftExecC2C(-1) planyn error=%d\n",cfrc);
         exit(1);
      }
   }
/* forward fourier transform */
   else if (isign > 0) {
/* perform y fft */
      cfrc = cufftExecC2C(planyn,(cufftComplex *)gn,(cufftComplex *)gn,
                          CUFFT_INVERSE);
/*    cudaThreadSynchronize(); */
      if (cfrc) {
         printf("cufftExecC2C(1) planyn error=%d\n",cfrc);
         exit(1);
      }
/* transpose g to f */
      crc = cudaGetLastError();
      gpuctpose4n<<<dimGridty,dimBlockt,ns>>>(gn,fn,ny,nxh1,ndim,nyd,
                                              nxh1d);
/*    cudaThreadSynchronize(); */
      crc = cudaGetLastError();
      if (crc) {
         printf("gpuctpose4n error=%d:%s\n",crc,
                cudaGetErrorString(crc));
         exit(1);
      }
/* perform x fft */
      cfrc = cufftExecC2R(planxrn,(cufftComplex *)fn,(cufftReal *)fn);
      cudaThreadSynchronize();
      if (cfrc) {
         printf("cufftExecC2R(1) planxrn error=%d\n",cfrc);
         exit(1);
      }
   }
   return;
}

/* Interfaces to Fortran */

/*--------------------------------------------------------------------*/	
extern "C" void gpufft2rrcuinit_(int *nx, int *ny, int *ndim) {
   gpufft2rrcuinit(*nx,*ny,*ndim);
   return;
}

/*--------------------------------------------------------------------*/
extern "C" void gpufft2cuinit_(int *nx, int *ny, int *ndim) {
   gpufft2cuinit(*nx,*ny,*ndim);
   return;
}

/*--------------------------------------------------------------------*/
extern "C" void gpufft2rrcudel_() {
   gpufft2rrcudel();
   return;
}

/*--------------------------------------------------------------------*/
extern "C" void gpufft2cudel_() {
   gpufft2cudel();
   return;
}

/*--------------------------------------------------------------------*/
extern "C" void gpufft2rrcu_(unsigned long *gp_f, unsigned long *gp_g,
                             int *isign, int *indx, int *indy,
                             int *nxh1d, int *nyd) {
   float2 *f, *g;
   f = (float2 *)*gp_f;
   g = (float2 *)*gp_g;
   gpufft2rrcu(f,g,*isign,*indx,*indy,*nxh1d,*nyd);
   return;
}

/*--------------------------------------------------------------------*/
extern "C" void gpufft2rrcun_(unsigned long *gp_fn,
                              unsigned long *gp_gn, int *isign,
                              int *indx, int *indy, int *ndim,
                              int *nxh1d, int *nyd) {
   float2 *fn, *gn;
   fn = (float2 *)*gp_fn;
   gn = (float2 *)*gp_gn;
   gpufft2rrcun(fn,gn,*isign,*indx,*indy,*ndim,*nxh1d,*nyd);
   return;
}

