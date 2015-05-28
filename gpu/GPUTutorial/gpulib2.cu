/*--------------------------------------------------------------------*/
/* CUDA Library for GPU Tutorial */
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

__global__ void emptyKernel() {}

/*--------------------------------------------------------------------*/
extern "C" void setgbsize(int nblock) {
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
         printf("    CUDA_GLOBAL_MEM_SIZE=%u(%f GB),Capability=%d\n",
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
/* terminate CUDA */
   crc = cudaThreadExit();
   if (crc) {
      printf("cudaThreadExit Error=%d:%s\n",crc,cudaGetErrorString(crc));
   }
   return;
}

/*--------------------------------------------------------------------*/
__global__ void gcopy1(float a[], float b[], int nx) {
/* 1d copy a = b                            */
/* one block of mx threads copies mx values */
/* ((nx-1)/mx+1) independent blocks         */
/* nx = size of arrays in x                 */
/* local data */
   int j, js, id, mx;
   mx = blockDim.x;
   j = threadIdx.x;
   id = blockIdx.x;

   js = j + mx*id;
   if (js < nx) a[js] = b[js];

   return;
}

/*--------------------------------------------------------------------*/
__global__ void gcopy2a(float a[], float b[], int nx, int ny) {
/* 2d copy a = b                            */
/* one block of mx threads copies mx values */
/* nbx*ny independent blocks                */
/* local data */
   int j, k, js, id, mx;
   mx = blockDim.x;
   j = threadIdx.x;
   id = blockIdx.x;
   k = blockIdx.y;

   js = j + mx*id;
   if ((js < nx) && (k < ny)) {
      a[js+nx*k] = b[js+nx*k];
   }

   return;
}

/*--------------------------------------------------------------------*/
__global__ void gcopy2b(float a[], float b[], int nx, int ny) {
/* 2d copy a = b                            */
/* one block of mx threads copies nx values */
/* ny independent blocks                    */
/* local data */
   int j, k, mx;
   mx = blockDim.x;
   k = blockIdx.x;

   j = threadIdx.x;
   while (j < nx) {
      if (k < ny)
         a[j+nx*k] = b[j+nx*k];
      j += mx;
   }

   return;
}

/*--------------------------------------------------------------------*/
__global__ void gsaxpy2(float a[], float b[], float s, int nx, int ny) {
/* 2d vector multiplye a = s*b + a          */
/* one block of mx threads copies nx values */
/* ny independent blocks                    */
/* local data */
   int j, k, mx;
   mx = blockDim.x;
   k = blockIdx.x;

   j = threadIdx.x;
   while (j < nx) {
      if (k < ny)
         a[j+nx*k] = s*b[j+nx*k] + a[j+nx*k];
      j += mx;
   }

   return;
}

/*--------------------------------------------------------------------*/
__global__ void gcopy3(float a[], float b[], int nx, int ny) {
/* 2d copy a = b                                  */
/* one block of mx*my threads copies mx*my values */
/* ((nx-1)/mx+1)*((ny-1)/my+1) independent blocks */
/* local data */
   int j, k, js, ks, idx, idy, mx, my;
   mx = blockDim.x; my = blockDim.y;
   j = threadIdx.x; k = threadIdx.y;
   idx = blockIdx.x; idy = blockIdx.y;

   ks = k + my*idy;
   js = j + mx*idx;
   if ((js < nx) && (ks < ny))
       a[js+nx*ks] = b[js+nx*ks];

   return;
}

/*--------------------------------------------------------------------*/
__global__ void  gtranspose2(float a[], float b[], int nx, int ny) {
/* a = transpose(b)                                   */
/* one block of mx*mx threads transposes mx*mx values */
/* ((nx-1)/mx+1)*((ny-1)/mx+1) independent blocks     */
/* local data */
   int j, k, js, ks, idx, idy, joff, koff, mx, mxv;
   extern __shared__ float s[];
   mx = blockDim.x; mxv = mx + 1;
   j = threadIdx.x; k = threadIdx.y;
   idx = blockIdx.x; idy = blockIdx.y;
   koff = mx*idy;
   joff = mx*idx;

   ks = k + koff;
   js = j + joff;
   if ((js < nx) && (ks < ny))
      s[j+mxv*k] = b[js+nx*ks];
/* synchronize threads */
   __syncthreads();
   js = k + joff;
   ks = j + koff;
   if ((js < nx) && (ks < ny))
      a[ks+ny*js] = s[k+mxv*j];

   return;
}

/*--------------------------------------------------------------------*/
__global__ void gsum1(float a[], float *sa, int nx) {
/* 1d serial sum reductions, each of length mx */
/* sa = sum(a)                                 */
/* local data */
   int j, js, jb, mx, joff, mxm;
   float t;
   extern __shared__ float s[];
   mx = blockDim.x;
   js = threadIdx.x;
   jb = blockIdx.x;
   joff = mx*jb;

   j = js + joff;
/* copy global data to shared memory */
   if (j < nx) s[js] = a[j];
/* synchronize to make sure each thread in block has the data */
   __syncthreads();
   if (js==0) {
      mxm = nx - joff;
      if (mxm > mx) mxm = mx;
/* perform serial local sum reduction: result in t */
      t = 0.0f;
      for (j = 0; j < mxm; j++) {
         t += s[j];
      }
/* accumulate results to global memory for each block */
/* for devices with compute capability 2.x            */
      atomicAdd(&sa[0],t);
   }

   return;
}

/*--------------------------------------------------------------------*/
__device__ void lsum2(float *sdata, int n) {
/* finds local sum of nths data items shared by threads */
/* using binary tree method. input is modified. */
/* local data */
   int l, k;
   float s;
   l = threadIdx.x;
   k = blockDim.x >> 1;
   s = 0.0f;

   if (l < n) s = sdata[l];
   while (k > 0) {
      if (l < k) {
         if ((l+k) < n) {
            s += sdata[l+k];
            sdata[l] = s;
         }
      }
      __syncthreads();
      k >>= 1;
   }

   return;
}

/*--------------------------------------------------------------------*/
__global__ void gsum2(float a[], float d[], int nx) {
/* segmented 1d sum reductions, each of length mx             */
/* forall (j = 1:nbx); d(j) = sum(a(1+mx*(j-1):min(nx,mx*j))) */
/* parallel summation                                         */
/* local data */
   int j, js, jb, mx, joff, mxm;
   extern __shared__ float s[];
   mx = blockDim.x;
   js = threadIdx.x;
   jb = blockIdx.x;
   joff = mx*jb;

   j = js + joff;
/* copy global data to shared memory */
   if (j < nx) s[js] = a[j];
/* synchronize to make sure each thread in block has the data */
   __syncthreads();
   mxm = nx - joff;
   if (mxm > mx) mxm = mx;
/* perform parallel local sum reduction: result in s[0] */
   lsum2(s,mxm);
/* write out result to global memory for each block */
   if (js==0) d[jb] = s[0];

   return;
}

/*--------------------------------------------------------------------*/
extern "C" void gpu_copy1(float *a, float *b, int mx, int nx) {
/* 2d copy of length nx, with block size mx */
/* one block of mx threads copies mx values */
/* ((nx-1)/mx+1) independent blocks         */
/* local data */
   int nbx;
   nbx = (nx - 1)/mx + 1;
   dim3 dimBlock(mx);
   dim3 dimGrid(nbx);
   crc = cudaGetLastError();
 
   gcopy1<<<dimGrid,dimBlock>>>(a,b,nx);
 
   cudaThreadSynchronize();
   crc = cudaGetLastError();
   if (crc) {
      printf("gcopy1 error=%d:%s\n",crc,cudaGetErrorString(crc));
      exit(1);
   }
   return;
}

/*--------------------------------------------------------------------*/
extern "C" void gpu_copy2a(float *a, float *b, int mx, int nx, int ny) {
/* 2d copy a = b                            */
/* one block of mx threads copies mx values */
/* nbx*ny independent blocks                */
/* local data */
   int nbx;
   nbx = (nx - 1)/mx + 1;
   dim3 dimBlock(mx);
   dim3 dimGrid(nbx,ny);
   crc = cudaGetLastError();

   gcopy2a<<<dimGrid,dimBlock>>>(a,b,nx,ny);
 
   cudaThreadSynchronize();
   crc = cudaGetLastError();
   if (crc) {
      printf("gcopy2a error=%d:%s\n",crc,cudaGetErrorString(crc));
      exit(1);
   }
   return;
}

/*--------------------------------------------------------------------*/
extern "C" void gpu_copy2b(float *a, float *b, int mx, int nx, int ny) {
/* 2d copy a = b                            */
/* one block of mx threads copies nx values */
/* ny independent blocks                    */
/* local data */
   dim3 dimBlock(mx);
   dim3 dimGrid(ny);
   crc = cudaGetLastError();

   gcopy2b<<<dimGrid,dimBlock>>>(a,b,nx,ny);
 
   cudaThreadSynchronize();
   crc = cudaGetLastError();
   if (crc) {
      printf("gcopy2b error=%d:%s\n",crc,cudaGetErrorString(crc));
      exit(1);
   }
   return;
}

/*--------------------------------------------------------------------*/
extern "C" void gpu_saxpy2(float *a, float *b, float s, int mx, int nx,
                           int ny) {
/* 2d vector multiply a = s*b + a           */
/* one block of mx threads copies nx values */
/* ny independent blocks                    */
/* local data */
   dim3 dimBlock(mx);
   dim3 dimGrid(ny);
   crc = cudaGetLastError();

   gsaxpy2<<<dimGrid,dimBlock>>>(a,b,s,nx,ny);
 
   cudaThreadSynchronize();
   crc = cudaGetLastError();
   if (crc) {
      printf("gsaxpy2 error=%d:%s\n",crc,cudaGetErrorString(crc));
      exit(1);
   }
   return;
}

/*--------------------------------------------------------------------*/
extern "C" void gpu_copy3(float *a, float *b, int mx, int my, int nx,
                          int ny) {
/* 2d copy a = b                                  */
/* one block of mx*my threads copies mx*my values */
/* ((nx-1)/mx+1)*((ny-1)/my+1) independent blocks */
/* local data */
   int nbx, nby;
   nbx = (nx - 1)/mx + 1; nby = (ny - 1)/my + 1;
   dim3 dimBlock(mx,my);
   dim3 dimGrid(nbx,nby);
   crc = cudaGetLastError();

   gcopy3<<<dimGrid,dimBlock>>>(a,b,nx,ny);
 
   cudaThreadSynchronize();
   crc = cudaGetLastError();
   if (crc) {
      printf("gcopy3 error=%d:%s\n",crc,cudaGetErrorString(crc));
      exit(1);
   }
   return;
}

/*--------------------------------------------------------------------*/
extern "C" void  gpu_transpose2(float *a, float *b, int mx, int nx,
                                int ny) {
/* 2d transpose of length nx, ny, with block size mx, mx */
/* one block of mx*mx threads transposes mx*mx values    */
/* ((nx-1)/mx+1)*((ny-1)/mx+1) independent blocks        */
/* local data */
   int nbx, nby, ns;
   nbx = (nx - 1)/mx + 1; nby = (ny - 1)/mx + 1;
   dim3 dimBlock(mx,mx);
   dim3 dimGrid(nbx,nby);
/* calculate size of shared memory */
   ns = (mx + 1)*mx*sizeof(float);
   crc = cudaGetLastError();

   gtranspose2<<<dimGrid,dimBlock,ns>>>(a,b,nx,ny);
 
   cudaThreadSynchronize();
   crc = cudaGetLastError();
   if (crc) {
      printf("gtranspose2 error=%d:%s\n",crc,cudaGetErrorString(crc));
      exit(1);
   }
   return;
}

/*--------------------------------------------------------------------*/
extern "C" void gpu_sum1(float *a, float *sa, int mx, int nx) {
/* 1d serial sum reductions, each of length mx */
/* one block of mx threads sums mx values      */
/* ((nx-1)/mx+1) independent blocks            */
/* local data */
   int nbx, ns;
   float t;
   nbx = (nx - 1)/mx + 1;
   dim3 dimBlock(mx);
   dim3 dimGrid(nbx);
   t = 0.0f;
   gpu_fcopyin(&t,sa,1);
/* calculate size of shared memory */
   ns = mx*sizeof(float);
   crc = cudaGetLastError();
 
   gsum1<<<dimGrid,dimBlock,ns>>>(a,sa,nx);
 
   cudaThreadSynchronize();
   crc = cudaGetLastError();
   if (crc) {
      printf("gsum1 error=%d:%s\n",crc,cudaGetErrorString(crc));
      exit(1);
   }
   return;
}

/*--------------------------------------------------------------------*/
extern "C" void gpu_sum2(float *a, float *d, int mx, int nx) {
/* segmented 1d parallel sum reductions, each of length mx */
/* one block of mx threads sums mx values                  */
/* ((nx-1)/mx+1) independent blocks                        */
/* local data */
   int nbx, ns;
   nbx = (nx - 1)/mx + 1;
   dim3 dimBlock(mx);
   dim3 dimGrid(nbx);
/* calculate size of shared memory */
   ns = mx*sizeof(float);
   crc = cudaGetLastError();
 
   gsum2<<<dimGrid,dimBlock,ns>>>(a,d,nx);
 
   cudaThreadSynchronize();
   crc = cudaGetLastError();
   if (crc) {
      printf("gsum2 error=%d:%s\n",crc,cudaGetErrorString(crc));
      exit(1);
   }
   return;
}

/*--------------------------------------------------------------------*/
extern "C" void gpu_sum3(float *a, float *d, float *sa, int mx,
                         int nx) {
/* segmented 1d parallel sum reductions   */
/* one block of mx threads sums mx values */
/* ((nx-1)/mx+1) independent blocks       */
/* local data */
   int nxs, nbx, n, ns;
   nxs = nx;
   nbx = (nxs - 1)/mx + 1;
   dim3 dimBlock(mx);
   dim3 dimGrid(nbx);
/* calculate size of shared memory */
   ns = mx*sizeof(float);
   crc = cudaGetLastError();
   gsum2<<<dimGrid,dimBlock,ns>>>(a,d,nxs);
/* cudaThreadSynchronize(); */
   crc = cudaGetLastError();
   if (crc) {
      printf("gsum2:0 error=%d:%s\n",crc,cudaGetErrorString(crc));
      exit(1);
   }
/* write out result */
   if (nbx==1) {
      dimGrid.x = 1;
      crc = cudaGetLastError();
      gcopy1<<<dimGrid,dimBlock>>>(sa,d,1);
      cudaThreadSynchronize();
      crc = cudaGetLastError();
      if (crc) {
         printf("gcopy1:0 error=%d:%s\n",crc,cudaGetErrorString(crc));
         exit(1);
      }
      return;
   }
/* reduce further if necessary */
   if (nbx > 1) {
      nxs = nbx;
      nbx = (nxs - 1)/mx + 1;
      dimGrid.x = nbx;
      crc = cudaGetLastError();
      gsum2<<<dimGrid,dimBlock,ns>>>(d,sa,nxs);
      if (nbx==1)
         cudaThreadSynchronize();
      crc = cudaGetLastError();
      if (crc) {
         printf("gsum2:1 error=%d:%s\n",crc,cudaGetErrorString(crc));
         exit(1);
      }
      n = 0;
   }
   if (nbx==1)
      return;
/* iterate if necessary */
   while (nbx > 1) {
      n += nbx;
      nxs = nbx;
      nbx = (nxs - 1)/mx + 1;
      dimGrid.x = nbx;
      crc = cudaGetLastError();
      gsum2<<<dimGrid,dimBlock,ns>>>(&sa[n-nxs],&sa[n],nxs);
/*    cudaThreadSynchronize(); */
      crc = cudaGetLastError();
      if (crc) {
         printf("gsum2:n error=%d:%s\n",crc,cudaGetErrorString(crc));
         exit(1);
      }
   }
/* write out result */
   dimGrid.x = 1;
   crc = cudaGetLastError();
   gcopy1<<<dimGrid,dimBlock>>>(sa,&sa[n],1);
   cudaThreadSynchronize();
   crc = cudaGetLastError();
   if (crc) {
      printf("gcopy1:n error=%d:%s\n",crc,cudaGetErrorString(crc));
      exit(1);
   }
   return;
}


/* Interfaces to Fortran */

/*--------------------------------------------------------------------*/
extern "C" void setgbsize_(int *nblock) {
   setgbsize(*nblock);
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
extern "C" void gpu_deallocate_(unsigned long *gp_d, int *irc) {
/* deallocate global memory on GPU, return pointer to Fortran */
   void *d;
   d = (void *)*gp_d;
   gpu_deallocate(d,irc);
   *gp_d = 0;
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
extern "C" void emptykernel_() {
   emptykernel();
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

/*--------------------------------------------------------------------*/
extern "C" void gpu_copy1_(unsigned long *gp_a, unsigned long *gp_b,
                           int *mx, int *nx) {
   float *a, *b;
   a = (float *)*gp_a;
   b = (float *)*gp_b;
   gpu_copy1(a,b,*mx,*nx);
   return;
}

/*--------------------------------------------------------------------*/
extern "C" void gpu_copy2a_(unsigned long *gp_a, unsigned long *gp_b,
                            int *mx, int *nx, int *ny) {
   float *a, *b;
   a = (float *)*gp_a;
   b = (float *)*gp_b;
   gpu_copy2a(a,b,*mx,*nx,*ny);
   return;
}

/*--------------------------------------------------------------------*/
extern "C" void gpu_copy2b_(unsigned long *gp_a, unsigned long *gp_b,
                            int *mx, int *nx, int *ny) {
   float *a, *b;
   a = (float *)*gp_a;
   b = (float *)*gp_b;
   gpu_copy2b(a,b,*mx,*nx,*ny);
   return;
}

/*--------------------------------------------------------------------*/
extern "C" void gpu_saxpy2_(unsigned long *gp_a, unsigned long *gp_b,
                            float *s, int *mx, int *nx, int *ny) {
   float *a, *b;
   a = (float *)*gp_a;
   b = (float *)*gp_b;
   gpu_saxpy2(a,b,*s,*mx,*nx,*ny);
   return;
}

/*--------------------------------------------------------------------*/
extern "C" void gpu_copy3_(unsigned long *gp_a, unsigned long *gp_b,
                           int *mx, int *my, int *nx, int *ny) {
   float *a, *b;
   a = (float *)*gp_a;
   b = (float *)*gp_b;
   gpu_copy3(a,b,*mx,*my,*nx,*ny);
   return;
}

/*--------------------------------------------------------------------*/
extern "C" void  gpu_transpose2_(unsigned long *gp_a,
                                 unsigned long *gp_b, int *mx, int *nx,
                                 int *ny) {
   float *a, *b;
   a = (float *)*gp_a;
   b = (float *)*gp_b;
   gpu_transpose2(a,b,*mx,*nx,*ny);
   return;
}

/*--------------------------------------------------------------------*/
extern "C" void gpu_sum1_(unsigned long *gp_a, unsigned long *gp_sa,
                          int *mx, int *nx) {
   float *a, *sa;
   a = (float *)*gp_a;
   sa = (float *)*gp_sa;
   gpu_sum1(a,sa,*mx,*nx);
   return;
}

/*--------------------------------------------------------------------*/
extern "C" void gpu_sum2_(unsigned long *gp_a, unsigned long *gp_d,
                          int *mx, int *nx) {
   float *a, *d;
   a = (float *)*gp_a;
   d = (float *)*gp_d;
   gpu_sum2(a,d,*mx,*nx);
   return;
}

/*--------------------------------------------------------------------*/
extern "C" void gpu_sum3_(unsigned long *gp_a, unsigned long *gp_d,
                          unsigned long *gp_sa,int *mx, int *nx) {
   float *a, *d, *sa;
   a = (float *)*gp_a;
   d = (float *)*gp_d;
   sa = (float *)*gp_sa;
   gpu_sum3(a,d,sa,*mx,*nx);
   return;
}
