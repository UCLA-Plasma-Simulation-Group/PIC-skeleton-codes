/* GPU vector add test program for OpenCL */
/* written by Viktor K. Decyk, UCLA */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
/* For Linux */
#include <CL/cl.h>
/* For Mac OS */
/* #include <OpenCL/opencl.h> */

#define NKERNELS                 2
#define NAMEL                    33

static int nblock_size = 64;

static cl_int crc, status;
static cl_uint kernels;
static char param[80];
static size_t psize;
static cl_device_id devid;
static cl_context cnt;
static cl_command_queue cmq;
static cl_program clp;
static cl_kernel clks[NKERNELS];
static char knames[NKERNELS][NAMEL];
static cl_event evnt;

/*--------------------------------------------------------------------*/
void setgbsize(int nblock) {
   nblock_size = nblock;
   return;
}

/*--------------------------------------------------------------------*/
void g_fallocate(float **g_f, int nsize, int *irc) {
/* allocate global float memory on GPU, return pointer to C */
/* requires a valid context cnt                             */
/* error code is modified only if there is an error         */
   cl_mem clm;
/* create a buffer object */
   clm = clCreateBuffer(cnt,CL_MEM_READ_WRITE,sizeof(float)*nsize,NULL,
                        &crc);
   if (crc != CL_SUCCESS) {
      printf("clCreateBuffer Error=%i\n",crc);
      *irc = 1;
   }
   *g_f = (float *)clm;
   return;
}

/*--------------------------------------------------------------------*/
void g_f4allocate(float **g_f, int nsize, int *irc) {
/* allocate global vector float memory on GPU, return pointer to C */
/* requires a valid context cnt                             */
/* error code is modified only if there is an error         */
   int ns;
   cl_mem clm;
/* create a buffer object */
   ns = (nsize - 1)/4 + 1;
   clm = clCreateBuffer(cnt,CL_MEM_READ_WRITE,sizeof(cl_float4)*ns,NULL,
                        &crc);
   if (crc != CL_SUCCESS) {
      printf("clCreateBuffer Vector Error=%i\n",crc);
      *irc = 1;
   }
   *g_f = (float *)clm;
   return;
}

/*--------------------------------------------------------------------*/
void g_deallocate(float **g_f, int *irc) {
/* deallocate global memory on GPU */
/* error code is modified only if there is an error       */
/* decrement memory object reference count */
   cl_mem clm;
   clm = (cl_mem) *g_f;
   crc = clReleaseMemObject(clm);
   if (crc != CL_SUCCESS) {
      printf("clReleaseMemObject Error=%i\n",crc);
      *irc = 1;
   }
   *g_f = NULL;
   return;
}

/*--------------------------------------------------------------------*/
void copyin_gmemptr(float *f, float *g_f, int nsize) {
/* copy float array from main memory to global GPU memory */
/* requires a valid command queue cmq                     */
/* error code is modified only if there is an error       */
   cl_mem clm;
   clm = (cl_mem) g_f;
   crc = clEnqueueWriteBuffer(cmq,clm,CL_TRUE,0,sizeof(float)*nsize,f,
                              0,NULL,NULL);
   if (crc != CL_SUCCESS) {
      printf("clEnqueueWriteBuffer Error=%i\n",crc);
      exit(1);
   }
   return;
}

/*--------------------------------------------------------------------*/
void copyout_gmemptr(float *f, float *g_f, int nsize) {
/* copy float array from global GPU memory to main memory */
/* requires a valid command queue cmq                     */
/* error code is modified only if there is an error       */
   cl_mem clm;
   clm = (cl_mem) g_f;
   crc = clEnqueueReadBuffer(cmq,clm,CL_TRUE,0,sizeof(float)*nsize,f,
                             0,NULL,NULL);
   if (crc != CL_SUCCESS) {
      printf("clEnqueueReadBuffer Error=%i\n",crc);
      exit(1);
   }
   return;
}

/*--------------------------------------------------------------------*/
const char* pSource[] = {
"__kernel void gadd(__global float a[], __global float b[],",
"                   __global float c[], int nx) {",
"   int j;",
"   j = get_local_id(0)+get_local_size(0)*get_group_id(0);",
"   if (j < nx)",
"      a[j] = b[j] + c[j];",
"   return;",
"}",
"__kernel void vadd(__global float4 a[], __global float4 b[],",
"                   __global float4 c[], int nx) {",
"   int j, nxv;",
"   nxv = (nx - 1)/4 + 1;",
"   j = get_local_id(0)+get_local_size(0)*get_group_id(0);",
"   if (j < nxv)",
"      a[j] = b[j] + c[j];",
"   return;",
"}"
};

/*--------------------------------------------------------------------*/
void init_cl(int platf, int dev, int *irc) {
/* initialize OpenCL with device dev or selects best GPU available   */
/* searches throughs NE platforms and devices, selects the device    */
/* with the most compute units, and saves the device id devid        */
/* creates context cnt, command queue cmq, and program object clp    */
/* compiles program and creates kernel objects clks                  */
/* saves names of kernels in kname and numbers of kernels in kernels */
/* if dev is a valid device, it is used, otherwise the GPU with the  */
/* most multi-processors is selected                                 */
/* error code is modified only if there is an error */
#define NE                 4
   int maxcpus = 0, jm = -1, km = -1;
   int j, k, n;
   double z;
   char *pLog;
   cl_uint pgml = sizeof(pSource)/sizeof(pSource[0]);
   cl_uint nplts, ndevs, maxunits;
   cl_ulong msize;
   size_t maxwks[3];
   cl_platform_id plids[NE];
   cl_device_id dvids[NE];
   cl_context_properties cps[3];
/* get list of platform ids */
   crc = clGetPlatformIDs(NE,plids,&nplts);
   if (crc != CL_SUCCESS) {
      printf("clGetPlatformIDs Error=%i\n",crc);
      *irc = 1;
      return;
   }
   if (nplts > NE) {
      printf("%i platform(s) found, only %i used\n",nplts,NE);
      nplts = NE;
   }
/* get information about platforms */
   for (k = 0; k < nplts; k++) {
      crc = clGetPlatformInfo(plids[k],CL_PLATFORM_NAME,sizeof(param),
                              param,&psize);
      if (crc != CL_SUCCESS) {
         printf("clGetPlatformInfo CL_PLATFORM_NAME Error=%i\n",crc);
         printf("current param size=%li\n",sizeof(param));
      }
      else
         printf("k=%i:CL_PLATFORM_NAME=%s\n",k,param);
      crc = clGetPlatformInfo(plids[k],CL_PLATFORM_VERSION,sizeof(param),
                              param,&psize);
      if (crc != CL_SUCCESS) {
         printf("clGetPlatformInfo CL_PLATFORM_VERSION Error=%i\n",crc);
         printf("current param size=%li\n",sizeof(param));
      }
      else 
         printf("    CL_PLATFORM_VERSION=%s\n",param);
/* get list of platform ids */
      crc = clGetDeviceIDs(plids[k],CL_DEVICE_TYPE_ALL,NE,dvids,&ndevs);
      if (crc != CL_SUCCESS) {
         printf("clGetDeviceIDs CL_DEVICE_TYPE_ALL Error=%i\n",crc);
         *irc = 1;
         return;
      }
      if (ndevs > NE) {
         printf("%i device(s) found, only %i used\n",ndevs,NE);
         ndevs = NE;
      }
/* get information about devices */
      for (j = 0; j < ndevs; j++) {
         crc = clGetDeviceInfo(dvids[j],CL_DEVICE_NAME,sizeof(param),
                               param,&psize);
         if (crc != CL_SUCCESS) {
            printf("clGetDeviceInfo CL_DEVICE_NAME Error=%i\n",crc);
            printf("current param size=%li\n",sizeof(param));
            param[0] = 0;
         }
         crc = clGetDeviceInfo(dvids[j],CL_DEVICE_MAX_COMPUTE_UNITS,
                               sizeof(cl_uint),&maxunits,&psize);
         if (crc != CL_SUCCESS) {
            printf("clGetDeviceInfo CL_DEVICE_MAX_COMPUTE_UNITS Error=%i\n",
                   crc);
         }
         printf("j=%i:CL_DEVICE_NAME=%s,CL_DEVICE_MAX_COMPUTE_UNITS=%i\n",
                j,param,maxunits);
         crc = clGetDeviceInfo(dvids[j],CL_DEVICE_GLOBAL_MEM_SIZE,
                               sizeof(cl_ulong),&msize,&psize);
         if (crc != CL_SUCCESS) {
            printf("clGetDeviceInfo CL_DEVICE_GLOBAL_MEM_SIZE Error=%i\n",
                   crc);
         }
         z = ((double) msize)/1073741824.0;
         printf("    CL_DEVICE_GLOBAL_MEM_SIZE=%lu(%f GB)\n",msize,
                (float) z);
         crc = clGetDeviceInfo(dvids[j],CL_DEVICE_MAX_WORK_ITEM_SIZES,
                               3*sizeof(size_t),maxwks,&psize);
         if (crc != CL_SUCCESS) {
            printf("clGetDeviceInfo CL_DEVICE_MAX_WORK_ITEM_SIZES Error=%i\n",
                   crc);
         }
         printf("    CL_DEVICE_MAX_WORK_ITEM_SIZES=(%lu,%lu,%lu)\n",
                maxwks[0],maxwks[1],maxwks[2]);
         if (maxunits > maxcpus) {
            maxcpus = maxunits;
            jm = j;
            km = k;
         }
      }
   }
   if ((platf >= 0) && (platf < nplts) && (dev >= 0) && (dev < ndevs)) {
      km = platf;
      jm = dev;
   }
   printf("context using platform k=%i, device j=%i\n",km,jm);
   devid = dvids[jm];
/* create OpenCL context */
   cps[0] = CL_CONTEXT_PLATFORM;
   cps[1] = (cl_context_properties)plids[km];
   cps[2] = 0;
   cnt = clCreateContext(cps,1,&devid,NULL,NULL,&crc);
   if (crc != CL_SUCCESS) {
      printf("clCreateContext Error=%i\n",crc);
      *irc = 1;
      return;
   }
/* create a command-queue on device */
   cmq = clCreateCommandQueue(cnt,dvids[jm],0,&crc);
   if (crc != CL_SUCCESS) {
      printf("clCreateCommandQueue Error=%i\n",crc);
      *irc = 1;
      return;
   }
/* create program object for a context */
   clp = clCreateProgramWithSource(cnt,pgml,pSource,NULL,&crc);
   if (crc != CL_SUCCESS) {
      printf("clCreateProgramWithSource Error=%i\n",crc);
      *irc = 1;
      return;
   }
/* compile and link a program executable */
   crc = clBuildProgram(clp,0,NULL,NULL,NULL,NULL);
   if (crc != CL_SUCCESS) {
      printf("clBuildProgram Error=%i\n",crc);
      n = 2048;
      pLog = (char *) malloc(n);
      crc = clGetProgramBuildInfo(clp,devid,CL_PROGRAM_BUILD_LOG,
                                  n,pLog,&psize);
      if (crc==CL_SUCCESS)
         printf("%s\n",pLog);
      else {
         printf("clGetProgramBuildInfo Error=%i\n",crc);
         printf("current pLog size=%i\n",n);
      }
      free(pLog);
      *irc = 1;
      return;
   }
/* creates kernel objects for all kernels in program */
   crc = clCreateKernelsInProgram(clp,NKERNELS*sizeof(cl_kernel),clks,
                                  &kernels);
   if (crc != CL_SUCCESS) {
      printf("clCreateKernelsInProgram Error=%i\n",crc);
      *irc = 1;
      return;
   }
   if (kernels > NKERNELS) {
      printf("%i kernel(s) created, only %i used\n",kernels,NKERNELS);
      kernels = NKERNELS;
   }
   for (n = 0; n < kernels; n++) {
/* return information about kernel object */
      crc = clGetKernelInfo(clks[n],CL_KERNEL_FUNCTION_NAME,
                            sizeof(knames[n]),knames[n],&psize);
      if (crc != CL_SUCCESS) {
         printf("GetKernelInfo Error=%i\n",crc);
         printf("current name size=%li\n",sizeof(knames[n]));
      *irc = 1;
      return;
      }
   }
   return;
#undef NE
}

/*--------------------------------------------------------------------*/
void gpadd(float *a, float *b, float *c, int nx) {
/* Vector Add Interface for C */
   static int kid = -1;
   int n;
   size_t lwks[1], gwks[1];
   cl_mem g_a, g_b, g_c;
/* find which kernel corresponds to function name */
   if (kid < 0) {
      for (n = 0; n < kernels; n++) {
         if (!strcmp(knames[n],"gadd")) {
            kid = n;
            break;
         }
      }
   }
   if (kid < 0) {
      printf("argument name %s not found\n","gadd");
      exit(1);
   }
   g_a = (cl_mem) a;
   g_b = (cl_mem) b;
   g_c = (cl_mem) c;
/* set argument value for kernel */
   crc = clSetKernelArg(clks[kid],0,sizeof(cl_mem),(void *)&g_a);
   if (crc != CL_SUCCESS) {
      printf("arg %i:clSetKernelArg Error=%i\n",0,crc);
      exit(1);
   }
   crc = clSetKernelArg(clks[kid],1,sizeof(cl_mem),(void *)&g_b);
   if (crc != CL_SUCCESS) {
      printf("arg %i:clSetKernelArg Error=%i\n",1,crc);
      exit(1);
   }
   crc = clSetKernelArg(clks[kid],2,sizeof(cl_mem),(void *)&g_c);
   if (crc != CL_SUCCESS) {
      printf("arg %i:clSetKernelArg Error=%i\n",2,crc);
      exit(1);
   }
   crc = clSetKernelArg(clks[kid],3,sizeof(int),(void *)&nx);
   if (crc != CL_SUCCESS) {
      printf("arg %i:clSetKernelArg Error=%i\n",3,crc);
      exit(1);
   }
/* set up work sizes */
   lwks[0] = nblock_size;
   gwks[0] = ((nx - 1)/nblock_size + 1)*nblock_size;
/* enqueue command to execute kernel */
   crc = clEnqueueNDRangeKernel(cmq,clks[kid],1,NULL,gwks,lwks,0,NULL,
                                &evnt);
   if (crc != CL_SUCCESS) {
      printf("clEnqueueNDRangeKernel Error=%i\n",crc);
      exit(1);
   }
/* wait for event to complete */
   crc = clWaitForEvents(1,&evnt);
   if (crc != CL_SUCCESS) {
      printf("clWaitForEvents Error=%i\n",crc);
      exit(1);
   }
/* returns information about event */
   crc = clGetEventInfo(evnt,CL_EVENT_COMMAND_EXECUTION_STATUS,
                        sizeof(cl_int),&status,&psize);
   if (crc != CL_SUCCESS) {
      printf("clGetEventInfo Error=%i\n",crc);
   }
   if (crc != CL_COMPLETE) {
      printf("event status=%i\n",status);
   }
}

/*--------------------------------------------------------------------*/
void vpadd(float *a, float *b, float *c, int nx) {
/* Vector Add Interface for C */
   static int kid = -1;
   int n;
   size_t lwks[1], gwks[1];
   cl_mem g_a, g_b, g_c;
/* find which kernel corresponds to function name */
   if (kid < 0) {
      for (n = 0; n < kernels; n++) {
         if (!strcmp(knames[n],"vadd")) {
            kid = n;
            break;
         }
      }
   }
   if (kid < 0) {
      printf("argument name %s not found\n","vadd");
      exit(1);
   }
   g_a = (cl_mem) a;
   g_b = (cl_mem) b;
   g_c = (cl_mem) c;
/* set argument value for kernel */
   crc = clSetKernelArg(clks[kid],0,sizeof(cl_mem),(void *)&g_a);
   if (crc != CL_SUCCESS) {
      printf("arg %i:clSetKernelArg Error=%i\n",0,crc);
      exit(1);
   }
   crc = clSetKernelArg(clks[kid],1,sizeof(cl_mem),(void *)&g_b);
   if (crc != CL_SUCCESS) {
      printf("arg %i:clSetKernelArg Error=%i\n",1,crc);
      exit(1);
   }
   crc = clSetKernelArg(clks[kid],2,sizeof(cl_mem),(void *)&g_c);
   if (crc != CL_SUCCESS) {
      printf("arg %i:clSetKernelArg Error=%i\n",2,crc);
      exit(1);
   }
   crc = clSetKernelArg(clks[kid],3,sizeof(int),(void *)&nx);
   if (crc != CL_SUCCESS) {
      printf("arg %i:clSetKernelArg Error=%i\n",3,crc);
      exit(1);
   }
/* set up work sizes */
   lwks[0] = nblock_size;
   gwks[0] = ((nx - 1)/(4*nblock_size) + 1)*nblock_size;
/* enqueue command to execute kernel */
   crc = clEnqueueNDRangeKernel(cmq,clks[kid],1,NULL,gwks,lwks,0,NULL,
                                &evnt);
   if (crc != CL_SUCCESS) {
      printf("clEnqueueNDRangeKernel Error=%i\n",crc);
      exit(1);
   }
/* wait for event to complete */
   crc = clWaitForEvents(1,&evnt);
   if (crc != CL_SUCCESS) {
      printf("clWaitForEvents Error=%i\n",crc);
      exit(1);
   }
/* returns information about event */
   crc = clGetEventInfo(evnt,CL_EVENT_COMMAND_EXECUTION_STATUS,
                        sizeof(cl_int),&status,&psize);
   if (crc != CL_SUCCESS) {
      printf("clGetEventInfo Error=%i\n",crc);
   }
   if (crc != CL_COMPLETE) {
      printf("event status=%i\n",status);
   }
}

/*--------------------------------------------------------------------*/
void end_cl(int *irc) {
/* error code is modified only if there is an error       */
   int n;
/* decrement kernel reference count */
   for (n = 0; n < kernels; n++) {
      crc = clReleaseKernel(clks[n]);
      if (crc != CL_SUCCESS) {
         printf("n=%i:clReleaseKernel Error=%i\n",n,crc);
         *irc = 1;
      }
   }
/* decrement program reference count */
   crc = clReleaseProgram(clp);
   if (crc != CL_SUCCESS) {
      printf("clReleaseProgram Error=%i\n",crc);
      *irc = 1;
   }
/* decrement command queue reference count */
   crc = clReleaseCommandQueue(cmq);
   if (crc != CL_SUCCESS) {
      printf("clReleaseCommandQueue Error=%i\n",crc);
      *irc = 1;
   }
/* decrement context reference count */
   crc = clReleaseContext(cnt);
   if (crc != CL_SUCCESS) {
      printf("clReleaseContext Error=%i\n",crc);
      *irc = 1;
   }
   return;
}

#undef NAMEL
#undef NKERNELS

/* Interfaces to Fortran */

/*--------------------------------------------------------------------*/
void setgbsize_(int *nblock) {
   setgbsize(*nblock);
   return;
}

/*--------------------------------------------------------------------*/
void g_fallocate_(unsigned long *gp_f, int *nsize, int *irc) {
/* allocate global float memory on GPU, return pointer to Fortran */
   float *fptr;
   g_fallocate(&fptr,*nsize,irc);
   *gp_f = (long )fptr;
   return;
}

/*--------------------------------------------------------------------*/
void g_f4allocate_(unsigned long *gp_f, int *nsize, int *irc) {
/* allocate global vector float memory on GPU, return pointer to Fortran */
   float *fptr;
   g_fallocate(&fptr,*nsize,irc);
   *gp_f = (long )fptr;
   return;
}

/*--------------------------------------------------------------------*/
void g_deallocate_(unsigned long *gp_f, int *irc) {
/* deallocate global memory on GPU, return pointer to Fortran */
   float *f;
   f = (float *)*gp_f;
   g_deallocate(&f,irc);
   *gp_f = 0;
   return;
}

/*--------------------------------------------------------------------*/
void copyin_gmemptr_(float *f, unsigned long *gp_f, int *nsize) {
/* copy float array from main memory to global GPU memory */
   float *g_f;
   g_f = (float *)*gp_f;
   copyin_gmemptr(f,g_f,*nsize);
   return;
}

/*--------------------------------------------------------------------*/
void copyout_gmemptr_(float *f, unsigned long *gp_f, int *nsize) {
/* copy float array from global GPU memory to main memory */
   float *g_f;
   g_f = (float *)*gp_f;
   copyout_gmemptr(f,g_f,*nsize);
   return;
}

/*--------------------------------------------------------------------*/
void init_cl_(int *platf, int *dev, int *irc) {
   init_cl(*platf,*dev,irc);
   return;
}

/*--------------------------------------------------------------------*/
void gpadd_(unsigned long *gp_a, unsigned long *gp_b,
            unsigned long *gp_c, int *nx) {
/* Vector Add Interface for Fortran */
   float *a, *b, *c;
   a = (float *)*gp_a;
   b = (float *)*gp_b;
   c = (float *)*gp_c;
   gpadd(a,b,c,*nx);
}

/*--------------------------------------------------------------------*/
void vpadd_(unsigned long *gp_a, unsigned long *gp_b,
            unsigned long *gp_c, int *nx) {
/* Vector Add Interface for Fortran */
   float *a, *b, *c;
   a = (float *)*gp_a;
   b = (float *)*gp_b;
   c = (float *)*gp_c;
   vpadd(a,b,c,*nx);
}

/*--------------------------------------------------------------------*/
void end_cl_(int *irc) {
   end_cl(irc);
   return;
}
