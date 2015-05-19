/* vector add test program for Pthreads */
/* written by Viktor K. Decyk, UCLA */

#include <stdlib.h>
#include <stdio.h>
#include "LnxMP.h"

#define MAXTHREADS                 16

static int nthreads = 1;
static int idtask[MAXTHREADS];

/*--------------------------------------------------------------------*/
void padd(float a[], float b[], float c[], int *nx) {
   int j, lx;
   lx = *nx;
   for (j = 0; j < lx; j++) {
      a[j] = b[j] + c[j];
   }
   return;
}

/*--------------------------------------------------------------------*/
void init_pt(int nth, int *irc) {
/* initialize multi-tasking library */
/* use nth threads if nth > 0; otherwise, use the number found */
/* error code is modified only if there is an error */
/* local data */
   int ncpus;
/* determine how many processors are available */
   MP_Init(&ncpus,irc);
   nthreads = ncpus;
   if (nthreads==0)
      nthreads = 1;
   printf("number of cpus found = %i\n",nthreads);
   if (nth > 0)
      nthreads = nth;
   printf("using %i cpu(s)\n",nthreads);
   return;
}

/*--------------------------------------------------------------------*/
void setnthsize(int nth) {
/* set number of threads */
   if (nth > 0)
      nthreads = nth;
   return;
}

/*--------------------------------------------------------------------*/
void copy_memptr(float *f, float *p_f, int nsize) {
/* copyin float array */
   int j;
   for (j = 0; j < nsize; j++) {
      p_f[j] = f[j];
   }
   return;
}

/*--------------------------------------------------------------------*/
void ptadd(float *a, float *b, float *c, int nx, int *irc) {
/*  multitasking vector add */
/*  irc = ierror indicator (0 = no error) */
/* local data */
   int i, nmt, nxp, nxo;
   int nargs;
   nargs = 4;
   nmt = nthreads - 1;
   nxp = (nx - 1)/nthreads + 1;
/* start add tasks */
   for (i = 0; i < nmt; i++) {
      nxo = nxp*i;
      MP_Taskstart(&idtask[i],&padd,&nargs,&a[nxo],&b[nxo],&c[nxo],&nxp);
/* check for errors */
      if (idtask[i]==0) {
         *irc = -1;
         return;
      }
   }
/* add remaining data */
   nxo = nxp*nmt;
   nxp = nx - nxp*nmt;
   padd(&a[nxo],&b[nxo],&c[nxo],&nxp);
/* wait for tasks to complete */
   for (i = 0; i < nmt; i++) {
      MP_Taskwait(&idtask[i]);
      if (idtask[i] != 0)
         *irc = -2;
   }
   return;
}

/*--------------------------------------------------------------------*/
void end_pt() {
/* terminate pthreads library */
   MP_End();
   return;
}
