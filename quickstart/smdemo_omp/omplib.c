/* vector add test program for OpenMP  */
/* written by Viktor K. Decyk, UCLA */

#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include "omplib.h"

static int nthreads = 1;

/*--------------------------------------------------------------------*/
void mpadd(float a[], float b[], float c[], int nx) {
   int j;
#pragma omp parallel for private(j)
   for (j = 0; j < nx; j++) {
      a[j] = b[j] + c[j];
   }
   return;
}

/*--------------------------------------------------------------------*/
void init_omp(int nth) {
/* initialize openmp library */
/* use nth threads if nth > 0; otherwise, use the number found */
/* local data */
   int ncpus;
/* determine how many processors are available */
   ncpus = omp_get_num_procs();
   printf("number of cpus found = %i\n",ncpus);
   nthreads = omp_get_max_threads();
   printf("maximum number of threads = %i\n",nthreads);
   if (nth > 0)
      nthreads = nth;
   omp_set_num_threads(nthreads);
   printf("using %i thread(s)\n",nthreads);
   return;
}

/*--------------------------------------------------------------------*/
void setnthsize(int nth) {
/* set number of threads */
   if (nth > 0)
      nthreads = nth;
   omp_set_num_threads(nthreads);
   return;
}
