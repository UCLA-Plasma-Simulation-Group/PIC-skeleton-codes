/* SSE2 utility Library */
/* written by Viktor K. Decyk, UCLA and Ricardo Fonseca, ISCTE */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <complex.h>
#include <emmintrin.h>
#include <mm_malloc.h>
#include "sselib3.h"

/* the following two functions work only with Intel icc compiler */
/*
int check_sse2() {
   int result;
   result = _may_i_use_cpu_feature(_FEATURE_SSE2);
   return result;
}

int check_avx() {
   int result;
   result = _may_i_use_cpu_feature(_FEATURE_AVX);
   return result;
}
*/

/*--------------------------------------------------------------------*/
void sse_fallocate(float **s_f, int nsize, int *irc) {
/* allocate aligned float memory on SSE return pointer to C */
/* size is padded to be a multiple of the alignment length */
/* local data */
/* NV = vector length for 32 bit data */
#define NV             4
   int ns;
   void *sptr = NULL;
   ns = NV*((nsize - 1)/NV + 1);
   sptr = _mm_malloc(ns*sizeof(float),4*NV);
   if (sptr==NULL) {
      printf("_mm_malloc float Error,len=%d\n",ns);
      *irc = 1;
   }
   *s_f = (float *)sptr;
   return;
#undef NV
}

/*--------------------------------------------------------------------*/
void sse_callocate(float complex **s_c, int nsize, int *irc) {
/* allocate aligned float complex memory on SSE return pointer to C */
/* size is padded to be a multiple of the alignment length          */
/* local data */
/* NV = vector length for 64 bit data */
#define NV             2
   int ns;
   void *sptr = NULL;
   ns = NV*((nsize - 1)/NV + 1);
   sptr = _mm_malloc(ns*sizeof(float complex),8*NV);
   if (sptr==NULL) {
      printf("_mm_malloc float complex Error,len=%d\n",ns);
      *irc = 1;
   }
   *s_c = (float complex *)sptr;
   return;
#undef NV
}

/*--------------------------------------------------------------------*/
void sse_iallocate(int **s_i, int nsize, int *irc) {
/* allocate aligned int memory on SSE, return pointer to C */
/* size is padded to be a multiple of the alignment length */
/* local data */
/* NV = vector length for 32 bit data */
#define NV             4
   int ns;
   void *sptr = NULL;
   ns = NV*((nsize - 1)/NV + 1);
   sptr = _mm_malloc(ns*sizeof(int),4*NV);
   if (sptr==NULL) {
      printf("_mm_malloc int Error,len=%d\n",ns);
      *irc = 1;
   }
   *s_i = (int *)sptr;
   return;
#undef NV
}

/*--------------------------------------------------------------------*/
void sse_deallocate(void *s_d) {
/* deallocate aligned memory on SSE */
   _mm_free(s_d);
   return;
}

/* Interfaces to Fortran */

/*--------------------------------------------------------------------*/
void sse_deallocate_(void *sp_d) {
/* pointer in Fortran should also be nullified */
   sse_deallocate(sp_d);
   return;
}

/* works only with Intel icc compiler */
/*
int check_sse2_() {
   return check_sse2();
}

int check_avx_() {
   return check_avx();
}
*/

void fcopyin_(float *f, float *g, int *n) {
   int j;
   for (j = 0; j < *n; j++) {
      f[j] = g[j];
   }
   return;
}

