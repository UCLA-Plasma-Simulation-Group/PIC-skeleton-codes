/* AVX512 utility Library */
/* written by Viktor K. Decyk, UCLA and Ricardo Fonseca, ISCTE */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <complex.h>
#include <immintrin.h>
#include <mm_malloc.h>
#include "avx512lib3.h"

/*--------------------------------------------------------------------*/
void avx512_fallocate(float **s_f, int nsize, int *irc) {
/* allocate aligned float memory on AVX512 return pointer to C */
/* size is padded to be a multiple of the alignment length     */
/* local data */
/* NV = vector length for 32 bit data */
#define NV             16
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
void avx512_callocate(float complex **s_c, int nsize, int *irc) {
/* allocate aligned float complex memory on AVX512 return pointer to C */
/* size is padded to be a multiple of the alignment length             */
/* local data */
/* NV = vector length for 64 bit data */
#define NV             8
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
void avx512_iallocate(int **s_i, int nsize, int *irc) {
/* allocate aligned int memory on AVX512, return pointer to C */
/* size is padded to be a multiple of the alignment length    */
/* local data */
/* NV = vector length for 32 bit data */
#define NV             16
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
void avx512_deallocate(void *s_d) {
/* deallocate aligned memory on AVX512 */
   _mm_free(s_d);
   return;
}

/* Interfaces to Fortran */

/*--------------------------------------------------------------------*/
void avx512_deallocate_(void *sp_d) {
/* pointer in Fortran should also be nullified */
   avx512_deallocate(sp_d);
   return;
}

void fcopyin_(float *f, float *g, int *n) {
   int j;
   for (j = 0; j < *n; j++) {
      f[j] = g[j];
   }
   return;
}
