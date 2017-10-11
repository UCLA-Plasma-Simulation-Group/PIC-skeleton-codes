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

/*--------------------------------------------------------------------*/
void cknciscan2(int *isdata, int nths) {
/* performs local prefix reduction of integer data shared by threads */
/* using binary tree method. */
/* requires KNC, isdata needs to be 64 byte aligned */
/* local data */
   int j, ns, isum, ist;
   __m512i v_m1, v_m2, v_it, v_is, v_ioff;
   ns = 16*(nths/16);
   v_m1 = _mm512_set_epi32(11,11,11,11,11,10,9,8,3,3,3,3,3,2,1,0);
   v_m2 = _mm512_set_epi32(7,7,7,7,7,7,7,7,7,6,5,4,3,2,1,0);
   isum = 0;
   v_ioff = _mm512_setzero_epi32();
/* vector loop over elements in blocks of 16 */
   for (j = 0; j < ns; j+=16) {
/* load data */
      v_it = _mm512_load_epi32(&isdata[j]);
/* first pass */
      v_is = _mm512_shuffle_epi32(v_it,177);
      v_it = _mm512_mask_add_epi32(v_it,_mm512_int2mask(43690),v_it,
             v_is);
/* second pass */
      v_is = _mm512_shuffle_epi32(v_it,80);
      v_it = _mm512_mask_add_epi32(v_it,_mm512_int2mask(52428),v_it,
             v_is);
/* third pass */
      v_is = _mm512_permutevar_epi32(v_m1,v_it);
      v_it = _mm512_mask_add_epi32(v_it,_mm512_int2mask(61680),v_it,
             v_is);
/* fourth pass */
      v_is = _mm512_permutevar_epi32(v_m2,v_it);
      v_it = _mm512_mask_add_epi32(v_it,_mm512_int2mask(65280),v_it,
             v_is);
/* add offset */
      v_it = _mm512_add_epi32(v_it,v_ioff);
/* next offset */
      v_ioff = _mm512_shuffle_epi32(v_it,255);
      v_ioff = _mm512_permute4f128_epi32(v_ioff,255);
/* write data */
      _mm512_store_epi32(&isdata[j],v_it);
   }
   if (ns > 0)
      isum = isdata[ns-1];
/* loop over remaining elements */
   for (j = ns; j < nths; j++) {
      ist = isdata[j];
      isum += ist;
      isdata[j] = isum;
   }
   return;
}

/* Interfaces to Fortran */

/*--------------------------------------------------------------------*/
void avx512_deallocate_(void *sp_d) {
/* pointer in Fortran should also be nullified */
   avx512_deallocate(sp_d);
   return;
}

/*--------------------------------------------------------------------*/
void cknciscan2_(int *isdata, int *nths) {
   cknciscan2(isdata,*nths);
   return;
}

void fcopyin_(float *f, float *g, int *n) {
   int j;
   for (j = 0; j < *n; j++) {
      f[j] = g[j];
   }
   return;
}
