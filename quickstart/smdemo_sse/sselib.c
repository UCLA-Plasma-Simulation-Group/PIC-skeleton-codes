/* SSE vector add test program */
/* written by Viktor K. Decyk, UCLA and Ricardo Fonseca, IST */

#include <stdlib.h>
#include <stdio.h>
#include <xmmintrin.h>
#include <mm_malloc.h>
#include "sselib.h"

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
void sse_deallocate(void *s_d) {
/* deallocate aligned memory on SSE */
   _mm_free(s_d);
   return;
}

/*--------------------------------------------------------------------*/
void ssadd(float s_a[], float s_b[], float s_c[], int nx) {
   int i;
    __m128 va, vb, vc;
   for (i = 0; i < nx; i+=4) {
      vb = _mm_load_ps(&s_b[i]);
      vc = _mm_load_ps(&s_c[i]);
      va = _mm_add_ps(vb,vc);
     _mm_store_ps(&s_a[i],va);
   }
   return;
}

/* Interfaces to Fortran77 */

/*--------------------------------------------------------------------*/
void sse_deallocate_(void *sp_d) {
/* pointer in Fortran should also be nullified */
   sse_deallocate(sp_d);
   return;
}

/*--------------------------------------------------------------------*/
/*
void ssadd_(unsigned long *gp_a, unsigned long *gp_b,
           unsigned long *gp_c, int *nx) {
*/
/* Vector Add Interface for Fortran */
/*
   float *g_a, *g_b, *g_c;
   g_a = (float *)*gp_a;
   g_b = (float *)*gp_b;
   g_c = (float *)*gp_c;
   ssadd(g_a,g_b,g_c,*nx);
}
*/

void ptrdump_(float s_f[]) {
   printf("dump s_f=%p\n",s_f);
   return;
}