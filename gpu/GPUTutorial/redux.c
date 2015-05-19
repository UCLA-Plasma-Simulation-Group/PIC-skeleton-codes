/*--------------------------------------------------------------------*/
/* C GPU Tutorial: Reduction        */
/* written by Viktor K. Decyk, UCLA */

#include "redux.h"

void sum0(float a[], float *sa, int nx) {
/* simple 1d sum reduction of length nx */
/* sa = sum(a)                           */
/* local data */
   int j;

   *sa = 0.0;
   for (j = 0; j < nx; j++) {
      *sa += a[j];
   }

   return;
}

void sum1(float a[], float *sa, int mx, int nx) {
/* 1d sum reductions, each of length mx */
/* sa = sum(a) */
/* local data */
   int j, js, jb, nbx, mxm;
   float t;
/* nbx = number of blocks */
   nbx = (nx - 1)/mx + 1;

   *sa = 0.0;
   for (jb = 0; jb < nbx; jb++) {
      mxm = mx < nx-mx*jb ? mx : nx-mx*jb;
      t = 0.0;
      for (js = 0; js < mxm; js++) {
         j = js + mx*jb;
         t += a[j];
      }
      *sa += t;
   }

   return;
}

void sum2(float a[], float d[], int mx, int nx) {
/* segmented 1d sum reductions, each of length mx             */
/* forall (j = 1:nbx); d(j) = sum(a(1+mx*(j-1):min(nx,mx*j))) */
/* local data */
   int j, js, jb, nbx, mxm;
   float t;
/* nbx = number of blocks */
   nbx = (nx - 1)/mx + 1;

   for (jb = 0; jb < nbx; jb++) {
      mxm = mx < nx-mx*jb ? mx : nx-mx*jb;
      t = 0.0;
      for (js = 0; js < mxm; js++) {
         j = js + mx*jb;
         t += a[j];
      }
      d[jb] = t;
   }

   return;
}

/* Interfaces to Fortran */

void csum0_(float *a, float *sa, int *nx) {
   sum0(a,sa,*nx);
   return;
}

void csum1_(float *a, float *sa, int *mx, int *nx) {
   sum1(a,sa,*mx,*nx);
   return;
}

void csum2_(float *a, float *d, int *mx, int *nx) {
   sum2(a,d,*mx,*nx);
   return;
}
