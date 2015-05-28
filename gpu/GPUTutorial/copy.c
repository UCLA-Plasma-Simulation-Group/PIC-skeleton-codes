/*--------------------------------------------------------------------*/
/* C GPU Tutorial: Copy             */
/* written by Viktor K. Decyk, UCLA */

#include "copy.h"

void copy0(float a[], float b[], int nx) {
/* simple 1d copy of length nx */
/* a = b                       */
/* nx = size of arrays in x    */
/* local data */
   int j;

   for (j = 0; j < nx; j++) {
      a[j] = b[j];
   }

   return;
}

void copy1(float a[], float b[], int mx, int nx) {
/* segmented 1d copy of length nx, with block size mx */
/* a = b                                              */
/* local data */
   int j, id, nbx, joff, mxm;
/* nbx = number of blocks */
   nbx = (nx - 1)/mx + 1;

   for (id = 0; id < nbx; id++) {
      joff = mx*id;
      mxm = mx < nx-joff ? mx : nx-joff;
      for (j = 0; j < mxm; j++) {
         a[j+joff] = b[j+joff];
      }
   }

   return;
}

void copy2(float a[], float b[], int mx, int nx, int ny) {
/* segmented 2d copy of length nx, ny, with block size mx */
/* a = b                                                  */
/* local data */
   int j, k, id, nbx, joff, mxm;
/* nbx = number of blocks in x */
   nbx = (nx - 1)/mx + 1;

   for (k = 0; k < ny; k++) {
      for (id = 0; id < nbx; id++) {
         joff = mx*id;
         mxm = mx < nx-joff ? mx : nx-joff;
         for (j = 0; j < mxm; j++) {
            a[j+joff+nx*k] = b[j+joff+nx*k];
         }
      }
   }

   return;
}

void saxpy2(float a[], float b[], float s, int mx, int nx, int ny) {
/* segmented 2d vector multiply of length nx, ny, with block size mx */
/* a = s*b + a                                                       */
/* local data */
   int j, k, id, nbx, joff, mxm;
/* nbx = number of blocks in x */
   nbx = (nx - 1)/mx + 1;

   for (k = 0; k < ny; k++) {
      for (id = 0; id < nbx; id++) {
         joff = mx*id;
         mxm = mx < nx-joff ? mx : nx-joff;
         for (j = 0; j < mxm; j++) {
            a[j+joff+nx*k] = s*b[j+joff+nx*k] + a[j+joff+nx*k];
         }
      }
   }

   return;
}

void copy3(float a[], float b[], int mx, int my, int nx, int ny) {
/* segmented 2d copy of length nx, ny, with block size mx, my */
/* a = b                                                      */
/* local data */
   int j, k, idx, idy, nbx, nby, joff, koff, mxm, mym;
/* nbx/nby = number of blocks in x/y */
   nbx = (nx - 1)/mx + 1; nby = (ny - 1)/my + 1;

   for (idy = 0; idy < nby; idy++) {
      koff = my*idy;
      mym = my < ny-koff ? my : ny-koff;
      for (idx = 0; idx < nbx; idx++) {
         joff = mx*idx;
         mxm = mx < nx-joff ? mx : nx-joff;
         for (k = 0; k < mym; k++) {
            for (j = 0; j < mxm; j++) {
               a[j+joff+nx*(k+koff)] = b[j+joff+nx*(k+koff)];
            }
         }
      }
   }

   return;
}

/* Interfaces to Fortran */

void ccopy0_(float *a, float *b, int *nx) {
   copy0(a,b,*nx);
   return;
}

void ccopy1_(float *a, float *b, int *mx, int *nx) {
   copy1(a,b,*mx,*nx);
   return;
}

void ccopy2_(float *a, float *b, int *mx, int *nx, int *ny) {
   copy2(a,b,*mx,*nx,*ny);
   return;
}

void csaxpy2_(float *a, float *b, float *s, int *mx, int *nx, int *ny) {
   saxpy2(a,b,*s,*mx,*nx,*ny);
   return;
}

void ccopy3_(float *a, float *b, int *mx, int *my, int *nx, int *ny) {
copy3(a,b,*mx,*my,*nx,*ny);
   return;
}
