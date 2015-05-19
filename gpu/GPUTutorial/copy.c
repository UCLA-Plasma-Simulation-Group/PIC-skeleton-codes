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
   int j, js, jb, nbx, mxm;
/* nbx = number of blocks */
   nbx = (nx - 1)/mx + 1;

   for (jb = 0; jb < nbx; jb++) {
      mxm = mx < nx-mx*jb ? mx : nx-mx*jb;
      for (js = 0; js < mxm; js++) {
         j = js + mx*jb;
         a[j] = b[j];
      }
   }

   return;
}

void copy2(float a[], float b[], int mx, int nx, int ny) {
/* segmented 2d copy of length nx, ny, with block size mx */
/* a = b                                                  */
/* local data */
   int j, k, js, jb, nbx, mxm;
/* nbx = number of blocks in x */
   nbx = (nx - 1)/mx + 1;

   for (k = 0; k < ny; k++) {
      for (jb = 0; jb < nbx; jb++) {
         mxm = mx < nx-mx*jb ? mx : nx-mx*jb;
         for (js = 0; js < mxm; js++) {
            j = js + mx*jb;
            a[j+nx*k] = b[j+nx*k];
         }
      }
   }

   return;
}

void saxpy2(float a[], float b[], float s, int mx, int nx, int ny) {
/* segmented 2d vector multiply of length nx, ny, with block size mx */
/* a = s*b + a                                                       */
/* local data */
   int j, k, js, jb, nbx, mxm;
/* nbx = number of blocks in x */
   nbx = (nx - 1)/mx + 1;

   for (k = 0; k < ny; k++) {
      for (jb = 0; jb < nbx; jb++) {
         mxm = mx < nx-mx*jb ? mx : nx-mx*jb;
         for (js = 0; js < mxm; js++) {
            j = js + mx*jb;
            a[j+nx*k] = s*b[j+nx*k] + a[j+nx*k];
         }
      }
   }

   return;
}

void copy3(float a[], float b[], int mx, int my, int nx, int ny) {
/* segmented 2d copy of length nx, ny, with block size mx, my */
/* a = b                                                      */
/* local data */
   int j, k, js, ks, jb, kb, nbx, nby, mxm, mym;
/* nbx/nby = number of blocks in x/y */
   nbx = (nx - 1)/mx + 1; nby = (ny - 1)/my + 1;

   for (kb = 0; kb < nby; kb++) {
      mym = my < ny-my*kb ? my : ny-my*kb;
      for (jb = 0; jb < nbx; jb++) {
         mxm = mx < nx-mx*jb ? mx : nx-mx*jb;
         for (ks = 0; ks < mym; ks++) {
            k = ks + my*kb;
            for (js = 0; js < mxm; js++) {
               j = js + mx*jb;
               a[j+nx*k] = b[j+nx*k];
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
