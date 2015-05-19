/*--------------------------------------------------------------------*/
/* C Library GPU Tutorial            */
/* written by Viktor K. Decyk, UCLA  */

#include "transpose.h"

void transpose0(float a[], float b[], int nx, int ny) {
/* simple 2d transpose of length nx, ny */
/* a = transpose(b)                     */
/* local data */
   int j, k;

   for (k = 0; k < ny; k++) {
      for (j = 0; j < nx; j++) {
         a[k+ny*j] = b[j+nx*k];
      }
   }

   return;
}

void  transpose2(float a[], float b[], int mx, int my, int nx, int ny) {
/* segmented 2d transpose of length nx, ny, with block size mx, my */
/* a = transpose(b)                                                */
/* local data */
   int j, k, js, ks, jb, kb, joff, koff, nbx, nby, mxv, mxm, mym;
/* scratch fast, local array */
   float s[(mx+1)*my];
/* nbx/nby = number of blocks in x/y */
   nbx = (nx - 1)/mx + 1; nby = (ny - 1)/my + 1;
   mxv = mx + 1;

   for (kb = 0; kb < nby; kb++) {
      koff = my*kb;
      mym = my < ny-koff ? my : ny-koff;
      for (jb = 0; jb < nbx; jb++) {
         joff = mx*jb;
         mxm = mx < nx-joff ? mx : nx-joff;

/* copy in from slow memory with stride 1 into block of fast memory */
         for (ks = 0; ks < mym; ks++) {
            k = ks + koff;
            for (js = 0; js < mxm; js++) {
               j = js + joff;
               s[js+mxv*ks] = b[j+nx*k];
            }
         }

/* copy out to slow memory with stride 1 from block of fast memory */
         for (js = 0; js < mxm; js++) {
            j = js + joff;
            for (ks = 0; ks < mym; ks++) {
               k = ks + koff;
               a[k+ny*j] = s[js+mxv*ks];
            }
         }

      }
   }

   return;
}

/* Interfaces to Fortran */

void ctranspose0_(float *a, float *b, int *nx, int *ny) {
   transpose0(a,b,*nx,*ny);
   return;
}

void ctranspose2_(float *a, float *b, int *mx, int *my, int *nx,
                  int *ny) {
   transpose2(a,b,*mx,*my,*nx,*ny);
   return;
}
