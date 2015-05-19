/* C Library for Skeleton 3D Electrostatic PIC Code field diagnostics */
/* Wrappers for calling the Fortran routines from a C main program */

#include <complex.h>

void potp3_(float complex *q, float complex *pot, float complex *ffc,
            float *we, int *nx, int *ny, int *nz, int *nxvh, int *nyv,
            int *nzv, int *nxhd, int *nyhd, int *nzhd);

void divf3_(float complex *f, float complex *df, int *nx, int *ny,
            int *nz, int *nxvh, int *nyv, int *nzv);

void gradf3_(float complex *df, float complex *f, int *nx, int *ny,
             int *nz, int *nxvh, int *nyv, int *nzv);

void smooth3_(float complex *q, float complex *qs, float complex *ffc,
              int *nx, int *ny, int *nz, int *nxvh, int *nyv, int *nzv,
              int *nxhd, int *nyhd, int *nzhd);

void rdmodes3_(float complex *pot, float complex *pott, int *nx,
               int *ny, int *nz, int *modesx, int *modesy, int *modesz,
               int *nxvh, int *nyv, int *nzv, int *modesxd,
               int *modesyd, int *modeszd);

void wrmodes3_(float complex *pot, float complex *pott, int *nx, 
               int *ny, int *nz, int *modesx, int *modesy, int *modesz,
               int *nxvh, int *nyv, int *nzv, int *modesxd,
               int *modesyd, int *modeszd);

/* Interfaces to C */

/*--------------------------------------------------------------------*/
void cpotp3(float complex q[], float complex pot[], float complex ffc[],
            float *we, int nx, int ny, int nz, int nxvh, int nyv,
            int nzv, int nxhd, int nyhd, int nzhd) {
   potp3_(q,pot,ffc,we,&nx,&ny,&nz,&nxvh,&nyv,&nzv,&nxhd,&nyhd,&nzhd);
   return;
}

/*--------------------------------------------------------------------*/
void cdivf3(float complex f[], float complex df[], int nx, int ny,
            int nz, int nxvh, int nyv, int nzv) {
   divf3_(f,df,&nx,&ny,&nz,&nxvh,&nyv,&nzv);
   return;
}

/*--------------------------------------------------------------------*/
void cgradf3(float complex df[], float complex f[], int nx, int ny,
             int nz, int nxvh, int nyv, int nzv) {
   gradf3_(df,f,&nx,&ny,&nz,&nxvh,&nyv,&nzv);
   return;
}

/*--------------------------------------------------------------------*/
void csmooth3(float complex q[], float complex qs[],
              float complex ffc[], int nx, int ny, int nz, int nxvh,
              int nyv, int nzv, int nxhd, int nyhd, int nzhd) {
   smooth3_(q,qs,ffc,&nx,&ny,&nz,&nxvh,&nyv,&nzv,&nxhd,&nyhd,&nzhd);
   return;
}

/*--------------------------------------------------------------------*/
void crdmodes3(float complex pot[], float complex pott[], int nx,
               int ny, int nz, int modesx, int modesy, int modesz,
               int nxvh, int nyv, int nzv, int modesxd, int modesyd,
               int modeszd) {
   rdmodes3_(pot,pott,&nx,&ny,&nz,&modesx,&modesy,&modesz,&nxvh,&nyv,
             &nzv,&modesxd,&modesyd,&modeszd);
   return;
}

/*--------------------------------------------------------------------*/
void cwrmodes3(float complex pot[], float complex pott[], int nx, 
               int ny, int nz, int modesx, int modesy, int modesz,
               int nxvh, int nyv, int nzv, int modesxd, int modesyd,
               int modeszd) {
   wrmodes3_(pot,pott,&nx,&ny,&nz,&modesx,&modesy,&modesz,&nxvh,&nyv,
             &nzv,&modesxd,&modesyd,&modeszd);
   return;
}
