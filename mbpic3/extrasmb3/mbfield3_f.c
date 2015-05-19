/* C Library for Skeleton 3D Electromagnetic OpenMP PIC Code field */
/* diagnostics                                                     */
/* Wrappers for calling the Fortran routines from a C main program */

#include <complex.h>

void mpotp3_(float complex *q, float complex *pot, float complex *ffc,
             float *we, int *nx, int *ny, int *nz, int *nxvh, int *nyv,
             int *nzv, int *nxhd, int *nyhd, int *nzhd);

void mdivf3_(float complex *f, float complex *df, int *nx, int *ny,
             int *nz, int *nxvh, int *nyv, int *nzv);

void mgradf3_(float complex *df, float complex *f, int *nx, int *ny,
              int *nz, int *nxvh, int *nyv, int *nzv);

void mcurlf3_(float complex *f, float complex *g, int *nx, int *ny,
              int *nz, int *nxvh, int *nyv, int *nzv);

void mavpot33_(float complex *bxyz, float complex *axyz, int *nx,
               int *ny, int *nz, int *nxvh, int *nyv, int *nzv);

void mavrpot33_(float complex *axyz, float complex *bxyz,
                float complex *ffc, float *ci, int *nx, int *ny,
                int *nz, int *nxvh, int *nyv, int *nzv, int *nxhd,
                int *nyhd, int *nzhd);

void msmooth3_(float complex *q, float complex *qs, float complex *ffc,
               int *nx, int *ny, int *nz, int *nxvh, int *nyv, int *nzv,
               int *nxhd, int *nyhd, int *nzhd);

void msmooth33_(float complex *cu, float complex *cus,
                float complex *ffc, int *nx, int *ny, int *nz,
                int *nxvh, int *nyv, int *nzv, int *nxhd, int *nyhd,
                int *nzhd);

void rdmodes3_(float complex *pot, float complex *pott, int *nx,
               int *ny, int *nz, int *modesx, int *modesy, int *modesz,
               int *nxvh, int *nyv, int *nzv, int *modesxd,
               int *modesyd, int *modeszd);

void wrmodes3_(float complex *pot, float complex *pott, int *nx, 
               int *ny, int *nz, int *modesx, int *modesy, int *modesz,
               int *nxvh, int *nyv, int *nzv, int *modesxd,
               int *modesyd, int *modeszd);

void rdvmodes3_(float complex *vpot, float complex *vpott, int *nx,
                int *ny, int *nz, int *modesx, int *modesy, int *modesz,
                int *ndim, int *nxvh, int *nyv, int *nzv, int *modesxd,
                int *modesyd, int *modeszd);

void wrvmodes3_(float complex *vpot, float complex *vpott, int *nx,
                int *ny, int *nz, int *modesx, int *modesy, int *modesz,
                int *ndim, int *nxvh, int *nyv, int *nzv, int *modesxd, 
                int *modesyd, int *modeszd);

/* Interfaces to C */

/*--------------------------------------------------------------------*/
void cmpotp3(float complex q[], float complex pot[],
             float complex ffc[], float *we, int nx, int ny, int nz,
             int nxvh, int nyv, int nzv, int nxhd, int nyhd, int nzhd) {
   mpotp3_(q,pot,ffc,we,&nx,&ny,&nz,&nxvh,&nyv,&nzv,&nxhd,&nyhd,&nzhd);
   return;
}

/*--------------------------------------------------------------------*/
void cmdivf3(float complex f[], float complex df[], int nx, int ny,
             int nz, int nxvh, int nyv, int nzv) {
   mdivf3_(f,df,&nx,&ny,&nz,&nxvh,&nyv,&nzv);
   return;
}

/*--------------------------------------------------------------------*/
void cgmradf3(float complex df[], float complex f[], int nx, int ny,
              int nz, int nxvh, int nyv, int nzv) {
   mgradf3_(df,f,&nx,&ny,&nz,&nxvh,&nyv,&nzv);
   return;
}

/*--------------------------------------------------------------------*/
void cmcurlf3(float complex f[], float complex g[], int nx, int ny,
              int nz, int nxvh, int nyv, int nzv) {
   mcurlf3_(f,g,&nx,&ny,&nz,&nxvh,&nyv,&nzv);
   return;
}

/*--------------------------------------------------------------------*/
void cmavpot33(float complex bxyz[], float complex axyz[], int nx,
               int ny, int nz, int nxvh, int nyv, int nzv) {
   mavpot33_(bxyz,axyz,&nx,&ny,&nz,&nxvh,&nyv,&nzv);
   return;
}

/*--------------------------------------------------------------------*/
void cmavrpot33(float complex axyz[], float complex bxyz[],
                float complex ffc[], float ci, int nx, int ny, int nz,
                int nxvh, int nyv, int nzv, int nxhd, int nyhd,
                int nzhd) {
   mavrpot33_(axyz,bxyz,ffc,&ci,&nx,&ny,&nz,&nxvh,&nyv,&nzv,&nxhd,&nyhd,
              &nzhd);
   return;
}

/*--------------------------------------------------------------------*/
void cmsmooth3(float complex q[], float complex qs[],
               float complex ffc[], int nx, int ny, int nz, int nxvh,
               int nyv, int nzv, int nxhd, int nyhd, int nzhd) {
   msmooth3_(q,qs,ffc,&nx,&ny,&nz,&nxvh,&nyv,&nzv,&nxhd,&nyhd,&nzhd);
   return;
}

/*--------------------------------------------------------------------*/
void cmsmooth33(float complex cu[], float complex cus[],
                float complex ffc[], int nx, int ny, int nz, int nxvh,
                int nyv, int nzv, int nxhd, int nyhd, int nzhd) {
   msmooth33_(cu,cus,ffc,&nx,&ny,&nz,&nxvh,&nyv,&nzv,&nxhd,&nyhd,&nzhd);
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

/*--------------------------------------------------------------------*/
void crdvmodes3(float complex vpot[], float complex vpott[], int nx,
                int ny, int nz, int modesx, int modesy, int modesz,
                int ndim, int nxvh, int nyv, int nzv, int modesxd,
                int modesyd, int modeszd) {
   rdvmodes3_(vpot,vpott,&nx,&ny,&nz,&modesx,&modesy,&modesz,&ndim,
              &nxvh,&nyv,&nzv,&modesxd,&modesyd,&modeszd);
   return;
}

/*--------------------------------------------------------------------*/
void cwrvmodes3(float complex vpot[], float complex vpott[], int nx,
                int ny, int nz, int modesx, int modesy, int modesz,
                int ndim, int nxvh, int nyv, int nzv, int modesxd, 
                int modesyd, int modeszd) {
   wrvmodes3_(vpot,vpott,&nx,&ny,&nz,&modesx,&modesy,&modesz,&ndim,
              &nxvh,&nyv,&nzv,&modesxd,&modesyd,&modeszd);
   return;
}

