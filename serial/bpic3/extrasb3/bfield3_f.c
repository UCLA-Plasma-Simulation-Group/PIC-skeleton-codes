/* C Library for Skeleton 3D Electromagnetic PIC Code field diagnostics */
/* Wrappers for calling the Fortran routines from a C main program */

#include <complex.h>

void potp3_(float complex *q, float complex *pot, float complex *ffc,
            float *we, int *nx, int *ny, int *nz, int *nxvh, int *nyv,
            int *nzv, int *nxhd, int *nyhd, int *nzhd);

void divf3_(float complex *f, float complex *df, int *nx, int *ny,
            int *nz, int *nxvh, int *nyv, int *nzv);

void gradf3_(float complex *df, float complex *f, int *nx, int *ny,
             int *nz, int *nxvh, int *nyv, int *nzv);

void curlf3_(float complex *f, float complex *g, int *nx, int *ny,
             int *nz, int *nxvh, int *nyv, int *nzv);

void avpot33_(float complex *bxyz, float complex *axyz, int *nx,
              int *ny, int *nz, int *nxvh, int *nyv, int *nzv);

void avrpot33_(float complex *axyz, float complex *bxyz,
               float complex *ffc, float *ci, int *nx, int *ny, int *nz,
               int *nxvh, int *nyv, int *nzv, int *nxhd, int *nyhd, 
               int *nzhd);

void smooth3_(float complex *q, float complex *qs, float complex *ffc,
              int *nx, int *ny, int *nz, int *nxvh, int *nyv, int *nzv,
              int *nxhd, int *nyhd, int *nzhd);

void smooth33_(float complex *cu, float complex *cus,
               float complex *ffc, int *nx, int *ny, int *nz, int *nxvh,
               int *nyv, int *nzv, int *nxhd, int *nyhd, int *nzhd);

void apotp33_(float complex *cu, float complex *axyz,
              float complex *ffc, float *ci, float *wm, int *nx,
              int *ny, int *nz, int *nxvh, int *nyv, int *nzv,
              int *nxhd, int *nyhd, int *nzhd);

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
void ccurlf3(float complex f[], float complex g[], int nx, int ny,
             int nz, int nxvh, int nyv, int nzv) {
   curlf3_(f,g,&nx,&ny,&nz,&nxvh,&nyv,&nzv);
   return;
}

/*--------------------------------------------------------------------*/
void cavpot33(float complex bxyz[], float complex axyz[], int nx,
              int ny, int nz, int nxvh, int nyv, int nzv) {
   avpot33_(bxyz,axyz,&nx,&ny,&nz,&nxvh,&nyv,&nzv);
   return;
}

/*--------------------------------------------------------------------*/
void cavrpot33(float complex axyz[], float complex bxyz[],
               float complex ffc[], float ci, int nx, int ny, int nz,
               int nxvh, int nyv, int nzv, int nxhd, int nyhd, 
               int nzhd) {
   avrpot33_(axyz,bxyz,ffc,&ci,&nx,&ny,&nz,&nxvh,&nyv,&nzv,&nxhd,&nyhd,
             &nzhd);
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
void csmooth33(float complex cu[], float complex cus[],
               float complex ffc[], int nx, int ny, int nz, int nxvh,
               int nyv, int nzv, int nxhd, int nyhd, int nzhd) {
   smooth33_(cu,cus,ffc,&nx,&ny,&nz,&nxvh,&nyv,&nzv,&nxhd,&nyhd,&nzhd);
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

