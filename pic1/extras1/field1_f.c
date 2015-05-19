/* C Library for Skeleton 1D Electrostatic PIC Code field diagnostics */
/* Wrappers for calling the Fortran routines from a C main program    */

#include <complex.h>

void potp1_(float complex *q, float complex *pot, float complex *ffc,
            float *we, int *nx, int *nxvh, int *nxhd);

void divf1_(float complex *f, float complex *df, int *nx, int *ndim,
            int *nxvh);

void gradf1_(float complex *df, float complex *f, int *nx, int *ndim,
             int *nxvh);

void smooth1_(float complex *q, float complex *qs, float complex *ffc,
              int *nx, int *nxvh, int *nxhd);

void rdmodes1_(float complex *pot, float complex *pott, int *nx,
               int *modesx, int *nxvh, int *modesxd);

void wrmodes1_(float complex *pot, float complex *pott, int *nx,
               int *modesx, int *nxvh, int *modesxd);

/* Interfaces to C */

/*--------------------------------------------------------------------*/
void cpotp1(float complex q[], float complex pot[], float complex ffc[],
            float *we, int nx, int nxvh, int nxhd) {
   potp1_(q,pot,ffc,we,&nx,&nxvh,&nxhd);
   return;
}

/*--------------------------------------------------------------------*/
void cdivf1(float complex f[], float complex df[], int nx, int ndim,
            int nxvh) {
   divf1_(f,df,&nx,&ndim,&nxvh);
   return;
}

/*--------------------------------------------------------------------*/
void cgradf1(float complex df[], float complex f[], int nx, int ndim,
             int nxvh) {
   gradf1_(df,f,&nx,&ndim,&nxvh);
   return;
}

/*--------------------------------------------------------------------*/
void csmooth1(float complex q[], float complex qs[],
              float complex ffc[], int nx, int nxvh, int nxhd) {
   smooth1_(q,qs,ffc,&nx,&nxvh,&nxhd);
   return;
}

/*--------------------------------------------------------------------*/
void crdmodes1(float complex pot[], float complex pott[], int nx,
               int modesx, int nxvh, int modesxd) {
   rdmodes1_(pot,pott,&nx,&modesx,&nxvh,&modesxd);
   return;
}

/*--------------------------------------------------------------------*/
void cwrmodes1(float complex pot[], float complex pott[], int nx,
               int modesx, int nxvh, int modesxd) {
   wrmodes1_(pot,pott,&nx,&modesx,&nxvh,&modesxd);
   return;
}

