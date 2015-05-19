/* C Library for Skeleton 2D Electrostatic PIC Code field diagnostics */
/* Wrappers for calling the Fortran routines from a C main program */

#include <complex.h>

void potp2_(float complex *q, float complex *pot, float complex *ffc,
            float *we, int *nx, int *ny, int *nxvh, int *nyv, int *nxhd,
            int *nyhd);

void divf2_(float complex *f, float complex *df, int *nx, int *ny,
            int *nxvh, int *nyv);

void gradf2_(float complex *df, float complex *f, int *nx, int *ny,
             int *nxvh, int *nyv);

void smooth2_(float complex *q, float complex *qs, float complex *ffc,
              int *nx, int *ny, int *nxvh, int *nyv, int *nxhd,
              int *nyhd);

void rdmodes2_(float complex *pot, float complex *pott, int *nx,
               int *ny, int *modesx, int *modesy, int *nxvh, int *nyv,
               int *modesxd, int *modesyd);

void wrmodes2_(float complex *pot, float complex *pott, int *nx, 
               int *ny, int *modesx, int *modesy, int *nxvh, int *nyv,
               int *modesxd, int *modesyd);

/* Interfaces to C */

/*--------------------------------------------------------------------*/
void cpotp2(float complex q[], float complex pot[], float complex ffc[],
            float *we, int nx, int ny, int nxvh, int nyv, int nxhd,
            int nyhd) {
   potp2_(q,pot,ffc,we,&nx,&ny,&nxvh,&nyv,&nxhd,&nyhd);
   return;
}

/*--------------------------------------------------------------------*/
void cdivf2(float complex f[], float complex df[], int nx, int ny,
            int nxvh, int nyv) {
   divf2_(f,df,&nx,&ny,&nxvh,&nyv);
   return;
}

/*--------------------------------------------------------------------*/
void cgradf2(float complex df[], float complex f[], int nx, int ny,
             int nxvh, int nyv) {
   gradf2_(df,f,&nx,&ny,&nxvh,&nyv);
   return;
}

/*--------------------------------------------------------------------*/
void csmooth2(float complex q[], float complex qs[],
              float complex ffc[], int nx, int ny, int nxvh, int nyv,
              int nxhd, int nyhd) {
   smooth2_(q,qs,ffc,&nx,&ny,&nxvh,&nyv,&nxhd,&nyhd);
   return;
}

/*--------------------------------------------------------------------*/
void crdmodes2(float complex pot[], float complex pott[], int nx,
               int ny, int modesx, int modesy, int nxvh, int nyv,
               int modesxd, int modesyd) {
   rdmodes2_(pot,pott,&nx,&ny,&modesx,&modesy,&nxvh,&nyv,&modesxd,
             &modesyd);
   return;
}

/*--------------------------------------------------------------------*/
void cwrmodes2(float complex pot[], float complex pott[], int nx, 
               int ny, int modesx, int modesy, int nxvh, int nyv,
               int modesxd, int modesyd) {
   wrmodes2_(pot,pott,&nx,&ny,&modesx,&modesy,&nxvh,&nyv,&modesxd,
             &modesyd);
   return;
}

