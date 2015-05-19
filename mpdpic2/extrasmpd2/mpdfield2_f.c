/* C Library for Skeleton 2-1/2D Darwin MPI/OpenMP PIC Code field  */
/* diagnostics                                                     */
/* Wrappers for calling the Fortran routines from a C main program */

#include <complex.h>

void mppotp2_(float complex *q, float complex *pot, float complex *ffc,
              float *we, int *nx, int *ny, int *kstrt, int *nyv,
              int *kxp, int *nyhd);

void mppdivf2_(float complex *f, float complex *df, int *nx, int *ny,
               int *kstrt, int *ndim, int *nyv, int *kxp);

void mppgradf2_(float complex *df, float complex *f, int *nx, int *ny,
                int *kstrt, int *ndim, int *nyv, int *kxp);

void mppcurlf2_(float complex *f, float complex *g, int *nx, int *ny,
                int *kstrt, int *nyv, int *kxp);

void mppapotp23_(float complex *cu, float complex *axy,
                 float complex *ffc, float *ci, float *wm, int *nx,
                 int *ny, int *kstrt, int *nyv, int *kxp, int *nyhd);

void mppsmooth2_(float complex *q, float complex *qs,
                 float complex *ffc, int *nx, int *ny, int *kstrt, 
                 int *nyv, int *kxp, int *nyhd);

void mppsmooth23_(float complex *cu, float complex *cus,
                  float complex *ffc, int *nx, int *ny, int *kstrt,
                  int *nyv, int *kxp, int *nyhd);

void pprdmodes2_(float complex *pot, float complex *pott, int *nx,
                 int *ny, int *modesx, int *modesy, int *kstrt, 
                 int *nyv, int *kxp, int *modesxpd, int *modesyd);

void ppwrmodes2_(float complex *pot, float complex *pott, int *nx,
                 int *ny, int *modesx, int *modesy, int *kstrt,
                 int *nyv, int *kxp, int *modesxpd, int *modesyd);

void pprdvmodes2_(float complex *vpot, float complex *vpott, int *nx,
                  int *ny, int *modesx, int *modesy, int *ndim,
                  int *kstrt, int *nyv, int *kxp, int *modesxpd,
                  int *modesyd);

void ppwrvmodes2_(float complex *vpot, float complex *vpott, int *nx,
                  int *ny, int *modesx, int *modesy, int *ndim,
                  int *kstrt, int *nyv, int *kxp, int *modesxpd, 
                  int *modesyd);

/* Interfaces to C */

/*--------------------------------------------------------------------*/
void cmppotp2(float complex q[], float complex pot[],
              float complex ffc[], float *we, int nx, int ny, int kstrt,
              int nyv, int kxp, int nyhd) {
   mppotp2_(q,pot,ffc,we,&nx,&ny,&kstrt,&nyv,&kxp,&nyhd);
   return;
}

/*--------------------------------------------------------------------*/
void cmppdivf2(float complex f[], float complex df[], int nx, int ny,
               int kstrt, int ndim, int nyv, int kxp) {
   mppdivf2_(f,df,&nx,&ny,&kstrt,&ndim,&nyv,&kxp);
   return;
}

/*--------------------------------------------------------------------*/
void cmppgradf2(float complex df[], float complex f[], int nx, int ny,
                int kstrt, int ndim, int nyv, int kxp) {
   mppgradf2_(df,f,&nx,&ny,&kstrt,&ndim,&nyv,&kxp);
   return;
}

/*--------------------------------------------------------------------*/
void cmppcurlf2(float complex f[], float complex g[], int nx, int ny,
                int kstrt, int nyv, int kxp) {
   mppcurlf2_(f,g,&nx,&ny,&kstrt,&nyv,&kxp);
   return;
}

/*--------------------------------------------------------------------*/
void cmppapotp23(float complex cu[], float complex axy[],
                 float complex ffc[], float ci, float *wm, int nx,
                 int ny, int kstrt, int nyv, int kxp, int nyhd) {
   mppapotp23_(cu,axy,ffc,&ci,wm,&nx,&ny,&kstrt,&nyv,&kxp,&nyhd);
   return;
}

/*--------------------------------------------------------------------*/
void cmppsmooth2(float complex q[], float complex qs[],
                 float complex ffc[], int nx, int ny, int kstrt, int nyv,
                 int kxp, int nyhd) {
   mppsmooth2_(q,qs,ffc,&nx,&ny,&kstrt,&nyv,&kxp, &nyhd);
   return;
}

/*--------------------------------------------------------------------*/
void cmppsmooth23(float complex cu[], float complex cus[],
                  float complex ffc[], int nx, int ny, int kstrt,
                  int nyv, int kxp, int nyhd) {
   mppsmooth23_(cu,cus,ffc,&nx,&ny,&kstrt,&nyv,&kxp,&nyhd);
   return;
}

/*--------------------------------------------------------------------*/
void cpprdmodes2(float complex pot[], float complex pott[], int nx,
                 int ny, int modesx, int modesy, int kstrt, int nyv,
                 int kxp, int modesxpd, int modesyd) {
   pprdmodes2_(pot,pott,&nx,&ny,&modesx,&modesy,&kstrt,&nyv,&kxp,
               &modesxpd,&modesyd);
   return;
}

/*--------------------------------------------------------------------*/
void cppwrmodes2(float complex pot[], float complex pott[], int nx,
                 int ny, int modesx, int modesy, int kstrt, int nyv,
                 int kxp, int modesxpd, int modesyd) {
   ppwrmodes2_(pot,pott,&nx,&ny,&modesx,&modesy,&kstrt,&nyv,&kxp,
               &modesxpd,&modesyd);
   return;
}

/*--------------------------------------------------------------------*/
void cpprdvmodes2(float complex vpot[], float complex vpott[], int nx,
                  int ny, int modesx, int modesy, int ndim, int kstrt, 
                  int nyv, int kxp, int modesxpd, int modesyd) {
   pprdvmodes2_(vpot,vpott,&nx,&ny,&modesx,&modesy,&ndim,&kstrt,&nyv,
                &kxp,&modesxpd,&modesyd);
   return;
}

/*--------------------------------------------------------------------*/
void cppwrvmodes2(float complex vpot[], float complex vpott[], int nx,
                  int ny, int modesx, int modesy, int ndim, int kstrt, 
                  int nyv, int kxp, int modesxpd, int modesyd) {
   ppwrvmodes2_(vpot,vpott,&nx,&ny,&modesx,&modesy,&ndim,&kstrt,&nyv,
                &kxp,&modesxpd,&modesyd);
   return;
}
