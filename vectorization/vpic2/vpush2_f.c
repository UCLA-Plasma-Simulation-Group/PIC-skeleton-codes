/* C Library for Skeleton 2D Electrostatic Vector PIC Code */
/* Wrappers for calling the Fortran routines from a C main program */

#include <complex.h>

double ranorm_();

double randum_();

void distr2_(float *part, float *vtx, float *vty, float *vdx, float *vdy,
             int *npx, int *npy, int *idimp, int *nop, int *nx, int *ny,
             int *ipbc);

void gpush2l_(float *part, float *fxy, float *qbm, float *dt, float *ek,
              int *idimp, int *nop, int *nx, int *ny, int *nxv,
              int *nyv, int *ipbc);

void gpush2lt_(float *part, float *fxy, float *qbm, float *dt, float *ek,
               int *idimp, int *nop, int *nx, int *ny, int *nxv,
               int *nyv, int *ipbc);

void gpost2l_(float *part, float *q, float *qm, int *nop, int *idimp,
              int *nxv, int *nyv);

void gpost2lt_(float *part, float *q, float *qm, int *nop, int *idimp,
               int *nxv, int *nyv);

void dsortp2yl_(float *parta, float *partb, int *npic, int *idimp,
                int *nop, int *ny1);

void cguard2l_(float *fxy, int *nx, int *ny, int *nxe, int *nye);

void aguard2l_(float *q, int *nx, int *ny, int *nxe, int *nye);

void pois22_(float complex *q, float complex *fxy, int *isign,
             float complex *ffc, float *ax, float *ay, float *affp,
             float *we, int *nx, int *ny, int *nxvh, int *nyv,
             int *nxhd, int *nyhd);

void wfft2rinit_(int *mixup, float complex *sct, int *indx, int *indy,
                 int *nxhyd, int *nxyhd);

void fft2rxx_(float complex *f, int *isign, int *mixup, float complex *sct,
              int *indx, int *indy, int *nyi, int *nyp, int *nxhd,
              int *nyd, int *nxhyd, int *nxyhd);

void fft2rxy_(float complex *f, int *isign, int *mixup, float complex *sct,
              int *indx, int *indy, int *nxi, int *nxp, int *nxhd,
              int *nyd, int *nxhyd, int *nxyhd);

void fft2r2x_(float complex *f, int *isign, int *mixup, float complex *sct,
              int *indx, int *indy, int *nyi, int *nyp, int *nxhd,
              int *nyd, int *nxhyd, int *nxyhd);

void fft2r2y_(float complex *f, int *isign, int *mixup, float complex *sct,
              int *indx, int *indy, int *nxi, int *nxp, int *nxhd,
              int *nyd, int *nxhyd, int *nxyhd);
     
void wfft2rx_(float complex *f, int *isign, int *mixup, float complex *sct,
              int *indx, int *indy, int *nxhd, int *nyd, int *nxhyd,
              int *nxyhd);

void wfft2r2_(float complex *f, int *isign, int *mixup, float complex *sct,
              int *indx, int *indy, int *nxhd, int *nyd, int *nxhyd,
              int *nxyhd);

/* Interfaces to C */

double ranorm() {
  return ranorm_();
}

double randum() {
   return randum_();
}

/*--------------------------------------------------------------------*/
void cdistr2(float part[], float vtx, float vty, float vdx, float vdy,
             int npx, int npy, int idimp, int nop, int nx, int ny,
             int ipbc) {
   distr2_(part,&vtx,&vty,&vdx,&vdy,&npx,&npy,&idimp,&nop,&nx,&ny,&ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cgpush2l(float part[], float fxy[], float qbm, float dt, float *ek,
              int idimp, int nop, int nx, int ny, int nxv, int nyv,
              int ipbc) {
   gpush2l_(part,fxy,&qbm,&dt,ek,&idimp,&nop,&nx,&ny,&nxv,&nyv,&ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cgpush2lt(float part[], float fxy[], float qbm, float dt,
               float *ek, int idimp, int nop, int nx, int ny, int nxv, 
               int nyv, int ipbc) {
   gpush2lt_(part,fxy,&qbm,&dt,ek,&idimp,&nop,&nx,&ny,&nxv,&nyv,&ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cgpost2l(float part[], float q[], float qm, int nop, int idimp,
              int nxv, int nyv) {
   gpost2l_(part,q,&qm,&nop,&idimp, &nxv,&nyv);
   return;
}

/*--------------------------------------------------------------------*/
void cgpost2lt(float part[], float q[], float qm, int nop, int idimp,
               int nxv, int nyv) {
   gpost2lt_(part,q,&qm,&nop,&idimp, &nxv,&nyv);
   return;
}

/*--------------------------------------------------------------------*/
void cdsortp2yl(float parta[], float partb[], int npic[], int idimp,
                int nop, int ny1) {
   dsortp2yl_(parta,partb,npic,&idimp,&nop,&ny1);
   return;
}

/*--------------------------------------------------------------------*/
void ccguard2l(float fxy[], int nx, int ny, int nxe, int nye) {
   cguard2l_(fxy,&nx,&ny,&nxe,&nye);
   return;
}

/*--------------------------------------------------------------------*/
void caguard2l(float q[], int nx, int ny, int nxe, int nye) {
   aguard2l_(q,&nx,&ny,&nxe,&nye);
   return;
}

/*--------------------------------------------------------------------*/
void cpois22(float complex q[], float complex fxy[], int isign,
             float complex ffc[], float ax, float ay, float affp,
             float *we, int nx, int ny, int nxvh, int nyv, int nxhd,
             int nyhd) {
   pois22_(q,fxy,&isign,ffc,&ax,&ay,&affp,we,&nx,&ny,&nxvh,&nyv,&nxhd,
           &nyhd);
   return;
}

/*--------------------------------------------------------------------*/
void cwfft2rinit(int mixup[], float complex sct[], int indx, int indy, 
                 int nxhyd, int nxyhd) {
   wfft2rinit_(mixup,sct,&indx,&indy,&nxhyd,&nxyhd);
   return;
}

/*--------------------------------------------------------------------*/
void cfft2rxx(float complex f[], int isign, int mixup[],
              float complex sct[], int indx, int indy, int nyi, int nyp,
              int nxhd, int nyd, int nxhyd, int nxyhd) {
   fft2rxx_(f,&isign,mixup,sct,&indx,&indy,&nyi,&nyp,&nxhd,&nyd,&nxhyd,
            &nxyhd);
   return;
}

/*--------------------------------------------------------------------*/
void cfft2rxy(float complex f[], int isign, int mixup[],
              float complex sct[], int indx, int indy, int nxi, int nxp,
              int nxhd, int nyd, int nxhyd, int nxyhd) {
   fft2rxy_(f,&isign,mixup,sct,&indx,&indy,&nxi,&nxp,&nxhd,&nyd,&nxhyd,
            &nxyhd);
   return;
}

/*--------------------------------------------------------------------*/
void cfft2r2x(float complex f[], int isign, int mixup[],
              float complex sct[], int indx, int indy, int nyi, int nyp,
              int nxhd, int nyd, int nxhyd, int nxyhd) {
   fft2r2x_(f,&isign,mixup,sct,&indx,&indy,&nyi,&nyp,&nxhd,&nyd,&nxhyd,
            &nxyhd);
   return;
}

/*--------------------------------------------------------------------*/
void cfft2r2y(float complex f[], int isign, int mixup[],
              float complex sct[], int indx, int indy, int nxi, int nxp,
              int nxhd, int nyd, int nxhyd, int nxyhd) {
   fft2r2y_(f,&isign,mixup,sct,&indx,&indy,&nxi,&nxp,&nxhd,&nyd,&nxhyd,
            &nxyhd);
   return;
}

/*--------------------------------------------------------------------*/
void cwfft2rx(float complex f[], int isign, int mixup[],
              float complex sct[], int indx, int indy, int nxhd, int nyd,
              int nxhyd, int nxyhd) {
   wfft2rx_(f,&isign,mixup,sct,&indx,&indy,&nxhd,&nyd,&nxhyd,&nxyhd);
   return;
}

/*--------------------------------------------------------------------*/
void cwfft2r2(float complex f[], int isign, int mixup[],
              float complex sct[], int indx, int indy, int nxhd, int nyd,
              int nxhyd, int nxyhd) {
   wfft2r2_(f,&isign,mixup,sct,&indx,&indy,&nxhd,&nyd,&nxhyd,&nxyhd);
   return;
}

