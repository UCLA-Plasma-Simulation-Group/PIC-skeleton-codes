/* C Library for Skeleton 2-1/2D Electromagnetic PIC Code */
/* Wrappers for calling the Fortran routines from a C main program */

#include <complex.h>

double ranorm_();

double randum_();

void distr2h_(float *part, float *vtx, float *vty, float *vtz,
              float *vdx, float *vdy, float *vdz, int *npx, int *npy,
              int *idimp, int *nop, int *nx, int *ny, int *ipbc);

void dgbpush23l_(double *part, double *fxy, double *bxy, double *qbm,
                 double *dt, double *dtc, double *ek, int *idimp,
                 int *nop, int *nx, int *ny, int *nxv, int *nyv,
                 int *ipbc);

void gbpush23l_(float *part, float *fxy, float *bxy, float *qbm,
                float *dt, float *dtc, float *ek, int *idimp, int *nop,
                int *nx, int *ny, int *nxv, int *nyv, int *ipbc);

void dgrbpush23l_(double *part, double *fxy, double *bxy, double *qbm,
                  double *dt, double *dtc, double *ci, double *ek,
                  int *idimp, int *nop, int *nx, int *ny, int *nxv,
                  int *nyv, int *ipbc);

void grbpush23l_(float *part, float *fxy, float *bxy, float *qbm,
                 float *dt, float *dtc, float *ci, float *ek,
                 int *idimp, int *nop, int *nx, int *ny, int *nxv,
                 int *nyv, int *ipbc);

void gpost2l_(float *part, float *q, float *qm, int *nop, int *idimp,
              int *nxv, int *nyv);

void gjpost2l_(float *part, float *cu, float *qm, float *dt, int *nop,
               int *idimp, int *nx, int *ny, int *nxv, int *nyv,
               int *ipbc);

void grjpost2l_(float *part, float *cu, float *qm, float *dt, float *ci,
                int *nop, int *idimp, int *nx, int *ny, int *nxv,
                int *nyv, int *ipbc);

void dsortp2yl_(float *parta, float *partb, int *npic, int *idimp,
                int *nop, int *ny1);

void bguard2l_(float *bxy, int *nx, int *ny, int *nxe, int *nye);

void acguard2l_(float *cu, int *nx, int *ny, int *nxe, int *nye);

void aguard2l_(float *q, int *nx, int *ny, int *nxe, int *nye);

void pois23_(float complex *q, float complex *fxy, int *isign,
             float complex *ffc, float *ax, float *ay, float *affp,
             float *we, int *nx, int *ny, int *nxvh, int *nyv,
             int *nxhd, int *nyhd);

void cuperp2_(float complex *cu, int *nx, int *ny, int *nxvh, int *nyv);

void ibpois23_(float complex *cu, float complex *bxy,
               float complex *ffc, float *ci, float *wm, int *nx,
               int *ny, int *nxvh, int *nyv, int *nxhd, int *nyhd);

void maxwel2_(float complex *exy, float complex *bxy,
              float complex *cu, float complex *ffc, float *ci,
              float *dt, float *wf, float *wm, int *nx, int *ny,
              int *nxvh, int *nyv, int *nxhd, int *nyhd);

void emfield2_(float complex *fxy, float complex *exy,
               float complex *ffc, int *isign, int *nx, int *ny,
               int *nxvh, int *nyv, int *nxhd, int *nyhd);

void wfft2rinit_(int *mixup, float complex *sct, int *indx, int *indy,
                 int *nxhyd, int *nxyhd);

void fft2rxx_(float complex *f, int *isign, int *mixup, float complex *sct,
              int *indx, int *indy, int *nyi, int *nyp, int *nxhd,
              int *nyd, int *nxhyd, int *nxyhd);

void fft2rxy_(float complex *f, int *isign, int *mixup, float complex *sct,
              int *indx, int *indy, int *nxi, int *nxp, int *nxhd,
              int *nyd, int *nxhyd, int *nxyhd);

void fft2r3x_(float complex *f, int *isign, int *mixup, float complex *sct,
              int *indx, int *indy, int *nyi, int *nyp, int *nxhd,
              int *nyd, int *nxhyd, int *nxyhd);

void fft2r3y_(float complex *f, int *isign, int *mixup, float complex *sct,
              int *indx, int *indy, int *nxi, int *nxp, int *nxhd,
              int *nyd, int *nxhyd, int *nxyhd);

void wfft2rx_(float complex *f, int *isign, int *mixup,
              float complex *sct, int *indx, int *indy, int *nxhd,
              int *nyd, int *nxhyd, int *nxyhd);

void wfft2r3_(float complex *f, int *isign, int *mixup,
              float complex *sct, int *indx, int *indy, int *nxhd,
              int *nyd, int *nxhyd, int *nxyhd);

/* Interfaces to C */

double ranorm() {
  return ranorm_();
}

double randum() {
   return randum_();
}

/*--------------------------------------------------------------------*/
void cdistr2h(float part[], float vtx, float vty, float vtz, float vdx,
              float vdy, float vdz, int npx, int npy, int idimp, int nop,
              int nx, int ny, int ipbc) {
   distr2h_(part,&vtx,&vty,&vtz,&vdx,&vdy,&vdz,&npx,&npy,&idimp,&nop,
            &nx,&ny,&ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cgbpush23l(float part[], float fxy[], float bxy[], float qbm,
                float dt, float dtc, float *ek, int idimp, int nop,
                int nx, int ny, int nxv, int nyv, int ipbc) {
   gbpush23l_(part,fxy,bxy,&qbm,&dt,&dtc,ek,&idimp,&nop,&nx,&ny,&nxv,
              &nyv,&ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cdgbpush23l(double part[], double fxy[], double bxy[], double qbm,
                 double dt, double dtc, double *ek, int idimp, int nop,
                 int nx, int ny, int nxv, int nyv, int ipbc) {
   dgbpush23l_(part,fxy,bxy,&qbm,&dt,&dtc,ek,&idimp,&nop,&nx,&ny,&nxv,
               &nyv,&ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cgrbpush23l(float part[], float fxy[], float bxy[], float qbm,
                 float dt, float dtc, float ci, float *ek, int idimp,
                 int nop, int nx, int ny, int nxv, int nyv, int ipbc) {
   grbpush23l_(part,fxy,bxy,&qbm,&dt,&dtc,&ci,ek,&idimp,&nop,&nx,&ny,
               &nxv,&nyv,&ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cdgrbpush23l(double part[], double fxy[], double bxy[], double qbm,
                  double dt, double dtc, double ci, double *ek, int idimp,
                  int nop, int nx, int ny, int nxv, int nyv, int ipbc) {
   dgrbpush23l_(part,fxy,bxy,&qbm,&dt,&dtc,&ci,ek,&idimp,&nop,&nx,&ny,
                &nxv,&nyv,&ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cgpost2l(float part[], float q[], float qm, int nop, int idimp,
              int nxv, int nyv) {
   gpost2l_(part,q,&qm,&nop,&idimp,&nxv,&nyv);
   return;
}

/*--------------------------------------------------------------------*/
void cgjpost2l(float part[], float cu[], float qm, float dt, int nop,
               int idimp, int nx, int ny, int nxv, int nyv, int ipbc) {
   gjpost2l_(part,cu,&qm,&dt,&nop,&idimp,&nx,&ny,&nxv,&nyv,&ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cgrjpost2l(float part[], float cu[], float qm, float dt, float ci,
                int nop, int idimp, int nx, int ny, int nxv, int nyv,
                int ipbc) {
   grjpost2l_(part,cu,&qm,&dt,&ci,&nop,&idimp,&nx,&ny,&nxv,&nyv,&ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cdsortp2yl(float parta[], float partb[], int npic[], int idimp,
                int nop, int ny1) {
   dsortp2yl_(parta,partb,npic,&idimp,&nop,&ny1);
   return;
}

/*--------------------------------------------------------------------*/
void cbguard2l(float bxy[], int nx, int ny, int nxe, int nye) {
   bguard2l_(bxy,&nx,&ny,&nxe,&nye);
   return;
}

/*--------------------------------------------------------------------*/
void cacguard2l(float cu[], int nx, int ny, int nxe, int nye) {
   acguard2l_(cu,&nx,&ny,&nxe,&nye);
   return;
}

/*--------------------------------------------------------------------*/
void caguard2l(float q[], int nx, int ny, int nxe, int nye) {
   aguard2l_(q,&nx,&ny,&nxe,&nye);
   return;
}

/*--------------------------------------------------------------------*/
void cpois23(float complex q[], float complex fxy[], int isign,
             float complex ffc[], float ax, float ay, float affp,
             float *we, int nx, int ny, int nxvh, int nyv, int nxhd,
             int nyhd) {
   pois23_(q,fxy,&isign,ffc,&ax,&ay,&affp,we,&nx,&ny,&nxvh,&nyv,&nxhd,
           &nyhd);
   return;
}

/*--------------------------------------------------------------------*/
void ccuperp2(float complex cu[], int nx, int ny, int nxvh, int nyv) {
   cuperp2_(cu,&nx,&ny,&nxvh,&nyv);
   return;
}

/*--------------------------------------------------------------------*/
void cibpois23(float complex cu[], float complex bxy[],
               float complex ffc[], float ci, float *wm, int nx, int ny,
               int nxvh, int nyv, int nxhd, int nyhd) {
   ibpois23_(cu,bxy,ffc,&ci,wm,&nx,&ny,&nxvh,&nyv,&nxhd,&nyhd);
   return;
}

/*--------------------------------------------------------------------*/
void cmaxwel2(float complex exy[], float complex bxy[],
              float complex cu[], float complex ffc[], float ci,
              float dt, float *wf, float *wm, int nx, int ny, int nxvh,
              int nyv, int nxhd, int nyhd) {
   maxwel2_(exy,bxy,cu,ffc,&ci,&dt,wf,wm,&nx,&ny,&nxvh,&nyv,&nxhd,&nyhd);
   return;
}

/*--------------------------------------------------------------------*/
void cemfield2(float complex fxy[], float complex exy[],
               float complex ffc[], int isign, int nx, int ny, int nxvh,
               int nyv, int nxhd, int nyhd) {
   emfield2_(fxy,exy,ffc,&isign,&nx,&ny,&nxvh,&nyv,&nxhd,&nyhd);
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
   fft2rxx_(f,&isign, mixup,sct,&indx,&indy,&nyi,&nyp,&nxhd,&nyd,&nxhyd,
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
void cfft2r3x(float complex f[], int isign, int mixup[],
              float complex sct[], int indx, int indy, int nyi, int nyp,
              int nxhd, int nyd, int nxhyd, int nxyhd) {
   fft2r3x_(f,&isign,mixup,sct,&indx,&indy,&nyi,&nyp,&nxhd,&nyd,&nxhyd,
            &nxyhd);
   return;
}

/*--------------------------------------------------------------------*/
void cfft2r3y(float complex f[], int isign, int mixup[],
              float complex sct[], int indx, int indy, int nxi, int nxp,
              int nxhd, int nyd, int nxhyd, int nxyhd) {
   fft2r3y_(f,&isign, mixup,sct,&indx,&indy,&nxi,&nxp,&nxhd,&nyd,&nxhyd,
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
void cwfft2r3(float complex f[],int isign, int mixup[],
              float complex sct[], int indx, int indy, int nxhd, int nyd,
              int nxhyd, int nxyhd) {
   wfft2r3_(f,&isign,mixup,sct,&indx,&indy,&nxhd,&nyd,&nxhyd,&nxyhd);
   return;
}
