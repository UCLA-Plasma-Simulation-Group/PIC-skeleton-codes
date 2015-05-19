/* C Library for Skeleton 2-1/2D Darwin PIC Code */
/* Wrappers for calling the Fortran routines from a C main program */

#include <complex.h>

double ranorm_();

double randum_();

void distr2h_(float *part, float *vtx, float *vty, float *vtz,
              float *vdx, float *vdy, float *vdz, int *npx, int *npy,
              int *idimp, int *nop, int *nx, int *ny, int *ipbc);

void gbpush23l_(float *part, float *fxy, float *bxy, float *qbm,
                float *dt, float *dtc, float *ek, int *idimp, int *nop,
                int *nx, int *ny, int *nxv, int *nyv, int *ipbc);

void gpost2l_(float *part, float *q, float *qm, int *nop, int *idimp,
              int *nxv, int *nyv);

void gjpost2l_(float *part, float *cu, float *qm, float *dt, int *nop,
               int *idimp, int *nx, int *ny, int *nxv, int *nyv,
               int *ipbc);

void gmjpost2l_(float *part, float *amu, float *qm, int *nop,
                int *idimp, int *nxv, int *nyv);

void gdjpost2l_(float *part, float *fxy, float *bxy, float *dcu,
                float *amu, float *qm, float *qbm, float *dt,
                int *idimp, int *nop, int *nxv, int *nyv);

void gdcjpost2l_(float *part, float *fxy, float *bxy, float *cu,
                 float *dcu, float *amu, float *qm, float *qbm,
                 float *dt, int *idimp, int *nop, int *nxv, int *nyv);

void dsortp2yl_(float *parta, float *partb, int *npic, int *idimp,
                int *nop, int *ny1);

void bguard2l_(float *bxy, int *nx, int *ny, int *nxe, int *nye);

void acguard2l_(float *cu, int *nx, int *ny, int *nxe, int *nye);

void aguard2l_(float *q, int *nx, int *ny, int *nxe, int *nye);

void amcguard2l_(float *amu, int *nx, int *ny, int *nxe, int *nye,
                 int *ndim);

void scfguard2l_(float *cus, float *cu, float *q2m0, int *nx, int *ny,
                 int *nxe, int *nye);

void ascfguard2l_(float *dcu, float *cus, float *q2m0, int *nx, int *ny,
                  int *nxe, int *nye);

void fwpminmx2_(float *qe, float *qbme, float *wpmax, float *wpmin,
                int *nx, int *ny, int *nxe, int *nye);

void pois23_(float complex *q, float complex *fxy, int *isign,
             float complex *ffc, float *ax, float *ay, float *affp,
             float *we, int *nx, int *ny, int *nxvh, int *nyv,
             int *nxhd, int *nyhd);

void cuperp2_(float complex *cu, int *nx, int *ny, int *nxvh, int *nyv);

void bbpois23_(float complex *cu, float complex *bxy,
               float complex *ffc, float *ci, float *wm, int *nx,
               int *ny, int *nxvh, int *nyv, int *nxhd, int *nyhd);

void baddext2_(float *bxy, float *omx, float *omy, float *omz, int *nx,
               int *ny, int *nxe, int *nye);

void dcuperp23_(float complex *dcu, float complex *amu, int *nx,
                int *ny, int *nxvh, int *nyv);

void adcuperp23_(float complex *dcu, float complex *amu, int *nx,
                 int *ny, int *nxvh, int *nyv);

void epois23_(float complex *dcu, float complex *exy, int *isign,
              float complex *ffe, float *ax, float *ay, float *affp,
              float *wp0, float *ci, float *wf, int *nx, int *ny,
              int *nxvh, int *nyv, int *nxhd, int *nyhd);

void addvrfield2_(float *a, float *b, float *c, int *ndim, int *nxe,
                  int *nye);

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

void fft2rnx_(float complex *f, float complex *ss, int *isign,
              int *mixup, float complex *sct, int *indx, int *indy,
              int *nyi, int *nyp, int *nxhd, int *nyd, int *ndim,
              int *nxhyd, int *nxyhd);

void fft2rny_(float complex *f, int *isign, int *mixup,
              float complex *sct, int *indx, int *indy, int *nxi,
              int *nxp, int *nxhd, int *nyd, int *ndim, int *nxhyd,
              int *nxyhd);

void wfft2rx_(float complex *f, int *isign, int *mixup,
              float complex *sct, int *indx, int *indy, int *nxhd,
              int *nyd, int *nxhyd, int *nxyhd);

void wfft2r3_(float complex *f, int *isign, int *mixup,
              float complex *sct, int *indx, int *indy, int *nxhd,
              int *nyd, int *nxhyd, int *nxyhd);

void wfft2rn_(float complex *f, float complex *ss, int *isign,
              int *mixup, float complex *sct, int *indx, int *indy,
              int *nxhd, int *nyd, int *ndim, int *nxhyd, int *nxyhd);

void swapc2n_(float *f, float *s, int *isign, int *nxh, int *nyi,
              int *nyt, int *nxhd, int *nyd, int *ndim);

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
void cgmjpost2l(float part[], float amu[], float qm, int nop, int idimp,
                int nxv, int nyv) {
   gmjpost2l_(part,amu,&qm,&nop,&idimp,&nxv,&nyv);
   return;
}

/*--------------------------------------------------------------------*/
void cgdjpost2l(float part[], float fxy[], float bxy[], float dcu[],
                float amu[], float qm, float qbm, float dt, int idimp,
                int nop, int nxv, int nyv) {
   gdjpost2l_(part,fxy,bxy,dcu,amu,&qm,&qbm,&dt,&idimp,&nop,&nxv,&nyv);
   return;
}

/*--------------------------------------------------------------------*/
void cgdcjpost2l(float part[], float fxy[], float bxy[], float cu[],
                 float dcu[], float amu[], float qm, float qbm,
                 float dt, int idimp, int nop, int nxv, int nyv) {
   gdcjpost2l_(part,fxy,bxy,cu,dcu,amu,&qm,&qbm,&dt,&idimp,&nop,&nxv,
               &nyv);
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
void camcguard2l(float amu[], int nx, int ny, int nxe, int nye,
                 int ndim) {
   amcguard2l_(amu,&nx,&ny,&nxe,&nye,&ndim);
   return;
}

/*--------------------------------------------------------------------*/
void cascfguard2l(float dcu[], float cus[], float q2m0, int nx, int ny,
                  int nxe, int nye) {
   ascfguard2l_(dcu,cus,&q2m0,&nx,&ny,&nxe,&nye);
   return;
}

/*--------------------------------------------------------------------*/
void cfwpminmx2(float qe[], float qbme, float *wpmax, float *wpmin,
                int nx, int ny, int nxe, int nye) {
   fwpminmx2_(qe,&qbme,wpmax,wpmin,&nx,&ny,&nxe,&nye);
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
void cbbpois23(float complex cu[], float complex bxy[],
               float complex ffc[], float ci, float *wm, int nx, int ny,
               int nxvh, int nyv, int nxhd, int nyhd) {
   bbpois23_(cu,bxy,ffc,&ci,wm,&nx,&ny,&nxvh,&nyv,&nxhd,&nyhd);
   return;
}

/*--------------------------------------------------------------------*/
void cbaddext2(float bxy[], float omx, float omy, float omz, int nx,
               int ny, int nxe, int nye) {
   baddext2_(bxy,&omx,&omy,&omz,&nx,&ny,&nxe,&nye);
   return;
}

/*--------------------------------------------------------------------*/
void cdcuperp23(float complex dcu[], float complex amu[], int nx,
                int ny, int nxvh, int nyv) {
   dcuperp23_(dcu,amu,&nx,&ny,&nxvh,&nyv);
   return;
}

/*--------------------------------------------------------------------*/
void cadcuperp23(float complex dcu[], float complex amu[], int nx,
                 int ny, int nxvh, int nyv) {
   adcuperp23_(dcu,amu,&nx,&ny,&nxvh,&nyv);
   return;
}

/*--------------------------------------------------------------------*/
void cepois23(float complex dcu[], float complex exy[], int isign,
              float complex ffe[], float ax, float ay, float affp,
              float wp0, float ci, float *wf, int nx, int ny, int nxvh,
              int nyv, int nxhd, int nyhd) {
   epois23_(dcu,exy,&isign,ffe,&ax,&ay,&affp,&wp0,&ci,wf,&nx,&ny,&nxvh,
            &nyv,&nxhd,&nyhd);
   return;
}

/*--------------------------------------------------------------------*/
void caddvrfield2(float a[], float b[], float c[], int ndim, int nxe,
                  int nye) {
   addvrfield2_(a,b,c,&ndim,&nxe,&nye);
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
void cfft2rnx(float complex f[], float complex ss[], int isign,
              int mixup[], float complex sct[], int indx, int indy,
              int nyi, int nyp, int nxhd, int nyd, int ndim, int nxhyd,
              int nxyhd) {
   fft2rnx_(f,ss,&isign,mixup,sct,&indx,&indy,&nyi,&nyp,&nxhd,&nyd,
            &ndim,&nxhyd,&nxyhd);
   return;
}

/*--------------------------------------------------------------------*/
void cfft2rny(float complex f[], int isign, int mixup[],
              float complex sct[], int indx, int indy, int nxi, int nxp,
              int nxhd, int nyd, int ndim, int nxhyd, int nxyhd) {
   fft2rny_(f,&isign,mixup,sct,&indx,&indy,&nxi,&nxp,&nxhd,&nyd,&ndim,
            &nxhyd,&nxyhd);
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

/*--------------------------------------------------------------------*/
void cwfft2rn(float complex f[], float complex ss[], int isign,
              int mixup[], float complex sct[], int indx, int indy,
              int nxhd, int nyd, int ndim, int nxhyd, int nxyhd) {
   wfft2rn_(f,ss,&isign,mixup,sct,&indx,&indy,&nxhd,&nyd,&ndim,&nxhyd,
            &nxyhd);
   return;
}

/*--------------------------------------------------------------------*/
void cswapc2n(float f[], float s[], int isign, int nxh, int nyi,
              int nyt, int nxhd, int nyd, int ndim) {
   swapc2n_(f,s,&isign,&nxh,&nyi,&nyt,&nxhd,&nyd,&ndim);
   return;
}

