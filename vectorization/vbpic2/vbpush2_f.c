/* C Library for Skeleton 2-1/2D Electromagnetic Vector PIC Code */
/* Wrappers for calling the Fortran routines from a C main program */

#include <complex.h>

double ranorm_();

double randum_();

void distr2ht_(float *part, float *vtx, float *vty, float *vtz,
               float *vdx, float *vdy, float *vdz, int *npx, int *npy,
               int *idimp, int *npe, int *nx, int *ny, int *ipbc);

void gbpush23lt_(float *part, float *fxy, float *bxy, float *qbm,
                 float *dt, float *dtc, float *ek, int *idimp,
                 int *nop, int *npe, int *nx, int *ny, int *nxv,
                 int *nyv, int *ipbc);

void grbpush23lt_(float *part, float *fxy, float *bxy, float *qbm,
                  float *dt, float *dtc, float *ci, float *ek,
                  int *idimp, int *nop, int *npe, int *nx, int *ny,
                  int *nxv, int *nyv, int *ipbc);

void vgbpush23lt_(float *part, float *fxy, float *bxy, float *qbm,
                  float *dt, float *dtc, float *ek, int *idimp, 
                  int *nop, int *npe, int *nx, int *ny, int *nxv,
                  int *nyv, int *ipbc);

void vgrbpush23lt_(float *part, float *fxy, float *bxy, float *qbm,
                   float *dt, float *dtc, float *ci, float *ek,
                   int *idimp, int *nop, int *npe, int *nx, int *ny,
                   int *nxv, int *nyv, int *ipbc);

void gpost2lt_(float *part, float *q, float *qm, int *nop, int *npe,
               int *idimp, int *nxv, int *nyv);

void vgpost2lt_(float *part, float *q, float *qm, int *nop, int *npe,
                int *idimp, int *nxv, int *nyv);

void gjpost2lt_(float *part, float *cu, float *qm, float *dt, int *nop,
                int *npe, int *idimp, int *nx, int *ny, int *nxv,
                int *nyv, int *ipbc);

void grjpost2lt_(float *part, float *cu, float *qm, float *dt,
                 float *ci, int *nop, int *npe, int *idimp, int *nx,
                 int *ny, int *nxv, int *nyv, int *ipbc);

void vgjpost2lt_(float *part, float *cu, float *qm, float *dt,
                 int *nop, int *npe, int *idimp, int *nx, int *ny,
                 int *nxv, int *nyv, int *ipbc);

void vgrjpost2lt_(float *part, float *cu, float *qm, float *dt,
                  float *ci, int *nop, int *npe, int *idimp, int *nx,
                  int *ny, int *nxv, int *nyv, int *ipbc);

void dsortp2ylt_(float *parta, float *partb, int *npic, int *idimp,
                 int *nop, int *npe, int *ny1);

void bguard2l_(float *bxy, int *nx, int *ny, int *nxe, int *nye);

void acguard2l_(float *cu, int *nx, int *ny, int *nxe, int *nye);

void aguard2l_(float *q, int *nx, int *ny, int *nxe, int *nye);

void vpois23_(float complex *q, float complex *fxy, int *isign,
              float complex *ffc, float *ax, float *ay, float *affp,
              float *we, int *nx, int *ny, int *nxvh, int *nyv,
              int *nxhd, int *nyhd);

void cuperp2_(float complex *cu, int *nx, int *ny, int *nxvh, int *nyv);

void vibpois23_(float complex *cu, float complex *bxy,
                float complex *ffc, float *ci, float *wm, int *nx,
                int *ny, int *nxvh, int *nyv, int *nxhd, int *nyhd);

void vmaxwel2_(float complex *exy, float complex *bxy,
               float complex *cu, float complex *ffc, float *ci,
               float *dt, float *wf, float *wm, int *nx, int *ny,
               int *nxvh, int *nyv, int *nxhd, int *nyhd);

void vemfield2_(float complex *fxy, float complex *exy,
                float complex *ffc, int *isign, int *nx, int *ny,
                int *nxvh, int *nyv, int *nxhd, int *nyhd);

void wfft2rinit_(int *mixup, float complex *sct, int *indx, int *indy,
                 int *nxhyd, int *nxyhd);

void wfft2rvx_(float complex *f, int *isign, int *mixup,
               float complex *sct, int *indx, int *indy, int *nxhd,
               int *nyd, int *nxhyd, int *nxyhd);

void wfft2rv3_(float complex *f, int *isign, int *mixup,
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
void cdistr2ht(float part[], float vtx, float vty, float vtz, float vdx,
               float vdy, float vdz, int npx, int npy, int idimp,
               int npe, int nx, int ny, int ipbc) {
   distr2ht_(part,&vtx,&vty,&vtz,&vdx,&vdy,&vdz,&npx,&npy,&idimp,&npe,
             &nx,&ny,&ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cgbpush23lt(float part[], float fxy[], float bxy[], float qbm,
                 float dt, float dtc, float *ek, int idimp, int nop,
                 int npe, int nx, int ny, int nxv, int nyv, int ipbc) {
   gbpush23lt_(part,fxy,bxy,&qbm,&dt,&dtc,ek,&idimp,&nop,&npe,&nx,&ny,
               &nxv,&nyv,&ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cgrbpush23lt(float part[], float fxy[], float bxy[], float qbm,
                  float dt, float dtc, float ci, float *ek, int idimp,
                  int nop, int npe, int nx, int ny, int nxv, int nyv,
                  int ipbc) {
   grbpush23lt_(part,fxy,bxy,&qbm,&dt,&dtc,&ci,ek,&idimp,&nop,&npe,&nx,
                &ny,&nxv,&nyv,&ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cvgbpush23lt(float part[], float fxy[], float bxy[], float qbm,
                  float dt, float dtc, float *ek, int idimp, int nop,
                  int npe, int nx, int ny, int nxv, int nyv, int ipbc) {
   vgbpush23lt_(part,fxy,bxy,&qbm,&dt,&dtc,ek,&idimp,&nop,&npe,&nx,&ny,
                &nxv,&nyv,&ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cvgrbpush23lt(float part[], float fxy[], float bxy[], float qbm,
                   float dt, float dtc, float ci, float *ek, int idimp,
                   int nop, int npe, int nx, int ny, int nxv, int nyv,
                   int ipbc) {
   vgrbpush23lt_(part,fxy,bxy,&qbm,&dt,&dtc,&ci,ek,&idimp,&nop,&npe,
                 &nx,&ny,&nxv,&nyv,&ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cgpost2lt(float part[], float q[], float qm, int nop, int npe,
               int idimp, int nxv, int nyv) {
   gpost2lt_(part,q,&qm,&nop,&npe,&idimp,&nxv,&nyv);
   return;
}

/*--------------------------------------------------------------------*/
void cvgpost2lt(float part[], float q[], float qm, int nop, int npe,
                int idimp, int nxv, int nyv) {
   vgpost2lt_(part,q,&qm,&nop,&npe,&idimp,&nxv,&nyv);
   return;
}

/*--------------------------------------------------------------------*/
void cgjpost2lt(float part[], float cu[], float qm, float dt, int nop,
                int npe, int idimp, int nx, int ny, int nxv, int nyv,
                int ipbc) {
   gjpost2lt_(part,cu,&qm,&dt,&nop,&npe,&idimp,&nx,&ny,&nxv,&nyv,
              &ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cgrjpost2lt(float part[], float cu[], float qm, float dt, float ci,
                int nop, int npe, int idimp, int nx, int ny, int nxv,
                int nyv, int ipbc) {
   grjpost2lt_(part,cu,&qm,&dt,&ci,&nop,&npe,&idimp,&nx,&ny,&nxv,&nyv,
               &ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cvgjpost2lt(float part[], float cu[], float qm, float dt, int nop,
                 int npe, int idimp, int nx, int ny, int nxv, int nyv,
                 int ipbc) {
   vgjpost2lt_(part,cu,&qm,&dt,&nop,&npe,&idimp,&nx,&ny,&nxv,&nyv,
               &ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cvgrjpost2lt(float part[], float cu[], float qm, float dt,
                  float ci, int nop, int npe, int idimp, int nx, int ny,
                  int nxv, int nyv, int ipbc) {
   vgrjpost2lt_(part,cu,&qm,&dt,&ci,&nop,&npe,&idimp,&nx,&ny,&nxv,&nyv,
                &ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cdsortp2ylt(float parta[], float partb[], int npic[], int idimp,
                 int nop, int npe, int ny1) {
   dsortp2ylt_(parta,partb,npic,&idimp,&nop,&npe,&ny1);
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
void cvpois23(float complex q[], float complex fxy[], int isign,
              float complex ffc[], float ax, float ay, float affp,
              float *we, int nx, int ny, int nxvh, int nyv, int nxhd,
              int nyhd) {
   vpois23_(q,fxy,&isign,ffc,&ax,&ay,&affp,we,&nx,&ny,&nxvh,&nyv,&nxhd,
            &nyhd);
   return;
}

/*--------------------------------------------------------------------*/
void ccuperp2(float complex cu[], int nx, int ny, int nxvh, int nyv) {
   cuperp2_(cu,&nx,&ny,&nxvh,&nyv);
   return;
}

/*--------------------------------------------------------------------*/
void cvibpois23(float complex cu[], float complex bxy[],
                float complex ffc[], float ci, float *wm, int nx, int ny,
                int nxvh, int nyv, int nxhd, int nyhd) {
   vibpois23_(cu,bxy,ffc,&ci,wm,&nx,&ny,&nxvh,&nyv,&nxhd,&nyhd);
   return;
}

/*--------------------------------------------------------------------*/
void cvmaxwel2(float complex exy[], float complex bxy[],
               float complex cu[], float complex ffc[], float ci,
               float dt, float *wf, float *wm, int nx, int ny, int nxvh,
               int nyv, int nxhd, int nyhd) {
   vmaxwel2_(exy,bxy,cu,ffc,&ci,&dt,wf,wm,&nx,&ny,&nxvh,&nyv,&nxhd,
             &nyhd);
   return;
}

/*--------------------------------------------------------------------*/
void cvemfield2(float complex fxy[], float complex exy[],
                float complex ffc[], int isign, int nx, int ny, int nxvh,
                int nyv, int nxhd, int nyhd) {
   vemfield2_(fxy,exy,ffc,&isign,&nx,&ny,&nxvh,&nyv,&nxhd,&nyhd);
   return;
}

/*--------------------------------------------------------------------*/
void cwfft2rinit(int mixup[], float complex sct[], int indx, int indy, 
                 int nxhyd, int nxyhd) {
   wfft2rinit_(mixup,sct,&indx,&indy,&nxhyd,&nxyhd);
   return;
}

/*--------------------------------------------------------------------*/
void cwfft2rvx(float complex f[], int isign, int mixup[],
               float complex sct[], int indx, int indy, int nxhd, int nyd,
               int nxhyd, int nxyhd) {
   wfft2rvx_(f,&isign,mixup,sct,&indx,&indy,&nxhd,&nyd,&nxhyd,&nxyhd);
   return;
}

/*--------------------------------------------------------------------*/
void cwfft2rv3(float complex f[],int isign, int mixup[],
               float complex sct[], int indx, int indy, int nxhd,
               int nyd, int nxhyd, int nxyhd) {
   wfft2rv3_(f,&isign,mixup,sct,&indx,&indy,&nxhd,&nyd,&nxhyd,&nxyhd);
   return;
}
