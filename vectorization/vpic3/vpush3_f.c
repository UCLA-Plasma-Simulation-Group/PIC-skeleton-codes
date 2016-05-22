/* C Library for Skeleton 3D Electrostatic Vector PIC Code */
/* Wrappers for calling the Fortran routines from a C main program */

#include <complex.h>
#include "vpush3.h"

double ranorm_();

double randum_();

void distr3t_(float *part, float *vtx, float *vty, float *vtz, 
              float *vdx, float *vdy, float *vdz, int *npx, int *npy,
              int *npz, int *idimp, int *npe, int *nx, int *ny, int *nz,
              int *ipbc);

void gpush3lt_(float *part, float *fxyz, float *qbm, float *dt,
               float *ek, int *idimp, int *nop, int *npe, int *nx,
               int *ny, int *nz, int *nxv, int *nyv, int *nzv,
               int *ipbc);

void vgpush3lt_(float *part, float *fxyz, float *qbm, float *dt,
                float *ek, int *idimp, int *nop, int *npe, int *nx,
                int *ny, int *nz, int *nxv, int *nyv, int *nzv,
                int *ipbc);

void gpost3lt_(float *part, float *q, float *qm, int *nop, int *npe,
               int *idimp, int *nxv, int *nyv, int *nzv);

void vgpost3lt_(float *part, float *q, float *qm, int *nop, int *npe,
                int *idimp, int *nxv, int *nyv, int *nzv);

void dsortp3yzlt_(float *parta, float *partb, int *npic, int *idimp,
                  int *nop, int *npe, int *ny1, int *nyz1);

void cguard3l_(float *fxyz, int *nx, int *ny, int *nz, int *nxe,
               int *nye, int *nze);

void aguard3l_(float *q, int *nx, int *ny, int *nz, int *nxe, int *nye,
               int *nze);

void vpois33_(float complex *q, float complex *fxyz, int *isign,
              float complex *ffc, float *ax, float *ay, float *az,
              float *affp, float *we, int *nx, int *ny, int *nz,
              int *nxvh, int *nyv, int *nzv, int *nxhd, int *nyhd,
              int *nzhd);

void wfft3rinit_(int *mixup, float complex *sct, int *indx, int *indy,
                 int *indz, int *nxhyzd, int *nxyzhd);

void fft3rvxy_(float complex *f, int *isign, int *mixup,
               float complex *sct, int *indx, int *indy, int *indz,
               int *nzi, int *nzp, int *nxhd, int *nyd, int *nzd,
               int *nxhyzd, int *nxyzhd);

void fft3rxz_(float complex *f, int *isign, int *mixup,
              float complex *sct, int *indx, int *indy, int *indz,
              int *nyi, int *nyp, int *nxhd, int *nyd, int *nzd,
              int *nxhyzd, int *nxyzhd);

void fft3rv3xy_(float complex *f, int *isign, int *mixup,
                float complex *sct, int *indx, int *indy, int *indz,
                int *nzi, int *nzp, int *nxhd, int *nyd, int *nzd,
                int *nxhyzd, int *nxyzhd);

void fft3rv3z_(float complex *f, int *isign, int *mixup,
               float complex *sct, int *indx, int *indy, int *indz,
               int *nyi, int *nyp, int *nxhd, int *nyd, int *nzd,
               int *nxhyzd, int *nxyzhd);

void wfft3rvx_(float complex *f, int *isign, int *mixup,
               float complex *sct, int *indx, int *indy, int *indz,
               int *nxhd, int *nyd, int *nzd, int *nxhyzd, int *nxyzhd);

void wfft3rv3_(float complex *f, int *isign, int *mixup,
               float complex *sct, int *indx, int *indy, int *indz,
               int *nxhd, int *nyd, int *nzd, int *nxhyzd, int *nxyzhd);

/* Interfaces to C */

double ranorm() {
  return ranorm_();
}

double randum() {
   return randum_();
}

/*--------------------------------------------------------------------*/
void cdistr3t(float part[], float vtx, float vty, float vtz, float vdx,
              float vdy, float vdz, int npx, int npy, int npz,
              int idimp, int npe, int nx, int ny, int nz, int ipbc) {
   distr3t_(part,&vtx,&vty,&vtz,&vdx,&vdy,&vdz,&npx,&npy,&npz,&idimp,
            &npe,&nx,&ny,&nz,&ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cgpush3lt(float part[], float fxyz[], float qbm, float dt,
               float *ek, int idimp, int nop, int npe, int nx, int ny,
               int nz, int nxv, int nyv, int nzv, int ipbc) {
   gpush3lt_(part,fxyz,&qbm,&dt,ek,&idimp,&nop,&npe,&nx,&ny,&nz,&nxv,
             &nyv,&nzv,&ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cvgpush3lt(float part[], float fxyz[], float qbm, float dt,
                float *ek, int idimp, int nop, int npe, int nx, int ny,
                int nz, int nxv, int nyv, int nzv, int ipbc) {
   vgpush3lt_(part,fxyz,&qbm,&dt,ek,&idimp,&nop,&npe,&nx,&ny,&nz,&nxv,
              &nyv,&nzv,&ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cv2gpush3lt(float part[], float fxyz[], float qbm, float dt,
                 float *ek, int idimp, int nop, int npe, int nx, int ny,
                int nz, int nxv, int nyv, int nzv, int ipbc) {
   v2gpush3lt_(part,fxyz,&qbm,&dt,ek,&idimp,&nop,&npe,&nx,&ny,&nz,&nxv,
               &nyv,&nzv,&ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cgpost3lt(float part[], float q[], float qm, int nop, int npe,
               int idimp, int nxv, int nyv, int nzv) {
   gpost3lt_(part,q,&qm,&nop,&npe,&idimp,&nxv,&nyv,&nzv);
   return;
}

/*--------------------------------------------------------------------*/
void cvgpost3lt(float part[], float q[], float qm, int nop, int npe,
                int idimp, int nxv, int nyv, int nzv) {
   vgpost3lt_(part,q,&qm,&nop,&npe,&idimp,&nxv,&nyv,&nzv);
   return;
}

/*--------------------------------------------------------------------*/
void cv2gpost3lt(float part[], float q[], float qm, int nop, int npe,
                 int idimp, int nxv, int nyv, int nzv) {
   v2gpost3lt_(part,q,&qm,&nop,&npe,&idimp,&nxv,&nyv,&nzv);
   return;
}

/*--------------------------------------------------------------------*/
void cdsortp3yzlt(float parta[], float partb[], int npic[], int idimp,
                  int nop, int npe, int ny1, int nyz1) {
   dsortp3yzlt_(parta,partb,npic,&idimp,&nop,&npe,&ny1,&nyz1);
   return;
}

/*--------------------------------------------------------------------*/
void ccguard3l(float fxyz[], int nx, int ny, int nz, int nxe, int nye,
               int nze) {
   cguard3l_(fxyz,&nx,&ny,&nz,&nxe,&nye,&nze);
   return;
}

/*--------------------------------------------------------------------*/
void caguard3l(float q[], int nx, int ny, int nz, int nxe, int nye,
               int nze) {
   aguard3l_(q,&nx,&ny,&nz,&nxe,&nye,&nze);
   return;
}

/*--------------------------------------------------------------------*/
void cvpois33(float complex q[], float complex fxyz[], int isign,
              float complex ffc[], float ax, float ay, float az,
              float affp, float *we, int nx, int ny, int nz, int nxvh,
              int nyv, int nzv, int nxhd, int nyhd, int nzhd) {
   vpois33_(q,fxyz,&isign,ffc,&ax,&ay,&az,&affp,we,&nx,&ny,&nz,&nxvh,
            &nyv,&nzv,&nxhd,&nyhd,&nzhd);
   return;
}

/*--------------------------------------------------------------------*/
void cwfft3rinit(int mixup[], float complex sct[], int indx, int indy,
                 int indz, int nxhyzd, int nxyzhd) {
   wfft3rinit_(mixup,sct,&indx,&indy,&indz,&nxhyzd,&nxyzhd);
   return;
}

/*--------------------------------------------------------------------*/
void cfft3rvxy(float complex f[], int isign, int mixup[],
               float complex sct[], int indx, int indy, int indz,
               int nzi, int nzp, int nxhd, int nyd, int nzd, int nxhyzd,
               int nxyzhd) {
   fft3rvxy_(f,&isign,mixup,sct,&indx,&indy,&indz,&nzi,&nzp,&nxhd,&nyd,
             &nzd,&nxhyzd,&nxyzhd);
   return;
}

/*--------------------------------------------------------------------*/
void cfft3rxz(float complex f[], int isign, int mixup[],
              float complex sct[], int indx, int indy, int indz,
              int nyi, int nyp, int nxhd, int nyd, int nzd, int nxhyzd,
              int nxyzhd) {
   fft3rxz_(f,&isign,mixup,sct,&indx,&indy,&indz,&nyi,&nyp,&nxhd,&nyd,
            &nzd,&nxhyzd,&nxyzhd);
   return;
}

/*--------------------------------------------------------------------*/
void cfft3rv3xy(float complex f[], int isign, int mixup[],
                float complex sct[], int indx, int indy, int indz,
                int nzi, int nzp, int nxhd, int nyd, int nzd,
                int nxhyzd, int nxyzhd) {
   fft3rv3xy_(f,&isign,mixup,sct,&indx,&indy,&indz,&nzi,&nzp,&nxhd,&nyd,
              &nzd,&nxhyzd,&nxyzhd);
   return;
}

/*--------------------------------------------------------------------*/
void cfft3rv3z(float complex f[], int isign, int mixup[],
               float complex sct[], int indx, int indy, int indz,
               int nyi, int nyp, int nxhd, int nyd, int nzd, int nxhyzd,
               int nxyzhd) {
   fft3rv3z_(f,&isign,mixup,sct,&indx,&indy,&indz,&nyi,&nyp,&nxhd,&nyd,
             &nzd,&nxhyzd,&nxyzhd);
   return;
}

/*--------------------------------------------------------------------*/
void cwfft3rvx(float complex f[], int isign, int mixup[],
               float complex sct[], int indx, int indy, int indz,
               int nxhd, int nyd, int nzd, int nxhyzd, int nxyzhd) {
   wfft3rvx_(f,&isign,mixup,sct,&indx,&indy,&indz,&nxhd,&nyd,&nzd,
             &nxhyzd,&nxyzhd);
   return;
}

/*--------------------------------------------------------------------*/
void cwfft3rv3(float complex f[], int isign, int mixup[],
               float complex sct[], int indx, int indy, int indz,
               int nxhd, int nyd, int nzd, int nxhyzd, int nxyzhd) {
   wfft3rv3_(f,&isign,mixup,sct,&indx,&indy,&indz,&nxhd,&nyd,&nzd,
             &nxhyzd,&nxyzhd);
   return;
}

