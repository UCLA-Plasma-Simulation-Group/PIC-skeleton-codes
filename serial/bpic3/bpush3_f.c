/* C Library for Skeleton 3D Electromagnetic PIC Code */
/* Wrappers for calling the Fortran routines from a C main program */

#include <complex.h>

double ranorm_();

double randum_();

void distr3_(float *part, float *vtx, float *vty, float *vtz,
             float *vdx, float *vdy, float *vdz, int *npx, int *npy,
             int *npz, int *idimp, int *nop, int *nx, int *ny, int *nz, 
             int *ipbc);

void gbpush3l_(float *part, float *fxyz, float *bxyz, float *qbm,
               float *dt, float *dtc, float *ek, int *idimp, int *nop,
               int *nx, int *ny, int *nz, int *nxv, int *nyv, int *nzv,
               int *ipbc);

void grbpush3l_(float *part, float *fxyz, float *bxyz, float *qbm,
                float *dt, float *dtc, float *ci, float *ek, int *idimp,
                int *nop, int *nx, int *ny, int *nz, int *nxv, int *nyv,
                int *nzv, int *ipbc);

void gpost3l_(float *part, float *q, float *qm, int *nop, int *idimp,
              int *nxv, int *nyv, int *nzv);

void gjpost3l_(float *part, float *cu, float *qm, float *dt, int *nop,
               int *idimp, int *nx, int *ny, int *nz, int *nxv,
               int *nyv, int *nzv, int *ipbc);

void grjpost3l_(float *part, float *cu, float *qm, float *dt, float *ci,
                int *nop, int *idimp, int *nx, int *ny, int *nz,
                int *nxv, int *nyv, int *nzv, int *ipbc);

void dsortp3yzl_(float *parta, float *partb, int *npic, int *idimp,
                 int *nop, int *ny1, int *nyz1);

void cguard3l_(float *fxyz, int *nx, int *ny, int *nz, int *nxe,
               int *nye, int *nze);

void acguard3l_(float *cu, int *nx, int *ny, int *nz, int *nxe,
                int *nye, int *nze);

void aguard3l_(float *q, int *nx, int *ny, int *nz, int *nxe, int *nye,
               int *nze);

void pois33_(float complex *q, float complex *fxyz, int *isign,
             float complex *ffc, float *ax, float *ay, float *az,
             float *affp, float *we, int *nx, int *ny, int *nz,
             int *nxvh, int *nyv, int *nzv, int *nxhd, int *nyhd,
             int *nzhd);

void cuperp3_(float complex *cu, int *nx, int *ny, int *nz, int *nxvh,
              int *nyv, int *nzv);

void ibpois33_(float complex *cu, float complex *bxyz,
               float complex *ffc, float *ci, float *wm, int *nx,
               int *ny, int *nz, int *nxvh, int *nyv, int *nzv,
               int *nxhd, int *nyhd, int *nzhd);

void maxwel3_(float complex *exyz, float complex *bxyz,
              float complex *cu, float complex *ffc, float *ci,
              float *dt, float *wf, float *wm, int *nx, int *ny,
              int *nz, int *nxvh, int *nyv, int *nzv, int *nxhd,
              int *nyhd, int *nzhd);

void emfield3_(float complex *fxyz, float complex *exyz,
               float complex *ffc, int *isign, int *nx, int *ny,
               int *nz, int *nxvh, int *nyv, int *nzv, int *nxhd,
               int *nyhd, int *nzhd);

void wfft3rinit_(int *mixup, float complex *sct, int *indx, int *indy,
                 int *indz, int *nxhyzd, int *nxyzhd);

void fft3rxy_(float complex *f, int *isign, int *mixup,
              float complex *sct, int *indx, int *indy, int *indz,
              int *nzi, int *nzp, int *nxhd, int *nyd, int *nzd,
              int *nxhyzd, int *nxyzhd);

void fft3rxz_(float complex *f, int *isign, int *mixup,
              float complex *sct, int *indx, int *indy, int *indz,
              int *nyi, int *nyp, int *nxhd, int *nyd, int *nzd,
              int *nxhyzd, int *nxyzhd);

void fft3r3xy_(float complex *f, int *isign, int *mixup,
               float complex *sct, int *indx, int *indy, int *indz,
               int *nzi, int *nzp, int *nxhd, int *nyd, int *nzd,
               int *nxhyzd, int *nxyzhd);

void fft3r3z_(float complex *f, int *isign, int *mixup,
              float complex *sct, int *indx, int *indy, int *indz,
              int *nyi, int *nyp, int *nxhd, int *nyd, int *nzd,
              int *nxhyzd,  int *nxyzhd);

void wfft3rx_(float complex *f, int *isign, int *mixup,
              float complex *sct, int *indx, int *indy, int *indz,
              int *nxhd, int *nyd, int *nzd, int *nxhyzd, int *nxyzhd);

void wfft3r3_(float complex *f, int *isign, int *mixup,
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
void cdistr3(float part[], float vtx, float vty, float vtz, float vdx,
             float vdy, float vdz, int npx, int npy, int npz, int idimp,
             int nop, int nx, int ny, int nz, int ipbc) {
   distr3_(part,&vtx,&vty,&vtz,&vdx,&vdy,&vdz,&npx,&npy,&npz,&idimp,
           &nop,&nx,&ny,&nz,&ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cgbpush3l(float part[], float fxyz[], float bxyz[], float qbm,
               float dt, float dtc, float *ek, int idimp, int nop,
               int nx, int ny, int nz, int nxv, int nyv, int nzv,
               int ipbc) {
   gbpush3l_(part,fxyz,bxyz,&qbm,&dt,&dtc,ek,&idimp,&nop,&nx,&ny,&nz,
             &nxv,&nyv,&nzv,&ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cgrbpush3l(float part[], float fxyz[], float bxyz[], float qbm,
                float dt, float dtc, float ci, float *ek, int idimp,
                int nop, int nx, int ny, int nz, int nxv, int nyv,
                int nzv, int ipbc) {
   grbpush3l_(part,fxyz,bxyz,&qbm,&dt,&dtc,&ci,ek,&idimp,&nop,&nx,&ny,
              &nz,&nxv,&nyv,&nzv,&ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cgpost3l(float part[], float q[], float qm, int nop, int idimp,
              int nxv, int nyv, int nzv) {
   gpost3l_(part,q,&qm,&nop,&idimp,&nxv,&nyv,&nzv);
   return;
}

/*--------------------------------------------------------------------*/
void cgjpost3l(float part[], float cu[], float qm, float dt, int nop,
               int idimp, int nx, int ny, int nz, int nxv, int nyv,
               int nzv, int ipbc) {
   gjpost3l_(part,cu,&qm,&dt,&nop,&idimp,&nx,&ny,&nz,&nxv,&nyv,&nzv,
             &ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cgrjpost3l(float part[], float cu[], float qm, float dt, float ci,
                int nop, int idimp, int nx, int ny, int nz, int nxv,
                int nyv, int nzv, int ipbc) {
   grjpost3l_(part,cu,&qm,&dt,&ci,&nop,&idimp,&nx,&ny,&nz,&nxv,&nyv,
              &nzv,&ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cdsortp3yzl(float parta[], float partb[], int npic[], int idimp,
                 int nop, int ny1, int nyz1) {
   dsortp3yzl_(parta,partb,npic,&idimp,&nop,&ny1,&nyz1);
   return;
}

/*--------------------------------------------------------------------*/
void ccguard3l(float fxyz[], int nx, int ny, int nz, int nxe, int nye,
               int nze) {
   cguard3l_(fxyz,&nx,&ny,&nz,&nxe,&nye,&nze);
   return;
}

/*--------------------------------------------------------------------*/
void cacguard3l(float cu[], int nx, int ny, int nz, int nxe, int nye,
                int nze) {
   acguard3l_(cu,&nx,&ny,&nz,&nxe,&nye,&nze);
   return;
}

/*--------------------------------------------------------------------*/
void caguard3l(float q[], int nx, int ny, int nz, int nxe, int nye,
               int nze) {
   aguard3l_(q,&nx,&ny, &nz,&nxe,&nye,&nze);
   return;
}

/*--------------------------------------------------------------------*/
void cpois33(float complex q[], float complex fxyz[], int isign,
             float complex ffc[], float ax, float ay, float az,
             float affp, float *we, int nx, int ny, int nz, int nxvh,
             int nyv, int nzv, int nxhd, int nyhd, int nzhd) {
   pois33_(q,fxyz,&isign,ffc,&ax,&ay,&az,&affp,we,&nx,&ny,&nz,&nxvh,
           &nyv,&nzv,&nxhd,&nyhd,&nzhd);
   return;
}

/*--------------------------------------------------------------------*/
void ccuperp3(float complex cu[], int nx, int ny, int nz, int nxvh,
              int nyv, int nzv) {
   cuperp3_(cu,&nx,&ny,&nz,&nxvh,&nyv,&nzv);
   return;
}

/*--------------------------------------------------------------------*/
void cibpois33(float complex cu[], float complex bxyz[],
               float complex ffc[], float ci, float *wm, int nx, int ny,
               int nz, int nxvh, int nyv, int nzv, int nxhd, int nyhd,
               int nzhd) {
   ibpois33_(cu,bxyz,ffc,&ci,wm,&nx,&ny,&nz,&nxvh,&nyv,&nzv,&nxhd,&nyhd,
             &nzhd);
   return;
}

/*--------------------------------------------------------------------*/
void cmaxwel3(float complex exyz[], float complex bxyz[],
              float complex cu[], float complex ffc[], float ci,
              float dt, float *wf, float *wm, int nx, int ny, int nz,
              int nxvh, int nyv, int nzv, int nxhd, int nyhd, int nzhd) {
   maxwel3_(exyz,bxyz,cu,ffc,&ci,&dt,wf,wm,&nx,&ny,&nz,&nxvh,&nyv,&nzv,
            &nxhd,&nyhd,&nzhd);
   return;
}

/*--------------------------------------------------------------------*/
void cemfield3(float complex fxyz[], float complex exyz[],
               float complex ffc[], int isign, int nx, int ny, int nz,
               int nxvh, int nyv, int nzv, int nxhd, int nyhd,
               int nzhd) {
   emfield3_(fxyz,exyz,ffc,&isign,&nx,&ny,&nz,&nxvh,&nyv,&nzv,&nxhd,
             &nyhd,&nzhd);
   return;
}

/*--------------------------------------------------------------------*/
void cwfft3rinit(int mixup[], float complex sct[], int indx, int indy,
                 int indz, int nxhyzd, int nxyzhd) {
   wfft3rinit_(mixup,sct,&indx,&indy,&indz,&nxhyzd,&nxyzhd);
   return;
}

/*--------------------------------------------------------------------*/
void cfft3rxy(float complex f[], int isign, int mixup[],
              float complex sct[], int indx, int indy, int indz,
              int nzi, int nzp, int nxhd, int nyd, int nzd, int nxhyzd,
              int nxyzhd) {
   fft3rxy_(f,&isign,mixup,sct,&indx,&indy,&indz,&nzi,&nzp,&nxhd,&nyd,
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
void cfft3r3xy(float complex f[], int isign, int mixup[],
               float complex sct[], int indx, int indy, int indz,
               int nzi, int nzp, int nxhd, int nyd, int nzd, int nxhyzd,
               int nxyzhd) {
   fft3r3xy_(f,&isign,mixup,sct,&indx,&indy,&indz,&nzi,&nzp,&nxhd,&nyd,
             &nzd,&nxhyzd,&nxyzhd);
   return;
}

/*--------------------------------------------------------------------*/
void cfft3r3z(float complex f[], int isign, int mixup[],
              float complex sct[], int indx, int indy, int indz,
              int nyi, int nyp, int nxhd, int nyd, int nzd, int nxhyzd, 
              int nxyzhd) {
   fft3r3z_(f,&isign,mixup,sct,&indx,&indy,&indz,&nyi,&nyp,&nxhd,&nyd,
            &nzd,&nxhyzd,&nxyzhd);
   return;
}

/*--------------------------------------------------------------------*/
void cwfft3rx(float complex f[], int isign, int mixup[],
              float complex sct[], int indx, int indy, int indz,
              int nxhd, int nyd, int nzd, int nxhyzd, int nxyzhd) {
   wfft3rx_(f,&isign,mixup,sct,&indx,&indy,&indz,&nxhd,&nyd,&nzd,
            &nxhyzd,&nxyzhd);
   return;
}

/*--------------------------------------------------------------------*/
void cwfft3r3(float complex f[], int isign, int mixup[],
              float complex sct[], int indx, int indy, int indz,
              int nxhd, int nyd, int nzd, int nxhyzd, int nxyzhd) {
   wfft3r3_(f,&isign,mixup,sct,&indx,&indy,&indz,&nxhd,&nyd,&nzd,
            &nxhyzd,&nxyzhd);
   return;
}

