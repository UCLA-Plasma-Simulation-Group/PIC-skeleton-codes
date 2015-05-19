/* C Library for Skeleton 1D Electromagnetic PIC Code */
/* Wrappers for calling the Fortran routines from a C main program */

#include <complex.h>

double ranorm_();

double randum_();

void distr1h_(float *part, float *vtx, float *vty, float *vtz,
              float *vdx, float *vdy, float *vdz, int *npx, int *idimp, 
              int *nop, int *nx, int *ipbc);

void gbpush13l_(float *part, float *fxyz, float *byz, float *omx,
                float *qbm, float *dt, float *dtc, float *ek,
                int *idimp, int *nop, int *nx, int *nxv, int *ipbc);

void grbpush13l_(float *part, float *fxyz, float *byz, float *omx,
                 float *qbm, float *dt, float *dtc, float *ci,
                 float *ek, int *idimp, int *nop, int *nx, int *nxv,
                 int *ipbc);

void gpost1l_(float *part, float *q, float *qm, int *nop, int *idimp,
              int *nxv);

void gjpost1l_(float *part, float *cu, float *qm, float *dt, int *nop,
               int *idimp, int *nx, int *nxv, int *ipbc);

void grjpost1l_(float *part, float *cu, float *qm, float *dt, float *ci,
                int *nop, int *idimp, int *nx, int *nxv, int *ipbc);

void dsortp1xl_(float *parta, float *partb, int *npic, int *idimp,
                int *nop, int *nx1);

void cguard1l_(float *byz, int *nx, int *nxe);

void bguard1l_(float *fxyz, int *nx, int *nxe);

void acguard1l_(float *cu, int *nx, int *nxe);

void aguard1l_(float *q, int *nx, int *nxe);

void pois1_(float complex *q, float complex *fx, int *isign,
            float complex *ffc, float *ax, float *affp, float *we,
            int *nx);

void ibpois13_(float complex *cu, float complex *byz,
               float complex *ffc, float *ci, float *wm, int *nx,
               int *nxvh, int *nxhd);

void maxwel1_(float complex *eyz, float complex *byz,
              float complex *cu, float complex *ffc, float *ci,
              float *dt, float *wf, float *wm, int *nx, int *nxvh,
              int *nxhd);

void emfield1_(float complex *fxyz, float complex *fx,
               float complex *eyz, float complex *ffc, int *nx,
               int *nxvh, int *nxhd);

void bmfield1_(float complex *fyz, float complex *eyz,
               float complex *ffc, int *nx, int *nxvh, int *nxhd);

void wfft1rinit_(int *mixup, float complex *sct, int *indx, int *nxhd);

void fft1rxx_(float complex *f, float complex *t, int *isign,
              int *mixup, float complex *sct, int *indx, int *nxd, 
              int *nxhd);

void fft1r2x_(float complex *f, float complex *t, int *isign,
              int *mixup, float complex *sct, int *indx, int *nxd,
              int *nxhd);

void fft1r3x_(float complex *f, float complex *t, int *isign,
              int *mixup, float complex *sct, int *indx, int *nxd,
              int *nxhd);


/* Interfaces to C */

double ranorm() {
  return ranorm_();
}

double randum() {
   return randum_();
}

/*--------------------------------------------------------------------*/
void cdistr1h(float part[], float vtx, float vty, float vtz, float vdx,
              float vdy, float vdz, int npx, int idimp, int nop, int nx,
              int ipbc) {
   distr1h_(part,&vtx,&vty,&vtz,&vdx,&vdy,&vdz,&npx,&idimp,&nop,&nx,
            &ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cgbpush13l(float part[], float fxyz[], float byz[], float omx,
                float qbm, float dt, float dtc, float *ek, int idimp,
                int nop, int nx, int nxv, int ipbc) {
   gbpush13l_(part,fxyz,byz,&omx,&qbm,&dt,&dtc,ek,&idimp,&nop,&nx,&nxv,
              &ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cgrbpush13l(float part[], float fxyz[], float byz[], float omx,
                 float qbm, float dt, float dtc, float ci, float *ek,
                 int idimp, int nop, int nx, int nxv, int ipbc) {
   grbpush13l_(part,fxyz,byz,&omx,&qbm,&dt,&dtc,&ci,ek,&idimp,&nop,&nx,
               &nxv,&ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cgpost1l(float part[], float q[], float qm, int nop, int idimp,
              int nxv) {
   gpost1l_(part,q,&qm,&nop,&idimp,&nxv);
   return;
}

/*--------------------------------------------------------------------*/
void cgjpost1l(float part[], float cu[], float qm, float dt, int nop,
               int idimp, int nx, int nxv, int ipbc) {
   gjpost1l_(part,cu,&qm,&dt,&nop,&idimp,&nx,&nxv,&ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cgrjpost1l(float part[], float cu[], float qm, float dt, float ci,
                int nop, int idimp, int nx, int nxv, int ipbc) {
   grjpost1l_(part,cu,&qm,&dt,&ci,&nop,&idimp,&nx,&nxv,&ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cdsortp1xl(float parta[], float partb[], int npic[], int idimp,
                int nop, int nx1) {
   dsortp1xl_(parta,partb,npic,&idimp,&nop,&nx1);
   return;
}

/*--------------------------------------------------------------------*/
void ccguard1l(float byz[], int nx, int nxe) {
   cguard1l_(byz,&nx,&nxe);
   return;
}

/*--------------------------------------------------------------------*/
void cbguard1l(float fxyz[], int nx, int nxe) {
   bguard1l_(fxyz,&nx,&nxe);
   return;
}

/*--------------------------------------------------------------------*/
void cacguard1l(float cu[], int nx, int nxe) {
   acguard1l_(cu,&nx,&nxe);
   return;
}

/*--------------------------------------------------------------------*/
void caguard1l(float q[], int nx, int nxe) {
   aguard1l_(q,&nx,&nxe);
   return;
}

/*--------------------------------------------------------------------*/
void cpois1(float complex q[], float complex fx[], int isign,
            float complex ffc[], float ax, float affp, float *we,
            int nx) {
   pois1_(q,fx,&isign,ffc,&ax,&affp,we,&nx);
   return;
}

/*--------------------------------------------------------------------*/
void cibpois13(float complex cu[], float complex byz[],
               float complex ffc[], float ci, float *wm, int nx,
               int nxvh, int nxhd) {
   ibpois13_(cu,byz,ffc,&ci,wm,&nx,&nxvh,&nxhd);
   return;
}

/*--------------------------------------------------------------------*/
void cmaxwel1(float complex eyz[], float complex byz[],
              float complex cu[], float complex ffc[], float ci,
              float dt, float *wf, float *wm, int nx, int nxvh,
              int nxhd) {
   maxwel1_(eyz,byz,cu,ffc,&ci,&dt,wf,wm,&nx,&nxvh,&nxhd);
   return;
}

/*--------------------------------------------------------------------*/
void cemfield1(float complex fxyz[], float complex fx[],
               float complex eyz[], float complex ffc[], int nx,
               int nxvh, int nxhd) {
   emfield1_(fxyz,fx,eyz,ffc,&nx,&nxvh,&nxhd);
   return;
}

/*--------------------------------------------------------------------*/
void cbmfield1(float complex fyz[], float complex eyz[],
               float complex ffc[], int nx, int nxvh, int nxhd) {
   bmfield1_(fyz,eyz,ffc,&nx,&nxvh,&nxhd);
   return;
}

/*--------------------------------------------------------------------*/
void cwfft1rinit(int mixup[], float complex sct[], int indx, int nxhd) {
   wfft1rinit_(mixup,sct,&indx,&nxhd);
   return;
}

/*--------------------------------------------------------------------*/
void cfft1rxx(float complex f[], float complex t[], int isign,
              int mixup[], float complex sct[], int indx, int nxd, 
              int nxhd) {
   fft1rxx_(f,t,&isign,mixup,sct,&indx,&nxd,&nxhd);
   return;
}

/*--------------------------------------------------------------------*/
void cfft1r2x(float complex f[], float complex t[], int isign,
              int mixup[], float complex sct[], int indx, int nxd,
              int nxhd) {
   fft1r2x_(f,t,&isign,mixup,sct,&indx,&nxd,&nxhd);
   return;
}

/*--------------------------------------------------------------------*/
void cfft1r3x(float complex f[], float complex t[], int isign,
              int mixup[], float complex sct[], int indx, int nxd,
              int nxhd) {
   fft1r3x_(f,t,&isign,mixup,sct,&indx,&nxd,&nxhd);
   return;
}
