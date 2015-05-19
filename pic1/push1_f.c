/* C Library for Skeleton 1D Electrostatic PIC Code */
/* Wrappers for calling the Fortran routines from a C main program */

#include <complex.h>

double ranorm_();

double randum_();

void distr1_(float *part, float *vtx, float *vdx, int *npx, int *idimp,
             int *nop, int *nx, int *ipbc);

void gpush1l_(float *part, float *fx, float *qbm, float *dt, float *ek,
              int *idimp, int *nop, int *nx, int *nxv, int *ipbc);

void gpost1l_(float *part, float *q, float *qm, int *nop, int *idimp,
              int *nxv);

void dsortp1xl_(float *parta, float *partb, int *npic, int *idimp,
                int *nop, int *nx1);

void cguard1l_(float *fx, int *nx, int *nxe);

void aguard1l_(float *q, int *nx, int *nxe);

void pois1_(float complex *q, float complex *fx, int *isign,
            float complex *ffc, float *ax, float *affp, float *we,
            int *nx);

void wfft1rinit_(int *mixup, float complex *sct, int *indx, int *nxhd);

void fft1rxx_(float complex *f, float complex *t, int *isign,
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
void cdistr1(float part[], float vtx, float vdx, int npx, int idimp,
             int nop, int nx, int ipbc) {
   distr1_(part,&vtx,&vdx,&npx,&idimp,&nop,&nx,&ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cgpush1l(float part[], float fx[], float qbm, float dt, float *ek,
              int idimp, int nop, int nx, int nxv, int ipbc) {
   gpush1l_(part,fx,&qbm,&dt,ek,&idimp,&nop,&nx,&nxv,&ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cgpost1l(float part[], float q[], float qm, int nop, int idimp,
              int nxv) {
   gpost1l_(part,q,&qm,&nop,&idimp,&nxv);
   return;
}

/*--------------------------------------------------------------------*/
void cdsortp1xl(float parta[], float partb[], int npic[], int idimp,
                int nop, int nx1) {
   dsortp1xl_(parta,partb,npic,&idimp,&nop,&nx1);
   return;
}

/*--------------------------------------------------------------------*/
void ccguard1l(float fx[], int nx, int nxe) {
   cguard1l_(fx,&nx,&nxe);
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
