/* C Library for Skeleton 1-2/2D Darwin PIC Code */
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

void gpost1l_(float *part, float *q, float *qm, int *nop, int *idimp,
              int *nxv);

void gjpost1l_(float *part, float *cu, float *qm, float *dt, int *nop,
               int *idimp, int *nx, int *nxv, int *ipbc);

void gmjpost1l_(float *part, float *amu, float *qm, int *nop,
                int *idimp, int *nxv);

void gdjpost1l_(float *part, float *fxyz, float *byz, float *dcu,
                float *amu, float *omx, float *qm, float *qbm,
                float *dt, int *idimp, int *nop, int *nxv);

void gdcjpost1l_(float *part, float *fxyz, float *byz, float *cu,
                 float *dcu, float *amu, float *omx, float *qm,
                 float *qbm, float *dt, int *idimp, int *nop,
                 int *nxv);

void dsortp1xl_(float *parta, float *partb, int *npic, int *idimp,
                int *nop, int *nx1);

void dguard1l_(float *fx, int *nx, int *nxe);

void cguard1l_(float *byz, int *nx, int *nxe);

void acguard1l_(float *cu, int *nx, int *nxe);

void aguard1l_(float *q, int *nx, int *nxe);

void ascfguard1l_(float *dcu, float *cus, float *q2m0, int *nx,
                  int *nxe);

void fwpminmx1_(float *qe, float *qbme, float *wpmax, float *wpmin,
                int *nx, int *nxe);

void pois1_(float complex *q, float complex *fx, int *isign,
            float complex *ffc, float *ax, float *affp, float *we,
            int *nx);

void bbpois13_(float complex *cu, float complex *byz,
               float complex *ffc, float *ci, float *wm, int *nx,
               int *nxvh, int *nxhd);

void baddext1_(float *byz, float *omy, float *omz, int *nx, int *nxe);

void dcuperp13_(float complex *dcu, float complex *amu, int *nx,
                int *nxvh);

void adcuperp13_(float complex *dcu, float complex *amu, int *nx,
                 int *nxvh);

void epois13_(float complex *dcu, float complex *eyz, int *isign,
              float complex *ffe, float *ax, float *affp, float *wp0,
              float *ci, float *wf, int *nx, int *nxvh, int *nxhd);

void addvrfield13_(float *fxyze, float *eyze, float *fxe, int *nxe);

void wfft1rinit_(int *mixup, float complex *sct, int *indx, int *nxhd);

void fft1rxx_(float complex *f, float complex *t, int *isign,
              int *mixup, float complex *sct, int *indx, int *nxd, 
              int *nxhd);

void fft1r2x_(float complex *f, float complex *t, int *isign,
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
void cgmjpost1l(float part[], float amu[], float qm, int nop, int idimp,
                int nxv) {
   gmjpost1l_(part,amu,&qm,&nop,&idimp,&nxv);
   return;
}

/*--------------------------------------------------------------------*/
void cgdjpost1l(float part[], float fxyz[], float byz[], float dcu[],
                float amu[], float omx, float qm, float qbm, float dt,
                int idimp, int nop, int nxv) {
   gdjpost1l_(part,fxyz,byz,dcu,amu,&omx,&qm,&qbm,&dt,&idimp,&nop,&nxv);
   return;
}

/*--------------------------------------------------------------------*/
void cgdcjpost1l(float part[], float fxyz[], float byz[], float cu[],
                 float dcu[], float amu[], float omx, float qm,
                 float qbm, float dt, int idimp, int nop, int nxv) {
   gdcjpost1l_(part,fxyz,byz,cu,dcu,amu,&omx,&qm,&qbm,&dt,&idimp,&nop,
               &nxv);
   return;
}

/*--------------------------------------------------------------------*/
void cdsortp1xl(float parta[], float partb[], int npic[], int idimp,
                int nop, int nx1) {
   dsortp1xl_(parta,partb,npic,&idimp,&nop,&nx1);
   return;
}

void cdguard1l(float fx[], int nx, int nxe){
   dguard1l_(fx,&nx,&nxe);
   return;
}

/*--------------------------------------------------------------------*/
void ccguard1l(float byz[], int nx, int nxe) {
   cguard1l_(byz,&nx,&nxe);
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
void cascfguard1l(float dcu[], float cus[], float q2m0, int nx,
                  int nxe) {
   ascfguard1l_(dcu,cus,&q2m0,&nx,&nxe);
   return;
}

/*--------------------------------------------------------------------*/
void cfwpminmx1(float qe[], float qbme, float *wpmax, float *wpmin,
                int nx, int nxe) {
   fwpminmx1_(qe,&qbme,wpmax,wpmin,&nx,&nxe);
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
void cbbpois13(float complex cu[], float complex byz[],
               float complex ffc[], float ci, float *wm, int nx,
               int nxvh, int nxhd) {
   bbpois13_(cu,byz,ffc,&ci,wm,&nx,&nxvh,&nxhd);
   return;
}

/*--------------------------------------------------------------------*/
void cbaddext1(float byz[], float omy, float omz, int nx, int nxe) {
   baddext1_(byz,&omy,&omz,&nx,&nxe);
   return;
}

/*--------------------------------------------------------------------*/
void cdcuperp13(float complex dcu[], float complex amu[], int nx,
                int nxvh) {
   dcuperp13_(dcu,amu,&nx,&nxvh);
   return;
}

/*--------------------------------------------------------------------*/
void cadcuperp13(float complex dcu[], float complex amu[], int nx,
                 int nxvh) {
   adcuperp13_(dcu,amu,&nx,&nxvh);
   return;
}

/*--------------------------------------------------------------------*/
void cepois13(float complex dcu[], float complex eyz[], int isign,
              float complex ffe[], float ax, float affp, float wp0,
              float ci, float *wf, int nx, int nxvh, int nxhd) {
   epois13_(dcu,eyz,&isign,ffe,&ax,&affp,&wp0,&ci,wf,&nx,&nxvh,&nxhd);
   return;
}

/*--------------------------------------------------------------------*/
void caddvrfield13(float fxyze[], float eyze[], float fxe[], int nxe) {
   addvrfield13_(fxyze,eyze,fxe,&nxe);
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
