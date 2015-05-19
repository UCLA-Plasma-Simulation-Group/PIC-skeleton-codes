/* C Library for Skeleton 1-2/2D Darwin PIC Code */
/* Wrappers for calling the Fortran routines from a C main program */

#include <complex.h>

double ranorm_();

void distr1h_(float *part, float *vtx, float *vty, float *vtz,
              float *vdx, float *vdy, float *vdz, int *npx, int *idimp, 
              int *nop, int *nx, int *ipbc);

void dblkp1l_(float *part, int *kpic, int *nppmx, int *idimp, int *nop,
              int *mx, int *mx1, int *irc);

void ppmovin1l_(float *part, float *ppart, int *kpic, int *nppmx,
                int *idimp, int *nop, int *mx, int *mx1, int *irc);

void ppcheck1l_(float *ppart, int *kpic, int *idimp, int *nppmx,
                int *nx, int *mx, int *mx1, int *irc);

void gbppush13l_(float *ppart, float *fxyz, float *byz, int *kpic,
                 float *omx, float *qbm, float *dt, float *dtc,
                 float *ek, int *idimp, int *nppmx, int *nx, int *mx,
                 int *nxv, int *mx1, int *ipbc);

void gbppushf13l_(float *ppart, float *fxyz, float *byz, int *kpic,
                  int *ncl, int *ihole, float *omx, float *qbm,
                  float *dt, float *dtc, float *ek, int *idimp,
                  int *nppmx, int *nx, int *mx, int *nxv, int *mx1,
                  int *ntmax, int *irc);

void gppost1l_(float *ppart, float *q, int *kpic, float *qm, int *nppmx,
               int *idimp, int *mx, int *nxv, int *mx1);

void gjppost1l_(float *ppart, float *cu, int *kpic, float *qm,
                float *dt, int *nppmx, int *idimp, int *nx, int *mx,
                int *nxv, int *mx1, int *ipbc);

void gmjppost1l_(float *ppart, float *amu, int *kpic, float *qm,
                 int *nppmx, int *idimp, int *mx, int *nxv, int *mx1);

void gdjppost1l_(float *ppart, float *fxyz, float *byz, float *dcu,
                 float *amu, int *kpic, float *omx, float *qm,
                 float *qbm, float *dt, int *idimp, int *nppmx, int *nx,
                 int *mx, int *nxv, int *mx1);

void gdcjppost1l_(float *ppart, float *fxyz, float *byz, float *cu,
                  float *dcu, float *amu, int *kpic, float *omx,
                  float *qm, float *qbm, float *dt, int *idimp,
                  int *nppmx, int *nx, int *mx, int *nxv, int *mx1);

void pporder1l_(float *ppart, float *ppbuff, int *kpic, int *ncl,
                int *ihole, int *idimp, int *nppmx, int *nx, int *mx,
                int *mx1, int *npbmx, int *ntmax, int *irc);

void pporderf1l_(float *ppart, float *ppbuff, int *kpic, int *ncl,
                 int *ihole, int *idimp, int *nppmx, int *mx1,
                 int *npbmx, int *ntmax, int *irc);

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

/*--------------------------------------------------------------------*/
void cdistr1h(float part[], float vtx, float vty, float vtz, float vdx,
              float vdy, float vdz, int npx, int idimp, int nop, int nx,
              int ipbc) {
   distr1h_(part,&vtx,&vty,&vtz,&vdx,&vdy,&vdz,&npx,&idimp,&nop,&nx,
            &ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cdblkp1l(float part[], int kpic[], int *nppmx, int idimp, int nop,
              int mx, int mx1, int *irc) {
   dblkp1l_(part,kpic,nppmx,&idimp,&nop,&mx,&mx1,irc);
   return;
}

/*--------------------------------------------------------------------*/
void cppmovin1l(float part[], float ppart[], int kpic[], int nppmx,
                int idimp, int nop, int mx, int mx1, int *irc) {
   ppmovin1l_(part,ppart,kpic,&nppmx,&idimp,&nop,&mx,&mx1,irc);
   return;
}

/*--------------------------------------------------------------------*/
void cppcheck1l(float ppart[], int kpic[], int idimp, int nppmx, int nx,
                int mx, int mx1, int *irc) {
   ppcheck1l_(ppart,kpic,&idimp,&nppmx,&nx,&mx,&mx1,irc);
   return;
}

/*--------------------------------------------------------------------*/
void cgbppush13l(float ppart[], float fxyz[], float byz[], int kpic[],
                 float omx, float qbm, float dt, float dtc, float *ek,
                 int idimp, int nppmx, int nx, int mx, int nxv, int mx1,
                 int ipbc) {
   gbppush13l_(ppart,fxyz,byz,kpic,&omx,&qbm,&dt,&dtc,ek,&idimp,&nppmx,
               &nx,&mx,&nxv,&mx1,&ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cgbppushf13l(float ppart[], float fxyz[], float byz[], int kpic[],
                  int ncl[], int ihole[], float omx, float qbm,
                  float dt, float dtc, float *ek, int idimp, int nppmx,
                  int nx, int mx, int nxv, int mx1, int ntmax,
                  int *irc) {
   gbppushf13l_(ppart,fxyz,byz,kpic,ncl,ihole,&omx,&qbm,&dt,&dtc,ek,
                &idimp,&nppmx,&nx,&mx,&nxv,&mx1,&ntmax,irc);
   return;
}

/*--------------------------------------------------------------------*/
void cgppost1l(float ppart[], float q[], int kpic[], float qm,
               int nppmx, int idimp, int mx, int nxv, int mx1) {
   gppost1l_(ppart,q,kpic,&qm,&nppmx,&idimp,&mx,&nxv,&mx1);
   return;
}

/*--------------------------------------------------------------------*/
void cgjppost1l(float ppart[], float cu[], int kpic[], float qm,
                float dt, int nppmx, int idimp, int nx, int mx, int nxv,
                int mx1, int ipbc) {
   gjppost1l_(ppart,cu,kpic,&qm,&dt,&nppmx,&idimp,&nx,&mx,&nxv,&mx1,
              &ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cgmjppost1l(float ppart[], float amu[], int kpic[], float qm,
                 int nppmx, int idimp, int mx, int nxv, int mx1) {
   gmjppost1l_(ppart,amu,kpic,&qm,&nppmx,&idimp,&mx,&nxv,&mx1);
   return;
}

/*--------------------------------------------------------------------*/
void cgdjppost1l(float ppart[], float fxyz[], float byz[], float dcu[],
                 float amu[], int kpic[], float omx, float qm,
                 float qbm, float dt, int idimp, int nppmx, int nx,
                 int mx, int nxv, int mx1) {
   gdjppost1l_(ppart,fxyz,byz,dcu,amu,kpic,&omx,&qm,&qbm,&dt,&idimp,
               &nppmx,&nx,&mx,&nxv,&mx1);
   return;
}

/*--------------------------------------------------------------------*/
void cgdcjppost1l(float ppart[], float fxyz[], float byz[], float cu[],
                  float dcu[], float amu[], int kpic[], float omx,
                  float qm, float qbm, float dt, int idimp, int nppmx,
                  int nx, int mx, int nxv, int mx1) {
   gdcjppost1l_(ppart,fxyz,byz,cu,dcu,amu,kpic,&omx,&qm,&qbm,&dt,&idimp,
                &nppmx,&nx,&mx,&nxv,&mx1);
   return;
}

/*--------------------------------------------------------------------*/
void cpporder1l(float ppart[], float ppbuff[], int kpic[], int ncl[],
                int ihole[], int idimp, int nppmx, int nx, int mx,
                int mx1, int npbmx, int ntmax, int *irc) {
   pporder1l_(ppart,ppbuff,kpic,ncl,ihole,&idimp,&nppmx,&nx,&mx,&mx1,
              &npbmx,&ntmax,irc);
   return;
}

/*--------------------------------------------------------------------*/
void cpporderf1l(float ppart[], float ppbuff[], int kpic[], int ncl[],
                 int ihole[], int idimp, int nppmx, int mx1, int npbmx,
                 int ntmax, int *irc) {
   pporderf1l_(ppart,ppbuff,kpic,ncl,ihole,&idimp,&nppmx,&mx1,&npbmx,
               &ntmax,irc);
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
