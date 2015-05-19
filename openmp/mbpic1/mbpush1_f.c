/* C Library for Skeleton 1-2/2D Electromagnetic PIC Code */
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

void grbppush13l_(float *ppart, float *fxyz, float *byz, int *kpic,
                  float *omx, float *qbm, float *dt, float *dtc,
                  float *ci, float *ek, int *idimp, int *nppmx, int *nx,
                  int *mx, int *nxv, int *mx1, int *ipbc);

void grbppushf13l_(float *ppart, float *fxyz, float *byz, int *kpic,
                   int *ncl, int *ihole, float *omx, float *qbm,
                   float *dt, float *dtc, float *ci, float *ek,
                   int *idimp, int *nppmx, int *nx, int *mx, int *nxv,
                   int *mx1, int *ntmax, int *irc);
                   
void gppost1l_(float *ppart, float *q, int *kpic, float *qm, int *nppmx,
               int *idimp, int *mx, int *nxv, int *mx1);

void gjppost1l_(float *ppart, float *cu, int *kpic, float *qm,
                float *dt, int *nppmx, int *idimp, int *nx, int *mx,
                int *nxv, int *mx1, int *ipbc);

void gjppostf1l_(float *ppart, float *cu, int *kpic, int *ncl,
                 int *ihole, float *qm, float *dt, int *nppmx,
                 int *idimp, int *nx, int *mx, int *nxv, int *mx1,
                 int *ntmax, int *irc);

void grjppost1l_(float *ppart, float *cu, int *kpic, float *qm,
                 float *dt, float *ci, int *nppmx, int *idimp,
                 int *nx, int *mx, int *nxv, int *mx1, int *ipbc);

void grjppostf1l_(float *ppart, float *cu, int *kpic, int *ncl,
                  int *ihole, float *qm, float *dt, float *ci,
                  int *nppmx, int *idimp, int *nx, int *mx, int *nxv,
                  int *mx1, int *ntmax, int *irc);

void pporder1l_(float *ppart, float *ppbuff, int *kpic, int *ncl,
                int *ihole, int *idimp, int *nppmx, int *nx, int *mx,
                int *mx1, int *npbmx, int *ntmax, int *irc);

void pporderf1l_(float *ppart, float *ppbuff, int *kpic, int *ncl,
                 int *ihole, int *idimp, int *nppmx, int *mx1,
                 int *npbmx, int *ntmax, int *irc);

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
void cgrbppush13l(float ppart[], float fxyz[], float byz[], int kpic[],
                  float omx, float qbm, float dt, float dtc, float ci,
                  float *ek, int idimp, int nppmx, int nx, int mx,
                  int nxv, int mx1, int ipbc) {
   grbppush13l_(ppart,fxyz,byz,kpic,&omx,&qbm,&dt,&dtc,&ci,ek,&idimp,
                &nppmx,&nx,&mx,&nxv,&mx1,&ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cgrbppushf13l(float ppart[], float fxyz[], float byz[], int kpic[],
                   int ncl[], int ihole[], float omx, float qbm,
                   float dt, float dtc, float ci, float *ek, int idimp,
                   int nppmx, int nx, int mx, int nxv, int mx1,
                   int ntmax, int *irc) {
   grbppushf13l_(ppart,fxyz,byz,kpic,ncl,ihole,&omx,&qbm,&dt,&dtc,&ci,
                 ek,&idimp,&nppmx,&nx,&mx,&nxv,&mx1,&ntmax,irc);
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
void cgjppostf1l(float ppart[], float cu[], int kpic[], int ncl[],
                 int ihole[], float qm, float dt, int nppmx, int idimp,
                 int nx, int mx, int nxv, int mx1, int ntmax,
                 int *irc) {
   gjppostf1l_(ppart,cu,kpic,ncl,ihole,&qm,&dt,&nppmx,&idimp,&nx,&mx,
               &nxv,&mx1,&ntmax,irc);
   return;
}

/*--------------------------------------------------------------------*/
void cgrjppost1l(float ppart[], float cu[], int kpic[], float qm,
                 float dt, float ci, int nppmx, int idimp, int nx,
                 int mx, int nxv, int mx1, int ipbc) {
   grjppost1l_(ppart,cu,kpic,&qm,&dt,&ci,&nppmx,&idimp,&nx,&mx,&nxv,
               &mx1,&ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cgrjppostf1l(float ppart[], float cu[], int kpic[], int ncl[],
                  int ihole[], float qm, float dt, float ci, int nppmx,
                  int idimp, int nx, int mx, int nxv, int mx1,
                  int ntmax, int *irc) {
   grjppostf1l_(ppart,cu,kpic,ncl,ihole,&qm,&dt,&ci,&nppmx,&idimp,&nx,
                &mx,&nxv,&mx1,&ntmax,irc);
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
