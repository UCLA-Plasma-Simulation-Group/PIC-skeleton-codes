/* C Library for Skeleton 2-1/2D Darwin OpenMP PIC Code */
/* Wrappers for calling the Fortran routines from a C main program */

#include <complex.h>

double ranorm_();

/* Interfaces to Fortran */

void distr2h_(float *part, float *vtx, float *vty, float *vtz,
              float *vdx, float *vdy, float *vdz, int *npx, int *npy,
              int *idimp, int *nop, int *nx, int *ny, int *ipbc);

void dblkp2l_(float *part, int *kpic, int *nppmx, int *idimp, int *nop,
              int *mx, int *my, int *mx1, int *mxy1, int *irc);

void ppmovin2l_(float *part, float *ppart, int *kpic, int *nppmx,
                int *idimp, int *nop, int *mx, int *my, int *mx1,
                int *mxy1, int *irc);

void ppcheck2l_(float *ppart, int *kpic, int *idimp, int *nppmx,
                int *nx, int *ny, int *mx, int *my, int *mx1, int *my1, 
                int *irc);

void gbppush23l_(float *ppart, float *fxy, float *bxy, int *kpic,
                 float *qbm, float *dt, float *dtc, float *ek,
                 int *idimp, int *nppmx, int *nx, int *ny, int *mx,
                 int *my, int *nxv, int *nyv, int *mx1, int *mxy1,
                 int *ipbc);

void gbppushf23l_(float *ppart, float *fxy, float *bxy, int *kpic,
                  int *ncl, int *ihole, float *qbm, float *dt,
                  float *dtc, float *ek, int *idimp, int *nppmx,
                  int *nx, int *ny, int *mx, int *my, int *nxv,
                  int *nyv, int *mx1, int *mxy1, int *ntmax, 
                  int *irc);

void gppost2l_(float *ppart, float *q, int *kpic, float *qm,
               int *nppmx, int *idimp, int *mx, int *my, int *nxv,
               int *nyv, int *mx1, int *mxy1);

void gjppost2l_(float *ppart, float *cu, int *kpic, float *qm,
                float *dt, int *nppmx, int *idimp, int *nx, int *ny, 
                int *mx, int *my, int *nxv, int *nyv, int *mx1,
                int *mxy1, int *ipbc);

void gmjppost2l_(float *ppart, float *amu, int *kpic, float *qm,
                 int *nppmx, int *idimp, int *mx, int *my, int *nxv,
                 int *nyv, int *mx1, int *mxy1);

void gdjppost2l_(float *ppart, float *fxy, float *bxy, float *dcu,
                 float *amu, int *kpic, float *qm, float *qbm,
                 float *dt, int *idimp, int *nppmx, int *nx, int *ny,
                 int *mx, int *my, int *nxv, int *nyv, int *mx1,
                 int *mxy1);

void gdcjppost2l_(float *ppart, float *fxy, float *bxy, float *cu,
                  float *dcu, float *amu, int *kpic, float *qm,
                  float *qbm, float *dt, int *idimp, int *nppmx,
                  int *nx, int *ny, int *mx, int *my, int *nxv,
                  int *nyv, int *mx1, int *mxy1);

void pporder2l_(float *ppart, float *ppbuff, int *kpic, int *ncl,
                int *ihole, int *idimp, int *nppmx, int *nx, int *ny,
                int *mx, int *my, int *mx1, int *my1, int *npbmx,
                int *ntmax, int *irc);

void pporderf2l_(float *ppart, float *ppbuff, int *kpic, int *ncl,
                 int *ihole, int *idimp, int *nppmx, int *mx1,
                 int *my1, int *npbmx, int *ntmax, int *irc);

void bguard2l_(float *bxy, int *nx, int *ny, int *nxe, int *nye);

void acguard2l_(float *cu, int *nx, int *ny, int *nxe, int *nye);

void aguard2l_(float *q, int *nx, int *ny, int *nxe, int *nye);

void amcguard2l_(float *amu, int *nx, int *ny, int *nxe, int *nye,
                 int *ndim);

void ascfguard2l_(float *dcu, float *cus, float *q2m0, int *nx, int *ny,
                  int *nxe, int *nye);

void fwpminmx2_(float *qe, float *qbme, float *wpmax, float *wpmin,
                int *nx, int *ny, int *nxe, int *nye);

void mpois23_(float complex *q, float complex *fxy, int *isign,
              float complex *ffc, float *ax, float *ay, float *affp,
              float *we, int *nx, int *ny, int *nxvh, int *nyv,
              int *nxhd, int *nyhd);

void mcuperp2_(float complex *cu, int *nx, int *ny, int *nxvh,
               int *nyv);

void mbbpois23_(float complex *cu, float complex *bxy,
                float complex *ffc, float *ci, float *wm, int *nx,
                int *ny, int *nxvh, int *nyv, int *nxhd, int *nyhd);

void baddext2_(float *bxy, float *omx, float *omy, float *omz, int *nx,
               int *ny, int *nxe, int *nye);

void mdcuperp23_(float complex *dcu, float complex *amu, int *nx,
                 int *ny, int *nxvh, int *nyv);

void madcuperp23_(float complex *dcu, float complex *amu, int *nx,
                  int *ny, int *nxvh, int *nyv);

void mepois23_(float complex *dcu, float complex *exy, int *isign,
               float complex *ffe, float *ax, float *ay, float *affp,
               float *wp0, float *ci, float *wf, int *nx, int *ny,
               int *nxvh, int *nyv, int *nxhd, int *nyhd);

void addvrfield2_(float *a, float *b, float *c, int *ndim, int *nxe,
                  int *nye);

void wfft2rinit_(int *mixup, float complex *sct, int *indx, int *indy,
                 int *nxhyd, int *nxyhd);

void fft2rmxx_(float complex *f, int *isign, int *mixup,
               float complex *sct, int *indx, int *indy, int *nyi,
               int *nyp, int *nxhd, int *nyd, int *nxhyd, int *nxyhd);

void fft2rmxy_(float complex *f, int *isign, int *mixup,
               float complex *sct, int *indx, int *indy, int *nxi,
               int *nxp, int *nxhd, int *nyd, int *nxhyd, int *nxyhd);

void fft2rm3x_(float complex *f, int *isign, int *mixup,
               float complex *sct, int *indx, int *indy, int *nyi,
               int *nyp, int *nxhd, int *nyd, int *nxhyd, int *nxyhd);

void fft2rm3y_(float complex *f, int *isign, int *mixup,
               float complex *sct, int *indx, int *indy, int *nxi,
               int *nxp, int *nxhd, int *nyd, int *nxhyd, int *nxyhd);

void fft2rmnx_(float complex *f, float complex *ss, int *isign,
               int *mixup, float complex *sct, int *indx, int *indy,
               int *nyi, int *nyp, int *nxhd, int *nyd, int *ndim,
               int *nxhyd, int *nxyhd);

void fft2rmny_(float complex *f, int *isign, int *mixup,
               float complex *sct, int *indx, int *indy, int *nxi,
               int *nxp, int *nxhd, int *nyd, int *ndim, int *nxhyd,
               int *nxyhd);

void wfft2rmx_(float complex *f, int *isign, int *mixup,
               float complex *sct, int *indx, int *indy, int *nxhd,
               int *nyd, int *nxhyd, int *nxyhd);

void wfft2rm3_(float complex *f, int *isign, int *mixup,
               float complex *sct, int *indx, int *indy, int *nxhd,
               int *nyd, int *nxhyd, int *nxyhd);

void wfft2rmn_(float complex *f, float complex *ss, int *isign,
               int *mixup, float complex *sct, int *indx, int *indy,
               int *nxhd, int *nyd, int *ndim, int *nxhyd, int *nxyhd);

void mswapc2n_(float *f, float *s, int *isign, int *nxh, int *nyi,
               int *nyt, int *nxhd, int *nyd, int *ndim);

/* Interfaces to C */

double ranorm() {
  return ranorm_();
}

/*--------------------------------------------------------------------*/
void cdistr2h(float part[], float vtx, float vty, float vtz, float vdx,
              float vdy, float vdz, int npx, int npy, int idimp,
              int nop, int nx, int ny, int ipbc) {
   distr2h_(part,&vtx,&vty, &vtz,&vdx,&vdy,&vdz,&npx,&npy,&idimp,&nop,
            &nx,&ny,&ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cdblkp2l(float part[], int kpic[], int *nppmx, int idimp, int nop,
              int mx, int my, int mx1, int mxy1, int *irc) {
   dblkp2l_(part,kpic,nppmx,&idimp,&nop,&mx,&my,&mx1,&mxy1,irc);
   return;
}

/*--------------------------------------------------------------------*/
void cppmovin2l(float part[], float ppart[], int kpic[], int nppmx,
                int idimp, int nop, int mx, int my, int mx1, int mxy1,
                int *irc) {
   ppmovin2l_(part,ppart,kpic,&nppmx,&idimp,&nop,&mx,&my,&mx1,&mxy1,
              irc);
   return;
}

/*--------------------------------------------------------------------*/
void cppcheck2l(float ppart[], int kpic[], int idimp, int nppmx, int nx,
                int ny, int mx, int my, int mx1, int my1, 
                int *irc) {
   ppcheck2l_(ppart,kpic,&idimp,&nppmx,&nx,&ny,&mx,&my,&mx1,&my1,irc);
   return;
}

/*--------------------------------------------------------------------*/
void cgbppush23l(float ppart[], float fxy[], float bxy[], int kpic[],
                 float qbm, float dt, float dtc, float *ek, int idimp,
                 int nppmx, int nx, int ny, int mx, int my, int nxv,
                 int nyv, int mx1, int mxy1, int ipbc) {
   gbppush23l_(ppart,fxy,bxy,kpic,&qbm,&dt,&dtc,ek,&idimp,&nppmx,&nx,
               &ny,&mx,&my,&nxv,&nyv,&mx1,&mxy1,&ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cgbppushf23l(float ppart[], float fxy[], float bxy[], int kpic[],
                  int ncl[], int ihole[], float qbm, float dt,
                  float dtc, float *ek, int idimp, int nppmx, int nx,
                  int ny, int mx, int my, int nxv, int nyv, int mx1,
                  int mxy1, int ntmax, int *irc) {
   gbppushf23l_(ppart,fxy,bxy,kpic,ncl,ihole,&qbm,&dt,&dtc,ek,&idimp,
                &nppmx,&nx,&ny,&mx,&my,&nxv,&nyv,&mx1,&mxy1,&ntmax,irc);
   return;
}

/*--------------------------------------------------------------------*/
void cgppost2l(float ppart[], float q[], int kpic[], float qm,
               int nppmx, int idimp, int mx, int my, int nxv, int nyv,
               int mx1, int mxy1) {
   gppost2l_(ppart,q,kpic,&qm,&nppmx,&idimp,&mx,&my,&nxv,&nyv,&mx1,
             &mxy1);
   return;
}

/*--------------------------------------------------------------------*/
void cgjppost2l(float ppart[], float cu[], int kpic[], float qm,
                float dt, int nppmx, int idimp, int nx, int ny, int mx,
                int my, int nxv, int nyv, int mx1, int mxy1, int ipbc) {
   gjppost2l_(ppart,cu,kpic,&qm,&dt,&nppmx,&idimp,&nx,&ny,&mx,&my,&nxv,
              &nyv,&mx1,&mxy1,&ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cgmjppost2l(float ppart[], float amu[], int kpic[], float qm,
                 int nppmx, int idimp, int mx, int my, int nxv, int nyv,
                 int mx1, int mxy1) {
   gmjppost2l_(ppart,amu,kpic,&qm,&nppmx,&idimp,&mx,&my,&nxv,&nyv,&mx1,
               &mxy1);
   return;
}

/*--------------------------------------------------------------------*/
void cgdjppost2l(float ppart[], float fxy[], float bxy[], float dcu[],
                 float amu[], int kpic[], float qm, float qbm, float dt,
                 int idimp, int nppmx, int nx, int ny, int mx, int my,
                 int nxv, int nyv, int mx1, int mxy1) {
   gdjppost2l_(ppart,fxy,bxy,dcu,amu,kpic,&qm,&qbm,&dt,&idimp,&nppmx,
               &nx,&ny,&mx,&my,&nxv,&nyv,&mx1,&mxy1);
   return;
}

/*--------------------------------------------------------------------*/
void cgdcjppost2l(float ppart[], float fxy[], float bxy[], float cu[],
                  float dcu[], float amu[], int kpic[], float qm,
                  float qbm, float dt, int idimp, int nppmx, int nx,
                  int ny, int mx, int my, int nxv, int nyv, int mx1,
                  int mxy1) {
   gdcjppost2l_(ppart,fxy,bxy,cu,dcu,amu,kpic,&qm,&qbm,&dt,&idimp,
                &nppmx,&nx,&ny,&mx,&my,&nxv,&nyv,&mx1,&mxy1);
   return;
}

/*--------------------------------------------------------------------*/
void cpporder2l(float ppart[], float ppbuff[], int kpic[], int ncl[],
                int ihole[], int idimp, int nppmx, int nx, int ny,
                int mx, int my, int mx1, int my1, int npbmx, int ntmax,
                int *irc) {
   pporder2l_(ppart,ppbuff,kpic,ncl,ihole,&idimp,&nppmx,&nx,&ny,&mx,&my,
              &mx1,&my1,&npbmx,&ntmax,irc);
   return;
}

/*--------------------------------------------------------------------*/
void cpporderf2l(float ppart[], float ppbuff[], int kpic[], int ncl[],
                 int ihole[], int idimp, int nppmx, int mx1, int my1,
                 int npbmx, int ntmax, int *irc) {
   pporderf2l_(ppart,ppbuff,kpic,ncl,ihole,&idimp,&nppmx,&mx1,&my1,
               &npbmx,&ntmax,irc);
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
void cmpois23(float complex q[], float complex fxy[], int isign,
              float complex ffc[], float ax, float ay, float affp,
              float *we, int nx, int ny, int nxvh, int nyv, int nxhd,
              int nyhd) {
   mpois23_(q,fxy,&isign,ffc,&ax,&ay,&affp,we,&nx,&ny,&nxvh,&nyv,&nxhd,
            &nyhd);
   return;
}

/*--------------------------------------------------------------------*/
void cmcuperp2(float complex cu[], int nx, int ny, int nxvh, int nyv) {
   mcuperp2_(cu,&nx,&ny,&nxvh,&nyv);
   return;
}

/*--------------------------------------------------------------------*/
void cmbbpois23(float complex cu[], float complex bxy[],
                float complex ffc[], float ci, float *wm, int nx,
                int ny, int nxvh, int nyv, int nxhd, int nyhd) {
   mbbpois23_(cu,bxy,ffc,&ci,wm,&nx,&ny,&nxvh,&nyv,&nxhd,&nyhd);
   return;
}

/*--------------------------------------------------------------------*/
void cbaddext2(float bxy[], float omx, float omy, float omz, int nx,
               int ny, int nxe, int nye) {
   baddext2_(bxy,&omx,&omy,&omz,&nx,&ny,&nxe,&nye);
   return;
}

/*--------------------------------------------------------------------*/
void cmdcuperp23(float complex dcu[], float complex amu[], int nx,
                 int ny, int nxvh, int nyv) {
   mdcuperp23_(dcu,amu,&nx,&ny,&nxvh,&nyv);
   return;
}

/*--------------------------------------------------------------------*/
void cmadcuperp23(float complex dcu[], float complex amu[], int nx,
                  int ny, int nxvh, int nyv) {
   madcuperp23_(dcu,amu,&nx,&ny,&nxvh,&nyv);
   return;
}

/*--------------------------------------------------------------------*/
void cmepois23(float complex dcu[], float complex exy[], int isign,
               float complex ffe[], float ax, float ay, float affp,
               float wp0, float ci, float *wf, int nx, int ny, int nxvh,
               int nyv, int nxhd, int nyhd) {
   mepois23_(dcu,exy,&isign,ffe,&ax,&ay,&affp,&wp0,&ci,wf,&nx,&ny,&nxvh,
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
void cfft2rmxx(float complex f[], int isign, int mixup[],
               float complex sct[], int indx, int indy, int nyi,
               int nyp, int nxhd, int nyd, int nxhyd, int nxyhd) {
   fft2rmxx_(f,&isign,mixup,sct,&indx,&indy,&nyi,&nyp,&nxhd,&nyd,&nxhyd,
             &nxyhd);
   return;
}

/*--------------------------------------------------------------------*/
void cfft2rmxy(float complex f[], int isign, int mixup[],
               float complex sct[], int indx, int indy, int nxi,
               int nxp, int nxhd, int nyd, int nxhyd, int nxyhd) {
   fft2rmxy_(f,&isign,mixup,sct,&indx,&indy,&nxi,&nxp,&nxhd,&nyd,&nxhyd,
             &nxyhd);
   return;
}

/*--------------------------------------------------------------------*/
void cfft2rm3x(float complex f[], int isign, int mixup[],
               float complex sct[], int indx, int indy, int nyi,
               int nyp, int nxhd, int nyd, int nxhyd, int nxyhd) {
   fft2rm3x_(f,&isign,mixup,sct,&indx,&indy,&nyi,&nyp,&nxhd,&nyd,&nxhyd,
             &nxyhd);
   return;
}

/*--------------------------------------------------------------------*/
void cfft2rm3y(float complex f[], int isign, int mixup[],
               float complex sct[], int indx, int indy, int nxi,
               int nxp, int nxhd, int nyd, int nxhyd, int nxyhd) {
   fft2rm3y_(f,&isign,mixup,sct,&indx,&indy,&nxi,&nxp,&nxhd,&nyd,&nxhyd,
             &nxyhd);
   return;
}

/*--------------------------------------------------------------------*/
void cfft2rmnx(float complex f[], float complex ss[], int isign,
               int mixup[], float complex sct[], int indx, int indy,
               int nyi, int nyp, int nxhd, int nyd, int ndim, int nxhyd,
               int nxyhd) {
   fft2rmnx_(f,ss,&isign,mixup,sct,&indx,&indy,&nyi,&nyp,&nxhd,&nyd,
             &ndim,&nxhyd,&nxyhd);
   return;
}

/*--------------------------------------------------------------------*/
void cfft2rmny(float complex f[], int isign, int mixup[],
               float complex sct[], int indx, int indy, int nxi,
               int nxp, int nxhd, int nyd, int ndim, int nxhyd,
               int nxyhd) {
/*--------------------------------------------------------------------*/
   fft2rmny_(f,&isign,mixup,sct,&indx,&indy,&nxi,&nxp,&nxhd,&nyd,&ndim,
             &nxhyd,&nxyhd);
   return;
}

/*--------------------------------------------------------------------*/
void cwfft2rmx(float complex f[], int isign, int mixup[],
               float complex sct[], int indx, int indy, int nxhd,
               int nyd, int nxhyd, int nxyhd) {
   wfft2rmx_(f,&isign,mixup,sct,&indx,&indy,&nxhd,&nyd,&nxhyd,&nxyhd);
   return;
}

/*--------------------------------------------------------------------*/
void cwfft2rm3(float complex f[], int isign, int mixup[],
               float complex sct[], int indx, int indy, int nxhd,
               int nyd, int nxhyd, int nxyhd) {
   wfft2rm3_(f,&isign,mixup,sct,&indx,&indy,&nxhd,&nyd,&nxhyd,&nxyhd);
   return;
}

/*--------------------------------------------------------------------*/
void cwfft2rmn(float complex f[], float complex ss[], int isign,
               int mixup[], float complex sct[], int indx, int indy,
               int nxhd, int nyd, int ndim, int nxhyd, int nxyhd) {
   wfft2rmn_(f,ss,&isign,mixup,sct,&indx,&indy,&nxhd,&nyd,&ndim,&nxhyd,
             &nxyhd);
   return;
}

/*--------------------------------------------------------------------*/
void cmswapc2n(float f[], float s[], int isign, int nxh, int nyi,
               int nyt, int nxhd, int nyd, int ndim) {
   mswapc2n_(f,s,&isign,&nxh,&nyi,&nyt,&nxhd,&nyd,&ndim);
   return;
}
