/* C Library for Skeleton 3D Darwin PIC Code */
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

void gpost3l_(float *part, float *q, float *qm, int *nop, int *idimp,
              int *nxv, int *nyv, int *nzv);

void gjpost3l_(float *part, float *cu, float *qm, float *dt, int *nop,
               int *idimp, int *nx, int *ny, int *nz, int *nxv,
               int *nyv, int *nzv, int *ipbc);

void gmjpost3l_(float *part, float *amu, float *qm, int *nop,
                int *idimp, int *nxv, int *nyv, int *nzv);

void gdjpost3l_(float *part, float *fxyz, float *bxyz, float *dcu,
                float *amu, float *qm, float *qbm, float *dt,
                int *idimp, int *nop, int *nxv, int *nyv, int *nzv);

void gdcjpost3l_(float *part, float *fxyz, float *bxyz, float *cu,
                 float *dcu, float *amu, float *qm, float *qbm,
                 float *dt, int *idimp, int *nop, int *nxv, int *nyv,
                 int *nzv);

void dsortp3yzl_(float *parta, float *partb, int *npic, int *idimp,
                 int *nop, int *ny1, int *nyz1);

void cguard3l_(float *fxyz, int *nx, int *ny, int *nz, int *nxe,
               int *nye, int *nze);

void acguard3l_(float *cu, int *nx, int *ny, int *nz, int *nxe,
                int *nye, int *nze);

void aguard3l_(float *q, int *nx, int *ny, int *nz, int *nxe, int *nye,
               int *nze);

void amcguard3l_(float *amu, int *nx, int *ny, int *nz, int *nxe,
                 int *nye, int *nze, int *ndim);

void ascfguard3l_(float *dcu, float *cus, float *q2m0, int *nx,
                  int *ny, int *nz, int *nxe, int *nye, int *nze);

void fwpminmx3_(float *qe, float *qbme, float *wpmax, float *wpmin,
                int *nx, int *ny, int *nz, int *nxe, int *nye, 
                int *nze);

void pois33_(float complex *q, float complex *fxyz, int *isign,
             float complex *ffc, float *ax, float *ay, float *az,
             float *affp, float *we, int *nx, int *ny, int *nz,
             int *nxvh, int *nyv, int *nzv, int *nxhd, int *nyhd,
             int *nzhd);

void cuperp3_(float complex *cu, int *nx, int *ny, int *nz, int *nxvh,
              int *nyv, int *nzv);

void bbpois33_(float complex *cu, float complex *bxyz,
               float complex *ffc, float *ci, float *wm, int *nx,
               int *ny, int *nz, int *nxvh, int *nyv, int *nzv,
               int *nxhd, int *nyhd, int *nzhd);

void baddext3_(float *bxyz, float *omx, float *omy, float *omz, int *nx,
               int *ny, int *nz, int *nxe, int *nye, int *nze);

void dcuperp3_(float complex *dcu, float complex *amu, int *nx, int *ny,
               int *nz, int *nxvh, int *nyv, int *nzv);

void adcuperp3_(float complex *dcu, float complex *amu, int *nx,
                int *ny, int *nz, int *nxvh, int *nyv, int *nzv);

void epois33_(float complex *dcu, float complex *exyz, int *isign,
              float complex *ffe, float *ax, float *ay, float *az,
              float *affp, float *wp0, float *ci, float *wf, int *nx,
              int *ny, int *nz, int *nxvh, int *nyv, int *nzv, int *nxhd,
              int *nyhd, int *nzhd);

void addvrfield3_(float *a, float *b, float *c, int *ndim, int *nxe,
                  int *nye, int *nze);

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

void fft3rnxy_(float complex *f, float complex *ss, int *isign,
               int *mixup, float complex *sct, int *indx, int *indy,
               int *indz, int *nzi, int *nzp, int *nxhd, int *nyd,
               int *nzd, int *ndim, int *nxhyzd, int *nxyzhd);

void fft3rnz_(float complex *f, int *isign, int *mixup,
              float complex *sct, int *indx, int *indy, int *indz,
              int *nyi, int *nyp, int *nxhd, int *nyd, int *nzd,
              int *ndim, int *nxhyzd, int *nxyzhd);

void wfft3rx_(float complex *f, int *isign, int *mixup,
              float complex *sct, int *indx, int *indy, int *indz,
              int *nxhd, int *nyd, int *nzd, int *nxhyzd, int *nxyzhd);

void wfft3r3_(float complex *f, int *isign, int *mixup,
              float complex *sct, int *indx, int *indy, int *indz,
              int *nxhd, int *nyd, int *nzd, int *nxhyzd, int *nxyzhd);

void wfft3rn_(float complex *f, float complex *ss, int *isign,
              int *mixup, float complex *sct, int *indx, int *indy,
              int *indz, int *nxhd, int *nyd, int *nzd, int *ndim, 
              int *nxhyzd, int *nxyzhd);

void swap3cn_(float *f, float *s, int *isign, int *nxh, int *ny,
              int *nzi, int *nzt, int *nxhd, int *nyd, int *nzd, 
              int *ndim);

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
void cgmjpost3l(float part[], float amu[], float qm, int nop, int idimp,
                int nxv, int nyv, int nzv) {
   gmjpost3l_(part,amu,&qm,&nop,&idimp,&nxv,&nyv,&nzv);
   return;
}

/*--------------------------------------------------------------------*/
void cgdjpost3l(float part[], float fxyz[], float bxyz[], float dcu[],
                float amu[], float qm, float qbm, float dt, int idimp,
                int nop, int nxv, int nyv, int nzv) {
   gdjpost3l_(part,fxyz,bxyz,dcu,amu,&qm,&qbm,&dt,&idimp,&nop,&nxv,&nyv,
              &nzv);
   return;
}

/*--------------------------------------------------------------------*/
void cgdcjpost3l(float part[], float fxyz[], float bxyz[], float cu[],
                 float dcu[], float amu[], float qm, float qbm,
                 float dt, int idimp, int nop, int nxv, int nyv,
                 int nzv) {
   gdcjpost3l_(part,fxyz,bxyz,cu,dcu,amu,&qm,&qbm,&dt,&idimp,&nop,&nxv,
               &nyv,&nzv);
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
void camcguard3l(float amu[], int nx, int ny, int nz, int nxe, int nye,
                 int nze, int ndim) {
   amcguard3l_(amu,&nx,&ny,&nz,&nxe,&nye,&nze,&ndim);
   return;
}

/*--------------------------------------------------------------------*/
void cascfguard3l(float dcu[], float cus[], float q2m0, int nx, int ny,
                  int nz, int nxe, int nye, int nze) {
   ascfguard3l_(dcu,cus,&q2m0,&nx,&ny,&nz,&nxe,&nye,&nze);
   return;
}

/*--------------------------------------------------------------------*/
void cfwpminmx3(float qe[], float qbme, float *wpmax, float *wpmin,
                int nx, int ny, int nz, int nxe, int nye, int nze) {
   fwpminmx3_(qe,&qbme,wpmax,wpmin,&nx,&ny,&nz,&nxe,&nye,&nze);
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
void cbbpois33(float complex cu[], float complex bxyz[],
               float complex ffc[], float ci, float *wm, int nx, int ny,
               int nz, int nxvh, int nyv, int nzv, int nxhd, int nyhd,
               int nzhd) {
   bbpois33_(cu,bxyz,ffc,&ci,wm,&nx,&ny,&nz,&nxvh,&nyv,&nzv,&nxhd,&nyhd,
             &nzhd);
   return;
}

/*--------------------------------------------------------------------*/
void cbaddext3(float bxyz[], float omx, float omy, float omz, int nx,
               int ny, int nz, int nxe, int nye, int nze) {
   baddext3_(bxyz,&omx,&omy,&omz,&nx,&ny,&nz,&nxe,&nye,&nze);
   return;
}

/*--------------------------------------------------------------------*/
void cdcuperp3(float complex dcu[], float complex amu[], int nx, int ny,
               int nz, int nxvh, int nyv, int nzv) {
   dcuperp3_(dcu,amu,&nx,&ny,&nz,&nxvh,&nyv,&nzv);
   return;
}

/*--------------------------------------------------------------------*/
void cadcuperp3(float complex dcu[], float complex amu[], int nx,
                int ny, int nz, int nxvh, int nyv, int nzv) {
   adcuperp3_(dcu,amu,&nx,&ny,&nz,&nxvh,&nyv,&nzv);
   return;
}

/*--------------------------------------------------------------------*/
void cepois33(float complex dcu[], float complex exyz[], int isign,
              float complex ffe[], float ax, float ay, float az,
              float affp, float wp0, float ci, float *wf, int nx,
              int ny, int nz, int nxvh, int nyv, int nzv, int nxhd,
              int nyhd, int nzhd) {
   epois33_(dcu,exyz,&isign,ffe,&ax,&ay,&az,&affp,&wp0,&ci,wf,&nx,&ny,
            &nz,&nxvh,&nyv,&nzv,&nxhd,&nyhd,&nzhd);
   return;
}

/*--------------------------------------------------------------------*/
void caddvrfield3(float a[], float b[], float c[], int ndim, int nxe,
                  int nye, int nze) {
   addvrfield3_(a,b,c,&ndim,&nxe,&nye,&nze);
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
void cfft3rnxy(float complex f[], float complex ss[], int isign,
               int mixup[], float complex sct[], int indx, int indy,
               int indz, int nzi, int nzp, int nxhd, int nyd, int nzd,
               int ndim, int nxhyzd, int nxyzhd) {
   fft3rnxy_(f,ss,&isign,mixup,sct,&indx,&indy,&indz,&nzi,&nzp,&nxhd,
             &nyd,&nzd,&ndim,&nxhyzd,&nxyzhd);
   return;
}

/*--------------------------------------------------------------------*/
void cfft3rnz(float complex f[], int isign, int mixup[],
              float complex sct[], int indx, int indy, int indz,
              int nyi, int nyp, int nxhd, int nyd, int nzd, int ndim,
              int nxhyzd, int nxyzhd) {
   fft3rnz_(f,&isign,mixup,sct,&indx,&indy,&indz,&nyi,&nyp,&nxhd,&nyd,
            &nzd,&ndim,&nxhyzd,&nxyzhd);
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

/*--------------------------------------------------------------------*/
void cwfft3rn(float complex f[], float complex ss[], int isign,
              int mixup[], float complex sct[], int indx, int indy,
              int indz, int nxhd, int nyd, int nzd, int ndim, 
              int nxhyzd, int nxyzhd) {
/*--------------------------------------------------------------------*/
   wfft3rn_(f,ss,&isign,mixup,sct,&indx,&indy,&indz,&nxhd,&nyd,&nzd,
            &ndim,&nxhyzd,&nxyzhd);
   return;
}

/*--------------------------------------------------------------------*/
void cswap3cn(float f[], float s[], int isign, int nxh, int ny, int nzi,
              int nzt, int nxhd, int nyd, int nzd, int ndim) {
   swap3cn_(f,s,&isign,&nxh,&ny,&nzi,&nzt,&nxhd,&nyd,&nzd,&ndim);
   return;
}
