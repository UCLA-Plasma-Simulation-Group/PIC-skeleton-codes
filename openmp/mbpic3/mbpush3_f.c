/* C Library for Skeleton 3D Electromagnetic OpenMP PIC Code */
/* Wrappers for calling the Fortran routines from a C main program */

#include <complex.h>

double ranorm_();

/* Interfaces to Fortran */

void distr3_(float *part, float *vtx, float *vty, float *vtz,
             float *vdx, float *vdy, float *vdz, int *npx, int *npy,
             int *npz, int *idimp, int *nop, int *nx, int *ny, int *nz,
             int *ipbc);

void dblkp3l_(float *part, int *kpic, int *nppmx, int *idimp, int *nop,
              int *mx, int *my, int *mz, int *mx1, int *my1, int *mxyz1,
              int *irc);

void ppmovin3l_(float *part, float *ppart, int *kpic, int *nppmx,
                int *idimp, int *nop, int *mx, int *my, int *mz,
                int *mx1, int *my1, int *mxyz1, int *irc);

void ppcheck3l_(float *ppart, int *kpic, int *idimp, int *nppmx,
                int *nx, int *ny, int *nz, int *mx, int *my, int *mz,
                int *mx1, int *my1, int *mz1, int *irc);

void gbppush3l_(float *ppart, float *fxyz, float *bxyz, int *kpic,
                float *qbm, float *dt, float *dtc, float *ek, int *idimp,
                int *nppmx, int *nx, int *ny, int *nz, int *mx, int *my,
                int *mz, int *nxv, int *nyv, int *nzv, int *mx1, int *my1,
                int *mxyz1, int *ipbc);

void gbppushf3l_(float *ppart, float *fxyz, float *bxyz, int *kpic,
                 int *ncl, int *ihole, float *qbm, float *dt,
                 float *dtc, float *ek, int *idimp, int *nppmx, int *nx,
                 int *ny, int *nz, int *mx, int *my, int *mz, int *nxv,
                 int *nyv, int *nzv, int *mx1, int *my1, int *mxyz1,
                 int *ntmax,
                 int *irc);

void grbppush3l_(float *ppart, float *fxyz, float *bxyz, int *kpic,
                 float *qbm, float *dt, float *dtc, float *ci,
                 float *ek, int *idimp, int *nppmx, int *nx, int *ny,
                 int *nz, int *mx, int *my, int *mz, int *nxv, int *nyv,
                 int *nzv, int *mx1, int *my1, int *mxyz1, int *ipbc);

void grbppushf3l_(float *ppart, float *fxyz, float *bxyz, int *kpic,
                  int *ncl, int *ihole, float *qbm, float *dt,
                  float *dtc, float *ci, float *ek, int *idimp,
                  int *nppmx, int *nx, int *ny, int *nz, int *mx,
                  int *my, int *mz, int *nxv, int *nyv, int *nzv,
                  int *mx1, int *my1, int *mxyz1, int *ntmax, int *irc);

void gppost3l_(float *ppart, float *q, int *kpic, float *qm, int *nppmx,
               int *idimp, int *mx, int *my, int *mz, int *nxv,
               int *nyv, int *nzv, int *mx1, int *my1, int *mxyz1);

void gjppost3l_(float *ppart, float *cu, int *kpic, float *qm,
                float *dt, int *nppmx, int *idimp, int *nx, int *ny,
                int *nz, int *mx, int *my, int *mz, int *nxv, int *nyv,
                int *nzv, int *mx1, int *my1, int *mxyz1, int *ipbc);

void gjppostf3l_(float *ppart, float *cu, int *kpic, int *ncl,
                 int *ihole, float *qm, float *dt, int *nppmx,
                 int *idimp, int *nx, int *ny, int *nz, int *mx,
                 int *my, int *mz, int *nxv, int *nyv, int *nzv,
                 int *mx1, int *my1, int *mxyz1, int *ntmax, int *irc);

void grjppost3l_(float *ppart, float *cu, int *kpic, float *qm,
                 float *dt, float *ci, int *nppmx, int *idimp, int *nx,
                 int *ny, int *nz, int *mx, int *my, int *mz, int *nxv,
                 int *nyv, int *nzv, int *mx1, int *my1, int *mxyz1,
                 int *ipbc);

void grjppostf3l_(float *ppart, float *cu, int *kpic, int *ncl,
                  int *ihole, float *qm, float *dt, float *ci,
                  int *nppmx, int *idimp, int *nx, int *ny, int *nz,
                  int *mx, int *my, int *mz, int *nxv, int *nyv,
                  int *nzv, int *mx1, int *my1, int *mxyz1, int *ntmax,
                  int *irc);

void pporder3l_(float *ppart, float *ppbuff, int *kpic, int *ncl,
                int *ihole, int *idimp, int *nppmx, int *nx, int *ny,
                int *nz, int *mx, int *my, int *mz, int *mx1, int *my1,
                int *mz1, int *npbmx, int *ntmax, int *irc);

void pporderf3l_(float *ppart, float *ppbuff, int *kpic, int *ncl,
                 int *ihole, int *idimp, int *nppmx, int *mx1, int *my1,
                 int *mz1, int *npbmx, int *ntmax, int *irc);

void cguard3l_(float *fxyz, int *nx, int *ny, int *nz, int *nxe,
               int *nye, int *nze);

void acguard3l_(float *cu, int *nx, int *ny, int *nz, int *nxe,
                int *nye, int *nze);

void aguard3l_(float *q, int *nx, int *ny, int *nz, int *nxe, int *nye,
               int *nze);

void mpois33_(float complex *q, float complex *fxyz, int *isign,
              float complex *ffc, float *ax, float *ay, float *az,
              float *affp, float *we, int *nx, int *ny, int *nz,
              int *nxvh, int *nyv, int *nzv, int *nxhd, int *nyhd,
              int *nzhd);

void mcuperp3_(float complex *cu, int *nx, int *ny, int *nz, int *nxvh,
               int *nyv, int *nzv);

void mibpois33_(float complex *cu, float complex *bxyz,
                float complex *ffc, float *ci, float *wm, int *nx,
                int *ny, int *nz, int *nxvh, int *nyv, int *nzv,
                int *nxhd, int *nyhd, int *nzhd);

void mmaxwel3_(float complex *exyz, float complex *bxyz,
               float complex *cu, float complex *ffc, float *ci,
               float *dt, float *wf, float *wm, int *nx, int *ny,
               int *nz, int *nxvh, int *nyv, int *nzv, int *nxhd,
               int *nyhd, int *nzhd);

void memfield3_(float complex *fxyz, float complex *exyz,
                float complex *ffc, int *isign, int *nx, int *ny,
                int *nz, int *nxvh, int *nyv, int *nzv, int *nxhd, 
                int *nyhd, int *nzhd);

void wfft3rinit_(int *mixup, float complex *sct, int *indx, int *indy,
                 int *indz, int *nxhyzd, int *nxyzhd);

void fft3rmxy_(float complex *f, int *isign, int *mixup,
               float complex *sct, int *indx, int *indy, int *indz,
               int *nzi, int *nzp, int *nxhd, int *nyd, int *nzd,
               int *nxhyzd, int *nxyzhd);

void fft3rmxz_(float complex *f, int *isign, int *mixup,
               float complex *sct, int *indx, int *indy, int *indz,
               int *nyi, int *nyp, int *nxhd, int *nyd, int *nzd,
               int *nxhyzd, int *nxyzhd);

void fft3rm3xy_(float complex *f, int *isign, int *mixup,
                float complex *sct, int *indx, int *indy, int *indz,
                int *nzi, int *nzp, int *nxhd, int *nyd, int *nzd,
                int *nxhyzd, int *nxyzhd);

void fft3rm3z_(float complex *f, int *isign, int *mixup,
               float complex *sct, int *indx, int *indy, int *indz,
               int *nyi, int *nyp, int *nxhd, int *nyd, int *nzd,
               int *nxhyzd, int *nxyzhd);

void wfft3rmx_(float complex *f, int *isign, int *mixup,
               float complex *sct, int *indx, int *indy, int *indz,
               int *nxhd, int *nyd, int *nzd, int *nxhyzd,
               int *nxyzhd);

void wfft3rm3_(float complex *f, int *isign, int *mixup,
               float complex *sct, int *indx, int *indy, int *indz,
               int *nxhd, int *nyd, int *nzd, int *nxhyzd,
               int *nxyzhd);

/* Interfaces to C */

double ranorm() {
  return ranorm_();
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
void cdblkp3l(float part[], int kpic[], int *nppmx, int idimp, int nop,
              int mx, int my, int mz, int mx1, int my1, int mxyz1,
              int *irc) {
   dblkp3l_(part,kpic,nppmx,&idimp,&nop,&mx,&my,&mz,&mx1,&my1,&mxyz1,
            irc);
   return;
}

/*--------------------------------------------------------------------*/
void cppmovin3l(float part[], float ppart[], int kpic[], int nppmx,
                int idimp, int nop, int mx, int my, int mz, int mx1,
                int my1, int mxyz1, int *irc) {
   ppmovin3l_(part,ppart,kpic,&nppmx,&idimp,&nop,&mx,&my,&mz,&mx1,&my1,
              &mxyz1,irc);
   return;
}

/*--------------------------------------------------------------------*/
void cppcheck3l(float ppart[], int kpic[], int idimp, int nppmx, int nx,
                int ny, int nz, int mx, int my, int mz, int mx1,
                int my1, int mz1, int *irc) {
   ppcheck3l_(ppart,kpic,&idimp,&nppmx,&nx,&ny,&nz,&mx,&my,&mz,&mx1,
              &my1,&mz1,irc);
   return;
}

/*--------------------------------------------------------------------*/
void cgbppush3l(float ppart[], float fxyz[], float bxyz[], int kpic[],
                float qbm, float dt, float dtc, float *ek, int idimp,
                int nppmx, int nx, int ny, int nz, int mx, int my,
                int mz, int nxv, int nyv, int nzv, int mx1, int my1,
                int mxyz1, int ipbc) {
   gbppush3l_(ppart,fxyz,bxyz, kpic,&qbm,&dt,&dtc,ek,&idimp,&nppmx,&nx,
              &ny,&nz,&mx,&my,&mz,&nxv,&nyv,&nzv,&mx1,&my1,&mxyz1,
              &ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cgbppushf3l(float ppart[], float fxyz[], float bxyz[], int kpic[],
                 int ncl[], int ihole[], float qbm, float dt, float dtc,
                 float *ek, int idimp, int nppmx, int nx, int ny,
                 int nz, int mx, int my, int mz, int nxv, int nyv,
                 int nzv, int mx1, int my1, int mxyz1, int ntmax,
                 int *irc) {
   gbppushf3l_(ppart,fxyz,bxyz,kpic,ncl,ihole,&qbm,&dt,&dtc,ek,&idimp,
               &nppmx,&nx,&ny,&nz,&mx,&my,&mz,&nxv,&nyv,&nzv,&mx1,&my1,
               &mxyz1,&ntmax,irc);
   return;
}

/*--------------------------------------------------------------------*/
void cgrbppush3l(float ppart[], float fxyz[], float bxyz[], int kpic[],
                 float qbm, float dt, float dtc, float ci, float *ek,
                 int idimp, int nppmx, int nx, int ny, int nz, int mx,
                 int my, int mz, int nxv, int nyv, int nzv, int mx1,
                 int my1, int mxyz1, int ipbc) {
   grbppush3l_(ppart,fxyz,bxyz,kpic,&qbm,&dt,&dtc,&ci,ek,&idimp,&nppmx,
               &nx,&ny,&nz,&mx,&my,&mz,&nxv,&nyv,&nzv,&mx1,&my1,&mxyz1,
               &ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cgrbppushf3l(float ppart[], float fxyz[], float bxyz[], int kpic[],
                  int ncl[], int ihole[], float qbm, float dt,
                  float dtc, float ci, float *ek, int idimp, int nppmx,
                  int nx, int ny, int nz, int mx, int my, int mz,
                  int nxv, int nyv, int nzv, int mx1, int my1,
                  int mxyz1, int ntmax, int *irc) {
   grbppushf3l_(ppart,fxyz,bxyz,kpic,ncl,ihole,&qbm,&dt,&dtc,&ci,ek,
                &idimp,&nppmx,&nx,&ny,&nz,&mx,&my,&mz,&nxv,&nyv,&nzv,
                &mx1,&my1,&mxyz1,&ntmax,irc);
   return;
}

/*--------------------------------------------------------------------*/
void cgppost3l(float ppart[], float q[], int kpic[], float qm,
               int nppmx, int idimp, int mx, int my, int mz, int nxv,
               int nyv, int nzv, int mx1, int my1, int mxyz1) {
   gppost3l_(ppart,q,kpic,&qm,&nppmx,&idimp,&mx,&my,&mz,&nxv,&nyv,&nzv,
             &mx1,&my1,&mxyz1);
   return;
}

/*--------------------------------------------------------------------*/
void cgjppost3l(float ppart[], float cu[], int kpic[], float qm,
                float dt, int nppmx, int idimp, int nx, int ny, int nz,
                int mx, int my, int mz, int nxv, int nyv, int nzv,
                int mx1, int my1, int mxyz1, int ipbc) {
   gjppost3l_(ppart,cu,kpic,&qm,&dt,&nppmx,&idimp,&nx,&ny,&nz,&mx,&my,
              &mz,&nxv,&nyv,&nzv,&mx1,&my1,&mxyz1,&ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cgjppostf3l(float ppart[], float cu[], int kpic[], int ncl[],
                 int ihole[], float qm, float dt, int nppmx, int idimp,
                 int nx, int ny, int nz, int mx, int my, int mz,
                 int nxv, int nyv, int nzv, int mx1, int my1, int mxyz1,
                 int ntmax, int *irc) {
   gjppostf3l_(ppart,cu,kpic,ncl,ihole,&qm,&dt,&nppmx,&idimp,&nx,&ny,
               &nz,&mx,&my,&mz,&nxv,&nyv,&nzv,&mx1,&my1,&mxyz1,&ntmax,
               irc);
   return;
}

/*--------------------------------------------------------------------*/
void cgrjppost3l(float ppart[], float cu[], int kpic[], float qm,
                 float dt, float ci, int nppmx, int idimp, int nx,
                 int ny, int nz, int mx, int my, int mz, int nxv,
                 int nyv, int nzv, int mx1, int my1, int mxyz1,
                 int ipbc) {
   grjppost3l_(ppart,cu,kpic,&qm,&dt,&ci,&nppmx,&idimp,&nx,&ny,&nz,&mx,
               &my,&mz,&nxv,&nyv,&nzv,&mx1,&my1,&mxyz1,&ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cgrjppostf3l(float ppart[], float cu[], int kpic[], int ncl[],
                  int ihole[], float qm, float dt, float ci, int nppmx,
                  int idimp, int nx, int ny, int nz, int mx, int my,
                  int mz, int nxv, int nyv, int nzv, int mx1, int my1,
                  int mxyz1, int ntmax, int *irc) {
   grjppostf3l_(ppart,cu,kpic,ncl,ihole,&qm,&dt,&ci,&nppmx,&idimp,&nx,
                &ny,&nz,&mx,&my,&mz,&nxv,&nyv,&nzv,&mx1,&my1,&mxyz1,
                &ntmax,irc);
   return;
}

/*--------------------------------------------------------------------*/
void cpporder3l(float ppart[], float ppbuff[], int kpic[], int ncl[],
                int ihole[], int idimp, int nppmx, int nx, int ny,
                int nz, int mx, int my, int mz, int mx1, int my1,
                int mz1, int npbmx, int ntmax, int *irc) {
   pporder3l_(ppart,ppbuff,kpic,ncl,ihole,&idimp,&nppmx,&nx,&ny,&nz,&mx,
              &my,&mz,&mx1,&my1,&mz1,&npbmx,&ntmax,irc);
   return;
}

/*--------------------------------------------------------------------*/
void cpporderf3l(float ppart[], float ppbuff[], int kpic[], int ncl[],
                 int ihole[], int idimp, int nppmx, int mx1, int my1,
                 int mz1, int npbmx, int ntmax, int *irc) {
   pporderf3l_(ppart,ppbuff,kpic,ncl,ihole,&idimp,&nppmx,&mx1,&my1,&mz1,
               &npbmx,&ntmax,irc);
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
   aguard3l_(q,&nx,&ny,&nz,&nxe,&nye,&nze);
   return;
}

/*--------------------------------------------------------------------*/
void cmpois33(float complex q[], float complex fxyz[], int isign,
              float complex ffc[], float ax, float ay, float az,
              float affp, float *we, int nx, int ny, int nz, int nxvh,
              int nyv, int nzv, int nxhd, int nyhd, int nzhd) {
   mpois33_(q,fxyz,&isign,ffc,&ax,&ay,&az,&affp,we,&nx,&ny,&nz,&nxvh,
            &nyv,&nzv,&nxhd,&nyhd,&nzhd);
   return;
}

/*--------------------------------------------------------------------*/
void cmcuperp3(float complex cu[], int nx, int ny, int nz, int nxvh,
               int nyv, int nzv) {
   mcuperp3_(cu,&nx,&ny,&nz,&nxvh,&nyv,&nzv);
   return;
}

/*--------------------------------------------------------------------*/
void cmibpois33(float complex cu[], float complex bxyz[],
                float complex ffc[], float ci, float *wm, int nx,
                int ny, int nz, int nxvh, int nyv, int nzv, int nxhd,
                int nyhd, int nzhd) {
   mibpois33_(cu,bxyz,ffc,&ci,wm,&nx,&ny,&nz,&nxvh,&nyv,&nzv,&nxhd,
              &nyhd,&nzhd);
   return;
}

/*--------------------------------------------------------------------*/
void cmmaxwel3(float complex exyz[], float complex bxyz[],
               float complex cu[], float complex ffc[], float ci,
               float dt, float *wf, float *wm, int nx, int ny, int nz,
               int nxvh, int nyv, int nzv, int nxhd, int nyhd, 
               int nzhd) {
   mmaxwel3_(exyz,bxyz,cu,ffc,&ci,&dt,wf,wm,&nx,&ny,&nz,&nxvh,&nyv,&nzv,
             &nxhd,&nyhd,&nzhd);
   return;
}

/*--------------------------------------------------------------------*/
void cmemfield3(float complex fxyz[], float complex exyz[],
                float complex ffc[], int isign, int nx, int ny, int nz,
                int nxvh, int nyv, int nzv, int nxhd, int nyhd,
                int nzhd) {
   memfield3_(fxyz,exyz,ffc,&isign,&nx,&ny,&nz,&nxvh,&nyv,&nzv,&nxhd,
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
void cfft3rmxy(float complex f[], int isign, int mixup[],
               float complex sct[], int indx, int indy, int indz,
               int nzi, int nzp, int nxhd, int nyd, int nzd, int nxhyzd,
               int nxyzhd) {
   fft3rmxy_(f,&isign,mixup,sct,&indx,&indy,&indz,&nzi,&nzp,&nxhd,&nyd,
             &nzd,&nxhyzd,&nxyzhd);

   return;
}

/*--------------------------------------------------------------------*/
void cfft3rmxz(float complex f[], int isign, int mixup[],
               float complex sct[], int indx, int indy, int indz,
               int nyi, int nyp, int nxhd, int nyd, int nzd, int nxhyzd,
               int nxyzhd) {
   fft3rmxz_(f,&isign,mixup,sct,&indx,&indy,&indz,&nyi,&nyp,&nxhd,&nyd,
             &nzd,&nxhyzd,&nxyzhd);
   return;
}

/*--------------------------------------------------------------------*/
void cfft3rm3xy(float complex f[], int isign, int mixup[],
                float complex sct[], int indx, int indy, int indz,
                int nzi, int nzp, int nxhd, int nyd, int nzd,
                int nxhyzd, int nxyzhd) {
   fft3rm3xy_(f,&isign,mixup,sct,&indx,&indy,&indz,&nzi,&nzp,&nxhd,&nyd,
              &nzd,&nxhyzd,&nxyzhd);
   return;
}

/*--------------------------------------------------------------------*/
void cfft3rm3z(float complex f[], int isign, int mixup[],
               float complex sct[], int indx, int indy, int indz,
               int nyi, int nyp, int nxhd, int nyd, int nzd, int nxhyzd,
               int nxyzhd) {
   fft3rm3z_(f,&isign,mixup,sct,&indx,&indy,&indz,&nyi,&nyp,&nxhd,&nyd,
             &nzd,&nxhyzd,&nxyzhd);
   return;
}

/*--------------------------------------------------------------------*/
void cwfft3rmx(float complex f[], int isign, int mixup[],
               float complex sct[], int indx, int indy, int indz,
               int nxhd, int nyd, int nzd, int nxhyzd, int nxyzhd) {
   wfft3rmx_(f,&isign,mixup,sct,&indx,&indy,&indz,&nxhd,&nyd,&nzd,
             &nxhyzd,&nxyzhd);
   return;
}

/*--------------------------------------------------------------------*/
void cwfft3rm3(float complex f[], int isign, int mixup[],
               float complex sct[], int indx, int indy, int indz,
               int nxhd, int nyd, int nzd, int nxhyzd, int nxyzhd) {
   wfft3rm3_(f,&isign,mixup,sct,&indx,&indy,&indz,&nxhd,&nyd,&nzd,
             &nxhyzd,&nxyzhd);
   return;
}
