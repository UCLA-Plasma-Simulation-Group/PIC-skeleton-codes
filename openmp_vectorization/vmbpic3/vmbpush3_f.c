/* C Library for Skeleton 3D Electromagnetic OpenMP/Vector PIC Code */
/* Wrappers for calling the Fortran routines from a C main program */
/* written by Viktor K. Decyk, UCLA */

#include <complex.h>

double ranorm_();

/* Interfaces to Fortran */

void distr3_(float *part, float *vtx, float *vty, float *vtz,
             float *vdx, float *vdy, float *vdz, int *npx, int *npy,
             int *npz, int *idimp, int *nop, int *nx, int *ny, int *nz, 
             int *ipbc);

void dblkp3l_(float *part, int *kpic, int *nppmx, int *idimp, int *nop,
              int *mx, int *my, int *mz, int *mx1, int *my1,
              int *mxyz1, int *irc);

void ppmovin3lt_(float *part, float *ppart, int *kpic, int *nppmx,
                 int *idimp, int *nop, int *mx, int *my, int *mz,
                 int *mx1, int *my1, int *mxyz1, int *irc);

void ppmovin3ltp_(float *part, float *ppart, int *kpic, int *kp,
                  int *nppmx, int *idimp, int *nop, int *mx, int *my,
                  int *mz, int *mx1, int *my1, int *mxyz1, int *irc);

void ppcheck3lt_(float *ppart, int *kpic, int *idimp, int *nppmx,
                 int *nx, int *ny, int *nz, int *mx, int *my, int *mz,
                 int *mx1, int *my1, int *mz1, int *irc);

void gbppush3lt_(float *ppart, float *fxyz, float *bxyz, int *kpic,
                 float *qbm, float *dt, float *dtc, float *ek,
                 int *idimp, int *nppmx, int *nx, int *ny, int *nz,
                 int *mx, int *my, int *mz, int *nxv, int *nyv,
                 int *nzv, int *mx1, int *my1, int *mxyz1, int *ipbc);

void gbppushf3lt_(float *ppart, float *fxyz, float *bxyz, int *kpic,
                  int *ncl, int *ihole, float *qbm, float *dt,
                  float *dtc, float *ek, int *idimp, int *nppmx,
                  int *nx, int *ny, int *nz, int *mx, int *my, int *mz,
                  int *nxv, int *nyv, int *nzv, int *mx1, int *my1, 
                  int *mxyz1, int *ntmax, int *irc);

void grbppush3lt_(float *ppart, float *fxyz, float *bxyz, int *kpic,
                  float *qbm, float *dt, float *dtc, float *ci,
                  float *ek, int *idimp, int *nppmx, int *nx, int *ny,
                  int *nz, int *mx, int *my, int *mz, int *nxv,
                  int *nyv, int *nzv, int *mx1, int *my1, int *mxyz1,
                  int *ipbc);

void grbppushf3lt_(float *ppart, float *fxyz, float *bxyz, int *kpic,
                   int *ncl, int *ihole, float *qbm, float *dt,
                   float *dtc, float *ci, float *ek, int *idimp,
                   int *nppmx, int *nx, int *ny, int *nz, int *mx,
                   int *my, int *mz, int *nxv, int *nyv, int *nzv,
                   int *mx1, int *my1, int *mxyz1, int *ntmax,
                   int *irc);

void vgbppush3lt_(float *ppart, float *fxyz, float *bxyz, int *kpic,
                  float *qbm, float *dt, float *dtc, float *ek,
                  int *idimp, int *nppmx, int *nx, int *ny, int *nz,
                  int *mx, int *my, int *mz, int *nxv, int *nyv,
                  int *nzv, int *mx1, int *my1, int *mxyz1,
                  int *ipbc);

void vgbppushf3lt_(float *ppart, float *fxyz, float *bxyz, int *kpic,
                   int *ncl, int *ihole, float *qbm, float *dt,
                   float *dtc, float *ek, int *idimp, int *nppmx,
                   int *nx, int *ny, int *nz, int *mx, int *my,
                   int *mz, int *nxv, int *nyv, int *nzv, int *mx1,
                   int *my1, int *mxyz1, int *ntmax, int *irc);

void vgrbppush3lt_(float *ppart, float *fxyz, float *bxyz, int *kpic,
                   float *qbm, float *dt, float *dtc, float *ci,
                   float *ek, int *idimp, int *nppmx, int *nx, int *ny,
                   int *nz, int *mx, int *my, int *mz, int *nxv,
                   int *nyv, int *nzv, int *mx1, int *my1, int *mxyz1, 
                   int *ipbc);

void vgrbppushf3lt_(float *ppart, float *fxyz, float *bxyz, int *kpic,
                    int *ncl, int *ihole, float *qbm, float *dt,
                    float *dtc, float *ci, float *ek, int *idimp,
                    int *nppmx, int *nx, int *ny, int *nz, int *mx,
                    int *my, int *mz, int *nxv, int *nyv, int *nzv,
                    int *mx1, int *my1, int *mxyz1, int *ntmax,
                    int *irc);

void v2gbppush3lt_(float *ppart, float *fxyz, float *bxyz, int *kpic,
                   float *qbm, float *dt, float *dtc, float *ek,
                   int *idimp, int *nppmx, int *nx, int *ny, int *nz,
                   int *mx, int *my, int *mz, int *nxv, int *nyv,
                   int *nzv, int *mx1, int *my1, int *mxyz1,
                   int *ipbc);

void v2gbppushf3lt_(float *ppart, float *fxyz, float *bxyz, int *kpic,
                    int *ncl, int *ihole, float *qbm, float *dt,
                    float *dtc, float *ek, int *idimp, int *nppmx,
                    int *nx, int *ny, int *nz, int *mx, int *my,
                    int *mz, int *nxv, int *nyv, int *nzv, int *mx1,
                    int *my1, int *mxyz1, int *ntmax, int *irc);

void v2grbppush3lt_(float *ppart, float *fxyz, float *bxyz, int *kpic,
                    float *qbm, float *dt, float *dtc, float *ci,
                    float *ek, int *idimp, int *nppmx, int *nx,
                    int *ny, int *nz, int *mx, int *my, int *mz,
                    int *nxv, int *nyv, int *nzv, int *mx1, int *my1,
                    int *mxyz1, int *ipbc);

void v2grbppushf3lt_(float *ppart, float *fxyz, float *bxyz, int *kpic,
                     int *ncl, int *ihole, float *qbm, float *dt,
                     float *dtc, float *ci, float *ek, int *idimp,
                     int *nppmx, int *nx, int *ny, int *nz, int *mx,
                     int *my, int *mz, int *nxv, int *nyv, int *nzv, 
                     int *mx1, int *my1, int *mxyz1, int *ntmax,
                     int *irc);

void gppost3lt_(float *ppart, float *q, int *kpic, float *qm,
                int *nppmx, int *idimp, int *mx, int *my, int *mz,
                int *nxv, int *nyv, int *nzv, int *mx1, int *my1,
                int *mxyz1);

void vgppost3lt_(float *ppart, float *q, int *kpic, float *qm,
                 int *nppmx, int *idimp, int *mx, int *my, int *mz,
                 int *nxv, int *nyv, int *nzv, int *mx1, int *my1,
                 int *mxyz1);

void gjppost3lt_(float *ppart, float *cu, int *kpic, float *qm,
                 float *dt, int *nppmx, int *idimp, int *nx, int *ny, 
                 int *nz, int *mx, int *my, int *mz, int *nxv,
                 int *nyv, int *nzv, int *mx1, int *my1, int *mxyz1,
                 int *ipbc);

void gjppostf3lt_(float *ppart, float *cu, int *kpic, int *ncl,
                  int *ihole, float *qm, float *dt, int *nppmx,
                  int *idimp, int *nx, int *ny, int *nz, int *mx,
                  int *my, int *mz, int *nxv, int *nyv, int *nzv,
                  int *mx1, int *my1,  int *mxyz1, int *ntmax,
                  int *irc);

void grjppost3lt_(float *ppart, float *cu, int *kpic, float *qm,
                  float *dt, float *ci, int *nppmx, int *idimp,
                  int *nx, int *ny, int *nz, int *mx, int *my, int *mz,
                  int *nxv, int *nyv, int *nzv, int *mx1, int *my1,
                  int *mxyz1, int *ipbc);

void grjppostf3lt_(float *ppart, float *cu, int *kpic, int *ncl,
                   int *ihole, float *qm, float *dt, float *ci,
                   int *nppmx, int *idimp, int *nx, int *ny, int *nz,
                   int *mx, int *my, int *mz, int *nxv, int *nyv,
                   int *nzv, int *mx1, int *my1, int *mxyz1,
                   int *ntmax, int *irc);

void vgjppost3lt_(float *ppart, float *cu, int *kpic, float *qm,
                  float *dt, int *nppmx, int *idimp, int *nx, int *ny,
                  int *nz, int *mx, int *my, int *mz, int *nxv,
                  int *nyv, int *nzv, int *mx1, int *my1, int *mxyz1,
                  int *ipbc);

void vgjppostf3lt_(float *ppart, float *cu, int *kpic, int *ncl,
                   int *ihole, float *qm, float *dt, int *nppmx, 
                   int *idimp, int *nx, int *ny, int *nz, int *mx,
                   int *my, int *mz, int *nxv, int *nyv, int *nzv,
                   int *mx1, int *my1, int *mxyz1, int *ntmax,
                   int *irc);

void vgrjppost3lt_(float *ppart, float *cu, int *kpic, float *qm,
                   float *dt, float *ci, int *nppmx, int *idimp,
                   int *nx, int *ny, int *nz, int *mx, int *my,
                   int *mz, int *nxv, int *nyv, int *nzv, int *mx1,
                   int *my1,int *mxyz1, int *ipbc);

void vgrjppostf3lt_(float *ppart, float *cu, int *kpic, int *ncl,
                    int *ihole, float *qm, float *dt, float *ci,
                    int *nppmx, int *idimp, int *nx, int *ny, int *nz,
                    int *mx, int *my, int *mz, int *nxv, int *nyv,
                    int *nzv, int *mx1, int *my1, int *mxyz1, 
                    int *ntmax, int *irc);

void viscan2_(int *isdata, int *mb, int *nths);

void pporder3lt_(float *ppart, float *ppbuff, int *kpic, int *ncl,
                 int *ihole, int *idimp, int *nppmx, int *nx, int *ny,
                 int *nz, int *mx, int *my, int *mz, int *mx1,
                 int *my1, int *mz1, int *npbmx, int *ntmax,
                 int *irc);

void pporderf3lt_(float *ppart, float *ppbuff, int *kpic, int *ncl,
                  int *ihole, int *idimp, int *nppmx, int *mx1,
                  int *my1, int *mz1, int *npbmx, int *ntmax,
                  int *irc);

void vpporder3lt_(float *ppart, float *ppbuff, int *kpic, int *ncl,
                  int *ihole, int *idimp, int *nppmx, int *nx, int *ny,
                  int *nz, int *mx, int *my, int *mz, int *mx1,
                  int *my1, int *mz1, int *npbmx, int *ntmax,
                  int *irc);

void vpporderf3lt_(float *ppart, float *ppbuff, int *kpic, int *ncl,
                   int *ihole, int *idimp, int *nppmx, int *mx1,
                   int *my1, int *mz1, int *npbmx, int *ntmax,
                   int *irc);

void v2pporderf3lt_(float *ppart, float *ppbuff, int *kpic, int *ncl,
                    int *ihole, int *idimp, int *nppmx, int *mx1, 
                    int *my1, int *mz1, int *npbmx, int *ntmax,
                    int *irc);

void cguard3l_(float *fxyz, int *nx, int *ny, int *nz, int *nxe,
               int *nye, int *nze);

void acguard3l_(float *cu, int *nx, int *ny, int *nz, int *nxe,
                int *nye, int *nze);

void aguard3l_(float *q, int *nx, int *ny, int *nz, int *nxe, int *nye,
               int *nze);

void vmpois33_(float complex *q, float complex *fxyz, int *isign,
               float complex *ffc, float *ax, float *ay, float *az,
               float *affp, float *we, int *nx, int *ny, int *nz,
               int *nxvh, int *nyv, int *nzv, int *nxhd, int *nyhd,
               int *nzhd);

void mcuperp3_(float complex *cu, int *nx, int *ny, int *nz, int *nxvh,
               int *nyv, int *nzv);

void vmibpois33_(float complex *cu, float complex* bxyz,
                 float complex *ffc, float *ci, float *wm, int *nx,
                 int *ny, int *nz, int *nxvh, int *nyv, int *nzv,
                 int *nxhd, int *nyhd, int *nzhd);

void vmmaxwel3_(float complex *exyz, float complex *bxyz,
                float complex *cu, float complex *ffc, float *ci,
                float *dt, float *wf, float *wm, int *nx, int *ny,
                int *nz, int *nxvh, int *nyv, int *nzv, int *nxhd, 
                int *nyhd, int *nzhd);

void vmemfield3_(float complex *fxyz, float complex *exyz,
                 float complex *ffc, int *isign, int *nx, int *ny,
                 int *nz, int *nxvh, int *nyv, int *nzv, int *nxhd,
                 int *nyhd, int *nzhd);

void wfft3rinit_(int *mixup, float complex *sct, int *indx, int *indy,
                 int *indz, int *nxhyzd, int *nxyzhd);

void fft3rvmxy_(float complex *f, int *isign, int *mixup,
                float complex *sct, int *indx, int *indy, int *indz,
                int *nzi, int *nzp, int *nxhd, int *nyd, int *nzd,
                int *nxhyzd, int *nxyzhd);

void fft3rvmxz_(float complex *f, int *isign, int *mixup,
                float complex *sct, int *indx, int *indy, int *indz,
                int *nyi, int *nyp, int *nxhd, int *nyd, int *nzd,
                int *nxhyzd, int *nxyzhd);

void fft3rvm3xy_(float complex *f, int *isign, int *mixup,
                 float complex *sct, int *indx, int *indy, int *indz,
                 int *nzi, int *nzp, int *nxhd, int *nyd, int *nzd,
                 int *nxhyzd, int *nxyzhd);

void fft3rvm3z_(float complex *f, int *isign, int *mixup,
                float complex *sct, int *indx, int *indy, int *indz,
                int *nyi, int *nyp, int *nxhd, int *nyd, int *nzd,
                int *nxhyzd, int *nxyzhd);

void wfft3rvmx_(float complex *f, int *isign, int *mixup,
                float complex *sct, int *indx, int *indy, int *indz,
                int *nxhd, int *nyd, int *nzd, int *nxhyzd,
                int *nxyzhd);

void wfft3rvm3_(float complex *f, int *isign, int *mixup,
                float complex *sct, int *indx, int *indy, int *indz,
                int *nxhd, int *nyd, int *nzd, int *nxhyzd,
                int *nxyzhd);

void set_szero3_(float *q, int *mx, int *my, int *mz, int *nxv,
                 int *nyv, int *nzv, int *mx1, int *my1, int *mxyz1);

void set_vzero3_(float *cu, int *mx, int *my, int *mz, int *ndim,
                 int *nxv, int *nyv, int *nzv, int *mx1, int *my1, 
                 int *mxyz1);

void set_cvzero3_(float complex *exyz, int *nx, int *ny, int *nz,
                  int *ndim, int *nxvh, int *nyv, int *nzv);

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
void cppmovin3lt(float part[], float ppart[], int kpic[], int nppmx,
                 int idimp, int nop, int mx, int my, int mz, int mx1,
                 int my1, int mxyz1, int *irc) {
   ppmovin3lt_(part,ppart,kpic,&nppmx,&idimp,&nop,&mx,&my,&mz,&mx1,&my1,
               &mxyz1,irc);
   return;
}

/*--------------------------------------------------------------------*/
void cppmovin3ltp(float part[], float ppart[], int kpic[], int kp[],
                  int nppmx, int idimp, int nop, int mx, int my, int mz,
                  int mx1, int my1, int mxyz1, int *irc) {
   ppmovin3ltp_(part,ppart,kpic,kp,&nppmx,&idimp,&nop,&mx,&my,&mz,&mx1,
                &my1,&mxyz1,irc);
   return;
}

/*--------------------------------------------------------------------*/
void cppcheck3lt(float ppart[], int kpic[], int idimp, int nppmx,
                 int nx, int ny, int nz, int mx, int my, int mz,
                 int mx1, int my1, int mz1, int *irc) {
   ppcheck3lt_(ppart,kpic,&idimp,&nppmx,&nx,&ny,&nz,&mx,&my,&mz,&mx1,
               &my1,&mz1,irc);
   return;
}

/*--------------------------------------------------------------------*/
void cgbppush3lt(float ppart[], float fxyz[], float bxyz[], int kpic[],
                 float qbm, float dt, float dtc, float *ek, int idimp,
                 int nppmx, int nx, int ny, int nz, int mx, int my,
                 int mz, int nxv, int nyv, int nzv, int mx1, int my1,
                 int mxyz1, int ipbc) {
   gbppush3lt_(ppart,fxyz,bxyz,kpic,&qbm,&dt,&dtc,ek,&idimp,&nppmx,&nx,
               &ny,&nz,&mx,&my,&mz,&nxv,&nyv,&nzv,&mx1,&my1,&mxyz1,
               &ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cgbppushf3lt(float ppart[], float fxyz[], float bxyz[], int kpic[],
                  int ncl[], int ihole[], float qbm, float dt, 
                  float dtc, float *ek, int idimp, int nppmx, int nx,
                  int ny, int nz, int mx, int my, int mz, int nxv,
                  int nyv, int nzv, int mx1, int my1, int mxyz1,
                  int ntmax, int *irc) {
   gbppushf3lt_(ppart,fxyz,bxyz,kpic,ncl,ihole,&qbm,&dt,&dtc,ek,&idimp,
                &nppmx,&nx,&ny,&nz,&mx,&my,&mz,&nxv,&nyv,&nzv,&mx1,&my1,
                &mxyz1,&ntmax,irc);
   return;
}

/*--------------------------------------------------------------------*/
void cgrbppush3lt(float ppart[], float fxyz[], float bxyz[], int kpic[],
                  float qbm, float dt, float dtc, float ci, float *ek,
                  int idimp, int nppmx, int nx, int ny, int nz, int mx,
                  int my, int mz, int nxv, int nyv, int nzv, int mx1,
                  int my1, int mxyz1, int ipbc) {
   grbppush3lt_(ppart,fxyz,bxyz,kpic,&qbm,&dt,&dtc,&ci,ek,&idimp,&nppmx,
                &nx,&ny,&nz,&mx,&my,&mz,&nxv,&nyv,&nzv,&mx1,&my1,&mxyz1,
                &ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cgrbppushf3lt(float ppart[], float fxyz[], float bxyz[],
                   int kpic[], int ncl[], int ihole[], float qbm,
                   float dt, float dtc, float ci, float *ek, int idimp,
                   int nppmx, int nx, int ny, int nz, int mx, int my,
                   int mz, int nxv, int nyv, int nzv, int mx1, int my1,
                   int mxyz1, int ntmax, int *irc) {
   grbppushf3lt_(ppart,fxyz,bxyz,kpic,ncl,ihole,&qbm,&dt,&dtc,&ci,ek,
                 &idimp,&nppmx,&nx,&ny,&nz,&mx,&my,&mz,&nxv,&nyv,&nzv,
                 &mx1,&my1,&mxyz1,&ntmax,irc);
   return;
}

/*--------------------------------------------------------------------*/
void cvgbppush3lt(float ppart[], float fxyz[], float bxyz[], int kpic[],
                  float qbm, float dt, float dtc, float *ek, int idimp,
                  int nppmx, int nx, int ny, int nz, int mx, int my,
                  int mz, int nxv, int nyv, int nzv, int mx1, int my1,
                  int mxyz1, int ipbc) {
   vgbppush3lt_(ppart,fxyz,bxyz,kpic,&qbm,&dt,&dtc,ek,&idimp,&nppmx,&nx,
                &ny,&nz,&mx,&my,&mz,&nxv,&nyv,&nzv,&mx1,&my1,&mxyz1,
                &ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cvgbppushf3lt(float ppart[], float fxyz[], float bxyz[], 
                   int kpic[], int ncl[], int ihole[], float qbm, 
                   float dt, float dtc, float *ek, int idimp, int nppmx,
                   int nx, int ny, int nz, int mx, int my, int mz,
                   int nxv, int nyv, int nzv, int mx1, int my1,
                   int mxyz1, int ntmax, int *irc) {
   vgbppushf3lt_(ppart,fxyz,bxyz,kpic,ncl,ihole,&qbm,&dt,&dtc,ek,&idimp,
                 &nppmx,&nx,&ny,&nz,&mx,&my,&mz,&nxv,&nyv,&nzv,&mx1,
                 &my1,&mxyz1,&ntmax,irc);
   return;
}

/*--------------------------------------------------------------------*/
void cvgrbppush3lt(float ppart[], float fxyz[], float bxyz[],
                   int kpic[], float qbm, float dt, float dtc, float ci,
                   float *ek, int idimp, int nppmx, int nx, int ny,
                   int nz, int mx, int my, int mz, int nxv, int nyv,
                   int nzv, int mx1, int my1, int mxyz1, int ipbc) {
   vgrbppush3lt_(ppart,fxyz,bxyz,kpic,&qbm,&dt,&dtc,&ci,ek,&idimp,
                 &nppmx,&nx,&ny,&nz,&mx,&my,&mz,&nxv,&nyv,&nzv,&mx1,
                 &my1,&mxyz1,&ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cvgrbppushf3lt(float ppart[], float fxyz[], float bxyz[],
                    int kpic[], int ncl[], int ihole[], float qbm,
                    float dt, float dtc, float ci, float *ek, int idimp,
                    int nppmx, int nx, int ny, int nz, int mx, int my,
                    int mz, int nxv, int nyv, int nzv, int mx1, int my1,
                    int mxyz1, int ntmax, int *irc) {
   vgrbppushf3lt_(ppart,fxyz,bxyz,kpic,ncl,ihole,&qbm,&dt,&dtc,&ci,ek,
                  &idimp,&nppmx,&nx,&ny,&nz,&mx,&my,&mz,&nxv,&nyv,&nzv,
                  &mx1,&my1,&mxyz1,&ntmax,irc);
   return;
}

/*--------------------------------------------------------------------*/
void cv2gbppush3lt(float ppart[], float fxyz[], float bxyz[], 
                   int kpic[], float qbm, float dt, float dtc,
                   float *ek, int idimp, int nppmx, int nx, int ny, 
                   int nz, int mx, int my, int mz, int nxv, int nyv,
                   int nzv, int mx1, int my1, int mxyz1, int ipbc) {
   v2gbppush3lt_(ppart,fxyz,bxyz,kpic,&qbm,&dt,&dtc,ek,&idimp,&nppmx,
                 &nx,&ny,&nz,&mx,&my,&mz,&nxv,&nyv,&nzv,&mx1,&my1,
                 &mxyz1,&ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cv2gbppushf3lt(float ppart[], float fxyz[], float bxyz[], 
                    int kpic[], int ncl[], int ihole[], float qbm, 
                    float dt, float dtc, float *ek, int idimp,
                    int nppmx, int nx, int ny, int nz, int mx, int my,
                    int mz, int nxv, int nyv, int nzv, int mx1, int my1,
                    int mxyz1, int ntmax, int *irc) {
   v2gbppushf3lt_(ppart,fxyz,bxyz,kpic,ncl,ihole,&qbm,&dt,&dtc,ek,
                  &idimp,&nppmx,&nx,&ny,&nz,&mx,&my,&mz,&nxv,&nyv,&nzv,
                  &mx1,&my1,&mxyz1,&ntmax,irc);
   return;
}

/*--------------------------------------------------------------------*/
void cv2grbppush3lt(float ppart[], float fxyz[], float bxyz[],
                    int kpic[], float qbm, float dt, float dtc,
                    float ci, float *ek, int idimp, int nppmx, int nx,
                    int ny, int nz, int mx, int my, int mz, int nxv,
                    int nyv, int nzv, int mx1, int my1, int mxyz1,
                    int ipbc) {
   v2grbppush3lt_(ppart,fxyz,bxyz,kpic,&qbm,&dt,&dtc,&ci,ek,&idimp,
                  &nppmx,&nx,&ny,&nz,&mx,&my,&mz,&nxv,&nyv,&nzv,&mx1,
                  &my1,&mxyz1,&ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cv2grbppushf3lt(float ppart[], float fxyz[], float bxyz[],
                     int kpic[], int ncl[], int ihole[], float qbm,
                     float dt, float dtc, float ci, float *ek,
                     int idimp, int nppmx, int nx, int ny, int nz,
                     int mx, int my, int mz, int nxv, int nyv, int nzv,
                     int mx1, int my1, int mxyz1, int ntmax, int *irc) {
   v2grbppushf3lt_(ppart,fxyz,bxyz,kpic,ncl,ihole,&qbm,&dt,&dtc,&ci,ek,
                   &idimp,&nppmx,&nx,&ny,&nz,&mx,&my,&mz,&nxv,&nyv,&nzv, 
                   &mx1,&my1,&mxyz1,&ntmax,irc);
   return;
}

/*--------------------------------------------------------------------*/
void cgppost3lt(float ppart[], float q[], int kpic[], float qm,
                int nppmx, int idimp, int mx, int my, int mz, int nxv,
                int nyv, int nzv, int mx1, int my1, int mxyz1) {
   gppost3lt_(ppart,q,kpic,&qm,&nppmx,&idimp,&mx,&my,&mz,&nxv,&nyv,&nzv,
              &mx1,&my1,&mxyz1);
   return;
}

/*--------------------------------------------------------------------*/
void cvgppost3lt(float ppart[], float q[], int kpic[], float qm,
                 int nppmx, int idimp, int mx, int my, int mz, int nxv,
                 int nyv, int nzv, int mx1, int my1, int mxyz1) {
   vgppost3lt_(ppart,q,kpic,&qm,&nppmx,&idimp,&mx,&my,&mz,&nxv,&nyv,
               &nzv,&mx1,&my1,&mxyz1);
   return;
}

/*--------------------------------------------------------------------*/
void cgjppost3lt(float ppart[], float cu[], int kpic[], float qm,
                 float dt, int nppmx, int idimp, int nx, int ny, int nz,
                 int mx, int my, int mz, int nxv, int nyv, int nzv,
                 int mx1, int my1, int mxyz1, int ipbc) {
   gjppost3lt_(ppart,cu,kpic,&qm,&dt,&nppmx,&idimp,&nx,&ny,&nz,&mx,&my,
               &mz,&nxv,&nyv,&nzv,&mx1,&my1,&mxyz1,&ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cgjppostf3lt(float ppart[], float cu[], int kpic[], int ncl[],
                  int ihole[], float qm, float dt, int nppmx, int idimp,
                  int nx, int ny, int nz, int mx, int my, int mz,
                  int nxv, int nyv, int nzv, int mx1, int my1, 
                  int mxyz1, int ntmax, int *irc) {
   gjppostf3lt_(ppart,cu,kpic,ncl,ihole,&qm,&dt,&nppmx,&idimp,&nx,&ny,
                &nz,&mx,&my,&mz,&nxv,&nyv,&nzv,&mx1,&my1,&mxyz1,&ntmax,
                irc);
   return;
}

/*--------------------------------------------------------------------*/
void cgrjppost3lt(float ppart[], float cu[], int kpic[], float qm,
                  float dt, float ci, int nppmx, int idimp, int nx,
                  int ny, int nz, int mx, int my, int mz, int nxv,
                  int nyv, int nzv, int mx1, int my1, int mxyz1,
                  int ipbc) {
   grjppost3lt_(ppart,cu,kpic,&qm,&dt,&ci,&nppmx,&idimp,&nx,&ny,&nz,&mx,
                &my,&mz,&nxv,&nyv,&nzv,&mx1,&my1,&mxyz1,&ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cgrjppostf3lt(float ppart[], float cu[], int kpic[], int ncl[],
                   int ihole[], float qm, float dt, float ci, int nppmx,
                   int idimp, int nx, int ny, int nz, int mx, int my,
                   int mz, int nxv, int nyv, int nzv, int mx1, int my1,
                   int mxyz1, int ntmax, int *irc) {
   grjppostf3lt_(ppart,cu,kpic,ncl,ihole,&qm,&dt,&ci,&nppmx,&idimp,&nx,
                 &ny,&nz,&mx,&my,&mz,&nxv,&nyv,&nzv,&mx1,&my1,&mxyz1,
                 &ntmax,irc);
   return;
}

/*--------------------------------------------------------------------*/
void cvgjppost3lt(float ppart[], float cu[], int kpic[], float qm,
                  float dt, int nppmx, int idimp, int nx, int ny,
                  int nz, int mx, int my, int mz, int nxv, int nyv,
                  int nzv, int mx1, int my1, int mxyz1, int ipbc) {
   vgjppost3lt_(ppart,cu,kpic,&qm,&dt,&nppmx,&idimp,&nx,&ny,&nz,&mx,&my,
                &mz,&nxv,&nyv,&nzv,&mx1,&my1,&mxyz1,&ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cvgjppostf3lt(float ppart[], float cu[], int kpic[], int ncl[],
                   int ihole[], float qm, float dt, int nppmx, 
                   int idimp, int nx, int ny, int nz, int mx, int my,
                   int mz, int nxv, int nyv, int nzv, int mx1, int my1, 
                   int mxyz1, int ntmax, int *irc) {
   vgjppostf3lt_(ppart,cu,kpic,ncl,ihole,&qm,&dt,&nppmx,&idimp,&nx,&ny,
                 &nz,&mx,&my,&mz,&nxv,&nyv,&nzv,&mx1,&my1,&mxyz1,&ntmax,
                 irc);
   return;
}

/*--------------------------------------------------------------------*/
void cvgrjppost3lt(float ppart[], float cu[], int kpic[], float qm,
                   float dt, float ci, int nppmx, int idimp, int nx,
                   int ny, int nz, int mx, int my, int mz, int nxv,
                   int nyv, int nzv, int mx1, int my1, int mxyz1,
                   int ipbc) {
   vgrjppost3lt_(ppart,cu,kpic,&qm,&dt,&ci,&nppmx,&idimp,&nx,&ny,&nz,
                 &mx,&my,&mz,&nxv,&nyv,&nzv,&mx1,&my1,&mxyz1,&ipbc);

   return;
}

/*--------------------------------------------------------------------*/
void cvgrjppostf3lt(float ppart[], float cu[], int kpic[], int ncl[],
                    int ihole[], float qm, float dt, float ci,
                    int nppmx, int idimp, int nx, int ny, int nz,
                    int mx, int my, int mz, int nxv, int nyv, int nzv,
                    int mx1, int my1, int mxyz1, int ntmax, int *irc) {
   vgrjppostf3lt_(ppart,cu,kpic,ncl,ihole,&qm,&dt,&ci,&nppmx,&idimp,&nx,
                  &ny,&nz,&mx,&my,&mz,&nxv,&nyv,&nzv,&mx1,&my1,&mxyz1,
                  &ntmax,irc);
   return;
}

/*--------------------------------------------------------------------*/
void cviscan2(int *isdata, int *mb, int nths) {
   viscan2_(isdata,mb,&nths);
   return;
}

/*--------------------------------------------------------------------*/
void cpporder3lt(float ppart[], float ppbuff[], int kpic[], int ncl[],
                 int ihole[], int idimp, int nppmx, int nx, int ny,
                 int nz, int mx, int my, int mz, int mx1, int my1,
                 int mz1, int npbmx, int ntmax, int *irc) {
   pporder3lt_(ppart,ppbuff,kpic,ncl,ihole,&idimp,&nppmx,&nx,&ny,&nz,
               &mx,&my,&mz,&mx1,&my1,&mz1,&npbmx,&ntmax,irc);
   return;
}

/*--------------------------------------------------------------------*/
void cpporderf3lt(float ppart[], float ppbuff[], int kpic[], int ncl[],
                  int ihole[], int idimp, int nppmx, int mx1, int my1,
                  int mz1, int npbmx, int ntmax, int *irc) {
   pporderf3lt_(ppart,ppbuff,kpic,ncl,ihole,&idimp,&nppmx,&mx1,&my1,
                &mz1,&npbmx,&ntmax,irc);
   return;
}

/*--------------------------------------------------------------------*/
void cvpporder3lt(float ppart[], float ppbuff[], int kpic[], int ncl[],
                  int ihole[], int idimp, int nppmx, int nx, int ny,
                  int nz, int mx, int my, int mz, int mx1, int my1,
                  int mz1, int npbmx, int ntmax, int *irc) {
   vpporder3lt_(ppart,ppbuff,kpic,ncl,ihole,&idimp,&nppmx,&nx,&ny,&nz,
                &mx,&my,&mz,&mx1,&my1,&mz1,&npbmx,&ntmax,irc);
   return;
}

/*--------------------------------------------------------------------*/
void cvpporderf3lt(float ppart[], float ppbuff[], int kpic[], int ncl[],
                   int ihole[], int idimp, int nppmx, int mx1, int my1,
                   int mz1, int npbmx, int ntmax, int *irc) {
   vpporderf3lt_(ppart,ppbuff,kpic,ncl,ihole,&idimp,&nppmx,&mx1,&my1,
                 &mz1,&npbmx,&ntmax,irc);
   return;
}

/*--------------------------------------------------------------------*/
void cv2pporderf3lt(float ppart[], float ppbuff[], int kpic[],
                    int ncl[], int ihole[], int idimp, int nppmx,
                    int mx1, int my1, int mz1, int npbmx, int ntmax,
                    int *irc) {
   v2pporderf3lt_(ppart,ppbuff,kpic,ncl,ihole,&idimp,&nppmx,&mx1,&my1,
                  &mz1,&npbmx,&ntmax,irc);
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
void cvmpois33(float complex q[], float complex fxyz[], int isign,
               float complex ffc[], float ax, float ay, float az,
               float affp, float *we, int nx, int ny, int nz, int nxvh,
               int nyv, int nzv, int nxhd, int nyhd, int nzhd) {
   vmpois33_(q,fxyz,&isign,ffc,&ax,&ay,&az,&affp,we,&nx,&ny,&nz,&nxvh,
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
void cvmibpois33(float complex cu[], float complex bxyz[],
                 float complex ffc[], float ci, float *wm, int nx,
                 int ny, int nz, int nxvh, int nyv, int nzv, int nxhd,
                 int nyhd, int nzhd) {
   vmibpois33_(cu,bxyz,ffc,&ci,wm,&nx,&ny,&nz,&nxvh,&nyv,&nzv,&nxhd,
               &nyhd,&nzhd);
   return;
}

/*--------------------------------------------------------------------*/
void cvmmaxwel3(float complex exyz[], float complex bxyz[],
                float complex cu[], float complex ffc[], float ci,
                float dt, float *wf, float *wm, int nx, int ny, int nz,
                int nxvh, int nyv, int nzv, int nxhd, int nyhd, 
                int nzhd) {
   vmmaxwel3_(exyz,bxyz,cu,ffc,&ci,&dt,wf,wm,&nx,&ny,&nz,&nxvh,&nyv,
              &nzv,&nxhd,&nyhd,&nzhd);
   return;
}

/*--------------------------------------------------------------------*/
void cvmemfield3(float complex fxyz[], float complex exyz[],
                 float complex ffc[], int isign, int nx, int ny, int nz,
                 int nxvh, int nyv, int nzv, int nxhd, int nyhd,
                 int nzhd) {
   vmemfield3_(fxyz,exyz,ffc,&isign,&nx,&ny,&nz,&nxvh,&nyv,&nzv,&nxhd,
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
void cfft3rvmxy(float complex f[], int isign, int mixup[],
                float complex sct[], int indx, int indy, int indz,
                int nzi, int nzp, int nxhd, int nyd, int nzd,
                int nxhyzd, int nxyzhd) {
   fft3rvmxy_(f,&isign,mixup,sct,&indx,&indy,&indz,&nzi,&nzp,&nxhd,&nyd,
              &nzd,&nxhyzd,&nxyzhd);
   return;
}

/*--------------------------------------------------------------------*/
void cfft3rvmxz(float complex f[], int isign, int mixup[],
                float complex sct[], int indx, int indy, int indz,
                int nyi, int nyp, int nxhd, int nyd, int nzd,
                int nxhyzd, int nxyzhd) {
   fft3rvmxz_(f,&isign,mixup,sct,&indx,&indy,&indz,&nyi,&nyp,&nxhd,&nyd,
              &nzd,&nxhyzd,&nxyzhd);
   return;
}

/*--------------------------------------------------------------------*/
void cfft3rvm3xy(float complex f[], int isign, int mixup[],
                 float complex sct[], int indx, int indy, int indz,
                 int nzi, int nzp, int nxhd, int nyd, int nzd,
                 int nxhyzd, int nxyzhd) {
   fft3rvm3xy_(f,&isign,mixup,sct,&indx,&indy,&indz,&nzi,&nzp,&nxhd,
               &nyd,&nzd,&nxhyzd,&nxyzhd);
   return;
}

/*--------------------------------------------------------------------*/
void cfft3rvm3z(float complex f[], int isign, int mixup[],
                float complex sct[], int indx, int indy, int indz,
                int nyi, int nyp, int nxhd, int nyd, int nzd,
                int nxhyzd, int nxyzhd) {
   fft3rvm3z_(f,&isign,mixup,sct,&indx,&indy,&indz,&nyi,&nyp,&nxhd,&nyd,
              &nzd,&nxhyzd,&nxyzhd);
   return;
}

/*--------------------------------------------------------------------*/
void cwfft3rvmx(float complex f[], int isign, int mixup[],
                float complex sct[], int indx, int indy, int indz,
                int nxhd, int nyd, int nzd, int nxhyzd, int nxyzhd) {
   wfft3rvmx_(f,&isign,mixup,sct,&indx,&indy,&indz,&nxhd,&nyd,&nzd,
              &nxhyzd,&nxyzhd);
   return;
}

/*--------------------------------------------------------------------*/
void cwfft3rvm3(float complex f[], int isign, int mixup[],
                float complex sct[], int indx, int indy, int indz,
                int nxhd, int nyd, int nzd, int nxhyzd, int nxyzhd) {
   wfft3rvm3_(f,&isign,mixup,sct,&indx,&indy,&indz,&nxhd,&nyd,&nzd,
              &nxhyzd,&nxyzhd);
   return;
}

/*--------------------------------------------------------------------*/
void cset_szero3(float q[], int mx, int my, int mz, int nxv, int nyv,
                 int nzv, int mx1, int my1, int mxyz1) {
   set_szero3_(q,&mx,&my,&mz,&nxv,&nyv,&nzv,&mx1,&my1,&mxyz1);
   return;
}

/*--------------------------------------------------------------------*/
void cset_vzero3(float cu[], int mx, int my, int mz, int ndim, int nxv,
                 int nyv, int nzv, int mx1, int my1, int mxyz1) {
   set_vzero3_(cu,&mx,&my,&mz,&ndim,&nxv,&nyv,&nzv,&mx1,&my1,&mxyz1);
   return;
}

/*--------------------------------------------------------------------*/
void cset_cvzero3(float complex exyz[], int nx, int ny, int nz,
                  int ndim, int nxvh, int nyv, int nzv) {
   set_cvzero3_(exyz,&nx,&ny,&nz,&ndim,&nxvh,&nyv,&nzv);
   return;
}

