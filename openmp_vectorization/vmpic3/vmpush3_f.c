/* C Library for Skeleton 3D Electrostatic OpenMP/Vector PIC Code */
/* Wrappers for calling the Fortran routines from a C main program */

#include <complex.h>

double ranorm_();

/*--------------------------------------------------------------------*/
void distr3_(float *part, float *vtx, float *vty, float *vtz,
             float *vdx, float *vdy, float *vdz, int *npx, int *npy,
             int *npz, int *idimp, int *nop, int *nx, int *ny, int *nz,
             int *ipbc);

/*--------------------------------------------------------------------*/
void dblkp3l_(float *part, int *kpic, int *nppmx, int *idimp, int *nop,
              int *mx, int *my, int *mz, int *mx1, int *my1, int *mxyz1,
              int *irc);

/*--------------------------------------------------------------------*/
void ppmovin3lt_(float *part, float *ppart, int *kpic, int *nppmx,
                 int *idimp, int *nop, int *mx, int *my, int *mz,
                 int *mx1, int *my1, int *mxyz1, int *irc);

/*--------------------------------------------------------------------*/
void ppmovin3ltp_(float *part, float *ppart, int *kpic, int *kp,
                  int *nppmx, int *idimp, int *nop, int *mx, int *my,
                  int *mz, int *mx1, int *my1, int *mxyz1, int *irc);

/*--------------------------------------------------------------------*/
void ppcheck3lt_(float *ppart, int *kpic, int *idimp, int *nppmx,
                 int *nx, int *ny, int *nz, int *mx, int *my, int *mz, 
                 int *mx1, int *my1, int *mz1, int *irc);

/*--------------------------------------------------------------------*/
void gppush3lt_(float *ppart, float *fxyz, int *kpic, float *qbm,
                float *dt, float *ek, int *idimp, int *nppmx, int *nx,
                int *ny, int *nz, int *mx, int *my, int *mz, int *nxv,
                int *nyv, int *nzv, int *mx1, int *my1, int *mxyz1,
                int *ipbc);

/*--------------------------------------------------------------------*/
void gppushf3lt_(float *ppart, float *fxyz, int *kpic, int *ncl,
                 int *ihole, float *qbm, float *dt, float *ek,
                 int *idimp, int *nppmx, int *nx, int *ny, int *nz,
                 int *mx, int *my, int *mz, int *nxv, int *nyv,
                 int *nzv, int *mx1, int *my1, int *mxyz1, int *ntmax,
                 int *irc);

/*--------------------------------------------------------------------*/
void vgppush3lt_(float *ppart, float *fxyz, int *kpic, float *qbm,
                 float *dt, float *ek, int *idimp, int *nppmx, int *nx,
                 int *ny, int *nz, int *mx, int *my, int *mz, int *nxv,
                 int *nyv, int *nzv, int *mx1, int *my1, int *mxyz1,
                 int *ipbc);

/*--------------------------------------------------------------------*/
void vgppushf3lt_(float *ppart, float *fxyz, int *kpic, int *ncl,
                  int *ihole, float *qbm, float *dt, float *ek,
                  int *idimp, int *nppmx, int *nx, int *ny, int *nz,
                  int *mx, int *my, int *mz, int *nxv, int *nyv,
                  int *nzv, int *mx1, int *my1, int *mxyz1, int *ntmax,
                  int *irc);

/*--------------------------------------------------------------------*/
void v2gppush3lt_(float *ppart, float *fxyz, int *kpic, float *qbm,
                  float *dt, float *ek, int *idimp, int *nppmx, int *nx,
                  int *ny, int *nz, int *mx, int *my, int *mz, int *nxv,
                  int *nyv, int *nzv, int *mx1, int *my1, int *mxyz1,
                  int *ipbc);

/*--------------------------------------------------------------------*/
void v2gppushf3lt_(float *ppart, float *fxyz, int *kpic, int *ncl,
                   int *ihole, float *qbm, float *dt, float *ek,
                   int *idimp, int *nppmx, int *nx, int *ny, int *nz,
                   int *mx, int *my, int *mz, int *nxv, int *nyv,
                   int *nzv, int *mx1, int *my1, int *mxyz1, int *ntmax,
                   int *irc);

/*--------------------------------------------------------------------*/
void gppost3lt_(float *ppart, float *q, int *kpic, float *qm,
                int *nppmx, int *idimp, int *mx, int *my, int *mz,
                int *nxv, int *nyv, int *nzv, int *mx1, int *my1,
                int *mxyz1);

/*--------------------------------------------------------------------*/
void vgppost3lt_(float *ppart, float *q, int *kpic, float *qm,
                 int *nppmx, int *idimp, int *mx, int *my, int *mz,
                 int *nxv, int *nyv, int *nzv, int *mx1, int *my1,
                 int *mxyz1);

/*--------------------------------------------------------------------*/
void viscan2_(int *isdata, int *mb, int *nths);

/*--------------------------------------------------------------------*/
void pporder3lt_(float *ppart, float *ppbuff, int *kpic, int *ncl,
                 int *ihole, int *idimp, int *nppmx, int *nx, int *ny,
                 int *nz, int *mx, int *my, int *mz, int *mx1, int *my1,
                 int *mz1, int *npbmx, int *ntmax, int *irc);

/*--------------------------------------------------------------------*/
void pporderf3lt_(float *ppart, float *ppbuff, int *kpic, int *ncl,
                  int *ihole, int *idimp, int *nppmx, int *mx1,
                  int *my1, int *mz1, int *npbmx, int *ntmax, int *irc);

/*--------------------------------------------------------------------*/
void vpporder3lt_(float *ppart, float *ppbuff, int *kpic, int *ncl,
                  int *ihole, int *idimp, int *nppmx, int *nx, int *ny,
                  int *nz, int *mx, int *my, int *mz, int *mx1,
                  int *my1, int *mz1, int *npbmx, int *ntmax, int *irc);

/*--------------------------------------------------------------------*/
void vpporderf3lt_(float *ppart, float *ppbuff, int *kpic, int *ncl,
                   int *ihole, int *idimp, int *nppmx, int *mx1, 
                   int *my1, int *mz1, int *npbmx, int *ntmax,
                   int *irc);

/*--------------------------------------------------------------------*/
void v2pporderf3lt_(float *ppart, float *ppbuff, int *kpic, int *ncl,
                    int *ihole, int *idimp, int *nppmx, int *mx1, 
                    int *my1, int *mz1, int *npbmx, int *ntmax,
                    int *irc);

/*--------------------------------------------------------------------*/
void cguard3l_(float *fxyz, int *nx, int *ny, int *nz, int *nxe,
               int *nye, int *nze);

/*--------------------------------------------------------------------*/
void aguard3l_(float *q, int *nx, int *ny, int *nz, int *nxe, int *nye,
               int *nze);

/*--------------------------------------------------------------------*/
void vmpois33_(float complex *q, float complex *fxyz, int *isign,
               float complex *ffc, float *ax, float *ay, float *az,
               float *affp, float *we, int *nx, int *ny, int *nz,
               int *nxvh, int *nyv, int *nzv, int *nxhd, int *nyhd,
               int *nzhd);

/*--------------------------------------------------------------------*/
void wfft3rinit_(int *mixup, float complex *sct, int *indx, int *indy,
                 int *indz, int *nxhyzd, int *nxyzhd);

/*--------------------------------------------------------------------*/
void fft3rvmxy_(float complex *f, int *isign, int *mixup,
                float complex *sct, int *indx, int *indy, int *indz,
                int *nzi, int *nzp, int *nxhd, int *nyd, int *nzd,
                int *nxhyzd, int *nxyzhd);

/*--------------------------------------------------------------------*/
void fft3rvmxz_(float complex *f, int *isign, int *mixup,
                float complex *sct, int *indx, int *indy, int *indz,
                int *nyi, int *nyp, int *nxhd, int *nyd, int *nzd,
                int *nxhyzd, int *nxyzhd);

/*--------------------------------------------------------------------*/
void fft3rvm3xy_(float complex *f, int *isign, int *mixup,
                 float complex *sct, int *indx, int *indy, int *indz,
                 int *nzi, int *nzp, int *nxhd, int *nyd, int *nzd,
                 int *nxhyzd, int *nxyzhd);

/*--------------------------------------------------------------------*/
void fft3rvm3z_(float complex *f, int *isign, int *mixup,
                float complex *sct, int *indx, int *indy, int *indz,
                int *nyi, int *nyp, int *nxhd, int *nyd, int *nzd,
                int *nxhyzd, int *nxyzhd);

/*--------------------------------------------------------------------*/
void wfft3rvmx_(float complex *f, int *isign, int *mixup,
                float complex *sct, int *indx, int *indy, int *indz,
                int *nxhd, int *nyd, int *nzd, int *nxhyzd,
                int *nxyzhd);

/*--------------------------------------------------------------------*/
void wfft3rvm3_(float complex *f, int *isign, int *mixup,
                float complex *sct, int *indx, int *indy, int *indz,
                int *nxhd, int *nyd, int *nzd, int *nxhyzd,
                int *nxyzhd);

/*--------------------------------------------------------------------*/
void set_szero3_(float *q, int *mx, int *my, int *mz, int *nxv,
                 int *nyv, int *nzv, int *mx1, int *my1, int *mxyz1);

/*--------------------------------------------------------------------*/
void set_vzero3_(float *cu, int *mx, int *my, int *mz, int *ndim,
                 int *nxv, int *nyv, int *nzv, int *mx1, int *my1, 
                 int *mxyz1);

/*--------------------------------------------------------------------*/
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
void cgppush3lt(float ppart[], float fxyz[], int kpic[], float qbm,
                float dt, float *ek, int idimp, int nppmx, int nx,
                int ny, int nz, int mx, int my, int mz, int nxv,
                int nyv, int nzv, int mx1, int my1, int mxyz1,
                int ipbc) {
   gppush3lt_(ppart,fxyz,kpic,&qbm,&dt,ek,&idimp,&nppmx,&nx,&ny,&nz,&mx,
              &my,&mz,&nxv,&nyv,&nzv,&mx1,&my1,&mxyz1,&ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cgppushf3lt(float ppart[], float fxyz[], int kpic[], int ncl[],
                 int ihole[], float qbm, float dt, float *ek, int idimp,  
                 int nppmx, int nx, int ny, int nz, int mx, int my,
                 int mz, int nxv, int nyv, int nzv, int mx1, int my1,
                 int mxyz1, int ntmax, int *irc) {
   gppushf3lt_(ppart,fxyz,kpic,ncl,ihole,&qbm,&dt,ek,&idimp,&nppmx,&nx,
               &ny,&nz,&mx,&my,&mz,&nxv,&nyv,&nzv,&mx1,&my1,&mxyz1,
               &ntmax,irc);
   return;
}

/*--------------------------------------------------------------------*/
void cvgppush3lt(float ppart[], float fxyz[], int kpic[], float qbm,
                 float dt, float *ek, int idimp, int nppmx, int nx,
                 int ny, int nz, int mx, int my, int mz, int nxv,
                 int nyv, int nzv, int mx1, int my1, int mxyz1,
                 int ipbc) {
   vgppush3lt_(ppart,fxyz,kpic,&qbm,&dt,ek,&idimp,&nppmx,&nx,&ny,&nz,
               &mx,&my,&mz,&nxv,&nyv,&nzv,&mx1,&my1,&mxyz1,&ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cvgppushf3lt(float ppart[], float fxyz[], int kpic[], int ncl[],
                  int ihole[], float qbm, float dt, float *ek,
                  int idimp, int nppmx, int nx, int ny, int nz, int mx,
                  int my, int mz, int nxv, int nyv, int nzv, int mx1,
                  int my1, int mxyz1, int ntmax, int *irc) {
   vgppushf3lt_(ppart,fxyz,kpic,ncl,ihole,&qbm,&dt,ek,&idimp,&nppmx,&nx,
                &ny,&nz,&mx,&my,&mz,&nxv,&nyv,&nzv,&mx1,&my1,&mxyz1,
                &ntmax,irc);
   return;
}

/*--------------------------------------------------------------------*/
void cv2gppush3lt(float ppart[], float fxyz[], int kpic[], float qbm,
                  float dt, float *ek, int idimp, int nppmx, int nx,
                  int ny, int nz, int mx, int my, int mz, int nxv,
                  int nyv, int nzv, int mx1, int my1, int mxyz1,
                  int ipbc) {
   v2gppush3lt_(ppart,fxyz,kpic,&qbm,&dt,ek,&idimp,&nppmx,&nx,&ny,&nz,
                &mx,&my,&mz,&nxv,&nyv,&nzv,&mx1,&my1,&mxyz1,&ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cv2gppushf3lt(float ppart[], float fxyz[], int kpic[], int ncl[],
                   int ihole[], float qbm, float dt, float *ek,
                   int idimp, int nppmx, int nx, int ny, int nz, int mx,
                   int my, int mz, int nxv, int nyv, int nzv, int mx1,
                   int my1, int mxyz1, int ntmax, int *irc) {
   v2gppushf3lt_(ppart,fxyz,kpic,ncl,ihole,&qbm,&dt,ek,&idimp,&nppmx,
                 &nx,&ny,&nz,&mx,&my,&mz,&nxv,&nyv,&nzv,&mx1,&my1,
                 &mxyz1,&ntmax,irc);
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
