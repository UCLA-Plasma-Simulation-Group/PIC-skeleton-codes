/* C Library for Skeleton 2-1/2D Electromagnetic OpenMP/Vector PIC */
/* Code                                                            */
/* Wrappers for calling the Fortran routines from a C main program */

#include <complex.h>

double ranorm_();

/*--------------------------------------------------------------------*/
void distr2h_(float *part, float *vtx, float *vty, float *vtz,
              float *vdx, float *vdy, float *vdz, int *npx, int *npy,
              int *idimp, int *nop, int *nx, int *ny, int *ipbc);

/*--------------------------------------------------------------------*/
void dblkp2l_(float *part, int *kpic, int *nppmx, int *idimp, int *nop,
              int *mx, int *my, int *mx1, int *mxy1, int *irc);

/*--------------------------------------------------------------------*/
void ppmovin2lt_(float *part, float *ppart, int *kpic, int *nppmx,
                 int *idimp, int *nop, int *mx, int *my, int *mx1,
                 int *mxy1, int *irc);

/*--------------------------------------------------------------------*/
void ppmovin2ltp_(float *part, float *ppart, int *kpic, int *kp,
                  int *nppmx, int *idimp, int *nop, int *mx, int *my,
                  int *mx1, int *mxy1, int *irc);

/*--------------------------------------------------------------------*/
void ppcheck2lt_(float *ppart, int *kpic, int *idimp, int *nppmx,
                 int *nx, int *ny, int *mx, int *my, int *mx1,
                 int *my1, int *irc);

/*--------------------------------------------------------------------*/
void gbppush23lt_(float *ppart, float *fxy, float *bxy, int *kpic,
                  float *qbm, float *dt, float *dtc, float *ek,
                  int *idimp, int *nppmx, int *nx, int *ny, int *mx,
                  int *my, int *nxv, int *nyv, int *mx1, int *mxy1,
                  int *ipbc);

/*--------------------------------------------------------------------*/
void gbppushf23lt_(float *ppart, float *fxy, float *bxy, int *kpic,
                   int *ncl, int *ihole, float *qbm, float *dt,
                   float *dtc, float *ek, int *idimp, int *nppmx,
                   int *nx, int *ny, int *mx, int *my, int *nxv,
                   int *nyv, int *mx1, int *mxy1, int *ntmax, 
                   int *irc);

/*--------------------------------------------------------------------*/
void grbppush23lt_(float *ppart, float *fxy, float *bxy, int *kpic,
                   float *qbm, float *dt, float *dtc, float *ci,
                   float *ek, int *idimp, int *nppmx, int *nx, int *ny,
                   int *mx, int *my, int *nxv, int *nyv, int *mx1,
                   int *mxy1, int *ipbc);

/*--------------------------------------------------------------------*/
void grbppushf23lt_(float *ppart, float *fxy, float *bxy, int *kpic,
                    int *ncl, int *ihole, float *qbm, float *dt,
                    float *dtc, float *ci, float *ek, int *idimp,
                    int *nppmx, int *nx, int *ny, int *mx, int *my,
                    int *nxv, int *nyv, int *mx1, int *mxy1,
                    int *ntmax, int *irc);

/*--------------------------------------------------------------------*/
void vgbppush23lt_(float *ppart, float *fxy, float *bxy, int *kpic,
                   float *qbm, float *dt, float *dtc, float *ek,
                   int *idimp, int *nppmx, int *nx, int *ny, int *mx,
                   int *my, int *nxv, int *nyv, int *mx1, int *mxy1,
                   int *ipbc);

/*--------------------------------------------------------------------*/
void vgbppushf23lt_(float *ppart, float *fxy, float *bxy, int *kpic,
                    int *ncl, int *ihole, float *qbm, float *dt,
                    float *dtc, float *ek, int *idimp, int *nppmx,
                    int *nx, int *ny, int *mx, int *my, int *nxv,
                    int *nyv, int *mx1, int *mxy1, int *ntmax, 
                    int *irc);

/*--------------------------------------------------------------------*/
void vgrbppush23lt_(float *ppart, float *fxy, float *bxy, int *kpic,
                    float *qbm, float *dt, float *dtc, float *ci,
                    float *ek, int *idimp, int *nppmx, int *nx,
                    int *ny, int *mx, int *my, int *nxv, int *nyv,
                    int *mx1, int *mxy1, int *ipbc);

/*--------------------------------------------------------------------*/
void vgrbppushf23lt_(float *ppart, float *fxy, float *bxy, int *kpic,
                     int *ncl, int *ihole, float *qbm, float *dt,
                     float *dtc, float *ci, float *ek, int *idimp,
                     int *nppmx, int *nx, int *ny, int *mx, int *my,
                     int *nxv, int *nyv, int *mx1, int *mxy1,
                     int *ntmax, int *irc);

/*--------------------------------------------------------------------*/
void gppost2lt_(float *ppart, float *q, int *kpic, float *qm,
                int *nppmx, int *idimp, int *mx, int *my, int *nxv,
                int *nyv, int *mx1, int *mxy1);

/*--------------------------------------------------------------------*/
void vgppost2lt_(float *ppart, float *q, int *kpic, float *qm,
                 int *nppmx, int *idimp, int *mx, int *my, int *nxv,
                 int *nyv, int *mx1, int *mxy1);

/*--------------------------------------------------------------------*/
void gjppost2lt_(float *ppart, float *cu, int *kpic, float *qm,
                 float *dt, int *nppmx, int *idimp, int *nx, int *ny, 
                 int *mx, int *my, int *nxv, int *nyv, int *mx1,
                 int *mxy1, int *ipbc);

/*--------------------------------------------------------------------*/
void gjppostf2lt_(float *ppart, float *cu, int *kpic, int *ncl,
                  int *ihole, float *qm, float *dt, int *nppmx, 
                  int *idimp, int *nx, int *ny, int *mx, int *my,
                  int *nxv, int *nyv, int *mx1, int *mxy1, int *ntmax,
                   int *irc);

/*--------------------------------------------------------------------*/
void grjppost2lt_(float *ppart, float *cu, int *kpic, float *qm,
                  float *dt, float *ci, int *nppmx, int *idimp,
                  int *nx, int *ny, int *mx, int *my, int *nxv,
                  int *nyv, int *mx1, int *mxy1, int *ipbc);

/*--------------------------------------------------------------------*/
void grjppostf2lt_(float *ppart, float *cu, int *kpic, int *ncl,
                   int *ihole, float *qm, float *dt, float *ci,
                   int *nppmx, int *idimp, int *nx, int *ny, int *mx,
                   int *my, int *nxv, int *nyv, int *mx1, int *mxy1,
                   int *ntmax, int *irc);

/*--------------------------------------------------------------------*/
void vgjppost2lt_(float *ppart, float *cu, int *kpic, float *qm,
                  float *dt, int *nppmx, int *idimp, int *nx, int *ny, 
                  int *mx, int *my, int *nxv, int *nyv, int *mx1,
                  int *mxy1, int *ipbc);

/*--------------------------------------------------------------------*/
void vgjppostf2lt_(float *ppart, float *cu, int *kpic, int *ncl,
                   int *ihole, float *qm, float *dt, int *nppmx, 
                   int *idimp, int *nx, int *ny, int *mx, int *my,
                   int *nxv, int *nyv, int *mx1, int *mxy1, int *ntmax,
                   int *irc);

/*--------------------------------------------------------------------*/
void vgrjppost2lt_(float *ppart, float *cu, int *kpic, float *qm,
                   float *dt, float *ci, int *nppmx, int *idimp,
                   int *nx, int *ny, int *mx, int *my, int *nxv,
                   int *nyv, int *mx1, int *mxy1, int *ipbc);

/*--------------------------------------------------------------------*/
void vgrjppostf2lt_(float *ppart, float *cu, int *kpic, int *ncl,
                    int *ihole, float *qm, float *dt, float *ci,
                    int *nppmx, int *idimp, int *nx, int *ny, int *mx,
                    int *my, int *nxv, int *nyv, int *mx1, int *mxy1,
                    int *ntmax, int *irc);

/*--------------------------------------------------------------------*/
void pporder2lt_(float *ppart, float *ppbuff, int *kpic, int *ncl,
                 int *ihole, int *idimp, int *nppmx, int *nx, int *ny,
                 int *mx, int *my, int *mx1, int *my1, int *npbmx,
                 int *ntmax, int *irc);

/*--------------------------------------------------------------------*/
void pporderf2lt_(float *ppart, float *ppbuff, int *kpic, int *ncl,
                  int *ihole, int *idimp, int *nppmx, int *mx1,
                  int *my1, int *npbmx, int *ntmax, int *irc);

/*--------------------------------------------------------------------*/
void vpporder2lt_(float *ppart, float *ppbuff, int *kpic, int *ncl,
                  int *ihole, int *idimp, int *nppmx, int *nx, int *ny,
                  int *mx, int *my, int *mx1, int *my1, int *npbmx,
                  int *ntmax, int *irc);

/*--------------------------------------------------------------------*/
void vpporderf2lt_(float *ppart, float *ppbuff, int *kpic, int *ncl,
                   int *ihole, int *idimp, int *nppmx, int *mx1,
                   int *my1, int *npbmx, int *ntmax, int *irc);

/*--------------------------------------------------------------------*/
void bguard2l_(float *bxy, int *nx, int *ny, int *nxe, int *nye);

/*--------------------------------------------------------------------*/
void acguard2l_(float *cu, int *nx, int *ny, int *nxe, int *nye);

/*--------------------------------------------------------------------*/
void aguard2l_(float *q, int *nx, int *ny, int *nxe, int *nye);

/*--------------------------------------------------------------------*/
void vmpois23_(float complex *q, float complex *fxy, int *isign,
               float complex *ffc, float *ax, float *ay, float *affp,
               float *we, int *nx, int *ny, int *nxvh, int *nyv,
               int *nxhd, int *nyhd);

/*--------------------------------------------------------------------*/
void mcuperp2_(float complex *cu, int *nx, int *ny, int *nxvh,
               int *nyv);

/*--------------------------------------------------------------------*/
void vmibpois23_(float complex *cu, float complex *bxy,
                 float complex *ffc, float *ci, float *wm, int *nx,
                 int *ny, int *nxvh, int *nyv, int *nxhd, int *nyhd);

/*--------------------------------------------------------------------*/
void vmmaxwel2_(float complex *exy, float complex *bxy,
                float complex *cu, float complex *ffc, float *ci,
                float *dt, float *wf, float *wm, int *nx, int *ny,
                int *nxvh, int *nyv, int *nxhd, int *nyhd);

/*--------------------------------------------------------------------*/
void vmemfield2_(float complex *fxy, float complex *exy,
                 float complex *ffc, int *isign, int *nx, int *ny,
                 int *nxvh, int *nyv, int *nxhd, int *nyhd);

/*--------------------------------------------------------------------*/
void wfft2rinit_(int *mixup, float complex *sct, int *indx, int *indy,
                 int *nxhyd, int *nxyhd);

/*--------------------------------------------------------------------*/
void fft2rvmxx_(float complex *f, int *isign, int *mixup,
                float complex *sct, int *indx, int *indy, int *nyi,
                int *nyp, int *nxhd, int *nyd, int *nxhyd,
                int *nxyhd);

/*--------------------------------------------------------------------*/
void fft2rmxy_(float complex *f, int *isign, int *mixup,
               float complex *sct, int *indx, int *indy, int *nxi,
               int *nxp, int *nxhd, int *nyd, int *nxhyd, int *nxyhd);

/*--------------------------------------------------------------------*/
void fft2rvm3x_(float complex *f, int *isign, int *mixup,
                float complex *sct, int *indx, int *indy, int *nyi,
                int *nyp, int *nxhd, int *nyd, int *nxhyd,
                int *nxyhd);

/*--------------------------------------------------------------------*/
void fft2rvm3y_(float complex *f, int *isign, int *mixup,
                float complex *sct, int *indx, int *indy, int *nxi,
                int *nxp, int *nxhd, int *nyd, int *nxhyd,
                int *nxyhd);

/*--------------------------------------------------------------------*/
void wfft2rvmx_(float complex *f, int *isign, int *mixup,
                float complex *sct, int *indx, int *indy, int *nxhd,
                int *nyd, int *nxhyd, int *nxyhd);

/*--------------------------------------------------------------------*/
void wfft2rvm3_(float complex *f, int *isign, int *mixup,
                float complex *sct, int *indx, int *indy, int *nxhd,
                int *nyd, int *nxhyd, int *nxyhd);

/* Interfaces to C */

double ranorm() {
  return ranorm_();
}

/*--------------------------------------------------------------------*/
void cdistr2h(float part[], float vtx, float vty, float vtz, float vdx,
              float vdy, float vdz, int npx, int npy, int idimp, int nop,
              int nx, int ny, int ipbc) {
   distr2h_(part,&vtx,&vty,&vtz,&vdx,&vdy,&vdz,&npx,&npy,&idimp,&nop,
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
void cppmovin2lt(float part[], float ppart[], int kpic[], int nppmx,
                 int idimp, int nop, int mx, int my, int mx1, int mxy1,
                 int *irc) {
   ppmovin2lt_(part,ppart,kpic,&nppmx,&idimp,&nop,&mx,&my,&mx1,&mxy1,
               irc);
   return;
}

/*--------------------------------------------------------------------*/
void cppmovin2ltp(float part[], float ppart[], int kpic[], int kp[],
                  int nppmx, int idimp, int nop, int mx, int my,
                  int mx1, int mxy1, int *irc) {
   ppmovin2ltp_(part,ppart,kpic,kp,&nppmx,&idimp,&nop,&mx,&my,&mx1,
                &mxy1,irc);
   return;
}

/*--------------------------------------------------------------------*/
void cppcheck2lt(float ppart[], int kpic[], int idimp, int nppmx, int nx,
                 int ny, int mx, int my, int mx1, int my1, 
                 int *irc) {
   ppcheck2lt_(ppart,kpic,&idimp,&nppmx,&nx,&ny,&mx,&my,&mx1,&my1,irc);
   return;
}

/*--------------------------------------------------------------------*/
void cgbppush23lt(float ppart[], float fxy[], float bxy[], int kpic[],
                  float qbm, float dt, float dtc, float *ek, int idimp,
                  int nppmx, int nx, int ny, int mx, int my, int nxv,
                  int nyv, int mx1, int mxy1, int ipbc) {
   gbppush23lt_(ppart,fxy,bxy,kpic,&qbm,&dt,&dtc,ek,&idimp,&nppmx,&nx,
                &ny,&mx,&my,&nxv,&nyv,&mx1,&mxy1,&ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cgbppushf23lt(float ppart[], float fxy[], float bxy[], int kpic[],
                   int ncl[], int ihole[], float qbm, float dt,
                   float dtc, float *ek, int idimp, int nppmx, int nx,
                   int ny, int mx, int my, int nxv, int nyv, int mx1,
                   int mxy1, int ntmax, int *irc) {
   gbppushf23lt_(ppart,fxy,bxy,kpic,ncl,ihole,&qbm,&dt,&dtc,ek,&idimp,
                 &nppmx,&nx,&ny,&mx,&my,&nxv,&nyv,&mx1,&mxy1,&ntmax,
                 irc);
   return;
}

/*--------------------------------------------------------------------*/
void cgrbppush23lt(float ppart[], float fxy[], float bxy[], int kpic[],
                   float qbm, float dt, float dtc, float ci, float *ek,
                   int idimp, int nppmx, int nx, int ny, int mx, int my,
                   int nxv, int nyv, int mx1, int mxy1, int ipbc) {
   grbppush23lt_(ppart,fxy,bxy,kpic,&qbm,&dt,&dtc,&ci,ek,&idimp,&nppmx,
                 &nx,&ny,&mx,&my,&nxv,&nyv,&mx1,&mxy1,&ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cgrbppushf23lt(float ppart[], float fxy[], float bxy[], int kpic[],
                    int ncl[], int ihole[], float qbm, float dt,
                    float dtc, float ci, float *ek, int idimp,
                    int nppmx, int nx, int ny, int mx, int my, int nxv,
                    int nyv, int mx1, int mxy1, int ntmax, int *irc) {
   grbppushf23lt_(ppart,fxy,bxy,kpic,ncl,ihole,&qbm,&dt,&dtc,&ci,ek,
                  &idimp,&nppmx,&nx,&ny,&mx,&my,&nxv,&nyv,&mx1,&mxy1,
                  &ntmax,irc);
   return;
}

/*--------------------------------------------------------------------*/
void cvgbppush23lt(float ppart[], float fxy[], float bxy[], int kpic[],
                   float qbm, float dt, float dtc, float *ek, int idimp,
                   int nppmx, int nx, int ny, int mx, int my, int nxv,
                   int nyv, int mx1, int mxy1, int ipbc) {
   vgbppush23lt_(ppart,fxy,bxy,kpic,&qbm,&dt,&dtc,ek,&idimp,&nppmx,&nx,
                 &ny,&mx,&my,&nxv,&nyv,&mx1,&mxy1,&ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cvgbppushf23lt(float ppart[], float fxy[], float bxy[], int kpic[],
                    int ncl[], int ihole[], float qbm, float dt,
                    float dtc, float *ek, int idimp, int nppmx, int nx,
                    int ny, int mx, int my, int nxv, int nyv, int mx1,
                    int mxy1, int ntmax, int *irc) {
   vgbppushf23lt_(ppart,fxy,bxy,kpic,ncl,ihole,&qbm,&dt,&dtc,ek,&idimp,
                  &nppmx,&nx,&ny,&mx,&my,&nxv,&nyv,&mx1,&mxy1,&ntmax,
                  irc);
   return;
}

/*--------------------------------------------------------------------*/
void cvgrbppush23lt(float ppart[], float fxy[], float bxy[], int kpic[],
                    float qbm, float dt, float dtc, float ci, float *ek,
                    int idimp, int nppmx, int nx, int ny, int mx,
                    int my, int nxv, int nyv, int mx1, int mxy1,
                    int ipbc) {
   vgrbppush23lt_(ppart,fxy,bxy,kpic,&qbm,&dt,&dtc,&ci,ek,&idimp,
                  &nppmx,&nx,&ny,&mx,&my,&nxv,&nyv,&mx1,&mxy1,&ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cvgrbppushf23lt(float ppart[], float fxy[], float bxy[],
                     int kpic[], int ncl[], int ihole[], float qbm,
                     float dt, float dtc, float ci, float *ek, 
                     int idimp, int nppmx, int nx, int ny, int mx,
                     int my, int nxv, int nyv, int mx1, int mxy1,
                     int ntmax, int *irc) {
   vgrbppushf23lt_(ppart,fxy,bxy,kpic,ncl,ihole,&qbm,&dt,&dtc,&ci,ek,
                   &idimp,&nppmx,&nx,&ny,&mx,&my,&nxv,&nyv,&mx1,&mxy1,
                   &ntmax,irc);
   return;
}

/*--------------------------------------------------------------------*/
void cgppost2lt(float ppart[], float q[], int kpic[], float qm,
                int nppmx, int idimp, int mx, int my, int nxv, int nyv,
                int mx1, int mxy1) {
   gppost2lt_(ppart,q,kpic,&qm,&nppmx,&idimp,&mx,&my,&nxv,&nyv,&mx1,
              &mxy1);
   return;
}

/*--------------------------------------------------------------------*/
void cvgppost2lt(float ppart[], float q[], int kpic[], float qm,
                 int nppmx, int idimp, int mx, int my, int nxv, int nyv,
                 int mx1, int mxy1) {
   vgppost2lt_(ppart,q,kpic,&qm,&nppmx,&idimp,&mx,&my,&nxv,&nyv,&mx1,
               &mxy1);
   return;
}

/*--------------------------------------------------------------------*/
void cgjppost2lt(float ppart[], float cu[], int kpic[], float qm,
                 float dt, int nppmx, int idimp, int nx, int ny, int mx,
                 int my, int nxv, int nyv, int mx1, int mxy1,
                 int ipbc) {
   gjppost2lt_(ppart,cu,kpic,&qm,&dt,&nppmx,&idimp,&nx,&ny, &mx,&my,
               &nxv,&nyv,&mx1,&mxy1,&ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cgjppostf2lt(float ppart[], float cu[], int kpic[], int ncl[],
                  int ihole[], float qm, float dt, int nppmx, int idimp,
                  int nx, int ny, int mx, int my, int nxv, int nyv,
                  int mx1, int mxy1, int ntmax, int *irc) {
   gjppostf2lt_(ppart,cu,kpic,ncl,ihole,&qm,&dt,&nppmx,&idimp,&nx,&ny,
                &mx,&my,&nxv,&nyv,&mx1,&mxy1,&ntmax,irc);
   return;
}

/*--------------------------------------------------------------------*/
void cgrjppost2lt(float ppart[], float cu[], int kpic[], float qm,
                  float dt, float ci, int nppmx, int idimp, int nx,
                  int ny, int mx, int my, int nxv, int nyv, int mx1,
                  int mxy1, int ipbc) {
   grjppost2lt_(ppart,cu,kpic,&qm,&dt,&ci,&nppmx,&idimp,&nx,&ny,&mx,
                &my,&nxv,&nyv,&mx1,&mxy1,&ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cgrjppostf2lt(float ppart[], float cu[], int kpic[], int ncl[],
                   int ihole[], float qm, float dt, float ci, int nppmx,
                   int idimp, int nx, int ny, int mx, int my, int nxv,
                   int nyv, int mx1, int mxy1, int ntmax, int *irc) {
   grjppostf2lt_(ppart,cu,kpic,ncl,ihole,&qm,&dt,&ci,&nppmx,&idimp,&nx,
                 &ny,&mx,&my,&nxv,&nyv,&mx1,&mxy1,&ntmax,irc);
   return;
}

/*--------------------------------------------------------------------*/
void cvgjppost2lt(float ppart[], float cu[], int kpic[], float qm,
                  float dt, int nppmx, int idimp, int nx, int ny,
                  int mx, int my, int nxv, int nyv, int mx1, int mxy1,
                  int ipbc) {
   vgjppost2lt_(ppart,cu,kpic,&qm,&dt,&nppmx,&idimp,&nx,&ny, &mx,&my,
                &nxv,&nyv,&mx1,&mxy1,&ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cvgjppostf2lt(float ppart[], float cu[], int kpic[], int ncl[],
                   int ihole[], float qm, float dt, int nppmx,
                   int idimp, int nx, int ny, int mx, int my, int nxv,
                   int nyv, int mx1, int mxy1, int ntmax, int *irc) {
   vgjppostf2lt_(ppart,cu,kpic,ncl,ihole,&qm,&dt,&nppmx,&idimp,&nx,&ny,
                 &mx,&my,&nxv,&nyv,&mx1,&mxy1,&ntmax,irc);
   return;
}

/*--------------------------------------------------------------------*/
void cvgrjppost2lt(float ppart[], float cu[], int kpic[], float qm,
                   float dt, float ci, int nppmx, int idimp, int nx,
                   int ny, int mx, int my, int nxv, int nyv, int mx1,
                   int mxy1, int ipbc) {
   vgrjppost2lt_(ppart,cu,kpic,&qm,&dt,&ci,&nppmx,&idimp,&nx,&ny,&mx,
                 &my,&nxv,&nyv,&mx1,&mxy1,&ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cvgrjppostf2lt(float ppart[], float cu[], int kpic[], int ncl[],
                    int ihole[], float qm, float dt, float ci,
                    int nppmx, int idimp, int nx, int ny, int mx,
                    int my, int nxv, int nyv, int mx1, int mxy1,
                    int ntmax, int *irc) {
   vgrjppostf2lt_(ppart,cu,kpic,ncl,ihole,&qm,&dt,&ci,&nppmx,&idimp,
                  &nx,&ny,&mx,&my,&nxv,&nyv,&mx1,&mxy1,&ntmax,irc);
   return;
}

/*--------------------------------------------------------------------*/
void cpporder2lt(float ppart[], float ppbuff[], int kpic[], int ncl[],
                 int ihole[], int idimp, int nppmx, int nx, int ny,
                 int mx, int my, int mx1, int my1, int npbmx, int ntmax,
                 int *irc) {
   pporder2lt_(ppart,ppbuff,kpic,ncl,ihole,&idimp,&nppmx,&nx,&ny,&mx,
               &my,&mx1,&my1,&npbmx,&ntmax,irc);
   return;
}

/*--------------------------------------------------------------------*/
void cpporderf2lt(float ppart[], float ppbuff[], int kpic[], int ncl[],
                  int ihole[], int idimp, int nppmx, int mx1, int my1,
                  int npbmx, int ntmax, int *irc) {
   pporderf2lt_(ppart,ppbuff,kpic,ncl,ihole,&idimp,&nppmx,&mx1,&my1,
                &npbmx,&ntmax,irc);
   return;
}

/*--------------------------------------------------------------------*/
void cvpporder2lt(float ppart[], float ppbuff[], int kpic[], int ncl[],
                  int ihole[], int idimp, int nppmx, int nx, int ny,
                  int mx, int my, int mx1, int my1, int npbmx,
                  int ntmax, int *irc) {
   vpporder2lt_(ppart,ppbuff,kpic,ncl,ihole,&idimp,&nppmx,&nx,&ny,&mx,
                &my,&mx1,&my1,&npbmx,&ntmax,irc);
   return;
}

/*--------------------------------------------------------------------*/
void cvpporderf2lt(float ppart[], float ppbuff[], int kpic[], int ncl[],
                   int ihole[], int idimp, int nppmx, int mx1, int my1,
                   int npbmx, int ntmax, int *irc) {
   vpporderf2lt_(ppart,ppbuff,kpic,ncl,ihole,&idimp,&nppmx,&mx1,&my1,
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
void cvmpois23(float complex q[], float complex fxy[], int isign,
               float complex ffc[], float ax, float ay, float affp,
               float *we, int nx, int ny, int nxvh, int nyv, int nxhd,
               int nyhd) {
   vmpois23_(q,fxy,&isign,ffc,&ax,&ay,&affp,we,&nx,&ny,&nxvh,&nyv,
             &nxhd,&nyhd);
   return;
}

/*--------------------------------------------------------------------*/
void cmcuperp2(float complex cu[], int nx, int ny, int nxvh, int nyv) {
   mcuperp2_(cu,&nx,&ny,&nxvh,&nyv);
   return;
}

/*--------------------------------------------------------------------*/
void cvmibpois23(float complex cu[], float complex bxy[],
                 float complex ffc[], float ci, float *wm, int nx,
                 int ny, int nxvh, int nyv, int nxhd, int nyhd) {
   vmibpois23_(cu,bxy,ffc,&ci,wm,&nx,&ny,&nxvh,&nyv,&nxhd,&nyhd);
   return;
}

/*--------------------------------------------------------------------*/
void cvmmaxwel2(float complex exy[], float complex bxy[],
                float complex cu[], float complex ffc[], float ci,
                float dt, float *wf, float *wm, int nx, int ny,
                int nxvh, int nyv, int nxhd, int nyhd) {
   vmmaxwel2_(exy,bxy,cu,ffc,&ci,&dt,wf,wm,&nx,&ny,&nxvh,&nyv,&nxhd,
              &nyhd);
   return;
}

/*--------------------------------------------------------------------*/
void cvmemfield2(float complex fxy[], float complex exy[],
                 float complex ffc[], int isign, int nx, int ny,
                 int nxvh, int nyv, int nxhd, int nyhd) {
   vmemfield2_(fxy,exy,ffc,&isign,&nx,&ny,&nxvh,&nyv,&nxhd,&nyhd);
   return;
}

/*--------------------------------------------------------------------*/
void cwfft2rinit(int mixup[], float complex sct[], int indx, int indy, 
                 int nxhyd, int nxyhd) {
   wfft2rinit_(mixup,sct,&indx,&indy,&nxhyd,&nxyhd);
   return;
}

/*--------------------------------------------------------------------*/
void cfft2rvmxx(float complex f[], int isign, int mixup[],
                float complex sct[], int indx, int indy, int nyi,
                int nyp, int nxhd, int nyd, int nxhyd, int nxyhd) {
   fft2rvmxx_(f,&isign,mixup,sct,&indx,&indy,&nyi,&nyp,&nxhd,&nyd,
              &nxhyd,&nxyhd);
   return;
}

/*--------------------------------------------------------------------*/
void cfft2rmxy(float complex f[], int isign, int mixup[],
               float complex sct[], int indx, int indy, int nxi,
               int nxp, int nxhd, int nyd, int nxhyd, int nxyhd) {
   fft2rmxy_(f,&isign,mixup,sct,&indx,&indy,&nxi,&nxp,&nxhd,&nyd,
             &nxhyd,&nxyhd);
   return;
}

/*--------------------------------------------------------------------*/
void cfft2rvm3x(float complex f[], int isign, int mixup[],
                float complex sct[], int indx, int indy, int nyi,
                int nyp, int nxhd, int nyd, int nxhyd, int nxyhd) {
   fft2rvm3x_(f,&isign,mixup,sct,&indx,&indy,&nyi,&nyp,&nxhd,&nyd,
              &nxhyd,&nxyhd);
   return;
}

/*--------------------------------------------------------------------*/
void cfft2rvm3y(float complex f[], int isign, int mixup[],
                float complex sct[], int indx, int indy, int nxi,
                int nxp, int nxhd, int nyd, int nxhyd, int nxyhd) {
   fft2rvm3y_(f,&isign,mixup,sct,&indx,&indy,&nxi,&nxp,&nxhd,&nyd,
              &nxhyd,&nxyhd);
   return;
}

/*--------------------------------------------------------------------*/
void cwfft2rvmx(float complex f[], int isign, int mixup[],
                float complex sct[], int indx, int indy, int nxhd,
                int nyd, int nxhyd, int nxyhd) {
   wfft2rvmx_(f,&isign,mixup,sct,&indx,&indy,&nxhd,&nyd,&nxhyd,&nxyhd);
   return;
}

/*--------------------------------------------------------------------*/
void cwfft2rvm3(float complex f[], int isign, int mixup[],
                float complex sct[], int indx, int indy, int nxhd,
                int nyd, int nxhyd, int nxyhd) {
   wfft2rvm3_(f,&isign,mixup,sct,&indx,&indy,&nxhd,&nyd,&nxhyd,&nxyhd);
   return;
}
