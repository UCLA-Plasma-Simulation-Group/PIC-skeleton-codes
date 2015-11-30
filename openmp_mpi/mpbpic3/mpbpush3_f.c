/* C Library for Skeleton 3D Electromagnetic MPI/OpenMP PIC Code */
/* Wrappers for calling the Fortran routines from a C main program */

#include <complex.h>
#include "mpbpush3.h"

double ranorm_();

void pdicomp32l_(float *edges, int *nyzp, int *noff, int *nypmx,
                 int *nzpmx, int *nypmn, int *nzpmn, int *ny, int *nz,
                 int *kstrt, int *nvpy, int *nvpz, int *idps,
                 int *idds);

void fcomp32_(int *nvp, int *nx, int *ny, int *nz, int *nvpy, int *nvpz,
              int *ierr);

void pdistr32_(float *part, float *edges, int *npp, int *nps,
               float *vtx, float *vty, float *vtz, float *vdx,
               float *vdy, float *vdz, int *npx, int *npy, int *npz,
               int *nx, int *ny, int *nz, int *idimp, int *npmax,
               int *idps, int *ipbc, int *ierr);

void ppdblkp3l_(float *part, int *kpic, int *npp, int *noff, int *nppmx,
                int *idimp, int *npmax, int *mx, int *my, int *mz,
                int *mx1, int *myp1, int *mxyzp1, int *idds, int *irc);

void pppmovin3l_(float *part, float *ppart, int *kpic, int *npp,
                 int *noff, int *nppmx, int *idimp, int *npmax, int *mx,
                 int *my, int *mz, int *mx1, int *myp1, int *mxyzp1,
                 int *idds, int *irc);

void pppcheck3l_(float *ppart, int *kpic, int *noff, int *nyzp,
                 int *idimp, int *nppmx, int *nx, int *mx, int *my,
                 int *mz, int *mx1, int *myp1, int *mzp1, int *idds,
                 int *irc);

void ppgbppush32l_(float *ppart, float *fxyz, float *bxyz, int *kpic,
                   int *noff, int *nyzp, float *qbm, float *dt,
                   float *dtc, float *ek, int *idimp, int *nppmx,
                   int *nx, int *ny, int *nz, int *mx, int *my, int *mz,
                   int *nxv, int *nypmx, int *nzpmx, int *mx1,
                   int *myp1, int *mxyzp1, int *idds, int *ipbc);

void ppgbppushf32l_(float *ppart, float *fxyz, float *bxyz, int *kpic,
                    int *ncl, int *ihole, int *noff, int *nyzp,
                    float *qbm, float *dt, float *dtc, float *ek,
                    int *idimp, int *nppmx, int *nx, int *ny, int *nz,
                    int *mx, int *my, int *mz, int *nxv, int *nypmx,
                    int *nzpmx, int *mx1, int *myp1, int *mxyzp1, 
                    int *ntmax, int *idds, int *irc);

void ppgrbppush32l_(float *ppart, float *fxyz, float *bxyz, int *kpic,
                    int *noff, int *nyzp, float *qbm, float *dt, 
                    float *dtc, float *ci, float *ek, int *idimp,
                    int *nppmx, int *nx, int *ny, int *nz, int *mx,
                    int *my, int *mz, int *nxv, int *nypmx, int *nzpmx,
                    int *mx1, int *myp1, int *mxyzp1, int *idds,
                    int *ipbc);

void ppgrbppushf32l_(float *ppart, float *fxyz, float *bxyz, int *kpic,
                     int *ncl, int *ihole, int *noff, int *nyzp,
                     float *qbm, float *dt, float *dtc, float *ci,
                     float *ek, int *idimp, int *nppmx, int *nx,
                     int *ny, int *nz, int *mx, int *my, int *mz,
                     int *nxv, int *nypmx, int *nzpmx, int *mx1,
                     int *myp1, int *mxyzp1, int *ntmax, int *idds,
                     int *irc);

void ppgppost32l_(float *ppart, float *q, int *kpic, int *noff,
                  float *qm, int *nppmx, int *idimp, int *mx, int *my,
                  int *mz, int *nxv, int *nypmx, int *nzpmx, int *mx1,
                  int *myp1, int *mxyzp1, int *idds);

void ppgjppost32l_(float *ppart, float *cu, int *kpic, int *noff,
                   float *qm, float *dt, int *nppmx, int *idimp,
                   int *nx, int *ny, int *nz, int *mx, int *my, int *mz,
                   int *nxv, int *nypmx, int *nzpmx, int *mx1, 
                   int *myp1, int *mxyzp1, int *idds, int *ipbc);

void ppgjppostf32l_(float *ppart, float *cu, int *kpic, int *ncl,
                    int *ihole, int *noff, int *nyzp, float *qm,
                    float *dt, int *nppmx, int *idimp, int *nx, int *ny,
                    int *nz, int *mx, int *my, int *mz, int *nxv,
                    int *nypmx, int *nzpmx, int *mx1, int *myp1,
                    int *mxyzp1, int *ntmax, int *idds, int *irc);

void ppgrjppost32l_(float *ppart, float *cu, int *kpic, int *noff,
                    float *qm, float *dt, float *ci, int *nppmx,
                    int *idimp, int *nx, int *ny, int *nz, int *mx,
                    int *my, int *mz, int *nxv, int *nypmx, int *nzpmx,
                    int *mx1, int *myp1, int *mxyzp1, int *idds,
                    int *ipbc);

void ppgrjppostf32l_(float *ppart, float *cu, int *kpic, int *ncl,
                     int *ihole, int *noff, int *nyzp, float *qm,
                     float *dt, float *ci, int *nppmx, int *idimp,
                     int *nx, int *ny, int *nz, int *mx, int *my,
                     int *mz, int *nxv, int *nypmx, int *nzpmx,
                     int *mx1, int *myp1, int *mxyzp1, int *ntmax,
                     int *idds, int *irc);

void ppporder32la_(float *ppart, float *ppbuff, float *sbufl,
                   float *sbufr, int *kpic, int *ncl, int *ihole,
                   int *ncll, int *nclr, int *noff, int *nyzp,
                   int *idimp, int *nppmx, int *nx, int *ny, int *nz,
                   int *mx, int *my, int *mz, int *mx1, int *myp1,
                   int *mzp1, int *mxzyp1, int *npbmx, int *ntmax,
                   int *nbmax, int *idds, int *irc);

void ppporderf32la_(float *ppart, float *ppbuff, float *sbufl,
                    float *sbufr, int *ncl, int *ihole, int *ncll,
                    int *nclr, int *idimp, int *nppmx, int *mx1,
                    int *myp1, int *mzp1, int *mxzyp1, int *npbmx,
                    int *ntmax, int *nbmax, int *irc);

void ppporder32lb_(float *ppart, float *ppbuff, float *rbufl,
                   float *rbufr, int *kpic, int *ncl, int *ihole,
                   int *mcll, int *mclr, int *mcls, int *idimp,
                   int *nppmx, int *mx1, int *myp1, int *mzp1,
                   int *mxzyp1, int *npbmx, int *ntmax, int *nbmax,
                   int *irc);

void ppcguard32xl_(float *fxyz, int *nyzp, int *nx, int *ndim, int *nxe,
                   int *nypmx, int *nzpmx, int *idds);

void ppaguard32xl_(float *q, int *nyzp, int *nx, int *nxe, int *nypmx,
                   int *nzpmx, int *idds);

void ppacguard32xl_(float *cu, int *nyzp, int *nx, int *ndim, int *nxe,
                    int *nypmx, int *nzpmx, int *idds);

void mppois332_(float complex *q, float complex *fxyz, int *isign,
                float complex *ffc, float *ax, float *ay, float *az,
                float *affp, float *we, int *nx, int *ny, int *nz,
                int *kstrt, int *nvpy, int *nvpz, int *nzv, int *kxyp,
                int *kyzp, int *nzhd);

void mppcuperp32_(float complex *cu, int *nx, int *ny, int *nz,
                  int *kstrt, int *nvpy, int *nvpz, int *nzv, int *kxyp,
                  int *kyzp);

void mippbpoisp332_(float complex *cu, float complex *bxyz,
                    float complex *ffc, float *ci, float *wm, int *nx,
                    int *ny, int *nz, int *kstrt, int *nvpy, int *nvpz,
                    int *nzv, int *kxyp, int *kyzp, int *nzhd);

void mppmaxwel32_(float complex *exyz, float complex *bxyz,
                  float complex *cu, float complex *ffc, float *affp,
                  float *ci, float *dt, float *wf, float *wm, int *nx,
                  int *ny, int *nz, int *kstrt, int *nvpy, int *nvpz,
                  int *nzv, int *kxyp, int *kyzp, int *nzhd);

void mppemfield32_(float complex *fxyz, float complex *exyz,
                   float complex *ffc, int *isign, int *nx, int *ny,
                   int *nz, int *kstrt, int *nvpy, int *nvpz, int *nzv,
                   int *kxyp, int *kyzp, int *nzhd);

void wpfft32rinit_(int *mixup, float complex *sct, int *indx, int *indy,
                   int *indz, int *nxhyzd, int *nxyzhd);

void ppfft32rmxx_(float complex *f, int *isign, int *mixup,
                  float complex *sct, int *indx, int *indy, int *indz,
                  int *kstrt, int *nvp, int *kypi, int *kypp, int *nxvh,
                  int *kzpp, int *kypd, int *kzpd, int *nxhyzd,
                  int *nxyzhd);

void ppfft32rmxy_(float complex *g, int *isign, int *mixup,
                  float complex *sct, int *indx, int *indy, int *indz,
                  int *kstrt, int *nvpy, int *nvpz, int *kxypi,
                  int *kxypp, int *nyv, int *kzpp, int *kxypd, 
                  int *kzpd, int *nxhyzd, int *nxyzhd);

void ppfft32rmxz_(float complex *h, int *isign, int *mixup,
                  float complex *sct, int *indx, int *indy, int *indz,
                  int *kstrt, int *nvpy, int *nvpz, int *kxypi,
                  int *kxypp, int *nzv, int *kyzp, int *kxypd,
                  int *kyzpd, int *nxhyzd, int *nxyzhd);

void ppfft32rm3xx_(float complex *f, int *isign, int *mixup,
                   float complex *sct, int *indx, int *indy, int *indz,
                   int *kstrt, int *nvp, int *kypi, int *kypp,
                   int *nxvh, int *kzpp, int *kypd, int *kzpd,
                   int *nxhyzd, int *nxyzhd);

void ppfft32rm3xy_(float complex *g, int *isign, int *mixup,
                   float complex *sct, int *indx, int *indy, int *indz,
                   int *kstrt, int *nvpy, int *nvpz, int *kxypi,
                   int *kxypp, int *nyv, int *kzpp, int *kxypd, 
                   int *kzpd, int *nxhyzd, int *nxyzhd);

void ppfft32rm3xz_(float complex *h, int *isign, int *mixup,
                   float complex *sct, int *indx, int *indy, int *indz,
                   int *kstrt, int *nvpy, int *nvpz, int *kxypi,
                   int *kxypp, int *nzv, int *kyzp, int *kxypd, 
                   int *kyzpd, int *nxhyzd, int *nxyzhd);

void wppfft32rm_(float complex *f, float complex *g, float complex *h,
                 float complex *bs, float complex *br, int *isign,
                 int *ntpose, int *mixup, float complex *sct,
                 float *ttp, int *indx, int *indy, int *indz,
                 int *kstrt, int *nvpy, int *nvpz, int *nxvh, int *nyv,
                 int *nzv, int *kxyp, int *kyp, int *kyzp, int *kzp,
                 int *kxypd, int *kypd, int *kyzpd, int *kzpd,
                 int *kzyp, int *nxhyzd, int *nxyzhd);

void wppfft32rm3_(float complex *f, float complex *g, float complex *h,
                  float complex *bs, float complex *br, int *isign,
                  int *ntpose, int *mixup, float complex *sct,
                  float *ttp, int *indx, int *indy, int *indz,
                  int *kstrt, int *nvpy, int *nvpz, int *nxvh, int *nyv,
                  int *nzv, int *kxyp, int *kyp, int *kyzp, int *kzp,
                  int *kxypd, int *kypd, int *kyzpd, int *kzpd,
                  int *kzyp, int *nxhyzd, int *nxyzhd);

/* Interfaces to C */

double ranorm() {
  return ranorm_();
}

/*--------------------------------------------------------------------*/
void cpdicomp32l(float edges[], int nyzp[], int noff[], int *nypmx,
                 int *nzpmx, int *nypmn, int *nzpmn, int ny, int nz,
                 int kstrt, int nvpy, int nvpz, int idps, int idds) {
   pdicomp32l_(edges,nyzp,noff,nypmx,nzpmx,nypmn,nzpmn,&ny,&nz,&kstrt,
               &nvpy,&nvpz,&idps,&idds);
   return;
}

/*--------------------------------------------------------------------*/
void cfcomp32(int nvp, int nx, int ny, int nz, int *nvpy, int *nvpz,
              int *ierr) {
   fcomp32_(&nvp,&nx,&ny,&nz,nvpy,nvpz,ierr);
   return;
}

/*--------------------------------------------------------------------*/
void cpdistr32(float part[], float edges[], int *npp, int nps,
               float vtx, float vty, float vtz, float vdx, float vdy,
               float vdz, int npx, int npy, int npz, int nx, int ny,
               int nz, int idimp, int npmax, int idps, int ipbc,
               int *ierr) {
   pdistr32_(part,edges,npp,&nps,&vtx,&vty,&vtz,&vdx,&vdy,&vdz,&npx,
             &npy,&npz,&nx,&ny,&nz,&idimp,&npmax,&idps,&ipbc,ierr);
   return;
}

/*--------------------------------------------------------------------*/
void cppdblkp3l(float part[], int kpic[], int npp, int noff[],
                int *nppmx, int idimp, int npmax, int mx, int my,
                int mz, int mx1, int myp1, int mxyzp1, int idds,
                int *irc) {
   ppdblkp3l_(part,kpic,&npp,noff,nppmx,&idimp,&npmax,&mx,&my,&mz,&mx1,
              &myp1,&mxyzp1,&idds,irc);
   return;
}

/*--------------------------------------------------------------------*/
void cpppmovin3l(float part[], float ppart[], int kpic[], int npp,
                 int noff[], int nppmx, int idimp, int npmax, int mx,
                 int my, int mz, int mx1, int myp1, int mxyzp1,
                 int idds, int *irc) {
   pppmovin3l_(part,ppart,kpic,&npp,noff,&nppmx,&idimp,&npmax,&mx,&my,
               &mz,&mx1,&myp1,&mxyzp1,&idds,irc);
   return;
}

/*--------------------------------------------------------------------*/
void cpppcheck3l(float ppart[], int kpic[], int noff[], int nyzp[],
                 int idimp, int nppmx, int nx, int mx, int my, int mz,
                 int mx1, int myp1, int mzp1, int idds, int *irc) {
   pppcheck3l_(ppart,kpic,noff,nyzp,&idimp,&nppmx,&nx,&mx,&my,&mz,&mx1,
               &myp1,&mzp1,&idds,irc);
   return;
}

/*--------------------------------------------------------------------*/
void cppgbppush32l(float ppart[], float fxyz[], float bxyz[],
                   int kpic[], int noff[], int nyzp[], float qbm,
                   float dt, float dtc, float *ek, int idimp, int nppmx,
                   int nx, int ny, int nz, int mx, int my, int mz, 
                   int nxv, int nypmx, int nzpmx, int mx1, int myp1,
                   int mxyzp1, int idds, int ipbc) {
   ppgbppush32l_(ppart,fxyz,bxyz,kpic,noff,nyzp,&qbm,&dt,&dtc,ek,&idimp,
                 &nppmx,&nx,&ny,&nz,&mx,&my,&mz,&nxv,&nypmx,&nzpmx,&mx1,
                 &myp1,&mxyzp1,&idds,&ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cppgbppushf32l(float ppart[], float fxyz[], float bxyz[],
                    int kpic[], int ncl[], int ihole[], int noff[],
                    int nyzp[], float qbm, float dt, float dtc, 
                    float *ek, int idimp, int nppmx, int nx, int ny,
                    int nz, int mx, int my, int mz, int nxv, int nypmx,
                    int nzpmx, int mx1, int myp1, int mxyzp1, int ntmax,
                    int idds, int *irc) {
   ppgbppushf32l_(ppart,fxyz,bxyz,kpic,ncl,ihole,noff,nyzp,&qbm,&dt,
                  &dtc,ek,&idimp,&nppmx,&nx,&ny,&nz,&mx,&my,&mz,&nxv,
                  &nypmx,&nzpmx,&mx1,&myp1,&mxyzp1, &ntmax,&idds,irc);
   return;
}

/*--------------------------------------------------------------------*/
void cppgrbppush32l(float ppart[], float fxyz[], float bxyz[],
                    int kpic[], int noff[], int nyzp[], float qbm,
                    float dt, float dtc, float ci, float *ek, int idimp,
                    int nppmx, int nx, int ny, int nz, int mx, int my,
                    int mz, int nxv, int nypmx, int nzpmx, int mx1,
                    int myp1, int mxyzp1, int idds, int ipbc) {
   ppgrbppush32l_(ppart,fxyz,bxyz,kpic,noff,nyzp,&qbm,&dt,&dtc,&ci,ek,
                  &idimp,&nppmx,&nx,&ny,&nz,&mx,&my,&mz,&nxv,&nypmx,
                  &nzpmx,&mx1,&myp1,&mxyzp1,&idds,&ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cppgrbppushf32l(float ppart[], float fxyz[], float bxyz[],
                     int kpic[], int ncl[], int ihole[], int noff[],
                     int nyzp[], float qbm, float dt, float dtc,
                     float ci, float *ek, int idimp, int nppmx, int nx,
                     int ny, int nz, int mx, int my, int mz, int nxv,
                     int nypmx, int nzpmx, int mx1, int myp1,
                     int mxyzp1, int ntmax, int idds, int *irc) {
   ppgrbppushf32l_(ppart,fxyz,bxyz,kpic,ncl,ihole,noff,nyzp,&qbm,&dt,
                   &dtc,&ci,ek,&idimp,&nppmx,&nx,&ny,&nz,&mx,&my,&mz,
                   &nxv,&nypmx,&nzpmx,&mx1,&myp1,&mxyzp1,&ntmax,&idds,
                   irc);
   return;
}

/*--------------------------------------------------------------------*/
void cppgppost32l(float ppart[], float q[], int kpic[], int noff[],
                  float qm, int nppmx, int idimp, int mx, int my,
                  int mz, int nxv, int nypmx, int nzpmx, int mx1,
                  int myp1, int mxyzp1, int idds) {
   ppgppost32l_(ppart,q,kpic,noff,&qm,&nppmx,&idimp,&mx,&my,&mz,&nxv,
                &nypmx,&nzpmx,&mx1,&myp1,&mxyzp1,&idds);
   return;
}

/*--------------------------------------------------------------------*/
void cppgjppost32l(float ppart[], float cu[], int kpic[], int noff[],
                   float qm, float dt, int nppmx, int idimp, int nx,
                   int ny, int nz, int mx, int my, int mz, int nxv,
                   int nypmx, int nzpmx, int mx1, int myp1, int mxyzp1,
                   int idds, int ipbc) {
   ppgjppost32l_(ppart,cu,kpic,noff,&qm,&dt,&nppmx,&idimp,&nx,&ny,&nz,
                 &mx,&my,&mz,&nxv,&nypmx,&nzpmx,&mx1,&myp1,&mxyzp1,
                 &idds,&ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cppgjppostf32l(float ppart[], float cu[], int kpic[], int ncl[],
                    int ihole[], int noff[], int nyzp[], float qm,
                    float dt, int nppmx, int idimp, int nx, int ny,
                    int nz, int mx, int my, int mz, int nxv, int nypmx,
                    int nzpmx, int mx1, int myp1, int mxyzp1, int ntmax,
                    int idds, int *irc) {
   ppgjppostf32l_(ppart,cu,kpic,ncl,ihole,noff,nyzp,&qm,&dt,&nppmx,
                  &idimp,&nx,&ny,&nz,&mx,&my,&mz,&nxv,&nypmx,&nzpmx,
                  &mx1,&myp1,&mxyzp1,&ntmax,&idds,irc);
   return;
}

/*--------------------------------------------------------------------*/
void cppgrjppost32l(float ppart[], float cu[], int kpic[], int noff[],
                    float qm, float dt, float ci, int nppmx, int idimp,
                    int nx, int ny, int nz, int mx, int my, int mz,
                    int nxv, int nypmx, int nzpmx, int mx1, int myp1,
                    int mxyzp1, int idds, int ipbc) {
   ppgrjppost32l_(ppart,cu,kpic,noff,&qm,&dt,&ci,&nppmx,&idimp,&nx,&ny,
                  &nz,&mx,&my,&mz,&nxv,&nypmx,&nzpmx,&mx1,&myp1,&mxyzp1,
                  &idds,&ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cppgrjppostf32l(float ppart[], float cu[], int kpic[], int ncl[],
                     int ihole[], int noff[], int nyzp[], float qm,
                     float dt, float ci, int nppmx, int idimp, int nx,
                     int ny, int nz, int mx, int my, int mz, int nxv,
                     int nypmx, int nzpmx, int mx1, int myp1,
                     int mxyzp1, int ntmax, int idds, int *irc) {
   ppgrjppostf32l_(ppart,cu,kpic,ncl,ihole,noff,nyzp,&qm,&dt,&ci,&nppmx,
                   &idimp,&nx,&ny,&nz,&mx,&my,&mz,&nxv,&nypmx,&nzpmx,
                   &mx1,&myp1,&mxyzp1,&ntmax,&idds,irc);
   return;
}

/*--------------------------------------------------------------------*/
void cppporder32la(float ppart[], float ppbuff[], float sbufl[],
                   float sbufr[], int kpic[], int ncl[], int ihole[],
                   int ncll[], int nclr[], int noff[], int nyzp[],
                   int idimp, int nppmx, int nx, int ny, int nz, int mx,
                   int my, int mz, int mx1, int myp1, int mzp1,
                   int mxzyp1, int npbmx, int ntmax, int nbmax,
                   int idds, int *irc) {
   ppporder32la_(ppart,ppbuff,sbufl,sbufr,kpic,ncl,ihole,ncll,nclr,noff,
                 nyzp,&idimp,&nppmx,&nx,&ny,&nz,&mx,&my,&mz,&mx1,&myp1,
                 &mzp1,&mxzyp1,&npbmx,&ntmax,&nbmax,&idds,irc);
   return;
}

/*--------------------------------------------------------------------*/
void cppporderf32la(float ppart[], float ppbuff[], float sbufl[],
                    float sbufr[], int ncl[], int ihole[], int ncll[],
                    int nclr[], int idimp, int nppmx, int mx1, int myp1,
                    int mzp1, int mxzyp1, int npbmx, int ntmax,
                    int nbmax, int *irc) {
   ppporderf32la_(ppart,ppbuff,sbufl,sbufr,ncl,ihole,ncll,nclr,&idimp,
                  &nppmx,&mx1,&myp1,&mzp1,&mxzyp1,&npbmx,&ntmax,&nbmax,
                  irc);
   return;
}

/*--------------------------------------------------------------------*/
void cppporder32lb(float ppart[], float ppbuff[], float rbufl[],
                   float rbufr[], int kpic[], int ncl[], int ihole[],
                   int mcll[], int mclr[], int mcls[], int idimp,
                   int nppmx, int mx1, int myp1, int mzp1, int mxzyp1,
                   int npbmx, int ntmax, int nbmax, int *irc) {
   ppporder32lb_(ppart,ppbuff,rbufl,rbufr,kpic,ncl,ihole,mcll,mclr,mcls,
                 &idimp,&nppmx,&mx1,&myp1,&mzp1,&mxzyp1,&npbmx,&ntmax,
                 &nbmax,irc);
   return;
}

/*--------------------------------------------------------------------*/
void cppcguard32xl(float fxyz[], int nyzp[], int nx, int ndim, int nxe,
                   int nypmx, int nzpmx, int idds) {
   ppcguard32xl_(fxyz,nyzp,&nx,&ndim,&nxe,&nypmx,&nzpmx,&idds);
   return;
}

/*--------------------------------------------------------------------*/
void cppaguard32xl(float q[], int nyzp[], int nx, int nxe, int nypmx,
                   int nzpmx, int idds) {
   ppaguard32xl_(q,nyzp,&nx,&nxe,&nypmx,&nzpmx,&idds);
   return;
}

/*--------------------------------------------------------------------*/
void cppacguard32xl(float cu[], int nyzp[], int nx, int ndim, int nxe,
                    int nypmx, int nzpmx, int idds) {
   ppacguard32xl_(cu,nyzp,&nx,&ndim,&nxe,&nypmx,&nzpmx,&idds);
   return;
}

/*--------------------------------------------------------------------*/
void cmppois332(float complex q[], float complex fxyz[], int isign,
                float complex ffc[], float ax, float ay, float az,
                float affp, float *we, int nx, int ny, int nz,
                int kstrt, int nvpy, int nvpz, int nzv, int kxyp,
                int kyzp, int nzhd) {
   mppois332_(q,fxyz,&isign,ffc,&ax,&ay,&az,&affp,we,&nx,&ny,&nz,&kstrt,
              &nvpy,&nvpz,&nzv,&kxyp,&kyzp,&nzhd);
   return;
}

/*--------------------------------------------------------------------*/
void cmppcuperp32(float complex cu[], int nx, int ny, int nz, int kstrt,
                  int nvpy, int nvpz, int nzv, int kxyp, int kyzp) {
   mppcuperp32_(cu,&nx,&ny,&nz,&kstrt,&nvpy,&nvpz,&nzv,&kxyp,&kyzp);
   return;
}

/*--------------------------------------------------------------------*/
void cmippbpoisp332(float complex cu[], float complex bxyz[],
                    float complex ffc[], float ci, float *wm, int nx,
                    int ny, int nz, int kstrt, int nvpy, int nvpz,
                    int nzv, int kxyp, int kyzp, int nzhd) {
   mippbpoisp332_(cu,bxyz,ffc,&ci,wm,&nx,&ny,&nz,&kstrt,&nvpy,&nvpz,&nzv,
                  &kxyp,&kyzp,&nzhd);
   return;
}

/*--------------------------------------------------------------------*/
void cmppmaxwel32(float complex exyz[], float complex bxyz[],
                  float complex cu[], float complex ffc[], float affp,
                  float ci, float dt, float *wf, float *wm, int nx,
                  int ny, int nz, int kstrt, int nvpy, int nvpz,
                  int nzv, int kxyp, int kyzp, int nzhd) {
   mppmaxwel32_(exyz,bxyz,cu,ffc,&affp,&ci,&dt,wf,wm,&nx,&ny,&nz,&kstrt,
                &nvpy,&nvpz,&nzv,&kxyp,&kyzp,&nzhd);
   return;
}

/*--------------------------------------------------------------------*/
void cmppemfield32(float complex fxyz[], float complex exyz[],
                   float complex ffc[], int isign, int nx, int ny,
                   int nz, int kstrt, int nvpy, int nvpz, int nzv,
                   int kxyp, int kyzp, int nzhd) {
   mppemfield32_(fxyz,exyz,ffc,&isign,&nx,&ny,&nz,&kstrt,&nvpy,&nvpz,
                 &nzv,&kxyp,&kyzp,&nzhd);
   return;
}

/*--------------------------------------------------------------------*/
void cwpfft32rinit(int mixup[], float complex sct[], int indx, int indy,
                   int indz, int nxhyzd, int nxyzhd) {
   wpfft32rinit_(mixup,sct,&indx,&indy,&indz,&nxhyzd,&nxyzhd);
   return;
}

/*--------------------------------------------------------------------*/
void cppfft32rmxx(float complex f[], int isign, int mixup[],
                  float complex sct[], int indx, int indy, int indz,
                  int kstrt, int nvp, int kypi, int kypp, int nxvh,
                  int kzpp, int kypd, int kzpd, int nxhyzd, int nxyzhd) {
   ppfft32rmxx_(f,&isign,mixup,sct,&indx,&indy,&indz,&kstrt,&nvp,&kypi,
                &kypp,&nxvh,&kzpp,&kypd,&kzpd,&nxhyzd,&nxyzhd);
   return;
}

/*--------------------------------------------------------------------*/
void cppfft32rmxy(float complex g[], int isign, int mixup[],
                  float complex sct[], int indx, int indy, int indz,
                  int kstrt, int nvpy, int nvpz, int kxypi, int kxypp,
                  int nyv, int kzpp, int kxypd, int kzpd, int nxhyzd,
                  int nxyzhd) {
   ppfft32rmxy_(g,&isign,mixup,sct,&indx,&indy,&indz,&kstrt,&nvpy,&nvpz,
                &kxypi,&kxypp,&nyv,&kzpp,&kxypd,&kzpd,&nxhyzd,&nxyzhd);
   return;
}

/*--------------------------------------------------------------------*/
void cppfft32rmxz(float complex h[], int isign, int mixup[],
                  float complex sct[], int indx, int indy, int indz,
                  int kstrt, int nvpy, int nvpz, int kxypi, int kxypp,
                  int nzv, int kyzp, int kxypd, int kyzpd, int nxhyzd,
                  int nxyzhd) {
   ppfft32rmxz_(h,&isign,mixup,sct,&indx,&indy,&indz,&kstrt,&nvpy,&nvpz,
                &kxypi,&kxypp,&nzv,&kyzp,&kxypd,&kyzpd,&nxhyzd,&nxyzhd);
   return;
}

/*--------------------------------------------------------------------*/
void cppfft32rm3xx(float complex f[], int isign, int mixup[],
                   float complex sct[], int indx, int indy, int indz,
                   int kstrt, int nvp, int kypi, int kypp, int nxvh,
                   int kzpp, int kypd, int kzpd, int nxhyzd, int nxyzhd) {
   ppfft32rm3xx_(f,&isign,mixup,sct,&indx,&indy,&indz,&kstrt,&nvp,&kypi,
                 &kypp,&nxvh,&kzpp,&kypd,&kzpd,&nxhyzd,&nxyzhd);
   return;
}

/*--------------------------------------------------------------------*/
void cppfft32rm3xy(float complex g[], int isign, int mixup[],
                   float complex sct[], int indx, int indy, int indz,
                   int kstrt, int nvpy, int nvpz, int kxypi, int kxypp,
                   int nyv, int kzpp, int kxypd, int kzpd, int nxhyzd,
                   int nxyzhd) {
   ppfft32rm3xy_(g,&isign,mixup,sct,&indx,&indy,&indz,&kstrt,&nvpy,
                 &nvpz,&kxypi,&kxypp,&nyv,&kzpp,&kxypd,&kzpd,&nxhyzd,
                 &nxyzhd);
   return;
}

/*--------------------------------------------------------------------*/
void cppfft32rm3xz(float complex h[], int isign, int mixup[],
                   float complex sct[], int indx, int indy, int indz,
                   int kstrt, int nvpy, int nvpz, int kxypi, int kxypp,
                   int nzv, int kyzp, int kxypd, int kyzpd, int nxhyzd,
                   int nxyzhd) {
   ppfft32rm3xz_(h,&isign,mixup,sct,&indx,&indy,&indz,&kstrt,&nvpy,
                 &nvpz,&kxypi,&kxypp,&nzv,&kyzp,&kxypd,&kyzpd,&nxhyzd,
                 &nxyzhd);
   return;
}

/*--------------------------------------------------------------------*/
void cwppfft32rm(float complex f[], float complex g[], 
                 float complex h[], float complex bs[],
                 float complex br[], int isign, int ntpose, int mixup[],
                 float complex sct[], float *ttp, int indx, int indy,
                 int indz, int kstrt, int nvpy, int nvpz, int nxvh,
                 int nyv, int nzv, int kxyp, int kyp, int kyzp, int kzp,
                 int kxypd, int kypd, int kyzpd, int kzpd, int kzyp,
                 int nxhyzd, int nxyzhd) {
   wppfft32rm_(f,g,h,bs,br,&isign,&ntpose,mixup,sct,ttp,&indx,&indy,
               &indz,&kstrt,&nvpy,&nvpz,&nxvh,&nyv,&nzv,&kxyp,&kyp,
               &kyzp,&kzp,&kxypd,&kypd,&kyzpd,&kzpd,&kzyp,&nxhyzd,
               &nxyzhd);
   return;
}

/*--------------------------------------------------------------------*/
void cwppfft32rm3(float complex f[], float complex g[],
                  float complex h[], float complex bs[],
                  float complex br[], int isign, int ntpose,
                  int mixup[], float complex sct[], float *ttp,
                  int indx, int indy, int indz, int kstrt, int nvpy,
                  int nvpz, int nxvh, int nyv, int nzv, int kxyp,
                  int kyp, int kyzp, int kzp, int kxypd, int kypd,
                  int kyzpd, int kzpd, int kzyp, int nxhyzd, int nxyzhd) {
   wppfft32rm3_(f,g,h,bs,br,&isign,&ntpose,mixup,sct,ttp,&indx,&indy,
                &indz,&kstrt,&nvpy,&nvpz,&nxvh,&nyv,&nzv,&kxyp,&kyp,
                &kyzp,&kzp,&kxypd,&kypd,&kyzpd,&kzpd,&kzyp,&nxhyzd,
                &nxyzhd);
   return;
}

