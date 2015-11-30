/* C Library for Skeleton 3D Darwin MPI/OpenMP PIC Code */
/* Wrappers for calling the Fortran routines from a C main program */

#include <complex.h>
#include "mpdpush3.h"

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

void ppgppost32l_(float *ppart, float *q, int *kpic, int *noff,
                  float *qm, int *nppmx, int *idimp, int *mx, int *my,
                  int *mz, int *nxv, int *nypmx, int *nzpmx, int *mx1,
                  int *myp1, int *mxyzp1, int *idds);

void ppgjppost32l_(float *ppart, float *cu, int *kpic, int *noff,
                   float *qm, float *dt, int *nppmx, int *idimp,
                   int *nx, int *ny, int *nz, int *mx, int *my, int *mz,
                   int *nxv, int *nypmx, int *nzpmx, int *mx1, 
                   int *myp1, int *mxyzp1, int *idds, int *ipbc);

void ppgmjppost32l_(float *ppart, float *amu, int *kpic, int *noff,
                    float *qm, int *nppmx, int *idimp, int *mx, int *my,
                    int *mz, int *nxv, int *nypmx, int *nzpmx, int *mx1,
                    int *myp1, int *mxyzp1, int *idds);

void ppgdjppost32l_(float *ppart, float *fxyz, float *bxyz, int *kpic,
                    int *noff, int *nyzp, float *dcu, float *amu,
                    float *qm, float *qbm, float *dt, int *idimp,
                    int *nppmx, int *nx, int *mx, int *my, int *mz,
                    int *nxv, int *nypmx, int *nzpmx, int *mx1,
                    int *myp1, int *mxyzp1, int *idds);

void ppgdcjppost32l_(float *ppart, float *fxyz, float *bxyz, int *kpic, 
                     int *noff, int *nyzp, float *cu, float *dcu,
                     float *amu, float *qm, float *qbm, float *dt,
                     int *idimp, int *nppmx, int *nx, int *mx, int *my,
                     int *mz, int *nxv, int *nypmx, int *nzpmx,
                     int *mx1, int *myp1, int *mxyzp1, int *idds);

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

void mppascfguard32l_(float *dcu, float *cus, int *nyzp, float *q2m0,
                      int *nx, int *nxe, int *nypmx, int *nzpmx,
                      int *idds);

void ppfwpminmx32_(float *qe, int *nyzp, float *qbme, float *wpmax,
                   float *wpmin, int *nx, int *nxe, int *nypmx,
                   int *nzpmx, int *idds);

void mppois332_(float complex *q, float complex *fxyz, int *isign,
                float complex *ffc, float *ax, float *ay, float *az,
                float *affp, float *we, int *nx, int *ny, int *nz,
                int *kstrt, int *nvpy, int *nvpz, int *nzv, int *kxyp,
                int *kyzp, int *nzhd);

void mppcuperp32_(float complex *cu, int *nx, int *ny, int *nz,
                  int *kstrt, int *nvpy, int *nvpz, int *nzv, int *kxyp,
                  int *kyzp);

void mppbbpoisp332_(float complex *cu, float complex *bxyz,
                    float complex *ffc, float *ci, float *wm, int *nx,
                    int *ny, int *nz, int *kstrt, int *nvpy, int *nvpz,
                    int *nzv, int *kxyp, int *kyzp, int *nzhd);

void mppbaddext32_(float *bxyz, int *nyzp, float *omx, float *omy,
                   float *omz, int *nx, int *nxe, int *nypmx,
                   int *nzpmx, int *idds);

void mppdcuperp32_(float complex *dcu, float complex *amu, int *nx,
                   int *ny, int *nz, int *kstrt, int *nvpy, int *nvpz,
                   int *nzv, int *kxyp, int *kyzp);

void mppadcuperp32_(float complex *dcu, float complex *amu, int *nx,
                    int *ny, int *nz, int *kstrt, int *nvpy, int *nvpz,
                    int *nzv, int *kxyp, int *kyzp);

void mppepoisp332_(float complex *dcu, float complex *exyz, int *isign,
                   float complex *ffe, float *ax, float *ay, float *az,
                   float *affp, float *wp0, float *ci, float *wf,
                   int *nx, int *ny, int *nz, int *kstrt, int *nvpy,
                   int *nvpz, int *nzv, int *kxyp, int *kyzp, 
                   int *nzhd);

void mppaddvrfield32_(float *a, float *b, float *c, int *ndim, int *nxe,
                      int *nypmx, int *nzpmx);

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

void ppfft32rmnxx_(float complex *f, float complex *ss, int *isign,
                   int *mixup, float complex *sct, int *indx, int *indy,
                   int *indz, int *kstrt, int *nvp, int *kypi,
                   int *kypp, int *nxvh, int *kzpp, int *kypd,
                   int *kzpd, int *ndim, int *nxhyzd, int *nxyzhd);

void ppfft32rmnxy_(float complex *g, int *isign, int *mixup,
                   float complex *sct, int *indx, int *indy, int *indz,
                   int *kstrt, int *nvpy, int *nvpz, int *kxypi,
                   int *kxypp, int *nyv, int *kzpp, int *kxypd,
                   int *kzpd, int *ndim, int *nxhyzd, int *nxyzhd);

void ppfft32rmnxz_(float complex *h, int *isign, int *mixup,
                   float complex *sct, int *indx, int *indy, int *indz,
                   int *kstrt, int *nvpy, int *nvpz, int *kxypi,
                   int *kxypp, int *nzv, int *kyzp, int *kxypd,
                   int *kyzpd, int *ndim, int *nxhyzd, int *nxyzhd);

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

void wppfft32rmn_(float complex *f, float complex *g, float complex *h,
                  float complex *bs, float complex *br,
                  float complex *ss, int *isign, int *ntpose,
                  int *mixup, float complex *sct, float *ttp, int *indx,
                  int *indy, int *indz, int *kstrt, int *nvpy,
                  int *nvpz, int *nxvh, int *nyv, int *nzv, int *kxyp,
                  int *kyp, int *kyzp, int *kzp, int *kxypd, int *kypd,
                  int *kyzpd, int *kzpd, int *kzyp, int *ndim,
                  int *nxhyzd, int *nxyzhd);

void mppswapc32n_(float *f, float *s, int *isign, int *nxh, int *kypi,
                  int *kypt, int *nxvh, int *kzpp, int *kypd, int *kzpd,
                  int *ndim);

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
void cppgmjppost32l(float ppart[], float amu[], int kpic[], int noff[],
                    float qm, int nppmx, int idimp, int mx, int my,
                    int mz, int nxv, int nypmx, int nzpmx, int mx1,
                    int myp1, int mxyzp1, int idds) {
   ppgmjppost32l_(ppart,amu,kpic,noff,&qm,&nppmx,&idimp,&mx,&my,&mz,
                  &nxv,&nypmx,&nzpmx,&mx1,&myp1,&mxyzp1,&idds);
   return;
}

/*--------------------------------------------------------------------*/
void cppgdjppost32l(float ppart[], float fxyz[], float bxyz[],
                    int kpic[], int noff[], int nyzp[], float dcu[],
                    float amu[], float qm, float qbm, float dt,
                    int idimp, int nppmx, int nx, int mx, int my,
                    int mz, int nxv, int nypmx, int nzpmx, int mx1,
                    int myp1, int mxyzp1, int idds) {
   ppgdjppost32l_(ppart,fxyz,bxyz,kpic,noff,nyzp,dcu,amu,&qm,&qbm,&dt,
                  &idimp,&nppmx,&nx,&mx,&my,&mz,&nxv,&nypmx,&nzpmx,&mx1,
                  &myp1,&mxyzp1,&idds);
   return;
}

/*--------------------------------------------------------------------*/
void cppgdcjppost32l(float ppart[], float fxyz[], float bxyz[],
                     int kpic[], int noff[], int nyzp[], float cu[],
                     float dcu[], float amu[], float qm, float qbm,
                     float dt, int idimp, int nppmx, int nx, int mx,
                     int my, int mz, int nxv, int nypmx, int nzpmx,
                     int mx1, int myp1, int mxyzp1, int idds) {
   ppgdcjppost32l_(ppart,fxyz,bxyz,kpic,noff,nyzp,cu,dcu,amu,&qm,&qbm,
                   &dt,&idimp,&nppmx,&nx,&mx,&my,&mz,&nxv,&nypmx,&nzpmx,
                   &mx1,&myp1,&mxyzp1,&idds);
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
void cmppascfguard32l(float dcu[], float cus[], int nyzp[], float q2m0,
                      int nx, int nxe, int nypmx, int nzpmx, int idds) {
   mppascfguard32l_(dcu,cus,nyzp,&q2m0,&nx,&nxe,&nypmx,&nzpmx,&idds);
   return;
}

/*--------------------------------------------------------------------*/
void cppfwpminmx32(float qe[], int nyzp[], float qbme, float *wpmax,
                   float *wpmin, int nx, int nxe, int nypmx, int nzpmx,
                   int idds) {
   ppfwpminmx32_(qe,nyzp,&qbme,wpmax,wpmin,&nx,&nxe,&nypmx,&nzpmx,
                 &idds);
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
void cmppbbpoisp332(float complex cu[], float complex bxyz[],
                    float complex ffc[], float ci, float *wm, int nx,
                    int ny, int nz, int kstrt, int nvpy, int nvpz,
                    int nzv, int kxyp, int kyzp, int nzhd) {
   mppbbpoisp332_(cu,bxyz,ffc,&ci,wm,&nx,&ny,&nz,&kstrt,&nvpy,&nvpz,
                  &nzv,&kxyp,&kyzp,&nzhd);
   return;
}

/*--------------------------------------------------------------------*/
void cmppbaddext32(float bxyz[], int nyzp[], float omx, float omy,
                   float omz, int nx, int nxe, int nypmx, int nzpmx,
                   int idds) {
   mppbaddext32_(bxyz,nyzp,&omx,&omy,&omz,&nx,&nxe,&nypmx,&nzpmx,&idds);
   return;
}

/*--------------------------------------------------------------------*/
void cmppdcuperp32(float complex dcu[], float complex amu[], int nx,
                   int ny, int nz, int kstrt, int nvpy, int nvpz,
                   int nzv, int kxyp, int kyzp) {
   mppdcuperp32_(dcu,amu,&nx,&ny,&nz,&kstrt,&nvpy,&nvpz,&nzv,&kxyp,
                 &kyzp);
   return;
}

/*--------------------------------------------------------------------*/
void cmppadcuperp32(float complex dcu[], float complex amu[], int nx,
                    int ny, int nz, int kstrt, int nvpy, int nvpz,
                    int nzv, int kxyp, int kyzp) {
   mppadcuperp32_(dcu,amu,&nx,&ny,&nz,&kstrt,&nvpy,&nvpz,&nzv,&kxyp,
                  &kyzp);
   return;
}

/*--------------------------------------------------------------------*/
void cmppepoisp332(float complex dcu[], float complex exyz[], int isign,
                   float complex ffe[], float ax, float ay, float az,
                   float affp, float wp0, float ci, float *wf, int nx,
                   int ny, int nz, int kstrt, int nvpy, int nvpz,
                   int nzv, int kxyp, int kyzp, int nzhd) {
   mppepoisp332_(dcu,exyz,&isign,ffe,&ax,&ay,&az,&affp,&wp0,&ci,wf,&nx,
                 &ny,&nz,&kstrt,&nvpy,&nvpz,&nzv,&kxyp,&kyzp,&nzhd);
   return;
}

/*--------------------------------------------------------------------*/
void cmppaddvrfield32(float a[], float b[], float c[], int ndim,
                      int nxe, int nypmx, int nzpmx) {
   mppaddvrfield32_(a,b,c,&ndim,&nxe,&nypmx,&nzpmx);
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
void cppfft32rmnxx(float complex f[], float complex ss[], int isign,
                   int mixup[], float complex sct[], int indx, int indy,
                   int indz, int kstrt, int nvp, int kypi, int kypp,
                   int nxvh, int kzpp, int kypd, int kzpd, int ndim,
                   int nxhyzd, int nxyzhd) {
   ppfft32rmnxx_(f,ss,&isign,mixup,sct,&indx,&indy,&indz,&kstrt,&nvp,
                 &kypi,&kypp,&nxvh,&kzpp,&kypd,&kzpd,&ndim,&nxhyzd,
                 &nxyzhd);
   return;
}

/*--------------------------------------------------------------------*/
void cppfft32rmnxy(float complex g[], int isign, int mixup[],
                   float complex sct[], int indx, int indy, int indz,
                   int kstrt, int nvpy, int nvpz, int kxypi, int kxypp,
                   int nyv, int kzpp, int kxypd, int kzpd, int ndim,
                   int nxhyzd, int nxyzhd) {
   ppfft32rmnxy_(g,&isign,mixup,sct,&indx,&indy,&indz,&kstrt,&nvpy,
                 &nvpz,&kxypi,&kxypp,&nyv,&kzpp,&kxypd,&kzpd,&ndim,
                 &nxhyzd,&nxyzhd);
   return;
}

/*--------------------------------------------------------------------*/
void cppfft32rmnxz(float complex h[], int isign, int mixup[],
                   float complex sct[], int indx, int indy, int indz,
                   int kstrt, int nvpy, int nvpz, int kxypi, int kxypp,
                   int nzv, int kyzp, int kxypd, int kyzpd, int ndim,
                   int  nxhyzd, int nxyzhd) {
   ppfft32rmnxz_(h,&isign,mixup,sct,&indx,&indy,&indz,&kstrt,&nvpy,
                 &nvpz,&kxypi,&kxypp,&nzv,&kyzp,&kxypd,&kyzpd,&ndim,
                 &nxhyzd,&nxyzhd);
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

/*--------------------------------------------------------------------*/
void cwppfft32rmn(float complex f[], float complex g[],
                  float complex h[], float complex bs[],
                  float complex br[], float complex ss[], int isign,
                  int ntpose, int mixup[], float complex sct[], 
                  float *ttp, int indx, int indy, int indz, int kstrt,
                  int nvpy, int nvpz, int nxvh, int nyv, int nzv,
                  int kxyp, int kyp, int kyzp, int kzp, int kxypd,
                  int kypd, int kyzpd, int kzpd, int kzyp, int ndim,
                  int nxhyzd, int nxyzhd) {
   wppfft32rmn_(f,g,h,bs,br,ss,&isign,&ntpose,mixup,sct,ttp,&indx,&indy,
                &indz,&kstrt,&nvpy,&nvpz,&nxvh,&nyv,&nzv,&kxyp,&kyp,
                &kyzp,&kzp,&kxypd,&kypd,&kyzpd,&kzpd,&kzyp,&ndim,
                &nxhyzd,&nxyzhd);
   return;
}

/*--------------------------------------------------------------------*/
void cmppswapc32n(float f[], float s[], int isign, int nxh, int kypi,
                  int kypt, int nxvh, int kzpp, int kypd, int kzpd,
                  int ndim) {
   mppswapc32n_(f,s,&isign,&nxh,&kypi,&kypt,&nxvh,&kzpp,&kypd,&kzpd,
                &ndim);
   return;
}

