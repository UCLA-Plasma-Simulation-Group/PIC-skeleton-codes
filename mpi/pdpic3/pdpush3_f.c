/* C Library for Skeleton 3D Darwin MPI PIC Code */
/* Wrappers for calling the Fortran routines from a C main program */

#include <complex.h>
#include "pdpush3.h"

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

void ppgbpush32l_(float *part, float *fxyz, float *bxyz, float *edges,
                  int *npp, int *noff, int *ihole, float *qbm,
                  float *dt, float *dtc, float *ek, int *nx, int *ny,
                  int *nz, int *idimp, int *npmax, int *nxv, int *nypmx,
                  int *nzpmx, int *idps, int *idds, int *ntmax,
                  int *ipbc);

void ppgpost32l_(float *part, float *q, int *npp, int *noff, float *qm,
                 int *idimp, int *npmax, int *nxv, int *nypmx,
                 int *nzpmx, int *idds);

void ppgjpost32l_(float *part, float *cu, float *edges, int *npp,
                  int *noff, int *ihole, float *qm, float *dt, int *nx,
                  int *ny, int *nz, int *idimp, int *npmax, int *nxv,
                  int *nypmx, int *nzpmx, int *idps, int *idds,
                  int *ntmax, int *ipbc);

void ppgmjpost32l_(float *part, float *amu, int *npp, int *noff,
                   float *qm, int *idimp, int *npmax, int *nxv,
                   int *nypmx, int *nzpmx, int *idds);

void ppgdjpost32l_(float *part, float *fxyz, float *bxyz, int *npp,
                   int *noff, float *dcu, float *amu, float *qm, 
                   float *qbm, float *dt, int *idimp, int *npmax,
                   int *nxv, int *nypmx, int *nzpmx, int *idds);

void ppgdcjpost32l_(float *part, float *fxyz, float *bxyz, int *npp,
                    int *noff, float *cu, float *dcu, float *amu,
                    float *qm, float *qbm, float *dt, int *idimp,
                    int *npmax, int *nxv, int *nypmx, int *nzpmx,
                    int *idds);

void ppdsortp32yzl_(float *parta, float *partb, int *npic, int *npp,
                    int *noff, int *nyzp, int *idimp, int *npmax,
                    int *nyzpm1, int *idds);

void ppcguard32xl_(float *fxyz, int *nyzp, int *nx, int *ndim, int *nxe,
                   int *nypmx, int *nzpmx, int *idds);

void ppaguard32xl_(float *q, int *nyzp, int *nx, int *nxe, int *nypmx,
                   int *nzpmx, int *idds);

void ppacguard32xl_(float *cu, int *nyzp, int *nx, int *ndim, int *nxe,
                    int *nypmx, int *nzpmx, int *idds);

void ppascfguard32l_(float *dcu, float *cus, int *nyzp, float *q2m0,
                     int *nx, int *nxe, int *nypmx, int *nzpmx,
                     int *idds);

void ppfwpminmx32_(float *qe, int *nyzp, float *qbme, float *wpmax,
                   float *wpmin, int *nx, int *nxe, int *nypmx,
                   int *nzpmx, int *idds);

void ppois332_(float complex *q, float complex *fxyz, int *isign,
               float complex *ffc, float *ax, float *ay, float *az,
               float *affp, float *we, int *nx, int *ny, int *nz,
               int *kstrt, int *nvpy, int *nvpz, int *nzv, int *kxyp,
               int *kyzp, int *nzhd);

void ppcuperp32_(float complex *cu, int *nx, int *ny, int *nz,
                 int *kstrt, int *nvpy, int *nvpz, int *nzv, int *kxyp,
                 int *kyzp);

void ppbbpoisp332_(float complex *cu, float complex *bxyz,
                   float complex *ffc, float *ci, float *wm, int *nx,
                   int *ny, int *nz, int *kstrt, int *nvpy, int *nvpz,
                   int *nzv, int *kxyp, int *kyzp, int *nzhd);

void ppbaddext32_(float *bxyz, int *nyzp, float *omx, float *omy,
                  float *omz, int *nx, int *nxe, int *nypmx, int *nzpmx,
                  int *idds);

void ppdcuperp32_(float complex *dcu, float complex *amu, int *nx,
                  int *ny, int *nz, int *kstrt, int *nvpy, int *nvpz,
                  int *nzv, int *kxyp, int *kyzp);

void ppadcuperp32_(float complex *dcu, float complex *amu, int *nx,
                   int *ny, int *nz, int *kstrt, int *nvpy, int *nvpz,
                   int *nzv, int *kxyp, int *kyzp);

void ppepoisp332_(float complex *dcu, float complex *exyz, int *isign,
                  float complex *ffe, float *ax, float *ay, float *az,
                  float *affp, float *wp0, float *ci, float *wf,
                  int *nx, int *ny, int *nz, int *kstrt, int *nvpy,
                  int *nvpz, int *nzv, int *kxyp, int *kyzp, int *nzhd);

void ppaddvrfield32_(float *a, float *b, float *c, int *ndim, int *nxe,
                     int *nypmx, int *nzpmx);

void wpfft32rinit_(int *mixup, float complex *sct, int *indx, int *indy,
                   int *indz, int *nxhyzd, int *nxyzhd);

void ppfft32rxx_(float complex *f, int *isign, int *mixup,
                 float complex *sct, int *indx, int *indy, int *indz,
                 int *kstrt, int *nvp, int *kypi, int *kypp, int *nxvh,
                 int *kzpp, int *kypd, int *kzpd, int *nxhyzd,
                 int *nxyzhd);

void ppfft32rxy_(float complex *g, int *isign, int *mixup,
                 float complex *sct, int *indx, int *indy, int *indz,
                 int *kstrt, int *nvpy, int *nvpz, int *kxypi,
                 int *kxypp, int *nyv, int *kzpp, int *kxypd, int *kzpd,
                 int *nxhyzd, int *nxyzhd);

void ppfft32rxz_(float complex *h, int *isign, int *mixup,
                 float complex *sct, int *indx, int *indy, int *indz,
                 int *kstrt, int *nvpy, int *nvpz, int *kxypi,
                 int *kxypp, int *nzv, int *kyzp, int *kxypd,
                 int *kyzpd, int *nxhyzd, int *nxyzhd);

void ppfft32r3xx_(float complex *f, int *isign, int *mixup,
                  float complex *sct, int *indx, int *indy, int *indz,
                  int *kstrt, int *nvp, int *kypi, int *kypp, int *nxvh,
                  int *kzpp, int *kypd, int *kzpd, int *nxhyzd,
                  int *nxyzhd);

void ppfft32r3xy_(float complex *g, int *isign, int *mixup,
                  float complex *sct, int *indx, int *indy, int *indz,
                  int *kstrt, int *nvpy, int *nvpz, int *kxypi,
                  int *kxypp, int *nyv, int *kzpp, int *kxypd, 
                  int *kzpd, int *nxhyzd, int *nxyzhd);

void ppfft32r3xz_(float complex *h, int *isign, int *mixup,
                  float complex *sct, int *indx, int *indy, int *indz,
                  int *kstrt, int *nvpy, int *nvpz, int *kxypi,
                  int *kxypp, int *nzv, int *kyzp, int *kxypd, 
                  int *kyzpd, int *nxhyzd, int *nxyzhd);

void ppfft32rnxx_(float complex *f, float complex *ss, int *isign,
                  int *mixup, float complex *sct, int *indx, int *indy,
                  int *indz, int *kstrt, int *nvp, int *kypi, int *kypp,
                  int *nxvh, int *kzpp, int *kypd, int *kzpd, int *ndim,
                  int *nxhyzd, int *nxyzhd);

void ppfft32rnxy_(float complex *g, int *isign, int *mixup,
                  float complex *sct, int *indx, int *indy, int *indz,
                  int *kstrt, int *nvpy, int *nvpz, int *kxypi,
                  int *kxypp, int *nyv, int *kzpp, int *kxypd,
                  int *kzpd, int *ndim, int *nxhyzd, int *nxyzhd);

void ppfft32rnxz_(float complex *h, int *isign, int *mixup,
                  float complex *sct, int *indx, int *indy, int *indz,
                  int *kstrt, int *nvpy, int *nvpz, int *kxypi,
                  int *kxypp, int *nzv, int *kyzp, int *kxypd,
                  int *kyzpd, int *ndim, int *nxhyzd, int *nxyzhd);

void wppfft32r_(float complex *f, float complex *g, float complex *h,
                float complex *bs, float complex *br, int *isign,
                int *ntpose, int *mixup, float complex *sct,
                float *ttp, int *indx, int *indy, int *indz, int *kstrt,
                int *nvpy, int *nvpz, int *nxvh, int *nyv, int *nzv,
                int *kxyp, int *kyp, int *kyzp, int *kzp, int *kxypd,
                int *kypd, int *kyzpd, int *kzpd, int *kzyp,
                int *nxhyzd, int *nxyzhd);

void wppfft32r3_(float complex *f, float complex *g, float complex *h,
                 float complex *bs, float complex *br, int *isign,
                 int *ntpose, int *mixup, float complex *sct,
                 float *ttp, int *indx, int *indy, int *indz,
                 int *kstrt, int *nvpy, int *nvpz, int *nxvh, int *nyv,
                 int *nzv, int *kxyp, int *kyp, int *kyzp, int *kzp,
                 int *kxypd, int *kypd, int *kyzpd, int *kzpd,
                 int *kzyp, int *nxhyzd, int *nxyzhd);

void wppfft32rn_(float complex *f, float complex *g, float complex *h,
                 float complex *bs, float complex *br,
                 float complex *ss, int *isign, int *ntpose,
                 int *mixup, float complex *sct, float *ttp, int *indx,
                 int *indy, int *indz, int *kstrt, int *nvpy,
                 int *nvpz, int *nxvh, int *nyv, int *nzv, int *kxyp,
                 int *kyp, int *kyzp, int *kzp, int *kxypd, int *kypd,
                 int *kyzpd, int *kzpd, int *kzyp, int *ndim,
                 int *nxhyzd, int *nxyzhd);

void ppswapc32n_(float *f, float *s, int *isign, int *nxh, int *kypi,
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
void cppgbpush32l(float part[], float fxyz[], float bxyz[],
                  float edges[], int npp, int noff[], int ihole[],
                  float qbm, float dt, float dtc, float *ek, int nx,
                  int ny, int nz, int idimp, int npmax, int nxv,
                  int nypmx, int nzpmx, int idps, int idds, int ntmax,
                  int ipbc) {
   ppgbpush32l_(part,fxyz,bxyz,edges,&npp,noff,ihole,&qbm,&dt,&dtc,ek,
                &nx,&ny,&nz,&idimp,&npmax,&nxv,&nypmx,&nzpmx,&idps,
                &idds,&ntmax,&ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cppgpost32l(float part[], float q[], int npp, int noff[], float qm,
                 int idimp, int npmax, int nxv, int nypmx, int nzpmx,
                 int idds) {
   ppgpost32l_(part,q,&npp,noff,&qm,&idimp,&npmax,&nxv,&nypmx,&nzpmx,
               &idds);
   return;
}

/*--------------------------------------------------------------------*/
void cppgjpost32l(float part[], float cu[], float edges[], int npp,
                  int noff[], int ihole[], float qm, float dt, int nx,
                  int ny, int nz, int idimp, int npmax, int nxv,
                  int nypmx, int nzpmx, int idps, int idds, int ntmax,
                  int ipbc) {
   ppgjpost32l_(part,cu,edges,&npp,noff,ihole,&qm,&dt,&nx,&ny,&nz,
                &idimp,&npmax,&nxv,&nypmx,&nzpmx,&idps,&idds,&ntmax,
                &ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cppgmjpost32l(float part[], float amu[], int npp, int noff[],
                   float qm, int idimp, int npmax, int nxv, int nypmx,
                   int nzpmx, int idds) {
   ppgmjpost32l_(part,amu,&npp,noff,&qm,&idimp,&npmax,&nxv,&nypmx,
                 &nzpmx,&idds);
   return;
}

/*--------------------------------------------------------------------*/
void cppgdjpost32l(float part[], float fxyz[], float bxyz[], int npp,
                   int noff[], float dcu[], float amu[], float qm, 
                   float qbm, float dt, int idimp, int npmax, int nxv,
                   int nypmx, int nzpmx, int idds) {
   ppgdjpost32l_(part,fxyz,bxyz,&npp,noff,dcu,amu,&qm,&qbm,&dt,&idimp,
                 &npmax,&nxv,&nypmx,&nzpmx,&idds);
   return;
}

/*--------------------------------------------------------------------*/
void cppgdcjpost32l(float part[], float fxyz[], float bxyz[], int npp,
                    int noff[], float cu[], float dcu[], float amu[],
                    float qm, float qbm, float dt, int idimp, int npmax,
                    int nxv, int nypmx, int nzpmx, int idds) {
   ppgdcjpost32l_(part,fxyz,bxyz,&npp,noff,cu,dcu,amu,&qm,&qbm,&dt,
                  &idimp,&npmax,&nxv,&nypmx,&nzpmx,&idds);
   return;
}

/*--------------------------------------------------------------------*/
void cppdsortp32yzl(float parta[], float partb[], int npic[], int npp,
                    int noff[], int nyzp[], int idimp, int npmax,
                    int nyzpm1, int idds) {
   ppdsortp32yzl_(parta,partb,npic,&npp,noff,nyzp,&idimp,&npmax,&nyzpm1,
                  &idds);
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
void cppascfguard32l(float dcu[], float cus[], int nyzp[], float q2m0,
                     int nx, int nxe, int nypmx, int nzpmx, int idds) {
   ppascfguard32l_(dcu,cus,nyzp,&q2m0,&nx,&nxe,&nypmx,&nzpmx,&idds);
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
void cppois332(float complex q[], float complex fxyz[], int isign,
               float complex ffc[], float ax, float ay, float az,
               float affp, float *we, int nx, int ny, int nz, int kstrt,
               int nvpy, int nvpz, int nzv, int kxyp, int kyzp,
               int nzhd) {
   ppois332_(q,fxyz,&isign,ffc,&ax,&ay,&az,&affp,we,&nx,&ny,&nz,&kstrt,
             &nvpy,&nvpz,&nzv,&kxyp,&kyzp,&nzhd);
   return;
}

/*--------------------------------------------------------------------*/
void cppcuperp32(float complex cu[], int nx, int ny, int nz, int kstrt,
                 int nvpy, int nvpz, int nzv, int kxyp, int kyzp) {
   ppcuperp32_(cu,&nx,&ny,&nz,&kstrt,&nvpy,&nvpz,&nzv,&kxyp,&kyzp);
   return;
}

/*--------------------------------------------------------------------*/
void cppbbpoisp332(float complex cu[], float complex bxyz[],
                   float complex ffc[], float ci, float *wm, int nx,
                   int ny, int nz, int kstrt, int nvpy, int nvpz,
                   int nzv, int kxyp, int kyzp, int nzhd) {
   ppbbpoisp332_(cu,bxyz,ffc,&ci,wm,&nx,&ny,&nz,&kstrt,&nvpy,&nvpz,&nzv,
                 &kxyp,&kyzp,&nzhd);
   return;
}

/*--------------------------------------------------------------------*/
void cppbaddext32(float bxyz[], int nyzp[], float omx, float omy,
                  float omz, int nx, int nxe, int nypmx, int nzpmx,
                  int idds) {
   ppbaddext32_(bxyz,nyzp,&omx,&omy,&omz,&nx,&nxe,&nypmx,&nzpmx,&idds);
   return;
}

/*--------------------------------------------------------------------*/
void cppdcuperp32(float complex dcu[], float complex amu[], int nx,
                  int ny, int nz, int kstrt, int nvpy, int nvpz,
                  int nzv, int kxyp, int kyzp) {
   ppdcuperp32_(dcu,amu,&nx,&ny,&nz,&kstrt,&nvpy,&nvpz,&nzv,&kxyp,
                &kyzp);
   return;
}

/*--------------------------------------------------------------------*/
void cppadcuperp32(float complex dcu[], float complex amu[], int nx,
                   int ny, int nz, int kstrt, int nvpy, int nvpz,
                   int nzv, int kxyp, int kyzp) {
   ppadcuperp32_(dcu,amu,&nx,&ny,&nz,&kstrt,&nvpy,&nvpz,&nzv,&kxyp,
                 &kyzp);
   return;
}

/*--------------------------------------------------------------------*/
void cppepoisp332(float complex dcu[], float complex exyz[], int isign,
                  float complex ffe[], float ax, float ay, float az,
                  float affp, float wp0, float ci, float *wf, int nx,
                  int ny, int nz, int kstrt, int nvpy, int nvpz,
                  int nzv, int kxyp, int kyzp, int nzhd) {
   ppepoisp332_(dcu,exyz,&isign,ffe,&ax,&ay,&az,&affp,&wp0,&ci,wf,&nx,
                &ny,&nz,&kstrt,&nvpy,&nvpz,&nzv,&kxyp,&kyzp,&nzhd);
   return;
}

/*--------------------------------------------------------------------*/
void cppaddvrfield32(float a[], float b[], float c[], int ndim, int nxe,
                     int nypmx, int nzpmx) {
   ppaddvrfield32_(a,b,c,&ndim,&nxe,&nypmx,&nzpmx);
   return;
}

/*--------------------------------------------------------------------*/
void cwpfft32rinit(int mixup[], float complex sct[], int indx, int indy,
                   int indz, int nxhyzd, int nxyzhd) {
   wpfft32rinit_(mixup,sct,&indx,&indy,&indz,&nxhyzd,&nxyzhd);
   return;
}

/*--------------------------------------------------------------------*/
void cppfft32rxx(float complex f[], int isign, int mixup[],
                 float complex sct[], int indx, int indy, int indz,
                 int kstrt, int nvp, int kypi, int kypp, int nxvh,
                 int kzpp, int kypd, int kzpd, int nxhyzd, int nxyzhd) {
   ppfft32rxx_(f,&isign,mixup,sct,&indx,&indy,&indz,&kstrt,&nvp,&kypi,
               &kypp,&nxvh,&kzpp,&kypd,&kzpd,&nxhyzd,&nxyzhd);
   return;
}

/*--------------------------------------------------------------------*/
void cppfft32rxy(float complex g[], int isign, int mixup[],
                 float complex sct[], int indx, int indy, int indz,
                 int kstrt, int nvpy, int nvpz, int kxypi, int kxypp,
                 int nyv, int kzpp, int kxypd, int kzpd, int nxhyzd,
                 int nxyzhd) {
   ppfft32rxy_(g,&isign,mixup,sct,&indx,&indy,&indz,&kstrt,&nvpy,&nvpz,
               &kxypi,&kxypp,&nyv,&kzpp,&kxypd,&kzpd,&nxhyzd,&nxyzhd);
   return;
}

/*--------------------------------------------------------------------*/
void cppfft32rxz(float complex h[], int isign, int mixup[],
                 float complex sct[], int indx, int indy, int indz,
                 int kstrt, int nvpy, int nvpz, int kxypi, int kxypp,
                 int nzv, int kyzp, int kxypd, int kyzpd, int nxhyzd,
                 int nxyzhd) {
   ppfft32rxz_(h,&isign,mixup,sct,&indx,&indy,&indz,&kstrt,&nvpy,&nvpz,
               &kxypi,&kxypp,&nzv,&kyzp,&kxypd,&kyzpd,&nxhyzd,&nxyzhd);
   return;
}

/*--------------------------------------------------------------------*/
void cppfft32r3xx(float complex f[], int isign, int mixup[],
                  float complex sct[], int indx, int indy, int indz,
                  int kstrt, int nvp, int kypi, int kypp, int nxvh,
                  int kzpp, int kypd, int kzpd, int nxhyzd, int nxyzhd) {
   ppfft32r3xx_(f,&isign,mixup,sct,&indx,&indy,&indz,&kstrt,&nvp,&kypi,
                &kypp,&nxvh,&kzpp,&kypd,&kzpd,&nxhyzd,&nxyzhd);
   return;
}

/*--------------------------------------------------------------------*/
void cppfft32r3xy(float complex g[], int isign, int mixup[],
                  float complex sct[], int indx, int indy, int indz,
                  int kstrt, int nvpy, int nvpz, int kxypi, int kxypp,
                  int nyv, int kzpp, int kxypd, int kzpd, int nxhyzd,
                  int nxyzhd) {
   ppfft32r3xy_(g,&isign,mixup,sct,&indx,&indy,&indz,&kstrt,&nvpy,&nvpz,
                &kxypi,&kxypp,&nyv,&kzpp,&kxypd,&kzpd,&nxhyzd,&nxyzhd);
   return;
}

/*--------------------------------------------------------------------*/
void cppfft32r3xz(float complex h[], int isign, int mixup[],
                  float complex sct[], int indx, int indy, int indz,
                  int kstrt, int nvpy, int nvpz, int kxypi, int kxypp,
                  int nzv, int kyzp, int kxypd, int kyzpd, int nxhyzd,
                  int nxyzhd) {
   ppfft32r3xz_(h,&isign,mixup,sct,&indx,&indy,&indz,&kstrt,&nvpy,&nvpz,
                &kxypi,&kxypp,&nzv,&kyzp,&kxypd,&kyzpd,&nxhyzd,&nxyzhd);
   return;
}
/*--------------------------------------------------------------------*/
void cppfft32rnxx(float complex f[], float complex ss[], int isign,
                  int mixup[], float complex sct[], int indx, int indy,
                  int indz, int kstrt, int nvp, int kypi, int kypp,
                  int nxvh, int kzpp, int kypd, int kzpd, int ndim,
                  int nxhyzd, int nxyzhd) {
   ppfft32rnxx_(f,ss,&isign,mixup,sct,&indx,&indy,&indz,&kstrt,&nvp,
                &kypi,&kypp,&nxvh,&kzpp,&kypd,&kzpd,&ndim,&nxhyzd,
                &nxyzhd);
   return;
}

/*--------------------------------------------------------------------*/
void cppfft32rnxy(float complex g[], int isign, int mixup[],
                  float complex sct[], int indx, int indy, int indz,
                  int kstrt, int nvpy, int nvpz, int kxypi, int kxypp,
                  int nyv, int kzpp, int kxypd, int kzpd, int ndim,
                  int nxhyzd, int nxyzhd) {
   ppfft32rnxy_(g,&isign,mixup,sct,&indx,&indy,&indz,&kstrt,&nvpy,&nvpz,
                &kxypi,&kxypp,&nyv,&kzpp,&kxypd,&kzpd,&ndim,&nxhyzd,
                &nxyzhd);
   return;
}

/*--------------------------------------------------------------------*/
void cppfft32rnxz(float complex h[], int isign, int mixup[],
                  float complex sct[], int indx, int indy, int indz,
                  int kstrt, int nvpy, int nvpz, int kxypi, int kxypp,
                  int nzv, int kyzp, int kxypd, int kyzpd, int ndim,
                  int nxhyzd, int nxyzhd) {
   ppfft32rnxz_(h,&isign,mixup,sct,&indx,&indy,&indz,&kstrt,&nvpy,&nvpz,
                &kxypi,&kxypp,&nzv,&kyzp,&kxypd,&kyzpd,&ndim,&nxhyzd,
                &nxyzhd);
   return;
}

/*--------------------------------------------------------------------*/
void cwppfft32r(float complex f[], float complex g[], float complex h[],
                float complex bs[], float complex br[], int isign,
                int ntpose, int mixup[], float complex sct[],
                float *ttp, int indx, int indy, int indz, int kstrt,
                int nvpy, int nvpz, int nxvh, int nyv, int nzv,
                int kxyp, int kyp, int kyzp, int kzp, int kxypd,
                int kypd, int kyzpd, int kzpd, int kzyp, int nxhyzd,
                int nxyzhd) {
   wppfft32r_(f,g,h,bs,br,&isign,&ntpose,mixup,sct,ttp,&indx,&indy,
              &indz,&kstrt,&nvpy,&nvpz,&nxvh,&nyv,&nzv,&kxyp,&kyp,&kyzp,
              &kzp,&kxypd,&kypd,&kyzpd,&kzpd,&kzyp,&nxhyzd,&nxyzhd);
   return;
}

/*--------------------------------------------------------------------*/
void cwppfft32r3(float complex f[], float complex g[],
                 float complex h[], float complex bs[],
                 float complex br[], int isign, int ntpose, int mixup[],
                 float complex sct[], float *ttp, int indx, int indy,
                 int indz, int kstrt, int nvpy, int nvpz, int nxvh,
                 int nyv, int nzv, int kxyp, int kyp, int kyzp, int kzp,
                 int kxypd, int kypd, int kyzpd, int kzpd, int kzyp,
                 int nxhyzd, int nxyzhd) {
   wppfft32r3_(f,g,h,bs,br,&isign,&ntpose,mixup,sct,ttp,&indx,&indy,
               &indz,&kstrt,&nvpy,&nvpz,&nxvh,&nyv,&nzv,&kxyp,&kyp,
               &kyzp,&kzp,&kxypd,&kypd,&kyzpd,&kzpd,&kzyp,&nxhyzd,
               &nxyzhd);
   return;
}

/*--------------------------------------------------------------------*/
void cwppfft32rn(float complex f[], float complex g[],
                 float complex h[], float complex bs[],
                 float complex br[], float complex ss[], int isign,
                 int ntpose, int mixup[], float complex sct[], 
                 float *ttp, int indx, int indy, int indz, int kstrt,
                 int nvpy, int nvpz, int nxvh, int nyv, int nzv,
                 int kxyp, int kyp, int kyzp, int kzp, int kxypd,
                 int kypd, int kyzpd, int kzpd, int kzyp, int ndim,
                 int nxhyzd, int nxyzhd) {
   wppfft32rn_(f,g,h,bs,br,ss,&isign,&ntpose,mixup,sct,ttp,&indx,&indy,
               &indz,&kstrt,&nvpy,&nvpz,&nxvh,&nyv,&nzv,&kxyp,&kyp,
               &kyzp,&kzp,&kxypd,&kypd,&kyzpd,&kzpd,&kzyp,&ndim,
               &nxhyzd,&nxyzhd);
   return;
}

/*--------------------------------------------------------------------*/
void cppswapc32n(float f[], float s[], int isign, int nxh, int kypi,
                 int kypt, int nxvh, int kzpp, int kypd, int kzpd,
                 int ndim) {
   ppswapc32n_(f,s,&isign,&nxh,&kypi,&kypt,&nxvh,&kzpp,&kypd,&kzpd,
               &ndim);
   return;
}

