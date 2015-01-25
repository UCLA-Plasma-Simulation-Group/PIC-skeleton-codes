/* C Library for Skeleton 2-1/2D Electromagnetic MPI PIC Code */
/* Wrappers for calling the Fortran routines from a C main program */

#include <complex.h>
#include "pbpush2.h"

double ranorm_();

void pdicomp2l_(float *edges, int *nyp, int *noff, int *nypmx,
                int *nypmn, int *ny, int *kstrt, int *nvp, int *idps);

void pdistr2h_(float *part, float *edges, int *npp, int *nps,
               float *vtx, float *vty, float *vtz, float *vdx,
               float *vdy, float *vdz, int *npx, int *npy, int *nx,
               int *ny, int *idimp, int *npmax, int *idps, int *ipbc,
               int *ierr);

void ppgbpush23l_(float *part, float *fxy, float *bxy, float *edges,
                  int *npp, int *noff, int *ihole, float *qbm,
                  float *dt, float *dtc, float *ek, int *nx, int *ny,
                  int *idimp, int *npmax, int *nxv, int *nypmx,
                  int *idps, int *ntmax, int *ipbc);

void ppgrbpush23l_(float *part, float *fxy, float *bxy, float *edges,
                   int *npp, int *noff, int *ihole, float *qbm,
                   float *dt, float *dtc, float *ci, float *ek,
                   int *nx, int *ny, int *idimp, int *npmax, int *nxv,
                   int *nypmx, int *idps, int *ntmax, int *ipbc);

void ppgpost2l_(float *part, float *q, int *npp, int *noff, float *qm,
                int *idimp, int *npmax, int *nxv, int *nypmx);

void ppgjpost2l_(float *part, float *cu, float *edges, int *npp,
                 int *noff, int *ihole, float *qm, float *dt, int *nx,
                 int *ny, int *idimp, int *npmax, int *nxv, int *nypmx,
                 int *idps, int *ntmax, int *ipbc);

void ppgrjpost2l_(float *part, float *cu, float *edges, int *npp,
                  int *noff, int *ihole, float *qm, float *dt,
                  float *ci, int *nx, int *ny, int *idimp, int *npmax,
                  int *nxv, int *nypmx, int *idps, int *ntmax,
                  int *ipbc);

void ppdsortp2yl_(float *parta, float *partb, int *npic, int *npp,
                  int *noff, int *nyp, int *idimp, int *npmax,
                  int *nypm1);

void ppcguard2xl_(float *fxy, int *myp, int *nx, int *ndim, int *nxe,
                  int *nypmx);

void ppaguard2xl_(float *q, int *myp, int *nx, int *nxe, int *nypmx);

void ppacguard2xl_(float *cu, int *myp, int *nx, int *ndim, int *nxe,
                   int *nypmx);

void ppois23_(float complex *q, float complex *fxy, int *isign,
              float complex *ffc, float *ax, float *ay, float *affp,
              float *we, int *nx, int *ny, int *kstrt, int *nyv,
              int *kxp, int *nyhd);

void ppcuperp2_(float complex *cu, int *nx, int *ny, int *kstrt,
                int *nyv, int *kxp);

void ippbpoisp23_(float complex *cu, float complex *bxy,
                  float complex *ffc, float *ci, float *wm, int *nx,
                  int *ny, int *kstrt, int *nyv, int *kxp, int *nyhd);

void ppmaxwel2_(float complex *exy, float complex *bxy,
                float complex *cu, float complex *ffc, float *affp,
                float *ci, float *dt, float *wf, float *wm, int *nx,
                int *ny, int *kstrt, int *nyv, int *kxp, int *nyhd);

void ppemfield2_(float complex *fxy, float complex *exy,
                 float complex *ffc, int *isign, int *nx, int *ny,
                 int *kstrt, int *nyv, int *kxp, int *nyhd);

void wpfft2rinit_(int *mixup, float complex *sct, int *indx, int *indy,
                  int *nxhyd, int *nxyhd);

void ppfft2rxx_(float complex *f, int *isign, int *mixup,
                float complex *sct, int *indx, int *indy, int *kstrt,
                int *kypi, int *kypp, int *nxvh, int *kypd, int *nxhyd,
                int *nxyhd);

void ppfft2rxy_(float complex *g, int *isign, int *mixup,
                float complex *sct, int *indx, int *indy, int *kstrt,
                int *kxpi, int *kxpp, int *nyv, int *kxp, int *nxhyd,
                int *nxyhd);

void ppfft2r3xx_(float complex *f, int *isign, int *mixup,
                 float complex *sct, int *indx, int *indy, int *kstrt,
                 int *kypi, int *kypp, int *nxvh, int *kypd, int *nxhyd,
                 int *nxyhd);

void ppfft2r3xy_(float complex *g, int *isign, int *mixup,
                 float complex *sct, int *indx, int *indy, int *kstrt,
                 int *kxpi, int *kxpp, int *nyv, int *kxp, int *nxhyd,
                 int *nxyhd);

void wppfft2r_(float complex *f, float complex *g, float complex *bs,
               float complex *br, int *isign, int *ntpose, int *mixup,
               float complex *sct, float *ttp, int *indx, int *indy,
               int *kstrt, int *nvp, int *nxvh, int *nyv, int *kxp,
               int *kyp, int *kypd, int *nxhyd, int *nxyhd);

void wppfft2r3_(float complex *f, float complex *g, float complex *bs,
                float complex *br, int *isign, int *ntpose, int *mixup,
                float complex *sct, float *ttp, int *indx, int *indy,
                int *kstrt, int *nvp, int *nxvh, int *nyv, int *kxp,
                int *kyp, int *kypd, int *nxhyd, int *nxyhd);

/* Interfaces to C */

double ranorm() {
  return ranorm_();
}

/*--------------------------------------------------------------------*/
void cpdicomp2l(float edges[], int *nyp, int *noff, int *nypmx,
                int *nypmn, int ny, int kstrt, int nvp, int idps) {
   pdicomp2l_(edges,nyp,noff,nypmx,nypmn,&ny,&kstrt,&nvp,&idps);
   return;
}

/*--------------------------------------------------------------------*/
void cpdistr2h(float part[], float edges[], int *npp, int nps,
               float vtx, float vty, float vtz, float vdx, float vdy,
               float vdz, int npx, int npy, int nx, int ny, int idimp,
               int npmax, int idps, int ipbc, int *ierr) {
   pdistr2h_(part,edges,npp,&nps,&vtx,&vty,&vtz,&vdx,&vdy,&vdz,&npx,
             &npy,&nx,&ny,&idimp,&npmax,&idps,&ipbc,ierr);
   return;
}

/*--------------------------------------------------------------------*/
void cppgbpush23l(float part[], float fxy[], float bxy[], float edges[],
                  int npp, int noff, int ihole[], float qbm, float dt, 
                  float dtc, float *ek, int nx, int ny, int idimp,
                  int npmax, int nxv, int nypmx, int idps, int ntmax,
                  int ipbc) {
   ppgbpush23l_(part,fxy,bxy,edges,&npp,&noff,ihole,&qbm,&dt,&dtc,ek,
                &nx,&ny,&idimp,&npmax,&nxv,&nypmx,&idps,&ntmax,&ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cppgrbpush23l(float part[], float fxy[], float bxy[],
                   float edges[], int npp, int noff, int ihole[],
                   float qbm, float dt, float dtc, float ci, float *ek,
                   int nx, int ny, int idimp, int npmax, int nxv,
                   int nypmx, int idps, int ntmax, int ipbc) {
   ppgrbpush23l_(part,fxy,bxy,edges,&npp,&noff,ihole,&qbm,&dt,&dtc,&ci,
                 ek,&nx,&ny,&idimp,&npmax,&nxv,&nypmx,&idps,&ntmax,
                 &ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cppgpost2l(float part[], float q[], int npp, int noff, float qm,
                int idimp, int npmax, int nxv, int nypmx) {
   ppgpost2l_(part,q,&npp,&noff,&qm,&idimp,&npmax,&nxv,&nypmx);
   return;
}

/*--------------------------------------------------------------------*/
void cppgjpost2l(float part[], float cu[], float edges[], int npp,
                 int noff, int ihole[], float qm, float dt, int nx,
                 int ny, int idimp, int npmax, int nxv, int nypmx,
                 int idps, int ntmax, int ipbc) {
   ppgjpost2l_(part,cu,edges,&npp,&noff,ihole,&qm,&dt,&nx,&ny,&idimp,
               &npmax,&nxv,&nypmx,&idps,&ntmax,&ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cppgrjpost2l(float part[], float cu[], float edges[], int npp,
                  int noff, int ihole[], float qm, float dt, float ci,
                  int nx, int ny, int idimp, int npmax, int nxv,
                  int nypmx, int idps, int ntmax, int ipbc) {
   ppgrjpost2l_(part,cu,edges,&npp,&noff,ihole,&qm,&dt,&ci,&nx,&ny,
                &idimp,&npmax,&nxv,&nypmx,&idps,&ntmax,&ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cppdsortp2yl(float parta[], float partb[], int npic[], int npp,
                  int noff, int nyp, int idimp, int npmax, int nypm1) {
   ppdsortp2yl_(parta,partb,npic,&npp,&noff,&nyp,&idimp,&npmax,&nypm1);
   return;
}

/*--------------------------------------------------------------------*/
void cppcguard2xl(float fxy[], int myp, int nx, int ndim, int nxe,
                  int nypmx) {
   ppcguard2xl_(fxy,&myp,&nx,&ndim,&nxe,&nypmx);
   return;
}

/*--------------------------------------------------------------------*/
void cppaguard2xl(float q[], int myp, int nx, int nxe, int nypmx) {
   ppaguard2xl_(q,&myp,&nx,&nxe,&nypmx);
   return;
}

/*--------------------------------------------------------------------*/
void cppacguard2xl(float cu[], int myp, int nx, int ndim, int nxe,
                   int nypmx) {
   ppacguard2xl_(cu,&myp,&nx,&ndim,&nxe,&nypmx);
   return;
}

/*--------------------------------------------------------------------*/
void cppois23(float complex q[], float complex fxy[], int isign,
              float complex ffc[], float ax, float ay, float affp,
              float *we, int nx, int ny, int kstrt, int nyv, int kxp,
              int nyhd) {
   ppois23_(q,fxy,&isign,ffc,&ax,&ay,&affp,we,&nx,&ny,&kstrt,&nyv,&kxp,
            &nyhd);
   return;
}

/*--------------------------------------------------------------------*/
void cppcuperp2(float complex cu[], int nx, int ny, int kstrt, int nyv,
                int kxp) {
   ppcuperp2_(cu,&nx,&ny,&kstrt,&nyv,&kxp);
   return;
}

/*--------------------------------------------------------------------*/
void cippbpoisp23(float complex cu[], float complex bxy[],
                  float complex ffc[], float ci, float *wm, int nx,
                  int ny, int kstrt, int nyv, int kxp, int nyhd) {
   ippbpoisp23_(cu,bxy,ffc,&ci,wm,&nx,&ny,&kstrt,&nyv,&kxp,&nyhd);
   return;
}

/*--------------------------------------------------------------------*/
void cppmaxwel2(float complex exy[], float complex bxy[],
                float complex cu[], float complex ffc[], float affp,
                float ci, float dt, float *wf, float *wm, int nx,
                int ny, int kstrt, int nyv, int kxp, int nyhd) {
   ppmaxwel2_(exy,bxy,cu,ffc,&affp,&ci,&dt,wf,wm,&nx,&ny,&kstrt,&nyv,
              &kxp,&nyhd);
   return;
}

/*--------------------------------------------------------------------*/
void cppemfield2(float complex fxy[], float complex exy[],
                 float complex ffc[], int isign, int nx, int ny,
                 int kstrt, int nyv, int kxp, int nyhd) {
   ppemfield2_(fxy,exy,ffc,&isign,&nx,&ny,&kstrt,&nyv,&kxp,&nyhd);
   return;
}

/*--------------------------------------------------------------------*/
void cwpfft2rinit(int mixup[], float complex sct[], int indx, int indy,
                  int nxhyd, int nxyhd) {
   wpfft2rinit_(mixup,sct,&indx,&indy,&nxhyd,&nxyhd);
   return;
}

/*--------------------------------------------------------------------*/
void cppfft2rxx(float complex f[], int isign, int mixup[],
                float complex sct[], int indx, int indy, int kstrt,
                int kypi, int kypp, int nxvh, int kypd, int nxhyd,
                int nxyhd) {
   ppfft2rxx_(f,&isign,mixup,sct,&indx,&indy,&kstrt,&kypi,&kypp,&nxvh,
              &kypd,&nxhyd,&nxyhd);
   return;
}

/*--------------------------------------------------------------------*/
void cppfft2rxy(float complex g[], int isign, int mixup[],
                float complex sct[], int indx, int indy, int kstrt,
                int kxpi, int kxpp, int nyv, int kxp, int nxhyd,
                int nxyhd) {
   ppfft2rxy_(g,&isign,mixup,sct,&indx,&indy,&kstrt,&kxpi,&kxpp,&nyv,
              &kxp,&nxhyd,&nxyhd);
   return;
}

/*--------------------------------------------------------------------*/
void cppfft2r3xx(float complex f[], int isign, int mixup[],
                 float complex sct[], int indx, int indy, int kstrt,
                 int kypi, int kypp, int nxvh, int kypd, int nxhyd,
                 int nxyhd) {
   ppfft2r3xx_(f,&isign,mixup,sct,&indx,&indy,&kstrt,&kypi,&kypp,&nxvh,
               &kypd,&nxhyd,&nxyhd);
   return;
}

/*--------------------------------------------------------------------*/
void cppfft2r3xy(float complex g[], int isign, int mixup[],
                 float complex sct[], int indx, int indy, int kstrt,
                 int kxpi, int kxpp, int nyv, int kxp, int nxhyd,
                 int nxyhd) {
   ppfft2r3xy_(g,&isign,mixup,sct,&indx,&indy,&kstrt,&kxpi,&kxpp,&nyv,
               &kxp,&nxhyd,&nxyhd);
   return;
}

/*--------------------------------------------------------------------*/
void cwppfft2r(float complex f[], float complex g[], float complex bs[],
               float complex br[], int isign, int ntpose, int mixup[],
               float complex sct[], float *ttp, int indx, int indy,
               int kstrt, int nvp, int nxvh, int nyv, int kxp, int kyp,
               int kypd, int nxhyd, int nxyhd) {
   wppfft2r_(f,g,bs,br,&isign,&ntpose,mixup,sct,ttp,&indx,&indy,&kstrt,
             &nvp,&nxvh,&nyv,&kxp,&kyp,&kypd,&nxhyd,&nxyhd);
   return;
}

/*--------------------------------------------------------------------*/
void cwppfft2r3(float complex f[], float complex g[], float complex bs[],
                float complex br[], int isign, int ntpose, int mixup[],
                float complex sct[], float *ttp, int indx, int indy,
                int kstrt, int nvp, int nxvh, int nyv, int kxp, int kyp,
                int kypd, int nxhyd, int nxyhd) {
   wppfft2r3_(f,g,bs,br,&isign,&ntpose,mixup,sct,ttp,&indx,&indy,&kstrt,
              &nvp,&nxvh,&nyv,&kxp,&kyp,&kypd,&nxhyd,&nxyhd);
   return;
}

