/* C Library for Skeleton 2D Electrostatic MPI PIC Code */
/* Wrappers for calling the Fortran routines from a C main program */

#include <complex.h>

double ranorm_();

void pdicomp2l_(float *edges, int *nyp, int *noff, int *nypmx,
                int *nypmn, int *ny, int *kstrt, int *nvp, int *idps);

void pdistr2_(float *part, float *edges, int *npp, int *nps, float *vtx,
              float *vty, float *vdx, float *vdy, int *npx, int *npy,
              int *nx, int *ny, int *idimp, int *npmax, int *idps,
              int *ipbc, int *ierr);

void ppgpush2l_(float *part, float *fxy, float *edges, int *npp,
                int *noff, int *ihole, float *qbm, float *dt, float *ek,
                int *nx, int *ny, int *idimp, int *npmax, int *nxv,
                int *nypmx, int *idps, int *ntmax, int *ipbc);

void ppgpost2l_(float *part, float *q, int *npp, int *noff, float *qm,
                int *idimp, int *npmax, int *nxv, int *nypmx);

void ppdsortp2yl_(float *parta, float *partb, int *npic, int *npp,
                  int *noff, int *nyp, int *idimp, int *npmax,
                  int *nypm1);

void ppcguard2xl_(float *fxy, int *nyp, int *nx, int *ndim, int *nxe,
                  int *nypmx);

void ppaguard2xl_(float *q, int *nyp, int *nx, int *nxe, int *nypmx);

void ppois22_(float complex *q, float complex *fxy, int *isign,
              float complex *ffc, float *ax, float *ay, float *affp,
              float *we, int *nx, int *ny, int *kstrt, int *nyv,
              int *kxp, int *nyhd);

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

void ppfft2r2xx_(float complex *f, int *isign, int *mixup,
                 float complex *sct, int *indx, int *indy, int *kstrt,
                 int *kypi, int *kypp, int *nxvh, int *kypd, int *nxhyd,
                 int *nxyhd);

void ppfft2r2xy_(float complex *g, int *isign, int *mixup,
                 float complex *sct, int *indx, int *indy, int *kstrt,
                 int *kxpi, int *kxpp, int *nyv, int *kxp, int *nxhyd,
                 int *nxyhd);

void wppfft2r_(float complex *f, float complex *g, float complex *bs,
               float complex *br, int *isign, int *ntpose, int *mixup,
               float complex *sct, float *ttp, int *indx, int *indy,
               int *kstrt, int *nvp, int *nxvh, int *nyv, int *kxp,
               int *kyp, int *kypd, int *nxhyd, int *nxyhd);

void wppfft2r2_(float complex *f, float complex *g, float complex *bs,
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
void cpdistr2(float part[], float edges[], int *npp, int nps, float vtx,
              float vty, float vdx, float vdy, int npx, int npy, int nx,
              int ny, int idimp, int npmax, int idps, int ipbc, int *ierr) {
   pdistr2_(part,edges,npp,&nps,&vtx,&vty,&vdx,&vdy,&npx,&npy,&nx,&ny,
            &idimp,&npmax,&idps,&ipbc,ierr);
   return;
}

/*--------------------------------------------------------------------*/
void cppgpush2l(float part[], float fxy[], float edges[], int npp,
                int noff, int ihole[], float qbm, float dt, float *ek,
                int nx, int ny, int idimp, int npmax, int nxv,
                int nypmx, int idps, int ntmax, int ipbc) {
   ppgpush2l_(part,fxy,edges,&npp,&noff,ihole,&qbm,&dt,ek,&nx,&ny,
              &idimp,&npmax,&nxv,&nypmx,&idps,&ntmax,&ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cppgpost2l(float part[], float q[], int npp, int noff, float qm,
                int idimp, int npmax, int nxv, int nypmx) {
   ppgpost2l_(part,q,&npp,&noff,&qm,&idimp,&npmax,&nxv,&nypmx);
   return;
}

/*--------------------------------------------------------------------*/
void cppdsortp2yl(float parta[], float partb[], int npic[], int npp,
                  int noff, int nyp, int idimp, int npmax, int nypm1) {
   ppdsortp2yl_(parta,partb,npic,&npp,&noff,&nyp,&idimp,&npmax,&nypm1);
   return;
}

/*--------------------------------------------------------------------*/
void cppcguard2xl(float fxy[], int nyp, int nx, int ndim, int nxe,
                  int nypmx) {
   ppcguard2xl_(fxy,&nyp,&nx,&ndim,&nxe,&nypmx);
   return;
}

/*--------------------------------------------------------------------*/
void cppaguard2xl(float q[], int nyp, int nx, int nxe, int nypmx) {
   ppaguard2xl_(q,&nyp,&nx,&nxe,&nypmx);
   return;
}

/*--------------------------------------------------------------------*/
void cppois22(float complex q[], float complex fxy[], int isign,
              float complex ffc[], float ax, float ay, float affp,
              float *we, int nx, int ny, int kstrt, int nyv, int kxp,
              int nyhd) {
   ppois22_(q,fxy,&isign,ffc,&ax,&ay,&affp,we,&nx,&ny,&kstrt,&nyv,&kxp,
            &nyhd);
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
void cppfft2r2xx(float complex f[], int isign, int mixup[],
                 float complex sct[], int indx, int indy, int kstrt,
                 int kypi, int kypp, int nxvh, int kypd, int nxhyd,
                 int nxyhd) {
   ppfft2r2xx_(f,&isign,mixup,sct,&indx,&indy,&kstrt,&kypi,&kypp,&nxvh,
               &kypd,&nxhyd,&nxyhd);
   return;
}

/*--------------------------------------------------------------------*/
void cppfft2r2xy(float complex g[], int isign, int mixup[],
                 float complex sct[], int indx, int indy, int kstrt,
                 int kxpi, int kxpp, int nyv, int kxp, int nxhyd,
                 int nxyhd) {
   ppfft2r2xy_(g,&isign,mixup,sct,&indx,&indy,&kstrt,&kxpi,&kxpp,&nyv,
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
void cwppfft2r2(float complex f[], float complex g[], float complex bs[],
                float complex br[], int isign, int ntpose, int mixup[],
                float complex sct[], float *ttp, int indx, int indy,
                int kstrt, int nvp, int nxvh, int nyv, int kxp, int kyp,
                int kypd, int nxhyd, int nxyhd) {
   wppfft2r2_(f,g,bs,br,&isign,&ntpose,mixup,sct,ttp,&indx,&indy,&kstrt,
              &nvp,&nxvh,&nyv,&kxp,&kyp,&kypd,&nxhyd,&nxyhd);
   return;
}
