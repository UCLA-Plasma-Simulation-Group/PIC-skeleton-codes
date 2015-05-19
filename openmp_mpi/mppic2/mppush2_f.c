/* C Library for Skeleton 2D Electrostatic MPI/OpenMP PIC Code */
/* Wrappers for calling the Fortran routines from a C main program */

#include <complex.h>

double ranorm_();

void pdicomp2l_(float *edges, int *nyp, int *noff, int *nypmx,
                int *nypmn, int *ny, int *kstrt, int *nvp, int *idps);

void pdistr2_(float *part, float *edges, int *npp, int *nps, float *vtx,
              float *vty, float *vdx, float *vdy, int *npx, int *npy,
              int *nx, int *ny, int *idimp, int *npmax, int *idps,
              int *ipbc, int *ierr);

void ppdblkp2l_(float *part, int *kpic, int *npp, int *noff, int *nppmx,
                int *idimp, int *npmax, int *mx, int *my, int *mx1,
                int *mxyp1, int *irc);

void pppmovin2l_(float *part, float *ppart, int *kpic, int *npp,
                 int *noff, int *nppmx, int *idimp, int *npmax, int *mx,
                 int *my, int *mx1, int *mxyp1, int *irc);

void pppcheck2l_(float *ppart, int *kpic, int *noff, int *nyp,
                 int *idimp, int *nppmx, int *nx, int *mx, int *my,
                 int *mx1, int *myp1, int *irc);

void ppgppush2l_(float *ppart, float *fxy, int *kpic, int *noff,
                 int *nyp, float *qbm, float *dt, float *ek, int *nx,
                 int *ny, int *mx, int *my, int *idimp, int *nppmx,
                 int *nxv, int *nypmx, int *mx1, int *mxyp1, 
                 int *ipbc);

void ppgppushf2l_(float *ppart, float *fxy, int *kpic, int *ncl,
                  int *ihole, int *noff, int *nyp, float *qbm,
                  float *dt, float *ek, int *nx, int *ny, int *mx,
                  int *my, int *idimp, int *nppmx, int *nxv, int *nypmx,
                  int *mx1, int *mxyp1, int *ntmax, int *irc);

void ppgppost2l_(float *ppart, float *q, int *kpic, int *noff,
                 float *qm, int *idimp, int *nppmx, int *mx, int *my,
                 int *nxv, int *nypmx, int *mx1, int *mxyp1);

void ppporder2la_(float *ppart, float *ppbuff, float *sbufl,
                  float *sbufr, int *kpic, int *ncl, int *ihole,
                  int *ncll, int *nclr, int *noff, int *nyp,
                  int *idimp, int *nppmx, int *nx, int *ny, int *mx,
                  int *my, int *mx1, int *myp1, int *npbmx, int *ntmax,
                  int *nbmax, int *irc);

void ppporderf2la_(float *ppart, float *ppbuff, float *sbufl,
                   float *sbufr, int *ncl, int *ihole, int *ncll,
                   int *nclr, int *idimp, int *nppmx, int *mx1,
                   int *myp1, int *npbmx, int *ntmax, int *nbmax,
                   int *irc);

void ppporder2lb_(float *ppart, float *ppbuff, float *rbufl,
                  float *rbufr, int *kpic, int *ncl, int *ihole,
                  int *mcll, int *mclr, int *idimp, int *nppmx,
                  int *mx1, int *myp1, int *npbmx, int *ntmax,
                  int *nbmax, int *irc);

void ppcguard2xl_(float *fxy, int *nyp, int *nx, int *ndim, int *nxe,
                  int *nypmx);

void ppaguard2xl_(float *q, int *nyp, int *nx, int *nxe, int *nypmx);

void mppois22_(float complex *q, float complex *fxy, int *isign,
               float complex *ffc, float *ax, float *ay, float *affp,
               float *we, int *nx, int *ny, int *kstrt, int *nyv,
               int *kxp, int *nyhd);

void wpfft2rinit_(int *mixup, float complex *sct, int *indx, int *indy,
                  int *nxhyd, int *nxyhd);

void ppfft2rmxx_(float complex *f, int *isign, int *mixup,
                 float complex *sct, int *indx, int *indy, int *kstrt,
                 int *kypi, int *kypp, int *nxvh, int *kypd, int *nxhyd,
                 int *nxyhd);

void ppfft2rmxy_(float complex *g, int *isign, int *mixup,
                 float complex *sct, int *indx, int *indy, int *kstrt,
                 int *kxpi, int *kxpp, int *nyv, int *kxp, int *nxhyd,
                 int *nxyhd);

void ppfft2rm2xx_(float complex *f, int *isign, int *mixup,
                  float complex *sct, int *indx, int *indy, int *kstrt,
                  int *kypi, int *kypp, int *nxvh, int *kypd,
                  int *nxhyd, int *nxyhd);

void ppfft2rm2xy_(float complex *g, int *isign, int *mixup,
                  float complex *sct, int *indx, int *indy, int *kstrt,
                  int *kxpi, int *kxpp, int *nyv, int *kxp, int *nxhyd,
                  int *nxyhd);

void wppfft2rm_(float complex *f, float complex *g, float complex *bs,
                float complex *br, int *isign, int *ntpose, int *mixup,
                float complex *sct, float *ttp, int *indx, int *indy,
                int *kstrt, int *nvp, int *nxvh, int *nyv, int *kxp,
                int *kyp, int *kypd, int *nxhyd, int *nxyhd);

void wppfft2rm2_(float complex *f, float complex *g, float complex *bs,
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
              int ny, int idimp, int npmax, int idps, int ipbc,
              int *ierr) {
   pdistr2_(part,edges,npp,&nps,&vtx,&vty,&vdx,&vdy,&npx,&npy,&nx,&ny,
            &idimp,&npmax,&idps,&ipbc,ierr);
   return;
}

/*--------------------------------------------------------------------*/
void cppdblkp2l(float part[], int kpic[], int npp, int noff, int *nppmx,
                int idimp, int npmax, int mx, int my, int mx1,
                int mxyp1, int *irc) {
   ppdblkp2l_(part,kpic,&npp,&noff,nppmx,&idimp,&npmax,&mx,&my,&mx1,
              &mxyp1,irc);
   return;
}

/*--------------------------------------------------------------------*/
void cpppmovin2l(float part[], float ppart[], int kpic[], int npp,
                 int noff, int nppmx, int idimp, int npmax, int mx,
                 int my, int mx1, int mxyp1, int *irc) {
   pppmovin2l_(part,ppart,kpic,&npp,&noff,&nppmx,&idimp,&npmax,&mx,&my,
               &mx1,&mxyp1,irc);
   return;
}

/*--------------------------------------------------------------------*/
void cpppcheck2l(float ppart[], int kpic[], int noff, int nyp,
                 int idimp, int nppmx, int nx, int mx, int my, int mx1,
                 int myp1, int *irc) {
   pppcheck2l_(ppart,kpic,&noff,&nyp,&idimp,&nppmx,&nx,&mx,&my,&mx1,
               &myp1,irc);
   return;
}

/*--------------------------------------------------------------------*/
void cppgppush2l(float ppart[], float fxy[], int kpic[], int noff,
                 int nyp, float qbm, float dt, float *ek, int nx,
                 int ny, int mx, int my, int idimp, int nppmx, int nxv,
                 int nypmx, int mx1, int mxyp1, int ipbc) {
   ppgppush2l_(ppart,fxy,kpic,&noff,&nyp,&qbm,&dt,ek,&nx,&ny,&mx,&my,
               &idimp,&nppmx,&nxv,&nypmx,&mx1,&mxyp1,&ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cppgppushf2l(float ppart[], float fxy[], int kpic[], int ncl[],
                  int ihole[], int noff, int nyp, float qbm, float dt,
                  float *ek, int nx, int ny, int mx, int my, int idimp,
                  int nppmx, int nxv, int nypmx, int mx1, int mxyp1,
                  int ntmax, int *irc) {
   ppgppushf2l_(ppart,fxy,kpic,ncl,ihole,&noff,&nyp,&qbm,&dt,ek,&nx,&ny,
                &mx,&my,&idimp,&nppmx,&nxv,&nypmx,&mx1,&mxyp1,&ntmax,
                irc);
   return;
}

/*--------------------------------------------------------------------*/
void cppgppost2l(float ppart[], float q[], int kpic[], int noff, 
                 float qm, int idimp, int nppmx, int mx, int my,
                 int nxv, int nypmx, int mx1, int mxyp1) {
   ppgppost2l_(ppart,q,kpic,&noff,&qm,&idimp,&nppmx,&mx,&my,&nxv,
               &nypmx,&mx1,&mxyp1);
   return;
}

/*--------------------------------------------------------------------*/
void cppporder2la(float ppart[], float ppbuff[], float sbufl[],
                  float sbufr[], int kpic[], int ncl[], int ihole[],
                  int ncll[], int nclr[], int noff, int nyp, int idimp,
                  int nppmx, int nx, int ny, int mx, int my, int mx1,
                  int myp1, int npbmx, int ntmax, int nbmax, int *irc) {
   ppporder2la_(ppart,ppbuff,sbufl,sbufr,kpic,ncl,ihole,ncll,nclr,&noff,
                &nyp,&idimp,&nppmx,&nx,&ny,&mx,&my,&mx1,&myp1,&npbmx,
                &ntmax,&nbmax,irc);
   return;
}

/*--------------------------------------------------------------------*/
void cppporderf2la(float ppart[], float ppbuff[], float sbufl[],
                   float sbufr[], int ncl[], int ihole[], int ncll[],
                   int nclr[], int idimp, int nppmx, int mx1, int myp1,
                   int npbmx, int ntmax, int nbmax, int *irc) {
   ppporderf2la_(ppart,ppbuff,sbufl,sbufr,ncl,ihole,ncll,nclr,&idimp,
                 &nppmx,&mx1,&myp1,&npbmx,&ntmax,&nbmax,irc);
   return;
}

/*--------------------------------------------------------------------*/
void cppporder2lb(float ppart[], float ppbuff[], float rbufl[],
                  float rbufr[], int kpic[], int ncl[], int ihole[],
                  int mcll[], int mclr[], int idimp, int nppmx, int mx1,
                  int myp1, int npbmx, int ntmax, int nbmax, int *irc) {
   ppporder2lb_(ppart,ppbuff,rbufl,rbufr,kpic,ncl,ihole,mcll,mclr,
                &idimp,&nppmx,&mx1,&myp1,&npbmx,&ntmax,&nbmax,irc);
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
void cmppois22(float complex q[], float complex fxy[], int isign,
               float complex ffc[], float ax, float ay, float affp,
               float *we, int nx, int ny, int kstrt, int nyv, int kxp,
               int nyhd) {
   mppois22_(q,fxy,&isign,ffc,&ax,&ay,&affp,we,&nx,&ny,&kstrt,&nyv,&kxp,
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
void cppfft2rmxx(float complex f[], int isign, int mixup[],
                 float complex sct[], int indx, int indy, int kstrt,
                 int kypi, int kypp, int nxvh, int kypd, int nxhyd,
                 int nxyhd) {
   ppfft2rmxx_(f,&isign,mixup,sct,&indx,&indy,&kstrt,&kypi,&kypp,&nxvh,
               &kypd,&nxhyd,&nxyhd);
   return;
}

/*--------------------------------------------------------------------*/
void cppfft2rmxy(float complex g[], int isign, int mixup[],
                 float complex sct[], int indx, int indy, int kstrt,
                 int kxpi, int kxpp, int nyv, int kxp, int nxhyd,
                 int nxyhd) {
   ppfft2rmxy_(g,&isign,mixup,sct,&indx,&indy,&kstrt,&kxpi,&kxpp,&nyv,
               &kxp,&nxhyd,&nxyhd);
   return;
}

/*--------------------------------------------------------------------*/
void cppfft2rm2xx(float complex f[], int isign, int mixup[],
                  float complex sct[], int indx, int indy, int kstrt,
                  int kypi, int kypp, int nxvh, int kypd, int nxhyd,
                  int nxyhd) {
   ppfft2rm2xx_(f,&isign,mixup,sct,&indx,&indy,&kstrt,&kypi,&kypp,&nxvh,
                &kypd,&nxhyd,&nxyhd);
   return;
}

/*--------------------------------------------------------------------*/
void cppfft2rm2xy(float complex g[], int isign, int mixup[],
                  float complex sct[], int indx, int indy, int kstrt,
                  int kxpi, int kxpp, int nyv, int kxp, int nxhyd,
                  int nxyhd) {
   ppfft2rm2xy_(g,&isign,mixup,sct,&indx,&indy,&kstrt,&kxpi,&kxpp,&nyv,
                &kxp,&nxhyd,&nxyhd);
   return;
}

/*--------------------------------------------------------------------*/
void cwppfft2rm(float complex f[], float complex g[],
                float complex bs[], float complex br[], int isign,
                int ntpose, int mixup[], float complex sct[],
                float *ttp, int indx, int indy, int kstrt, int nvp,
                int nxvh, int nyv, int kxp, int kyp, int kypd,
                int nxhyd, int nxyhd) {
   wppfft2rm_(f,g,bs,br,&isign,&ntpose,mixup,sct,ttp,&indx,&indy,&kstrt,
              &nvp,&nxvh,&nyv,&kxp,&kyp,&kypd,&nxhyd,&nxyhd);
   return;
}

/*--------------------------------------------------------------------*/
void cwppfft2rm2(float complex f[], float complex g[],
                 float complex bs[], float complex br[], int isign,
                 int ntpose, int mixup[], float complex sct[],
                 float *ttp, int indx, int indy, int kstrt, int nvp,
                 int nxvh, int nyv, int kxp, int kyp, int kypd,
                 int nxhyd, int nxyhd) {
   wppfft2rm2_(f,g,bs,br,&isign,&ntpose,mixup,sct,ttp,&indx,&indy,&kstrt,
               &nvp,&nxvh,&nyv,&kxp,&kyp,&kypd,&nxhyd,&nxyhd);
   return;
}
