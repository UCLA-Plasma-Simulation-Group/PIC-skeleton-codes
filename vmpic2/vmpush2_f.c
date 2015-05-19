/* C Library for Skeleton 2D Electrostatic OpenMP/Vector PIC Code */
/* Wrappers for calling the Fortran routines from a C main program */

#include <complex.h>

double ranorm_();

/*--------------------------------------------------------------------*/
void distr2_(float *part, float *vtx, float *vty, float *vdx,
             float *vdy, int *npx, int *npy, int *idimp, int *nop,
             int *nx, int *ny, int *ipbc);

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
                 int *nx, int *ny, int *mx, int *my, int *mx1, int *my1,
                 int *irc);

/*--------------------------------------------------------------------*/
void gppush2lt_(float *ppart, float *fxy, int *kpic, float *qbm,
                float *dt, float *ek, int *idimp, int *nppmx, int *nx,
                int *ny, int *mx, int *my, int *nxv, int *nyv, int *mx1,
                int *mxy1, int *ipbc);

/*--------------------------------------------------------------------*/
void gppushf2lt_(float *ppart, float *fxy, int *kpic, int *ncl,
                 int *ihole, float *qbm, float *dt, float *ek,
                 int *idimp, int *nppmx, int *nx, int *ny, int *mx,
                 int *my, int *nxv, int *nyv, int *mx1, int *mxy1,
                 int *ntmax, int *irc);

/*--------------------------------------------------------------------*/
void vgppush2lt_(float *ppart, float *fxy, int *kpic, float *qbm,
                 float *dt, float *ek, int *idimp, int *nppmx, int *nx,
                 int *ny, int *mx, int *my, int *nxv, int *nyv,
                 int *mx1, int *mxy1, int *ipbc);

/*--------------------------------------------------------------------*/
void vgppushf2lt_(float *ppart, float *fxy, int *kpic, int *ncl,
                  int *ihole, float *qbm, float *dt, float *ek,
                  int *idimp, int *nppmx, int *nx, int *ny, int *mx,
                  int *my, int *nxv, int *nyv, int *mx1, int *mxy1,
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
void viscan2_(int *isdata, int *mb, int *nths);

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
void cguard2l_(float *fxy, int *nx, int *ny, int *nxe, int *nye);

/*--------------------------------------------------------------------*/
void aguard2l_(float *q, int *nx, int *ny, int *nxe, int *nye);

/*--------------------------------------------------------------------*/
void vmpois22_(float complex *q, float complex *fxy, int *isign,
               float complex *ffc, float *ax, float *ay, float *affp,
               float *we, int *nx, int *ny, int *nxvh, int *nyv,
               int *nxhd, int *nyhd);

/*--------------------------------------------------------------------*/
void wfft2rinit_(int *mixup, float complex *sct, int *indx, int *indy,
                 int *nxhyd, int *nxyhd);

/*--------------------------------------------------------------------*/
void fft2rmxx_(float complex *f, int *isign, int *mixup,
               float complex *sct, int *indx, int *indy, int *nyi,
               int *nyp, int *nxhd, int *nyd, int *nxhyd, int *nxyhd);

/*--------------------------------------------------------------------*/
void fft2rmxy_(float complex *f, int *isign, int *mixup,
               float complex *sct, int *indx, int *indy, int *nxi,
               int *nxp, int *nxhd, int *nyd, int *nxhyd, int *nxyhd);

/*--------------------------------------------------------------------*/
void fft2rm2x_(float complex *f, int *isign, int *mixup,
               float complex *sct, int *indx, int *indy, int *nyi,
               int *nyp, int *nxhd, int *nyd, int *nxhyd, int *nxyhd);

/*--------------------------------------------------------------------*/
void fft2rm2y_(float complex *f, int *isign, int *mixup,
               float complex *sct, int *indx, int *indy, int *nxi,
               int *nxp, int *nxhd, int *nyd, int *nxhyd, int *nxyhd);

/*--------------------------------------------------------------------*/
void wfft2rmx_(float complex *f, int *isign, int *mixup,
               float complex *sct, int *indx, int *indy, int *nxhd,
               int *nyd, int *nxhyd, int *nxyhd);

/*--------------------------------------------------------------------*/
void wfft2rm2_(float complex *f, int *isign, int *mixup,
               float complex *sct, int *indx, int *indy, int *nxhd,
               int *nyd, int *nxhyd, int *nxyhd);

/* Interfaces to C */

double ranorm() {
  return ranorm_();
}

/*--------------------------------------------------------------------*/
void cdistr2(float part[], float vtx, float vty, float vdx, float vdy,
             int npx, int npy, int idimp, int nop, int nx, int ny,
             int ipbc) {
   distr2_(part,&vtx,&vty,&vdx,&vdy,&npx,&npy,&idimp,&nop,&nx,&ny,
           &ipbc);
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
void cgppush2lt(float ppart[], float fxy[], int kpic[], float qbm,
                float dt, float *ek, int idimp, int nppmx, int nx,
                int ny, int mx, int my, int nxv, int nyv, int mx1,
                int mxy1, int ipbc) {
   gppush2lt_(ppart,fxy,kpic,&qbm,&dt,ek,&idimp,&nppmx,&nx,&ny,&mx,&my,
              &nxv,&nyv,&mx1,&mxy1,&ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cgppushf2lt(float ppart[], float fxy[], int kpic[], int ncl[],
                 int ihole[], float qbm, float dt, float *ek, int idimp,
                 int nppmx, int nx, int ny, int mx, int my, int nxv,
                 int nyv, int mx1, int mxy1, int ntmax, int *irc) {
   gppushf2lt_(ppart,fxy,kpic,ncl,ihole,&qbm,&dt,ek,&idimp,&nppmx,&nx,
               &ny,&mx,&my,&nxv,&nyv,&mx1,&mxy1,&ntmax,irc);
   return;
}

/*--------------------------------------------------------------------*/
void cvgppush2lt(float ppart[], float fxy[], int kpic[], float qbm,
                 float dt, float *ek, int idimp, int nppmx, int nx,
                 int ny, int mx, int my, int nxv, int nyv, int mx1,
                 int mxy1, int ipbc) {
   vgppush2lt_(ppart,fxy,kpic,&qbm,&dt,ek,&idimp,&nppmx,&nx,&ny,&mx,&my,
               &nxv,&nyv,&mx1,&mxy1,&ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cvgppushf2lt(float ppart[], float fxy[], int kpic[], int ncl[],
                  int ihole[], float qbm, float dt, float *ek,
                  int idimp, int nppmx, int nx, int ny, int mx, int my,
                  int nxv, int nyv, int mx1, int mxy1, int ntmax,
                  int *irc) {
   vgppushf2lt_(ppart,fxy,kpic,ncl,ihole,&qbm,&dt,ek,&idimp,&nppmx,&nx,
                &ny,&mx,&my,&nxv,&nyv,&mx1,&mxy1,&ntmax,irc);
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
void cviscan2(int *isdata, int *mb, int nths) {
   viscan2_(isdata,mb,&nths);
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
void ccguard2l(float fxy[], int nx, int ny, int nxe, int nye) {
   cguard2l_(fxy,&nx,&ny,&nxe,&nye);
   return;
}

/*--------------------------------------------------------------------*/
void caguard2l(float q[], int nx, int ny, int nxe, int nye) {
   aguard2l_(q,&nx,&ny,&nxe,&nye);
   return;
}

/*--------------------------------------------------------------------*/
void cvmpois22(float complex q[], float complex fxy[], int isign,
               float complex ffc[], float ax, float ay, float affp,
               float *we, int nx, int ny, int nxvh, int nyv, int nxhd,
               int nyhd) {
   vmpois22_(q,fxy,&isign,ffc,&ax,&ay,&affp,we,&nx,&ny,&nxvh,&nyv,&nxhd,
             &nyhd);
   return;
}

/*--------------------------------------------------------------------*/
void cwfft2rinit(int mixup[], float complex sct[], int indx, int indy, 
                 int nxhyd, int nxyhd) {
   wfft2rinit_(mixup,sct,&indx,&indy,&nxhyd,&nxyhd);
   return;
}

/*--------------------------------------------------------------------*/
void cfft2rmxx(float complex f[], int isign, int mixup[],
               float complex sct[], int indx, int indy, int nyi,
               int nyp, int nxhd, int nyd, int nxhyd, int nxyhd) {
   fft2rmxx_(f,&isign,mixup,sct,&indx,&indy,&nyi,&nyp,&nxhd,&nyd,&nxhyd,
             &nxyhd);
   return;
}

/*--------------------------------------------------------------------*/
void cfft2rmxy(float complex f[], int isign, int mixup[],
               float complex sct[], int indx, int indy, int nxi,
               int nxp, int nxhd, int nyd, int nxhyd, int nxyhd) {
   fft2rmxy_(f,&isign,mixup,sct,&indx,&indy,&nxi,&nxp,&nxhd,&nyd,&nxhyd,
             &nxyhd);
   return;
}

/*--------------------------------------------------------------------*/
void cfft2rm2x(float complex f[], int isign, int mixup[],
               float complex sct[], int indx, int indy, int nyi,
               int nyp, int nxhd, int nyd, int nxhyd, int nxyhd) {
   fft2rm2x_(f,&isign,mixup,sct,&indx,&indy,&nyi,&nyp,&nxhd,&nyd,&nxhyd,
             &nxyhd);
   return;
}

/*--------------------------------------------------------------------*/
void cfft2rm2y(float complex f[], int isign, int mixup[],
               float complex sct[], int indx, int indy, int nxi,
               int nxp, int nxhd, int nyd, int nxhyd, int nxyhd) {
   fft2rm2y_(f,&isign,mixup,sct,&indx,&indy,&nxi,&nxp,&nxhd,&nyd,&nxhyd,
             &nxyhd);
   return;
}

/*--------------------------------------------------------------------*/
void cwfft2rmx(float complex f[], int isign, int mixup[],
               float complex sct[], int indx, int indy, int nxhd,
               int nyd, int nxhyd, int nxyhd) {
   wfft2rmx_(f,&isign,mixup,sct,&indx,&indy,&nxhd,&nyd,&nxhyd,&nxyhd);
   return;
}

/*--------------------------------------------------------------------*/
void cwfft2rm2(float complex f[], int isign, int mixup[],
               float complex sct[], int indx, int indy, int nxhd,
               int nyd, int nxhyd, int nxyhd) {
   wfft2rm2_(f,&isign,mixup,sct,&indx,&indy,&nxhd,&nyd,&nxhyd,&nxyhd);
   return;
}
