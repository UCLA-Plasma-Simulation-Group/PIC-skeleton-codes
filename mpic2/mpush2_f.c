/* C Library for Skeleton 2D Electrostatic OpenMP PIC Code */
/* Wrappers for calling the Fortran routines from a C main program */

#include <complex.h>

double ranorm_();

/* Interfaces to Fortran */

void distr2_(float *part, float *vtx, float *vty, float *vdx, float *vdy,
             int *npx, int *npy, int *idimp, int *nop, int *nx, int *ny,
             int *ipbc);

void dblkp2l_(float *part, int *kpic, int *nppmx, int *idimp, int *nop,
              int *mx, int *my, int *mx1, int *mxy1, int *irc);

void ppmovin2l_(float *part, float *ppart, int *kpic, int *nppmx,
                int *idimp, int *nop, int *mx, int *my, int *mx1,
                int *mxy1, int *irc);

void ppcheck2l_(float *ppart, int *kpic, int *idimp, int *nppmx,
                int *nx, int *ny, int *mx, int *my, int *mx1, int *my1, 
                int *irc);

void gppush2l_(float *ppart, float *fxy, int *kpic, float *qbm,
               float *dt, float *ek, int *idimp, int *nppmx, int *nx,
               int *ny, int *mx, int *my, int *nxv, int *nyv, int *mx1,
               int *mxy1, int *ipbc);

void gppushf2l_(float *ppart, float *fxy, int *kpic, int *ncl,
                int *ihole, float *qbm, float *dt, float *ek,
                int *idimp, int *nppmx, int *nx, int *ny, int *mx,
                int *my, int *nxv, int *nyv, int *mx1, int *mxy1,
                int *ntmax, int *irc);

void gppost2l_(float *ppart, float *q, int *kpic, float *qm,
               int *nppmx, int *idimp, int *mx, int *my, int *nxv,
               int *nyv, int *mx1, int *mxy1);

void cguard2l_(float *fxy, int *nx, int *ny, int *nxe, int *nye);

void aguard2l_(float *q, int *nx, int *ny, int *nxe, int *nye);

void pporder2l_(float *ppart, float *ppbuff, int *kpic, int *ncl,
                int *ihole, int *idimp, int *nppmx, int *nx, int *ny,
                int *mx, int *my, int *mx1, int *my1, int *npbmx,
                int *ntmax, int *irc);

void pporderf2l_(float *ppart, float *ppbuff, int *kpic, int *ncl,
                 int *ihole, int *idimp, int *nppmx, int *mx1,
                 int *my1, int *npbmx, int *ntmax, int *irc);

void mpois22_(float complex *q, float complex *fxy, int *isign,
              float complex *ffc, float *ax, float *ay, float *affp,
              float *we, int *nx, int *ny, int *nxvh, int *nyv,
              int *nxhd, int *nyhd);

void wfft2rinit_(int *mixup, float complex *sct, int *indx, int *indy,
                 int *nxhyd, int *nxyhd);

void fft2rmxx_(float complex *f, int *isign, int *mixup,
               float complex *sct, int *indx, int *indy, int *nyi,
               int *nyp, int *nxhd, int *nyd, int *nxhyd, int *nxyhd);

void fft2rmxy_(float complex *f, int *isign, int *mixup,
               float complex *sct, int *indx, int *indy, int *nxi,
               int *nxp, int *nxhd, int *nyd, int *nxhyd, int *nxyhd);

void fft2rm2x_(float complex *f, int *isign, int *mixup,
               float complex *sct, int *indx, int *indy, int *nyi,
               int *nyp, int *nxhd, int *nyd, int *nxhyd, int *nxyhd);

void fft2rm2y_(float complex *f, int *isign, int *mixup,
               float complex *sct, int *indx, int *indy, int *nxi,
               int *nxp, int *nxhd, int *nyd, int *nxhyd, int *nxyhd);

void wfft2rmx_(float complex *f, int *isign, int *mixup,
               float complex *sct, int *indx, int *indy, int *nxhd,
               int *nyd, int *nxhyd, int *nxyhd);

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
void cppmovin2l(float part[], float ppart[], int kpic[], int nppmx,
                int idimp, int nop, int mx, int my, int mx1, int mxy1,
                int *irc) {
   ppmovin2l_(part,ppart,kpic,&nppmx,&idimp,&nop,&mx,&my,&mx1,&mxy1,
              irc);
   return;
}

/*--------------------------------------------------------------------*/
void cppcheck2l(float ppart[], int kpic[], int idimp, int nppmx, int nx,
                int ny, int mx, int my, int mx1, int my1, 
                int *irc) {
   ppcheck2l_(ppart,kpic,&idimp,&nppmx,&nx,&ny,&mx,&my,&mx1,&my1,irc);
   return;
}

/*--------------------------------------------------------------------*/
void cgppush2l(float ppart[], float fxy[], int kpic[], float qbm,
               float dt, float *ek, int idimp, int nppmx, int nx,
               int ny, int mx, int my, int nxv, int nyv, int mx1,
               int mxy1, int ipbc) {
   gppush2l_(ppart,fxy,kpic,&qbm,&dt,ek,&idimp,&nppmx,&nx,&ny,&mx,&my,
             &nxv,&nyv,&mx1,&mxy1,&ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cgppushf2l(float ppart[], float fxy[], int kpic[], int ncl[],
                int ihole[], float qbm, float dt, float *ek, int idimp,
                int nppmx, int nx, int ny, int mx, int my, int nxv,
                int nyv, int mx1, int mxy1, int ntmax, int *irc) {
   gppushf2l_(ppart,fxy,kpic,ncl,ihole,&qbm,&dt,ek,&idimp,&nppmx,&nx,
              &ny,&mx,&my,&nxv,&nyv,&mx1,&mxy1,&ntmax,irc);
   return;
}

/*--------------------------------------------------------------------*/
void cgppost2l(float ppart[], float q[], int kpic[], float qm,
               int nppmx, int idimp, int mx, int my, int nxv, int nyv,
               int mx1, int mxy1) {
   gppost2l_(ppart,q,kpic,&qm,&nppmx,&idimp,&mx,&my,&nxv,&nyv,&mx1,
             &mxy1);
   return;
}

/*--------------------------------------------------------------------*/
void cpporder2l(float ppart[], float ppbuff[], int kpic[], int ncl[],
                int ihole[], int idimp, int nppmx, int nx, int ny,
                int mx, int my, int mx1, int my1, int npbmx, int ntmax,
                int *irc) {
   pporder2l_(ppart,ppbuff,kpic,ncl,ihole,&idimp,&nppmx,&nx,&ny,&mx,&my,
              &mx1,&my1,&npbmx,&ntmax,irc);
   return;
}

/*--------------------------------------------------------------------*/
void cpporderf2l(float ppart[], float ppbuff[], int kpic[], int ncl[],
                 int ihole[], int idimp, int nppmx, int mx1, int my1,
                 int npbmx, int ntmax, int *irc) {
   pporderf2l_(ppart,ppbuff,kpic,ncl,ihole,&idimp,&nppmx,&mx1,&my1,
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
void cmpois22(float complex q[], float complex fxy[], int isign,
              float complex ffc[], float ax, float ay, float affp,
              float *we, int nx, int ny, int nxvh, int nyv, int nxhd,
              int nyhd) {
   mpois22_(q,fxy,&isign,ffc,&ax,&ay,&affp,we,&nx,&ny,&nxvh,&nyv,&nxhd,
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
