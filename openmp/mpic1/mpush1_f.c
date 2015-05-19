/* C Library for Skeleton 1D Electrostatic OpenMP PIC Code */
/* Wrappers for calling the Fortran routines from a C main program */

#include <complex.h>

double ranorm_();

void distr1_(float *part, float *vtx, float *vdx, int *npx, int *idimp,
             int *nop, int *nx, int *ipbc);

void dblkp1l_(float *part, int *kpic, int *nppmx, int *idimp, int *nop,
              int *mx, int *mx1, int *irc);

void ppmovin1l_(float *part, float *ppart, int *kpic, int *nppmx,
                int *idimp, int *nop, int *mx, int *mx1, int *irc);

void ppcheck1l_(float *ppart, int *kpic, int *idimp, int *nppmx,
                int *nx, int *mx, int *mx1, int *irc);

void gppush1l_(float *ppart, float *fx, int *kpic, float *qbm,
               float *dt, float *ek, int *idimp, int *nppmx, int *nx,
               int *mx, int *nxv, int *mx1, int *ipbc);

void gppushf1l_(float *ppart, float *fx, int *kpic, int *ncl,
                int *ihole, float *qbm, float *dt, float *ek,
                int *idimp, int *nppmx, int *nx, int *mx, int *nxv,
                int *mx1, int *ntmax, int *irc);

void gppost1l_(float *ppart, float *q, int *kpic, float *qm, int *nppmx,
               int *idimp, int *mx, int *nxv, int *mx1);

void pporder1l_(float *ppart, float *ppbuff, int *kpic, int *ncl,
                int *ihole, int *idimp, int *nppmx, int *nx, int *mx,
                int *mx1, int *npbmx, int *ntmax, int *irc);

void pporderf1l_(float *ppart, float *ppbuff, int *kpic, int *ncl,
                 int *ihole, int *idimp, int *nppmx, int *mx1,
                 int *npbmx, int *ntmax, int *irc);

void cguard1l_(float *fx, int *nx, int *nxe);

void aguard1l_(float *q, int *nx, int *nxe);

void pois1_(float complex *q, float complex *fx, int *isign,
            float complex *ffc, float *ax, float *affp, float *we,
            int *nx);

void wfft1rinit_(int *mixup, float complex *sct, int *indx, int *nxhd);

void fft1rxx_(float complex *f, float complex *t, int *isign,
              int *mixup, float complex *sct, int *indx, int *nxd, 
              int *nxhd);

/* Interfaces to C */

double ranorm() {
  return ranorm_();
}


/*--------------------------------------------------------------------*/
void cdistr1(float part[], float vtx, float vdx, int npx, int idimp,
             int nop, int nx, int ipbc) {
   distr1_(part,&vtx,&vdx,&npx,&idimp,&nop,&nx,&ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cdblkp1l(float part[], int kpic[], int *nppmx, int idimp, int nop,
              int mx, int mx1, int *irc) {
   dblkp1l_(part,kpic,nppmx,&idimp,&nop,&mx,&mx1,irc);
   return;
}

/*--------------------------------------------------------------------*/
void cppmovin1l(float part[], float ppart[], int kpic[], int nppmx,
                int idimp, int nop, int mx, int mx1, int *irc) {
   ppmovin1l_(part,ppart,kpic,&nppmx,&idimp,&nop,&mx,&mx1,irc);
   return;
}

/*--------------------------------------------------------------------*/
void cppcheck1l(float ppart[], int kpic[], int idimp, int nppmx, int nx,
                int mx, int mx1, int *irc) {
   ppcheck1l_(ppart,kpic,&idimp,&nppmx,&nx,&mx,&mx1,irc);
   return;
}

/*--------------------------------------------------------------------*/
void cgppush1l(float ppart[], float fx[], int kpic[], float qbm,
               float dt, float *ek, int idimp, int nppmx, int nx,
               int mx, int nxv, int mx1, int ipbc) {
   gppush1l_(ppart,fx,kpic,&qbm,&dt,ek,&idimp,&nppmx,&nx,&mx,&nxv,&mx1,
             &ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cgppushf1l(float ppart[], float fx[], int kpic[], int ncl[],
                int ihole[], float qbm, float dt, float *ek, int idimp,
                int nppmx, int nx, int mx, int nxv, int mx1, int ntmax,
                int *irc) {
   gppushf1l_(ppart,fx,kpic,ncl,ihole,&qbm,&dt,ek,&idimp,&nppmx,&nx,&mx,
              &nxv,&mx1,&ntmax,irc);
   return;
}

/*--------------------------------------------------------------------*/
void cgppost1l(float ppart[], float q[], int kpic[], float qm,
               int nppmx, int idimp, int mx, int nxv, int mx1) {
   gppost1l_(ppart,q,kpic,&qm,&nppmx,&idimp,&mx,&nxv,&mx1);
   return;
}

/*--------------------------------------------------------------------*/
void cpporder1l(float ppart[], float ppbuff[], int kpic[], int ncl[],
                int ihole[], int idimp, int nppmx, int nx, int mx,
                int mx1, int npbmx, int ntmax, int *irc) {
   pporder1l_(ppart,ppbuff,kpic,ncl,ihole,&idimp,&nppmx,&nx,&mx,&mx1,
              &npbmx,&ntmax,irc);
   return;
}

/*--------------------------------------------------------------------*/
void cpporderf1l(float ppart[], float ppbuff[], int kpic[], int ncl[],
                 int ihole[], int idimp, int nppmx, int mx1, int npbmx,
                 int ntmax, int *irc) {
   pporderf1l_(ppart,ppbuff,kpic,ncl,ihole,&idimp,&nppmx,&mx1,&npbmx,
               &ntmax,irc);
   return;
}

/*--------------------------------------------------------------------*/
void ccguard1l(float fx[], int nx, int nxe) {
   cguard1l_(fx,&nx,&nxe);
   return;
}

/*--------------------------------------------------------------------*/
void caguard1l(float q[], int nx, int nxe) {
   aguard1l_(q,&nx,&nxe);
   return;
}

/*--------------------------------------------------------------------*/
void cpois1(float complex q[], float complex fx[], int isign,
            float complex ffc[], float ax, float affp, float *we,
            int nx) {
   pois1_(q,fx,&isign,ffc,&ax,&affp,we,&nx);
   return;
}

/*--------------------------------------------------------------------*/
void cwfft1rinit(int mixup[], float complex sct[], int indx, int nxhd) {
   wfft1rinit_(mixup,sct,&indx,&nxhd);
   return;
}

/*--------------------------------------------------------------------*/
void cfft1rxx(float complex f[], float complex t[], int isign,
              int mixup[], float complex sct[], int indx, int nxd, 
              int nxhd) {
   fft1rxx_(f,t,&isign,mixup,sct,&indx,&nxd,&nxhd);
   return;
}
