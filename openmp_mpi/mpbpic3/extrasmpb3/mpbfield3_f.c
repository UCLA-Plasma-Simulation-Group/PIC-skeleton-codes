/* C Library for Skeleton 3D Electromagnetic MPI PIC Code field    */
/* diagnostics */
/* Wrappers for calling the Fortran routines from a C main program */

#include <complex.h>

void mppotp32_(float complex *q, float complex *pot, float complex *ffc,
               float *we, int *nx, int *ny, int *nz, int *kstrt,
               int *nvpy, int *nvpz, int *nzv, int *kxyp, int *kyzp,
               int *nzhd);

void mppdivf32_(float complex *f, float complex *df, int *nx, int *ny,
                int *nz, int *kstrt, int *nvpy, int *nvpz, int *nzv,
                int *kxyp, int *kyzp);

void mppgradf32_(float complex *df, float complex *f, int *nx, int *ny,
                 int *nz, int *kstrt, int *nvpy, int *nvpz, int *nzv,
                 int *kxyp, int *kyzp);

void mppcurlf32_(float complex *f, float complex *g, int *nx, int *ny,
                 int *nz, int *kstrt, int *nvpy, int *nvpz, int *nzv,
                 int *kxyp, int *kyzp);

void mppavpot332_(float complex *bxyz, float complex *axyz, int *nx,
                  int *ny, int *nz, int *kstrt, int *nvpy, int *nvpz,
                  int *nzv, int *kxyp, int *kyzp);

void mppavrpot332_(float complex *axyz, float complex *bxyz,
                   float complex *ffc, float *affp, float *ci, int *nx,
                   int *ny, int *nz, int *kstrt, int *nvpy, int *nvpz,
                   int *nzv, int *kxyp, int *kyzp, int *nzhd);

void mppsmooth32_(float complex *q, float complex *qs,
                  float complex *ffc, int *nx, int *ny, int *nz,
                  int *kstrt, int *nvpy, int *nvpz, int *nzv, int *kxyp, 
                  int *kyzp, int *nzhd);

void mppsmooth332_(float complex *cu, float complex *cus,
                   float complex *ffc, int *nx, int *ny, int *nz,
                   int *kstrt, int *nvpy, int *nvpz, int *nzv,
                   int *kxyp, int *kyzp, int *nzhd);

void pprdmodes32_(float complex *pot, float complex *pott, int *nx,
                  int *ny, int *nz, int *modesx, int *modesy, 
                  int *modesz, int *kstrt, int *nvpy, int *nvpz,
                  int *nzv, int *kxyp, int *kyzp, int *modesxpd,
                  int *modesypd, int *modeszd);

void ppwrmodes32_(float complex *pot, float complex *pott, int *nx,
                  int *ny, int *nz, int *modesx, int *modesy,
                  int *modesz, int *kstrt, int *nvpy, int *nvpz,
                  int *nzv, int *kxyp, int *kyzp, int *modesxpd,
                  int *modesypd, int *modeszd);

void pprdvmodes32_(float complex *vpot, float complex *vpott, int *nx,
                   int *ny, int *nz, int *modesx, int *modesy,
                   int *modesz, int *ndim, int *kstrt, int *nvpy,
                   int *nvpz, int *nzv, int *kxyp, int *kyzp,
                   int *modesxpd, int *modesypd, int *modeszd);

void ppwrvmodes32_(float complex *vpot, float complex *vpott, int *nx,
                   int *ny, int *nz, int *modesx, int *modesy,
                   int *modesz, int *ndim, int *kstrt, int *nvpy,
                   int *nvpz, int *nzv, int *kxyp, int *kyzp,
                   int *modesxpd, int *modesypd, int *modeszd);

/* Interfaces to C */

/*--------------------------------------------------------------------*/
void cmppotp32(float complex q[], float complex pot[],
               float complex ffc[], float *we, int nx, int ny, int nz,
               int kstrt, int nvpy, int nvpz, int nzv, int kxyp,
               int kyzp, int nzhd) {
   mppotp32_(q,pot,ffc,we,&nx,&ny,&nz,&kstrt,&nvpy,&nvpz,&nzv,&kxyp,
             &kyzp,&nzhd);
   return;
}

/*--------------------------------------------------------------------*/
void cmppdivf32(float complex f[], float complex df[], int nx, int ny,
                int nz, int kstrt, int nvpy, int nvpz, int nzv,
                int kxyp, int kyzp) {
   mppdivf32_(f,df,&nx,&ny,&nz,&kstrt,&nvpy,&nvpz,&nzv,&kxyp,&kyzp);
   return;
}

/*--------------------------------------------------------------------*/
void cmppgradf32(float complex df[], float complex f[], int nx, int ny,
                 int nz, int kstrt, int nvpy, int nvpz, int nzv,
                 int kxyp, int kyzp) {
   mppgradf32_(df,f,&nx,&ny,&nz,&kstrt,&nvpy,&nvpz,&nzv,&kxyp,&kyzp);
   return;
}

/*--------------------------------------------------------------------*/
void cmppcurlf32(float complex f[], float complex g[], int nx, int ny,
                 int nz, int kstrt, int nvpy, int nvpz, int nzv,
                 int kxyp, int kyzp) {
   mppcurlf32_(f,g,&nx,&ny,&nz,&kstrt,&nvpy,&nvpz,&nzv,&kxyp,&kyzp);
   return;
}

/*--------------------------------------------------------------------*/
void cmppavpot332(float complex bxyz[], float complex axyz[], int nx,
                  int ny, int nz, int kstrt, int nvpy, int nvpz,
                  int nzv, int kxyp, int kyzp) {
   mppavpot332_(bxyz,axyz,&nx,&ny,&nz,&kstrt,&nvpy,&nvpz,&nzv,&kxyp,
               &kyzp);
   return;
}

/*--------------------------------------------------------------------*/
void cmppavrpot332(float complex axyz[], float complex bxyz[],
                   float complex ffc[], float affp, float ci, int nx,
                   int ny, int nz, int kstrt, int nvpy, int nvpz,
                   int nzv, int kxyp, int kyzp, int nzhd) {
   mppavrpot332_(axyz,bxyz,ffc,&affp,&ci,&nx,&ny,&nz,&kstrt,&nvpy,&nvpz,
                &nzv,&kxyp,&kyzp,&nzhd);
   return;
}

/*--------------------------------------------------------------------*/
void cmppsmooth32(float complex q[], float complex qs[],
                  float complex ffc[], int nx, int ny, int nz,
                  int kstrt, int nvpy, int nvpz, int nzv, int kxyp,
                  int kyzp, int nzhd) {
   mppsmooth32_(q,qs,ffc,&nx,&ny,&nz,&kstrt,&nvpy,&nvpz,&nzv,&kxyp,
                &kyzp,&nzhd);
   return;
}

/*--------------------------------------------------------------------*/
void cmppsmooth332(float complex cu[], float complex cus[],
                   float complex ffc[], int nx, int ny, int nz,
                   int kstrt, int nvpy, int nvpz, int nzv, int kxyp,
                   int kyzp, int nzhd) {
   mppsmooth332_(cu,cus,ffc,&nx,&ny,&nz,&kstrt,&nvpy,&nvpz,&nzv,&kxyp,
                 &kyzp,&nzhd);
   return;
}

/*--------------------------------------------------------------------*/
void cpprdmodes32(float complex pot[], float complex pott[], int nx,
                  int ny, int nz, int modesx, int modesy, int modesz,    
                  int kstrt, int nvpy, int nvpz, int nzv, int kxyp,
                  int kyzp, int modesxpd, int modesypd, int modeszd) {
   pprdmodes32_(pot,pott,&nx,&ny,&nz,&modesx,&modesy,&modesz,&kstrt,
                &nvpy,&nvpz,&nzv,&kxyp,&kyzp,&modesxpd,&modesypd,
                &modeszd);
   return;
}

/*--------------------------------------------------------------------*/
void cppwrmodes32(float complex pot[], float complex pott[], int nx,
                  int ny, int nz, int modesx, int modesy, int modesz,
                  int kstrt, int nvpy, int nvpz, int nzv, int kxyp,
                  int kyzp, int modesxpd, int modesypd, int modeszd) {
   ppwrmodes32_(pot,pott,&nx,&ny,&nz,&modesx,&modesy,&modesz,&kstrt,
                &nvpy,&nvpz,&nzv,&kxyp,&kyzp,&modesxpd,&modesypd,
                &modeszd);
   return;
}

/*--------------------------------------------------------------------*/
void cpprdvmodes32(float complex vpot[], float complex vpott[], int nx,
                   int ny, int nz, int modesx, int modesy, int modesz, 
                   int ndim, int kstrt, int nvpy, int nvpz, int nzv,
                   int kxyp, int kyzp, int modesxpd, int modesypd, 
                   int modeszd) {
   pprdvmodes32_(vpot,vpott,&nx,&ny,&nz,&modesx,&modesy,&modesz,&ndim,
                 &kstrt,&nvpy,&nvpz,&nzv,&kxyp,&kyzp,&modesxpd,
                 &modesypd,&modeszd);
   return;
}

/*--------------------------------------------------------------------*/
void cppwrvmodes32(float complex vpot[], float complex vpott[], int nx,
                   int ny, int nz, int modesx, int modesy, int modesz, 
                   int ndim, int kstrt, int nvpy, int nvpz, int nzv,
                   int kxyp, int kyzp, int modesxpd, int modesypd,
                   int modeszd) {
   ppwrvmodes32_(vpot,vpott,&nx,&ny,&nz,&modesx,&modesy,&modesz,&ndim,
                 &kstrt,&nvpy,&nvpz,&nzv,&kxyp,&kyzp,&modesxpd,
                 &modesypd,&modeszd);
   return;
}
