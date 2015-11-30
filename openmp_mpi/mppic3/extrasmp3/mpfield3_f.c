/* C Library for Skeleton 3D Electrostatic MPI PIC Code field    */
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

void mppsmooth32_(float complex *q, float complex *qs,
                  float complex *ffc, int *nx, int *ny, int *nz,
                  int *kstrt, int *nvpy, int *nvpz, int *nzv, int *kxyp, 
                  int *kyzp, int *nzhd);

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
void cmppsmooth32(float complex q[], float complex qs[],
                  float complex ffc[], int nx, int ny, int nz,
                  int kstrt, int nvpy, int nvpz, int nzv, int kxyp,
                  int kyzp, int nzhd) {
   mppsmooth32_(q,qs,ffc,&nx,&ny,&nz,&kstrt,&nvpy,&nvpz,&nzv,&kxyp,
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
