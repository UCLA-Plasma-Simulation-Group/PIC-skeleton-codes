/* Basic parallel PIC library for MPI communications */
/* Wrappers for calling the Fortran routines from a C main program */

#include <complex.h>

void ppinit2_(int *idproc, int *nvp, int *argc, char *argv[]);

void ppexit_();

void ppabort_();

void pwtimera_(int *icntrl, float *time, double *dtime);

void ppsum_(float *f, float *g, int *nxp);

void ppdsum_(double *f, double *g, int *nxp);

void ppimax_(int *f, int *g, int *nxp);

void ppdmax_(double *f, double *g, int *nxp);

void ppncguard2l_(float *f, int *nyp, int *kstrt, int *nvp, int *nxv,
                  int *nypmx);

void ppnaguard2l_(float *f, float *scr, int *nyp, int *nx, int *kstrt,
                  int *nvp, int *nxv, int *nypmx);

void ppnacguard2l_(float *f, float *scr, int *nyp, int *nx, int *ndim,
                   int *kstrt, int *nvp, int *nxv, int *nypmx);

void pptpose_(float complex *f, float complex *g, float complex *s,
              float complex *t, int *nx, int *ny, int *kxp, int *kyp,
              int *kstrt, int *nvp, int *nxv, int *nyv, int *kxpd,
              int *kypd);

void ppntpose_(float complex *f, float complex *g, float complex *s,
               float complex *t, int *nx, int *ny, int *kxp, int *kyp,
               int *kstrt, int *nvp, int *ndim, int *nxv, int *nyv,
               int *kxpd, int *kypd);

void ppmove2_(float *part, float *edges, int *npp, float *sbufr,
              float *sbufl, float *rbufr, float *rbufl, int *ihole,
              int *ny, int *kstrt, int *nvp, int *idimp, int *npmax,
              int *idps, int *nbmax, int *ntmax, int *info);

/* Interfaces to C */

/*--------------------------------------------------------------------*/
void cppinit2(int *idproc, int *nvp, int argc, char *argv[]) {
   ppinit2_(idproc,nvp,&argc,argv);
   return;
}

void cppexit() {
   ppexit_();
   return;
}

void cppabort() {
   ppabort_();
   return;
}

/*--------------------------------------------------------------------*/
void cpwtimera(int icntrl, float *time, double *dtime) {
   pwtimera_(&icntrl,time,dtime);
   return;
}

/*--------------------------------------------------------------------*/
void cppsum(float f[], float g[], int nxp) {
   ppsum_(f,g,&nxp);
   return;
}

/*--------------------------------------------------------------------*/
void cppdsum(double f[], double g[], int nxp) {
   ppdsum_(f,g,&nxp);
   return;
}

/*--------------------------------------------------------------------*/
void cppimax(int f[], int g[], int nxp) {
   ppimax_(f,g,&nxp);
   return;
}

/*--------------------------------------------------------------------*/
void cppdmax(double f[], double g[], int nxp) {
   ppdmax_(f,g,&nxp);
   return;
}

/*--------------------------------------------------------------------*/
void cppncguard2l(float f[], int nyp, int kstrt, int nvp, int nxv,
                  int nypmx) {
   ppncguard2l_(f,&nyp,&kstrt,&nvp,&nxv,&nypmx);
   return;
}

/*--------------------------------------------------------------------*/
void cppnaguard2l(float f[], float scr[], int nyp, int nx, int kstrt,
                  int nvp, int nxv, int nypmx) {
   ppnaguard2l_(f,scr,&nyp,&nx,&kstrt,&nvp,&nxv,&nypmx);
   return;
}

/*--------------------------------------------------------------------*/
void cppnacguard2l(float f[], float scr[], int nyp, int nx, int ndim,
                   int kstrt, int nvp, int nxv, int nypmx) {
   ppnacguard2l_(f,scr,&nyp,&nx,&ndim,&kstrt,&nvp,&nxv,&nypmx);
                   
   return;
}

/*--------------------------------------------------------------------*/
void cpptpose(float complex f[], float complex g[], float complex s[],
              float complex t[], int nx, int ny, int kxp, int kyp,
              int kstrt, int nvp, int nxv, int nyv, int kxpd, int kypd) {
   pptpose_(f,g,s,t,&nx,&ny,&kxp,&kyp,&kstrt,&nvp,&nxv,&nyv,&kxpd,
            &kypd);
   return;
}

/*--------------------------------------------------------------------*/
void cppntpose(float complex f[], float complex g[], float complex s[],
               float complex t[], int nx, int ny, int kxp, int kyp,
               int kstrt, int nvp, int ndim, int nxv, int nyv, int kxpd,
               int kypd) {
   ppntpose_(f,g,s,t,&nx,&ny,&kxp,&kyp,&kstrt,&nvp,&ndim,&nxv,&nyv,
             &kxpd,&kypd);
   return;
}

/*--------------------------------------------------------------------*/
void cppmove2(float part[], float edges[], int *npp, float sbufr[],
              float sbufl[], float rbufr[], float rbufl[], int ihole[],
              int ny, int kstrt, int nvp, int idimp, int npmax, int idps,
              int nbmax, int ntmax, int info[]) {
   ppmove2_(part,edges,npp,sbufr,sbufl,rbufr,rbufl,ihole,&ny,&kstrt,
            &nvp,&idimp,&npmax,&idps,&nbmax,&ntmax,info);
   return;
}
