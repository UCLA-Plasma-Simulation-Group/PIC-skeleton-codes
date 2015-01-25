/* Basic parallel PIC library for MPI communications */
/* Wrappers for calling the Fortran routines from a C main program */

#include <complex.h>

void ppinit2_(int *idproc, int *nvp, int *argc, char *argv[]);

void ppfndgrp_(int *locl, int *kstrt, int *nvp, int *idev, int *ndev);

void ppexit_();

void ppabort_();

void pwtimera_(int *icntrl, float *time, double *dtime);

void ppsum_(float *f, float *g, int *nxp);

void ppdsum_(double *f, double *g, int *nxp);

void ppimax_(int *f, int *g, int *nxp);

void pppcncguard2l_(float complex *scs, float complex *scr, int *kstrt,
                    int *nvp, int *nxvh);

void pppcnaguard2l_(float complex *scs, float complex *scr, int *kstrt,
                    int *nvp, int *nxvh);

void ppptpose_(float complex *sm, float complex *tm, int *nx, int *ny,
               int *kxp, int *kyp, int *kstrt, int *nvp);

void ppptposen_(float complex *sm, float complex *tm, int *nx, int *ny,
                int *kxp, int *kyp, int *kstrt, int *nvp, int *ndim);

void acsndrec_(float complex *stm, int *idproc, int *nsize, int *ntag,
               int *mode);

void pppmove2_(float *sbufr, float *sbufl, float *rbufr, float *rbufl,
               int *ncll, int *nclr, int *mcll, int *mclr, int *kstrt,
               int *nvp, int *idimp, int *nbmax, int *mx1);

/* Interfaces to C */

/*--------------------------------------------------------------------*/
void cppinit2(int *idproc, int *nvp, int argc, char *argv[]) {
   ppinit2_(idproc,nvp,&argc,argv);
   return;
}

/*--------------------------------------------------------------------*/
void cppfndgrp(int locl[], int kstrt, int nvp, int *idev, int *ndev) {
   ppfndgrp_(locl,&kstrt,&nvp,idev,ndev);
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
void cpppcncguard2l(float complex scs[], float complex scr[], int kstrt,
                   int nvp, int nxvh) {
   pppcncguard2l_(scs,scr,&kstrt,&nvp,&nxvh);
   return;
}

/*--------------------------------------------------------------------*/
void cpppcnaguard2l(float complex scs[], float complex scr[], int kstrt,
                    int nvp, int nxvh) {
   pppcnaguard2l_(scs,scr,&kstrt,&nvp,&nxvh);
   return;
}

/*--------------------------------------------------------------------*/
void cppptpose(float complex sm[], float complex tm[], int nx, int ny,
               int kxp, int kyp, int kstrt, int nvp) {
   ppptpose_(sm,tm,&nx,&ny,&kxp,&kyp,&kstrt,&nvp);
   return;
}

/*--------------------------------------------------------------------*/
void cppptposen(float complex sm[], float complex tm[], int nx, int ny,
                int kxp, int kyp, int kstrt, int nvp, int ndim) {
   ppptposen_(sm,tm,&nx,&ny,&kxp,&kyp,&kstrt,&nvp,&ndim);
   return;
}

/*--------------------------------------------------------------------*/
void cacsndrec(float complex stm[], int idproc, int nsize, int ntag,
               int mode) {
acsndrec_(stm,&idproc,&nsize,&ntag,&mode);
   return;
}

/*--------------------------------------------------------------------*/
void cpppmove2(float sbufr[], float sbufl[], float rbufr[], 
               float rbufl[], int ncll[], int nclr[], int mcll[],
               int mclr[], int kstrt, int nvp, int idimp, int nbmax,
               int mx1) {
   pppmove2_(sbufr,sbufl,rbufr,rbufl,ncll,nclr,mcll,mclr,&kstrt,&nvp,
             &idimp,&nbmax,&mx1);
   return;
}
