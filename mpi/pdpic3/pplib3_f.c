/* Basic parallel PIC library for MPI communications */
/* Wrappers for calling the Fortran routines from a C main program */

#include <complex.h>

void ppinit2_(int *idproc, int *nvp);

void ppexit_();

void ppabort_();

void pwtimera_(int *icntrl, float *time, double *dtime);

void ppsum_(float *f, float *g, int *nxp);

void ppdsum_(double *f, double *g, int *nxp);

void ppimax_(int *f, int *g, int *nxp);

void ppdmax_(double *f, double *g, int *nxp);

void ppncguard32l_(float *f, float *scs, int *nyzp, int *kstrt,
                   int *nvpy, int *nvpz, int *nxv, int *nypmx, 
                   int *nzpmx, int *idds);

void ppnaguard32l_(float *f, float *scs, float *scr, int *nyzp,
                   int *kstrt, int *nvpy, int *nvpz, int *nx, int *nxv,
                   int *nypmx, int *nzpmx, int *idds);

void ppnacguard32l_(float *f, float *scs, float *scr, int *nyzp,
                    int *ndim, int *kstrt, int *nvpy, int *nvpz,
                    int *nx, int *nxv, int *nypmx, int *nzpmx,
                    int *idds);

void pptpos3a_(float complex *f, float complex *g, float complex *s, 
               float complex *t, int *nx, int *ny, int *nz, int *kxyp,
               int *kyp, int *kzp, int *kstrt, int *nvpy, int *nxv,
               int *nyv, int *kxypd, int *kypd, int *kzpd);

void pptpos3b_(float complex *g, float complex *h, float complex *s,
               float complex *t, int *nx, int *ny, int *nz, int *kxyp,
               int *kyzp, int *kzp, int *kstrt, int *nvpy, int *nvpz,
               int *nyv, int *nzv, int *kxypd, int *kyzpd, int *kzpd);

void ppntpos3a_(float complex *f, float complex *g, float complex *s,
                float complex *t, int *nx, int *ny, int *nz, int *kxyp,
                int *kyp, int *kzp, int *kstrt, int *nvpy, int *ndim,
                int *nxv, int *nyv, int *kxypd, int *kypd, int *kzpd);

void ppntpos3b_(float complex *g, float complex *h, float complex *s,
                float complex *t, int* nx, int *ny, int *nz, int *kxyp,
                int *kyzp, int *kzp, int *kstrt, int *nvpy, int *nvpz,
                int *ndim, int *nyv, int *nzv, int *kxypd, int *kyzpd,
                int *kzpd);

void ppmove32_(float *part, float *edges, int *npp, float *sbufr,
               float *sbufl, float *rbufr, float *rbufl, int *ihole, 
               int *ny, int *nz, int *kstrt, int *nvpy, int *nvpz,
               int *idimp, int *npmax, int *idps, int *nbmax,
               int *ntmax, int *info);

void ppmoveg32_(float *part, float *edges, int *npp, float *sbufr,
                float *sbufl, float *rbufr, float *rbufl, int *ihole, 
                int *ny, int *nz, int *kstrt, int *nvpy, int *nvpz,
                int *idimp, int *npmax, int *idps, int *nbmax,
                int *ntmax, int *info);

/* Interfaces to C */

/*--------------------------------------------------------------------*/
void cppinit2(int *idproc, int *nvp, int argc, char *argv[]) {
   ppinit2_(idproc,nvp);
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
void cppncguard32l(float f[], float scs[], int nyzp[], int kstrt,
                   int nvpy, int nvpz, int nxv, int nypmx, int nzpmx,
                   int idds) {
   ppncguard32l_(f,scs,nyzp,&kstrt,&nvpy,&nvpz,&nxv,&nypmx,&nzpmx,
                 &idds);
   return;
}

/*--------------------------------------------------------------------*/
void cppnaguard32l(float f[], float scs[], float scr[], int nyzp[],
                   int kstrt, int nvpy, int nvpz, int nx, int nxv,
                   int nypmx, int nzpmx, int idds) {
   ppnaguard32l_(f,scs,scr,nyzp,&kstrt,&nvpy,&nvpz,&nx,&nxv,&nypmx,
                 &nzpmx,&idds);
   return;
}

/*--------------------------------------------------------------------*/
void cppnacguard32l(float f[], float scs[], float scr[], int nyzp[],
                    int ndim, int kstrt, int nvpy, int nvpz, int nx,
                    int nxv, int nypmx, int nzpmx, int idds) {
   ppnacguard32l_(f,scs,scr,nyzp,&ndim,&kstrt,&nvpy,&nvpz,&nx,&nxv,
                  &nypmx,&nzpmx,&idds);
   return;
}

/*--------------------------------------------------------------------*/
void cpptpos3a(float complex f[], float complex g[], float complex s[], 
               float complex t[], int nx, int ny, int nz, int kxyp,
               int kyp, int kzp, int kstrt, int nvpy, int nxv, int nyv,
               int kxypd, int kypd, int kzpd) {
   pptpos3a_(f,g,s,t,&nx,&ny,&nz,&kxyp,&kyp,&kzp,&kstrt,&nvpy,&nxv,&nyv,
             &kxypd,&kypd,&kzpd);
   return;
}

/*--------------------------------------------------------------------*/
void cpptpos3b(float complex g[], float complex h[], float complex s[],
               float complex t[], int nx, int ny, int nz, int kxyp,
               int kyzp, int kzp, int kstrt, int nvpy, int nvpz,
               int nyv, int nzv, int kxypd, int kyzpd, int kzpd) {
   pptpos3b_(g,h,s,t,&nx,&ny,&nz,&kxyp,&kyzp,&kzp,&kstrt,&nvpy,&nvpz,
             &nyv,&nzv,&kxypd,&kyzpd,&kzpd);
   return;
}

/*--------------------------------------------------------------------*/
void cppntpos3a(float complex f[], float complex g[], float complex s[],
                float complex t[], int nx, int ny, int nz, int kxyp,
                int kyp, int kzp, int kstrt, int nvpy, int ndim,
                int nxv, int nyv, int kxypd, int kypd, int kzpd) {
   ppntpos3a_(f,g,s,t,&nx,&ny,&nz,&kxyp,&kyp,&kzp,&kstrt,&nvpy,&ndim,
              &nxv,&nyv,&kxypd,&kypd,&kzpd);
   return;
}

/*--------------------------------------------------------------------*/
void cppntpos3b(float complex g[], float complex h[], float complex s[],
                float complex t[], int nx, int ny, int nz, int kxyp,
                int kyzp, int kzp, int kstrt, int nvpy, int nvpz,
                int ndim, int nyv, int nzv, int kxypd, int kyzpd,
                int kzpd) {
   ppntpos3b_(g,h,s,t,&nx,&ny,&nz,&kxyp,&kyzp,&kzp,&kstrt,&nvpy,&nvpz,
                &ndim,&nyv,&nzv,&kxypd,&kyzpd,&kzpd);
   return;
}

/*--------------------------------------------------------------------*/
void cppmove32(float part[], float edges[], int *npp, float sbufr[],
               float sbufl[], float rbufr[], float rbufl[], int ihole[], 
               int ny, int nz, int kstrt, int nvpy, int nvpz, int idimp,
               int npmax, int idps, int nbmax, int ntmax, int info[]) {
   ppmove32_(part,edges,npp,sbufr,sbufl,rbufr,rbufl,ihole,&ny,&nz,
             &kstrt,&nvpy,&nvpz,&idimp,&npmax,&idps,&nbmax,&ntmax,info);
   return;
}

/*--------------------------------------------------------------------*/
void cppmoveg32(float part[], float edges[], int *npp, float sbufr[],
               float sbufl[], float rbufr[], float rbufl[], int ihole[], 
               int ny, int nz, int kstrt, int nvpy, int nvpz, int idimp,
               int npmax, int idps, int nbmax, int ntmax, int info[]) {
   ppmoveg32_(part,edges,npp,sbufr,sbufl,rbufr,rbufl,ihole,&ny,&nz,
              &kstrt,&nvpy,&nvpz,&idimp,&npmax,&idps,&nbmax,&ntmax,
              info);
   return;
}
