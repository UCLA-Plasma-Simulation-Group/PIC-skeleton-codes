/* header file for mpplib2.c */

void cppinit2(int *idproc, int *nvp, int argc, char *argv[]);

void cppexit();

void cppabort();

void cpwtimera(int icntrl, float *time, double *dtime);

void cppsum(float f[], float g[], int nxp);

void cppdsum(double f[], double g[], int nxp);

void cppimax(int f[], int g[], int nxp);

void cppdmax(double f[], double g[], int nxp);

void cppncguard2l(float f[], int nyp, int kstrt, int nvp, int nxv,
                  int nypmx);

void cppnaguard2l(float f[], float scr[], int nyp, int nx, int kstrt,
                  int nvp, int nxv, int nypmx);

void cppnacguard2l(float f[], float scr[], int nyp, int nx, int ndim,
                   int kstrt, int nvp, int nxv, int nypmx);

void cpptpose(float complex f[], float complex g[], float complex s[],
              float complex t[], int nx, int ny, int kxp, int kyp,
              int kstrt, int nvp, int nxv, int nyv, int kxpd, int kypd);

void cppntpose(float complex f[], float complex g[], float complex s[],
               float complex t[], int nx, int ny, int kxp, int kyp,
               int kstrt, int nvp, int ndim, int nxv, int nyv, int kxpd,
               int kypd);

void cpppmove2(float sbufr[], float sbufl[], float rbufr[], 
               float rbufl[], int ncll[], int nclr[], int mcll[],
               int mclr[], int kstrt, int nvp, int idimp, int nbmax,
               int mx1);
