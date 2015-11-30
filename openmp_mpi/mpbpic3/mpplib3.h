/* header file for mpplib3.c */

void cppinit2(int *idproc, int *nvp, int argc, char *argv[]);

void cppexit();

void cppabort();

void cpwtimera(int icntrl, float *time, double *dtime);

void cppsum(float f[], float g[], int nxp);

void cppdsum(double f[], double g[], int nxp);

void cppimax(int f[], int g[], int nxp);

void cppdmax(double f[], double g[], int nxp);

void cppncguard32l(float f[], float scs[], int nyzp[], int kstrt,
                   int nvpy, int nvpz, int nxv, int nypmx, int nzpmx,
                   int idds);

void cppnaguard32l(float f[], float scs[], float scr[], int nyzp[],
                   int kstrt, int nvpy, int nvpz, int nx, int nxv,
                   int nypmx, int nzpmx, int idds);

void cppnacguard32l(float f[], float scs[], float scr[], int nyzp[],
                    int ndim, int kstrt, int nvpy, int nvpz, int nx,
                    int nxv, int nypmx, int nzpmx, int idds);

void cpptpos3a(float complex f[], float complex g[], float complex s[], 
               float complex t[], int nx, int ny, int nz, int kxyp,
               int kyp, int kzp, int kstrt, int nvpy, int nxv, int nyv,
               int kxypd, int kypd, int kzpd);

void cpptpos3b(float complex g[], float complex h[], float complex s[],
               float complex t[], int nx, int ny, int nz, int kxyp,
               int kyzp, int kzp, int kstrt, int nvpy, int nvpz,
               int nyv, int nzv, int kxypd, int kyzpd, int kzpd);

void cppntpos3a(float complex f[], float complex g[], float complex s[],
                float complex t[], int nx, int ny, int nz, int kxyp,
                int kyp, int kzp, int kstrt, int nvpy, int ndim,
                int nxv, int nyv, int kxypd, int kypd, int kzpd);

void cppntpos3b(float complex g[], float complex h[], float complex s[],
                float complex t[], int nx, int ny, int nz, int kxyp,
                int kyzp, int kzp, int kstrt, int nvpy, int nvpz,
                int ndim, int nyv, int nzv, int kxypd, int kyzpd,
                int kzpd);

void cpppmove32(float sbufr[], float sbufl[], float rbufr[], 
                float rbufl[], int ncll[], int nclr[], int mcll[], 
                int mclr[], int mcls[], int kstrt, int nvpy, int nvpz,
                int idimp, int nbmax, int mx1, int myp1, int mzp1,
                int mxzyp1, int *irc);
