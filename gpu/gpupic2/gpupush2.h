/* C header file for gpupush2.cu */

#include <complex.h>

void cgpuppush2l(float ppart[], float fxy[], int kpic[], float qbm,
                 float dt, float *ek, int idimp, int nppmx, int nx,
                 int ny, int mx, int my, int nxv, int nyv, int mx1,
                 int mxy1, int ipbc);

void cgpuppushf2l(float ppart[], float fxy[], int kpic[], int ncl[],
                  int ihole[], float qbm, float dt, float *ek,
                  int idimp, int nppmx, int nx, int ny, int mx, int my,
                  int nxv, int nyv, int mx1, int mxy1, int ntmax, 
                  int *irc);

void cgpu2ppost2l(float ppart[], float q[], int kpic[], float qm,
                  int nppmx, int idimp, int mx, int my, int nxv,
                  int nyv, int mx1, int mxy1);

void cgpucaguard2l(float complex qc[], float q[], int nx, int ny, 
                   int nxe, int nye, int nxvh, int nyv);

void cgpuccguard2l(float complex fxyc[], float fxy[], int nx, int ny,
                   int nxe, int nye, int nxvh, int nyv);

void cgpuppord2l(float ppart[], float ppbuff[], int kpic[], int ncl[],
                 int ihole[], int idimp, int nppmx, int nx, int ny,
                 int mx, int my, int mx1, int my1, int npbmx, int ntmax,
                 int *irc);

void cgpuppordf2l(float ppart[], float ppbuff[], int kpic[], int ncl[],
                  int ihole[], int idimp, int nppmx, int mx1, int my1,
                  int npbmx, int ntmax, int *irc);

void cgpupois22t(float complex qt[], float complex fxyt[],
                 float complex ffct[], float *we, int nx, int ny,
                 int nxvh, int nyv, int nxhd, int nyhd);

void cgpuwfft2rcs(float complex f[], float complex g[], int isign,
                  int mixup[], float complex sct[], int indx, int indy,
                  int nxhd, int nyd, int nxhyd, int nxyhd);

void cgpuwfft2rcsn(float complex fn[], float complex gn[], int isign,
                   int mixup[], float complex sct[], int indx, int indy,
                   int ndim, int nxhd, int nyd, int nxhyd, int nxyhd);

void cgpusum2(float a[], float *sa, int nx);
