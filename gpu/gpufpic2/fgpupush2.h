/* C header file for gpupush2.cuf */

#include <complex.h>

void fgpuppush2l(float ppart[], float fxy[], int kpic[], float qbm,
                 float dt, float *ek, int idimp, int nppmx, int nx,
                 int ny, int mx, int my, int nxv, int nyv, int mx1,
                 int mxy1, int ipbc);

void fgpuppushf2l(float ppart[], float fxy[], int kpic[], int ncl[],
                  int ihole[], float qbm, float dt, float *ek,
                  int idimp, int nppmx, int nx, int ny, int mx, int my,
                  int nxv, int nyv, int mx1, int mxy1, int ntmax, 
                  int *irc);

void fgpu2ppost2l(float ppart[], float q[], int kpic[], float qm,
                  int nppmx, int idimp, int mx, int my, int nxv,
                  int nyv, int mx1, int mxy1);

void fgpucaguard2l(float complex qc[], float q[], int nx, int ny, 
                   int nxe, int nye, int nxvh, int nyv);

void fgpuccguard2l(float complex fxyc[], float fxy[], int nx, int ny,
                   int nxe, int nye, int nxvh, int nyv);

void fgpuppord2l(float ppart[], float ppbuff[], int kpic[], int ncl[],
                 int ihole[], int idimp, int nppmx, int nx, int ny,
                 int mx, int my, int mx1, int my1, int npbmx, int ntmax,
                 int *irc);

void fgpuppordf2l(float ppart[], float ppbuff[], int kpic[], int ncl[],
                  int ihole[], int idimp, int nppmx, int mx1, int my1,
                  int npbmx, int ntmax, int *irc);

void fgpupois22t(float complex qt[], float complex fxyt[],
                 float complex ffct[], float *we, int nx, int ny,
                 int nxvh, int nyv, int nxhd, int nyhd);

void fgpuwfft2rcs(float complex f[], float complex g[], int isign,
                  int mixup[], float complex sct[], int indx, int indy,
                  int nxhd, int nyd, int nxhyd, int nxyhd);

void fgpuwfft2rcsn(float complex fn[], float complex gn[], int isign,
                   int mixup[], float complex sct[], int indx, int indy,
                   int ndim, int nxhd, int nyd, int nxhyd, int nxyhd);

void fgpusum2(float a[], float *sa, int nx);
