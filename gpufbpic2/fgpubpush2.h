/* C header file for gpubpush2.cuf */

#include <complex.h>

void fgpubppush23l(float ppart[], float fxy[], float bxy[], int kpic[],
                   float qbm, float dt, float dtc, float *ek, int idimp,
                   int nppmx, int nx, int ny, int mx, int my, int nxv,
                   int nyv, int mx1, int mxy1, int ipbc);

void fgpubppushf23l(float ppart[], float fxy[], float bxy[], int kpic[],
                    int ncl[], int ihole[], float qbm, float dt,
                    float dtc, float *ek, int idimp, int nppmx, int nx,
                    int ny, int mx, int my, int nxv, int nyv, int mx1,
                    int mxy1, int ntmax, int *irc);

void fgpurbppush23l(float ppart[], float fxy[], float bxy[], int kpic[],
                    float qbm, float dt, float dtc, float ci, float *ek,
                    int idimp, int nppmx, int nx, int ny, int mx,
                    int my, int nxv, int nyv, int mx1, int mxy1,
                    int ipbc);

void fgpurbppushf23l(float ppart[], float fxy[], float bxy[],
                     int kpic[], int ncl[], int ihole[], float qbm,
                     float dt, float dtc, float ci, float *ek,
                     int idimp, int nppmx, int nx, int ny, int mx,
                     int my, int nxv, int nyv, int mx1, int mxy1,
                     int ntmax, int *irc);

void fgpu2ppost2l(float ppart[], float q[], int kpic[], float qm,
                  int nppmx, int idimp, int mx, int my, int nxv,
                  int nyv, int mx1, int mxy1);

void fgpu2jppost2l(float ppart[], float cu[], int kpic[], float qm, 
                   float dt, int nppmx, int idimp, int nx, int ny,
                   int mx, int my, int nxv, int nyv, int mx1, int mxy1,
                   int ipbc);

void fgpu2jppostf2l(float ppart[], float cu[], int kpic[], int ncl[],
                    int ihole[], float qm, float dt, int nppmx,
                    int idimp, int nx, int ny, int mx, int my, int nxv,
                    int nyv, int mx1, int mxy1, int ntmax, int *irc);

void fgpu2rjppost2l(float ppart[], float cu[], int kpic[], float qm,
                    float dt, float ci, int nppmx, int idimp, int nx,
                    int ny, int mx, int my, int nxv, int nyv, int mx1,
                    int mxy1, int ipbc);

void fgpu2rjppostf2l(float ppart[], float cu[], int kpic[], int ncl[],
                     int ihole[], float qm, float dt, float ci,
                     int nppmx, int idimp, int nx, int ny, int mx,
                     int my, int nxv, int nyv, int mx1, int mxy1,
                     int ntmax, int *irc);

void fgpucaguard2l(float complex qc[], float q[], int nx, int ny,
                   int nxe, int nye, int nxvh, int nyv);

void fgpucacguard2l(float complex cuc[], float cu[], int nx, int ny,
                    int nxe, int nye, int nxvh, int nyv);

void fgpucbguard2l(float complex bxyc[], float bxy[], int nx, int ny,
                   int nxe, int nye, int nxvh, int nyv);

void fgpuppord2l(float ppart[], float ppbuff[], int kpic[], int ncl[],
                 int ihole[], int idimp, int nppmx, int nx, int ny,
                 int mx, int my, int mx1, int my1, int npbmx, int ntmax,
                 int *irc);

void fgpuppordf2l(float ppart[], float ppbuff[], int kpic[], int ncl[],
                  int ihole[], int idimp, int nppmx, int mx1, int my1,
                  int npbmx, int ntmax, int *irc);

void fgpupois23t(float complex qt[], float complex fxyt[],
                 float complex ffct[], float *we, int nx, int ny,
                 int nxvh, int nyv, int nxhd, int nyhd);

void fgpucuperp2t(float complex cut[], int nx, int ny, int nxvh,
                  int nyv);

void fgpuibpois23t(float complex cut[], float complex bxyt[],
                   float complex ffct[], float ci, float *wm, int nx,
                   int ny, int nxvh, int nyv, int nxhd, int nyhd);

void fgpumaxwel2t(float complex exyt[], float complex bxyt[],
                  float complex cut[], float complex ffct[], float ci,
                  float dt, float *wf, float *wm, int nx, int ny,
                  int nxvh, int nyv, int nxhd, int nyhd);

void fgpuemfield2t(float complex fxyt[], float complex exyt[],
                   float complex ffct[], int isign, int nx, int ny,
                   int nxvh, int nyv, int nxhd, int nyhd);

void fgpuwfft2rcs(float complex f[], float complex g[], int isign,
                  int mixup[], float complex sct[], int indx, int indy,
                  int nxhd, int nyd, int nxhyd, int nxyhd);

void fgpuwfft2rcsn(float complex fn[], float complex gn[], int isign,
                   int mixup[], float complex sct[], int indx, int indy,
                   int ndim, int nxhd, int nyd, int nxhyd, int nxyhd);

void fgpusum2(float a[], float *sa, int nx);
