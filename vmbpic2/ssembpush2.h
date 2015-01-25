/* header file for ssembpush2.c */

void csse2gbppush23lt(float ppart[], float fxy[], float bxy[],
                      int kpic[], float qbm, float dt, float dtc,
                      float *ek, int idimp, int nppmx, int nx, int ny,
                      int mx, int my, int nxv, int nyv, int mx1,
                      int mxy1, int ipbc);

void csse2gbppushf23lt(float ppart[], float fxy[], float bxy[],
                       int kpic[], int ncl[], int ihole[], float qbm,
                       float dt, float dtc, float *ek, int idimp,
                       int nppmx, int nx, int ny, int mx, int my,
                       int nxv, int nyv, int mx1, int mxy1, int ntmax,
                       int *irc);

void csse2grbppush23lt(float ppart[], float fxy[], float bxy[],
                       int kpic[], float qbm, float dt, float dtc, 
                       float ci, float *ek, int idimp, int nppmx,
                       int nx, int ny, int mx, int my, int nxv, int nyv,
                       int mx1, int mxy1, int ipbc);

void csse2grbppushf23lt(float ppart[], float fxy[], float bxy[],
                        int kpic[], int ncl[], int ihole[], float qbm,
                        float dt, float dtc, float ci, float *ek,
                        int idimp, int nppmx, int nx, int ny, int mx,
                        int my, int nxv, int nyv, int mx1, int mxy1,
                        int ntmax, int *irc);

void csse2gppost2lt(float ppart[], float q[], int kpic[], float qm,
                    int nppmx, int idimp, int mx, int my, int nxv,
                    int nyv, int mx1, int mxy1);

void csse2gjppost2lt(float ppart[], float cu[], int kpic[], float qm,
                    float dt, int nppmx, int idimp, int nx, int ny,
                    int mx, int my, int nxv, int nyv, int mx1, int mxy1,
                    int ipbc);

void csse2gjppostf2lt(float ppart[], float cu[], int kpic[], int ncl[],
                      int ihole[], float qm, float dt, int nppmx,
                      int idimp, int nx, int ny, int mx, int my,
                      int nxv, int nyv, int mx1, int mxy1, int ntmax,
                      int *irc);

void csse2grjppost2lt(float ppart[], float cu[], int kpic[], float qm,
                      float dt, float ci, int nppmx, int idimp, int nx,
                      int ny, int mx, int my, int nxv, int nyv, int mx1,
                      int mxy1, int ipbc);

void csse2grjppostf2lt(float ppart[], float cu[], int kpic[], int ncl[],
                   int ihole[], float qm, float dt, float ci, int nppmx,
                   int idimp, int nx, int ny, int mx, int my, int nxv,
                   int nyv, int mx1, int mxy1, int ntmax, int *irc);

void csse2pporder2lt(float ppart[], float ppbuff[], int kpic[],
                     int ncl[], int ihole[], int idimp, int nppmx,
                     int nx, int ny, int mx, int my, int mx1, int my1,
                     int npbmx, int ntmax, int *irc);

void csse2pporderf2lt(float ppart[], float ppbuff[], int kpic[],
                      int ncl[], int ihole[], int idimp, int nppmx,
                      int mx1, int my1, int npbmx, int ntmax,
                      int *irc);

void csse2bguard2l(float bxy[], int nx, int ny, int nxe, int nye);

void csse2acguard2l(float cu[], int nx, int ny, int nxe, int nye);

void csse2aguard2l(float q[], int nx, int ny, int nxe, int nye);

void csse2mpois23(float complex q[], float complex fxy[], int isign,
                  float complex ffc[], float ax, float ay, float affp,
                  float *we, int nx, int ny, int nxvh, int nyv,
                  int nxhd, int nyhd);

void csse2mcuperp2(float complex cu[], int nx, int ny, int nxvh,
                   int nyv);

void csse2mibpois23(float complex cu[], float complex bxy[],
                    float complex ffc[], float ci, float *wm, int nx,
                    int ny, int nxvh, int nyv, int nxhd, int nyhd);

void csse2mmaxwel2(float complex exy[], float complex bxy[],
                   float complex cu[], float complex ffc[], float ci,
                   float dt, float *wf, float *wm, int nx, int ny,
                   int nxvh, int nyv, int nxhd, int nyhd);

void csse2memfield2(float complex fxy[], float complex exy[],
                    float complex ffc[], int isign, int nx, int ny,
                    int nxvh, int nyv, int nxhd, int nyhd);

void csse2fft2rmxx(float complex f[], int isign, int mixup[],
                   float complex sct[], int indx, int indy, int nyi,
                   int nyp, int nxhd, int nyd, int nxhyd, int nxyhd);

void csse2fft2rmxy(float complex f[], int isign, int mixup[],
                   float complex sct[], int indx, int indy, int nxi,
                   int nxp, int nxhd, int nyd, int nxhyd, int nxyhd);

void csse2fft2rm3x(float complex f[], int isign, int mixup[],
                   float complex sct[], int indx, int indy, int nyi,
                   int nyp, int nxhd, int nyd, int nxhyd, int nxyhd);

void csse2fft2rm3y(float complex f[], int isign, int mixup[],
                   float complex sct[], int indx, int indy, int nxi,
                   int nxp, int nxhd, int nyd, int nxhyd, int nxyhd);

void csse2wfft2rmx(float complex f[], int isign, int mixup[],
                   float complex sct[], int indx, int indy, int nxhd,
                   int nyd, int nxhyd, int nxyhd);

void csse2wfft2rm3(float complex f[], int isign, int mixup[],
                   float complex sct[], int indx, int indy, int nxhd,
                   int nyd, int nxhyd, int nxyhd);
