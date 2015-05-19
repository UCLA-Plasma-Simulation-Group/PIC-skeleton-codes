/* header file for ssempush2.c */

void csse2gppush2lt(float ppart[], float fxy[], int kpic[], float qbm,
                    float dt, float *ek, int idimp, int nppmx, int nx,
                    int ny, int mx, int my, int nxv, int nyv, int mx1,
                    int mxy1, int ipbc);

void csse2gppushf2lt(float ppart[], float fxy[], int kpic[], int ncl[],
                     int ihole[], float qbm, float dt, float *ek, 
                     int idimp, int nppmx, int nx, int ny, int mx,
                     int my, int nxv, int nyv, int mx1, int mxy1,
                     int ntmax, int *irc);

void csse2gppost2lt(float ppart[], float q[], int kpic[], float qm,
                    int nppmx, int idimp, int mx, int my, int nxv,
                    int nyv, int mx1, int mxy1);

void csse2pporder2lt(float ppart[], float ppbuff[], int kpic[],
                     int ncl[], int ihole[], int idimp, int nppmx,
                     int nx, int ny, int mx, int my, int mx1, int my1,
                     int npbmx, int ntmax, int *irc);

void csse2pporderf2lt(float ppart[], float ppbuff[], int kpic[],
                      int ncl[], int ihole[], int idimp, int nppmx,
                      int mx1, int my1, int npbmx, int ntmax,
                      int *irc);

void csse2cguard2l(float fxy[], int nx, int ny, int nxe, int nye);

void csse2aguard2l(float q[], int nx, int ny, int nxe, int nye);

void csse2mpois22(float complex q[], float complex fxy[], int isign,
                  float complex ffc[], float ax, float ay, float affp,
                  float *we, int nx, int ny, int nxvh, int nyv,
                  int nxhd, int nyhd);

void csse2fft2rmxx(float complex f[], int isign, int mixup[],
                   float complex sct[], int indx, int indy, int nyi,
                   int nyp, int nxhd, int nyd, int nxhyd, int nxyhd);

void csse2fft2rmxy(float complex f[], int isign, int mixup[],
                   float complex sct[], int indx, int indy, int nxi,
                   int nxp, int nxhd, int nyd, int nxhyd, int nxyhd);

void csse2fft2rm2x(float complex f[], int isign, int mixup[],
                   float complex sct[], int indx, int indy, int nyi,
                   int nyp, int nxhd, int nyd, int nxhyd, int nxyhd);

void csse2fft2rm2y(float complex f[], int isign, int mixup[],
                   float complex sct[], int indx, int indy, int nxi,
                   int nxp, int nxhd, int nyd, int nxhyd, int nxyhd);

void csse2wfft2rmx(float complex f[], int isign, int mixup[],
                   float complex sct[], int indx, int indy, int nxhd,
                   int nyd, int nxhyd, int nxyhd);

void csse2wfft2rm2(float complex f[], int isign, int mixup[],
                   float complex sct[], int indx, int indy, int nxhd,
                   int nyd, int nxhyd, int nxyhd);
