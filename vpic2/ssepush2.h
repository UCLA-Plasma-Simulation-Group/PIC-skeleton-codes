/* header file for ssepush2.c */

void csse2xiscan2(int *isdata, int nths);

void csse2gpush2lt(float part[], float fxy[], float qbm, float dt,
                   float *ek, int idimp, int nop, int npe, int nx,
                   int ny, int nxv, int nyv, int ipbc);

void csse2gpost2lt(float part[], float q[], float qm, int nop, int npe,
                   int idimp, int nxv, int nyv);

void csse2dsortp2ylt(float parta[], float partb[], int npic[],
                     int idimp, int nop, int npe, int ny1);

void csse2cguard2l(float fxy[], int nx, int ny, int nxe, int nye);

void csse2aguard2l(float q[], int nx, int ny, int nxe, int nye);

void csse2pois22(float complex q[], float complex fxy[], int isign,
                 float complex ffc[], float ax, float ay, float affp,
                 float *we, int nx, int ny, int nxvh, int nyv, int nxhd,
                 int nyhd);

void csse2fft2rxx(float complex f[], int isign, int mixup[],
                  float complex sct[], int indx, int indy, int nyi,
                  int nyp, int nxhd, int nyd, int nxhyd, int nxyhd);

void csse2fft2rxy(float complex f[], int isign, int mixup[],
                 float complex sct[], int indx, int indy, int nxi, 
                 int nxp, int nxhd, int nyd, int nxhyd, int nxyhd);

void csse2fft2r2x(float complex f[], int isign, int mixup[],
                  float complex sct[], int indx, int indy, int nyi,
                  int nyp, int nxhd, int nyd, int nxhyd, int nxyhd);

void csse2fft2r2y(float complex f[], int isign, int mixup[],
                  float complex sct[], int indx, int indy, int nxi,
                  int nxp, int nxhd, int nyd, int nxhyd, int nxyhd);

void csse2wfft2rx(float complex f[], int isign, int mixup[],
                  float complex sct[], int indx, int indy, int nxhd,
                  int nyd, int nxhyd, int nxyhd);

void csse2wfft2r2(float complex f[], int isign, int mixup[],
                  float complex sct[], int indx, int indy, int nxhd,
                  int nyd, int nxhyd, int nxyhd);
