/* header file for ssebpush2.c */

void csse2xiscan2(int *isdata, int nths);

void csse2gbpush23lt(float part[], float fxy[], float bxy[], float qbm,
                     float dt, float dtc, float *ek, int idimp, int nop,
                     int npe, int nx, int ny, int nxv, int nyv,
                     int ipbc);

void csse2grbpush23lt(float part[], float fxy[], float bxy[], float qbm,
                      float dt, float dtc, float ci, float *ek,
                      int idimp, int nop, int npe, int nx, int ny,
                      int nxv, int nyv, int ipbc);

void csse2gpost2lt(float part[], float q[], float qm, int nop, int npe,
                   int idimp, int nxv, int nyv);

void csse2gjpost2lt(float part[], float cu[], float qm, float dt,
                    int nop, int npe, int idimp, int nx, int ny, 
                    int nxv, int nyv, int ipbc);

void csse2grjpost2lt(float part[], float cu[], float qm, float dt,
                     float ci, int nop, int npe, int idimp, int nx,
                     int ny, int nxv, int nyv, int ipbc);

void csse2dsortp2ylt(float parta[], float partb[], int npic[],
                     int idimp, int nop, int npe, int ny1);

void csse2bguard2l(float bxy[], int nx, int ny, int nxe, int nye);

void csse2acguard2l(float cu[], int nx, int ny, int nxe, int nye);

void csse2aguard2l(float q[], int nx, int ny, int nxe, int nye);

void csse2pois23(float complex q[], float complex fxy[], int isign,
                 float complex ffc[], float ax, float ay, float affp,
                 float *we, int nx, int ny, int nxvh, int nyv, int nxhd,
                 int nyhd);

void csse2cuperp2(float complex cu[], int nx, int ny, int nxvh,
                  int nyv);

void csse2ibpois23(float complex cu[], float complex bxy[],
                   float complex ffc[], float ci, float *wm, int nx,
                   int ny, int nxvh, int nyv, int nxhd, int nyhd);

void csse2maxwel2(float complex exy[], float complex bxy[],
                  float complex cu[], float complex ffc[], float ci,
                  float dt, float *wf, float *wm, int nx, int ny,
                  int nxvh, int nyv, int nxhd, int nyhd);

void csse2emfield2(float complex fxy[], float complex exy[],
                   float complex ffc[], int isign, int nx, int ny,
                   int nxvh, int nyv, int nxhd, int nyhd);

void csse2fft2rxx(float complex f[], int isign, int mixup[],
                  float complex sct[], int indx, int indy, int nyi,
                  int nyp, int nxhd, int nyd, int nxhyd, int nxyhd);

void csse2fft2rxy(float complex f[], int isign, int mixup[],
                 float complex sct[], int indx, int indy, int nxi, 
                 int nxp, int nxhd, int nyd, int nxhyd, int nxyhd);

void csse2fft2r3x(float complex f[], int isign, int mixup[],
                  float complex sct[], int indx, int indy, int nyi,
                  int nyp, int nxhd, int nyd, int nxhyd, int nxyhd);

void csse2fft2r3y(float complex f[], int isign, int mixup[],
                  float complex sct[], int indx, int indy, int nxi,
                  int nxp, int nxhd, int nyd, int nxhyd, int nxyhd);

void csse2wfft2rx(float complex f[], int isign, int mixup[],
                  float complex sct[], int indx, int indy, int nxhd,
                  int nyd, int nxhyd, int nxyhd);

void csse2wfft2r3(float complex f[],int isign, int mixup[],
                  float complex sct[], int indx, int indy, int nxhd,
                  int nyd, int nxhyd, int nxyhd);
