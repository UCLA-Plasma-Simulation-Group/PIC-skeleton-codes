/* header file for vpush2.c */

double ranorm();

double randum();

void cdistr2t(float part[], float vtx, float vty, float vdx, float vdy,
              int npx, int npy, int idimp, int npe, int nx, int ny,
              int ipbc);

void cgpush2lt(float part[], float fxy[], float qbm, float dt,
               float *ek, int idimp, int nop, int npe, int nx, int ny,
               int nxv, int nyv, int ipbc);

void cvgpush2lt(float part[], float fxy[], float qbm, float dt,
                float *ek, int idimp, int nop, int npe, int nx, int ny,
                int nxv, int nyv, int ipbc);

void cgpost2lt(float part[], float q[], float qm, int nop, int npe,
               int idimp, int nxv, int nyv);

void cvgpost2lt(float part[], float q[], float qm, int nop, int npe,
                int idimp, int nxv, int nyv);

void cdsortp2ylt(float parta[], float partb[], int npic[], int idimp,
                 int nop, int npe, int ny1);

void ccguard2l(float fxy[], int nx, int ny, int nxe, int nye);

void caguard2l(float q[], int nx, int ny, int nxe, int nye);

void cvpois22(float complex q[], float complex fxy[], int isign,
              float complex ffc[], float ax, float ay, float affp,
              float *we, int nx, int ny, int nxvh, int nyv, int nxhd,
              int nyhd);

void cwfft2rinit(int mixup[], float complex sct[], int indx, int indy, 
                 int nxhyd, int nxyhd);

void cfft2rvxx(float complex f[], int isign, int mixup[],
               float complex sct[], int indx, int indy, int nyi,
               int nyp, int nxhd, int nyd, int nxhyd, int nxyhd);

void cfft2rxy(float complex f[], int isign, int mixup[],
              float complex sct[], int indx, int indy, int nxi, int nxp,
              int nxhd, int nyd, int nxhyd, int nxyhd);

void cfft2rv2x(float complex f[], int isign, int mixup[],
               float complex sct[], int indx, int indy, int nyi,
               int nyp, int nxhd, int nyd, int nxhyd, int nxyhd);

void cfft2r2y(float complex f[], int isign, int mixup[],
              float complex sct[], int indx, int indy, int nxi, int nxp,
              int nxhd, int nyd, int nxhyd, int nxyhd);

void cwfft2rvx(float complex f[], int isign, int mixup[],
               float complex sct[], int indx, int indy, int nxhd, int nyd,
               int nxhyd, int nxyhd);

void cwfft2rv2(float complex f[], int isign, int mixup[],
               float complex sct[], int indx, int indy, int nxhd,
               int nyd, int nxhyd, int nxyhd);

