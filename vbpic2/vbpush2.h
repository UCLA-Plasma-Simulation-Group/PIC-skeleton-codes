/* header file for vbpush2.c */

double ranorm();

double randum();

void cdistr2ht(float part[], float vtx, float vty, float vtz, float vdx,
               float vdy, float vdz, int npx, int npy, int idimp,
               int npe, int nx, int ny, int ipbc);

void cgbpush23lt(float part[], float fxy[], float bxy[], float qbm,
                 float dt, float dtc, float *ek, int idimp, int nop,
                 int npe, int nx, int ny, int nxv, int nyv, int ipbc);

void cgrbpush23lt(float part[], float fxy[], float bxy[], float qbm,
                  float dt, float dtc, float ci, float *ek, int idimp,
                  int nop, int npe, int nx, int ny, int nxv, int nyv,
                  int ipbc);

void cvgbpush23lt(float part[], float fxy[], float bxy[], float qbm,
                  float dt, float dtc, float *ek, int idimp, int nop,
                  int npe, int nx, int ny, int nxv, int nyv, int ipbc);

void cvgrbpush23lt(float part[], float fxy[], float bxy[], float qbm,
                   float dt, float dtc, float ci, float *ek, int idimp,
                   int nop, int npe, int nx, int ny, int nxv, int nyv,
                   int ipbc);

void cgpost2lt(float part[], float q[], float qm, int nop, int npe,
               int idimp, int nxv, int nyv);

void cvgpost2lt(float part[], float q[], float qm, int nop, int npe,
                int idimp, int nxv, int nyv);

void cgjpost2lt(float part[], float cu[], float qm, float dt, int nop,
                int npe, int idimp, int nx, int ny, int nxv, int nyv,
                int ipbc);

void cgrjpost2lt(float part[], float cu[], float qm, float dt, float ci,
                int nop, int npe, int idimp, int nx, int ny, int nxv,
                int nyv, int ipbc);

void cvgjpost2lt(float part[], float cu[], float qm, float dt, int nop,
                 int npe, int idimp, int nx, int ny, int nxv, int nyv,
                 int ipbc);

void cvgrjpost2lt(float part[], float cu[], float qm, float dt,
                  float ci, int nop, int npe, int idimp, int nx, int ny,
                  int nxv, int nyv, int ipbc);

void cdsortp2ylt(float parta[], float partb[], int npic[], int idimp,
                 int nop, int npe, int ny1);

void cbguard2l(float bxy[], int nx, int ny, int nxe, int nye);

void cacguard2l(float cu[], int nx, int ny, int nxe, int nye);

void caguard2l(float q[], int nx, int ny, int nxe, int nye);

void cvpois23(float complex q[], float complex fxy[], int isign,
              float complex ffc[], float ax, float ay, float affp,
              float *we, int nx, int ny, int nxvh, int nyv, int nxhd,
              int nyhd);

void ccuperp2(float complex cu[], int nx, int ny, int nxvh, int nyv);

void cvibpois23(float complex cu[], float complex bxy[],
                float complex ffc[], float ci, float *wm, int nx, int ny,
                int nxvh, int nyv, int nxhd, int nyhd);

void cvmaxwel2(float complex exy[], float complex bxy[],
               float complex cu[], float complex ffc[], float ci,
               float dt, float *wf, float *wm, int nx, int ny, int nxvh,
               int nyv, int nxhd, int nyhd);

void cvemfield2(float complex fxy[], float complex exy[],
                float complex ffc[], int isign, int nx, int ny, int nxvh,
                int nyv, int nxhd, int nyhd);

void cwfft2rinit(int mixup[], float complex sct[], int indx, int indy, 
                 int nxhyd, int nxyhd);

void cfft2rvxx(float complex f[], int isign, int mixup[],
               float complex sct[], int indx, int indy, int nyi,
               int nyp, int nxhd, int nyd, int nxhyd, int nxyhd);

void cfft2rxy(float complex f[], int isign, int mixup[],
              float complex sct[], int indx, int indy, int nxi, int nxp,
              int nxhd, int nyd, int nxhyd, int nxyhd);

void cfft2rv3x(float complex f[], int isign, int mixup[],
               float complex sct[], int indx, int indy, int nyi,
               int nyp, int nxhd, int nyd, int nxhyd, int nxyhd);

void cfft2rv3y(float complex f[], int isign, int mixup[],
               float complex sct[], int indx, int indy, int nxi,
               int nxp, int nxhd, int nyd, int nxhyd, int nxyhd);

void cwfft2rvx(float complex f[], int isign, int mixup[],
               float complex sct[], int indx, int indy, int nxhd, int nyd,
               int nxhyd, int nxyhd);

void cwfft2rv3(float complex f[],int isign, int mixup[],
               float complex sct[], int indx, int indy, int nxhd,
               int nyd, int nxhyd, int nxyhd);
