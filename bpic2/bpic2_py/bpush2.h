/* header file for bpush2.c */

double ranorm();

double randum();

void cdistr2h(float part[], float vtx, float vty, float vtz, float vdx,
              float vdy, float vdz, int npx, int npy, int idimp, int nop,
              int nx, int ny, int ipbc);

void cgbpush23l(float part[], float fxy[], float bxy[], float qbm,
                float dt, float dtc, float *ek, int idimp, int nop,
                int nx, int ny, int nxv, int nyv, int ipbc);

void cdgbpush23l(double part[], double fxy[], double bxy[], double qbm,
                 double dt, double dtc, double *ek, int idimp, int nop,
                 int nx, int ny, int nxv, int nyv, int ipbc);

void cgrbpush23l(float part[], float fxy[], float bxy[], float qbm,
                 float dt, float dtc, float ci, float *ek, int idimp,
                 int nop, int nx, int ny, int nxv, int nyv, int ipbc);

void cdgrbpush23l(double part[], double fxy[], double bxy[], double qbm,
                  double dt, double dtc, double ci, double *ek, int idimp,
                  int nop, int nx, int ny, int nxv, int nyv, int ipbc);

void cgpost2l(float part[], float q[], float qm, int nop, int idimp,
              int nxv, int nyv);

void cgjpost2l(float part[], float cu[], float qm, float dt, int nop,
               int idimp, int nx, int ny, int nxv, int nyv, int ipbc);

void cgrjpost2l(float part[], float cu[], float qm, float dt, float ci,
                int nop, int idimp, int nx, int ny, int nxv, int nyv,
                int ipbc);

void cdsortp2yl(float parta[], float partb[], int npic[], int idimp,
                int nop, int ny1);

void cbguard2l(float bxy[], int nx, int ny, int nxe, int nye);

void cacguard2l(float cu[], int nx, int ny, int nxe, int nye);

void caguard2l(float q[], int nx, int ny, int nxe, int nye);

void cpois23(float complex q[], float complex fxy[], int isign,
             float complex ffc[], float ax, float ay, float affp,
             float *we, int nx, int ny, int nxvh, int nyv, int nxhd,
             int nyhd);

void ccuperp2(float complex cu[], int nx, int ny, int nxvh, int nyv);

void cibpois23(float complex cu[], float complex bxy[],
               float complex ffc[], float ci, float *wm, int nx, int ny,
               int nxvh, int nyv, int nxhd, int nyhd);

void cmaxwel2(float complex exy[], float complex bxy[],
              float complex cu[], float complex ffc[], float ci,
              float dt, float *wf, float *wm, int nx, int ny, int nxvh,
              int nyv, int nxhd, int nyhd);

void cemfield2(float complex fxy[], float complex exy[],
               float complex ffc[], int isign, int nx, int ny, int nxvh,
               int nyv, int nxhd, int nyhd);

void cwfft2rinit(int mixup[], float complex sct[], int indx, int indy, 
                 int nxhyd, int nxyhd);

void cfft2rxx(float complex f[], int isign, int mixup[],
              float complex sct[], int indx, int indy, int nyi, int nyp,
              int nxhd, int nyd, int nxhyd, int nxyhd);

void cfft2rxy(float complex f[], int isign, int mixup[],
              float complex sct[], int indx, int indy, int nxi, int nxp,
              int nxhd, int nyd, int nxhyd, int nxyhd);

void cfft2r3x(float complex f[], int isign, int mixup[],
              float complex sct[], int indx, int indy, int nyi, int nyp,
              int nxhd, int nyd, int nxhyd, int nxyhd);

void cfft2r3y(float complex f[], int isign, int mixup[],
              float complex sct[], int indx, int indy, int nxi, int nxp,
              int nxhd, int nyd, int nxhyd, int nxyhd);

void cwfft2rx(float complex f[], int isign, int mixup[],
              float complex sct[], int indx, int indy, int nxhd, int nyd,
              int nxhyd, int nxyhd);

void cwfft2r3(float complex f[],int isign, int mixup[],
              float complex sct[], int indx, int indy, int nxhd, int nyd,
              int nxhyd, int nxyhd);
