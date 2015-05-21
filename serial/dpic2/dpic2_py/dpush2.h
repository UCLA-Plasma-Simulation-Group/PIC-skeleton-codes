/* header file for dpush2.c */

double ranorm();

double randum();

void cdistr2h(float part[], float vtx, float vty, float vtz, float vdx,
              float vdy, float vdz, int npx, int npy, int idimp, int nop,
              int nx, int ny, int ipbc);

void cgbpush23l(float part[], float fxy[], float bxy[], float qbm,
                float dt, float dtc, float *ek, int idimp, int nop,
                int nx, int ny, int nxv, int nyv, int ipbc);

void cgpost2l(float part[], float q[], float qm, int nop, int idimp,
              int nxv, int nyv);

void cgjpost2l(float part[], float cu[], float qm, float dt, int nop,
               int idimp, int nx, int ny, int nxv, int nyv, int ipbc);

void cgmjpost2l(float part[], float amu[], float qm, int nop, int idimp,
                int nxv, int nyv);

void cgdjpost2l(float part[], float fxy[], float bxy[], float dcu[],
                float amu[], float qm, float qbm, float dt, int idimp,
                int nop, int nxv, int nyv);

void cgdcjpost2l(float part[], float fxy[], float bxy[], float cu[],
                 float dcu[], float amu[], float qm, float qbm,
                 float dt, int idimp, int nop, int nxv, int nyv);

void cdsortp2yl(float parta[], float partb[], int npic[], int idimp,
                int nop, int ny1);

void cbguard2l(float bxy[], int nx, int ny, int nxe, int nye);

void cacguard2l(float cu[], int nx, int ny, int nxe, int nye);

void caguard2l(float q[], int nx, int ny, int nxe, int nye);

void camcguard2l(float amu[], int nx, int ny, int nxe, int nye,
                 int ndim);

void cascfguard2l(float dcu[], float cus[], float q2m0, int nx, int ny,
                  int nxe, int nye);

void cfwpminmx2(float qe[], float qbme, float *wpmax, float *wpmin,
                int nx, int ny, int nxe, int nye);

void cpois23(float complex q[], float complex fxy[], int isign,
             float complex ffc[], float ax, float ay, float affp,
             float *we, int nx, int ny, int nxvh, int nyv, int nxhd,
             int nyhd);

void ccuperp2(float complex cu[], int nx, int ny, int nxvh, int nyv);

void cbbpois23(float complex cu[], float complex bxy[],
               float complex ffc[], float ci, float *wm, int nx, int ny,
               int nxvh, int nyv, int nxhd, int nyhd);

void cbaddext2(float bxy[], float omx, float omy, float omz, int nx,
               int ny, int nxe, int nye);

void cdcuperp23(float complex dcu[], float complex amu[], int nx,
                int ny, int nxvh, int nyv);

void cadcuperp23(float complex dcu[], float complex amu[], int nx,
                 int ny, int nxvh, int nyv);

void cepois23(float complex dcu[], float complex exy[], int isign,
              float complex ffe[], float ax, float ay, float affp,
              float wp0, float ci, float *wf, int nx, int ny, int nxvh,
              int nyv, int nxhd, int nyhd);

void caddvrfield2(float a[], float b[], float c[], int ndim, int nxe,
                  int nye);

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

void cfft2rnx(float complex f[], float complex ss[], int isign,
              int mixup[], float complex sct[], int indx, int indy,
              int nyi, int nyp, int nxhd, int nyd, int ndim, int nxhyd,
              int nxyhd);

void cfft2rny(float complex f[], int isign, int mixup[],
              float complex sct[], int indx, int indy, int nxi, int nxp,
              int nxhd, int nyd, int ndim, int nxhyd, int nxyhd);

void cwfft2rx(float complex f[], int isign, int mixup[],
              float complex sct[], int indx, int indy, int nxhd, int nyd,
              int nxhyd, int nxyhd);

void cwfft2r3(float complex f[], int isign, int mixup[],
              float complex sct[], int indx, int indy, int nxhd, int nyd,
              int nxhyd, int nxyhd);

void cwfft2rn(float complex f[], float complex ss[], int isign,
              int mixup[], float complex sct[], int indx, int indy,
              int nxhd, int nyd, int ndim, int nxhyd, int nxyhd);

void cswapc2n(float f[], float s[], int isign, int nxh, int nyi,
              int nyt, int nxhd, int nyd, int ndim);
