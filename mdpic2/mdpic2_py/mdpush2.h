/* header file for mdpush2.c */

double ranorm();

void cdistr2h(float part[], float vtx, float vty, float vtz, float vdx,
              float vdy, float vdz, int npx, int npy, int idimp, int nop,
              int nx, int ny, int ipbc);

void cdblkp2l(float part[], int kpic[], int *nppmx, int idimp, int nop,
              int mx, int my, int mx1, int mxy1, int *irc);

void cppmovin2l(float part[], float ppart[], int kpic[], int nppmx,
                int idimp, int nop, int mx, int my, int mx1, int mxy1,
                int *irc);

void cppcheck2l(float ppart[], int kpic[], int idimp, int nppmx, int nx,
                int ny, int mx, int my, int mx1, int my1, 
                int *irc);

void cgbppush23l(float ppart[], float fxy[], float bxy[], int kpic[],
                 float qbm, float dt, float dtc, float *ek, int idimp,
                 int nppmx, int nx, int ny, int mx, int my, int nxv,
                 int nyv, int mx1, int mxy1, int ipbc);

void cgbppushf23l(float ppart[], float fxy[], float bxy[], int kpic[],
                  int ncl[], int ihole[], float qbm, float dt,
                  float dtc, float *ek, int idimp, int nppmx, int nx,
                  int ny, int mx, int my, int nxv, int nyv, int mx1,
                  int mxy1, int ntmax, int *irc);

void cgppost2l(float ppart[], float q[], int kpic[], float qm,
               int nppmx, int idimp, int mx, int my, int nxv, int nyv,
               int mx1, int mxy1);

void cgjppost2l(float ppart[], float cu[], int kpic[], float qm,
                float dt, int nppmx, int idimp, int nx, int ny, int mx,
                int my, int nxv, int nyv, int mx1, int mxy1, int ipbc);

void cgmjppost2l(float ppart[], float amu[], int kpic[], float qm,
                 int nppmx, int idimp, int mx, int my, int nxv, int nyv,
                 int mx1, int mxy1);

void cgdjppost2l(float ppart[], float fxy[], float bxy[], float dcu[],
                 float amu[], int kpic[], float qm, float qbm, float dt,
                 int idimp, int nppmx, int nx, int ny, int mx, int my,
                 int nxv, int nyv, int mx1, int mxy1);

void cgdcjppost2l(float ppart[], float fxy[], float bxy[], float cu[],
                  float dcu[], float amu[], int kpic[], float qm,
                  float qbm, float dt, int idimp, int nppmx, int nx,
                  int ny, int mx, int my, int nxv, int nyv, int mx1,
                  int mxy1);

void cpporder2l(float ppart[], float ppbuff[], int kpic[], int ncl[],
                int ihole[], int idimp, int nppmx, int nx, int ny,
                int mx, int my, int mx1, int my1, int npbmx, int ntmax,
                int *irc);

void cpporderf2l(float ppart[], float ppbuff[], int kpic[], int ncl[],
                 int ihole[], int idimp, int nppmx, int mx1, int my1,
                 int npbmx, int ntmax, int *irc);

void cbguard2l(float bxy[], int nx, int ny, int nxe, int nye);

void cacguard2l(float cu[], int nx, int ny, int nxe, int nye);

void caguard2l(float q[], int nx, int ny, int nxe, int nye);

void camcguard2l(float amu[], int nx, int ny, int nxe, int nye,
                 int ndim);

void cascfguard2l(float dcu[], float cus[], float q2m0, int nx, int ny,
                  int nxe, int nye);

void cfwpminmx2(float qe[], float qbme, float *wpmax, float *wpmin,
                int nx, int ny, int nxe, int nye);

void cmpois23(float complex q[], float complex fxy[], int isign,
              float complex ffc[], float ax, float ay, float affp,
              float *we, int nx, int ny, int nxvh, int nyv, int nxhd,
              int nyhd);

void cmcuperp2(float complex cu[], int nx, int ny, int nxvh, int nyv);

void cmbbpois23(float complex cu[], float complex bxy[],
                float complex ffc[], float ci, float *wm, int nx,
                int ny, int nxvh, int nyv, int nxhd, int nyhd);

void cbaddext2(float bxy[], float omx, float omy, float omz, int nx,
               int ny, int nxe, int nye);

void cmdcuperp23(float complex dcu[], float complex amu[], int nx,
                 int ny, int nxvh, int nyv);

void cmadcuperp23(float complex dcu[], float complex amu[], int nx,
                  int ny, int nxvh, int nyv);

void cmepois23(float complex dcu[], float complex exy[], int isign,
               float complex ffe[], float ax, float ay, float affp,
               float wp0, float ci, float *wf, int nx, int ny, int nxvh,
               int nyv, int nxhd, int nyhd);

void caddvrfield2(float a[], float b[], float c[], int ndim, int nxe,
                  int nye);

void cwfft2rinit(int mixup[], float complex sct[], int indx, int indy, 
                 int nxhyd, int nxyhd);

void cfft2rmxx(float complex f[], int isign, int mixup[],
               float complex sct[], int indx, int indy, int nyi,
               int nyp, int nxhd, int nyd, int nxhyd, int nxyhd);

void cfft2rmxy(float complex f[], int isign, int mixup[],
               float complex sct[], int indx, int indy, int nxi,
               int nxp, int nxhd, int nyd, int nxhyd, int nxyhd);

void cfft2rm3x(float complex f[], int isign, int mixup[],
               float complex sct[], int indx, int indy, int nyi,
               int nyp, int nxhd, int nyd, int nxhyd, int nxyhd);

void cfft2rm3y(float complex f[], int isign, int mixup[],
               float complex sct[], int indx, int indy, int nxi,
               int nxp, int nxhd, int nyd, int nxhyd, int nxyhd);

void cfft2rmnx(float complex f[], float complex ss[], int isign,
               int mixup[], float complex sct[], int indx, int indy,
               int nyi, int nyp, int nxhd, int nyd, int ndim, int nxhyd,
               int nxyhd);

void cfft2rmny(float complex f[], int isign, int mixup[],
               float complex sct[], int indx, int indy, int nxi,
               int nxp, int nxhd, int nyd, int ndim, int nxhyd,
               int nxyhd);

void cfft2rnx(float complex f[], float complex ss[], int isign,
              int mixup[], float complex sct[], int indx, int indy,
              int nyi, int nyp, int nxhd, int nyd, int ndim, int nxhyd,
              int nxyhd);

void cfft2rny(float complex f[], int isign, int mixup[],
              float complex sct[], int indx, int indy, int nxi, int nxp,
              int nxhd, int nyd, int ndim, int nxhyd, int nxyhd);

void cwfft2rmx(float complex f[], int isign, int mixup[],
               float complex sct[], int indx, int indy, int nxhd,
               int nyd, int nxhyd, int nxyhd);

void cwfft2rm3(float complex f[], int isign, int mixup[],
               float complex sct[], int indx, int indy, int nxhd,
               int nyd, int nxhyd, int nxyhd);

void cwfft2rmn(float complex f[], float complex ss[], int isign,
               int mixup[], float complex sct[], int indx, int indy,
               int nxhd, int nyd, int ndim, int nxhyd, int nxyhd);

void cmswapc2n(float f[], float s[], int isign, int nxh, int nyi,
               int nyt, int nxhd, int nyd, int ndim);


