/* header file for mpush2.c */

double ranorm();

void cdistr2(float part[], float vtx, float vty, float vdx, float vdy,
             int npx, int npy, int idimp, int nop, int nx, int ny,
             int ipbc);

void cdblkp2l(float part[], int kpic[], int *nppmx, int idimp, int nop,
              int mx, int my, int mx1, int mxy1, int *irc);

void cppmovin2l(float part[], float ppart[], int kpic[], int nppmx,
                int idimp, int nop, int mx, int my, int mx1, int mxy1,
                int *irc);

void cppcheck2l(float ppart[], int kpic[], int idimp, int nppmx, int nx,
                int ny, int mx, int my, int mx1, int my1, 
                int *irc);

void cgppush2l(float ppart[], float fxy[], int kpic[], float qbm,
               float dt, float *ek, int idimp, int nppmx, int nx,
               int ny, int mx, int my, int nxv, int nyv, int mx1,
               int mxy1, int ipbc);

void cgppushf2l(float ppart[], float fxy[], int kpic[], int ncl[],
                int ihole[], float qbm, float dt, float *ek, int idimp,
                int nppmx, int nx, int ny, int mx, int my, int nxv,
                int nyv, int mx1, int mxy1, int ntmax, int *irc);

void cgppost2l(float ppart[], float q[], int kpic[], float qm,
               int nppmx, int idimp, int mx, int my, int nxv, int nyv,
               int mx1, int mxy1);

void cpporder2l(float ppart[], float ppbuff[], int kpic[], int ncl[],
                int ihole[], int idimp, int nppmx, int nx, int ny,
                int mx, int my, int mx1, int my1, int npbmx, int ntmax,
                int *irc);

void cpporderf2l(float ppart[], float ppbuff[], int kpic[], int ncl[],
                 int ihole[], int idimp, int nppmx, int mx1, int my1,
                 int npbmx, int ntmax, int *irc);

void ccguard2l(float fxy[], int nx, int ny, int nxe, int nye);

void caguard2l(float q[], int nx, int ny, int nxe, int nye);

void cppcheck2l(float ppart[], int kpic[], int idimp, int nppmx, int nx,
                int ny, int mx, int my, int mx1, int my1, 
                int *irc);

void cmpois22(float complex q[], float complex fxy[], int isign,
              float complex ffc[], float ax, float ay, float affp,
              float *we, int nx, int ny, int nxvh, int nyv, int nxhd,
              int nyhd);

void cwfft2rinit(int mixup[], float complex sct[], int indx, int indy, 
                 int nxhyd, int nxyhd);

void cfft2rmxx(float complex f[], int isign, int mixup[],
               float complex sct[], int indx, int indy, int nyi,
               int nyp, int nxhd, int nyd, int nxhyd, int nxyhd);

void cfft2rmxy(float complex f[], int isign, int mixup[],
               float complex sct[], int indx, int indy, int nxi,
               int nxp, int nxhd, int nyd, int nxhyd, int nxyhd);

void cfft2r2x(float complex f[], int isign, int mixup[],
              float complex sct[], int indx, int indy, int nyi, int nyp,
              int nxhd, int nyd, int nxhyd, int nxyhd);

void cfft2r2y(float complex f[], int isign, int mixup[],
              float complex sct[], int indx, int indy, int nxi, int nxp,
              int nxhd, int nyd, int nxhyd, int nxyhd);

void cwfft2rmx(float complex f[], int isign, int mixup[],
               float complex sct[], int indx, int indy, int nxhd,
               int nyd, int nxhyd, int nxyhd);

void cwfft2rm2(float complex f[], int isign, int mixup[],
               float complex sct[], int indx, int indy, int nxhd,
               int nyd, int nxhyd, int nxyhd);
