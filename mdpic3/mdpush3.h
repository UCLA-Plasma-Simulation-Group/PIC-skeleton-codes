/* header file for mdpush3.c */

double ranorm();

void cdistr3(float part[], float vtx, float vty, float vtz, float vdx,
             float vdy, float vdz, int npx, int npy, int npz, int idimp,
             int nop, int nx, int ny, int nz, int ipbc);

void cdblkp3l(float part[], int kpic[], int *nppmx, int idimp, int nop,
              int mx, int my, int mz, int mx1, int my1, int mxyz1,
              int *irc);

void cppmovin3l(float part[], float ppart[], int kpic[], int nppmx,
                int idimp, int nop, int mx, int my, int mz, int mx1,
                int my1, int mxyz1, int *irc);

void cppcheck3l(float ppart[], int kpic[], int idimp, int nppmx, int nx,
                int ny, int nz, int mx, int my, int mz, int mx1,
                int my1, int mz1, int *irc);

void cgbppush3l(float ppart[], float fxyz[], float bxyz[], int kpic[],
                float qbm, float dt, float dtc, float *ek, int idimp,
                int nppmx, int nx, int ny, int nz, int mx, int my,
                int mz, int nxv, int nyv, int nzv, int mx1, int my1,
                int mxyz1, int ipbc);

void cgbppushf3l(float ppart[], float fxyz[], float bxyz[], int kpic[],
                 int ncl[], int ihole[], float qbm, float dt, float dtc,
                 float *ek, int idimp, int nppmx, int nx, int ny,
                 int nz, int mx, int my, int mz, int nxv, int nyv,
                 int nzv, int mx1, int my1, int mxyz1, int ntmax,
                 int *irc);

void cgppost3l(float ppart[], float q[], int kpic[], float qm,
               int nppmx, int idimp, int mx, int my, int mz, int nxv,
               int nyv, int nzv, int mx1, int my1, int mxyz1);

void cgjppost3l(float ppart[], float cu[], int kpic[], float qm,
                float dt, int nppmx, int idimp, int nx, int ny, int nz,
                int mx, int my, int mz, int nxv, int nyv, int nzv,
                int mx1, int my1, int mxyz1, int ipbc);

void cgmjppost3l(float ppart[], float amu[], int kpic[], float qm,
                 int nppmx, int idimp, int mx, int my, int mz, int nxv,
                 int nyv, int nzv, int mx1, int my1, int mxyz1);

void cgdjppost3l(float ppart[], float fxyz[], float bxyz[], int kpic[], 
                 float dcu[], float amu[], float qm, float qbm,
                 float dt, int idimp, int nppmx, int nx, int ny, int nz,
                 int mx, int my, int mz, int nxv, int nyv, int nzv,
                 int mx1, int my1, int mxyz1);

void cgdcjppost3l(float ppart[], float fxyz[], float bxyz[], int kpic[],
                  float cu[], float dcu[], float amu[], float qm,
                  float qbm, float dt, int idimp, int nppmx, int nx,
                  int ny, int nz, int mx, int my, int mz, int nxv,
                  int nyv, int nzv, int mx1, int my1, int mxyz1);

void cpporder3l(float ppart[], float ppbuff[], int kpic[], int ncl[],
                int ihole[], int idimp, int nppmx, int nx, int ny,
                int nz, int mx, int my, int mz, int mx1, int my1,
                int mz1, int npbmx, int ntmax, int *irc);

void cpporderf3l(float ppart[], float ppbuff[], int kpic[], int ncl[],
                 int ihole[], int idimp, int nppmx, int mx1, int my1,
                 int mz1, int npbmx, int ntmax, int *irc);

void ccguard3l(float fxyz[], int nx, int ny, int nz, int nxe, int nye,
               int nze);

void cacguard3l(float cu[], int nx, int ny, int nz, int nxe, int nye,
                int nze);

void caguard3l(float q[], int nx, int ny, int nz, int nxe, int nye,
               int nze);

void camcguard3l(float amu[], int nx, int ny, int nz, int nxe, int nye,
                 int nze, int ndim);

void cascfguard3l(float dcu[], float cus[], float q2m0, int nx, int ny,
                  int nz, int nxe, int nye, int nze);

void cfwpminmx3(float qe[], float qbme, float *wpmax, float *wpmin,
                int nx, int ny, int nz, int nxe, int nye, int nze);

void cmpois33(float complex q[], float complex fxyz[], int isign,
              float complex ffc[], float ax, float ay, float az,
              float affp, float *we, int nx, int ny, int nz, int nxvh,
              int nyv, int nzv, int nxhd, int nyhd, int nzhd);

void cmcuperp3(float complex cu[], int nx, int ny, int nz, int nxvh,
               int nyv, int nzv);

void cmbbpois33(float complex cu[], float complex bxyz[],
                float complex ffc[], float ci, float *wm, int nx,
                int ny, int nz, int nxvh, int nyv, int nzv, int nxhd,
                int nyhd, int nzhd);

void cbaddext3(float bxyz[], float omx, float omy, float omz, int nx,
               int ny, int nz, int nxe, int nye, int nze);

void cmdcuperp3(float complex dcu[], float complex amu[], int nx,
                int ny, int nz, int nxvh, int nyv, int nzv);

void cmadcuperp3(float complex dcu[], float complex amu[], int nx,
                 int ny, int nz, int nxvh, int nyv, int nzv);

void cmepois33(float complex dcu[], float complex exyz[], int isign,
               float complex ffe[], float ax, float ay, float az,
               float affp, float wp0, float ci, float *wf, int nx,
               int ny, int nz, int nxvh, int nyv, int nzv, int nxhd,
               int nyhd, int nzhd);

void caddvrfield3(float a[], float b[], float c[], int ndim, int nxe,
                  int nye, int nze);

void cwfft3rinit(int mixup[], float complex sct[], int indx, int indy,
                 int indz, int nxhyzd, int nxyzhd);

void cfft3rmxy(float complex f[], int isign, int mixup[],
               float complex sct[], int indx, int indy, int indz,
               int nzi, int nzp, int nxhd, int nyd, int nzd, int nxhyzd,
               int nxyzhd);

void cfft3rmxz(float complex f[], int isign, int mixup[],
               float complex sct[], int indx, int indy, int indz,
               int nyi, int nyp, int nxhd, int nyd, int nzd, int nxhyzd,
               int nxyzhd);

void cfft3rm3xy(float complex f[], int isign, int mixup[],
                float complex sct[], int indx, int indy, int indz,
                int nzi, int nzp, int nxhd, int nyd, int nzd,
                int nxhyzd, int nxyzhd);

void cfft3rm3z(float complex f[], int isign, int mixup[],
               float complex sct[], int indx, int indy, int indz,
               int nyi, int nyp, int nxhd, int nyd, int nzd, int nxhyzd,
               int nxyzhd);

void cfft3rmnxy(float complex f[], float complex ss[], int isign,
                int mixup[], float complex sct[], int indx, int indy,
                int indz, int nzi, int nzp, int nxhd, int nyd, int nzd,
                int ndim, int nxhyzd, int nxyzhd);

void cfft3rmnz(float complex f[], int isign, int mixup[],
               float complex sct[], int indx, int indy, int indz,
               int nyi, int nyp, int nxhd, int nyd, int nzd, int ndim,
               int nxhyzd, int nxyzhd);

void cwfft3rmx(float complex f[], int isign, int mixup[],
               float complex sct[], int indx, int indy, int indz,
               int nxhd, int nyd, int nzd, int nxhyzd, int nxyzhd);

void cwfft3rm3(float complex f[], int isign, int mixup[],
               float complex sct[], int indx, int indy, int indz,
               int nxhd, int nyd, int nzd, int nxhyzd, int nxyzhd);

void cwfft3rmn(float complex f[], float complex ss[], int isign,
               int mixup[], float complex sct[], int indx, int indy,
               int indz, int nxhd, int nyd, int nzd, int ndim, 
               int nxhyzd, int nxyzhd);

void cmswap3cn(float f[], float s[], int isign, int nxh, int ny,
               int nzi, int nzt, int nxhd, int nyd, int nzd, int ndim);
