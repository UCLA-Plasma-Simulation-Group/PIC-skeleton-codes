/* header file for vpush3.c */

double ranorm();

double randum();

void cdistr3t(float part[], float vtx, float vty, float vtz, float vdx,
              float vdy, float vdz, int npx, int npy, int npz,
              int idimp, int npe, int nx, int ny, int nz, int ipbc);

void cgpush3lt(float part[], float fxyz[], float qbm, float dt,
               float *ek, int idimp, int nop, int npe, int nx, int ny,
               int nz, int nxv, int nyv, int nzv, int ipbc);

void cvgpush3lt(float part[], float fxyz[], float qbm, float dt,
                float *ek, int idimp, int nop, int npe, int nx, int ny,
                int nz, int nxv, int nyv, int nzv, int ipbc);

void cgpost3lt(float part[], float q[], float qm, int nop, int npe,
               int idimp, int nxv, int nyv, int nzv);

void cvgpost3lt(float part[], float q[], float qm, int nop, int npe,
                int idimp, int nxv, int nyv, int nzv);

void cdsortp3yzlt(float parta[], float partb[], int npic[], int idimp,
                  int nop, int npe, int ny1, int nyz1);

void ccguard3l(float fxyz[], int nx, int ny, int nz, int nxe, int nye,
               int nze);

void caguard3l(float q[], int nx, int ny, int nz, int nxe, int nye,
               int nze);

void cvpois33(float complex q[], float complex fxyz[], int isign,
              float complex ffc[], float ax, float ay, float az,
              float affp, float *we, int nx, int ny, int nz, int nxvh,
              int nyv, int nzv, int nxhd, int nyhd, int nzhd);

void cwfft3rinit(int mixup[], float complex sct[], int indx, int indy,
                 int indz, int nxhyzd, int nxyzhd);

void cfft3rvxy(float complex f[], int isign, int mixup[],
               float complex sct[], int indx, int indy, int indz,
               int nzi, int nzp, int nxhd, int nyd, int nzd, int nxhyzd,
               int nxyzhd);

void cfft3rxz(float complex f[], int isign, int mixup[],
              float complex sct[], int indx, int indy, int indz,
              int nyi, int nyp, int nxhd, int nyd, int nzd, int nxhyzd,
              int nxyzhd);

void cfft3rv3xy(float complex f[], int isign, int mixup[],
                float complex sct[], int indx, int indy, int indz,
                int nzi, int nzp, int nxhd, int nyd, int nzd,
                int nxhyzd, int nxyzhd);

void cfft3rv3z(float complex f[], int isign, int mixup[],
               float complex sct[], int indx, int indy, int indz,
               int nyi, int nyp, int nxhd, int nyd, int nzd, int nxhyzd,
               int nxyzhd);

void cwfft3rvx(float complex f[], int isign, int mixup[],
               float complex sct[], int indx, int indy, int indz,
               int nxhd, int nyd, int nzd, int nxhyzd, int nxyzhd);

void cwfft3rv3(float complex f[], int isign, int mixup[],
               float complex sct[], int indx, int indy, int indz,
               int nxhd, int nyd, int nzd, int nxhyzd, int nxyzhd);
