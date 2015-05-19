/* header file for push3.c */

double ranorm();

double randum();

void cdistr3(float part[], float vtx, float vty, float vtz, float vdx,
             float vdy, float vdz, int npx, int npy, int npz, int idimp,
             int nop, int nx, int ny, int nz, int ipbc);

void cgpush3l(float part[], float fxyz[], float qbm, float dt,
              float *ek, int idimp, int nop, int nx, int ny, int nz,
              int nxv, int nyv, int nzv, int ipbc);

void cgpost3l(float part[], float q[], float qm, int nop, int idimp,
              int nxv, int nyv, int nzv);

void cdsortp3yzl(float parta[], float partb[], int npic[], int idimp,
                 int nop, int ny1, int nyz1);

void ccguard3l(float fxyz[], int nx, int ny, int nz, int nxe, int nye,
               int nze);

void caguard3l(float q[], int nx, int ny, int nz, int nxe, int nye,
               int nze);

void cpois33(float complex q[], float complex fxyz[], int isign,
             float complex ffc[], float ax, float ay, float az,
             float affp, float *we, int nx, int ny, int nz, int nxvh,
             int nyv, int nzv, int nxhd, int nyhd, int nzhd);

void cwfft3rinit(int mixup[], float complex sct[], int indx, int indy,
                 int indz, int nxhyzd, int nxyzhd);

void cfft3rxy(float complex f[], int isign, int mixup[],
              float complex sct[], int indx, int indy, int indz,
              int nzi, int nzp, int nxhd, int nyd, int nzd, int nxhyzd,
              int nxyzhd);

void cfft3rxz(float complex f[], int isign, int mixup[],
              float complex sct[], int indx, int indy, int indz,
              int nyi, int nyp, int nxhd, int nyd, int nzd, int nxhyzd,
              int nxyzhd);

void cfft3r3xy(float complex f[], int isign, int mixup[],
               float complex sct[], int indx, int indy, int indz,
               int nzi, int nzp, int nxhd, int nyd, int nzd, int nxhyzd,
               int nxyzhd);

void cfft3r3z(float complex f[], int isign, int mixup[],
              float complex sct[], int indx, int indy, int indz,
              int nyi, int nyp, int nxhd, int nyd, int nzd, int nxhyzd, 
              int nxyzhd);

void cwfft3rx(float complex f[], int isign, int mixup[],
              float complex sct[], int indx, int indy, int indz,
              int nxhd, int nyd, int nzd, int nxhyzd, int nxyzhd);

void cwfft3r3(float complex f[], int isign, int mixup[],
              float complex sct[], int indx, int indy, int indz,
              int nxhd, int nyd, int nzd, int nxhyzd, int nxyzhd);
