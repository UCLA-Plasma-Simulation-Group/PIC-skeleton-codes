/* header file for dpush3.c */

double ranorm();

double randum();

void cdistr3(float part[], float vtx, float vty, float vtz, float vdx,
             float vdy, float vdz, int npx, int npy, int npz, int idimp,
             int nop, int nx, int ny, int nz, int ipbc);

void cgbpush3l(float part[], float fxyz[], float bxyz[], float qbm,
               float dt, float dtc, float *ek, int idimp, int nop,
               int nx, int ny, int nz, int nxv, int nyv, int nzv,
               int ipbc);

void cgpost3l(float part[], float q[], float qm, int nop, int idimp,
              int nxv, int nyv, int nzv);

void cgjpost3l(float part[], float cu[], float qm, float dt, int nop,
               int idimp, int nx, int ny, int nz, int nxv, int nyv,
               int nzv, int ipbc);

void cgmjpost3l(float part[], float amu[], float qm, int nop, int idimp,
                int nxv, int nyv, int nzv);

void cgdjpost3l(float part[], float fxyz[], float bxyz[], float dcu[],
                float amu[], float qm, float qbm, float dt, int idimp,
                int nop, int nxv, int nyv, int nzv);

void cgdcjpost3l(float part[], float fxyz[], float bxyz[], float cu[],
                 float dcu[], float amu[], float qm, float qbm,
                 float dt, int idimp, int nop, int nxv, int nyv,
                 int nzv);

void cdsortp3yzl(float parta[], float partb[], int npic[], int idimp,
                 int nop, int ny1, int nyz1);

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

void cpois33(float complex q[], float complex fxyz[], int isign,
             float complex ffc[], float ax, float ay, float az,
             float affp, float *we, int nx, int ny, int nz, int nxvh,
             int nyv, int nzv, int nxhd, int nyhd, int nzhd);

void ccuperp3(float complex cu[], int nx, int ny, int nz, int nxvh,
              int nyv, int nzv);

void cbbpois33(float complex cu[], float complex bxyz[],
               float complex ffc[], float ci, float *wm, int nx, int ny,
               int nz, int nxvh, int nyv, int nzv, int nxhd, int nyhd,
               int nzhd);

void cbaddext3(float bxyz[], float omx, float omy, float omz, int nx,
               int ny, int nz, int nxe, int nye, int nze);

void cdcuperp3(float complex dcu[], float complex amu[], int nx, int ny,
               int nz, int nxvh, int nyv, int nzv);

void cadcuperp3(float complex dcu[], float complex amu[], int nx,
                int ny, int nz, int nxvh, int nyv, int nzv);

void cepois33(float complex dcu[], float complex exyz[], int isign,
              float complex ffe[], float ax, float ay, float az,
              float affp, float wp0, float ci, float *wf, int nx, int ny,
              int nz, int nxvh, int nyv, int nzv, int nxhd, int nyhd,
              int nzhd);

void caddvrfield3(float a[], float b[], float c[], int ndim, int nxe,
                  int nye, int nze);

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

void cfft3rnxy(float complex f[], float complex ss[], int isign,
               int mixup[], float complex sct[], int indx, int indy,
               int indz, int nzi, int nzp, int nxhd, int nyd, int nzd,
               int ndim, int nxhyzd, int nxyzhd);

void cfft3rnz(float complex f[], int isign, int mixup[],
              float complex sct[], int indx, int indy, int indz,
              int nyi, int nyp, int nxhd, int nyd, int nzd, int ndim,
              int nxhyzd, int nxyzhd);

void cwfft3rx(float complex f[], int isign, int mixup[],
              float complex sct[], int indx, int indy, int indz,
              int nxhd, int nyd, int nzd, int nxhyzd, int nxyzhd);

void cwfft3r3(float complex f[], int isign, int mixup[],
              float complex sct[], int indx, int indy, int indz,
              int nxhd, int nyd, int nzd, int nxhyzd, int nxyzhd);

void cwfft3rn(float complex f[], float complex ss[], int isign,
              int mixup[], float complex sct[], int indx, int indy,
              int indz, int nxhd, int nyd, int nzd, int ndim, 
              int nxhyzd, int nxyzhd);

void cswap3cn(float f[], float s[], int isign, int nxh, int ny, int nzi,
              int nzt, int nxhd, int nyd, int nzd, int ndim);
