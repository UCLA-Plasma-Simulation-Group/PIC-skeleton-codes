/* header file for kncpush3.c */

void ckncxiscan2(int *isdata, int nths);

void ckncgpush3lt(float part[], float fxyz[], float qbm, float dt,
                  float *ek, int idimp, int nop, int npe, int nx,
                  int ny, int nz, int nxv, int nyv, int nzv, int ipbc);

void ckncgpost3lt(float part[], float q[], float qm, int nop, int npe,
                  int idimp, int nxv, int nyv, int nzv);

void cknc2gpost3lt(float part[], float q[], float qm, int nop, int npe,
                   int idimp, int nxv, int nyv, int nzv);

void ckncdsortp3yzlt(float parta[], float partb[], int npic[],
                     int idimp, int nop, int npe, int ny1, int nyz1);

void cknccguard3l(float fxyz[], int nx, int ny, int nz, int nxe,
                  int nye, int nze);

void ckncaguard3l(float q[], int nx, int ny, int nz, int nxe, int nye,
                  int nze);

void ckncpois33(float complex q[], float complex fxyz[], int isign,
                float complex ffc[], float ax, float ay, float az,
                float affp, float *we, int nx, int ny, int nz, int nxvh,
                int nyv, int nzv, int nxhd, int nyhd, int nzhd);

void ckncfft3rvxy(float complex f[], int isign, int mixup[],
                  float complex sct[], int indx, int indy, int indz,
                  int nzi, int nzp, int nxhd, int nyd, int nzd,
                  int nxhyzd, int nxyzhd);

void ckncfft3rxz(float complex f[], int isign, int mixup[],
                 float complex sct[], int indx, int indy, int indz,
                 int nyi, int nyp, int nxhd, int nyd, int nzd,
                 int nxhyzd, int nxyzhd);

void ckncfft3rv3xy(float complex f[], int isign, int mixup[],
                   float complex sct[], int indx, int indy, int indz,
                   int nzi, int nzp, int nxhd, int nyd, int nzd,
                   int nxhyzd, int nxyzhd);

void ckncfft3rv3z(float complex f[], int isign, int mixup[],
                  float complex sct[], int indx, int indy, int indz,
                  int nyi, int nyp, int nxhd, int nyd, int nzd,
                  int nxhyzd, int nxyzhd);

void ckncwfft3rvx(float complex f[], int isign, int mixup[],
                  float complex sct[], int indx, int indy, int indz,
                  int nxhd, int nyd, int nzd, int nxhyzd, int nxyzhd) ;

void ckncwfft3rv3(float complex f[], int isign, int mixup[],
                  float complex sct[], int indx, int indy, int indz,
                  int nxhd, int nyd, int nzd, int nxhyzd, int nxyzhd);
