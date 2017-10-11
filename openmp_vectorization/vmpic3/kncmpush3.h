/* header file for kncmpush3.c */

void ckncgppush3lt(float ppart[], float fxyz[], int kpic[], float qbm,
                   float dt, float *ek, int idimp, int nppmx, int nx,
                   int ny, int nz, int mx, int my, int mz, int nxv,
                   int nyv, int nzv, int mx1, int my1, int mxyz1,
                   int ipbc);

void ckncgppushf3lt(float ppart[], float fxyz[], int kpic[], int ncl[],
                    int ihole[], float qbm, float dt, float *ek,
                    int idimp, int nppmx, int nx, int ny, int nz,
                    int mx, int my, int mz, int nxv, int nyv, int nzv,
                    int mx1, int my1, int mxyz1, int ntmax, int *irc);

void ckncgppost3lt(float ppart[], float q[], int kpic[], float qm,
                   int nppmx, int idimp, int mx, int my, int mz,
                   int nxv, int nyv, int nzv, int mx1, int my1,
                   int mxyz1);

void cknc2gppost3lt(float ppart[], float q[], int kpic[], float qm,
                    int nppmx, int idimp, int mx, int my, int mz,
                    int nxv, int nyv, int nzv, int mx1, int my1,
                    int mxyz1);

void ckncpporder3lt(float ppart[], float ppbuff[], int kpic[],
                    int ncl[], int ihole[], int idimp, int nppmx, 
                    int nx, int ny, int nz, int mx, int my, int mz,
                    int mx1, int my1, int mz1, int npbmx, int ntmax,
                    int *irc);

void ckncpporderf3lt(float ppart[], float ppbuff[], int kpic[],
                     int ncl[], int ihole[], int idimp, int nppmx,
                     int mx1, int my1, int mz1, int npbmx, int ntmax,
                     int *irc);

void cknccguard3l(float fxyz[], int nx, int ny, int nz, int nxe,
                  int nye, int nze);

void ckncaguard3l(float q[], int nx, int ny, int nz, int nxe, int nye,
                  int nze);

void ckncmpois33(float complex q[], float complex fxyz[], int isign,
                 float complex ffc[], float ax, float ay, float az,
                 float affp, float *we, int nx, int ny, int nz, int nxvh,
                 int nyv, int nzv, int nxhd, int nyhd, int nzhd);

void ckncfft3rmxy(float complex f[], int isign, int mixup[],
                  float complex sct[], int indx, int indy, int indz,
                  int nzi, int nzp, int nxhd, int nyd, int nzd,
                  int nxhyzd, int nxyzhd);

void ckncfft3rmz(float complex f[], int isign, int mixup[],
                 float complex sct[], int indx, int indy, int indz,
                 int nyi, int nyp, int nxhd, int nyd, int nzd,
                 int nxhyzd, int nxyzhd);

void ckncfft3rm3xy(float complex f[], int isign, int mixup[],
                   float complex sct[], int indx, int indy, int indz,
                   int nzi, int nzp, int nxhd, int nyd, int nzd,
                   int nxhyzd, int nxyzhd);

void ckncfft3rm3z(float complex f[], int isign, int mixup[],
                  float complex sct[], int indx, int indy, int indz,
                  int nyi, int nyp, int nxhd, int nyd, int nzd,
                  int nxhyzd, int nxyzhd);

void ckncwfft3rmx(float complex f[], int isign, int mixup[],
                  float complex sct[], int indx, int indy, int indz,
                  int nxhd, int nyd, int nzd, int nxhyzd, int nxyzhd) ;

void ckncwfft3rm3(float complex f[], int isign, int mixup[],
                  float complex sct[], int indx, int indy, int indz,
                  int nxhd, int nyd, int nzd, int nxhyzd, int nxyzhd);
