/* header file for vmpush3.c */

double ranorm();

/*--------------------------------------------------------------------*/
void cdistr3(float part[], float vtx, float vty, float vtz, float vdx,
             float vdy, float vdz, int npx, int npy, int npz, int idimp, 
             int nop, int nx, int ny, int nz, int ipbc);

/*--------------------------------------------------------------------*/
void cdblkp3l(float part[], int kpic[], int *nppmx, int idimp, int nop,
              int mx, int my, int mz, int mx1, int my1, int mxyz1,
              int *irc);

/*--------------------------------------------------------------------*/
void cppmovin3lt(float part[], float ppart[], int kpic[], int nppmx,
                 int idimp, int nop, int mx, int my, int mz, int mx1,
                 int my1, int mxyz1, int *irc);

/*--------------------------------------------------------------------*/
void cppmovin3ltp(float part[], float ppart[], int kpic[], int kp[],
                  int nppmx, int idimp, int nop, int mx, int my, int mz,
                  int mx1, int my1, int mxyz1, int *irc);

/*--------------------------------------------------------------------*/
void cppcheck3lt(float ppart[], int kpic[], int idimp, int nppmx,
                 int nx, int ny, int nz, int mx, int my, int mz,
                 int mx1, int my1, int mz1, int *irc);

/*--------------------------------------------------------------------*/
void cgppush3lt(float ppart[], float fxyz[], int kpic[], float qbm,
                float dt, float *ek, int idimp, int nppmx, int nx,
                int ny, int nz, int mx, int my, int mz, int nxv,
                int nyv, int nzv, int mx1, int my1, int mxyz1,
                int ipbc);

/*--------------------------------------------------------------------*/
void cgppushf3lt(float ppart[], float fxyz[], int kpic[], int ncl[],
                 int ihole[], float qbm, float dt, float *ek, int idimp,  
                 int nppmx, int nx, int ny, int nz, int mx, int my,
                 int mz, int nxv, int nyv, int nzv, int mx1, int my1,
                 int mxyz1, int ntmax, int *irc);

/*--------------------------------------------------------------------*/
void cvgppush3lt(float ppart[], float fxyz[], int kpic[], float qbm,
                 float dt, float *ek, int idimp, int nppmx, int nx,
                 int ny, int nz, int mx, int my, int mz, int nxv,
                 int nyv, int nzv, int mx1, int my1, int mxyz1,
                 int ipbc);

/*--------------------------------------------------------------------*/
void cvgppushf3lt(float ppart[], float fxyz[], int kpic[], int ncl[],
                  int ihole[], float qbm, float dt, float *ek,
                  int idimp, int nppmx, int nx, int ny, int nz, int mx,
                  int my, int mz, int nxv, int nyv, int nzv, int mx1,
                  int my1, int mxyz1, int ntmax, int *irc);

/*--------------------------------------------------------------------*/
void cv2gppush3lt(float ppart[], float fxyz[], int kpic[], float qbm,
                  float dt, float *ek, int idimp, int nppmx, int nx,
                  int ny, int nz, int mx, int my, int mz, int nxv,
                  int nyv, int nzv, int mx1, int my1, int mxyz1,
                  int ipbc);

/*--------------------------------------------------------------------*/
void cv2gppushf3lt(float ppart[], float fxyz[], int kpic[], int ncl[],
                   int ihole[], float qbm, float dt, float *ek,
                   int idimp, int nppmx, int nx, int ny, int nz, int mx,
                   int my, int mz, int nxv, int nyv, int nzv, int mx1,
                   int my1, int mxyz1, int ntmax, int *irc);

/*--------------------------------------------------------------------*/
void cgppost3lt(float ppart[], float q[], int kpic[], float qm,
                int nppmx, int idimp, int mx, int my, int mz, int nxv,
                int nyv, int nzv, int mx1, int my1, int mxyz1);

/*--------------------------------------------------------------------*/
void cvgppost3lt(float ppart[], float q[], int kpic[], float qm,
                 int nppmx, int idimp, int mx, int my, int mz, int nxv,
                 int nyv, int nzv, int mx1, int my1, int mxyz1);

/*--------------------------------------------------------------------*/
void cviscan2(int *isdata, int *mb, int nths);

/*--------------------------------------------------------------------*/
void cpporder3lt(float ppart[], float ppbuff[], int kpic[], int ncl[],
                 int ihole[], int idimp, int nppmx, int nx, int ny,
                 int nz, int mx, int my, int mz, int mx1, int my1,
                 int mz1, int npbmx, int ntmax, int *irc);

/*--------------------------------------------------------------------*/
void cpporderf3lt(float ppart[], float ppbuff[], int kpic[], int ncl[],
                  int ihole[], int idimp, int nppmx, int mx1, int my1,
                  int mz1, int npbmx, int ntmax, int *irc);

/*--------------------------------------------------------------------*/
void cvpporder3lt(float ppart[], float ppbuff[], int kpic[], int ncl[],
                  int ihole[], int idimp, int nppmx, int nx, int ny,
                  int nz, int mx, int my, int mz, int mx1, int my1,
                  int mz1, int npbmx, int ntmax, int *irc);

/*--------------------------------------------------------------------*/
void cvpporderf3lt(float ppart[], float ppbuff[], int kpic[], int ncl[],
                   int ihole[], int idimp, int nppmx, int mx1, int my1,
                   int mz1, int npbmx, int ntmax, int *irc);

/*--------------------------------------------------------------------*/
void cv2pporderf3lt(float ppart[], float ppbuff[], int kpic[],
                    int ncl[], int ihole[], int idimp, int nppmx,
                    int mx1, int my1, int mz1, int npbmx, int ntmax,
                    int *irc);

/*--------------------------------------------------------------------*/
void ccguard3l(float fxyz[], int nx, int ny, int nz, int nxe, int nye,
               int nze);

/*--------------------------------------------------------------------*/
void caguard3l(float q[], int nx, int ny, int nz, int nxe, int nye,
               int nze);

/*--------------------------------------------------------------------*/
void cvmpois33(float complex q[], float complex fxyz[], int isign,
               float complex ffc[], float ax, float ay, float az,
               float affp, float *we, int nx, int ny, int nz, int nxvh,
               int nyv, int nzv, int nxhd, int nyhd, int nzhd);

/*--------------------------------------------------------------------*/
void cwfft3rinit(int mixup[], float complex sct[], int indx, int indy,
                 int indz, int nxhyzd, int nxyzhd);

/*--------------------------------------------------------------------*/
void cfft3rvmxy(float complex f[], int isign, int mixup[],
                float complex sct[], int indx, int indy, int indz,
                int nzi, int nzp, int nxhd, int nyd, int nzd,
                int nxhyzd, int nxyzhd);

/*--------------------------------------------------------------------*/
void cfft3rvmxz(float complex f[], int isign, int mixup[],
                float complex sct[], int indx, int indy, int indz,
                int nyi, int nyp, int nxhd, int nyd, int nzd,
                int nxhyzd, int nxyzhd);

/*--------------------------------------------------------------------*/
void cfft3rvm3xy(float complex f[], int isign, int mixup[],
                 float complex sct[], int indx, int indy, int indz,
                 int nzi, int nzp, int nxhd, int nyd, int nzd,
                 int nxhyzd, int nxyzhd);

/*--------------------------------------------------------------------*/
void cfft3rvm3z(float complex f[], int isign, int mixup[],
                float complex sct[], int indx, int indy, int indz,
                int nyi, int nyp, int nxhd, int nyd, int nzd,
                int nxhyzd, int nxyzhd);

/*--------------------------------------------------------------------*/
void cwfft3rvmx(float complex f[], int isign, int mixup[],
                float complex sct[], int indx, int indy, int indz,
                int nxhd, int nyd, int nzd, int nxhyzd, int nxyzhd);

/*--------------------------------------------------------------------*/
void cwfft3rvm3(float complex f[], int isign, int mixup[],
                float complex sct[], int indx, int indy, int indz,
                int nxhd, int nyd, int nzd, int nxhyzd, int nxyzhd);

/*--------------------------------------------------------------------*/
void cset_szero3(float q[], int mx, int my, int mz, int nxv, int nyv,
                 int nzv, int mx1, int my1, int mxyz1);

/*--------------------------------------------------------------------*/
void cset_vzero3(float cu[], int mx, int my, int mz, int ndim, int nxv,
                 int nyv, int nzv, int mx1, int my1, int mxyz1);

/*--------------------------------------------------------------------*/
void cset_cvzero3(float complex exyz[], int nx, int ny, int nz,
                  int ndim, int nxvh, int nyv, int nzv);
