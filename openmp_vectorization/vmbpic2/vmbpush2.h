/* header file for vmbpush2.c */

double ranorm();

void cdistr2h(float part[], float vtx, float vty, float vtz, float vdx,
              float vdy, float vdz, int npx, int npy, int idimp, int nop,
              int nx, int ny, int ipbc);

void cdblkp2l(float part[], int kpic[], int *nppmx, int idimp, int nop,
              int mx, int my, int mx1, int mxy1, int *irc);

void cppmovin2lt(float part[], float ppart[], int kpic[], int nppmx,
                 int idimp, int nop, int mx, int my, int mx1, int mxy1,
                 int *irc);

void cppmovin2ltp(float part[], float ppart[], int kpic[], int kp[],
                  int nppmx, int idimp, int nop, int mx, int my,
                  int mx1, int mxy1, int *irc);

void cppcheck2lt(float ppart[], int kpic[], int idimp, int nppmx, int nx,
                 int ny, int mx, int my, int mx1, int my1, 
                 int *irc);

void cgbppush23lt(float ppart[], float fxy[], float bxy[], int kpic[],
                  float qbm, float dt, float dtc, float *ek, int idimp,
                  int nppmx, int nx, int ny, int mx, int my, int nxv,
                  int nyv, int mx1, int mxy1, int ipbc);

void cgbppushf23lt(float ppart[], float fxy[], float bxy[], int kpic[],
                   int ncl[], int ihole[], float qbm, float dt,
                   float dtc, float *ek, int idimp, int nppmx, int nx,
                   int ny, int mx, int my, int nxv, int nyv, int mx1,
                   int mxy1, int ntmax, int *irc);

void cgrbppush23lt(float ppart[], float fxy[], float bxy[], int kpic[],
                   float qbm, float dt, float dtc, float ci, float *ek,
                   int idimp, int nppmx, int nx, int ny, int mx, int my,
                   int nxv, int nyv, int mx1, int mxy1, int ipbc);

void cgrbppushf23lt(float ppart[], float fxy[], float bxy[], int kpic[],
                    int ncl[], int ihole[], float qbm, float dt,
                    float dtc, float ci, float *ek, int idimp,
                    int nppmx, int nx, int ny, int mx, int my, int nxv,
                    int nyv, int mx1, int mxy1, int ntmax, int *irc);

void cvgbppush23lt(float ppart[], float fxy[], float bxy[], int kpic[],
                   float qbm, float dt, float dtc, float *ek, int idimp,
                   int nppmx, int nx, int ny, int mx, int my, int nxv,
                   int nyv, int mx1, int mxy1, int ipbc);

void cvgbppushf23lt(float ppart[], float fxy[], float bxy[], int kpic[],
                    int ncl[], int ihole[], float qbm, float dt,
                    float dtc, float *ek, int idimp, int nppmx, int nx,
                    int ny, int mx, int my, int nxv, int nyv, int mx1,
                    int mxy1, int ntmax, int *irc);

void cvgrbppush23lt(float ppart[], float fxy[], float bxy[], int kpic[],
                    float qbm, float dt, float dtc, float ci, float *ek,
                    int idimp, int nppmx, int nx, int ny, int mx,
                    int my, int nxv, int nyv, int mx1, int mxy1,
                    int ipbc);

void cvgrbppushf23lt(float ppart[], float fxy[], float bxy[],
                     int kpic[], int ncl[], int ihole[], float qbm,
                     float dt, float dtc, float ci, float *ek, 
                     int idimp, int nppmx, int nx, int ny, int mx,
                     int my, int nxv, int nyv, int mx1, int mxy1,
                     int ntmax, int *irc);

void cgppost2lt(float ppart[], float q[], int kpic[], float qm,
                int nppmx, int idimp, int mx, int my, int nxv, int nyv,
                int mx1, int mxy1);

void cvgppost2lt(float ppart[], float q[], int kpic[], float qm,
                 int nppmx, int idimp, int mx, int my, int nxv, int nyv,
                 int mx1, int mxy1);

void cviscan2(int *isdata, int *mb, int nths);

void cgjppost2lt(float ppart[], float cu[], int kpic[], float qm,
                 float dt, int nppmx, int idimp, int nx, int ny, int mx,
                 int my, int nxv, int nyv, int mx1, int mxy1,
                 int ipbc);

void cgjppostf2lt(float ppart[], float cu[], int kpic[], int ncl[],
                  int ihole[], float qm, float dt, int nppmx, int idimp,
                  int nx, int ny, int mx, int my, int nxv, int nyv,
                  int mx1, int mxy1, int ntmax, int *irc);

void cgrjppost2lt(float ppart[], float cu[], int kpic[], float qm,
                  float dt, float ci, int nppmx, int idimp, int nx,
                  int ny, int mx, int my, int nxv, int nyv, int mx1,
                  int mxy1, int ipbc);

void cgrjppostf2lt(float ppart[], float cu[], int kpic[], int ncl[],
                   int ihole[], float qm, float dt, float ci, int nppmx,
                   int idimp, int nx, int ny, int mx, int my, int nxv,
                   int nyv, int mx1, int mxy1, int ntmax, int *irc);

void cvgjppost2lt(float ppart[], float cu[], int kpic[], float qm,
                  float dt, int nppmx, int idimp, int nx, int ny,
                  int mx, int my, int nxv, int nyv, int mx1, int mxy1,
                  int ipbc);

void cvgjppostf2lt(float ppart[], float cu[], int kpic[], int ncl[],
                   int ihole[], float qm, float dt, int nppmx,
                   int idimp, int nx, int ny, int mx, int my, int nxv,
                   int nyv, int mx1, int mxy1, int ntmax, int *irc);

void cvgrjppost2lt(float ppart[], float cu[], int kpic[], float qm,
                   float dt, float ci, int nppmx, int idimp, int nx,
                   int ny, int mx, int my, int nxv, int nyv, int mx1,
                   int mxy1, int ipbc);

void cvgrjppostf2lt(float ppart[], float cu[], int kpic[], int ncl[],
                    int ihole[], float qm, float dt, float ci,
                    int nppmx, int idimp, int nx, int ny, int mx,
                    int my, int nxv, int nyv, int mx1, int mxy1,
                    int ntmax, int *irc);

void cpporder2lt(float ppart[], float ppbuff[], int kpic[], int ncl[],
                 int ihole[], int idimp, int nppmx, int nx, int ny,
                 int mx, int my, int mx1, int my1, int npbmx, int ntmax,
                 int *irc);

void cpporderf2lt(float ppart[], float ppbuff[], int kpic[], int ncl[],
                  int ihole[], int idimp, int nppmx, int mx1, int my1,
                  int npbmx, int ntmax, int *irc);

void cvpporder2lt(float ppart[], float ppbuff[], int kpic[], int ncl[],
                  int ihole[], int idimp, int nppmx, int nx, int ny,
                  int mx, int my, int mx1, int my1, int npbmx,
                  int ntmax, int *irc);

void cvpporderf2lt(float ppart[], float ppbuff[], int kpic[], int ncl[],
                   int ihole[], int idimp, int nppmx, int mx1, int my1,
                   int npbmx, int ntmax, int *irc);

void cbguard2l(float bxy[], int nx, int ny, int nxe, int nye);

void cacguard2l(float cu[], int nx, int ny, int nxe, int nye);

void caguard2l(float q[], int nx, int ny, int nxe, int nye);

void cvmpois23(float complex q[], float complex fxy[], int isign,
               float complex ffc[], float ax, float ay, float affp,
               float *we, int nx, int ny, int nxvh, int nyv, int nxhd,
               int nyhd);

void cmcuperp2(float complex cu[], int nx, int ny, int nxvh, int nyv);

void cvmibpois23(float complex cu[], float complex bxy[],
                 float complex ffc[], float ci, float *wm, int nx,
                 int ny, int nxvh, int nyv, int nxhd, int nyhd);

void cvmmaxwel2(float complex exy[], float complex bxy[],
                float complex cu[], float complex ffc[], float ci,
                float dt, float *wf, float *wm, int nx, int ny,
                int nxvh, int nyv, int nxhd, int nyhd);

void cvmemfield2(float complex fxy[], float complex exy[],
                 float complex ffc[], int isign, int nx, int ny,
                 int nxvh, int nyv, int nxhd, int nyhd);

void cwfft2rinit(int mixup[], float complex sct[], int indx, int indy, 
                 int nxhyd, int nxyhd);

void cfft2rvmxx(float complex f[], int isign, int mixup[],
                float complex sct[], int indx, int indy, int nyi,
                int nyp, int nxhd, int nyd, int nxhyd, int nxyhd);

void cfft2rmxy(float complex f[], int isign, int mixup[],
               float complex sct[], int indx, int indy, int nxi,
               int nxp, int nxhd, int nyd, int nxhyd, int nxyhd);

void cfft2rvm3x(float complex f[], int isign, int mixup[],
                float complex sct[], int indx, int indy, int nyi,
                int nyp, int nxhd, int nyd, int nxhyd, int nxyhd);

void cfft2rvm3y(float complex f[], int isign, int mixup[],
                float complex sct[], int indx, int indy, int nxi,
                int nxp, int nxhd, int nyd, int nxhyd, int nxyhd);

void cwfft2rvmx(float complex f[], int isign, int mixup[],
                float complex sct[], int indx, int indy, int nxhd,
                int nyd, int nxhyd, int nxyhd);

void cwfft2rvm3(float complex f[], int isign, int mixup[],
                float complex sct[], int indx, int indy, int nxhd,
                int nyd, int nxhyd, int nxyhd);

