/* header file for mpbpush2.c */

double ranorm();

void cpdicomp2l(float edges[], int *nyp, int *noff, int *nypmx,
                int *nypmn, int ny, int kstrt, int nvp, int idps);

void cpdistr2h(float part[], float edges[], int *npp, int nps,
               float vtx, float vty, float vtz, float vdx, float vdy,
               float vdz, int npx, int npy, int nx, int ny, int idimp,
               int npmax, int idps, int ipbc, int *ierr);

void cppdblkp2l(float part[], int kpic[], int npp, int noff, int *nppmx,
                int idimp, int npmax, int mx, int my, int mx1,
                int mxyp1, int *irc);

void cpppmovin2l(float part[], float ppart[], int kpic[], int npp,
                 int noff, int nppmx, int idimp, int npmax, int mx,
                 int my, int mx1, int mxyp1, int *irc);

void cpppcheck2l(float ppart[], int kpic[], int noff, int nyp,
                 int idimp, int nppmx, int nx, int mx, int my, int mx1,
                 int myp1, int *irc);

void cppgbppush23l(float ppart[], float fxy[], float bxy[], int kpic[],
                   int noff, int nyp, float qbm, float dt, float dtc,
                   float *ek, int idimp, int nppmx, int nx, int ny,
                   int mx, int my, int nxv, int nypmx, int mx1,
                   int mxyp1, int ipbc);

void cppgbppushf23l(float ppart[], float fxy[], float bxy[], int kpic[],
                    int ncl[], int ihole[], int noff, int nyp,
                    float qbm, float dt, float dtc, float *ek,
                    int idimp, int nppmx, int nx, int ny,
                    int mx, int my, int nxv, int nypmx, int mx1,
                    int mxyp1, int ntmax, int *irc);

void cppgrbppush23l(float ppart[], float fxy[], float bxy[], int kpic[],
                    int noff, int nyp, float qbm, float dt, float dtc,
                    float ci, float *ek, int idimp, int nppmx, int nx,
                    int ny, int mx, int my, int nxv, int nypmx, int mx1,
                    int mxyp1, int ipbc);

void cppgrbppushf23l(float ppart[], float fxy[], float bxy[],
                     int kpic[], int ncl[], int ihole[], int noff,
                     int nyp, float qbm, float dt, float dtc, float ci,
                     float *ek, int idimp, int nppmx, int nx, int ny,
                     int mx, int my, int nxv, int nypmx, int mx1,
                     int mxyp1, int ntmax, int *irc);

void cppgppost2l(float ppart[], float q[], int kpic[], int noff, 
                 float qm, int idimp, int nppmx, int mx, int my,
                 int nxv, int nypmx, int mx1, int mxyp1);

void cppgjppost2l(float ppart[], float cu[], int kpic[], int noff,
                  float qm, float dt, int nppmx, int idimp, int nx,
                  int ny, int mx, int my, int nxv, int nypmx, int mx1,
                  int mxyp1, int ipbc);

void cppgjppostf2l(float ppart[], float cu[], int kpic[], int ncl[],
                   int ihole[], int noff, int nyp, float qm, float dt,
                   int nppmx, int idimp, int nx, int ny, int mx, int my,
                   int nxv, int nypmx, int mx1, int mxyp1, int ntmax,
                   int *irc);

void cppgrjppost2l(float ppart[], float cu[], int kpic[], int noff,
                   float qm, float dt, float ci, int nppmx, int idimp,
                   int nx, int ny, int mx, int my, int nxv, int nypmx,
                   int mx1, int mxyp1, int ipbc);

void cppgrjppostf2l(float ppart[], float cu[], int kpic[], int ncl[],
                    int ihole[], int noff, int nyp, float qm, float dt,
                    float ci, int nppmx, int idimp, int nx, int ny,
                    int mx, int my, int nxv, int nypmx, int mx1,
                    int mxyp1, int ntmax, int *irc);

void cppporder2la(float ppart[], float ppbuff[], float sbufl[],
                  float sbufr[], int kpic[], int ncl[], int ihole[],
                  int ncll[], int nclr[], int noff, int nyp, int idimp,
                  int nppmx, int nx, int ny, int mx, int my, int mx1,
                  int myp1, int npbmx, int ntmax, int nbmax, int *irc);

void cppporderf2la(float ppart[], float ppbuff[], float sbufl[],
                   float sbufr[], int ncl[], int ihole[], int ncll[],
                   int nclr[], int idimp, int nppmx, int mx1, int myp1,
                   int npbmx, int ntmax, int nbmax, int *irc);

void cppporder2lb(float ppart[], float ppbuff[], float rbufl[],
                  float rbufr[], int kpic[], int ncl[], int ihole[],
                  int mcll[], int mclr[], int idimp, int nppmx, int mx1,
                  int myp1, int npbmx, int ntmax, int nbmax, int *irc);

void cppcguard2xl(float fxy[], int myp, int nx, int ndim, int nxe,
                  int nypmx);

void cppaguard2xl(float q[], int myp, int nx, int nxe, int nypmx);

void cppacguard2xl(float cu[], int myp, int nx, int ndim, int nxe,
                   int nypmx);

void cmppois23(float complex q[], float complex fxy[], int isign,
               float complex ffc[], float ax, float ay, float affp,
               float *we, int nx, int ny, int kstrt, int nyv, int kxp,
               int nyhd);

void cmppcuperp2(float complex cu[], int nx, int ny, int kstrt, int nyv,
                 int kxp);

void cmippbpoisp23(float complex cu[], float complex bxy[],
                   float complex ffc[], float ci, float *wm, int nx,
                   int ny, int kstrt, int nyv, int kxp, int nyhd);

void cmppmaxwel2(float complex exy[], float complex bxy[],
                 float complex cu[], float complex ffc[], float affp,
                 float ci, float dt, float *wf, float *wm, int nx,
                 int ny, int kstrt, int nyv, int kxp, int nyhd);

void cmppemfield2(float complex fxy[], float complex exy[],
                  float complex ffc[], int isign, int nx, int ny,
                  int kstrt, int nyv, int kxp, int nyhd);

void cwpfft2rinit(int mixup[], float complex sct[], int indx, int indy,
                  int nxhyd, int nxyhd);

void cppfft2rmxx(float complex f[], int isign, int mixup[],
                 float complex sct[], int indx, int indy, int kstrt,
                 int kypi, int kypp, int nxvh, int kypd, int nxhyd,
                 int nxyhd);

void cppfft2rmxy(float complex g[], int isign, int mixup[],
                 float complex sct[], int indx, int indy, int kstrt,
                 int kxpi, int kxpp, int nyv, int kxp, int nxhyd,
                 int nxyhd);

void cppfft2rm3xx(float complex f[], int isign, int mixup[],
                  float complex sct[], int indx, int indy, int kstrt,
                  int kypi, int kypp, int nxvh, int kypd, int nxhyd,
                  int nxyhd);

void cppfft2rm3xy(float complex g[], int isign, int mixup[],
                  float complex sct[], int indx, int indy, int kstrt,
                  int kxpi, int kxpp, int nyv, int kxp, int nxhyd,
                  int nxyhd);

void cwppfft2rm(float complex f[], float complex g[],
                float complex bs[], float complex br[], int isign,
                int ntpose, int mixup[], float complex sct[],
                float *ttp, int indx, int indy, int kstrt, int nvp,
                int nxvh, int nyv, int kxp, int kyp, int kypd,
                int nxhyd, int nxyhd);

void cwppfft2rm3(float complex f[], float complex g[],
                 float complex bs[], float complex br[], int isign,
                 int ntpose, int mixup[], float complex sct[],
                 float *ttp, int indx, int indy, int kstrt, int nvp,
                 int nxvh, int nyv, int kxp, int kyp, int kypd,
                 int nxhyd, int nxyhd);
