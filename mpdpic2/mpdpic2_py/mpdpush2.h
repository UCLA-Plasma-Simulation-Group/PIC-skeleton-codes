/* header file for mpdpush2.c */

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

void cppgppost2l(float ppart[], float q[], int kpic[], int noff, 
                 float qm, int idimp, int nppmx, int mx, int my,
                 int nxv, int nypmx, int mx1, int mxyp1);

void cppgjppost2l(float ppart[], float cu[], int kpic[], int noff,
                  float qm, float dt, int nppmx, int idimp, int nx,
                  int ny, int mx, int my, int nxv, int nypmx, int mx1,
                  int mxyp1, int ipbc);

void cppgmjppost2l(float ppart[], float amu[], int kpic[], int noff,
                   float qm, int nppmx, int idimp, int mx, int my,
                   int nxv, int nypmx, int mx1, int mxyp1);

void cppgdjppost2l(float ppart[], float fxy[], float bxy[], float dcu[],
                   float amu[], int kpic[], int noff, int nyp, float qm,
                   float qbm, float dt, int idimp, int nppmx, int nx,
                   int mx, int my, int nxv, int nypmx, int mx1,
                   int mxyp1);

void cppgdcjppost2l(float ppart[], float fxy[], float bxy[], float cu[],
                    float dcu[], float amu[], int kpic[], int noff,
                    int nyp, float qm, float qbm, float dt, int idimp,
                    int nppmx, int nx, int mx, int my, int nxv,
                    int nypmx, int mx1, int mxyp1);

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

void cppascfguard2l(float dcu[], float cus[], int nyp, float q2m0,
                    int nx, int nxe, int nypmx);

void cppfwpminmx2(float qe[], int nyp, float qbme, float *wpmax,
                  float *wpmin, int nx, int nxe, int nypmx);

void cmppois23(float complex q[], float complex fxy[], int isign,
               float complex ffc[], float ax, float ay, float affp,
               float *we, int nx, int ny, int kstrt, int nyv, int kxp,
               int nyhd);

void cmppcuperp2(float complex cu[], int nx, int ny, int kstrt, int nyv,
                 int kxp);

void cmppbbpoisp23(float complex cu[], float complex bxy[],
                   float complex ffc[], float ci, float *wm, int nx,
                   int ny, int kstrt, int nyv, int kxp, int nyhd);

void cppbaddext2(float bxy[], int nyp, float omx, float omy, float omz,
                 int nx, int nxe, int nypmx);

void cmppdcuperp23(float complex dcu[], float complex amu[], int nx,
                   int ny, int kstrt, int nyv, int kxp);

void cmppadcuperp23(float complex dcu[], float complex amu[], int nx,
                    int ny, int kstrt, int nyv, int kxp);

void cmppepoisp23(float complex dcu[], float complex exy[], int isign,
                  float complex ffe[], float ax, float ay, float affp,
                  float wp0, float ci, float *wf, int nx, int ny,
                  int kstrt, int nyv, int kxp, int nyhd);

void cppaddvrfield2(float a[], float b[], float c[], int ndim, int nxe,
                    int nypmx);

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

void cppfft2rmnxx(float complex f[], float complex ss[], int isign,
                  int mixup[], float complex sct[], int indx, int indy,
                  int kstrt, int kypi, int kypp, int nxvh, int kypd,
                  int ndim, int nxhyd, int nxyhd);

void cppfft2rmnxy(float complex g[], int isign, int mixup[],
                  float complex sct[], int indx, int indy, int kstrt,
                  int kxpi, int kxpp, int nyv, int kxp, int ndim,
                  int nxhyd, int nxyhd);

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

void cwppfft2rmn(float complex f[], float complex g[],
                 float complex bs[], float complex br[],
                 float complex ss[], int isign, int ntpose, int mixup[],
                 float complex sct[], float *ttp, int indx, int indy,
                 int kstrt, int nvp, int nxvh, int nyv, int kxp,
                 int kyp, int kypd, int ndim, int nxhyd, int nxyhd);

void cmppswapc2n(float f[], float s[], int isign, int nxh, int kypi,
                 int kypt, int nxvh, int kypd, int ndim);
