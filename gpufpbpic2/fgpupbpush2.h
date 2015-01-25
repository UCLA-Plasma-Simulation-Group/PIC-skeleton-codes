/* C header file for gpupbpush2.cuf */

void fgpuppgbppush23l(float ppart[], float fxy[], float bxy[],
                      int kpic[], int noff, int nyp, float qbm,
                      float dt, float dtc, float *ek, int idimp,
                      int nppmx, int nx, int ny, int mx, int my, 
                      int nxv, int nypmx, int mx1, int mxyp1, int ipbc);

void fgpuppgrbppush23l(float ppart[], float fxy[], float bxy[],
                       int kpic[], int noff, int nyp, float qbm,
                       float dt, float dtc, float ci, float *ek,
                       int idimp, int nppmx, int nx, int ny, int mx,
                       int my, int nxv, int nypmx, int mx1, int mxyp1,
                       int ipbc);

void fgpu2ppgppost2l(float ppart[], float q[], int kpic[], int noff,
                     float qm, int idimp, int nppmx, int mx, int my,
                     int nxv, int nypmx, int mx1, int mxyp1);

void fgpu2ppjppost2l(float ppart[], float cu[], int kpic[], int noff,
                     float qm, float dt, int nppmx, int idimp, int nx,
                     int ny, int mx, int my, int nxv, int nypmx,
                     int mx1, int mxyp1, int ipbc);

void fgpu2pprjppost2l(float ppart[], float cu[], int kpic[], int noff,
                      float qm, float dt, float ci, int nppmx,
                      int idimp, int nx, int ny, int mx, int my,
                      int nxv, int nypmx, int mx1, int mxyp1, int ipbc);

void fgpuppcaguard2xl(float complex qc[], float scs[], float q[],
                      int nyp, int nx, int nxe, int nypmx, int nxvh,
                      int kypd);

void fgpuppcaguard2yl(float complex fc[], float scr[], int nx, int nxvh,
                      int kypd);

void fgpuppcacguard2xl(float complex cuc[], float scs[], float cu[],
                       int nyp, int nx, int nxe, int nypmx, int nxvh,
                       int kypd);

void fgpuppcacguard2yl(float complex fvc[], float scr[], int nx,
                       int nxvh, int kypd);

void fgpuppcbguard2xl(float complex fxyc[], float scs[], float fxy[],
                      int nyp, int nx, int nxe, int nypmx, int nxvh,
                      int kypd);

void fgpuppcbguard2yl(float fxy[], float scr[], int nyp, int nx,
                      int nxe, int nxvh, int nypmx);

void fgpupppord2la(float ppart[], float ppbuff[], float sbufl[],
                   float sbufr[], int kpic[], int ncl[], int ihole[], 
                   int ncll[], int nclr[], int noff, int nyp, int idimp,
                   int nppmx, int nx, int ny, int mx, int my, int mx1,
                   int myp1, int npbmx, int ntmax, int nbmax, int *irc);

void fgpupppord2lb(float ppart[], float ppbuff[], float rbufl[],
                   float rbufr[], int kpic[], int ncl[], int ihole[],
                   int mcll[], int mclr[], int idimp, int nppmx, int mx1,
                   int myp1, int npbmx, int ntmax, int nbmax, int *irc);

void fgpuppois23t(float complex qt[], float complex fxyt[],
                  float complex ffct[], float *we, int nx, int ny,
                  int kstrt, int nyv, int kxp1, int nyhd);

void fgpuppcuperp2t(float complex cut[], int nx, int ny, int kstrt,
                    int nyv, int kxp1);

void fgpuippbpoisp23t(float complex cut[], float complex bxyt[],
                      float complex ffct[], float ci, float *wm, int nx,
                      int ny, int kstrt, int nyv, int kxp1, int nyhd);

void fgpuppmaxwel2t(float complex exyt[], float complex bxyt[],
                    float complex cut[], float complex ffct[],
                    float affp, float ci, float dt, float *wf,
                    float *wm, int nx, int ny, int kstrt, int nyv,
                    int kxp1, int nyhd);

void fgpuppemfield2t(float complex fxyt[], float complex exyt[],
                     float complex ffct[], int isign, int nx, int ny,
                     int kstrt, int nyv, int kxp1, int nyhd);

void fgpuwppfft2rcsx(float complex f[], float complex bsm[], int isign,
                     int mixup[], float complex sct[], int indx,
                     int indy, int kstrt, int nvp, int kxp1, int kyp,
                     int nxhd, int kypd, int nxhyd, int nxyhd);

void fgpuwppfft2rcsy(float complex g[], float complex brm[], int isign,
                     int mixup[], float complex sct[], int indx,
                     int indy, int kstrt, int nvp, int kxp1, int kyp,
                     int nyd, int nxhyd, int nxyhd);

void fgpuwppfft2rcsxn(float complex fn[], float complex bsm[], int isign,
                      int mixup[], float complex sct[], int indx,
                      int indy, int ndim, int kstrt, int nvp, int kxp1,
                      int kyp, int nxhd, int kypd, int nxhyd, int nxyhd);

void fgpuwppfft2rcsyn(float complex gn[], float complex brm[], int isign,
                      int mixup[], float complex sct[], int indx,
                      int indy, int ndim, int kstrt, int nvp, int kxp1,
                      int kyp, int nyd, int nxhyd, int nxyhd);

void fgpuppltpose(float complex f[], float complex g[], int nx, int ny,
                  int kxp, int kyp, int kstrt, int nxv, int nyv);

void fgpuppltposen(float complex fn[], float complex gn[], int nx,
                   int ny, int kxp, int kyp, int kstrt, int ndim,
                   int nxv, int nyv);

void fgpusum2(float a[], float *sa, int nx);

