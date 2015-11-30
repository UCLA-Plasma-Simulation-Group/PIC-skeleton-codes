/* header file for mppush3.c */

double ranorm();

void cpdicomp32l(float edges[], int nyzp[], int noff[], int *nypmx,
                 int *nzpmx, int *nypmn, int *nzpmn, int ny, int nz,
                 int kstrt, int nvpy, int nvpz, int idps, int idds);

void cfcomp32(int nvp, int nx, int ny, int nz, int *nvpy, int *nvpz,
              int *ierr);

void cpdistr32(float part[], float edges[], int *npp, int nps,
               float vtx, float vty, float vtz, float vdx, float vdy,
               float vdz, int npx, int npy, int npz, int nx, int ny,
               int nz, int idimp, int npmax, int idps, int ipbc,
               int *ierr);

void cppdblkp3l(float part[], int kpic[], int npp, int noff[],
                int *nppmx, int idimp, int npmax, int mx, int my,
                int mz, int mx1, int myp1, int mxyzp1, int idds,
                int *irc);

void cpppmovin3l(float part[], float ppart[], int kpic[], int npp,
                 int noff[], int nppmx, int idimp, int npmax, int mx,
                 int my, int mz, int mx1, int myp1, int mxyzp1,
                 int idds, int *irc);

void cpppcheck3l(float ppart[], int kpic[], int noff[], int nyzp[],
                 int idimp, int nppmx, int nx, int mx, int my, int mz,
                 int mx1, int myp1, int mzp1, int idds, int *irc);

void cppgppush32l(float ppart[], float fxyz[], int kpic[], int noff[],
                  int nyzp[], float qbm, float dt, float *ek, int idimp, 
                  int nppmx, int nx, int ny, int nz, int mx, int my,
                  int mz, int nxv, int nypmx, int nzpmx, int mx1,
                  int myp1, int mxyzp1, int idds, int ipbc);

void cppgppushf32l(float ppart[], float fxyz[], int kpic[], int ncl[],
                   int ihole[], int noff[], int nyzp[], float qbm,
                   float dt, float *ek, int idimp, int nppmx, int nx,
                   int ny, int nz, int mx, int my, int mz, int nxv,
                   int nypmx, int nzpmx, int mx1, int myp1, int mxyzp1,
                   int ntmax, int idds, int *irc);

void cppgppost32l(float ppart[], float q[], int kpic[], int noff[],
                  float qm, int nppmx, int idimp, int mx, int my,
                  int mz, int nxv, int nypmx, int nzpmx, int mx1,
                  int myp1, int mxyzp1, int idds);

void cppporder32la(float ppart[], float ppbuff[], float sbufl[],
                   float sbufr[], int kpic[], int ncl[], int ihole[],
                   int ncll[], int nclr[], int noff[], int nyzp[],
                   int idimp, int nppmx, int nx, int ny, int nz, int mx,
                   int my, int mz, int mx1, int myp1, int mzp1,
                   int mxzyp1, int npbmx, int ntmax, int nbmax,
                   int idds, int *irc);

void cppporderf32la(float ppart[], float ppbuff[], float sbufl[],
                    float sbufr[], int ncl[], int ihole[], int ncll[],
                    int nclr[], int idimp, int nppmx, int mx1, int myp1,
                    int mzp1, int mxzyp1, int npbmx, int ntmax,
                    int nbmax, int *irc);

void cppporder32lb(float ppart[], float ppbuff[], float rbufl[],
                   float rbufr[], int kpic[], int ncl[], int ihole[],
                   int mcll[], int mclr[], int mcls[], int idimp,
                   int nppmx, int mx1, int myp1, int mzp1, int mxzyp1,
                   int npbmx, int ntmax, int nbmax, int *irc);

void cppcguard32xl(float fxyz[], int nyzp[], int nx, int ndim, int nxe,
                   int nypmx, int nzpmx, int idds);

void cppaguard32xl(float q[], int nyzp[], int nx, int nxe, int nypmx,
                   int nzpmx, int idds);

void cmppois332(float complex q[], float complex fxyz[], int isign,
                float complex ffc[], float ax, float ay, float az,
                float affp, float *we, int nx, int ny, int nz,
                int kstrt, int nvpy, int nvpz, int nzv, int kxyp,
                int kyzp, int nzhd);

void cwpfft32rinit(int mixup[], float complex sct[], int indx, int indy,
                   int indz, int nxhyzd, int nxyzhd);

void cppfft32rmxx(float complex f[], int isign, int mixup[],
                  float complex sct[], int indx, int indy, int indz,
                  int kstrt, int nvp, int kypi, int kypp, int nxvh,
                  int kzpp, int kypd, int kzpd, int nxhyzd, int nxyzhd);

void cppfft32rmxy(float complex g[], int isign, int mixup[],
                  float complex sct[], int indx, int indy, int indz,
                  int kstrt, int nvpy, int nvpz, int kxypi, int kxypp,
                  int nyv, int kzpp, int kxypd, int kzpd, int nxhyzd,
                  int nxyzhd);

void cppfft32rmxz(float complex h[], int isign, int mixup[],
                  float complex sct[], int indx, int indy, int indz,
                  int kstrt, int nvpy, int nvpz, int kxypi, int kxypp,
                  int nzv, int kyzp, int kxypd, int kyzpd, int nxhyzd,
                  int nxyzhd);

void cppfft32rm3xx(float complex f[], int isign, int mixup[],
                   float complex sct[], int indx, int indy, int indz,
                   int kstrt, int nvp, int kypi, int kypp, int nxvh,
                   int kzpp, int kypd, int kzpd, int nxhyzd,
                   int nxyzhd);

void cppfft32rm3xy(float complex g[], int isign, int mixup[],
                   float complex sct[], int indx, int indy, int indz,
                   int kstrt, int nvpy, int nvpz, int kxypi, int kxypp,
                   int nyv, int kzpp, int kxypd, int kzpd, int nxhyzd,
                   int nxyzhd);

void cppfft32rm3xz(float complex h[], int isign, int mixup[],
                   float complex sct[], int indx, int indy, int indz,
                   int kstrt, int nvpy, int nvpz, int kxypi, int kxypp,
                   int nzv, int kyzp, int kxypd, int kyzpd, int nxhyzd,
                   int nxyzhd);

void cwppfft32rm(float complex f[], float complex g[],
                 float complex h[], float complex bs[], 
                 float complex br[], int isign, int ntpose, int mixup[],
                 float complex sct[], float *ttp, int indx, int indy,
                 int indz, int kstrt, int nvpy, int nvpz, int nxvh,
                 int nyv, int nzv, int kxyp, int kyp, int kyzp, int kzp,
                 int kxypd, int kypd, int kyzpd, int kzpd, int kzyp,
                 int nxhyzd, int nxyzhd);

void cwppfft32rm3(float complex f[], float complex g[],
                  float complex h[], float complex bs[],
                  float complex br[], int isign, int ntpose,
                  int mixup[], float complex sct[], float *ttp, 
                  int indx, int indy, int indz, int kstrt, int nvpy,
                  int nvpz, int nxvh, int nyv, int nzv, int kxyp, 
                  int kyp, int kyzp, int kzp, int kxypd, int kypd,
                  int kyzpd, int kzpd, int kzyp, int nxhyzd, 
                  int nxyzhd);
