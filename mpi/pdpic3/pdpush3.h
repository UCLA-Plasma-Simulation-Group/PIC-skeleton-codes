/* header file for pdpush3.c */

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

void cppgbpush32l(float part[], float fxyz[], float bxyz[],
                  float edges[], int npp, int noff[], int ihole[],
                  float qbm, float dt, float dtc, float *ek, int nx,
                  int ny, int nz, int idimp, int npmax, int nxv,
                  int nypmx, int nzpmx, int idps, int idds, int ntmax,
                  int ipbc);

void cppgpost32l(float part[], float q[], int npp, int noff[], float qm,
                 int idimp, int npmax, int nxv, int nypmx, int nzpmx,
                 int idds);

void cppgjpost32l(float part[], float cu[], float edges[], int npp,
                  int noff[], int ihole[], float qm, float dt, int nx,
                  int ny, int nz, int idimp, int npmax, int nxv,
                  int nypmx, int nzpmx, int idps, int idds, int ntmax,
                  int ipbc);

void cppgmjpost32l(float part[], float amu[], int npp, int noff[],
                   float qm, int idimp, int npmax, int nxv, int nypmx,
                   int nzpmx, int idds);

void cppgdjpost32l(float part[], float fxyz[], float bxyz[], int npp,
                   int noff[], float dcu[], float amu[], float qm, 
                   float qbm, float dt, int idimp, int npmax, int nxv,
                   int nypmx, int nzpmx, int idds);

void cppgdcjpost32l(float part[], float fxyz[], float bxyz[], int npp,
                    int noff[], float cu[], float dcu[], float amu[],
                    float qm, float qbm, float dt, int idimp, int npmax,
                    int nxv, int nypmx, int nzpmx, int idds);

void cppdsortp32yzl(float parta[], float partb[], int npic[], int npp,
                    int noff[], int nyzp[], int idimp, int npmax,
                    int nyzpm1, int idds);

void cppcguard32xl(float fxyz[], int nyzp[], int nx, int ndim, int nxe,
                   int nypmx, int nzpmx, int idds);

void cppaguard32xl(float q[], int nyzp[], int nx, int nxe, int nypmx,
                   int nzpmx, int idds);

void cppacguard32xl(float cu[], int nyzp[], int nx, int ndim, int nxe,
                    int nypmx, int nzpmx, int idds);

void cppascfguard32l(float dcu[], float cus[], int nyzp[], float q2m0,
                     int nx, int nxe, int nypmx, int nzpmx, int idds);

void cppfwpminmx32(float qe[], int nyzp[], float qbme, float *wpmax,
                   float *wpmin, int nx, int nxe, int nypmx, int nzpmx,
                   int idds);

void cppois332(float complex q[], float complex fxyz[], int isign,
               float complex ffc[], float ax, float ay, float az,
               float affp, float *we, int nx, int ny, int nz, int kstrt,
               int nvpy, int nvpz, int nzv, int kxyp, int kyzp,
               int nzhd);

void cppcuperp32(float complex cu[], int nx, int ny, int nz, int kstrt,
                 int nvpy, int nvpz, int nzv, int kxyp, int kyzp);

void cppbbpoisp332(float complex cu[], float complex bxyz[],
                   float complex ffc[], float ci, float *wm, int nx,
                   int ny, int nz, int kstrt, int nvpy, int nvpz,
                   int nzv, int kxyp, int kyzp, int nzhd);

void cppbaddext32(float bxyz[], int nyzp[], float omx, float omy,
                  float omz, int nx, int nxe, int nypmx, int nzpmx,
                  int idds);

void cppdcuperp32(float complex dcu[], float complex amu[], int nx,
                  int ny, int nz, int kstrt, int nvpy, int nvpz,
                  int nzv, int kxyp, int kyzp);

void cppadcuperp32(float complex dcu[], float complex amu[], int nx,
                   int ny, int nz, int kstrt, int nvpy, int nvpz,
                   int nzv, int kxyp, int kyzp);

void cppepoisp332(float complex dcu[], float complex exyz[], int isign,
                  float complex ffe[], float ax, float ay, float az,
                  float affp, float wp0, float ci, float *wf, int nx,
                  int ny, int nz, int kstrt, int nvpy, int nvpz,
                  int nzv, int kxyp, int kyzp, int nzhd);

void cppaddvrfield32(float a[], float b[], float c[], int ndim, int nxe,
                     int nypmx, int nzpmx);

void cwpfft32rinit(int mixup[], float complex sct[], int indx, int indy,
                   int indz, int nxhyzd, int nxyzhd);

void cppfft32rxx(float complex f[], int isign, int mixup[],
                 float complex sct[], int indx, int indy, int indz,
                 int kstrt, int nvp, int kypi, int kypp, int nxvh,
                 int kzpp, int kypd, int kzpd, int nxhyzd, int nxyzhd);

void cppfft32rxy(float complex g[], int isign, int mixup[],
                 float complex sct[], int indx, int indy, int indz,
                 int kstrt, int nvpy, int nvpz, int kxypi, int kxypp,
                 int nyv, int kzpp, int kxypd, int kzpd, int nxhyzd,
                 int nxyzhd);

void cppfft32rxz(float complex h[], int isign, int mixup[],
                 float complex sct[], int indx, int indy, int indz,
                 int kstrt, int nvpy, int nvpz, int kxypi, int kxypp,
                 int nzv, int kyzp, int kxypd, int kyzpd, int nxhyzd,
                 int nxyzhd);

void cppfft32r3xx(float complex f[], int isign, int mixup[],
                  float complex sct[], int indx, int indy, int indz,
                  int kstrt, int nvp, int kypi, int kypp, int nxvh,
                  int kzpp, int kypd, int kzpd, int nxhyzd, int nxyzhd);

void cppfft32r3xy(float complex g[], int isign, int mixup[],
                  float complex sct[], int indx, int indy, int indz,
                  int kstrt, int nvpy, int nvpz, int kxypi, int kxypp,
                  int nyv, int kzpp, int kxypd, int kzpd, int nxhyzd,
                  int nxyzhd);

void cppfft32r3xz(float complex h[], int isign, int mixup[],
                  float complex sct[], int indx, int indy, int indz,
                  int kstrt, int nvpy, int nvpz, int kxypi, int kxypp,
                  int nzv, int kyzp, int kxypd, int kyzpd, int nxhyzd,
                  int nxyzhd);

void cppfft32rnxx(float complex f[], float complex ss[], int isign,
                  int mixup[], float complex sct[], int indx, int indy,
                  int indz, int kstrt, int nvp, int kypi, int kypp,
                  int nxvh, int kzpp, int kypd, int kzpd, int ndim,
                  int nxhyzd, int nxyzhd);

void cppfft32rnxy(float complex g[], int isign, int mixup[],
                  float complex sct[], int indx, int indy, int indz,
                  int kstrt, int nvpy, int nvpz, int kxypi, int kxypp,
                  int nyv, int kzpp, int kxypd, int kzpd, int ndim,
                  int nxhyzd, int nxyzhd);

void cppfft32rnxz(float complex h[], int isign, int mixup[],
                  float complex sct[], int indx, int indy, int indz,
                  int kstrt, int nvpy, int nvpz, int kxypi, int kxypp,
                  int nzv, int kyzp, int kxypd, int kyzpd, int ndim,
                  int nxhyzd, int nxyzhd);

void cwppfft32r(float complex f[], float complex g[], float complex h[],
                float complex bs[], float complex br[], int isign,
                int ntpose, int mixup[], float complex sct[],
                float *ttp, int indx, int indy, int indz, int kstrt,
                int nvpy, int nvpz, int nxvh, int nyv, int nzv,
                int kxyp, int kyp, int kyzp, int kzp, int kxypd,
                int kypd, int kyzpd, int kzpd, int kzyp, int nxhyzd,
                int nxyzhd);

void cwppfft32r3(float complex f[], float complex g[],
                 float complex h[], float complex bs[],
                 float complex br[], int isign, int ntpose, int mixup[],
                 float complex sct[], float *ttp, int indx, int indy,
                 int indz, int kstrt, int nvpy, int nvpz, int nxvh,
                 int nyv, int nzv, int kxyp, int kyp, int kyzp, int kzp,
                 int kxypd, int kypd, int kyzpd, int kzpd, int kzyp,
                 int nxhyzd, int nxyzhd);

void cwppfft32rn(float complex f[], float complex g[],
                 float complex h[], float complex bs[],
                 float complex br[], float complex ss[], int isign,
                 int ntpose, int mixup[], float complex sct[], 
                 float *ttp, int indx, int indy, int indz, int kstrt,
                 int nvpy, int nvpz, int nxvh, int nyv, int nzv,
                 int kxyp, int kyp, int kyzp, int kzp, int kxypd,
                 int kypd, int kyzpd, int kzpd, int kzyp, int ndim,
                 int nxhyzd, int nxyzhd);

void cppswapc32n(float f[], float s[], int isign, int nxh, int kypi,
                 int kypt, int nxvh, int kzpp, int kypd, int kzpd,
                 int ndim);

