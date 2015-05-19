/* header file for pdpush2.c */

double ranorm();

void cpdicomp2l(float edges[], int *nyp, int *noff, int *nypmx,
                int *nypmn, int ny, int kstrt, int nvp, int idps);

void cpdistr2h(float part[], float edges[], int *npp, int nps,
               float vtx, float vty, float vtz, float vdx, float vdy,
               float vdz, int npx, int npy, int nx, int ny, int idimp,
               int npmax, int idps, int ipbc, int *ierr);

void cppgbpush23l(float part[], float fxy[], float bxy[], float edges[],
                  int npp, int noff, int ihole[], float qbm, float dt, 
                  float dtc, float *ek, int nx, int ny, int idimp,
                  int npmax, int nxv, int nypmx, int idps, int ntmax,
                  int ipbc);

void cppgpost2l(float part[], float q[], int npp, int noff, float qm,
                int idimp, int npmax, int nxv, int nypmx);

void cppgjpost2l(float part[], float cu[], float edges[], int npp,
                 int noff, int ihole[], float qm, float dt, int nx,
                 int ny, int idimp, int npmax, int nxv, int nypmx,
                 int idps, int ntmax, int ipbc);

void cppgmjpost2l(float part[], float amu[], int npp, int noff,
                  float qm, int idimp, int npmax, int nxv, int nypmx);

void cppgdjpost2l(float part[], float fxy[], float bxy[], int npp,
                  int noff, float dcu[], float amu[], float qm,
                  float qbm, float dt, int idimp, int npmax, int nxv,
                  int nypmx);

void cppgdcjpost2l(float part[], float fxy[], float bxy[], int npp,
                   int noff, float cu[], float dcu[], float amu[],
                   float qm, float qbm, float dt, int idimp, int npmax,
                   int nxv, int nypmx);

void cppdsortp2yl(float parta[], float partb[], int npic[], int npp,
                  int noff, int nyp, int idimp, int npmax, int nypm1);

void cppcguard2xl(float fxy[], int myp, int nx, int ndim, int nxe,
                  int nypmx);

void cppaguard2xl(float q[], int myp, int nx, int nxe, int nypmx);

void cppacguard2xl(float cu[], int myp, int nx, int ndim, int nxe,
                   int nypmx);

void cppascfguard2l(float dcu[], float cus[], int nyp, float q2m0,
                    int nx, int nxe, int nypmx);

void cppfwpminmx2(float qe[], int nyp, float qbme, float *wpmax,
                  float *wpmin, int nx, int nxe, int nypmx);

void cppois23(float complex q[], float complex fxy[], int isign,
              float complex ffc[], float ax, float ay, float affp,
              float *we, int nx, int ny, int kstrt, int nyv, int kxp,
              int nyhd);

void cppcuperp2(float complex cu[], int nx, int ny, int kstrt, int nyv,
                int kxp);

void cppbbpoisp23(float complex cu[], float complex bxy[],
                  float complex ffc[], float ci, float *wm, int nx,
                  int ny, int kstrt, int nyv, int kxp, int nyhd);

void cppbaddext2(float bxy[], int nyp, float omx, float omy, float omz,
                 int nx, int nxe, int nypmx);

void cppdcuperp23(float complex dcu[], float complex amu[], int nx,
                  int ny, int kstrt, int nyv, int kxp);

void cppadcuperp23(float complex dcu[], float complex amu[], int nx,
                   int ny, int kstrt, int nyv, int kxp);

void cppepoisp23(float complex dcu[], float complex exy[], int isign,
                 float complex ffe[], float ax, float ay, float affp,
                 float wp0, float ci, float *wf, int nx, int ny,
                 int kstrt, int nyv, int kxp, int nyhd);

void cppaddvrfield2(float a[], float b[], float c[], int ndim, int nxe,
                    int nypmx);

void cwpfft2rinit(int mixup[], float complex sct[], int indx, int indy,
                  int nxhyd, int nxyhd);

void cppfft2rxx(float complex f[], int isign, int mixup[],
                float complex sct[], int indx, int indy, int kstrt,
                int kypi, int kypp, int nxvh, int kypd, int nxhyd,
                int nxyhd);

void cppfft2rxy(float complex g[], int isign, int mixup[],
                float complex sct[], int indx, int indy, int kstrt,
                int kxpi, int kxpp, int nyv, int kxp, int nxhyd,
                int nxyhd);

void cppfft2r3xx(float complex f[], int isign, int mixup[],
                 float complex sct[], int indx, int indy, int kstrt,
                 int kypi, int kypp, int nxvh, int kypd, int nxhyd,
                 int nxyhd);

void cppfft2r3xy(float complex g[], int isign, int mixup[],
                 float complex sct[], int indx, int indy, int kstrt,
                 int kxpi, int kxpp, int nyv, int kxp, int nxhyd,
                 int nxyhd);

void cppfft2rnxx(float complex f[], float complex ss[], int isign,
                 int mixup[], float complex sct[], int indx, int indy,
                 int kstrt, int kypi, int kypp, int nxvh, int kypd,
                 int ndim, int nxhyd, int nxyhd);

void cppfft2rnxy(float complex g[], int isign, int mixup[],
                 float complex sct[], int indx, int indy, int kstrt,
                 int kxpi, int kxpp, int nyv, int kxp, int ndim,
                 int nxhyd, int nxyhd);

void cwppfft2r(float complex f[], float complex g[], float complex bs[],
               float complex br[], int isign, int ntpose, int mixup[],
               float complex sct[], float *ttp, int indx, int indy,
               int kstrt, int nvp, int nxvh, int nyv, int kxp, int kyp,
               int kypd, int nxhyd, int nxyhd);

void cwppfft2r3(float complex f[], float complex g[], float complex bs[],
                float complex br[], int isign, int ntpose, int mixup[],
                float complex sct[], float *ttp, int indx, int indy,
                int kstrt, int nvp, int nxvh, int nyv, int kxp, int kyp,
                int kypd, int nxhyd, int nxyhd);

void cwppfft2rn(float complex f[], float complex g[],
                float complex bs[], float complex br[],
                float complex ss[], int isign, int ntpose, int mixup[],
                float complex sct[], float *ttp, int indx, int indy,
                int kstrt, int nvp, int nxvh, int nyv, int kxp, int kyp,
                int kypd, int ndim, int nxhyd, int nxyhd);

void cppswapc2n(float f[], float s[], int isign, int nxh, int kypi,
                int kypt, int nxvh, int kypd, int ndim);
                
                
