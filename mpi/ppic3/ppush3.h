/* header file for ppush3.c */

#include <complex.h>

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

void cppgpush32l(float part[], float fxyz[], float edges[], int npp,
                 int noff[], int ihole[], float qbm, float dt,
                 float *ek, int nx, int ny, int nz, int idimp,
                 int npmax, int nxv, int nypmx, int nzpmx, int idps,
                 int idds, int ntmax, int ipbc);

void cppgpost32l(float part[], float q[], int npp, int noff[], float qm,
                 int idimp, int npmax, int nxv, int nypmx, int nzpmx,
                 int idds);

void cppdsortp32yzl(float parta[], float partb[], int npic[], int npp,
                    int noff[], int nyzp[], int idimp, int npmax,
                    int nyzpm1, int idds);

void cppcguard32xl(float fxyz[], int nyzp[], int nx, int ndim, int nxe,
                   int nypmx, int nzpmx, int idds);

void cppaguard32xl(float q[], int nyzp[], int nx, int nxe, int nypmx,
                   int nzpmx, int idds);

void cppois332(float complex q[], float complex fxyz[], int isign,
               float complex ffc[], float ax, float ay, float az,
               float affp, float *we, int nx, int ny, int nz, int kstrt,
               int nvpy, int nvpz, int nzv, int kxyp, int kyzp,
               int nzhd);

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
