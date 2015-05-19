/* header file for mdpush1.c */

double ranorm();

void cdistr1h(float part[], float vtx, float vty, float vtz, float vdx,
              float vdy, float vdz, int npx, int idimp, int nop, int nx,
              int ipbc);

void cdblkp1l(float part[], int kpic[], int *nppmx, int idimp, int nop,
              int mx, int mx1, int *irc);

void cppmovin1l(float part[], float ppart[], int kpic[], int nppmx,
                int idimp, int nop, int mx, int mx1, int *irc);

void cppcheck1l(float ppart[], int kpic[], int idimp, int nppmx, int nx,
                int mx, int mx1, int *irc);

void cgbppush13l(float ppart[], float fxyz[], float byz[], int kpic[],
                 float omx, float qbm, float dt, float dtc, float *ek,
                 int idimp, int nppmx, int nx, int mx, int nxv, int mx1,
                 int ipbc);

void cgbppushf13l(float ppart[], float fxyz[], float byz[], int kpic[],
                  int ncl[], int ihole[], float omx, float qbm,
                  float dt, float dtc, float *ek, int idimp, int nppmx,
                  int nx, int mx, int nxv, int mx1, int ntmax,
                  int *irc);

void cgppost1l(float ppart[], float q[], int kpic[], float qm,
               int nppmx, int idimp, int mx, int nxv, int mx1);

void cgjppost1l(float ppart[], float cu[], int kpic[], float qm,
                float dt, int nppmx, int idimp, int nx, int mx, int nxv,
                int mx1, int ipbc);

void cgmjppost1l(float ppart[], float amu[], int kpic[], float qm,
                 int nppmx, int idimp, int mx, int nxv, int mx1);

void cgdjppost1l(float ppart[], float fxyz[], float byz[], float dcu[],
                 float amu[], int kpic[], float omx, float qm,
                 float qbm, float dt, int idimp, int nppmx, int nx,
                 int mx, int nxv, int mx1);

void cgdcjppost1l(float ppart[], float fxyz[], float byz[], float cu[],
                  float dcu[], float amu[], int kpic[], float omx,
                  float qm, float qbm, float dt, int idimp, int nppmx,
                  int nx, int mx, int nxv, int mx1);

void cgmjpost1l(float part[], float amu[], float qm, int nop, int idimp,
                int nxv);

void cgdjpost1l(float part[], float fxyz[], float byz[], float dcu[],
                float amu[], float omx, float qm, float qbm, float dt,
                int idimp, int nop, int nxv);

void cgdcjpost1l(float part[], float fxyz[], float byz[], float cu[],
                 float dcu[], float amu[], float omx, float qm,
                 float qbm, float dt, int idimp, int nop, int nxv);

void cpporder1l(float ppart[], float ppbuff[], int kpic[], int ncl[],
                int ihole[], int idimp, int nppmx, int nx, int mx,
                int mx1, int npbmx, int ntmax, int *irc);

void cpporderf1l(float ppart[], float ppbuff[], int kpic[], int ncl[],
                 int ihole[], int idimp, int nppmx, int mx1, int npbmx,
                 int ntmax, int *irc);

void cdguard1l(float fx[], int nx, int nxe);

void ccguard1l(float byz[], int nx, int nxe);

void cacguard1l(float cu[], int nx, int nxe);

void caguard1l(float q[], int nx, int nxe);

void cascfguard1l(float dcu[], float cus[], float q2m0, int nx,
                  int nxe);

void cfwpminmx1(float qe[], float qbme, float *wpmax, float *wpmin,
                int nx, int nxe);

void cpois1(float complex q[], float complex fx[], int isign,
            float complex ffc[], float ax, float affp, float *we,
            int nx);

void cbbpois13(float complex cu[], float complex byz[],
               float complex ffc[], float ci, float *wm, int nx,
               int nxvh, int nxhd);

void cbaddext1(float byz[], float omy, float omz, int nx, int nxe);

void cdcuperp13(float complex dcu[], float complex amu[], int nx,
                int nxvh);

void cadcuperp13(float complex dcu[], float complex amu[], int nx,
                 int nxvh);

void cepois13(float complex dcu[], float complex eyz[], int isign,
              float complex ffe[], float ax, float affp, float wp0,
              float ci, float *wf, int nx, int nxvh, int nxhd);

void caddvrfield13(float fxyze[], float eyze[], float fxe[], int nxe);

void cwfft1rinit(int mixup[], float complex sct[], int indx, int nxhd);

void cfft1rxx(float complex f[], float complex t[], int isign,
              int mixup[], float complex sct[], int indx, int nxd, 
              int nxhd);

void cfft1r2x(float complex f[], float complex t[], int isign,
              int mixup[], float complex sct[], int indx, int nxd,
              int nxhd);

