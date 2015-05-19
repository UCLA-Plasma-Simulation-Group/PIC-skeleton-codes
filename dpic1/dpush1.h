/* header file for dpush1.c */

double ranorm();

double randum();

void cdistr1h(float part[], float vtx, float vty, float vtz, float vdx,
              float vdy, float vdz, int npx, int idimp, int nop, int nx,
              int ipbc);

void cgbpush13l(float part[], float fxyz[], float byz[], float omx,
                float qbm, float dt, float dtc, float *ek, int idimp,
                int nop, int nx, int nxv, int ipbc);

void cgpost1l(float part[], float q[], float qm, int nop, int idimp,
              int nxv);

void cgjpost1l(float part[], float cu[], float qm, float dt, int nop,
               int idimp, int nx, int nxv, int ipbc);

void cgmjpost1l(float part[], float amu[], float qm, int nop, int idimp,
                int nxv);

void cgdjpost1l(float part[], float fxyz[], float byz[], float dcu[],
                float amu[], float omx, float qm, float qbm, float dt,
                int idimp, int nop, int nxv);

void cgdcjpost1l(float part[], float fxyz[], float byz[], float cu[],
                 float dcu[], float amu[], float omx, float qm,
                 float qbm, float dt, int idimp, int nop, int nxv);

void cdsortp1xl(float parta[], float partb[], int npic[], int idimp,
                int nop, int nx1);

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

