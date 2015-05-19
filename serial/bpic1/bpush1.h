/* header file for bpush1.c */

double ranorm();

double randum();

void cdistr1h(float part[], float vtx, float vty, float vtz, float vdx,
              float vdy, float vdz, int npx, int idimp, int nop, int nx,
              int ipbc);

void cgbpush13l(float part[], float fxyz[], float byz[], float omx,
                float qbm, float dt, float dtc, float *ek, int idimp,
                int nop, int nx, int nxv, int ipbc);

void cgrbpush13l(float part[], float fxyz[], float byz[], float omx,
                 float qbm, float dt, float dtc, float ci, float *ek,
                 int idimp, int nop, int nx, int nxv, int ipbc);

void cgpost1l(float part[], float q[], float qm, int nop, int idimp,
              int nxv);

void cgjpost1l(float part[], float cu[], float qm, float dt, int nop,
               int idimp, int nx, int nxv, int ipbc);

void cgrjpost1l(float part[], float cu[], float qm, float dt, float ci,
                int nop, int idimp, int nx, int nxv, int ipbc);

void cdsortp1xl(float parta[], float partb[], int npic[], int idimp,
                int nop, int nx1);

void ccguard1l(float byz[], int nx, int nxe);

void cbguard1l(float fxyz[], int nx, int nxe);

void cacguard1l(float cu[], int nx, int nxe);

void caguard1l(float q[], int nx, int nxe);

void cpois1(float complex q[], float complex fx[], int isign,
            float complex ffc[], float ax, float affp, float *we,
            int nx);

void cibpois13(float complex cu[], float complex byz[],
               float complex ffc[], float ci, float *wm, int nx,
               int nxvh, int nxhd);

void cmaxwel1(float complex eyz[], float complex byz[],
              float complex cu[], float complex ffc[], float ci,
              float dt, float *wf, float *wm, int nx, int nxvh,
              int nxhd);

void cemfield1(float complex fxyz[], float complex fx[],
               float complex eyz[], float complex ffc[], int nx,
               int nxvh, int nxhd);

void cbmfield1(float complex fyz[], float complex eyz[],
               float complex ffc[], int nx, int nxvh, int nxhd);

void cwfft1rinit(int mixup[], float complex sct[], int indx, int nxhd);

void cfft1rxx(float complex f[], float complex t[], int isign,
              int mixup[], float complex sct[], int indx, int nxd, 
              int nxhd);

void cfft1r2x(float complex f[], float complex t[], int isign,
              int mixup[], float complex sct[], int indx, int nxd,
              int nxhd);

void cfft1r3x(float complex f[], float complex t[], int isign,
              int mixup[], float complex sct[], int indx, int nxd,
              int nxhd);
