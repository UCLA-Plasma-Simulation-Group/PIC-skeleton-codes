/* header file for push1.c */

double ranorm();

double randum();

void cdistr1(float part[], float vtx, float vdx, int npx, int idimp,
             int nop, int nx, int ipbc);

void cgpush1l(float part[], float fx[], float qbm, float dt, float *ek,
              int idimp, int nop, int nx, int nxv, int ipbc);

void cgpost1l(float part[], float q[], float qm, int nop, int idimp,
              int nxv);

void cdsortp1xl(float parta[], float partb[], int npic[], int idimp,
                int nop, int nx1);

void ccguard1l(float fx[], int nx, int nxe);

void caguard1l(float q[], int nx, int nxe);

void cpois1(float complex q[], float complex fx[], int isign,
            float complex ffc[], float ax, float affp, float *we,
            int nx);

void cwfft1rinit(int mixup[], float complex sct[], int indx, int nxhd);

void cfft1rxx(float complex f[], float complex t[], int isign,
              int mixup[], float complex sct[], int indx, int nxd, 
              int nxhd);
