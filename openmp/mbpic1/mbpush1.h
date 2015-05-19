/* header file for mbpush1.c */

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

void cgrbppush13l(float ppart[], float fxyz[], float byz[], int kpic[],
                  float omx, float qbm, float dt, float dtc, float ci,
                  float *ek, int idimp, int nppmx, int nx, int mx,
                  int nxv, int mx1, int ipbc);
     
void cgrbppushf13l(float ppart[], float fxyz[], float byz[], int kpic[],
                   int ncl[], int ihole[], float omx, float qbm,
                   float dt, float dtc, float ci, float *ek, int idimp,
                   int nppmx, int nx, int mx, int nxv, int mx1,
                   int ntmax, int *irc);

void cgppost1l(float ppart[], float q[], int kpic[], float qm,
               int nppmx, int idimp, int mx, int nxv, int mx1);

void cgjppost1l(float ppart[], float cu[], int kpic[], float qm,
                float dt, int nppmx, int idimp, int nx, int mx, int nxv,
                int mx1, int ipbc);

void cgjppostf1l(float ppart[], float cu[], int kpic[], int ncl[],
                 int ihole[], float qm, float dt, int nppmx, int idimp,
                 int nx, int mx, int nxv, int mx1, int ntmax, int *irc);

void cgrjppost1l(float ppart[], float cu[], int kpic[], float qm,
                 float dt, float ci, int nppmx, int idimp, int nx,
                 int mx, int nxv, int mx1, int ipbc);

void cgrjppostf1l(float ppart[], float cu[], int kpic[], int ncl[],
                  int ihole[], float qm, float dt, float ci, int nppmx,
                  int idimp, int nx, int mx, int nxv, int mx1,
                  int ntmax, int *irc);

void cpporder1l(float ppart[], float ppbuff[], int kpic[], int ncl[],
                int ihole[], int idimp, int nppmx, int nx, int mx,
                int mx1, int npbmx, int ntmax, int *irc);

void cpporderf1l(float ppart[], float ppbuff[], int kpic[], int ncl[],
                 int ihole[], int idimp, int nppmx, int mx1, int npbmx,
                 int ntmax, int *irc);

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
