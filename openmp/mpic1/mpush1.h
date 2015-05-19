/* header file for mpush1.c */

double ranorm();

void cdistr1(float part[], float vtx, float vdx, int npx, int idimp,
             int nop, int nx, int ipbc);

void cdblkp1l(float part[], int kpic[], int *nppmx, int idimp, int nop,
              int mx, int mx1, int *irc);

void cppmovin1l(float part[], float ppart[], int kpic[], int nppmx,
                int idimp, int nop, int mx, int mx1, int *irc);

void cppcheck1l(float ppart[], int kpic[], int idimp, int nppmx, int nx,
                int mx, int mx1, int *irc);

void cgppush1l(float ppart[], float fx[], int kpic[], float qbm,
               float dt, float *ek, int idimp, int nppmx, int nx,
               int mx, int nxv, int mx1, int ipbc);

void cgppushf1l(float ppart[], float fx[], int kpic[], int ncl[],
                int ihole[], float qbm, float dt, float *ek, int idimp,
                int nppmx, int nx, int mx, int nxv, int mx1, int ntmax,
                int *irc);

void cgppost1l(float ppart[], float q[], int kpic[], float qm,
               int nppmx, int idimp, int mx, int nxv, int mx1);

void cpporder1l(float ppart[], float ppbuff[], int kpic[], int ncl[],
                int ihole[], int idimp, int nppmx, int nx, int mx,
                int mx1, int npbmx, int ntmax, int *irc);

void cpporderf1l(float ppart[], float ppbuff[], int kpic[], int ncl[],
                 int ihole[], int idimp, int nppmx, int mx1, int npbmx,
                 int ntmax, int *irc);

void ccguard1l(float fx[], int nx, int nxe);

void caguard1l(float q[], int nx, int nxe);

void cpois1(float complex q[], float complex fx[], int isign,
            float complex ffc[], float ax, float affp, float *we,
            int nx);

void cwfft1rinit(int mixup[], float complex sct[], int indx, int nxhd);

void cfft1rxx(float complex f[], float complex t[], int isign,
              int mixup[], float complex sct[], int indx, int nxd, 
              int nxhd);
