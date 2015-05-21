/* C header file for gpplib2.cuf */

void fgppcaguard2l(float complex g_q[], float g_qe[], float g_scs[],
                   float scs[], float scr[], int nx, int nyp, int kstrt,
                   int nvp, int nxe, int nypmx, int nxvh, int kypd);

void fgppccguard2l(float complex g_fxy[], float g_fxye[], float g_scs[],
                   float scs[], float scr[], int nx, int nyp, int kstrt,
                   int nvp, int ndim, int nxe, int nypmx, int nxvh,
                   int kypd);

void fwappfft2rcs(float complex g_f[], float complex g_g[],
                  float complex g_bsm[], float complex g_brm[],
                  float complex bsm[], float complex brm[], int isign,
                  int g_mixup[], float complex g_sct[], float ttp[],
                  int indx, int indy, int kstrt, int nvp, int kxpd,
                  int kyp, int nxhd, int nyd, int kypd, int nxhyd,
                  int nxyhd);

void fwappfft2rcsn(float complex g_fn[], float complex g_gn[],
                   float complex g_bsm[], float complex g_brm[],
                   float complex bsm[], float complex brm[], int isign,
                   int g_mixup[], float complex g_sct[], float ttp[],
                   int indx, int indy, int kstrt, int nvp, int ndim,
                   int kxpd, int kyp, int nxhd, int nyd, int kypd,
                   int nxhyd, int nxyhd);

void fgpuppfft2rrcu(float complex g_f[], float complex g_g[],
                    float complex g_bsm[], float complex g_brm[],
                    float complex bsm[], float complex brm[], int isign,
                    float ttp[], int indx, int indy, int kstrt, int nvp,
                    int kxpd, int kyp, int nxhd, int nyd, int kypd);

void fgpuppfft2rrcun(float complex g_fn[], float complex g_gn[],
                     float complex g_bsm[], float complex g_brm[],
                     float complex bsm[], float complex brm[], int isign,
                     float ttp[], int indx, int indy, int kstrt, int nvp,
                     int ndim, int kxpd, int kyp, int nxhd, int nyd,
                     int kypd);

void fgpporder2l(float g_ppart[], float g_ppbuff[], float g_sbufl[],
                 float g_sbufr[], int g_kpic[], int g_ncl[],
                 int g_ihole[], int g_ncll[], int g_nclr[], 
                 float sbufl[], float sbufr[], float rbufl[], 
                 float rbufr[], int ncll[], int nclr[], int mcll[],
                 int mclr[], float ttp[], int noff, int nyp, int kstrt,
                 int nvp, int idimp, int nppmx, int nx, int ny, int mx,
                 int my, int mx1, int myp1, int npbmx, int ntmax,
                 int nbmax, int *g_irc);

