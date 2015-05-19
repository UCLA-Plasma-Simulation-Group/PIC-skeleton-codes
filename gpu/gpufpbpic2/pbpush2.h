/* header file for pbpush2.c */

double ranorm();

void cpdicomp2l(float edges[], int *nyp, int *noff, int *nypmx,
                int *nypmn, int ny, int kstrt, int nvp, int idps);

void cpdistr2h(float part[], float edges[], int *npp, int nps,
               float vtx, float vty, float vtz, float vdx, float vdy,
               float vdz, int npx, int npy, int nx, int ny, int idimp,
               int npmax, int idps, int ipbc, int *ierr);

void cppdblkp2l(float part[], int kpic[], int npp, int noff, int *nppmx,
                int idimp, int npmax, int mx, int my, int mx1,
                int mxyp1, int *irc);

void cpppmovin2lt(float part[], float ppart[], int kpic[], int npp,
                  int noff, int nppmx, int idimp, int npmax, int mx,
                  int my, int mx1, int mxyp1, int *irc);

void cpppcheck2lt(float ppart[], int kpic[], int noff, int nyp,
                  int idimp, int nppmx, int nx, int mx, int my, int mx1,
                  int myp1, int *irc);

void cppois23t(float complex qt[], float complex fxyt[], int isign,
               float complex ffct[], float ax, float ay, float affp,
               float *we, int nx, int ny, int kstrt, int nyv, int kxp1,
               int nyhd);

void cwpfft2rinit(int mixup[], float complex sct[], int indx, int indy,
                  int nxhyd, int nxyhd);
