/* C header file for gpupfft2.cu */

void gpupfft2rrcuinit(int nx, int kypp, int ndim);

void gpupfft2cuinit(int kxpp, int ny, int ndim);

void gpupfft2rrcudel();

void gpupfft2cudel();

void gpupfft2rrcux(float complex f[], float complex bsm[], int isign,
                   int indx, int indy, int kstrt, int nvp, int kxp1,
                   int kyp, int nxh1d, int kypd);

void gpupfft2rrcuy(float complex g[], float complex brm[], int isign,
                   int indx, int indy, int kstrt, int nvp, int kxp1,
                   int kyp, int nyd);

void gpupfft2rrcuxn(float complex fn[], float complex bsm[], int isign,
                    int indx, int indy, int ndim, int kstrt, int nvp,
                    int kxp1, int kyp, int nxh1d, int kypd);

void gpupfft2rrcuyn(float complex gn[], float complex brm[], int isign,
                    int indx, int indy, int ndim, int kstrt, int nvp,
                    int kxp1, int kyp, int nyd);

void cgpuppsltpose(float complex f[], float complex g[], float ani,
                   int nx, int ny, int kxp, int kyp, int kstrt, int nxv,
                   int nyv);

void cgpuppsltposen(float complex fn[], float complex gn[], float ani,
                    int nx, int ny, int kxp, int kyp, int kstrt,
                    int ndim, int nxv, int nyv);
