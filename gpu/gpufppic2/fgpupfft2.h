/* C header file for gpupfft2.cuf */

#include <complex.h>

void fgpupfft2rrcuinit(int nx, int kypp, int ndim);

void fgpupfft2cuinit(int kxpp, int ny, int ndim);

void fgpupfft2rrcudel();

void fgpupfft2cudel();

void fgpupfft2rrcux(float complex f[], float complex bsm[], int isign,
                    int indx, int indy, int kstrt, int nvp, int kxp1,
                    int kyp, int nxh1d, int kypd);

void fgpupfft2rrcuy(float complex g[], float complex brm[], int isign,
                    int indx, int indy, int kstrt, int nvp, int kxp1,
                    int kyp, int nyd);

void fgpupfft2rrcuxn(float complex fn[], float complex bsm[], int isign,
                     int indx, int indy, int ndim, int kstrt, int nvp,
                     int kxp1, int kyp, int nxh1d, int kypd);

void fgpupfft2rrcuyn(float complex gn[], float complex brm[], int isign,
                     int indx, int indy, int ndim, int kstrt, int nvp,
                     int kxp1, int kyp, int nyd);

void fgpuppsltpose(float complex f[], float complex g[], float ani,
                   int nx, int ny, int kxp, int kyp, int kstrt, int nxv,
                   int nyv);

void fgpuppsltposen(float complex fn[], float complex gn[], float ani,
                    int nx, int ny, int kxp, int kyp, int kstrt,
                    int ndim, int nxv, int nyv);
