/* C header file for gpufft2.cuf */

#include <complex.h>

void fgpufft2rrcuinit(int nx, int ny, int ndim);

void fgpufft2cuinit(int nx, int ny, int ndim);

void fgpufft2rrcudel();

void fgpufft2cudel();

void fgpufft2rrcu(float complex f[], float complex g[], int isign,
                 int indx, int indy, int nxh1d, int nyd);

void fgpufft2rrcun(float complex fn[], float complex gn[], int isign,
                  int indx, int indy, int ndim, int nxh1d, int nyd);
