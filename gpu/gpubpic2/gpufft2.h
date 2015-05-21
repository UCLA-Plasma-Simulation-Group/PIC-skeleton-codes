/* C header file for gpufft2.cu */

#include <complex.h>

void gpufft2rrcuinit(int nx, int ny, int ndim);

void gpufft2cuinit(int nx, int ny, int ndim);

void gpufft2rrcudel();

void gpufft2cudel();

void gpufft2rrcu(float complex f[], float complex g[], int isign,
                 int indx, int indy, int nxh1d, int nyd);

void gpufft2rrcun(float complex fn[], float complex gn[], int isign,
                  int indx, int indy, int ndim, int nxh1d, int nyd);
