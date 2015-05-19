/* header file for bfield3.c */

void cpotp3(float complex q[], float complex pot[], float complex ffc[],
            float *we, int nx, int ny, int nz, int nxvh, int nyv,
            int nzv, int nxhd, int nyhd, int nzhd);

void cdivf3(float complex f[], float complex df[], int nx, int ny,
            int nz, int nxvh, int nyv, int nzv);

void cgradf3(float complex df[], float complex f[], int nx, int ny,
             int nz, int nxvh, int nyv, int nzv);

void ccurlf3(float complex f[], float complex g[], int nx, int ny,
             int nz, int nxvh, int nyv, int nzv);

void cavpot33(float complex bxyz[], float complex axyz[], int nx,
              int ny, int nz, int nxvh, int nyv, int nzv);

void cavrpot33(float complex axyz[], float complex bxyz[],
               float complex ffc[], float ci, int nx, int ny, int nz,
               int nxvh, int nyv, int nzv, int nxhd, int nyhd, 
               int nzhd);

void csmooth3(float complex q[], float complex qs[],
              float complex ffc[], int nx, int ny, int nz, int nxvh,
              int nyv, int nzv, int nxhd, int nyhd, int nzhd);

void csmooth33(float complex cu[], float complex cus[],
               float complex ffc[], int nx, int ny, int nz, int nxvh,
               int nyv, int nzv, int nxhd, int nyhd, int nzhd);

void crdmodes3(float complex pot[], float complex pott[], int nx,
               int ny, int nz, int modesx, int modesy, int modesz,
               int nxvh, int nyv, int nzv, int modesxd, int modesyd,
               int modeszd);

void cwrmodes3(float complex pot[], float complex pott[], int nx, 
               int ny, int nz, int modesx, int modesy, int modesz,
               int nxvh, int nyv, int nzv, int modesxd, int modesyd,
               int modeszd);

void crdvmodes3(float complex vpot[], float complex vpott[], int nx,
                int ny, int nz, int modesx, int modesy, int modesz,
                int ndim, int nxvh, int nyv, int nzv, int modesxd,
                int modesyd, int modeszd);

void cwrvmodes3(float complex vpot[], float complex vpott[], int nx,
                int ny, int nz, int modesx, int modesy, int modesz,
                int ndim, int nxvh, int nyv, int nzv, int modesxd, 
                int modesyd, int modeszd);
