/* header file for bfield2.c */

void cpotp2(float complex q[], float complex pot[], float complex ffc[],
            float *we, int nx, int ny, int nxvh, int nyv, int nxhd,
            int nyhd);

void cdivf2(float complex f[], float complex df[], int nx, int ny,
            int nxvh, int nyv);

void cgradf2(float complex df[], float complex f[], int nx, int ny,
             int nxvh, int nyv);

void ccurlf2(float complex f[], float complex g[], int nx, int ny,
             int nxvh, int nyv);

void cavpot23(float complex bxy[], float complex axy[], int nx, int ny,
              int nxvh, int nyv);

void cavrpot23(float complex axy[], float complex bxy[],
               float complex ffc[], float ci, int nx, int ny, int nxvh,
               int nyv, int nxhd, int nyhd);

void csmooth2(float complex q[], float complex qs[],
              float complex ffc[], int nx, int ny, int nxvh, int nyv,
              int nxhd, int nyhd);

void csmooth23(float complex cu[], float complex cus[],
               float complex ffc[], int nx, int ny, int nxvh, int nyv,
               int nxhd, int nyhd);

void crdmodes2(float complex pot[], float complex pott[], int nx,
               int ny, int modesx, int modesy, int nxvh, int nyv,
               int modesxd, int modesyd);

void cwrmodes2(float complex pot[], float complex pott[], int nx, 
               int ny, int modesx, int modesy, int nxvh, int nyv, 
               int modesxd, int modesyd);

void crdvmodes2(float complex vpot[], float complex vpott[], int nx,
                int ny, int modesx, int modesy, int ndim, int nxvh,
                int nyv, int modesxd, int modesyd);

void cwrvmodes2(float complex vpot[], float complex vpott[], int nx,
                int ny, int modesx, int modesy, int ndim, int nxvh,
                int nyv, int modesxd, int modesyd);
