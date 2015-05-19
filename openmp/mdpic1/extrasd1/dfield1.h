/* header file for dfield1.c */

void cpotp1(float complex q[], float complex pot[], float complex ffc[],
            float *we, int nx, int nxvh, int nxhd);

void cdivf1(float complex f[], float complex df[], int nx, int ndim,
            int nxvh);

void cgradf1(float complex df[], float complex f[], int nx, int ndim,
             int nxvh);

void ccurlf1(float complex f[], float complex g[], int nx, int nxvh);

void capotp13(float complex cu[], float complex ayz[],
              float complex ffc[], float ci, float *wm, int nx,
              int nxvh, int nxhd);

void csmooth1(float complex q[], float complex qs[],
              float complex ffc[], int nx, int nxvh, int nxhd);

void csmooth13(float complex cu[], float complex cus[], 
               float complex ffc[], int nx, int nxvh, int nxhd);

void crdmodes1(float complex pot[], float complex pott[], int nx,
               int modesx, int nxvh, int modesxd);

void cwrmodes1(float complex pot[], float complex pott[], int nx,
               int modesx, int nxvh, int modesxd);

void crdvmodes1(float complex vpot[], float complex vpott[], int nx,
                int modesx, int ndim, int nxvh, int modesxd);

void cwrvmodes1(float complex vpot[], float complex vpott[], int nx,
                int modesx, int ndim, int nxvh, int modesxd);

