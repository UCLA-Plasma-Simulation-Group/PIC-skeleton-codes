/* header file for field1.c */

void cpotp1(float complex q[], float complex pot[], float complex ffc[],
            float *we, int nx, int nxvh, int nxhd);

void cdivf1(float complex f[], float complex df[], int nx, int ndim,
            int nxvh);

void cgradf1(float complex df[], float complex f[], int nx, int ndim,
             int nxvh);

void csmooth1(float complex q[], float complex qs[],
              float complex ffc[], int nx, int nxvh, int nxhd);

void crdmodes1(float complex pot[], float complex pott[], int nx,
               int modesx, int nxvh, int modesxd);

void cwrmodes1(float complex pot[], float complex pott[], int nx,
               int modesx, int nxvh, int modesxd);
