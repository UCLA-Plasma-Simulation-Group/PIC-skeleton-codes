/* header file for field2.c */

void cpotp2(float complex q[], float complex pot[], float complex ffc[],
            float *we, int nx, int ny, int nxvh, int nyv, int nxhd,
            int nyhd);

void cdivf2(float complex f[], float complex df[], int nx, int ny,
            int nxvh, int nyv);

void cgradf2(float complex df[], float complex f[], int nx, int ny,
             int nxvh, int nyv);

void csmooth2(float complex q[], float complex qs[],
              float complex ffc[], int nx, int ny, int nxvh, int nyv,
              int nxhd, int nyhd);

void crdmodes2(float complex pot[], float complex pott[], int nx,
               int ny, int modesx, int modesy, int nxvh, int nyv,
               int modesxd, int modesyd);

void cwrmodes2(float complex pot[], float complex pott[], int nx, 
               int ny, int modesx, int modesy, int nxvh, int nyv, 
               int modesxd, int modesyd);
