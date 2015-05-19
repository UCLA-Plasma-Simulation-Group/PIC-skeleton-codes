/* header file for field3.c */

void cpotp3(float complex q[], float complex pot[], float complex ffc[],
            float *we, int nx, int ny, int nz, int nxvh, int nyv,
            int nzv, int nxhd, int nyhd, int nzhd);

void cdivf3(float complex f[], float complex df[], int nx, int ny,
            int nz, int nxvh, int nyv, int nzv);

void cgradf3(float complex df[], float complex f[], int nx, int ny,
             int nz, int nxvh, int nyv, int nzv);

void csmooth3(float complex q[], float complex qs[],
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
