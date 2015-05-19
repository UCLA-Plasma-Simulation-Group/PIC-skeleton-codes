/* header file for mfield2.c */

void cmpotp2(float complex q[], float complex pot[],
             float complex ffc[], float *we, int nx, int ny, int nxvh,
             int nyv, int nxhd, int nyhd);

void cmdivf2(float complex f[], float complex df[], int nx, int ny,
             int nxvh, int nyv);

void cmgradf2(float complex df[], float complex f[], int nx, int ny,
              int nxvh, int nyv);

void cmsmooth2(float complex q[], float complex qs[],
               float complex ffc[], int nx, int ny, int nxvh, int nyv,
               int nxhd, int nyhd);

void crdmodes2(float complex pot[], float complex pott[], int nx,
               int ny, int modesx, int modesy, int nxvh, int nyv,
               int modesxd, int modesyd);

void cwrmodes2(float complex pot[], float complex pott[], int nx, 
               int ny, int modesx, int modesy, int nxvh, int nyv, 
               int modesxd, int modesyd);
