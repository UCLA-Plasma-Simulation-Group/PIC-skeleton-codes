/* header file for mbfield2.c */

void cmpotp2(float complex q[], float complex pot[],
             float complex ffc[], float *we, int nx, int ny, int nxvh,
             int nyv, int nxhd, int nyhd);

void cmdivf2(float complex f[], float complex df[], int nx, int ny,
             int nxvh, int nyv);

void cmgradf2(float complex df[], float complex f[], int nx, int ny,
              int nxvh, int nyv);

void cmcurlf2(float complex f[], float complex g[], int nx, int ny,
              int nxvh, int nyv);

void cmavpot23(float complex bxy[], float complex axy[], int nx, int ny,
               int nxvh, int nyv);

void cmavrpot23(float complex axy[], float complex bxy[],
                float complex ffc[], float ci, int nx, int ny, int nxvh,
                int nyv, int nxhd, int nyhd);

void cmsmooth2(float complex q[], float complex qs[],
               float complex ffc[], int nx, int ny, int nxvh, int nyv,
               int nxhd, int nyhd);

void cmsmooth23(float complex cu[], float complex cus[],
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
