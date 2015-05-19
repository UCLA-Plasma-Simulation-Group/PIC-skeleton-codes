/* header file for mdfield3.c */

void cmpotp3(float complex q[], float complex pot[],
             float complex ffc[], float *we, int nx, int ny, int nz,
             int nxvh, int nyv, int nzv, int nxhd, int nyhd, int nzhd);

void cmdivf3(float complex f[], float complex df[], int nx, int ny,
             int nz, int nxvh, int nyv, int nzv);

void cmgradf3(float complex df[], float complex f[], int nx, int ny,
              int nz, int nxvh, int nyv, int nzv);

void cmcurlf3(float complex f[], float complex g[], int nx, int ny,
              int nz, int nxvh, int nyv, int nzv);

void cmapotp33(float complex cu[], float complex axyz[],
               float complex ffc[], float ci, float *wm, int nx, int ny,
               int nz, int nxvh, int nyv, int nzv, int nxhd, int nyhd,
               int nzhd);

void cmsmooth3(float complex q[], float complex qs[],
               float complex ffc[], int nx, int ny, int nz, int nxvh,
               int nyv, int nzv, int nxhd, int nyhd, int nzhd);

void cmsmooth33(float complex cu[], float complex cus[],
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
