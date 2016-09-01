/* header file for dfield2.c */

void cpotp2(float complex q[], float complex pot[], float complex ffc[],
            float *we, int nx, int ny, int nxvh, int nyv, int nxhd,
            int nyhd);

void cdivf2(float complex f[], float complex df[], int nx, int ny,
            int nxvh, int nyv);

void cgradf2(float complex df[], float complex f[], int nx, int ny,
             int nxvh, int nyv);

void ccurlf2(float complex f[], float complex g[], int nx, int ny,
             int nxvh, int nyv);

void capotp23(float complex cu[], float complex axy[],
              float complex ffc[], float ci, float *wm, int nx, int ny,
              int nxvh, int nyv, int nxhd, int nyhd);

void cetfield23(float complex dcu[], float complex exy[],
                float complex ffe[], float ci, float *wf, int nx,
                int ny, int nxvh, int nyv, int nxhd, int nyhd);

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
