/* header file for mpdfield2.c */

void cmppotp2(float complex q[], float complex pot[],
              float complex ffc[], float *we, int nx, int ny, int kstrt,
              int nyv, int kxp, int nyhd);

void cmppdivf2(float complex f[], float complex df[], int nx, int ny,
               int kstrt, int ndim, int nyv, int kxp);

void cmppgradf2(float complex df[], float complex f[], int nx, int ny,
                int kstrt, int ndim, int nyv, int kxp);

void cmppcurlf2(float complex f[], float complex g[], int nx, int ny,
                int kstrt, int nyv, int kxp);

void cmppapotp23(float complex cu[], float complex axy[],
                 float complex ffc[], float ci, float *wm, int nx,
                 int ny, int kstrt, int nyv, int kxp, int nyhd);

void cmppsmooth2(float complex q[], float complex qs[],
                 float complex ffc[], int nx, int ny, int kstrt,
                 int nyv, int kxp, int nyhd);

void cmppsmooth23(float complex cu[], float complex cus[],
                  float complex ffc[], int nx, int ny, int kstrt,
                  int nyv, int kxp, int nyhd);

void cpprdmodes2(float complex pot[], float complex pott[], int nx,
                 int ny, int modesx, int modesy, int kstrt, int nyv,
                 int kxp, int modesxpd, int modesyd);

void cppwrmodes2(float complex pot[], float complex pott[], int nx,
                 int ny, int modesx, int modesy, int kstrt, int nyv,
                 int kxp, int modesxpd, int modesyd);

void cpprdvmodes2(float complex vpot[], float complex vpott[], int nx,
                  int ny, int modesx, int modesy, int ndim, int kstrt, 
                  int nyv, int kxp, int modesxpd, int modesyd);

void cppwrvmodes2(float complex vpot[], float complex vpott[], int nx,
                  int ny, int modesx, int modesy, int ndim, int kstrt, 
                  int nyv, int kxp, int modesxpd, int modesyd);
