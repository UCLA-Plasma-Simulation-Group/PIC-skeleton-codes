/* header file for pbfield2.c */

void cppotp2(float complex q[], float complex pot[],
             float complex ffc[], float *we, int nx, int ny, int kstrt,
             int nyv, int kxp, int nyhd);

void cppdivf2(float complex f[], float complex df[], int nx, int ny,
              int kstrt, int ndim, int nyv, int kxp);

void cppgradf2(float complex df[], float complex f[], int nx, int ny,
               int kstrt, int ndim, int nyv, int kxp);

void cppcurlf2(float complex f[], float complex g[], int nx, int ny,
               int kstrt, int nyv, int kxp);

void cppavpot23(float complex bxy[], float complex axy[], int nx,
                int ny, int kstrt, int nyv, int kxp);

void cppavrpot23(float complex axy[], float complex bxy[],
                 float complex ffc[], int affp, int ci, int nx, int ny,
                 int kstrt, int nyv, int kxp, int nyhd);

void cppsmooth2(float complex q[], float complex qs[],
                float complex ffc[], int nx, int ny, int kstrt, int nyv,
                int kxp, int nyhd);

void cppsmooth23(float complex cu[], float complex cus[],
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
