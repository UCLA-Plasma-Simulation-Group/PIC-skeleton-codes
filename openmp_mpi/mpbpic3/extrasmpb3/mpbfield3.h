/* header file for mpbfield3.c */

void cmppotp32(float complex q[], float complex pot[],
               float complex ffc[], float *we, int nx, int ny, int nz,
               int kstrt, int nvpy, int nvpz, int nzv, int kxyp,
               int kyzp, int nzhd);

void cmppdivf32(float complex f[], float complex df[], int nx, int ny,
                int nz, int kstrt, int nvpy, int nvpz, int nzv,
                int kxyp, int kyzp);

void cmppgradf32(float complex df[], float complex f[], int nx, int ny,
                 int nz, int kstrt, int nvpy, int nvpz, int nzv,
                 int kxyp, int kyzp);

void cmppcurlf32(float complex f[], float complex g[], int nx, int ny,
                 int nz, int kstrt, int nvpy, int nvpz, int nzv,
                 int kxyp, int kyzp);

void cmppavpot332(float complex bxyz[], float complex axyz[], int nx,
                  int ny, int nz, int kstrt, int nvpy, int nvpz,
                  int nzv, int kxyp, int kyzp);

void cmppavrpot332(float complex axyz[], float complex bxyz[],
                   float complex ffc[], float affp, float ci, int nx,
                   int ny, int nz, int kstrt, int nvpy, int nvpz,
                   int nzv, int kxyp, int kyzp, int nzhd);

void cmppsmooth32(float complex q[], float complex qs[],
                  float complex ffc[], int nx, int ny, int nz,
                  int kstrt, int nvpy, int nvpz, int nzv, int kxyp,
                  int kyzp, int nzhd);

void cmppsmooth332(float complex cu[], float complex cus[],
                   float complex ffc[], int nx, int ny, int nz,
                   int kstrt, int nvpy, int nvpz, int nzv, int kxyp,
                   int kyzp, int nzhd);

void cpprdmodes32(float complex pot[], float complex pott[], int nx,
                  int ny, int nz, int modesx, int modesy, int modesz,    
                  int kstrt, int nvpy, int nvpz, int nzv, int kxyp,
                  int kyzp, int modesxpd, int modesypd, int modeszd);

void cppwrmodes32(float complex pot[], float complex pott[], int nx,
                  int ny, int nz, int modesx, int modesy, int modesz,
                  int kstrt, int nvpy, int nvpz, int nzv, int kxyp,
                  int kyzp, int modesxpd, int modesypd, int modeszd);

void cpprdvmodes32(float complex vpot[], float complex vpott[], int nx,
                   int ny, int nz, int modesx, int modesy, int modesz, 
                   int ndim, int kstrt, int nvpy, int nvpz, int nzv,
                   int kxyp, int kyzp, int modesxpd, int modesypd, 
                   int modeszd);

void cppwrvmodes32(float complex vpot[], float complex vpott[], int nx,
                   int ny, int nz, int modesx, int modesy, int modesz, 
                   int ndim, int kstrt, int nvpy, int nvpz, int nzv,
                   int kxyp, int kyzp, int modesxpd, int modesypd,
                   int modeszd);
