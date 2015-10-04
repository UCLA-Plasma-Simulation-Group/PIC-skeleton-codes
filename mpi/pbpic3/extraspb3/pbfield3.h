/* header file for pbfield3.c */

void cppotp32(float complex q[], float complex pot[],
              float complex ffc[], float *we, int nx, int ny, int nz,
              int kstrt, int nvpy, int nvpz, int nzv, int kxyp,
              int kyzp, int nzhd);

void cppdivf32(float complex f[], float complex df[], int nx, int ny,
               int nz, int kstrt, int nvpy, int nvpz, int nzv, int kxyp,
               int kyzp);

void cppgradf32(float complex df[], float complex f[], int nx, int ny,
                int nz, int kstrt, int nvpy, int nvpz, int nzv,
                int kxyp, int kyzp);

void cppcurlf32(float complex f[], float complex g[], int nx, int ny,
                int nz, int kstrt, int nvpy, int nvpz, int nzv,
                int kxyp, int kyzp);

void cppavpot332(float complex bxyz[], float complex axyz[], int nx,
                 int ny, int nz, int kstrt, int nvpy, int nvpz, int nzv,
                 int kxyp, int kyzp);

void cppavrpot332(float complex axyz[], float complex bxyz[],
                  float complex ffc[], float affp, float ci, int nx,
                  int ny, int nz, int kstrt, int nvpy, int nvpz,
                  int nzv, int kxyp, int kyzp, int nzhd);

void cppsmooth32(float complex q[], float complex qs[],
                 float complex ffc[], int nx, int ny, int nz, int kstrt,
                 int nvpy, int nvpz, int nzv, int kxyp, int kyzp,
                 int nzhd);

void cppsmooth332(float complex cu[], float complex cus[],
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
