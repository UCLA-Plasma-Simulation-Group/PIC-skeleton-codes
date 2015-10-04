/* header file for pfield3.c */

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

void cppsmooth32(float complex q[], float complex qs[],
                 float complex ffc[], int nx, int ny, int nz, int kstrt,
                 int nvpy, int nvpz, int nzv, int kxyp, int kyzp,
                 int nzhd);

void cpprdmodes32(float complex pot[], float complex pott[], int nx,
                  int ny, int nz, int modesx, int modesy, int modesz,    
                  int kstrt, int nvpy, int nvpz, int nzv, int kxyp,
                  int kyzp, int modesxpd, int modesypd, int modeszd);

void cppwrmodes32(float complex pot[], float complex pott[], int nx,
                  int ny, int nz, int modesx, int modesy, int modesz,
                  int kstrt, int nvpy, int nvpz, int nzv, int kxyp,
                  int kyzp, int modesxpd, int modesypd, int modeszd);
