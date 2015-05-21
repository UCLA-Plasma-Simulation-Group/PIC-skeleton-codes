void mpadd(float a[], float b[], float c[], int nx);

void init_mpi(int *idproc, int *nvp, int *irc, int argc, char *argv[]);

void vscatter(float f[],float g[], int nx, int nxp);

void vgather(float f[],float g[], int nx, int nxp);

void end_mpi();
