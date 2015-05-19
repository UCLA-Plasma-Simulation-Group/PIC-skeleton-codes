/* header file for pplib2.c */

void cppinit2(int *idproc, int *nvp, int argc, char *argv[]);

void cppfndgrp(int locl[], int kstrt, int nvp, int *idev, int *ndev);

void cppexit();

void cppabort();

void cpwtimera(int icntrl, float *time, double *dtime);

void cppsum(float f[], float g[], int nxp);

void cppdsum(double f[], double g[], int nxp);

void cppimax(int f[], int g[], int nxp);

void cpppcncguard2l(float complex scs[], float complex scr[], int kstrt,
                    int nvp, int nxvh);

void cpppcnaguard2l(float complex scs[], float complex scr[], int kstrt,
                    int nvp, int nxvh);

void cppptpose(float complex sm[], float complex tm[], int nx, int ny,
               int kxp, int kyp, int kstrt, int nvp);

void cppptposen(float complex sm[], float complex tm[], int nx, int ny,
                int kxp, int kyp, int kstrt, int nvp, int ndim);

void cacsndrec(float complex stm[], int idproc, int nsize, int ntag,
               int mode);

void cpppmove2(float sbufr[], float sbufl[], float rbufr[], 
               float rbufl[], int ncll[], int nclr[], int mcll[],
               int mclr[], int kstrt, int nvp, int idimp, int nbmax,
               int mx1);
