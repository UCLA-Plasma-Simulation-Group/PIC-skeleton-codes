/* OpenMP utility library */
/* Wrappers for calling the Fortran routines from a C main program */

void init_omp_(int *nth);

void setnthsize_(int *nth);

int getnthsize_();

/* Interfaces to C */

void cinit_omp(int nth) {
   init_omp_(&nth);
   return;
}

void csetnthsize(int nth) {
   setnthsize_(&nth);
   return;
}

int cgetnthsize() {
   return getnthsize_();
}
