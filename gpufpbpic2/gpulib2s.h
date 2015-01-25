/* C header file for CUDA special utility Library */
/* written by Viktor K. Decyk, UCLA */

#include <complex.h>

void gpu_fallocate(float **g_f, int nsize, int *irc);

void gpu_iallocate(int **g_i, int nsize, int *irc);

void gpu_callocate(float complex **g_c, int nsize, int *irc);

void gpu_deallocate(void *g_d, int *irc);

void hpl_fallocate(float **h_f, int nsize, int *irc);

void hpl_callocate(float complex **h_c, int nsize, int *irc);

void hpl_deallocate(void *h_d, int *irc);
