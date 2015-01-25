/* C header file for CUDA utility Library */
/* written by Viktor K. Decyk, UCLA */

#include <complex.h>

void gpu_setgbsize(int nblock);

int getmmcc();

void gpu_fallocate(float **g_f, int nsize, int *irc);

void gpu_iallocate(int **g_i, int nsize, int *irc);

void gpu_callocate(float complex **g_c, int nsize, int *irc);

void gpu_deallocate(void *g_d, int *irc);

void gpu_fcopyin(float f[], float g_f[], int nsize);

void gpu_fcopyout(float f[], float g_f[], int nsize);

void gpu_icopyin(int f[], int g_f[], int nsize);

void gpu_icopyout(int f[], int g_f[], int nsize);

void gpu_ccopyin(float complex f[], float complex g_f[], int nsize);

void gpu_ccopyout(float complex f[], float complex g_f[], int nsize);

void gpu_zfmem(float g_f[], int nsize);

void gpu_set_cache_size(int nscache);

void emptykernel();

void init_cu(int dev, int *irc);

void end_cu();
