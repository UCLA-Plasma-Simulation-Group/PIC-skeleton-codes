/* C header file for CUDA Library for GPU Tutorial */
/* written by Viktor K. Decyk, UCLA */

void setgbsize(int nblock);

int getmmcc();

void gpu_fallocate(float **g_f, int nsize, int *irc);

void gpu_iallocate(int **g_i, int nsize, int *irc);

void gpu_deallocate(void *g_d, int *irc);

void gpu_fcopyin(float *f, float *g_f, int nsize);

void gpu_fcopyout(float *f, float *g_f, int nsize);

void emptykernel();

void init_cu(int dev, int *irc);

void end_cu();

void gpu_copy1(float *a, float *b, int mx, int nx);

void gpu_copy2a(float *a, float *b, int mx, int nx, int ny);

void gpu_copy2b(float *a, float *b, int mx, int nx, int ny);

void gpu_saxpy2(float *a, float *b, float s, int mx, int nx, int ny);

void gpu_copy3(float *a, float *b, int mx, int my, int nx, int ny);

void  gpu_transpose2(float *a, float *b, int mx, int nx, int ny);

void gpu_sum1(float *a, float *sa, int mx, int nx);

void gpu_sum2(float *a, float *d, int mx, int nx);

void gpu_sum3(float *a, float *d, float *sa, int mx, int nx);