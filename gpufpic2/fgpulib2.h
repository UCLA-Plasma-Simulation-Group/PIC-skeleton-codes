/* C header file for CUDA Fortran utility Library */
/* written by Viktor K. Decyk, UCLA */

void fgpu_setgbsize(int nblock);

int fgetmmcc();

void fgpu_zfmem(float g_f[], int nsize);

void fgpu_set_cache_size(int nscache);

void init_cuf(int dev, int *irc);

void end_cuf();

void fgpu_fcopyin(float f[], float g_f[], int nsize);

void fgpu_fcopyout(float f[], float g_f[], int nsize);

void fgpu_icopyin(int f[], int g_f[], int nsize);

void fgpu_icopyout(int f[], int g_f[], int nsize);

void fgpu_ccopyin(float complex f[], float complex g_f[], int nsize);

void fgpu_ccopyout(float complex f[], float complex g_f[], int nsize);
