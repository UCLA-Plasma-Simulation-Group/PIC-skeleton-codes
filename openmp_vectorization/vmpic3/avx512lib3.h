/* header file for avx512lib3.c */

void avx512_fallocate(float **s_f, int nsize, int *irc);

void avx512_callocate(float complex **s_c, int nsize, int *irc);

void avx512_iallocate(int **s_i, int nsize, int *irc);

void avx512_deallocate(void *s_d);

void cknciscan2(int *isdata, int nths);
