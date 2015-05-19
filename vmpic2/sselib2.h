/* header file for sselib2.c */

void sse_fallocate(float **s_f, int nsize, int *irc);

void sse_callocate(float complex **s_c, int nsize, int *irc);

void sse_iallocate(int **s_i, int nsize, int *irc);

void sse_deallocate(void *s_d);

void csse2iscan2(int *isdata, int nths);

int check_sse2();

int check_avx();

