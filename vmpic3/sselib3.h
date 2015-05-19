/* header file for sselib3.c */

void sse_fallocate(float **s_f, int nsize, int *irc);

void sse_callocate(float complex **s_c, int nsize, int *irc);

void sse_iallocate(int **s_i, int nsize, int *irc);

void sse_deallocate(void *s_d);

int check_sse2();

int check_avx();
