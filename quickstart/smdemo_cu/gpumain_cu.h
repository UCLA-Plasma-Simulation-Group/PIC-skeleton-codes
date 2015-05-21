void setgbsize(int nblock);

void g_fallocate(float **g_f, int nsize, int *irc);

void g_deallocate(float **g_f, int *irc);

void copyin_gmemptr(float *f, float *g_f, int nsize);

void copyout_gmemptr(float *f, float *g_f, int nsize);

void init_cu(int dev, int *irc);

void gpadd(float *a, float *b, float *c, int nx);

void end_cu();

    
