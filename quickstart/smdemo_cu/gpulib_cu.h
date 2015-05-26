void gpadd(float *a, float *b, float *c, int nx);

void init_cu(int dev, int *irc);

void end_cu();

void setgbsize(int nblock);

void gpu_fallocate(float **g_f, int nsize, int *irc);

void gpu_deallocate(float **g_f, int *irc);

void gpu_fcopyin(float *f, float *g_f, int nsize);

void gpu_fcopyout(float *f, float *g_f, int nsize);



    
