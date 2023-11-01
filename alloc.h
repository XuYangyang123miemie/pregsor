#ifndef ALLOC_H_
#define ALLOC_H_
void *alloc1(int n1,int size);
void *realloc1(void *v,int n1,int size);
void **alloc2(int n1,int n2,int size);
void ***alloc3(int n1,int n2,int n3,int size);
void free1 (void *p);
void free2 (void **p);
void free3 (void ***p);

int *alloc1int(int n1);
void free1int(int *p);
int *realloc1int(int *v, int n1);

int **alloc2int(int n1, int n2);
void free2int(int **p);

int ***alloc3int(int n1, int n2, int n3);
void free3int(int ***p);

float *alloc1float(int n1);
float *realloc1float(float *v, int n1);
void free1float(float *p);

float **alloc2float(int n1, int n2);
void free2float(float **p);

float ***alloc3float(int n1, int n2, int n3);
void free3float(float ***p);

double *alloc1double(int n1);
double *realloc1double(double *v, int n1);
void free1double(double *p);

double **alloc2double(int n1, int n2);
void free2double(double **p);

double ***alloc3double(int n1, int n2, int n3);
void free3double(double ***p);

complex *alloc1complex(int n1);

void free1complex(complex *p);

complex **alloc2complex(int n1, int n2);

void free2complex(complex **p);

complex ***alloc3complex(int n1, int n2, int n3);

void free3complex(complex ***p);

dcomplex *alloc1dcomplex (size_t n1);
dcomplex **alloc2dcomplex (size_t n1, size_t n2);
dcomplex ***alloc3dcomplex (size_t n1, size_t n2, size_t n3);

void free1dcomplex (dcomplex *p);
void free2dcomplex (dcomplex **p);
void free3dcomplex (dcomplex ***p);

#endif

