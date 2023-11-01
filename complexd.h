

#ifndef dcomplex
typedef struct _dcomplexStruct { /* double-precision complex number */
        double r,i;
}  dcomplex;
#endif


/* double complex */
dcomplex dcadd (dcomplex a, dcomplex b);
dcomplex dcsub (dcomplex a, dcomplex b);
dcomplex dcmul (dcomplex a, dcomplex b);
dcomplex dcdiv (dcomplex a, dcomplex b);
double drcabs (dcomplex z);
dcomplex dcmplx (double re, double im);
dcomplex dconjg (dcomplex z);
dcomplex dcneg (dcomplex z);
dcomplex dcinv (dcomplex z);
dcomplex dcsqrt (dcomplex z);
dcomplex dcexp (dcomplex z);
dcomplex dcrmul (dcomplex a, double x);





