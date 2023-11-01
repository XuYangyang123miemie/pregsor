#ifndef COMPLEX_H_
#define COMPLEX_H_
typedef struct _complexStruct {
	double r;
	double i;
} complex;

/* complex a plus complex b and return a complex structure c=a+b*/
complex cadd (complex a, complex b);

/* complex a sub complex b and return a complex structure c=a-b*/
complex csub (complex a, complex b);

/* complex a multiply complex b and return a complex structure c=a*b*/
complex cmul (complex a, complex b);

/* complex a divide complex b and return a complex structure c=a/b*/
complex cdiv (complex a, complex b);

/* *************************************************
cmplx	make a complex number from two real numbers
conjg	complex conjugate of a complex number 
cneg	negate a complex number
cinv	invert a complex number
csqrt	complex square root of a complex number
cexp	complex exponential of a complex number
crmul	multiply a complex number by a real number 
rcabs	real magnitude of a complex number(
************************************************* */

complex cmplx (double re, double im);

complex conjg (complex z);

complex cneg (complex z);

complex cinv (complex z);

complex ccsqrt (complex z);

complex ccexp (complex z);

complex crmul (complex a, double x);

double rcabs (complex z);
#endif




