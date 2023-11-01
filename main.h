#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <errno.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <limits.h>
#include <time.h>
#include"complex.c"
#include"complexd.c"


#include "cbesselxuyy.h"

#define NINT(x) ((int)((x)>0.0?(x)+0.5:(x)-0.5))
#define	MAX(x,y) ((x) > (y) ? (x) : (y))
#define	MIN(x,y) ((x) < (y) ? (x) : (y))
#define PI (3.141592653589793)
#ifndef ABS
#define ABS(x) ((x) < 0 ? -(x) : (x))
#endif
#define EPS 1.0e-6
#define EPS2 1.0e-3
#define INITIAL_T 99999

#define SU_NFLTS	65535   /* Arbitrary limit on data array size*/
#define FSIZE sizeof(double)

#undef off_t
#define off_t long long

#define LOOKFAC	2	/* Look ahead factor for npfaro	  */
#define PFA_MAX	720720	/* Largest allowed ntfft	          */


/*ricker wavelet*/
void ricker1_wavelet (int nt, double dt, double fpeak, double *wavelet);
double *ricker(double Freq,double dt,int *Npoint);
/*hanker*/
double Bessel_first_zero ( double x );
double Bessel_first_first ( double x );
double Hanker_first_zero_imag(double x);
double Hanker_first_zero_real(double x);

/*Green*/
double Greenr(double r);
double Greeni(double r);
double GreenCWr(dcomplex r);
double GreenCWi(dcomplex r);

/*compute the limit of integration*/
void ss2(double fxr,double fzr,double h,int vi,int vj,int j,int n, double lim[]);
void ss4(double fxr,double fzr,double h,int vi,int vj,int vi1,int vj1,int j,int n, double lim[]);
void ss2_CW(double fxr,double fzr,double h,int vi,int vj,int j,int n, double lim[]);
void ss4_CW(double fxr,double fzr,double h,int vi,int vj,int vi1,int vj1,int j,int n, double lim[]);
/* integrand function*/
double fre2(int n,double x,double z,double kk,double var[]);
double fim2(int n,double x,double z,double kk,double var[]);
double fre4(int n,double kk,double var[]);
double fim4(int n,double kk,double var[]);
double fre2_CW(int n,double x,double z,dcomplex kk,double var[]);
double fim2_CW(int n,double x,double z,dcomplex kk,double var[]);
double fre4_CW(int n,dcomplex kk,double var[]);
double fim4_CW(int n,dcomplex kk,double var[]);
/*Gauss numerical integration*/
double gaus2(double x,double z,double fxr,double fzr,double h,int vi,int vj, double kk,int n,int js[],void (*ss)(),double (*f)());
double gaus4(double fxr,double fzr,double h,int vi,int vj,int vi1,int vj1, double kk,int n,int js[],void (*ss)(),double (*f)());
/*Complex wavenumber Gauss numerical integration*/
double gaus2_CW(double x,double z,double fxr,double fzr,double h,int vi,int vj, dcomplex kk,int n,int js[],void (*ss)(),double (*f)());
double gaus4_CW(double fxr,double fzr,double h,int vi,int vj,int vi1,int vj1, dcomplex kk,int n,int js[],void (*ss)(),double (*f)());

/*generate a real number diag matrix*/
void Rdiag(double *a,int r,double **matri);

/*2dim matrix A*B*/
void matmul_0(double *a,int r1,int c1,double *b,double *c);
void matmul_1(double **a,int r1,int c1,double *b,double *c);
void matmul_2(double **a,int r1,int c1,double **b,int c2,double **c);

/*matrix transp*/
complex *mattransp(complex *p,int m,int n);
dcomplex *mattransp_d(dcomplex *p,int m,int n);

/*Block Toeplitz matrix G(r,r') mutiply vector y.___FFT1*/
void BTmatrixmutip(complex **G,complex **y,int Nx,int Nz,complex **multip);
void BTmatrixmutip_d(dcomplex **G,dcomplex **y,int Nx,int Nz,dcomplex **multip);

/*a Nx*Nz order BTTB matrix G multiply vector y.__FFT2*/
void BTTB_y(complex **GT,complex **y,int Nx,int Nz,complex **multip_GT);
void BTTB_y_d(dcomplex **GT,dcomplex **y,int Nx,int Nz,dcomplex **multip_GT);

/*BCCB matrix GC multiply vector y,or the inverse of preconditon BCCB matrix multiply vector.___FFT2*/
/* GC*y ,or M(-1)*b */
void BCCB_y(complex **GC,complex **yC,int Nx,int Nz,complex **multip_GC,int isign);
void BCCB_y_d(dcomplex **GC,dcomplex **yC,int Nx,int Nz,dcomplex **multip_GC,int isign);

/*make a Strang BCCB precondition matrix from a BTTB matrix*/
void BCCB_S(complex **GT,int Nx,int Nz,complex **Pre_S);
void BCCB_S_d(dcomplex **GT,int Nx,int Nz,dcomplex **Pre_S);

/*make a T.Chen BCCB precondition matrix from a BTTB matrix*/
void BCCB_TC(complex **GT,int Nx,int Nz,complex **Pre_C);
void BCCB_TC_d(dcomplex **GT,int Nx,int Nz,dcomplex **Pre_C);

/*matrix interpolate*/
void interp2d(float *original_data, float *processed_data,int W,int H);

/*conjugate gradient method*/
/*void sprsin(double **a,int n,double thresh,int nmax,double sa[],int ija[]);
void asolve(int n,double b[],double x[],double sa[],int ija[],int itrnsp);
void atimes(int n,double x[],double r[],double sa[],int ija[],int itrnsp);
double snrm(int n,double sx[],int itol);
void linbcg(int n,double b[],double x[],double sa[],int ija[],int itol,double tol,int itmax,int *iter,double *err);*/

/*complex coefficients conjugate gradient method*/
void asolve(int nx,int nz,complex **b,complex **x,complex **M);
double snrm(int nx,int nz,complex **sx,int itol);
void linbcg(double **diag,complex **G,complex **M,int nx,int nz,complex **b,complex **x,
             int itol,double tol,int itmax,int *iter,double *err);
void linbcg_CW(complex **diag,complex **G,complex **M,int nx,int nz,complex **b,complex **x,
             int itol,double tol,int itmax,int *iter,double *err);
void linbcgstab(double **diag,complex **G,complex **M,int nx,int nz,complex **b,complex **x,
             int itol,double tol,int itmax,int *iter,double *err);
void linbcg_psc(double **diag,complex **G,complex **M,int nx,int nz,complex **b,complex **x,
             int itol,double tol,int itmax,int *iter,double *err);
void linbcgstab_psc(double **diag,complex **G,complex **M,int nx,int nz,complex **b,complex **x,
             int itol,double tol,int itmax,int *iter,double *err);
void linbcgstab_sor_psc(double **diag,complex **G,complex **M,int nx,int nz,complex **b,complex **x,
             int itol,double tol,int itmax,int *iter,double *err);
			 
/*GSOR*/			 
void linssor(complex **g,double **diag,complex **s,int nx,int nz,complex **u0,int itmax);
void linssor_d(dcomplex **g,double **diag,dcomplex **s,int nx,int nz,dcomplex **u0,int itmax);

void linssor_CW_PreO(complex **g,complex **diag,complex **M,complex **s,int nx,int nz,complex **u0,int itmax);
void linssor_CW_PreO_d(dcomplex **g,dcomplex **diag,dcomplex **M,dcomplex **s,int nx,int nz,dcomplex **u0,int itmax);

void linssor_CW(complex **g,complex **diag,complex **s,int nx,int nz,complex **u0,int itmax);
void linssor_preOCW(complex **g,complex **oor,complex **M,complex **s,int nx,int nz,complex **u0,int itmax);
//void linssor_CW_PreO(complex **g,complex **diag,complex **M,complex **s,int nx,int nz,complex **u0,int itmax);
void linssor_CW_PreM(complex **g,complex **diag,complex **M,complex **s,int nx,int nz,complex **u0,int itmax);
void linssor_CW_PrePL(complex **g,complex **diag,complex **M,complex **N,complex **s,int nx,int nz,complex **u0,int itmax);
void linssor_CW_PreP(complex **g,complex **diag,complex **s,int nx,int nz,complex **u0,int itmax);