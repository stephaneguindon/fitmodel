#ifndef EIGEN_H
#define EIGEN_H

typedef struct { fit_double re, im; } complex;
#define csize(a) ((fit_double)fabs(a.re)+(fit_double)fabs(a.im))


int eigenQREV (fit_double Q[], fit_double pi[], int n, 
               fit_double Root[], fit_double U[], fit_double V[], fit_double spacesqrtpi[]);
void HouseholderRealSym(fit_double a[], int n, fit_double d[], fit_double e[]);
int EigenTridagQLImplicit(fit_double d[], fit_double e[], int n, fit_double z[]);
void EigenSort(fit_double d[], fit_double U[], int n);
int eigenRealSym(fit_double A[], int n, fit_double Root[], fit_double work[]);



int eigen(int job,fit_double *A,int n,fit_double *rr,fit_double *ri,fit_double *vr,fit_double *vi,fit_double *work);
complex compl(fit_double re,fit_double im);
complex _conj(complex a);
complex cplus(complex a,complex b);
complex cminus(complex a,complex b);
complex cby(complex a,complex b);
complex cdiv(complex a,complex b);
complex my_cexp(complex a);
complex cfactor(complex x,fit_double a);
int cxtoy(complex *x,complex *y,int n);
int cmatby(complex *a,complex *b,complex *c,int n,int m,int k);
int cmatinv(complex *x,int n,int m,fit_double *space);
void balance(fit_double *mat,int n,int *low,int *hi,fit_double *scale);
void unbalance(int n,fit_double *vr,fit_double *vi,int low,int hi,fit_double *scale);
void elemhess(int job,fit_double *mat,int n,int low,int hi,fit_double *vr,fit_double *vi,int *work);
int realeig(int job,fit_double *mat,int n,int low,int hi,fit_double *valr,fit_double *vali,fit_double *vr,fit_double *vi);


#endif

