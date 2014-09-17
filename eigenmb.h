#ifndef EIGENMB_H
#define EIGENMB_H
#include <math.h>
#include <float.h>

#define RC_COMPLEX_EVAL 2	/* code that complex eigenvalue obtained */

int      InvertMatrix (fit_double **a, int n, fit_double *col, int *indx, fit_double **a_inv);
int      LUDecompose (fit_double **a, int n, fit_double *vv, int *indx, fit_double *pd);
int      EigenRealGeneral (int n, fit_double **a, fit_double *v, fit_double *vi, fit_double **u, int *iwork, fit_double *work);
void     LUBackSubst (fit_double **a, int n, int *indx, fit_double *b);
int      EigenRG (int n, fit_double **a, fit_double *wr, fit_double *wi, fit_double **z, int *iv1, fit_double *fv1);
void     Balanc (int n, fit_double **a, int *pLow, int *pHigh, fit_double *scale);
void     Exchange (int j, int k, int l, int m, int n, fit_double **a, fit_double *scale);
void     ElmHes (int n, int low, int high, fit_double **a, int *intchg);
void     ElTran (int n, int low, int high, fit_double **a, int *intchg, fit_double **z);
int      Hqr2 (int n, int low, int high, fit_double **h, fit_double *wr, fit_double *wi, fit_double **z);
void     BalBak (int n, int low, int high, fit_double *scale, int m, fit_double **z);
void     CDiv (fit_double ar, fit_double ai, fit_double br, fit_double bi, fit_double *cr, fit_double *ci);
fit_double   D_sign (fit_double a, fit_double b);


#define TINY		1.0e-20
#if !defined(MAX)
#	define MAX(a,b)	(((a) > (b)) ? (a) : (b))
#endif
#if !defined(MIN)
#	define MIN(a,b)	(((a) < (b)) ? (a) : (b))
#endif

#define FALSE					0
#define TRUE					1
#define NO_ERROR				0
#define ERROR					1

#endif
