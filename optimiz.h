/*

FitModeL :  a program that Fits substitution models to the data
using Maximum Likelihood.  

Copyright (C) Stephane Guindon. jul. 2004 onward

All parts of  the source except where indicated  are distributed under
the GNU public licence.  See http://www.opensource.org for details.

*/


#ifndef OPTIMIZ_H
#define OPTIMIZ_H

fit_double Br_Len_Golden(fit_double ax,fit_double bx,fit_double cx,fit_double tol,fit_double *xmin,edge *b_fcus,arbre *tree);
int Generic_Brak(fit_double *param,fit_double *ax,fit_double *bx,fit_double *cx,fit_double *fa,fit_double *fb,fit_double *fc,fit_double lim_inf,fit_double lim_sup,arbre *tree);
fit_double Generic_Golden(fit_double *param,fit_double ax,fit_double bx,fit_double cx,fit_double tol,fit_double lim_inf,fit_double lim_sup,fit_double *xmin,arbre *tree,int n_iter_max);
fit_double Generic_Brent(fit_double *param, 
			 fit_double ax, fit_double bx, fit_double cx, fit_double tol, 
			 fit_double lim_inf,
			 fit_double lim_sup,
			 arbre *tree, int n_iter_max);
fit_double RRparam_GTR_Golden(fit_double ax,fit_double bx,fit_double cx,fit_double tol,fit_double *xmin,arbre *tree,allseq *alldata,fit_double *param,int n_iter_max);
int Br_Len_Brak(fit_double *ax,fit_double *bx,fit_double *cx,fit_double *fa,fit_double *fb,fit_double *fc,edge *b_fcus,arbre *tree);
fit_double Br_Len_Brent(fit_double ax,fit_double bx,fit_double cx,fit_double tol,fit_double *xmin,edge *b_fcus,arbre *tree,int n_iter_max);
void Optimize_Single_Param_Generic(arbre *tree,fit_double *param,fit_double start,fit_double lim_inf,fit_double lim_sup,int n_max_iter);
void Round_Optimize(arbre *tree,allseq *data);
void Optimize_Br_Len_Serie(node *a,node *d,edge *b_fcus,arbre *tree,allseq *alldata);
void Optimize_Tpos_Serie(node *a,node *d,edge *b_fcus,arbre *tree,allseq *alldata,int n_passes);
void Print_Lk_Progress(arbre *tree,fit_double lk_new,fit_double lk_old,int n_iter);
void Optimiz_All_Free_Param(arbre *tree,int verbose);
void BFGS(arbre *tree, 
	  fit_double *p, 
	  int n, 
	  fit_double gtol, 
	  fit_double step_size,
	  fit_double(*func)(arbre *tree), 
	  int(*dfunc)(arbre *tree,fit_double *param,int n_param,fit_double stepsize,fit_double(*func)(arbre *tree),fit_double *derivatives), 
	  int(*lnsrch)(arbre *tree, int n, fit_double *xold, fit_double fold,fit_double *g, fit_double *p, fit_double *x,fit_double *f, fit_double stpmax, int *check),
	  int *failed);
void Lnsrch_GTR(arbre *tree,int n,fit_double *xold,fit_double fold,fit_double *g,fit_double *p,fit_double *x,fit_double *f,fit_double stpmax,int *check, fit_double (*func)());
void Lnsrch_Omega(arbre *tree,int n,fit_double *xold,fit_double fold,fit_double *g,fit_double *p,fit_double *x,fit_double *f,fit_double stpmax,int *check, fit_double (*func)());
void Lnsrch_Qmat_Proba(arbre *tree,int n,fit_double *xold,fit_double fold,fit_double *g,fit_double *p,fit_double *x,fit_double *f,fit_double stpmax,int *check, fit_double (*func)());
void Lnsrch_Omega_Proba(arbre *tree,int n,fit_double *xold,fit_double fold,fit_double *g,fit_double *p,fit_double *x,fit_double *f,fit_double stpmax,int *check, fit_double (*func)());
void Lnsrch_Nucleotide_Frequencies(arbre *tree,int n,fit_double *xold,fit_double fold,fit_double *g,fit_double *p,fit_double *x,fit_double *f,fit_double stpmax,int *check, fit_double (*func)());
void Lnsrch_Theta(arbre *tree,int n,fit_double *xold,fit_double fold,fit_double *g,fit_double *p,fit_double *x,fit_double *f,fit_double stpmax,int *check, fit_double (*func)());
void Lnsrch_Selec_Reg_Freq(arbre *tree, int n, fit_double *xold, fit_double fold, fit_double *g, fit_double *p, fit_double *x,
			   fit_double *f, fit_double stpmax, int *check, fit_double (*func)());
void Lnsrch_Trans_Omega_Proba(arbre *tree, int n, fit_double *xold, fit_double fold, fit_double *g, fit_double *p, fit_double *x,
                              fit_double *f, fit_double stpmax, int *check, fit_double (*func)());

void Lnsrch_Trans_Qmat_Proba(arbre *tree, int n, fit_double *xold, fit_double fold, fit_double *g, fit_double *p, fit_double *x,
                             fit_double *f, fit_double stpmax, int *check, fit_double (*func)());

int Lnsrch(arbre *tree, int n, fit_double *xold, fit_double fold, fit_double *g, fit_double *p, fit_double *x,
	   fit_double *f, fit_double stpmax, int *check);

fit_double Linmin_On_A_Direction(fit_double *x, 
                             fit_double *dir, 
                             fit_double n, 
                             fit_double *grad,
                             fit_double f0,
                             fit_double (*func)(),
                             arbre *tree);

fit_double Vect_Prod(fit_double *a, 
                 fit_double *b, 
                 int n);

void Conjugate_Gradients(fit_double *x, 
                         int n, 
                         arbre *tree);

#endif

