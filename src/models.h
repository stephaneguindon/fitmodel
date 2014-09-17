/*

FitModeL :  a program that Fits substitution models to the data
using Maximum Likelihood.  

Copyright (C) Stephane Guindon. jul. 2004 onward

All parts of  the source except where indicated  are distributed under
the GNU public licence.  See http://www.opensource.org for details.

*/
#ifndef MODELS_H
#define MODELS_H

int PMat_TN93(fit_double l,model *mod,fit_double **Pij);
void dPMat_TN93(fit_double l,fit_double **dPij,model *mod,fit_double rr);
void d2PMat_TN93(fit_double l,fit_double **d2Pij,model *mod,fit_double rr);
int Matinv(fit_double *x,int n,int m,fit_double *space);
int PMat(fit_double l,model *mod,fit_double **Pij,qmat *qmat_struct,arbre *tree);
void dPMat(fit_double l,fit_double rr,model *mod,fit_double **dPij);
void d2PMat(fit_double l,fit_double rr,model *mod,fit_double **d2Pij);
int Init_Qmat_Dayhoff(fit_double *daa,fit_double *pi);
int Init_Qmat_DCMut(fit_double *daa,fit_double *pi);
int Init_Qmat_WAG(fit_double *daa,fit_double *pi);
int Init_Qmat_JTT(fit_double *daa,fit_double *pi);
int Init_Qmat_MtREV(fit_double *daa,fit_double *pi);
model *Init_Model(allseq *data,option *input);
int PMat_Numeric(fit_double l,int ns,fit_double **Pij,qmat *qmat_struct,arbre *tree);
void Update_Qmat_GTR(model *mod);
void Set_Model_Parameters(arbre *tree);
void Update_Qmat_Codon(model *mod,qmat *qmat_struct);
void Update_Eigen(int dim,fit_double *Q,fit_double *U,fit_double *V,fit_double *Root, model *mod);
void Init_T_Omega(model *mod);
void Init_T_Omega_Proba(model *mod);

#endif
