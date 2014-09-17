/*

FitModeL :  a program that Fits substitution models to the data
using Maximum Likelihood.  

Copyright (C) Stephane Guindon. jul. 2004 onward

All parts of  the source except where indicated  are distributed under
the GNU public licence.  See http://www.opensource.org for details.

*/


#ifndef LK_H
#define LK_H

void Init_Tips_At_One_Site_Nucleotides(char state,fit_double **p_lk);
void Init_Tips_At_One_Site_AA(char aa,fit_double **p_lk);
void Get_All_Partial_Lk(arbre *tree,edge *b_fcus,node *a,node *d);
void Get_All_Partial_Lk_Scale(arbre *tree,edge *b_fcus,node *a,node *d);
void Pre_Order_Lk(node *pere,node *fils,arbre *tree);
void Post_Order_Lk(node *pere,node *fils,arbre *tree);
void Lk(arbre *tree);
fit_double Lk_At_Given_Edge(arbre *tree,edge *b_fcus);
void Site_Lk(arbre *tree,allseq *alldata);
void Site_Lk_Codon(arbre *tree,allseq *alldata);
fit_double Lk_At_Given_Edge_Codon(arbre *tree,edge *b_fcus);
void Site_Lk_Nt_AA(arbre *tree,allseq *alldata);
fit_double Lk_At_Given_Edge_Nt_AA(arbre *tree,edge *b_fcus);
void Update_P(arbre *tree,int t_edge_num);
fit_double Return_Lk(arbre *tree);
fit_double Return_Abs_Lk(arbre *tree);
fit_double ****Get_Partial_Lk_Struct(arbre *tree,int len,int n_catg,int n_catq);
void Update_P_Lk(arbre *tree,edge *b_fcus,node *n);
void Make_Tree_4_Lk(arbre *tree,allseq *alldata,int n_site);
void Get_Infered_States(edge *b_fcus,node *a,node *d,arbre *tree,int site_num);
void Get_PMat(arbre *tree,fit_double len,fit_double ****Pij,qmat *qmat_struct);
void Get_All_PMat_Post(node *a,node *d,edge *b,arbre *tree);
void Proba_Omega_On_One_Edge(edge *b,qmat *qmat_struct,arbre *tree);
fit_double Unscaled_Dwell_S1(int beg,int mid,int end,edge *b,arbre *tree);
void Proba_Omega_At_One_Site(arbre *tree);
void Init_Tips(arbre *tree);
fit_double Unscaled_Dwell_S2(int beg, int mid, int end, edge *b, arbre *tree);
void Compute_Pmat_S2(fit_double l, fit_double **Pij,qmat *qmat_struct, arbre *tree);
void Eigen_S2(qmat *qmat_struct, arbre *tree);
void Compute_All_Dwell_Probs_S2(edge *b_fcus,qmat *qmat_struct,arbre *tree);
void Compute_Proba_Omega_On_Edges(arbre *tree);
void Compute_Proba_Omega_On_Every_Edge_At_One_Site(arbre *tree);
void Init_Tips_At_One_Site_Codon(char state1, char state2, char state3, fit_double **p_lk, code *c_code, int ns, int *stop);

fit_double Lk_Arg(fit_double *param_val, 
              arbre *tree);
void Integral_Term_On_One_Edge(edge *b, arbre *tree);
void Compute_Proba_Omega_On_One_Edge_At_One_Site(edge *b, fit_double *p, arbre *tree);
fit_double  Select_Regime_Change_Proba_S1(int beg, int end, edge *b, fit_double length, arbre *tree);
fit_double Compute_Proba_Omega_At_One_Spot_At_One_Site(edge *b, fit_double pos, int wclass, int site, arbre *tree);

#endif






