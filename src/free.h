/*

FitModeL :  a program that Fits substitution models to the data
using Maximum Likelihood.  

Copyright (C) Stephane Guindon. jul. 2004 onward

All parts of  the source except where indicated  are distributed under
the GNU public licence.  See http://www.opensource.org for details.

*/


#ifndef FREE_H
#define FREE_H

void Free_All_Nodes_Light(arbre *tree);
void Free_All_Edges_Light(arbre *tree);
void Free_Partial_Lk(fit_double ****p_lk,int len,int n_catg);
void Free_Tree(arbre *tree);
void Free_Edge(edge *b);
void Free_Node(node *n);
void Free_Cseq(allseq *data);
void Free_Seq(seq **d,int n_otu);
void Free_All(seq **d,allseq *alldata,arbre *tree);
void Free_SubTree(edge *b_fcus,node *a,node *d,arbre *tree);
void Free_Tree_Ins_Tar(arbre *tree);
void Free_Tree_Lk(arbre *tree);
void Free_dPij(arbre *tree);
void Free_Edge_P_Lk_Struct(edge *b,arbre *tree);
void Free_Node_Lk(node *n);
void Free_Edge_Lk(arbre *tree,edge *b);
void Free_Model(model *mod);
void Free(void *p);
void Free_Input(option *input);
void Free_Code(code *c_code);
void Free_Qmat(qmat *this1,model *mod);
void Free_NNI(nni *this1);
void Free_Qmat_Struct(model *mod);
void Free_Numerical_Param(numpar *this_one);
void Free_Integral_Term_On_One_Edge(edge *b, arbre *tree);
void Free_Edge_Labels(edge *b);

#endif
