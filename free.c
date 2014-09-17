/*

FitModeL :  a program that Fits substitution models to the data
using Maximum Likelihood.  

Copyright (C) Stephane Guindon. jul. 2004 onward

All parts of  the source except where indicated  are distributed under
the GNU public licence.  See http://www.opensource.org for details.

*/


#include "config.h"
#include "utilities.h"
#include "free.h"

/*********************************************************/

void Free_All_Nodes_Light(arbre *tree)
{
  int i;
  For(i,2*tree->n_otu-2) 
    Free_Node(tree->noeud[i]);
}

/*********************************************************/

void Free_All_Edges_Light(arbre *tree)
{
  int i;
  For(i,2*tree->n_otu-3) 
    Free_Edge(tree->t_edges[i]);
}


/*********************************************************/

void Free_Partial_Lk(fit_double ****p_lk, int len, int n_catg)
{
  int i,j;

  For(i,len)
    {
      For(j,n_catg) Free((*p_lk)[i][j]);
      Free((*p_lk)[i]);
    }
  Free((*p_lk));
  (*p_lk) = NULL;
}

/*********************************************************/

void Free_Tree(arbre *tree)
{
    int i,j,k;
    edge *b;
    node *n;
    

/*     if(tree->n_root) */
/*       { */
/* 	Free_Node(tree->noeud[2*tree->n_otu-1]); */
/* 	Free_Node(tree->noeud[2*tree->n_otu-2]); */
/* 	Free_Edge(tree->t_edges[2*tree->n_otu-2]); */
/* 	Free_Edge(tree->t_edges[2*tree->n_otu-3]); */
/*       } */


    if(tree->has_bip)
      {
	For(i,2*tree->n_otu-2)
	  {
	    Free(tree->noeud[i]->bip_size);
	    For(j,3)
	      {
		Free(tree->noeud[i]->bip_node[j]);
		For(k,tree->n_otu) Free(tree->noeud[i]->bip_name[j][k]);
		Free(tree->noeud[i]->bip_name[j]);
	      }
	    Free(tree->noeud[i]->bip_node);
	    Free(tree->noeud[i]->bip_name);
	  }
      }
    
    For(i,2*tree->n_otu-2)
      {
	b = tree->t_edges[i];
	Free_Edge(b);
      }
    Free(tree->t_edges);
    
    
    For(i,2*tree->n_otu-2)
      {
	n = tree->noeud[i];
	Free_Node_Lk(n);
      }


    For(i,2*tree->n_otu-1)
      {
	n = tree->noeud[i];
	Free_Node(n);
      }

    Free(tree->noeud);
    
    
    Free(tree);
}

/*********************************************************/

void Free_Edge(edge *b)
{
  Free_Edge_Labels(b);
/*   Free_NNI(b->s_nni); */
  Free(b);
}

/*********************************************************/

void Free_Node(node *n)
{
  Free(n->b);
  Free(n->v);
  Free(n->l);
  Free(n->l_prime);
  Free(n->name);
  Free(n->tpos_min);
  Free(n->post_w_ave);
  Free(n);
}

/*********************************************************/

void Free_Edge_Labels(edge *b)
{
  int i;
  For(i,b->n_labels+b->n_labels%BLOCK_LABELS) Free(b->labels[i]);
  Free(b->labels);
  b->labels = NULL;
}

/*********************************************************/

void Free_Cseq(allseq *data)
{
  int i;
  
  Free(data->selclass);
  Free(data->invar);
  Free(data->wght);
  Free(data->ambigu);
  Free(data->b_frq);
  For(i,data->init_len) Free(data->pospatt[i]);
  Free(data->pospatt);
  For(i,data->n_otu)
    {
      Free(data->c_seq[i]->name);
      Free(data->c_seq[i]->state);
      Free(data->c_seq[i]);
    }
  Free(data->c_seq);
  Free(data);
}

/*********************************************************/

void Free_Seq(seq **d, int n_otu)
{
  int i;
  For(i,n_otu)
    {
      Free(d[i]->name);
      Free(d[i]->state);
      Free(d[i]);
    }
  Free(d);
}


/*********************************************************/

void Free_All(seq **d, allseq *alldata, arbre *tree)
{
  Free_Cseq(alldata);
  Free_Seq(d,tree->n_otu);
  Free_Tree(tree);
}      

/*********************************************************/
void Free_SubTree(edge *b_fcus, node *a, node *d, arbre *tree)
{
  int i;

  if(d->tax) return;
  else
    {
      For(i,3)
	{
	  if(d->v[i] != a)
	    {
	      Free_SubTree(d->b[i],d,d->v[i],tree);
	      Free_Edge(d->b[i]);
	      Free_Node(d->v[i]);
	    }
	}
    }
}

/*********************************************************/
void Free_Tree_Ins_Tar(arbre *tree)
{
  return;
}

/*********************************************************/

void Free_Tree_Lk(arbre *tree)
{
  int i;
  edge *b;
  
  b = NULL;
  
  Free(tree->site_loglk_sorted);
  Free(tree->site_lk);
  Free(tree->w_patt);
  Free(tree->n_changes_per_site);
  Free(tree->loglk_cat);

  For(i,MMAX(tree->mod->n_omega,tree->mod->n_catq))
      Free(tree->sel_regime_prob[i]);
  Free(tree->sel_regime_prob);

  For(i,2*tree->n_otu-3)
    {
      b = tree->t_edges[i];      
      Free_Edge_Lk(tree,b);
    }
  Free_Numerical_Param(tree->numerical_param);
}


/*********************************************************/

void Free_Numerical_Param(numpar *this_one)
{
    Free(this_one->param_val);
    Free(this_one->param_size);
    Free(this_one->param_beg);
    Free(this_one);
}

/*********************************************************/

void Free_dPij(arbre *tree)
{
}

/*********************************************************/

void Free_Edge_P_Lk_Struct(edge *b, arbre *tree)
{
  int i,j;

  if(b->p_lk_left) 
    {
      For(i,tree->n_pattern)
	{
	  For(j,tree->mod->n_catg) Free(b->p_lk_left[i][j]);
	  Free(b->p_lk_left[i]);
	}
      Free(b->p_lk_left);
      b->p_lk_left = NULL;
    }

  if(b->p_lk_rght) 
    {
      For(i,tree->n_pattern)
	{
	  For(j,tree->mod->n_catg) Free(b->p_lk_rght[i][j]);
	  Free(b->p_lk_rght[i]);
	}
      Free(b->p_lk_rght);
      b->p_lk_rght = NULL;
    }
}

/*********************************************************/

void Free_Node_Lk(node *n)
{
    if(n->seq) 
        Free(n->seq);
}

/*********************************************************/

void Free_Edge_Lk(arbre *tree, edge *b)
{
  int j;
  int catg, catq,site;

  if(b->p_lk_left)
    {
      For(site,tree->n_pattern)
	{
	  For(catq,tree->mod->n_catq)
	    {
	      For(catg,tree->mod->n_catg)
		{
		  Free(b->p_lk_left[site][catq][catg]);
		}
	      Free(b->p_lk_left[site][catq]);
	    }
	  
	  Free(b->p_lk_left[site]);
	}
      Free(b->p_lk_left);
      Free(b->sum_scale_f_left);
    }

  if(b->p_lk_rght)
    {
      For(site,tree->n_pattern)
	{
	  For(catq,tree->mod->n_catq)
	    {
	      For(catg,tree->mod->n_catg)
		{
		  Free(b->p_lk_rght[site][catq][catg]);
		}
	      Free(b->p_lk_rght[site][catq]);
	    }
	  Free(b->p_lk_rght[site]);
	}
      Free(b->p_lk_rght);
      Free(b->sum_scale_f_rght);
    }

  For(catq,tree->mod->n_catq)
    {
      For(catg,tree->mod->n_catg)
	{
	  For(j,tree->mod->ns)
	    {
	      Free(b->Pij_rr[catq][catg][j]);
	    }
	  
	  Free(b->Pij_rr[catq][catg]);
}
      Free(b->Pij_rr[catq]);
    }

  Free(b->Pij_rr);
}
  
/*********************************************************/

void Free_Model(model *mod)
{
  Free(mod->omega);
  Free(mod->omega_min);
  Free(mod->omega_max);
  Free(mod->omega_proba);
  Free(mod->selec_reg_freq);
  Free(mod->rsubst_rate);

  Free(mod->pi);
  Free(mod->r_proba);
  Free(mod->rr_mixturem);
  Free(mod->gtr_param);

  Free(mod->s_opt);
  Free_Code(mod->c_code);

  Free(mod);
}

/*********************************************************/

void Free(void *p)
{
  free(p);
}

/*********************************************************/

void Free_Input(option *input)
{
  Free(input->use_default_ssvs_switches);
  Free(input->output_tree_file);
  Free(input->output_stat_file);
  Free(input->seqfile);
  Free(input->input_tree_file);
  Free(input->use_default_dnds);
  Free(input);
}

/*********************************************************/

void Free_Code(code *c_code)
{
  int i;

  Free(c_code->n_diff_b_2_codons);
  Free(c_code->tstvtable);
  Free(c_code->sense_c);
  Free(c_code->name);
  Free(c_code->from_64_2_61);
  Free(c_code->aa);
  Free(c_code->gtr_4_codon);
  For(i,c_code->number_of_codes)
    Free(c_code->genetc[i]);
  Free(c_code->genetc);

  Free(c_code);
}

/*********************************************************/

void Free_Qmat(qmat *this1, model *mod)
{
  Free(this1->u_mat);
  Free(this1->v_mat);
  Free(this1->root_vct);
  Free(this1->expD_mr_vct);
  Free(this1->qmat);
  Free(this1->qmat_proba);
  if(mod->switch_modelname != NO_SWITCH) Free(this1->trans_qmat_proba);
  Free(this1->omega);
  Free(this1->omega_proba);
  if(mod->switch_modelname == NO_SWITCH) Free(this1->trans_omega_proba);
  Free(this1->theta);
  Free(this1->t_omega);
  Free(this1->t_omega_proba);
  Free(this1);
  
}

/*********************************************************/

void Free_NNI(nni *this1)
{
  Free(this1->bl_info);
  Free(this1);
}

/*********************************************************/

void Free_Qmat_Struct(model *mod)
{
  Free_Qmat(mod->qmat_struct[0],mod);
  Free(mod->qmat_struct);
}

/*********************************************************/

void Free_Integral_Term_On_One_Edge(edge *b, arbre *tree)
{
  int i,j;

  For(i,tree->mod->ns)
    {
      For(j,tree->mod->ns)
	{
	  Free(b->integral[i][j]);
	}
      Free(b->integral[i]);
    }
  Free(b->integral);
}

/*********************************************************/

