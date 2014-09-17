/*

FitModeL :  a program that Fits substitution models to the data
using Maximum Likelihood.  

Copyright (C) Stephane Guindon. jul. 2004 onward

All parts of  the source except where indicated  are distributed under
the GNU public licence.  See http://www.opensource.org for details.

*/

#include "config.h"
#include "utilities.h"
#include "lk.h"
#include "optimiz.h"
#include "models.h"
#include "free.h"


/*********************************************************/

void Init_Tips_At_One_Site_Nucleotides(char state, fit_double **p_lk)
{
  /* Initialize a vector of partial likelihood at a tip of the tree.
     'state' is the nucleotide observed at this tip */


  switch(state){
  case 'A' : (*p_lk)[0]=1.; (*p_lk)[1]=(*p_lk)[2]=(*p_lk)[3]=.0;
    break;
  case 'C' : (*p_lk)[1]=1.; (*p_lk)[0]=(*p_lk)[2]=(*p_lk)[3]=.0;
    break;
  case 'G' : (*p_lk)[2]=1.; (*p_lk)[1]=(*p_lk)[0]=(*p_lk)[3]=.0;
    break;
  case 'T' : (*p_lk)[3]=1.; (*p_lk)[1]=(*p_lk)[2]=(*p_lk)[0]=.0;
    break;
  case 'U' : (*p_lk)[3]=1.; (*p_lk)[1]=(*p_lk)[2]=(*p_lk)[0]=.0;
    break;
  case 'M' : (*p_lk)[0]=(*p_lk)[1]=1.; (*p_lk)[2]=(*p_lk)[3]=.0;
    break;
  case 'R' : (*p_lk)[0]=(*p_lk)[2]=1.; (*p_lk)[1]=(*p_lk)[3]=.0;
    break;
  case 'W' : (*p_lk)[0]=(*p_lk)[3]=1.; (*p_lk)[1]=(*p_lk)[2]=.0;
    break;
  case 'S' : (*p_lk)[1]=(*p_lk)[2]=1.; (*p_lk)[0]=(*p_lk)[3]=.0;
    break;
  case 'Y' : (*p_lk)[1]=(*p_lk)[3]=1.; (*p_lk)[0]=(*p_lk)[2]=.0;
    break;
  case 'K' : (*p_lk)[2]=(*p_lk)[3]=1.; (*p_lk)[0]=(*p_lk)[1]=.0;
    break;
  case 'B' : (*p_lk)[1]=(*p_lk)[2]=(*p_lk)[3]=1.; (*p_lk)[0]=.0;
    break;
  case 'D' : (*p_lk)[0]=(*p_lk)[2]=(*p_lk)[3]=1.; (*p_lk)[1]=.0;
    break;
  case 'H' : (*p_lk)[0]=(*p_lk)[1]=(*p_lk)[3]=1.; (*p_lk)[2]=.0;
    break;
  case 'V' : (*p_lk)[0]=(*p_lk)[1]=(*p_lk)[2]=1.; (*p_lk)[3]=.0;
    break;
  case 'N' : case 'X' : case '?' : case 'O' : case '-' : 
    (*p_lk)[0]=(*p_lk)[1]=(*p_lk)[2]=(*p_lk)[3]=1.;break;
  default : 
    {
      printf("\n. Unknown character state : %c\n",state);
      Exit("\n. Init failed (check the data type)\n"); 
      break;
    }
  }
}

/*********************************************************/

void Init_Tips_At_One_Site_AA(char aa, fit_double **p_lk)
{
  /* Initialize a vector of partial likelihood at a tip of the tree.
     'state' is the amino acid observed at this tip */

  int i;

  For(i,20) (*p_lk)[i] = .0;

  switch(aa){
  case 'A' : (*p_lk)[0]= 1.; break;/* Alanine */
  case 'R' : (*p_lk)[1]= 1.; break;/* Arginine */
  case 'N' : (*p_lk)[2]= 1.; break;/* Asparagine */
  case 'D' : (*p_lk)[3]= 1.; break;/* Aspartic acid */
  case 'C' : (*p_lk)[4]= 1.; break;/* Cysteine */
  case 'Q' : (*p_lk)[5]= 1.; break;/* Glutamine */
  case 'E' : (*p_lk)[6]= 1.; break;/* Glutamic acid */ 
  case 'G' : (*p_lk)[7]= 1.; break;/* Glycine */
  case 'H' : (*p_lk)[8]= 1.; break;/* Histidine */ 
  case 'I' : (*p_lk)[9]= 1.; break;/* Isoleucine */
  case 'L' : (*p_lk)[10]=1.; break;/* Leucine */
  case 'K' : (*p_lk)[11]=1.; break;/* Lysine */
  case 'M' : (*p_lk)[12]=1.; break;/* Methionine */
  case 'F' : (*p_lk)[13]=1.; break;/* Phenylalanin */
  case 'P' : (*p_lk)[14]=1.; break;/* Proline */
  case 'S' : (*p_lk)[15]=1.; break;/* Serine */
  case 'T' : (*p_lk)[16]=1.; break;/* Threonine */
  case 'W' : (*p_lk)[17]=1.; break;/* Tryptophan */
  case 'Y' : (*p_lk)[18]=1.; break;/* Tyrosine */
  case 'V' : (*p_lk)[19]=1.; break;/* Valine */
    
  case 'B' : (*p_lk)[2]= 1.; break;/* Asparagine */
  case 'Z' : (*p_lk)[5]= 1.; break;/* Glutamine */

  case 'X' : case '?' : case '-' : For(i,20) (*p_lk)[i] = 1.; break;
  default : 
    {
      printf("\n. Unknown character state : %c\n",aa);
      Exit("\n. Init failed (check the data type)\n");  
      break;
    }
  }
}

/*********************************************************/

void Init_Tips_At_One_Site_Codon(char state1, char state2, char state3, fit_double **p_lk, code *c_code, int ns, int *stop)
{
  /* Initialize a vector of partial likelihood at a tip of the tree.
     'state' is the codon observed at this tip */

  fit_double *p_lk_state1, *p_lk_state2, *p_lk_state3;
  int i,j,k,l;
  int num_codon;

  p_lk_state1 = (fit_double *)mCalloc(4,sizeof(fit_double));
  p_lk_state2 = (fit_double *)mCalloc(4,sizeof(fit_double));
  p_lk_state3 = (fit_double *)mCalloc(4,sizeof(fit_double));

  Init_Tips_At_One_Site_Nucleotides(state1,&p_lk_state1);
  Init_Tips_At_One_Site_Nucleotides(state2,&p_lk_state2);
  Init_Tips_At_One_Site_Nucleotides(state3,&p_lk_state3);
  
  *stop = 0;

  For(i,4)
    {
      if(p_lk_state1[i])
	{
	  For(j,4)
	    {
	      if(p_lk_state2[j])
		{
		  For(k,4)
		    {
		      if(p_lk_state3[k])
			{
			  num_codon = i*16+j*4+k;
			  if(c_code->from_64_2_61[num_codon] > -1) 
			    {
			      l = 0;
			      do
				{
				  (*p_lk)[c_code->from_64_2_61[num_codon]+
					  l*c_code->n_sense_c] = 1.0;
				  l++;
				}
			      while(l*c_code->n_sense_c < ns);
			    }
			  else
			    {
			      if(
				 (!Is_Ambigu(&state1,NT,1)) &&
				 (!Is_Ambigu(&state2,NT,1)) &&
				 (!Is_Ambigu(&state3,NT,1))
				 )
				{
				  int i;
                                  *stop = 1;
				  l = 0;
				  do
				    {
				      For(i,61) (*p_lk)[i+l*c_code->n_sense_c] = 1.0;
				      l++;
				    }
				  while(l*c_code->n_sense_c < ns);
				}
			    }
			}
		    }
		}
	    }
	}
    }
  Free(p_lk_state1);
  Free(p_lk_state2);
  Free(p_lk_state3);
}

/*********************************************************/

void Get_All_Partial_Lk(arbre *tree, edge *b_fcus, node *a, node *d)
{
  /* Core of Felsenstein's algorithm (without scaling to avoid underflows) */

  int i,j;
  fit_double p1_lk1,p2_lk2;
  fit_double ****p_lk,****p_lk_v1,****p_lk_v2;
  int catg,catq,site;
  fit_double ****Pij1,****Pij2;
  
  
  if(d->tax) return;
  else
    {
      int dir1,dir2;
      
      dir1=dir2=-1;
      For(i,3) if(d->v[i] != a) (dir1<0)?(dir1=i):(dir2=i);

      p_lk = 
	(d == b_fcus->left)?
	(b_fcus->p_lk_left):
	(b_fcus->p_lk_rght);
      
      p_lk_v1 = 
	(d == d->b[dir1]->left)?
	(d->b[dir1]->p_lk_rght):
	(d->b[dir1]->p_lk_left);
      
      p_lk_v2 = 
	(d == d->b[dir2]->left)?
	(d->b[dir2]->p_lk_rght):
	(d->b[dir2]->p_lk_left);
      
      Pij1 = d->b[dir1]->Pij_rr;
      Pij2 = d->b[dir2]->Pij_rr;
      

      For(catq,tree->mod->n_catq)
	{
	  For(catg,tree->mod->n_catg)
	    {	  
	      For(site,tree->n_pattern)
		{
		  For(i,tree->mod->ns)
		    {
		      p1_lk1 = p2_lk2 = .0;
		      For(j,tree->mod->ns)
			{
			  p1_lk1 += (Pij1[catq][catg][i][j] * p_lk_v1[site][catq][catg][j]);
			  p2_lk2 += (Pij2[catq][catg][i][j] * p_lk_v2[site][catq][catg][j]);
			}
		      
		      p_lk[site][catq][catg][i] = p1_lk1*p2_lk2;


		      if(p_lk[site][catq][catg][i] < LIM_SCALE_VAL)
			{
			  printf("\n. WARNING : scaling is required at site %d\n",site);
			  printf("\n. Alpha = %f\n",tree->mod->alpha);
			  Exit("");
			}
		    }
                }
	    }
	}
    }
}

/*********************************************************/

void Get_All_Partial_Lk_Scale(arbre *tree, edge *b_fcus, node *a, node *d)
{
  /* Core of Felsenstein's algorithm (with scaling to avoid underflows) */

  int i,j;
  fit_double p1_lk1,p2_lk2;
  fit_double ****p_lk,****p_lk_v1,****p_lk_v2;
  int catq,catg,site;
  fit_double ****Pij1,****Pij2;
  fit_double max_p_lk;
  fit_double sum_scale_d1,sum_scale_d2;
  fit_double try;
  
  
  p1_lk1 = p2_lk2 = .0;
  if(d->tax) return;
  else
    {
      int dir1,dir2;
      
      dir1=dir2=-1;
      For(i,3) if(d->v[i] != a) (dir1<0)?(dir1=i):(dir2=i);

/*       if(b_fcus->l < BL_MIN) b_fcus->l = BL_MIN; */

      p_lk = 
	(d == b_fcus->left)?
	(b_fcus->p_lk_left):
	(b_fcus->p_lk_rght);
      
      p_lk_v1 = 
	(d == d->b[dir1]->left)?
	(d->b[dir1]->p_lk_rght):
	(d->b[dir1]->p_lk_left);
      
      p_lk_v2 = 
	(d == d->b[dir2]->left)?
	(d->b[dir2]->p_lk_rght):
	(d->b[dir2]->p_lk_left);
      
      Pij1 = d->b[dir1]->Pij_rr;
      Pij2 = d->b[dir2]->Pij_rr;
      

      For(site,tree->n_pattern)
	{	  
	  sum_scale_d1 = sum_scale_d2 = .0;

	  (d == d->b[dir1]->left)?
	    (sum_scale_d1 = d->b[dir1]->sum_scale_f_rght[site]):
	    (sum_scale_d1 = d->b[dir1]->sum_scale_f_left[site]);
	  
	  (d == d->b[dir2]->left)?
	    (sum_scale_d2 = d->b[dir2]->sum_scale_f_rght[site]):
	    (sum_scale_d2 = d->b[dir2]->sum_scale_f_left[site]);

	  (d == b_fcus->left)?
	    (b_fcus->sum_scale_f_left[site] = sum_scale_d1 + sum_scale_d2):
	    (b_fcus->sum_scale_f_rght[site] = sum_scale_d1 + sum_scale_d2);
	  
	  max_p_lk = -MDBL_MAX;
	  
	  For(catq,tree->mod->n_catq)
	    {
	      For(catg,tree->mod->n_catg)
		{
		  For(i,tree->mod->ns)
		    {
		      p1_lk1 = p2_lk2 = .0;
		      For(j,tree->mod->ns)
			{
			  p1_lk1 += (Pij1[catq][catg][i][j] * p_lk_v1[site][catq][catg][j]);
			  p2_lk2 += (Pij2[catq][catg][i][j] * p_lk_v2[site][catq][catg][j]);
			}
		      
		      try = p1_lk1*p2_lk2;
		      
		      p_lk[site][catq][catg][i]=try;
		      
/* 		      printf("SITE %3d p_lk = %E (i=%3d) (CATG=%3d)\n", */
/* 			     site, */
/* 			     p_lk[site][catq][catg][i], */
/* 			     i,catg); */

		      if(p_lk[site][catq][catg][i] > max_p_lk) max_p_lk = p_lk[site][catq][catg][i];
		    }
		}
	    }
	      


	  if(max_p_lk < LIM_SCALE_VAL)
	    {

	      if(max_p_lk < 10.*MDBL_MIN) 
		{
		  Exit("\n. Numerical underflow !\n");
		}
	      For(catq,tree->mod->n_catq)
		{
		  For(catg,tree->mod->n_catg)
		    {
		      For(i,tree->mod->ns)
			{
			  p_lk[site][catq][catg][i] /= max_p_lk;
			  
/* 			  if(p_lk[site][catq][catg][i] > MDBL_MAX) */
/* 			    { */
/* 			      printf("\n. Numerical overflow ! \n"); */
/* 			      /\*  			  p_lk[site][catg][i] = p_lk[site][catg-1][i] ; *\/ */
/* 			    } */
			}
		    }
		}
		      
	      (d == b_fcus->left)?
		(b_fcus->sum_scale_f_left[site] += (fit_double)log(max_p_lk)):
		(b_fcus->sum_scale_f_rght[site] += (fit_double)log(max_p_lk));
	    }

	  if(max_p_lk > (1./LIM_SCALE_VAL))
	    {
	      For(catq,tree->mod->n_catq)
		{
		  For(catg,tree->mod->n_catg)
		    {
		      For(i,tree->mod->ns)
			{
			  p_lk[site][catq][catg][i] /= max_p_lk;
			  			  
/* 			  if(p_lk[site][catq][catg][i] < MDBL_MIN) */
/* 			    { */
/* 			      printf("\n. Numerical overflow ! \n"); */
/* 			      /\*  			  p_lk[site][catg][i] = p_lk[site][catg-1][i] ; *\/ */
/* 			    } */
			}
		    }
		}
		      
	      (d == b_fcus->left)?
		(b_fcus->sum_scale_f_left[site] += (fit_double)log(max_p_lk)):
		(b_fcus->sum_scale_f_rght[site] += (fit_double)log(max_p_lk));
	    }
	}
    }
}

/*********************************************************/

void Pre_Order_Lk(node *pere, node *fils, arbre *tree)
{
  /* Pre-order traversal for the computation of vectors of partial likelihood */
 
  int i,dir1,dir2,dir3;

  dir1 = dir2 = dir3 = -1;

  if(fils->tax) return;
  else
    {
      For(i,3)
	{
	  if(fils->v[i] != pere)
	    {
	      Pre_Order_Lk(fils,fils->v[i],tree);
	      if(dir1 < 0) dir1 = i;
	      else dir2 = i;
	    }
	  else dir3 = i;
	}

      (tree->n_otu > LIM_SCALE)?
	(Get_All_Partial_Lk_Scale(tree,fils->b[dir3],pere,fils)):
	(Get_All_Partial_Lk(tree,fils->b[dir3],pere,fils));
    }
}

/*********************************************************/

void Post_Order_Lk(node *pere, node *fils, arbre *tree)
{
  /* Post-order traversal for the computation of vectors of partial likelihood */

  int i,j,dir1,dir2;

  dir1 = dir2 = -1;
  
  if(fils->tax) return;
  else
    {
      For(i,3)
	{
	  if(fils->v[i] != pere)
	    {
	      For(j,3)
		{
		  if(j != i)
		    {
		      if(dir1 < 0) dir1 = j;
		      else dir2 = j;
		    }
		}
	      
/* 	      Partial_Lk(fils,dir1,dir2,tree); */
	      (tree->n_otu > LIM_SCALE)?
	      (Get_All_Partial_Lk_Scale(tree,fils->b[i],fils->v[i],fils)):
	      (Get_All_Partial_Lk(tree,fils->b[i],fils->v[i],fils));
	      dir1 = dir2 = -1;
	      Post_Order_Lk(fils,fils->v[i],tree);
	    }
	}
    }
}

/*********************************************************/

void Lk(arbre *tree)
{
  /* Compute the log likelihood of a model given the sequences */ 

  int br,site;

  if(tree->numerical_param->replace_num_param) Replace_Num_Parameters(tree);
  
  Set_Model_Parameters(tree);

  Get_All_PMat_Post(tree->noeud[0],
		    tree->noeud[0]->v[0],
		    tree->noeud[0]->b[0],
		    tree);

  Pre_Order_Lk(tree->noeud[0],tree->noeud[0]->v[0],tree);
  if(tree->both_sides)
    Post_Order_Lk(tree->noeud[0],
		  tree->noeud[0]->v[0],
		  tree);


  tree->tot_loglk = .0;    
  tree->curr_site =  0;

  //
  tree->loglk_cat[0] = .0;
  tree->loglk_cat[1] = .0;
  tree->loglk_cat[2] = .0;
  //


  For(site,tree->n_pattern)
    {
      tree->curr_site=site;
      Site_Lk(tree,tree->data);
    }  

/*   qksort(tree->site_loglk_sorted, 0, tree->n_pattern-1); */
  tree->tot_loglk = .0;
  For(site, tree->n_pattern)
    tree->tot_loglk += tree->site_loglk_sorted[site];
                                                                                                 
  For(br,2*tree->n_otu-3)
    {
      if(tree->t_edges[br]->get_p_lk_left) tree->t_edges[br]->ud_p_lk_left = 1;
      if(tree->t_edges[br]->get_p_lk_rght) tree->t_edges[br]->ud_p_lk_rght = 1;
    }

  if(tree->mod->s_opt->print) Print_Param(stdout,tree);

  if(tree->numerical_param->replace_num_param) Replace_Num_Parameters(tree);

}

/*********************************************************/

fit_double Lk_At_Given_Edge(arbre *tree, edge *b_fcus)
{
  (tree->mod->model_number >= 20)?
    (Lk_At_Given_Edge_Codon(tree,b_fcus)):
    (Lk_At_Given_Edge_Nt_AA(tree,b_fcus));
  return tree->tot_loglk;
}

/*********************************************************/

void Site_Lk(arbre *tree, allseq *alldata)
{
  (tree->mod->model_number >= 20)?
    (Site_Lk_Codon(tree,alldata)):
    (Site_Lk_Nt_AA(tree,alldata));
}

/*********************************************************/

void Site_Lk_Codon(arbre *tree, allseq *alldata)
{
  /* Compute the likelihood at a given site (codon-based models) */

  int k,l;
  fit_double log_site_lk,site_lk;
  fit_double site_lk_level_1,site_lk_level_2, site_lk_level_3,site_lk_level_4;
  edge *eroot;
  int catq, catg;
  int state_elsewhere;
  fit_double site_lk_cat[30];

  For(k,30) site_lk_cat[k] = 0.0;

  if(!tree->w_patt[tree->curr_site]) return;

  eroot = tree->noeud[0]->b[0];
    
  
  state_elsewhere = tree->data->invar[tree->curr_site];
  
  log_site_lk = site_lk = .0;      

  site_lk_level_1 = .0;
  For(catq,tree->mod->n_catq)
    {
      site_lk_level_2 = .0;
      For(catg,tree->mod->n_catg)
	{
	  site_lk_level_3 = .0;
	  	  
	  For(k,tree->mod->ns)
	    {
	      site_lk_level_4 = .0;
	      For(l,tree->mod->ns)
		{
		  site_lk_level_4 += 
		    eroot->Pij_rr[catq][catg][k][l] *		    
		    eroot->p_lk_rght[tree->curr_site][catq][catg][l]; 
		}
	      
	      site_lk_level_3 +=
		site_lk_level_4 *
		tree->mod->pi[k%tree->mod->ns_codon] * 
		(*(eroot->qmat_struct->t_omega_proba[catq*eroot->qmat_struct->n_omega+
						     (int)k/tree->mod->ns_codon])) *
		eroot->p_lk_left[tree->curr_site][catq][catg][k];
	    }
	  site_lk_level_2 += site_lk_level_3 * tree->mod->r_proba[catg];
	}
      
      site_lk_level_1 += site_lk_level_2 * eroot->qmat_struct->qmat_proba[catq];
      site_lk_cat[catq] =  site_lk_level_2;
    }
  
/*   printf("%E (%f)  %E (%f)  %E (%f)\n", */
/* 	 site_lk_cat[0], */
/* 	 eroot->qmat_struct->qmat_proba[0], */
/* 	 site_lk_cat[1], */
/* 	 eroot->qmat_struct->qmat_proba[1], */
/* 	 site_lk_cat[2], */
/* 	 eroot->qmat_struct->qmat_proba[2]); */


  site_lk = site_lk_level_1;

/*   For(k,tree->mod->n_omega) */
#if defined POSSELSITEID

  if((tree->input->pos_sel_sites_id_method == ADAPTIVE_CONTROL) ||
     (tree->input->pos_sel_sites_id_method == PARAMETRIC_BOOTSTRAP))
      {
	fit_double p0,p1,p2,w0,w1,w2;

	w0 = tree->mod->qmat_struct[0]->omega[0];
	w1 = tree->mod->qmat_struct[0]->omega[1];
	w2 = tree->mod->qmat_struct[0]->omega[2];

	p0 = 
	  site_lk_cat[0] * eroot->qmat_struct->qmat_proba[0] / 
	  site_lk; 

	p1 = 
	  site_lk_cat[1] * eroot->qmat_struct->qmat_proba[1] / 
	  site_lk; 

	p2 = 
	  site_lk_cat[2] * eroot->qmat_struct->qmat_proba[2] / 
	  site_lk; 

          tree->sel_regime_prob[2][tree->curr_site] = p2;
      }
  else if((tree->input->pos_sel_sites_id_method == NEWTON_ETAL)  || 
          (tree->input->pos_sel_sites_id_method == FIXED_REGION))
      {
          
          tree->sel_regime_prob[2][tree->curr_site] = 
              site_lk_cat[2] * eroot->qmat_struct->qmat_proba[2] /
              site_lk;
      }

#endif

#if defined FITMODEL

  if(tree->mod->switch_modelname == NO_SWITCH)
    {
      tree->sel_regime_prob[0][tree->curr_site] =
	site_lk_cat[0] * eroot->qmat_struct->qmat_proba[0] /
	site_lk;
      
      if((tree->mod->n_catq >= 2) || (tree->mod->n_omega >= 2))
	tree->sel_regime_prob[1][tree->curr_site] =
          site_lk_cat[1] * eroot->qmat_struct->qmat_proba[1] /
          site_lk;
      
      if((tree->mod->n_catq >= 3) || (tree->mod->n_omega >= 3))
	tree->sel_regime_prob[2][tree->curr_site] =
	  site_lk_cat[2] * eroot->qmat_struct->qmat_proba[2] /
	  site_lk;
    }

  
#endif



/*       site_lk_cat[k] *  */
/*       eroot->qmat_struct->qmat_proba[k] / site_lk; */
  	   
/*   fprintf(tree->input->fp_output_stat, */
/*   printf( */
/*           "%d %E %E %E %E %E\n", */
/*           tree->curr_site, */
/*           site_lk_cat[2] * eroot->qmat_struct->qmat_proba[2] / site_lk, */
/*           site_lk_cat[0] * eroot->qmat_struct->qmat_proba[0] / site_lk, */
/*           site_lk_cat[1] * eroot->qmat_struct->qmat_proba[1] / site_lk, */
/*           site_lk_cat[2] * eroot->qmat_struct->qmat_proba[2] / site_lk, */
/*           log(site_lk) */
/*           ); */
 
  /* take into account scaling factors and invariants */
  if (!tree->mod->invar)
    {
      log_site_lk = (fit_double)log(site_lk) + 
	            eroot->sum_scale_f_left[tree->curr_site] + 
	            eroot->sum_scale_f_rght[tree->curr_site];
    }
  else
    {
      if ((fit_double)tree->data->invar[tree->curr_site] > -0.5)
        {
          if (!(eroot->sum_scale_f_left[tree->curr_site] + 
		eroot->sum_scale_f_rght[tree->curr_site]==0.0))
            site_lk *= 
	      (fit_double)exp(eroot->sum_scale_f_left[tree->curr_site] + 
		  eroot->sum_scale_f_rght[tree->curr_site]);
	  
	  log_site_lk = (fit_double)log(
			    site_lk*(1.0-tree->mod->pinvar) + 
			    tree->mod->pinvar*tree->mod->pi[state_elsewhere]
			    );
        }
      else
        {
          log_site_lk = 
	    (fit_double)log(site_lk*(1.0-tree->mod->pinvar)) + 
	    eroot->sum_scale_f_left[tree->curr_site] + 
	    eroot->sum_scale_f_rght[tree->curr_site];
        }
    }


  if(log_site_lk < -MDBL_MAX)
    {      
      printf("\n. site_lk = %E\n",site_lk);
      printf("\n. log_site_lk < -MDBL_MAX at site %d \n",tree->curr_site);
      Exit("\n");
/*       log_site_lk = log(1.E-200); */
    }

  tree->site_lk[tree->curr_site] = log_site_lk;

  tree->loglk_cat[0] += tree->w_patt[tree->curr_site]*log(site_lk_cat[0]);
  tree->loglk_cat[1] += tree->w_patt[tree->curr_site]*log(site_lk_cat[1]);
  tree->loglk_cat[2] += tree->w_patt[tree->curr_site]*log(site_lk_cat[2]);

  tree->site_loglk_sorted[tree->curr_site] = tree->w_patt[tree->curr_site]*log_site_lk;
}

/*********************************************************/

fit_double Lk_At_Given_Edge_Codon(arbre *tree, edge *b_fcus)
{
  /* Compute the likelihood using the vector of partial likelihood
     at both extremities of a given branch (codon-based models) */

  int site,catq,catg,k,l;
  fit_double site_lk,log_site_lk;
  fit_double site_lk_level_1,site_lk_level_2, site_lk_level_3,site_lk_level_4;
  fit_double rr,len;
  
  tree->tot_loglk = log_site_lk = .0;

  if(b_fcus->l < BL_MIN) b_fcus->l = BL_MIN;
  if(b_fcus->l > BL_MAX) b_fcus->l = BL_MAX;
  
  len = rr = -1.0;
  For(catg,tree->mod->n_catg)
    {
      len = b_fcus->l*tree->mod->rr_mixturem[catg];
      if(len < BL_MIN) len = BL_MIN;
      if(len > BL_MAX) len = BL_MAX;
     
      For(catq,tree->mod->n_catq)
	{
	  if(tree->mod->analytical)
	    {
	      PMat(len,
		   tree->mod,
		   b_fcus->Pij_rr[catq][catg],
		   NULL,tree);
	    }
	  else
	    {
	      b_fcus->qmat_struct->curr_qmat_cat = catq;
	      PMat(len,
		   tree->mod,
		   b_fcus->Pij_rr[catq][catg],
		   b_fcus->qmat_struct,tree);
	    }
	}
    }



  For(site,tree->n_pattern)
    {
      if(tree->w_patt[site])
	{
	  log_site_lk = site_lk = .0;      
	
	  site_lk_level_1 = .0;
	  For(catq,tree->mod->n_catq)
	    {
	      site_lk_level_2 = .0;
	      For(catg,tree->mod->n_catg)
		{
		  site_lk_level_3 = .0;
		  For(k,tree->mod->ns)
		    {
		      site_lk_level_4 = .0;
		      For(l,tree->mod->ns)
			{  
			  site_lk_level_4 +=
			    b_fcus->p_lk_rght[site][catq][catg][l] *
			    b_fcus->Pij_rr[catq][catg][k][l];
			  
/* 			  if(site == 15)  */
/* 			    { */
/* 			      printf("site_lk_level_4=%E\n",site_lk_level_4); */
/* 			    } */
			}

		      site_lk_level_3 += 
			site_lk_level_4 *
 			tree->mod->pi[k%tree->mod->ns_codon] *
			(*(b_fcus->qmat_struct->t_omega_proba[catq*
							      b_fcus->qmat_struct->n_omega+
							      (int)k/tree->mod->ns_codon])) *
			b_fcus->p_lk_left[site][catq][catg][k];

		    }
		  site_lk_level_2 +=
		    site_lk_level_3 *
		    tree->mod->r_proba[catg];    
		}
	      site_lk_level_1 += 
		site_lk_level_2 *		
		b_fcus->qmat_struct->qmat_proba[catq];
	    }

	  site_lk = site_lk_level_1;

	  /* take into account scaling factors and invariants */
	  if (!tree->mod->invar)
	    {
	      log_site_lk = (fit_double)log(site_lk) + 
		b_fcus->sum_scale_f_left[site] + 
		b_fcus->sum_scale_f_rght[site];
	    }
	  else
	    {
	      if ((fit_double)tree->data->invar[site] > -0.5)
		{
		  if (!(b_fcus->sum_scale_f_left[site] + 
			b_fcus->sum_scale_f_rght[site]==0.0))
		    site_lk *= 
		      (fit_double)exp(b_fcus->sum_scale_f_left[site] + 
			  b_fcus->sum_scale_f_rght[site]);
		  
		  log_site_lk = (fit_double)log(
				    site_lk*(1.0-tree->mod->pinvar) + 
				    tree->mod->pinvar*tree->mod->pi[tree->data->invar[site]]
				    );
		}
	      else
		{
		  log_site_lk = 
		    (fit_double)log(site_lk*(1.0-tree->mod->pinvar)) + 
		    b_fcus->sum_scale_f_left[site] + 
		    b_fcus->sum_scale_f_rght[site];
		}
	    }

          if(log_site_lk < -MDBL_MAX) Exit("\nlog_site_lk < -MDBL_MAX\n");

          tree->site_lk[site] = log_site_lk;

/* 	  printf("site %d %f\n",site,log_site_lk); */

          tree->site_loglk_sorted[site] = tree->w_patt[site]*log_site_lk;
	}
    }

  /* Sort and add likelihood at each site from the smallest value to the biggest.
     This is done to limit round-off errors */
/*   qksort(tree->site_loglk_sorted, 0, tree->n_pattern-1); */
  tree->tot_loglk = .0;
  For(site, tree->n_pattern) 
    {
      tree->tot_loglk += tree->site_loglk_sorted[site];
    }
  return tree->tot_loglk;
}
  
/*********************************************************/

void Site_Lk_Nt_AA(arbre *tree, allseq *alldata)
{
  /* Compute the likelihood at a given site (nucleotide or amino acid based models) */

  int k,l;
  fit_double log_site_lk,site_lk;
  fit_double site_lk_level_1,site_lk_level_2, site_lk_level_3,site_lk_level_4;
  int left;
  edge *eroot;
  int is_ambigu;
  int catq, catg;
  int state_root, state_elsewhere;


  if(!tree->w_patt[tree->curr_site]) return;

  eroot = tree->noeud[0]->b[0];
  (eroot->rght->tax)?(left=1):(left=0);
  
  is_ambigu = alldata->ambigu[tree->curr_site];
  
  state_root = -1;
  state_root = Assign_State(alldata->c_seq[eroot->rght->num]->state+tree->curr_site*tree->mod->stepsize,
			    tree->mod->datatype,
			    tree->mod->stepsize,
			    tree->mod->c_code);

  state_elsewhere = -1;
  state_elsewhere = tree->data->invar[tree->curr_site];

  /**/
/*    is_ambigu = 1; */
  /**/
  

  log_site_lk = site_lk = .0;      
  
  if(is_ambigu)
    {
      site_lk_level_1 = .0;
      For(catq,tree->mod->n_catq)
	{
	  site_lk_level_2 = .0;
	  For(catg,tree->mod->n_catg)
	    {
	      site_lk_level_3 = .0;
              
	      For(k,tree->mod->ns)
		{
		  site_lk_level_4 = .0;
		  For(l,tree->mod->ns)
		    {
		      site_lk_level_4 += 
			eroot->Pij_rr[catq][catg][k][l] *		    
			eroot->p_lk_rght[tree->curr_site][catq][catg][l]; 
		    }
		  
		  site_lk_level_3 +=
		    site_lk_level_4 *
		    tree->mod->pi[k] * 
		    eroot->p_lk_left[tree->curr_site][catq][catg][k];
		}
	      site_lk_level_2 +=
		site_lk_level_3 *
		tree->mod->r_proba[catg];
	    }
	  
	  site_lk_level_1 +=
	    site_lk_level_2 *
	    eroot->qmat_struct->qmat_proba[catq];
	}
    }
  else
    {
      
      site_lk_level_1 = .0;
      For(catq,tree->mod->n_catq)
	{
	  site_lk_level_2 = .0;
	  For(catg,tree->mod->n_catg)
	    {
	      site_lk_level_3 = .0;                          
	      For(k,tree->mod->ns)
		{
		  site_lk_level_3 += 
		    eroot->Pij_rr[catq][catg][k][state_root] *		    
		    eroot->p_lk_rght[tree->curr_site][catq][catg][state_root] *
		    tree->mod->pi[k] * 
		    eroot->p_lk_left[tree->curr_site][catq][catg][k];
		  
		}
	      site_lk_level_2 +=
		site_lk_level_3 *
		tree->mod->r_proba[catg];
	    }
	  site_lk_level_1 +=
	    site_lk_level_2 *
	    eroot->qmat_struct->qmat_proba[catq];
	}
    }
  
  site_lk = site_lk_level_1;

  /* Take into account scaling factors and invariants */
  if (!tree->mod->invar)
    {
      log_site_lk = (fit_double)log(site_lk) + 
	eroot->sum_scale_f_left[tree->curr_site] + 
	eroot->sum_scale_f_rght[tree->curr_site];
    }
  else
    {
      if ((fit_double)tree->data->invar[tree->curr_site] > -0.5)
        {
          if (!(eroot->sum_scale_f_left[tree->curr_site] + 
		eroot->sum_scale_f_rght[tree->curr_site]==0.0))
            site_lk *= 
	      (fit_double)exp(eroot->sum_scale_f_left[tree->curr_site] + 
			      eroot->sum_scale_f_rght[tree->curr_site]);
	  
	  log_site_lk = (fit_double)log(
					site_lk*(1.0-tree->mod->pinvar) + 
					tree->mod->pinvar*tree->mod->pi[state_elsewhere]
					);
        }
      else
        {
          log_site_lk = 
	    (fit_double)log(site_lk*(1.0-tree->mod->pinvar)) + 
	    eroot->sum_scale_f_left[tree->curr_site] + 
	    eroot->sum_scale_f_rght[tree->curr_site];
        }
    }
  
  
  if(log_site_lk < -MDBL_MAX)
    {
      printf("CURR_SITE = %d\n",tree->curr_site);
      Exit("\nlog_site_lk < -MDBL_MAX\n");
      /*       log_site_lk = log(1.E-200); */
    }
  
  tree->site_lk[tree->curr_site] = log_site_lk;
  
  tree->site_loglk_sorted[tree->curr_site] = tree->w_patt[tree->curr_site]*log_site_lk;
}

/*********************************************************/

fit_double Lk_At_Given_Edge_Nt_AA(arbre *tree, edge *b_fcus)
{
  /* Compute the likelihood using the vector of partial likelihood
     at both extremities of a given branch (amino acid and nucleotide
     based models) 
  */
  
  int site,catq,catg,k,l;
  fit_double site_lk,log_site_lk;
  fit_double site_lk_level_1,site_lk_level_2, site_lk_level_3,site_lk_level_4;
  fit_double rr,len;
    
  tree->tot_loglk = log_site_lk = .0;
    
  if(b_fcus->l < BL_MIN) b_fcus->l = BL_MIN;
  if(b_fcus->l > BL_MAX) b_fcus->l = BL_MAX;
  
  
  len = rr = -1.0;
  For(catg,tree->mod->n_catg)
    {
      len = b_fcus->l*tree->mod->rr_mixturem[catg];
      if(len < BL_MIN) len = BL_MIN;
      if(len > BL_MAX) len = BL_MAX;
     
      For(catq,tree->mod->n_catq)
	{
	  if(tree->mod->analytical)
	    {
	      PMat(len,
		   tree->mod,
		   b_fcus->Pij_rr[catq][catg],
		   NULL,tree);
	    }
	  else
	    {
	      b_fcus->qmat_struct->curr_qmat_cat = catq;
	      PMat(len,
		   tree->mod,
		   b_fcus->Pij_rr[catq][catg],
		   b_fcus->qmat_struct,tree);
	    }
	}
    }



  For(site,tree->n_pattern)
    {
      if(tree->w_patt[site])
	{
	  log_site_lk = site_lk = .0;      
	
	  site_lk_level_1 = .0;
	  For(catq,tree->mod->n_catq)
	    {
	      site_lk_level_2 = .0;
	      For(catg,tree->mod->n_catg)
		{
		  site_lk_level_3 = .0;
		  For(k,tree->mod->ns)
		    {
		      site_lk_level_4 = .0;
		      For(l,tree->mod->ns)
			{  
			  site_lk_level_4 +=
                              b_fcus->p_lk_rght[site][catq][catg][l] *
                              b_fcus->Pij_rr[catq][catg][k][l];

/*                           site_lk +=  */
/*                               tree->mod->r_proba[catg] * */
/*                               tree->mod->pi[k] *  */
/*                               b_fcus->p_lk_left[site][catq][catg][k] * */
/*                               b_fcus->p_lk_rght[site][catq][catg][l] * */
/*                               b_fcus->Pij_rr[catq][catg][k][l]; */
			}

		      site_lk_level_3 += 
			site_lk_level_4 *
			tree->mod->pi[k] *
			b_fcus->p_lk_left[site][catq][catg][k];
		    }
		  site_lk_level_2 +=
		    site_lk_level_3 *
		    tree->mod->r_proba[catg];    
		}
	      site_lk_level_1 += 
		site_lk_level_2 *		
		b_fcus->qmat_struct->qmat_proba[catq];
	    }

	  site_lk = site_lk_level_1;


	  /* Take into account scaling factors and invariants */
	  if (!tree->mod->invar)
	    {
	      log_site_lk = (fit_double)log(site_lk) + 
		b_fcus->sum_scale_f_left[site] + 
		b_fcus->sum_scale_f_rght[site];
	    }
	  else
	    {
	      if ((fit_double)tree->data->invar[site] > -0.5)
		{
		  if (!(b_fcus->sum_scale_f_left[site] + 
			b_fcus->sum_scale_f_rght[site]==0.0))
		    site_lk *= 
		      (fit_double)exp(b_fcus->sum_scale_f_left[site] + 
			  b_fcus->sum_scale_f_rght[site]);
		  
		  log_site_lk = (fit_double)log(
				    site_lk*(1.0-tree->mod->pinvar) + 
				    tree->mod->pinvar*tree->mod->pi[tree->data->invar[site]]
				    );
		}
	      else
		{
		  log_site_lk = 
		    (fit_double)log(site_lk*(1.0-tree->mod->pinvar)) + 
		    b_fcus->sum_scale_f_left[site] + 
		    b_fcus->sum_scale_f_rght[site];
		}
	    }

          if(log_site_lk < -MDBL_MAX) Exit("\nlog_site_lk < -MDBL_MAX\n");

          tree->site_lk[site] = log_site_lk;

          tree->site_loglk_sorted[site] = /* code 2.3 */
            tree->w_patt[site]*log_site_lk;
	}
    }

  /* Sort and add likelihood at each site from the smallest value to the biggest.
     This is done to limit round-off errors */
/*   qksort(tree->site_loglk_sorted, 0, tree->n_pattern-1); */
  tree->tot_loglk = .0;
  For(site, tree->n_pattern) tree->tot_loglk += tree->site_loglk_sorted[site];
  
  return tree->tot_loglk;
}

/*********************************************************/

void Update_P(arbre *tree, int t_edge_num)
{
  /* Update a matrix of change probabilities */
  
  int catq,catg;
  fit_double len;

  len = -1.0;
  For(catg,tree->mod->n_catg)
    {
      len = tree->t_edges[t_edge_num]->l*tree->mod->rr_mixturem[catg];
      if(len < BL_MIN) len = BL_MIN;
      if(len > BL_MAX) len = BL_MAX;
      
      if(tree->mod->analytical)
	{
	  For(catq,tree->mod->n_catq)
	    {
	      PMat(len,
		   tree->mod,
		   tree->t_edges[t_edge_num]->Pij_rr[catq][catg],
		   NULL,tree);
	    }
	}
      else
	{
	  For(catq,tree->mod->n_catq)
	    {
	      tree->t_edges[t_edge_num]->qmat_struct->curr_qmat_cat = catq;
	      PMat(len,
		   tree->mod,
		   tree->t_edges[t_edge_num]->Pij_rr[catq][catg],
		   tree->t_edges[t_edge_num]->qmat_struct,tree);
	    }
	}
    }
}

/*********************************************************/

fit_double Return_Lk(arbre *tree)
{
  Lk(tree);
  return tree->tot_loglk;
}

/*********************************************************/

fit_double Return_Abs_Lk(arbre *tree)
{
  Lk(tree);
  return (fit_double)fabs(tree->tot_loglk);
}

/*********************************************************/

fit_double ****Get_Partial_Lk_Struct(arbre *tree, int len, int n_catg, int n_catq)
{
  /* Allocate memory for vectors of partial likelihood */

  fit_double ****p_lk;
  int site,catg,catq;

  p_lk = (fit_double ****)mCalloc(len,sizeof(fit_double ***)); 
  For(site,len)
    {
      p_lk[site] = (fit_double ***)mCalloc(n_catq,sizeof(fit_double **));
      For(catq,n_catq) 
	{
	  p_lk[site][catq] = (fit_double **)mCalloc(n_catg,sizeof(fit_double *));
	  For(catg,n_catg)
	    p_lk[site][catq][catg] = (fit_double *)mCalloc(tree->mod->ns,sizeof(fit_double ));
	}
    }

  return p_lk;
}

/*********************************************************/

void Update_P_Lk(arbre *tree, edge *b_fcus, node *n)
{
  /* Update a vector of partial likelihood using the two vectors
     of partial likelihood underneath */

  int k,l;
  int site, catg, catq;
  fit_double ****p_lk, ****p_lk_v1, ****p_lk_v2;
  fit_double **Pij1, **Pij2;
  fit_double *n_scale_f, *d1_scale_f, *d2_scale_f;
  fit_double p1_lk1,p2_lk2;
  fit_double max_p_lk;
  edge *b1, *b2;

  
  b1 = b2  = NULL;
  p_lk = p_lk_v1 = p_lk_v2 = NULL;
  max_p_lk = MDBL_MIN;


  if(n == b_fcus->left)
    {
      p_lk = b_fcus->p_lk_left;
      
      p_lk_v1 = 
      (n == n->b[b_fcus->l_v1]->left)?
      (n->b[b_fcus->l_v1]->p_lk_rght):
      (n->b[b_fcus->l_v1]->p_lk_left);

      p_lk_v2 = 
      (n == n->b[b_fcus->l_v2]->left)?
      (n->b[b_fcus->l_v2]->p_lk_rght):
      (n->b[b_fcus->l_v2]->p_lk_left);

      n_scale_f = b_fcus->sum_scale_f_left;

      d1_scale_f = 
      (n == n->b[b_fcus->l_v1]->left)?
      (n->b[b_fcus->l_v1]->sum_scale_f_rght):
      (n->b[b_fcus->l_v1]->sum_scale_f_left);

      d2_scale_f = 
      (n == n->b[b_fcus->l_v2]->left)?
      (n->b[b_fcus->l_v2]->sum_scale_f_rght):
      (n->b[b_fcus->l_v2]->sum_scale_f_left);
    
      b_fcus->get_p_lk_left = 1;
      b_fcus->ud_p_lk_left  = 1;
    }
  else
    {
      p_lk = b_fcus->p_lk_rght;
      
      p_lk_v1 = 
      (n == n->b[b_fcus->r_v1]->left)?
      (n->b[b_fcus->r_v1]->p_lk_rght):
      (n->b[b_fcus->r_v1]->p_lk_left);

      p_lk_v2 = 
      (n == n->b[b_fcus->r_v2]->left)?
      (n->b[b_fcus->r_v2]->p_lk_rght):
      (n->b[b_fcus->r_v2]->p_lk_left);

      n_scale_f = b_fcus->sum_scale_f_rght;

      d1_scale_f = 
      (n == n->b[b_fcus->r_v1]->left)?
      (n->b[b_fcus->r_v1]->sum_scale_f_rght):
      (n->b[b_fcus->r_v1]->sum_scale_f_left);

      d2_scale_f = 
      (n == n->b[b_fcus->r_v2]->left)?
      (n->b[b_fcus->r_v2]->sum_scale_f_rght):
      (n->b[b_fcus->r_v2]->sum_scale_f_left);

      b_fcus->get_p_lk_rght = 1;
      b_fcus->ud_p_lk_rght  = 1;
    }
 

  if(n == b_fcus->left) 
    {
      b1 = n->b[b_fcus->l_v1];
      b2 = n->b[b_fcus->l_v2];
    }
  else
    {
      b1 = n->b[b_fcus->r_v1];
      b2 = n->b[b_fcus->r_v2];
    }


  if(tree->n_otu <= LIM_SCALE)
    {
      /* Without scaling */

      For(site,tree->n_pattern)
	{
	  if(tree->w_patt[site])
	    {
	      For(catq,tree->mod->n_catq)
		{
		  For(catg,tree->mod->n_catg)
		    {
		      
		      Pij1 = b1->Pij_rr[catq][catg];
		      Pij2 = b2->Pij_rr[catq][catg];
		  		  
		      For(k,tree->mod->ns)
			{
			  p1_lk1 = p2_lk2 = .0;
			  
			  For(l,tree->mod->ns)
			    {
			      p1_lk1 += Pij1[k][l] * p_lk_v1[site][catq][catg][l];
			      p2_lk2 += Pij2[k][l] * p_lk_v2[site][catq][catg][l];
			    }
			  p_lk[site][catq][catg][k] = p1_lk1 * p2_lk2;
			}
		    }
		}
	    }
	}
    }
  else
    {
      /* With scaling */

      For(site,tree->n_pattern)
	{
	  if(tree->w_patt[site])
	    {
	      For(catq,tree->mod->n_catq)
		{
		  For(catg,tree->mod->n_catg)
		    {
		  
		      Pij1 = b1->Pij_rr[catq][catg];
		      Pij2 = b2->Pij_rr[catq][catg];
		  
		      if(!catg) 
			{
			  n_scale_f[site] = d1_scale_f[site] + d2_scale_f[site];
			  max_p_lk = -MDBL_MAX;
			}
		  
		      For(k,tree->mod->ns)
			{
			  p_lk[site][catq][catg][k] = .0;
			  p1_lk1 = p2_lk2     = .0;
			  
			  For(l,tree->mod->ns)
			    {
			      p1_lk1 += Pij1[k][l] * p_lk_v1[site][catq][catg][l];
			      p2_lk2 += Pij2[k][l] * p_lk_v2[site][catq][catg][l];
			    }
			  p_lk[site][catq][catg][k] = p1_lk1 * p2_lk2;
			  
			  if(p_lk[site][catq][catg][k] > max_p_lk) max_p_lk = p_lk[site][catq][catg][k];
			}
		    }
		}

	      if(max_p_lk < LIM_SCALE_VAL)
		{
		  For(catq,tree->mod->n_catq)
		    {
		      For(catg,tree->mod->n_catg)
			{
			  For(k,tree->mod->ns) 
			    {
			      p_lk[site][catq][catg][k] /= max_p_lk;
				
/* 			      if(p_lk[site][catq][catg][k] > (1./LIM_SCALE_VAL)) */
/* 				{ */
/* 				  printf("\n. Numerical overflow ! \n"); */
/* 				  /\*  			  p_lk[site][catg][i] = p_lk[site][catg-1][i] ; *\/ */
/* 				}			   */
			    }
			}
		    }
		  n_scale_f[site] += (fit_double)log(max_p_lk);
		}
 
	      if(max_p_lk > (1./LIM_SCALE_VAL))
		{
		  For(catq,tree->mod->n_catq)
		    {
		      For(catg,tree->mod->n_catg)
			{
			  For(k,tree->mod->ns) 
			    {
			      p_lk[site][catq][catg][k] /= max_p_lk;
				
/* 			      if(p_lk[site][catq][catg][k] < LIM_SCALE_VAL) */
/* 				{ */
/* 				  printf("\n. Numerical overflow ! \n"); */
/* 				  /\*  			  p_lk[site][catg][i] = p_lk[site][catg-1][i] ; *\/ */
/* 				}			   */
			    }
			}
		    }
		  n_scale_f[site] += (fit_double)log(max_p_lk);
		}
	    }
	}
    }
}


/*********************************************************/

void Make_Tree_4_Lk(arbre *tree, allseq *alldata, int n_site)
{
  /* Prepare a tree for a likelihood analysis */

  int i,n;
/*   int k,curr_site; */
/*   int catq, catg; */

  tree->n_pattern         = tree->data->crunch_len/tree->mod->stepsize;

  if(!tree->site_loglk_sorted) tree->site_loglk_sorted = (fit_double *)mCalloc(tree->n_pattern, sizeof(fit_double));
  if(!tree->site_lk)           tree->site_lk           = (fit_double *)mCalloc(alldata->crunch_len,sizeof(fit_double));
  if(!tree->w_patt)            tree->w_patt            = (int *)mCalloc(alldata->crunch_len,sizeof(int));
  
  tree->loglk_cat = (fit_double *)mCalloc(3,sizeof(fit_double ));
  

  tree->sel_regime_prob = (fit_double **)mCalloc(MMAX(tree->mod->n_omega,tree->mod->n_catq),sizeof(fit_double *));
  For(i,MMAX(tree->mod->n_omega,tree->mod->n_catq))
      tree->sel_regime_prob[i] = (fit_double *)mCalloc(tree->n_pattern,sizeof(fit_double ));

  tree->n_changes_per_site = (int *)mCalloc(tree->n_pattern,sizeof(int));


  Check_Memory_Amount(tree);

  printf("\n. Allocating memory...\n");

  For(i,2*tree->n_otu-3)
    {
      Make_Edge_Lk(tree->t_edges[i]->left,
		   tree->t_edges[i]->rght,
		   tree->t_edges[i],
		   tree);
    }


  For(i,2*tree->n_otu-2)
    Make_Node_Lk(tree,tree->noeud[i]);


  Alloc_All_P_Lk(tree);

  n=0;
#if defined(POSSELSITEID)
  Fors(i,tree->data->crunch_len,tree->mod->stepsize)
    {
      tree->w_patt[n] = 1;
      n++;
    }
#else
  Fors(i,tree->data->crunch_len,tree->mod->stepsize)
    {
      tree->w_patt[n] = (int)tree->data->wght[i];
      n++;
    }
#endif

#if defined(FITMODEL)

  Init_Tips(tree);

#endif

  Allocate_Num_Parameters(tree);
}

/*********************************************************/

void Init_Tips(arbre *tree)
{
  int curr_site,i,j,k,catg,catq,stop;
  
  stop = 0;
  Fors(curr_site,tree->data->crunch_len,tree->mod->stepsize)
    {
      For(i,tree->n_otu)
	{
	  For(j,tree->mod->ns) 
	    {
              /* printf("\n. i=%d j=%d %p %p %p", */
              /*        i,j, */
              /*        tree->noeud[i], */
              /*        tree->noeud[i]->b[0], */
              /*        tree->noeud[i]->b[0]->p_lk_rght); */
              /* fflush(NULL); */
              tree->noeud[i]->b[0]->p_lk_rght[curr_site/tree->mod->stepsize][0][0][j] = 0.0;
            }

	  switch(tree->mod->datatype)
	    {
	    case NT : 
	      {
		switch(tree->mod->model_applies_to)
		  {
		  case NT :
		    {
		      Init_Tips_At_One_Site_Nucleotides(tree->data->c_seq[i]->state[curr_site],
							&tree->noeud[i]->b[0]->p_lk_rght[curr_site][0][0]);
		      break;
		    }
		  case CODONS :
		    {			  
		      Init_Tips_At_One_Site_Codon(tree->data->c_seq[i]->state[curr_site],
						  tree->data->c_seq[i]->state[curr_site+1],
						  tree->data->c_seq[i]->state[curr_site+2],
						  &tree->noeud[i]->b[0]->p_lk_rght[curr_site/tree->mod->stepsize][0][0],
						  tree->mod->c_code,
						  tree->mod->ns,&stop);
		      if(stop)
			{
			  printf("\n. Stop codon found in sequence %s at position %d\n",
				 tree->data->c_seq[i]->name,
				 curr_site*tree->mod->stepsize);
			  /* 					    tree->data->c_seq[i]->state[curr_site] = '-'; */
			  /* 					    tree->data->c_seq[i]->state[curr_site+1] = '-'; */
			  /* 					    tree->data->c_seq[i]->state[curr_site+2] = '-'; */
			}
		      break;
		    }
		  default : 
		    {
		      printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
		      Exit("");
		    }
		  }
		break;
	      }
	    case AA : 
	      {
		Init_Tips_At_One_Site_AA(tree->data->c_seq[i]->state[curr_site],
					 &tree->noeud[i]->b[0]->p_lk_rght[curr_site][0][0]);
		
		
		break;
	      }
	    default : 
	      {
		printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
		Exit("");
	      }	      
	    }

	  if((tree->noeud[i]->b[0]->p_lk_rght) && (tree->noeud[i]->b[0]->get_p_lk_rght))
	    {
	      For(catq,tree->mod->n_catq)
		{
		  For(catg,tree->mod->n_catg)
		    {
		      For(k,tree->mod->ns)
			{
			  tree->noeud[i]->b[0]->p_lk_rght[curr_site/tree->mod->stepsize][catq][catg][k]=
			    tree->noeud[i]->b[0]->p_lk_rght[curr_site/tree->mod->stepsize][0][0][k];
			}
		    }
		}	      
	    }
	  else
	    {
	      For(catq,tree->mod->n_catq)
		{
		  For(catg,tree->mod->n_catg)
		    For(k,tree->mod->ns)
		    {
		      tree->noeud[i]->b[0]->p_lk_rght[curr_site/tree->mod->stepsize][catq][catg][k]=
			tree->noeud[i]->b[0]->p_lk_rght[curr_site/tree->mod->stepsize][0][0][k];
		    }
		}
	    }
	}
    }	    
}


/*********************************************************/

void Get_Infered_States(edge *b_fcus, node *a, node *d, arbre *tree, int site_num)
{
  /* Estimate posterior probabilities of states at a given node */

  int i,k,l;
  int catq,catg,state_proba;
  fit_double ****p_lk,****p_lk_opposite,buff,*post_proba;
  fit_double joint_lk,lk,proba_max;
  fit_double sum;
  int ns;

  if(d->tax) 
    {
      p_lk = (d == b_fcus->left)?(b_fcus->p_lk_left):(b_fcus->p_lk_rght);
      For(i,tree->mod->ns)
	{
	  if(p_lk[site_num][0][0][i] > 0.0) 
	    {
	      d->infered_state = i;
	      break;
	    }
	}
      return;
    }
  else
    {      
      if(a->tax)
	{
	  p_lk = (a == b_fcus->left)?(b_fcus->p_lk_left):(b_fcus->p_lk_rght);
	  For(i,tree->mod->ns)
	    {
	      if(p_lk[site_num][0][0][i] > 0.0) 
		{
		  a->infered_state = i;
		  break;
		}
	    }
	}
      For(i,3) if(d->v[i] != a) Get_Infered_States(d->b[i],d,d->v[i],tree,site_num);

      if(d == b_fcus->left)
	{
	  p_lk = b_fcus->p_lk_left;
	  p_lk_opposite = b_fcus->p_lk_rght;
	}
      else
	{
	  p_lk = b_fcus->p_lk_rght;
	  p_lk_opposite = b_fcus->p_lk_left;
	}

      post_proba = (fit_double *)mCalloc(tree->mod->ns_codon,sizeof(fit_double));
 
      lk = buff = .0;
      
      For(catq,tree->mod->n_catq) 
	{
	  For(catg,tree->mod->n_catg)
	    {
	      For(k,tree->mod->ns)
		{
		  For(l,tree->mod->ns)
		    {			  
		      lk += 
			b_fcus->qmat_struct->qmat_proba[(int)k/tree->mod->ns_codon+catq] *
			tree->mod->r_proba[catg] *
			tree->mod->pi[k%tree->mod->ns_codon] * 
			p_lk[site_num][catq][catg][k] *
			p_lk_opposite[site_num][catq][catg][l] *
			b_fcus->Pij_rr[catq][catg][k][l];				   
		    }
		}
	    }
	}


      buff = (fit_double)log(lk);
      buff += b_fcus->sum_scale_f_left[site_num];
      buff += b_fcus->sum_scale_f_rght[site_num];
      lk = (fit_double)exp(buff);
      
      
      sum = 0.0;
      joint_lk = 0.0;
      For(state_proba,tree->mod->ns)
	{
	  joint_lk = .0;
	  For(catq,tree->mod->n_catq) 
	    {
	      For(catg,tree->mod->n_catg)
		{
		  For(l,tree->mod->ns)
		    {			  
		      joint_lk += 
			b_fcus->qmat_struct->qmat_proba[(int)state_proba/tree->mod->ns_codon+catq] *
			tree->mod->r_proba[catg] *
			p_lk[site_num][catq][catg][state_proba] *
			p_lk_opposite[site_num][catq][catg][l] *
			b_fcus->Pij_rr[catq][catg][state_proba][l];
		    }
		}
	    }

	  buff = (fit_double)log(joint_lk);
	  buff += b_fcus->sum_scale_f_left[site_num];
	  buff += b_fcus->sum_scale_f_rght[site_num];
	  joint_lk = (fit_double)exp(buff);
	  

	  if(tree->mod->model_number > 22)
	    post_proba[state_proba%tree->mod->ns_codon] += 
	      joint_lk * 
	      tree->mod->pi[state_proba%tree->mod->ns_codon] /
	      lk;        
	  else
	    post_proba[state_proba%tree->mod->ns_codon] += 
	      joint_lk * 
	      tree->mod->pi[state_proba%tree->mod->ns_codon] /
	      lk;        
	}

      proba_max = 0.0;
      sum = 0.0;
      ns = (tree->mod->model_number < 20)?(tree->mod->ns):(tree->mod->ns_codon);
      For(state_proba,ns)
	{
	  sum += post_proba[state_proba];
	  if(post_proba[state_proba] > proba_max)
	    {
	      d->infered_state = state_proba;
	      proba_max = post_proba[state_proba];
	    }
	}
      
      if((sum < 1-1.E-06) || (sum > 1+1.E-06)) Exit("\n. Problem in Get_Infered_States\n");
      Free(post_proba);
    }
}

/*********************************************************/

void Get_PMat(arbre *tree, fit_double len, fit_double ****Pij, qmat *qmat_struct)
{
  /* Update a single Pij */
  
  int catg, catq;
  fit_double init_len;

  init_len = len;

  For(catg,tree->mod->n_catg)
    {
      len = init_len*tree->mod->rr_mixturem[catg];
	
      if(len < BL_MIN) len = BL_MIN;
      if(len > BL_MAX) len = BL_MAX;
      
      For(catq,tree->mod->n_catq)
	{
	  if(tree->mod->analytical)
	    {
	      PMat(len,
		   tree->mod,
		   Pij[catq][catg],
		   NULL,tree);
	    }
	  else
	    {
	      qmat_struct->curr_qmat_cat = catq;
	      PMat(len,
		   tree->mod,
		   Pij[catq][catg],
		   qmat_struct,tree);
	    }
	}
    }
}

/*********************************************************/

void Get_All_PMat_Post(node *a, node *d, edge *b, arbre *tree)
{
  
  /* Update Pij on branch 'b' */

  int i;
      
  For(i,tree->n_pattern)
    {
      b->sum_scale_f_rght[i] = .0;
      b->sum_scale_f_left[i] = .0;
    }

  if(b->l < BL_MIN) b->l = BL_MIN;
  if(b->l > BL_MAX) b->l = BL_MAX;
  
  Get_PMat(tree,
	   b->l,
	   b->Pij_rr,
	   b->qmat_struct);
  

  if(d->tax) return;
  else
    For(i,3) if(d->v[i] != a) Get_All_PMat_Post(d,d->v[i],d->b[i],tree);
}

/*********************************************************/


fit_double  Select_Regime_Change_Proba_S1(int beg, int end, edge *b, fit_double length, arbre *tree)
{
  /* Return the probability of change from selection regime 'beg' to 'end'
     on branch 'b'. Switches between selection regimes follow S1 
     (see Guindon et al., 2004, PNAS) 
  */


  if(beg == end)
    {
      return
	b->qmat_struct->omega_proba[beg] * 
	(1. - (fit_double)exp(-(length*b->qmat_struct->theta[0]))) +
	 (fit_double)exp(-(length*b->qmat_struct->theta[0]));
    }
  else
    {
      return
	b->qmat_struct->omega_proba[end] * 
	(1. - (fit_double)exp(-(length*b->qmat_struct->theta[0])));
    }
}

/*********************************************************/

fit_double Unscaled_Dwell_S1(int beg, int mid, int end, edge *b, arbre *tree)
{
  /* Calculation of (dwell times in 'mid') x (proba of change from 'beg' to 'end'). 
     Switches between selection regimes follow S1 (see Guindon et al., 2004, PNAS) 
  */

  fit_double tta;
  fit_double bl;
  fit_double pmid, pend;
  
  bl = b->l;
  tta = b->qmat_struct->theta[0];
  pmid = b->qmat_struct->omega_proba[mid];
  pend = b->qmat_struct->omega_proba[end];
  if((beg == mid) && (mid == end))
    {
      return
	(pmid*(fit_double)exp(tta*bl) + pmid*pend*bl*tta*(fit_double)exp(tta*bl) - 2.*pmid*pend*(fit_double)exp(tta*bl) - 
	 pmid*bl*tta - pend + pmid*pend*tta*bl + 2.*pmid*pend - pend*bl*tta + tta*bl - 
	 pmid + pend*(fit_double)exp(tta*bl))/(tta*(fit_double)exp(tta*bl));
    }

  if((beg == mid) && (mid != end))
    {
      return
	pend*(2.*pmid + pmid*bl*tta*(fit_double)exp(tta*bl) - 2.*pmid*(fit_double)exp(tta*bl) + pmid*bl*tta - 1. 
	      - tta*bl + (fit_double)exp(tta*bl))/(tta*(fit_double)exp(tta*bl));
    }

  if((beg != mid) && (mid == end))
    {

      return
	pmid*((fit_double)exp(tta*bl) + pend*bl*tta*(fit_double)exp(tta*bl) - 2.*pend*(fit_double)exp(tta*bl) - tta*bl + 
	      pend*bl*tta + 2.*pend - 1.)/(tta*(fit_double)exp(tta*bl));
    }
  
  if((beg != mid) && (mid != end))
    {
      return
	(pmid*pend*(2. + bl*tta*(fit_double)exp(tta*bl) - 2.*(fit_double)exp(tta*bl) + tta*bl))/(tta*(fit_double)exp(tta*bl));
    }
  else Exit(". Unexamined failure !\n");

  return -1;
}

/*********************************************************/

void Compute_All_Dwell_Probs_S2(edge *b_fcus,qmat *qmat_struct,arbre *tree)
{
  int step;
  int i,j,k,l;
  fit_double **P1,**P2;

  P1 = (fit_double **)mCalloc(tree->mod->n_omega,sizeof(fit_double *));
  For(i,tree->mod->n_omega)
      P1[i] = (fit_double *)mCalloc(tree->mod->n_omega,sizeof(fit_double ));
  P2 = (fit_double **)mCalloc(tree->mod->n_omega,sizeof(fit_double *));
  For(i,tree->mod->n_omega)
      P2[i] = (fit_double *)mCalloc(tree->mod->n_omega,sizeof(fit_double ));
  
  For(j,tree->mod->n_omega)
      {
          For(k,tree->mod->n_omega)
              {
                  For(l,tree->mod->n_omega)
                      {
                          b_fcus->dwell_probs[j][k][l] = 0.0;
                      }
              }
      }


  step = 10000;

  For(i,step)
      {
          Compute_Pmat_S2(i*b_fcus->l/(fit_double)step,P1,qmat_struct,tree);
          Compute_Pmat_S2((step-i)*b_fcus->l/(fit_double)step,P2,qmat_struct,tree);
          
          For(j,tree->mod->n_omega)
              {
                  For(k,tree->mod->n_omega)
                      {
                          For(l,tree->mod->n_omega)
                              {
                                  b_fcus->dwell_probs[j][k][l] += P1[j][k]*P2[k][l] * b_fcus->l / (fit_double)step;
                              }
                      }
              }
      }

  Compute_Pmat_S2(b_fcus->l,P1,qmat_struct,tree);
  
  For(j,tree->mod->n_omega)
      {
          For(k,tree->mod->n_omega)
              {
                  For(l,tree->mod->n_omega)
                      {
                          b_fcus->dwell_probs[j][k][l] /= P1[j][l];
                      }
              }
      }

  For(i,tree->mod->n_omega) Free(P1[i]); Free(P1);
  For(i,tree->mod->n_omega) Free(P2[i]); Free(P2);
}

/*********************************************************/

void Compute_Pmat_S2(fit_double l, fit_double **Pij, qmat *qmat_struct, arbre *tree)
{
    int i,j,k;
    int n,neg;
    fit_double *U, *V, *R;
    fit_double *expt, *uexpt;

    n = tree->mod->n_omega;
    
    expt  = (fit_double *)mCalloc(n,  sizeof(fit_double));
    uexpt = (fit_double *)mCalloc(n*n,sizeof(fit_double));
    

    U = qmat_struct->u_mat;
    V = qmat_struct->v_mat;
    R = qmat_struct->expD_mr_vct;


    /* Compute pow((fit_double)exp(D/mr),l) into mat_eDmrl */
    For(i,n)  expt[i] = (fit_double)pow(R[i],l);

    /* Multiply Vr*pow((fit_double)exp(D/mr),l)*Vi into Pij */
    For(i,n) For(j,n) uexpt[i*n+j] = U[i*n+j] * expt[j];
  
    neg = 0;
    For(i,n)
        {      
            For(j,n)
                {
                    Pij[i][j] = .0;
                    
                    For(k,n) 
                        Pij[i][j] += uexpt[i*n+k] * V[k*n+j];

                    if(Pij[i][j] < MDBL_MIN)
                        {
                            Pij[i][j] = 1.E-30;
                        }
                    if(Pij[i][j] < 0.0)
                        {
                            if(Pij[i][j] < 0.0-1.E-5) neg = 1;
                            Pij[i][j] = 1.E-60;
                        }
                }
        }
    
    
    Free(expt); Free(uexpt);
    
    if(neg) 
        {
            printf(". Err. in 'fit_double **Get_Pmat_S2(fit_double l, qmat *qmat_struct, arbre *tree)'\n");
        }

}

/*********************************************************/

/* void Eigen_S2(qmat *qmat_struct, arbre *tree) */
/* { */
/*     fit_double p0,p1,p2,alpha,beta,delta; */
/*     fit_double *Q, *U, *V, *Root; */
/*     int i; */

/*     Q    = qmat_struct->qmat; */
/*     U    = qmat_struct->u_mat; */
/*     V    = qmat_struct->v_mat; */
/*     Root = qmat_struct->root_vct; */

/*     p0 = tree->mod->qmat_struct[0]->omega_proba[0]; */
/*     p1 = tree->mod->qmat_struct[0]->omega_proba[1]; */
/*     p2 = tree->mod->qmat_struct[0]->omega_proba[2]; */

/*     alpha = tree->mod->qmat_struct[0]->theta[1]/tree->mod->qmat_struct[0]->theta[0]; */
/*     beta  = tree->mod->qmat_struct[0]->theta[2]/tree->mod->qmat_struct[0]->theta[0]; */
/*     delta = tree->mod->qmat_struct[0]->theta[0]; */
       
/*     Q[0*tree->mod->n_omega+0] = -p1 * alpha * delta - p2 * beta * delta; */
/*     Q[0*tree->mod->n_omega+1] =  p1 * alpha * delta; */
/*     Q[0*tree->mod->n_omega+2] =  p2 * beta  * delta; */
/*     Q[1*tree->mod->n_omega+0] =  p0 * alpha * delta; */
/*     Q[1*tree->mod->n_omega+1] = -p0 * alpha * delta - p2 * delta; */
/*     Q[1*tree->mod->n_omega+2] =  p2 * delta; */
/*     Q[2*tree->mod->n_omega+0] =  p0 * beta  * delta; */
/*     Q[2*tree->mod->n_omega+1] =  p1 * delta; */
/*     Q[2*tree->mod->n_omega+2] = -p0 * beta  * delta - p1 * delta; */
    
/*     Update_Eigen(tree->mod->n_omega,Q,U,V,Root,tree->mod); */
/* /\*     PMat_Numeric(b->l,tree->mod,Pij,qmat_struct); *\/ */

/*     For(i,tree->mod->n_omega) */
/*         qmat_struct->expD_mr_vct[i] = (fit_double)exp(Root[i]); */
/* } */

/*********************************************************/

void Compute_Proba_Omega_On_Edges(arbre *tree)
{
  int edgenum,selregime;
  char *s;
  fit_double ***probs,*dwell,*largest_prob;
  int *best_class;
  fit_double tbl;
  fit_double sum;
  FILE **fp_w, *fp_best, *fp;
  int i;

  fp = NULL;
  
  s = (char *)mCalloc(T_MAX_FILE,sizeof(char));
  best_class = (int *)mCalloc(2*tree->n_otu-3,sizeof(int));
  largest_prob = (fit_double *)mCalloc(2*tree->n_otu-3,sizeof(fit_double));
  fp_w = (FILE **)mCalloc(tree->mod->n_omega,sizeof(FILE *));


  For(i,tree->mod->n_omega)
    {
      strcpy(s,tree->input->seqfile);
      strcat(s,"_trees_w");
      sprintf(s+strlen(s),"%d",i+1);
      fp_w[i] = Openfile(s,1);
    }

  strcpy(s,tree->input->seqfile);
  strcat(s,"_trees_wbest");
  fp_best = Openfile(s,1);

  Free(s);

/*   qmat_struct = NULL; */
/*   if(tree->mod->s_opt->opt_theta == 3) */
/*     { */
/*       qmat_struct = Make_Qmat_Struct(tree->mod->n_omega,1,tree->mod->n_omega); */
/*       Eigen_S2(qmat_struct,tree); */
/*       For(i,2*tree->n_otu-3) */
/* 	{ */
/* 	  Compute_All_Dwell_Probs_S2(tree->t_edges[i],qmat_struct,tree); */
/* 	} */
/*       Free_Qmat(qmat_struct,tree->mod); */
/*     } */


  dwell = (fit_double *)mCalloc(tree->mod->n_omega,sizeof(fit_double));
  probs = (fit_double ***)mCalloc(2*tree->n_otu-3,sizeof(fit_double **));

  For(edgenum,2*tree->n_otu-3)
    {
      probs[edgenum] = (fit_double **)mCalloc(tree->n_pattern,sizeof(fit_double *));
      For(tree->curr_site,tree->n_pattern)
	probs[edgenum][tree->curr_site] = (fit_double *)mCalloc(tree->mod->n_omega,sizeof(fit_double));
    }


  tree->both_sides =1;
  Lk(tree);

  // Compute posterior probabilities on each edge, at each site
  For(edgenum,2*tree->n_otu-3) 
    {
      printf("\n. Edge %4d/%4d",edgenum+1,2*tree->n_otu-3);

      Integral_Term_On_One_Edge(tree->t_edges[edgenum],tree);

      For(tree->curr_site,tree->n_pattern)
	Compute_Proba_Omega_On_One_Edge_At_One_Site(tree->t_edges[edgenum],
						    probs[edgenum][tree->curr_site],
						    tree);

      Free_Integral_Term_On_One_Edge(tree->t_edges[edgenum],tree);
    }
    
  For(edgenum,2*tree->n_otu-3) 
    {
      Make_New_Edge_Label(tree->t_edges[edgenum]);
      tree->t_edges[edgenum]->n_labels++;
    }

  // Print posterior probabilities on each edge, at each site
  For(tree->curr_site,tree->n_pattern) 
    {                  
      For(edgenum,2*tree->n_otu-3) 
	{
	  largest_prob[edgenum] = -1.0;
	  best_class[edgenum] = -1;
	}

      For(selregime,tree->mod->n_omega)
      /* For(selregime,3) */
	{
	  For(edgenum,2*tree->n_otu-3) 
	    {
	      sprintf(tree->t_edges[edgenum]->labels[0],":%f",probs[edgenum][tree->curr_site][selregime]);
	      if(probs[edgenum][tree->curr_site][selregime] > largest_prob[edgenum])
		{
		  largest_prob[edgenum] = probs[edgenum][tree->curr_site][selregime];
		  best_class[edgenum] = selregime;
		}
	    }

	  fp = fp_w[selregime];

	  fprintf(fp,"tree %d = ",tree->curr_site+1);
	  s = Write_Tree(tree);
	  fprintf(fp,"%s\n",s);
	  Free(s);
	}
      
      For(edgenum,2*tree->n_otu-3) 
	{
	  switch(best_class[edgenum])
	    {
	    case 0:
	      {
		sprintf(tree->t_edges[edgenum]->labels[0],":%f",0.0);
		break;
	      }
	    case 1:
	      {
		sprintf(tree->t_edges[edgenum]->labels[0],":%f",0.5);
		break;
	      }
	    default:
	      {
		sprintf(tree->t_edges[edgenum]->labels[0],":%f",1.0);
		break;
	      }
	    }
	}

      fprintf(fp_best,"tree %d = ",tree->curr_site+1);      
      s = Write_Tree(tree);
      fprintf(fp_best,"%s\n",s);
      Free(s);
      
      For(edgenum,2*tree->n_otu-3) 
	{
	  sum = .0;
	  For(selregime,tree->mod->n_omega) sum += probs[edgenum][tree->curr_site][selregime];
	  
	  if((sum < 0.9999) || (sum > 1.0001))
	    {
	      printf("\n. sum = %f\n",sum);
	      printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	      Exit("\n");
	    }
	}

      fflush(NULL);
    }


  // Compute the posterior probabilities at each site
  tbl = .0;
/*   For(edgenum,2*tree->n_otu-3) tbl += tree->t_edges[edgenum]->l; */
  For(edgenum,2*tree->n_otu-3) tbl += 1.0;
  For(tree->curr_site,tree->n_pattern) 
    {    
      dwell[0] = dwell[1] = dwell[2] = .0;

      For(edgenum,2*tree->n_otu-3) 
	For(selregime,tree->mod->n_omega)
	dwell[selregime] += probs[edgenum][tree->curr_site][selregime];
      
      For(selregime,tree->mod->n_omega)
	tree->sel_regime_prob[selregime][tree->curr_site] = dwell[selregime]/tbl;      
    }

  Free(dwell);

  printf("\n. Computed probs at each site\n");

  For(edgenum,2*tree->n_otu-3)
    {
      For(tree->curr_site,tree->n_pattern)
	Free(probs[edgenum][tree->curr_site]);
      Free(probs[edgenum]);
    }
  Free(probs);
  Free(best_class);
  Free(largest_prob);


  For(i,tree->mod->n_omega)
    fclose(fp_w[i]);
  Free(fp_w);
  fclose(fp_best);

  fflush(NULL);
}

/*********************************************************/

/* void Compute_Proba_Omega_On_Every_Edge_At_One_Site(arbre *tree) */
/* {  */
/*     int i; */
/*     fit_double tbl; */
/*     fit_double dwell[3]; */
    
    
/*     dwell[0] = dwell[1] = dwell[2] = .0; */
/*     tbl = .0; */
/*     For(i,2*tree->n_otu-3) */
/*         { */
/*             Compute_Proba_Omega_On_One_Edge_At_One_Site(tree->t_edges[i],tree); */
            
/*             tbl += tree->t_edges[i]->l; */
            
/*             dwell[0] += tree->t_edges[i]->dwell[0]; */
/*             dwell[1] += tree->t_edges[i]->dwell[1]; */
/*             dwell[2] += tree->t_edges[i]->dwell[2]; */
/*         } */
    
    
/*     For(i,tree->mod->n_omega) */
/*       tree->sel_regime_prob[i][tree->curr_site] = dwell[i]/tbl; */
    
    
    /*   printf("%3d  %10f  %10E  %10E  %10E\n\n", */
    /* 	 tree->curr_site+1, */
    /* 	 tree->site_lk[tree->curr_site], */
    /* 	 dwell[0]/tbl, */
    /* 	 dwell[1]/tbl, */
    /* 	 dwell[2]/tbl); */
/* } */

/*********************************************************/

/* void Compute_Proba_Omega_On_One_Edge_At_One_Site(edge *b, fit_double *proba, arbre *tree) */
/* { */
/*     /\* Calculation of the expected frequencies of each of */
/*        the selection regime at a given site, based on dwell times */
/*     *\/ */
    
/*     int ofcus, oleft, orght; */
/*     int sleft, srght; */
/*     fit_double dwell_time; */
/*     fit_double change_proba; */
/*     fit_double **lk_oleft_orght; */
/*     fit_double sum_lk_oleft_orght; */
/*     int i; */
    
/*     if(!tree->w_patt[tree->curr_site]) return; */
    
/*     lk_oleft_orght = (fit_double **)mCalloc(3,sizeof(fit_double *)); */
/*     For(i,3) lk_oleft_orght[i] = (fit_double *)mCalloc(3,sizeof(fit_double)); */
    
/*     /\* WARNING : p_lks' MUST NOT BE SCALED !!! CHECK INIT_CONSTANTS *\/ */
    
    
/*     sum_lk_oleft_orght = .0; */
/*     For(sleft,tree->mod->ns) */
/*         { */
/*             For(srght,tree->mod->ns) */
/*                 { */
/*                     oleft = (int)(sleft/tree->mod->ns_codon); */
/*                     orght = (int)(srght/tree->mod->ns_codon); */
                    
/*                     lk_oleft_orght[oleft][orght] += */
/*                         b->qmat_struct->omega_proba[oleft] * */
/*                         tree->mod->pi[sleft%tree->mod->ns_codon] * */
/*                         b->p_lk_left[tree->curr_site][0][0][sleft] * */
/*                         b->p_lk_rght[tree->curr_site][0][0][srght] * */
/*                         b->Pij_rr[0][0][sleft][srght]; */
                    
/*                 } */
/*         } */
    
/*     For(oleft,3) For(orght,3) */
/*         sum_lk_oleft_orght += */
/*         lk_oleft_orght[oleft][orght]; */
    
    
/*     dwell_time = .0; */
/*     change_proba = .0; */
    
/*     switch(tree->mod->s_opt->opt_theta) */
/*       { */
/*       case 1 : */
/* 	{ */
/* 	  For(ofcus,tree->mod->n_omega) */
/* 	    { */

/* 	      proba[ofcus] = .0; */

/* 	      For(oleft,tree->mod->n_omega) */
/* 		{ */
/* 		  For(orght,tree->mod->n_omega) */
/* 		    { */
/* 		      dwell_time = Unscaled_Dwell_S1(oleft,ofcus,orght,b,tree); */
		      
/* 		      dwell_time /= Select_Regime_Change_Proba_S1(oleft,orght,b,b->l,tree); */
		      
/* 		      dwell_time /= b->l; */
		      
/* 		      proba[ofcus] += */
/* 			dwell_time * */
/* 			lk_oleft_orght[oleft][orght] / */
/* 			sum_lk_oleft_orght; */
		      
/* 		    } */
/* 		} */
/* 	    } */
/* 	  break; */
/* 	} */
/*       case 3 : */
/* 	{ */
/* 	  For(ofcus,tree->mod->n_omega) */
/* 	    { */
/* 	      For(oleft,tree->mod->n_omega) */
/* 		{ */
/* 		  For(orght,tree->mod->n_omega) */
/* 		    { */
/* 		      dwell_time = b->dwell_probs[oleft][ofcus][orght]; */
		      
/* 		      dwell_time /= b->l; */
		      
/* 		      proba[ofcus] += */
/* 			dwell_time * */
/* 			lk_oleft_orght[oleft][orght] / */
/* 			sum_lk_oleft_orght; */
		      
/* 		    } */
/* 		} */
/* 	    } */
/* 	  break; */
/* 	} */
/*       default : Exit("\n. Err. in 'void Proba_Omega_On_One_Edge_S1(edge *b, arbre *tree)'"); */
/*       } */
    
    
/*     For(i,3) */
/*       b->dwell[i] = proba[i] * b->l; */
    
    
/*     sprintf(b->labels[0],":%f",proba[2]); */
    

/*   if(((proba[0] + proba[1] + proba[2]) > 1.0+1.E-5) || */
/*      ((proba[0] + proba[1] + proba[2]) < 1.0-1.E-5)) */
/*      Exit("\n. Err in Compute_Proba_Omega_On_One_Edge_At_One_Site\n"); */

/* /\*       printf("%f %f %f\n", *\/ */
/* /\*              proba[0], *\/ */
/* /\*              proba[1], *\/ */
/* /\*              proba[2]); *\/ */
    
/*     For(i,3) Free(lk_oleft_orght[i]); */
/*     Free(lk_oleft_orght); */
/* } */

/*********************************************************/

void Integral_Term_On_One_Edge(edge *b, arbre *tree)
{

  fit_double ***integral,**P1,**P2;  
  int ns;
  int i,j,k,l;
  int step;


  ns = tree->mod->ns;


  P1 = (fit_double **)mCalloc(ns,sizeof(fit_double *));
  For(i,ns) P1[i] = (fit_double *)mCalloc(ns,sizeof(fit_double));

  P2 = (fit_double **)mCalloc(ns,sizeof(fit_double *));
  For(i,ns) P2[i] = (fit_double *)mCalloc(ns,sizeof(fit_double));



  integral = (fit_double ***)mCalloc(ns,sizeof(fit_double **));
  For(i,ns)
    {
      integral[i] = (fit_double **)mCalloc(ns,sizeof(fit_double *));
      For(j,ns) integral[i][j] = (fit_double *)mCalloc(ns,sizeof(fit_double));
    }

  b->integral = integral;

  //Integral calculation 
  step = 10;

  printf("\n. [");
  For(i,step)
    {
      PMat(((fit_double)(i+0.5)/step)*b->l,
	   tree->mod,
	   P1,
	   b->qmat_struct,tree);
      
      PMat(((fit_double)(step-i-0.5)/step)*b->l,
	   tree->mod,
	   P2,
	   b->qmat_struct,tree);
      

      For(j,ns)
	{
	  For(k,ns)
	    {
	      For(l,ns)
		{
		  integral[j][k][l] += P1[j][k] * P2[j][l]  / ((fit_double)(step));
		}
	    }
	}      

      printf("."); fflush(NULL);
    }
  printf("]\n");

  For(i,ns) Free(P1[i]); 
  Free(P1);

  For(i,ns) Free(P2[i]);
  Free(P2);

}

/*********************************************************/

void Compute_Proba_Omega_On_One_Edge_At_One_Site(edge *b, fit_double *postprob, arbre *tree)
{
    /* Calculation of the expected frequencies of each of
       the selection regime at a given site, based on dwell times
    */

  fit_double site_lk;
  int i,j,k,l;
  int nomega;
  fit_double sum;

  nomega = tree->mod->n_omega;

/*   For(i,nomega) postprob[i] = Compute_Proba_Omega_At_One_Spot_At_One_Site(b,0.5,i,tree->curr_site,tree); */

  site_lk = (fit_double)exp(tree->site_lk[tree->curr_site]);

  sum = .0;
  For(i,nomega)
    {
      postprob[i] = .0;
      For(j,tree->mod->ns_codon)
	{
	  For(k,tree->mod->ns)
	    {
	      For(l,tree->mod->ns)
		{
		    postprob[i] +=
		      tree->mod->pi[j] *
		      b->p_lk_left[tree->curr_site][0][0][k] *
		      b->p_lk_rght[tree->curr_site][0][0][l] *
		      b->integral[i*tree->mod->ns_codon+j][k][l];
		}
	    }
	}
      sum += postprob[i];
    }

  For(i,nomega) postprob[i] *= b->qmat_struct->omega_proba[i] / site_lk;
  For(i,nomega) postprob[i] *= exp(b->sum_scale_f_left[tree->curr_site] + b->sum_scale_f_rght[tree->curr_site]);

/*   printf("\n. sum = %f lk = %f",sum,site_lk); */


  For(i,tree->mod->n_omega) 
    if((postprob[i] < 0.0) || (postprob[i] > 1.0))
      {
	printf("\n. Cat : %d Prob : %f\n",i,postprob[i]);
	printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	Exit("\n");
      }

  sum = 0.0;
  For(i,tree->mod->n_omega) sum += postprob[i];


  if((sum > 1.0+1.E-5) ||
     (sum < 1.0-1.E-5))
    {
      For(i,nomega) printf("\n. %f %f %f",postprob[i],b->sum_scale_f_left[tree->curr_site],b->sum_scale_f_rght[tree->curr_site]);
      printf("\n. Sum = %f\n",sum);
      printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Exit("\n");
    }


  return;
}

/*********************************************************/
/* pos=0 -> spot is the left node; pos=1 -> spot is the right node */
fit_double Compute_Proba_Omega_At_One_Spot_At_One_Site(edge *b, fit_double pos, int wclass, int site, arbre *tree)
{
  int i,j,k;
  fit_double **P1,**P2;  
  fit_double site_lk;
  fit_double p;

  P1 = (fit_double **)mCalloc(tree->mod->ns,sizeof(fit_double *));
  For(i,tree->mod->ns) P1[i] = (fit_double *)mCalloc(tree->mod->ns,sizeof(fit_double));

  P2 = (fit_double **)mCalloc(tree->mod->ns,sizeof(fit_double *));
  For(i,tree->mod->ns) P2[i] = (fit_double *)mCalloc(tree->mod->ns,sizeof(fit_double));

  PMat(pos*b->l,     tree->mod,P1,b->qmat_struct,tree);
  PMat((1.-pos)*b->l,tree->mod,P2,b->qmat_struct,tree);

  p = 0.0;
  For(i,tree->mod->ns)
    {
      For(j,tree->mod->ns)
	{
	  For(k,tree->mod->ns)
	    {
	      if((i / tree->mod->ns_codon) == wclass)
		{
		  p +=
		    tree->mod->pi[i%tree->mod->ns_codon] *
		    b->p_lk_left[site][0][0][j] *
		    b->p_lk_rght[site][0][0][k] *
		    P1[i][j] * P2[i][k];
		}
	    }
	}
    }


  p *= b->qmat_struct->omega_proba[wclass];

  site_lk = (fit_double)exp(tree->site_lk[site]);

  p /= site_lk;

  For(i,tree->mod->ns) Free(P1[i]); 
  Free(P1);

  For(i,tree->mod->ns) Free(P2[i]);
  Free(P2);

  return p;
}

/*********************************************************/

fit_double Lk_Arg(fit_double *param_val, arbre *tree)
{    
    int i,j;
    int currently_optimized;
    int beg, end;

    tree->numerical_param->replace_num_param = 1;

    currently_optimized = tree->numerical_param->currently_opt;
            
    if(currently_optimized < 0) Exit("\n. Err. in Lk_Param_Update\n");
            
    j = 0;

    beg = tree->numerical_param->param_beg[currently_optimized];
    end = 
        tree->numerical_param->param_beg[currently_optimized] +
        tree->numerical_param->param_size[currently_optimized];
    
    for(i=beg;i<end;i++)
        {
            tree->numerical_param->param_val[i] = param_val[j];
            j++;
        }

    Lk(tree);

    tree->numerical_param->replace_num_param = 0;

    return tree->tot_loglk;
}

/*********************************************************/


