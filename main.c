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
#include "options.h"
#include "draw.h"


#ifdef FITMODEL

int main(int argc, char **argv)
{
  seq **data;
  allseq *alldata;
  option *input;
  char *s,*buff;
  arbre *tree;
  int n_data_sets;
  model *mod;
  time_t t_beg,t_end;
  div_t hour,min;


  tree = NULL;
  mod  = NULL;


  fflush(stdout);

  buff = (char *)mCalloc(T_MAX_LINE,sizeof(char));

  /* Menu and allocate memory for the substitution model */
  input = (option *)Get_Input(argc,argv);
  
  Make_Model_Complete(input->mod);
  /* The substitution model is now operational */


  /* Output files (tree and statistics) */
  input->fp_output_stat = Openfile(input->output_stat_file,input->stat_file_open_mode);
  input->fp_output_tree = Openfile(input->output_tree_file,input->tree_file_open_mode);

  n_data_sets = 0;
  do
    {
        n_data_sets++;
      
        time(&t_beg);
        
        if(n_data_sets > input->n_data_sets) 
            data = NULL;
        else                                 
            data = Get_Seq(input,0);               /* Read sequences */


        if(n_data_sets < input->start_at_data_set)
	  if(!fgets(buff,T_MAX_LINE,input->fp_tree))
	    {
	      printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	      Exit("");
	    }

        if((data) && (n_data_sets >= input->start_at_data_set))
            {
                if(n_data_sets > 1) printf(". Data set [#%d]\n",n_data_sets);
                
                printf("\n. Reading sequences...\n");
                
                alldata = Compact_Seq(data,input); /* Compress sequences */
                
                Free_Seq(data,alldata->n_otu);
                
                Make_All_Qmat_Struct(input->mod); /* Allocate memory for the 
                                                     instantaneous rate matrix 
                                                  */
                
                mod = Init_Model(alldata,input); /* Inititalise the substitution 
                                                    model 
                                                 */
                
                Check_Ambiguities(alldata,input->mod->datatype,input->mod->stepsize);
                /* identify sites with ambiguous character states */
                
                
                printf("\n. Reading user tree...\n");

                tree = Read_Tree_File(input->fp_tree); /* Read input tree */

                if(!tree) 
                    {
                        printf("\n. Missing tree for data set #%d\n",n_data_sets);
                        printf("  This data set is not analyzed.\n");
                        data = NULL;
                    }
                
                if(!tree) continue;
                
                tree->mod        = mod;
                tree->data       = alldata;
                tree->both_sides = 1;
                tree->input      = input;
                
                
                
                Order_Tree_CSeq(tree,alldata); /* Sequences are ordered with respect to 
                                                  the order in which taxa appear in the 
                                                  input tree
                                               */                

                Make_Tree_4_Lk(tree,alldata,alldata->init_len); /* Prepare the tree for 
                                                                   a likelihood analysis 
                                                                */
                
                tree->both_sides                 = 1;
                tree->mod->s_opt->opt_subst_rate = 0;
                tree->mod->tpos_ols              = 0;
                tree->mod->update_bl_using_tpos  = 0;
                tree->n_root                     = NULL;
                

                /*           fprintf(input->fp_output_stat,"Site p2 lk0 lk1 lk2 lk\n"); */
                Init_Num_Parameters(tree);
				
		Print_Model_Param(tree);
                

                if(tree->mod->s_opt->opt_param) Round_Optimize(tree,tree->data); /* Optimise parameters */
                else                            Return_Lk(tree); /* Just compute the likelihood */


                fflush(NULL);
                
                Update_BrLen_Invar(tree); /* Update branch lengths by taking into 
                                             account invariable sites */



                if(tree->mod->switch_modelname != NO_SWITCH)
                    {
                        /* Computation of the posterior probabilities of selection regimes */
                        printf("\n. Calculation of the posterior probabilities of selection regimes\n");

                        tree->both_sides = 1; 
                        Lk(tree);

			rewind(tree->input->fp_output_tree);
                        fprintf(tree->input->fp_output_tree,"begin basstrees;\n");
                        fprintf(tree->input->fp_output_tree,"dimensions ntax=%d ntrees=%d;\n",
				tree->n_otu,
				tree->n_pattern);

                        Compute_Proba_Omega_On_Edges(tree);

                        fprintf(tree->input->fp_output_tree,"end;\n");
                    }
		else
		  {
		    s = Write_Tree(tree); /* Output tree */	  
		    if(input->n_data_sets == 1) rewind(input->fp_output_tree);
		    fprintf(input->fp_output_tree,"%s\n",s);
		    Free(s);
		  }
                
                time(&t_end);
                
                hour = div(t_end-t_beg,3600);
                min  = div(t_end-t_beg,60  );
                
                min.quot -= hour.quot*60;
                  
		Print_Param(tree->input->fp_output_stat,tree);
                /* Print_Fp_Out(tree->input->fp_output_stat, t_beg, t_end, tree, input, n_data_sets); */
                
/*                 Print_CSeq(tree->input->fp_output_stat,tree->data); */

                printf("\n. Time used %dh%dm%ds\n", hour.quot,min.quot,(int)(t_end-t_beg)%60);
                printf("\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n\n");
                
                fflush(NULL);
                
                Free_Tree_Lk(tree);
                
                Free_Cseq(alldata);
                
                Free_Qmat_Struct(mod);
                
                Free_Tree(tree);
            }
    }while(data);
  

  Free_Model(mod);

  
  if(input->fp_seq )        fclose(input->fp_seq );
  if(input->fp_tree)        fclose(input->fp_tree);
  if(input->fp_output_stat) fclose(input->fp_output_stat);
  if(input->fp_output_tree) fclose(input->fp_output_tree);
  
  Free_Input(input);

  Free(buff);
    
  return 0;
}

#elif EVOLVE

int main(int argc, char **argv)
{
  seq **ref_data;
  allseq *ref_alldata,*sim_seq;
  option *input;
  arbre *tree;
  int n_data_sets;
  model *mod;
  time_t t_beg,t_end;
  div_t hour,min;
  
  tree = NULL;
  mod  = NULL;
  
  srand(time(NULL));
  
  
  fflush(stdout);
  
  /* Menu and allocate memory for the substitution model */
  input = (option *)Get_Input(argc,argv);
  
  Make_Model_Complete(input->mod);
  /* The substitution model is now operational */
    
  /* Output files (simulated sequences and statistics) */
  input->fp_output_stat = Openfile(input->output_stat_file,input->stat_file_open_mode);
  input->fp_output_tree = Openfile(input->output_tree_file,input->tree_file_open_mode);
  
  
  ref_data = Get_Seq(input,0); /* Read sequences */
  
  ref_alldata = Compact_Seq(ref_data,input); /* Compress sequences */
  
  Free_Seq(ref_data,ref_alldata->n_otu);
  
  Make_All_Qmat_Struct(input->mod); /* Allocate memory for the 
				       instantaneous rate matrix 
				    */
  
  mod = Init_Model(ref_alldata,input); /* Inititalise the substitution 
					  model 
				       */
  
  Check_Ambiguities(ref_alldata,
		    input->mod->datatype,
		    input->mod->stepsize);
  /* identify sites with ambiguous character states */
  
  printf("\n. Reading user tree...\n");
  
  tree = Read_Tree_File(input->fp_tree); /* Read input tree */
  
  tree->mod        = mod;
  tree->data       = ref_alldata;
  tree->both_sides = 1;
  tree->input      = input;
  
  if(input->user_len != -1) tree->data->init_len = tree->data->crunch_len = input->user_len;
  else                      tree->data->init_len = tree->data->crunch_len = ref_alldata->init_len;
  
  Order_Tree_CSeq(tree,ref_alldata); /* Sequences are ordered with respect to 
					the order in which taxa appear in the 
					input tree
				     */
  
  Make_Tree_4_Lk(tree,ref_alldata,ref_alldata->init_len); /* Prepare the tree for 
							     a likelihood analysis 
							  */
  
  tree->both_sides                 = 1;
  tree->mod->s_opt->opt_subst_rate = 0;
  tree->mod->tpos_ols              = 0;
  tree->mod->update_bl_using_tpos  = 0;
  tree->n_root                     = tree->noeud[0];
  
  
  n_data_sets = 0;
  do
    {
      n_data_sets++;
      
      time(&t_beg);
      
                  /* fprintf(input->fp_output_stat,"Site trueclass\n"); */
      
      sim_seq = Evolve(tree);
      
      Print_CSeq(input->fp_output_tree,sim_seq);
      
      /*             Print_Ancestral_Seq(tree->noeud[0],tree->noeud[0]->v[0],tree); */
      
      fflush(NULL);
      
      time(&t_end);
      
      hour = div(t_end-t_beg,3600);
      min  = div(t_end-t_beg,60  );
      
      min.quot -= hour.quot*60;
      
      Print_Fp_Out(input->fp_output_stat, t_beg, t_end, tree, input, n_data_sets);
      
      printf("\n. Time used %dh%dm%ds\n", hour.quot,min.quot,(int)(t_end-t_beg)%60);
      printf("\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n\n");
      
      fflush(NULL);
      
      
    }
  while(n_data_sets < input->n_data_set_asked);
  
  Free_Cseq(ref_alldata);
  Free_Tree_Lk(tree);
  Free_Qmat_Struct(mod);
  Free_Tree(tree);
  Free_Model(mod);
  
  if(input->fp_seq )        fclose(input->fp_seq );
  if(input->fp_tree)        fclose(input->fp_tree);
  if(input->fp_output_tree) fclose(input->fp_output_tree);
  if(input->fp_output_stat) fclose(input->fp_output_stat);
  
  Free_Input(input);
  
  return 0;
}

#elif ANCESTRAL

int main(int argc, char **argv)
{
    seq **ref_data;
    allseq *ref_alldata;
    option *input;
    arbre *tree;
    int n_data_sets;
    model *mod;
    time_t t_beg,t_end;
    div_t hour,min;
    
    tree = NULL;
    mod  = NULL;
    
    srand(time(NULL));

    
    fflush(stdout);
    
    /* Menu and allocate memory for the substitution model */
    input = (option *)Get_Input(argc,argv);
    
    Make_Model_Complete(input->mod);
    /* The substitution model is now operational */
    
    /* Output files (simulated sequences and statistics) */
    input->fp_output_stat = Openfile(input->output_stat_file,input->stat_file_open_mode);
    input->fp_output_tree = Openfile(input->output_tree_file,input->tree_file_open_mode);
    
    ref_data = Get_Seq(input,0); /* Read sequences */
    ref_alldata = Compact_Seq(ref_data,input); /* Compress sequences */
    
    Free_Seq(ref_data,ref_alldata->n_otu);
    Make_All_Qmat_Struct(input->mod); /* Allocate memory for the 
                                         instantaneous rate matrix 
                                      */
    
    mod = Init_Model(ref_alldata,input); /* Inititalise the substitution 
                                            model 
                                         */


    Check_Ambiguities(ref_alldata,
                      input->mod->datatype,
                      input->mod->stepsize);
    /* identify sites with ambiguous character states */
    
    printf("\n. Reading user tree...\n");
    
    tree = Read_Tree_File(input->fp_tree); /* Read input tree */
    
       
    tree->mod        = mod;
    tree->data       = ref_alldata;
    tree->both_sides = 1;
    tree->input      = input;

    if(input->user_len != -1) tree->data->init_len = tree->data->crunch_len = input->user_len;
    else                      tree->data->init_len = tree->data->crunch_len = ref_alldata->init_len;
    
    Order_Tree_CSeq(tree,ref_alldata); /* Sequences are ordered with respect to 
                                          the order in which taxa appear in the 
                                          input tree
                                       */
    
    Make_Tree_4_Lk(tree,ref_alldata,ref_alldata->init_len); /* Prepare the tree for 
                                                               a likelihood analysis 
                                                            */
    
    tree->both_sides                 = 1;
    tree->mod->s_opt->opt_subst_rate = 0;
    tree->mod->tpos_ols              = 0;
    tree->mod->update_bl_using_tpos  = 0;
    tree->root                       = tree->noeud[0];
    
    n_data_sets = 0;
    do
        {
            n_data_sets++;
            
            time(&t_beg);
            
/*             fprintf(input->fp_output_stat,"Site trueclass\n"); */

            Estimate_Ancestral_States(tree);
            
            Print_CSeq(input->fp_output_tree,tree->data);
            
            Print_Ancestral_Seq(tree->noeud[0],tree->noeud[0]->v[0],tree);
            
            fflush(NULL);
            
            time(&t_end);
      
            hour = div(t_end-t_beg,3600);
            min  = div(t_end-t_beg,60  );
      
            min.quot -= hour.quot*60;
      
/*             Print_Fp_Out(input->fp_output_stat, t_beg, t_end, tree, input, n_data_sets); */
            
            printf("\n. Time used %dh%dm%ds\n", hour.quot,min.quot,(int)(t_end-t_beg)%60);
            printf("\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n\n");
      
            fflush(NULL);
      
            
        }
    while(n_data_sets < input->n_data_set_asked);

 
    Free_Cseq(ref_alldata);
    
    Free_Tree_Lk(tree);

    Free_Qmat_Struct(mod);
    
    Free_Tree(tree);
            
    Free_Model(mod);
  
  
    if(input->fp_seq )        fclose(input->fp_seq );
    if(input->fp_tree)        fclose(input->fp_tree);
    if(input->fp_output_tree) fclose(input->fp_output_tree);
    if(input->fp_output_stat) fclose(input->fp_output_stat);

    Free_Input(input);
    
    return 0;
}


#elif RF


int main(int argc, char **argv)
{
    arbre *tree1,*tree2;
    FILE *fp_tree1, *fp_tree2;
    int n_tree1, n_tree2;
    
    srand(time(NULL));

        
    fp_tree1 = Openfile(argv[1],0);
    
    fp_tree2 = Openfile(argv[2],0);
        
    n_tree1 = n_tree2 = 0;

    do
        {
            n_tree1++;
            n_tree2 = 0;
            rewind(fp_tree2);

            tree1 = Read_Tree_File(fp_tree1); 
                        
            if(!tree1) break;


            tree1->root = tree1->noeud[0];

            Alloc_Bip(tree1);
  
            Get_Bip(tree1->noeud[0],
                    tree1->noeud[0]->v[0],
                    tree1);

            do
                {
                    
                    n_tree2++;

                    tree2 = Read_Tree_File(fp_tree2);                     

                    if(!tree2) break;

                    tree2->root = tree2->noeud[0];
                    
                    Alloc_Bip(tree2);
                    
                    Get_Bip(tree2->noeud[0],
                            tree2->noeud[0]->v[0],
                            tree2);

                    printf("Tree1 num %4d, tree2 num %4d ; ",n_tree1,n_tree2);

                    printf("RF=%f\n",Rf(tree1,tree2) );

                    Free_Tree(tree2);

                }while(1);
 

            Free_Tree(tree1);
        
        }while(1);
    

    fclose(fp_tree1);
    fclose(fp_tree2);

    return 0;

}

#elif POSSELSITEID

int main(int argc, char **argv)
{
   seq **ref_data;
   allseq *ref_alldata;
   option *input;
   arbre *tree;
   model *mod;
   int n_sites,site,i,j;
   fit_double *sortedprob,buff,*orig_site_cat_prob;
   int orig_n_sites;
   fit_double xtop, xbot, thresh;
   fit_double p_thresh;
   fit_double *fdr,fixed_fdr;
   int *siteclass;
   fit_double *orig_pi;

   tree = NULL;
   mod  = NULL;

   srand(time(NULL));

   #undef  LIM_SCALE
   #define LIM_SCALE     300
   #undef  LIM_SCALE_VAL
   #define LIM_SCALE_VAL 1.E-500


   fflush(stdout);

   /* Menu and allocate memory for the substitution model */
   input = (option *)Get_Input(argc,argv);

   Make_Model_Complete(input->mod);
   /* The substitution model is now operational */

   /* Output files (simulated sequences and statistics) */
   input->fp_output_stat = Openfile(input->output_stat_file,input->stat_file_open_mode);
   input->fp_output_tree = Openfile(input->output_tree_file,input->tree_file_open_mode);

   ref_data = Get_Seq(input,0);               /* Read sequences     */
   ref_alldata = Compact_Seq(ref_data,input); /* Compress sequences */

   Free_Seq(ref_data,ref_alldata->n_otu);
   Make_All_Qmat_Struct(input->mod); /* Allocate memory for the 
   instantaneous rate matrix 
   */

   mod = Init_Model(ref_alldata,input); /* Inititalise the substitution 
   model 
   */

   Check_Ambiguities(ref_alldata,
   input->mod->datatype,
   input->mod->stepsize);  /* identify sites with ambiguous character states */

   printf("\n. Reading user tree...\n");

   tree = Read_Tree_File(input->fp_tree); /* Read input tree */

   tree->mod        = mod;
   tree->data       = ref_alldata;
   tree->both_sides = 1;
   tree->input      = input;

   orig_n_sites = tree->data->crunch_len/tree->mod->stepsize;

   if(input->user_len != -1)
       {
           tree->data->init_len =
               (input->user_len > ref_alldata->init_len)?
               (input->user_len):
               (ref_alldata->init_len);
           
           tree->data->crunch_len = 
               (input->user_len > ref_alldata->init_len)?
               (input->user_len):
               (ref_alldata->init_len);
       }
   else
       {
           tree->data->init_len   = ref_alldata->init_len;
           tree->data->crunch_len = ref_alldata->init_len;
       }
   
   
   Order_Tree_CSeq(tree,ref_alldata); /* Sequences are ordered with respect to 
                                         the order in which taxa appear in the 
                                         input tree
                                      */
   
   Make_Tree_4_Lk(tree,ref_alldata,ref_alldata->init_len); /* Prepare the tree for 
                                                              a likelihood analysis 
                                                           */
   
   tree->both_sides                 = 1;
   tree->mod->s_opt->opt_subst_rate = 0;
   tree->mod->tpos_ols              = 0;
   tree->mod->update_bl_using_tpos  = 0;
   tree->root                       = tree->noeud[0];

   orig_site_cat_prob = (fit_double *)mCalloc(orig_n_sites,sizeof(fit_double));
   orig_pi            = (fit_double *)mCalloc(tree->mod->ns,sizeof(fit_double));

   tree->n_pattern        = orig_n_sites;
   tree->data->crunch_len = orig_n_sites * tree->mod->stepsize;

   For(i,tree->mod->ns) orig_pi[i] = tree->mod->pi[i];

   Init_Tips(tree);

   Return_Lk(tree);

   For(site,orig_n_sites)
       {
           orig_site_cat_prob[site] = tree->sel_regime_prob[2][site];
       }


   if(mod->omega[2] < 1.0) 
       {
           For(site,orig_n_sites) 
               { 
                   fprintf(input->fp_output_stat,"%15d     0\n",0); 
                   fprintf(input->fp_output_tree,"0\n");
               }
           return 0;
       }

   if((input->pos_sel_sites_id_method == ADAPTIVE_CONTROL) ||
      (input->pos_sel_sites_id_method == PARAMETRIC_BOOTSTRAP))
       {
           if(input->user_len != -1)
               {
                   tree->n_pattern        = input->user_len/tree->mod->stepsize;
                   tree->data->crunch_len = input->user_len;
               }
           else
               {
                   tree->n_pattern        = ref_alldata->init_len/tree->mod->stepsize;
                   tree->data->crunch_len = ref_alldata->init_len;
               }
         
       }

   sortedprob = (fit_double *)mCalloc(tree->n_pattern,sizeof(fit_double));
   siteclass  = (int *)mCalloc(tree->n_pattern,sizeof(int));

   n_sites = -1;

   if(input->pos_sel_sites_id_method == FIXED_REGION)
       {

           fdr = (fit_double *)mCalloc(orig_n_sites,sizeof(fit_double));

           /* classify sites */
           
           For(site,orig_n_sites) 
               { 
                   if(orig_site_cat_prob[site] < mod->thresholdp2)
                       fprintf(input->fp_output_stat,"%15f     0\n",
                               orig_site_cat_prob[site]); 
                   else 
                       fprintf(input->fp_output_stat,"%15f     1\n",
                               orig_site_cat_prob[site]); 
               }
           
           
           /* estimate the FDR */
           /* Bubble sort of the site/cat proba values in decreasing order */

           For(site,orig_n_sites)
               {
                   sortedprob[site] = tree->sel_regime_prob[2][site];
                   siteclass[site]  = tree->data->selclass[site];
               }

           For(i,orig_n_sites-1)
               {
                   for(j=i+1;j<orig_n_sites;j++)
                       {
                           if(sortedprob[i] < sortedprob[j])
                               {
                                   buff          = sortedprob[j];
                                   sortedprob[j] = sortedprob[i];
                                   sortedprob[i] = buff;
                               }
                       }
               }
           
           fdr[0] = (1.-sortedprob[0]);
           for(i=1;i<orig_n_sites;i++) 
               {
                   fdr[i] = (fit_double)((fit_double)fdr[i-1]*(i) + (1.-sortedprob[i]))/(i+1);
               }

           For(i,orig_n_sites) if(sortedprob[i] < mod->thresholdp2) break;
           
           fixed_fdr = -1.;
           if(i == 0) fixed_fdr = fdr[0];
           else
               {

                   fixed_fdr = 
                       (mod->thresholdp2 - sortedprob[i]) / (sortedprob[i-1] - sortedprob[i]) * 
                       (fdr[i-1] - fdr[i]) + fdr[i];
               }
           
           fprintf(input->fp_output_tree,"%15f\n",fixed_fdr);
           Free(fdr);
       }


   if(input->pos_sel_sites_id_method == ADAPTIVE_CONTROL)
       {
           
           Evolve_Under_H0(tree);
           
           Print_CSeq(input->fp_output_tree,tree->data);
           
           Init_Tips(tree);
           
           Init_Model(tree->data, input);
           
           Return_Lk(tree);           
       
           n_sites =  tree->n_pattern;
           
           For(site,n_sites)
               {
                   sortedprob[site] = tree->sel_regime_prob[2][site];
                   siteclass[site]  = tree->data->selclass[site];
               }
           
           thresh = MDBL_MAX;
           
           /* Bubble sort of the site/cat proba values in increasing order */
           
           For(i,n_sites-1)
               {
                   for(j=i+1;j<n_sites;j++)
                       {
                           if(sortedprob[i] > sortedprob[j])
                               {
                                   buff          = sortedprob[j];
                                   sortedprob[j] = sortedprob[i];
                                   sortedprob[i] = buff;
                                   
                               }
                       }
               }
           
           
           p_thresh = 
               tree->mod->fdr *
               (mod->qmat_struct[0]->qmat_proba[0]+
                mod->qmat_struct[0]->qmat_proba[1]);
           
           xtop = (fit_double)ceil (tree->n_pattern * p_thresh);
           xbot = (fit_double)floor(tree->n_pattern * p_thresh);
           
           /* linear extrapolation */
           thresh = 
               sortedprob[(int)xbot] +
               (sortedprob[(int)xtop] - sortedprob[(int)xbot]) *
               (tree->n_pattern * p_thresh - (int)xbot)/((int)xtop - (int)xbot);
           
           
           For(site,orig_n_sites)
               fprintf(input->fp_output_stat,"%15E     1\n",
                       orig_site_cat_prob[site]);
           
           For(site,n_sites)
               fprintf(input->fp_output_stat,"%15E     0\n",
                       sortedprob[site]);
           
       }
   else if(input->pos_sel_sites_id_method == PARAMETRIC_BOOTSTRAP)
       {
           fit_double num_pfp,denom_pfp;
           fit_double *pfp,*thresharray;
           int n_run,n_tot_run;

           n_tot_run = 10;
           n_sites = orig_n_sites;

           pfp = (fit_double *)mCalloc(n_sites,sizeof(fit_double));
           thresharray = (fit_double *)mCalloc(n_tot_run,sizeof(fit_double));

           n_run = 0;
           while(n_run < n_tot_run)
               {
                   For(i,tree->mod->ns) tree->mod->pi[i] = orig_pi[i];
                   
                   Evolve(tree);
                                      
                   Init_Tips(tree);
           
                   Init_Model(tree->data, input);
           
                   Return_Lk(tree);
                   
                   
                   /* Bubble sort of the site/cat proba values in increasing order */
           
                   For(site,n_sites)
                       {
                           sortedprob[site] = tree->sel_regime_prob[2][site];
                           siteclass[site]  = tree->data->selclass[site];
                       }

                   For(i,n_sites-1)
                       {
                           for(j=i+1;j<n_sites;j++)
                               {
                                   if(sortedprob[i] > sortedprob[j])
                                       {
                                           buff          = sortedprob[j];
                                           sortedprob[j] = sortedprob[i];
                                           sortedprob[i] = buff;

                                           buff          = siteclass[j];
                                           siteclass[j]  = siteclass[i];
                                           siteclass[i]  = (int) buff;
                                       }
                               }
                       }


                   num_pfp = denom_pfp = 0.0;
                   For(i,n_sites)
                       {
                           denom_pfp = n_sites - i;
			   num_pfp = 0;
			   for(j=i;j<n_sites;j++)
			     {
			       if(siteclass[j] != 2) num_pfp += 1;
			     }
                           pfp[i] = (fit_double)num_pfp/(fit_double)denom_pfp; 
                       }

                   For(i,n_sites)
                       {
                           if(pfp[n_sites-1-i] > tree->mod->fdr) break;
                       }
                   if(i) i-=1;
                   thresharray[n_run] = sortedprob[n_sites-1-i];

/*                    For(i,n_sites) */
/*                        { */
/*                            if(pfp[i] < tree->mod->fdr) break; */
/*                        } */
/*                    if(i==n_sites) i-=1; */
/*                    thresharray[n_run] = sortedprob[i]; */
           
                   fprintf(stderr,"RUN %3d THRESH %10f %d\n",n_run+1,thresharray[n_run],i);
                   n_run++;
               }


           // Bubble sort of the threshold values

	   thresh = 0.0;
           For(i,n_tot_run-1)
               {
                   for(j=i+1;j<n_tot_run;j++)
                       {
                           if(thresharray[i] > thresharray[j])
                               {
                                   buff           = thresharray[j];
                                   thresharray[j] = thresharray[i];
                                   thresharray[i] = buff;
                               }
                       }
               }
           
	   thresh = 0.0;
           For(i,n_tot_run)
               {
                   thresh += thresharray[i];
               }

           // thresh = thresharray[(n_tot_run-1)/2];
	   thresh = thresh / (fit_double)(n_tot_run);
	   fprintf(stderr,"THRESH = %10f\n",thresh);


           For(site,orig_n_sites)
               {
                   if(orig_site_cat_prob[site] > thresh) 
                       fprintf(input->fp_output_stat,"%15E     1\n",
                               orig_site_cat_prob[site]);
                   else 
                       fprintf(input->fp_output_stat,"%15E     0\n",
                               orig_site_cat_prob[site]);
               }

           Free(pfp);
           Free(thresharray);
       }

   else if(input->pos_sel_sites_id_method == NEWTON_ETAL)
       {
           n_sites = orig_n_sites;

           fdr = (fit_double *)mCalloc(n_sites,sizeof(fit_double));
           
           /* Bubble sort of the site/cat proba values in decreasing order */
           
           For(site,n_sites)
               {
                   sortedprob[site] = tree->sel_regime_prob[2][site];
                   siteclass[site]  = tree->data->selclass[site];
               }

           For(i,n_sites-1)
               {
                   for(j=i+1;j<n_sites;j++)
                       {
                           if(sortedprob[i] < sortedprob[j])
                               {
                                   buff          = sortedprob[j];
                                   sortedprob[j] = sortedprob[i];
                                   sortedprob[i] = buff;
                               }
                       }
               }
           
           fdr[0] = (1.-sortedprob[0]);
           for(i=1;i<n_sites;i++)
               {
                   fdr[i] = (fit_double)((fit_double)fdr[i-1]*(i) + (1.-sortedprob[i]))/(i+1);
                   if((fit_double)fdr[i] > tree->mod->fdr) break;
               }
           
           /* linear extrapolation */
           if((i == 1) || (i == n_sites))
               thresh = 1.0;
           else
               {
                   thresh =
                       (1.-sortedprob[i-1]) +
                       ((1.-sortedprob[i]) - (1.-sortedprob[i-1])) *
                       (tree->mod->fdr - (fit_double)fdr[i-1])/((fit_double)fdr[i] - (fit_double)fdr[i-1]);
                   thresh = 1. - thresh;
               }
           
           printf("THRESH = %f\n",thresh);
           
           For(site,orig_n_sites)
               {
                   if(orig_site_cat_prob[site] > thresh) 
                       fprintf(input->fp_output_stat,"%15E     1\n",
                               orig_site_cat_prob[site]);
                   else 
                       fprintf(input->fp_output_stat,"%15E     0\n",
                               orig_site_cat_prob[site]);
               }

           Free(fdr);
       }
   
   
   
   Free(sortedprob);

   Free(siteclass);
   
   Free_Cseq(ref_alldata);
   
   Free_Tree_Lk(tree);
   
   Free_Qmat_Struct(mod);
   
   Free_Tree(tree);
   
   Free_Model(mod);
   
   Free(orig_pi);

   if(input->fp_seq )        fclose(input->fp_seq );
   if(input->fp_tree)        fclose(input->fp_tree);
   if(input->fp_output_tree) fclose(input->fp_output_tree);
   if(input->fp_output_stat) fclose(input->fp_output_stat);
   
   Free_Input(input);
   
   return 0;
}

#elif TOOL


int main(int argc, char **argv)
{
    FILE *fp_probs;
    int curr_site;
    int i,j;
    int n_data_sets, n_codon_sites;
    fit_double prob,*postprob,*sum,*sortedprob;
    fit_double buff;
    fit_double fdr,thresh;

    srand(time(NULL));
        
    fp_probs      = Openfile(argv[1],0);
    n_data_sets   = (int)atoi(argv[2]);    
    n_codon_sites = (int)atoi(argv[3]);
    fdr           = (fit_double)atof(argv[4]);

    postprob   = (fit_double *)mCalloc(n_codon_sites,sizeof(fit_double));
    sum        = (fit_double *)mCalloc(n_codon_sites,sizeof(fit_double));
    sortedprob = (fit_double *)mCalloc(n_codon_sites,sizeof(fit_double));

    prob = 0.0;
    thresh = -1.0;
    do
        {
            curr_site = 0;
            do
                {
                    fscanf(fp_probs,"%lf",&prob);
                    fscanf(fp_probs,"%lf",&prob);
                    fscanf(fp_probs,"%lf",&prob);
                    postprob[curr_site] = prob;
                    curr_site++;
                }
            while(curr_site < n_codon_sites);
            

           
           /* Bubble sort of the site/cat proba values in decreasing order */
           
           For(curr_site,n_codon_sites)
               {
                   sortedprob[curr_site] = postprob[curr_site];
               }

           For(i,n_codon_sites-1)
               {
                   for(j=i+1;j<n_codon_sites;j++)
                       {
                           if(sortedprob[i] < sortedprob[j])
                               {
                                   buff          = sortedprob[j];
                                   sortedprob[j] = sortedprob[i];
                                   sortedprob[i] = buff;
                               }
                       }
               }

           sum[0] = (1.-sortedprob[0]);
           for(i=1;i<n_codon_sites;i++)
               {
                   sum[i] = (fit_double)((fit_double)sum[i-1]*(i) + (1.-sortedprob[i]))/(i+1);
                   if((fit_double)sum[i] > fdr) break;
               }
           
           /* linear extrapolation */
           if((i == 1) || (i == n_codon_sites))
               thresh = 1.0;
           else
               {
                   thresh =
                       (1.-sortedprob[i-1]) +
                       ((1.-sortedprob[i]) - (1.-sortedprob[i-1])) *
                       (fdr - (fit_double)sum[i-1])/((fit_double)sum[i] - (fit_double)sum[i-1]);
                   thresh = 1. - thresh;
               }
           
/*            printf("THRESH = %f\n",thresh); */
           
           For(curr_site,n_codon_sites)
               {
                   if(postprob[curr_site] > thresh) 
                       fprintf(stdout,"%15E     1\n",
                               postprob[curr_site]);
                   else 
                       fprintf(stdout,"%15E     0\n",
                               postprob[curr_site]);
               }
                        
           n_data_sets--;

        }
    while(n_data_sets);

    Free(postprob);
    Free(sum);

    return 0;

}

#elif LUMP

int main(int argc, char **argv)
{
    int i,j,k;
    int u_acgt,v_acgt,w_acgt,x_acgt,y_acgt,z_acgt;
    int u_ry, v_ry, w_ry, x_ry, y_ry, z_ry;

    fit_double prr_xw, prr_xv, prr_xu, prr_wv, prr_wu, prr_vu;
    fit_double pyy_xw, pyy_xv, pyy_xu, pyy_wv, pyy_wu, pyy_vu;
    fit_double pry_xw, pry_xv, pry_xu, pry_wv, pry_wu, pry_vu;
    fit_double pyr_xw, pyr_xv, pyr_xu, pyr_wv, pyr_wu, pyr_vu;

    fit_double prr_uw, prr_wx, prr_vx, prr_ux, prr_vw, prr_uv;
    fit_double pyy_uw, pyy_wx, pyy_vx, pyy_ux, pyy_vw, pyy_uv;
    fit_double pry_uw, pry_wx, pry_vx, pry_ux, pry_vw, pry_uv;
    fit_double pyr_uw, pyr_wx, pyr_vx, pyr_ux, pyr_vw, pyr_uv;

    fit_double d_xw, d_xv, d_xu, d_wv, d_wu, d_vu;
    fit_double d_wx, d_vx, d_ux, d_vw, d_uw, d_uv;
    fit_double *p_patt_acgt,*p_patt_ry,*p_patt_ry_crunch;
    int patt_ry, patt_acgt;
    fit_double **Puy,**Pyv,**Pyz,**Pzw,**Pzx,**Pyw,**Pzv,**Pyx,**Puv, **Pwx, **Puw, **Pxv, **Pvu;
    fit_double luy,lyv,lyz,lzw,lzx;
    fit_double luy_rr,lyv_rr,lyz_rr,lzw_rr,lzx_rr,lyw_rr,lzv_rr,lyx_rr;
    fit_double luy_rr_loop,lyv_rr_loop,lyz_rr_loop,lzw_rr_loop,lzx_rr_loop,lyw_rr_loop,lzv_rr_loop,lyx_rr_loop;
    fit_double *pi_acgt,*pi_ry;
    qmat *qmat_struct;
    fit_double mr;
    fit_double a,b,c,d,e,f;
    fit_double A,B,C,D;
    fit_double R,Y;
    fit_double sum;
    fit_double **hadamard,**inv_hadamard;
    fit_double *br_len_spectrum,*buff;
    fit_double ry_support_uv,ry_support_uw,ry_support_ux;
    fit_double acgt_support_uv,acgt_support_uw,acgt_support_ux;
    fit_double lk,site_lk,lk_max;
    fit_double puw_bound, pwx_bound, pxv_bound, pvu_bound;

    p_patt_acgt      = (fit_double *)mCalloc((int)pow(4,4),sizeof(fit_double));
    p_patt_ry        = (fit_double *)mCalloc((int)pow(2,4),sizeof(fit_double));
    p_patt_ry_crunch = (fit_double *)mCalloc((int)pow(2,3),sizeof(fit_double));

    Puy = (fit_double **)mCalloc(4,sizeof(fit_double *));
    For(i,4) Puy[i] = (fit_double *)mCalloc(4,sizeof(fit_double ));

    Pyv = (fit_double **)mCalloc(4,sizeof(fit_double *));
    For(i,4) Pyv[i] = (fit_double *)mCalloc(4,sizeof(fit_double ));

    Pyz = (fit_double **)mCalloc(4,sizeof(fit_double *));
    For(i,4) Pyz[i] = (fit_double *)mCalloc(4,sizeof(fit_double ));

    Pzw = (fit_double **)mCalloc(4,sizeof(fit_double *));
    For(i,4) Pzw[i] = (fit_double *)mCalloc(4,sizeof(fit_double ));

    Pzx = (fit_double **)mCalloc(4,sizeof(fit_double *));
    For(i,4) Pzx[i] = (fit_double *)mCalloc(4,sizeof(fit_double ));

    Pyw = (fit_double **)mCalloc(4,sizeof(fit_double *));
    For(i,4) Pyw[i] = (fit_double *)mCalloc(4,sizeof(fit_double ));

    Pzv = (fit_double **)mCalloc(4,sizeof(fit_double *));
    For(i,4) Pzv[i] = (fit_double *)mCalloc(4,sizeof(fit_double ));

    Pyx = (fit_double **)mCalloc(4,sizeof(fit_double *));
    For(i,4) Pyx[i] = (fit_double *)mCalloc(4,sizeof(fit_double ));

    Puv = (fit_double **)mCalloc(4,sizeof(fit_double *));
    For(i,4) Puv[i] = (fit_double *)mCalloc(4,sizeof(fit_double ));

    Pwx = (fit_double **)mCalloc(4,sizeof(fit_double *));
    For(i,4) Pwx[i] = (fit_double *)mCalloc(4,sizeof(fit_double ));

    Puw = (fit_double **)mCalloc(4,sizeof(fit_double *));
    For(i,4) Puw[i] = (fit_double *)mCalloc(4,sizeof(fit_double ));

    Pxv = (fit_double **)mCalloc(4,sizeof(fit_double *));
    For(i,4) Pxv[i] = (fit_double *)mCalloc(4,sizeof(fit_double ));

    Pvu = (fit_double **)mCalloc(4,sizeof(fit_double *));
    For(i,4) Pvu[i] = (fit_double *)mCalloc(4,sizeof(fit_double ));

    
    pi_acgt = (fit_double *)mCalloc(4,sizeof(fit_double));
    pi_ry   = (fit_double *)mCalloc(2,sizeof(fit_double));

    hadamard = (fit_double **)mCalloc((int)pow(2,3),sizeof(fit_double *));
    For(i,(int)pow(2,3))
        hadamard[i] = (fit_double *)mCalloc((int)pow(2,3),sizeof(fit_double));

    inv_hadamard = (fit_double **)mCalloc((int)pow(2,3),sizeof(fit_double *));
    For(i,(int)pow(2,3))
        inv_hadamard[i] = (fit_double *)mCalloc((int)pow(2,3),sizeof(fit_double));

    br_len_spectrum = (fit_double *)mCalloc((int)pow(2,3),sizeof(fit_double));
    buff            = (fit_double *)mCalloc((int)pow(2,3),sizeof(fit_double));

    qmat_struct = Make_Qmat_Struct(4,1,1);
    
    pi_acgt[0] = 1.;     /* A */
    pi_acgt[1] = 1.;     /* C */
    pi_acgt[2] = 1.;     /* G */
    pi_acgt[3] = 1.;     /* T */
    
    sum = 0.0;
    For(i,4) sum += pi_acgt[i];
    For(i,4) pi_acgt[i] /= sum;

/*     printf("pi_A = %f\n",pi_acgt[0]); */
/*     printf("pi_C = %f\n",pi_acgt[1]); */
/*     printf("pi_G = %f\n",pi_acgt[2]); */
/*     printf("pi_T = %f\n",pi_acgt[3]); */
/*     printf("\n"); */

    R = pi_acgt[0] + pi_acgt[2];
    Y = pi_acgt[1] + pi_acgt[3];

    pi_ry[0] = R;
    pi_ry[1] = Y;


    /*

    u      w
     \y__z/
     /    \
    v      x

    */

    /* BEG Set branch lengths */
    luy = 2.0;
    lzw = 2.0;
    lyz = 0.0001;
    lyv = 0.5;
    lzx = 0.0001;
    /* END Set branch lengths */
    


    /* BEG Set substitution rate matrix */
    a = 40.0;   /* A <-> C */
    b = 1.0;   /* A <-> G */
    c = 1.0;   /* A <-> T */
    d = 1.0;   /* C <-> G */
    e = 1.0;   /* C <-> T */
    f = 1.0;   /* G <-> T */

    qmat_struct->qmat[0*4+1] = a*pi_acgt[1]; qmat_struct->qmat[0*4+2] = b*pi_acgt[2]; 
    qmat_struct->qmat[0*4+3] = c*pi_acgt[3];
    qmat_struct->qmat[1*4+0] = a*pi_acgt[0]; qmat_struct->qmat[1*4+2] = d*pi_acgt[2]; 
    qmat_struct->qmat[1*4+3] = e*pi_acgt[3];
    qmat_struct->qmat[2*4+0] = b*pi_acgt[0]; qmat_struct->qmat[2*4+1] = d*pi_acgt[1]; 
    qmat_struct->qmat[2*4+3] = f*pi_acgt[3];
    qmat_struct->qmat[3*4+0] = c*pi_acgt[0]; qmat_struct->qmat[3*4+1] = e*pi_acgt[1]; 
    qmat_struct->qmat[3*4+2] = f*pi_acgt[2];
 
    qmat_struct->qmat[0*4+0] = -(a*pi_acgt[1]+b*pi_acgt[2]+c*pi_acgt[3]); 
    qmat_struct->qmat[1*4+1] = -(a*pi_acgt[0]+d*pi_acgt[2]+e*pi_acgt[3]);
    qmat_struct->qmat[2*4+2] = -(b*pi_acgt[0]+d*pi_acgt[1]+f*pi_acgt[3]); 
    qmat_struct->qmat[3*4+3] = -(c*pi_acgt[0]+e*pi_acgt[1]+f*pi_acgt[2]);
    
    mr = .0;
    For(i,4) mr -= pi_acgt[i] * (qmat_struct->qmat[i*4+i]);
    For(i,4*4) qmat_struct->qmat[i] /= mr;
    /* END Set substitution rate matrix */
    
    

    /* BEG Compute transition probabilities */
    Update_Eigen(4,qmat_struct->qmat,qmat_struct->u_mat,qmat_struct->v_mat,qmat_struct->root_vct,NULL);
    
    For(i,4) qmat_struct->expD_mr_vct[i] = (fit_double)exp(qmat_struct->root_vct[i]);
		
    PMat_Numeric(luy,4,Puy,qmat_struct); PMat_Numeric(lyv,4,Pyv,qmat_struct);
    PMat_Numeric(lyz,4,Pyz,qmat_struct); PMat_Numeric(lzw,4,Pzw,qmat_struct);
    PMat_Numeric(lzx,4,Pzx,qmat_struct);
    /* END Compute transition probabilities */
		



    /* BEG Used for maximum parsimony */
    acgt_support_uv = acgt_support_uw = acgt_support_ux = 0.0;
    ry_support_uv = ry_support_uw = ry_support_ux = 0.0;
    /* END Used for maximum parsimony */
    
    For(patt_ry,pow(2,4)) p_patt_ry[patt_ry] = .0;
		
    lk = .0;
    For(patt_acgt,pow(4,4))
      {
	/* BEG ACGT state at every leaf */
	x_acgt = patt_acgt%4;
	w_acgt = ((int)(patt_acgt/4.))%4;
	v_acgt = ((int)(patt_acgt/(4.*4.)))%4;
	u_acgt = ((int)(patt_acgt/(4.*4.*4.)))%4;
	/* END ACGT state at every leaf */
	
	

	/* BEG Probabilities of the different ACGT site patterns */
	p_patt_acgt[patt_acgt] = 0.;
	For(y_acgt,4)
	  {
	    For(z_acgt,4)
	      {
		p_patt_acgt[patt_acgt] += 
		  pi_acgt[u_acgt] *
		  Puy[u_acgt][y_acgt] * Pyv[y_acgt][v_acgt] *
		  Pyz[y_acgt][z_acgt] * Pzw[z_acgt][w_acgt] *
		  Pzx[z_acgt][x_acgt];                                
	      }
	  }

	lk += (fit_double)log(p_patt_acgt[patt_acgt]) * p_patt_acgt[patt_acgt];


	/* END Probabilities of the different ACGT site patterns */
	
	/*             printf("%d%d%d%d %f\n", */
	/*                    x_acgt, */
	/*                    w_acgt, */
	/*                    v_acgt, */
	/*                    u_acgt, */
	/*                    p_patt_acgt[patt_acgt]); */
	
	

	/* BEG Corresponding RY states */ 
	x_ry = ACGT_to_RY(x_acgt);
	w_ry = ACGT_to_RY(w_acgt);
	v_ry = ACGT_to_RY(v_acgt);
	u_ry = ACGT_to_RY(u_acgt);
	/* END Corresponding RY states */ 
	
	

	/* BEG Corresponding RY site pattern number */
	patt_ry =  pow(2,4)/pow(2,1)*x_ry + pow(2,4)/pow(2,2)*w_ry + 
	  pow(2,4)/pow(2,3)*v_ry + pow(2,4)/pow(2,4)*u_ry ;
	/* END Corresponding RY site pattern number */
	
	

	/* BEG Corresponding RY site pattern frequency */
	p_patt_ry[patt_ry] += p_patt_acgt[patt_acgt];
	/* END Corresponding RY site pattern frequency */
		    
	

	/* BEG Used for maximum parsimony */
	if((u_acgt == v_acgt) && (w_acgt == x_acgt) && (u_acgt != w_acgt))
	  acgt_support_uv += p_patt_acgt[patt_acgt];
	if((u_acgt == w_acgt) && (v_acgt == x_acgt) && (u_acgt != v_acgt))
	  acgt_support_uw += p_patt_acgt[patt_acgt];
	if((u_acgt == x_acgt) && (w_acgt == v_acgt) && (u_acgt != w_acgt))
	  acgt_support_ux += p_patt_acgt[patt_acgt];
	
	if((u_ry == v_ry) && (w_ry == x_ry) && (u_ry != w_ry))
	  ry_support_uv += p_patt_acgt[patt_acgt];
	if((u_ry == w_ry) && (v_ry == x_ry) && (u_ry != v_ry))
	  ry_support_uw += p_patt_acgt[patt_acgt];
	if((u_ry == x_ry) && (w_ry == v_ry) && (u_ry != w_ry))
	  ry_support_ux += p_patt_acgt[patt_acgt];
	/* END Used for maximum parsimony */
      }
    
    /*     printf("{uv}/{wx} -> %f %f\n",acgt_support_uv,ry_support_uv); */
    /*     printf("{uw}/{vx} -> %f %f\n",acgt_support_uw,ry_support_uw); */
    /*     printf("{ux}/{wv} -> %f %f\n",acgt_support_ux,ry_support_ux); */

    printf("LnL acgt = %f\n",lk);

    /* BEG Site pattern frequencies under a RY model where freq(R)=freq(Y) */ 
    For(i,pow(2,3)) p_patt_ry_crunch[i] = p_patt_ry[i] + p_patt_ry[(int)pow(2,4)-i-1];
    /* END Site pattern frequencies under a RY model where freq(R)=freq(Y) */ 

    

    /* BEG Check that the previous calculation went fine */ 
    sum = .0;
    For(patt_ry,pow(2,4)) sum+=p_patt_ry[patt_ry];
    if((sum < 1.-1.E-5) || (sum > 1.+1.E-5)) Exit("\n. Error in sum 1\n");
    sum = .0;
    For(patt_ry,pow(2,3)) sum+=p_patt_ry_crunch[patt_ry];
    if((sum < 1.-1.E-5) || (sum > 1.+1.E-5)) Exit("\n. Error in sum 2\n");
    /* END Check that the previous calculation went fine */ 
    
    
    
    

    /* BEG Compute the expected pairwise RY site pattern frequencies */ 
    prr_xw = prr_xv = prr_xu = prr_wv = prr_wu = prr_vu = .0;
    pyy_xw = pyy_xv = pyy_xu = pyy_wv = pyy_wu = pyy_vu = .0;
    pry_xw = pry_xv = pry_xu = pry_wv = pry_wu = pry_vu = .0;
    pyr_xw = pyr_xv = pyr_xu = pyr_wv = pyr_wu = pyr_vu = .0;
    For(patt_ry,pow(2,4))
      {
	u_ry = patt_ry%2;
	v_ry = ((int)(patt_ry/2.))%2;
	w_ry = ((int)(patt_ry/(2.*2.)))%2;
	x_ry = ((int)(patt_ry/(2.*2.*2.)))%2;
	
	if(     (x_ry == 0) && (w_ry == 0)) prr_xw += p_patt_ry[patt_ry];
	else if((x_ry == 1) && (w_ry == 1)) pyy_xw += p_patt_ry[patt_ry];
	if(     (x_ry == 0) && (w_ry == 1)) pry_xw += p_patt_ry[patt_ry];
	else if((x_ry == 1) && (w_ry == 0)) pyr_xw += p_patt_ry[patt_ry];
	
	if(     (x_ry == 0) && (v_ry == 0)) prr_xv += p_patt_ry[patt_ry];
	else if((x_ry == 1) && (v_ry == 1)) pyy_xv += p_patt_ry[patt_ry];
	if(     (x_ry == 0) && (v_ry == 1)) pry_xv += p_patt_ry[patt_ry];
	else if((x_ry == 1) && (v_ry == 0)) pyr_xv += p_patt_ry[patt_ry];
	
	if(     (x_ry == 0) && (u_ry == 0)) prr_xu += p_patt_ry[patt_ry];
	else if((x_ry == 1) && (u_ry == 1)) pyy_xu += p_patt_ry[patt_ry];
	if(     (x_ry == 0) && (u_ry == 1)) pry_xu += p_patt_ry[patt_ry];
	else if((x_ry == 1) && (u_ry == 0)) pyr_xu += p_patt_ry[patt_ry];
	
	if(     (w_ry == 0) && (v_ry == 0)) prr_wv += p_patt_ry[patt_ry];
	else if((w_ry == 1) && (v_ry == 1)) pyy_wv += p_patt_ry[patt_ry];
	if(     (w_ry == 0) && (v_ry == 1)) pry_wv += p_patt_ry[patt_ry];
	else if((w_ry == 1) && (v_ry == 0)) pyr_wv += p_patt_ry[patt_ry];
	
	if(     (w_ry == 0) && (u_ry == 0)) prr_wu += p_patt_ry[patt_ry];
	else if((w_ry == 1) && (u_ry == 1)) pyy_wu += p_patt_ry[patt_ry];
	if(     (w_ry == 0) && (u_ry == 1)) pry_wu += p_patt_ry[patt_ry];
	else if((w_ry == 1) && (u_ry == 0)) pyr_wu += p_patt_ry[patt_ry];
		    
	if(     (v_ry == 0) && (u_ry == 0)) prr_vu += p_patt_ry[patt_ry];
	else if((v_ry == 1) && (u_ry == 1)) pyy_vu += p_patt_ry[patt_ry];
	if(     (v_ry == 0) && (u_ry == 1)) pry_vu += p_patt_ry[patt_ry];
	else if((v_ry == 1) && (u_ry == 0)) pyr_vu += p_patt_ry[patt_ry];
      }
    /* END Compute the expected pairwise RY site pattern frequencies */ 
    
		
    prr_uw = prr_wu; pry_uw = pry_wu; pyr_uw = pyr_wu; pyy_uw = pyy_wu;
    prr_wx = prr_xw; pry_wx = pry_xw; pyr_wx = pyr_xw; pyy_wx = pyy_xw;
    prr_vx = prr_xv; pry_vx = pry_xv; pyr_vx = pyr_xv; pyy_vx = pyy_xv;
    prr_ux = prr_xu; pry_ux = pry_xu; pyr_ux = pyr_xu; pyy_ux = pyy_xu;
    prr_vw = prr_wv; pry_vw = pry_wv; pyr_vw = pyr_wv; pyy_vw = pyy_wv;
    prr_uv = prr_vu; pry_uv = pry_vu; pyr_uv = pyr_vu; pyy_uv = pyy_vu;



    

    /* BEG Compute the maximum likelihood RY pairwise distances */
    A = R * Y;
    
    B = R*R + Y*Y - prr_xw*R - pyy_xw*Y; C = -(prr_xw*Y + pyy_xw*R - R*Y); D = B*B-4.*A*C;
    d_xw = -(2.*R*Y)*(fit_double)log((-B+sqrt(D))/(2.*A));
    
    B = R*R + Y*Y - prr_xv*R - pyy_xv*Y; C = -(prr_xv*Y + pyy_xv*R - R*Y); D = B*B-4.*A*C;
    d_xv = -(2.*R*Y)*(fit_double)log((-B+sqrt(D))/(2.*A));
    
    B = R*R + Y*Y - prr_xu*R - pyy_xu*Y; C = -(prr_xu*Y + pyy_xu*R - R*Y); D = B*B-4.*A*C;
    d_xu = -(2.*R*Y)*(fit_double)log((-B+sqrt(D))/(2.*A));
    
    B = R*R + Y*Y - prr_wv*R - pyy_wv*Y; C = -(prr_wv*Y + pyy_wv*R - R*Y); D = B*B-4.*A*C;
    d_wv = -(2.*R*Y)*(fit_double)log((-B+sqrt(D))/(2.*A));
    
    B = R*R + Y*Y - prr_wu*R - pyy_wu*Y; C = -(prr_wu*Y + pyy_wu*R - R*Y); D = B*B-4.*A*C;
    d_wu = -(2.*R*Y)*(fit_double)log((-B+sqrt(D))/(2.*A));
    
    B = R*R + Y*Y - prr_vu*R - pyy_vu*Y; C = -(prr_vu*Y + pyy_vu*R - R*Y); D = B*B-4.*A*C;
    d_vu = -(2.*R*Y)*(fit_double)log((-B+sqrt(D))/(2.*A));

    d_wx = d_xw; d_vx = d_xv; d_ux = d_xu; d_vw = d_wv; d_uw = d_wu; d_uv = d_vu;

    /*     printf("d_vu = %f\n",d_vu); */
    
    /*     printf("\n"); */
    /*     printf("d_vu + d_xw = %f\n",d_vu+d_xw); */
    /*     printf("d_wu + d_xv = %f\n",d_wu+d_xv); */
    /*     printf("d_xu + d_wv = %f\n",d_xu+d_wv); */
    
    /* END Compute the maximum likelihood RY pairwise distances */

    /*    
	  TI
	 u      w
	  \y__z/
	  /    \
	 v      x

	 TII
	 u      v
	  \y__z/
	  /    \
	 w      x

	 TIII
	 u      w
	  \y__z/
	  /    \
	 x      v
    */

    /* BEG Least square branch lengths estimates for topology TI */
    luy_rr = (1./4.)*(d_uw + d_ux - d_vw - d_vx + 2.*d_uv); if(luy_rr < .0) luy_rr = .0;
    lyv_rr = (1./4.)*(d_vw + d_vx - d_uw - d_ux + 2.*d_uv); if(lyv_rr < .0) lyv_rr = .0;
    lzw_rr = (1./4.)*(d_uw + d_vw - d_ux - d_vx + 2.*d_wx); if(lzw_rr < .0) lzw_rr = .0;
    lzx_rr = (1./4.)*(d_ux + d_vx - d_uw - d_vw + 2.*d_wx); if(lzx_rr < .0) lzx_rr = .0;
    lyz_rr = (1./4.)*(d_uw + d_vx + d_ux + d_vw - 2.*d_uv - 2.*d_wx); if(lyz_rr < .0) lyz_rr = 1.E-6;
    /* END Least square branch lengths estimates for topology TI */
 

    printf("luy = %f\n",luy_rr);
    printf("lyv = %f\n",lyv_rr);
    printf("lzw = %f\n",lzw_rr);
    printf("lzx = %f\n",lzx_rr);
    printf("lyz = %f\n",lyz_rr);

    lk_max = -1.E10;
    for(luy_rr_loop = luy_rr/2.;luy_rr_loop < luy_rr*5.; luy_rr_loop += (luy_rr*2. - luy_rr/2.)/5.)
      {
	for(lyv_rr_loop = lyv_rr/2.;lyv_rr_loop < lyv_rr*5.; lyv_rr_loop += (lyv_rr*2. - lyv_rr/2.)/5.)
	  {
	    for(lzw_rr_loop = lzw_rr/2.;lzw_rr_loop < lzw_rr*5.; lzw_rr_loop += (lzw_rr*2. - lzw_rr/2.)/5.)
	      {
		for(lzx_rr_loop = lzx_rr/2.;lzx_rr_loop < lzx_rr*5.; lzx_rr_loop += (lzx_rr*2. - lzx_rr/2.)/5.)
		  {
		    for(lyz_rr_loop = lyz_rr/2.;lyz_rr_loop < lyz_rr*5.; lyz_rr_loop += (lyz_rr*2. - lyz_rr/2.)/5.)
		      {
			Puy[0][0]=R+Y*(fit_double)exp(-luy_rr_loop/(2.*R*Y)); Puy[0][1]=Y-Y*(fit_double)exp(-luy_rr_loop/(2.*R*Y)); 
			Puy[1][0]=R-R*(fit_double)exp(-luy_rr_loop/(2.*R*Y)); Puy[1][1]=Y+R*(fit_double)exp(-luy_rr_loop/(2.*R*Y));
			Pyv[0][0]=R+Y*(fit_double)exp(-lyv_rr_loop/(2.*R*Y)); Pyv[0][1]=Y-Y*(fit_double)exp(-lyv_rr_loop/(2.*R*Y)); 
			Pyv[1][0]=R-R*(fit_double)exp(-lyv_rr_loop/(2.*R*Y)); Pyv[1][1]=Y+R*(fit_double)exp(-lyv_rr_loop/(2.*R*Y));
			Pzw[0][0]=R+Y*(fit_double)exp(-lzw_rr_loop/(2.*R*Y)); Pzw[0][1]=Y-Y*(fit_double)exp(-lzw_rr_loop/(2.*R*Y)); 
			Pzw[1][0]=R-R*(fit_double)exp(-lzw_rr_loop/(2.*R*Y)); Pzw[1][1]=Y+R*(fit_double)exp(-lzw_rr_loop/(2.*R*Y));
			Pzx[0][0]=R+Y*(fit_double)exp(-lzx_rr_loop/(2.*R*Y)); Pzx[0][1]=Y-Y*(fit_double)exp(-lzx_rr_loop/(2.*R*Y)); 
			Pzx[1][0]=R-R*(fit_double)exp(-lzx_rr_loop/(2.*R*Y)); Pzx[1][1]=Y+R*(fit_double)exp(-lzx_rr_loop/(2.*R*Y));
			Pyz[0][0]=R+Y*(fit_double)exp(-lyz_rr_loop/(2.*R*Y)); Pyz[0][1]=Y-Y*(fit_double)exp(-lyz_rr_loop/(2.*R*Y)); 
			Pyz[1][0]=R-R*(fit_double)exp(-lyz_rr_loop/(2.*R*Y)); Pyz[1][1]=Y+R*(fit_double)exp(-lyz_rr_loop/(2.*R*Y));

			lk = .0;
			For(patt_ry,pow(2,4))
			  {
			    /* BEG RY state at every leaf */
			    u_ry = patt_ry%2;
			    v_ry = ((int)(patt_ry/2.))%2;
			    w_ry = ((int)(patt_ry/(2.*2.)))%2;
			    x_ry = ((int)(patt_ry/(2.*2.*2.)))%2;
			    /* END RY state at every leaf */
			    
			    /* BEG Likelihood calculation */
			    site_lk = .0;
			    For(y_ry,2) For(z_ry,2)
			      site_lk += pi_ry[u_ry] *
			      Puy[u_ry][y_ry] * Pyv[y_ry][v_ry] *
			      Pyz[y_ry][z_ry] * Pzw[z_ry][w_ry] *
			      Pzx[z_ry][x_ry];                                
			    /* END Likelihood calculation */
			    lk += p_patt_ry[patt_ry] * (fit_double)log(site_lk);
			  }
			if(lk > lk_max) lk_max = lk;
		      }
		  }
	      }
	  }
      }

    printf("TI, lk_max = %f\n",lk_max);
    

    /* END Least square branch lengths estimates for topology TII */
    luy_rr = (1./4.)*(d_uv + d_ux - d_vw - d_wx + 2.*d_uw); if(luy_rr < .0) luy_rr = .0;
    lyw_rr = (1./4.)*(d_vw + d_wx - d_uv - d_ux + 2.*d_uw); if(lyw_rr < .0) lyw_rr = .0;
    lzv_rr = (1./4.)*(d_vu + d_vw - d_ux - d_wx + 2.*d_vx); if(lzv_rr < .0) lzv_rr = .0;
    lzx_rr = (1./4.)*(d_xu + d_xw - d_vw - d_uv + 2.*d_vx); if(lzx_rr < .0) lzx_rr = .0;
    lyz_rr = (1./4.)*(d_uv + d_wx + d_ux + d_wv - 2.*d_wu - 2.*d_vx); if(lyz_rr < .0) lyz_rr = 1.E-6;
    /* END Least square branch lengths estimates for topology TII */

    printf("luy = %f\n",luy_rr);
    printf("lyw = %f\n",lyw_rr);
    printf("lzv = %f\n",lzv_rr);
    printf("lzx = %f\n",lzx_rr);
    printf("lyz = %f\n",lyz_rr);

    lk_max = -1.E10;
    for(luy_rr_loop = luy_rr/2.;luy_rr_loop < luy_rr*5.; luy_rr_loop += (luy_rr*2. - luy_rr/2.)/5.)
      {
	for(lyw_rr_loop = lyw_rr/2.;lyw_rr_loop < lyw_rr*5.; lyw_rr_loop += (lyw_rr*2. - lyw_rr/2.)/5.)
	  {
	    for(lzv_rr_loop = lzv_rr/2.;lzv_rr_loop < lzv_rr*5.; lzv_rr_loop += (lzv_rr*2. - lzv_rr/2.)/5.)
	      {
		for(lzx_rr_loop = lzx_rr/2.;lzx_rr_loop < lzx_rr*5.; lzx_rr_loop += (lzx_rr*2. - lzx_rr/2.)/5.)
		  {
		    for(lyz_rr_loop = lyz_rr/2.;lyz_rr_loop < lyz_rr*5.; lyz_rr_loop += (lyz_rr*2. - lyz_rr/2.)/5.)
		      {
			Puy[0][0]=R+Y*(fit_double)exp(-luy_rr_loop/(2.*R*Y)); Puy[0][1]=Y-Y*(fit_double)exp(-luy_rr_loop/(2.*R*Y)); 
			Puy[1][0]=R-R*(fit_double)exp(-luy_rr_loop/(2.*R*Y)); Puy[1][1]=Y+R*(fit_double)exp(-luy_rr_loop/(2.*R*Y));
			Pyw[0][0]=R+Y*(fit_double)exp(-lyw_rr_loop/(2.*R*Y)); Pyw[0][1]=Y-Y*(fit_double)exp(-lyw_rr_loop/(2.*R*Y)); 
			Pyw[1][0]=R-R*(fit_double)exp(-lyw_rr_loop/(2.*R*Y)); Pyw[1][1]=Y+R*(fit_double)exp(-lyw_rr_loop/(2.*R*Y));
			Pzv[0][0]=R+Y*(fit_double)exp(-lzv_rr_loop/(2.*R*Y)); Pzv[0][1]=Y-Y*(fit_double)exp(-lzv_rr_loop/(2.*R*Y)); 
			Pzv[1][0]=R-R*(fit_double)exp(-lzv_rr_loop/(2.*R*Y)); Pzv[1][1]=Y+R*(fit_double)exp(-lzv_rr_loop/(2.*R*Y));
			Pzx[0][0]=R+Y*(fit_double)exp(-lzx_rr_loop/(2.*R*Y)); Pzx[0][1]=Y-Y*(fit_double)exp(-lzx_rr_loop/(2.*R*Y)); 
			Pzx[1][0]=R-R*(fit_double)exp(-lzx_rr_loop/(2.*R*Y)); Pzx[1][1]=Y+R*(fit_double)exp(-lzx_rr_loop/(2.*R*Y));
			Pyz[0][0]=R+Y*(fit_double)exp(-lyz_rr_loop/(2.*R*Y)); Pyz[0][1]=Y-Y*(fit_double)exp(-lyz_rr_loop/(2.*R*Y)); 
			Pyz[1][0]=R-R*(fit_double)exp(-lyz_rr_loop/(2.*R*Y)); Pyz[1][1]=Y+R*(fit_double)exp(-lyz_rr_loop/(2.*R*Y));
			
			lk = .0;
			For(patt_ry,pow(2,4))
			  {
			    /* BEG RY state at every leaf */
			    u_ry = patt_ry%2;
			    v_ry = ((int)(patt_ry/2.))%2;
			    w_ry = ((int)(patt_ry/(2.*2.)))%2;
			    x_ry = ((int)(patt_ry/(2.*2.*2.)))%2;
			    /* END RY state at every leaf */
			    
			    /* BEG Likelihood calculation */
			    site_lk = .0;
			    For(y_ry,2) For(z_ry,2)
			      site_lk += pi_ry[u_ry] *
			      Puy[u_ry][y_ry] * Pyw[y_ry][w_ry] *
			      Pyz[y_ry][z_ry] * Pzv[z_ry][v_ry] *
			      Pzx[z_ry][x_ry];                                
			    /* END Likelihood calculation */
			    lk += p_patt_ry[patt_ry] * (fit_double)log(site_lk);
			  }
			if(lk > lk_max) lk_max = lk;
		      }
		  }
	      }
	  }
      }

    printf("TII, lk_max = %f\n",lk_max);

    /* END Least square branch lengths estimates for topology TIII */
    luy_rr = (1./4.)*(d_uw + d_uv - d_xw - d_xv + 2.*d_ux); if(luy_rr < .0) luy_rr = .0;
    lyx_rr = (1./4.)*(d_xw + d_xv - d_uw - d_uv + 2.*d_ux); if(lyx_rr < .0) lyx_rr = .0;
    lzw_rr = (1./4.)*(d_xw + d_uw - d_uv - d_xv + 2.*d_vw); if(lzw_rr < .0) lzw_rr = .0;
    lzv_rr = (1./4.)*(d_uv + d_xv - d_xw - d_uw + 2.*d_vw); if(lzv_rr < .0) lzv_rr = .0;
    lyz_rr = (1./4.)*(d_uw + d_xv + d_uv + d_xw - 2.*d_ux - 2.*d_vw); if(lyz_rr < .0) lyz_rr = 1.E-6;
    /* END Least square branch lengths estimates for topology TIII */
		
    printf("luy = %f\n",luy_rr);
    printf("lyx = %f\n",lyx_rr);
    printf("lzw = %f\n",lzw_rr);
    printf("lzv = %f\n",lzv_rr);
    printf("lyz = %f\n",lyz_rr);


    lk_max = -1.E10;
    for(luy_rr_loop = luy_rr/2.;luy_rr_loop < luy_rr*5.; luy_rr_loop += (luy_rr*2. - luy_rr/2.)/5.)
      {
	for(lyx_rr_loop = lyx_rr/2.;lyx_rr_loop < lyx_rr*5.; lyx_rr_loop += (lyx_rr*2. - lyx_rr/2.)/5.)
	  {
	    for(lzw_rr_loop = lzw_rr/2.;lzw_rr_loop < lzw_rr*5.; lzw_rr_loop += (lzw_rr*2. - lzw_rr/2.)/5.)
	      {
		for(lzv_rr_loop = lzv_rr/2.;lzv_rr_loop < lzv_rr*5.; lzv_rr_loop += (lzv_rr*2. - lzx_rr/2.)/5.)
		  {
		    for(lyz_rr_loop = lyz_rr/2.;lyz_rr_loop < lyz_rr*5.; lyz_rr_loop += (lyz_rr*2. - lyz_rr/2.)/5.)
		      {
			
			Puy[0][0]=R+Y*(fit_double)exp(-luy_rr_loop/(2.*R*Y)); Puy[0][1]=Y-Y*(fit_double)exp(-luy_rr_loop/(2.*R*Y)); 
			Puy[1][0]=R-R*(fit_double)exp(-luy_rr_loop/(2.*R*Y)); Puy[1][1]=Y+R*(fit_double)exp(-luy_rr_loop/(2.*R*Y));
			Pyx[0][0]=R+Y*(fit_double)exp(-lyx_rr_loop/(2.*R*Y)); Pyx[0][1]=Y-Y*(fit_double)exp(-lyx_rr_loop/(2.*R*Y)); 
			Pyx[1][0]=R-R*(fit_double)exp(-lyx_rr_loop/(2.*R*Y)); Pyx[1][1]=Y+R*(fit_double)exp(-lyx_rr_loop/(2.*R*Y));
			Pzw[0][0]=R+Y*(fit_double)exp(-lzw_rr_loop/(2.*R*Y)); Pzw[0][1]=Y-Y*(fit_double)exp(-lzw_rr_loop/(2.*R*Y)); 
			Pzw[1][0]=R-R*(fit_double)exp(-lzw_rr_loop/(2.*R*Y)); Pzw[1][1]=Y+R*(fit_double)exp(-lzw_rr_loop/(2.*R*Y));
			Pzv[0][0]=R+Y*(fit_double)exp(-lzv_rr_loop/(2.*R*Y)); Pzv[0][1]=Y-Y*(fit_double)exp(-lzv_rr_loop/(2.*R*Y)); 
			Pzv[1][0]=R-R*(fit_double)exp(-lzv_rr_loop/(2.*R*Y)); Pzv[1][1]=Y+R*(fit_double)exp(-lzv_rr_loop/(2.*R*Y));
			Pzv[0][0]=R+Y*(fit_double)exp(-lzv_rr_loop/(2.*R*Y)); Pzv[0][1]=Y-Y*(fit_double)exp(-lzv_rr_loop/(2.*R*Y)); 
			Pzv[1][0]=R-R*(fit_double)exp(-lzv_rr_loop/(2.*R*Y)); Pzv[1][1]=Y+R*(fit_double)exp(-lzv_rr_loop/(2.*R*Y));
			
			lk = .0;
			For(patt_ry,pow(2,4))
			  {
			    /* BEG RY state at every leaf */
			    u_ry = patt_ry%2;
			    v_ry = ((int)(patt_ry/2.))%2;
			    w_ry = ((int)(patt_ry/(2.*2.)))%2;
			    x_ry = ((int)(patt_ry/(2.*2.*2.)))%2;
			    /* END RY state at every leaf */
			    
			    /* BEG Likelihood calculation */
			    site_lk = .0;
			    For(y_ry,2) For(z_ry,2)
			      site_lk += pi_ry[u_ry] *
			      Puy[u_ry][y_ry] * Pyx[y_ry][x_ry] *
			      Pyz[y_ry][z_ry] * Pzw[z_ry][w_ry] *
			      Pzv[z_ry][v_ry];                                
			    /* END Likelihood calculation */
			    
			    lk += p_patt_ry[patt_ry] * (fit_double)log(site_lk);
			  }
			if(lk > lk_max) lk_max = lk;
		      }
		  }
	      }
	  }
      }
    
    printf("TIII lk = %f\n",lk_max);
			

    /* BEG Upper bound */
    Puv[0][0]=R+Y*(fit_double)exp(-d_uv/(2.*R*Y)); Puv[0][1]=Y-Y*(fit_double)exp(-d_uv/(2.*R*Y)); 
    Puv[1][0]=R-R*(fit_double)exp(-d_uv/(2.*R*Y)); Puv[1][1]=Y+R*(fit_double)exp(-d_uv/(2.*R*Y));
    Pwx[0][0]=R+Y*(fit_double)exp(-d_wx/(2.*R*Y)); Pwx[0][1]=Y-Y*(fit_double)exp(-d_wx/(2.*R*Y)); 
    Pwx[1][0]=R-R*(fit_double)exp(-d_wx/(2.*R*Y)); Pwx[1][1]=Y+R*(fit_double)exp(-d_wx/(2.*R*Y));
    /* END Upper bound */
    
    
    printf("Upper bound (Hendy & Holland) = %f\n",
	   -
	   (prr_uv * (fit_double)log(R*Puv[0][0]) +
	    pry_uv * (fit_double)log(R*Puv[0][1]) +
	    pyr_uv * (fit_double)log(Y*Puv[1][0]) +
	    pyy_uv * (fit_double)log(Y*Puv[1][1])) *
	   (prr_wx * (fit_double)log(R*Pwx[0][0]) +
	    pry_wx * (fit_double)log(R*Pwx[0][1]) +
	    pyr_wx * (fit_double)log(Y*Pwx[1][0]) +
	    pyy_wx * (fit_double)log(Y*Pwx[1][1])));




    /* BEG Upper bound */
    Puw[0][0]=R+Y*(fit_double)exp(-d_uw/(2.*R*Y)); Puw[0][1]=Y-Y*(fit_double)exp(-d_uw/(2.*R*Y));
    Puw[1][0]=R-R*(fit_double)exp(-d_uw/(2.*R*Y)); Puw[1][1]=Y+R*(fit_double)exp(-d_uw/(2.*R*Y));
    puw_bound = pow(Puw[0][0],prr_uw)*pow(Puw[0][1],pry_uw)*pow(Puw[1][0],pyr_uw)*pow(Puw[1][1],pyy_uw);

    Pwx[0][0]=R+Y*(fit_double)exp(-d_wx/(2.*R*Y)); Pwx[0][1]=Y-Y*(fit_double)exp(-d_wx/(2.*R*Y));
    Pwx[1][0]=R-R*(fit_double)exp(-d_wx/(2.*R*Y)); Pwx[1][1]=Y+R*(fit_double)exp(-d_wx/(2.*R*Y));
    pwx_bound = pow(Pwx[0][0],prr_wx)*pow(Pwx[0][1],pry_wx)*pow(Pwx[1][0],pyr_wx)*pow(Pwx[1][1],pyy_wx);

    Pxv[0][0]=R+Y*(fit_double)exp(-d_xv/(2.*R*Y)); Pxv[0][1]=Y-Y*(fit_double)exp(-d_xv/(2.*R*Y));
    Pxv[1][0]=R-R*(fit_double)exp(-d_xv/(2.*R*Y)); Pxv[1][1]=Y+R*(fit_double)exp(-d_xv/(2.*R*Y));
    pxv_bound = pow(Pxv[0][0],prr_xv)*pow(Pxv[0][1],pry_xv)*pow(Pxv[1][0],pyr_xv)*pow(Pxv[1][1],pyy_xv);

    Pvu[0][0]=R+Y*(fit_double)exp(-d_vu/(2.*R*Y)); Pvu[0][1]=Y-Y*(fit_double)exp(-d_vu/(2.*R*Y));
    Pvu[1][0]=R-R*(fit_double)exp(-d_vu/(2.*R*Y)); Pvu[1][1]=Y+R*(fit_double)exp(-d_vu/(2.*R*Y));
    pvu_bound = pow(Pvu[0][0],prr_vu)*pow(Pvu[0][1],pry_vu)*pow(Pvu[1][0],pyr_vu)*pow(Pvu[1][1],pyy_vu);
/*     END Upper bound */
    
    
    printf("Upper bound (Han & Bryant) = %f\n",
	   (fit_double)log(sqrt(puw_bound*pwx_bound*pxv_bound*pvu_bound)));



    /* BEG Hadamard conjugation */
    For(i,pow(2,3)) For(j,pow(2,3)) hadamard[i][j] = 1.0;
    
    For(i,3)
      {
	for(j=0;j<pow(2,i+1);j++) for(k=pow(2,i);k<pow(2,i+1);k++)
	  hadamard[j][k] = hadamard[j][k-(int)pow(2,i)];
	
	for(j=pow(2,i);j<pow(2,i+1);j++) for(k=0;k<pow(2,i+1);k++)
	  hadamard[j][k] = hadamard[j-(int)pow(2,i)][k];
	
	for(j=pow(2,i);j<pow(2,i+1);j++) for(k=pow(2,i);k<pow(2,i+1);k++)
	  hadamard[j][k] = -hadamard[j-(int)pow(2,i)][k-(int)pow(2,i)];
      }
    
    For(i,pow(2,3)) For(j,pow(2,3)) inv_hadamard[i][j] = hadamard[i][j] / pow(2,3);
    
    /*     For(i,pow(2,3)) printf("patt -> %f\n",p_patt_ry_crunch[i]); */
    /*     p_patt_ry_crunch[0] = 0.408712; */
    /*     p_patt_ry_crunch[1] = 0.177096; */
    /*     p_patt_ry_crunch[2] = 0.041928; */
    /*     p_patt_ry_crunch[3] = 0.048264; */
    /*     p_patt_ry_crunch[4] = 0.177096; */
    /*     p_patt_ry_crunch[5] = 0.077832; */
    /*     p_patt_ry_crunch[6] = 0.027144; */
    /*     p_patt_ry_crunch[7] = 0.041928; */
    
    For(i,pow(2,3)) For(j,pow(2,3)) br_len_spectrum[i] += hadamard[i][j] * p_patt_ry_crunch[j];
    For(i,pow(2,3)) br_len_spectrum[i] = (fit_double)log(br_len_spectrum[i]);
    For(i,pow(2,3)) buff[i] = br_len_spectrum[i];
    
    For(i,pow(2,3))
      {
	a = .0;
	For(j,pow(2,3)) a += inv_hadamard[i][j] * buff[j];
	br_len_spectrum[i] = a;
      }
    
    /*     printf("====================\n"); */
    /*     printf("0,{u,v,w,x} -> %f\n",br_len_spectrum[0]); */
    /*     printf("{u},{v,w,x} -> %f\n",br_len_spectrum[1]); */
    /*     printf("{v},{u,w,x} -> %f\n",br_len_spectrum[2]); */
    /*     printf("{u,v},{w,x} -> %f\n",br_len_spectrum[3]); */
    /*     printf("{w},{u,v,x} -> %f\n",br_len_spectrum[4]); */
    /*     printf("{u,w},{v,x} -> %f\n",br_len_spectrum[5]); */
    /*     printf("{v,w},{u,x} -> %f\n",br_len_spectrum[6]); */
    /*     printf("{u,v,w},{x} -> %f\n",br_len_spectrum[7]); */
    /*     printf("====================\n"); */
    
    /* END Hadamard conjugation */
    
    
    Free(p_patt_acgt);
    Free(p_patt_ry);
    For(i,4) Free(Puy[i]); Free(Puy);
    For(i,4) Free(Pyv[i]); Free(Pyv);
    For(i,4) Free(Pyz[i]); Free(Pyz);
    For(i,4) Free(Pzw[i]); Free(Pzw);
    For(i,4) Free(Pzx[i]); Free(Pzx);
    return -1;
}


#elif COLTREE

int main(int argc, char **argv) 
{
  FILE *fp_tree,*fp_tree_ps;
  tdraw *w;
  int tree_num;
  arbre *tree;
  edge *b_root;
  char c;
  int i;
  int b_root_num;

  tree   = NULL;
  b_root = NULL;
  b_root_num = 0;

  fp_tree    = Openfile(argv[1],0);
  fp_tree_ps = Openfile(argv[2],1);
  
  tree_num = 0;
  do
    {
      c=fgetc(fp_tree);
      if(c==';') tree_num++;
      else if(c==EOF) break;
    }while(1);
  fclose(fp_tree);

  fp_tree = Openfile(argv[1],0);
  Print_Postscript_Header(tree_num,fp_tree_ps);
  tree_num = 0;
  do
    {
      tree = Read_Tree_File(fp_tree);
      
      if(!tree) break;
      else
	{	  
	  if(!tree_num)
	    {
	      For(i,2*tree->n_otu-2)
		{
		  if(tree->noeud[i]->check_branch > -1)
		    {
		      b_root_num = tree->noeud[i]->b[tree->noeud[i]->check_branch]->num;
		      break;
		    }
		}	      
	    }
	  
	  b_root = tree->t_edges[b_root_num];
	  tree->fp_tree = fp_tree;
	  w = Make_Tdraw_Struct(tree);
/* 	  printf("%s\n",Write_Tree(tree)); */
	  Draw_Tree(b_root,w,tree);
	  Print_Tree_Postscript(b_root,fp_tree_ps,tree_num,w,tree);
	  tree_num++;
	}
      
    }while(1);

  Print_Postscript_EOF(fp_tree_ps);
  printf("\n\n");
  
  return -1;
}


#elif KARIN

int main(int argc, char **argv) 
{
  char **tax_set,**out_tax;
  char *s_tree;
  int i;
  arbre *annot_tree;
  arbre *fit_tree;
  edge *root_edge;
  int w,s;
  FILE *fit_fp;
  char **node_labels;


  tax_set = (char **)mCalloc(5,sizeof(char *));
  For(i,5) tax_set[i] = (char *)mCalloc(20,sizeof(char));

  out_tax = (char **)mCalloc(5,sizeof(char *));
  For(i,5) out_tax[i] = (char *)mCalloc(20,sizeof(char));

  strcpy(tax_set[0],"ENSDARP\0");
  strcpy(tax_set[1],"ENSGACP\0");
  strcpy(tax_set[2],"ENSORLP\0");
  strcpy(tax_set[3],"GSTENP\0");
  strcpy(tax_set[4],"SINFRUP\0");

  strcpy(out_tax[0],"ENSMUSP\0");
  strcpy(out_tax[1],"ENSP\0");

  KARIN_Analyze_Output(argv[1],argv[2]);

/*   root_edge = KARIN_Find_Root_Edge(out_tax,2,tree); */
/*   if(root_edge)  */
/*     { */
/*       Add_Root(root_edge,tree); */
/*       KARIN_Find_Duplication_Node(tree->n_root,tree->n_root->v[0],tax_set,5,tree); */
/*       KARIN_Find_Duplication_Node(tree->n_root,tree->n_root->v[1],tax_set,5,tree); */
/*       w = s = 0; */
/*       For(i,2*tree->n_otu-2) if(!strcmp(tree->noeud[i]->n_label,"[W]\0")) w=1; */
/*       For(i,2*tree->n_otu-2) if(!strcmp(tree->noeud[i]->n_label,"[S]\0")) s=1; */
/*       if(w+s==2) printf("AMBIGUOUS "); */
/*       else if(w==1) printf("YES "); */
/*       else printf("NO "); */
/*     } */
/*   else */
/*     { */
/*       printf("CHECK OUTGROUP "); */
/*     } */

  

  return -1;
}

/*********************************************************/


#endif


