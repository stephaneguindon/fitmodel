/*

FitModeL :  a program that Fits substitution models to the data
using Maximum Likelihood.  

Copyright (C) Stephane Guindon. jul. 2004 onward

All parts of  the source except where indicated  are distributed under
the GNU public licence.  See http://www.opensource.org for details.

*/
#include <config.h>


#ifndef UTILITIES_H
#define UTILITIES_H


#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include <time.h>
#include <limits.h>
#include <errno.h>
#include <float.h>
#include <stdbool.h>


#define YES 1
#define NO 0

#define LOG log
#define POW pow
#define EXP exp
#define FABS fabs
#define SQRT sqrt
#define CEIL ceil
#define FLOOR floor
#define RINT rint
#define ROUND round
#define TRUNC trunc
#define COS cos
#define SIN sin
#define TAN tan
#define SMALL DBL_MIN
#define BIG  DBL_MAX
#define SMALL_PIJ 1.E-20
#define LOGBIG 690.
#define LOGSMALL -690.

#ifdef FITMODEL
#define PROGNAME "FitModeL"
#elif EVOLVE
#define PROGNAME "Evolve"
#elif RF
#define PROGNAME "RF"
#elif TOOL
#define PROGNAME "tool"
#elif POSSELSITEID
#define PROGNAME "posselsiteid"
#elif LUMP
#define PROGNAME "lump"
#elif COLTREE
#define PROGNAME "coltree"
#else
#define PROGNAME "UNKNOWN"
#endif

#define RELEASE "v0.5.3"

typedef	double fit_double;
typedef	float  fit_float;

#define ADAPTIVE_CONTROL     0
#define NEWTON_ETAL          1  /* Newton, Noueiry, Sarkar and Ahlquist's method for multiple testing.
                                   Biostatistic, 2004, 5 p155-176 */

#define SMALL DBL_MIN

#define FIXED_REGION         2
#define PARAMETRIC_BOOTSTRAP 3
#define NODE_DEG_MAX        50

#define MAX(a,b)                     ((a)>(b)?(a):(b))
#define MIN(a,b)                     ((a)<(b)?(a):(b))

#define JC69       1
#define K80        2
#define F81        3
#define HKY85      4
#define F84        5
#define TN93       6
#define GTR        7

#define DAYHOFF   11
#define JTT       12
#define MTREV     13
#define WAG       14
#define DCMUT     15

#define M0        21
#define M1        22
#define M1a       23
#define M2        24
#define M2a       25
#define M3        26
#define MX        27

#define NO_SWITCH 31
#define SWITCH_S1 32
#define SWITCH_S2 33

#define NT        1
#define AA        2
#define CODONS    3

#define base_A    0
#define base_C    1
#define base_G    2
#define base_T    3

#define N_MAX_CATQ  20
#define N_MAX_OMEGA 10

#define  BRENT_ITMAX          100            
#define  BRENT_CGOLD    0.3819660            
#define  BRENT_ZEPS        1.e-10            
#define  MNBRAK_GOLD     1.618034            
#define  MNBRAK_GLIMIT      100.0            
#define  MNBRAK_TINY       1.e-20            
#define  ALPHA_MIN           0.04            
#define  ALPHA_MAX            100            
/* #define  BL_MIN            1.e-10 */
#define  BL_MIN            1.e-5
#define  BL_START          1.e-03            
#define  BL_MAX             100.0
/* #define  MIN_DIFF_LK       1.e-04             */
#define  MIN_DIFF_LK       1.e-03            
#define  GOLDEN_R      0.61803399            
#define  GOLDEN_C  (1.0-GOLDEN_R)            
#define  T_MAX_FILE           200            
#define  T_MAX_LINE       1000000            
#define  T_MAX_NAME           100            
#define  T_MAX_MODEL_NAME     100            
#define  T_MAX_SEQ        1000000            
#define  N_MAX_INSERT          20            
#define  N_MAX_OTU           1000            
#define  UNLIKELY          -1.e10            
#define  NJ_SEUIL             0.1            
#define  ROUND_MAX            100            
#define  DIST_MAX              4.            
#define  AROUND_LK           50.0            
#define  PROP_STEP            1.0            
#define  T_MAX_ALPHABET      1000            
#define  MDBL_MIN   2.225074E-308            
#define  MDBL_MAX   1.797693E+308            
#define  POWELL_ITMAX         200            
#define  LINMIN_TOL       2.0E-04            
#define  LIM_SCALE              3
#define  LIM_SCALE_VAL     1.E-50
/* #define  LIM_SCALE            3000 */
/* #define  LIM_SCALE_VAL     1.E-200 */
#define  T_MAX_LABEL           10
#define  N_MAX_LABEL           10
#define  BLOCK_LABELS         100


#define For(i,n)                     for(i=0; i<n; i++)
#define Fors(i,n,s)                  for(i=0; i<n; i+=s)
#define PointGamma(prob,alpha,beta)  PointChi2(prob,2.0*(alpha))/(2.0*(beta))
#define SHFT2(a,b,c)                 (a)=(b);(b)=(c);
#define SHFT3(a,b,c,d)               (a)=(b);(b)=(c);(c)=(d);
#define MMAX(a,b)                    ((a)>(b)?(a):(b))
#define MMIN(a,b)                    ((a)<(b)?(a):(b))
#define SIGN(a,b)                    ((b) > 0.0 ? fabs(a) : -fabs(a))
#define SHFT(a,b,c,d)                (a)=(b);(b)=(c);(c)=(d);

#ifdef COLTREE
struct __Arbre *tree;
#endif

/*********************************************************/

typedef struct __Node {

  struct __Node                **v; /* neighbours */
  struct __Edge                **b; /* adjacent edges */
  struct __Node            ***sons; /* three lists of pointer to descendant nodes. Each list corresponds to one direction */
  struct __Node        ***bip_node; /* three lists of pointer to tip nodes. One list for each direction */  
  struct __Node               *anc;

  fit_double                    *l; /* lenghts of adjacent edges */
  fit_double              *l_prime; /* another branch label */
  fit_double                  tpos; /* time position */
  fit_double             *tpos_min; /* min of time positions in each direction */
  fit_double              tpos_old; /* old time position */
  fit_double        tpos_lim_below; /* tpos_lim_below and t_pos_lim_above both define */
  fit_double        tpos_lim_above; /* the upper and lower bound for tpos */
  fit_double          dist_max_tip; /* distance to the most distant  tip */
  fit_double          dist_to_root; /* distance to root */
  fit_double               coord_x;
  fit_double               coord_y;
  fit_double           *post_w_ave; /* lenghts of adjacent edges */
  

  int                          num; /* node number */
  int                          tax; /* tax=1 -> this node is a tip */
  int                     ni,agglo; /* these values are used during the estimation of a tree using BIONJ */
  int                infered_state; /* ancestral state estimate */
  int                      *n_sons; /* size of each of the three lists in sons */
  int                    *bip_size; /* Size of each of the three lists from bip_node */ 
  int                 check_branch;
  int                non_root_dir1;
  int                non_root_dir2;
  int                     root_dir;


  char                       *name; /* name='NULL' is internal node, name='corresponding sequence name', otherwise */
  char                 ***bip_name; /* three lists of tip node names. One list for each direction */
  char                        *seq; /* (ancestral) sequence */  
  char                    *n_label; /* node label */
}node;


/*********************************************************/

typedef struct __Edge {
  struct __Node           *left,*rght; /* nodes located at the two extremities */ 
  struct __Qmat          *qmat_struct; /* rate matrix associated to this edge */
  struct __NNI                 *s_nni; /* used to perform NNI */
  
  fit_double                            l; /* edge length */
  fit_double                      l_prime; /* edge length (used for probabilities of selection classes) */
  fit_double                          l_old; /* edge length */
  fit_double                         l_ml; /* ML edge length */
  fit_double                      diff_lk; /* used in tree topology search based on NNI */
  fit_double  ****p_lk_left,****p_lk_rght; /* partial likelihoods on the left/right of this edge, at each site */
  fit_double                   ****Pij_rr; /* probabilities of change */
  fit_double        site_sum_scale_f_left; /* sum of the scaling factors on the left, at a given site */
  fit_double        site_sum_scale_f_rght; /* sum of the scaling factors on the right */
  fit_double            site_scale_f_left; /* scaling factor on the left */
  fit_double            site_scale_f_rght; /* scaling factor on the right */
  fit_double            *sum_scale_f_left; /* sum of the scaling factors on the left over the whole set of sites */
  fit_double            *sum_scale_f_rght; /* sum of the scaling factors on the left over the whole set of sites */
  
  int     l_r,r_l,l_v1,l_v2,r_v1,r_v2; /* directions */ 
  int                             num; /* edge number */
  int                  check_this_one; /* check this edge ? */
  int                       best_conf; /* used in tree topology search based on NNI */
  int                     num_st_left; /* number of the subtree on the left */
  int                     num_st_rght; /* number of the subtree on the right */
  int                   get_p_lk_left; /* compute the partial likelihood on the left of this edge ? */
  int                   get_p_lk_rght; /* compute the partial likelihood on the right of this edge ? */
  int                    ud_p_lk_left; /* update the partial likelihood on the left of this edge ? */
  int                    ud_p_lk_rght; /* update the partial likelihood on the right of this edge ? */
  int               substitution_type; /* -1 : a synonymous substitution may have occured on this edge */
  /*  0 : no substitution */
  /* +1 : a nonsynonymous substitution may have occured on this edge */
  int               n_nonsyno_on_left; /* number of nonsynonymous substitutions on the left */
  int               n_nonsyno_on_rght; /* number of nonsynonymous substitutions on the right */
  int                       bip_score;
  
  
  char                       **labels; /* string of characters that labels the corresponding edge */
  int                        n_labels; /* number of labels */
  
  fit_double                     dwell[3]; /* dell times under three different selection regimes */
  fit_double         dwell_probs[3][3][3]; /* matrices of dwell times in the three selection regimes 
					      along a given branch */
  fit_double              prob_sel_regime;
  fit_double            ***integral;
}edge;

/*********************************************************/

typedef struct __Arbre { /* 'Arbre' -> 'tree' in french... */
  struct __Node                      *n_root; /* if != NULL -> tree is rooted */
  struct __Node                    **noeud; /* list of pointers to nodes */
  struct __Edge                  **t_edges; /* list of pointers to edges */
  struct __Model                      *mod; /* substitution model */
  struct __AllSeq                    *data; /* compressed sequences */
  struct __Option                   *input;
  struct __Num_Parameters *numerical_param; /* numerical parameters that are used to compute the likelihood */
  struct __Edge                    *e_root;
 
  fit_double          tot_loglk; /* log likelihood */
  fit_double           *site_lk; /* log likelihood site by site */
  fit_double         *loglk_cat;
  fit_double *site_loglk_sorted; /* log likelihood site by site (sorted in increasing order) */
  fit_double    unconstraint_lk; /* multinomial model */
  fit_double  **sel_regime_prob; /* posterior probabilities or of the different selection classes
				    or likelihood of each selection regime, at each site */
  fit_double         n_root_pos;
  
  int                n_swap; /* used in tree topology search algorithm */
  int             n_pattern; /* number of site patterns (can be nucleotide, aa or codon patterns) */
  int               *w_patt; /* pattern weights */
  int            both_sides; /* if 1 -> compute the partial likelihoods of every subtree */
  int                 n_otu; /* number of taxa */
  int             curr_site; /* current site */
  int               has_bip;
  int   *n_changes_per_site; /* (minimum) number of substitutions per site */
  int            param_size; 
  int             param_num;
  int               verbose;
  int                    font;
  int        render_tax_names;
  int      render_edge_colors;
  int             tree_number;
  int      has_branch_lengths;
  int num_curr_branch_available;

  FILE *fp_tree;
  
}arbre;


/*********************************************************/

typedef struct __Seq { /* uncompressed sequences */
  char  *name; /* sequence ID */
  int     len; /* sequence length */
  char *state; /* sequence itself */
}seq;

/*********************************************************/


typedef struct __AllSeq { /* compressed sequences */
    struct __Model  *mod; /* substitution model */
    struct __Seq **c_seq; /* crunched sequences      */

    fit_double    *b_frq; /* state frequencies */

    int            *wght; /* weights of every site pattern */
    short int     *invar; /* 1 -> no substitutions are observed at this site */
    int            n_otu; /* number of sequences */ 
    int        clean_len; /* uncompressed sequence lenghts without gaps */
    int       crunch_len; /* compressed sequence lengths */
    int         init_len; /* uncompressed sequence length */
    short int    *ambigu; /* ambiguity characters */
    int        **pospatt; /* record of the positions of the different site patterns in the uncompressed data set */ 
    int        *selclass; /* selection class (estimated or true) at each site */
}allseq;

/*********************************************************/

typedef struct __Model {
  struct __Code          *c_code; /* genetic code structure (see below) */
  struct __Qmat    **qmat_struct; /* array of pointers to instantaneous rate matrices structure (see below) */
  struct __Optimiz        *s_opt; /* pointer to optimisation structure (see below) */

  fit_double                 *pi; /* state frequencies */
  fit_double               kappa; /* transition/transversion rate */
  fit_double              lambda; /* lambda = (ts/tv ratio btw purines) / (ts/tv ratio btw pyrimidines) */
  fit_double               alpha; /* gamma shape parameter */
  fit_double        *rr_mixturem; /* relative rates of substitution for mixture models (e.g. discrete gamma distribution) */
  fit_double            *r_proba; /* prior probabilities of the relative rates of substitution */
  fit_double              pinvar; /* proportion of invariant sites */
  fit_double          *gtr_param; /* relative rate parameters of the GTR model */
  fit_double                  mr; /* mean rate  of substitution */
  fit_double              *omega; /* dn/ds ratio values */
  fit_double        *omega_proba; /* equilibrium frequencies of the different dn/ds ratio values */
  fit_double          *omega_min; /* minimum of dn/ds ratio values */
  fit_double          *omega_max; /* maximum of dn/ds ratio values */

  fit_double         thresholdp2; /* threshold value for classification of sites based on estimated selection classes */ 
  fit_double                 fdr; /* false detection rate */
  fit_double        delta_switch; /* value of the switching parameter */
  fit_double     io_delta_switch; /* user-defined value of the switching parameter */
  fit_double        *rsubst_rate; /* relative rate of substitution */

  int              n_rsubst_rate; /* number of relative rate classes */
  int               model_number; /* model number 1 => JC69 2 => K2P 3 => F81 
				     4 => HKY85 5 => F84 6 => TN93 7 => GTR 
				     11 => Dayhoff 12 => JTT 13 => MtREV 
				     20 => M2/M3 21 => M2+S1/M2+S2/M3+S1/M3+S2 */
  int            subst_modelname;
  int           switch_modelname;
  int           model_applies_to; /* NT, AA or CODONS */
  int                         ns; /* number of states (4 for ADN, 20 for AA, 60 to 63 for codons (depending on the genetic code) */
  int                   ns_codon; /* number of sense codons */
  int                      ns_nt; /* number of nucleotides */
  int                   datatype; /* 0->DNA, 1->AA */
  int                     n_catg; /* number of classes of the discrete gamma distribution */
  int                      invar; /* invar=1 -> no substitutions are observed at this site */  
  int               update_eigen; /* update_eigen=1 -> the instantaneous rate matrix must be inverted */
  int                      n_otu; /* number of txta */
  int                   stepsize; /* stepsize=1 for DNA and amino acid models, stepsize=0 for codon models */ 
  int              n_qmat_struct; /* number of instantaneous rate matrices to allocate */
  int                     n_catq; /* number of distinct instantaneous rate matrices */
  int                    n_omega; /* number of distinct dn/ds ratio values */
  int                 analytical; /* analytical=1 -> change probabilities (Pij(t)) are computed analytically */
  int                 codon_freq; /* 0 -> 1/61  ; 1 -> Freq(XYZ) = Freq(X on first position)
				     x Freq(Y on secnd position)
				     x Freq(Z on third position) */
  
  fit_double         subst_rate;  /* deprecated */
  int                  tpos_ols;  /* deprecated */
  int      update_bl_using_tpos;  /* deprecated */
  
  fit_double    *selec_reg_freq;
  int	          n_catq_negsel;
  int	          n_catq_possel;
  int	         n_catq_neutral;

}model;

/*********************************************************/

typedef struct __Option {
    struct __Model      *mod; /* pointer to the model structure */
    struct __Arbre     *tree; /* pointer to the tree structure */
    struct __Seq      **data; /* pointer to the sequences */
    struct __AllSeq *alldata; /* pointer to the compressed sequences */
    
    FILE             *fp_seq; /* pointer to the sequence file */
    FILE            *fp_tree; /* pointer to the tree file */
    FILE     *fp_output_tree; /* pointer to the output tree file */
    FILE     *fp_output_stat; /* pointer to the output stats file */
    FILE   *fp_ancestral_seq;

    int             interleaved; /* interleaved=1=>PHYLIP interleaved format, 0=>sequential format */ 
    int               inputtree; /* is there any tree as input ?... 1=>Yes */
    int     stat_file_open_mode; /* 1=>replace file, 2=>append to file */
    int     tree_file_open_mode; /* 1=>replace file, 2=>append to file */  
    int             n_data_sets; /* number of data sets to analyse */
    int        n_data_set_asked; /* deprecated (only used if EVOLVE=1) */
    int                   n_otu; /* number of taxa */
    int                user_len; /* deprecated (only used if EVOLVE=1) */
    int       start_at_data_set;
    int pos_sel_sites_id_method; /* positively selected sites identification method */
    int        sel_class_target;


    char                   *seqfile; /* sequence file name */
    char           *input_tree_file; /* input tree file name */ 
    char          *output_tree_file; /* output tree file name */
    char          *output_stat_file; /* output stat file name */
    char          *use_default_dnds;  
    char *use_default_ssvs_switches;
}option;

/*********************************************************/

typedef struct __Optimiz {
    int             print; /* 1=>verbose mode */
    int         opt_alpha; /* optimise gamma shape parameter ?   */
    int         opt_kappa; /*    "     kappa                     */
    int        opt_lambda; /*    "     lambda                    */
    int         opt_omega; /*    "     dn/ds ratio               */
    int         opt_theta; /*    "     theta                     */
    int        opt_pinvar; /*    "     pinvar                    */
    int            opt_bl; /*    "     branch lengths            */     
    int         opt_bfreq; /*  Get ML estimates of base frequencies ? */
    int         opt_param; /* optimise substitution parameters ? */
    int          opt_tpos; /* deprecated */
    int    opt_subst_rate; /*     "      */
    int       opt_p_omega; /* optimise the switching rate beta   */
    int        sort_omega;

    fit_double     last_alpha; /* last value of the step in the linmin search */
}optimiz;

/*********************************************************/

typedef struct __Code {
  char               *aa; /* table of amino acid symbols */
  char             *name; /* name of that code */ 
  
  int           **genetc; /* map that relates codons to amino acids characters */
  int *n_diff_b_2_codons; /* number of differences between two codons */
  int         *tstvtable; /* tstvtable[i*61+j]=2 if at least one transversional 
			     change occured between codons i & j, 
			     =1 if transition, =0 otherwise */
  int      num_curr_code; /* genetic code number */
  int           *sense_c; /* array of sense codons */
  int          n_sense_c; /* number of sense codons */
  int      *from_64_2_61; /* this table is used to get rid of non-sense codons */
  fit_double   **gtr_4_codon; /* table used for the GTR version of codon models */
  int    number_of_codes;
}code;

/*********************************************************/

typedef struct __Qmat{
  fit_double           *u_mat; /* right eigen vectors */
  fit_double           *v_mat; /* left eigen vectors = inv(u_mat) */
  fit_double        *root_vct; /* eigen values */  
  fit_double     *expD_mr_vct; /* exp(root_vect) */
  fit_double            *qmat; /* instantaneous rate matrix */
  
  fit_double      *qmat_proba; /* instantaneous rate matrix */
  int       curr_qmat_cat; /* current instantaneous rate matrix */
  int              n_qmat; /* number of instantaneous rate matrices */
  
  int             n_omega; /* number of dn/ds ratio values */
  fit_double           *omega; /* dn/ds ratio values */
  fit_double       *omega_min; /* dn/ds ratio values */
  fit_double       *omega_max; /* dn/ds ratio values */
  fit_double        **t_omega; /* dn/ds ratio values used to fill the instantaneous rate matrices */
  fit_double     *omega_proba; /* equilibrium frequencies of dn/ds ratio values */
  fit_double  **t_omega_proba; /* equilibrium frequencies of dn/ds ratio values used to fill instantaneous rate matrices */
  fit_double           *theta; /* switching parameters */
  
  fit_double  *trans_qmat_proba;
  fit_double *trans_omega_proba;
}qmat;

/*********************************************************/

typedef struct __NNI{/* deprecated */
  /*
     v1              v1                 v1                    
     |		     |         	        |         
     |		     |         	        |         
     v2		     v2        	        v2        
    /  \	    /  \       	       /  \       
   /    \	   /    \      	      /    \      
  v3     \	  v5     \     	     v6     \     
         v4	         v4    	            v4    
        /  \	        /  \   	           /  \   
       /    \	       /    \  	          /    \  
      v5    v6        v3    v6 	         v5    v3 
  */    

  struct __Edge                 *b_fcus;
  struct __Node *v1,*v2,*v3,*v4,*v5,*v6;
  int                         best_conf;
  fit_double                       *bl_info;
  fit_double                   bl_info_init;
}nni;


/*********************************************************/

typedef struct __Num_Parameters {

    fit_double     *param_val;
    int       *param_size;
    int        *param_beg;
    fit_double *one_param_val;

    int        br_len_num;
    int            pi_num;
    int         kappa_num;
    int        lambda_num;
    int         alpha_num;
    int        pinvar_num;
    int         omega_num;
    int         theta_num;
    int       p_omega_num;
    int        p_qmat_num;
    int           gtr_num;

    int     currently_opt;
    int replace_num_param;
}numpar;

/*********************************************************/
/*********************************************************/
/*********************************************************/

fit_double bico(int n,int k);
fit_double factln(int n);
fit_double gammln(fit_double xx);
fit_double Pbinom(int N,int ni,fit_double p);
void Plim_Binom(fit_double pH0,int N,fit_double *pinf,fit_double *psup);
fit_double LnGamma(fit_double alpha);
fit_double IncompleteGamma(fit_double x,fit_double alpha,fit_double ln_gamma_alpha);
fit_double PointChi2(fit_double prob,fit_double v);
fit_double PointNormal(fit_double prob);
int DiscreteGamma(fit_double freqK[],fit_double rK[],fit_double alfa,fit_double beta,int K,int median);
arbre *Read_Tree_File(FILE *fp_tree);
arbre *Read_Tree(char *s_tree);
void Make_All_Edges_Light(node *a,node *d);
void Make_All_Edges_Lk(node *a, node *d, edge *b, arbre *tree);
void R_rtree(char *s_tree_a, char *s_tree_d, node *a, arbre *tree, int *n_int, int *n_ext);
char **Sub_Trees(char *tree, int *degree);
int Next_Par(char *s,int pos);
void Test_Tree(char *s);
char *Write_Tree(arbre *tree);
void R_wtree(node *pere, node *fils, char *s_tree, arbre *tree);
void Init_Tree(arbre *tree, int n_otu);
edge *Make_Edge_Light(node *a, node *d, int num);
void Make_NNI(nni *this1);
void Init_Edge_Light(edge *b, int num);
void Make_Edge_Dirs(edge *b,node *a,node *d);
void Make_Edge_Lk(node *a, node *d, edge *b, arbre *tree);
node *Make_Node_Light();
void Init_Node_Light(node *n);
void Make_Node_Lk(arbre *tree, node *n);
seq **Get_Seq(option *input,int rw);
seq **Read_Seq_Sequential(FILE *in,int *n_otu);
seq **Read_Seq_Interleaved(FILE *in,int *n_otu);
int Read_One_Line_Seq(seq ***data,int num_otu,FILE *in);
void Uppercase(char *ch);
allseq *Compact_Seq(seq **data,option *input);
allseq *Compact_CSeq(allseq *data,model *mod);
void Get_Codon_Freqs(allseq *data,model *mod);
void Get_Base_Freqs(allseq *data,int position,int step);
void Get_AA_Freqs(allseq *data);
void Init_Tree_Edges(node *a,node *d,arbre *tree,int *cur);
void Exit(char *message);
void *mCalloc(int nb,size_t size);
void *mRealloc(void *p,int nb,size_t size);
arbre *Make_Light_Tree_Struct(int n_otu);
int Sort_Float_Decrease(const void *a,const void *b);
void Print_Site(allseq *alldata,int num,int n_otu,char *sep,int stepsize);
void Print_CSeq(FILE *fp,allseq *alldata);
void Order_Tree_Seq(arbre *tree,seq **data);
void Order_Tree_CSeq(arbre *tree,allseq *data);
int Sort_String(const void *a,const void *b);
void Get_SubTrees_Numbers(arbre *tree);
optimiz *Alloc_Optimiz();
void Init_Optimiz(optimiz *s_opt);
int Filexists(char *filename);
FILE *Openfile(char *filename,int mode);
void Print_Fp_Out(FILE *fp_out,time_t t_beg,time_t t_end,arbre *tree,option *input,int n_data_set);
void Alloc_All_P_Lk(arbre *tree);
int Is_Ambigu(char *state,int datatype,int stepsize);
int Is_Completely_Ambigu(char *state,int datatype,int stepsize);
void Check_Ambiguities(allseq *data,int datatype,int stepsize);
int Assign_State(char *c,int datatype,int stepsize,code *c_code);
void Update_BrLen_Invar(arbre *tree);
void Getstring_Stdin(char *file_name);
allseq *Evolve(arbre *tree);
void Generate_One_Site(node *a,node *d,int anc_state,edge *b,int rrate,int qmat_cat,allseq *simu_seq,int *n_switches, arbre *tree);
int Generate_One_State(int ns,fit_double *p);
void Return_State(arbre *tree,int num_state,char *state);
int Get_AA(code *c_code,int num_codon);
code *Make_Code(int num_code);
void Init_Nb_Diff_Between_Two_Codons(code *c_code,int ns);
void Init_Tstv_Table(code *c_code,int ns);
int Compare_Two_States(char *state1,char *state2,int state_size);
void Copy_One_State(char *from,char *to,int state_size);
qmat *Make_Qmat_Struct(int ns,int n_qmat,int n_omega,model *mod);
void Copy_Qmat_Struct(qmat *from,qmat *to);
char *Get_Codon(int num_codon);
char Get_Base(int num_base);
void Init_Codons(code *c_code);
char *Print_State(int num,arbre *tree);
void Infer_Synonymous_Substitutions(edge *b_fcus,node *a,node *d,arbre *tree);
void Infer_NonSynonymous_Substitutions(edge *b_fcus,node *a,node *d,arbre *tree);
void Make_Sons(arbre *tree);
void Make_Sons_Post(node *a,node *d,arbre *tree);
void Make_Sons_Pre(node *a,node *d,arbre *tree);
int Is_Son(node *a,node *x,int dir);
void Distance_Between_2_nodes(edge *b_fcus,node *a,node *d,node *final,fit_double *distance);
void qksort(fit_double *A,int ilo,int ihi);
void Update_Tpos_Min(arbre *tree);
void Update_Tpos_Min_Post(node *a,node *d);
void Update_Tpos_Min_Pre(node *a,node *d);
void Tpos_OLS(arbre *tree,node *a,node *d);
void Update_Bl_Using_Tpos(arbre *tree,node *a,node *d);
void Set_Defaults_Model(model *mod);
void Init_NNI(node *a,node *d,edge *b_fcus,arbre *tree);
void Check_Tpos(node *a,node *d);
void Add_Root(edge *target, arbre *tree);
void Remove_Root(arbre *tree);
model *Make_Model_Basic();
void Make_Model_Complete(model *mod);
void Make_All_Qmat_Struct(model *mod);
void Print_Param(FILE *fp, arbre *tree);
void Check_Memory_Amount(arbre *tree);
void Hide_Ambiguities(allseq *data);
int Are_Compatible(char *statea, char *stateb, int stepsize, int datatype);
void Alloc_Bip(arbre *tree);
void Compare_Bip(arbre *tree1, arbre *tree2);
void Get_Bip(node *a, node *d, arbre *tree);
fit_double Get_Proportion_Of_False_Positives(fit_double threshold, fit_double *postp2, int *selclass, int size, int *num, int *denom, int print);
fit_double Get_Proportion_Of_False_Negatives(fit_double threshold, fit_double *postp2, int *selclass, int size, int *num, int *denom, int print);
void Print_Ancestral_Seq(node *a,node *d,arbre *tree);
void Estimate_Ancestral_States(arbre *tree);
void Estimate_Ancestral_States_Post_Codon(node *a, node *d, edge *b, arbre *tree);
void Estimate_Ancestral_States_Post_NtAA(node *a, node *d, edge *b, arbre *tree);
void Print_Tree_Nodes(node *a, node *d, arbre *tree);
fit_double Rf(arbre *tree1, arbre *tree2);
void Init_GTR_4_Codon_Table(model *mod);
void Evolve_Under_H0(arbre *tree);
void Print_Seq(seq **data, int n_otu, FILE *fp);
void Allocate_Num_Parameters(arbre *tree);
void Init_Num_Parameters(arbre *tree);
void Replace_Num_Parameters(arbre *tree);
void Init_Lagrange_Qmat_Probs(arbre *tree);
void Allocate_Lagrange(arbre *tree);
void Init_Lagrange_Qmat_Probs(arbre *tree);
void Update_Ancestors(node *a, node *d);
void Unroot_Tree(char **subtrees);
void Clean_Multifurcation(char **subtrees, int current_deg, int end_deg);

fit_double Num_Derivatives_One_Param(fit_double (*func)(arbre *tree), arbre *tree,
				 fit_double f0, fit_double *param, fit_double stepsize,
                                     fit_double *err, int precise);

int Num_Derivative_Several_Param(arbre *tree, fit_double *param, int n_param, fit_double stepsize,
                                 fit_double (*func)(arbre *tree), fit_double *derivatives);

fit_double Lk_Arg_Edge_Len(fit_double *param_val,
                       int br_num,
                       arbre *tree);

void Untransform_Probs(fit_double *t_probs, /* transformed probabilities */ 
                       fit_double *probs,   /* 'real' probabilities      */ 
                       int n);

int ACGT_to_RY(int base_acgt);
void Print_Translated_Alignment_DNA_to_AA(FILE *fp, allseq *alldata, code *c_code);
int Get_Base_Num(char base);
void ChangeSize(int w, int h);
void OrientMe(float ang);
void MoveMeFlat(int i);
void Print_Model_Param(arbre *tree);
void Sort_Categories(arbre *tree);
void Get_Score_Mat(code *c_code);
void Get_Rid_Of_Prefix(char delim, arbre *tree);
int Is_Duplication_Node(node *a, node *d, char **tax_set, int n_tax, arbre *tree);
void KARIN_Find_Duplication_Node(node *a, node *d, char **tax_set, int n_tax, arbre *tree);
edge *KARIN_Find_Root_Edge(char **out_tax, int n_out_tax, arbre *tree);
void Make_All_Tree_Edges(arbre *tree);
void Make_All_Tree_Nodes(arbre *tree);
void Connect_One_Edge_To_Two_Nodes(node *a, node *d, edge *b, arbre *tree);
arbre *Make_Tree(int n_otu);
void Read_Branch_Length(char *s_d, char *s_a, arbre *tree);
void Read_Branch_Label(char *s_d, char *s_a, edge *b);
void Read_Node_Name(node *d, char *s_tree_d, arbre *tree);
void Make_New_Edge_Label(edge *b);
char **KARIN_Branch_Labels_To_Node_Labels(arbre *tree);
int KARIN_Find_Match_Nodes(node *n, arbre *tree1, arbre *tree2);
void KARIN_Analyze_Output(char *annot_tree_file, char *fit_tree_file);
void KARIN_Return_DnDs_Subtree(node *n, double *sum_w, int *n_edges, arbre *tree);
void Add_Ambiguities(allseq *ori_seq, allseq *simu_seq, int datatype);
allseq *Copy_Cseq(allseq *ori, int len, int ns);
allseq *Make_Cseq(int n_otu, int crunch_len, int init_len, char **sp_names);


#endif
