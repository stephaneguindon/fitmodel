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

/* #define ERANGE 1 */
extern int errno;

/*********************************************************/
/* NUMERICAL RECIPES ROUTINES FOR COMPUTING C(n,k)       */
/*********************************************************/

fit_double bico(int n, int k)
{
  return (fit_double)floor(0.5+(fit_double)exp(factln(n)-factln(k)-factln(n-k)));
}

/*********************************************************/

fit_double factln(int n)
{
   static fit_double a[101];

   if (n < 0){ Exit("Err: negative factorial in routine FACTLN"); }
   if (n <= 1) return 0.0;
   if (n <= 100) return a[n] ? a[n] : (a[n]=gammln(n+1.0));
   else return gammln(n+1.0);
}

/*********************************************************/

fit_double gammln(fit_double xx)
{
   fit_double x,tmp,ser;
   static fit_double cof[6]={76.18009173,-86.50532033,24.01409822,
      -1.231739516,0.120858003e-2,-0.536382e-5};
   int j;

   x=xx-1.0;
   tmp=x+5.5;
   tmp -= (x+0.5)*(fit_double)log(tmp);
   ser=1.0;
   for (j=0;j<=5;j++) {
      x += 1.0;
      ser += cof[j]/x;
   }
   return -tmp+(fit_double)log(2.50662827465*ser);
}

/*********************************************************/
/*          END OF NUMERICAL RECIPES ROUTINES            */
/*********************************************************/

fit_double Pbinom(int N, int ni, fit_double p)
{
  return bico(N,ni)*(fit_double)pow(p,ni)*(fit_double)pow(1-p,N-ni);
}

/*********************************************************/

void Plim_Binom(fit_double pH0, int N, fit_double *pinf, fit_double *psup)
{
  *pinf = pH0 - 1.64*(fit_double)sqrt(pH0*(1-pH0)/(fit_double)N);
  if(*pinf < 0) *pinf = .0;
  *psup = pH0 + 1.64*(fit_double)sqrt(pH0*(1-pH0)/(fit_double)N);
}

/*********************************************************/

fit_double LnGamma (fit_double alpha)
{
/* COURTESY OF ZIHENG YANG
   returns ln(gamma(alpha)) for alpha>0, accurate to 10 decimal places.  
   Stirling's formula is used for the central polynomial part of the procedure.
   Pike MC & Hill ID (1966) Algorithm 291: Logarithm of the gamma function.
   Communications of the Association for Computing Machinery, 9:684
*/
   fit_double x=alpha, f=0, z;

   if (x<7) {
      f=1;  z=x-1;
      while (++z<7)  f*=z;
      x=z;   f=-(fit_double)log(f);
   }
   z = 1/(x*x);
   return  f + (x-0.5)*(fit_double)log(x) - x + .918938533204673 
	  + (((-.000595238095238*z+.000793650793651)*z-.002777777777778)*z
	       +.083333333333333)/x;  
}

/*********************************************************/

fit_double IncompleteGamma (fit_double x, fit_double alpha, fit_double ln_gamma_alpha)
{
/* COURTESY OF ZIHENG YANG
   returns the incomplete gamma ratio I(x,alpha) where x is the upper 
	   limit of the integration and alpha is the shape parameter.
   returns (-1) if in error
   ln_gamma_alpha = ln(Gamma(alpha)), is almost redundant.
   (1) series expansion     if (alpha>x || x<=1)
   (2) continued fraction   otherwise
   RATNEST FORTRAN by
   Bhattacharjee GP (1970) The incomplete gamma integral.  Applied Statistics,
   19: 285-287 (AS32)
*/
   int i;
   fit_double p=alpha, g=ln_gamma_alpha;
   fit_double accurate=1e-8, overflow=1e30;
   fit_double factor, gin=0, rn=0, a=0,b=0,an=0,dif=0, term=0, pn[6];

   if (x==0) return (0);
   if (x<0 || p<=0) return (-1);

   factor=(fit_double)exp(p*(fit_double)log(x)-x-g);   
   if (x>1 && x>=p) goto l30;
   /* (1) series expansion */
   gin=1;  term=1;  rn=p;
 l20:
   rn++;
   term*=x/rn;   gin+=term;

   if (term > accurate) goto l20;
   gin*=factor/p;
   goto l50;
 l30:
   /* (2) continued fraction */
   a=1-p;   b=a+x+1;  term=0;
   pn[0]=1;  pn[1]=x;  pn[2]=x+1;  pn[3]=x*b;
   gin=pn[2]/pn[3];
 l32:
   a++;  b+=2;  term++;   an=a*term;
   for (i=0; i<2; i++) pn[i+4]=b*pn[i+2]-an*pn[i];
   if (pn[5] == 0) goto l35;
   rn=pn[4]/pn[5];   dif=(fit_double)fabs(gin-rn);
   if (dif>accurate) goto l34;
   if (dif<=accurate*rn) goto l42;
 l34:
   gin=rn;
 l35:
   for (i=0; i<4; i++) pn[i]=pn[i+2];
   if ((fit_double)fabs(pn[4]) < overflow) goto l32;
   for (i=0; i<4; i++) pn[i]/=overflow;
   goto l32;
 l42:
   gin=1-factor*gin;

 l50:
   return (gin);
}


/*********************************************************/

fit_double PointChi2 (fit_double prob, fit_double v)
{
/* COURTESY OF ZIHENG YANG
   returns z so that Prob{x<z}=prob where x is Chi2 distributed with df=v
   returns -1 if in error.   0.000002<prob<0.999998
   RATNEST FORTRAN by
       Best DJ & Roberts DE (1975) The percentage points of the 
       Chi2 distribution.  Applied Statistics 24: 385-388.  (AS91)
   Converted into C by Ziheng Yang, Oct. 1993.
*/
   fit_double e=.5e-6, aa=.6931471805, p=prob, g;
   fit_double xx, c, ch, a=0,q=0,p1=0,p2=0,t=0,x=0,b=0,s1,s2,s3,s4,s5,s6;

   if (p<.000002 || p>.999998 || v<=0) return (-1);

   g = LnGamma (v/2);
   xx=v/2;   c=xx-1;
   if (v >= -1.24*(fit_double)log(p)) goto l1;

   ch=(fit_double)pow((p*xx*(fit_double)exp(g+xx*aa)), 1/xx);
   if (ch-e<0) return (ch);
   goto l4;
l1:
   if (v>.32) goto l3;
   ch=0.4;   a=(fit_double)log(1-p);
l2:
   q=ch;  p1=1+ch*(4.67+ch);  p2=ch*(6.73+ch*(6.66+ch));
   t=-0.5+(4.67+2*ch)/p1 - (6.73+ch*(13.32+3*ch))/p2;
   ch-=(1-(fit_double)exp(a+g+.5*ch+c*aa)*p2/p1)/t;
   if ((fit_double)fabs(q/ch-1)-.01 <= 0) goto l4;
   else                       goto l2;
  
l3: 
   x=PointNormal (p);
   p1=0.222222/v;   ch=v*(fit_double)pow((x*(fit_double)sqrt(p1)+1-p1), 3.0);
   if (ch>2.2*v+6)  ch=-2*((fit_double)log(1-p)-c*(fit_double)log(.5*ch)+g);
l4:
   q=ch;   p1=.5*ch;
   if ((t=IncompleteGamma (p1, xx, g))<0) {
      printf ("\nerr IncompleteGamma");
      return (-1);
   }
   p2=p-t;
   t=p2*(fit_double)exp(xx*aa+g+p1-c*(fit_double)log(ch));   
   b=t/ch;  a=0.5*t-b*c;

   s1=(210+a*(140+a*(105+a*(84+a*(70+60*a))))) / 420;
   s2=(420+a*(735+a*(966+a*(1141+1278*a))))/2520;
   s3=(210+a*(462+a*(707+932*a)))/2520;
   s4=(252+a*(672+1182*a)+c*(294+a*(889+1740*a)))/5040;
   s5=(84+264*a+c*(175+606*a))/2520;
   s6=(120+c*(346+127*c))/5040;
   ch+=t*(1+0.5*t*s1-b*c*(s1-b*(s2-b*(s3-b*(s4-b*(s5-b*s6))))));
   if ((fit_double)fabs(q/ch-1) > e) goto l4;

   return (ch);
}

/*********************************************************/

fit_double PointNormal (fit_double prob)
{
/* COURTESY OF ZIHENG YANG
   returns z so that Prob{x<z}=prob where x ~ N(0,1) and (1e-12)<prob<1-(1e-12)
   returns (-9999) if in error
   Odeh RE & Evans JO (1974) The percentage points of the normal distribution.
   Applied Statistics 22: 96-97 (AS70)

   Newer methods:
     Wichura MJ (1988) Algorithm AS 241: the percentage points of the
       normal distribution.  37: 477-484.
     Beasley JD & Springer SG  (1977).  Algorithm AS 111: the percentage 
       points of the normal distribution.  26: 118-121.

*/
   fit_double a0=-.322232431088, a1=-1, a2=-.342242088547, a3=-.0204231210245;
   fit_double a4=-.453642210148e-4, b0=.0993484626060, b1=.588581570495;
   fit_double b2=.531103462366, b3=.103537752850, b4=.0038560700634;
   fit_double y, z=0, p=prob, p1;

   p1 = (p<0.5 ? p : 1-p);
   if (p1<1e-20) return (-9999);

   y = sqrt ((fit_double)log(1/(p1*p1)));   
   z = y + ((((y*a4+a3)*y+a2)*y+a1)*y+a0) / ((((y*b4+b3)*y+b2)*y+b1)*y+b0);
   return (p<0.5 ? -z : z);
}
/*********************************************************/

int DiscreteGamma (fit_double freqK[], fit_double rK[], 
    fit_double alfa, fit_double beta, int K, int median)
{
/* COURTESY OF ZIHENG YANG
   discretization of gamma distribution with equal proportions in each 
   category
*/
   int i;
   fit_double gap05=1.0/(2.0*K), t, factor=alfa/beta*K, lnga1;

   if(K==1) 
     {
       rK[0] = 1.0;
       return 0;
     }

   if (median) {
      for (i=0; i<K; i++) rK[i]=PointGamma((i*2.0+1)*gap05, alfa, beta);
      for (i=0,t=0; i<K; i++) t+=rK[i];
      for (i=0; i<K; i++)     rK[i]*=factor/t;
   }
   else {
      lnga1=LnGamma(alfa+1);
      for (i=0; i<K-1; i++)
	 freqK[i]=PointGamma((i+1.0)/K, alfa, beta);
      for (i=0; i<K-1; i++)
	 freqK[i]=IncompleteGamma(freqK[i]*beta, alfa+1, lnga1);
      rK[0] = freqK[0]*factor;
      rK[K-1] = (1-freqK[K-2])*factor;
      for (i=1; i<K-1; i++)  rK[i] = (freqK[i]-freqK[i-1])*factor;
   }
   for (i=0; i<K; i++) freqK[i]=1.0/K;

   return (0);
}

/*********************************************************/

arbre *Read_Tree_File(FILE *fp_tree)
{
  /* Read a tree in a previously opened file */

  char *line;
  arbre *tree;
  int i;
  char c;

  line = (char *)mCalloc(T_MAX_LINE,sizeof(char));

  do
    c=fgetc(fp_tree);
  while((c != '(') && (c != EOF));

  if(c==EOF) 
      {
	Free(line);
	return NULL;
      }

  i=0;
  for(;;)
    {
        if((c==' ') || (c=='\n') || (c=='\r'))
	{
	  c=fgetc(fp_tree); 
	  if(c==EOF) break;
	  else continue;
	}
      
      line[i]=c;
      i++;
      c=fgetc(fp_tree);
      if((c==EOF) || (c==';'))  break;
    }
  
  tree = Read_Tree(line);
  Free(line);
  return tree;
}

/*********************************************************/

arbre *Read_Tree(char *s_tree)
{
  char **subs;
  int i,n_ext,n_int,n_otu;
  arbre *tree;
  int degree;


  n_otu=0;
  For(i,(int)strlen(s_tree)) if(s_tree[i] == ',') n_otu++;
  n_otu+=1;

  tree = (arbre *)Make_Tree(n_otu);
  Init_Tree(tree,n_otu);
  Make_All_Tree_Nodes(tree);
  Make_All_Tree_Edges(tree);

  tree->noeud[n_otu]->num = n_otu;
  tree->noeud[n_otu]->tax = 0;
  
  subs = Sub_Trees(s_tree,&degree);
  Clean_Multifurcation(subs,degree,3);
  if(degree == 2) 
    {
      tree->n_root = tree->noeud[2*tree->n_otu-2];
      tree->n_root->num = 2*tree->n_otu-2;
      tree->n_root->tax = 0;

      printf("\n. The input tree must be unrooted!");
      printf("\n. Err. in file %s at line %d\n",__FILE__,__LINE__);
      Exit("");
    }

  tree->has_branch_lengths = 0;
  tree->num_curr_branch_available = 0;
  n_ext = 0;
  if(!tree->n_root)
    { 
      n_int = 0;
      For(i,degree) R_rtree(s_tree,subs[i],tree->noeud[n_otu],tree,&n_int,&n_ext);
    }
  else
    {
      n_int = -1;
      For(i,degree) R_rtree(s_tree,subs[i],tree->n_root,tree,&n_int,&n_ext);
      Update_Ancestors(tree->n_root,tree->n_root->v[0]);
      Update_Ancestors(tree->n_root,tree->n_root->v[1]);
      tree->n_root->anc = NULL;
    }

  For(i,2*tree->n_otu-3) 
    {
      tree->t_edges[i]->prob_sel_regime = tree->t_edges[i]->l_prime;
    }


  For(i,NODE_DEG_MAX) Free(subs[i]);
  Free(subs);
  return tree;
}

/*********************************************************/

void Connect_One_Edge_To_Two_Nodes(node *a, node *d, edge *b, arbre *tree)
{
  int i,dir_a_d;

  dir_a_d = -1;
  For(i,3) if(a->v[i] == d) {dir_a_d = i; break;}

  if(dir_a_d == -1)
    {
      printf("\n. a->num = %d d->num = %d",a->num,d->num);
      printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Exit("");
    }


  a->b[dir_a_d] = b;
  b->num        = tree->num_curr_branch_available;
  b->left       = a;
  b->rght       = d;
  if(a->tax) {b->rght = a; b->left = d;} /* root */
  /* a tip is necessary on the right side of the edge */

  (b->left == a)?
    (Make_Edge_Dirs(b,a,d)):
    (Make_Edge_Dirs(b,d,a));

  b->l                    = a->l[b->l_r];
  if(a->tax) b->l         = a->l[b->r_l];
  if(b->l < BL_MIN)  b->l = BL_MIN;
  else if(b->l > BL_MAX) b->l = BL_MAX;
  b->l_old                = b->l;
}

/*********************************************************/

void Make_All_Edges_Lk(node *a, node *d, edge *b, arbre *tree)
{
  int i;

  Make_Edge_Lk(a,d,b,tree);
  if(d->tax) return;
  else
    {
      For(i,3)
	{
	  if(d->v[i] != a)
	    Make_All_Edges_Lk(d,d->v[i],d->b[i],tree);
	}
    }
}
	    
/*********************************************************/

void R_rtree(char *s_tree_a, char *s_tree_d, node *a, arbre *tree, int *n_int, int *n_ext)
{
  int i;
  node *d;
  int n_otu = tree->n_otu;


  if(strstr(s_tree_a," ")) Exit("\n Err : tree must not contain a ' ' character\n");

  if(s_tree_d[0] == '(')
    {
      char **subs;
      int degree;

      (*n_int)+=1;
      d      = tree->noeud[n_otu+*n_int];
      d->num = n_otu+*n_int;
      d->tax = 0;

  
      Read_Branch_Label(s_tree_d,s_tree_a,tree->t_edges[tree->num_curr_branch_available]);
      Read_Branch_Length(s_tree_d,s_tree_a,tree);

      For(i,3)
       {
	 if(!a->v[i])
	   {
	     a->v[i]=d;
	     d->l[0]=a->l[i]=tree->t_edges[tree->num_curr_branch_available]->l;
	     d->l_prime[0]=a->l_prime[i]=tree->t_edges[tree->num_curr_branch_available]->l_prime;
	     break;
	   }
       }
      d->v[0]=a;

      Connect_One_Edge_To_Two_Nodes(a,d,tree->t_edges[tree->num_curr_branch_available],tree);
      tree->num_curr_branch_available++;

      subs=Sub_Trees(s_tree_d,&degree);
      Clean_Multifurcation(subs,degree,2);
      R_rtree(s_tree_d,subs[0],d,tree,n_int,n_ext);
      R_rtree(s_tree_d,subs[1],d,tree,n_int,n_ext);
      For(i,NODE_DEG_MAX) Free(subs[i]);
      Free(subs);
    }

  else
    {
      int i;

      d      = tree->noeud[*n_ext];
      d->tax = 1;
  
      Read_Branch_Label(s_tree_d,s_tree_a,tree->t_edges[tree->num_curr_branch_available]); 
      Read_Branch_Length(s_tree_d,s_tree_a,tree);
      Read_Node_Name(d,s_tree_d,tree);
      
      For(i,3)
	{
	 if(!a->v[i])
	   {
	     a->v[i]=d;
	     d->l[0]=a->l[i]=tree->t_edges[tree->num_curr_branch_available]->l;
	     d->l_prime[0]=a->l_prime[i]=tree->t_edges[tree->num_curr_branch_available]->l_prime;
	     break;
	   }
	}
      d->v[0]=a;

      Connect_One_Edge_To_Two_Nodes(a,d,tree->t_edges[tree->num_curr_branch_available],tree);
      tree->num_curr_branch_available++;
      
      d->num=*n_ext;
      (*n_ext)+=1;
    }
}

/*********************************************************/

char **Sub_Trees(char *tree, int *degree)
{
  char **subs;
  int posbeg,posend;
  int i;

  if(tree[0] != '(') {*degree = 1; return NULL;}

  subs=(char **)mCalloc(NODE_DEG_MAX,sizeof(char *));

  For(i,NODE_DEG_MAX) subs[i]=(char *)mCalloc(strlen(tree)+1,sizeof(char));

  
  posbeg=posend=1;
  (*degree)=0;
  do
    {
      posbeg = posend;
      if(tree[posend] != '(')
	{
	  while((tree[posend] != ',' ) &&
		(tree[posend] != ':' ) &&
		(tree[posend] != '#' ) &&
		(tree[posend] != ')' )) 
	    {
	      posend++ ;
	    }
	  posend -= 1;
	}
      else posend=Next_Par(tree,posend);

      while((tree[posend+1] != ',') &&
	    (tree[posend+1] != ':') &&
	    (tree[posend+1] != '#') &&
	    (tree[posend+1] != ')')) {posend++;}


      strncpy(subs[(*degree)],tree+posbeg,posend-posbeg+1);
      strcat(subs[(*degree)],"\0");

      posend += 1;
      while((tree[posend] != ',') &&
	    (tree[posend] != ')')) {posend++;}
      posend+=1;


      (*degree)++;
      if((*degree) == NODE_DEG_MAX)
	{
	  For(i,(*degree))
	    printf("\n. Subtree %d : %s\n",i+1,subs[i]);

	  printf("\n. The degree of a node cannot be greater than %d\n",NODE_DEG_MAX);
	  Exit("\n");
	}
    }
  while(tree[posend-1] != ')');

  return subs;
}

/*********************************************************/

int Next_Par(char *s, int pos)
{
  /* go to the next parenthesis */

  int curr;
  curr=pos+1;
  
  while(*(s+curr) != ')')
    {
      if(*(s+curr) == '(') curr=Next_Par(s,curr);
      curr++;
    }

  return curr; 
}

/*********************************************************/


void Test_Tree(char *s)
{
  errno = 0;
  if(*s != ',' && *s != ')' && *s != ':' && *s != '*' && (!(fit_double)strtod(s,(char **)NULL) && errno))
    {
      printf("\n ! Tree bad format !\n");
      exit(1);
    }
  return;
}

/*********************************************************/

char *Write_Tree(arbre *tree)
{
  char *s;
  int i;

  s=(char *)mCalloc(T_MAX_LINE,sizeof(char));

  s[0]='(';
  
  #ifdef PHYML 
  tree->n_root = NULL;
  tree->e_root = NULL;
  #endif
  

  if(!tree->n_root)
    {
      i = 0;
      while((!tree->noeud[tree->n_otu+i]->v[0]) ||
	    (!tree->noeud[tree->n_otu+i]->v[1]) ||
	    (!tree->noeud[tree->n_otu+i]->v[2])) i++;
      
      R_wtree(tree->noeud[tree->n_otu+i],tree->noeud[tree->n_otu+i]->v[0],s,tree);
      R_wtree(tree->noeud[tree->n_otu+i],tree->noeud[tree->n_otu+i]->v[1],s,tree);
      R_wtree(tree->noeud[tree->n_otu+i],tree->noeud[tree->n_otu+i]->v[2],s,tree);
    }
  else
    {
      R_wtree(tree->n_root,tree->n_root->v[0],s,tree);
      R_wtree(tree->n_root,tree->n_root->v[1],s,tree);
    }

  s[(int)strlen(s)-1]=')';
  s[(int)strlen(s)]=';';

  return s;
}

/*********************************************************/

void R_wtree(node *a, node *d, char *s_tree, arbre *tree)
{
  int i,p;

  p = -1;
  if(d->tax)
    {
      sprintf(s_tree+(int)strlen(s_tree),"%s",d->n_label);
      sprintf(s_tree+(int)strlen(s_tree),"%s",d->name);

      if((d->b) && (d->b[0]) && (d->b[0]->l != -1))
	{
	  strcat(s_tree,":");

	  if(a != tree->n_root)
	    {
	      sprintf(s_tree+(int)strlen(s_tree),"%f",d->b[0]->l);
	      strcat(s_tree,":");
	      if(d->b[0]->labels) sprintf(s_tree+(int)strlen(s_tree),"%s",d->b[0]->labels[0]);
	    }
	  else
	    {
	      if(tree->n_root->v[0] == d)
		{
		  sprintf(s_tree+(int)strlen(s_tree),"%f",tree->n_root->l[0]);
		}
	      else
		{
		  sprintf(s_tree+(int)strlen(s_tree),"%f",tree->n_root->l[1]);
		}
	    }
	}
      sprintf(s_tree+(int)strlen(s_tree),",");
   }
  else
    {
      s_tree[(int)strlen(s_tree)]='(';

      if(tree->n_root)
	{
	  For(i,3)
	    {
	      if((d->v[i] != a) && (d->b[i] != tree->e_root))
		R_wtree(d,d->v[i],s_tree,tree);
	      else p=i;
	    }
	}
      else
	{
	  For(i,3)
	    {
	      if(d->v[i] != a)
		R_wtree(d,d->v[i],s_tree,tree);
	      else p=i;
	    }
	}

      s_tree[(int)strlen(s_tree)-1]=')';
      if((d->b) && (d->b[0]->l != -1))
	{
	  sprintf(s_tree+(int)strlen(s_tree),"%s",d->n_label);

	  strcat(s_tree,":");

	  if(a != tree->n_root)
	    {
	      sprintf(s_tree+(int)strlen(s_tree),"%f",d->b[p]->l);
	      strcat(s_tree,":");
	      if(d->b[0]->labels) sprintf(s_tree+(int)strlen(s_tree),"%s",d->b[p]->labels[0]);
	    }
	  else
	    {
	      if(tree->n_root->v[0] == d)
		{
		  sprintf(s_tree+(int)strlen(s_tree),"%f",tree->n_root->l[0]);
		}
	      else
		{
		  sprintf(s_tree+(int)strlen(s_tree),"%f",tree->n_root->l[1]);
		}
	    }



	}
      strcat(s_tree,",");
    }
}

/*********************************************************/

void Init_Tree(arbre *tree, int n_otu)
{
  tree->n_otu              = n_otu;
  tree->noeud              = (node **)mCalloc(2*tree->n_otu-2,sizeof(node *));
  tree->t_edges            = (edge **)mCalloc(2*tree->n_otu-3,sizeof(edge *));
  
  tree->e_root             = NULL;
  tree->n_root             = NULL;
  tree->has_bip            = 0;

  tree->tot_loglk          = UNLIKELY;
  tree->n_swap             = 0;

  tree->n_pattern          = -1;

  tree->render_tax_names   = 0;
  tree->render_edge_colors = 0;
  tree->tree_number        = 0;
  tree->n_root_pos         = .5;
}

/*********************************************************/

edge *Make_Edge_Light(node *a, node *d, int num)
{
  edge *b;

  b = (edge *)mCalloc(1,sizeof(edge));

  Init_Edge_Light(b,num);

  if(a && b)
    {
      b->left = a;  b->rght = d;
      if(a->tax) {b->rght = a; b->left = d;} /* root */
      /* a tip is necessary on the right side of the edge */

      (b->left == a)?
	(Make_Edge_Dirs(b,a,d)):
	(Make_Edge_Dirs(b,d,a));

      b->l                    = a->l[b->l_r];
      b->prob_sel_regime      = a->l_prime[b->l_r];

      if(a->tax) {b->l        = a->l[b->r_l]; b->prob_sel_regime = a->l_prime[b->r_l];}
      if(b->l < BL_MIN)  b->l = BL_MIN;
      else if(b->l > BL_MAX) b->l = BL_MAX;
      b->l_old                = b->l;
    }
  else
    {
      b->left = NULL;
      b->rght = NULL;
    }

  return b;

}

/*********************************************************/

void Make_NNI(nni *this1)
{
  /* deprecated */
  this1->v1 = 
    this1->v2 = 
    this1->v3 = 
    this1->v4 = 
    this1->v5 = 
    this1->v6 = NULL;
  this1->best_conf = -1;
  this1->b_fcus = NULL;
  this1->bl_info = (fit_double *)mCalloc(3,sizeof(fit_double));
  this1->bl_info_init = -1.;
}

/*********************************************************/

void Init_Edge_Light(edge *b, int num)
{
  b->num               = num;
  b->substitution_type = 0;
  b->l_ml              = 1.E+10;

  b->p_lk_left         = NULL;
  b->p_lk_rght         = NULL;
  b->Pij_rr            = NULL;
  b->bip_score         = 0;

}

/*********************************************************/

void Make_Edge_Dirs(edge *b, node *a, node *d)
{
  /* Set the directions (0, 1 or 2) around the nodes at the two extremities of a branch */

  int i;

  b->l_r = b->r_l = -1;
  For(i,3)
    {
      if((a->v[i]) && (a->v[i] == d)) 
	{
	  b->l_r  = i;
	  a->b[i] = b;
	}
      if((d->v[i]) && (d->v[i] == a)) 
	{
	  b->r_l  = i;
	  d->b[i] = b;
	}
    }

  if(a->tax) {b->r_l = 0; For(i,3) if(d->v[i]==a) {b->l_r = i; break;}}


  b->l_v1 = b->l_v2 = b->r_v1 = b->r_v2 = -1;
  For(i,3)
    {
      if(b->left->v[i] != b->rght)
	{
	  if(b->l_v1 < 0) b->l_v1 = i;
	  else            b->l_v2 = i;
	}
     
      if(b->rght->v[i] != b->left)
	{
	  if(b->r_v1 < 0) b->r_v1 = i;
	  else            b->r_v2 = i;
	}
    }
}

/*********************************************************/

void Make_Edge_Lk(node *a, node *d, edge *b, arbre *tree)
{
  /* Allocate the matrix of change probabilities in the edge connecting nodes a & d */

  int i,j,k;
  int n_catq, n_catg;

  b->diff_lk   = 0.0;

  n_catq = tree->mod->n_catq;
  n_catg = tree->mod->n_catg;

  b->qmat_struct = tree->mod->qmat_struct[b->num];

  if(!b->Pij_rr)
    {
      b->best_conf = 1;

      b->Pij_rr    = (fit_double ****)mCalloc(n_catq,sizeof(fit_double ***));
      
      For(i,n_catq)
	{
	  b->Pij_rr[i] = (fit_double ***)mCalloc(n_catg,sizeof(fit_double **));
	  
	  For(j,n_catg)
	    {
	      b->Pij_rr[i][j]   = (fit_double **)mCalloc(tree->mod->ns,sizeof(fit_double *));
	      	      
	      For(k,tree->mod->ns)
		{
		  b->Pij_rr[i][j][k]   = (fit_double *)mCalloc(tree->mod->ns,sizeof(fit_double ));
		}
	    }
	}

  
      b->p_lk_left   = NULL;
      b->p_lk_rght  = NULL;
      
      b->sum_scale_f_rght = NULL;
      b->sum_scale_f_left  = NULL;

      b->get_p_lk_left = 0;
      b->get_p_lk_rght = 0;
      b->ud_p_lk_left  = 0;
      b->ud_p_lk_rght  = 0;
    }
}

/*********************************************************/

node *Make_Node_Light()
{
  node *n;
  n             = (node *)mCalloc(1,sizeof(node));
  n->v          = (node **)mCalloc(3,sizeof(node *));
  n->l          = (fit_double *)mCalloc(3,sizeof(fit_double));
  n->l_prime    = (fit_double *)mCalloc(3,sizeof(fit_double));
  n->b          = (edge **)mCalloc(3,sizeof(edge *));
  n->name       = (char *)mCalloc(T_MAX_NAME,sizeof(char));
  n->tpos_min   = (fit_double *)mCalloc(3,sizeof(fit_double));
  n->n_label    = (char *)mCalloc(T_MAX_NAME,sizeof(char));
  n->post_w_ave = (fit_double *)mCalloc(3,sizeof(fit_double));
  Init_Node_Light(n);  
  return n;
}

/*********************************************************/

void Init_Node_Light(node *n)
{
  int i;

  n->seq            = NULL;
  n->ni = n->agglo  = 0;
  n->tpos           = -1.;
  n->tpos_lim_above = -1.;
  n->tpos_lim_below = -1.;
  n->check_branch   = -1;

  For(i,3)
    {
      n->v[i]=NULL;
      n->b[i]=NULL;
      n->l[i]=-1;
      n->l_prime[i]=-1;
      n->tpos_min[i]=-1.;
    }
}

/*********************************************************/

void Make_Node_Lk(arbre *tree, node *n)
{
  n->seq = (char *)mCalloc(tree->n_pattern*tree->mod->stepsize+1,sizeof(char));
  return;
}

/*********************************************************/

seq **Get_Seq(option *input, int rw)
{
  /* Read sequences */

  seq **data;
  int i,j;
  char **buff;
  int n_unkn,n_removed,pos;
  int *remove;


/*   rewind(fp_seq); */

  if(input->interleaved) data = Read_Seq_Interleaved(input->fp_seq,&(input->mod->n_otu));
  else                   data = Read_Seq_Sequential(input->fp_seq,&(input->mod->n_otu));

  if(data)
    {
      buff = (char **)mCalloc(input->mod->n_otu,sizeof(char *));
      For(i,input->mod->n_otu) buff[i] = (char *)mCalloc(data[0]->len,sizeof(char));
      remove = (int *)mCalloc(data[0]->len,sizeof(int));
  
      n_removed = 0;
      For(i,data[0]->len)
	{
	  For(j,input->mod->n_otu)
	    {
	      if((data[j]->state[i] == '?') ||
		 (data[j]->state[i] == '-') ||
		 (data[j]->state[i] == 'X')) data[j]->state[i] = 'X';
	      if((input->mod->datatype == NT) && (data[j]->state[i] == 'N')) data[j]->state[i] = 'X';
	    }
	  
	  n_unkn = 0;
	  For(j,input->mod->n_otu) if(data[j]->state[i] == 'X') n_unkn++; 
	  
	  if(n_unkn == input->mod->n_otu)
	    {
	      remove[i] = 1;
	      n_removed++;
	    }
	  
	  For(j,input->mod->n_otu) buff[j][i] = data[j]->state[i];
	}
      
      if(n_removed > 0) 
	{
	  if(input->mod->datatype == NT)
	    printf("\n. %d sites are made from completely undetermined states ('X', '-', '?' or 'N')...\n",n_removed);
	  else
	    printf("\n. %d sites are made from completely undetermined states ('X', '-', '?')...\n",n_removed);
	}

      pos = 0;
      For(i,data[0]->len)
	{
/* 	  if(!remove[i]) */
/* 	    { */
	      For(j,input->mod->n_otu) data[j]->state[pos] = buff[j][i];
	      pos++;
/* 	    } */
	}

      For(i,input->mod->n_otu) data[i]->len = pos;
      For(i,input->mod->n_otu) Free(buff[i]);
      Free(buff);
      Free(remove);
    }
  return data;
}

/*********************************************************/
   
seq **Read_Seq_Sequential(FILE *in, int *n_otu)
{
  /* Read sequences in 'sequential' format */

  int i;
  char *line;
  int len,readok;
  seq **data;
  char c;
  char *format;

  format = (char *)mCalloc(T_MAX_NAME, sizeof(char));
  line = (char *)mCalloc(T_MAX_LINE,sizeof(char));

  readok = len = 0;
  do
    {
      if(fscanf(in,"%s",line) == EOF)
	{
	  Free(line); 
	  Free(format);
	  return NULL;
	}
      else
	{
	  if(strcmp(line,"\n") && strcmp(line,"\n") && strcmp(line,"\t"))
	    {
	      *n_otu = atoi(line);
	      if(*n_otu <= 0) Exit("\n. Problem with sequence format\n");
	      data = (seq **)mCalloc(*n_otu,sizeof(seq *));
	      if(!fscanf(in,"%s",line))
		{
		  printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
		  Exit("");
		}

	      len = atoi(line);
	      if(len <= 0) Exit("\n. Problem with sequence format\n");
	      else readok = 1;
	    }
	}
    }while(!readok);
  
  while(((c=fgetc(in))!='\n')); // Go to the end of the line
  while(((c=fgetc(in)) == ' ') || (c == '\r') || (c == '\t'));
  fseek(in,-1*sizeof(char),SEEK_CUR);

  For(i,*n_otu)
    {
      data[i] = (seq *)mCalloc(1,sizeof(seq));
      data[i]->len = 0;
      data[i]->name = (char *)mCalloc(T_MAX_NAME,sizeof(char));
      data[i]->state = (char *)mCalloc(T_MAX_SEQ,sizeof(char));
      sprintf(format, "%%%ds", T_MAX_NAME);
/*       sprintf(format, "%%%ds", 10); */
      if(!fscanf(in,format,data[i]->name))
	{
	  printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	  Exit("");
	}
/*       while(data[i]->len < len) */
      Read_One_Line_Seq(&data,i,in);
	
      if(data[i]->len != len) 
	{
	  printf("\n. Err: Problem with species %s's sequence (check the format)\n",data[i]->name);
	  Exit("");
	}
    }

  /*   fgets(line,T_MAX_LINE,in);  */
  /* inter data sets */
  
  Free(format);
  Free(line);
  return data;
}

/*********************************************************/

seq **Read_Seq_Interleaved(FILE *in, int *n_otu)
{
  /* Read sequences in 'interleaved' format */

  int i,end,num_block;
  char *line;
  int len,readok;
  seq **data;
  char c;
  char *format;

  line = (char *)mCalloc(T_MAX_LINE,sizeof(char));
  format = (char *)mCalloc(T_MAX_NAME, sizeof(char));

  readok = len = 0;
  do
    {
      if(fscanf(in,"%s",line) == EOF)
	{
	  Free(format);
	  Free(line); 
	  return NULL;
	}
      else
	{
	  if(strcmp(line,"\n") && strcmp(line,"\r") && strcmp(line,"\t"))
	    {
	      *n_otu = atoi(line);
	      if(*n_otu <= 0) Exit("\n. Problem with sequence format\n");
	      data = (seq **)mCalloc(*n_otu,sizeof(seq *));
	      if(!fscanf(in,"%s",line))
		{
		  printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
		  Exit("");
		}

	      len = atoi(line);
	      if(len <= 0) Exit("\n. Problem with sequence format\n");
	      else readok = 1;
 
 	    }
	}
    }while(!readok);

  while(((c=fgetc(in))!='\n')); // Go to the end of the line
  while(((c=fgetc(in)) == ' ') || (c == '\r') || (c == '\t'));
  fseek(in,-1*sizeof(char),SEEK_CUR);

  end = 0;
  For(i,*n_otu)
    {
      data[i] = (seq *)mCalloc(1,sizeof(seq));      
      data[i]->len = 0;
      data[i]->name = (char *)mCalloc(T_MAX_NAME,sizeof(char));
      data[i]->state = (char *)mCalloc(T_MAX_SEQ,sizeof(char));
      sprintf(format, "%%%ds", T_MAX_NAME);
      /* sprintf(format, "%%%ds", 10); */
      if(!fscanf(in,format,data[i]->name))
	{
	  printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	  Exit("");
	}
      if(!Read_One_Line_Seq(&data,i,in)) 
	{
	  end = 1;
	  if(i != *n_otu) 
	    {
	      printf("\n. Err: Problem with species %s's sequence\n",data[i]->name);
	      Exit("");
	    }
	  break;
	}
    }

  if(data[0]->len == len) end = 1;

  if(!end)
    {
      end = 0;

      num_block = 1;
      do
	{
	  num_block++;
	  
	  /* interblock */
	  if(!fgets(line,T_MAX_LINE,in)) break;
	  
	  if(line[0] != 13 && line[0] != 10) 
	    {
                printf("\n. One or more missing sequences in block %d\n",num_block-1);
                Exit("");
	    }

	  For(i,*n_otu)
	    if(data[i]->len != len)
	      break;

	  if(i == *n_otu) break;

	  
	  For(i,*n_otu)
	    {
	      if(data[i]->len > len) 
		{
		  printf("\n. Err: Problem with species %s's sequence\n",data[i]->name);
		  Exit("");
		}
	      else if(!Read_One_Line_Seq(&data,i,in)) 
		{
		  end = 1;
		  if(i != *n_otu) 
		    {
		      printf("\n. Err: Problem with species %s's sequence\n",data[i]->name);
		      Exit("");
		    }
		  break;
		}
	    }
	}while(!end);
    }

  For(i,*n_otu)
    if(data[i]->len != len)
      Exit("\n. Problem with sequence length...\n");

  Free(format);
  Free(line);
  return data;
}
 
/*********************************************************/

int Read_One_Line_Seq(seq ***data, int num_otu, FILE *in)
{
  char c = ' ';
  int nchar = 0;

  while(1)
    {
/*       if((c == EOF) || (c == '\n') || (c == '\r')) break; */

      if((c == 13) || (c == 10))
        {
          /* 	  PhyML_Printf("[%d %d]\n",c,nchar); fflush(NULL); */
          if(!nchar)
            {
              c=(char)fgetc(in);
              continue;
            }
          else
            {
              /* 	      PhyML_Printf("break\n");  */
              break;
            }
        }
      else if(c == EOF)
        {
          /* 	  PhyML_Printf("EOL\n"); */
          break;
        }
      else if((c == ' ') || (c == '\t') || (c == 32))
        {
          /* 	  PhyML_Printf("[%d]",c); */
          c=(char)fgetc(in);
          continue;
        }
      
      nchar++;
      Uppercase(&c);
 
      if (strchr("ABCDEFGHIKLMNOPQRSTUVWXY?-.", c) == NULL)
	{
	  printf("\n. Err: bad symbol: \"%c\" at position %d of species %s\n",
		 c,(*data)[num_otu]->len,(*data)[num_otu]->name);
	  Exit("");
	}

     if(c == '.')
        {
          c = (*data)[0]->state[(*data)[num_otu]->len];
          if(!num_otu)
            Exit("\n== Err: Symbol \".\" should not appear in the first sequence\n");
        }
      (*data)[num_otu]->state[(*data)[num_otu]->len]=c;
      (*data)[num_otu]->len++;
      c = (char)fgetc(in);
      /* PhyML_Printf("[%c %d]",c,c); */
      if(c == ';') break;
    }
  
  /* printf("\n. Exit nchar: %d [%d]\n",nchar,c==EOF); */
  if(c == EOF) return 0;
  else return 1;
}
   
/*********************************************************/

void Uppercase(char *ch)
{
  /* convert ch to upper case -- either ASCII or EBCDIC */
   *ch = isupper((int)*ch) ? *ch : toupper((int)*ch);
}

/*********************************************************/

allseq *Compact_Seq(seq **data, option *input)
{
  /* Compress sequences. For example, two sites that show the same pattern 
     are considered as one wite with weight 1 */
  
  allseq *alldata;
  int i,j,k,diff,site;
  int n_patt,which_patt,n_invar;
  char **sp_names;
  int n_otu;
  int len;


  n_otu = input->mod->n_otu;
  
  sp_names = (char **)mCalloc(n_otu,sizeof(char *));
  For(i,n_otu) 
    {
      sp_names[i] = (char *)mCalloc(T_MAX_NAME,sizeof(char));
      strcpy(sp_names[i],data[i]->name);
    }

  /*   alldata = Make_Seq(n_otu, */
  /* 		     (input->mod->seq_len == -1)? */
  /* 		     (data[0]->len): */
  /* 		     (input->mod->seq_len),sp_names); */
  
  if(data[0]->len > input->user_len)
    len = data[0]->len;
  else
    len = input->user_len;
  
  alldata = Make_Cseq(n_otu,len,data[0]->len,sp_names);

  
  For(i,n_otu) Free(sp_names[i]);
  Free(sp_names);
  
  n_patt = which_patt = 0;
  diff = -1;
  
  
  if(len%input->mod->stepsize) 
    {
      printf("\n. Sequence length (%d) is not a multiple of %d\n",len,input->mod->stepsize);
      Exit("");
    }
  
  
  Fors(site,data[0]->len,input->mod->stepsize) 
    { 
      Fors(k,n_patt,input->mod->stepsize)
	{
	  For(j,n_otu) /* Compare two patterns */
	    {
	      diff = 0;
	      if(!Compare_Two_States(alldata->c_seq[j]->state+k,
				     data[j]->state+site,
				     input->mod->stepsize)) diff = 1;
	      if(diff) break;
	    }
	  
	  if(j == n_otu)
	    {
	      which_patt = k;
	      break;
	    }
	}
      
      
#if defined(EVOLVE) || defined(POSSELSITEID)  || defined(FITMODEL) 
      /* sites are not compressed */
      k = n_patt;
#endif
      
      if(k == n_patt) /* k == n_patt iif current site is a new one */
	{
	  For(j,n_otu) 
	    Copy_One_State(data[j]->state+site,
			   alldata->c_seq[j]->state+n_patt,
			   input->mod->stepsize);
	  
          
	  for(j=0;j<n_otu;j++) 
	    {
	      if(!(Are_Compatible(alldata->c_seq[j]->state+n_patt,
				  alldata->c_seq[0]->state+n_patt,
				  input->mod->stepsize,
				  input->mod->datatype)))
		break; /* this site not invariable or this site does not 
			  contain only ambiguous character states */
	    }
	  
	  if(j==n_otu) /* invariable site */
	    {
	      For(j,n_otu)
		{
		  alldata->invar[n_patt] = Assign_State(alldata->c_seq[j]->state+n_patt,
							input->mod->datatype,
							input->mod->stepsize,
							input->mod->c_code);
		  break;
		}
	    }
	  else    
	    alldata->invar[n_patt] = -1; /* note : completely ambiguous site are not considered as
					    invariant sites */
	  
          
	  alldata->pospatt[n_patt] = (int *)mRealloc(alldata->pospatt[n_patt],(int)alldata->wght[n_patt]+1,sizeof(int));
	  alldata->pospatt[n_patt][(int)alldata->wght[n_patt]] = site;
	  
	  For(j,input->mod->stepsize) alldata->wght[n_patt+j] += 1;
          
	  n_patt+=input->mod->stepsize;
	}
      else /* current site has already been encountered before */
	{
	  alldata->pospatt[which_patt] = (int *)mRealloc(alldata->pospatt[which_patt],
							 (int)alldata->wght[which_patt]+1,sizeof(int));
	  alldata->pospatt[which_patt][(int)alldata->wght[which_patt]] = site;
	  
	  For(j,input->mod->stepsize) alldata->wght[which_patt+j] += 1;
	}
    }
  
  
  alldata->init_len = len;
  alldata->crunch_len = n_patt;
  For(i,n_otu) alldata->c_seq[i]->len = n_patt;
  
  n_invar=0;
  For(i,alldata->crunch_len) if(alldata->invar[i]>-1) n_invar+=(int)alldata->wght[i];
  
  if(input->mod->datatype == NT)
    Get_Base_Freqs(alldata,0,1);
  else
    Get_AA_Freqs(alldata);
  
  return alldata;
}

/*********************************************************/

allseq *Compact_CSeq(allseq *data, model *mod)
{
  /* This is the same function as 'Compact_Seq', except that it takes
     already compressed sequences as input. This function might be usefull
     when you want to compress subset of compressed sequences */


  allseq *alldata;
  int i,j,k,site;
  int n_patt,which_patt;
  int n_otu;

  n_otu = data->n_otu;

  alldata = (allseq *)mCalloc(1,sizeof(allseq));
  alldata->n_otu=n_otu;
  alldata->c_seq = (seq **)mCalloc(n_otu,sizeof(seq *));
  alldata->wght = (int *)mCalloc(data->crunch_len,sizeof(int));
  alldata->b_frq = (fit_double *)mCalloc(mod->ns,sizeof(fit_double));
  alldata->ambigu = (short int *)mCalloc(data->crunch_len,sizeof(short int));
  alldata->pospatt = (int **)mCalloc(data->init_len,sizeof(int *));
  alldata->invar = (short int *)mCalloc(data->crunch_len,sizeof(short int));

  alldata->crunch_len = alldata->init_len = -1;
  For(j,n_otu)
    {
      alldata->c_seq[j] = (seq *)mCalloc(1,sizeof(seq));
      alldata->c_seq[j]->name = (char *)mCalloc(T_MAX_NAME,sizeof(char));
      strcpy(alldata->c_seq[j]->name,data->c_seq[j]->name);
      alldata->c_seq[j]->state = (char *)mCalloc(data->crunch_len,sizeof(char));
      alldata->c_seq[j]->state[0] = data->c_seq[j]->state[0];
    }
  
  n_patt = which_patt =  0;


  Fors(site,data->crunch_len,mod->stepsize) 
    {
      Fors(k,n_patt,mod->stepsize)
	{
	  For(j,n_otu)
	    {
	      if(!Compare_Two_States(alldata->c_seq[j]->state+k,
				     data->c_seq[j]->state+site,
				     mod->stepsize))
		break;
	    }

	  if(j == n_otu)
	    {
	      which_patt = k;
	      break;
	    }
	}
      
      
      if(k == n_patt)
	{
	  For(j,n_otu) Copy_One_State(data->c_seq[j]->state+site,
				      alldata->c_seq[j]->state+n_patt,
				      mod->stepsize);
	  
	  for(j=1;j<n_otu;j++) 
	    if(!Compare_Two_States(alldata->c_seq[j]->state+n_patt,
				   alldata->c_seq[j-1]->state+n_patt,
				   mod->stepsize)) break;
	  
	  if(j==n_otu) alldata->invar[n_patt] = 1;
	  alldata->wght[n_patt] += data->wght[site];
	  n_patt+=mod->stepsize;
	}
      else alldata->wght[which_patt] += data->wght[site];


    }

  alldata->init_len = data->crunch_len;
  alldata->crunch_len = n_patt;
  For(i,n_otu) alldata->c_seq[i]->len = n_patt;


  return alldata;
}

/*********************************************************/

void Get_Codon_Freqs(allseq *data, model *mod)
{
  /* Compute codon frequencies using base frequencies at the three coding 
     positions (F3X4 in CODEML) */

  int i;

  if(!mod->codon_freq) 
    {
      For(i,mod->c_code->n_sense_c) 
	data->b_frq[i] = 1./mod->c_code->n_sense_c;
    }
  else
    {
      int j,k,num_codon;
      fit_double **base_freqs_pos;
      fit_double sum;

      base_freqs_pos = (fit_double **)mCalloc(4,sizeof(fit_double *));
      For(i,4) base_freqs_pos[i] = (fit_double *)mCalloc(3,sizeof(fit_double));

      For(i,3) 
	{
	  Get_Base_Freqs(data,i,3);
	  For(j,4) 
	    {
	      base_freqs_pos[j][i] = data->b_frq[j];
	    }
	}
      sum = .0;
      For(i,4)
	{
	  For(j,4)
	    {
	      For(k,4)
		{
		  num_codon = i*16+j*4+k;
		  if(mod->c_code->from_64_2_61[num_codon] > -1) 
		    {
		      data->b_frq[mod->c_code->from_64_2_61[num_codon]] = 
			base_freqs_pos[i][0] *
			base_freqs_pos[j][1] *
			base_freqs_pos[k][2];
		      sum += data->b_frq[mod->c_code->from_64_2_61[num_codon]];
		    }
		}
	    }
	}

      For(i,mod->c_code->n_sense_c) data->b_frq[i] /= sum;

      For(i,4) Free(base_freqs_pos[i]);
      Free(base_freqs_pos);
    }
}

/*********************************************************/

void Get_Base_Freqs(allseq *data, int position, int step)
{
  /* Compute observed base frequencies */

  int otu,site,k;
  fit_double A,C,G,T;
  fit_double fA,fC,fG,fT;
  fit_double w;

  fA = fC = fG = fT = .25;

  For(k,8) /* this loop allows to takes into account ambiguous characters in the computation
	      of base frequencies */
    {
      A = C = G = T = .0;
      For(otu,data->n_otu)
	{
	  for(site=position;site < data->crunch_len;site+=step)
	    {
	      w = data->wght[site];
	      if(w)
		{
		  switch(data->c_seq[otu]->state[site]){
		  case 'A' : A+=w;
		    break;
		  case 'C' : C+=w;
		    break;
		  case 'G' : G+=w;
		    break;
		  case 'T' : T+=w;
		    break;
		  case 'U' : T+=w;
		    break;
		  case 'M' : C+=w*fC/(fC+fA); A+=w*fA/(fA+fC);
		    break;
		  case 'R' : G+=w*fG/(fA+fG); A+=w*fA/(fA+fG);
		    break;
		  case 'W' : T+=w*fT/(fA+fT); A+=w*fA/(fA+fT);
		    break;
		  case 'S' : C+=w*fC/(fC+fG); G+=w*fG/(fC+fG);
		    break;
		  case 'Y' : C+=w*fC/(fC+fT); T+=w*fT/(fT+fC);
		    break;
		  case 'K' : G+=w*fG/(fG+fT); T+=w*fT/(fT+fG);
		    break;
		  case 'B' : C+=w*fC/(fC+fG+fT); G+=w*fG/(fC+fG+fT); T+=w*fT/(fC+fG+fT);
		    break;
		  case 'D' : A+=w*fA/(fA+fG+fT); G+=w*fG/(fA+fG+fT); T+=w*fT/(fA+fG+fT);
		    break;
		  case 'H' : A+=w*fA/(fA+fC+fT); C+=w*fC/(fA+fC+fT); T+=w*fT/(fA+fC+fT);
		    break;
		  case 'V' : A+=w*fA/(fA+fC+fG); C+=w*fC/(fA+fC+fG); G+=w*fG/(fA+fC+fG);
		    break;
		  case 'N' : case 'X' : case '?' : case 'O' : case '-' : 
		    A+=w*fA; C+=w*fC; G+=w*fG; T+=w*fT; break;
		  default : break;
		  }
		}
	    }  
	}
      fA = A/(A+C+G+T);
      fC = C/(A+C+G+T);
      fG = G/(A+C+G+T);
      fT = T/(A+C+G+T);
    }

  data->b_frq[0] = fA;
  data->b_frq[1] = fC;
  data->b_frq[2] = fG;
  data->b_frq[3] = fT;

}

/*********************************************************/

void Get_AA_Freqs(allseq *data)
{
  /* Compute amino acid observed frequencies */ 
  int i,j,k;
  fit_double A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y;
  fit_double fA,fC,fD,fE,fF,fG,fH,fI,fK,fL,fM,fN,fP,fQ,fR,fS,fT,fV,fW,fY;
  fit_double w;
  fit_double sum;


  fA = fC = fD = fE = fF = fG = fH = fI = fK = fL = 
  fM = fN = fP = fQ = fR = fS = fT = fV = fW = fY = 1./20.;

  For(k,8) /* this loop allows to takes into account ambiguous characters in the computation
	      of aa frequencies */
    {

      A = C = D = E = F = G = H = I = K = L = 
      M = N = P = Q = R = S = T = V = W = Y = .0;
      
      For(i,data->n_otu)
	{
	  For(j,data->crunch_len)
	    {
	      w = data->wght[j];
	      if(w)
		{
		  switch(data->c_seq[i]->state[j]){
		  case 'A' : A+=w;		break;
		  case 'C' : C+=w;		break;
		  case 'D' : D+=w;		break;
		  case 'E' : E+=w;		break;
		  case 'F' : F+=w;		break;
		  case 'G' : G+=w;		break;
		  case 'H' : H+=w;		break;
		  case 'I' : I+=w;		break;
		  case 'K' : K+=w;		break;
		  case 'L' : L+=w;		break;
		  case 'M' : M+=w;		break;
		  case 'N' : N+=w;		break;
		  case 'P' : P+=w;		break;
		  case 'Q' : Q+=w;		break;
		  case 'R' : R+=w;		break;
		  case 'S' : S+=w;		break;
		  case 'T' : T+=w;		break;
		  case 'V' : V+=w;		break;
		  case 'W' : W+=w;		break;
		  case 'Y' : Y+=w;		break;
		  case 'X' : case '?' : case 'O' : case '-' : 
		    A+=w*fA;
		    C+=w*fC; 
		    D+=w*fD; 
		    E+=w*fE; 
		    F+=w*fF; 
		    G+=w*fG; 
		    H+=w*fH; 
		    I+=w*fI; 
		    K+=w*fK; 
		    L+=w*fL; 
		    M+=w*fM; 
		    N+=w*fN; 
		    P+=w*fP; 
		    Q+=w*fQ; 
		    R+=w*fR; 
		    S+=w*fS; 
		    T+=w*fT; 
		    V+=w*fV; 
		    W+=w*fW; 
		    Y+=w*fY; 
		    break;
		  default : break;
		  }
		}
	    }  
	}
      sum = (A+C+D+E+F+G+H+I+K+L+M+N+P+Q+R+S+T+V+W+Y);
      fA = A/sum;      fC = C/sum;      fD = D/sum;      fE = E/sum;
      fF = F/sum;      fG = G/sum;      fH = H/sum;      fI = I/sum;
      fK = K/sum;      fL = L/sum;      fM = M/sum;      fN = N/sum;
      fP = P/sum;      fQ = Q/sum;      fR = R/sum;      fS = S/sum;
      fT = T/sum;      fV = V/sum;      fW = W/sum;      fY = Y/sum;
    }

  data->b_frq[0]  = fA;  data->b_frq[1]  = fR;  data->b_frq[2]  = fN;  data->b_frq[3]  = fD;
  data->b_frq[4]  = fC;  data->b_frq[5]  = fQ;  data->b_frq[6]  = fE;  data->b_frq[7]  = fG;
  data->b_frq[8]  = fH;  data->b_frq[9]  = fI;  data->b_frq[10] = fL;  data->b_frq[11] = fK;
  data->b_frq[12] = fM;  data->b_frq[13] = fF;  data->b_frq[14] = fP;  data->b_frq[15] = fS;
  data->b_frq[16] = fT;  data->b_frq[17] = fW;  data->b_frq[18] = fY;  data->b_frq[19] = fV;
}

/*********************************************************/

void Init_Tree_Edges(node *a, node *d, arbre *tree, int *cur)
{
  int i,dir_a_d;

  dir_a_d = -1;
  For(i,3) if(a->v[i] == d) {dir_a_d = i; break;}

  
  tree->t_edges[*cur] = a->b[dir_a_d];
  tree->t_edges[*cur]->num = *cur;
  *cur = *cur + 1;

  if(d->tax) return;
  else
    {
      For(i,3)
	{
	  if(d->v[i] != a)
	    Init_Tree_Edges(d,d->v[i],tree,cur);
	}
    }
}

/*********************************************************/

void Exit(char *message)
{
  fflush(NULL);
  fprintf(stderr,"%s",message);
  exit(1);
}

/*********************************************************/

void *mCalloc(int nb, size_t size)
{
  void *allocated;


  if((allocated = calloc((size_t)nb,(size_t)size)) != NULL)
    {
      return allocated;
    }
  else 
    Exit("\n. Err: low memory\n");

  return NULL;
}

/*********************************************************/

void *mRealloc(void *p,int nb, size_t size)
{

  if((p = realloc(p,(size_t)nb*size)) != NULL)
	return p;
  else
    Exit("\n. Err: low memory\n");
  
  return NULL;
}

/*********************************************************/


/*********************************************************/

int Sort_Float_Decrease(const void *a, const void *b)
{
  if((*(fit_double *)(a)) >= (*(fit_double *)(b))) return -1;
  else return 1;
}

/*********************************************************/

void Print_Site(allseq *alldata, int num, int n_otu, char *sep, int stepsize)
{
  /* Print a site of a compressed sequence data set */

  int i,j;
  For(i,n_otu) 
    {
      printf("%20s   ",alldata->c_seq[i]->name);
      For(j,stepsize)
	printf("%c",alldata->c_seq[i]->state[num+j]);
      printf("%s",sep);
    }
  fprintf(stderr,"%s",sep);
}

/*********************************************************/

void Print_Seq(seq **data, int n_otu, FILE *fp)
{
  /* Print uncompressed sequences */
  
  int i,j;
  
  fprintf(fp,"%d\t%d\n",n_otu,data[0]->len-3);
  For(i,n_otu)
    {
      For(j,23)
	{
	  if(j<(int)strlen(data[i]->name))
              putc(data[i]->name[j],fp);
	  else putc(' ',fp);
	}
/*       For(j,data[i]->len) */
      for(j=2;j<data[i]->len-1;j++)
	{
            fprintf(fp,"%c",data[i]->state[j]);
	}
      fprintf(fp,"\n");
    }
}

/*********************************************************/

void Print_CSeq(FILE *fp, allseq *alldata)
{
  /* Print compressed sequences */

  int i,j,k;
  int n_otu;
  
  n_otu = alldata->n_otu;
  fprintf(fp,"%d\t%d\n",n_otu,alldata->crunch_len);
  For(i,n_otu)
    {
      For(j,23)
	{
	  if(j<(int)strlen(alldata->c_seq[i]->name))
	    fputc(alldata->c_seq[i]->name[j],fp);
	  else fputc(' ',fp);
	}
      
      For(j,alldata->crunch_len)
	{	  
	  For(k,alldata->wght[j])
	    fprintf(fp,"%c",alldata->c_seq[i]->state[j]);
	}
      fprintf(fp,"\n");
    }
  fprintf(fp,"\n");
  
  /*   printf("\t");   */
  /*   For(j,alldata->crunch_len) */
  /*     printf("%.0f",alldata->wght[j]); */
  /*   printf("\n"); */
}

/*********************************************************/

void Print_Translated_Alignment_DNA_to_AA(FILE *fp, allseq *alldata, code *c_code)
{
    int i,j,k;
    int n_otu;
    char codon[3];
    int codon_num,corrected_codon_num;
    
    n_otu = alldata->n_otu;
    fprintf(fp,"%d\t%d\n",n_otu,alldata->crunch_len/3);
    For(i,n_otu)
        {
            For(j,23)
                {
                    if(j<(int)strlen(alldata->c_seq[i]->name))
                        fputc(alldata->c_seq[i]->name[j],fp);
                    else fputc(' ',fp);
                }
            
            for(j=0;j<alldata->crunch_len;j+=3) // read codon by codon
	      {	  
		For(k,3) 
		  {
		    codon[k] = alldata->c_seq[i]->state[j+k];
		    if((codon[k] != 'A') && 
		       (codon[k] != 'C') &&
		       (codon[k] != 'G') &&
		       (codon[k] != 'T')) break;
		  }

		if(k != 3) 
		  {
		    fprintf(fp,"-");
		  }
		else
		  {
		    codon_num = Get_Base_Num(codon[0])*16 + Get_Base_Num(codon[1])*4 + Get_Base_Num(codon[2]);
		    corrected_codon_num  = c_code->sense_c[codon_num];
		    fprintf(fp,"%c",c_code->aa[c_code->genetc[c_code->num_curr_code][corrected_codon_num]]);
		  }
	      }
            fprintf(fp,"\n");
        }
    fprintf(fp,"\n");
    
    /*   printf("\t");   */
    /*   For(j,alldata->crunch_len) */
    /*     printf("%.0f",alldata->wght[j]); */
    /*   printf("\n"); */
}


/*********************************************************/

void Order_Tree_Seq(arbre *tree, seq **data)
{
  /* Put the uncompressed sequences in the same order as in the tree */

  int i,j,n_otu;
  seq *buff;

  n_otu = tree->n_otu;

  For(i,n_otu)
    {
      For(j,n_otu)
	{
            if(!strcmp(tree->noeud[i]->name,data[j]->name))
                break;
	}
      buff = data[j];
      data[j] = data[i];
      data[i] = buff;
    }
}

/*********************************************************/

void Order_Tree_CSeq(arbre *tree, allseq *data)
{
  /* Put the compressed sequences in the same order as in the tree */

  int i,j,n_otu_tree,n_otu_seq;
  seq *buff;
      
      
  n_otu_tree = tree->n_otu;
  n_otu_seq  = data->n_otu;
  
  if(n_otu_tree != n_otu_seq) 
    Exit("\n. The number of tips in the tree is not the same as the number of sequences\n");

  For(i,MMAX(n_otu_tree,n_otu_seq))
    {
      For(j,MMIN(n_otu_tree,n_otu_seq))
	{
            if(!strcmp(tree->noeud[i]->name,data->c_seq[j]->name))
                break;
	}

      if(j==MMIN(n_otu_tree,n_otu_seq))
	{
	  printf("\n. Err: %s is not found in sequences data set\n",
		 tree->noeud[i]->name);
	  Exit("");
	}

      buff = data->c_seq[j];
      data->c_seq[j] = data->c_seq[i];
      data->c_seq[i] = buff;

    }
}

/*********************************************************/

int Sort_String(const void *a, const void *b)
{
  return(strcmp((*(const char **)(a)), (*(const char **)(b))));
}

/*********************************************************/

void Get_SubTrees_Numbers(arbre *tree)
{
  int i;
  int curr_num;

  curr_num = 0;
  For(i,2*tree->n_otu-3)
    {
      tree->t_edges[i]->num_st_left = curr_num;
      curr_num++;
      tree->t_edges[i]->num_st_rght = curr_num;
      curr_num++;
    }
}

/*********************************************************/

optimiz *Alloc_Optimiz()
{
  /* Allocate memory for the optimisation structure */

  optimiz *s_opt;
  s_opt = (optimiz *)mCalloc(1,sizeof(optimiz));
  return s_opt;
}

/*********************************************************/

void Init_Optimiz(optimiz *s_opt)
{
  /* Initialize the optimisation structure */

  s_opt->sort_omega      = 1;
  s_opt->print           = 1;
  s_opt->opt_alpha       = 0;
  s_opt->opt_bfreq       = 0;
  s_opt->opt_kappa       = 0;
  s_opt->opt_bl          = 0;
  s_opt->opt_pinvar      = 0;
  s_opt->opt_param       = 1;
  s_opt->opt_subst_rate  = 0;
  s_opt->opt_tpos        = 0;
  s_opt->opt_theta       = 1;
  s_opt->opt_p_omega     = 1;

  s_opt->last_alpha      = 1.E-6;
}

/*********************************************************/
	
int Filexists(char *filename)
{ 
  /* Test if a file named 'filename' already exists in the
     current directory */

  FILE *fp;

  fp =fopen(filename,"r");
  if (fp) {
    fclose(fp);
    return 1;
  } else
    return 0;
}

/*********************************************************/

FILE *Openfile(char *filename, int mode)
{
  /* mode = 0 -> read */
  /* mode = 1 -> write */
  /* mode = 2 -> append */

  FILE *fp;
  char *s;
  int try;


/*   s = (char *)mCalloc(T_MAX_FILE,sizeof(char)); */

/*   strcpy(s,filename); */

  s = filename;

  try = 0;
  while(!(int)strlen(s))
    {
      printf("\nCan't open file '', enter a new name : ");
      Getstring_Stdin(s);
      fflush(NULL);
      if(++try == 10) Exit("\n");
    };


  fp = NULL;

  switch(mode)
    {
    case 0 :
      {
	try = 0;
	while(!(fp = (FILE *)fopen(s,"r")))
	  {
	    printf("\nCan't open file '%s', enter a new name : ",s);
	    Getstring_Stdin(s);
	    fflush(NULL);
	    if(++try == 10) Exit("\n");
	  }
	break;
      }
    case 1 :
      {
	fp = (FILE *)fopen(s,"w");
	break;
      }
    case 2 :
      {
	fp = (FILE *)fopen(s,"a");
	break;
      }
 
    default : break;
    
    }

/*   Free(s); */

  return fp;
}

/*********************************************************/

void Print_Fp_Out(FILE *fp_out, time_t t_beg, time_t t_end, arbre *tree, option *input, int n_data_set)
{

  /* Print the statistics */

  char *s,*subst_modelname,*switch_modelname;
  div_t hour,min;

  s = (char *)mCalloc(T_MAX_LINE,sizeof(char));

  fprintf(fp_out,". Sequence file : [%s]\n\n", input->seqfile);
  
  fprintf(fp_out,". Data set [#%d]\n\n",n_data_set);

  subst_modelname  = Return_Subst_Model_Name(input->mod);
  switch_modelname = Return_Switch_Model_Name(input->mod);

  (tree->mod->datatype == NT)?
    ((tree->mod->model_number >= 20)?
     ((tree->mod->switch_modelname == NO_SWITCH)?
      (fprintf(fp_out,". Model of codon substitution : %s\n\n",subst_modelname)):
      (fprintf(fp_out,". Model of codon substitution : %s+%s\n\n",subst_modelname,switch_modelname))):
     (fprintf(fp_out,". Model of nucleotide substitution : %s\n\n",subst_modelname))):
    (fprintf(fp_out,". Model of amino acid substitution : %s\n\n",subst_modelname));

  if(tree->mod->model_number <= 5)
    {
      fprintf(fp_out,". Transition/transversion ratio : %.3f\n\n",tree->mod->kappa);
    }
  else if(tree->mod->model_number == 6)
    {
      fprintf(fp_out,". Transition/transversion ratio for purines :     %.3f\n",
	      tree->mod->kappa*2.*tree->mod->lambda/(1.+tree->mod->lambda));
      fprintf(fp_out,". Transition/transversion ratio for pyrimidines : %.3f\n\n",
	      tree->mod->kappa*2./(1.+tree->mod->lambda));
    }

  fprintf(fp_out,". Discrete gamma model : %s\n",
	  (tree->mod->n_catg>1)?("Yes"):("No\n"));
  if(tree->mod->n_catg > 1)
    {
      fprintf(fp_out,"  - Number of categories : %d\n",tree->mod->n_catg);
      fprintf(fp_out,"  - Gamma shape parameter : %.3f\n\n",tree->mod->alpha);
    }
  
  if(tree->mod->invar)
    fprintf(fp_out,". Proportion of invariant : %.3f\n\n",tree->mod->pinvar);

  if(tree->mod->model_applies_to == NT)
    {
      fprintf(fp_out,". Nucleotides frequencies :\n\n");
      fprintf(fp_out,"  - f(A)=%8.5f\n",tree->mod->pi[0]);
      fprintf(fp_out,"  - f(C)=%8.5f\n",tree->mod->pi[1]);
      fprintf(fp_out,"  - f(G)=%8.5f\n",tree->mod->pi[2]);
      fprintf(fp_out,"  - f(T)=%8.5f\n\n",tree->mod->pi[3]);
    }


  if(tree->mod->subst_modelname == GTR)
    {
      int i,j;
      
      printf("\n");
      fprintf(fp_out,". GTR rate parameters : \n\n");
      fprintf(fp_out," A <-> C   %8.5f\n",tree->mod->gtr_param[0]);
      fprintf(fp_out," A <-> G   %8.5f\n",tree->mod->gtr_param[1]);
      fprintf(fp_out," A <-> T   %8.5f\n",tree->mod->gtr_param[2]);
      fprintf(fp_out," C <-> G   %8.5f\n",tree->mod->gtr_param[3]);
      fprintf(fp_out," C <-> T   %8.5f\n",tree->mod->gtr_param[4]);
      fprintf(fp_out," G <-> T   1.0 (fixed)\n");
      
      fprintf(fp_out,"\n. Instantaneous rate matrix : \n");
      fprintf(fp_out,"\n[A---------C---------G---------T------]\n");
      For(i,4)
	{
	  For(j,4)
	    fprintf(fp_out,"%8.5f  ",tree->mod->qmat_struct[0]->qmat[i*4+j]);
	  fprintf(fp_out,"\n");
	}
      fprintf(fp_out,"\n");
      fprintf(fp_out,"eg., the instantaneous rate of change from 'C' to 'A' is %8.5f x %8.5f = %8.5f\n\n",
	      tree->mod->pi[0],
	      tree->mod->gtr_param[0],
	      tree->mod->qmat_struct[0]->qmat[1*4+0]);
    }
  
  
  if(tree->mod->model_applies_to == CODONS)
    {
      if(tree->mod->switch_modelname == NO_SWITCH)
	{
          int i,j;
	  fprintf(fp_out,". Transition/transversion ratio : %.3f\n\n",tree->mod->kappa);
	  
	  For(i,tree->mod->n_catq)
            fprintf(fp_out,". p%d =\t\t%f\n",
                    i+1,
                    tree->mod->qmat_struct[0]->qmat_proba[i]);
	  For(i,tree->mod->n_catq)
            fprintf(fp_out,". w%d =\t\t%f\n",
                    i+1,
                    tree->mod->qmat_struct[0]->omega[i]);

          
          fprintf(fp_out,"\n. Posterior probabilities of selection regimes at individual sites\n");
          For(i,tree->n_pattern)
	    {
	      fprintf(fp_out,"%5d\t",i+1);
	      For(j,tree->mod->n_catq)
		{
		  fprintf(fp_out,"%10f\t",tree->sel_regime_prob[j][i]);
		}
	      
	      if((tree->input->mod->subst_modelname != M0) && 
		 (tree->input->mod->subst_modelname != M1) &&                     
		 (tree->input->mod->subst_modelname != M1a)) 
		{
		  if(tree->sel_regime_prob[2][i] > 0.99) fprintf(fp_out,"*");
		  if(tree->sel_regime_prob[2][i] > 0.95) fprintf(fp_out,"*");
		}
	      fprintf(fp_out,"\n");
	    }
	}

      if((tree->mod->switch_modelname == SWITCH_S1) || (tree->mod->switch_modelname == SWITCH_S2))
	{
          int i,j;
	  
	  fprintf(fp_out,". Transition/transversion ratio : %.3f\n\n",tree->mod->kappa);
	  
	  fprintf(fp_out,". Codon model parameters : \n\n");
	  
	  /* fprintf(fp_out,". p1 =\t\t%f\n. p2 =\t\t%f\n. p3 =\t\t%f\n. w1 =\t\t%f\n. w2 =\t\t%f\n. w3 =\t\t%f\n. delta =\t%f\n. R01 =\t%f\n. R02 =\t%f\n. R12 =\t%f\n", */
	  /* 	  tree->mod->qmat_struct[0]->omega_proba[0], */
	  /* 	  tree->mod->qmat_struct[0]->omega_proba[1], */
	  /* 	  tree->mod->qmat_struct[0]->omega_proba[2], */
	  /* 	  tree->mod->qmat_struct[0]->omega[0], */
	  /* 	  tree->mod->qmat_struct[0]->omega[1], */
	  /* 	  tree->mod->qmat_struct[0]->omega[2], */
	  /* 	  tree->mod->delta_switch, */
	  /* 	  tree->mod->qmat_struct[0]->theta[0]/(tree->mod->delta_switch), */
	  /* 	  tree->mod->qmat_struct[0]->theta[1]/(tree->mod->delta_switch), */
	  /* 	  tree->mod->qmat_struct[0]->theta[2]/(tree->mod->delta_switch)); */
	  
	  For(i,tree->mod->n_omega)
	    {
	      fprintf(fp_out,"\n. p%d = %f | w%d = %f ",
	  	      i,tree->mod->qmat_struct[0]->omega_proba[i],
	  	      i,tree->mod->qmat_struct[0]->omega[i]);
	    }
	  fprintf(fp_out,"\n. delta = %f",tree->mod->delta_switch);

	  For(i,tree->mod->n_omega-1)
	    for(j=i+1;j<tree->mod->n_omega;j++)
	      {
	  	fprintf(fp_out,"\n. R%d%d = %f",i,j,
			tree->mod->qmat_struct[0]->theta[MMIN(i,j) * tree->mod->n_omega + MMAX(i,j) -
							 (MMIN(i,j)+1+(int)pow(MMIN(i,j)+1,2))/2]);
	      }
          
          fprintf(fp_out,"\n. Posterior probabilities of selection regimes at individual sites\n");
          For(i,tree->n_pattern)
	    {
	      fprintf(fp_out,"%5d\t",i+1);
	      For(j,tree->mod->n_omega)
		{
		  fprintf(fp_out,"%10f\t",tree->sel_regime_prob[j][i]);
		}
	      
	      
	      if((tree->input->mod->subst_modelname != M0) && 
		 (tree->input->mod->subst_modelname != M1) &&
		 (tree->input->mod->subst_modelname != M1a) && 
		 (tree->input->mod->subst_modelname != MX)) 
		{
		  if(tree->sel_regime_prob[2][i] > 0.99) fprintf(fp_out,"*");
		  if(tree->sel_regime_prob[2][i] > 0.95) fprintf(fp_out,"*");
		}
	      fprintf(fp_out,"\n");
	    }
	}
    }
      
  fprintf(fp_out,"\n. log likelihood = %f\n\n",tree->tot_loglk);
  
  /*   For(site,tree->n_pattern) */
  /*     { */
  /*       if((tree->sel_regime_prob[2][site] > 0.95) && (tree->mod->qmat_struct[0]->omega[2] > 1.0)) */
  /* 	{ */
  /* 	  fprintf(input->fp_output_stat,". Site %3d p=%f %f %f\n", */
  /* 		  site+1, */
  /* 		  tree->sel_regime_prob[2][site], */
  /* 		  tree->sel_regime_prob[1][site], */
  /* 		  tree->sel_regime_prob[0][site]); */
  
  /* 	} */
  /*     } */
  /*   fprintf(input->fp_output_stat,"\n"); */
  
  
  hour = div(t_end-t_beg,3600);
  min  = div(t_end-t_beg,60  );
  
  min.quot -= hour.quot*60;
  
  fprintf(fp_out,". Time used %dh%dm%ds\n", hour.quot,min.quot,(int)(t_end-t_beg)%60);
  if(t_end-t_beg > 60)
    fprintf(fp_out,". -> %d seconds\n",(int)(t_end-t_beg));
  
  fprintf(fp_out,"\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n\n");
  
  fflush(NULL);
  
  Free(subst_modelname);
  Free(switch_modelname);
  Free(s);
}
  
/*********************************************************/

void Alloc_All_P_Lk(arbre *tree)
{
  /* Allocate memory for the vectors of partial likelihood */
  
  int i,j,k,l;


  For(i,2*tree->n_otu-3)
    {
      tree->t_edges[i]->get_p_lk_left = 1; 
      tree->t_edges[i]->get_p_lk_rght = 1;

      /* There is one vector on the left side and one on the right side
	 of each edge */

      tree->t_edges[i]->p_lk_left = 
      (fit_double ****)mCalloc(tree->n_pattern,sizeof(fit_double ***));
      
      tree->t_edges[i]->p_lk_rght = 
      (fit_double ****)mCalloc(tree->n_pattern,sizeof(fit_double ***));
      
      /* printf("\n. MAKE %p %p", */
      /*        tree->t_edges[i], */
      /*        tree->t_edges[i]->p_lk_rght); */
      /* fflush(NULL); */

      For(j,tree->n_pattern)
	{
	  tree->t_edges[i]->p_lk_left[j] = 
	  (fit_double ***)mCalloc(tree->mod->n_catq,sizeof(fit_double **));
	  
	  For(k,tree->mod->n_catq)
	    {
	      tree->t_edges[i]->p_lk_left[j][k] = 
		(fit_double **)mCalloc(tree->mod->n_catg,sizeof(fit_double *));

	      For(l,tree->mod->n_catg)
		tree->t_edges[i]->p_lk_left[j][k][l] = 
		(fit_double *)mCalloc(tree->mod->ns,sizeof(fit_double ));
	    }


	  tree->t_edges[i]->p_lk_rght[j] = 
	  (fit_double ***)mCalloc(tree->mod->n_catq,sizeof(fit_double **));
	  
	  For(k,tree->mod->n_catq)
	    {
	      tree->t_edges[i]->p_lk_rght[j][k] = 
		(fit_double **)mCalloc(tree->mod->n_catg,sizeof(fit_double *));

	      For(l,tree->mod->n_catg)
		tree->t_edges[i]->p_lk_rght[j][k][l] = 
		(fit_double *)mCalloc(tree->mod->ns,sizeof(fit_double ));
	    }
	}

      tree->t_edges[i]->sum_scale_f_left = 
      (fit_double *)mCalloc(tree->n_pattern,sizeof(fit_double ));
      
      tree->t_edges[i]->sum_scale_f_rght = 
      (fit_double *)mCalloc(tree->n_pattern,sizeof(fit_double ));
      
    }
}

/*********************************************************/

int Is_Ambigu(char *state, int datatype, int stepsize)
{
  /* Is the state 'state' an ambiguous character ? */

  int i;

  if(datatype == NT) 
    {
      For(i,stepsize)
	{
	  if(strchr("MRWSYKBDHVNXO?-.",state[i]))
	    return 1;
	}
    }
  else if(datatype == AA)
    {
      if(strchr("X?-.",state[0])) return 1;       
    }
  else 
    {
      printf("\n. datatype : %d\n",datatype);
      printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Exit("");
    }
  return 0;
}

/*********************************************************/

int Is_Completely_Ambigu(char *state, int datatype, int stepsize)
{
  /* Is the state 'state' a completely ambiguous character (e.g. 'X'
     is a completely ambiguous character, 'R' (for DNA sequences) 
     is not) */
    
  int i;

  if(datatype == NT) 
    {
      For(i,stepsize)
	{
	  if(strchr("NXO?-.",state[i]))
	    return 1;
	}
    }
  else if(datatype == AA)
    {
      if(strchr("X?-.",state[0])) return 1;       
    }
  else
    {
      printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Exit("");
    }
  return 0;
}

/*********************************************************/

void Check_Ambiguities(allseq *data, int datatype, int stepsize)
{
  /* Identify ambiguous sites */

  int i,j;
 
  Fors(j,data->crunch_len,stepsize) For(i,data->n_otu)
    {
      if(Is_Ambigu(data->c_seq[i]->state+j,
		   datatype,
		   stepsize))
	{
	  data->ambigu[j] = 1;
	  break;
	}
    }
}

/*********************************************************/

int Assign_State(char *c, int datatype, int stepsize, code *c_code)
{
  /* Get state from 'c' */

  int state[3];
  int i;

  state[0] = -1;
  if(datatype == NT) /* nucleotides */
    {	  
      For(i,stepsize)
	{
	  switch(c[i])
	    {
	    case 'A' : state[i]=0; break;
	    case 'C' : state[i]=1; break;
	    case 'G' : state[i]=2; break;
	    case 'T' : state[i]=3; break;
	    case 'U' : state[i]=3; break;
	    default  : state[i]=-1;
	      break;
	    }
	}


      return (stepsize>1)?(c_code->from_64_2_61[state[0]*16+state[1]*4+state[2]]):(state[0]);
    }
  else if(datatype == AA)/* amino acids */
    {
      switch(c[0]){
      case 'A' : state[0]=0;  break;
      case 'R' : state[0]=1;  break;
      case 'N' : state[0]=2;  break;
      case 'D' : state[0]=3;  break;
      case 'C' : state[0]=4;  break;
      case 'Q' : state[0]=5;  break;
      case 'E' : state[0]=6;  break;
      case 'G' : state[0]=7;  break;
      case 'H' : state[0]=8;  break;
      case 'I' : state[0]=9;  break;
      case 'L' : state[0]=10; break;
      case 'K' : state[0]=11; break;
      case 'M' : state[0]=12; break;
      case 'F' : state[0]=13; break;
      case 'P' : state[0]=14; break;
      case 'S' : state[0]=15; break;
      case 'T' : state[0]=16; break;
      case 'W' : state[0]=17; break;
      case 'Y' : state[0]=18; break;
      case 'V' : state[0]=19; break;
	
      case 'B' : state[0] = 2; break;
      case 'Z' : state[0] = 5; break;
	
      default : state[0]=-1;
	break;
      }
      return state[0];
    }
  else
    {
      printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Exit("");
    }
  return -1;
}

/*********************************************************/

void Update_BrLen_Invar(arbre *tree)
{
  /* Adjust branch length by taking into account a proportion of invariable sites */

  int i;
  For(i,2*tree->n_otu-3) tree->t_edges[i]->l*=(1.0-tree->mod->pinvar);
}

/*********************************************************/

void Getstring_Stdin(char *file_name)
{ 
  /* Read a string of characters from standard input */
 
  if(!fgets(file_name,T_MAX_LINE,stdin))
    {
      printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Exit("");
    }
  if (strchr(file_name, '\n') != NULL)
    *strchr(file_name, '\n') = '\0';
  if (strchr(file_name, ' ') != NULL)
    *strchr(file_name, ' ') = '\0';
}

/*********************************************************/

allseq *Evolve(arbre *tree)
{
  /* Simulate sequences along a phylogenetic tree given a model of substitution */

  int i,j;
  int root_state,rrate_cat,qmat_cat;
  int br;
  fit_double *stationnary_p,sum;
  edge *eroot;
  allseq *sim_seq;
  int n_switches;

  stationnary_p = (fit_double *)mCalloc(tree->mod->ns,sizeof(fit_double));

  eroot = tree->noeud[0]->b[0];
  
  if(tree->mod->model_number > 20)
    {
      sim_seq = Copy_Cseq(tree->data, tree->data->crunch_len, tree->mod->ns_codon);

      /* init stationnary frequencies for switching models */
      if(tree->mod->switch_modelname != NO_SWITCH)
	{
	  For(i,tree->mod->n_catq)
	    {
	      For(j,tree->mod->ns)
		{
		  stationnary_p[i*tree->mod->n_catq + j] = 
		    tree->mod->pi[j%tree->mod->ns_codon] * 
		    (*(eroot->qmat_struct->t_omega_proba[i*eroot->qmat_struct->n_omega+
							 (int)j/tree->mod->ns_codon]));
		}
	    }
	}
      else /* non-switching model */
	{
	  For(j,tree->mod->ns) stationnary_p[j] = tree->mod->pi[j];
	}
      

      sum = 0.0;
      For(i,tree->mod->ns) sum += stationnary_p[i];
      /* For(i,tree->mod->ns) printf("\n. State freq %3d: %15f",i,stationnary_p[i]); */
      /* printf("\n. Sum of freq = %f\n",sum); */
    }
  else
    {
      sim_seq = Copy_Cseq(tree->data, tree->data->crunch_len, tree->mod->ns);
      For(i,tree->mod->ns) stationnary_p[i] = tree->mod->pi[i];
    }
        
  tree->curr_site = 0;
    
  DiscreteGamma (tree->mod->r_proba, tree->mod->rr_mixturem, tree->mod->alpha,
		 tree->mod->alpha,tree->mod->n_catg,0);
  
  Print_Param(tree->input->fp_output_stat,tree);
  Print_Param(stderr,tree);

  For(br,2*tree->n_otu-3) 
    {
      Get_PMat(tree,
	       tree->t_edges[br]->l,
	       tree->t_edges[br]->Pij_rr,
	       tree->t_edges[br]->qmat_struct);
    }
    

  while(tree->curr_site < tree->n_pattern)
    {
      
      For(br,2*tree->n_otu-3)
	{
	  tree->t_edges[br]->sum_scale_f_rght[tree->curr_site] = .0;
	  tree->t_edges[br]->sum_scale_f_left[tree->curr_site] = .0;
	}
      
      root_state = Generate_One_State(tree->mod->ns,stationnary_p);
      
      /* the root is also a tip of the tree */
      Return_State(tree,
		   root_state,
		   sim_seq->c_seq[tree->n_root->num]->state+tree->curr_site*tree->mod->stepsize);
           
      rrate_cat = Generate_One_State(tree->mod->n_catg,tree->mod->r_proba);
      qmat_cat  = Generate_One_State(tree->mod->n_catq,tree->mod->qmat_struct[0]->qmat_proba);
      
#if defined(EVOLVE)
      /* fprintf(tree->input->fp_output_stat, */
      /*         "%3d %3d %10f\n", */
      /*         tree->curr_site, */
      /*         qmat_cat,tree->mod->qmat_struct[0]->omega[qmat_cat]); */
#endif
      
      n_switches = 0;
      Generate_One_Site(tree->noeud[0],
			tree->noeud[0]->v[0],
			root_state,
			tree->noeud[0]->b[0],
			rrate_cat,qmat_cat,sim_seq,&n_switches,tree);

      printf("\n. Site %4d # of selection class switches: %4d",tree->curr_site,n_switches);

      For(i,tree->mod->stepsize) sim_seq->wght[tree->curr_site*tree->mod->stepsize+i] = 1;
      
      tree->curr_site+=1;      
    }
  
  /* Add_Ambiguities(tree->data,sim_seq,tree->mod->datatype); */

  Free(stationnary_p);

  return sim_seq;
}

/*********************************************************/

void Evolve_Under_H0(arbre *tree)
{
  /* Simulate sequences along a phylogenetic tree given a model of substitution */

    int i,j;
    int root_state,rrate_cat,qmat_cat;
    int br;
    fit_double *stationnary_p;
    allseq *sim_seq;

    stationnary_p = (fit_double *)mCalloc(tree->mod->ns,sizeof(fit_double));

    if(tree->mod->model_number > 20)
        {
	  sim_seq = Copy_Cseq(tree->data,tree->data->init_len,tree->mod->ns_codon);

	  /* init stationnary frequencies for switching models */
	  For(i,tree->mod->n_catq)
	    For(j,tree->mod->ns_codon)
	    stationnary_p[i*tree->mod->n_catq + j] *= tree->mod->omega_proba[j];
        }
    else
        {
	  sim_seq = Copy_Cseq(tree->data,tree->data->init_len,tree->mod->ns);
	  For(i,tree->mod->ns) stationnary_p[i] = tree->mod->pi[i];
        }
    
    tree->curr_site = 0;
    DiscreteGamma (tree->mod->r_proba, tree->mod->rr_mixturem, tree->mod->alpha,
                   tree->mod->alpha,tree->mod->n_catg,0);
    
    Print_Param(tree->input->fp_output_stat,tree);
    
    For(br,2*tree->n_otu-3) 
        {
            Get_PMat(tree,
                     tree->t_edges[br]->l,
                     tree->t_edges[br]->Pij_rr,
                     tree->t_edges[br]->qmat_struct);
        }
    
    
    while(tree->curr_site < tree->n_pattern)
        {
            For(br,2*tree->n_otu-3) 
                {
                    tree->t_edges[br]->sum_scale_f_rght[tree->curr_site] = .0;
                    tree->t_edges[br]->sum_scale_f_left[tree->curr_site] = .0;
                }
            
            root_state = Generate_One_State(tree->mod->ns,stationnary_p);
            
            
            /* the root is also a tip of the tree */
            Return_State(tree,root_state,sim_seq->c_seq[tree->n_root->num]->state+tree->curr_site*tree->mod->stepsize);
            
            rrate_cat = Generate_One_State(tree->mod->n_catg,tree->mod->r_proba);
            
            qmat_cat  = Generate_One_State(tree->mod->n_catq,tree->mod->qmat_struct[0]->qmat_proba);
            
            if(qmat_cat != 2) // Generate under H0
                {
                    
                    tree->data->selclass[tree->curr_site] = qmat_cat;
                    
#if defined(EVOLVE)
                    fprintf(tree->input->fp_output_stat,
                            "%d %d\n",
                            tree->curr_site,
                            qmat_cat);
#endif

                    Generate_One_Site(tree->noeud[0],
                                      tree->noeud[0]->v[0],
                                      root_state,
                                      tree->noeud[0]->b[0],
                                      rrate_cat,qmat_cat,sim_seq,NULL,tree);

                    For(i,tree->mod->stepsize) sim_seq->wght[tree->curr_site*tree->mod->stepsize+i] = 1;
                    tree->curr_site+=1;
                }
        }
    Free(stationnary_p);
}

/*********************************************************/

void Generate_One_Site(node *a, node *d, int anc_state, edge  *b, int rrate, int qmat_cat, allseq *sim_seq, int *n_switches, arbre *tree)
{
  /* Simulate one site */
  
  int d_state,i;
  char *s;
  
  s = (char *)mCalloc(4,sizeof(char));
  
  d_state = -1;
  
  d_state = Generate_One_State(tree->mod->ns,b->Pij_rr[qmat_cat][rrate][anc_state]);
  
  if((tree->mod->model_applies_to == CODONS) &&  (d_state/tree->mod->ns_codon != anc_state/tree->mod->ns_codon))
    {
      if(n_switches) (*n_switches)++;
    }

  if(d_state != anc_state) tree->n_changes_per_site[tree->curr_site] += 1;
  
  Return_State(tree,d_state,s);

  s[3] = '\0';

#if defined(EVOLVE)
  strcat(d->seq+tree->curr_site*tree->mod->stepsize,s);
#endif
  
  if(d->tax)
    {
      Return_State(tree,
		   d_state,
		   sim_seq->c_seq[d->num]->state+tree->curr_site*tree->mod->stepsize);
      
      /*       printf("%s\n",s);  */
      Free(s);
      return;
    }
  else
    {
      For(i,3)
	if(d->v[i] != a)
	  Generate_One_Site(d,d->v[i],d_state,d->b[i],rrate,qmat_cat,sim_seq,n_switches,tree);
    }

  Free(s);
  return;
}			

/*********************************************************/

int Generate_One_State(int ns, fit_double *p)
{
  /* Simulate one state */

  int state;

  state = 0;
  do
    {
/*       state = (int)(ns*(random()/(RAND_MAX+1.0))); */
/*       if((fit_double)random()/(RAND_MAX+1.0) < p[state]) break; */
      state = (int)(ns*(rand()/(RAND_MAX+1.0)));
      if((fit_double)rand()/(RAND_MAX+1.0) < p[state]) break;
    }
  while(1);

  return state;
}

/*********************************************************/

void Return_State(arbre *tree, int num_state, char *state)
{
  /* Convert state number to character */

  switch(tree->mod->datatype)
    {
    case NT : 
      {
	if(tree->mod->model_applies_to == NT)
	  {
	    /* nucleotides */
	    switch(num_state)
	      {
	      case 0: {state[0] = 'A'; break;}
	      case 1: {state[0] = 'C'; break;}
	      case 2: {state[0] = 'G'; break;}
	      case 3: {state[0] = 'T'; break;}
	      case 4: {state[0] = 'U'; break;}
	      case 5: {state[0] = 'M'; break;}
	      case 6: {state[0] = 'R'; break;}
	      case 7: {state[0] = 'W'; break;}
	      case 8: {state[0] = 'S'; break;}
	      case 9: {state[0] = 'Y'; break;}
	      case 10:{state[0] = 'K'; break;}
	      case 11:{state[0] = 'B'; break;}
	      case 12:{state[0] = 'D'; break;}
	      case 13:{state[0] = 'H'; break;}
	      case 14:{state[0] = 'V'; break;}
	      case 15:{state[0] = 'N'; break;}
	      default : {break;}
	      }
	  }
	else if(tree->mod->model_applies_to == CODONS)
	  {
	    char *tmp_state;
	    int i;
	   
	    /* codons */
	    num_state = tree->mod->c_code->sense_c[num_state%tree->mod->ns_codon];
	    tmp_state = Get_Codon(num_state);
	    For(i,3) state[i] = tmp_state[i];
	    Free(tmp_state);
	  }
	else
	  {
	    printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	    Exit("");
	  }
	break;
      }
    case AA :
      {
	/* aa */
	switch(num_state)
	  {
	  case 0: {state[0] =  'A'; break;}
	  case 1: {state[0] =  'C'; break;}
	  case 2: {state[0] =  'D'; break;}
	  case 3: {state[0] =  'E'; break;}
	  case 4: {state[0] =  'F'; break;}
	  case 5: {state[0] =  'G'; break;}
	  case 6: {state[0] =  'H'; break;}
	  case 7: {state[0] =  'I'; break;}
	  case 8: {state[0] =  'K'; break;}
	  case 9: {state[0] =  'L'; break;}
	  case 10:{state[0] =  'M'; break;}
	  case 11:{state[0] =  'N'; break;}
	  case 12:{state[0] =  'P'; break;}
	  case 13:{state[0] =  'Q'; break;}
	  case 14:{state[0] =  'R'; break;}
	  case 15:{state[0] =  'S'; break;}
	  case 16:{state[0] =  'T'; break;}
	  case 17:{state[0] =  'V'; break;}
	  case 18:{state[0] =  'W'; break;}
	  case 19:{state[0] =  'Y'; break;}
	  case 20:{state[0] =  'X'; break;}
	  default : {break;}
	  }
	break;
      }
    default : break;
    }
}

/*********************************************************/

fit_double Num_Derivatives_One_Param(fit_double (*func)(arbre *tree), arbre *tree,
				 fit_double f0, fit_double *param, fit_double stepsize,
				 fit_double *err, int precise)
{
  int i,j;
  fit_double errt,fac,hh,**a,ans;
  int n_iter;
  a = (fit_double **)mCalloc(11,sizeof(fit_double *));
  For(i,11) a[i] = (fit_double *)mCalloc(11,sizeof(fit_double));


  n_iter = 10; /* */

  ans  = .0;

  if(stepsize < SMALL) Exit("\n. h must be nonzero in Dfridr.");

  hh=stepsize;

  if(!precise)
    {
      *param   = *param+hh;
      a[0][0]  = (*func)(tree);
/*       printf("\n. f0=%f f1=%f hh=%G",f0,a[0][0],hh); */
      a[0][0]  -= f0;
      a[0][0]  /= hh;
      *param   = *param-hh;

      ans =  a[0][0];

    }
  else
    {
      *param   = *param+hh;
      a[0][0]  = (*func)(tree);
      /*   *param   = *param-2*hh; */
      /*   a[0][0] -= (*func)(tree); */
      /*   a[0][0] /= (2.0*hh); */
      /*   *param   = *param+hh; */
      a[0][0]  -= f0;
      a[0][0]  /= hh;
      *param   = *param-hh;

      *err=1e30;
      for(i=1;i<n_iter;i++)
	{
	  hh /= 1.4;

	  /*       *param   = *param+hh; */
	  /*       a[0][i]  = (*func)(tree); */
	  /*       *param   = *param-2*hh; */
	  /*       a[0][i] -= (*func)(tree); */
	  /*       a[0][i] /= (2.0*hh); */
	  /*       *param   = *param+hh; */


	  *param   = *param+hh;
	  a[0][i]  = (*func)(tree);
	  /*   *param   = *param-2*hh; */
	  /*   a[0][i] -= (*func)(tree); */
	  /*   a[0][i] /= (2.0*hh); */
	  /*   *param   = *param+hh; */
	  a[0][i]  -= f0;
	  a[0][i]  /= hh;
	  *param   = *param-hh;


	  fac=1.4*1.4;
	  for (j=1;j<=i;j++)
	    {
	      a[j][i]=(a[j-1][i]*fac-a[j-1][i-1])/(fac-1.0);
	      fac=1.4*1.4*fac;

	      errt=MAX(FABS(a[j][i]-a[j-1][i]),FABS(a[j][i]-a[j-1][i-1]));

	      if (errt <= *err)
		{
		  *err=errt;
		  ans=a[j][i];
		}
	    }

	  if(FABS(a[i][i]-a[i-1][i-1]) >= 2.0*(*err)) break;
	}
    }
  For(i,11) Free(a[i]);
  Free(a);

  return ans;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

int Num_Derivative_Several_Param(arbre *tree, fit_double *param, int n_param, fit_double stepsize,
				  fit_double (*func)(arbre *tree), fit_double *derivatives)
{
  int i;
  fit_double err,f0;

  f0 = (*func)(tree);

  For(i,n_param)
    {
      derivatives[i] = Num_Derivatives_One_Param(func,
						 tree,
						 f0,
						 param+i,
						 stepsize,
						 &err,
						 0
						 );
    }
  return 1;
}

/*********************************************************/

int Get_AA(code *c_code, int num_codon)
{
  /* return the amino acid number that correspond to codon 'num_codon' */
  return c_code->genetc[c_code->num_curr_code][num_codon];

}

/*********************************************************/

code *Make_Code(int num_code)
{
  /* Genetic codes 
     http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi?mode=c#thetop
  */

  int i,j;
  code *c_code;
  int number_of_codes;

  number_of_codes = 18;

  int genetc[18][64] ={
 
     {11,2,11,2,16,16,16,16,1,15,1,15,9,9,12,9,5,8,5,
      8,14,14,14,14,1,1,1,1,10,10,10,10,6,3,6,3,0,0,0,
      0,7,7,7,7,19,19,19,19,-1,18,-1,18,15,15,15,15,-1,
      4,17,4,10,13,10,13},    /* 0: universal */

     {11,2,11,2,16,16,16,16,-1,15,-1,15,12,9,12,9,5,8,
      5,8,14,14,14,14,1,1,1,1,10,10,10,10,6,3,6,3,0,0,
      0,0,7,7,7,7,19,19,19,19,-1,18,-1,18,15,15,15,15,
      17,4,17,4,10,13,10,13}, /* 1:vertebrate mt.*/
     
     {11,2,11,2,16,16,16,16,1,15,1,15,12,9,12,9,5,8,5,
      8,14,14,14,14,1,1,1,1,16,16,16,16,6,3,6,3,0,0,0,
      0,7,7,7,7,19,19,19,19,-1,18,-1,18,15,15,15,15,17,
      4,17,4,10,13,10,13},    /* 2: yeast mt. */

     {11,2,11,2,16,16,16,16,1,15,1,15,9,9,12,9,5,8,5,8,
      14,14,14,14,1,1,1,1,10,10,10,10,6,3,6,3,0,0,0,0,7,
      7,7,7,19,19,19,19,-1,18,-1,18,15,15,15,15,17,4,17,
      4,10,13,10,13},         /* 3: mold mt. */

     {11,2,11,2,16,16,16,16,15,15,15,15,12,9,12,9,5,8,
      5,8,14,14,14,14,1,1,1,1,10,10,10,10,6,3,6,3,0,0,
      0,0,7,7,7,7,19,19,19,19,-1,18,-1,18,15,15,15,15,
      17,4,17,4,10,13,10,13}, /* 4: invertebrate mt. */
     
     {11,2,11,2,16,16,16,16,1,15,1,15,9,9,12,9,5,8,5,8,
      14,14,14,14,1,1,1,1,10,10,10,10,6,3,6,3,0,0,0,0,
      7,7,7,7,19,19,19,19,5,18,5,18,15,15,15,15,-1,4,
      17,4,10,13,10,13},      /* 5: ciliate nu.*/
     
     {2,2,11,2,16,16,16,16,15,15,15,15,9,9,12,9,5,8,5,8,
      14,14,14,14,1,1,1,1,10,10,10,10,6,3,6,3,0,0,0,0,7,
      7,7,7,19,19,19,19,-1,18,-1,18,15,15,15,15,17,4,17,
      4,10,13,10,13},         /* 6: echinoderm mt.*/
     
     {11,2,11,2,16,16,16,16,1,15,1,15,9,9,12,9,5,8,5,8,
      14,14,14,14,1,1,1,1,10,10,10,10,6,3,6,3,0,0,0,0,7,
      7,7,7,19,19,19,19,-1,18,-1,18,15,15,15,15,4,4,17,
      4,10,13,10,13},         /* 7: euplotid mt. */
     
     {11,2,11,2,16,16,16,16,1,15,1,15,9,9,12,9,5,8,5,8,
      14,14,14,14,1,1,1,1,10,10,10,10,6,3,6,3,0,0,0,0,
      7,7,7,7,19,19,19,19,-1,18,-1,18,15,15,15,15,-1,4,
      17,4,10,13,10,13},      /* 8: bacterial and plant plastid */

     {11,2,11,2,16,16,16,16,1,15,1,15,9,9,12,9,5,8,5,8,
      14,14,14,14,1,1,1,1,10,10,15,10,6,3,6,3,0,0,0,0,
      7,7,7,7,19,19,19,19,-1,18,-1,18,15,15,15,15,-1,
      4,17,4,10,13,10,13},    /* 9: alternative yeast nu.*/

     {11,2,11,2,16,16,16,16,7,15,7,15,12,9,12,9,5,8,5,
      8,14,14,14,14,1,1,1,1,10,10,10,10,6,3,6,3,0,0,0,
      0,7,7,7,7,19,19,19,19,-1,18,-1,18,15,15,15,15,17,
      4,17,4,10,13,10,13},    /* 10: ascidian mt. */
     
     {2,2,11,2,16,16,16,16,15,15,15,15,9,9,12,9,5,8,5,
      8,14,14,14,14,1,1,1,1,10,10,10,10,6,3,6,3,0,0,0,
      0,7,7,7,7,19,19,19,19,18,18,-1,18,15,15,15,15,17,
      4,17,4,10,13,10,13},    /* 11: alternative flatworm mt. */

     {11,2,11,2,16,16,16,16,1,15,1,15,9,9,12,9,5,8,5,8,
      14,14,14,14,1,1,1,1,10,10,10,10,6,3,6,3,0,0,0,0,
      7,7,7,7,19,19,19,19,-1,18,5,18,15,15,15,15,-1,4,
      17,4,10,13,10,13},      /* 12: blepharisma nu. */
     
     {11,2,11,2,16,16,16,16,1,15,1,15,9,9,12,9,5,8,5,8,
      14,14,14,14,1,1,1,1,10,10,10,10,6,3,6,3,0,0,0,0,
      7,7,7,7,19,19,19,19,-1,18,10,18,15,15,15,15,-1,4,
      17,4,10,13,10,13},      /* 13: chlorophycean mt. */
     
     {2,2,11,2,16,16,16,16,15,15,15,15,12,9,12,9,5,8,5,
      8,14,14,14,14,1,1,1,1,10,10,10,10,6,3,6,3,0,0,0,
      0,7,7,7,7,19,19,19,19,-1,18,-1,18,15,15,15,15,17,
      4,17,4,10,13,10,13},    /* 14: trematode mt. */

     {11,2,11,2,16,16,16,16,1,15,1,15,9,9,12,9,5,8,5,8,
      14,14,14,14,1,1,1,1,10,10,10,10,6,3,6,3,0,0,0,0,7,
      7,7,7,19,19,19,19,-1,18,10,18,-1,15,15,15,-1,4,17,
      4,10,13,10,13},         /* 15: scenedesmus obliquus mt. */

     {11,2,11,2,16,16,16,16,1,15,1,15,9,9,12,9,5,8,5,8,
      14,14,14,14,1,1,1,1,10,10,10,10,6,3,6,3,0,0,0,0,
      7,7,7,7,19,19,19,19,-1,18,-1,18,15,15,15,15,-1,4,
      17,4,-1,13,10,13},       /* 16: thraustochytrium mt. */

     {0,0,0,0,1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4,5,5,5,5,6,6,6,6,
      7,7,7,7,8,8,8,8,9,9,9,9,10,10,10,10,11,11,11,11,12,12,12,12,
      13,13,13,13,14,14,14,14,15,15,15,15}       /* 17: Ziheng Yang's regular code */
    };

  c_code                    = (code *)mCalloc(1,sizeof(code));
  c_code->n_diff_b_2_codons = (int *)mCalloc(64*64,sizeof(int));
  c_code->tstvtable         = (int *)mCalloc(64*64,sizeof(int));
  c_code->gtr_4_codon       = (fit_double **)mCalloc(64*64,sizeof(fit_double *));
  c_code->sense_c           = (int *)mCalloc(64,sizeof(int));
  c_code->from_64_2_61      = (int *)mCalloc(64,sizeof(int));
  c_code->name              = (char *)mCalloc(100,sizeof(char));
  
  c_code->number_of_codes = number_of_codes;

  switch(num_code)
    {
    case 0 :
      {
	c_code->num_curr_code = 0;
	strcpy(c_code->name,"universal");
	c_code->n_sense_c = 61;
	break;
      }
    case 1 :
      {
	c_code->num_curr_code = 1;
	strcpy(c_code->name,"vertebrate mt.");
	c_code->n_sense_c = 60;
	break;
      }
    case 2 :
      {
	c_code->num_curr_code = 2;
	strcpy(c_code->name,"yeast mt.");
	c_code->n_sense_c = 62;
	break;
      }
    case 3 :
      {
	c_code->num_curr_code = 3;
	strcpy(c_code->name,"mold mt.");
	c_code->n_sense_c = 62;
	break;
      }
    case 4 :
      {
	c_code->num_curr_code = 4;
	strcpy(c_code->name,"invertebrate mt.");
	c_code->n_sense_c = 62;
	break;
      }
    case 5 :
      {
	c_code->num_curr_code = 5;
	strcpy(c_code->name,"ciliate nu.");
	c_code->n_sense_c = 63;
	break;
      }
    case 6 :
      {
	c_code->num_curr_code = 6;
	strcpy(c_code->name,"echinoderm nu.");
	c_code->n_sense_c = 62;
	break;
      }
    case 7 :
      {
	c_code->num_curr_code = 7;
	strcpy(c_code->name,"euplotid nu.");
	c_code->n_sense_c = 62;
	break;
      }
    case 8 :
      {
	c_code->num_curr_code = 8;
	strcpy(c_code->name,"bacterial & plant plastid");
	c_code->n_sense_c = 61;
	break;
      }
    case 9 :
      {
	c_code->num_curr_code = 9;
	strcpy(c_code->name,"alternative yeast nu.");
	c_code->n_sense_c = 61;
	break;
      }
    case 10 :
      {
	c_code->num_curr_code = 10;
	strcpy(c_code->name,"ascidian mt.");
	c_code->n_sense_c = 62;
	break;
      }
    case 11 :
      {
	c_code->num_curr_code = 11;
	strcpy(c_code->name,"alternative flatworm mt.");
	c_code->n_sense_c = 63;
	break;
      }
    case 12 :
      {
	c_code->num_curr_code = 12;
	strcpy(c_code->name,"blepharisma nu.");
	c_code->n_sense_c = 62;
	break;
      }
    case 13 :
      {
	c_code->num_curr_code = 13;
	strcpy(c_code->name,"chlorophycean mt.");
	c_code->n_sense_c = 62;
	break;
      }
    case 14 :
      {
	c_code->num_curr_code = 14;
	strcpy(c_code->name,"trematode mt.");
	c_code->n_sense_c = 62;
	break;
      }
    case 15 :
      {
	c_code->num_curr_code = 15;
	strcpy(c_code->name,"scenedesmus obliquus mt.");
	c_code->n_sense_c = 61;
	break;
      }
    case 16 :
      {
	c_code->num_curr_code = 16;
	strcpy(c_code->name,"thraustochytrium mt.");
	c_code->n_sense_c = 60;
	break;
      }
    case 17 :
      {
	c_code->num_curr_code = 17;
	strcpy(c_code->name,"Ziheng Yang's regular code");
	c_code->n_sense_c = 64;
	break;
      }
    default : {Exit("\n. Err. in Make_Code (unknown genetic code)\n");}
    }


  c_code->aa = (char *)mCalloc(23,sizeof(char));
  c_code->genetc = (int **)mCalloc(number_of_codes,sizeof(int *));
  For(i,c_code->number_of_codes) c_code->genetc[i] = (int *)mCalloc(64,sizeof(int));


  strcpy(c_code->aa,"ARNDCQEGHILKMFPSTWYV?-");
  
  For(i,number_of_codes) For(j,64) c_code->genetc[i][j] = genetc[i][j];
  
  
  Init_Codons(c_code);
  Init_Nb_Diff_Between_Two_Codons(c_code,c_code->n_sense_c);
  Init_Tstv_Table(c_code,c_code->n_sense_c);

  return c_code;
}

/*********************************************************/

void Init_Nb_Diff_Between_Two_Codons(code *c_code, int ns)
{
  /* Get the number of differences for each pair of codon */

  char *codon1,*codon2;
  int i,j,k;
  int diff;
  int num_codon1, num_codon2;

  codon1 = (char *)mCalloc(3,sizeof(char));
  codon2 = (char *)mCalloc(3,sizeof(char));

  num_codon1 = num_codon2 = -1;
  For(i,ns)
    {      
      num_codon1 = c_code->sense_c[i];

      codon1[0] = (num_codon1/16)%4; 
      codon1[1] = (num_codon1/4)%4; 
      codon1[2] = num_codon1%4;

      
      For(j,ns)
	{	  
	  num_codon2 = c_code->sense_c[j];

	  codon2[0] = (num_codon2/16)%4; 
	  codon2[1] = (num_codon2/4)%4; 
	  codon2[2] = num_codon2%4;	  
	  
	  diff = 0;
	  For(k,3) if(codon1[k] != codon2[k]) diff++;
	  c_code->n_diff_b_2_codons[i*ns+j] = diff;

	}      
    }

  Free(codon1);
  Free(codon2);

  return;
}

/*********************************************************/

void Init_Tstv_Table(code *c_code, int ns)
{
  /* If two codons differ by only one change, determine
     if this change is a transversion or a transition */


  int *codon1,*codon2;
  int i,j,k;
  int diff,pos;
  int num_codon1, num_codon2;

  codon1 = (int *)mCalloc(3,sizeof(int));
  codon2 = (int *)mCalloc(3,sizeof(int));

  num_codon1 = num_codon2 = -1;
  For(i,ns)
    {
      num_codon1 = c_code->sense_c[i];
            
      codon1[0] = (num_codon1/16)%4; 
      codon1[1] = (num_codon1/4)%4; 
      codon1[2] = num_codon1%4;

      
      For(j,ns)
	{
	  num_codon2 = c_code->sense_c[j];

	  codon2[0] = (num_codon2/16)%4; 
	  codon2[1] = (num_codon2/4)%4; 
	  codon2[2] = num_codon2%4;	  
	  
	  diff = 0;
	  pos = -1;
	  For(k,3) if(codon1[k] != codon2[k]) 
	    {
	      diff++;
	      pos = k;
	    }
	  if(diff > 1) c_code->tstvtable[i*ns+j] = -1;
	  else 
	    {
	      if(diff == 1)
		{
		  if((((codon1[pos] == 0) || (codon1[pos] == 2))  &&
		      ((codon2[pos] == 1) || (codon2[pos] == 3))) ||
		     (((codon1[pos] == 1) || (codon1[pos] == 3))  &&
		      ((codon2[pos] == 0) || (codon2[pos] == 2))))
		    c_code->tstvtable[i*ns+j] = 2; /* tv */
		  else
		    c_code->tstvtable[i*ns+j] = 1; /* ts */
		}
	      else  c_code->tstvtable[i*ns+j] = 0;	  
	    }
	}
    }

  Free(codon1);
  Free(codon2);
  return;
}

/*********************************************************/

void Init_GTR_4_Codon_Table(model *mod)
{
  int *codon1,*codon2;
  int i,j,k;
  int diff,pos;
  int num_codon1, num_codon2;
  int ns;

  codon1 = (int *)mCalloc(3,sizeof(int));
  codon2 = (int *)mCalloc(3,sizeof(int));

  ns = mod->c_code->n_sense_c;

  For(i,ns*ns) mod->c_code->gtr_4_codon[i] = NULL;

  num_codon1 = num_codon2 = -1;
  For(i,ns)
    {
      num_codon1 = mod->c_code->sense_c[i];
            
      codon1[0] = (num_codon1/16)%4; 
      codon1[1] = (num_codon1/4)%4; 
      codon1[2] = num_codon1%4;

      
      For(j,ns)
	{
	  num_codon2 = mod->c_code->sense_c[j];

	  codon2[0] = (num_codon2/16)%4; 
	  codon2[1] = (num_codon2/4)%4; 
	  codon2[2] = num_codon2%4;	  
	  



	  diff = 0;
	  pos = -1;
	  For(k,3) if(codon1[k] != codon2[k]) 
	    {
	      diff++;
	      pos = k;
	    }
	  if(diff > 1) mod->c_code->gtr_4_codon[i*ns+j] = NULL;
	  else 
	    {
                if(diff == 1)
                    {
                        switch(codon1[pos])
                            {
                            case 0 : /* A */
                                {
                                    switch(codon2[pos])
                                        {
                                        case 1 : {mod->c_code->gtr_4_codon[i*ns+j] = &(mod->gtr_param[0]); break;}/* A->C */
                                        case 2 : {mod->c_code->gtr_4_codon[i*ns+j] = &(mod->gtr_param[1]); break;}/* A->G */
                                        case 3 : {mod->c_code->gtr_4_codon[i*ns+j] = &(mod->gtr_param[2]); break;}/* A->T */
                                        default : {Exit("\n. Err. in Init_GTR_4_Codon_Table\n");}
                                        }
                                    break;
                                }
                                
                            case 1 : /* C */
                                {
                                    switch(codon2[pos])
                                        {
                                        case 0 : {mod->c_code->gtr_4_codon[i*ns+j] = &(mod->gtr_param[0]); break;}/* C->A */
                                        case 2 : {mod->c_code->gtr_4_codon[i*ns+j] = &(mod->gtr_param[3]); break;}/* C->G */
                                        case 3 : {mod->c_code->gtr_4_codon[i*ns+j] = &(mod->gtr_param[4]); break;}/* C->T */
                                        default : {Exit("\n. Err. in Init_GTR_4_Codon_Table\n");}
                                        }
                                    break;
                                }
                            case 2 : /* G */
                                {
                                    switch(codon2[pos])
                                        {
                                        case 0 : {mod->c_code->gtr_4_codon[i*ns+j] = &(mod->gtr_param[1]); break;}/* G->A */
                                        case 1 : {mod->c_code->gtr_4_codon[i*ns+j] = &(mod->gtr_param[3]); break;}/* G->C */
                                        case 3 : {mod->c_code->gtr_4_codon[i*ns+j] = &(mod->gtr_param[5]); break;}/* G->T */
                                        default : {Exit("\n. Err. in Init_GTR_4_Codon_Table\n");}
                                        }
                                    break;
                                }
                            case 3 : /* T */
                                {
                                    switch(codon2[pos])
                                        {
                                        case 0 : {mod->c_code->gtr_4_codon[i*ns+j] = &(mod->gtr_param[2]); break;}/* T->A */
                                        case 1 : {mod->c_code->gtr_4_codon[i*ns+j] = &(mod->gtr_param[4]); break;}/* T->C */
                                        case 2 : {mod->c_code->gtr_4_codon[i*ns+j] = &(mod->gtr_param[5]); break;}/* T->G */
                                        default : {Exit("\n. Err. in Init_GTR_4_Codon_Table\n");}
                                        }
                                    break;
                                }
                            default : Exit("\n. Err. in Init_GTR_4_Codon_Table\n");
                            }
                    }
	    }
	}
    }

  Free(codon1);
  Free(codon2);
  return;
}

/*********************************************************/

int Compare_Two_States(char *state1, char *state2, int state_size)
{
  /* 1 the two states are identical */
  /* 0 the two states are different */
  int i;

  For(i,state_size) if(state1[i] != state2[i]) break;
  
  return (i==state_size)?(1):(0);
}

/*********************************************************/

void Copy_One_State(char *from, char *to, int state_size)
{
  int i;
  For(i,state_size) to[i] = from[i];
}

/*********************************************************/

qmat *Make_Qmat_Struct(int ns, int n_qmat, int n_omega, model *mod)
{
  /* Set up the instantaneous rate matrix */

  qmat *qmat_struct;

  qmat_struct  = (qmat *)mCalloc(1,sizeof(qmat));

  qmat_struct->n_qmat  = n_qmat;
  qmat_struct->n_omega = n_omega;

  qmat_struct->u_mat             = (fit_double *) mCalloc(ns*ns*n_qmat,              sizeof(fit_double ));
  qmat_struct->v_mat             = (fit_double *) mCalloc(ns*ns*n_qmat,              sizeof(fit_double ));
  qmat_struct->root_vct          = (fit_double *) mCalloc(ns*n_qmat,                 sizeof(fit_double ));
  qmat_struct->expD_mr_vct       = (fit_double *) mCalloc(ns*n_qmat,                 sizeof(fit_double ));
  qmat_struct->qmat              = (fit_double *) mCalloc(ns*ns*n_qmat,              sizeof(fit_double )); 
  qmat_struct->theta             = (fit_double *) mCalloc((int)n_omega*(n_omega-1)/2,sizeof(fit_double )); 
  qmat_struct->qmat_proba        = (fit_double *) mCalloc(n_qmat,                    sizeof(fit_double )); 

  fit_double *p;
  
  p = (fit_double *) mCalloc(2*MAX(n_omega,n_qmat),sizeof(fit_double ));
  
  qmat_struct->omega = p;

  if(mod->switch_modelname == NO_SWITCH)
    {
      qmat_struct->trans_qmat_proba  = p + n_qmat;
      qmat_struct->trans_omega_proba = (fit_double *) mCalloc(n_omega,                   sizeof(fit_double )); 
    }
  else
    {
      qmat_struct->trans_omega_proba = p + n_omega;
      qmat_struct->trans_qmat_proba = (fit_double *) mCalloc(n_qmat,                   sizeof(fit_double )); 
    }

  qmat_struct->omega_min         = (fit_double *) mCalloc(n_omega,                   sizeof(fit_double )); 
  qmat_struct->omega_max         = (fit_double *) mCalloc(n_omega,                   sizeof(fit_double )); 
  qmat_struct->t_omega           = (fit_double **)mCalloc(n_omega*n_qmat,            sizeof(fit_double *)); 
  qmat_struct->omega_proba       = (fit_double *) mCalloc(n_omega,                   sizeof(fit_double )); 
  qmat_struct->t_omega_proba     = (fit_double **)mCalloc(n_omega*n_qmat,            sizeof(fit_double *)); 
  
  return qmat_struct;
}

/*********************************************************/

void Copy_Qmat_Struct(qmat *from, qmat *to)
{
  to = from;
}

/*********************************************************/

char *Get_Codon(int num_codon)
{
  char *s;
  s = (char *)mCalloc(3,sizeof(char));
  s[0] = Get_Base((num_codon/16)%4);
  s[1] = Get_Base((num_codon/4)%4);
  s[2] = Get_Base((num_codon)%4);
  return s;
}

/*********************************************************/

char Get_Base(int num_base)
{
  switch(num_base)
    {
    case 0 : return 'A'; break;
    case 1 : return 'C'; break;
    case 2 : return 'G'; break;
    case 3 : return 'T'; break;
    default : break;  
    }
  return -1;
}

/*********************************************************/

int Get_Base_Num(char base)
{
  switch(base)
    {
    case 'A' : return 0; break;
    case 'C' : return 1; break;
    case 'G' : return 2; break;
    case 'T' : return 3; break;
    default : break;  
    }
  return -1;
}

/*********************************************************/

void Init_Codons(code *c_code)
{
  int i,n;

  n = 0;
  For(i,64) 
    {
      if(c_code->genetc[c_code->num_curr_code][i] > -1) 
	{
	  c_code->from_64_2_61[i] = n;
	  c_code->sense_c[n] = i;
	  n++;
	}
      else
	  c_code->from_64_2_61[i] = -1;
    }
  if(n != c_code->n_sense_c) Exit("\n. Err. in Init_Codons !\n");
}

/*********************************************************/

char *Print_State(int num, arbre *tree)
{
  char *s;
  
  s = (char *)mCalloc(T_MAX_LINE,sizeof(char));

  Return_State(tree,num,s);

/*   printf("%s\n",s); */

  return s;
}

/*********************************************************/

void Infer_Synonymous_Substitutions(edge *b_fcus, node *a, node *d, arbre *tree)
{
  /* Put a label on the edge if the two codons at both extremities of this edge
     differ by at least one synonymous substitution */

  int i;
  int num_codon_a, num_codon_d;
  char *sa, *sd;

  b_fcus->substitution_type = 0;

  num_codon_a = tree->mod->c_code->sense_c[a->infered_state];
  num_codon_d = tree->mod->c_code->sense_c[d->infered_state];
  
  if((a->infered_state != d->infered_state) &&
     (Get_AA(tree->mod->c_code,num_codon_a) ==
      Get_AA(tree->mod->c_code,num_codon_d)))
    {
      b_fcus->check_this_one = 1;
      b_fcus->substitution_type = -1;
      sa = Print_State(a->infered_state,tree);
      sd = Print_State(d->infered_state,tree);
      strcat(b_fcus->labels[0],"s");
      strcat(b_fcus->labels[0],sa);
      strcat(b_fcus->labels[0],"<->");
      strcat(b_fcus->labels[0],sd);
      printf("%s <-> %s\n",sa,sd);
      Free(sa); Free(sd);
    }
      

  if(d->tax) return;
  else
    For(i,3) if(d->v[i] != a) 
      Infer_Synonymous_Substitutions(d->b[i],d,d->v[i],tree);

}

/*********************************************************/

void Infer_NonSynonymous_Substitutions(edge *b_fcus, node *a, node *d, arbre *tree)
{
  /* Put a label on the edge if the two codons at both extremities of this edge
     differ by at least one synonymous substitution */

  int i;
  int num_codon_a, num_codon_d;
  char *sa,*sd;

  b_fcus->substitution_type = 0;

  num_codon_a = tree->mod->c_code->sense_c[a->infered_state];
  num_codon_d = tree->mod->c_code->sense_c[d->infered_state];

  
  if((a->infered_state != d->infered_state) &&
     (Get_AA(tree->mod->c_code,num_codon_a) !=
      Get_AA(tree->mod->c_code,num_codon_d)))
    {
      b_fcus->check_this_one = 1;
      b_fcus->substitution_type = 1;
      sa = Print_State(a->infered_state,tree);
      sd = Print_State(d->infered_state,tree);
      strcat(b_fcus->labels[0],"ns");
      strcat(b_fcus->labels[0],sa);
      strcat(b_fcus->labels[0],"<->");
      strcat(b_fcus->labels[0],sd);
      printf("\n%s <-> %s",sa,sd);
      Free(sa); Free(sd);
    }
  
  if(d->tax) return;
  else
    For(i,3) if(d->v[i] != a) 
      Infer_NonSynonymous_Substitutions(d->b[i],d,d->v[i],tree);
}

/*********************************************************/

void Make_Sons(arbre *tree)
{
  /* Fill the lists of nodes that can be reached from a given
     node in each of the three directions */
  
  int i,j;

  For(i,2*tree->n_otu-2)
    {
      tree->noeud[i]->n_sons = (int *)mCalloc(3,sizeof(int));;
      tree->noeud[i]->sons = (node ***)mCalloc(3,sizeof(node **));
      For(j,3)
	tree->noeud[i]->sons[j] = (node **)mCalloc(2*tree->n_otu-2,sizeof(node *));
    }

  Make_Sons_Post(tree->noeud[0],tree->noeud[0]->v[0],tree);
  Make_Sons_Pre(tree->noeud[0],tree->noeud[0]->v[0],tree);

}

/*********************************************************/

void Make_Sons_Post(node *a, node *d, arbre *tree)
{
  /* Post-oder traversal for Make_Sons */

  int i;

  if (d->tax) return;
  else
    {
      int j;
      int dir_d_a;
      node **ddsons;
      int ddnsons;

      dir_d_a = -1;
      For(i,3)
	{
	  if(d->v[i] != a)
	    Make_Sons_Post(d,d->v[i],tree);
	  else
	    dir_d_a = i;
	}

      For(i,3)
	if(d->v[i] != a)
	  {
	    ddsons = (d == d->b[i]->left)?
	      (d->b[i]->rght->sons[d->b[i]->r_l]):
	      (d->b[i]->left->sons[d->b[i]->l_r]);

	    ddnsons = (d == d->b[i]->left)?
	      (d->b[i]->rght->n_sons[d->b[i]->r_l]):
	      (d->b[i]->left->n_sons[d->b[i]->l_r]);
	    
	    For(j,ddnsons)
	      {
		d->sons[dir_d_a][d->n_sons[dir_d_a]] = ddsons[j];
		d->n_sons[dir_d_a] += 1;
	      }
	    d->sons[dir_d_a][d->n_sons[dir_d_a]] = d->v[i];
	    d->n_sons[dir_d_a] += 1;
	  }
    }
}

/*********************************************************/

void Make_Sons_Pre(node *a, node *d, arbre *tree)
{
  
  /* Pre-oder traversal for Make_Sons */

  if(d->tax) return;
  else
    {
      int i,j,k;
      int ddnsons;
      node **ddsons;

      For(i,3)
	{
	  if(d->v[i] != a)
	    {
	      For(j,3)
		{
		  if(j != i)
		    {
		      ddsons = (d == d->b[j]->left)?
			(d->b[j]->rght->sons[d->b[j]->r_l]):
			(d->b[j]->left->sons[d->b[j]->l_r]);
		      
		      ddnsons = (d == d->b[j]->left)?
			(d->b[j]->rght->n_sons[d->b[j]->r_l]):
			(d->b[j]->left->n_sons[d->b[j]->l_r]);
		      
		      For(k,ddnsons)
			{
			  d->sons[i][d->n_sons[i]] = ddsons[k];
			  d->n_sons[i] += 1;
			}
		      d->sons[i][d->n_sons[i]] = d->v[j];
		      d->n_sons[i] += 1;
		    }
		}
	    }
	}

      For(i,3)
	if(d->v[i] != a)
	  Make_Sons_Pre(d,d->v[i],tree);
    }
}

/*********************************************************/

int Is_Son(node *a, node *x, int dir)
{
  /* Is 'x' a son of 'a' in direction 'dir' ? */

  int i;
  For(i,a->n_sons[dir]) 
    if(a->sons[dir][i] == x) 
      return 1;
  return 0;
}

/*********************************************************/

void Distance_Between_2_nodes(edge *b_fcus, node *a, node *d, node *final, fit_double *distance)
{  
  /* Distance (expected number of substitutions per site) between two nodes in the tree */

  if(a == final) return;
  if(d == final) 
    {
      (*distance)+=b_fcus->l;
      return;
    }
  else
    {
      int i,j;

      For(i,3)
	{
	  if((d->v[i] == a) && (Is_Son(d,final,i)))
	    {
	      (*distance)+=b_fcus->l;
	      For(j,3)
		{
		  if(d->v[j] != a)
		    Distance_Between_2_nodes(d->b[j],d,d->v[j],final,distance);
		}
	    }
	}
    }
}

/*********************************************************/

void qksort(fit_double* A, int ilo, int ihi)
{
  fit_double pivot;       // pivot value for partitioning array
  int ulo, uhi;       // indices at ends of unpartitioned region
  int ieq;            // least index of array entry with value equal to pivot
  fit_double tempEntry;   // temporary entry used for swapping
                                                                                                                                                                                                   
    if (ilo >= ihi) {
        return;
    }
    // Select a pivot value.
    pivot = A[(ilo + ihi)/2];
    // Initialize ends of unpartitioned region and least index of entry
    // with value equal to pivot.
    ieq = ulo = ilo;
    uhi = ihi;
    // While the unpartitioned region is not empty, try to reduce its size.
    while (ulo <= uhi) {
        if (A[uhi] > pivot) {
	  // Here, we can reduce the size of the unpartitioned region and
	  // try again.
            uhi--;
        } else {
	  // Here, A[uhi] <= pivot, so swap entries at indices ulo and
	  // uhi.
            tempEntry = A[ulo];
            A[ulo] = A[uhi];
            A[uhi] = tempEntry;
            // After the swap, A[ulo] <= pivot.
            if (A[ulo] < pivot) {
	      // Swap entries at indices ieq and ulo.
                tempEntry = A[ieq];
                A[ieq] = A[ulo];
                A[ulo] = tempEntry;
                // After the swap, A[ieq] < pivot, so we need to change
                // ieq.
                ieq++;
                // We also need to change ulo, but we also need to do
                // that when A[ulo] = pivot, so we do it after this if
                // statement.
            }
            // Once again, we can reduce the size of the unpartitioned
            // region and try again.
            ulo++;
        }
    }
    // Now, all entries from index ilo to ieq - 1 are less than the pivot
    // and all entries from index uhi to ihi + 1 are greater than the
    // pivot.  So we have two regions of the array that can be sorted
    // recursively to put all of the entries in order.
    qksort(A, ilo, ieq - 1);
    qksort(A, uhi + 1, ihi);
}

/*********************************************************/

void Update_Tpos_Min(arbre *tree)
{
  /* Only used with serial sample data */

  int i;
  node *root;

  root = tree->n_root;

  For(i,3)
    {
      if(root->v[i])
	{
	  Update_Tpos_Min_Post(root,root->v[i]);
	  Update_Tpos_Min_Pre(root,root->v[i]);
	}
    }
}

/*********************************************************/

void Update_Tpos_Min_Post(node *a, node *d)
{ 
  /* Only used with serial sample data */

 if(d->tax) 
    {
      /* tpos never changes at a tip node */
      d->tpos_min[0] = d->tpos;
      d->tpos_lim_below = d->tpos;
      d->tpos_lim_above = d->tpos;
      return;
    }
  else
    {
      int i,j;
      int d_v1,d_v2,d_a;
      fit_double min1, min2;
      fit_double tpos_below1, tpos_below2;


      d_v1 = d_v2 = d_a = -1;
      For(i,3)
	{
	  if(d->v[i] != a)
	    {
	      Update_Tpos_Min_Post(d,d->v[i]);
	      if(d_v1 < 0)  d_v1 = i;
	      else          d_v2 = i;
	    }
	  else d_a = i;
	}


      min1 = min2 = tpos_below1 = tpos_below2 = -1.;
      For(i,3)
	{
	  if(d->v[i] != a)
	    {
	      For(j,3) if(d->v[i]->v[j] == d) 
		{
		  if(min1 < .0) 
		    {
		      min1 = d->v[i]->tpos_min[j];
		      tpos_below1 = d->v[i]->tpos;
		    }
		  else          
		    {
		      min2 = d->v[i]->tpos_min[j];
		      tpos_below2 = d->v[i]->tpos;
		    }
		}
	    }
	}

      d->tpos_min[d_a] = MMIN(min1,min2);      
      d->tpos_lim_below = MMIN(tpos_below1, tpos_below2);
      d->tpos_lim_above = a->tpos;


      if(d->tpos_lim_above > d->tpos_lim_below) 
	{
	  printf("\nnode %d -> %f ; below %f ; above (%d) %f\n",
		 d->num,
		 d->tpos,
		 d->tpos_lim_below,
		 a->num,
		 d->tpos_lim_above);
	  printf("\n. Warning in Update_Tpos_Min_Post\n");
	}
    }
}

/*********************************************************/

void Update_Tpos_Min_Pre(node *a, node *d)
{
  /* Only used with serial sample data */

  if(d->tax) return;
  else
    {
      int i,j,k;
      fit_double min1, min2;

      For(i,3)
	{
	  if(d->v[i] != a)
	    {
	      min1 = min2 = -1.;
	      For(j,3)
		{
		  if(d->v[j] != d->v[i])
		    {
		      For(k,3)
			{
			  if(d->v[j]->v[k] == d)
			    {
			      if(min1 < .0) min1 = d->v[j]->tpos_min[k];
			      else          min2 = d->v[j]->tpos_min[k];
				
			    }
			}
		    }
		}
	      d->tpos_min[i] = MMIN(min1,min2);
	    }	  
	}
      
      For(i,3)
	if(d->v[i] != a)
	  Update_Tpos_Min_Pre(d,d->v[i]);
    }
}

/*********************************************************/

void Tpos_OLS(arbre *tree, node *a, node *d)
{
  /* Only used with serial sample data */

  if(d->tax) return;
  else
    {
      int i,d_v1,d_v2;

      d_v1 = d_v2 = -1;
      For(i,3)
	{
	  if(d->v[i] != a)
	    {
	      Tpos_OLS(tree,d,d->v[i]);
	      if(d_v1 < 0) d_v1 = i;
	      else         d_v2 = i;
	    }
	}
      
      d->tpos = 
	(d->v[d_v1]->tpos*tree->mod->subst_rate - d->b[d_v1]->l_ml + 
	 d->v[d_v2]->tpos*tree->mod->subst_rate - d->b[d_v2]->l_ml)/(2.*tree->mod->subst_rate);

/*       printf("Node %d -> %f (%f %f) (%f %f)\n", */
/* 	     d->num,d->tpos, */
/* 	     d->v[d_v1]->tpos,d->v[d_v2]->tpos, */
/* 	     d->b[d_v1]->l_ml,d->b[d_v2]->l_ml); */

      if(d->tpos > MMIN(d->v[d_v1]->tpos,d->v[d_v2]->tpos))
	d->tpos = MMIN(d->v[d_v1]->tpos,d->v[d_v2]->tpos);

    }
}

/*********************************************************/

void Update_Bl_Using_Tpos(arbre *tree, node *a, node *d)
{
  /* Only used with serial sample data */

  if(d->tax) return;
  else
    {
      int i;

      For(i,3)
	if(d->v[i] != a)
	  Update_Bl_Using_Tpos(tree,d,d->v[i]);
      
      if(tree->n_root->tpos > a->tpos) tree->n_root->tpos = a->tpos;

      For(i,3)
	if(d->v[i] != a)
	  {
	    d->b[i]->l = (d->v[i]->tpos - d->tpos)*tree->mod->subst_rate;
	    if(d->v[i]->tpos < (d->tpos - 1.E-5)) 
	      {
		printf("%.10f %.10f\n",d->v[i]->tpos,d->tpos);
		Exit("\n. Err in Update_Bl_Using_Tpos\n");
	      }
	  }
    }
}
/*********************************************************/

void Set_Defaults_Model(model *mod)
{
  /* Default values of substitution model parameters */
  int i;

  For(i,MMAX(N_MAX_CATQ,N_MAX_OMEGA)) mod->omega_min[i] = 0.0;
  For(i,MMAX(N_MAX_CATQ,N_MAX_OMEGA)) mod->omega_max[i] = 20.;
  
  mod->ns                   = 4;
  mod->tpos_ols             = 0;
  mod->update_bl_using_tpos = 0;
  mod->subst_rate           = .0;
  mod->model_number         = 4;
  mod->n_catg               = 1;
  mod->n_catq               = 1;
  mod->kappa                = 4.0;
  mod->alpha                = 2.0;
  mod->lambda               = 1.0;
  mod->pinvar               = 0.0;
  mod->s_opt->opt_alpha     = 0;
  mod->s_opt->opt_kappa     = 0;
  mod->s_opt->opt_lambda    = 0;
  mod->s_opt->opt_bl        = 0;
  mod->s_opt->opt_pinvar    = 0;
  mod->s_opt->opt_omega     = 1;
  mod->s_opt->opt_theta     = 1;
  mod->invar                = 0;
  mod->stepsize             = 1;
  mod->codon_freq           = 1;
  mod->omega[0]             = 0.0;
  mod->omega[1]             = 1.0;
  mod->omega[2]             = 4.0;
  mod->omega_proba[0]       = 0.6;
  mod->omega_proba[1]       = 0.3;
  mod->omega_proba[2]       = 0.1;
  mod->delta_switch         = 0.0;
  mod->io_delta_switch      = 0.0;
  mod->selec_reg_freq[0]    = -1.;
  mod->selec_reg_freq[1]    = -1.;
  mod->selec_reg_freq[2]    = -1.;
  mod->n_catq_negsel        = -1;
  mod->n_catq_possel        = -1;
  mod->thresholdp2          = 0.95;
  mod->fdr                  = 0.05;
  mod->datatype             = NT;
  mod->subst_modelname      = HKY85;
  mod->switch_modelname     = NO_SWITCH;
  mod->model_number         = 4;
  mod->model_applies_to     = NT;
  mod->n_rsubst_rate        = 2;
  
}

/*********************************************************/


void Init_NNI(node *a, node *d, edge *b_fcus, arbre *tree)
{
  /* Deprecated */

  int i;
  int first;
  int d_v1, d_v2;


  b_fcus->s_nni->bl_info_init = -1.;

  b_fcus->s_nni->bl_info[0]   = d->tpos;
  b_fcus->s_nni->bl_info[1]   = -1.;
  b_fcus->s_nni->bl_info[2]   = -1.;

  b_fcus->s_nni->best_conf    = 1;

  b_fcus->s_nni->b_fcus       = b_fcus;

  b_fcus->s_nni->v2           = a;
  b_fcus->s_nni->v4           = d;
  

  first = 0; 
  d_v1 = d_v2 = -1;
  For(i,3) 
    if(d->v[i] != a)
      {
	if(!first) 
	  {
	    d_v1 = i;
	    b_fcus->s_nni->v5 = d->v[i];
	    first = 1;
	  }
	else
	  {
	    d_v2 = i;
	    b_fcus->s_nni->v6 = d->v[i];
	  }
      }
  
  


  if(d->tax) return;
  else
    {
      d->b[d_v1]->s_nni->v1 = a;
      d->b[d_v1]->s_nni->v3 = d->v[d_v2];
      
      d->b[d_v2]->s_nni->v1 = a;
      d->b[d_v2]->s_nni->v3 = d->v[d_v1];
      
      For(i,3)
	if(d->v[i] != a) Init_NNI(d,d->v[i],d->b[i],tree);
    }
}


/*********************************************************/

void Check_Tpos(node *a, node *d)
{
  /* Only used with serial sample data */

  if(d->tpos < a->tpos) 
    {
      printf("\n.tpos incompatibility repaired (shift=%f)",
	     (d->tpos-a->tpos)/2.);
      a->tpos = a->tpos - (fit_double)fabs(d->tpos-a->tpos)/2.;
      d->tpos = a->tpos;
      fflush(NULL);
    }
  if(d->tax) return;
  else
    {
      int i;
      For(i,3)
	if(d->v[i] != a) Check_Tpos(d,d->v[i]);
    }
}

/*********************************************************/

void Update_Ancestors(node *a, node *d)
{
  d->anc = a;
  if(d->tax) return;
  else
    {
      int i;
      For(i,3) if(d->v[i] != d->anc) Update_Ancestors(d,d->v[i]);
    }
}

/*********************************************************/

void Add_Root(edge *target, arbre *tree)
{
/*   printf("\n. Add root on edge %d left = %d right = %d",target->num,target->left->num,target->rght->num); fflush(NULL); */
  tree->e_root = target;

  /* Create the root node if it does not exist yet */
  if((!tree->n_root) || (tree->n_root->num != 2*tree->n_otu-2))
    {      
      tree->n_root = (node *)Make_Node_Light(2*tree->n_otu-2);
    }

  tree->n_root->tax = 0;

  /* Set the position of the root */
  tree->n_root->v[0] = tree->e_root->left;
  tree->n_root->v[1] = tree->e_root->rght;

  tree->n_root->b[0] = tree->e_root;
  tree->n_root->b[1] = tree->e_root;

  if(tree->n_root_pos > -1.0)
    {
      if(tree->n_root_pos < 1.E-6 &&  tree->n_root_pos > -1.E-6)
	printf("\n. WARNING: you put the root at a weird position...");

/*       tree->n_root->l[0] = tree->e_root->l * (tree->n_root_pos/(1.+tree->n_root_pos)); */
/*       tree->n_root->l[1] = tree->e_root->l - tree->n_root->l[0]; */

      tree->n_root->l[0] = tree->e_root->l * tree->n_root_pos;
      tree->n_root->l[1] = tree->e_root->l * (1. - tree->n_root_pos);
    }
  else
    {
      tree->n_root->l[0] = tree->e_root->l / 2.;
      tree->n_root->l[1] = tree->e_root->l / 2.;
      tree->n_root_pos = 0.5;
    }
  
  Update_Ancestors(tree->n_root,tree->n_root->v[0]);
  Update_Ancestors(tree->n_root,tree->n_root->v[1]);
  tree->n_root->anc = NULL;
}

/*********************************************************/

void Remove_Root(arbre *tree)
{
  /* Remove the root of a tree */

  int i;

  For(i,3)
    {
      if(tree->n_root->v[0]->v[1]->v[i] == tree->n_root->v[0])
	tree->n_root->v[0]->v[1]->v[i] = tree->n_root->v[0]->v[2];

      if(tree->n_root->v[0]->v[2]->v[i] == tree->n_root->v[0])
	tree->n_root->v[0]->v[2]->v[i] = tree->n_root->v[0]->v[1];	  
    }

  (tree->n_root->v[0]->b[1]->left == tree->n_root->v[0]->v[1])?
    (tree->n_root->v[0]->b[1]->rght = tree->n_root->v[0]->v[2]):
    (tree->n_root->v[0]->b[1]->left = tree->n_root->v[0]->v[2]);
}

/*********************************************************/

model *Make_Model_Basic()
{
  /* Start the allocation of memory for the substitution model */

  model *mod;

  mod                 = (model *)mCalloc(1,sizeof(model));  
  mod->s_opt          = (optimiz *)Alloc_Optimiz();
  mod->gtr_param      = (fit_double *)mCalloc(6,sizeof(fit_double));
  mod->omega          = (fit_double *)mCalloc(MMAX(N_MAX_CATQ,N_MAX_OMEGA),sizeof(fit_double));
  mod->omega_proba    = (fit_double *)mCalloc(MMAX(N_MAX_CATQ,N_MAX_OMEGA),sizeof(fit_double));
  mod->omega_min      = (fit_double *)mCalloc(MMAX(N_MAX_CATQ,N_MAX_OMEGA),sizeof(fit_double));
  mod->omega_max      = (fit_double *)mCalloc(MMAX(N_MAX_CATQ,N_MAX_OMEGA),sizeof(fit_double));
  mod->selec_reg_freq = (fit_double *)mCalloc(MMAX(N_MAX_CATQ,N_MAX_OMEGA),sizeof(fit_double));

  return mod;
}

/*********************************************************/

void Make_Model_Complete(model *mod)
{
  /* Continue the allocation of memory for the substitution model... */
  mod->pi             = (fit_double *)mCalloc(mod->ns+1,                     sizeof(fit_double));
  mod->r_proba        = (fit_double *)mCalloc(mod->n_catg,                   sizeof(fit_double));
  mod->rr_mixturem    = (fit_double *)mCalloc(mod->n_catg,                   sizeof(fit_double));
  mod->rsubst_rate    = (fit_double *)mCalloc(mod->n_rsubst_rate,            sizeof(fit_double));
}

/*********************************************************/

void Make_All_Qmat_Struct(model *mod)
{
  /* Allocate memory of the instantaneous rate matrices structure */

  int i;

  if((mod->model_applies_to == NT) || (mod->model_applies_to == AA))
    {
      mod->qmat_struct    = (qmat **)mCalloc(2*mod->n_otu-3,sizeof(qmat *));
      mod->qmat_struct[0] = Make_Qmat_Struct(mod->ns,1,1,mod);		
      for(i=1;i<2*mod->n_otu-3;i++) mod->qmat_struct[i] = mod->qmat_struct[0];
      mod->qmat_struct[0]->qmat_proba[0] = 1.0;
      mod->qmat_struct[0]->curr_qmat_cat = 0;
    }
  else if(mod->switch_modelname == NO_SWITCH)
    {
      mod->qmat_struct = (qmat **)mCalloc(2*mod->n_otu-3,sizeof(qmat *));
      mod->qmat_struct[0] = Make_Qmat_Struct(mod->ns,mod->n_catq,mod->n_omega,mod);		
      for(i=1;i<2*mod->n_otu-3;i++) mod->qmat_struct[i] =  mod->qmat_struct[0];
    }
  else if(mod->switch_modelname != NO_SWITCH)
    {      
      mod->qmat_struct = (qmat **)mCalloc(2*mod->n_otu-3,sizeof(qmat *));
      mod->qmat_struct[0] = Make_Qmat_Struct(mod->ns,1,mod->n_omega,mod);		
      for(i=1;i<2*mod->n_otu-3;i++) mod->qmat_struct[i] =  mod->qmat_struct[0];
    }
  else
    {
      printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Exit("");
    }
}

/*********************************************************/

void Print_Param(FILE *fp, arbre *tree)
{
  /* Print the values of the substitution model parameters */
  int i,j;

  /* fprintf(fp,". "); */

  /* if(tree->mod->s_opt->opt_kappa)  fprintf(fp,"kappa=%f "  ,tree->mod->kappa); */
  /* if(tree->mod->s_opt->opt_alpha)  fprintf(fp,"alpha=%f "  ,tree->mod->alpha); */
  /* if(tree->mod->s_opt->opt_lambda) fprintf(fp,"lambda=%f " ,tree->mod->lambda); */
  /* if(tree->mod->s_opt->opt_pinvar) fprintf(fp,"pinv=%f "   ,tree->mod->pinvar); */
  
  /* if((tree->mod->model_applies_to == CODONS) && (tree->mod->switch_modelname == NO_SWITCH)) */
  /*     { */
  /*       For(i,tree->mod->n_omega) */
  /*         fprintf(fp,"p%d=%f ",i,tree->mod->qmat_struct[0]->qmat_proba[i]); */

  /*       For(i,tree->mod->n_omega) */
  /*         fprintf(fp,"w%d=%f ",i,tree->mod->qmat_struct[0]->omega[i]); */
  /*     } */
  /* else if((tree->mod->model_applies_to == CODONS) && */
  /*         ((tree->mod->switch_modelname == SWITCH_S1) ||(tree->mod->switch_modelname == SWITCH_S2))) */
  /*     { */

  /*       For(i,tree->mod->n_omega) */
  /*         fprintf(fp,"p%d=%f ",i,tree->mod->qmat_struct[0]->omega_proba[i]); */

  /*       For(i,tree->mod->n_omega) */
  /*         fprintf(fp,"w%d=%f ",i,tree->mod->qmat_struct[0]->omega[i]); */

  /*       fprintf(fp,"delta=%f ",tree->mod->qmat_struct[0]->theta[0]); */
  /*       fprintf(fp,"alpha=%f ",tree->mod->qmat_struct[0]->theta[1]/tree->mod->qmat_struct[0]->theta[0]); */
  /*       fprintf(fp,"beta=%f ",tree->mod->qmat_struct[0]->theta[2]/tree->mod->qmat_struct[0]->theta[0]); */
  /*       fprintf(fp,"meanswitch=%f ",tree->mod->delta_switch); */

  /*       For(i,tree->mod->n_omega-1) */
  /*         for(j=i+1;j<tree->mod->n_omega;j++) */
  /*           fprintf(fp,"R%d%d=%f ", */
  /*       	   i,j, */
  /*       	   tree->mod->qmat_struct[0]->theta[MMIN(i,j) * tree->mod->n_omega + MMAX(i,j) - */
  /*       					   (MMIN(i,j)+1+(int)pow(MMIN(i,j)+1,2))/2]/(tree->mod->delta_switch)); */
  /*     } */

  /* fprintf(fp,"[%d] ",tree->numerical_param->replace_num_param); */
  /* fprintf(fp,"[lnL=%f]\n",tree->tot_loglk); */
  /* fflush(NULL); */


  if(tree->mod->s_opt->opt_kappa)  fprintf(fp,"kappa=%4.2f "  ,tree->mod->kappa);
  if(tree->mod->s_opt->opt_alpha)  fprintf(fp,"alpha=%4.2f "  ,tree->mod->alpha);
  if(tree->mod->s_opt->opt_lambda) fprintf(fp,"lambda=%4.2f " ,tree->mod->lambda);
  if(tree->mod->s_opt->opt_pinvar) fprintf(fp,"pinv=%4.2f "   ,tree->mod->pinvar);
  
  if((tree->mod->model_applies_to == CODONS) && (tree->mod->switch_modelname == NO_SWITCH))
      {
        For(i,tree->mod->n_omega)
          fprintf(fp,"p%d=%3.4f ",i,tree->mod->qmat_struct[0]->qmat_proba[i]);

        For(i,tree->mod->n_omega)
          fprintf(fp,"w%d=%5.4f ",i,tree->mod->qmat_struct[0]->omega[i]);
      }
  else if((tree->mod->model_applies_to == CODONS) &&
          ((tree->mod->switch_modelname == SWITCH_S1) ||(tree->mod->switch_modelname == SWITCH_S2)))
      {

        For(i,tree->mod->n_omega)
          fprintf(fp,"p%d=%3.4f ",i,tree->mod->qmat_struct[0]->omega_proba[i]);

        For(i,tree->mod->n_omega)
          fprintf(fp,"w%d=%5.4f ",i,tree->mod->qmat_struct[0]->omega[i]);

        fprintf(fp,"delta=%5.2f ",tree->mod->qmat_struct[0]->theta[0]);
        fprintf(fp,"alpha=%5.2f ",tree->mod->qmat_struct[0]->theta[1]/tree->mod->qmat_struct[0]->theta[0]);
        fprintf(fp,"beta=%5.2f ",tree->mod->qmat_struct[0]->theta[2]/tree->mod->qmat_struct[0]->theta[0]);
        fprintf(fp,"meanswitch=%5.2f ",tree->mod->delta_switch);

        For(i,tree->mod->n_omega-1)
          for(j=i+1;j<tree->mod->n_omega;j++)
            fprintf(fp,"R%d%d=%5.2f ",
        	   i,j,
        	   tree->mod->qmat_struct[0]->theta[MMIN(i,j) * tree->mod->n_omega + MMAX(i,j) -
        					   (MMIN(i,j)+1+(int)pow(MMIN(i,j)+1,2))/2]/(tree->mod->delta_switch));
      }

  fprintf(fp,"[%d] ",tree->numerical_param->replace_num_param);
  fprintf(fp,"[lnL=%f]\n",tree->tot_loglk);
  fflush(NULL);
}

/*********************************************************/

void Check_Memory_Amount(arbre *tree)
{
  /* Rough estimate of the amount of memory that has to be used */

  int nbytes;
  model *mod;

  mod = tree->mod;

  nbytes = 0;


  /* Pmat */
  nbytes += (2*mod->n_otu-3) * mod->n_catg * mod->n_catq * mod->ns * mod->ns * sizeof(fit_double);
  nbytes += (2*mod->n_otu-3) * mod->n_catg * mod->n_catq * mod->ns * sizeof(fit_double *);
  nbytes += (2*mod->n_otu-3) * mod->n_catg * mod->n_catq * sizeof(fit_double **);
  nbytes += (2*mod->n_otu-3) * mod->n_catg * sizeof(fit_double ***);
  
  /* Partial Lk */
  nbytes += (2*mod->n_otu-3) * 2 * tree->n_pattern * mod->n_catq * mod->n_catg * mod->ns * sizeof(fit_double);
  nbytes += (2*mod->n_otu-3) * 2 * tree->n_pattern * mod->n_catq * mod->n_catg * sizeof(fit_double *);
  nbytes += (2*mod->n_otu-3) * 2 * tree->n_pattern * mod->n_catq * sizeof(fit_double **);
  nbytes += (2*mod->n_otu-3) * 2 * tree->n_pattern * sizeof(fit_double ***);
  nbytes += (2*mod->n_otu-3) * 2 * tree->n_pattern * sizeof(fit_double);

  if(((fit_double)nbytes/(1.E+06)) > 256.)
    {
      char answer;
      printf("\n. WARNING: this analysis requires at least %.0fMo of memory space.\n",(fit_double)nbytes/(1.E+06));
      printf("  Do you really want to proceed ? [Y/n] ");
      if(!scanf("%c", &answer))
	{
	  printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	  Exit("");
	}
      if(answer == '\n') answer = 'Y';
      else if(answer == 'n' || answer == 'N') Exit("\n\n");
      else getchar();
    }
  else if(((fit_double)nbytes/(1.E+06)) > 100.)
      printf("\n. WARNING: this analysis will use at least %.0fMo of memory space...\n",(fit_double)nbytes/(1.E+06));
      
}

/*********************************************************/

int Are_Compatible(char *statea, char *stateb, int stepsize, int datatype)
{
  int i,j;
  char a,b;


  if(datatype == NT) 
    {
      For(i,stepsize)
	{
	  a = statea[i];
	  For(j,stepsize)
	    {
	      b = stateb[j];

	      switch(a)
		{
		case 'A':
		  {
		    switch(b)
		      {
		      case 'A' : 
		      case 'M' : 
		      case 'R' : 
		      case 'W' : 
		      case 'D' : 
		      case 'H' : 
		      case 'V' : 
		      case 'X' : {b=b; break;}
		      default : return 0;
		      }
		    break;
		  }
		case 'G':
		  {
		    switch(b)
		      {
		      case 'G' : 
		      case 'R' : 
		      case 'S' : 
		      case 'K' : 
		      case 'B' : 
		      case 'D' : 
		      case 'V' : 
		      case 'X' : {b=b; break;}
		      default : return 0;
		      }
		    break;
		  }
		case 'C':
		  {
		    switch(b)
		      {
		      case 'C' : 
		      case 'M' : 
		      case 'S' : 
		      case 'Y' : 
		      case 'B' : 
		      case 'H' : 
		      case 'V' : 
		      case 'X' : {b=b; break;}
		      default : return 0;
		      }
		    break;
		  }
		case 'T':
		  {
		    switch(b)
		      {
		      case 'T' : 
		      case 'W' : 
		      case 'Y' : 
		      case 'K' : 
		      case 'B' : 
		      case 'D' : 
		      case 'H' : 
		      case 'X' : 
			{b=b; break;}
		      default : return 0;
		      }
		    break;
		  }
		case 'M' : 
		  {
		    switch(b)
		      {
		      case 'M' : 
		      case 'A' :
		      case 'C' :
		      case 'R' : 
		      case 'W' : 
		      case 'S' : 
		      case 'Y' : 
		      case 'B' : 
		      case 'D' : 
		      case 'H' : 
		      case 'V' : 
		      case 'X' :
			{b=b; break;}
		      default : return 0;
		      }	
		    break;
		  }
		case 'R' :
		  {
		    switch(b)
		      {
		      case 'R' :
		      case 'A' :
		      case 'G' :
		      case 'M' :
		      case 'W' :
		      case 'S' :
		      case 'K' :
		      case 'B' :
		      case 'D' :
		      case 'H' :
		      case 'V' :
		      case 'X' : {b=b; break;}
		      default : return 0;
		      }
		    break;
		  }
		  
		case 'W' :
		  {
		    switch(b)
		      {
		      case 'W' :
		      case 'A' :
		      case 'T' :
		      case 'M' :
		      case 'R' :
		      case 'Y' :
		      case 'K' :
		      case 'B' :
		      case 'D' :
		      case 'H' :
		      case 'V' :
		      case 'X' : {b=b; break;}
		      default : return 0;
		      }
		    break;
		  }
		  
		case 'S' :
		  {
		    switch(b)
		      {
		      case 'S' :
		      case 'C' :
		      case 'G' :
		      case 'M' :
		      case 'R' :
		      case 'Y' :
		      case 'K' :
		      case 'B' :
		      case 'D' :
		      case 'H' :
		      case 'V' :
		      case 'X' : {b=b; break;}
		      default : return 0;
		      }
		    break;
		  }
		  
		case 'Y' :
		  {
		    switch(b)
		      {
		      case 'Y' :
		      case 'C' :
		      case 'T' :
		      case 'M' :
		      case 'W' :
		      case 'S' :
		      case 'K' :
		      case 'B' :
		      case 'D' :
		      case 'H' :
		      case 'V' :
		      case 'X' : {b=b; break;}
		      default : return 0;
		      }
		    break;
		  }
		  
		case 'K' :
		  {
		    switch(b)
		      {
		      case 'K' :
		      case 'G' :
		      case 'T' :
		      case 'R' :
		      case 'W' :
		      case 'S' :
		      case 'Y' :
		      case 'B' :
		      case 'D' :
		      case 'H' :
		      case 'V' :
		      case 'X' : {b=b; break;}
		      default : return 0;
		      }
		    break;
		  }
		case 'B' :
		  {
		    switch(b)
		      {
		      case 'B' :
		      case 'C' :
		      case 'G' :
		      case 'T' :
		      case 'M' :
		      case 'R' :
		      case 'W' :
		      case 'S' :
		      case 'Y' :
		      case 'K' :
		      case 'D' :
		      case 'H' :
		      case 'V' :
		      case 'X' : {b=b; break;}
		      default : return 0;
		      }
		    break;
		  }
		case 'D' :
		  {
		    switch(b)
		      {
		      case 'D' :
		      case 'A' :
		      case 'G' :
		      case 'T' :
		      case 'M' :
		      case 'R' :
		      case 'W' :
		      case 'S' :
		      case 'Y' :
		      case 'K' :
		      case 'B' :
		      case 'H' :
		      case 'V' :
		      case 'X' : {b=b; break;}
		      default : return 0;
		      }
		    break;
		  }
		case 'H' :
		  {
		    switch(b)
		      {
		      case 'H' :
		      case 'A' :
		      case 'C' :
		      case 'T' :
		      case 'M' :
		      case 'R' :
		      case 'W' :
		      case 'S' :
		      case 'Y' :
		      case 'K' :
		      case 'B' :
		      case 'D' :
		      case 'V' :
		      case 'X' : {b=b; break;}
		      default : return 0;
		      }
		    break;
		  }
		case 'V' :
		  {
		    switch(b)
		      {
		      case 'V' :
		      case 'A' :
		      case 'C' :
		      case 'G' :
		      case 'M' :
		      case 'R' :
		      case 'W' :
		      case 'S' :
		      case 'Y' :
		      case 'K' :
		      case 'B' :
		      case 'D' :
		      case 'H' :
		      case 'X' : {b=b; break;}
		      default : return 0;
		      }
		    break;
		  }
		case 'X' :
		  {
		    switch(b)
		      {
		      case 'X' :
		      case 'A' :
		      case 'C' :
		      case 'G' :
		      case 'T' :
		      case 'M' :
		      case 'R' :
		      case 'W' :
		      case 'S' :
		      case 'Y' :
		      case 'K' :
		      case 'B' :
		      case 'D' :
		      case 'H' :
		      case 'V' : {b=b; break;}
		      default : return 0;
		      }
		    break;
		  }
		default : 
		  {
                      printf("\n. a = %c, b = %c\n",a,b);
		    Exit("\n. Err. in Are_Compatible\n");
		    return 0;
		  }
		}
	    }
	}
    }
  else
    {
      a = statea[0]; b = stateb[0];
      switch(a)
	{
	case 'A' :
	  {
	    switch(b)
	      {
	      case 'A' :
	      case 'X' : {b=b; break;}
	      default : return 0;
	      }
	    break;
	  }
	case 'R' :
	  {
	    switch(b)
	      {
	      case 'R' :
	      case 'X' : {b=b; break;}
	      default : return 0;
	      }
	    break;
	  }
	case 'N' :
	  {
	    switch(b)
	      {
	      case 'N' :
	      case 'B' :
	      case 'X' : {b=b; break;}
	      default : return 0;
	      }
	    break;
	  }
	case 'B' :
	  {
	    switch(b)
	      {
	      case 'N' :
	      case 'B' :
	      case 'X' : {b=b; break;}
	      default : return 0;
	      }
	    break;
	  }
	case 'D' :
	  {
	    switch(b)
	      {
	      case 'D' :
	      case 'X' : {b=b; break;}
	      default : return 0;
	      }
	    break;
	  }
	case 'C' :
	  {
	    switch(b)
	      {
	      case 'C' :
	      case 'X' : {b=b; break;}
	      default : return 0;
	      }
	    break;
	  }
	case 'Q' :
	  {
	    switch(b)
	      {
	      case 'Q' :
	      case 'Z' :
	      case 'X' : {b=b; break;}
	      default : return 0;
	      }
	    break;
	  }
	case 'Z' :
	  {
	    switch(b)
	      {
	      case 'Q' :
	      case 'Z' :
	      case 'X' : {b=b; break;}
	      default : return 0;
	      }
	    break;
	  }
	case 'E' :
	  {
	    switch(b)
	      {
	      case 'E' :
	      case 'X' : {b=b; break;}
	      default : return 0;
	      }
	    break;
	  }
	case 'G' :
	  {
	    switch(b)
	      {
	      case 'G' :
	      case 'X' : {b=b; break;}
	      default : return 0;
	      }
	    break;
	  }
	case 'H' :
	  {
	    switch(b)
	      {
	      case 'H' :
	      case 'X' : {b=b; break;}
	      default : return 0;
	      }
	    break;
	  }
	case 'I' :
	  {
	    switch(b)
	      {
	      case 'I' :
	      case 'X' : {b=b; break;}
	      default : return 0;
	      }
	    break;
	  }
	case 'L' :
	  {
	    switch(b)
	      {
	      case 'L' :
	      case 'X' : {b=b; break;}
	      default : return 0;
	      }
	    break;
	  }
	case 'K' :
	  {
	    switch(b)
	      {
	      case 'K' :
	      case 'X' : {b=b; break;}
	      default : return 0;
	      }
	    break;
	  }
	case 'M' :
	  {
	    switch(b)
	      {
	      case 'M' :
	      case 'X' : {b=b; break;}
	      default : return 0;
	      }
	    break;
	  }
	case 'F' :
	  {
	    switch(b)
	      {
	      case 'F' :
	      case 'X' : {b=b; break;}
	      default : return 0;
	      }
	    break;
	  }
	case 'P' :
	  {
	    switch(b)
	      {
	      case 'P' :
	      case 'X' : {b=b; break;}
	      default : return 0;
	      }
	    break;
	  }
	case 'S' :
	  {
	    switch(b)
	      {
	      case 'S' :
	      case 'X' : {b=b; break;}
	      default : return 0;
	      }
	    break;
	  }
	case 'T' :
	  {
	    switch(b)
	      {
	      case 'T' :
	      case 'X' : {b=b; break;}
	      default : return 0;
	      }
	    break;
	  }
	case 'W' :
	  {
	    switch(b)
	      {
	      case 'W' :
	      case 'X' : {b=b; break;}
	      default : return 0;
	      }
	    break;
	  }
	case 'Y' :
	  {
	    switch(b)
	      {
	      case 'Y' :
	      case 'X' : {b=b; break;}
	      default : return 0;
	      }
	    break;
	  }
	case 'V' :
	  {
	    switch(b)
	      {
	      case 'V' :
	      case 'X' : {b=b; break;}
	      default : return 0;
	      }
	    break;
	  }
	case 'X' :
	  {
	    switch(b)
	      {
	      case 'A':case 'R':case 'N' :case 'B' :case 'D' :
	      case 'C':case 'Q':case 'Z' :case 'E' :case 'G' :
	      case 'H':case 'I':case 'L' :case 'K' :case 'M' :
	      case 'F':case 'P':case 'S' :case 'T' :case 'W' :
	      case 'Y':case 'V': case 'X' : {b=b; break;}
	      default : return 0;
	      }
	    break;
	  }
	default : 
	  {
	    Exit("\n. Err. in Are_Compatible\n");
	    return 0;
	  }
	}
    }
  return 1;
}

/*********************************************************/

void Hide_Ambiguities(allseq *data)
{
  int i;

  For(i,data->crunch_len)
    {
      if(data->ambigu[i]) 
	{
	  data->wght[i] = 0.0;
	}
    }
}

/*********************************************************/

fit_double Rf(arbre *tree1, arbre *tree2)
{
    int i;
    fit_double rf;


    if(tree1->n_otu != tree2->n_otu) Exit("\n. Both trees should have the same number of tips\n");

    For(i,2*tree1->n_otu-3)
        {
            tree1->t_edges[i]->bip_score = 0;
            tree2->t_edges[i]->bip_score = 0;
        }
  
    Compare_Bip(tree1,tree2);

    rf = 0.0;
    For(i,2*tree1->n_otu-3)
        {
            if((!tree1->t_edges[i]->left->tax) &&
               (!tree1->t_edges[i]->rght->tax))
                {
                    rf += tree1->t_edges[i]->bip_score;
                }
        }
    
    rf = (tree1->n_otu-3) - rf;
    rf /= (fit_double)(tree1->n_otu-3);

    return rf;
}

/*********************************************************/


void Alloc_Bip(arbre *tree)
{
  int i,j,k;

  tree->has_bip = 1;
  
  For(i,2*tree->n_otu-2)
    {
      tree->noeud[i]->bip_size = (int *)mCalloc(3,sizeof(int));
      tree->noeud[i]->bip_node = (node ***)mCalloc(3,sizeof(node **));
      tree->noeud[i]->bip_name = (char ***)mCalloc(3,sizeof(char **));

      For(j,3)
	{
	  tree->noeud[i]->bip_node[j] = 
	    (node **)mCalloc(tree->n_otu,sizeof(node *));

	  tree->noeud[i]->bip_name[j] = 
	    (char **)mCalloc(tree->n_otu,sizeof(char *));
	  
	  For(k,tree->n_otu)
	    tree->noeud[i]->bip_name[j][k] = 
	    (char *)mCalloc(T_MAX_NAME,sizeof(char ));	  
	}
    }
}

/*********************************************************/

void Compare_Bip(arbre *tree1, arbre *tree2)
{
  int i,j,k;
  edge *b1,*b2;
  char **bip1,**bip2;
  int bip_size;
  
  
  For(i,2*tree1->n_otu-3)
    {
      if((!tree1->t_edges[i]->left->tax) &&
	 (!tree1->t_edges[i]->rght->tax))
	{
	  b1 = tree1->t_edges[i];
          
	  For(j,2*tree2->n_otu-3)
	    {
	      if((!tree2->t_edges[j]->left->tax) &&
		 (!tree2->t_edges[j]->rght->tax))
		{
		  b2 = tree2->t_edges[j];
                  
		  if(MMIN(b1->left->bip_size[b1->l_r],b1->rght->bip_size[b1->r_l]) ==
		     MMIN(b2->left->bip_size[b2->l_r],b2->rght->bip_size[b2->r_l]))
		    {
		      bip_size = MMIN(b1->left->bip_size[b1->l_r],b1->rght->bip_size[b1->r_l]);
		      
		      if(b1->left->bip_size[b1->l_r] == b1->rght->bip_size[b1->r_l])
			{
			  if(b1->left->bip_name[b1->l_r][0][0] < b1->rght->bip_name[b1->r_l][0][0])
			    {
			      bip1 = b1->left->bip_name[b1->l_r];
			    }
			  else
			    {
			      bip1 = b1->rght->bip_name[b1->r_l];
			    }
			}
		      else if(b1->left->bip_size[b1->l_r] < b1->rght->bip_size[b1->r_l])
			{
			  bip1 = b1->left->bip_name[b1->l_r];
			}
		      else
			{
			  bip1 = b1->rght->bip_name[b1->r_l];
			}
		      
		      
		      if(b2->left->bip_size[b2->l_r] == b2->rght->bip_size[b2->r_l])
			{
			  if(b2->left->bip_name[b2->l_r][0][0] < b2->rght->bip_name[b2->r_l][0][0])
			    {
			      bip2 = b2->left->bip_name[b2->l_r];
			    }
			  else
			    {
			      bip2 = b2->rght->bip_name[b2->r_l];
			    }
			}
		      else if(b2->left->bip_size[b2->l_r] < b2->rght->bip_size[b2->r_l])
			{
			  bip2 = b2->left->bip_name[b2->l_r];
			}
		      else
			{
			  bip2 = b2->rght->bip_name[b2->r_l];
			}
		      
		      if(bip_size == 1) Exit("\n. Problem in Compare_Bip\n");
                      
                      
		      For(k,bip_size) 
			{
			  if(strcmp(bip1[k],bip2[k])) break;
			}
		      
		      if(k == bip_size)
			{
			  b1->bip_score++;
			  b2->bip_score++;
			  break;
			}
		    }
		}
	    }
	}
    }
}

/*********************************************************/

void Get_Bip(node *a, node *d, arbre *tree)
{
  int i,j;
  
  if(d->tax)
    {      
      d->bip_node[0][0] = d;
      d->bip_size[0]    = 1;
      strcpy(d->bip_name[0][0],d->name);
      
      For(i,3)
	{
	  if(a->v[i] == d)
	    {
	      a->bip_size[i] = 0;
	      For(j,tree->n_otu)
		{
		  if(strcmp(tree->noeud[j]->name,d->name))
		    {
		      a->bip_node[i][a->bip_size[i]] = d;
		      strcpy(a->bip_name[i][a->bip_size[i]],tree->noeud[j]->name);
		      a->bip_size[i]++;
		    }
		}
	      qsort(a->bip_name[i],a->bip_size[i],sizeof(char *),Sort_String);
	      break;
	    }
	}
      return;
    }
  else
    {
      int k;
      int d_a;
      
      d_a = -1;
      
      For(i,3)
	{
	  if(d->v[i] != a) Get_Bip(d,d->v[i],tree);
	  else d_a = i;
	}
      
      d->bip_size[d_a] = 0;
      For(i,3)
	if(d->v[i] != a)
	  {
	    For(j,3)
	      {
		if(d->v[i]->v[j] == d)
		  {
		    For(k,d->v[i]->bip_size[j])
		      {
			d->bip_node[d_a][d->bip_size[d_a]] = d->v[i]->bip_node[j][k];
			strcpy(d->bip_name[d_a][d->bip_size[d_a]],d->v[i]->bip_node[j][k]->name);
			d->bip_size[d_a]++;
		      }
		    break;
		  }
	      }
	  }
      
      qsort(d->bip_name[d_a],d->bip_size[d_a],sizeof(char *),Sort_String);
      
      if(a != tree->n_root)
	{
	  For(i,3)
	    if(a->v[i] == d)
	      {
		a->bip_size[i] = 0;
		For(j,tree->n_otu)
		  {
		    For(k,d->bip_size[d_a])
		      {
			if(d->bip_node[d_a][k] == tree->noeud[j])
			  break;
		      }
		    
		    if(k == d->bip_size[d_a])
 		      {
			a->bip_node[i][a->bip_size[i]] = tree->noeud[j];
			strcpy(a->bip_name[i][a->bip_size[i]],tree->noeud[j]->name);
			a->bip_size[i]++;
		      }
		  }
		
		qsort(a->bip_name[i],a->bip_size[i],sizeof(char *),Sort_String);
		
		break;
	      }
	}
    }
}

/*********************************************************/


fit_double Get_Proportion_Of_False_Positives(fit_double threshold, fit_double *postp2, int *selclass, int size, int *num, int *denom, int print)
{
  int site;
  int numerator, denominator;
    
  numerator = denominator = 0;
  
  For(site, size)
    {
      if((postp2[site] > threshold) && 
	 (selclass[site] != 2)) /* False positive */
	numerator++;
      if(postp2[site] > threshold) 
	{
	  denominator++;
	}
    }
  
  if(num)     *num =   numerator;
  if(denom) *denom = denominator;
  
  if(print) printf("Num=%d Denom=%d\n",numerator,denominator);
  if(!denominator) return -1.0;
  else             return (fit_double)numerator/denominator;
}

/*********************************************************/

fit_double Get_Proportion_Of_False_Negatives(fit_double threshold, fit_double *postp2, int *selclass, int size, int *num, int *denom, int print)
{
    int site;
    int numerator, denominator;
    
    numerator = denominator = 0;
    For(site, size)
        {
            if((postp2[site] < threshold) && 
               (selclass[site] == 2)) /* False negative */
                numerator++;
            if(postp2[site] < threshold) denominator++;
        }

    if(num)     *num =   numerator;
    if(denom) *denom = denominator;

    if(print) printf("Num=%d Denom=%d\n",numerator,denominator);
    if(!denominator) return -1.0;
    else             return (fit_double)numerator/denominator;
}



/*********************************************************/

void Print_Ancestral_Seq(node *a, node *d, arbre *tree)
{
    int i;

/*     printf("d->num=%d\n",d->num); */

    if(!d->tax)
        fprintf(tree->input->fp_output_tree,"%s\n",d->seq);
    else return;
    For(i,3)
        if(d->v[i] != a) Print_Ancestral_Seq(d,d->v[i],tree);
}

/*********************************************************/

void Estimate_Ancestral_States(arbre *tree)
{

#undef  LIM_SCALE
#define LIM_SCALE     300
#undef  LIM_SCALE_VAL
#define LIM_SCALE_VAL 1.E-500

    tree->both_sides = 1;
    Lk(tree);

    Estimate_Ancestral_States_Post_Codon(tree->n_root,
                                         tree->n_root->v[0],
                                         tree->n_root->b[0],
                                         tree);
}

/*********************************************************/

void Estimate_Ancestral_States_Post_Codon(node *a, node *d, edge *b, arbre *tree)
{
    if(d->tax) return;
    else
        {
            int site,catq,catg,k,l;
            fit_double *site_lk;
            fit_double ****p_lk_a,****p_lk_d;
            fit_double max_p;
            
            site_lk = (fit_double *)mCalloc(tree->mod->ns,sizeof(fit_double));

            p_lk_a = (b->left == a)?(b->p_lk_left):(b->p_lk_rght);
            p_lk_d = (b->left == d)?(b->p_lk_left):(b->p_lk_rght);

            For(site,tree->n_pattern)
                {
                    For(catq,tree->mod->n_catq)
                        {
                            For(catg,tree->mod->n_catg)
                                {
                                    For(k,tree->mod->ns) /* d node */
                                        {
                                            For(l,tree->mod->ns) /* a node */
                                                {
                                                    site_lk[k] += 
                                                        b->Pij_rr[catq][catg][k][l] *		    
                                                        p_lk_a[site][catq][catg][l] *  
                                                        p_lk_d[site][catq][catg][k] * 
                                                        (*(b->qmat_struct->t_omega_proba[catq*
                                                                                         b->qmat_struct->n_omega+
                                                                                         (int)k/tree->mod->ns_codon])) *
                                                        tree->mod->r_proba[catg] *
                                                        b->qmat_struct->qmat_proba[catq];
                                                    
                                                }
                                        }
                                }
                        }
                    
                    max_p      = -1.0;
                    For(k,tree->mod->ns)
                        {
                            if(site_lk[k] * tree->mod->pi[k%tree->mod->ns_codon] > max_p)
                                {
                                    max_p = site_lk[k] * tree->mod->pi[k%tree->mod->ns_codon];
                                }
                        }
                }

            Free(site_lk);

            For(k,3)
                if(d->v[k] != a) 
                    Estimate_Ancestral_States_Post_Codon(d,d->v[k],d->b[k],tree);
        }
}

/*********************************************************/

void Estimate_Ancestral_States_Post_NtAA(node *a, node *d, edge *b, arbre *tree)
{
    if(d->tax) return;
    else
        {
            int site,catq,catg,k,l;
            fit_double *site_lk;
            fit_double ****p_lk_a,****p_lk_d;
            fit_double max_p;

            site_lk = (fit_double *)mCalloc(tree->mod->ns,sizeof(fit_double));

            p_lk_a = (b->left == a)?(b->p_lk_left):(b->p_lk_rght);
            p_lk_d = (b->left == d)?(b->p_lk_left):(b->p_lk_rght);

            For(site,tree->n_pattern)
                {
                    For(catq,tree->mod->n_catq)
                        {
                            For(catg,tree->mod->n_catg)
                                {
                                    For(k,tree->mod->ns) /* d node */
                                        {
                                            For(l,tree->mod->ns) /* a node */
                                                {
                                                    site_lk[k] += 
                                                        b->Pij_rr[catq][catg][k][l] *		    
                                                        p_lk_a[site][catq][catg][l] *  
                                                        p_lk_d[site][catq][catg][k] * 
                                                        tree->mod->r_proba[catg] *
                                                        b->qmat_struct->qmat_proba[catq];
                                                    
                                                }
                                        }
                                }
                        }
                    
                    max_p      = -1.0;
                    For(k,tree->mod->ns)
                        {
                            if(site_lk[k] * tree->mod->pi[k] > max_p)
                                {
                                    max_p = site_lk[k] * tree->mod->pi[k];
                                }
                        }
                }

            Free(site_lk);

            For(k,3)
                if(d->v[k] != a) 
                    Estimate_Ancestral_States_Post_NtAA(d,d->v[k],d->b[k],tree);
        }
}

/*********************************************************/

void Print_Tree_Nodes(node *a, node *d, arbre *tree)
{

    if(a->tax) printf("%d..%d\n",d->num+1,a->num+1);
    else printf("%d..%d\n",a->num+1,d->num+1);

    if(d->tax) return;
    else
        {
            int i;
            For(i,3)
                if(d->v[i] != a) Print_Tree_Nodes(d,d->v[i],tree);
        }
}

/*********************************************************/

void Untransform_Probs(fit_double *t_probs, /* transformed probabilities */ 
                       fit_double *probs,   /* 'real' probabilities      */ 
                       int n)   
{
    fit_double sum1,sum2;
    int i;
    
/*     For(i,n) if((fit_double)fabs(t_probs[i] < 1.E-2)) t_probs[i] = (t_probs[i] < 0.)?(-1.E-2):(1.E-2);  */
/*     For(i,n) if((fit_double)fabs(t_probs[i] > 1.E+2)) t_probs[i] = (t_probs[i] < 0.)?(-1.E+2):(1.E+2);  */

    sum1 = sum2 = .0;
    For(i,n) sum1 += (fit_double)fabs(t_probs[i]);

    For(i,n) 
      {
        probs[i] = (fit_double)fabs(t_probs[i])/sum1;
        /* if(probs[i] < 1.E-3) probs[i] = 1.E-3; */
        sum2 += probs[i];
      }
    For(i,n) 
      {
	probs[i] /= sum2;
      } 
}

/*********************************************************/

void Allocate_Num_Parameters(arbre *tree)
{

    tree->numerical_param = (numpar *)mCalloc(1,sizeof(numpar));
    
    tree->numerical_param->replace_num_param = 0;

    tree->numerical_param->param_val = (fit_double *)mCalloc(
                                                         2*tree->n_otu-3+
                                                         tree->mod->ns+
                                                         tree->mod->n_omega+
                                                         7+
                                                         tree->mod->n_omega+1+
                                                         tree->mod->n_catq+1+
							 6,
                                                         sizeof(fit_double ));           

    tree->numerical_param->param_size = (int *)mCalloc(
                                                       2*tree->n_otu-3+
                                                       tree->mod->ns+
                                                       tree->mod->n_omega+
                                                       7+
                                                       tree->mod->n_omega+1+
                                                       tree->mod->n_catq+1+
						       6,
                                                       sizeof(int ));           
    
    tree->numerical_param->param_beg = (int *)mCalloc(
                                                      2*tree->n_otu-3+
                                                      tree->mod->ns+
                                                      tree->mod->n_omega+
                                                      7+
                                                      tree->mod->n_omega+1+
                                                      tree->mod->n_catq+1+
						      6,
                                                      sizeof(int )); 
}

/*********************************************************/

void Replace_Num_Parameters(arbre *tree)
{
    int i;
    int tot_size;
    fit_double tmp;

    tot_size = 0;

    tree->numerical_param->br_len_num  = 0;
    tree->numerical_param->param_size[0] = 2*tree->n_otu-3;
    For(i,tree->numerical_param->param_size[0]) 
        {
            tmp = tree->numerical_param->param_val[i+tot_size];
            tree->numerical_param->param_val[i+tot_size] = tree->t_edges[i]->l;
            tree->t_edges[i]->l = tmp;
/*             printf("l > %f\n",tree->t_edges[i]->l); */
        }
    tree->numerical_param->param_beg[0] = tot_size;
    tot_size += tree->numerical_param->param_size[0];
/*     printf("\n"); */

    tree->numerical_param->pi_num  = 1;
    tree->numerical_param->param_size[1] = tree->mod->ns;
    For(i,tree->numerical_param->param_size[1]) 
        {            
            tmp = tree->numerical_param->param_val[i+tot_size];
            tree->numerical_param->param_val[i+tot_size] = tree->mod->pi[i];
            tree->mod->pi[i] = tmp;
/*             printf("pi > %f\n",tree->mod->pi[i]); */
        }
    tree->numerical_param->param_beg[1] = tot_size;
    tot_size += tree->numerical_param->param_size[1];
/*     printf("\n"); */

    tree->numerical_param->kappa_num  = 2;
    tree->numerical_param->param_size[2] = 1;
    For(i,tree->numerical_param->param_size[2]) 
        {
            tmp = tree->numerical_param->param_val[i+tot_size];
            tree->numerical_param->param_val[i+tot_size] = tree->mod->kappa;
            tree->mod->kappa = tmp;
/*             printf("k > %f\n",tree->mod->kappa); */
        }
    tree->numerical_param->param_beg[2] = tot_size;
    tot_size += tree->numerical_param->param_size[2];
/*     printf("\n"); */

    tree->numerical_param->lambda_num  = 3;
    tree->numerical_param->param_size[3] = 1;
    For(i,tree->numerical_param->param_size[3]) 
        {
            tmp = tree->numerical_param->param_val[i+tot_size];
            tree->numerical_param->param_val[i+tot_size] = tree->mod->lambda;
            tree->mod->lambda = tmp;
/*             printf("lbda > %f\n",tree->mod->lambda); */
        }
    tree->numerical_param->param_beg[3] = tot_size;
    tot_size += tree->numerical_param->param_size[3];
/*     printf("\n"); */

    tree->numerical_param->alpha_num  = 4;
    tree->numerical_param->param_size[4] = 1;
    For(i,tree->numerical_param->param_size[4]) 
        {
            tmp = tree->numerical_param->param_val[i+tot_size];
            tree->numerical_param->param_val[i+tot_size] = tree->mod->alpha;
            tree->mod->alpha = tmp;
/*             printf("alpha > %f\n",tree->mod->alpha); */
        }
    tree->numerical_param->param_beg[4] = tot_size;
    tot_size += tree->numerical_param->param_size[4];
/*     printf("\n"); */

    tree->numerical_param->pinvar_num  = 5;
    tree->numerical_param->param_size[5] = 1;
    For(i,tree->numerical_param->param_size[5]) 
        {
            tmp = tree->numerical_param->param_val[i+tot_size];
            tree->numerical_param->param_val[i+tot_size] = tree->mod->pinvar;
            tree->mod->pinvar = tmp;
/*             printf("pinv > %f\n",tree->mod->pinvar); */
        }
    tree->numerical_param->param_beg[5] = tot_size;
    tot_size += tree->numerical_param->param_size[5];
/*     printf("\n"); */

    tree->numerical_param->omega_num  = 6;
    tree->numerical_param->param_size[6] = tree->mod->n_omega;
    For(i,tree->numerical_param->param_size[6]) 
        {
            tmp = tree->numerical_param->param_val[i+tot_size];
            tree->numerical_param->param_val[i+tot_size] = tree->mod->qmat_struct[0]->omega[i];
            tree->mod->qmat_struct[0]->omega[i] = tmp;
/*             printf("omega > %f\n",tree->mod->qmat_struct[0]->omega[i]); */
        }
    tree->numerical_param->param_beg[6] = tot_size;
    tot_size += tree->numerical_param->param_size[6];
/*     printf("\n"); */

    tree->numerical_param->theta_num  = 7;
    tree->numerical_param->param_size[7] = tree->mod->n_omega;
    For(i,tree->numerical_param->param_size[7]) 
        {
            tmp = tree->numerical_param->param_val[i+tot_size];
            tree->numerical_param->param_val[i+tot_size] = tree->mod->qmat_struct[0]->theta[i];
            tree->mod->qmat_struct[0]->theta[i] = tmp;
/*             printf("theta > %f\n",tree->mod->qmat_struct[0]->theta[i]); */
        }
    tree->numerical_param->param_beg[7] = tot_size;
    tot_size += tree->numerical_param->param_size[7];
/*     printf("\n"); */

    tree->numerical_param->p_omega_num  = 8;
    tree->numerical_param->param_size[8] = tree->mod->n_omega;
    For(i,tree->numerical_param->param_size[8]) 
        {
            tmp = tree->numerical_param->param_val[i+tot_size];
            tree->numerical_param->param_val[i+tot_size] = tree->mod->qmat_struct[0]->trans_omega_proba[i];
            tree->mod->qmat_struct[0]->trans_omega_proba[i] = tmp;
/*             printf("omega_p > %f\n",tree->mod->qmat_struct[0]->omega_proba[i]); */
        }
    tree->numerical_param->param_beg[8] = tot_size;
    tot_size += tree->numerical_param->param_size[8];
/*     printf("\n"); */


    tree->numerical_param->p_qmat_num  = 9;
    tree->numerical_param->param_size[9] = tree->mod->n_catq;
    For(i,tree->numerical_param->param_size[9]) 
        {
            tmp = tree->numerical_param->param_val[i+tot_size];
            tree->numerical_param->param_val[i+tot_size] = tree->mod->qmat_struct[0]->trans_qmat_proba[i];
            tree->mod->qmat_struct[0]->trans_qmat_proba[i] = tmp;
/*             printf("qmat_p > %f\n",tree->mod->qmat_struct[0]->qmat_proba[i]); */
        }
    tree->numerical_param->param_beg[9] = tot_size;
    tot_size += tree->numerical_param->param_size[9];
/*     printf("\n"); */

    tree->numerical_param->gtr_num  = 10;
    tree->numerical_param->param_size[10] = 6;
    For(i,tree->numerical_param->param_size[10]) 
        {
            tmp = tree->numerical_param->param_val[i+tot_size];
            tree->numerical_param->param_val[i+tot_size] = tree->mod->gtr_param[i];
            tree->mod->gtr_param[i] = tmp;
        }
    tree->numerical_param->param_beg[10] = tot_size;
    tot_size += tree->numerical_param->param_size[10];
/*     printf("\n"); */


}

/*********************************************************/

void Init_Num_Parameters(arbre *tree)
{
    int i;
    int tot_size;

    tot_size = 0;

    tree->numerical_param->br_len_num  = 0;
    tree->numerical_param->param_size[0] = 2*tree->n_otu-3;
    For(i,tree->numerical_param->param_size[0]) 
        {
            tree->numerical_param->param_val[i+tot_size] = tree->t_edges[i]->l;
        }
    tree->numerical_param->param_beg[0] = tot_size;
    tot_size += tree->numerical_param->param_size[0];


    tree->numerical_param->pi_num  = 1;
    tree->numerical_param->param_size[1] = tree->mod->ns;
    For(i,tree->numerical_param->param_size[1]) 
        {            
            tree->numerical_param->param_val[i+tot_size] = tree->mod->pi[i];
        }
    tree->numerical_param->param_beg[1] = tot_size;
    tot_size += tree->numerical_param->param_size[1];
    
    
    tree->numerical_param->kappa_num  = 2;
    tree->numerical_param->param_size[2] = 1;
    For(i,tree->numerical_param->param_size[2]) 
        {
            tree->numerical_param->param_val[i+tot_size] = tree->mod->kappa;
        }
    tree->numerical_param->param_beg[2] = tot_size;
    tot_size += tree->numerical_param->param_size[2];
    
    
    tree->numerical_param->lambda_num  = 3;
    tree->numerical_param->param_size[3] = 1;
    For(i,tree->numerical_param->param_size[3]) 
        {
            tree->numerical_param->param_val[i+tot_size] = tree->mod->lambda;
        }
    tree->numerical_param->param_beg[3] = tot_size;
    tot_size += tree->numerical_param->param_size[3];
    
    
    tree->numerical_param->alpha_num  = 4;
    tree->numerical_param->param_size[4] = 1;
    For(i,tree->numerical_param->param_size[4]) 
        {
            tree->numerical_param->param_val[i+tot_size] = tree->mod->alpha;
        }
    tree->numerical_param->param_beg[4] = tot_size;
    tot_size += tree->numerical_param->param_size[4];
    
    
    tree->numerical_param->pinvar_num  = 5;
    tree->numerical_param->param_size[5] = 1;
    For(i,tree->numerical_param->param_size[5]) 
        {
            tree->numerical_param->param_val[i+tot_size] = tree->mod->pinvar;
        }
    tree->numerical_param->param_beg[5] = tot_size;
    tot_size += tree->numerical_param->param_size[5];
    

    tree->numerical_param->omega_num  = 6;
    tree->numerical_param->param_size[6] = tree->mod->n_omega;
    For(i,tree->numerical_param->param_size[6]) 
        {
            tree->numerical_param->param_val[i+tot_size] = tree->mod->qmat_struct[0]->omega[i];
        }
    tree->numerical_param->param_beg[6] = tot_size;
    tot_size += tree->numerical_param->param_size[6];
    
    
    tree->numerical_param->theta_num  = 7;
    tree->numerical_param->param_size[7] = tree->mod->n_omega;
    For(i,tree->numerical_param->param_size[7]) 
        {
            tree->numerical_param->param_val[i+tot_size] = tree->mod->qmat_struct[0]->theta[i];
        }
    tree->numerical_param->param_beg[7] = tot_size;
    tot_size += tree->numerical_param->param_size[7];


    tree->numerical_param->p_omega_num  = 8;
    tree->numerical_param->param_size[8] = tree->mod->n_omega;
    For(i,tree->numerical_param->param_size[8]) 
        {
            tree->numerical_param->param_val[i+tot_size] = tree->mod->qmat_struct[0]->omega_proba[i];
        }
    tree->numerical_param->param_beg[8] = tot_size;
    tot_size += tree->numerical_param->param_size[8];

    tree->numerical_param->p_qmat_num  = 9;
    tree->numerical_param->param_size[9] = tree->mod->n_catq;
    For(i,tree->numerical_param->param_size[9]) 
        {
            tree->numerical_param->param_val[i+tot_size] = tree->mod->qmat_struct[0]->qmat_proba[i];
        }
    tree->numerical_param->param_beg[9] = tot_size;
    tot_size += tree->numerical_param->param_size[9];

    tree->numerical_param->gtr_num  = 10;
    tree->numerical_param->param_size[10] = 6;
    For(i,tree->numerical_param->param_size[10]) 
        {
            tree->numerical_param->param_val[i+tot_size] = tree->mod->gtr_param[i];
        }
    tree->numerical_param->param_beg[10] = tot_size;
    tot_size += tree->numerical_param->param_size[10];
}

/*********************************************************/

int ACGT_to_RY(int base_acgt)
{
    int base_ry;

    base_ry = -1;

    switch(base_acgt)
        {
        case 0 : {base_ry = 0; break;} /* A <-> R */
        case 1 : {base_ry = 1; break;} /* C <-> Y */
        case 2 : {base_ry = 0; break;} /* G <-> R */
        case 3 : {base_ry = 1; break;} /* T <-> Y */
        default : Exit("\n. Unkown base used in ACGT_to_RY\n");
        }

    return base_ry;
}

/*********************************************************/


	/* 	if(a->v[i] == d) */
	/* 	  { */
	/* 	    if(a->b[i]->prob_sel_regime <= 0.1) */
	/* 	      glColor3f(.0,.0,1.); */
	/* 	    else if(a->b[i]->prob_sel_regime > 0.1 && a->b[i]->prob_sel_regime <= 0.2) */
	/* 	      glColor3f(.0,.5,1.); */
	/* 	    else if(a->b[i]->prob_sel_regime > 0.2 && a->b[i]->prob_sel_regime <= 0.3) */
	/* 	      glColor3f(.0,1.,1.); */
	/* 	    else if(a->b[i]->prob_sel_regime > 0.3 && a->b[i]->prob_sel_regime <= 0.4) */
	/* 	      glColor3f(.0,1.,.5); */
	/* 	    else if(a->b[i]->prob_sel_regime > 0.4 && a->b[i]->prob_sel_regime <= 0.5) */
	/* 	      glColor3f(.0,1.,.0); */
	/* 	    else if(a->b[i]->prob_sel_regime > 0.5 && a->b[i]->prob_sel_regime <= 0.6) */
	/* 	      glColor3f(.5,1.,.0); */
	/* 	    else if(a->b[i]->prob_sel_regime > 0.6 && a->b[i]->prob_sel_regime <= 0.7) */
	/* 	      glColor3f(1.,1.,0.); */
	/* 	    else if(a->b[i]->prob_sel_regime > 0.7 && a->b[i]->prob_sel_regime <= 0.8) */
	/* 	      glColor3f(1.,.5,0.); */
	/* 	    else if(a->b[i]->prob_sel_regime > 0.8 && a->b[i]->prob_sel_regime <= 0.9) */
	/* 	      glColor3f(1.,0.,0.); */
	/* 	    else if(a->b[i]->prob_sel_regime > 0.9) */
	/* 	      glColor3f(1.,.0,.0); */
	/* 	    break; */
	/* 	  }  */

/* 	if(a->v[i] == d) */
/* 	  { */
/* 	    if(a->b[i]->l < 1.E-3) render_branch = 0; */
/* 	    glColor3f(0.0,0.0,0.0); */
/* 	    if(a->b[i]->prob_sel_regime <= 0.1) */
/* 	      { */
/* 		glLineWidth(10*0.1); */
/* 		gl2psLineWidth(10*0.1); */
/* 	      } */
/* 	    else if(a->b[i]->prob_sel_regime > 0.1 && a->b[i]->prob_sel_regime <= 0.2) */
/* 	      { */
/* 		glLineWidth(10*0.15); */
/* 		gl2psLineWidth(10*0.1); */
/* 	      } */
/* 	    else if(a->b[i]->prob_sel_regime > 0.2 && a->b[i]->prob_sel_regime <= 0.3) */
/* 	      { */
/* 		glLineWidth(10*0.25); */
/* 		gl2psLineWidth(10*0.1); */
/* 	      } */
/* 	    else if(a->b[i]->prob_sel_regime > 0.3 && a->b[i]->prob_sel_regime <= 0.4) */
/* 	      { */
/* 		glLineWidth(10*0.35); */
/* 		gl2psLineWidth(10*0.1); */
/* 	      } */
/* 	    else if(a->b[i]->prob_sel_regime > 0.4 && a->b[i]->prob_sel_regime <= 0.5) */
/* 	      { */
/* 		glLineWidth(10*0.45); */
/* 		gl2psLineWidth(10*0.4); */
/* 	      } */
/* 	    else if(a->b[i]->prob_sel_regime > 0.5 && a->b[i]->prob_sel_regime <= 0.6) */
/* 	      { */
/* 		glLineWidth(10*0.55); */
/* 		gl2psLineWidth(10*0.4); */
/* 	      } */
/* 	    else if(a->b[i]->prob_sel_regime > 0.6 && a->b[i]->prob_sel_regime <= 0.7) */
/* 	      { */
/* 		glLineWidth(10*0.65); */
/* 		gl2psLineWidth(10*0.4); */
/* 	      } */
/* 	    else if(a->b[i]->prob_sel_regime > 0.7 && a->b[i]->prob_sel_regime <= 0.8) */
/* 	      { */
/* 		glLineWidth(10*0.75); */
/* 		gl2psLineWidth(10*0.7); */
/* 	      } */
/* 	    else if(a->b[i]->prob_sel_regime > 0.8 && a->b[i]->prob_sel_regime <= 0.9) */
/* 	      { */
/* 		glLineWidth(10*0.85); */
/* 		gl2psLineWidth(10*0.7); */
/* 	      } */
/* 	    else if(a->b[i]->prob_sel_regime > 0.9) */
/* 	      { */
/* 		glLineWidth(10*0.95); */
/* 		gl2psLineWidth(10*0.7); */
/* 	      } */
/* 	    break; */
/* 	  } */
/*      } */


void Print_Model_Param(arbre *tree)
{
  model *mod;
  char *s,*line;
  int i;
  mod = tree->mod;
  
  line = (char *)mCalloc(50,sizeof(char));

  
  printf("\n\n ");
  For(i,39) printf("_"); printf("\n");

  For(i,42) line[i] = ' ';
  line[0] = '|';
  sprintf(line+1,"# taxa : %d",tree->n_otu); 
  line[strlen(line)] = ' ';
  line[40] = '|';
  line[41] = '\n';
  printf("%s",line); fflush(NULL);

  For(i,42) line[i] = ' ';
  line[0] = '|';
  sprintf(line+1,"# distinct patterns : %d",tree->n_pattern); 
  line[strlen(line)] = ' ';
  line[40] = '|';
  line[41] = '\n';
  printf("%s",line); fflush(NULL);

  For(i,42) line[i] = ' ';
  line[0] = '|';
  s = Return_Model_Applies_To(mod);
  sprintf(line+1,"Model applies to : %s",s); Free(s);
  line[strlen(line)] = ' ';
  line[40] = '|';
  line[41] = '\n';
  printf("%s",line); fflush(NULL);


  For(i,42) line[i] = ' ';
  line[0] = '|';
  s = Return_Subst_Model_Name(mod);
  sprintf(line+1,"Subst. model name : %s",s); Free(s);
  line[strlen(line)] = ' ';
  line[40] = '|';
  line[41] = '\n';
  printf("%s",line); fflush(NULL);

  For(i,42) line[i] = ' ';
  line[0] = '|';
  s = Return_Switch_Model_Name(mod);
  sprintf(line+1,"Switch model name : %s",s); Free(s);
  line[strlen(line)] = ' ';
  line[40] = '|';
  line[41] = '\n';
  printf("%s",line); fflush(NULL);

  printf("|");
  For(i,39) printf("_");
  printf("|\n");

  Free(line);
}

/*********************************************************/

void Sort_Categories(arbre *tree)
{
  int i,j,ii,jj;
  double tmp;

  if(tree->mod->model_applies_to == CODONS)
    {      
      if(tree->mod->s_opt->sort_omega)
	{
	  if((tree->mod->switch_modelname == SWITCH_S1) || (tree->mod->switch_modelname == SWITCH_S2))
	    {
              Lk(tree);
              printf("\n. Log-likelihood before swap: %f\n",tree->tot_loglk);

	      For(i,tree->mod->n_omega-1)
		{
		  for(j=i+1;j<tree->mod->n_omega;j++)
		    {
		      if(tree->mod->qmat_struct[0]->omega[j] <
			 tree->mod->qmat_struct[0]->omega[i])
			{
			  printf("\n. Swap %f and %f",
				 tree->mod->qmat_struct[0]->omega[j],
				 tree->mod->qmat_struct[0]->omega[i]);

			  tmp = tree->mod->qmat_struct[0]->omega[j];
			  tree->mod->qmat_struct[0]->omega[j] = 
			    tree->mod->qmat_struct[0]->omega[i];
			  tree->mod->qmat_struct[0]->omega[i] = tmp;
			  
			  tmp = tree->mod->qmat_struct[0]->omega_proba[j];
			  tree->mod->qmat_struct[0]->omega_proba[j] = 
			    tree->mod->qmat_struct[0]->omega_proba[i];
			  tree->mod->qmat_struct[0]->omega_proba[i] = tmp;
			  
			  tmp = tree->mod->qmat_struct[0]->trans_omega_proba[j];
			  tree->mod->qmat_struct[0]->trans_omega_proba[j] = 
			    tree->mod->qmat_struct[0]->trans_omega_proba[i];
			  tree->mod->qmat_struct[0]->trans_omega_proba[i] = tmp;

			  tmp = tree->mod->qmat_struct[0]->omega_min[j];
			  tree->mod->qmat_struct[0]->omega_min[j] = 
			    tree->mod->qmat_struct[0]->omega_min[i];
			  tree->mod->qmat_struct[0]->omega_min[i] = tmp;

			  tmp = tree->mod->qmat_struct[0]->omega_max[j];
			  tree->mod->qmat_struct[0]->omega_max[j] = 
			    tree->mod->qmat_struct[0]->omega_max[i];
			  tree->mod->qmat_struct[0]->omega_max[i] = tmp;

                          
                          if(tree->mod->switch_modelname == SWITCH_S2)
                            {
                              double *buff_theta;
                              
                              buff_theta = (double *)mCalloc(tree->mod->n_omega,sizeof(double));
                              
                              For(ii,tree->mod->n_omega)
                                {
                                  For(jj,tree->mod->n_omega)
                                    {
                                      if(ii != jj)
                                        {
                                          /* printf("\n. ii: %d jj: %d theta: %f",ii,jj, */
                                          /*        tree->mod->qmat_struct[0]->theta[MMIN(ii,jj) * tree->mod->n_omega + MMAX(ii,jj) - */
                                          /*                                         (MMIN(ii,jj)+1+(int)pow(MMIN(ii,jj)+1,2))/2]); */
                                          
                                          buff_theta[MMIN(ii,jj) * tree->mod->n_omega + MMAX(ii,jj) -
                                                     (MMIN(ii,jj)+1+(int)pow(MMIN(ii,jj)+1,2))/2] = 
                                            tree->mod->qmat_struct[0]->theta[MMIN(ii,jj) * tree->mod->n_omega + MMAX(ii,jj) -
                                                                             (MMIN(ii,jj)+1+(int)pow(MMIN(ii,jj)+1,2))/2];
                                        }
                                    }
                                }

                              For(ii,tree->mod->n_omega)
                                {
                                  For(jj,tree->mod->n_omega)
                                    {
                                      if(ii != jj)
                                        {
                                          /* printf("\n. ii: %d jj: %d theta: %f",ii,jj, */
                                          /*        tree->mod->qmat_struct[0]->theta[MMIN(ii,jj) * tree->mod->n_omega + MMAX(ii,jj) - */
                                          /*                                         (MMIN(ii,jj)+1+(int)pow(MMIN(ii,jj)+1,2))/2]); */

                                          if(ii == i && jj != j)
                                            {
                                              tree->mod->qmat_struct[0]->theta[MMIN(ii,jj) * tree->mod->n_omega + MMAX(ii,jj) -
                                                                               (MMIN(ii,jj)+1+(int)pow(MMIN(ii,jj)+1,2))/2] = 
                                                buff_theta[MMIN(j,jj) * tree->mod->n_omega + MMAX(j,jj) -
                                                           (MMIN(j,jj)+1+(int)pow(MMIN(j,jj)+1,2))/2] ;
                                            }
                                          
                                          else if(ii != i && jj == j)
                                            {                                      
                                              tree->mod->qmat_struct[0]->theta[MMIN(ii,jj) * tree->mod->n_omega + MMAX(ii,jj) -
                                                                               (MMIN(ii,jj)+1+(int)pow(MMIN(ii,jj)+1,2))/2] =
                                                buff_theta[MMIN(ii,i) * tree->mod->n_omega + MMAX(ii,i) -
                                                           (MMIN(ii,i)+1+(int)pow(MMIN(ii,i)+1,2))/2] ;
                                            }
                                          /* printf("\n. ii: %d jj: %d theta: %f",ii,jj, */
                                          /*        tree->mod->qmat_struct[0]->theta[MMIN(ii,jj) * tree->mod->n_omega + MMAX(ii,jj) - */
                                          /*                                         (MMIN(ii,jj)+1+(int)pow(MMIN(ii,jj)+1,2))/2]); */
                                        }
                                    }
                                }
                              free(buff_theta);
                            }
			}
		    }
		}
              Lk(tree);
              printf("\n. Log-likelihood after swap: %f\n",tree->tot_loglk);
	    }
	  else if(tree->mod->switch_modelname == NO_SWITCH)
	    {
	      For(i,tree->mod->n_catq-1)
		{
		  for(j=i+1;j<tree->mod->n_catq;j++)
		    {
		      if(tree->mod->qmat_struct[0]->omega[j] <
			 tree->mod->qmat_struct[0]->omega[i])
			{
			  printf("\n. Swap %f and %f",
				 tree->mod->qmat_struct[0]->omega[j],
				 tree->mod->qmat_struct[0]->omega[i]);

			  tmp = tree->mod->qmat_struct[0]->omega[j];
			  tree->mod->qmat_struct[0]->omega[j] = 
			    tree->mod->qmat_struct[0]->omega[i];
			  tree->mod->qmat_struct[0]->omega[i] = tmp;
			  
			  tmp = tree->mod->qmat_struct[0]->qmat_proba[j];
			  tree->mod->qmat_struct[0]->qmat_proba[j] = 
			    tree->mod->qmat_struct[0]->qmat_proba[i];
			  tree->mod->qmat_struct[0]->qmat_proba[i] = tmp;
			  
			  tmp = tree->mod->qmat_struct[0]->trans_qmat_proba[j];
			  tree->mod->qmat_struct[0]->trans_qmat_proba[j] = 
			    tree->mod->qmat_struct[0]->trans_qmat_proba[i];
			  tree->mod->qmat_struct[0]->trans_qmat_proba[i] = tmp;

			  tmp = tree->mod->qmat_struct[0]->omega_min[j];
			  tree->mod->qmat_struct[0]->omega_min[j] = 
			    tree->mod->qmat_struct[0]->omega_min[i];
			  tree->mod->qmat_struct[0]->omega_min[i] = tmp;

			  tmp = tree->mod->qmat_struct[0]->omega_max[j];
			  tree->mod->qmat_struct[0]->omega_max[j] = 
			    tree->mod->qmat_struct[0]->omega_max[i];
			  tree->mod->qmat_struct[0]->omega_max[i] = tmp;
			}
		    }
		}
	    }
	  else
	    {
	      printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	      Exit("");		    
	    }
	}     
    }
}

/*********************************************************/

void Get_Score_Mat(code *c_code)
{
  double score[20][20];
  double size[20][20];
  int i,j;
  int codi,codj;
  int aai, aaj;

  For(i,20) For(j,20) score[i][j] = 0.0;
  For(i,20) For(j,20) size[i][j]  = 0.0;

  For(i,61)
    {
      For(j,61)
	{
	  codi = c_code->sense_c[i];
	  codj = c_code->sense_c[j];

	  aai  = Get_AA(c_code,codi);
	  aaj  = Get_AA(c_code,codj);

	  score[aai][aaj] += c_code->n_diff_b_2_codons[i*61+j];
	  size[aai][aaj]  += 1.0;;
	}
    }

  For(i,20)
    {
      For(j,20)
	{
	  score[i][j] /= size[i][j];
 	}      
      score[i][i] = 0.0;
    }

  For(i,20)
    {
      For(j,20)
	{
	  printf("%4d ",(int)rint(score[i][j]));
 	}
      printf("\n");
    }

  For(i,20)
    {
      For(j,20)
	{
	  printf("tree->step_mat[%2d*tree->mod->ns+%2d] = %4d ;",i,j,(int)rint(score[i][j]));
	  printf("\n");
 	}
    }
  fflush(NULL);
}

/*********************************************************/

void KARIN_Analyze_Output(char *annot_tree_file, char *fit_tree_file)
{
  arbre *annot_tree,*fit_tree_w0,*fit_tree_w1,*fit_tree_w2;
  int i,j,nsites;
  FILE *fit_fp_w0,*fit_fp_w1,*fit_fp_w2;
  FILE *fit_fp_stats;
  double sum;
  node *n0,*n1,*n2,*match_node;
  char **node_labels;
  char *line;
  double p0,p1,p2;
  double w0,w1,w2;
  int root_dir,non_root_dir1,non_root_dir2;
  double wb[3];
  double mean,var;

  line = (char *)mCalloc(T_MAX_LINE,sizeof(char));

  annot_tree  = Read_Tree_File(fopen(annot_tree_file,"r"));
  node_labels = KARIN_Branch_Labels_To_Node_Labels(annot_tree);  
  Alloc_Bip(annot_tree);
  Get_Bip(annot_tree->n_root,annot_tree->n_root->v[0],annot_tree);
  Get_Bip(annot_tree->n_root,annot_tree->n_root->v[1],annot_tree);

  strcpy(line,fit_tree_file);
  fit_fp_w0    = fopen(strcat(line,"_trees_w1"),"r");
  strcpy(line,fit_tree_file);
  fit_fp_w1    = fopen(strcat(line,"_trees_w2"),"r");
  strcpy(line,fit_tree_file);
  fit_fp_w2    = fopen(strcat(line,"_trees_w3"),"r");
  strcpy(line,fit_tree_file);
  fit_fp_stats = fopen(strcat(line,"_fitmodel_stats"),"r");


  p0 = p1 = p2 = w0 = w1 = w2 = -1.0;

  For(i,12) 
    if(!fgets(line,T_MAX_LINE,fit_fp_stats))
      {
	printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	Exit("");
      }
  if(!fgets(line,T_MAX_LINE,fit_fp_stats))
    {
      printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Exit("");
    }
  sscanf(line+6,"%lf",&p0);
  if(!fgets(line,T_MAX_LINE,fit_fp_stats))
    {
      printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Exit("");
    }
  sscanf(line+6,"%lf",&p1);
  if(!fgets(line,T_MAX_LINE,fit_fp_stats))
    {
      printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Exit("");
    }
  sscanf(line+6,"%lf",&p2);
  if(!fgets(line,T_MAX_LINE,fit_fp_stats))
    {
      printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Exit("");
    }
  sscanf(line+6,"%lf",&w0);
  if(!fgets(line,T_MAX_LINE,fit_fp_stats))
    {
      printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Exit("");
    }
  sscanf(line+6,"%lf",&w1);
  if(!fgets(line,T_MAX_LINE,fit_fp_stats))
    {
      printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Exit("");
    }
  sscanf(line+6,"%lf",&w2);

  printf("# p0=%f p1=%f p2=%f w0=%f w1=%f w2=%f\n",p0,p1,p2,w0,w1,w2);
  printf("# SITE\tNODE\tDUP\tA\tB\tC\tVAR\n");

  nsites = 0;
  while(1)
    {
      fit_tree_w0 = Read_Tree_File(fit_fp_w0);
      fit_tree_w1 = Read_Tree_File(fit_fp_w1);
      fit_tree_w2 = Read_Tree_File(fit_fp_w2);
      
      if(!fit_tree_w0) break;

      Alloc_Bip(fit_tree_w0);
      Get_Bip(fit_tree_w0->noeud[0],fit_tree_w0->noeud[0]->v[0],fit_tree_w0);
      
      For(j,2*fit_tree_w0->n_otu-3)
	{
	  Make_New_Edge_Label(fit_tree_w0->t_edges[j]);
	  fit_tree_w0->t_edges[j]->n_labels++;
	}

      nsites++;

      For(i,2*fit_tree_w0->n_otu-2)
	{
	  n0 = fit_tree_w0->noeud[i];
	  n1 = fit_tree_w1->noeud[i];
	  n2 = fit_tree_w2->noeud[i];

	  if(!n0->tax)
	    {
	      /* Posterior means of the dn/ds ratios in the three directions */
	      wb[0] = n0->l_prime[0] * w0 + n1->l_prime[0] * w1 + n2->l_prime[0] * w2;
	      wb[1] = n0->l_prime[1] * w0 + n1->l_prime[1] * w1 + n2->l_prime[1] * w2;
	      wb[2] = n0->l_prime[2] * w0 + n1->l_prime[2] * w1 + n2->l_prime[2] * w2;
	      
	      sprintf(n0->b[0]->labels[0],"%f",wb[0]);
	      sprintf(n0->b[1]->labels[0],"%f",wb[1]);
	      sprintf(n0->b[2]->labels[0],"%f",wb[2]);

	      sum = n0->l_prime[0]+n1->l_prime[0]+n2->l_prime[0];
	      if((sum < 0.99) || (sum > 1.01)) { printf("\n(i) sum = %f",sum); Exit("\n"); }

	      sum = n0->l_prime[1]+n1->l_prime[1]+n2->l_prime[1];
	      if((sum < 0.99) || (sum > 1.01)) { printf("\n(ii) sum = %f",sum); Exit("\n"); }

	      sum = n0->l_prime[2]+n1->l_prime[2]+n2->l_prime[2];
	      if((sum < 0.99) || (sum > 1.01)) { printf("\n(iii) sum = %f",sum); Exit("\n"); }

	      match_node = annot_tree->noeud[KARIN_Find_Match_Nodes(fit_tree_w0->noeud[i],fit_tree_w0,annot_tree)];
	      
	      n0 = annot_tree->noeud[KARIN_Find_Match_Nodes(fit_tree_w0->noeud[i]->v[0],fit_tree_w0,annot_tree)];
	      n1 = annot_tree->noeud[KARIN_Find_Match_Nodes(fit_tree_w0->noeud[i]->v[1],fit_tree_w0,annot_tree)];
	      n2 = annot_tree->noeud[KARIN_Find_Match_Nodes(fit_tree_w0->noeud[i]->v[2],fit_tree_w0,annot_tree)];
	      
	      root_dir = non_root_dir1 = non_root_dir2 = -1;
	      For(j,3)
		{
		  if(match_node->v[j] != match_node->anc)
		    {
		      if(non_root_dir1 < 0) non_root_dir1 = j;
		      else                  non_root_dir2 = j;
		    }
		}
	      
	      if(match_node->v[non_root_dir1] == n0) non_root_dir1 = 0;
	      if(match_node->v[non_root_dir1] == n1) non_root_dir1 = 1;
	      if(match_node->v[non_root_dir1] == n2) non_root_dir1 = 2;

	      if(match_node->v[non_root_dir2] == n0) non_root_dir2 = 0;
	      if(match_node->v[non_root_dir2] == n1) non_root_dir2 = 1;
	      if(match_node->v[non_root_dir2] == n2) non_root_dir2 = 2;
	      
	      sum = non_root_dir1 + non_root_dir2;
	      switch((int)sum)
		{
		case 1 : { root_dir = 2; break; }
		case 2 : { root_dir = 1; break; }
		case 3 : { root_dir = 0; break; }
		}
	      
	      fit_tree_w0->noeud[i]->root_dir      = root_dir;
	      fit_tree_w0->noeud[i]->non_root_dir1 = non_root_dir1;
	      fit_tree_w0->noeud[i]->non_root_dir2 = non_root_dir2;

	      fit_tree_w0->noeud[i]->post_w_ave[root_dir]      = wb[root_dir];
	      fit_tree_w0->noeud[i]->post_w_ave[non_root_dir1] = wb[non_root_dir1];
	      fit_tree_w0->noeud[i]->post_w_ave[non_root_dir2] = wb[non_root_dir2];

	      mean = (wb[root_dir] + wb[non_root_dir1] + wb[non_root_dir2])/3.;
	      var  = 
		pow(wb[root_dir]-mean,2) +
		pow(wb[non_root_dir1]-mean,2) +
		pow(wb[non_root_dir2]-mean,2) ;
	      
	      line = node_labels[KARIN_Find_Match_Nodes(fit_tree_w0->noeud[i],fit_tree_w0,annot_tree)];
	      if(!strcmp(line,""))
		{
		  printf("%5d %5d N %10f %10f %10f %10f\n",
			 nsites,
			 i,
			 wb[root_dir],
			 wb[non_root_dir1],
			 wb[non_root_dir2],
			 var);
		}
	      else
		{
/* 		  printf("single %5d %5d %s %10f %10f %10f %10f\n", */
/* 			 nsites, */
/* 			 i, */
/* 			 line, */
/* 			 wb[root_dir], */
/* 			 wb[non_root_dir1], */
/* 			 wb[non_root_dir2], */
/* 			 var); */
		  printf("%5d %5d %s %10f %10f %10f %10f\n",
			 nsites,
			 i,
			 line,
			 wb[root_dir],
			 wb[non_root_dir1],
			 wb[non_root_dir2],
			 var);
		}
	    }
	}

/*       double sum_w0, sum_w1, sum_w2; */
/*       int n_edges0, n_edges1, n_edges2; */
/*       For(i,2*fit_tree_w0->n_otu-2) */
/* 	{ */
/* 	  if(!fit_tree_w0->noeud[i]->tax) */
/* 	    { */
/* 	      line = node_labels[KARIN_Find_Match_Nodes(fit_tree_w0->noeud[i],fit_tree_w0,annot_tree)]; */
	      	      
/* 	      sum_w0   = .0; */
/* 	      n_edges0 =  0; */
/* 	      KARIN_Return_DnDs_Subtree(fit_tree_w0->noeud[i]->v[fit_tree_w0->noeud[i]->root_dir],&sum_w0,&n_edges0,fit_tree_w0); */

/* 	      sum_w1   = .0; */
/* 	      n_edges1 =  0; */
/* 	      KARIN_Return_DnDs_Subtree(fit_tree_w0->noeud[i]->v[fit_tree_w0->noeud[i]->non_root_dir1],&sum_w1,&n_edges1,fit_tree_w0); */
	      
/* 	      sum_w2   = .0; */
/* 	      n_edges2 =  0; */
/* 	      KARIN_Return_DnDs_Subtree(fit_tree_w0->noeud[i]->v[fit_tree_w0->noeud[i]->non_root_dir2],&sum_w2,&n_edges2,fit_tree_w0); */
	      
	      
/* 	      if(!strcmp(line,"")) */
/* 		{ */
/* 		  printf("multip %5d %5d N %10f %10f %10f\n", */
/* 			 nsites, */
/* 			 i, */
/* 			 (sum_w0 + fit_tree_w0->noeud[i]->post_w_ave[fit_tree_w0->noeud[i]->root_dir])/(n_edges0+1.), */
/* 			 (sum_w1 + fit_tree_w0->noeud[i]->post_w_ave[fit_tree_w0->noeud[i]->non_root_dir1])/(n_edges1+1.), */
/* 			 (sum_w2 + fit_tree_w0->noeud[i]->post_w_ave[fit_tree_w0->noeud[i]->non_root_dir2])/(n_edges2+1.));   */
/* 		} */
/* 	      else */
/* 		{ */
/* 		  printf("multip %5d %5d %s %10f %10f %10f\n", */
/* 			 nsites, */
/* 			 i, */
/* 			 line, */
/* 			 (sum_w0 + fit_tree_w0->noeud[i]->post_w_ave[fit_tree_w0->noeud[i]->root_dir])/(n_edges0+1.), */
/* 			 (sum_w1 + fit_tree_w0->noeud[i]->post_w_ave[fit_tree_w0->noeud[i]->non_root_dir1])/(n_edges1+1.), */
/* 			 (sum_w2 + fit_tree_w0->noeud[i]->post_w_ave[fit_tree_w0->noeud[i]->non_root_dir2])/(n_edges2+1.));	   */
/* 		} */
/* 	    } */
/* 	} */
      
/*       line = Write_Tree(fit_tree_w0); */
/*       printf("\ntree %3d %s\n",nsites,line); */

      Free_Tree(fit_tree_w0);
      Free_Tree(fit_tree_w1);
      Free_Tree(fit_tree_w2);
    }



  Free(line);
  fclose(fit_fp_w0);
  fclose(fit_fp_w1);
  fclose(fit_fp_w2);
  fclose(fit_fp_stats);




/*   arbre *annot_tree,*fit_tree; */
/*   int i,j,nsites,n_otu; */
/*   FILE *fit_fp; */
/*   double val; */
/*   int class[3],sum; */
/*   double **mean; */
/*   node *n; */
/*   char **node_labels; */


/*   fit_fp = fopen(fit_tree_file,"r"); */
/*   fit_tree = Read_Tree_File(fit_fp); */

/*   n_otu = fit_tree->n_otu; */

/*   mean = (double **)mCalloc(3,sizeof(double *)); */
/*   For(i,3) mean[i] = (double *)mCalloc(2*fit_tree->n_otu-1,sizeof(double)); */

/*   Free_Tree(fit_tree); */
/*   rewind(fit_fp); */

/*   For(i,3) For(j,2*n_otu-2) mean[i][j] = 0.0; */

/*   nsites = 0; */
/*   while(1) */
/*     { */
/*       fit_tree = Read_Tree_File(fit_fp); */
      
/*       if(!fit_tree) break; */
/*       nsites++; */

/*       For(i,2*fit_tree->n_otu-2) */
/* 	{ */
/* 	  n = fit_tree->noeud[i]; */
/* 	  if(!n->tax) */
/* 	    { */
/* 	      For(j,3) class[j] = 0; */
/* 	      For(j,3) */
/* 		{ */
/* 		  val = n->l_prime[j]; */
/* 		  if(val < 0.01)                    class[0] = 1; */
/* 		  else if(val < 0.51 && val > 0.49) class[1] = 1; */
/* 		  else                              class[2] = 1; */
/* 		} */

/* 	      sum = class[0] + class[1] + class[2]; */
/* 	      switch(sum) */
/* 		{ */
/* 		case 1 : { mean[0][i] += 1.0; break; } */
/* 		case 2 : { mean[1][i] += 1.0; break; } */
/* 		case 3 : { mean[2][i] += 1.0; break; } */
/* 		} */
/* 	    } */
/* 	} */
/*       Free_Tree(fit_tree); */
/*     } */

/*   annot_tree  = Read_Tree_File(fopen(annot_tree_file,"r")); */
/*   node_labels = KARIN_Branch_Labels_To_Node_Labels(annot_tree); */
/*   Alloc_Bip(annot_tree); */
/*   Get_Bip(annot_tree->n_root,annot_tree->n_root->v[0],annot_tree); */
/*   Get_Bip(annot_tree->n_root,annot_tree->n_root->v[1],annot_tree); */
  
/*   rewind(fit_fp); */
/*   fit_tree = Read_Tree_File(fit_fp); */
/*   Alloc_Bip(fit_tree); */
/*   Get_Bip(fit_tree->noeud[0],fit_tree->noeud[0]->v[0],fit_tree); */

  
/*   For(i,2*fit_tree->n_otu-2) */
/*     { */
/*       if(!fit_tree->noeud[i]->tax) */
/* 	{ */
/* 	  mean[0][i] /= (double)nsites; */
/* 	  mean[1][i] /= (double)nsites; */
/* 	  mean[2][i] /= (double)nsites; */
	  	  
/* 	  printf("%s\t[%s]\t%f\t%f\t%f", */
/* 		 fit_tree_file, */
/* 		 node_labels[KARIN_Find_Match_Nodes(fit_tree->noeud[i],fit_tree,annot_tree)], */
/* 		 mean[0][i], */
/* 		 mean[1][i], */
/* 		 mean[2][i]); */
	  
/* 	  printf("\n"); */
/* 	} */
/*     } */

/*   Free(mean); */
/*   fclose(fit_fp); */


}

/*********************************************************/

char **KARIN_Branch_Labels_To_Node_Labels(arbre *tree)
{
  int i,num;
  char **node_labels;

  node_labels = (char **)mCalloc(2*tree->n_otu-2,sizeof(char *));
  For(i,2*tree->n_otu-2) node_labels[i] = (char *)mCalloc(T_MAX_LABEL,sizeof(char ));
 

  num = -1;
  For(i,2*tree->n_otu-2)
    {
      if(tree->t_edges[i]->labels)
	{
	  num = 
	    (tree->t_edges[i]->left == tree->t_edges[i]->rght->anc)?
	    (tree->t_edges[i]->rght->num):
	    (tree->t_edges[i]->left->num);
	  
	  strcpy(node_labels[num],tree->t_edges[i]->labels[0]);
	}
    }

  return node_labels;
}

/*********************************************************/

edge *KARIN_Find_Root_Edge(char **out_tax, int n_out_tax, arbre *tree)
{
  int i,j,k,l;
  node *n;
  int n_matches;
  
  Alloc_Bip(tree);
  Get_Bip(tree->noeud[0],tree->noeud[0]->v[0],tree);
  
  n = NULL;
  for(i=tree->n_otu;i<2*tree->n_otu-2;i++)
    {
      n = tree->noeud[i];

      For(j,3)
	{
	  if(n->bip_size[j] == n_out_tax)
	    {
	      n_matches = 0;
	      For(k,n_out_tax)
		{		  
		  For(l,n->bip_size[j])
		    {
		      if(strstr(n->bip_name[j][l],out_tax[k])) n_matches++;
		    }
		}

	      if(n_matches == n_out_tax)
		{		  
		  return n->b[j];
		}
	    }
	}
    }
  return NULL;
}

/*********************************************************/

void KARIN_Find_Duplication_Node(node *a, node *d, char **tax_set, int n_tax, arbre *tree)
{
  int i;
  int score;

  if(d->tax) return;
  else
    {
      score = Is_Duplication_Node(a,d,tax_set,n_tax,tree);
      if((score > 0) && (a == tree->n_root)) strcpy(d->n_label,"[W]");
      else if(score > 0) strcpy(d->n_label,"[S]");
      
    }
  
  For(i,3)
    if((d->v[i] != a) && (d->b[i] != tree->e_root))
      KARIN_Find_Duplication_Node(d,d->v[i],tax_set,n_tax,tree);  
}

/*********************************************************/
/* Find the node in tree2 that matches node n in tree1 */
int KARIN_Find_Match_Nodes(node *n, arbre *tree1, arbre *tree2)
{
  int i,j,k;
  int match_0, match_1, match_2;

/*   printf("\n. n->num = %d (%d %d %d)",n->num,n->bip_size[0],n->bip_size[1],n->bip_size[2]); */

  match_0 = match_1 = match_2 = -1;
  For(i,2*tree2->n_otu-2)
    {
/*       printf("\nx x->num = %d (%d %d %d)",i,tree2->noeud[i]->bip_size[0],tree2->noeud[i]->bip_size[1],tree2->noeud[i]->bip_size[2]); */

      For(j,3)
	{
	  if(n->bip_size[0] != tree2->noeud[i]->bip_size[j]) match_0 = 0;
	  else
	    {	      
	      For(k,n->bip_size[0])
		{
		  if(strcmp(n->bip_name[0][k],tree2->noeud[i]->bip_name[j][k])) break;
		}
	      if(k == n->bip_size[0])
		{
		  match_0 = 1;
		}
	    }
	  if(match_0 == 1) break;
	}


      For(j,3)
	{
	  if(n->bip_size[1] != tree2->noeud[i]->bip_size[j]) match_1 = 0;
	  else
	    {
	      For(k,n->bip_size[1])
		{
		  if(strcmp(n->bip_name[1][k],tree2->noeud[i]->bip_name[j][k])) break;
		}
	      if(k == n->bip_size[1])
		{
		  match_1 = 1;
		}
	    }
	  if(match_1 == 1) break;
	}


      For(j,3) 
	{
	  if(n->bip_size[2] != tree2->noeud[i]->bip_size[j]) match_2 = 0;
	  else
	    {
	      For(k,n->bip_size[2])
		{
		  if(strcmp(n->bip_name[2][k],tree2->noeud[i]->bip_name[j][k])) break;
		}
	      if(k == n->bip_size[2])
		{
		  match_2 = 1;
		}
	    }
	  if(match_2 == 1) break;
	}


      if((match_0 == 1) && (match_1 == 1) && (match_2 == 1)) break;
    }
      
  if((match_0 + match_1 + match_2) != 3)
    {
/*       printf("\n. No match found"); */
      return -1;
    }
  else
    {
/*       printf("\n. node %d in tree %p matches node %d in tree %p",n->num,tree1,i,tree2); */
      return i;
    }
}

/*********************************************************/

int Is_Duplication_Node(node *a, node *d, char **tax_set, int n_tax, arbre *tree)
{
  int match;
  int i,j,k,l;  
  int **scores;
  int s_l, s_r;


  scores = (int **)mCalloc(3,sizeof(int *));
  For(i,3) scores[i] = (int *)mCalloc(n_tax,sizeof(int));

  For(j,3) For(k,n_tax) scores[j][k] = 0;
  match = 0;

  For(i,3)
    {
      if((d->v[i] != a) && (d->b[i] != tree->e_root))
	{
	  For(j,3)
	    {
	      if(d->v[i]->v[j] == d)
		{
		  For(k,d->v[i]->bip_size[j])
		    {
		      For(l,n_tax) 
			if(strstr(d->v[i]->bip_name[j][k],tax_set[l])) 
			  {
			    scores[i][l]++;
			    break;
			  }
		    }
		  break;
		}
	    }
	}
    }
  
  
  match = 0;
  For(l,n_tax)
    {
      s_l = s_r = -1;
      For(i,3)
	{
	  if((d->v[i] != a) && (d->b[i] != tree->e_root))
	    {
	      if(s_l < 0) s_l = scores[i][l];
	      else        s_r = scores[i][l];
	    }	  
	}
      if((s_l > 0) && (s_r > 0)) match = 1;
      if(((s_l > 1) && (s_r == 0)) || ((s_r > 1) && (s_l == 0))) {match = 0; break;}
    }
  

  For(i,3) Free(scores[i]);
  Free(scores);

  return match;
}

/*********************************************************/

void Get_Rid_Of_Prefix(char delim, arbre *tree)
{
  char *name;
  int i,j,k,start;

  name = (char *)mCalloc(T_MAX_NAME,sizeof(char));
  
  For(i,tree->n_otu)
    {
      k = 0;
      start = NO;
      For(j,strlen(tree->noeud[i]->name)+1)
	{	  
	  if(start) name[k++] = tree->noeud[i]->name[j];
	  if(tree->noeud[i]->name[j] == delim) start = YES;	  
	}
      strcpy(tree->noeud[i]->name,name);
    }

  Free(name);
}

/*********************************************************/

void Unroot_Tree(char **subtrees)
{
  char **tmp_sub;
  int degree,i,j;

  printf("\n. Removing the root...\n");
  
  tmp_sub = Sub_Trees(subtrees[0],&degree);
  if(degree >= 2)
    {
      strcpy(subtrees[2],subtrees[1]);
      Clean_Multifurcation(tmp_sub,degree,2);
      For(j,2) strcpy(subtrees[j],tmp_sub[j]);
    }
  else
    {
      tmp_sub = Sub_Trees(subtrees[1],&degree);
      strcpy(subtrees[2],subtrees[0]);
      Clean_Multifurcation(tmp_sub,degree,2);
      For(j,2) strcpy(subtrees[j],tmp_sub[j]);
    }

  For(i,degree) Free(tmp_sub[i]);
  Free(tmp_sub);
}

/*********************************************************/

void Clean_Multifurcation(char **subtrees, int current_deg, int end_deg)
{

  if(current_deg <= end_deg) return;
  else
    {
      char *s_tmp;
      int i;

      s_tmp = (char *)mCalloc(T_MAX_LINE,sizeof(char));

      strcat(s_tmp,"(\0");
      strcat(s_tmp,subtrees[0]);
      strcat(s_tmp,",\0");
      strcat(s_tmp,subtrees[1]);
      strcat(s_tmp,")\0");
      Free(subtrees[0]);
      subtrees[0] = s_tmp;

      for(i=1;i<current_deg-1;i++) strcpy(subtrees[i],subtrees[i+1]);

      Clean_Multifurcation(subtrees,current_deg-1,end_deg);
    }
}

/*********************************************************/

void Make_All_Tree_Edges(arbre *tree)
{
  int i;

  tree->t_edges = (edge **)mCalloc(2*tree->n_otu-2,sizeof(edge *));

  For(i,2*tree->n_otu-2) tree->t_edges[i] = (edge *)Make_Edge_Light(NULL,NULL,i);
}

/*********************************************************/

void Make_All_Tree_Nodes(arbre *tree)
{
  int i;

  tree->noeud = (node **)mCalloc(2*tree->n_otu-1,sizeof(node *));

  For(i,2*tree->n_otu-1)
    {
      tree->noeud[i] = (node *)Make_Node_Light(i);
      if(i < tree->n_otu) tree->noeud[i]->tax = 1;
      else                tree->noeud[i]->tax = 0;
    }
}

/*********************************************************/

arbre *Make_Tree(int n_otu)
{
  arbre *tree;
  tree = (arbre *)mCalloc(1,sizeof(arbre ));
  Init_Tree(tree,n_otu);
  return tree;
}

/*********************************************************/

void Read_Branch_Label(char *s_d, char *s_a, edge *b)
{
  char *sub_tp;
  char *p;
  int i,pos;

  sub_tp = (char *)mCalloc(T_MAX_LINE,sizeof(char));

  strcpy(sub_tp,s_d);
  strcat(sub_tp,"#");
  p = strstr(s_a,sub_tp);
  i = 0;
  b->n_labels = 0;
  if(p)
    {

      if(!(b->n_labels%BLOCK_LABELS)) Make_New_Edge_Label(b);
      b->n_labels++;
      
      pos = 0;
      do 
	{
	  b->labels[b->n_labels-1][pos] = p[i+strlen(s_d)+1];
	  i++;
	  pos++;
	  if(p[i+strlen(s_d)+1] == '#') 
	    { 
	      b->labels[b->n_labels-1][pos] = '\0';
	      b->n_labels++;
	      if(!(b->n_labels%BLOCK_LABELS)) Make_New_Edge_Label(b);
	      i++;
	      pos=0;
	    }
	}
      while((p[i+strlen(s_d)+1] != ':') && 
	    (p[i+strlen(s_d)+1] != ',') && 
	    (p[i+strlen(s_d)+1] != ')') && 
	    (p[i+strlen(s_d)+1] != '('));

      b->labels[b->n_labels-1][pos] = '\0';
    }

/*   if(p) */
/*     { */
/*       if(b->n_labels == 1) */
/* 	printf("\n. Found label '%s' on edge %3d.\n",b->labels[0],b->num); */
/*       else */
/* 	{ */
/* 	  printf("\n. Found labels "); */
/* 	  For(i,b->n_labels) printf("'%s' ",b->labels[i]); */
/* 	  printf("on edge %3d.",b->num); */
/* 	} */
/*     } */
/*   else */
/*     { */
/*       if(!(b->n_labels%BLOCK_LABELS)) Make_New_Edge_Label(b); */
/*       b->n_labels++; */
/*       strcpy(b->labels[0],"N"); */
/*       printf("\n. Added label '%s' on edge %3d.\n",b->labels[0],b->num); */
/*     } */

  Free(sub_tp);
}

/*********************************************************/

void Read_Branch_Length(char *s_d, char *s_a, arbre *tree)
{
  char *sub_tp;
  char *p;
  edge *b;
  int i;

  b = tree->t_edges[tree->num_curr_branch_available];

  sub_tp = (char *)mCalloc(T_MAX_LINE,sizeof(char));

  For(i,b->n_labels) 
    {
      strcat(s_d,"#");
      strcat(s_d,b->labels[i]);
    }

  strcpy(sub_tp,s_d);
  strcat(sub_tp,":");
  p = strstr(s_a,sub_tp);
  if(p) 
    {
      b->l = atof((char *)p+(int)strlen(sub_tp));
      tree->has_branch_lengths = 1;
    }
  
  sprintf(sub_tp+strlen(sub_tp),"%f",b->l); // Use %f and not %.xf otherwise strstr(s_a,sub_tp) will always return null
  strcat(sub_tp,"::");
  p = strstr(s_a,sub_tp);
  if(p) 
    {
      b->l_prime = atof((char *)p+(int)strlen(sub_tp));
    }

  Free(sub_tp);
}

/*********************************************************/

void Read_Node_Name(node *d, char *s_tree_d, arbre *tree)
{
  int i;

  if(!tree->t_edges[tree->num_curr_branch_available]->n_labels)
    {
      strcpy(d->name,s_tree_d);
    }
  else
    {
      i = 0;
      do
	{
	  d->name[i] = s_tree_d[i];
	  i++;
	}
      while(s_tree_d[i] != '#');
      d->name[i] = '\0';
    }
}

/*********************************************************/

void Make_New_Edge_Label(edge *b)
{
  int i;

  b->labels = (char **)realloc(b->labels,(b->n_labels+BLOCK_LABELS)*sizeof(char *));

  if(!b->labels)
    {
      printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Exit("");
    }
  else
    {
      for(i=b->n_labels;i<b->n_labels+100;i++) b->labels[i] = (char *)mCalloc(T_MAX_LABEL,sizeof(char));
    }
}

/*********************************************************/

void KARIN_Return_DnDs_Subtree(node *n, double *sum_w, int *n_edges, arbre *tree)
{
  if(n->tax) 
    {
      return;
    }
  else
    {
      (*sum_w) += n->post_w_ave[n->non_root_dir1];
      (*sum_w) += n->post_w_ave[n->non_root_dir2];
      (*n_edges) += 2;
      KARIN_Return_DnDs_Subtree(n->v[n->non_root_dir1],sum_w,n_edges,tree);
      KARIN_Return_DnDs_Subtree(n->v[n->non_root_dir2],sum_w,n_edges,tree);
    }
}

/*********************************************************/

void Add_Ambiguities(allseq *ori_seq, allseq *simu_seq, int datatype)
 {
   int i,j;

   For(i,ori_seq->n_otu)
     {
       printf("\n. Adding ambiguities to %s",ori_seq->c_seq[i]->name);
       For(j,ori_seq->crunch_len)
	 {
	   if(Is_Ambigu(ori_seq->c_seq[i]->state+j,datatype,1))
	     {
	       simu_seq->c_seq[i]->state[j] = '?';
	     }
	 }
     }
 }

/*********************************************************/

allseq *Copy_Cseq(allseq *ori, int len, int ns)
{
  allseq *new;
  int i,j,n_otu;
  char **sp_names;

  n_otu = ori->n_otu;

  sp_names = (char **)mCalloc(n_otu,sizeof(char *));
  For(i,n_otu)
    {
      sp_names[i] = (char *)mCalloc(strlen(ori->c_seq[i]->name)+1,sizeof(char));
      strcpy(sp_names[i],ori->c_seq[i]->name);
    }

  new = Make_Cseq(n_otu,len,ori->init_len,sp_names);

  For(j,ori->crunch_len)
    {
      For(i,ori->n_otu) 
	{
	  new->c_seq[i]->state[j]     = ori->c_seq[i]->state[j];
	}

      new->wght[j]   = ori->wght[j];
      new->ambigu[j] = ori->ambigu[j];
      new->invar[j]  = ori->invar[j];
    }

  For(i,ori->n_otu)
    {
      new->c_seq[i]->len = ori->c_seq[i]->len;
      strcpy(new->c_seq[i]->name,ori->c_seq[i]->name);
    }

  new->init_len           = ori->init_len;
  new->clean_len          = ori->clean_len;
  new->crunch_len         = ori->crunch_len;
  For(i,ns) new->b_frq[i] = ori->b_frq[i];
  new->n_otu              = ori->n_otu;

  For(i,n_otu) Free(sp_names[i]);
  Free(sp_names);

  return new;
}

/*********************************************************/

allseq *Make_Cseq(int n_otu, int crunch_len, int init_len, char **sp_names)
{
  allseq *alldata;
  int j;

  alldata                        = (allseq *)mCalloc(1,sizeof(allseq));
  alldata->n_otu                 = n_otu;
  alldata->c_seq                 = (seq **)mCalloc(n_otu,sizeof(seq *));
  alldata->b_frq                 = (fit_double *)mCalloc(T_MAX_ALPHABET,sizeof(fit_double));
  alldata->wght                  = (int *)mCalloc(crunch_len,sizeof(int));
  alldata->ambigu                = (short int *)mCalloc(crunch_len,sizeof(short int));
  alldata->invar                 = (short int *)mCalloc(crunch_len,sizeof(short int));
  alldata->pospatt               = (int **)mCalloc(crunch_len,sizeof(int *));

  alldata->crunch_len = crunch_len;
  alldata->init_len   = init_len;

  For(j,n_otu)
    {
      alldata->c_seq[j]            = (seq *)mCalloc(1,sizeof(seq));
      alldata->c_seq[j]->name      = (char *)mCalloc((int)(strlen(sp_names[j])+1),sizeof(char));
      strcpy(alldata->c_seq[j]->name,sp_names[j]);
      alldata->c_seq[j]->state     = (char *)mCalloc(crunch_len,sizeof(char));
    }

  return alldata;
}

/*********************************************************/

