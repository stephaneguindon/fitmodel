/*

FitModeL :  a program that Fits substitution models to the data
using Maximum Likelihood.  

Copyright (C) Stephane Guindon. jul. 2004 onward

All parts of  the source except where indicated  are distributed under
the GNU public licence.  See http://www.opensource.org for details.

*/


#include "config.h"
#include "utilities.h"
#include "optimiz.h"
#include "lk.h"
#include "free.h"


/*********************************************************/

fit_double Br_Len_Golden(fit_double ax, fit_double bx, fit_double cx, fit_double tol, 
			 fit_double *xmin, edge *b_fcus, arbre *tree)
{
  /* 
     Adapted from 'Numerical Recipes in C, Second Edition, Press et al., 1988'
     
     Optimise a branch length using the Golden section search algorithm. 
  */

   fit_double f1,f2,x0,x1,x2,x3;

   x0=ax;
   x3=cx;
   if (fabs(cx-bx) > fabs(bx-ax)) 
     {
       x1=bx;
       x2=bx+GOLDEN_C*(cx-bx);
     } 
   else 
     {
       x2=bx;
       x1=bx-GOLDEN_C*(bx-ax);
     }
   
   b_fcus->l=fabs(x1);
   f1 = -Lk_At_Given_Edge(tree,b_fcus);
   b_fcus->l=fabs(x2);
   f2 = -Lk_At_Given_Edge(tree,b_fcus);
   while ((fit_double)fabs(x3-x0) > tol*((fit_double)fabs(x1)+(fit_double)fabs(x2))) 
     {
       if (f2 < f1) 
	 {
	   SHFT3(x0,x1,x2,GOLDEN_R*x1+GOLDEN_C*x3)
	   b_fcus->l=(fit_double)fabs(x2);
	   SHFT2(f1,f2,-Lk_At_Given_Edge(tree,b_fcus))
	 } 
       else 
	 {
	   SHFT3(x3,x2,x1,GOLDEN_R*x2+GOLDEN_C*x0)
	   b_fcus->l=(fit_double)fabs(x1);
	   SHFT2(f2,f1,-Lk_At_Given_Edge(tree,b_fcus))
	 }
/*        printf("l=%f lnL=%f\n",b_fcus->l,tree->tot_loglk); */
     }
   if (f1 < f2) 
     {
       *xmin=(fit_double)fabs(x1);
       return -f1;
     } 
   else 
     {
       *xmin=(fit_double)fabs(x2);
       return -f2;
     }
}

/*********************************************************/

int Generic_Brak(fit_double *param,
		 fit_double *ax, fit_double *bx, fit_double *cx, 
		 fit_double *fa, fit_double *fb, fit_double *fc,
		 fit_double lim_inf, fit_double lim_sup,
		 arbre *tree)
{
  /* 
     Adapted from 'Numerical Recipes in C, Second Edition, Press et al., 1988'
     
     Bracket a maximum of the likelihood function with respect to parameter 'param'.
  */

   fit_double ulim,u,r,q,fu,dum;

   u = 0.0;
   *param = *ax;

   if(*param > lim_sup) *param = lim_sup;
   if(*param < lim_inf) *param = lim_inf;
   *fa=-Return_Lk(tree);
   *param = *bx;
   if(*param > lim_sup) *param = lim_sup;
   if(*param < lim_inf) *param = lim_inf;
   *fb=-Return_Lk(tree);
   if (*fb > *fa) {
      SHFT(dum,*ax,*bx,dum)
      SHFT(dum,*fb,*fa,dum)
   }
   *cx=(*bx)+MNBRAK_GOLD*(*bx-*ax);
   *param = (fit_double)fabs(*cx);
   if(*param > lim_sup) *param = lim_sup;
   if(*param < lim_inf) *param = lim_inf;
   *fc=-Return_Lk(tree); 
   while (*fb > *fc) 
     {
        
       if(*ax > lim_sup) *ax = lim_sup;
       if(*ax < lim_inf) *ax = lim_inf;
       if(*bx > lim_sup) *bx = lim_sup;
       if(*bx < lim_inf) *bx = lim_inf;
       if(*cx > lim_sup) *cx = lim_sup;
       if(*cx < lim_inf) *cx = lim_inf;
       if(u   > lim_sup) u   = lim_sup;
       if(u   < lim_inf) u   = lim_inf;

       r=(*bx-*ax)*(*fb-*fc);
       q=(*bx-*cx)*(*fb-*fa);
       u=(*bx)-((*bx-*cx)*q-(*bx-*ax)*r)/
               (2.0*SIGN(MMAX((fit_double)fabs(q-r),MNBRAK_TINY),q-r));
       ulim=(*bx)+MNBRAK_GLIMIT*(*cx-*bx);
       
       if ((*bx-u)*(u-*cx) > lim_inf) 
	 {
	   *param = (fit_double)fabs(u);
	   if(*param > lim_sup) {*param = u = lim_sup;}
	   if(*param < lim_inf) {*param = u = lim_inf;}
	   fu=-Return_Lk(tree);
	   if (fu < *fc) 
	     {
	       *ax=(*bx);
	       *bx=u;
	       *fa=(*fb);
	       *fb=fu;
	       (*ax)=(fit_double)fabs(*ax);
	       (*bx)=(fit_double)fabs(*bx);
	       (*cx)=(fit_double)fabs(*cx);
	       return(0);
	     } 
	   else if (fu > *fb) 
	     {
	       *cx=u;
	       *fc=fu;	
	       (*ax)=(fit_double)fabs(*ax);
	       (*bx)=(fit_double)fabs(*bx);
	       (*cx)=(fit_double)fabs(*cx);
	       return(0);
	     }
	   u=(*cx)+MNBRAK_GOLD*(*cx-*bx);
	   *param = (fit_double)fabs(u);
	   if(*param > lim_sup) {*param = u = lim_sup;}
	   if(*param < lim_inf) {*param = u = lim_inf;}
	   fu=-Return_Lk(tree);
	 } 
       else if ((*cx-u)*(u-ulim) > lim_inf) 
	 {
	   *param = (fit_double)fabs(u);
	   if(*param > lim_sup) {*param = u = lim_sup;}
	   if(*param < lim_inf) {*param = u = lim_inf;}
	   fu=-Return_Lk(tree);
	   if (fu < *fc) 
	     {
	       SHFT(*bx,*cx,u,*cx+MNBRAK_GOLD*(*cx-*bx))
	       *param = (fit_double)fabs(u); 
	       SHFT(*fb,*fc,fu,-Return_Lk(tree))
	     }
	 } 
       else if ((u-ulim)*(ulim-*cx) >= lim_inf) 
	 {
	   u=ulim;
	   *param = (fit_double)fabs(u);
	   if(*param > lim_sup) {*param = u = lim_sup;}
	   if(*param < lim_inf) {*param = u = lim_inf;}
	   fu=-Return_Lk(tree);
	 } 
       else 
	 {
	   u=(*cx)+MNBRAK_GOLD*(*cx-*bx);
	   *param = (fit_double)fabs(u);
	   if(*param > lim_sup) {*param = u = lim_sup;}
	   if(*param < lim_inf) {*param = u = lim_inf;}
	   fu=-Return_Lk(tree);
	 }
       SHFT(*ax,*bx,*cx,u)
       SHFT(*fa,*fb,*fc,fu)


     }
   (*ax)=(fit_double)fabs(*ax);
   (*bx)=(fit_double)fabs(*bx);
   (*cx)=(fit_double)fabs(*cx);
   return(0);
}

/*********************************************************/

fit_double Generic_Golden(fit_double *param, 
			  fit_double ax, fit_double bx, fit_double cx, fit_double tol,
			  fit_double lim_inf, fit_double lim_sup,
			  fit_double *xmin, arbre *tree, int n_iter_max)
{
  /* 
     Adapted from 'Numerical Recipes in C, Second Edition, Press et al., 1988'
     
     Golden section search for finding a maximum of the likelihood function with 
     respect to parameter 'param'.
  */

   fit_double f1,f2,x0,x1,x2,x3;
   int n_iter;
   fit_double init_lk;
   
   x0=ax;
   x3=cx;
   if ((fit_double)fabs(cx-bx) > (fit_double)fabs(bx-ax)) 
     {
       x1=bx;
       x2=bx+GOLDEN_C*(cx-bx);
     } 
   else 
     {
       x2=bx;
       x1=bx-GOLDEN_C*(bx-ax);
     }

   *param=x1;
   if(*param > lim_sup) *param = lim_sup;
   if(*param < lim_inf) *param = lim_inf;

   Lk(tree);

   f1=-tree->tot_loglk;

   *param=x2;

   if(*param > lim_sup) *param = lim_sup;
   if(*param < lim_inf) *param = lim_inf;

   Lk(tree);

   f2=-tree->tot_loglk;

   init_lk = MMAX(-f1,-f2);


   n_iter = 0;
   while ((fit_double)fabs(x3-x0) > tol*((fit_double)fabs(x1)+(fit_double)fabs(x2))) 
     {
       if(x3 > lim_sup) x3 = lim_sup;
       if(x3 < lim_inf) x3 = lim_inf;
       if(x2 > lim_sup) x2 = lim_sup;
       if(x2 < lim_inf) x2 = lim_inf;
       if(x1 > lim_sup) x1 = lim_sup;
       if(x1 < lim_inf) x1 = lim_inf;
       if(x0 > lim_sup) x0 = lim_sup;
       if(x0 < lim_inf) x0 = lim_inf;

       if (f2 < f1) 
	 {
	   SHFT3(x0,x1,x2,GOLDEN_R*x1+GOLDEN_C*x3)
	   *param=x2;
	   if(*param > lim_sup) *param = lim_sup;
	   if(*param < lim_inf) *param = lim_inf;

	   Lk(tree);

	   SHFT2(f1,f2,-tree->tot_loglk)
	 } 
       else 
	 {
	   SHFT3(x3,x2,x1,GOLDEN_R*x2+GOLDEN_C*x0)
	   *param=x1;
	   if(*param > lim_sup) *param = lim_sup;
	   if(*param < lim_inf) *param = lim_inf;

	   Lk(tree);

	   SHFT2(f2,f1,-tree->tot_loglk)
	 }
       
       if((n_iter++ > n_iter_max) && (tree->tot_loglk > init_lk)) break;
       
/*        printf("param=%E %f %d\n",*param,tree->tot_loglk,n_iter); */
     }

   if(tree->tot_loglk < init_lk-MIN_DIFF_LK) 
     {
       printf("\n. %f %f\n",
	      tree->tot_loglk,
	      init_lk);

       printf("\n. WARNING : optimisation failed !\n");
     }

   if (f1 < f2) 
    {
       *xmin=x1;
       return f1;
     } 
   else 
     {
       *xmin=x2;
       return f2;
     }
}

/*********************************************************/

fit_double Generic_Brent(fit_double *param, 
			 fit_double ax, fit_double bx, fit_double cx, fit_double tol, 
			 fit_double lim_inf,
			 fit_double lim_sup,
			 arbre *tree, int n_iter_max)
{
  /* 
     Adapted from 'Numerical Recipes in C, Second Edition, Press et al., 1988'
     
     Brent's method for finding a maximum of the likelihood function with 
     respect to parameter 'param'.
  */

  int iter;
  fit_double a,b,d,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
  fit_double e=0.0;
  fit_double init_loglk,max_loglk;
  fit_double bestx;

  d=0.0;
  a=((ax < cx) ? ax : cx);
  b=((ax > cx) ? ax : cx);
  x=w=v=bx;
  bestx = *param;
  Lk(tree);
  fw=fv=fx=-tree->tot_loglk;
  init_loglk = tree->tot_loglk;
  max_loglk = tree->tot_loglk;

  /* printf("init param=%f loglk=%f [%f %f]\n",*param,tree->tot_loglk,lim_inf,lim_sup); */


  for(iter=1;iter<=BRENT_ITMAX;iter++) 
    {
      xm=0.5*(a+b);
      tol2=2.0*(tol1=tol*(fit_double)fabs(x)+BRENT_ZEPS);
      if((fit_double)fabs(x-xm) <= (tol2-0.5*(b-a))) 
	{
	  if(tree->tot_loglk < init_loglk - MIN_DIFF_LK)
              {
                  printf("\n. WARNING : Brent failed\n");
/* 		  printf("\n\n= Debug info below."); */
/* 		  printf("\n= Go back to normal ?\n"); */
/* 		  *param = initx; */
/* 		  Lk(tree); */
/* 		  printf("param=%f loglk=%f\n",*param,tree->tot_loglk); */
/* 		  Exit("\n"); */
              }
          /* printf("MIN_DIFF param=%f loglk=%f [%f %f]\n",*param,tree->tot_loglk,lim_inf,lim_sup); */

          *param = bestx;
          Lk(tree);
          return -fx;
	}
      
      if((fit_double)fabs(e) > tol1) 
	{
	  r=(x-w)*(fx-fv);
	  q=(x-v)*(fx-fw);
	  p=(x-v)*q-(x-w)*r;
	  q=2.0*(q-r);
	  if(q > 0.0) p = -p;
	  q=(fit_double)fabs(q);
	  etemp=e;
	  e=d;
	  if((fit_double)fabs(p) >= (fit_double)fabs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
	    d=BRENT_CGOLD*(e=(x >= xm ? a-x : b-x));
	  else{
	    d=p/q;
	    u=x+d;
	    if (u-a < tol2 || b-u < tol2)
	      d=SIGN(tol1,xm-x);
	  }
	}
      else
	{
	  d=BRENT_CGOLD*(e=(x >= xm ? a-x : b-x));
	}
      u=((fit_double)fabs(d) >= tol1 ? x+d : x+SIGN(tol1,d));
      *param=u;
      /* if(u < lim_inf) u = lim_inf; */
      /* if(u > lim_sup) u = lim_sup; */
      Lk(tree);
      fu=-tree->tot_loglk;

      if(tree->tot_loglk > max_loglk)
	{
	  max_loglk = tree->tot_loglk;
	  bestx = *param;
	}

      /* printf("param=%f loglk=%f [%f %f]\n",*param,tree->tot_loglk,lim_inf,lim_sup); */
      /*       printf("param=%f loglk=%f\n",*param,tree->tot_loglk); */

      if(fu <= fx)
	{
	  if(iter > n_iter_max) 
	    {
	      if(tree->tot_loglk < init_loglk - MIN_DIFF_LK)
		{
                  printf("\n. WARNING : Brent failed\n");
                }

              /* printf("ITER MAX param=%f loglk=%f initlk=%f [%f %f]\n",*param,tree->tot_loglk,init_loglk,lim_inf,lim_sup); */

	      *param = bestx;
	      Lk(tree);
	      return tree->tot_loglk;
	    }
	  if(u >= x) a=x; else b=x;
	  SHFT(v,w,x,u)
	    SHFT(fv,fw,fx,fu)
	    } 
      else
	{
	  if (u < x) a=u; else b=u;
	  if (fu <= fw || w == x) 
	    {
	      v=w;
	      w=u;
	      fv=fw;
	      fw=fu;
	    } 
	  else if (fu <= fv || v == x || v == w) {
            v=u;
            fv=fu;
	  }
	}
    }
  printf("\n Brent ended \n");
  printf("\n. Too many iterations in BRENT !");
  *param = bestx;
  Lk(tree);
  return(-1);
}

/*********************************************************/

fit_double RRparam_GTR_Golden(fit_double ax, fit_double bx, fit_double cx, fit_double tol, 
			  fit_double *xmin, arbre *tree, allseq *alldata, fit_double *param, int n_iter_max)
{
  /* 
     Adapted from 'Numerical Recipes in C, Second Edition, Press et al., 1988'
     
     Golden section search for finding a maximum of the likelihood function with 
     respect to one of the relative rates of substitution of the GTR model.
  */

   fit_double f1,f2,x0,x1,x2,x3;
   int n_iter;


   x0=ax;
   x3=cx;
   if ((fit_double)fabs(cx-bx) > (fit_double)fabs(bx-ax)) 
     {
       x1=bx;
       x2=bx+GOLDEN_C*(cx-bx);
     } 
   else 
     {
       x2=bx;
       x1=bx-GOLDEN_C*(bx-ax);
     }
   (*param)=x1;

   Lk(tree);
   f1=-tree->tot_loglk;
   (*param)=x2;

   Lk(tree);
   f2=-tree->tot_loglk;

   n_iter = 0;
   while ((fit_double)fabs(x3-x0) > tol*((fit_double)fabs(x1)+(fit_double)fabs(x2))) 
     {

       if (f2 < f1) 
	 {
	   SHFT3(x0,x1,x2,GOLDEN_R*x1+GOLDEN_C*x3)
	   (*param)=x2;
	   Lk(tree);
	   SHFT2(f1,f2,-tree->tot_loglk)
	 } 
       else 
	 {
	   SHFT3(x3,x2,x1,GOLDEN_R*x2+GOLDEN_C*x0)
	   (*param)=x1;
	   Lk(tree);
	   SHFT2(f2,f1,-tree->tot_loglk)
	 }
       
       if(n_iter++ > n_iter_max) break;
       
/*        printf("p=%E %f\n",(*param),tree->tot_loglk); */
     }
   if (f1 < f2) 
    {
       *xmin=x1;
       return f1;
     } 
   else 
     {
       *xmin=x2;
       return f2;
     }
}

/*********************************************************/

int Br_Len_Brak(fit_double *ax, fit_double *bx, fit_double *cx,
		fit_double *fa, fit_double *fb, fit_double *fc,
		edge *b_fcus, arbre *tree)
{
  /* 
     Adapted from 'Numerical Recipes in C, Second Edition, Press et al., 1988'
     
     Bracket a maximum of the likelihood function with respect to the length
     of branch 'b_fcus'.
  */

   fit_double ulim,u,r,q,fu,dum;

   b_fcus->l = (fit_double)fabs(*ax);
   *fa=-Lk_At_Given_Edge(tree,b_fcus);
   b_fcus->l = (fit_double)fabs(*bx);
   *fb=-Lk_At_Given_Edge(tree,b_fcus);
   if (*fb > *fa) {
      SHFT(dum,*ax,*bx,dum)
      SHFT(dum,*fb,*fa,dum)
   }
   *cx=(*bx)+MNBRAK_GOLD*(*bx-*ax);
   b_fcus->l = (fit_double)fabs(*cx);
   *fc=-Lk_At_Given_Edge(tree,b_fcus);
   while (*fb > *fc)
     {

       r=(*bx-*ax)*(*fb-*fc);
       q=(*bx-*cx)*(*fb-*fa);
       u=(*bx)-((*bx-*cx)*q-(*bx-*ax)*r)/
               (2.0*SIGN(MMAX((fit_double)fabs(q-r),MNBRAK_TINY),q-r));
       ulim=(*bx)+MNBRAK_GLIMIT*(*cx-*bx);
       
       if ((*bx-u)*(u-*cx) > 0.0)
	 {
	   b_fcus->l = (fit_double)fabs(u);
	   fu=-Lk_At_Given_Edge(tree,b_fcus);
	   if (fu < *fc)
	     {
	       *ax=(*bx);
	       *bx=u;
	       *fa=(*fb);
	       *fb=fu;
	       (*ax)=(fit_double)fabs(*ax);
	       (*bx)=(fit_double)fabs(*bx);
	       (*cx)=(fit_double)fabs(*cx);
	       return(0);
	     }
	   else if (fu > *fb)
	     {
	       *cx=u;
	       *fc=fu;
	       (*ax)=(fit_double)fabs(*ax);
	       (*bx)=(fit_double)fabs(*bx);
	       (*cx)=(fit_double)fabs(*cx);
	       return(0);
	     }
	   u=(*cx)+MNBRAK_GOLD*(*cx-*bx);
	   b_fcus->l = (fit_double)fabs(u);
	   fu=-Lk_At_Given_Edge(tree,b_fcus);
	 }
       else if ((*cx-u)*(u-ulim) > 0.0)
	 {
	   b_fcus->l = (fit_double)fabs(u);
	   fu=-Lk_At_Given_Edge(tree,b_fcus);
	   if (fu < *fc)
	     {
	       SHFT(*bx,*cx,u,*cx+MNBRAK_GOLD*(*cx-*bx))
	       b_fcus->l = (fit_double)fabs(u);
	       SHFT(*fb,*fc,fu,-Lk_At_Given_Edge(tree,b_fcus))
	     }
	 }
       else if ((u-ulim)*(ulim-*cx) >= 0.0)
	 {
	   u=ulim;
	   b_fcus->l = (fit_double)fabs(u);
	   fu=-Lk_At_Given_Edge(tree,b_fcus);
	 }
       else
	 {
	   u=(*cx)+MNBRAK_GOLD*(*cx-*bx);
	   b_fcus->l = (fit_double)fabs(u);
	   fu=-Lk_At_Given_Edge(tree,b_fcus);
	 }
       SHFT(*ax,*bx,*cx,u)
       SHFT(*fa,*fb,*fc,fu)
      }
   (*ax)=(fit_double)fabs(*ax);
   (*bx)=(fit_double)fabs(*bx);
   (*cx)=(fit_double)fabs(*cx);
   return(0);
}

/*********************************************************/

fit_double Br_Len_Brent(fit_double ax, fit_double bx, fit_double cx, fit_double tol,
		    fit_double *xmin, edge *b_fcus, arbre *tree, int n_iter_max)
{

  /* 
     Adapted from 'Numerical Recipes in C, Second Edition, Press et al., 1988'
     
     Brent's method for finding a maximum of the likelihood function with 
     respect to the length of branch 'b_fcus'.
  */


  int iter;
  fit_double a,b,d,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
  fit_double e=0.0;
  fit_double init_loglk, max_loglk;
  fit_double bestx;

  d=0.0;
  a=((ax < cx) ? ax : cx);
  b=((ax > cx) ? ax : cx);
  x=w=v=bx;
  b_fcus->l = (fit_double)fabs(bx);
  fw=fv=fx=-Lk_At_Given_Edge(tree,b_fcus);
  init_loglk = -fw;
  max_loglk = UNLIKELY;
  bestx = bx;

/*   printf("BE %E %f\n",b_fcus->l,Lk_At_Given_Edge(tree,b_fcus)); */

  for(iter=1;iter<=BRENT_ITMAX;iter++)
    {
      xm=0.5*(a+b);
      tol2=2.0*(tol1=tol*(fit_double)fabs(x)+BRENT_ZEPS);
      if((fit_double)fabs(x-xm) <= (tol2-0.5*(b-a)))
	{
	  *xmin=x;
	  Lk_At_Given_Edge(tree,b_fcus);

          if(tree->tot_loglk < init_loglk - MIN_DIFF_LK)
	    {
	      printf("\n. WARNING : Brent (branch lengths) failed (diff=%E,%E)\n",init_loglk,tree->tot_loglk);
	      b_fcus->l = bestx;
	      Lk_At_Given_Edge(tree,b_fcus);
	    }
	  return -fx;
	}
      
      if((fit_double)fabs(e) > tol1)
	{
	  r=(x-w)*(fx-fv);
	  q=(x-v)*(fx-fw);
	  p=(x-v)*q-(x-w)*r;
	  q=2.0*(q-r);
	  if(q > 0.0) p = -p;
	  q=(fit_double)fabs(q);
	  etemp=e;
	  e=d;
	  if((fit_double)fabs(p) >= (fit_double)fabs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
	    d=BRENT_CGOLD*(e=(x >= xm ? a-x : b-x));
	  else{
	    d=p/q;
	    u=x+d;
	    if (u-a < tol2 || b-u < tol2)
	      d=SIGN(tol1,xm-x);
	  }
	}
      else
	{
	  d=BRENT_CGOLD*(e=(x >= xm ? a-x : b-x));
	}
      u=((fit_double)fabs(d) >= tol1 ? x+d : x+SIGN(tol1,d));
      if(u<BL_MIN) u = BL_MIN;
      b_fcus->l=(fit_double)fabs(u);
      fu=-Lk_At_Given_Edge(tree,b_fcus);

      if(tree->tot_loglk > max_loglk)
	{
	  max_loglk = tree->tot_loglk;
	  bestx = (fit_double)fabs(u);
	}

/*       printf("edge %d l=%f lnL=%f\n",b_fcus->num,b_fcus->l,fu); */
      if(fu <= fx)
	{
	  if(iter > n_iter_max) 
	    {
	      printf("\n. WARNING : too many iterations in Brent\n");
	      b_fcus->l = (fit_double)fabs(bx);
	      Lk_At_Given_Edge(tree,b_fcus);
	      return tree->tot_loglk;
	    }
	  if(u >= x) a=x; else b=x;
	  SHFT(v,w,x,u)
	    SHFT(fv,fw,fx,fu)
	    }
      else
	{
	  if (u < x) a=u; else b=u;
	  if (fu <= fw || w == x)
	    {
	      v=w;
	      w=u;
	      fv=fw;
            fw=fu;
	    }
	  else if (fu <= fv || v == x || v == w) {
            v=u;
            fv=fu;
	  }
	}
    }
  printf("Too many iterations in BRENT");
  return(-1);
}

/*********************************************************/

void Optimize_Single_Param_Generic(arbre *tree, fit_double *param, 
				   fit_double start, 
				   fit_double lim_inf, fit_double lim_sup,
				   int n_max_iter)
{
  /* Call the optimisation routine */
  
  fit_double ax,bx,cx;
  fit_double lk_init;

  lk_init = tree->tot_loglk;

  tree->mod->s_opt->opt_bl = 
  tree->both_sides         = 0;

  ax =  lim_inf;
  if((*param < lim_inf) ||
     (*param > lim_sup)) bx = lim_inf + (lim_sup-lim_inf)/2.;
  else if(((*param < lim_inf + 1.E-10) && (*param > lim_inf - 1.E-10)) ||
	  ((*param < lim_sup + 1.E-10) && (*param > lim_sup - 1.E-10)))
    {
      printf("\n. WARNING : one parameter on the boundary ! (value=%f)\n",*param);
      bx = lim_inf+(lim_sup-lim_inf)/2.;
    }
  else bx = start;
  cx = lim_sup;
  
/*   Generic_Brak(param, */
/* 	       &ax,&bx,&cx, */
/* 	       &fa,&fb,&fc, */
/* 	       lim_inf, lim_sup, */
/* 	       tree); */

/*   Generic_Golden(param, */
/* 		 ax,bx,cx,1.e-5, */
/* 		 lim_inf,lim_sup, */
/* 		 param,tree,n_max_iter); */
  
  Generic_Brent(param,
		ax,bx,cx,1.E-5,
		lim_inf,lim_sup,
		tree,n_max_iter);

  if(tree->tot_loglk < lk_init-MIN_DIFF_LK) 
    {
      printf("\n. curr = %f, init = %f\n",tree->tot_loglk,lk_init);
      Exit("\n. Optimisation failed !\n");
    }
}

/*********************************************************/

void Round_Optimize(arbre *tree, allseq *data)
{

  /* Core of the optimisation procedure.
     Branch lengths are first optimised one by one. 
     Parameters of the substitution model are adjusted
     thereafter. These two steps are repeated until 
     convergence to a maximum of the likelihood function.
  */
     

  int n_round,each,n_passes;
  fit_double lk_old, lk_new, lk_after_br_lengths;
  int end_of_br_lengths, end_of_params;
  node *root;
  fit_double diff_lk;
  char *s;


  root = tree->n_root;
  lk_new = tree->tot_loglk;
  lk_old = lk_after_br_lengths = UNLIKELY;
  n_round = 0;
  each = 0;
  diff_lk = 1.e-5;
  end_of_br_lengths = end_of_params = 0;

  tree->mod->tpos_ols             = 0;
  tree->mod->update_bl_using_tpos = tree->mod->s_opt->opt_tpos;
  tree->mod->s_opt->opt_bl        = 0;
  tree->both_sides                = 1;
  tree->mod->s_opt->print         = 0;
  Lk(tree);
  tree->mod->s_opt->print         = 1;


  while(n_round < ROUND_MAX)
    {
      For(n_passes,1)
	{
	  if(!tree->mod->s_opt->opt_tpos)
	    {
	      (!((n_round+2)%2))?(root=tree->noeud[0]):(root=tree->noeud[tree->n_otu-1]);
	      printf("\n. [Optimisation step] => branch lengths...\n");
	      
	      Optimize_Br_Len_Serie(root,
				    root->v[0],
				    root->b[0],
				    tree,
				    data);
	    }
	  else
	    {
	      Optimize_Tpos_Serie(root,
				  root->v[0],
				  root->b[0],
				  tree,
				  data,
				  1);
	    }

	  tree->mod->s_opt->opt_bl = 0;
	  tree->both_sides         = 1;
	  Lk(tree);
	}
      
      lk_after_br_lengths = tree->tot_loglk;
      
      end_of_params     = 0;
      end_of_br_lengths = 0;
      if((fit_double)fabs(lk_after_br_lengths-lk_old) < diff_lk) end_of_br_lengths = 1;


      if(!each)
	{
	  each = 1;
	  Optimiz_All_Free_Param(tree,1);
	  tree->mod->s_opt->opt_bl        = 0;
	  tree->both_sides                = 1;
	  Lk(tree);
	  if((fit_double)fabs(tree->tot_loglk - lk_after_br_lengths) < diff_lk) end_of_params = 1;
	}
      
      if(tree->input->n_data_sets == 1)
	{
	  s = Write_Tree(tree); /* Output tree */	  
	  rewind(tree->input->fp_output_tree);
	  fprintf(tree->input->fp_output_tree,"WARNING : this tree is not the maximum likelihood one\n");
	  fprintf(tree->input->fp_output_tree,"%s\n",s);
	  Free(s);
	}

      lk_new = tree->tot_loglk;
      if(end_of_br_lengths && end_of_params) break;
      lk_old  = lk_new;
      n_round++;
      each--;
    }
}

/*********************************************************/

void Optimize_Br_Len_Serie(node *a, node *d, edge *b_fcus,
			  arbre *tree,allseq *alldata)
{
  /* Pre-order traversal for the optimisation of branch lengths
     one by one.
  */

  int i;
  fit_double l_infa,l_max,l_infb;
  fit_double lk_init;
  
  lk_init = tree->tot_loglk;

  l_infa = 2.*b_fcus->l;
  l_max  = b_fcus->l;
  l_infb = 0.5*b_fcus->l;
 
/*   printf("AV %E %f %f\n",b_fcus->l,Lk_At_Given_Edge(tree,b_fcus),lk_init); */
  Br_Len_Brent(l_infa,l_max,l_infb,
	       1.e-5,
	       &(b_fcus->l),
	       b_fcus,tree,1000);

  /* Golden method is generally slower than Brent */
  /*   Br_Len_Golden(l_infa,l_max,l_infb, */
  /* 		1.e-5, */
  /* 		&(b_fcus->l), */
  /* 		b_fcus,tree); */
  /*   printf("Edge %d -> %20f\n",b_fcus->num,tree->tot_loglk); */

  if(tree->tot_loglk < lk_init - MIN_DIFF_LK)
    {
      printf("\n. %f %f %f %f",l_infa,l_max,l_infb,b_fcus->l);
      printf("\n. %f -- %f \n",lk_init,tree->tot_loglk);
      printf("\n. Optimisation failed. File %s, line %d\n",__FILE__,__LINE__);
      Exit("\n");
    }

/*   printf("Edge %3d -> %f %f\n", */
/* 	 b_fcus->num, */
/* 	 tree->tot_loglk, */
/* 	 b_fcus->l); fflush(NULL); */
    
  
  if(d->tax) return;
  else For(i,3) if(d->v[i] != a)
    {
      Update_P_Lk(tree,d->b[i],d);
      Optimize_Br_Len_Serie(d,d->v[i],d->b[i],tree,alldata);
    }
  For(i,3) if((d->v[i] == a) && !(d->v[i]->tax)) Update_P_Lk(tree,d->b[i],d);
}

/*********************************************************/

void Optimize_Tpos_Serie(node *a, node *d, edge *b_fcus, 
			 arbre *tree,allseq *alldata, int n_passes)
{    
  /* Only used with serial sample data */

  if(d->tax) return;
  else 
    {
      int i;
      fit_double infa,max,infb;      

      infa = d->tpos_lim_above;
      max  = d->tpos;
      infb = d->tpos_lim_below;
      
      Generic_Brent(&(d->tpos),
		    infa,max,infb,1.e-10,
		    1.e-10,1.e+10,
		    tree,n_passes);

      
      printf("Node %4d -> %f ; %5.2f < %5.2f < %5.2f \n",
	     d->num,tree->tot_loglk,
	     d->tpos_lim_above,
	     d->tpos,
	     d->tpos_lim_below);
      
      For(i,3) if(d->v[i] != a)
	Optimize_Tpos_Serie(d,d->v[i],d->b[i],tree,alldata,n_passes);
    }
}

/*********************************************************/

void Print_Lk_Progress(arbre *tree, fit_double lk_new, fit_double lk_old, int n_iter)
{
  if(!n_iter)
    printf("\n. Log(lk) :               * -> %15.6f ",lk_new);
  else
    printf("\n. Log(lk) : %15.6f -> %15.6f ",lk_old,lk_new);
  fflush(stdout);
}

/*********************************************************/

void Optimiz_All_Free_Param(arbre *tree, int verbose)
{
  /* Optimisation of the parameters of the substitution model */

    int  init_both_sides, init_derivatives;
    int failed;
    fit_double *init_values;
    int i;
    

    init_both_sides     = tree->both_sides;
    init_derivatives    = tree->mod->s_opt->opt_bl;
    tree->both_sides    = 0;
    tree->mod->s_opt->opt_bl = 0;
    
    init_values = (fit_double *)mCalloc(100,sizeof(fit_double));

    
    if(tree->mod->model_number == 7)
        {
            failed = 0;
            if(verbose) printf("\n. [Optimisation step] => GTR parameters...\n");
            
            For(i,6) init_values[i] = tree->mod->gtr_param[i];
            
            tree->mod->update_eigen = 1;

            /* BFGS(tree,tree->mod->gtr_param,6,1.e-5,1.e-7,&Return_Abs_Lk,&Num_Derivative_Several_Param,&Lnsrch_GTR,&failed); */
            BFGS(tree,tree->mod->gtr_param,
                 6,
                 1.e-5,1.e-3,
                 &Return_Abs_Lk,
                 &Num_Derivative_Several_Param,
                 &Lnsrch,&failed);

            if(failed)
                {
                    
                    printf("\n. Optimising one-by-one...\n");
                    
                    For(i,6) tree->mod->gtr_param[i] = init_values[i];
                    Lk(tree);
                    
                    For(i,5) 
                        Optimize_Single_Param_Generic(tree,&(tree->mod->gtr_param[i]),tree->mod->gtr_param[i],1.E-10,1.E+10,1000);
                }
            tree->mod->update_eigen = 0;
        }
    
    if(tree->mod->s_opt->opt_pinvar)
        {
            if(verbose) printf("\n. [Optimisation step] => proportion of invariable sites (pinv)...\n"); fflush(stdout);
            tree->mod->pinvar = 0.5;
            Optimize_Single_Param_Generic(tree,&(tree->mod->pinvar),tree->mod->pinvar,0.0001,0.9999,100);
        }
    
    if(tree->mod->s_opt->opt_kappa) 
        {
            if(verbose) printf("\n. [Optimisation step] => ts/tv ratio (kappa)...\n"); fflush(stdout);      
            Optimize_Single_Param_Generic(tree,&(tree->mod->kappa),tree->mod->kappa,0.1,1000.,100);
        }
    
    if(tree->mod->s_opt->opt_lambda) 
        {
            Optimize_Single_Param_Generic(tree,&(tree->mod->lambda),tree->mod->lambda,0.01,1000.,100);
        }
    
    if(tree->mod->s_opt->opt_alpha) 
        { 
            if(verbose) printf("\n. [Optimisation step] => gamma shape parameter (alpha)...\n"); fflush(stdout);
            Optimize_Single_Param_Generic(tree,&(tree->mod->alpha),tree->mod->alpha,0.001,100.,100);
        }
    
    if(tree->mod->s_opt->opt_bfreq)
        {
            
            failed = 0;
            tree->mod->update_eigen = 1;
            if(verbose) printf("\n. [Optimisation step] => nucleotide frequencies...\n");
            
            For(i,4) init_values[i] = tree->mod->pi[i];
            
            BFGS(tree,tree->mod->pi,
                 4,
                 1.e-5,1.e-3,
                 &Return_Abs_Lk,
                 &Num_Derivative_Several_Param,
                 &Lnsrch,&failed);

            if(failed)
                {
                    printf("\n. Optimising one-by-one...\n");
                    
                    For(i,4) tree->mod->pi[i] = init_values[i];
                    Lk(tree);
                    
                    For(i,4) 
                        Optimize_Single_Param_Generic(tree,&(tree->mod->pi[i]),tree->mod->pi[i],0.0001,0.99999,100);
                }
            
            tree->mod->update_eigen = 0;
        }
    

    if(tree->mod->model_number > 19)
      {
	switch(tree->mod->switch_modelname) 
	  {
	  case NO_SWITCH :
	    {
	      switch(tree->mod->subst_modelname)
		{
		case M2 : case M2a : case M3 : case MX :
		  {
		    printf("\n. [Optimisation step] => equilibrium frequencies of selection regimes (p0, p1, p2)...\n");
		    

		    For(i,tree->mod->n_catq) 
		      init_values[i] = tree->mod->qmat_struct[0]->qmat_proba[i];

		    /* failed = NO; */
                    /* BFGS(tree,tree->mod->qmat_struct[0]->trans_qmat_proba, */
                    /*      tree->mod->n_catq, */
                    /*      1.e-5,1.e-3, */
                    /*      &Return_Abs_Lk, */
                    /*      &Num_Derivative_Several_Param, */
                    /*      &Lnsrch,&failed); */
                    
		    failed = YES;
                    if(failed)
                      {
                        For(i,tree->mod->n_catq) 
                          Optimize_Single_Param_Generic(tree,
                                                        &(tree->mod->qmat_struct[0]->trans_qmat_proba[i]),
                                                        tree->mod->qmat_struct[0]->trans_qmat_proba[i],
                                                        1.E-2,
                                                        1.E+2,
                                                        20);
                      }

		    break;
		  }
		case M1 : case M1a : 
		  {
		    printf("\n. [Optimisation step] => equilibrium frequencies of selection regimes (p0, p1)...\n");
		    
		    Optimize_Single_Param_Generic(tree,
						  &(tree->mod->qmat_struct[0]->trans_qmat_proba[0]),
						  tree->mod->qmat_struct[0]->trans_qmat_proba[0],
						  1.E-2,
						  1.E+2,
						  20);
		    break;
		  }
		case M0 : break;
		default : 
		  {
		    printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
		    Exit("");
		  }
		}
	      
	      if(tree->mod->s_opt->opt_omega)
		{
		  switch(tree->mod->subst_modelname)
		    {
		    case M0 :
		      {
			printf("\n. [Optimisation step] => dn/ds ratio, class of free w... \n");
			Optimize_Single_Param_Generic(tree,
						      &(tree->mod->qmat_struct[0]->omega[0]),
						      tree->mod->qmat_struct[0]->omega[0],
						      1.E-3,
						      20.0,
						      20);
			break;
		      }
		    case M1a :
		      {
			printf("\n. [Optimisation step] => dn/ds ratio, class of free w... \n");
			Optimize_Single_Param_Generic(tree,
						      &(tree->mod->qmat_struct[0]->omega[0]),
						      tree->mod->qmat_struct[0]->omega[0],
						      1.E-3,
						      1.00,
						      20);
			break;
		      }
		      
		    case M2 :
		      {
			printf("\n. [Optimisation step] => dn/ds ratio, class of free w... \n");
			Optimize_Single_Param_Generic(tree,
						      &(tree->mod->qmat_struct[0]->omega[2]),
						      tree->mod->qmat_struct[0]->omega[2],
						      1.E-3,
						      20.0,
						      20);
/* 			Sort_Categories(tree); */
			break;
		      }
		      
		    case M2a :
		      {
			
			printf("\n. [Optimisation step] => dn/ds ratios (w0, w1, w2)... \n");
			
			For(i,3)
			  if((tree->mod->qmat_struct[0]->omega[i] < 1.0 - MDBL_MIN) ||
			     (tree->mod->qmat_struct[0]->omega[i] > 1.0 + MDBL_MIN))
			    Optimize_Single_Param_Generic(tree,
							  &(tree->mod->qmat_struct[0]->omega[i]),
							  tree->mod->qmat_struct[0]->omega[i],
							  1.E-3,
							  20.0,
							  20);
			Sort_Categories(tree);
			break;
		      }
		      
		    case M3 : case MX : 
		      {
			tree->mod->s_opt->sort_omega = 0;
			failed = 0;
			
			printf("\n. [Optimisation step] => dn/ds ratios... \n");
			
			For(i,tree->mod->n_catq) 
			  init_values[i] = tree->mod->qmat_struct[0]->omega[i];
			
                        /* failed = NO; */
                        /* BFGS(tree,tree->mod->qmat_struct[0]->omega, */
                        /*      tree->mod->n_catq, */
                        /*      1.e-5,1.e-3, */
                        /*      &Return_Abs_Lk, */
                        /*      &Num_Derivative_Several_Param, */
                        /*      &Lnsrch,&failed); */
 
                        failed = YES;
                        if(failed)
                          {
                            For(i,tree->mod->n_catq) tree->mod->qmat_struct[0]->omega[i] = init_values[i];
                            Lk(tree);

                            For(i,tree->mod->n_catq)
                              {
                                Optimize_Single_Param_Generic(tree,
                                                              &(tree->mod->qmat_struct[0]->omega[i]),
                                                              tree->mod->qmat_struct[0]->omega[i],
                                                              /* tree->mod->qmat_struct[0]->omega_min[i], */
                                                              /* tree->mod->qmat_struct[0]->omega_max[i], */
                                                              0.001,
                                                              20.,
                                                              20);
                              }
                          }

/* 			tree->numerical_param->currently_opt = tree->numerical_param->omega_num; */
/* 			Conjugate_Gradients(tree->mod->qmat_struct[0]->omega, */
/* 					    3, */
/* 					    tree); */
						

                        if(tree->mod->s_opt->opt_p_omega && tree->mod->s_opt->opt_omega)
                          {
                            printf("\n. Multidimensional optimisation step (dN/dS and corresponding frequencies)...\n");                  
                            BFGS(tree,tree->mod->qmat_struct[0]->omega,
                                 2*tree->mod->n_catq,
                                 1.e-5,1.e-3,
                                 &Return_Abs_Lk,
                                 &Num_Derivative_Several_Param,
                                 &Lnsrch,&failed);
                          }
                        
			tree->mod->s_opt->sort_omega = 1;
			Sort_Categories(tree);
			break;
		      }
		    default : 
		      {
			printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
			Exit("");		  
		      }
		    }
		}
	      break;
	    }
	  case SWITCH_S1 : case SWITCH_S2 :
	    {
	      if(tree->mod->s_opt->opt_p_omega)
		{                  

		  printf("\n. [Optimisation step] => equilibrium frequencies of selection regimes (p0, p1, p2)...\n");
		  
		  For(i,tree->mod->n_omega) 
		    init_values[i] = tree->mod->qmat_struct[0]->trans_omega_proba[i];
		  
                  /* failed = NO; */
                  /* BFGS(tree,tree->mod->qmat_struct[0]->trans_omega_proba, */
                  /*      tree->mod->n_omega, */
                  /*      1.e-5,1.e-3, */
                  /*      &Return_Abs_Lk, */
                  /*      &Num_Derivative_Several_Param, */
                  /*      &Lnsrch,&failed); */
                  
                  failed = YES;
                  if(failed)
                    {
                      For(i,tree->mod->n_omega) tree->mod->qmat_struct[0]->trans_omega_proba[i] = init_values[i];
                      Lk(tree);

                      For(i,tree->mod->n_omega) 
                        Optimize_Single_Param_Generic(tree,
                                                      &(tree->mod->qmat_struct[0]->trans_omega_proba[i]),
                                                      tree->mod->qmat_struct[0]->trans_omega_proba[i],
                                                      1.E-1,
                                                      1.E+1,
                                                      20);                                            		    
                    }
                }

                  
	      if(tree->mod->s_opt->opt_omega)
		{
		  switch(tree->mod->subst_modelname)
		    {
		    case M2 :
		      {
			printf("\n. [Optimisation step] => dn/ds ratio, class of free w... \n");
			Optimize_Single_Param_Generic(tree,
						      &(tree->mod->qmat_struct[0]->omega[2]),
						      tree->mod->qmat_struct[0]->omega[2],
						      1.E-3,
						      20.,
						      20);
/* 			Sort_Categories(tree); */
			break;
		      }
		    case M2a :
		      {
			printf("\n. [Optimisation step] => dn/ds ratios (w0, w1, w2)... \n");
			
			For(i,3)
			  if((tree->mod->qmat_struct[0]->omega[i] < 1.0 - MDBL_MIN) ||
			     (tree->mod->qmat_struct[0]->omega[i] > 1.0 + MDBL_MIN))
			    Optimize_Single_Param_Generic(tree,
							  &(tree->mod->qmat_struct[0]->omega[i]),
							  tree->mod->qmat_struct[0]->omega[i],
							  1.E-3,
							  20.,
							  20);
			Sort_Categories(tree);
			break;
		      }
		      
		    case M3 : case MX : 
		      {			
			tree->mod->s_opt->sort_omega = 0;
			
			printf("\n. [Optimisation step] => dn/ds ratios... \n");

			For(i,tree->mod->n_omega) init_values[i] = tree->mod->qmat_struct[0]->omega[i];

                        /* failed = NO; */
                        /* BFGS(tree,tree->mod->qmat_struct[0]->omega, */
                        /*      tree->mod->n_omega, */
                        /*      1.e-5,1.e-3, */
                        /*      &Return_Abs_Lk, */
                        /*      &Num_Derivative_Several_Param, */
                        /*      &Lnsrch,&failed); */
                        
                        failed = YES;
                        if(failed)
                          {
                            For(i,tree->mod->n_omega) tree->mod->qmat_struct[0]->omega[i] = init_values[i];
                            Lk(tree);

                            For(i,tree->mod->n_omega)
                              {
                                Optimize_Single_Param_Generic(tree,
                                                              &(tree->mod->qmat_struct[0]->omega[i]),
                                                              tree->mod->qmat_struct[0]->omega[i],
                                                              /* tree->mod->qmat_struct[0]->omega_min[i], */
                                                              /* tree->mod->qmat_struct[0]->omega_max[i], */
                                                              0.001,
                                                              20.,
                                                              20);
                              }
                          }

                        if(tree->mod->s_opt->opt_p_omega && tree->mod->s_opt->opt_omega)
                          {
                            printf("\n. Multidimensional optimisation step (dN/dS and corresponding frequencies)...\n");
                            BFGS(tree,tree->mod->qmat_struct[0]->omega,
                                 2*tree->mod->n_omega,
                                 1.e-5,1.e-3,
                                 &Return_Abs_Lk,
                                 &Num_Derivative_Several_Param,
                                 &Lnsrch,&failed);
                          }

			tree->mod->s_opt->sort_omega = 1;
			Sort_Categories(tree);       
			break;
		      }
		      
		    case M0 :
		      {
			printf("\n. [Optimisation step] => dn/ds ratio, class of free w... \n");
			Optimize_Single_Param_Generic(tree,
						      &(tree->mod->qmat_struct[0]->omega[0]),
						      tree->mod->qmat_struct[0]->omega[0],
						      1.E-3,
						      20.,
						      20);
			break;
		      }
		    default : Exit("\n= Err. in switch case subst_modelname 3\n"); 
		    }
		  
		  if(tree->mod->s_opt->opt_theta == 1)
		    {
		      printf("\n. [Optimisation step] => switching rate parameter (delta)... \n");
		      
		      Optimize_Single_Param_Generic(tree,
						    &(tree->mod->qmat_struct[0]->theta[0]),
						    tree->mod->qmat_struct[0]->theta[0],
						    0.01,
						    500,
						    20);
		    }
		  else if(tree->mod->s_opt->opt_theta > 1)
		    {
		      
		      printf("\n. [Optimisation step] => switching rate parameters... \n");
		      
		      For(i,tree->mod->n_omega*(tree->mod->n_omega-1)/2) init_values[i] = tree->mod->qmat_struct[0]->theta[i];
		      
		      For(i,tree->mod->n_omega*(tree->mod->n_omega-1)/2)
			Optimize_Single_Param_Generic(tree,
						      &(tree->mod->qmat_struct[0]->theta[i]),
						      tree->mod->qmat_struct[0]->theta[i],
						      1.E-3,
						      500.,
						      20);
		      
		    }
		}
              break;
	    }
	  default : Exit("\n= Err. in switch case switch_modelname\n");
	  }
      }
    if(tree->mod->s_opt->opt_subst_rate)
      {
	printf("\n. Optimization of the substitution rate...\n");
	Optimize_Single_Param_Generic(tree,
				      &(tree->mod->subst_rate),
				      tree->mod->subst_rate,
				      1.E-10,
				      100.,
				      2000);
	printf("\n. Substitution rate = %f\n",tree->mod->subst_rate);
      }
    
    tree->both_sides         = init_both_sides;
    tree->mod->s_opt->opt_bl = init_derivatives;
    
    Free(init_values);
}

/*********************************************************/


#define ITMAX 200
#define EPS   3.0e-8
#define TOLX (4*EPS)
#define STPMX 100.0
static fit_double sqrarg;
#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)

void BFGS(arbre *tree, 
	  fit_double *p, 
	  int n, 
	  fit_double gtol, 
	  fit_double step_size,
	  fit_double(*func)(arbre *tree), 
	  int(*dfunc)(arbre *tree,fit_double *param,int n_param,fit_double stepsize,fit_double(*func)(arbre *tree),fit_double *derivatives), 
	  int(*lnsrch)(arbre *tree, int n, fit_double *xold, fit_double fold,fit_double *g, fit_double *p, fit_double *x,fit_double *f, fit_double stpmax, int *check),
	  int *failed)
{

  int check,i,its,j;
  fit_double den,fac,fad,fae,fp,stpmax,sum=0.0,sumdg,sumxi,temp,test,fret;
  fit_double *dg,*g,*hdg,**hessin,*pnew,*xi;
  hessin = (fit_double **)mCalloc(n,sizeof(fit_double *));
  For(i,n) hessin[i] = (fit_double *)mCalloc(n,sizeof(fit_double));
  dg   = (fit_double *)mCalloc(n,sizeof(fit_double ));
  g    = (fit_double *)mCalloc(n,sizeof(fit_double ));
  pnew = (fit_double *)mCalloc(n,sizeof(fit_double ));
  hdg  = (fit_double *)mCalloc(n,sizeof(fit_double ));
  xi   = (fit_double *)mCalloc(n,sizeof(fit_double ));
  

  /* PhyML_Printf("\n. ENTER BFGS WITH: %f\n",Lk(NULL,tree)); */

  fp=(*func)(tree);
  (*dfunc)(tree,p,n,step_size,func,g);

  for (i=0;i<n;i++) 
    {
      for (j=0;j<n;j++) hessin[i][j]=0.0;
      hessin[i][i]=1.0;
      xi[i] = -g[i];
      sum += p[i]*p[i];
    }

  stpmax=STPMX*MAX(SQRT(sum),(fit_double)n);

  for(its=1;its<=ITMAX;its++) 
    {
      lnsrch(tree,n,p,fp,g,xi,pnew,&fret,stpmax,&check);

      fp = fret;
      
      for (i=0;i<n;i++) 
	{
	  xi[i]=pnew[i]-p[i];
	  p[i]=pnew[i];
	}

      test=0.0;
      for (i=0;i<n;i++) 
	{
	  temp=FABS(xi[i])/MAX(FABS(p[i]),1.0);
	  if (temp > test) test=temp;
	}

      if (test < TOLX) 
	{
	  (*func)(tree);
	  For(i,n) Free(hessin[i]);
	  free(hessin);
	  free(xi);
	  free(pnew);
	  free(hdg);
	  free(g);
	  free(dg);   
	  return;
	}

      for (i=0;i<n;i++) dg[i]=g[i];

      (*dfunc)(tree,p,n,step_size,func,g);

      test=0.0;
      den=MAX(fret,1.0);
      for (i=0;i<n;i++) 
	{
	  temp=FABS(g[i])*MAX(FABS(p[i]),1.0)/den;
	  if (temp > test) test=temp;
	}


      if (test < gtol) 
	{
	  (*func)(tree);
	  For(i,n) Free(hessin[i]);
	  free(hessin);
	  free(xi);
	  free(pnew);
	  free(hdg);
	  free(g);
	  free(dg);   
	  return;
	}

    for (i=0;i<n;i++) dg[i]=g[i]-dg[i];

    for (i=0;i<n;i++) 
      {
	hdg[i]=0.0;
	for (j=0;j<n;j++) hdg[i] += hessin[i][j]*dg[j];
      }

    fac=fae=sumdg=sumxi=0.0;
    for (i=0;i<n;i++) 
      {
	fac += dg[i]*xi[i];
	fae += dg[i]*hdg[i];
	sumdg += SQR(dg[i]);
	sumxi += SQR(xi[i]);
      }
    
    if(fac*fac > EPS*sumdg*sumxi) 
      {
	fac=1.0/fac;
	fad=1.0/fae;
	for (i=0;i<n;i++) dg[i]=fac*xi[i]-fad*hdg[i];
	for (i=0;i<n;i++) 
	  {
	    for (j=0;j<n;j++) 
	      {
		hessin[i][j] += fac*xi[i]*xi[j]
		  -fad*hdg[i]*hdg[j]+fae*dg[i]*dg[j];
	      }
	  }
      }
    for (i=0;i<n;i++) 
      {
	xi[i]=0.0;
	for (j=0;j<n;j++) xi[i] -= hessin[i][j]*g[j];
      }
    }

  printf("\n. Too many iterations in BFGS...\n");
  *failed = 1;
  For(i,n) Free(hessin[i]);
  free(hessin);
  free(xi);
  free(pnew);
  free(hdg);
  free(g);
  free(dg);   
}

#undef ITMAX
#undef EPS
#undef TOLX
#undef STPMX

/*********************************************************/

#define ALF 1.0e-4
#define TOLX 1.0e-7

int Lnsrch(arbre *tree, int n, fit_double *xold, fit_double fold, fit_double *g, fit_double *p, fit_double *x,
	   fit_double *f, fit_double stpmax, int *check)
{
    //     Adapted from 'Numerical Recipes in C, Second Edition, Press et al., 1988'
    //     Line search for the optimisation of the GTR model parameters.
 

  int i;
  fit_double a,alam,alam2,alamin,b,disc,f2,fold2,rhs1,rhs2,slope,sum,temp,test,tmplam;
  fit_double *local_xold;

  local_xold = (fit_double *)mCalloc(n,sizeof(fit_double));
  For(i,n) local_xold[i] = xold[i];

  alam = alam2 = f2 = fold2 = tmplam = .0;

  *check=0;
  for(sum=0.0,i=0;i<n;i++) sum += p[i]*p[i];
  sum=(fit_double)sqrt(sum);
  if(sum > stpmax)
    for(i=0;i<n;i++) p[i] *= stpmax/sum;
  for(slope=0.0,i=0;i<n;i++)
    slope += g[i]*p[i];
  test=0.0;
  for(i=0;i<n;i++) 
    {
      temp=(fit_double)fabs(p[i])/MMAX((fit_double)fabs(xold[i]),1.0);
      if (temp > test) test=temp;
    }
  alamin=TOLX/test;
  alam=1.0;
  for (;;) 
    {
      for(i=0;i<n;i++) 
          {
              /* x[i]=xold[i]+alam*p[i]; */

            x[i]=FABS(local_xold[i]+alam*p[i]);
            //
            xold[i] = x[i];
          }
      if(i==n) 
	{
	  *f=Return_Abs_Lk(tree);
	}
      else     
          *f=1.+fold+ALF*alam*slope;
      if (alam < alamin)
	{
	  *check=1;
	  For(i,n) xold[i] = local_xold[i];
	  Free(local_xold);
	  return 0;
	} 
      else if (*f <= fold+ALF*alam*slope) 
          {
            For(i,n) xold[i] = local_xold[i];
            Free(local_xold); 
            return 0;
          }
      else 
	{
            /* if (alam == 1.0) */
	  if ((alam < 1.0+SMALL) && (alam > 1.0-SMALL))
                tmplam = -slope/(2.0*(*f-fold-slope));
            else 
                {
                    rhs1 = *f-fold-alam*slope;
                    rhs2=f2-fold2-alam2*slope;
                    a=(rhs1/(alam*alam)-rhs2/(alam2*alam2))/(alam-alam2);
                    b=(-alam2*rhs1/(alam*alam)+alam*rhs2/(alam2*alam2))/(alam-alam2);
                    if (a == 0.0) tmplam = -slope/(2.0*b);
                    else 
                        {
                            disc=b*b-3.0*a*slope;
                            if (disc<0.0) 
                                {
                                    disc=b*b-3.0*a*slope;
                                    if (disc<0.0) tmplam = 0.5*alam;
                                    else if(b <= 0.0) tmplam=(-b+(fit_double)sqrt(disc))/(3.0*a);
                                    else tmplam = -slope/(b+(fit_double)sqrt(disc));
                                }
                            else tmplam=(-b+(fit_double)sqrt(disc))/(3.0*a);
                        }
                    if (tmplam>0.5*alam) tmplam=0.5*alam;
                }
	}
      alam2=alam;
      f2 = *f;
      fold2=fold;
      alam=MMAX(tmplam,0.1*alam);
    }
}

#undef ALF
#undef TOLX
#undef NRANSI

/*********************************************************/

/* #define ITMAX 200 */
/* #define EPS   3.0e-8 */
/* #define TOLX (4*EPS) */
/* #define STPMX 100.0 */
/* static fit_double sqrarg; */
/* #define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg) */

/* void BFGS(arbre *tree, fit_double *p, int n, fit_double gtol, fit_double step_size, */
/* 	  fit_double(*func)(), void (*dfunc)(), void (*lnsrch)(),int *failed) */
/* { */

  
/* //   Adapted from 'Numerical Recipes in C, Second Edition, Press et al., 1988' */

/* //   Broyden-Fletcher-Golfarb-Shanno quasi-Newton method for optimisation */
/* //   in multidimensions. */
  

/*   int check,i,its,j; */
/*   fit_double den,fac,fad,fae,fp,stpmax,sum=0.0,sumdg,sumxi,temp,test,fret; */
/*   fit_double *dg,*g,*hdg,**hessin,*pnew,*xi; */
/*   fit_double curr_abs_lk, init_abs_lk; */

/*   hessin = (fit_double **)mCalloc(n,sizeof(fit_double *)); */
/*   For(i,n) hessin[i] = (fit_double *)mCalloc(n,sizeof(fit_double)); */
/*   dg   = (fit_double *)mCalloc(n,sizeof(fit_double )); */
/*   g    = (fit_double *)mCalloc(n,sizeof(fit_double )); */
/*   pnew = (fit_double *)mCalloc(n,sizeof(fit_double )); */
/*   hdg  = (fit_double *)mCalloc(n,sizeof(fit_double )); */
/*   xi   = (fit_double *)mCalloc(n,sizeof(fit_double )); */
  
/*   fp=(*func)(tree); */
/*   curr_abs_lk = init_abs_lk = fp; */
/*   (*dfunc)(tree,p,n,step_size,func,g); */

/*   for (i=0;i<n;i++)  */
/*     { */
/*       for (j=0;j<n;j++) hessin[i][j]=0.0; */
/*       hessin[i][i]=1.0; */
/*       xi[i] = -g[i]; */
/*       sum += p[i]*p[i]; */
/*     } */

/*   stpmax=STPMX*MMAX((fit_double)sqrt(sum),(fit_double)n); */

/*   for(its=1;its<=ITMAX;its++)  */
/*     { */
/*         lnsrch(tree,n,p,fp,g,xi,pnew,&fret,stpmax,&check,func); */

/*       fp = fret; */
      
/*       for (i=0;i<n;i++)  */
/* 	{ */
/* 	  xi[i]=pnew[i]-p[i]; */
/* 	  p[i]=pnew[i]; */
/* 	} */

/*       test=0.0; */
/*       for (i=0;i<n;i++)  */
/* 	{ */
/* 	  temp=(fit_double)fabs(xi[i])/MMAX((fit_double)fabs(p[i]),1.0); */
/* 	  if (temp > test) test=temp; */
/* 	} */
/*       if (test < TOLX)  */
/* 	{ */
/* 	  curr_abs_lk = (*func)(tree); */
/* 	  For(i,n) Free(hessin[i]); */
/* 	  free(hessin); */
/* 	  free(xi); */
/* 	  free(pnew); */
/* 	  free(hdg); */
/* 	  free(g); */
/* 	  free(dg); */
/* 	  if((curr_abs_lk > init_abs_lk-MIN_DIFF_LK) || (its == 1))  */
/*               { */
/*                   printf("\n. WARNING : BFGS failed ! \n"); */
/*                   *failed = 1; */
/*               } */
/* 	  return; */
/* 	} */

/*       for (i=0;i<n;i++) dg[i]=g[i]; */

/*       (*dfunc)(tree,p,n,step_size,func,g); */

/*       test=0.0; */
/*       den=MMAX(fret,1.0); */
/*       for (i=0;i<n;i++)  */
/* 	{ */
/* 	  temp=(fit_double)fabs(g[i])*MMAX((fit_double)fabs(p[i]),1.0)/den; */
/* 	  if (temp > test) test=temp; */
/* 	} */
/*       if (test < gtol)  */
/* 	{ */
/*             curr_abs_lk = (*func)(tree); */
/*             For(i,n) Free(hessin[i]); */
/*             free(hessin); */
/*             free(xi); */
/*             free(pnew); */
/*             free(hdg); */
/*             free(g); */
/*             free(dg); */
/*             if(curr_abs_lk > init_abs_lk-MIN_DIFF_LK) */
/*                 printf("\n. WARNING : BFGS failed ! \n"); */
/*             return; */
/* 	} */

/*     for (i=0;i<n;i++) dg[i]=g[i]-dg[i]; */

/*     for (i=0;i<n;i++)  */
/*       { */
/* 	hdg[i]=0.0; */
/* 	for (j=0;j<n;j++) hdg[i] += hessin[i][j]*dg[j]; */
/*       } */

/*     fac=fae=sumdg=sumxi=0.0; */
/*     for (i=0;i<n;i++)  */
/*       { */
/* 	fac += dg[i]*xi[i]; */
/* 	fae += dg[i]*hdg[i]; */
/* 	sumdg += SQR(dg[i]); */
/* 	sumxi += SQR(xi[i]); */
/*       } */
    
/*     if(fac*fac > EPS*sumdg*sumxi)  */
/*       { */
/* 	fac=1.0/fac; */
/* 	fad=1.0/fae; */
/* 	for (i=0;i<n;i++) dg[i]=fac*xi[i]-fad*hdg[i]; */
/* 	for (i=0;i<n;i++)  */
/* 	  { */
/* 	    for (j=0;j<n;j++)  */
/* 	      { */
/* 		hessin[i][j] += fac*xi[i]*xi[j] */
/* 		  -fad*hdg[i]*hdg[j]+fae*dg[i]*dg[j]; */
/* 	      } */
/* 	  } */
/*       } */
/*     for (i=0;i<n;i++)  */
/*       { */
/* 	xi[i]=0.0; */
/* 	for (j=0;j<n;j++) xi[i] -= hessin[i][j]*g[j]; */
/*       } */
/*     } */

/*   (*func)(tree); */
/*   printf("\n. Too many iterations in BFGS...\n"); */
/*   *failed = 1; */
/*   For(i,n) Free(hessin[i]); */
/*   free(hessin); */
/*   free(xi); */
/*   free(pnew); */
/*   free(hdg); */
/*   free(g); */
/*   free(dg);    */
/* } */

/* #undef ITMAX */
/* #undef EPS */
/* #undef TOLX */
/* #undef STPMX */

/*********************************************************/

#define ALF 1.0e-4
#define TOLX 1.0e-7

void Lnsrch_GTR(arbre *tree, int n, fit_double *xold, fit_double fold, fit_double *g, fit_double *p, fit_double *x,
		fit_double *f, fit_double stpmax, int *check, fit_double(*func)())
{
  /* 
     Adapted from 'Numerical Recipes in C, Second Edition, Press et al., 1988'

     Line search for the optimisation of the GTR model parameters.
  */

  int i;
  fit_double a,alam,alam2,alamin,b,disc,f2,fold2,rhs1,rhs2,slope,sum,temp,test,tmplam;
  fit_double *local_xold;

  alam = alam2 = f2 = fold2 = tmplam = .0;

  local_xold = (fit_double *)mCalloc(n,sizeof(fit_double));
  For(i,n) local_xold[i] = xold[i];


  *check=0;
  for(sum=0.0,i=0;i<n;i++) sum += p[i]*p[i];
  sum=(fit_double)sqrt(sum);
  if(sum > stpmax)
    for(i=0;i<n;i++) p[i] *= stpmax/sum;
  for(slope=0.0,i=0;i<n;i++)
    slope += g[i]*p[i];
  test=0.0;
  for(i=0;i<n;i++) 
    {
      temp=(fit_double)fabs(p[i])/MMAX((fit_double)fabs(local_xold[i]),1.0);
      if (temp > test) test=temp;
    }
  alamin=TOLX/test;
  alam=1.0;
  for (;;) 
    {
      for(i=0;i<n;i++) x[i]=local_xold[i]+alam*p[i];
      /**/      
      for(i=0;i<n;i++) 
	{
	  tree->mod->gtr_param[i]=local_xold[i]+alam*p[i];
	  if(tree->mod->gtr_param[i] < 0.0) break;
	}
      /**/
      if(i==n) 
	{
            *f=(*func)(tree);
/* 	  printf("loglk = %f\n",*f); */
	}
      else     *f=1.+fold+ALF*alam*slope;
      if (alam < alamin)
	{
	  for (i=0;i<n;i++) {x[i]=local_xold[i]; if(x[i] < .0) break;}

	  /**/      
	  for(i=0;i<n;i++)
	    {
	      tree->mod->gtr_param[i]=local_xold[i]+alam*p[i];
	      if(tree->mod->gtr_param[i] < 0.0) 
		tree->mod->gtr_param[i] = 0.0;
	    }
	  /**/


	  *check=1;
	  For(i,n) xold[i] = local_xold[i];
	  Free(local_xold);
	  return;
	} 
      else if (*f <= fold+ALF*alam*slope) 
	{
	  For(i,n) xold[i] = local_xold[i];
	  Free(local_xold); 
	  return;
	}
      else 
	{
	  if (alam == 1.0)
	    tmplam = -slope/(2.0*(*f-fold-slope));
	  else 
	    {
	      rhs1 = *f-fold-alam*slope;
	      rhs2=f2-fold2-alam2*slope;
	      a=(rhs1/(alam*alam)-rhs2/(alam2*alam2))/(alam-alam2);
	      b=(-alam2*rhs1/(alam*alam)+alam*rhs2/(alam2*alam2))/(alam-alam2);
	      if (a == 0.0) tmplam = -slope/(2.0*b);
	      else 
		{
		  disc=b*b-3.0*a*slope;
		  if (disc<0.0) 
		    {
		      disc=b*b-3.0*a*slope;
		      if (disc<0.0) tmplam = 0.5*alam;
		      else if(b <= 0.0) tmplam=(-b+(fit_double)sqrt(disc))/(3.0*a);
		      else tmplam = -slope/(b+(fit_double)sqrt(disc));
		    }
		  else tmplam=(-b+(fit_double)sqrt(disc))/(3.0*a);
		}
	      if (tmplam>0.5*alam) tmplam=0.5*alam;
	    }
	}
      alam2=alam;
      f2 = *f;
      fold2=fold;
      alam=MMAX(tmplam,0.1*alam);
    }
  Free(local_xold);
}

#undef ALF
#undef TOLX
#undef NRANSI

/*********************************************************/

#define ALF 1.0e-4
#define TOLX 1.0e-7

void Lnsrch_Omega(arbre *tree, int n, fit_double *xold, fit_double fold, fit_double *g, fit_double *p, fit_double *x,
		  fit_double *f, fit_double stpmax, int *check, fit_double (*func)())
{
  /* 
     Adapted from 'Numerical Recipes in C, Second Edition, Press et al., 1988'
     
     Line search for the optimisation of dn/ds ratios.
  */

  int i;
  fit_double a,alam,alam2,alamin,b,disc,f2,fold2,rhs1,rhs2,slope,sum,temp,test,tmplam;
  fit_double *local_xold;

  alam = alam2 = f2 = fold2 = tmplam = .0;

  local_xold = (fit_double *)mCalloc(n,sizeof(fit_double));
  For(i,n) local_xold[i] = xold[i];


  *check=0;
  for(sum=0.0,i=0;i<n;i++) sum += p[i]*p[i];
  sum=(fit_double)sqrt(sum);
  if(sum > stpmax)
    for(i=0;i<n;i++) p[i] *= stpmax/sum;
  for(slope=0.0,i=0;i<n;i++)
    slope += g[i]*p[i];
  test=0.0;
  for(i=0;i<n;i++) 
    {
      temp=(fit_double)fabs(p[i])/MMAX((fit_double)fabs(local_xold[i]),1.0);
      if (temp > test) test=temp;
    }
  alamin=TOLX/test;
  alam=1.0;
  for (;;) 
    {
      for(i=0;i<n;i++) x[i]=(fit_double)fabs(local_xold[i]+alam*p[i]);
      /**/      
      for(i=0;i<n;i++) 
	{
	  tree->mod->qmat_struct[0]->omega[i]=(fit_double)fabs(local_xold[i]+alam*p[i]);
/* 	  printf("%f\n",tree->mod->qmat_struct[0]->omega[i]); */
/* 	  if(tree->mod->qmat_struct[0]->omega[i] < 0.0) break; */
	}
      /**/
      if(i==n) 
	{
	  *f=(*func)(tree);
	}
      else     *f=1.+fold+ALF*alam*slope;
      if (alam < alamin)
	{
	  for (i=0;i<n;i++) x[i]=local_xold[i];
	  for (i=0;i<n;i++) tree->mod->qmat_struct[0]->omega[i]=(fit_double)fabs(local_xold[i]);
	  *check=1;
	  For(i,n) xold[i] = local_xold[i];
	  Free(local_xold);
	  return;
	} 
      else if (*f <= fold+ALF*alam*slope) 
	{
	  For(i,n) xold[i] = local_xold[i];
	  Free(local_xold); 
	  return;
	}
      else 
	{
	  if (alam == 1.0)
	    tmplam = -slope/(2.0*(*f-fold-slope));
	  else 
	    {
	      rhs1 = *f-fold-alam*slope;
	      rhs2=f2-fold2-alam2*slope;
	      a=(rhs1/(alam*alam)-rhs2/(alam2*alam2))/(alam-alam2);
	      b=(-alam2*rhs1/(alam*alam)+alam*rhs2/(alam2*alam2))/(alam-alam2);
	      if (a == 0.0) tmplam = -slope/(2.0*b);
	      else 
		{
		  disc=b*b-3.0*a*slope;
		  if (disc<0.0) 
		    {
		      disc=b*b-3.0*a*slope;
		      if (disc<0.0) tmplam = 0.5*alam;
		      else if(b <= 0.0) tmplam=(-b+(fit_double)sqrt(disc))/(3.0*a);
		      else tmplam = -slope/(b+(fit_double)sqrt(disc));
		    }
		  else tmplam=(-b+(fit_double)sqrt(disc))/(3.0*a);
		}
	      if (tmplam>0.5*alam) tmplam=0.5*alam;
	    }
	}
      alam2=alam;
      f2 = *f;
      fold2=fold;
      alam=MMAX(tmplam,0.1*alam);
    }
  Free(local_xold);
}

#undef ALF
#undef TOLX
#undef NRANSI

/*********************************************************/

#define ALF 1.0e-4
#define TOLX 1.0e-7

void Lnsrch_Qmat_Proba(arbre *tree, int n, fit_double *xold, fit_double fold, fit_double *g, fit_double *p, fit_double *x,
		       fit_double *f, fit_double stpmax, int *check, fit_double (*func)())
{
  /* 
     Adapted from 'Numerical Recipes in C, Second Edition, Press et al., 1988'
     
     Line search for the optimisation of the prior probabilities of the dn/ds ratios
     (without switches).
  */
  
  fit_double a,alam,alam2,alamin,b,disc,f2,fold2,rhs1,rhs2,slope,sum,temp,test,tmplam;
  fit_double *local_xold;
  int i;
  
  alam = alam2 = f2 = fold2 = tmplam = .0;

  local_xold = (fit_double *)mCalloc(n,sizeof(fit_double));
  For(i,n) local_xold[i] = xold[i];
  
  
  *check=0;
  for(sum=0.0,i=0;i<n;i++) sum += p[i]*p[i];
  sum=(fit_double)sqrt(sum);
  if(sum > stpmax)
    for(i=0;i<n;i++) p[i] *= stpmax/sum;
  for(slope=0.0,i=0;i<n;i++)
    slope += g[i]*p[i];
  test=0.0;
  for(i=0;i<n;i++) 
    {
      temp=(fit_double)fabs(p[i])/MMAX((fit_double)fabs(local_xold[i]),1.0);
      if (temp > test) test=temp;
    }
  alamin=TOLX/test;
  alam=1.0;
  for (;;) 
    {
      sum = .0;
      for(i=0;i<n;i++) x[i]=(fit_double)fabs(local_xold[i]+alam*p[i]);
/*       for(i=0;i<n;i++) sum += x[i]; */
/*       for(i=0;i<n;i++) x[i]/=sum; */
      for(i=0;i<n;i++) tree->mod->qmat_struct[0]->qmat_proba[i] = x[i];

      if(i==n) 
	{
	  *f=(*func)(tree);
	}
      else     *f=1.+fold+ALF*alam*slope;
      if (alam < alamin)
	{
	  for (i=0;i<n;i++) x[i]=local_xold[i];
	  for (i=0;i<n;i++) tree->mod->qmat_struct[0]->qmat_proba[i]=(fit_double)fabs(local_xold[i]);
	  *check=1;
	  For(i,n) xold[i] = local_xold[i];
	  Free(local_xold);
	  return;
	} 
      else if (*f <= fold+ALF*alam*slope) 
	{
	  For(i,n) xold[i] = local_xold[i];
	  Free(local_xold); 
	  return;
	}
      else 
	{
	  if (alam == 1.0)
	    tmplam = -slope/(2.0*(*f-fold-slope));
	  else 
	    {
	      rhs1 = *f-fold-alam*slope;
	      rhs2=f2-fold2-alam2*slope;
	      a=(rhs1/(alam*alam)-rhs2/(alam2*alam2))/(alam-alam2);
	      b=(-alam2*rhs1/(alam*alam)+alam*rhs2/(alam2*alam2))/(alam-alam2);
	      if (a == 0.0) tmplam = -slope/(2.0*b);
	      else 
		{
		  disc=b*b-3.0*a*slope;
		  if (disc<0.0) 
		    {
		      disc=b*b-3.0*a*slope;
		      if (disc<0.0) tmplam = 0.5*alam;
		      else if(b <= 0.0) tmplam=(-b+(fit_double)sqrt(disc))/(3.0*a);
		      else tmplam = -slope/(b+(fit_double)sqrt(disc));
		    }
		  else tmplam=(-b+(fit_double)sqrt(disc))/(3.0*a);
		}
	      if (tmplam>0.5*alam) tmplam=0.5*alam;
	    }
	}
      alam2=alam;
      f2 = *f;
      fold2=fold;
      alam=MMAX(tmplam,0.1*alam);
    }
  Free(local_xold);
}

#undef ALF
#undef TOLX
#undef NRANSI

/*********************************************************/

#define ALF 1.0e-4
#define TOLX 1.0e-7

void Lnsrch_Omega_Proba(arbre *tree, int n, fit_double *xold, fit_double fold, fit_double *g, fit_double *p, fit_double *x,
			fit_double *f, fit_double stpmax, int *check, fit_double (*func)())
{
  /* 
     Adapted from 'Numerical Recipes in C, Second Edition, Press et al., 1988'
     
     Line search for the optimisation of the prior probabilities of the dn/ds ratios
     (with switches).
  */

  int i;
  fit_double a,alam,alam2,alamin,b,disc,f2,fold2,rhs1,rhs2,slope,sum,temp,test,tmplam;
  fit_double *local_xold;

  alam = alam2 = f2 = fold2 = tmplam = .0;

  local_xold = (fit_double *)mCalloc(n,sizeof(fit_double));
  For(i,n) local_xold[i] = xold[i];


  *check=0;
  for(sum=0.0,i=0;i<n;i++) sum += p[i]*p[i];
  sum=(fit_double)sqrt(sum);
  if(sum > stpmax)
    for(i=0;i<n;i++) p[i] *= stpmax/sum;
  for(slope=0.0,i=0;i<n;i++)
    slope += g[i]*p[i];
  test=0.0;
  for(i=0;i<n;i++) 
    {
      temp=(fit_double)fabs(p[i])/MMAX((fit_double)fabs(local_xold[i]),1.0);
      if (temp > test) test=temp;
    }
  alamin=TOLX/test;
  alam=1.0;
  for (;;) 
    {
        sum = .0;
        for(i=0;i<n;i++) x[i]=(fit_double)fabs(local_xold[i]+alam*p[i]);
/*         for(i=0;i<n;i++) sum += x[i]; */
/*         for(i=0;i<n;i++) x[i]/=sum; */
        for(i=0;i<n;i++) tree->mod->qmat_struct[0]->omega_proba[i] = x[i];


      /**/
      if(i==n) 
	  *f=(*func)(tree);
      else     
          *f=1.+fold+ALF*alam*slope;

      if (alam < alamin)
          {
              for (i=0;i<n;i++) 
                  {
                      x[i]=local_xold[i];
                      tree->mod->qmat_struct[0]->omega_proba[i] = x[i];
                  }
              *check=1;
              For(i,n) xold[i] = local_xold[i];
              Free(local_xold);
              return;
          }

      else if (*f <= fold+ALF*alam*slope) 
	{
	  For(i,n) xold[i] = local_xold[i];
	  Free(local_xold); 
	  return;
	}
      else 
	{
	  if (alam == 1.0)
	    tmplam = -slope/(2.0*(*f-fold-slope));
	  else 
	    {
                rhs1 = *f-fold-alam*slope;
                rhs2=f2-fold2-alam2*slope;
                a=(rhs1/(alam*alam)-rhs2/(alam2*alam2))/(alam-alam2);
                b=(-alam2*rhs1/(alam*alam)+alam*rhs2/(alam2*alam2))/(alam-alam2);
                if (a == 0.0) tmplam = -slope/(2.0*b);
                else 
                    {
                        disc=b*b-3.0*a*slope;
                        if (disc<0.0) 
                            {
                                disc=b*b-3.0*a*slope;
                                if (disc<0.0) tmplam = 0.5*alam;
                                else if(b <= 0.0) tmplam=(-b+(fit_double)sqrt(disc))/(3.0*a);
                                else tmplam = -slope/(b+(fit_double)sqrt(disc));
                            }
                        else tmplam=(-b+(fit_double)sqrt(disc))/(3.0*a);
                    }
                if (tmplam>0.5*alam) tmplam=0.5*alam;
	    }
	}
      alam2=alam;
      f2 = *f;
      fold2=fold;
      alam=MMAX(tmplam,0.1*alam);
    }
  Free(local_xold);
}

#undef ALF
#undef TOLX
#undef NRANSI


/*********************************************************/

#define ALF 1.0e-4
#define TOLX 1.0e-7

void Lnsrch_Trans_Omega_Proba(arbre *tree, int n, fit_double *xold, fit_double fold, fit_double *g, fit_double *p, fit_double *x,
                              fit_double *f, fit_double stpmax, int *check, fit_double (*func)())
{
  /* 
     Adapted from 'Numerical Recipes in C, Second Edition, Press et al., 1988'
     
     Line search for the optimisation of the prior probabilities of the dn/ds ratios
     (with switches).
  */

  int i;
  fit_double a,alam,alam2,alamin,b,disc,f2,fold2,rhs1,rhs2,slope,sum,temp,test,tmplam;
  fit_double *local_xold;

  alam = alam2 = f2 = fold2 = tmplam = .0;

  local_xold = (fit_double *)mCalloc(n,sizeof(fit_double));
  For(i,n) local_xold[i] = xold[i];


  *check=0;
  for(sum=0.0,i=0;i<n;i++) sum += p[i]*p[i];
  sum=(fit_double)sqrt(sum);
  if(sum > stpmax)
    for(i=0;i<n;i++) p[i] *= stpmax/sum;
  for(slope=0.0,i=0;i<n;i++)
    slope += g[i]*p[i];
  test=0.0;
  for(i=0;i<n;i++) 
    {
      temp=(fit_double)fabs(p[i])/MMAX((fit_double)fabs(local_xold[i]),1.0);
      if (temp > test) test=temp;
    }
  alamin=TOLX/test;
  alam=1.0;
  for (;;) 
    {
        sum = .0;
        for(i=0;i<n;i++) x[i]=(fit_double)fabs(local_xold[i]+alam*p[i]);
/*         for(i=0;i<n;i++) sum += x[i]; */
/*         for(i=0;i<n;i++) x[i]/=sum; */
        for(i=0;i<n;i++) tree->mod->qmat_struct[0]->trans_omega_proba[i] = x[i];


      /**/
      if(i==n) 
	  *f=(*func)(tree);
      else     
          *f=1.+fold+ALF*alam*slope;

      if (alam < alamin)
          {
              for (i=0;i<n;i++) 
                  {
                      x[i]=local_xold[i];
                      tree->mod->qmat_struct[0]->trans_omega_proba[i] = x[i];
                  }
              *check=1;
              For(i,n) xold[i] = local_xold[i];
              Free(local_xold);
              return;
          }

      else if (*f <= fold+ALF*alam*slope) 
	{
	  For(i,n) xold[i] = local_xold[i];
	  Free(local_xold); 
	  return;
	}
      else 
	{
	  if (alam == 1.0)
	    tmplam = -slope/(2.0*(*f-fold-slope));
	  else 
	    {
                rhs1 = *f-fold-alam*slope;
                rhs2=f2-fold2-alam2*slope;
                a=(rhs1/(alam*alam)-rhs2/(alam2*alam2))/(alam-alam2);
                b=(-alam2*rhs1/(alam*alam)+alam*rhs2/(alam2*alam2))/(alam-alam2);
                if (a == 0.0) tmplam = -slope/(2.0*b);
                else 
                    {
                        disc=b*b-3.0*a*slope;
                        if (disc<0.0) 
                            {
                                disc=b*b-3.0*a*slope;
                                if (disc<0.0) tmplam = 0.5*alam;
                                else if(b <= 0.0) tmplam=(-b+(fit_double)sqrt(disc))/(3.0*a);
                                else tmplam = -slope/(b+(fit_double)sqrt(disc));
                            }
                        else tmplam=(-b+(fit_double)sqrt(disc))/(3.0*a);
                    }
                if (tmplam>0.5*alam) tmplam=0.5*alam;
	    }
	}
      alam2=alam;
      f2 = *f;
      fold2=fold;
      alam=MMAX(tmplam,0.1*alam);
    }
  Free(local_xold);
}

#undef ALF
#undef TOLX
#undef NRANSI


/*********************************************************/

#define ALF 1.0e-4
#define TOLX 1.0e-7

void Lnsrch_Trans_Qmat_Proba(arbre *tree, int n, fit_double *xold, fit_double fold, fit_double *g, fit_double *p, fit_double *x,
                             fit_double *f, fit_double stpmax, int *check, fit_double (*func)())
{
  /* 
     Adapted from 'Numerical Recipes in C, Second Edition, Press et al., 1988'
     
     Line search for the optimisation of the prior probabilities of the dn/ds ratios
     (with switches).
  */

  int i;
  fit_double a,alam,alam2,alamin,b,disc,f2,fold2,rhs1,rhs2,slope,sum,temp,test,tmplam;
  fit_double *local_xold;

  alam = alam2 = f2 = fold2 = tmplam = .0;

  local_xold = (fit_double *)mCalloc(n,sizeof(fit_double));
  For(i,n) local_xold[i] = xold[i];


  *check=0;
  for(sum=0.0,i=0;i<n;i++) sum += p[i]*p[i];
  sum=(fit_double)sqrt(sum);
  if(sum > stpmax)
    for(i=0;i<n;i++) p[i] *= stpmax/sum;
  for(slope=0.0,i=0;i<n;i++)
    slope += g[i]*p[i];
  test=0.0;
  for(i=0;i<n;i++) 
    {
      temp=(fit_double)fabs(p[i])/MMAX((fit_double)fabs(local_xold[i]),1.0);
      if (temp > test) test=temp;
    }
  alamin=TOLX/test;
  alam=1.0;
  for (;;) 
    {
        sum = .0;
        for(i=0;i<n;i++) x[i]=(fit_double)fabs(local_xold[i]+alam*p[i]);
        for(i=0;i<n;i++) tree->mod->qmat_struct[0]->trans_qmat_proba[i] = x[i];


      /**/
      if(i==n) 
	  *f=(*func)(tree);
      else     
          *f=1.+fold+ALF*alam*slope;

      if (alam < alamin)
          {
              for (i=0;i<n;i++) 
                  {
                      x[i]=local_xold[i];
                      tree->mod->qmat_struct[0]->trans_qmat_proba[i] = x[i];
                  }
              *check=1;
              For(i,n) xold[i] = local_xold[i];
              Free(local_xold);
              return;
          }

      else if (*f <= fold+ALF*alam*slope) 
	{
	  For(i,n) xold[i] = local_xold[i];
	  Free(local_xold); 
	  return;
	}
      else 
	{
	  if (alam == 1.0)
	    tmplam = -slope/(2.0*(*f-fold-slope));
	  else 
	    {
                rhs1 = *f-fold-alam*slope;
                rhs2=f2-fold2-alam2*slope;
                a=(rhs1/(alam*alam)-rhs2/(alam2*alam2))/(alam-alam2);
                b=(-alam2*rhs1/(alam*alam)+alam*rhs2/(alam2*alam2))/(alam-alam2);
                if (a == 0.0) tmplam = -slope/(2.0*b);
                else 
                    {
                        disc=b*b-3.0*a*slope;
                        if (disc<0.0) 
                            {
                                disc=b*b-3.0*a*slope;
                                if (disc<0.0) tmplam = 0.5*alam;
                                else if(b <= 0.0) tmplam=(-b+(fit_double)sqrt(disc))/(3.0*a);
                                else tmplam = -slope/(b+(fit_double)sqrt(disc));
                            }
                        else tmplam=(-b+(fit_double)sqrt(disc))/(3.0*a);
                    }
                if (tmplam>0.5*alam) tmplam=0.5*alam;
	    }
	}
      alam2=alam;
      f2 = *f;
      fold2=fold;
      alam=MMAX(tmplam,0.1*alam);
    }
  Free(local_xold);
}

#undef ALF
#undef TOLX
#undef NRANSI


/*********************************************************/


#define ALF 1.0e-4
#define TOLX 1.0e-7

void Lnsrch_Nucleotide_Frequencies(arbre *tree, int n, fit_double *xold, fit_double fold, fit_double *g, fit_double *p, fit_double *x,
				   fit_double *f, fit_double stpmax, int *check, fit_double (*func)())
{
  /* 
     Adapted from 'Numerical Recipes in C, Second Edition, Press et al., 1988'
     
     Line search for the optimisation of the nucleotide frequencies.
  */

  int i;
  fit_double a,alam,alam2,alamin,b,disc,f2,fold2,rhs1,rhs2,slope,sum,temp,test,tmplam;
  fit_double *local_xold;

  alam = alam2 = f2 = fold2 = tmplam = .0;

  local_xold = (fit_double *)mCalloc(n,sizeof(fit_double));
  For(i,n) local_xold[i] = xold[i];


  *check=0;
  for(sum=0.0,i=0;i<n;i++) sum += p[i]*p[i];
  sum=(fit_double)sqrt(sum);
  if(sum > stpmax)
    for(i=0;i<n;i++) p[i] *= stpmax/sum;
  for(slope=0.0,i=0;i<n;i++)
    slope += g[i]*p[i];
  test=0.0;
  for(i=0;i<n;i++) 
    {
      temp=(fit_double)fabs(p[i])/MMAX((fit_double)fabs(local_xold[i]),1.0);
      if (temp > test) test=temp;
    }
  alamin=TOLX/test;
  alam=1.0;
  for (;;) 
    {
      for(i=0;i<n;i++) x[i]=(fit_double)fabs(local_xold[i]+alam*p[i]);
      /**/      
      for(i=0;i<n;i++) 
          tree->mod->pi[i]=x[i];
      /**/
      if(i==n) 
	{
	  *f=(*func)(tree);
	}
      else     *f=1.+fold+ALF*alam*slope;
      if (alam < alamin)
	{
	  for (i=0;i<n;i++) x[i]=local_xold[i];
	  for (i=0;i<n;i++) tree->mod->pi[i]=(fit_double)fabs(local_xold[i]);
	  *check=1;
	  For(i,n) xold[i] = local_xold[i];
	  Free(local_xold);
	  return;
	} 
      else if (*f <= fold+ALF*alam*slope) 
	{
	  For(i,n) xold[i] = local_xold[i];
	  Free(local_xold); 
	  return;
	}
      else 
	{
	  if (alam == 1.0)
	    tmplam = -slope/(2.0*(*f-fold-slope));
	  else 
	    {
	      rhs1 = *f-fold-alam*slope;
	      rhs2=f2-fold2-alam2*slope;
	      a=(rhs1/(alam*alam)-rhs2/(alam2*alam2))/(alam-alam2);
	      b=(-alam2*rhs1/(alam*alam)+alam*rhs2/(alam2*alam2))/(alam-alam2);
	      if (a == 0.0) tmplam = -slope/(2.0*b);
	      else 
		{
		  disc=b*b-3.0*a*slope;
		  if (disc<0.0) 
		    {
		      disc=b*b-3.0*a*slope;
		      if (disc<0.0) tmplam = 0.5*alam;
		      else if(b <= 0.0) tmplam=(-b+(fit_double)sqrt(disc))/(3.0*a);
		      else tmplam = -slope/(b+(fit_double)sqrt(disc));
		    }
		  else tmplam=(-b+(fit_double)sqrt(disc))/(3.0*a);
		}
	      if (tmplam>0.5*alam) tmplam=0.5*alam;
	    }
	}
      alam2=alam;
      f2 = *f;
      fold2=fold;
      alam=MMAX(tmplam,0.1*alam);
    }
  Free(local_xold);
}

#undef ALF
#undef TOLX
#undef NRANSI

/*********************************************************/

#define ALF 1.0e-4
#define TOLX 1.0e-7

void Lnsrch_Theta(arbre *tree, int n, fit_double *xold, fit_double fold, fit_double *g, fit_double *p, fit_double *x,
		  fit_double *f, fit_double stpmax, int *check, fit_double (*func)())
{
  /* 
     Adapted from 'Numerical Recipes in C, Second Edition, Press et al., 1988'
     
     Line search for the optimisation of the switching parameters.
  */

  int i;
  fit_double a,alam,alam2,alamin,b,disc,f2,fold2,rhs1,rhs2,slope,sum,temp,test,tmplam;
  fit_double *local_xold;

  alam = alam2 = f2 = fold2 = tmplam = .0;

  local_xold = (fit_double *)mCalloc(n,sizeof(fit_double));
  For(i,n) local_xold[i] = xold[i];



  *check=0;
  for(sum=0.0,i=0;i<n;i++) sum += p[i]*p[i];
  sum=(fit_double)sqrt(sum);
  if(sum > stpmax)
    for(i=0;i<n;i++) p[i] *= stpmax/sum;
  for(slope=0.0,i=0;i<n;i++)
    slope += g[i]*p[i];
  test=0.0;
  for(i=0;i<n;i++) 
    {
      temp=(fit_double)fabs(p[i])/MMAX((fit_double)fabs(local_xold[i]),1.0);
      if (temp > test) test=temp;
    }
  alamin=TOLX/test;
  alam=1.0;
  for (;;) 
    {
      for(i=0;i<n;i++) x[i]=(fit_double)fabs(local_xold[i]+alam*p[i]);
      /**/      
      for(i=0;i<n;i++) 
	{
	  tree->mod->qmat_struct[0]->theta[i]=(fit_double)fabs(local_xold[i]+alam*p[i]);
/* 	  if(tree->mod->qmat_struct[0]->theta[i] < 0.0) break; */
	}
      /**/
      if(i==n) 
	{
	  *f=(*func)(tree);
	}
      else     *f=1.+fold+ALF*alam*slope;
      if (alam < alamin)
	{
	  for (i=0;i<n;i++) x[i]=local_xold[i];
	  for (i=0;i<n;i++) tree->mod->qmat_struct[0]->theta[i]=(fit_double)fabs(local_xold[i]);
	  *check=1;
	  For(i,n) xold[i] = local_xold[i];
	  Free(local_xold);
	  return;
	} 
      else if (*f <= fold+ALF*alam*slope) 
	{
	  For(i,n) xold[i] = local_xold[i];
	  Free(local_xold); 
	  return;
	}
      else 
	{
	  if (alam == 1.0)
	    tmplam = -slope/(2.0*(*f-fold-slope));
	  else 
	    {
	      rhs1 = *f-fold-alam*slope;
	      rhs2=f2-fold2-alam2*slope;
	      a=(rhs1/(alam*alam)-rhs2/(alam2*alam2))/(alam-alam2);
	      b=(-alam2*rhs1/(alam*alam)+alam*rhs2/(alam2*alam2))/(alam-alam2);
	      if (a == 0.0) tmplam = -slope/(2.0*b);
	      else 
		{
		  disc=b*b-3.0*a*slope;
		  if (disc<0.0) 
		    {
		      disc=b*b-3.0*a*slope;
		      if (disc<0.0) tmplam = 0.5*alam;
		      else if(b <= 0.0) tmplam=(-b+(fit_double)sqrt(disc))/(3.0*a);
		      else tmplam = -slope/(b+(fit_double)sqrt(disc));
		    }
		  else tmplam=(-b+(fit_double)sqrt(disc))/(3.0*a);
		}
	      if (tmplam>0.5*alam) tmplam=0.5*alam;
	    }
	}
      alam2=alam;
      f2 = *f;
      fold2=fold;
      alam=MMAX(tmplam,0.1*alam);
    }
  Free(local_xold);
}

#undef ALF
#undef TOLX
#undef NRANSI

/*********************************************************/
/*********************************************************/
/*********************************************************/


/*********************************************************/
/*********************************************************/
/*********************************************************/

fit_double Vect_Prod(fit_double *a, 
                 fit_double *b, 
                 int n)
{
    fit_double prod;
    int i;

    prod = 0.0;    
    For(i,n) prod += a[i]*b[i];

    return prod;

}

/*********************************************************/
/*********************************************************/
/*********************************************************/



fit_double Linmin_On_A_Direction(fit_double *x, 
                             fit_double *dir, 
                             fit_double n, 
                             fit_double *grad,
                             fit_double f0,
                             fit_double (*func)(),
                             arbre *tree)
{

    fit_double a,b,c,d;
    fit_double *derivatives;
    fit_double *xa,*xb,*xc,*xd;
    fit_double fxa,fxb,fxc,fxd,f0_old;
    fit_double golden;
    fit_double alpha;
    int i, n_iter;
    fit_double *init_x;

    if(tree->verbose) printf("\n. Entrering Linmin...\n");

    xa = (fit_double *)mCalloc(n,sizeof(fit_double));
    xb = (fit_double *)mCalloc(n,sizeof(fit_double));
    xc = (fit_double *)mCalloc(n,sizeof(fit_double));
    xd = (fit_double *)mCalloc(n,sizeof(fit_double));

    init_x = (fit_double *)mCalloc(n,sizeof(fit_double));
    For(i,n) init_x[i] = x[i];

    /* bracket a maximum */
    
    derivatives    = grad;
    alpha          = -1.0;    

/*     if(tree->verbose)  */
/*         { */
            For(i,n) printf("\n. Derivative[%3d]=%10g",i,derivatives[i]);
            printf("\n");
/*         } */


    alpha = 1.E-10;

    n_iter = 0;
    f0 = (*func)(x,tree);
    do
        {
            f0_old = f0;
            
            For(i,n) 
	      {
		x[i] = init_x[i] + alpha * dir[i];
		printf("x[%d]=%f\n",i,x[i]);
	      }

            f0 = (*func)(x,tree);

/*             if(tree->verbose)  */
	      printf("\n. alpha=%g f0=%.10f f0_old=%g\n",alpha,f0,f0_old);
            
            alpha *= 10.;

	    if(++n_iter > 30) break;

        }while((f0 > f0_old) || ((fit_double)fabs(f0-f0_old) < 1.E-10));

    alpha /= 10.;

    b = alpha/10.;
    d = alpha;
    a = (alpha/10.)/10.;


    For(i,n) xa[i] = init_x[i] + a * dir[i];
    For(i,n) xd[i] = init_x[i] + d * dir[i];
    For(i,n) xb[i] = init_x[i] + b * dir[i];

    fxa = (*func)(xa,tree);
    fxb = (*func)(xb,tree);
    fxd = (*func)(xd,tree);

/*     j = 1; */
/*     b = 1.0; */
/*     do */
/*         { */
/*             if(fxa>fxd) b = a+(d-a)*(1.E-3*j); */
/*             else        b = d-(d-a)*(1.E-3*j);   */

/*             For(i,n) xb[i] = init_x[i] + b * dir[i]; */
/*             if(tree->verbose) printf("\n. b=%g\n",b); */
/*             fxb = (*func)(xb,tree); */
/*             j++; */
/*             if(j > 100) break; */
/*         }while((fxa>fxb) || (fxd>fxb)); */



    golden = (3. - (fit_double)sqrt(5)) / 2.; 



    if((fxb < fxa) || (fxb < fxd)) 
        {
            printf("\n. WARNING : failed to bracket a maximum !");
            alpha = (fxa > fxb)?(a):(b);
        }
    else
        {
            /*
              
            bracket a maximum
            
            
            GOLDEN SECTION SEARCH (A FEW TIPS) 
            
            |--------|----|--------| 
            a        b    c        d 
            
            (b-a) / (d-a) = golden; 
            c-a = d-b; 
            (c-b) / (c-a) = golden 
            (c-d) / (d-b) = golden 
            
            */
            
            
            
            c = b + (d-b)*0.5;
            For(i,n) xc[i] = init_x[i] + c * dir[i];
            fxc = (*func)(xc,tree);
            
            do
                {
                    
/*                     if(tree->verbose)  */
/*                         { */
                            printf("\n. a=%g, b=%g, c=%g, d=%g\n",a,b,c,d);
                            printf("\n. fxa=%g, fxb=%g, fxc=%g, fxd=%g\n",fxa,fxb,fxc,fxd);
/*                         } */
                    
                    if(fxc > fxb)
                        {
/*                             if(tree->verbose)  */
			      printf("\n. Choosing [b,c,d]\n");
                            a = b;
                            b = c;
                            c = golden * (d-b) + b;
                            
                            fxa = fxb;
                            fxb = fxc;
                            
                            For(i,n) xa[i] = init_x[i] + a * dir[i];
                            For(i,n) xb[i] = init_x[i] + b * dir[i];
                            For(i,n) xc[i] = init_x[i] + c * dir[i];
                            For(i,n) xd[i] = init_x[i] + d * dir[i];
                            
                            fxc = (*func)(xc,tree);
                        }
                    else
                        {
/*                             if(tree->verbose)  */
			      printf("\n. Choosing [a,b,c]\n");
                            d = c;
                            c = b;
                            b = golden * (a-c) + c;
                            
                            fxd = fxc;
                            fxc = fxb;
                            
                            For(i,n) xa[i] = init_x[i] + a * dir[i];
                            For(i,n) xb[i] = init_x[i] + b * dir[i];
                            For(i,n) xc[i] = init_x[i] + c * dir[i];
                            For(i,n) xd[i] = init_x[i] + d * dir[i];
                            
                            fxb = (*func)(xb,tree);
                            
                        }
                    
                }while((fit_double)fabs(fxa-fxd)>MIN_DIFF_LK);
        
            alpha = (fxc > fxb)?(c):(b);
        }

    For(i,n) x[i] = init_x[i] + alpha * dir[i];

    tree->mod->s_opt->last_alpha = alpha;

    Free(xa);
    Free(xb);
    Free(xc);
    Free(xd);
    Free(init_x);

    if(tree->verbose) printf("\n. Exiting Linmin...\n");

    return alpha;
}

/*********************************************************/
/*********************************************************/
/*********************************************************/
