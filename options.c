/*

FitModeL :  a program that Fits substitution models to the data
using Maximum Likelihood.  

Copyright (C) Stephane Guindon. jul. 2004 onward

All parts of  the source except where indicated  are distributed under
the GNU public licence.  See http://www.opensource.org for details.

*/

#include "config.h"
#include "utilities.h"
#include "options.h"
#include "free.h"
#include <sys/types.h>
#include <unistd.h>



/*********************************************************/

void Usage()
{
  Exit("Usage: just type 'fitmodel'\n");
}

/*********************************************************/

option *Get_Input(int argc, char **argv)
{
  option *input;
  char choix;
  char *s, *buff;
  int n_trial;
  
  input = NULL;

  s     = (char *)mCalloc(T_MAX_LINE,sizeof(char));
  buff  = (char *)mCalloc(T_MAX_LINE,sizeof(char));


  putchar('\n');
  input                            = (option *)mCalloc(1,sizeof(option));
  input->fp_seq                    = NULL;
  input->fp_tree                   = NULL;
  input->mod                       = Make_Model_Basic();
  input->seqfile                   = (char *)mCalloc(T_MAX_FILE,sizeof(char));
  input->use_default_dnds          = (char *)mCalloc(10,sizeof(char));
  input->use_default_ssvs_switches = (char *)mCalloc(10,sizeof(char));
  input->input_tree_file           = (char *)mCalloc(T_MAX_FILE,sizeof(char));
  input->output_tree_file          = (char *)mCalloc(T_MAX_FILE,sizeof(char));
  input->output_stat_file          = (char *)mCalloc(T_MAX_FILE,sizeof(char));
  input->stat_file_open_mode       = 1;
  input->tree_file_open_mode       = 1;
  input->n_data_set_asked          = -1;
  input->user_len                  = -1;
  input->sel_class_target          = 2;


  Set_Defaults_Model(input->mod);

  Init_Optimiz(input->mod->s_opt);

  input->mod->c_code = Make_Code(0);

  Init_GTR_4_Codon_Table(input->mod);


  if(argc != 1) Command_Line(input,argv,argc);
  else
    {
#if defined(EVOLVE)
      char *n_data_sets;

      printf("Enter the tree file name > "); fflush(NULL);
      Getstring_Stdin(input->input_tree_file);
      input->fp_tree = Openfile(input->input_tree_file,0);
      printf("\n");
      
      printf("Enter the reference sequence file name > "); fflush(NULL);
      Getstring_Stdin(input->seqfile);
      input->fp_seq = Openfile(input->seqfile,0);
      printf("\n");
      
      printf("Number of data sets > ");
      n_data_sets = (char *)mCalloc(T_MAX_LINE,sizeof(char));
      Getstring_Stdin(n_data_sets);
      n_trial = 0;
      while((!atoi(n_data_sets)) || (atoi(n_data_sets) < 0))
	{
	  if(++n_trial > 10) Exit("\n. Err : the number of sets must be a positive integer\n");
	  printf("\nThe number of sets must be a positive integer\n");
	  printf("Enter a new value > ");
	  Getstring_Stdin(n_data_sets);
	}
      input->n_data_set_asked = atoi(n_data_sets);
      Free(n_data_sets);
      
#elif defined(FITMODEL) || defined(POSSELSITEID)
      
      printf("Enter the tree file name > "); fflush(NULL);
      Getstring_Stdin(input->input_tree_file);
      input->fp_tree = Openfile(input->input_tree_file,0);
      printf("\n");
      
      printf("Enter the reference sequence file name > "); fflush(NULL);
      Getstring_Stdin(input->seqfile);
      input->fp_seq = Openfile(input->seqfile,0);
      printf("\n");

#elif PHYML
      printf("Enter the sequence file name > "); fflush(NULL);
      Getstring_Stdin(input->seqfile);
      input->fp_seq = Openfile(input->seqfile,0);
#endif
            
      Open_Output_Files(input,input->seqfile);
      
      if(Filexists(input->output_stat_file))
	{
	  printf("\n");
	  printf("A file '%s' already exists\n",input->output_stat_file);
	  printf("Do you want to Replace or Append to it ?\n");
	  n_trial = 0;
	  do
	    {
	      printf("Please type R or A > ");
	      if(!scanf("%c",&choix))
		{
		  printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
		  Exit("");
		}
	      if(choix == '\n') choix = 'r'; 
	      else getchar();
	      if(++n_trial>10) Exit("\n");
	      Uppercase(&choix);
	    }
	  while((choix != 'R') && (choix != 'A'));
	  
	  switch(choix)
	    {
	    case 'R' :
	      {
		input->stat_file_open_mode = 1;
		break;
	      }
	    case 'A' :
	      {
		input->stat_file_open_mode = 2;
		break;
	      }
	    default : 
	      {
		printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
		Exit("");
		break;
	      }
	    }
	}

      
      if(Filexists(input->output_tree_file)) 
	{
	  printf("\n");
	  printf("A file '%s' already exists\n",input->output_tree_file);
	  printf("Do you want to Replace or Append to it ?\n");
	  n_trial = 0;
	  do
	    {
	      printf("Please type R or A > ");
	      if(!scanf("%c",&choix))
		{
		  printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
		  Exit("");
		}
	      if(choix == '\n') choix = 'r'; 
	      else getchar();
	      Uppercase(&choix);
	      if(++n_trial>10) Exit("\n");
	    }
	  while((choix != 'R') && (choix != 'A'));

	  
	  switch(choix)
	    {
	    case 'R' :
	      {
		input->stat_file_open_mode = 1;
		break;
	      }
	    case 'A' :
	      {
		input->stat_file_open_mode = 2;
		break;
	      }
	    default : 
	      {
		printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
		Exit("");
		break;
	      }
	    }
	}
    
      /* Set default values */
      choix                          = 0;
      input->inputtree               = 1;
      input->tree                    = NULL;
      input->interleaved             = 1;
      input->n_data_sets             = 1;
      input->start_at_data_set       = 1;
      strcpy(input->use_default_dnds,"Yes");
      strcpy(input->use_default_ssvs_switches,"Yes");
      input->pos_sel_sites_id_method = 0;

      do
	{
#ifdef WIN32
	  system("cls");
#elif UNIX
	  printf("\033[2J\033[H");
#endif
	  
	  
	  printf("\n - %s %s - \n\n\n",PROGNAME,RELEASE);
	  
	  printf("Settings for this run:\n\n");
	  
	  switch(input->mod->model_applies_to)
	    {
	      case NT : 
		{
		  choix = Print_Menu_Nt_AA(input,s);
		  break;
		}
	    case CODONS :
	      {
		choix = Print_Menu_Codon(input,s);
		break;
	      }
	    default : 
	      {
		printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
		Exit("");
	      }
	    }

	  if (choix == 'Y') 
	    {
	      break;
	    }
	  
	  if(choix != -1)
	    {
	      switch(choix)
		{
		  
#ifdef PHYML
		case 'U' :
		  {
		    if(!input->inputtree) input->inputtree = 1;
		    else                  
		      {
			input->inputtree = 0;
			printf("Enter the name of the tree file > ");
			Getstring_Stdin(input->input_tree_file);
			input->fp_tree = Openfile(input->input_tree_file,0);
		      }
		    break;
		  }
#endif
		  
		case 'N' :
		  {
		    if(input->mod->datatype == AA) 
		      Exit("\n. 'N' is not a valid choice with amino acids sequences.\n");
		    else
		      {
			switch(input->mod->model_applies_to)
			  {
			  case NT : 
			    {
			      input->mod->model_applies_to = CODONS;
			      input->mod->stepsize         = 3;
			      input->mod->ns               = input->mod->c_code->n_sense_c;
			      input->mod->ns_codon         = input->mod->c_code->n_sense_c;
			      input->mod->n_catq           = 3;
			      input->mod->n_omega          = 3;
			      input->mod->n_catg           = 1;
			      input->mod->invar            = 0;
			      input->mod->pinvar           = .0;
			      input->mod->model_number     = 24;
			      input->mod->subst_modelname  = M2;
			      break;
			    }
			  case CODONS:
			    {
			      input->mod->model_applies_to = NT;
			      input->mod->stepsize         = 1;
			      input->mod->ns               = 4;
			      input->mod->n_catq           = 1;
			      input->mod->n_omega          = 0;
			      input->mod->subst_modelname  = HKY85;
			      input->mod->model_number     = 4;
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
		  
		case 'D' :
		  {
		    if(input->mod->datatype == NT)
		      {
			input->mod->datatype         = AA;
			input->mod->model_applies_to = AA;
			input->mod->stepsize         = 1;
			input->mod->model_number     = 11;
			input->mod->subst_modelname  = DAYHOFF;
		      }
		    else
		      {
			input->mod->datatype         = NT;
			input->mod->model_applies_to = NT;
			input->mod->stepsize         = 1;
			input->mod->model_number     = 4;
			input->mod->subst_modelname  = HKY85;
		      }
		    break;
		  }
		  
		  
		  
		case 'T' :
		  {
		    char answer;
		    
		    if((input->mod->datatype == AA)   || 
		       (input->mod->subst_modelname == JC69) ||
		       (input->mod->subst_modelname == F81) ||
		       (input->mod->subst_modelname == GTR)) 
		      Exit("\n 'T' is not a valid choice for this model\n");
		    
		    switch(input->mod->s_opt->opt_kappa)
		      {
		      case 0 : 
			{
			  printf("Optimize ts/tv ratio ? [Y/n] ");
			  if(!scanf("%c", &answer))
			    {
			      printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
			      Exit("");
			    }
			  if(answer == '\n') answer = 'Y';
			  else getchar();
			  break;
			}
		      case 1 : 
			{
			  printf("Optimize ts/tv ratio ? [N/y] ");
			  if(!scanf("%c", &answer))
			    {
			      printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
			      Exit("");
			    }
			  if(answer == '\n') answer = 'N';
			  else getchar();
			  break;
			}
		      default : Exit("\n");
		      }
		    
		    n_trial = 0;
		    while((answer != 'Y') && (answer != 'y') &&
			  (answer != 'N') && (answer != 'n'))  
		      {
			if(++n_trial > 10) Exit("\nErr : wrong answers !");
			printf("Optimize ts/tv ratio ? [N/y] ");
			if(!scanf("%c", &answer))
			  {
			    printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
			    Exit("");
			  }
			if(answer == '\n') answer = 'N';
			else getchar();
		      }
		    
		    switch(answer)
		      {
		      case 'Y' : case 'y' : 
			{
			  input->mod->kappa = 4.0;
			  input->mod->s_opt->opt_kappa = 1; 
			  if(input->mod->model_number == 6) 
			    input->mod->s_opt->opt_lambda = 1;
			  break;
			}
		      case 'N' : case 'n' : 
			{
			  char *t;
			  t = (char *)mCalloc(T_MAX_LINE,sizeof(char));
			  input->mod->s_opt->opt_kappa = 0; 
			  input->mod->kappa = 4.0;
			  
			  Read_Param("kappa",4.0,.0001,100.,t,&(input->mod->kappa));

			  input->mod->s_opt->opt_kappa  = 0;
			  input->mod->s_opt->opt_lambda = 0;
			  Free(t);		
			  break;
			}
		      }	 
		    break;
		  }	  
		  
		case 'I' : 
		  {
		    if(input->interleaved)
		      input->interleaved = 0;
		    else input->interleaved = 1;
		    break;
		  }
		  
		case 'O' : 
		  {
		    input->mod->s_opt->opt_param = (input->mod->s_opt->opt_param)?(0):(1);
		    break;
		  }
		  
#if defined(EVOLVE) || defined(POSSELSITEID) 
		case 'L' :
		  {
		    char *len;
                    if((input->pos_sel_sites_id_method != 0) &&
                       (input->pos_sel_sites_id_method != 3)) Exit("\n. Wrong option !\n");
		    len = (char *)mCalloc(T_MAX_LINE,sizeof(char));
		    printf("Enter the sequence length > ");
		    Getstring_Stdin(len);
		    n_trial = 0;
		    while((!atof(len)) || (atof(len) < 0.0-MDBL_MIN))
		      {
			if(++n_trial > 10)
			  Exit("\nErr : sequence length must be a positive integer \n");
			printf("Sequence length must be a positive integer \n");
			printf("Enter a new value > ");
			Getstring_Stdin(len);
		      }
		    input->user_len = atoi(len);	    
		    Free(len);
		    break;
		  }
#endif
		  
		default : 
		  {
		    break;
		  }
		}
	    }
	}while(1);
      
      if((input->mod->model_number == 1) || (input->mod->model_number == 3))
	{
	  input->mod->s_opt->opt_kappa  = 0;
	  input->mod->s_opt->opt_lambda = 0;
	}
      if(input->mod->subst_modelname != TN93) input->mod->s_opt->opt_lambda = 0;           
    }

  if(input->mod->subst_modelname == MX)
    {
      char *s;
      int i;
      fit_double w;

      s = (char *)mCalloc(T_MAX_LINE,sizeof(char));

      if(!strcmp(input->use_default_dnds,"Yes"))
	{
	  printf("Number of rate matrices > ");
	  Getstring_Stdin(s);
	  n_trial = 0;
	  while((!atof(s)) || (atof(s) < 0.0-MDBL_MIN) || (atof(s) > N_MAX_CATQ))
	    {
	      if(++n_trial > 10)
		{
		  printf("\nErr : number of rate matrices must be a positive integer smaller than %d.\n",N_MAX_CATQ);
		  Exit("\n");
		}
	      printf("The number of rate matrices must be a positive integer smaller than %d.\n",N_MAX_CATQ);
	      printf("Enter a new value > ");
	      Getstring_Stdin(s);
	    }
	  
	  input->mod->n_omega              = (fit_double)atoi(s);
	}
      
      if(input->mod->switch_modelname == NO_SWITCH)
	{
	  input->mod->n_catq           = input->mod->n_omega;	    
	  input->mod->ns               = input->mod->c_code->n_sense_c;
	}
      else
	{
	  input->mod->n_catq           = 1;	    
	  input->mod->ns               = input->mod->n_omega*input->mod->c_code->n_sense_c;
	}
	  
      input->mod->s_opt->opt_omega   = input->mod->n_catq;	    

#ifndef EVOLVE
      if(input->mod->s_opt->opt_param == 1)
#endif
	{
	  For(i,input->mod->n_omega) 
	    {
	      printf("Minimum value that w%d can take > ",i);
	      Getstring_Stdin(s);
	      w = atof(s);
	      if(w < 0.0)
		{
		  printf("\nErr: min must be greater than 0.\n");
		  Exit("\n");
		}
	      input->mod->omega_min[i] = w;
              if(input->mod->omega_min[i] < 0.001) 
                {
                  printf("\n. WARNING: min value of omega set to 0.001 so as to avoid numerical precision issues");
                  input->mod->omega_min[i] = 0.001;
                }
	    }
	  
	  For(i,input->mod->n_omega) 
	    {
	      printf("Maximum value that w%d can take > ",i);
	      Getstring_Stdin(s);
	      w = atof(s);
	      if(w < input->mod->omega_min[i])
		{
		  printf("\nErr: max must be greater than %.2f\n",input->mod->omega_min[i]);
		  Exit("\n");
		}
	      input->mod->omega_max[i] = w;

              if(input->mod->omega_max[i] > 20.) 
                {
                  printf("\n. WARNING: max value of omega set to 20. so as to avoid numerical precision issues");
                  input->mod->omega_max[i] = 20.;
                }

	    }
	  For(i,input->mod->n_omega) 
	    {
	      if(input->mod->omega[i] < input->mod->omega_min[i] ||
		 input->mod->omega[i] > input->mod->omega_max[i])
		{
		  input->mod->omega[i] = input->mod->omega_min[i] + (input->mod->omega_max[i]-input->mod->omega_min[i])/2.;
		}
	    }
          
          For(i,input->mod->n_omega) 
            input->mod->omega_proba[i] = (double)1./input->mod->n_omega;


	}
      Free(s);
    }
  
  Free(buff);
  Free(s);

  Check_Inconsistencies(input);

  return input;
}

/*********************************************************/
 
char Print_Menu_Nt_AA(option *input, char *s)
{
  char choix;
  int n_trial;

  choix = -1;



  printf("  D "
	 "                                Data type (DNA/AA) "
	 " %-15s \n",
	 (input->mod->datatype == AA)?("AA"):("DNA"));

  if((input->mod->model_number < 10) && (input->mod->model_number > 3))
    printf("  F "
	   "           Base frequency estimates (empirical/ML) "
	   " %-15s \n",
	   (input->mod->s_opt->opt_bfreq)?("ML"):("empirical"));

  
#ifdef PHYML
  printf("  U "
	 "                      Input tree (BIONJ/user tree) "
	 " %-15s \n",
	 (input->inputtree)?("BIONJ"):("user tree"));
#endif
  
  if (input->mod->datatype == NT)
    {
      char *s;
 
      s = Return_Model_Applies_To(input->mod);
      printf("  N  "
	     "                DNA substitution model applies to  "
	     "%-15s \n",s);
      Free(s);

      s = Return_Subst_Model_Name(input->mod);
      printf("  M  "
	     "                Model of nucleotides substitution "
	     " %-15s \n",s);
      Free(s);
    }
  else
    {
      char *s;
        
      s = Return_Subst_Model_Name(input->mod);
      printf("  M  "
	     "                Model of amino-acids substitution "
	     " %-15s \n",s);
      Free(s);
    }
  
  if(input->mod->s_opt->opt_pinvar)
    {
      strcpy(s,"estimated");      
    }
  else
    {
      strcpy(s,"fixed");
      strcat(s," (p-invar = ");
      sprintf(s+(int)strlen(s),"%3.2f)",input->mod->pinvar);
    }
  
  printf("  V  "
	 " Proportion of invariable sites (fixed/estimated)"
	 "  %-15s \n",s);
  
  if ((input->mod->datatype == NT) && 
      (input->mod->subst_modelname != JC69) && 
      (input->mod->subst_modelname != F81) && 
      (input->mod->subst_modelname != GTR))
    {	        
      if(input->mod->s_opt->opt_kappa)
	{
	  strcpy(s,"estimated");
	}
      else
	{
	  strcpy(s,"fixed");
	  strcat(s," (ts/tv = ");
	  sprintf(s+(int)strlen(s),"%3.2f)",input->mod->kappa);
	}

      printf("  T "
	     "                     Ts/tv ratio (fixed/estimated) "
	     " %-15s \n",s);
    }
  
  printf("  R "
	 "       One category of substitution rates (yes/no) "
	 " %-15s \n",
	 (input->mod->n_catg > 1)?("no"):("yes"));
  
  
  if(input->mod->n_catg > 1)
    {
      printf("  C "
	     "           Number of substitution rates categories "
	     " %-15d \n",
	     input->mod->n_catg);
    }
  
  if(input->mod->n_catg > 1)
    {
      if(input->mod->s_opt->opt_alpha)
	{
	  strcpy(s,"estimated");
	}
      else
	{
	  strcpy(s,"fixed");
	  strcat(s,"(alpha = ");
	  sprintf(s+(int)strlen(s),"%3.2f)",input->mod->alpha);
	}

      printf("  A "
	     "                           Alpha (fixed/estimated) "
	     " %-15s \n",s);
    }
  
  
  strcpy(s,"");
  
  /*       sprintf(s," (%d sets)",input->n_data_sets); */
  /*       strcpy(buff,(input->n_data_sets > 1)?("yes"):("no")); */
  /*       buff=strcat(buff,(input->n_data_sets > 1)?(s):("\0")); */
  
  strcpy(s,(input->n_data_sets > 1)?("yes"):("no"));
  if(input->n_data_sets > 1)
    sprintf(s+(int)strlen(s)+1," (%d sets)",input->n_data_sets);
  else
    strcat(s,"\0");
  
  
  
  
  printf("  S "
	 "                        Analyse multiple data sets "
	 " %-15s \n",s);
  
  
  printf("  I "
	 "       Input sequences interleaved (or sequential) "
	 " %-15s \n",
	 (input->interleaved)?("interleaved"):("sequential"));
  

  printf("  O "
	 "                  Optimise substitution parameters "
	 " %-15s \n",
	 (input->mod->s_opt->opt_param)?("yes"):("no"));



#if defined(EVOLVE) || defined(POSSELSITEID)
  strcpy(s,"");
  if(input->user_len==-1) strcpy(s,"Reference data set length");
  else sprintf(s,"l = %d",input->user_len);
 
  printf("  L "
	 "                                   Sequence length "
	 " %-15s \n",s);
#endif
  
  
  printf("\n");


  printf("\nAre these settings correct? "
	 "(type  Y  or letter for one to change)  ");
  
  if(!scanf("%c",&choix))
    {
      printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Exit("");
    }
  if(choix == '\n') choix = 'X'; 
  else getchar(); /* \n */
  
  Uppercase(&choix);
    
  switch(choix)
    {

    case 'F' :
      {
	if((input->mod->model_number > 10) ||
	   (input->mod->model_number < 4)) Exit("\n. Invalid choice...\n");
	input->mod->s_opt->opt_bfreq = (input->mod->s_opt->opt_bfreq)?(0):(1);
	break;
      }


    case 'M' :
      {
	if(input->mod->datatype == NT)
	  {
	    if(input->mod->subst_modelname == JC69)
	      {
		input->mod->model_number    = 2;
		input->mod->subst_modelname = K80;
	      }
	    else if(input->mod->subst_modelname == K80)
	      {
		input->mod->model_number     = 3;
		input->mod->subst_modelname  = F81;
		input->mod->s_opt->opt_kappa = 0;
	      }
	    else if(input->mod->subst_modelname == F81)
	      {
		input->mod->model_number    = 4;
		input->mod->subst_modelname = HKY85;
	      }
	    else if(input->mod->subst_modelname == HKY85)
	      {
		input->mod->model_number    = 5;
		input->mod->subst_modelname = F84;
	      }
	    else if(input->mod->subst_modelname == F84)
	      {
		input->mod->model_number    = 6;
		input->mod->subst_modelname = TN93;
		if(input->mod->s_opt->opt_kappa) input->mod->s_opt->opt_lambda = 1;
	      }
	    else if(input->mod->subst_modelname == TN93)
	      {
		input->mod->model_number     = 7;
		input->mod->subst_modelname  = GTR;
		input->mod->s_opt->opt_kappa = 0;
	      }
	    else if(input->mod->subst_modelname == GTR)
	      {
		input->mod->model_number     = 1;
		input->mod->subst_modelname  = JC69;
		input->mod->s_opt->opt_kappa = 0;
	      }
	  }
	else
	  {
	    if(input->mod->subst_modelname == DAYHOFF)
	      {
		input->mod->model_number    = 12;
		input->mod->subst_modelname = JTT;
	      }
	    else if(input->mod->subst_modelname == JTT)
	      {
		input->mod->model_number    = 13;
		input->mod->subst_modelname = MTREV;
	      }
	    else if(input->mod->subst_modelname == MTREV)
	      {
		input->mod->model_number    = 14;
		input->mod->subst_modelname = WAG;
	      }

	    else if(input->mod->subst_modelname == WAG)
	      {
		input->mod->model_number    = 15;
		input->mod->subst_modelname = DCMUT;
	      }

	    else if(input->mod->subst_modelname == DCMUT)
	      {
		input->mod->model_number    = 11;
		input->mod->subst_modelname = DAYHOFF;
	      }
	  }
	break;
      }

    case 'R' :
      {
	(input->mod->n_catg == 1)?(input->mod->n_catg = 4):(input->mod->n_catg = 1);
	break;
      }

    case 'C' :
      {
	char *c;
	printf("Enter your number of categories > ");
	c = (char *)mCalloc(T_MAX_LINE,sizeof(char));
	Getstring_Stdin(c);
	n_trial = 0;
	while((!atoi(c)) || (atoi(c) < 0))
	  {
	    if(++n_trial > 10) Exit("\nErr : the number of categories must be a positive integer\n");
	    printf("\nThe number of categories must be a positive integer\n");
	    printf("Enter a new value > ");
	    Getstring_Stdin(c);
	  }
	input->mod->n_catg = atoi(c);
	Free(c);
	break;
      }
    case 'A' :
      {
	char answer;
	
	switch(input->mod->s_opt->opt_alpha)
	  {
	  case 0 : 
	    {
	      printf("Optimize alpha ? [Y/n] ");
	      if(!scanf("%c",&answer))
		{
		  printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
		  Exit("");
		}
	      if(answer == '\n') answer = 'Y';
	      else getchar();
	      break;
	    }
	  case 1 : 
	    {
	      printf("Optimize alpha ? [N/y] ");
	      if(!scanf("%c",&answer))
		{
		  printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
		  Exit("");
		}
	      if(answer == '\n') answer = 'N';
	      else getchar();
	      break;
	    }
	  default : Exit("\n");
	  }
	
	n_trial = 0;
	while((answer != 'Y') && (answer != 'y') &&
	      (answer != 'N') && (answer != 'n'))  
	  {
	    if(++n_trial > 10) Exit("\nErr : wrong answers !");
	    printf("Optimize alpha ? [N/y] ");
	    if(!scanf("%c",&answer))
	      {
		printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
		Exit("");
	      }
	    if(answer == '\n') answer = 'N';
	    else getchar();
	  }
	
	switch(answer)
	  {
	  case 'Y' : case 'y' : 
	    {
	      input->mod->s_opt->opt_alpha = 1; 
	      break;
	    }
	  case 'N' : case 'n' : 
	    {
	      char *a;
	      a = (char *)mCalloc(T_MAX_LINE,sizeof(char));
	      input->mod->alpha = 1.0;
	      input->mod->s_opt->opt_alpha = 0; 
	      
	      Read_Param("alpha",1.0,.01,100.,a,&(input->mod->alpha));
	      
	      Free(a);
	      input->mod->s_opt->opt_alpha  = 0;
	      break;
	    }
	  }	 
	break;
      }
    
    case 'V' : 
      {
	char answer;
	
	switch(input->mod->s_opt->opt_pinvar)
	  {
	  case 0 : 
	    {
	      printf("Optimize p-invar ? [Y/n] ");
	      if(!scanf("%c", &answer))
		{
		  printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
		  Exit("");
		}
	      if(answer == '\n') answer = 'Y';
	      else getchar();
	      break;
	    }
	  case 1 : 
	    {
	      printf("Optimize p-invar ? [N/y] ");
	      if(!scanf("%c", &answer))
		{
		  printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
		  Exit("");
		}
	      if(answer == '\n') answer = 'N';
	      else getchar();
	      break;
	    }
	  default : Exit("\n");
	  }
	
	n_trial = 0;
	while((answer != 'Y') && (answer != 'y') &&
	      (answer != 'N') && (answer != 'n'))  
	  {
	    if(++n_trial > 10) Exit("\nErr : wrong answers !");
	    printf("Optimize p-invar ? [N/y] ");
	    if(!scanf("%c", &answer))
	      {
		printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
		Exit("");
	      }
	    if(answer == '\n') answer = 'N';
	    else getchar();
	  }
	
	switch(answer)
	  {
	  case 'Y' : case 'y' : 
	    {
	      input->mod->s_opt->opt_pinvar = 1; 
	      input->mod->pinvar            = 0.5;
	      input->mod->invar             = 1;
	      break;
	    }
	  case 'N' : case 'n' : 
	    {
	      char *p;
	      p = (char *)mCalloc(T_MAX_LINE,sizeof(char));
	      
	      Read_Param("p-invar",0.3,.0,1.,p,&(input->mod->pinvar));
	      
	      if(input->mod->pinvar > 0.0+MDBL_MIN) input->mod->invar = 1;
	      else                                  input->mod->invar = 0;
	      
	      Free(p);
	      
	      input->mod->s_opt->opt_pinvar = 0;
	      break;
	    }
	  }	 
	break;
      }
    case 'S' :
      {
	char *c;
	int n_trial;
	
	printf("How many data sets > ");
	c = (char *)mCalloc(T_MAX_LINE,sizeof(char));
	Getstring_Stdin(c);
	
	n_trial = 0;
	while((!atoi(c)) || (atoi(c) < 0))
	  {
	    if(++n_trial > 10) Exit("\nErr : The number of data sets must be a positive integer\n");
	    printf("\nThe number of data sets must be a positive integer\n");
	    printf("Enter a new value > ");
	    Getstring_Stdin(c);
	  }
	input->n_data_sets = atoi(c);
	
        printf("Start at data set [default=1] > ");
	Getstring_Stdin(c);
	
	n_trial = 0;
	while(atoi(c) < 0)
	  {
	    if(++n_trial > 10) Exit("\nErr : The number of data sets must be a positive integer\n");
	    printf("\nThe number of data sets must be a positive integer\n");
	    printf("Enter a new value > ");
	    Getstring_Stdin(c);
	  }
	input->start_at_data_set = atoi(c);

	input->n_data_sets += input->start_at_data_set-1;


	Free(c);
	break;
      }
    }
  return choix;
}

/*********************************************************/

char Print_Menu_Codon(option *input, char *s)
{
  
  char choix;

  choix = -1;


  printf("  D "
	 "                                Data type (DNA/AA) "
	 " %-15s \n",
	 (input->mod->datatype == AA)?("AA"):("DNA"));

#ifdef PHYML
  printf("  U "
	 "                      Input tree (BIONJ/user tree) "
	 " %-15s \n",
	 (input->inputtree)?("BIONJ"):("user tree"));
#endif
  
  if (input->mod->datatype == NT)
    {
      char *s;

      s = Return_Model_Applies_To(input->mod);
      printf("  N  "
	     "                DNA substitution model applies to  "
	     "%-15s \n", s);
      Free(s);

      printf("  G "
	     "                                      Genetic code  "
	     "%-15s \n",input->mod->c_code->name);
        
      s = Return_Subst_Model_Name(input->mod);
      printf("  C  "
	     "                      Model of codon substitution "
	     " %-15s \n",s);
      Free(s);
        
        
      s = Return_Switch_Model_Name(input->mod);
      printf("  M  "
	     "      Model of switches between selection regimes "
	     " %-15s \n",s);
      Free(s);


    }
    

  printf("  W "
	 "         Start with default dn/ds ratio parameters  "
	 "%-15s \n",input->use_default_dnds);
  
  
  if(input->mod->switch_modelname != NO_SWITCH)
    {
      printf("  R "
	     "     Start with default values for switching rates  "
	     "%-15s \n",input->use_default_ssvs_switches);
      
/*       s = Return_Target_Selection_Class(input); */
/*       printf("  B " */
/* 	     "                       Selection class to focus on  " */
/* 	     "%-15s \n",s); */
/*       Free(s); */
    }
  
  if(input->mod->s_opt->opt_kappa)
    {
      strcpy(s,"estimated");
    }
  else
    {
      strcpy(s,"fixed");
      strcat(s," (ts/tv = ");
      sprintf(s+(int)strlen(s),"%3.2f)",input->mod->kappa);
    }
  
  printf("  T "
	 "                     Ts/tv ratio (fixed/estimated) "
	 " %-15s \n",s);
      
  
  printf("  F "
	 "                Codon frequency (uniform, or F3X4) "
	 " %-15s \n",(input->mod->codon_freq)?("F3X4"):("uniform"));

    

  strcpy(s,(input->n_data_sets > 1)?("yes"):("no"));
  if(input->n_data_sets > 1)
    sprintf(s+(int)strlen(s)+1," (%d sets)",input->n_data_sets);
  else
    strcat(s,"\0");
  

  printf("  S "
	 "                        Analyze multiple data sets "
	 " %-15s \n",s);
  
  
  printf("  I "
	 "       Input sequences interleaved (or sequential) "
	 " %-15s \n",
	 (input->interleaved)?("interleaved"):("sequential"));
  
  printf("  O "
	 "                  Optimise substitution parameters "
	 " %-15s \n",
	 (input->mod->s_opt->opt_param)?("yes"):("no"));


#if defined(EVOLVE) || defined(POSSELSITEID)
  if(
     (input->pos_sel_sites_id_method == ADAPTIVE_CONTROL) ||
     (input->pos_sel_sites_id_method == PARAMETRIC_BOOTSTRAP)
     )
      
    {
      strcpy(s,"");
      if(input->user_len==-1)
	strcpy(s,"Reference data set length");
      else sprintf(s,"l = %d",input->user_len);
          
      printf("  L "
	     "                                   Sequence length "
	     " %-15s \n",s);
    }

  switch(input->pos_sel_sites_id_method)
    {
    case 0 :  {strcpy(s,"Adaptive control");     break;}
    case 1 :  {strcpy(s,"Newton et al. (2004)"); break;}
    case 2 :  {strcpy(s,"Fixed region");         break;}
    case 3 :  {strcpy(s,"Parametric bootstrap"); break;}
    default : {Exit("\n");                       break;}
    }

  printf("  P "
	 "            Detection of positively selected sites "
         " %-15s \n",s);

  if(input->pos_sel_sites_id_method == FIXED_REGION)
    {
      printf("  V "
	     "                                   Threshold value "
	     " %.3f \n",input->mod->thresholdp2);
    }

  if(input->pos_sel_sites_id_method != FIXED_REGION)
    {
      printf("  X "
	     "                        False detection rate value "
	     " %.3f \n",input->mod->fdr);
    }
    
#endif
  
  
  printf("\n");
  
  printf("\nAre these settings correct? "
	 "(type  Y  or letter for one to change)  ");
  
  if(!scanf("%c",&choix))
    {
      printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Exit("");
    }
  if(choix == '\n') choix = 'X'; 
  else getchar(); /* \n */
  
  Uppercase(&choix);
  
  
  switch(choix)
    {
#if defined(EVOLVE) || defined(POSSELSITEID)
    case 'X' :
      {
	char *t;

	t = (char *)mCalloc(T_MAX_LINE,sizeof(char));
	Read_Param("fdr",input->mod->fdr,.0,1.,t,&(input->mod->fdr));
	Free(t);
	break;
      }
#endif

    case 'B' :
      {
	input->sel_class_target++;
	if(input->sel_class_target > 2) 
	  input->sel_class_target = 0;	
	break;
      }
    case 'V' :
      {
	char *t;

	t = (char *)mCalloc(T_MAX_LINE,sizeof(char));
	Read_Param("threshold",input->mod->thresholdp2,.0,1.,t,&(input->mod->thresholdp2));
	Free(t);
	break;
      }
    case 'P' :
      {
	switch(input->pos_sel_sites_id_method)
	  {
	  case 0 : {input->pos_sel_sites_id_method = 1; break;}
	  case 1 : {input->pos_sel_sites_id_method = 2; break;}
	  case 2 : {input->pos_sel_sites_id_method = 3; break;}
	  case 3 : {input->pos_sel_sites_id_method = 0; break;}
	  default : {Exit("\n"); break;}
	  }
	break;
      }
    case 'G' :
      {
	int code_num;

	code_num = input->mod->c_code->num_curr_code;
	Free_Code(input->mod->c_code);

	if(code_num == 17)
	  code_num = 0;
	else code_num++;
	input->mod->c_code   = Make_Code(code_num);
	input->mod->ns       = input->mod->c_code->n_sense_c;
	input->mod->ns_codon = input->mod->c_code->n_sense_c;

	break;
      }

    case 'W' :
      {
	char *t;
	
	t = (char *)mCalloc(T_MAX_LINE,sizeof(char));
	
	if(!strcmp(input->use_default_dnds,"Yes"))
	  {
	    if(input->mod->subst_modelname == M3)
	      {
		strcpy(input->use_default_dnds,"No");
		Read_Param("w1",input->mod->omega[0],.0,10000.,t,&(input->mod->omega[0]));
		strcpy(t,"");
		Read_Param("w2",input->mod->omega[1],.0,10000.,t,&(input->mod->omega[1]));
		strcpy(t,"");
		Read_Param("w3",input->mod->omega[2],.0,10000.,t,&(input->mod->omega[2]));
		strcpy(t,"");		
		Read_Param("p1",input->mod->omega_proba[0],.0,1.,t,&(input->mod->omega_proba[0]));
		strcpy(t,"");		
		Read_Param("p2",input->mod->omega_proba[1],.0,1.,t,&(input->mod->omega_proba[1]));
		strcpy(t,"");		
		Read_Param("p3",input->mod->omega_proba[2],.0,1.,t,&(input->mod->omega_proba[2]));
		strcpy(t,"");
	      }

	    else if(input->mod->subst_modelname == M2a)
	      {
		strcpy(input->use_default_dnds,"No");
		Read_Param("w1",input->mod->omega[0],.0,1.,t,&(input->mod->omega[0]));
		strcpy(t,"");
		Read_Param("w3",input->mod->omega[2],.0,10000.,t,&(input->mod->omega[2]));
		strcpy(t,"");		
		Read_Param("p1",input->mod->omega_proba[0],.0,1.,t,&(input->mod->omega_proba[0]));
		strcpy(t,"");		
		Read_Param("p2",input->mod->omega_proba[1],.0,1.,t,&(input->mod->omega_proba[1]));
		strcpy(t,"");		
		Read_Param("p3",input->mod->omega_proba[2],.0,1.,t,&(input->mod->omega_proba[2]));
		strcpy(t,"");
	      }

	    else if(input->mod->subst_modelname == MX)
	      {
		int i,n_trial;
		char *s;
		s = (char *)mCalloc(T_MAX_LINE,sizeof(char));
		
		strcpy(input->use_default_dnds,"No");
      
		printf("Number of rate matrices > ");
		Getstring_Stdin(s);
		n_trial = 0;
		while((!atof(s)) || (atof(s) < 0.0-MDBL_MIN) || (atof(s) > N_MAX_CATQ))
		  {
		    if(++n_trial > 10)
		      {
			printf("\nErr : number of rate matrices must be a positive integer smaller than %d.\n",N_MAX_CATQ);
			Exit("\n");
		      }
		    printf("The number of rate matrices must be a positive integer smaller than %d.\n",N_MAX_CATQ);
		    printf("Enter a new value > ");
		    Getstring_Stdin(s);
		  }		
		input->mod->n_omega  = (fit_double)atoi(s);

		For(i,input->mod->n_omega)
		  {
		    strcpy(s,"w");
		    sprintf(s+1,"%d",i+1);
		    Read_Param(s,input->mod->omega[i],.0,100.,t,&(input->mod->omega[i]));
		  }

		For(i,input->mod->n_omega)
		  {
		    strcpy(s,"p");
		    sprintf(s+1,"%d",i+1);
		    Read_Param(s,input->mod->omega_proba[i],.0,1.,t,&(input->mod->omega_proba[i]));
		  }
		Free(s);
	      }

	    else
	      {
		input->mod->omega[0] = 0.0;
		input->mod->omega[1] = 1.0;
	      }

	    
	    input->mod->omega_proba[0] /= 
	      input->mod->omega_proba[0]+
	      input->mod->omega_proba[1]+
	      input->mod->omega_proba[2];

	    input->mod->omega_proba[1] /= 
	      input->mod->omega_proba[0]+
	      input->mod->omega_proba[1]+
	      input->mod->omega_proba[2];

	    input->mod->omega_proba[2] /= 
	      input->mod->omega_proba[0]+
	      input->mod->omega_proba[1]+
	      input->mod->omega_proba[2];

	    Free(t);		
	  }
	else
	  {
	    strcpy(input->use_default_dnds,"Yes");
	    
	    input->mod->omega[0] = 0.0;
	    input->mod->omega[1] = 1.0;
	    input->mod->omega[2] = 4.0;
	    
	    input->mod->omega_proba[0] = 0.6;
	    input->mod->omega_proba[1] = 0.3;
	    input->mod->omega_proba[2] = 0.1;
	  }
	break;
      }
      
    case 'M' : 
      {
	if(input->mod->switch_modelname == NO_SWITCH)
	  {
	    input->mod->switch_modelname   = SWITCH_S1;
            input->mod->subst_modelname    = M2;
	    input->mod->model_number       = 24;
	    input->mod->n_catq             = 1;
	    input->mod->s_opt->opt_theta   = 1;
            input->mod->s_opt->opt_p_omega = 1;
            input->mod->n_omega            = 3;
            input->mod->s_opt->opt_omega   = 1;
            input->mod->omega[0]           = 0.0;
            input->mod->omega[1]           = 1.0;
            input->mod->omega[2]           = 4.0;
            input->mod->omega_proba[0]     = 0.6;
            input->mod->omega_proba[1]     = 0.3;
            input->mod->omega_proba[2]     = 0.1;
	    input->mod->ns                 = input->mod->n_omega*input->mod->c_code->n_sense_c;
	  }
	else if(input->mod->switch_modelname == SWITCH_S1)
	  {
	    input->mod->switch_modelname   = SWITCH_S2;
	    input->mod->s_opt->opt_theta   = 3;
	  }
	else
	  {
	    input->mod->switch_modelname   = NO_SWITCH;
            input->mod->subst_modelname    = M2;
	    input->mod->model_number       = 24;
	    input->mod->ns                 = input->mod->c_code->n_sense_c;
	    input->mod->n_catq             = 3;
	    input->mod->s_opt->opt_theta   = 1;
            input->mod->n_catq             = 3;
            input->mod->s_opt->opt_p_omega = 1;
            input->mod->n_omega            = 3;
            input->mod->s_opt->opt_omega   = 1;
            input->mod->omega[0]           = 0.0;
            input->mod->omega[1]           = 1.0;
            input->mod->omega[2]           = 4.0;
            input->mod->omega_proba[0]     = 0.6;
            input->mod->omega_proba[1]     = 0.3;
            input->mod->omega_proba[2]     = 0.1;
	  }
	break;
      }

    case 'R' :
      {
    	char *t;
	
    	t = (char *)mCalloc(T_MAX_LINE,sizeof(char));
	
    	if(input->mod->switch_modelname == NO_SWITCH) Exit("\n. Invalid choice...\n");

    	if(!strcmp(input->use_default_ssvs_switches,"Yes"))
    	  {
    	    strcpy(input->use_default_ssvs_switches,"No");
	    
    	    Read_Param("delta",1.0,0.001,100.,t,&(input->mod->io_delta_switch));
    	    strcpy(t,"");

            /* if(input->mod->switch_modelname == SWITCH_S2) */
    	    /*   { */
    	    /* 	Read_Param("alpha",1.0,0.001,100.,t,&(input->mod->qmat_struct[0]->theta[1])); */
    	    /* 	input->mod->qmat_struct[0]->theta[1] /= input->mod->qmat_struct[0]->theta[0];  */
    	    /* 	strcpy(t,""); */
                    
    	    /* 	Read_Param("beta", 1.0,0.001,100.,t,&(input->mod->qmat_struct[0]->theta[2])); */
    	    /* 	input->mod->qmat_struct[0]->theta[2] /= input->mod->qmat_struct[0]->theta[0];  */
    	    /* 	strcpy(t,""); */
    	    /*   } */
            /* else if(input->mod->switch_modelname == SWITCH_S1) */
    	    /*   { */
    	    /* 	/\* input->mod->alpha_sw = 1.0; *\/ */
    	    /* 	/\* input->mod->beta_sw  = 1.0; *\/ */
    	    /* 	input->mod->qmat_struct[0]->theta[1] = 1.0 / input->mod->qmat_struct[0]->theta[0]; */
    	    /* 	input->mod->qmat_struct[0]->theta[2] = 1.0 / input->mod->qmat_struct[0]->theta[0]; */
    	    /*   } */
    	  }
    	else
    	  {
    	    strcpy(input->use_default_ssvs_switches,"Yes");
	    	    
    	    /* input->mod->qmat_struct[0]->theta[0] = 1.0; */
    	    /* input->mod->qmat_struct[0]->theta[1] = 1.0; */
    	    /* input->mod->qmat_struct[0]->theta[2] = 1.0; */

    	    /* input->mod->delta_sw = 1.0; */
    	    /* input->mod->alpha_sw = 1.0; */
    	    /* input->mod->beta_sw  = 1.0; */
    	  }
    	Free(t);
    	break;
      }

    case 'F' :
      {
	(input->mod->codon_freq)?
	  (input->mod->codon_freq = 0):
	  (input->mod->codon_freq = 1);
	break;
      }
    case 'C' :
      {
	if(input->mod->switch_modelname == NO_SWITCH)
	  {
	    if(input->mod->subst_modelname == M2)
	      {
		input->mod->subst_modelname    = M2a;
		input->mod->model_number       = 25;
		input->mod->s_opt->opt_p_omega = 1;
		input->mod->n_catq             = 3;
		input->mod->n_omega            = 3;
		input->mod->s_opt->opt_omega   = 2;
		input->mod->omega[0]           = 0.3;
		input->mod->omega[1]           = 1.0;
		input->mod->omega[2]           = 4.0;
		input->mod->omega_proba[0]     = 0.6;
		input->mod->omega_proba[1]     = 0.3;
		input->mod->omega_proba[2]     = 0.1;
		input->mod->ns                 = input->mod->c_code->n_sense_c;
	      }
	    else if(input->mod->subst_modelname == M2a)
	      { 
		input->mod->subst_modelname    = M3;
		input->mod->model_number       = 26;
		input->mod->s_opt->opt_p_omega = 1;
		input->mod->n_catq             = 3;
		input->mod->n_omega            = 3;
		input->mod->s_opt->opt_omega   = 3;
		input->mod->omega[0]           = 0.3;
		input->mod->omega[1]           = 1.0;
		input->mod->omega[2]           = 4.0;
		input->mod->omega_proba[0]     = 0.6;
		input->mod->omega_proba[1]     = 0.3;
		input->mod->omega_proba[2]     = 0.1;
		input->mod->ns                 = input->mod->c_code->n_sense_c;
	      }
	    else if(input->mod->subst_modelname == M3)
	      {	
		input->mod->subst_modelname    = MX;
		input->mod->model_number       = 27;
		input->mod->s_opt->opt_p_omega = 1;
		input->mod->n_catq             = 3;
		input->mod->n_omega            = 3;
		input->mod->s_opt->opt_omega   = 3;
		input->mod->omega[0]           = 0.3;
		input->mod->omega[1]           = 1.0;
		input->mod->omega[2]           = 4.0;
		input->mod->omega_proba[0]     = 0.6;
		input->mod->omega_proba[1]     = 0.3;
		input->mod->omega_proba[2]     = 0.1;
		input->mod->ns                 = input->mod->c_code->n_sense_c;
	      }
	    else if(input->mod->subst_modelname == MX)
	      {
		input->mod->subst_modelname    = M0;
		input->mod->model_number       = 21;
		input->mod->n_catq             = 1;
		input->mod->n_omega            = 1;
		input->mod->s_opt->opt_omega   = 1;
		input->mod->s_opt->opt_p_omega = 0;
		input->mod->omega[0]           = 1.0;
		input->mod->omega_proba[0]     = 1.0;
		input->mod->ns                 = input->mod->c_code->n_sense_c;
	      }

	    else if(input->mod->subst_modelname == M0)
	      {
		input->mod->subst_modelname    = M1;
		input->mod->model_number       = 22;
		input->mod->n_catq             = 2;
		input->mod->n_omega            = 2;
		input->mod->s_opt->opt_omega   = 0;
		input->mod->s_opt->opt_p_omega = 1;
		input->mod->omega[0]           = 0.0;
		input->mod->omega[1]           = 1.0;
                input->mod->omega_min[0]       = 0.0;
                input->mod->omega_max[0]       = 20.0;
                input->mod->omega_min[1]       = 0.0;
                input->mod->omega_max[1]       = 20.0;
		input->mod->omega_proba[0]     = 0.5;
		input->mod->omega_proba[1]     = 0.5;
		input->mod->ns                 = input->mod->c_code->n_sense_c;
	      }
	    else if(input->mod->subst_modelname == M1)
	      {
		input->mod->subst_modelname    = M1a;
		input->mod->model_number       = 23;
		input->mod->n_catq             = 2;
		input->mod->n_omega            = 2;
		input->mod->s_opt->opt_omega   = 1;
		input->mod->s_opt->opt_p_omega = 1;
		input->mod->omega[0]           = 0.0;
		input->mod->omega[1]           = 1.0;
		input->mod->omega_proba[0]     = 0.5;
		input->mod->omega_proba[1]     = 0.5;
		input->mod->ns                 = input->mod->c_code->n_sense_c;
	      }
	    else if(input->mod->subst_modelname == M1a)
	      {
		input->mod->subst_modelname    = M2;
		input->mod->model_number       = 24;
		input->mod->n_catq             = 3;
		input->mod->s_opt->opt_p_omega = 1;
		input->mod->n_omega            = 3;
		input->mod->s_opt->opt_omega   = 1;
		input->mod->omega[0]           = 0.0;
		input->mod->omega[1]           = 1.0;
		input->mod->omega[2]           = 4.0;
                input->mod->omega_min[0]       = 0.0;
                input->mod->omega_max[0]       = 20.0;
                input->mod->omega_min[1]       = 0.0;
                input->mod->omega_max[1]       = 20.0;
                input->mod->omega_min[2]       = 0.0;
                input->mod->omega_max[2]       = 20.0;
		input->mod->omega_proba[0]     = 0.6;
		input->mod->omega_proba[1]     = 0.3;
		input->mod->omega_proba[2]     = 0.1;
		input->mod->ns                 = input->mod->c_code->n_sense_c;
	      }
	  }

	else if((input->mod->switch_modelname == SWITCH_S1) || 
		(input->mod->switch_modelname == SWITCH_S2))
	  {
            
	    if(input->mod->subst_modelname == M2)
	      {
		input->mod->subst_modelname    = M2a;
		input->mod->model_number       = 25;
		input->mod->s_opt->opt_p_omega = 1;
		input->mod->n_catq             = 1;
		input->mod->n_omega            = 3;
		input->mod->s_opt->opt_omega   = 2;
		input->mod->omega[0]           = 0.3;
		input->mod->omega[1]           = 1.0;
		input->mod->omega[2]           = 4.0;
                input->mod->omega_min[0]       = 0.0;
                input->mod->omega_max[0]       = 20.0;
                input->mod->omega_min[1]       = 0.0;
                input->mod->omega_max[1]       = 20.0;
                input->mod->omega_min[2]       = 0.0;
                input->mod->omega_max[2]       = 20.0;
		input->mod->omega_proba[0]     = 0.6;
		input->mod->omega_proba[1]     = 0.3;
		input->mod->omega_proba[2]     = 0.1;
		input->mod->ns                 = input->mod->n_omega*input->mod->c_code->n_sense_c;
	      }
	    else if(input->mod->subst_modelname == M2a)
	      { 
		input->mod->subst_modelname    = M3;
		input->mod->model_number       = 26;
		input->mod->s_opt->opt_p_omega = 1;
		input->mod->n_catq             = 1;
		input->mod->n_omega            = 3;
		input->mod->s_opt->opt_omega   = 3;
		input->mod->omega[0]           = 0.3;
		input->mod->omega[1]           = 1.0;
		input->mod->omega[2]           = 4.0;
                input->mod->omega_min[0]       = 0.0;
                input->mod->omega_max[0]       = 20.0;
                input->mod->omega_min[1]       = 0.0;
                input->mod->omega_max[1]       = 20.0;
                input->mod->omega_min[2]       = 0.0;
                input->mod->omega_max[2]       = 20.0;
		input->mod->omega_proba[0]     = 0.6;
		input->mod->omega_proba[1]     = 0.3;
		input->mod->omega_proba[2]     = 0.1;
		input->mod->ns                 = input->mod->n_omega*input->mod->c_code->n_sense_c;
	      }
	    else if(input->mod->subst_modelname == M3)
	      {
		input->mod->subst_modelname    = MX;
		input->mod->model_number       = 27;
		input->mod->s_opt->opt_p_omega = 1;
		input->mod->n_catq             = 1;
		input->mod->n_omega            = 3;
		input->mod->s_opt->opt_omega   = 3;
		input->mod->omega[0]           = 0.3;
		input->mod->omega[1]           = 1.0;
		input->mod->omega[2]           = 4.0;
                input->mod->omega_min[0]       = 0.0;
                input->mod->omega_max[0]       = 20.0;
                input->mod->omega_min[1]       = 0.0;
                input->mod->omega_max[1]       = 20.0;
                input->mod->omega_min[2]       = 0.0;
                input->mod->omega_max[2]       = 20.0;
		input->mod->omega_proba[0]     = 0.6;
		input->mod->omega_proba[1]     = 0.3;
		input->mod->omega_proba[2]     = 0.1;
		input->mod->ns                 = input->mod->n_omega*input->mod->c_code->n_sense_c;
	      }
	    else if(input->mod->subst_modelname == MX)
	      {
		input->mod->subst_modelname    = M2;
		input->mod->model_number       = 24;
		input->mod->n_catq             = 1;
		input->mod->s_opt->opt_p_omega = 1;
		input->mod->n_omega            = 3;
		input->mod->s_opt->opt_omega   = 1;
		input->mod->omega[0]           = 0.0;
		input->mod->omega[1]           = 1.0;
		input->mod->omega[2]           = 4.0;
		input->mod->omega_proba[0]     = 0.6;
		input->mod->omega_proba[1]     = 0.3;
		input->mod->omega_proba[2]     = 0.1;
		input->mod->ns                 = input->mod->n_omega*input->mod->c_code->n_sense_c;
	      }
	  }
            
	break;
      }

    case 'S' :
      {
	char *c;
	int n_trial;
          
	printf("How many data sets > ");
	c = (char *)mCalloc(T_MAX_LINE,sizeof(char));
	Getstring_Stdin(c);
          
	n_trial = 0;
	while((!atoi(c)) || (atoi(c) < 0))
	  {
	    if(++n_trial > 10) Exit("\nErr : The number of data sets must be a positive integer\n");
	    printf("\nThe number of data sets must be a positive integer\n");
	    printf("Enter a new value > ");
	    Getstring_Stdin(c);
	  }
	input->n_data_sets = atoi(c);
          

	printf("Start at data set [default=1] > ");
	Getstring_Stdin(c);
          
	n_trial = 0;
	while(atoi(c) < 0)
	  {
	    if(++n_trial > 10) Exit("\nErr : The number of data sets must be a positive integer\n");
	    printf("\nThe number of data sets must be a positive integer\n");
	    printf("Enter a new value > ");
	    Getstring_Stdin(c);
	  }
	input->start_at_data_set = atoi(c);
          
	input->n_data_sets += input->start_at_data_set-1;

	Free(c);
	break;
      }
	
    }
  return choix;
}

/*********************************************************/

void Read_Param(char *param_name, fit_double defaut, fit_double b_inf, fit_double b_sup, char *s, fit_double *update_param)
{
  printf("Enter the value for %s [default %s=%.2f] > ",
	 param_name,
	 param_name,
	 defaut);
  Getstring_Stdin(s);
  Check_Param_Value(s,param_name,b_inf,b_sup,update_param);
}

/*********************************************************/

void Check_Param_Value(char *param_val, char *param_name, fit_double b_inf, fit_double b_sup, fit_double *update_param)
{
  int n_trial;

  if(strcmp(param_val,"\0"))
    {
      n_trial = 0;
      while((atof(param_val) < b_inf) || (atof(param_val) > b_sup))
	{
	  if(++n_trial > 10) 
	    {
	      printf("\nErr : %s must be a greater than %f and smaller than %.2f !\n",
		     param_name, b_inf, b_sup);
	      Exit("");
	    }

	  printf("\n%s must be a greater than %f and smaller than %.2f !\n",
		 param_name, b_inf, b_sup);
	  printf("Enter a new value > ");
	  Getstring_Stdin(param_val);
	}
      *(update_param) = (fit_double)atof(param_val);
    }
}

/*********************************************************/

void Command_Line(option *input, char **argv, int argc)
{
  int optionnum;

  optionnum = Search_Option("-help",argv,argc);
  if(optionnum > -1) Help();

  optionnum = Search_Option("-treefile",argv,argc);
  if(optionnum < 0) Exit("\n. You must specify a tree file name (option -treefile)\n\n");
  else
    {
      Read_Cmdl_Option(input->input_tree_file,argv[optionnum+1]);
      input->fp_tree = Openfile(input->input_tree_file,0);
      
    }

  optionnum = Search_Option("-seqfile",argv,argc);
  if(optionnum < 0) Exit("\n. You must specify a sequence file name (option -seqfile)\n\n");
  else
    {
      Read_Cmdl_Option(input->seqfile,argv[optionnum+1]);
      input->fp_seq = Openfile(input->seqfile,0);
    }

  Open_Output_Files(input,input->seqfile);

  /* Set default values */
  input->mod->datatype          = NT;
  input->inputtree              = 1;
  input->tree                   = NULL;
  input->mod->subst_modelname   = HKY85;
  input->mod->switch_modelname  = NO_SWITCH;
  input->mod->model_number      = 4;
  input->interleaved            = 1;
  input->n_data_sets            = 1;
  strcpy(input->use_default_dnds,"Yes");
  strcpy(input->use_default_ssvs_switches,"Yes");
  

  
  optionnum = Search_Option("-type",argv,argc);
  if(optionnum > -1)
    {
      char *s;
      s = (char *)mCalloc(T_MAX_LINE,sizeof(char));
      Read_Cmdl_Option(s,argv[optionnum+1]);
      if(!(strcmp(s,"nt")))      
	{
	  input->mod->datatype = NT;	
	}
      else if(!(strcmp(s,"aa"))) 
	{
	  input->mod->datatype = AA;
	}
      else Exit("\n. Unknown data type (must be 'nt' or 'aa')\n");
      Free(s);
    }
  
  optionnum = Search_Option("-freq",argv,argc);
  if(optionnum > -1)
    {
      char *s;
      s = (char *)mCalloc(T_MAX_LINE,sizeof(char));
      Read_Cmdl_Option(s,argv[optionnum+1]);
      if(!(strcmp(s,"empirical"))) input->mod->s_opt->opt_bfreq = 0;	
      else if(!(strcmp(s,"ml")))   input->mod->s_opt->opt_bfreq = 1;
      else if(!(strcmp(s,"uniform"))) 
	{
	  input->mod->s_opt->opt_bfreq = 0;
	  input->mod->codon_freq       = 0;
	}
      else if(!(strcmp(s,"F3X4"))) 
	{
	  input->mod->s_opt->opt_bfreq = 0;
	  input->mod->codon_freq       = 1;
	}
      else Exit("\n. Unknown frequence option (must be 'empirical' or 'ml' or 'uniform' or 'F3X4')\n");
      Free(s);
    }
  
  optionnum = Search_Option("-codon",argv,argc);
  if(optionnum > -1)
    {
      char *s;
      s = (char *)mCalloc(T_MAX_LINE,sizeof(char));
      Read_Cmdl_Option(s,argv[optionnum+1]);
      if(!(strcmp(s,"yes")))     input->mod->model_applies_to = CODONS;	
      else if(!(strcmp(s,"no"))) input->mod->model_applies_to = NT;
      else Exit("\n. Unknown option (must be 'yes' or 'no')\n");
      Free(s);
    }

  optionnum = Search_Option("-model",argv,argc);
  if(optionnum > -1)
    {
      char *s;
      s = (char *)mCalloc(T_MAX_LINE,sizeof(char));
      Read_Cmdl_Option(s,argv[optionnum+1]);
      if(!(strcmp(s,"JC69")))     
	{
	  input->mod->n_omega          = 0;
	  input->mod->ns               = 4;
	  input->mod->model_number     = 1;
	  input->mod->subst_modelname  = JC69;
	  input->mod->s_opt->opt_kappa = 0;
	}
      else if(!(strcmp(s,"K80")))     
	{
	  input->mod->n_omega         = 0;
	  input->mod->ns              = 4;
	  input->mod->model_number    = 2;
	  input->mod->subst_modelname = K80;
	}
      else if(!(strcmp(s,"F81")))     
	{
	  input->mod->n_omega          = 0;
	  input->mod->ns               = 4;
	  input->mod->model_number     = 3;
	  input->mod->subst_modelname  = F81;
	  input->mod->s_opt->opt_kappa = 0;
	}
      else if(!(strcmp(s,"HKY85")))     
	{
	  input->mod->n_omega         = 0;
	  input->mod->ns              = 4;
	  input->mod->model_number    = 4;
	  input->mod->subst_modelname = HKY85;
	}
      else if(!(strcmp(s,"F84")))     
	{
	  input->mod->n_omega         = 0;
	  input->mod->ns              = 4;
	  input->mod->model_number    = 5;
	  input->mod->subst_modelname = F84;
	}
      else if(!(strcmp(s,"TN93")))     
	{
	  input->mod->n_omega         = 0;
	  input->mod->ns              = 4;
	  input->mod->model_number    = 6;
	  input->mod->subst_modelname = TN93;
	  if(input->mod->s_opt->opt_kappa) input->mod->s_opt->opt_lambda = 1;
	}
      else if(!(strcmp(s,"GTR")))     
	{
	  input->mod->n_omega          = 0;
	  input->mod->ns               = 4;
	  input->mod->model_number     = 7;
	  input->mod->subst_modelname  = GTR;
	  input->mod->s_opt->opt_kappa = 0;
	}
      else if(!(strcmp(s,"Dayhoff")))     
	{
	  input->mod->n_omega         = 0;
	  input->mod->ns              = 20;
	  input->mod->model_number    = 11;
	  input->mod->subst_modelname = DAYHOFF;

	}
      else if(!(strcmp(s,"JTT")))     
	{
	  input->mod->n_omega         = 0;
	  input->mod->ns              = 20;	  
	  input->mod->model_number    = 12;
	  input->mod->subst_modelname = JTT;
	}
      else if(!(strcmp(s,"MtREV")))     
	{
	  input->mod->n_omega         = 0;
	  input->mod->ns              = 20;	  
	  input->mod->model_number    = 13;
	  input->mod->subst_modelname = MTREV;
	}
      else if(!(strcmp(s,"WAG")))     
	{
	  input->mod->n_omega         = 0;
	  input->mod->ns              = 20;	  
	  input->mod->model_number    = 14;
	  input->mod->subst_modelname = WAG;
	}
      else if(!(strcmp(s,"DCMut")))     
	{
	  input->mod->n_omega         = 0;
	  input->mod->ns              = 20;	  
	  input->mod->model_number    = 15;
	  input->mod->subst_modelname = WAG;
	}

      else if(!(strcmp(s,"M0")))     
	{
	  input->mod->s_opt->opt_p_omega = 0;
	  input->mod->ns                 = input->mod->c_code->n_sense_c;
	  input->mod->ns_codon           = input->mod->c_code->n_sense_c;
	  input->mod->stepsize           = 3;
	  input->mod->n_catq             = 1;
	  input->mod->n_omega            = 1;
	  input->mod->n_catg             = 1;
	  input->mod->invar              = 0;
	  input->mod->pinvar             = .0;
	  optionnum = Search_Option("-opt_omega",argv,argc);
	  if(optionnum == -1) input->mod->s_opt->opt_omega = 1;
	  input->mod->subst_modelname    = M0;
	  input->mod->model_number       = 21;
	}
      else if(!(strcmp(s,"M1")))     
	{
	  input->mod->s_opt->opt_p_omega = 1;
	  input->mod->ns                 = input->mod->c_code->n_sense_c;
	  input->mod->ns_codon           = input->mod->c_code->n_sense_c;
	  input->mod->stepsize           = 3;
	  input->mod->n_catq             = 2;
	  input->mod->n_omega            = 2;
	  input->mod->n_catg             = 1;
	  input->mod->invar              = 0;
	  input->mod->pinvar             = .0;
	  optionnum = Search_Option("-opt_omega",argv,argc);
	  if(optionnum == -1) input->mod->s_opt->opt_omega = 1;
	  input->mod->subst_modelname    = M1;
	  input->mod->model_number       = 22;
	}
      else if(!(strcmp(s,"M2")))     
	{
	  input->mod->ns               = input->mod->c_code->n_sense_c;
	  input->mod->ns_codon         = input->mod->c_code->n_sense_c;
	  input->mod->stepsize         = 3;
	  input->mod->n_catq           = 3;
	  input->mod->n_omega          = 3;
	  input->mod->n_catg           = 1;
	  input->mod->invar            = 0;
	  input->mod->pinvar           = .0;
	  optionnum = Search_Option("-opt_omega",argv,argc);
	  if(optionnum == -1) input->mod->s_opt->opt_omega = 1;
	  input->mod->subst_modelname  = M2;
	  input->mod->model_number     = 24;
	}
      else if(!(strcmp(s,"M2a")))     
	{
	  input->mod->ns               = input->mod->c_code->n_sense_c;
	  input->mod->ns_codon         = input->mod->c_code->n_sense_c;
	  input->mod->stepsize         = 3;
	  input->mod->n_catq           = 3;
	  input->mod->n_omega          = 3;
	  input->mod->n_catg           = 1;
	  input->mod->invar            = 0;
	  input->mod->pinvar           = .0;
	  optionnum = Search_Option("-opt_omega",argv,argc);
	  if(optionnum == -1) input->mod->s_opt->opt_omega = 2;
	  input->mod->subst_modelname  = M2a;
	  input->mod->model_number     = 25;
	}
      else if(!(strcmp(s,"M3")))     
	{
	  input->mod->ns               = input->mod->c_code->n_sense_c;
	  input->mod->ns_codon         = input->mod->c_code->n_sense_c;
	  input->mod->stepsize         = 3;
	  input->mod->n_catq           = 3;
	  input->mod->n_omega          = 3;
	  input->mod->n_catg           = 1;
	  input->mod->invar            = 0;
	  input->mod->pinvar           = .0;
	  optionnum = Search_Option("-opt_omega",argv,argc);
	  if(optionnum == -1)  input->mod->s_opt->opt_omega = 3;
	  input->mod->subst_modelname  = M3;
	  input->mod->model_number     = 26;
	}
      else
	Exit("\n. Unknown model name (type fitmodel -h for more informations)\n");
      
      Free(s);
    }
  
  optionnum = Search_Option("-opt_omega",argv,argc);
  if(optionnum > -1)
    {
      char *s;
      s = (char *)mCalloc(T_MAX_LINE,sizeof(char));
      Read_Cmdl_Option(s,argv[optionnum+1]);
      if(!(strcmp(s,"yes")))      {s=s;}	
      else if(!(strcmp(s,"no")))  input->mod->s_opt->opt_omega = 0;
      else Exit("\n. Unknown option (must be 'yes' or 'no')\n");
      Free(s);

    }

  optionnum = Search_Option("-optpw",argv,argc);
  if(optionnum > -1)
    {
      char *s;
      s = (char *)mCalloc(T_MAX_LINE,sizeof(char));
      Read_Cmdl_Option(s,argv[optionnum+1]);
      if(!(strcmp(s,"yes")))      input->mod->s_opt->opt_p_omega = 1;	
      else if(!(strcmp(s,"no")))  input->mod->s_opt->opt_p_omega = 0;
      else Exit("\n. Unknown option (must be 'yes' or 'no')\n");
      Free(s);

    }

  optionnum = Search_Option("-pinvar",argv,argc);
  if(optionnum > -1)
    {
      char *s;
      s = (char *)mCalloc(T_MAX_LINE,sizeof(char));
      Read_Cmdl_Option(s,argv[optionnum+1]);
      Check_Param_Value(s,"pinvar",.0,1.,&(input->mod->pinvar));
      Free(s);
    }

  optionnum = Search_Option("-optpinvar",argv,argc);
  if(optionnum > -1)
    {
      char *s;
      s = (char *)mCalloc(T_MAX_LINE,sizeof(char));
      Read_Cmdl_Option(s,argv[optionnum+1]);
      if(!(strcmp(s,"yes")))      
	{
	  input->mod->s_opt->opt_pinvar = 1;
	  input->mod->invar             = 1;
	}
      else if(!(strcmp(s,"no")))  
	{
	  input->mod->s_opt->opt_pinvar = 0;
	  input->mod->invar             = 0;
	}
      else Exit("\n. Unknown option (must be 'yes' or 'no')\n");
      Free(s);
    }

  optionnum = Search_Option("-kappa",argv,argc);
  if(optionnum > -1)
    {
      char *s;
      s = (char *)mCalloc(T_MAX_LINE,sizeof(char));
      Read_Cmdl_Option(s,argv[optionnum+1]);
      Check_Param_Value(s,"kappa",.00001,100.,&(input->mod->kappa));
      Free(s);
    }

  optionnum = Search_Option("-optkappa",argv,argc);
  if(optionnum > -1)
    {
      char *s;
      s = (char *)mCalloc(T_MAX_LINE,sizeof(char));
      Read_Cmdl_Option(s,argv[optionnum+1]);
      if(!(strcmp(s,"yes")))      input->mod->s_opt->opt_kappa = 1;	
      else if(!(strcmp(s,"no")))  input->mod->s_opt->opt_kappa = 0;
      else Exit("\n. Unknown option (must be 'yes' or 'no')\n");
      Free(s);
    }

  optionnum = Search_Option("-ncatg",argv,argc);
  if(optionnum > -1)
    {
      char *s;
      int n_trial;

      s = (char *)mCalloc(T_MAX_LINE,sizeof(char));
      Read_Cmdl_Option(s,argv[optionnum+1]);
      n_trial = 0;
      while((!atoi(s)) || (atoi(s) < 0))
	{
	  if(++n_trial > 10) Exit("\nErr : the number of categories must be a positive integer\n");
	  printf("\nThe number of categories must be a positive integer\n");
	  printf("Enter a new value > ");
	  Getstring_Stdin(s);
	}
      input->mod->n_catg = atoi(s);
      Free(s);
    }


  optionnum = Search_Option("-alpha",argv,argc);
  if(optionnum > -1)
    {
      char *s;
      s = (char *)mCalloc(T_MAX_LINE,sizeof(char));
      Read_Cmdl_Option(s,argv[optionnum+1]);
      Check_Param_Value(s,"alpha",.01,100.,&(input->mod->alpha));
      Free(s);
    }

  optionnum = Search_Option("-optalpha",argv,argc);
  if(optionnum > -1)
    {
      char *s;
      s = (char *)mCalloc(T_MAX_LINE,sizeof(char));
      Read_Cmdl_Option(s,argv[optionnum+1]);
      if(!(strcmp(s,"yes")))      input->mod->s_opt->opt_alpha = 1;	
      else if(!(strcmp(s,"no")))  input->mod->s_opt->opt_alpha = 0;
      else Exit("\n. Unknown option (must be 'yes' or 'no')\n");
      Free(s);
    }

  optionnum = Search_Option("-multiple",argv,argc);
  if(optionnum > -1)
    {
      char *s;
      int n_trial;

      s = (char *)mCalloc(T_MAX_LINE,sizeof(char));
      Read_Cmdl_Option(s,argv[optionnum+1]);
     
      n_trial = 0;
      while((!atoi(s)) || (atoi(s) < 0))
	{
	  if(++n_trial > 10) Exit("\nErr : The number of data sets must be a positive integer\n");
	  printf("\nThe number of data sets must be a positive integer\n");
	  printf("Enter a new value > ");
	  Getstring_Stdin(s);
	}
      input->n_data_sets = atoi(s);
      Free(s);
    }

  optionnum = Search_Option("-optall",argv,argc);
  if(optionnum > -1)
    {
      char *s;
      s = (char *)mCalloc(T_MAX_LINE,sizeof(char));
      Read_Cmdl_Option(s,argv[optionnum+1]);
      if(!(strcmp(s,"yes")))      input->mod->s_opt->opt_param = 1;	
      else if(!(strcmp(s,"no")))  input->mod->s_opt->opt_param = 0;
      else Exit("\n. Unknown option (must be 'yes' or 'no')\n");
      Free(s);
    }
  
  optionnum = Search_Option("-code",argv,argc);
  if(optionnum > -1)
    {
      char *s;
      int code,n_trial;
      int code_num;

      s = (char *)mCalloc(T_MAX_LINE,sizeof(char));
      Read_Cmdl_Option(s,argv[optionnum+1]);
      n_trial = 0;
      while((!atoi(s)) || (atoi(s) < 0))
	{
	  if(++n_trial > 10) Exit("\nErr : The number of the genetic code must be a positive integer\n");
	  printf("\nThe number of the genetic code must be a positive integer\n");
	  printf("Enter a new value > ");
	  Getstring_Stdin(s);
	}
      code = atoi(s);

      /* code = transl_table (see http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi?mode=c#thetop) */

      switch(code)
	{
	case 1 : 
	  {
	    input->mod->c_code->num_curr_code = 0;
	    break;
	  }
	case 2 : 
	  {
	    input->mod->c_code->num_curr_code = 1;
	    break;
	  }
	case 3 : 
	  {
	    input->mod->c_code->num_curr_code = 2;
	    break;
	  }
	case 4 : 
	  {
	    input->mod->c_code->num_curr_code = 3;
	    break;
	  }
	case 5 : 
	  {
	    input->mod->c_code->num_curr_code = 4;
	    break;
	  }
	case 6 : 
	  {
	    input->mod->c_code->num_curr_code = 5;
	    break;
	  }
	case 9 : 
	  {
	    input->mod->c_code->num_curr_code = 6;
	    break;
	  }
	case 10 : 
	  {
	    input->mod->c_code->num_curr_code = 7;
	    break;
	  }
	case 11 : 
	  {
	    input->mod->c_code->num_curr_code = 8;
	    break;
	  }
	case 12 : 
	  {
	    input->mod->c_code->num_curr_code = 9;
	    break;
	  }
	case 13 : 
	  {
	    input->mod->c_code->num_curr_code = 10;
	    break;
	  }
	case 14 : 
	  {
	    input->mod->c_code->num_curr_code = 11;
	    break;
	  }
	case 15 : 
	  {
	    input->mod->c_code->num_curr_code = 12;
	    break;
	  }
	case 16 : 
	  {
	    input->mod->c_code->num_curr_code = 13;
	    break;
	  }
	case 21 : 
	  {
	    input->mod->c_code->num_curr_code = 14;
	    break;
	  }
	case 22 : 
	  {
	    input->mod->c_code->num_curr_code = 15;
	    break;
	  }
	case 23 : 
	  {
	    input->mod->c_code->num_curr_code = 16;
	    break;
	  }

	default : {Exit("\n. Unknown genetic code.\n");}
	}

      code_num = input->mod->c_code->num_curr_code;
      Free_Code(input->mod->c_code);
      input->mod->c_code   = Make_Code(code_num);
      input->mod->ns       = input->mod->c_code->n_sense_c;
      input->mod->ns_codon = input->mod->c_code->n_sense_c;

    }
  
  optionnum = Search_Option("-switches",argv,argc);
  if(optionnum > -1)
    {
      char *s;

      s = (char *)mCalloc(T_MAX_LINE,sizeof(char));
      Read_Cmdl_Option(s,argv[optionnum+1]);     
      if(!(strcmp(s,"no")))       
	{
	  printf("\n. Fix needed in file %s at line %d\n",__FILE__,__LINE__);
	  Exit("");
		
	  input->mod->switch_modelname = NO_SWITCH;
	  input->mod->ns               = input->mod->c_code->n_sense_c;
	  input->mod->model_number     = 20;
	  input->mod->n_catq           = 3;
	}
      else if(!(strcmp(s,"S1")))
	{
	  printf("\n. Fix needed in file %s at line %d\n",__FILE__,__LINE__);
	  Exit("");

	  input->mod->switch_modelname = SWITCH_S1;
	  input->mod->n_omega          = 3;
	  input->mod->ns               = input->mod->n_omega*input->mod->c_code->n_sense_c;
	  input->mod->model_number     = 21;
	  input->mod->n_catq           = 1;
	  input->mod->s_opt->opt_theta = 1;
	}

      else if(!(strcmp(s,"S2")))
	{
	  printf("\n. Fix needed in file %s at line %d\n",__FILE__,__LINE__);
	  Exit("");

	  input->mod->switch_modelname = SWITCH_S2;
	  input->mod->n_omega          = 3;
	  input->mod->ns               = input->mod->n_omega*input->mod->c_code->n_sense_c;
	  input->mod->model_number     = 21;
	  input->mod->n_catq           = 1;
	  input->mod->s_opt->opt_theta = 3;
	}
      else Exit("\n. Unknown option (must be 'no' or 'S1' or 'S2')\n");
      Free(s);
    }
  
  optionnum = Search_Option("-optdelta_sw",argv,argc);
  if(optionnum > -1)
    {
      char *s;
      s = (char *)mCalloc(T_MAX_LINE,sizeof(char));

      Read_Cmdl_Option(s,argv[optionnum+1]);

      if(!(strcmp(s,"yes")))      
	{
	  if(input->mod->switch_modelname == SWITCH_S1)
	    input->mod->s_opt->opt_theta = 1;
	  else if(input->mod->switch_modelname == SWITCH_S2)
	    input->mod->s_opt->opt_theta = 3;
	}

      else if(!(strcmp(s,"no"))) 
	{
	  input->mod->s_opt->opt_theta = 0;
	}

      else Exit("\n. Unknown option (must be 'yes' or 'no')\n");
      Free(s);
    }

  /*   optionnum = Search_Option("-optalpha_sw",argv,argc); */
  /*   if(optionnum > -1) */
  /*     { */
  /*       char *s; */
  /*       s = (char *)mCalloc(T_MAX_LINE,sizeof(char)); */
  /*       Read_Cmdl_Option(s,argv[optionnum+1]); */
  /*       if(!(strcmp(s,"yes")))      input->mod->s_opt->opt_alpha_sw = 1;	 */
  /*       else if(!(strcmp(s,"no")))  input->mod->s_opt->opt_alpha_sw = 0; */
  /*       else Exit("\n. Unknown option (must be 'yes' or 'no')\n"); */
  /*       Free(s); */
  /*     } */

  /*   optionnum = Search_Option("-optbeta_sw",argv,argc); */
  /*   if(optionnum > -1) */
  /*     { */
  /*       char *s; */
  /*       s = (char *)mCalloc(T_MAX_LINE,sizeof(char)); */
  /*       Read_Cmdl_Option(s,argv[optionnum+1]); */
  /*       if(!(strcmp(s,"yes")))      input->mod->s_opt->opt_beta_sw = 1;	 */
  /*       else if(!(strcmp(s,"no")))  input->mod->s_opt->opt_beta_sw = 0; */
  /*       else Exit("\n. Unknown option (must be 'yes' or 'no')\n"); */
  /*       Free(s); */
  /*     } */
  

  optionnum = Search_Option("-w1",argv,argc);
  if(optionnum > -1)
    {
      char *s;
      s = (char *)mCalloc(T_MAX_LINE,sizeof(char));
      Read_Cmdl_Option(s,argv[optionnum+1]);
      Check_Param_Value(s,"w1",.00000001,10000000.,&(input->mod->omega[0]));
      Free(s);
    }

  optionnum = Search_Option("-w2",argv,argc);
  if(optionnum > -1)
    {
      char *s;
      s = (char *)mCalloc(T_MAX_LINE,sizeof(char));
      Read_Cmdl_Option(s,argv[optionnum+1]);
      Check_Param_Value(s,"w2",.00000001,10000000.,&(input->mod->omega[1]));
      Free(s);
    }

  optionnum = Search_Option("-w3",argv,argc);
  if(optionnum > -1)
    {
      char *s;
      s = (char *)mCalloc(T_MAX_LINE,sizeof(char));
      Read_Cmdl_Option(s,argv[optionnum+1]);
      Check_Param_Value(s,"w3",.00000001,10000000.,&(input->mod->omega[2]));
      Free(s);
    }

  optionnum = Search_Option("-p1",argv,argc);
  if(optionnum > -1)
    {
      char *s;
      s = (char *)mCalloc(T_MAX_LINE,sizeof(char));
      Read_Cmdl_Option(s,argv[optionnum+1]);
      Check_Param_Value(s,"p1",.0,1.,&(input->mod->omega_proba[0]));
      Free(s);
    }

  optionnum = Search_Option("-p2",argv,argc);
  if(optionnum > -1)
    {
      char *s;
      s = (char *)mCalloc(T_MAX_LINE,sizeof(char));
      Read_Cmdl_Option(s,argv[optionnum+1]);
      Check_Param_Value(s,"p2",.0,1.,&(input->mod->omega_proba[1]));
      Free(s);
    }

  optionnum = Search_Option("-p3",argv,argc);
  if(optionnum > -1)
    {
      char *s;
      s = (char *)mCalloc(T_MAX_LINE,sizeof(char));
      Read_Cmdl_Option(s,argv[optionnum+1]);
      Check_Param_Value(s,"p3",.0,1.,&(input->mod->omega_proba[2]));
      Free(s);
    }

  Check_Inconsistencies(input);
}

/*********************************************************/

int Search_Option(char *option_name, char **argv, int argc)
{
  int i;
  For(i,argc) if(strstr(argv[i],option_name)) return i;
  return -1;
}

/*********************************************************/

void Read_Cmdl_Option(char *optionname, char *cmdl)
{
  char *line;
  int i,pos;

  line = (char *)mCalloc((int)strlen(cmdl)+1,sizeof(char));
  
  pos = 0;
  For(i,(int)strlen(cmdl))
    {
      if(cmdl[i] != ' ')
	{
	  line[pos] = cmdl[i];
	  pos++;
	}
    }
  line[(int)strlen(line)]='\0';
  strcpy(optionname,line);
  Free(line);
}

/*********************************************************/
 
void Open_Output_Files(option *input,char *seqfile)
{
  
  char answer, *file_name;
  int n_trial;

  file_name = (char *)mCalloc(T_MAX_FILE,sizeof(char));


  strcpy(input->output_stat_file,input->seqfile);
  strcpy(input->output_tree_file,input->seqfile);
  
#ifdef PHYML
  strcat(input->output_stat_file,"_phyml_stats");
  strcat(input->output_tree_file,"_phyml_tree");
#endif
    
#ifdef FITMODEL
  strcat(input->output_stat_file,"_fitmodel_stats");
  strcat(input->output_tree_file,"_fitmodel_tree");
#endif
  
#ifdef EVOLVE
  int pid;
  pid = (int)getpid();
  sprintf(input->output_stat_file+strlen(input->output_stat_file),".%d",pid);
  sprintf(input->output_tree_file+strlen(input->output_tree_file),".%d",pid);
  strcat(input->output_stat_file,".evolve_stats");
  strcat(input->output_tree_file,".evolve_seq");
#endif

#ifdef RF
  strcat(input->output_stat_file,"_rf_stats");
  strcat(input->output_tree_file,"_rf_tree");
#endif

#ifdef POSSELSITEID
  strcat(input->output_stat_file,"_ps_stats");
  strcat(input->output_tree_file,"_ps_tree");
#endif

#ifdef WIN32
  strcat(input->output_stat_file,".txt");
  strcat(input->output_tree_file,".txt");
#endif

  printf("The first output file will be called '%s'\n",input->output_stat_file);
  printf("Do you want to keep that name ? [Y/n] ");
  if(!scanf("%c", &answer))
    {
      printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Exit("");
    }
  if(answer == '\n') answer = 'Y';
  else getchar();
  n_trial = 0;
  while((answer != 'Y') && (answer != 'y') &&
	(answer != 'N') && (answer != 'n'))  
    {
      if(++n_trial > 10) Exit("\nErr : wrong answers !");
      printf("Do you want to keep that name ? [Y/n] ");
      if(!scanf("%c", &answer))
	{
	  printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	  Exit("");
	}
      if(answer == '\n') answer = 'Y';
      else getchar();
    }
      
  if((answer == 'N') || (answer == 'n'))
    {
      printf("Enter the name you want then > ");
      Getstring_Stdin(file_name);
      strcpy(input->output_stat_file,file_name);
    }
      
  printf("\n");
  printf("The second output file will be called '%s'\n",input->output_tree_file);
  printf("Do you want to keep that name ? [Y/n] ");
  if(!scanf("%c", &answer))
    {
      printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Exit("");
    }
  if(answer == '\n') answer = 'Y';
  else getchar();
  n_trial = 0;
  while((answer != 'Y') && (answer != 'y') &&
	(answer != 'N') && (answer != 'n'))  
    {
      if(++n_trial > 10) Exit("\nErr : wrong answers !");
      printf("Do you want to keep that name ? [Y/n] ");
      if(!scanf("%c", &answer))
	{
	  printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	  Exit("");
	}
      if(answer == '\n') answer = 'Y';
      else getchar();
    }
      
  if((answer == 'N') || (answer == 'n'))
    {
      printf("Enter the name you want then > ");
      Getstring_Stdin(file_name);
      strcpy(input->output_tree_file,file_name);
    }
            
  Free(file_name);
      
}

/*********************************************************/

void Check_Inconsistencies(option *input)
{


  if(input->mod->datatype == AA)
    {
      if(input->mod->s_opt->opt_bfreq)
	Exit("\n. Inconsistency in your settings (data type vs. nucleotide freq.)\n");

      if(input->mod->model_number < 10)
	Exit("\n. Inconsistency in your settings (data type vs. model) \n");

      if(input->mod->s_opt->opt_kappa)
	Exit("\n. Inconsistency in your settings (data type vs. optimise kappa) \n");

      if(input->mod->model_applies_to == CODONS)
	Exit("\n. Inconsistency in your settings (data type vs. codon-based model) \n");

    }

  if(input->mod->s_opt->opt_bfreq)
    {
      if((input->mod->model_number > 10) || ((input->mod->model_number < 3)))
	Exit("\n. Inconsistency in your settings (nucleotide freq. vs. model) \n");

      if(input->mod->datatype == AA)
	Exit("\n. Inconsistency in your settings (nucleotide freq. vs. data type) \n");

      if(input->mod->model_applies_to == CODONS)
	Exit("\n. Inconsistency in your settings (nucleotide freq. vs. codon-based model) \n");

    }


  if(input->mod->model_applies_to == CODONS)
    {
      if(input->mod->datatype == AA)
	Exit("\n. Inconsistency in your settings (codon-based model vs. data type) \n");
	
      if(input->mod->s_opt->opt_bfreq)
	Exit("\n. Inconsistency in your settings (codon-based model vs. nucleotide freq.)\n");

      if(input->mod->model_number < 20)
	Exit("\n. Inconsistency in your settings (codon-based model vs. model) \n");

      if(input->mod->pinvar > 0.0)
	Exit("\n. Inconsistency in your settings (codon-based model vs. pinvar) \n");

      if(input->mod->s_opt->opt_pinvar)
	Exit("\n. Inconsistency in your settings (codon-based model vs. optimise pinvar) \n");

      if(input->mod->s_opt->opt_alpha)
	Exit("\n. Inconsistency in your settings (codon-based model vs. optimise alpha) \n");



      int i;
      For(i,input->mod->n_omega)
	{
	  if(input->mod->omega[i] < input->mod->omega_min[i] ||
	     input->mod->omega[i] > input->mod->omega_max[i])
	    {
	      printf("\n. w%d=%f min=%f max=%f",i+1,input->mod->omega[i],input->mod->omega_min[i],input->mod->omega_max[i]);
	      Exit("\n. Inconsistency in your settings (value of omega is out of its user-defined boundaries) \n");
	    }
	}
    }


  switch(input->mod->model_number)
    {
    case 1 : case 7 : case 11 : case 12 : case 13 : case 14 : case 15 :
      {
	if(input->mod->s_opt->opt_kappa)
	  Exit("\n. Inconsistency in your settings (model vs. optimise kappa) \n");
      }
    default : break;
    }



  if(input->mod->model_number >= 20)
    {
      if(input->mod->pinvar > 0.0)
	Exit("\n. Inconsistency in your settings (model vs. pinvar > 0.0) \n");

      if(input->mod->s_opt->opt_pinvar)
	Exit("\n. Inconsistency in your settings (model vs. optimise pinvar) \n");

      if(input->mod->n_catg > 1)
	Exit("\n. Inconsistency in your settings (model vs. optimise n_catg > 1) \n");

      if(input->mod->s_opt->opt_alpha)
	Exit("\n. Inconsistency in your settings (model vs. optimise alpha) \n");
    }
}

/*********************************************************/

#define BOLD "\033[00;01m"
#define FLAT "\033[00;00m"
#define LINE "\033[00;04m"

void Help()
{
  printf(BOLD"NAME\n");
  printf(FLAT"\tFitModel %s\n\n",RELEASE);

  printf(BOLD"SYNOPSIS\n");
  printf(FLAT"\tfitmodel -treefile "LINE"treefilename"FLAT" -seqfile "LINE"seqfilename"FLAT" ["LINE"options"FLAT"]\n\n");

  printf(BOLD"DESCRIPTION\n");
  printf(FLAT"\tFitmodel fits substitution models to sequence data sets (given in "LINE"seqfilename"FLAT")\n\ton a fixed topology (given in "LINE"treefilename"FLAT") using Maximum Likelihood.\n\n");

  printf(BOLD"COMMAND LINE OPTIONS\n");
  printf(BOLD"\t-type "LINE"nt"FLAT" or "LINE"aa"FLAT" (default=nt)\n");
  printf("\n");

  printf(BOLD"\t-freq "LINE"empirical"FLAT" or "LINE"ml"FLAT" or "LINE"uniform"FLAT" or "LINE"F3X4"FLAT" (defaults=empirical or F3X4)\n");
  printf("\n");

  printf(BOLD"\t-codon "LINE"no"FLAT" or "LINE"yes"FLAT" (defaults=no)\n");
  printf("\n");

  printf(BOLD"\t-model "LINE"JC69"FLAT", "LINE"K80"FLAT", "LINE"F81"FLAT", "LINE"HKY85"FLAT", "LINE"F84"FLAT", "LINE"TN93"FLAT", "LINE"GTR"FLAT", "LINE"Dayhoff"FLAT", "LINE"JTT"FLAT", "LINE"MtREV"FLAT", "LINE"WAG"FLAT", "LINE"DCMut"FLAT", "LINE"M2"FLAT" or "LINE"M3"FLAT" (default=HKY85)\n");
  printf("\n");

  printf(BOLD"\t-pinvar "LINE"[0.0;1.0]\n");
  printf("\n");
  
  printf(BOLD"\t-optpinvar "LINE"no"FLAT" or "LINE"yes"FLAT" (defaults=no)\n");
  printf("\n");

  printf(BOLD"\t-kappa "LINE"[0.01;100.0]\n");
  printf("\n");
  
  printf(BOLD"\t-optkappa "LINE"no"FLAT" or "LINE"yes"FLAT" (default=no)\n");
  printf("\n");

  printf(BOLD"\t-ncatg "LINE"integer > 0\n");
  printf("\n");

  printf(BOLD"\t-alpha "LINE"[0.01;100.0]\n");
  printf("\n");
  
  printf(BOLD"\t-optalpha "LINE"no"FLAT" or "LINE"yes"FLAT" (default=no)\n");
  printf("\n");

  printf("\n");

  printf(BOLD"\t-code "LINE"1,2,3,4,5,6,9,10,11,12,13,14,15,16,21,22,23"FLAT" (see NCBI Taxonomy webpage)"FLAT" (default=yes)\n");
  printf("\n");

  printf(BOLD"\t-p1 "LINE"[0.0;1.0]\n");
  printf("\n");

  printf(BOLD"\t-p2 "LINE"[0.0;1.0]\n");
  printf("\n");

  printf(BOLD"\t-p3 "LINE"[0.0;1.0]\n");
  printf("\n");

  printf(BOLD"\t-w1 "LINE"[1E-7;1E+7]\n");
  printf("\n");

  printf(BOLD"\t-w2 "LINE"[1E-7;1E+7]\n");
  printf("\n");

  printf(BOLD"\t-w3 "LINE"[1E-7;1E+7]\n");
  printf("\n");

  printf(BOLD"\t-switches "LINE"no"FLAT" or "LINE"S1"FLAT" or "LINE"S2"FLAT" (default=no)\n");
  printf("\n");


  printf(BOLD"\t-optpw "LINE"yes"FLAT" or "LINE"no"FLAT" (default=yes)\n");
  printf("\n");



  printf("\n");
  printf(BOLD"\t-multiple "LINE"integer > 0\n");
  printf("\n");

  printf(BOLD"\t-interleaved "LINE"yes"FLAT" or "LINE"no"FLAT" (default=yes)\n");
  printf("\n");
  
  printf(BOLD"\t-optall "LINE"yes"FLAT" or "LINE"no"FLAT" (default=yes)\n");
  printf("\n");
  
  Exit("\n");
}

/*********************************************************/
 
char *Return_Subst_Model_Name(model *mod)
{
  char *s;

  s = (char *)mCalloc(T_MAX_MODEL_NAME,sizeof(char));

  switch(mod->subst_modelname)
    {
    case JC69 : 
      {
	strcpy(s,"JC69");
	break;
      }
    case K80 : 
      {
	strcpy(s,"K80");
	break;
      }
    case F81 : 
      {
	strcpy(s,"F81");
	break;
      }
    case HKY85 : 
      {
	strcpy(s,"HKY85");
	break;
      }
    case F84 : 
      {
	strcpy(s,"F84");
	break;
      }
    case TN93 : 
      {
	strcpy(s,"TN93");
	break;
      }
    case GTR : 
      {
	strcpy(s,"GTR");
	break;
      }
    case DAYHOFF : 
      {
	strcpy(s,"Dayhoff");
	break;
      }
    case JTT : 
      {
	strcpy(s,"JTT");
	break;
      }
    case MTREV : 
      {
	strcpy(s,"MtREV");
	break;
      }
    case WAG : 
      {
	strcpy(s,"WAG");
	break;
      }
    case DCMUT : 
      {
	strcpy(s,"DCMut");
	break;
      }
    case M0 : 
      {
	strcpy(s,"M0");
	break;
      }
    case M1 : 
      {
	strcpy(s,"M1");
	break;
      }
    case M1a : 
      {
	strcpy(s,"M1a");
	break;
      }
    case M2 : 
      {
	strcpy(s,"M2");
	break;
      }
    case M2a : 
      {
	strcpy(s,"M2a");
	break;
      }
    case M3 : 
      {
	strcpy(s,"M3");
	break;
      } 
    case MX : 
      {
	strcpy(s,"MX");
	break;
      }
    case NO_SWITCH : 
      {
	strcpy(s,"No switch");
	break;
      }
    case SWITCH_S1 : 
      {
	strcpy(s,"S1");
	break;
      }
    case SWITCH_S2 : 
      {
	strcpy(s,"S2");
	break;
      }
    default : 
      {
	printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	Exit("");
      }
    }
  return s;
}

/*********************************************************/

char *Return_Target_Selection_Class(option *input)
{
  char *s;
  
  s = (char *)mCalloc(T_MAX_MODEL_NAME,sizeof(char));

  switch(input->sel_class_target)
    {
    case 0 :
      {
	strcpy(s,"first\0");
	break;
      }
    case 1 :
      {
	strcpy(s,"second\0");
	break;
      }
    case 2 :
      {
	strcpy(s,"third\0");
	break;
      }
    default : break;
    }

  return s;
}

/*********************************************************/

char *Return_Switch_Model_Name(model *mod)
{
  char *s;

  s = (char *)mCalloc(T_MAX_MODEL_NAME,sizeof(char));

  switch(mod->switch_modelname)
    {
    case NO_SWITCH :
      {
	strcpy(s,"No switch");
	break;
      }
    case SWITCH_S1 : 
      {
	strcpy(s,"S1");
	break;
      }
    case SWITCH_S2 : 
      {
	strcpy(s,"S2");
	break;
      }
    default : 
      {
	printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	Exit("");
      }
    }
  return s;
}

/*********************************************************/

char *Return_Model_Applies_To(model *mod)
{
  char *s;

  s = (char *)mCalloc(T_MAX_MODEL_NAME,sizeof(char));

  switch(mod->model_applies_to)
    {
    case AA :
      {
	strcpy(s,"amino-acids");
	break;
      }
    case CODONS : 
      {
	strcpy(s,"codons");
	break;
      }
    case NT : 
      {
	strcpy(s,"nucleotides");
	break;
      }
    default : 
      {
	printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	Exit("");
      }
    }
  return s;
}

/*********************************************************/
/*********************************************************/
/*********************************************************/
/*********************************************************/
/*********************************************************/
/*********************************************************/
/*********************************************************/
/*********************************************************/
/*********************************************************/
