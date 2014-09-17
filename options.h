/*

FitModeL :  a program that Fits substitution models to the data
using Maximum Likelihood.  

Copyright (C) Stephane Guindon. jul. 2004 onward

All parts of  the source except where indicated  are distributed under
the GNU public licence.  See http://www.opensource.org for details.

*/

#ifndef OPTIONS_H
#define OPTIONS_H

void Usage();
option *Get_Input(int argc,char **argv);
char Print_Menu_Nt_AA(option *input,char *s);
char Print_Menu_Codon(option *input,char *s);
void Read_Param(char *name_param, fit_double defaut, fit_double b_inf, fit_double b_sup, char *s, fit_double *update_param);
void Command_Line(option *input,char **argv, int argc);
int Search_Option(char *option_name, char **argv, int argc);
void Read_Cmdl_Option(char *optionname, char *cmdl);
void Open_Output_Files(option *input,char *seqfile);
void Check_Param_Value(char *param_val, char *param_name, fit_double b_inf, fit_double b_sup, fit_double *update_param);
void Check_Inconsistencies(option *input);
void Help();
char *Return_Switch_Model_Name(model *mod);
char *Return_Subst_Model_Name(model *mod);
char *Return_Model_Applies_To(model *mod);
char *Return_Target_Selection_Class(option *input);

#endif
