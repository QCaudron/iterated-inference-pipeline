/**************************************************************************
 *    This file is part of plom.
 *
 *    plom is free software: you can redistribute it and/or modify it
 *    under the terms of the GNU General Public License as published
 *    by the Free Software Foundation, either version 3 of the
 *    License, or (at your option) any later version.
 *
 *    plom is distributed in the hope that it will be useful, but
 *    WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public
 *    License along with plom.  If not, see
 *    <http://www.gnu.org/licenses/>.
 *************************************************************************/

#include "plom.h"


unsigned int n2d(char *filename, int Np)
{
  /*return number of lines of filename */
  int n,p;
  n=1;
  p=0;
  double tmp=0.0;
  FILE* myfile = NULL;
  myfile = fopen(filename, "r");

  if( myfile != NULL )
    {
      while( (fscanf(myfile,"%lf",&tmp)!=EOF) )
        {
          if(p<Np)
            {
              p++;
            }
          else
            {
              n++;
              p=0;
              p++;
            }
        }
      fclose(myfile);
      myfile=NULL;
    }
  else
    {
      fprintf(stderr, "failed to read %s, the program will now quit\n", filename);
      exit(EXIT_FAILURE);
    }

  return n;
}

unsigned int n2d_nan(char *filename, int Np)
{
  /*return number of lines of filename allowing possible NaN in the file */
  int n,p;
  n=1;
  p=0;
  char tmp[255];
  FILE* myfile = NULL;
  myfile = fopen(filename, "r");

  if( myfile != NULL )
    {
      while( (fscanf(myfile,"%s",tmp)!=EOF) )
        {
          if(p<Np)
            {
              p++;
            }
          else
            {
              n++;
              p=0;
              p++;
            }
        }
      fclose(myfile);
      myfile=NULL;
    }
  else
    {
      fprintf(stderr, "failed to read %s, the program will now quit\n", filename);
      exit(EXIT_FAILURE);
    }

  return n;
}



void load1u(unsigned int *tab, char *filename)
{
  int i;
  unsigned int tmp=0;

  FILE* myfile = NULL;
  
  i=0;
  myfile = fopen(filename, "r");
  if (myfile != NULL)
    {
      while(fscanf(myfile,"%u",&tmp)!=EOF)
        {
          tab[i++]=tmp;
        }
      fclose(myfile);
      myfile=NULL;
    }
  else
    {
      fprintf(stderr, "failed to read %s\n", filename);
    }
}


void load1d(double *tab, char *filename)
{
  int i;
  double tmp=0.0;

  FILE* myfile = NULL;

  i=0;
  myfile = fopen(filename, "r");
  if (myfile != NULL)
    {
      while(fscanf(myfile,"%lf",&tmp)!=EOF)
        {
          tab[i++]=tmp;
        }
      fclose(myfile);
      myfile=NULL;
    }
  else
    {
      fprintf(stderr, "failed to read %s\n", filename);
    }

}


void load2u(unsigned int **tab, char *filename, int Np)
{
  int n,p;
  n=0;
  p=0;
  unsigned int tmp=0;

  FILE* myfile = NULL;
  myfile = fopen(filename, "r");

  if( myfile != NULL )
    {
      while( (fscanf(myfile,"%d",&tmp)!=EOF) )
        {
          if(p<Np)
            {
              tab[n][p++]=tmp;
            }
          else
            {
              n++;
              p=0;
              tab[n][p++]=tmp;
            }
        }
      fclose(myfile);
      myfile=NULL;
    }
  else
    {
      fprintf(stderr, "failed to read %s\n", filename);
    }
}


void load2d(double **tab, char *filename, int Np)
{

  int n,p;
  n=0;
  p=0;
  double tmp=0.0;
  FILE* myfile = NULL;
  myfile = fopen(filename, "r");

  if( myfile != NULL )
    {
      while( (fscanf(myfile,"%lf",&tmp)!=EOF) )
        {
          if(p<Np)
            {
              tab[n][p++]=tmp;
            }
          else
            {
              n++;
              p=0;
              tab[n][p++]=tmp;
            }
        }
      fclose(myfile);
      myfile=NULL;
    }
  else
    {
      fprintf(stderr, "failed to read %s\n", filename);
    }
}



void load2d_nan(double **tab, char *filename, int Np)
{
  int n,p;
  n=0;
  p=0;

  char tmp[255];
  char *end;

  FILE* myfile = NULL;
  myfile = fopen(filename, "r");

  if( myfile != NULL )
    {
      while( (fscanf(myfile,"%s", tmp)!=EOF) )
        {
          if(p<Np)
            {
              tab[n][p++]=strtod(tmp, &end);
            }
          else
            {
              n++;
              p=0;
              tab[n][p++]=strtod(tmp, &end);
            }
        }
      fclose(myfile);
      myfile=NULL;
    }
  else
    {
      fprintf(stderr, "failed to read %s\n", filename);
    }

}





/**
 *load constants defined as global variable
 */
json_t *load_settings(const char *path)
{
    json_error_t settings_error;
    json_t *settings = json_load_file(path, 0, &settings_error);
    if(!settings) {
        print_err(settings_error.text);
        exit(EXIT_FAILURE);
    }

#if FLAG_VERBOSE
    print_log("load plom settings...");
#endif

    POP_SIZE_EQ_SUM_SV = fast_get_json_boolean(settings, "POP_SIZE_EQ_SUM_SV");

    json_t *cst = fast_get_json_object(settings, "cst");

    /*dimensions parameters*/
    N_C = fast_get_json_integer(cst, "N_C");
    N_AC = fast_get_json_integer(cst, "N_AC");
    N_CAC = N_C*N_AC;
    N_PAR_PROC = fast_get_json_integer(cst, "N_PAR_PROC");
    N_PAR_OBS = fast_get_json_integer(cst, "N_PAR_OBS");
    N_PAR_SV = fast_get_json_integer(cst, "N_PAR_SV");
    N_PAR_FIXED = fast_get_json_integer(cst, "N_PAR_FIXED");
    N_TS = fast_get_json_integer(cst, "N_TS");
    N_TS_INC = fast_get_json_integer(cst, "N_TS_INC");
    N_TS_INC_UNIQUE = fast_get_json_integer(cst, "N_TS_INC_UNIQUE");
    N_DATA = fast_get_json_integer(cst, "N_DATA");
    //N_DATA_NONAN is computed and assigned in build_data
    N_OBS_ALL = fast_get_json_integer(cst, "N_OBS_ALL");
    N_OBS_INC = fast_get_json_integer(cst, "N_OBS_INC");
    N_OBS_PREV = fast_get_json_integer(cst, "N_OBS_PREV");

    N_DRIFT = fast_get_json_integer(cst, "N_DRIFT");

    IS_SCHOOL_TERMS = fast_get_json_integer(cst, "IS_SCHOOL_TERMS");
    
    return settings;
}

/**
 * integrate best data from the webApp
 */
void load_best(struct s_best *p_best, struct s_data *p_data, json_t *theta,  int update_guess)
{
    int i, k, g, offset;

    struct s_router **routers = p_data->routers;
    json_t *parameters = fast_get_json_object(theta, "parameter");
    struct s_iterator *p_it_all = p_data->p_it_all;

    for(i=0; i<p_it_all->length; i++) {
        const char *par_key = routers[i]->name;

        json_t *par = fast_get_json_object(parameters, par_key);

	if(par){
	    json_t *follow = json_object_get(par, "follow");
	    if(follow){
		const char *par_follow_key = fast_get_json_string_from_object(par, "follow");
		par = fast_get_json_object(parameters, par_follow_key);
	    }
	    json_t *groups = fast_get_json_object(par, "group");

	    for (g=0; g < routers[i]->n_gp; g++) {
		const char *group_key = routers[i]->group_name[g];
		json_t *group = fast_get_json_object(groups, group_key);
		offset = p_it_all->offset[i]+g;

		const char *prior = fast_get_json_string_from_object(fast_get_json_object(group, "prior"), "value");
	    
		if (strcmp(prior, "normal") == 0) {
		    p_best->prior[offset] = &normal_prior;
		} else if (strcmp(prior, "pseudo_uniform") == 0) {
		    p_best->prior[offset] = &pseudo_unif_prior;
		} else {
		    p_best->prior[offset] = &gsl_ran_flat_pdf;
		}
	    
		if (update_guess) {
		    gsl_vector_set(p_best->mean, offset, fast_get_json_real_from_object(fast_get_json_object(group, "guess"), "value"));
		}

		gsl_matrix_set(p_best->var, offset, offset, (follow) ? 0.0 : fast_get_json_real_from_object(fast_get_json_object(group, "sd_transf"), "value"));

		p_best->par_prior[offset][0] = fast_get_json_real_from_object(fast_get_json_object(group, "min"), "value");
		p_best->par_prior[offset][1] = fast_get_json_real_from_object(fast_get_json_object(group, "max"), "value");	   
	    }
	} else { //theta remainder
	    
	    for (g=0; g < routers[i]->n_gp; g++) {
		offset = p_it_all->offset[i]+g;
		p_best->prior[offset] = &gsl_ran_flat_pdf;
		if (update_guess) {
		    gsl_vector_set(p_best->mean, offset, 0.0);
		}
		gsl_matrix_set(p_best->var, offset, offset, 0.0);
		p_best->par_prior[offset][0] = 0.0;
		p_best->par_prior[offset][1] = 1.0;	   
	    }
	}
    }

    if (json_object_get(theta, "covariance")) {
        load_covariance(p_best->var, fast_get_json_array(theta, "covariance"), p_data);
    }

    //zero terms of covariance matrix corresponding to env noises parameters
    if(p_data->noises_off & PLOM_NO_ENV_STO){
        struct s_iterator *p_it_noise = p_data->p_it_noise;
        for(i=0; i<p_it_noise->length; i++){
            for(g=0; g< routers[ p_it_noise->ind[i] ]->n_gp; g++){
		offset = p_it_noise->offset[i]+g;
		for(k=0; k< p_best->var->size1; k++){
		    gsl_matrix_set(p_best->var, offset, k, 0.0); // set row to 0
		    gsl_matrix_set(p_best->var, k, offset, 0.0); // set column to 0
		}
            }
        }
    }

    //zero terms of covariance matrix corresponding to the volatilities
    if(p_data->noises_off & PLOM_NO_DRIFT){
        for(i=0; i< N_DRIFT; i++) { //N_DRIFT as iterator could have been edited in case the user or a method imposed PLOM_NO_DRIFT (--no_drift)
            int ind_volatility = p_data->drift[i]->ind_volatility_Xdrift;
            for(g=0; g< routers[ind_volatility]->n_gp; g++) {
		offset = p_it_all->offset[ind_volatility]+g;
		for(k=0; k< p_best->var->size1; k++){
		    gsl_matrix_set(p_best->var, offset, k, 0.0); // set row to 0
		    gsl_matrix_set(p_best->var, k, offset, 0.0); // set column to 0
		}
            }
        }
    }
}


/**
 * covjson do not contain entry for theta remainder but covariance
 * does => pad with zeros where needed
 */
void load_covariance(gsl_matrix *covariance, json_t *covjson, struct s_data *p_data)
{
    char str[STR_BUFFSIZE];
    int i, k, offset, ii, kk, ooffset, indcovi, indcovk;
    struct s_iterator *p_it = p_data->p_it_all;
    struct s_router **routers = p_data->routers;

    int ind_theta_remainder = (p_data->p_it_theta_remainder->length) ? p_data->p_it_theta_remainder->ind[0]: -1;

    indcovi = 1;  //1 to skip header
    for (i=0; i< p_it->length; i++) {
	for(k=0; k< routers[ p_it->ind[i] ]->n_gp; k++) {
	    offset = p_it->offset[i]+k;

	    json_t *array_i;
	    if(i != ind_theta_remainder){
		array_i = json_array_get(covjson, indcovi);
		if (!json_is_array(array_i)) {
		    sprintf(str, "error: covariance[%d] is not an array\n", i);
		    print_err(str);
		    exit(EXIT_FAILURE);
		}
		indcovi++;
	    }

	    indcovk = 0;
	    for (ii=0; ii< p_it->length; ii++) {
		for(kk=0; kk< routers[ p_it->ind[ii] ]->n_gp; kk++) {
		    ooffset = p_it->offset[ii]+kk;
		
		    json_t *value_ik;

		    if((i != ind_theta_remainder) && (ii != ind_theta_remainder)){
			value_ik = json_array_get(array_i, indcovk);

			if (json_is_number(value_ik)) {
			    gsl_matrix_set(covariance, offset, ooffset, json_number_value(value_ik));
			} else {
			    sprintf(str, "error: covariance[%d][%d] is not a number\n", i, k);
			    print_err(str);
			    exit(EXIT_FAILURE);
			}
			indcovk++;
		    } else {
			gsl_matrix_set(covariance, offset, ooffset, 0.0);
		    }
		}
	    }
	}
    }
}
