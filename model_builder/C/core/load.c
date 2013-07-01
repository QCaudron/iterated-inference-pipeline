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



void load2u_var(unsigned int **tab, char *filename, unsigned int *Np)
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
          if(p<Np[n])
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




void load3d_var(double ***tab, int n, unsigned int *colbreaks1, unsigned int **colbreaks2, char *filename)
{
  /*we have to divided lines into different sections*/

  FILE* myfile =NULL;
  myfile=fopen(filename, "r");

  /*nesting: (i|c|ac) i encompass c that encompass ac*/
  int i=0;
  int c=0;
  int ac=0;
  int column=0;
  double temp;

  unsigned int *column_size = init1u_set0(n);

  for(i=0;i<n;i++)
    {
      column_size[i]=0;
      for(c=0;c<colbreaks1[i];c++)
        {
          column_size[i]+=colbreaks2[i][c];
        }
    }

  i=0;
  c=0;
  ac=0;

  if(myfile != NULL)
    {
      while(fscanf(myfile,"%lf",&temp)!=EOF)
        {
          if(column<column_size[i])
            {
              column++;
              if(ac<colbreaks2[i][c])
                {
                  tab[i][c][ac++]=temp;
                }
              else
                {
                  ac=0;
                  tab[i][++c][ac++]=temp;
                }
            }
          else
            {
              column=1;
              ac=0;
              c=0;
              tab[++i][c][ac++]=temp;
            }
        }/*fin du while*/

      fclose(myfile);
      myfile=NULL;
    }
  else
    {
      fprintf(stderr, "failed to read %s\n", filename);
    }

  FREE(column_size);
}


void load3u_var(unsigned int ***tab, int n, unsigned int *colbreaks1, unsigned int **colbreaks2, char *filename)
{
  /*we have to divided lines into different sections*/
  FILE* myfile =NULL;
  myfile=fopen(filename, "r");
  /*nesting: (i|c|ac) i encompass c that encompass ac*/
  int i=0;
  int c=0;
  int ac=0;
  int column=0;
  unsigned int temp;

  unsigned int *column_size = init1u_set0(n);

  for(i=0;i<n;i++)
    {
      column_size[i]=0;
      for(c=0;c<colbreaks1[i];c++)
        {
          column_size[i]+=colbreaks2[i][c];
        }
    }

  i=0;
  c=0;
  ac=0;

  if(myfile != NULL)
    {
      while(fscanf(myfile,"%d",&temp)!=EOF)
        {
          if(column<column_size[i])
            {
              column++;
              if(ac<colbreaks2[i][c])
                {
                  tab[i][c][ac++]=temp;
                }
              else
                {
                  ac=0;
                  tab[i][++c][ac++]=temp;
                }
            }
          else
            {
              column=1;
              ac=0;
              c=0;
              tab[++i][c][ac++]=temp;
            }
        }/*fin du while*/

      fclose(myfile);
      myfile=NULL;
    }
  else
    {
      fprintf(stderr, "failed to read %s\n", filename);
    }

  FREE(column_size);
}


void load3u_varp1(unsigned int ***tab, int n, unsigned int *colbreaks1, unsigned int colbreaks2, char *filename)
{
  int i, j;

  unsigned int **temp = init2u_var_set0(n, colbreaks1); /*we want to use load3u_var() so we crete a temporary array*/

  for(i=0; i<n; i++)
    for(j=0; j<colbreaks1[i]; j++)
      temp[i][j]=colbreaks2;

  load3u_var(tab, n, colbreaks1, temp, filename);

  clean2u(temp, n);
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

    offset = 0;

    for(i=0; i<p_it_all->length; i++) {
        const char *par_key = routers[i]->name;

        json_t *par = fast_get_json_object(parameters, par_key);
        json_t *follow = json_object_get(par, "follow");
        if(follow){
            const char *par_follow_key = fast_get_json_string_from_object(par, "follow");
            par = fast_get_json_object(parameters, par_follow_key);
        }
	json_t *groups = fast_get_json_object(par, "group");

        for (g=0; g < routers[i]->n_gp; g++) {
	    const char *group_key = routers[i]->group_name[g];
	    json_t *group = fast_get_json_object(groups, group_key);


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
	    
            offset++;
	}
    }

    if (json_object_get(theta, "covariance")) {
        load_covariance(p_best->var, fast_get_json_array(theta, "covariance"));
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


void load_covariance(gsl_matrix *covariance, json_t *array2d)
{
    char str[STR_BUFFSIZE];
    int i, k;

    for (i=1; i< json_array_size(array2d); i++) {//start at 1 to skip header
        json_t *array_i;
        array_i = json_array_get(array2d, i);
        if (!json_is_array(array_i)) {
            sprintf(str, "error: covariance[%d] is not an array\n", i);
            print_err(str);

            exit(EXIT_FAILURE);
        }

        for (k=0; k< json_array_size(array_i); k++) {
            json_t *value_ik;
            value_ik = json_array_get(array_i, k);

            if (json_is_number(value_ik)) {
                gsl_matrix_set(covariance, i-1, k, json_number_value(value_ik)); //i-1 due to header
            } else {
                sprintf(str, "error: covariance[%d][%d] is not a number\n", i, k);
                print_err(str);
                exit(EXIT_FAILURE);
            }
        }
    }
}
