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


//these functions provide syntaxic shortcuts to work with jansson

json_t *fast_get_json_object(const json_t *container, const char *obj_name)
{
  json_t *object;
  object = json_object_get(container, obj_name);
  if(!json_is_object(object))
    {
      char str[STR_BUFFSIZE];
      sprintf(str, "error: %s is not an object\n", obj_name);
      print_err(str);
      exit(EXIT_FAILURE);
    }

  return object;
}


double fast_get_json_real_from_object(json_t *object, const char *key)
{
    json_t *my_real;
    my_real = json_object_get(object, key);
    if(!json_is_number(my_real))
        {
            char str[STR_BUFFSIZE];
            sprintf(str, "error: %s is not a number\n", key);
            print_err(str);
            exit(EXIT_FAILURE);
        }

    return json_number_value(my_real);
}


const char *fast_get_json_string_from_object(const json_t *object, const char *key)
{
    json_t *string_value;
    string_value = json_object_get(object, key);

    if (!json_is_string(string_value) ) {
        char str[STR_BUFFSIZE];
        sprintf(str, "error: %s is not a string\n", key);
        print_err(str);
        exit(EXIT_FAILURE);
    }

    return json_string_value(string_value);
}



json_t *fast_get_json_array(const json_t *container, const char *array_name)
{
  json_t *array;
  array = json_object_get(container, array_name);
  if(!json_is_array(array))
    {
      char str[STR_BUFFSIZE];
      sprintf(str, "error: %s is not an array\n", array_name);
      print_err(str);
      exit(EXIT_FAILURE);
    }

  return array;
}


double fast_get_json_real_from_array(json_t *array, int i, const char *array_name)
{
  json_t *my_real;
  my_real = json_array_get(array, i);
  if(!json_is_number(my_real))
    {
      char str[STR_BUFFSIZE];
      sprintf(str, "error: %s[%d] is not a number\n", array_name, i);
      print_err(str);
      exit(EXIT_FAILURE);
    }

  return json_number_value(my_real);
}

const char *fast_get_json_string_from_array(json_t *array, int i, const char *array_name)
{
    json_t *key = json_array_get(array, i);
    if (!json_is_string(key)) {
        char str[STR_BUFFSIZE];
        sprintf(str, "error: %s[%d] is not a string\n", array_name, i);
        print_err(str);
        exit(EXIT_FAILURE);
    }

    return json_string_value(key);
}


int fast_get_json_boolean(json_t *container, char *obj_name)
{
  json_t *tmp;
  tmp = json_object_get(container, obj_name);
  if(!json_is_boolean(tmp))
    {
      char str[STR_BUFFSIZE];
      sprintf(str, "error: %s is not a boolean\n", obj_name);
      print_err(str);
      exit(EXIT_FAILURE);
    }

  return (int) json_is_true(tmp);
}



int fast_get_json_integer(json_t *container, char *obj_name)
{
  json_t *tmp;
  tmp = json_object_get(container, obj_name);
  if(!json_is_number(tmp))
    {
      char str[STR_BUFFSIZE];
      sprintf(str, "error: %s is not a number\n", obj_name);
      print_err(str);
      exit(EXIT_FAILURE);
    }

  return (int) json_integer_value(tmp);
}



double fast_get_json_real(json_t *container, char *obj_name)
{
  json_t *tmp;
  tmp = json_object_get(container, obj_name);
  if(!json_is_number(tmp))
    {
      char str[STR_BUFFSIZE];
      sprintf(str, "error: %s is not a number\n", obj_name);
      print_err(str);
      exit(EXIT_FAILURE);
    }

  return json_number_value(tmp);
}


unsigned int *fast_load_fill_json_1u(json_t *array, char *array_name)
{
    char str[STR_BUFFSIZE];

    int i;
    unsigned int *tab = malloc(json_array_size(array) * sizeof(unsigned int));
    if(tab==NULL) {
        sprintf(str, "Allocation impossible in file :%s line : %d",__FILE__,__LINE__);
        print_err(str);
        exit(EXIT_FAILURE);
    }

    for(i=0; i< json_array_size(array); i++) {
        json_t *array_i;
        array_i = json_array_get(array, i);

        if(!json_is_number(array_i)){
            sprintf(str, "error: %s[%d] is not a number\n", array_name, i);
            print_err(str);
            exit(EXIT_FAILURE);
        }
        tab[i] = (int) json_integer_value(array_i);

    }

    return tab;
}


unsigned int **fast_load_fill_json_2u(json_t *array, char *array_name)
{
    char str[STR_BUFFSIZE];
    int i;
    char array_name_i[STR_BUFFSIZE];

    unsigned int **tab = malloc(json_array_size(array) * sizeof(unsigned int *));
    if(tab==NULL) {
        sprintf(str, "Allocation impossible in file :%s line : %d",__FILE__,__LINE__);
        print_err(str);
        exit(EXIT_FAILURE);
    }

    for(i=0; i< json_array_size(array); i++) {
        json_t *array_i;
        array_i = json_array_get(array, i);
        if(!json_is_array(array_i)) {
            sprintf(str, "error: %s[%d] is not an array\n", array_name, i);
            print_err(str);
            exit(EXIT_FAILURE);
        }
        sprintf(array_name_i, "%s[%d]", array_name, i);
        tab[i] = fast_load_fill_json_1u(array_i, array_name_i);
    }

    return tab;

}

unsigned int ***fast_load_fill_json_3u(json_t *array, char *array_name)
{
    char str[STR_BUFFSIZE];
    int i;
    char array_name_i[STR_BUFFSIZE];

    unsigned int ***tab = malloc(json_array_size(array) * sizeof(unsigned int **));
    if(tab==NULL) {
        sprintf(str, "Allocation impossible in file :%s line : %d",__FILE__,__LINE__);
        print_err(str);
        exit(EXIT_FAILURE);
    }


    for(i=0; i< json_array_size(array); i++) {
        json_t *array_i;
        array_i = json_array_get(array, i);
        if(!json_is_array(array_i)) {
            sprintf(str, "error: %s[%d] is not an array\n", array_name, i);
            print_err(str);
            exit(EXIT_FAILURE);
        }
        sprintf(array_name_i, "%s[%d]", array_name, i);
        tab[i] = fast_load_fill_json_2u(array_i, array_name_i);
    }

    return tab;

}



double *fast_load_fill_json_1d(json_t *array, char *array_name)
{
    /* handles properly null values */

    char str[STR_BUFFSIZE];
    int i;
    double *tab = malloc(json_array_size(array) * sizeof(double));
    if(tab==NULL){
        sprintf(str, "Allocation impossible in file :%s line : %d",__FILE__,__LINE__);
        print_err(str);
        exit(EXIT_FAILURE);
    }

    for(i=0; i< json_array_size(array); i++) {
        json_t *array_i;
        array_i = json_array_get(array, i);

        if(json_is_number(array_i)) {
            tab[i] = json_number_value(array_i);
        }
        else if(json_is_null(array_i)) {
            tab[i] = NAN;
        }
        else {
            sprintf(str, "error: %s[%d] is not a number nor null\n", array_name, i);
            print_err(str);
            exit(EXIT_FAILURE);
        }
    }

    return tab;

}

double **fast_load_fill_json_2d(json_t *array, char *array_name)
{
    /* handles properly null values */
    char str[STR_BUFFSIZE];
    int i;
    char array_name_i[STR_BUFFSIZE];

    double **tab = malloc(json_array_size(array) * sizeof(double *));
    if(tab==NULL) {
        sprintf(str, "Allocation impossible in file :%s line : %d",__FILE__,__LINE__);
        print_err(str);
        exit(EXIT_FAILURE);
    }

    for(i=0; i< json_array_size(array); i++){
        json_t *array_i;
        array_i = json_array_get(array, i);
        if(!json_is_array(array_i)) {
            sprintf(str, "error: %s[%d] is not an array\n", array_name, i);
            print_err(str);

            exit(EXIT_FAILURE);
        }
        sprintf(array_name_i, "%s[%d]", array_name, i);
        tab[i] = fast_load_fill_json_1d(array_i, array_name_i);
    }

    return tab;

}





double ***fast_load_fill_json_3d(json_t *array, char *array_name)
{
    /* handles properly null values */
    char str[STR_BUFFSIZE];
    int i;
    char array_name_i[STR_BUFFSIZE];

    double ***tab = malloc(json_array_size(array) * sizeof(double **));
    if(tab==NULL) {
        sprintf(str, "Allocation impossible in file :%s line : %d",__FILE__,__LINE__);
        print_err(str);
        exit(EXIT_FAILURE);
    }


    for(i=0; i< json_array_size(array); i++) {
        json_t *array_i;
        array_i = json_array_get(array, i);
        if(!json_is_array(array_i)) {
            sprintf(str, "error: %s[%d] is not an array\n", array_name, i);
            print_err(str);
            exit(EXIT_FAILURE);
        }
        sprintf(array_name_i, "%s[%d]", array_name, i);
        tab[i] = fast_load_fill_json_2d(array_i, array_name_i);
    }

    return tab;

}


json_t *load_json(void)
{
    char *buffer;
    buffer = malloc(BUFFER_SIZE* sizeof(char));

    fgets(buffer, BUFFER_SIZE, stdin);
//    printf("%s\n", buffer);

    json_t *root;
    json_error_t error;

    root = json_loads(buffer, 0, &error);
    if(!root) {
        char str[STR_BUFFSIZE];
        sprintf(str, "could not parse theta.json\nerror: on line %d: %s\n", error.line, error.text);
        print_err(str);
        exit(EXIT_FAILURE);
    }

    //  buffer = json_dumps(root, 0);

    free(buffer);

    return root;
}



/**
 * return index of element in the json array of string array
 */
int index_of_json_array(const json_t *array, const char *element){

    int k = 0;
    while((k < json_array_size(array)) && strcmp(json_string_value(json_array_get(array, k)), element)) {
        k++;
    }

    if(k == json_array_size(array)) {
        return -1;
    }

    return k;
}


json_t *plom_theta_remainder_new(json_t* theta){
    int g;
    
    //create group object
    json_t *par_group = json_object();
    json_t *partition = fast_get_json_object(fast_get_json_object(theta, "partition"), "variable_population");
    json_t *group = fast_get_json_array(partition, "group");

    int n_gp = json_array_size(group);

    for (g=0; g < n_gp; g++) { //for each group
        json_t *my_group = json_array_get(group, g); 
        const char *my_group_id = fast_get_json_string_from_object(my_group, "id");

	json_t *values = json_pack("{s:{s:f}, s:{s:f}, s:{s:f}, s:{s:f}, s:{s:s}}", "min", "value", 0.0, "max", "value", 1.0, "guess", "value", 0.0, "sd_transf", "value", 0.0, "prior", "value", "uniform");
	json_object_set_new(par_group, my_group_id, values);
    }

    return json_pack("{s:s, s:s, s:o}", "partition_id", "variable_population", "transformation", "logit", "group", par_group);
}
