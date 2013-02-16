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

struct s_iterator *build_iterator(json_t *settings, struct s_router **routers, struct s_drift **drift, char *it_type, const enum plom_noises_off noises_off)
{
    int i, k;

    struct s_iterator *p_it;
    p_it = malloc(sizeof(struct s_iterator));
    if (p_it==NULL) {
        char str[STR_BUFFSIZE];
        snprintf(str, STR_BUFFSIZE, "Allocation impossible in file :%s line : %d",__FILE__,__LINE__);
        print_err(str);
        exit(EXIT_FAILURE);
    }

    int is_drift = ! (noises_off & PLOM_NO_DRIFT);

    //length and ind
    if (strcmp(it_type, "all") == 0) {

        p_it->length = N_PAR_SV + N_PAR_PROC + N_PAR_OBS;

        if (p_it->length) {
            p_it->ind = init1u_set0(p_it->length);

            for (i=0; i<p_it->length; i++) {
                p_it->ind[i] = i;
            }
        }

    } else if (strcmp(it_type, "only_drift") == 0) {
	
	p_it->length = (is_drift)? N_DRIFT: 0;

        if (p_it->length) {
            p_it->ind = init1u_set0(p_it->length);

            for (i=0; i<p_it->length; i++) {
                p_it->ind[i] = drift[i]->ind_par_Xdrift_applied;
            }
        }

    } else if (strcmp(it_type, "par_sv") == 0) {

        p_it->length = N_PAR_SV;

        if (p_it->length) {
            p_it->ind = init1u_set0(p_it->length);

            for (i=0; i<p_it->length; i++) {
                p_it->ind[i] = i;
            }
        }


    } else if (strcmp(it_type, "all_no_drift") == 0) {

        p_it->length = N_PAR_SV + N_PAR_PROC + N_PAR_OBS;
	if(is_drift){
	    p_it->length -= N_DRIFT;
	}
	
        if (p_it->length) {
            p_it->ind = init1u_set0(p_it->length);

            k = 0;
            for (i=0; i < N_PAR_SV + N_PAR_PROC + N_PAR_OBS; i++) {
                if(! is_drift) {
                    p_it->ind[k] = i;
		    k++;
		} else {
		    if(! in_drift(i, drift)) {
			p_it->ind[k] = i;
			k++;
		    }
		}
            }

        }

    } else if (strcmp(it_type, "par_proc_par_obs_no_drift") == 0) {

        p_it->length = N_PAR_PROC + N_PAR_OBS;
	if(is_drift){
	    p_it->length -= N_DRIFT;
	}

        if (p_it->length) {
            p_it->ind = init1u_set0(p_it->length);
            k = 0;
            for (i=N_PAR_SV; i< N_PAR_SV + N_PAR_PROC + N_PAR_OBS; i++) {
                if(! is_drift) {
                    p_it->ind[k] = i;
		    k++;
		} else {
		    if(! in_drift(i, drift)) {
			p_it->ind[k] = i;
			k++;
		    }
		}
            }
        }

    } else if (strcmp(it_type, "par_sv_and_drift") == 0) {

        p_it->length = N_PAR_SV;
	if(is_drift){
	    p_it->length += N_DRIFT;
	}

        if (p_it->length) {
            p_it->ind = init1u_set0(p_it->length);

            for (i=0; i< N_PAR_SV; i++) {
                p_it->ind[i] = i;
            }

	    if(is_drift){
		for (i=0; i< N_DRIFT; i++) {
		    p_it->ind[N_PAR_SV+i] = drift[i]->ind_par_Xdrift_applied;
		}
	    }
        }

    } else if (strcmp(it_type, "noise") == 0) {

        json_t *ind_noise_sd = fast_get_json_array(settings, "ind_noise_sd");
        p_it->length = json_array_size(ind_noise_sd);
        if (p_it->length) {
            p_it->ind = fast_load_fill_json_1u(ind_noise_sd, "ind_noise_sd");
        }
    }

    //alloc ind and offset and assign nbtot
    p_it->nbtot = 0; // will be incremented later if p_it->length but need to be assigned systematically!

    if (p_it->length) {
        unsigned int *all_offset = init1u_set0(N_PAR_SV+N_PAR_PROC+N_PAR_OBS);
        for(i=1; i< (N_PAR_SV+N_PAR_PROC+N_PAR_OBS); i++) {
            all_offset[i] = all_offset[i-1] + routers[ i-1 ]->n_gp; //cumsum
        }

        p_it->offset = init1u_set0(p_it->length);
        for(i=0; i<p_it->length; i++) {
            p_it->offset[i] = all_offset[ p_it->ind[i] ];
            p_it->nbtot += routers[ p_it->ind[i] ]->n_gp;
        }
        FREE(all_offset);
    }

    return p_it;
}


void clean_iterator(struct s_iterator *p_it)
{
    if (p_it->length) {
        FREE(p_it->ind);
        FREE(p_it->offset);
    }

    FREE(p_it);
}


struct s_router *build_router(const json_t *par, const char *par_key, const json_t *partition, const json_t *order, const char *link_key, const char *u_data, int is_bayesian)
{
    char str[STR_BUFFSIZE];
    int g, k;

    struct s_router *p_router;
    p_router = malloc(sizeof(struct s_router));
    if (p_router == NULL) {
        snprintf(str, STR_BUFFSIZE, "Allocation impossible in file :%s line : %d",__FILE__,__LINE__);
        print_err(str);
        exit(EXIT_FAILURE);
    }

    //store the name
    p_router->name = init1c(strlen(par_key) + 1);
    strcpy(p_router->name, par_key);

    json_t *groups = fast_get_json_array(partition, "group");
    p_router->p = json_array_size(order);

    p_router->n_gp = json_array_size(groups);

    p_router->map = init1u_set0(p_router->p);
    p_router->min = init1d_set0(p_router->n_gp);
    p_router->max = init1d_set0(p_router->n_gp);
    p_router->min_z = init1d_set0(p_router->n_gp);
    p_router->max_z = init1d_set0(p_router->n_gp);

    //alloc for group_name
    p_router->group_name = malloc(p_router->n_gp * sizeof(char *));
    if (p_router->group_name == NULL) {
        snprintf(str, STR_BUFFSIZE, "Allocation impossible in file :%s line : %d",__FILE__,__LINE__);
        print_err(str);
        exit(EXIT_FAILURE);
    }

    json_t *par_min = fast_get_json_object(par, "min");
    json_t *par_max = fast_get_json_object(par, "max");

    for (g=0; g < p_router->n_gp; g++) { //for each group

        json_t *my_group = json_array_get(groups, g); //TODO check if we get an object...

        const char *my_group_id = fast_get_json_string_from_object(my_group, "id");

        p_router->group_name[g] = init1c(strlen(my_group_id) + 1);
        strcpy(p_router->group_name[g], my_group_id);

        p_router->min[g] = fast_get_json_real_from_object(par_min, my_group_id);
        p_router->max[g] = fast_get_json_real_from_object(par_max, my_group_id);

        json_t *compo = fast_get_json_array(my_group, link_key); //link_key is population_id or time_series_id
        for (k = 0; k< json_array_size(compo); k++) {
            const char * element = fast_get_json_string_from_array(compo, k, link_key);
            int index = index_of_json_array(order, element);
            if(index != -1){
                p_router->map[index] = g;
            } else {
                char str[STR_BUFFSIZE];
                snprintf(str, STR_BUFFSIZE, "could not find element '%s' in array: '%s'", element, link_key);
                print_err(str);
                exit(EXIT_FAILURE);
            }
        }
    }

    set_f_trans(p_router, par, u_data, is_bayesian);
    set_ab_z(p_router);

    return p_router;
}


void clean_router(struct s_router *p_router)
{
    int g;
    for (g=0; g < p_router->n_gp; g++) {
        FREE(p_router->group_name[g]);
    }
    FREE(p_router->group_name);
    FREE(p_router->name);

    FREE(p_router->min);
    FREE(p_router->max);
    FREE(p_router->min_z);
    FREE(p_router->max_z);
    FREE(p_router->map);
    FREE(p_router);
}

struct s_router **build_routers(json_t *settings, json_t *theta, int is_bayesian)
{
    int i, j, offset;
    json_t *parameters = fast_get_json_object(theta, "value");
    json_t *partitions = fast_get_json_object(theta, "partition");
    json_t *orders = fast_get_json_object(settings, "orders");
    const char *frequency = fast_get_json_string_from_object(fast_get_json_object(settings, "cst"), "FREQUENCY");

    const char par_types[][10] = { "par_sv", "par_proc", "par_obs" };
    const char pop_ts_types[][20] = { "cac_id", "cac_id", "ts_id" };
    const char link_types[][20] = { "population_id", "population_id", "time_series_id" };

    struct s_router **routers;
    routers = malloc((N_PAR_SV+N_PAR_PROC+N_PAR_OBS)*sizeof(struct s_router *));
    if (routers==NULL) {
        char str[STR_BUFFSIZE];
        snprintf(str, STR_BUFFSIZE, "Allocation impossible in file :%s line : %d",__FILE__,__LINE__);
        print_err(str);
        exit(EXIT_FAILURE);
    }

    offset = 0;
    for(i=0; i<3; i++) {
        json_t *my_par_list = fast_get_json_array(orders, par_types[i]);

        for(j=0; j< json_array_size(my_par_list); j++) {

            const char *par_key = fast_get_json_string_from_array(my_par_list, j, par_types[i]);
            json_t *par = fast_get_json_object(parameters, par_key);

            //overwrite par if follower
            if(json_object_get(par, "follow")){
                const char *par_follow_key = fast_get_json_string_from_object(par, "follow");

                //check that parameter types of follower and follow are compatibles
                int is_par = ((index_of_json_array(fast_get_json_array(orders, "par_sv"), par_follow_key) != -1) || (index_of_json_array(fast_get_json_array(orders, "par_proc"), par_follow_key) != -1));
                int is_obs = (index_of_json_array(fast_get_json_array(orders, "par_obs"), par_follow_key) != -1);

                if( (!is_par && !is_obs) || ((i==0) && is_obs) || ((i==1) && is_obs) || ((i==2) && is_par) ){
                    char str[STR_BUFFSIZE];
                    snprintf(str, STR_BUFFSIZE, "invalid parameter type between %s (follower) and %s (follow)", par_key, par_follow_key);
                    print_err(str);
                    exit(EXIT_FAILURE);
                } else {
                    par = fast_get_json_object(parameters, par_follow_key);
                }
            }

            const char *partition_key = fast_get_json_string_from_object(par, "partition_id");
            routers[offset] = build_router(par, par_key,
                                           fast_get_json_object(partitions, partition_key),
                                           fast_get_json_array(orders, pop_ts_types[i]),
                                           link_types[i],
                                           frequency,
                                           is_bayesian);
            offset++;
        }
    }

    return routers;
}



void clean_routers(struct s_router **routers)
{
    int i;
    for(i=0; i< (N_PAR_SV+N_PAR_PROC+N_PAR_OBS); i++) {
        clean_router(routers[i]);
    }

    FREE(routers);
}

struct s_par *build_par(struct s_data *p_data)
{
    char str[STR_BUFFSIZE];
    int i;
    struct s_par *p_par;
    p_par = malloc(sizeof(struct s_par));
    if (p_par==NULL) {
        char str[STR_BUFFSIZE];
        snprintf(str, STR_BUFFSIZE, "Allocation impossible in file :%s line : %d",__FILE__,__LINE__);
        print_err(str);
        exit(EXIT_FAILURE);
    }

    p_par->size_natural = N_PAR_SV+N_PAR_PROC+N_PAR_OBS;

    p_par->natural = malloc(p_par->size_natural* sizeof (double *));
    if (p_par->natural==NULL) {
        snprintf(str, STR_BUFFSIZE, "Allocation impossible in file :%s line : %d",__FILE__,__LINE__);
        print_err(str);
        exit(EXIT_FAILURE);
    }

    for(i=0; i<p_par->size_natural; i++) {
        p_par->natural[i] = init1d_set0(p_data->routers[i]->n_gp);
    }

    return p_par;
}


void clean_par(struct s_par *p_par)
{
    clean2d(p_par->natural, p_par->size_natural);

    FREE(p_par);
}


struct s_par **build_J_p_par(struct s_data *p_data)
{
    int j;

    struct s_par **J_p_par;
    J_p_par = malloc(J * sizeof(struct s_par *));
    if (J_p_par==NULL) {
        char str[STR_BUFFSIZE];
        snprintf(str, STR_BUFFSIZE, "Allocation impossible in file :%s line : %d",__FILE__,__LINE__);
        print_err(str);
        exit(EXIT_FAILURE);
    }

    for(j=0; j<J; j++) {
        J_p_par[j] = build_par(p_data);
    }

    return J_p_par;
}

void clean_J_p_par(struct s_par **J_p_par)
{
    int j;
    for(j=0; j<J; j++) {
        clean_par(J_p_par[j]);
    }

    FREE(J_p_par);
}


struct s_obs2ts **build_obs2ts(json_t *json_obs2ts)
{
    int o;

    struct s_obs2ts **obs2ts;
    obs2ts = malloc(N_OBS_ALL *sizeof(struct s_obs2ts *));
    if (obs2ts==NULL) {
        char str[STR_BUFFSIZE];
        sprintf(str, "Allocation impossible in file :%s line : %d",__FILE__,__LINE__);
        print_err(str);
        exit(EXIT_FAILURE);
    }

    for(o=0; o<N_OBS_ALL; o++) {
        obs2ts[o] = malloc(sizeof(struct s_obs2ts));
        if (obs2ts[o]==NULL) {
            char str[STR_BUFFSIZE];
            snprintf(str, STR_BUFFSIZE, "Allocation impossible in file :%s line : %d",__FILE__,__LINE__);
            print_err(str);
            exit(EXIT_FAILURE);
        }
    }

    for(o=0; o<N_OBS_ALL; o++) {
        json_t *json_obs2ts_o = json_array_get(json_obs2ts, o);

        obs2ts[o]->n_ts_unique = fast_get_json_integer(json_obs2ts_o, "n_ts_unique");
        obs2ts[o]->n_stream = fast_load_fill_json_1u(fast_get_json_array(json_obs2ts_o, "n_stream"), "n_stream");
        obs2ts[o]->n_cac = fast_load_fill_json_1u(fast_get_json_array(json_obs2ts_o, "n_cac"), "n_cac");
        obs2ts[o]->cac = fast_load_fill_json_3u(fast_get_json_array(json_obs2ts_o, "cac"), "cac");
        obs2ts[o]->offset = fast_get_json_integer(json_obs2ts_o, "offset");
    }

    return obs2ts;
}

void clean_obs2ts(struct s_obs2ts **obs2ts)
{
    int o;

    for(o=0; o<N_OBS_ALL; o++) {
        clean3u_var(obs2ts[o]->cac, obs2ts[o]->n_ts_unique, obs2ts[o]->n_cac);
        FREE(obs2ts[o]->n_stream);
        FREE(obs2ts[o]->n_cac);
        FREE(obs2ts[o]);
    }

    FREE(obs2ts);
}


struct s_drift **build_drift(json_t *json_drift, struct s_router **routers)
{
    struct s_drift **drift;
    drift = malloc(sizeof(struct s_drift *));
    if (drift==NULL) {
        char str[STR_BUFFSIZE];
        snprintf(str, STR_BUFFSIZE, "Allocation impossible in file :%s line : %d",__FILE__,__LINE__);
        print_err(str);
        exit(EXIT_FAILURE);
    }

    if(N_DRIFT){
        json_t *ind = fast_get_json_array(json_drift, "ind_par_Xdrift_applied");
        json_t *vol = fast_get_json_array(json_drift, "ind_volatility_Xdrift");

        int i;
        for(i=0; i<N_DRIFT; i++) {
            drift[i] = malloc(sizeof(struct s_drift));
            if (drift[i]==NULL) {
                char str[STR_BUFFSIZE];
                snprintf(str, STR_BUFFSIZE, "Allocation impossible in file :%s line : %d",__FILE__,__LINE__);
                print_err(str);
                exit(EXIT_FAILURE);
            }


            drift[i]->ind_par_Xdrift_applied = (int) fast_get_json_real_from_array(ind, i, "ind_par_Xdrift_applied");
            drift[i]->ind_volatility_Xdrift = (int) fast_get_json_real_from_array(vol, i, "ind_volatility_Xdrift");
            if(i == 0){
                drift[i]->offset = N_PAR_SV*N_CAC;
            }else{
                drift[i]->offset = drift[i-1]->offset + routers[drift[i-1]->ind_par_Xdrift_applied]->n_gp;
            }
        }
    }

    return drift;
}


void clean_drift(struct s_drift **drift)
{

    if (N_DRIFT) {
        int i;
        for(i=0; i<N_DRIFT; i++) {
            FREE(drift[i]);
        }
    }

    FREE(drift);
}

struct s_data *build_data(json_t *settings, json_t *theta, enum plom_implementations implementation, enum plom_noises_off noises_off, int is_bayesian)
{
    int n, ts, cac, k;
    int tmp_n_data_nonan, count_n_nan;

    json_t *json_data = fast_get_json_object(settings, "data");

    struct s_data *p_data;
    p_data = malloc(sizeof(struct s_data));
    if (p_data==NULL) {
        char str[STR_BUFFSIZE];
        snprintf(str, STR_BUFFSIZE, "Allocation impossible in file :%s line : %d",__FILE__,__LINE__);
        print_err(str);
        exit(EXIT_FAILURE);
    }

    p_data->implementation = implementation;
    p_data->noises_off = noises_off;

    //always present
    p_data->obs2ts = build_obs2ts(fast_get_json_array(json_data, "obs2ts"));
    p_data->routers = build_routers(settings, theta, is_bayesian);

    //ts names
    json_t *ts_name = fast_get_json_array(fast_get_json_object(settings, "orders"), "ts_id");
    p_data->ts_name = malloc(N_TS * sizeof(char *));
    if (p_data->ts_name == NULL) {
        char str[STR_BUFFSIZE];
        snprintf(str, STR_BUFFSIZE, "Allocation impossible in file :%s line : %d",__FILE__,__LINE__);
        print_err(str);
        exit(EXIT_FAILURE);
    }

    for (ts = 0; ts < N_TS; ts++ ) {
        const char *ts_key = fast_get_json_string_from_array(ts_name, ts, "ts_id");
        p_data->ts_name[ts] = init1c(strlen(ts_key) + 1);
        strcpy(p_data->ts_name[ts], ts_key);
    }

    //cac names
    json_t *cac_name = fast_get_json_array(fast_get_json_object(settings, "orders"), "cac_id");
    p_data->cac_name = malloc(N_CAC * sizeof(char *));
    if (p_data->cac_name == NULL) {
        char str[STR_BUFFSIZE];
        snprintf(str, STR_BUFFSIZE, "Allocation impossible in file :%s line : %d",__FILE__,__LINE__);
        print_err(str);
        exit(EXIT_FAILURE);
    }

    for (cac = 0; cac < N_CAC; cac++ ) {
        const char *cac_key = fast_get_json_string_from_array(cac_name, cac, "cac_id");
        p_data->cac_name[cac] = init1c(strlen(cac_key) + 1);
        strcpy(p_data->cac_name[cac], cac_key);
    }

    /* drift (if any) */
    p_data->drift = build_drift(fast_get_json_object(settings, "drift"), p_data->routers);

    /* iterators */
    p_data->p_it_all = build_iterator(settings, p_data->routers, p_data->drift, "all", p_data->noises_off);
    p_data->p_it_par_sv = build_iterator(settings, p_data->routers, p_data->drift, "par_sv", p_data->noises_off);
    p_data->p_it_all_no_drift = build_iterator(settings, p_data->routers, p_data->drift, "all_no_drift", p_data->noises_off);
    p_data->p_it_par_proc_par_obs_no_drift = build_iterator(settings, p_data->routers, p_data->drift, "par_proc_par_obs_no_drift", p_data->noises_off);
    p_data->p_it_par_sv_and_drift = build_iterator(settings, p_data->routers, p_data->drift, "par_sv_and_drift", p_data->noises_off);
    p_data->p_it_only_drift = build_iterator(settings, p_data->routers, p_data->drift, "only_drift", p_data->noises_off);
    p_data->p_it_noise = build_iterator(settings, p_data->routers, p_data->drift, "noise", p_data->noises_off);

    p_data->rep1 = fast_load_fill_json_2d(fast_get_json_array(json_data, "rep1"), "rep1");
    p_data->pop_size_t0 = fast_load_fill_json_1d(fast_get_json_array(json_data, "pop_size_t0"), "pop_size_t0");

    //the following is optional (for instance it is non needed for simulation models)
    if (N_DATA) {
        /*mandatory non fitted parameters and data*/
        p_data->data = fast_load_fill_json_2d(fast_get_json_array(json_data, "data"), "data");

        /*get N_DATA_NONAN and times */
        tmp_n_data_nonan=0;
        for(n=0; n<N_DATA; n++) {
            count_n_nan = 0;
            for(ts=0; ts<N_TS; ts++) //are all the lines composed only of NaN
                if (isnan(p_data->data[n][ts])) count_n_nan++;

            if(count_n_nan < N_TS) {
                tmp_n_data_nonan++;
                if (tmp_n_data_nonan == 1) {
                    p_data->times = init1u_set0(1);
                } else {
                    unsigned int *tmp;
                    tmp = realloc(p_data->times, tmp_n_data_nonan * sizeof(double ) );
                    if ( tmp == NULL ) {
                        print_err("Reallocation impossible");
                        FREE(p_data->times);
                        exit(EXIT_FAILURE);
                    }
                    else {
                        p_data->times = tmp;
                    }
                }
                p_data->times[tmp_n_data_nonan-1] = n+1;
            }
        }
        N_DATA_NONAN = tmp_n_data_nonan;

        /*data_ind*/
        struct s_data_ind **data_ind;
        data_ind = malloc(N_DATA_NONAN * sizeof(struct s_data_ind *));
        if (data_ind==NULL) {
            char str[STR_BUFFSIZE];
            snprintf(str, STR_BUFFSIZE, "Allocation impossible in file :%s line : %d",__FILE__,__LINE__);
            print_err(str);
            exit(EXIT_FAILURE);
        }

        for(n=0; n<N_DATA_NONAN; n++) {
            data_ind[n] = malloc(sizeof(struct s_data_ind));
            if (data_ind[n]==NULL) {
                char str[STR_BUFFSIZE];
                snprintf(str, STR_BUFFSIZE, "Allocation impossible in file :%s line : %d",__FILE__,__LINE__);
                print_err(str);
                exit(EXIT_FAILURE);
            }

            data_ind[n]->n_nonan=0;
            for(ts=0; ts<N_TS; ts++)
                if (!isnan(p_data->data[p_data->times[n]-1][ts]))
                    data_ind[n]->n_nonan += 1;

            if (data_ind[n]->n_nonan) {
                data_ind[n]->ind_nonan = init1u_set0(data_ind[n]->n_nonan);
                k=0;
                for(ts=0; ts<N_TS; ts++)
                    if (!isnan(p_data->data[p_data->times[n]-1][ts]))
                        data_ind[n]->ind_nonan[k++] = ts;
            } else {
                data_ind[n]->ind_nonan = NULL;
            }
        }

        p_data->data_ind = data_ind;
    }

    /*covariates (non fitted parameters) (forcing parameters a.k.a par_fixed_values)*/
    if (N_PAR_FIXED) {
        json_t *json_par_fixed = fast_get_json_array(fast_get_json_object(settings, "orders"), "par_fixed");
        json_t *json_par_fixed_values = fast_get_json_object(json_data, "par_fixed_values");

        p_data->par_fixed = malloc(N_PAR_FIXED * sizeof(double **));
        if (p_data->par_fixed==NULL) {
            char str[STR_BUFFSIZE];
            snprintf(str, STR_BUFFSIZE, "Allocation impossible in file :%s line : %d",__FILE__,__LINE__);
            print_err(str);
            exit(EXIT_FAILURE);
        }

        for(k=0; k< json_array_size(json_par_fixed); k++) {
            char par_fixed_name[255];
            json_t *tmp_str = json_array_get(json_par_fixed, k);
            json_t *json_my_par_fixed_values;

            if (!json_is_string(tmp_str)) {
                char str[STR_BUFFSIZE];
                snprintf(str, STR_BUFFSIZE, "error: par_fixed[%d] is not a string\n", k);
                print_err(str);
                exit(EXIT_FAILURE);
            }
            strcpy(par_fixed_name, json_string_value(tmp_str));
            json_my_par_fixed_values = fast_get_json_array(json_par_fixed_values, par_fixed_name);
            p_data->par_fixed[k] = fast_load_fill_json_2d(json_my_par_fixed_values, par_fixed_name);
        }
    }

    if (IS_SCHOOL_TERMS) {
        json_t *json_school = fast_get_json_array(json_data, "school_terms");
        int cac;

        p_data->n_terms = init1u_set0(N_CAC);

        for(cac=0 ; cac< N_CAC ; cac++) {
            json_t *json_school_cac = json_array_get(json_school, cac);
            if (!json_is_array(json_school_cac)) {
                char str[STR_BUFFSIZE];
                snprintf(str, STR_BUFFSIZE, "error: json_school[%d] is not an array\n", cac);
                print_err(str);
                exit(EXIT_FAILURE);
            }
            p_data->n_terms[cac] = json_array_size(json_school_cac) / 2; //we divide by 2 because begining and end of terms are on the same line
        }

        //to fill p_data->school_terms we use a temporary array:
        double **tmp_sch = fast_load_fill_json_2d(json_school, "school_terms"); //the temporary one
        p_data->school_terms = init3d_varp1_set0(N_CAC, p_data->n_terms, 2); //the true one
        for(cac=0 ; cac< N_CAC ; cac++) {
            for(k=0 ; k< p_data->n_terms[cac] ; k++) {
                p_data->school_terms[cac][k][0] = tmp_sch[cac][k*2];
                p_data->school_terms[cac][k][1] = tmp_sch[cac][k*2+1];
            }
        }
        clean2d(tmp_sch, N_CAC); //no longer needed

        /* compute proportion of year taken up by school !!!shool terms are in years!!! : */
        p_data->prop_school = init1d_set0(N_CAC);
        for(cac=0 ; cac< N_CAC ; cac++) {
            for(k=0 ; k<p_data->n_terms[cac]; k++) {
                p_data->prop_school[cac]+= (p_data->school_terms[cac][k][1]-p_data->school_terms[cac][k][0]);
            }
            p_data->prop_school[cac] = (1.0 - p_data->prop_school[cac]);
        }
    }

    if (N_AC > 1) {
        p_data->waifw = fast_load_fill_json_2d(fast_get_json_array(json_data, "waifw"), "waifw");
    }

    //p_data->mat_d = init2d_set0(N_C, N_C);
    //load2d(p_data->mat_d, "data/mat_d.data", N_C);

    return p_data;

}


void clean_data(struct s_data *p_data)
{
    int n;

    clean2c(p_data->ts_name, N_TS);
    clean2c(p_data->cac_name, N_CAC);

    clean_iterator(p_data->p_it_all);
    clean_iterator(p_data->p_it_par_sv);
    clean_iterator(p_data->p_it_all_no_drift);
    clean_iterator(p_data->p_it_par_proc_par_obs_no_drift);
    clean_iterator(p_data->p_it_par_sv_and_drift);
    clean_iterator(p_data->p_it_only_drift);
    clean_iterator(p_data->p_it_noise);

    clean_obs2ts(p_data->obs2ts);
    clean_routers(p_data->routers);
    clean_drift(p_data->drift);

    clean2d(p_data->rep1, (N_DATA == 0) ? 1 : N_DATA); //for simulation models
    FREE(p_data->pop_size_t0);

    if (N_DATA) {
        clean2d(p_data->data, N_DATA);
        FREE(p_data->times);

        /*data_ind*/
        for(n=0; n<N_DATA_NONAN; n++) {
            FREE(p_data->data_ind[n]->ind_nonan);
            FREE(p_data->data_ind[n]);
        }
        FREE(p_data->data_ind);
    }

    /*extra non fitted parameters*/
    if (N_PAR_FIXED) {
        clean3d(p_data->par_fixed, N_PAR_FIXED, N_DATA_PAR_FIXED);
    }

    /* school terms */
    if (IS_SCHOOL_TERMS) {
        clean3d_var(p_data->school_terms, N_CAC, p_data->n_terms);
        FREE(p_data->prop_school);
        FREE(p_data->n_terms);
    }

    if (N_AC > 1) {
        clean2d(p_data->waifw, N_AC);
    }

    //clean2d(p_data->mat_d, N_C);
    FREE(p_data);
}



struct s_calc **build_calc(int *n_threads, int general_id, double eps_abs, double eps_rel, int J, int dim_ode, int (*func_step_ode) (double, const double *, double *, void *), struct s_data *p_data)


{
    char str[STR_BUFFSIZE];
    int nt;


    *n_threads = sanitize_n_threads(*n_threads, J);
    omp_set_num_threads(*n_threads);

    struct s_calc **calc;
    calc=malloc(*n_threads*sizeof (struct s_calc *));
    if (calc==NULL) {
        snprintf(str, STR_BUFFSIZE, "Allocation impossible in file :%s line : %d",__FILE__,__LINE__);
        print_err(str);
        exit(EXIT_FAILURE);
    }

    for(nt=0; nt< *n_threads; nt++) {
        calc[nt]=malloc(sizeof (struct s_calc));
        if (calc[nt]==NULL) {
            snprintf(str, STR_BUFFSIZE, "Allocation impossible in file :%s line : %d",__FILE__,__LINE__);
            print_err(str);
            exit(EXIT_FAILURE);
        }
    }

    unsigned long int seed; /*SEED*/
#if FLAG_JSON
    seed = (unsigned) time(NULL);
#else
    seed=2;
#endif
    seed += general_id; /*we ensure uniqueness of seed in case of parrallel runs*/

    snprintf(str, STR_BUFFSIZE, "seed=%ld", seed);
    print_log(str);

    /* we create as many rng as parallel threads *but* note that for the operations not prarallelized, we always use cacl[0].randgsl */
    for (nt=0; nt< *n_threads; nt++) {
        calc[nt] = build_p_calc(*n_threads, nt, seed, eps_abs, eps_rel, dim_ode, func_step_ode, p_data);
    }

    return calc;
}


/**
 * eps_abs, eps_rel: absolute and relative error 
 */

struct s_calc *build_p_calc(int n_threads, int thread_id, int seed, double eps_abs, double eps_rel, int dim_ode, int (*func_step_ode) (double, const double *, double *, void *), struct s_data *p_data)
{
    struct s_calc *p_calc = malloc(sizeof(struct s_calc));
    if (p_calc==NULL) {
	char str[STR_BUFFSIZE];
        snprintf(str, STR_BUFFSIZE, "Allocation impossible in file :%s line : %d",__FILE__,__LINE__);
        print_err(str);
        exit(EXIT_FAILURE);
    }


    p_calc->n_threads = n_threads;

    p_calc->current_n = 0;
    p_calc->current_nn = 0;

    /* ref */
    p_calc->p_data = p_data;

    /* thread_id */
    p_calc->thread_id = thread_id;


    /*random numbers...*/
    /*random number generator and parallel MC simulations*/
    /*
      idea using one different seed per thread but is it realy uncorelated ???
      Should I go through the trouble of changing from GSL to SPRNG????
      answer:
      I would recommend using ranlxd.  The seeds should give 2^31
      effectively independent streams of length 10^171.  A discussion of the
      seeding procedure can be found in the file notes.ps at
      http://www.briangough.ukfsn.org/ranlux_2.2/
      --
      Brian Gough
    */

    const gsl_rng_type *Type; /*type de generateur aleatoire*/
    if (n_threads == 1){ //we don't need a rng supporting parallel computing, we use mt19937 that is way faster than ranlxs0 (1754 k ints/sec vs 565 k ints/sec)
        Type = gsl_rng_mt19937; /*MT19937 generator of Makoto Matsumoto and Takuji Nishimura*/
    } else {
        Type = gsl_rng_ranlxs0; //gsl_rng_ranlxs2 is better than gsl_rng_ranlxs0 but 2 times slower
    }

    /*we create as many rng as parallel threads *but* note that for the operations not prarallelized, we always use cacl[0].randgsl*/

    /* rng */
    p_calc->randgsl=gsl_rng_alloc(Type);
    gsl_rng_set(p_calc->randgsl, seed + thread_id);

    if (p_data->implementation == PLOM_ODE){

        p_calc->T = gsl_odeiv2_step_rkf45;
        p_calc->control = gsl_odeiv2_control_y_new(eps_abs, eps_rel);
        p_calc->step = gsl_odeiv2_step_alloc(p_calc->T, dim_ode);
        p_calc->evolve = gsl_odeiv2_evolve_alloc(dim_ode);
        (p_calc->sys).function = func_step_ode;
        (p_calc->sys).jacobian = jac;
        (p_calc->sys).dimension=(dim_ode);
        (p_calc->sys).params= p_calc;

        p_calc->yerr = init1d_set0(dim_ode);
    } else if (p_data->implementation == PLOM_SDE){
	p_calc->y_pred = init1d_set0(dim_ode);
    } else if (p_data->implementation == PLOM_PSR){
        build_psr(p_calc);
    }

    //multi-threaded sorting
    p_calc->to_be_sorted = init1d_set0(J);
    p_calc->index_sorted = init1st_set0(J);

    return p_calc;
}

void clean_p_calc(struct s_calc *p_calc, struct s_data *p_data)
{
    gsl_rng_free(p_calc->randgsl);

    if (p_data->implementation == PLOM_ODE){

        gsl_odeiv2_step_free(p_calc->step);
        gsl_odeiv2_evolve_free(p_calc->evolve);
        gsl_odeiv2_control_free(p_calc->control);

        FREE(p_calc->yerr);

    } else if (p_data->implementation == PLOM_SDE){
	FREE(p_calc->y_pred);

    } else if (p_data->implementation == PLOM_PSR){

        clean2d(p_calc->prob, N_PAR_SV+2); //+2 for U and DU of the universes
        clean3u(p_calc->inc, N_PAR_SV+2, N_CAC); //+2 for U and DU of the universes
    }

    FREE(p_calc->to_be_sorted);
    FREE(p_calc->index_sorted);

    FREE(p_calc);
}


void clean_calc(struct s_calc **calc, struct s_data *p_data)
{
    int nt;
    int n_threads = calc[0]->n_threads;

    /*clean calc*/
    for(nt=0; nt<n_threads; nt++) {
        clean_p_calc(calc[nt], p_data);
    }

    FREE(calc);
}



/**
 * dt: integration time step. If <=0.0, default to 0.25/365.0 * ONE_YEAR_IN_DATA_UNIT
 */

struct s_X *build_X(int size_proj, int size_obs, struct s_data *p_data, double dt)
{
    struct s_X *p_X;
    p_X = malloc(sizeof(struct s_X));
    if (p_X==NULL) {
        char str[STR_BUFFSIZE];
        snprintf(str, STR_BUFFSIZE, "Allocation impossible in file :%s line : %d",__FILE__,__LINE__);
        print_err(str);
        exit(EXIT_FAILURE);
    }

    //integration time step
    if (dt <= 0.0){ //default to 0.25/365.0 * ONE_YEAR_IN_DATA_UNIT
	dt = 0.25/365.0 * ONE_YEAR_IN_DATA_UNIT;
    }
    //IMPORTANT: we ensure an integer multiple of dt in between 2 data points    
    p_X->dt = 1.0/ ((double) round(1.0/dt));

    p_X->proj = init1d_set0(size_proj);
    p_X->obs = init1d_set0(size_obs);

    return p_X;
}

void clean_X(struct s_X *p_X)
{
    FREE(p_X->proj);
    FREE(p_X->obs);

    FREE(p_X);
}


struct s_X **build_J_p_X(int size_proj, int size_obs, struct s_data *p_data, double dt)
{
    int j;

    struct s_X **J_p_X;
    J_p_X = malloc(J *sizeof(struct s_X *));
    if (J_p_X==NULL) {
        char str[STR_BUFFSIZE];
        snprintf(str, STR_BUFFSIZE, "Allocation impossible in file :%s line : %d",__FILE__,__LINE__);
        print_err(str);
        exit(EXIT_FAILURE);
    }

    for(j=0; j<J; j++){
        J_p_X[j] = build_X(size_proj, size_obs, p_data, dt);
    }

    return J_p_X;
}

void clean_J_p_X(struct s_X **J_p_X)
{
    int j;
    for(j=0; j<J; j++) {
        clean_X(J_p_X[j]);
    }

    FREE(J_p_X);
}


struct s_X ***build_D_J_p_X(int size_proj, int size_obs, struct s_data *p_data, double dt)
{
    int n;

    struct s_X ***D_J_p_X;
    D_J_p_X = malloc((N_DATA+1) *sizeof(struct s_X **));
    if (D_J_p_X==NULL){
        char str[STR_BUFFSIZE];
        snprintf(str, STR_BUFFSIZE, "Allocation impossible in file :%s line : %d",__FILE__,__LINE__);
        print_err(str);
        exit(EXIT_FAILURE);
    }

    for(n=0; n<(N_DATA+1); n++) {
        D_J_p_X[n] = build_J_p_X(size_proj, size_obs, p_data, dt);
    }

    return D_J_p_X;
}


void clean_D_J_p_X(struct s_X ***D_J_p_X)
{
    int n;
    for(n=0; n<(N_DATA+1); n++) {
        clean_J_p_X(D_J_p_X[n]);
    }

    FREE(D_J_p_X);
}


struct s_hat *build_hat(struct s_data *p_data)
{
    struct s_hat *p_hat = malloc(sizeof (struct s_hat));
    if (p_hat==NULL) {
        char str[STR_BUFFSIZE];
        snprintf(str, STR_BUFFSIZE, "Allocation impossible in file :%s line : %d",__FILE__,__LINE__);
        print_err(str);
        exit(EXIT_FAILURE);
    }

    p_hat->state = init1d_set0(N_PAR_SV*N_CAC);
    p_hat->state_95 = init2d_set0(N_PAR_SV*N_CAC, 2);

    p_hat->obs = init1d_set0(N_TS);
    p_hat->obs_95 = init2d_set0(N_TS, 2);

    if(p_data->p_it_only_drift->nbtot) {
        p_hat->drift = init1d_set0(p_data->p_it_only_drift->nbtot);
        p_hat->drift_95 = init2d_set0(p_data->p_it_only_drift->nbtot, 2);
    }

    return p_hat;
}


void clean_hat(struct s_hat *p_hat, struct s_data *p_data)
{
    FREE(p_hat->state);
    clean2d(p_hat->state_95, N_PAR_SV*N_CAC);

    FREE(p_hat->obs);
    clean2d(p_hat->obs_95, N_TS);

    if(p_data->p_it_only_drift->nbtot) {
        FREE(p_hat->drift);
        clean2d(p_hat->drift_95, p_data->p_it_only_drift->nbtot);
    }

    FREE(p_hat);
}


struct s_hat **build_D_p_hat(struct s_data *p_data)
{
    int n;

    struct s_hat **D_p_hat = malloc(N_DATA * sizeof (struct s_hat *));
    if (D_p_hat==NULL) {
        char str[STR_BUFFSIZE];
        snprintf(str, STR_BUFFSIZE, "Allocation impossible in file :%s line : %d",__FILE__,__LINE__);
        print_err(str);
        exit(EXIT_FAILURE);
    }

    for(n=0; n<N_DATA; n++) {
        D_p_hat[n] = build_hat(p_data);
    }

    return D_p_hat;
}


void clean_D_p_hat(struct s_hat **D_p_hat, struct s_data *p_data)
{
    int n;

    for(n=0; n<N_DATA; n++) {
        clean_hat(D_p_hat[n], p_data);
    }

    FREE(D_p_hat);
}




struct s_likelihood *build_likelihood(void)
{
    struct s_likelihood *p_like = malloc(sizeof (struct s_likelihood));
    if (p_like==NULL) {
        char str[STR_BUFFSIZE];
        snprintf(str, STR_BUFFSIZE, "Allocation impossible in file :%s line : %d",__FILE__,__LINE__);
        print_err(str);
        exit(EXIT_FAILURE);
    }

    p_like->ess_n = 0.0;
    p_like->Llike_best_n = 0.0;
    p_like->Llike_best = 0.0;
    p_like->weights = init1d_set0(J);

    p_like->select = init2u_set0(N_DATA_NONAN, J);

    p_like->n_all_fail = 0;

    /* for bayesian methods */
    p_like->Llike_prev = 0.0;
    p_like->Llike_new = 0.0;

    return p_like;
}

void clean_likelihood(struct s_likelihood *p_like)
{
    FREE(p_like->weights);
    clean2u(p_like->select, N_DATA_NONAN);

    FREE(p_like);
}

struct s_best *build_best(struct s_data *p_data, json_t *theta, int update_covariance)
{
    int i,j,k;
    struct s_router **routers = p_data->routers;

    struct s_best *p_best;
    p_best = malloc(sizeof(struct s_best));
    if (p_best==NULL) {
        char str[STR_BUFFSIZE];
        snprintf(str, STR_BUFFSIZE, "Allocation impossible in file :%s line : %d",__FILE__,__LINE__);
        print_err(str);
        exit(EXIT_FAILURE);
    }

    p_best->length = p_data->p_it_all->nbtot;

    p_best->mean = gsl_vector_calloc(p_best->length);
    p_best->proposed = gsl_vector_calloc(p_best->length);

    p_best->var = gsl_matrix_calloc(p_best->length, p_best->length);

    p_best->par_prior = init2d_set0(p_best->length, 2);
    p_best->prior = malloc(p_best->length * sizeof(double (*) (double)) );

    if (p_best->prior==NULL) {
        char str[STR_BUFFSIZE];
        snprintf(str, STR_BUFFSIZE, "Allocation impossible in file :%s line : %d",__FILE__,__LINE__);
        print_err(str);
        exit(EXIT_FAILURE);
    }

    p_best->to_be_estimated = init1u_set0(p_best->length);

    p_best->is_follower = init1u_set0(p_best->length);
    p_best->n_follow = 0;
    for(i=0; i<p_data->p_it_all->length; i++) {
        const char *par_key = routers[i]->name;
        json_t *par = fast_get_json_object(fast_get_json_object(theta, "value"), par_key);

        if(json_object_get(par, "follow")) {
            const char *par_follow_key = fast_get_json_string_from_object(par, "follow");

            p_best->n_follow++;
            if (p_best->n_follow == 1) {
                p_best->follower = init1u_set0(1);
                p_best->follow = init1u_set0(1);
            } else {
                unsigned int *tmp, *tmp2;
                tmp = realloc(p_best->follower, p_best->n_follow * sizeof(unsigned int) );
                tmp2 = realloc(p_best->follow, p_best->n_follow * sizeof(unsigned int) );
                if ( tmp == NULL || tmp2 == NULL ) {
                    print_err("Reallocation impossible");
                    FREE(p_best->follower); FREE(p_best->follow);
                    exit(EXIT_FAILURE);
                } else {
                    p_best->follower = tmp;
                    p_best->follow = tmp2;
                }
            }
            p_best->follower[p_best->n_follow-1] = i;

            for(j=0; j< p_data->p_it_all->length; j++){
                if(strcmp(par_follow_key, routers[j]->name) == 0){
                    p_best->follow[p_best->n_follow-1] = j;
                    break;
                }
            }
            if(j==p_data->p_it_all->length){
                char str[STR_BUFFSIZE];
                snprintf(str, STR_BUFFSIZE, "parameters to follow could not be found (%s -> %s)", routers[i]->name, par_follow_key);
                print_err(str);
                exit(EXIT_FAILURE);
            }

            for (k=0; k<routers[i]->n_gp; k++) {
                p_best->is_follower[p_data->p_it_all->offset[i] + k] = 1;
            }

        } else {

            for (k=0; k<routers[i]->n_gp; k++) {
                p_best->is_follower[p_data->p_it_all->offset[i] + k] = 0;
            }
        }
    }

    load_best(p_best, p_data, theta, 1, update_covariance);

    gsl_vector_memcpy(p_best->proposed, p_best->mean);

    update_to_be_estimated(p_best);

    p_best->mean_sampling = init1d_set0(p_best->length);
    p_best->var_sampling = gsl_matrix_calloc(p_best->length, p_best->length);

    return p_best;
}


void clean_best(struct s_best *p_best)
{
    gsl_vector_free(p_best->mean);
    gsl_vector_free(p_best->proposed);
    gsl_matrix_free(p_best->var);

    clean2d(p_best->par_prior, p_best->length);
    FREE(p_best->prior);

    FREE(p_best->mean_sampling);
    gsl_matrix_free(p_best->var_sampling);

    FREE(p_best->to_be_estimated);

    FREE(p_best->is_follower);
    if(p_best->n_follow){
        FREE(p_best->follower);
        FREE(p_best->follow);
    }

    FREE(p_best);
}
