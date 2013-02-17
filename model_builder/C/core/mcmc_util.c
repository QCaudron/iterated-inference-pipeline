#include "plom.h"
#include "mcmc_util.h"


struct s_mcmc_calc_data *build_mcmc_calc_data(struct s_best *p_best, const double a, const int m_switch, const int m_epsilon, const double epsilon_max, const int is_smoothed_tunning, const double alpha)
{
    /* which parameters have to be estimated (parameters with jump_size > 0.0) */
    char str[STR_BUFFSIZE];
    int k;

    struct s_mcmc_calc_data *p;
    p = malloc(sizeof(struct s_mcmc_calc_data));
    if(p==NULL) {
	snprintf(str, STR_BUFFSIZE, "Allocation impossible in file :%s line : %d",__FILE__,__LINE__);
        print_err(str);
        exit(EXIT_FAILURE);
    }

    p->n_acceptance_rates = p_best->length;
    p->acceptance_rates = init1d_set0(p->n_acceptance_rates);
    p->smoothed_acceptance_rates = init1d_set0(p->n_acceptance_rates);

    //start with 1.0 (for non fitted parameters it will stay at 1.0 (in a way, non fitted parameters are always accepted)
    for(k=0; k< (p->n_acceptance_rates); k++) {
        p->acceptance_rates[k] = 1.0;
        p->smoothed_acceptance_rates[k] = 1.0;
    }

    p->global_acceptance_rate = 1.0;
    p->smoothed_global_acceptance_rate = 1.0;

    p->has_cycled = 1;
    p->m_full_iteration = 0;
    p->cycle_id = p_best->n_to_be_estimated -1;

    p->epsilon = 1.0;
    p->epsilon_max = epsilon_max;
    p->a = a;

    p->alpha = alpha;
    p->is_smoothed_tunning = is_smoothed_tunning;

    // iteration to swith between initial and empirical covariances
    p->m_switch = m_switch;

    int min_switch = 5*p_best->n_to_be_estimated*p_best->n_to_be_estimated;
    if (m_switch < 0) {
        p->m_switch = min_switch;
    } else if (p->m_switch < min_switch) {
	snprintf(str, STR_BUFFSIZE, "attention: covariance switching iteration (%i) is smaller than proposed one (%i)\n", m_switch, min_switch);
        print_warning(str);
    }

    p->m_epsilon = m_epsilon;

    return (p);
}

void clean_mcmc_calc_data(struct s_mcmc_calc_data *p)
{
    FREE(p->acceptance_rates);
    FREE(p->smoothed_acceptance_rates);
    FREE(p);
}


/**
 * IN FULL UPDATE CASE:
 *
 * return the empirical covariance matrix or the initial one
 * (depending on the iteration value and options) and the evaluated
 * tuning factor sd_fac
 * 
 *
 * IN SEQUENTIAL UPDATE CASE:
 *
 * Generate a random value for a component "k" of proposed chosen
 * randomly.  Every number of parameters with jump_size > 0.0
 * (n_to_be_estimated) we randomly shuffle the index of the n_to_be
 * estimated parameter index. In between these shuffling event, we
 * cycle throufh the shuffled indexes. Traces are printed every
 * n_to_be_estimated iterations. This should mimate block update of
 * the n_to_be_estimated component of theta while increasing the
 * acceptance ratio
 *
 * return the initial covariance matrix
 */
gsl_matrix * get_var_and_sd_fac(double *sd_fac, struct s_best *p_best, struct s_mcmc_calc_data *p, struct s_calc *p_calc, int m)
{
    if (OPTION_FULL_UPDATE) {
        //////////////////////
        // FULL UPDATE CASE //
        //////////////////////

        // evaluate epsilon(m) = epsilon(m-1) * exp(a^(m-1) * (acceptance_rate(m-1) - 0.234))

        if ( (m > p->m_epsilon) && (m*p->global_acceptance_rate < p->m_switch) ) {

	    double ar = (p->is_smoothed_tunning) ? p->smoothed_global_acceptance_rate : p->global_acceptance_rate;
            p->epsilon *=  exp(pow(p->a, (double)(m-1)) * (ar - 0.234));

        } else {
	    // after switching epsilon is set back to 1
            p->epsilon = 1.0;
        }

	p->epsilon = GSL_MIN(p->epsilon, p->epsilon_max);
	
#if FLAG_VERBOSE
	char str[STR_BUFFSIZE];
	snprintf(str, STR_BUFFSIZE, "epsilon = %g", p->epsilon);
	print_log(str);
#endif

        // evaluate tuning factor sd_fac = epsilon * 2.38/sqrt(n_to_be_estimated)
        *sd_fac = p->epsilon * 2.38/sqrt(p_best->n_to_be_estimated);
	
        if( (m * p->global_acceptance_rate) >= p->m_switch) {
	    return p_best->var_sampling;	    
	} else {
	    return p_best->var;
	}

    } else {

        ////////////////////////////
        // SEQUENTIAL UPDATE CASE //
        ////////////////////////////

        *sd_fac = 1.0;

        if(p_best->n_to_be_estimated > 0) { //due to the webApp all jump size can be 0.0...
            if(p->has_cycled) {
                gsl_ran_shuffle(p_calc->randgsl, p_best->to_be_estimated, p_best->n_to_be_estimated, sizeof (unsigned int));
            }
        }

	return p_best->var;
    }
}




/**
 * Acceptance rate(s): we use an average filter for the local version.
 * For the webApp, we use a low pass filter (a.k.a exponential
 * smoothing or exponential moving average) to reflect live tunning of
 * the walk rates
 */
void compute_acceptance_rates(struct s_best *p_best, struct s_mcmc_calc_data *p, double is_accepted, int m)
{

    if (!OPTION_FULL_UPDATE) {
        if(p_best->n_to_be_estimated > 0) { //due to the webApp all jump size can be 0.0...
            int mm = p_best->to_be_estimated[ p->cycle_id ];
            p->smoothed_acceptance_rates[mm] = (1.0 - p->alpha) * p->smoothed_acceptance_rates[mm] + p->alpha * is_accepted;
            p->acceptance_rates[mm] += ((is_accepted - p->acceptance_rates[mm])/ ((double) p->m_full_iteration));
        }
    }

    p->smoothed_global_acceptance_rate = (1.0 - p->alpha) * p->smoothed_global_acceptance_rate + p->alpha * is_accepted;
    p->global_acceptance_rate += ( (is_accepted - p->global_acceptance_rate) / ((double) m) ); //note that we divide by m and not p->m_full_iteration

}



void increment_iteration_counters(struct s_mcmc_calc_data *p, struct s_best *p_best, const int is_full_update)
{

    if(is_full_update) {
        p->has_cycled = 1;
        p->m_full_iteration ++;
        p->cycle_id = p_best->n_to_be_estimated;
    } else {

        if (p->cycle_id >= (p_best->n_to_be_estimated -1)) { // >= instead of  == because due to the webApp all jump size can be 0.0...
            p->has_cycled = 1;
            p->m_full_iteration += 1;
            p->cycle_id = 0;
        } else {
            p->has_cycled = 0;
            p->cycle_id += 1;
        }

    }
}




void ran_proposal_sequential(gsl_vector *proposed, struct s_best *p_best, gsl_matrix *var, double sd_fac, struct s_calc *p_calc)
{
    int k;

    struct s_mcmc_calc_data *p =  (struct s_mcmc_calc_data *) p_calc->method_specific_shared_data;

    if (p_best->n_to_be_estimated > 0) { //due to the webApp all jump size can be 0.0...
        k = p_best->to_be_estimated[ p->cycle_id ];
        gsl_vector_set(proposed, k, gsl_vector_get(p_best->mean, k) + gsl_ran_gaussian(p_calc->randgsl, sd_fac*sqrt(gsl_matrix_get(var, k, k))));
    }

}



void header_acceptance_rates(FILE *p_file, struct s_data *p_data)
{
    int i, g;
    struct s_router **routers = p_data->routers;

    fprintf(p_file, "index,epsilon,global_ar");

    if (!OPTION_FULL_UPDATE) {
	for(i=0; i<p_data->p_it_all->length; i++) {
	    const char *name = routers[i]->name;
	    for(g=0; g<p_data->routers[i]->n_gp; g++) {
		const char *group = routers[i]->group_name[g];
		fprintf(p_file, ",%s:%s", name, group);
	    }
	}
    }

    fprintf(p_file, "\n");
}


/**
 * for the webApp, we print the smoothed values
 */
void print_acceptance_rates(FILE *p_file, struct s_mcmc_calc_data *p, int m_full_iteration)
{
    int k;

#if FLAG_JSON
    json_t *root;
    json_t *j_print = json_array();
#endif

#if FLAG_JSON
    json_array_append_new(j_print, json_integer(m_full_iteration));
    json_array_append_new(j_print, json_real(p->epsilon));
    json_array_append_new(j_print, json_real(p->smoothed_global_acceptance_rate));
#else
    fprintf(p_file, "%d,%g,%g", m_full_iteration, p->epsilon, p->global_acceptance_rate);
#endif

    /* parameter specific acceptance rates */
    if (!OPTION_FULL_UPDATE) {
        for(k=0; k< (p->n_acceptance_rates); k++) {
#if FLAG_JSON
            json_array_append_new(j_print, json_real(p->smoothed_acceptance_rates[k]));
#else
	    fprintf(p_file, ",%g", p->acceptance_rates[k]);
#endif
        }
    }

#if FLAG_JSON
    root = json_pack("{s,s,s,o}", "flag", "pmcmc", "msg", j_print);
    json_dumpf(root, stdout, JSON_COMPACT); printf("\n");
    fflush(stdout);
    json_decref(root);
#else
    fprintf(p_file, "\n");
#endif

}


/**
 * print empirical matrix of variance covariance
 */
void print_covariance(FILE *p_file_cov, gsl_matrix *covariance)
{

    int row, col;
    double x;
#if FLAG_JSON
    json_t *root;
    json_t *json_print = json_array();
    json_t *json_print_n;
#endif


    for(row=0; row<covariance->size1; row++) {

#if FLAG_JSON
        json_print_n = json_array();
#endif
        for(col=0; col<covariance->size2; col++) {
            x = gsl_matrix_get(covariance, row, col);
#if FLAG_JSON
            json_array_append_new(json_print_n, json_real(x));
#else
            fprintf(p_file_cov,"%g%s", x, (col < (covariance->size2-1)) ? ",": "");
#endif
        }

#if FLAG_JSON
        json_array_append_new(json_print, json_print_n);
#else
        fprintf(p_file_cov,"\n");
#endif
    }


#if FLAG_JSON
    root = json_pack("{s,s,s,o}", "flag", "cov", "msg", json_print);
    json_dumpf(root, stdout, JSON_COMPACT); printf("\n");
    fflush(stdout);
    json_decref(root);
#endif

}
