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

#include "pmcmc.h"

/**
 * propose new p_best->proposed from p_best->mean
 *
 * IN FULL UPDATE CASE:
 *
 * new p_best->proposed is sampled by a MVN generator using the covariance matrix
 * 2.38²*epsilon²/n_to_be_estimated * var where:
 * - epsilon is dynamicaly tuned following the iterative law:
 *   epsilon(m) = epsilon(m-1) * exp(a^m * (acceptance rate - 0.234))
 * - var is the initial covariance if the number of accepted particles
 *   is lower than SWITCH, and the empirical one otherwise
 *
 * For the calculation, the collapsed Cholesky decomposition of var is used,
 * so the tuning factor is 2.38*epsilon/sqrt(n_to_be_estimated).
 *
 * return the empirical covariance matrix or the initial one depending on if(m*global_acceptance_rate >= m_switch)
 * 
 *
 * IN SEQUENTIAL UPDATE CASE:
 *
 * One component of p_best->proposed is sampled at each iteration
 * following a gaussian law. The covariances are the diagonal terms of
 * the initial covariance matrix (loaded from a covariace.output file
 * or filled with jump sizes. The components sampling order is shuffled
 * each time all components have been sampled (p_pMCMC_specific_data->has_cycled).
 *
 * return the initial covariance matrix
 */
gsl_matrix * propose_new_theta_and_load_X0(double *sd_fac, struct s_best *p_best, struct s_X *p_X, struct s_par *p_par, struct s_data *p_data, struct s_mcmc_calc_data *p, struct s_calc *p_calc, int m)
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
	    //p_best->proposed ~ MVN(p_best->mean, p_best->var_sampling)
	    propose_safe_theta_and_load_X0(p_best->proposed, p_best, p_best->var_sampling, *sd_fac, p_par, p_X, p_data, p_calc, plom_rmvnorm);
	
	    return p_best->var_sampling;
	
	} else {
	    //p_best->proposed ~ MVN(p_best->mean, p_best->var)
	    propose_safe_theta_and_load_X0(p_best->proposed, p_best, p_best->var, *sd_fac, p_par, p_X, p_data, p_calc, plom_rmvnorm);

	    return p_best->var;
	}

    } else {

        ////////////////////////////
        // SEQUENTIAL UPDATE CASE //
        ////////////////////////////

        /*generate a random value for a component "k" of proposed chosen
          randomly.  every number of parameters with jump_size > 0.0
          (n_to_be_estimated) we randomly shuffle the index of the n_to_be
          estimated parameter index. In between these shuffling event, we
          cycle throufh the shuffled indexes. Traces are printed every
          n_to_be_estimated iterations. This should mimate block update of
          the n_to_be_estimated component of theta while increasing the
          acceptance ratio
        */
        *sd_fac = 1.0;

        if(p_best->n_to_be_estimated > 0) { //due to the webApp all jump size can be 0.0...
            if(p->has_cycled) {
                gsl_ran_shuffle(p_calc->randgsl, p_best->to_be_estimated, p_best->n_to_be_estimated, sizeof (unsigned int));
            }
        }
        propose_safe_theta_and_load_X0(p_best->proposed, p_best, p_best->var, 1.0, p_par, p_X, p_data, p_calc, ran_proposal_sequential);

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


/**
 * recursive expression for the average so that it can be used in
 * real time by print_hat (zmq and co...)
 */
void compute_best_traj(struct s_hat **D_p_hat_best, struct s_hat **D_p_hat_prev, struct s_hat **D_p_hat_new, struct s_data *p_data, double alpha, double m)
{
    /* recursive expression for the average so that it can be used in
       real time by print_hat (zmq and co...) */

    int n, i, ts;

    for(n=0; n<N_DATA; n++) {
        //sv
        for(i=0 ; i<N_PAR_SV*N_CAC ; i++) {
            D_p_hat_best[n]->state[i] = ((m-1.0)/m)*D_p_hat_best[n]->state[i] + (1.0/m)*(alpha*D_p_hat_new[n]->state[i] + (1.0-alpha)*D_p_hat_prev[n]->state[i]);
            D_p_hat_best[n]->state_95[i][0] = ((m-1.0)/m)*D_p_hat_best[n]->state_95[i][0] + (1.0/m)*(alpha*D_p_hat_new[n]->state_95[i][0] + (1.0-alpha)*D_p_hat_prev[n]->state_95[i][0]);
            D_p_hat_best[n]->state_95[i][1] = ((m-1.0)/m)*D_p_hat_best[n]->state_95[i][1] + (1.0/m)*(alpha*D_p_hat_new[n]->state_95[i][1] + (1.0-alpha)*D_p_hat_prev[n]->state_95[i][1]);
        }

        //ts
        for(ts=0; ts< N_TS; ts++) {
            D_p_hat_best[n]->obs[ts] = ((m-1.0)/m)*D_p_hat_best[n]->obs[ts] + (1.0/m)*(alpha*D_p_hat_new[n]->obs[ts] + (1.0-alpha)*D_p_hat_prev[n]->obs[ts]);
            D_p_hat_best[n]->obs_95[ts][0] = ((m-1.0)/m)*D_p_hat_best[n]->obs_95[ts][0] + (1.0/m)*(alpha*D_p_hat_new[n]->obs_95[ts][0] + (1.0-alpha)*D_p_hat_prev[n]->obs_95[ts][0]);
            D_p_hat_best[n]->obs_95[ts][1] = ((m-1.0)/m)*D_p_hat_best[n]->obs_95[ts][1] + (1.0/m)*(alpha*D_p_hat_new[n]->obs_95[ts][1] + (1.0-alpha)*D_p_hat_prev[n]->obs_95[ts][1]);
        }

        //drift
        for(i=0; i< p_data->p_it_only_drift->nbtot; i++) {
            D_p_hat_best[n]->drift[i] = ((m-1.0)/m)*D_p_hat_best[n]->drift[i] + (1.0/m)*(alpha*D_p_hat_new[n]->drift[i] + (1.0-alpha)*D_p_hat_prev[n]->drift[i]);
            D_p_hat_best[n]->drift_95[i][0] = ((m-1.0)/m)*D_p_hat_best[n]->drift_95[i][0] + (1.0/m)*(alpha*D_p_hat_new[n]->drift_95[i][0] + (1.0-alpha)*D_p_hat_prev[n]->drift_95[i][0]);
            D_p_hat_best[n]->drift_95[i][1] = ((m-1.0)/m)*D_p_hat_best[n]->drift_95[i][1] + (1.0/m)*(alpha*D_p_hat_new[n]->drift_95[i][1] + (1.0-alpha)*D_p_hat_prev[n]->drift_95[i][1]);
        }
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
