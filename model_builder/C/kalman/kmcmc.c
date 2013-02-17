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

#include "kalman.h"

/**
 * run KMCMC see pmcmc/pmcmc.c for doc
 */
void kmcmc(struct s_kalman *p_kalman, struct s_likelihood *p_like, struct s_mcmc_calc_data *p_mcmc_calc_data, plom_f_pred_t f_pred,  int OPTION_ACC)
{

    int m;              // iteration index
    int is_accepted;    // boolean
    double alpha;       // acceptance rate
    double sd_fac;
    int accept = 0; // number of proposed values of the parameters accepted by Metropolis Hastings */

    gsl_matrix *var; 

    // syntactical shortcut
    struct s_X *p_X = p_kalman->p_X;
    struct s_par *p_par = p_kalman->p_par;
    struct s_data *p_data = p_kalman->p_data;
    struct s_calc **calc = p_kalman->calc;
    struct s_best *p_best = p_kalman->p_best;

    // initialize time to calculate the computational time of a pMCMC iteration
#if FLAG_VERBOSE
    char str[255];
    int64_t time_pmcmc_begin, time_pmcmc_end; /* to calculate the computational time of a pMCMC iteration */
    time_pmcmc_begin = s_clock();
#endif

    // open output files
    FILE *p_file_best = sfr_fopen(SFR_PATH, GENERAL_ID, "best", "w", header_best, p_data);
    FILE *p_file_X = NULL;
    if (OPTION_TRAJ) {
        p_file_X = sfr_fopen(SFR_PATH, GENERAL_ID, "X", "w", header_X, p_data);
    }

    FILE *p_file_acc = NULL;
    if (OPTION_ACC){
        p_file_acc = sfr_fopen(SFR_PATH, GENERAL_ID, "acc", "w", header_acceptance_rates, p_data);
    }

    /////////////////////////
    // initialization step //
    /////////////////////////

    m=0;

    back_transform_theta2par(p_par, p_best->proposed, p_data->p_it_all, p_data);
    linearize_and_repeat(p_X, p_par, p_data, p_data->p_it_par_sv);
    prop2Xpop_size(p_X, p_data);
    theta_driftIC2Xdrift(p_X, p_best->proposed, p_data);

    //run Kalman
    p_like->Llike_best = run_kalman(p_X, p_best, p_par, p_kalman->p_kalman_update, p_data, calc, f_pred, p_file_X, m);
    p_like->Llike_new = p_like->Llike_best;

    //the initial iteration is "accepted"
    p_like->Llike_prev = p_like->Llike_new;
    //store the accepted value in p_best->mean
    gsl_vector_memcpy(p_best->mean, p_best->proposed);

    print_best(p_file_best, 0, p_best, p_data, p_like->Llike_best);

#if FLAG_VERBOSE
    time_pmcmc_end = s_clock();
    struct s_duration t_exec = time_exec(time_pmcmc_begin, time_pmcmc_end);
    sprintf(str, "iteration number:%d\t logV: %g\t accepted:%d computed in:= %dd %dh %dm %gs", m, p_like->Llike_best, accept, t_exec.d, t_exec.h, t_exec.m, t_exec.s);
    print_log(str);
#endif

    ////////////////
    // iterations //
    ////////////////

    for(m=1; m<M; m++) {
#if FLAG_VERBOSE
        time_pmcmc_begin = s_clock();
#endif

        increment_iteration_counters(p_mcmc_calc_data, p_best, OPTION_FULL_UPDATE);

        // web interface
#if FLAG_JSON
        if (p_mcmc_calc_data->has_cycled) {
            update_walk_rates(p_best, NULL, NULL, p_data);
            update_to_be_estimated(p_best);
        }
#endif

        // generate new theta
	var = get_var_and_sd_fac(&sd_fac, p_best, p_mcmc_calc_data, calc[0], m);
	propose_safe_theta_and_load_X0(p_best->proposed, p_best, var, sd_fac, p_par, p_X, p_data, calc[0],  (OPTION_FULL_UPDATE) ? plom_rmvnorm : ran_proposal_sequential);

        back_transform_theta2par(p_par, p_best->proposed, p_data->p_it_par_proc_par_obs_no_drift, p_data);

        //run Kalman
        p_like->Llike_best = run_kalman(p_X, p_best, p_par, p_kalman->p_kalman_update, p_data, calc, f_pred, p_file_X, m);

        p_like->Llike_new = p_like->Llike_best;

        // acceptance
        is_accepted = metropolis_hastings(p_best, p_like, &alpha, p_data, calc[0], var, sd_fac, OPTION_FULL_UPDATE);

        if (is_accepted) {
            p_like->Llike_prev = p_like->Llike_new;
            gsl_vector_memcpy(p_best->mean, p_best->proposed);
	    accept++;
        } else if(!OPTION_FULL_UPDATE) {
            //required if sequential update:
            gsl_vector_set(p_best->proposed,
                           p_best->to_be_estimated[ p_mcmc_calc_data->cycle_id ],
                           gsl_vector_get(p_best->mean, p_best->to_be_estimated[ p_mcmc_calc_data->cycle_id ]));
        }

        compute_acceptance_rates(p_best, p_mcmc_calc_data, (double) is_accepted, m);

#if FLAG_VERBOSE
        time_pmcmc_end = s_clock();
        struct s_duration t_exec = time_exec(time_pmcmc_begin, time_pmcmc_end);
        sprintf(str, "iteration number: %d (%d / %d)\t logV: %g (previous was %g) accepted: %d computed in:= %dd %dh %dm %gs", p_mcmc_calc_data->m_full_iteration, p_mcmc_calc_data->cycle_id, p_best->n_to_be_estimated, p_like->Llike_best, p_like->Llike_prev, accept, t_exec.d, t_exec.h, t_exec.m, t_exec.s);
        print_log(str);
#endif

        // append output files
        if (p_mcmc_calc_data->cycle_id >= (p_best->n_to_be_estimated -1)) { // >= instead of  == because due to the webApp all jump size can be 0.0...

	    // evaluate empirical covariance
            eval_var_emp(p_best, (double) p_mcmc_calc_data->m_full_iteration);

            print_best(p_file_best, p_mcmc_calc_data->m_full_iteration, p_best, p_data, p_like->Llike_prev);

            if (OPTION_ACC) {
		print_acceptance_rates(p_file_acc, p_mcmc_calc_data, p_mcmc_calc_data->m_full_iteration);
	    }

#if FLAG_VERBOSE
	    snprintf(str, STR_BUFFSIZE, "acceptance rate(s) at iteration %d: %g (smoothed: %g)", p_mcmc_calc_data->m_full_iteration, p_mcmc_calc_data->global_acceptance_rate, p_mcmc_calc_data->smoothed_global_acceptance_rate);
	    print_log(str);
#endif
        }

    }

    /////////////////
    // terminating //
    /////////////////

    sfr_fclose(p_file_best);
    if (OPTION_TRAJ) {
        sfr_fclose(p_file_X);
    }
    if (OPTION_ACC) {
        sfr_fclose(p_file_acc);
    }

}
