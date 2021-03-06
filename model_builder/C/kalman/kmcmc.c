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
void kmcmc(struct s_kalman *p_kalman, struct s_likelihood *p_like, struct s_mcmc_calc_data *p_mcmc_calc_data, plom_f_pred_t f_pred, enum plom_print print_opt, int thin_traj)
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

    char str[255];
    int64_t time_kmcmc_begin, time_kmcmc_end; /* to calculate the computational time of a pMCMC iteration */
    if (!(print_opt & PLOM_QUIET)) {
	time_kmcmc_begin = s_clock();
    }

    // open output files
    FILE *p_file_trace = plom_fopen(SFR_PATH, GENERAL_ID, "trace", "w", header_trace, p_data);
    FILE *p_file_X = NULL;
    int is_print_x = (print_opt & PLOM_PRINT_X);
    if (is_print_x) {
        p_file_X = plom_fopen(SFR_PATH, GENERAL_ID, "X", "w", header_X, p_data);
    }

    FILE *p_file_acc = NULL;
    if (print_opt & PLOM_PRINT_ACC){
        p_file_acc = plom_fopen(SFR_PATH, GENERAL_ID, "acc", "w", header_acceptance_rates, p_data);
    }

    /////////////////////////
    // initialization step //
    /////////////////////////

    m=0;

    back_transform_theta2par(p_par, p_best->proposed, p_data->p_it_all, p_data);
    linearize_and_repeat(p_X, p_par, p_data, p_data->p_it_par_sv);
    prop2Xpop_size(p_X, p_data, calc[0]);
    theta_driftIC2Xdrift(p_X, p_best->proposed, p_data);

    //run Kalman    
    p_like->Llike_best = run_kalman(p_X, p_best, p_par, p_kalman->p_kalman_update, p_data, calc, f_pred, m, p_file_X, NULL, NULL, print_opt);

    p_like->Llike_new = p_like->Llike_best;

    //the initial iteration is "accepted"
    p_like->Llike_prev = p_like->Llike_new;
    //store the accepted value in p_best->mean
    gsl_vector_memcpy(p_best->mean, p_best->proposed);

    print_trace(p_file_trace, 0, p_best, p_data, p_like->Llike_best);

    if (!(print_opt & PLOM_QUIET)) {
	time_kmcmc_end = s_clock();
	struct s_duration t_exec = time_exec(time_kmcmc_begin, time_kmcmc_end);
	sprintf(str, "iteration number:%d\t logV: %g\t accepted:%d computed in:= %dd %dh %dm %gs", m, p_like->Llike_best, accept, t_exec.d, t_exec.h, t_exec.m, t_exec.s);
	print_log(str);
    }

    ////////////////
    // iterations //
    ////////////////

    for(m=1; m<M; m++) {
	if (!(print_opt & PLOM_QUIET)) {
	    time_kmcmc_begin = s_clock();
	}

        increment_iteration_counters(p_mcmc_calc_data, p_best, OPTION_FULL_UPDATE);

        // web interface
#if FLAG_JSON
        if (p_mcmc_calc_data->has_cycled) {
            update_walk_rates(p_best, p_data);
            update_to_be_estimated(p_best);
        }
#endif

        // generate new theta
	var = get_var_and_sd_fac(&sd_fac, p_best, p_mcmc_calc_data, calc[0], m);
	propose_safe_theta_and_load_X0(p_best->proposed, p_best, var, sd_fac, p_par, p_X, p_data, calc[0],  (OPTION_FULL_UPDATE) ? plom_rmvnorm : ran_proposal_sequential);

        back_transform_theta2par(p_par, p_best->proposed, p_data->p_it_par_proc_par_obs_no_drift, p_data);

        //run Kalman

	//overwrite print_opt
	if ( (is_print_x) && ( (p_mcmc_calc_data->m_full_iteration % thin_traj) == 0)  && p_data->nb_obs && (p_mcmc_calc_data->cycle_id >= (p_best->n_to_be_estimated -1)) ) {
	    print_opt |= PLOM_PRINT_X;
	} else {
	    print_opt &= ~PLOM_PRINT_X;
	}

        p_like->Llike_best = run_kalman(p_X, p_best, p_par, p_kalman->p_kalman_update, p_data, calc, f_pred, m, p_file_X, NULL, NULL, print_opt);

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

	if (!(print_opt & PLOM_QUIET)) {
	    time_kmcmc_end = s_clock();
	    struct s_duration t_exec = time_exec(time_kmcmc_begin, time_kmcmc_end);
	    sprintf(str, "iteration number: %d (%d / %d)\t logV: %g (previous was %g) accepted: %d computed in:= %dd %dh %dm %gs", p_mcmc_calc_data->m_full_iteration, p_mcmc_calc_data->cycle_id, p_best->n_to_be_estimated, p_like->Llike_best, p_like->Llike_prev, accept, t_exec.d, t_exec.h, t_exec.m, t_exec.s);
	    print_log(str);
	}

        // append output files
        if (p_mcmc_calc_data->cycle_id >= (p_best->n_to_be_estimated -1)) { // >= instead of  == because due to the webApp all jump size can be 0.0...

	    // evaluate empirical covariance
            eval_var_emp(p_best, (double) p_mcmc_calc_data->m_full_iteration);

            print_trace(p_file_trace, p_mcmc_calc_data->m_full_iteration, p_best, p_data, p_like->Llike_prev);

	    if (print_opt & PLOM_PRINT_ACC) {
		print_acceptance_rates(p_file_acc, p_mcmc_calc_data, p_mcmc_calc_data->m_full_iteration);
	    }

	    if (!(print_opt & PLOM_QUIET)) {
		snprintf(str, STR_BUFFSIZE, "acceptance rate(s) at iteration %d: %g (smoothed: %g)", p_mcmc_calc_data->m_full_iteration, p_mcmc_calc_data->global_acceptance_rate, p_mcmc_calc_data->smoothed_global_acceptance_rate);
		print_log(str);
	    }
        }
    }

    /////////////////
    // terminating //
    /////////////////

    plom_fclose(p_file_trace);

    if (is_print_x) {
        plom_fclose(p_file_X);
    }
    if (print_opt & PLOM_PRINT_ACC) {
        plom_fclose(p_file_acc);
    }

}
