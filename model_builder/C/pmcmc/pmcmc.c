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
 * run SMC
 */
void run_propag(struct s_X ***D_J_p_X, struct s_X ***D_J_p_X_tmp, struct s_par *p_par, struct s_hat ***D_p_hat_new,
                struct s_likelihood *p_like, struct s_data *p_data, struct s_calc **calc, plom_f_pred_t f_pred,
                void *sender, void *receiver, void *controller, const enum plom_print print_opt)
{
    if (OPTION_PIPELINE) {
	run_SMC_zmq(D_J_p_X, D_J_p_X_tmp, p_par, *D_p_hat_new, p_like, p_data, calc, f_pred, JCHUNK, print_opt, sender, receiver, controller);
    } else {
#if FLAG_OMP
	run_SMC(D_J_p_X, D_J_p_X_tmp, p_par, *D_p_hat_new, p_like, p_data, calc, f_pred, 1, NULL, NULL, NULL, print_opt);
#else       
	run_SMC_zmq_inproc(D_J_p_X, D_J_p_X_tmp, p_par, *D_p_hat_new, p_like, p_data, calc, f_pred, 1, NULL, NULL, NULL, print_opt, sender, receiver, controller);
#endif
    }
}


void pmcmc(struct s_best *p_best, struct s_X ***D_J_p_X, struct s_X ***D_J_p_X_tmp, struct s_par *p_par, struct s_hat ***D_p_hat_prev, struct s_hat ***D_p_hat_new, struct s_hat **D_p_hat_best, struct s_likelihood *p_like, struct s_data *p_data, struct s_calc **calc, plom_f_pred_t f_pred, const enum plom_print print_opt, int thin_traj)
{
    int m;              // iteration index
    int is_accepted;    // boolean
    double alpha;       // acceptance rate
    double sd_fac;
    int accept = 0; // number of proposed values of the parameters accepted by Metropolis Hastings */

    gsl_matrix *var; 

    // syntactical shortcut
    struct s_mcmc_calc_data *p_mcmc_calc_data =  (struct s_mcmc_calc_data *) calc[0]->method_specific_shared_data;

    char str[STR_BUFFSIZE];

    // initialize time to calculate the computational time of a pMCMC iteration
#if FLAG_VERBOSE
    int64_t time_pmcmc_begin, time_pmcmc_end; /* to calculate the computational time of a pMCMC iteration */
    time_pmcmc_begin = s_clock();
#endif

    // open output files
    FILE *p_file_best = sfr_fopen(SFR_PATH, GENERAL_ID, "best", "w", header_best, p_data);
    FILE *p_file_X = NULL;
    if (print_opt & PLOM_PRINT_X_SMOOTH) {
        p_file_X = sfr_fopen(SFR_PATH, GENERAL_ID, "X", "w", header_X, p_data);
    }

    FILE *p_file_acc = NULL;
    if (print_opt & PLOM_PRINT_ACC){
        p_file_acc = sfr_fopen(SFR_PATH, GENERAL_ID, "acc", "w", header_acceptance_rates, p_data);
    }

    void *context = NULL;
    void *sender = NULL;
    void *receiver = NULL;
    void *controller = NULL;

#if !FLAG_OMP
    struct s_thread_smc *p_thread_smc = NULL;
    pthread_t *worker = NULL;
    int nt;
#endif

    if (OPTION_PIPELINE) {

#if FLAG_VERBOSE
        print_log("setting up zmq sockets...");
#endif
        context = zmq_ctx_new ();

        //  Socket to send messages on
        sender = zmq_socket (context, ZMQ_PUSH);
        zmq_bind (sender, "tcp://*:5557");

        //  Socket to receive messages on
        receiver = zmq_socket (context, ZMQ_PULL);
        zmq_bind (receiver, "tcp://*:5558");

        //  Socket for worker control
        controller = zmq_socket (context, ZMQ_PUB);
        zmq_bind (controller, "tcp://*:5559");

    } else {

#if !FLAG_OMP

	context = zmq_ctx_new ();
	sender = zmq_socket (context, ZMQ_PUSH);
	zmq_bind (sender, "inproc://server_sender");

	receiver = zmq_socket (context, ZMQ_PULL);
	zmq_bind (receiver, "inproc://server_receiver");

	controller = zmq_socket (context, ZMQ_PUB);
	zmq_bind (controller, "inproc://server_controller");

	p_thread_smc = malloc(calc[0]->n_threads*sizeof(struct s_thread_smc));
	worker = malloc(calc[0]->n_threads*sizeof(pthread_t));
	int id;
	int J_chunk = J/calc[0]->n_threads;
	
	for (nt = 0; nt < calc[0]->n_threads; nt++) {
	    p_thread_smc[nt].thread_id = nt;       	    	    
	    p_thread_smc[nt].J_chunk = J_chunk;
	    p_thread_smc[nt].J = J;
	    p_thread_smc[nt].p_data = p_data;
	    p_thread_smc[nt].p_par = p_par;
	    p_thread_smc[nt].D_J_p_X = D_J_p_X;
	    p_thread_smc[nt].p_calc = calc[nt];	
	    p_thread_smc[nt].p_like = p_like;
	    p_thread_smc[nt].context = context;
	    pthread_create (&worker[nt], NULL, worker_routine_smc_inproc, (void*) &p_thread_smc[nt]);
	    snprintf(str, STR_BUFFSIZE, "worker %d started", nt);
	    print_log(str);
	}

	//wait that all worker are connected
	for (nt = 0; nt < calc[0]->n_threads; nt++) {
	    zmq_recv(receiver, &id, sizeof (int), 0);
	    snprintf(str, STR_BUFFSIZE, "worker %d connected", id);
	    print_log(str);
	}
	       
#endif

    }

    /////////////////////////
    // initialization step //
    /////////////////////////

    m=0;

    // initialize SMC arguments (particle 0)
    back_transform_theta2par(p_par, p_best->proposed, p_data->p_it_all, p_data);
    linearize_and_repeat(D_J_p_X[0][0], p_par, p_data, p_data->p_it_par_sv);
    prop2Xpop_size(D_J_p_X[0][0], p_data, calc[0]);
    theta_driftIC2Xdrift(D_J_p_X[0][0], p_best->proposed, p_data);

    //load X_0 for the J-1 other particles
    replicate_J_p_X_0(D_J_p_X[0], p_data);

    //run SMC
    run_propag(D_J_p_X, D_J_p_X_tmp, p_par, D_p_hat_new, p_like, p_data, calc, f_pred, sender, receiver, controller, print_opt);

    p_like->Llike_new = p_like->Llike_best;

    if ( ( print_opt & PLOM_PRINT_X_SMOOTH ) && p_data->nb_obs ) {
        sample_traj_and_print(p_file_X, D_J_p_X, p_par, p_data, p_like, calc[0], 0);
    }

    //the initial iteration is "accepted"
    p_like->Llike_prev = p_like->Llike_new;
    //store the accepted value in p_best->mean
    gsl_vector_memcpy(p_best->mean, p_best->proposed);

    print_best(p_file_best, 0, p_best, p_data, p_like->Llike_best);

    // print iteration info
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
            update_walk_rates(p_best, p_data);
            update_to_be_estimated(p_best);
        }
#endif

	//currently depreciated (waiting for compute quantile online (https://github.com/plom-io/plom-sfi/issues/9))
        //swap_D_p_hat(D_p_hat_prev, D_p_hat_new);

        // generate new theta
	var = get_var_and_sd_fac(&sd_fac, p_best, p_mcmc_calc_data, calc[0], m);
	propose_safe_theta_and_load_X0(p_best->proposed, p_best, var, sd_fac, p_par, D_J_p_X[0][0], p_data, calc[0],  (OPTION_FULL_UPDATE) ? plom_rmvnorm : ran_proposal_sequential);

        //load X_0 for the J-1 other particles
        replicate_J_p_X_0(D_J_p_X[0], p_data);

        back_transform_theta2par(p_par, p_best->proposed, p_data->p_it_par_proc_par_obs_no_drift, p_data);

        //run SMC
        run_propag(D_J_p_X, D_J_p_X_tmp, p_par, D_p_hat_new, p_like, p_data, calc, f_pred, sender, receiver, controller, print_opt);

        p_like->Llike_new = p_like->Llike_best;

        // acceptance
        is_accepted = metropolis_hastings(p_best, p_like, &alpha, p_data, calc[0], var, sd_fac, OPTION_FULL_UPDATE);
	//currently depreciated (waiting for compute quantile online (https://github.com/plom-io/plom-sfi/issues/9))
        //compute_best_traj(D_p_hat_best, *D_p_hat_prev, *D_p_hat_new, p_data, (alpha>1.0) ? 1.0 : alpha, (double) m);

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
	snprintf(str, STR_BUFFSIZE, "iteration number: %d (%d / %d)\t logV: %g (previous was %g) accepted: %d computed in:= %dd %dh %dm %gs", p_mcmc_calc_data->m_full_iteration, p_mcmc_calc_data->cycle_id, p_best->n_to_be_estimated, p_like->Llike_best, p_like->Llike_prev, accept, t_exec.d, t_exec.h, t_exec.m, t_exec.s);
	print_log(str);
#endif


        if (p_mcmc_calc_data->cycle_id >= (p_best->n_to_be_estimated -1)) { // >= instead of  == because due to the webApp all jump size can be 0.0...
            // evaluate empirical covariance
            eval_var_emp(p_best, (double) p_mcmc_calc_data->m_full_iteration);

	    if ( (print_opt & PLOM_PRINT_X_SMOOTH) && ( (p_mcmc_calc_data->m_full_iteration % thin_traj) == 0)  && p_data->nb_obs ) {
		sample_traj_and_print(p_file_X, D_J_p_X, p_par, p_data, p_like, calc[0], p_mcmc_calc_data->m_full_iteration);
	    }

            print_best(p_file_best, p_mcmc_calc_data->m_full_iteration, p_best, p_data, p_like->Llike_prev);

	    if (print_opt & PLOM_PRINT_ACC) {
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

    if (print_opt & PLOM_PRINT_X_SMOOTH) {
        sfr_fclose(p_file_X);
    }
    if (print_opt & PLOM_PRINT_ACC) {
        sfr_fclose(p_file_acc);
    }


    if (OPTION_PIPELINE){

#if FLAG_VERBOSE
        print_log("killing the workers...");
#endif
	zmq_send (controller, "KILL", 5, 0);        

#if FLAG_VERBOSE
        print_log("closing zmq sockets...");
#endif
        zmq_close (sender);
        zmq_close (receiver);
        zmq_close (controller);

        zmq_ctx_destroy (context);

    } else {

#if !FLAG_OMP

	zmq_send (controller, "KILL", 5, 0);        
	zmq_close (sender);
	zmq_close (receiver);
	zmq_close (controller);

	for(nt = 0; nt < calc[0]->n_threads; nt++){
	    pthread_join(worker[nt], NULL);
	}

	free(worker);
	free(p_thread_smc);

	zmq_ctx_destroy (context);	       
#endif

    }

}
