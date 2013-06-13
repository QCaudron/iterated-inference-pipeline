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

#include "mif.h"

void mif(struct s_calc **calc, struct s_data *p_data, struct s_best *p_best, struct s_X ***J_p_X, struct s_X ***J_p_X_tmp, struct s_par **J_p_par, struct s_likelihood *p_like, gsl_vector **J_theta, gsl_vector **J_theta_tmp, double **D_theta_bart, double **D_theta_Vt, plom_f_pred_t f_pred, int is_mvn)
{
    int i, j, k;
    int m, n, nn; /*m:filtering iteration index, n: indeces N_DATA and nn N_DATA_NONAN */
    char str[STR_BUFFSIZE];

#if FLAG_VERBOSE
    int64_t time_mif_begin, time_mif_end; /* to calculate the computational time of the MIF iteration */
#endif
    int t0, t1;

#if FLAG_OMP

    int thread_id, nnp1;

#else

    void *context = zmq_ctx_new();

    void *sender = zmq_socket (context, ZMQ_PUSH);
    zmq_bind (sender, "inproc://server_sender");

    void *receiver = zmq_socket (context, ZMQ_PULL);
    zmq_bind (receiver, "inproc://server_receiver");

    void *controller = zmq_socket (context, ZMQ_PUB);
    zmq_bind (controller, "inproc://server_controller");

    struct s_thread_mif *p_thread_mif = malloc(calc[0]->n_threads*sizeof(struct s_thread_mif));
    pthread_t *worker = malloc(calc[0]->n_threads*sizeof(pthread_t));
    int nt, id, the_nt;
    int J_chunk = J/calc[0]->n_threads;
	
    for (nt = 0; nt < calc[0]->n_threads; nt++) {
	p_thread_mif[nt].thread_id = nt;       	    	    
	p_thread_mif[nt].J_chunk = J_chunk;
	p_thread_mif[nt].J = J;
	p_thread_mif[nt].p_data = p_data;
	p_thread_mif[nt].J_p_par = J_p_par;
	p_thread_mif[nt].J_p_X = J_p_X;
	p_thread_mif[nt].p_calc = calc[nt];	
	p_thread_mif[nt].p_like = p_like;
	p_thread_mif[nt].context = context;
	pthread_create (&worker[nt], NULL, worker_routine_mif_inproc, (void*) &p_thread_mif[nt]);

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

    int n_max;
    if(OPTION_IC_ONLY) {
        n_max = ((L+1) > N_DATA_NONAN) ? N_DATA_NONAN: L+1;
    } else {
        n_max = N_DATA_NONAN;
    }

    FILE *p_file_best = sfr_fopen(SFR_PATH, GENERAL_ID, "best", "w", header_best, p_data);

    FILE *p_file_mif = NULL;
    if (OPTION_TRAJ) {
        p_file_mif = sfr_fopen(SFR_PATH, GENERAL_ID, "mif", "w", header_mean_var_theoretical_mif, p_data);
    }

    print_best(p_file_best, 0, p_best, p_data, NAN);
#if ! FLAG_JSON
    fflush(p_file_best);
#endif

    gsl_matrix *var_fitted = NULL; 
    if(is_mvn){
	var_fitted = gsl_matrix_calloc(p_best->n_to_be_estimated, p_best->n_to_be_estimated);
	for(i=0; i<p_best->n_to_be_estimated; i++) {
	    for(k=0; k<p_best->n_to_be_estimated; k++) {
		gsl_matrix_set(var_fitted, i, k, gsl_matrix_get(p_best->var, p_best->to_be_estimated[i], p_best->to_be_estimated[k]));
	    }
	}
    }
    gsl_matrix *chol = (gsl_matrix *) calc[0]->method_specific_shared_data;
    void (*my_ran_proposal) (theta_t *proposed, struct s_best *p_best, gsl_matrix *var, double sd_fac, struct s_calc *p_calc);

    for(m=1; m<=M; m++) {
#if FLAG_VERBOSE
        time_mif_begin = s_clock();
#endif

        fill_theta_bart_and_Vt_mif(D_theta_bart, D_theta_Vt, p_best, p_data, m);

	if(is_mvn){
	    //compute chol only once
	    gsl_matrix_memcpy(chol, var_fitted);	    
	    gsl_matrix_scale(chol, pow(MIF_b*FREEZE, 2));

	    int status = gsl_linalg_cholesky_decomp(chol);
	    if(status == GSL_EDOM) {
		// error: matrix not positive
		print_err("Covariance matrix is not positive definite");
	    }
	   
	    //the chol computed above will be used (and passed by calc in ran_proposal_mif)
	    my_ran_proposal = &ran_proposal_chol;
	} else {
	    my_ran_proposal = &ran_proposal;
	}

        for(j=0; j<J; j++) {
            propose_safe_theta_and_load_X0(J_theta[j], p_best, p_best->var, MIF_b*FREEZE, J_p_par[j], (*J_p_X)[j], p_data, calc[0], my_ran_proposal);
        }

        t0=0; p_like->Llike_best = 0.0; p_like->n_all_fail = 0;

        for(n=0; n<n_max; n++) {

#if FLAG_JSON //for the webApp, we block at every iterations to prevent the client to be saturated with msg
            if(n % 10 == 0){
                block();
            }
#endif

            t1=p_data->times[n];

            for(j=0; j<J; j++) {
		back_transform_theta2par(J_p_par[j], J_theta[j], p_data->p_it_par_proc_par_obs_no_drift, p_data);
            }

            /*we have to use this subloop to mimate equaly spaced time step and hence set the incidence to 0 every time unit...*/
            for(nn=t0 ; nn<t1 ; nn++) {
                store_state_current_n_nn(calc, n, nn);
#if FLAG_OMP
                nnp1 = nn+1;
#pragma omp parallel for private(thread_id)
                for(j=0;j<J;j++) {
		    thread_id = omp_get_thread_num();

                    reset_inc((*J_p_X)[j], p_data);
		    f_pred((*J_p_X)[j], nn, nnp1, J_p_par[j], p_data, calc[thread_id]);

                    if(nnp1 == t1) {
                        proj2obs((*J_p_X)[j], p_data);
                        p_like->weights[j] = exp(get_log_likelihood((*J_p_X)[j], J_p_par[j], p_data, calc[thread_id]));
                    }
                }
#else 
		//send work           
		for (nt=0; nt<calc[0]->n_threads; nt++) {
		    zmq_send(sender, &nt, sizeof (int), 0);
		}

		//get results from the workers
		for (nt=0; nt<calc[0]->n_threads; nt++) {
		    zmq_recv(receiver, &the_nt, sizeof (int), 0);	       
		}

#endif



            } /* end for on nn */

            if (OPTION_PRIOR) {
                patch_likelihood_prior(p_like, p_best, J_theta, p_data, n_max, n, L);
            }

            int success = weight(p_like, n);
	    mean_var_theta_theoretical_mif(D_theta_bart[n+1], D_theta_Vt[n+1], J_theta, p_like, p_data, p_best, m, ((double) (t1-t0)), OPTION_TRAJ);

            if (OPTION_TRAJ) {
		print_mean_var_theta_theoretical_mif(p_file_mif, D_theta_bart[n+1], D_theta_Vt[n+1], p_like, p_data, m, t1);
            }

            if(success) {
		systematic_sampling(p_like, calc[0], n);
            }

	    resample_and_mut_theta_mif(p_like->select[n], J_theta, J_theta_tmp, calc, p_data, p_best, FREEZE*sqrt(((double) (t1-t0))), var_fitted, is_mvn);
	    resample_X(p_like->select[n], J_p_X, J_p_X_tmp, p_data);

	    if(n==L){
		update_fixed_lag_smoothing(p_best, p_like, J_theta, p_data);
	    }

            t0=t1;

        } /*end of for loop on the time (n)*/

        /* update theta_best */
	(m<=SWITCH) ? update_theta_best_stable_mif(p_best, D_theta_bart, p_data) : update_theta_best_king_mif(p_best, D_theta_bart, D_theta_Vt, p_data, m);

#if FLAG_VERBOSE
        time_mif_end = s_clock();
        struct s_duration t_exec = time_exec(time_mif_begin, time_mif_end);
        sprintf(str, "iteration number:%d\t logV: %g\t n_all_fail: %d\t computed in:= %dd %dh %dm %gs", m, p_like->Llike_best, p_like->n_all_fail, t_exec.d, t_exec.h, t_exec.m, t_exec.s);
        print_log(str);
#endif

        print_best(p_file_best, m, p_best, p_data, p_like->Llike_best);
#if ! FLAG_JSON
        fflush(p_file_best);
#endif

    } /*end for on m*/
    
    if(is_mvn){
        gsl_matrix_free(var_fitted);
    }

    sfr_fclose(p_file_best);

    if (OPTION_TRAJ) {
        sfr_fclose(p_file_mif);
    }

#if !FLAG_OMP
    zmq_send (controller, "KILL", 5, 0);        
    zmq_close (sender);
    zmq_close (receiver);
    zmq_close (controller);

    for(nt = 0; nt < calc[0]->n_threads; nt++){
	pthread_join(worker[nt], NULL);
    }

    free(worker);
    free(p_thread_mif);

    zmq_ctx_destroy (context);	       
#endif

}
