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

/**
* Computes the weight of the particles.
* Note that p_like->weights already contains the likelihood
* @return the sucess status (sucess if some particles have a likelihood > LIKE_MIN)
*/

int weight(struct s_likelihood *p_like, int n)
{

#if FLAG_WARNING
    char str[100];
#endif

    int j;

    double like_tot_n = 0.0;
    int nfailure_n = 0;
    int success = 1;

    p_like->ess_n = 0.0;

    for(j=0;j<J;j++) {
        /*compute first part of weights (non divided by sum likelihood)*/
        if (p_like->weights[j] <= pow(LIKE_MIN, N_TS)) {
            p_like->weights[j] = 0.0;
            nfailure_n += 1;
        } else {
            like_tot_n += p_like->weights[j]; //note that like_tot_n contains only like of part having a like>LIKE_MIN
            p_like->ess_n += p_like->weights[j]*p_like->weights[j]; //first part of ess computation (sum of square)
        }

    } /*end of for on j*/

    /*compute second part of weights (divided by sum likelihood)*/
    if(nfailure_n == J) {
        success = 0;
        p_like->n_all_fail += 1;
#if FLAG_WARNING
        sprintf(str,"warning: nfailure = %d, at n=%d we keep all particles and assign equal weights", nfailure_n , n);
        print_warning(str);
#endif
        p_like->Llike_best_n = LOG_LIKE_MIN*N_TS;

        double invJ=1.0/ ((double) J);
        for(j=0;j<J;j++) {
            p_like->weights[j]= invJ;
            p_like->select[n][j]= j;
        }
        p_like->ess_n = 0.0;

    } else {
        for(j=0;j<J;j++) {
            p_like->weights[j] /= like_tot_n;
        }

        p_like->Llike_best_n = log(like_tot_n / ((double) J));
        p_like->ess_n = (like_tot_n*like_tot_n)/p_like->ess_n;
    }

    p_like->Llike_best += p_like->Llike_best_n;

    return success;
}


/**
   Systematic sampling.  Systematic sampling is faster than
   multinomial sampling and introduces less monte carlo variability
*/
void systematic_sampling(struct s_likelihood *p_like, struct s_calc *p_calc, int n)
{
    unsigned int *select = p_like->select[n];
    double *prob = p_like->weights;

    int i,j;

    double ran;
    double inc = 1.0/((double) J);

    ran = gsl_ran_flat(p_calc->randgsl, 0.0, inc);
    i = 0;
    double weight_cum = prob[0];

    for(j=0; j<J; j++) {
        while(ran > weight_cum) {
            i++;
            weight_cum += prob[i];
        }
        select[j] = i;
        ran += inc;
    }
}

/**
   Multinomial sampling.
   @see systematic_sampling
*/

void multinomial_sampling(struct s_likelihood *p_like, struct s_calc *p_calc, int n)
{
    unsigned int *select = p_like->select[n];
    double *prob = p_like->weights;

    int i,j;
    unsigned int *count = init1u_set0(J);

    gsl_ran_multinomial(p_calc->randgsl, J, J, prob, count);

    i=0;
    for(j=0;j<J;j++) {
        while(count[j] >0) {
            count[j] -= 1;
            select[i++] = j;
        }
    }

    FREE(count);
}



void resample_X(unsigned int *select, struct s_X ***J_p_X, struct s_X ***J_p_X_tmp, struct s_data *p_data)
{
    int k, j;
    int size_resample_proj = N_PAR_SV*N_CAC + p_data->p_it_only_drift->nbtot;

    for(j=0;j<J;j++) {

        (*J_p_X_tmp)[j]->dt = (*J_p_X)[select[j]]->dt;

        for(k=0;k<size_resample_proj;k++) { //we don't need N_TS_INC_UNIQUE as they are present in obs
            (*J_p_X_tmp)[j]->proj[k] = (*J_p_X)[select[j]]->proj[k];
        }

        for(k=0;k<N_TS;k++) {
            (*J_p_X_tmp)[j]->obs[k] = (*J_p_X)[select[j]]->obs[k];
        }
    }

    swap_X(J_p_X, J_p_X_tmp);
}


/**
 * load X_0 for the J-1 other particles
 */
void replicate_J_p_X_0(struct s_X **J_p_X, struct s_data *p_data)
{
    int j;
    int size_proj = N_PAR_SV*N_CAC + p_data->p_it_only_drift->nbtot + N_TS_INC_UNIQUE;

    for(j=1; j<J; j++) {
        memcpy(J_p_X[j]->proj, J_p_X[0]->proj, size_proj * sizeof(double));
        J_p_X[j]->dt = J_p_X[0]->dt;
    }
}



/**
   run an SMC algorithm.

   D_J_p_X is [N_DATA+1][J]. The "+1" is for initial condition (one
   time step before first data) N_DATA+1 is indexed by nn,
   D_J_p_X[nn=0] must have been filled with the initial
   conditions.

   Note that this function assume that p_par are the same for every
   particles. This might change in the future if this function is
   extended to work for the MIF as well. See weight() on how to
   handle this case.

   After running this function:

   -p_like is filled with the likelihood values (llike at each n (n
   index N_DATA_NONAN), ess...)

   -D_J_p_X contains the J resampled particles at every nn.

   -p_hat contain the estimation of the states and obs for every nn
   and the 95%CI
*/

void run_SMC(struct s_X ***D_J_p_X, struct s_X ***D_J_p_X_tmp,
             struct s_par *p_par, struct s_hat **D_p_hat, struct s_likelihood *p_like,
             struct s_data *p_data, struct s_calc **calc,
             plom_f_pred_t f_pred,
             int option_filter, FILE *p_file_X, FILE *p_file_hat, FILE *p_file_pred_res,
             const enum plom_print print_opt
             )
{
    int j, n, np1, t0,t1;
    int thread_id;
    int size_proj = N_PAR_SV*N_CAC + p_data->p_it_only_drift->nbtot + N_TS_INC_UNIQUE;

    p_like->Llike_best = 0.0;

    //print t0 (n=-1)
    if ((print_opt & PLOM_PRINT_X) || (print_opt & PLOM_PRINT_HAT)) {
	for(j=0;j<J;j++) {
	    reset_inc(D_J_p_X[0][j], p_data);
	    proj2obs(D_J_p_X[0][j], p_data);
	}       

	if (print_opt & PLOM_PRINT_HAT) {
	    compute_hat_nn(D_J_p_X[0], &p_par, p_data, calc, D_p_hat[0], 1, -1, 0.0);
	    print_p_hat(p_file_hat, NULL, D_p_hat[0], p_data, 0.0);
	}

	if (print_opt & PLOM_PRINT_X) {
	    print_X(p_file_X, &p_par, D_J_p_X[0], p_data, calc[0], 1, 0, 0, -1, 0.0);
	}
    }

    for(n=0; n<p_data->nb_obs; n++) {

#if FLAG_JSON //for the webApp, we block at every iterations to prevent the client to be saturated with msg
        if (print_opt & PLOM_PRINT_HAT) {
            if(n % 10 == 0){
                block();
            }
        }
#endif

	np1 = n+1;
	t0 = p_data->times[n];
	t1 = p_data->times[np1];
		    
	//we are going to overwrite the content of the [np1] pointer: initialise it with values from [n]
	for(j=0;j<J;j++) {
	    memcpy(D_J_p_X[np1][j]->proj, D_J_p_X[n][j]->proj, size_proj * sizeof(double));
	    D_J_p_X[np1][j]->dt = D_J_p_X[n][j]->dt;
	}

#if FLAG_OMP
#pragma omp parallel for private(thread_id)
#endif
	for(j=0;j<J;j++) {
#if FLAG_OMP
	    thread_id = omp_get_thread_num();
#else
	    thread_id = 0;
#endif

	    reset_inc(D_J_p_X[np1][j], p_data);
	    (*f_pred)(D_J_p_X[np1][j], t0, t1, p_par, p_data, calc[thread_id]);
	    //round_inc(D_J_p_X[nnp1][j]);

	    proj2obs(D_J_p_X[np1][j], p_data);

	    if(p_data->data_ind[n]->n_nonan) {
		p_like->weights[j] = exp(get_log_likelihood(D_J_p_X[np1][j], p_par, p_data, calc[thread_id], n, t1));
	    }
	}

        if(option_filter && p_data->data_ind[n]->n_nonan) {

            if(weight(p_like, n)) {
                systematic_sampling(p_like, calc[0], n);
            }

	    if (print_opt & PLOM_PRINT_HAT) {
		compute_hat(D_J_p_X[np1], p_par, p_data, calc, D_p_hat[np1], p_like->weights, n, t1);
	    }

            resample_X(p_like->select[n], &(D_J_p_X[np1]), &(D_J_p_X_tmp[np1]), p_data);

            if (print_opt & PLOM_PRINT_PRED_RES) {
                print_prediction_residuals(p_file_pred_res, &p_par, p_data, calc[0], D_J_p_X[np1], p_like->Llike_best_n, p_like->ess_n, 1, n, t1);
            }

        } else if (print_opt & PLOM_PRINT_HAT) { 
            //we do not filter or all data ara NaN (no info).
	    compute_hat_nn(D_J_p_X[np1], &p_par, p_data, calc, D_p_hat[np1], 1, n, t1);
        }

        if (print_opt & PLOM_PRINT_HAT) {
            print_p_hat(p_file_hat, NULL, D_p_hat[np1], p_data, t1);
        }

	if (print_opt & PLOM_PRINT_X) {
	    print_X(p_file_X, &p_par, D_J_p_X[np1], p_data, calc[0], 1, 0, 0, n, t1);
	}

    }

}



/**
 *  same as run_SMC but deleguates work to plom workers. Each
 *  worker receive Jchunk particles
 */
void run_SMC_zmq(struct s_X ***D_J_p_X, struct s_X ***D_J_p_X_tmp, struct s_par *p_par, struct s_hat **D_p_hat, struct s_likelihood *p_like, struct s_data *p_data, struct s_calc **calc, plom_f_pred_t f_pred, int Jchunk, const enum plom_print print_opt, void *sender, void *receiver, void *controller)
{
    int j, n, np1;
    int the_j;

    p_like->Llike_best = 0.0; p_like->n_all_fail = 0;

    if (print_opt & PLOM_PRINT_HAT) {
	for(j=0;j<J;j++) {
	    for(j=0;j<J;j++) {
		reset_inc(D_J_p_X[0][j], p_data);
		proj2obs(D_J_p_X[0][j], p_data);
	    }
	}
	compute_hat(D_J_p_X[0], p_par, p_data, calc, D_p_hat[0], p_like->weights, -1, 0.0);
    }


    for (n=0; n < p_data->nb_obs; n++) {

	np1 = n+1;

	//send work
	//we are going to overwrite the content of the [nnp1] pointer: initialise it with values from [nn]
	for (j=0;j<J;j++) {
	    //we try to minimize the number of times we send the parameters...
	    if ( (j % Jchunk) == 0) {
		zmq_send(sender, &n, sizeof (int), ZMQ_SNDMORE);
		send_par(sender, p_par, p_data, ZMQ_SNDMORE);
	    }

	    zmq_send(sender, &j, sizeof (int), ZMQ_SNDMORE);                   	       	       
	    send_X(sender, D_J_p_X[n][j], p_data, ( ((j+1) % Jchunk) == 0) ? 0: ZMQ_SNDMORE);
	    //printf("part %d sent %d\n", j, ( ((j+1) % Jchunk) == 0) ? 0: ZMQ_SNDMORE);
	}

	//get results from the workers
	for (j=0; j<J; j++) {
	    zmq_recv(receiver, &the_j, sizeof (int), 0);
	    recv_X(D_J_p_X[np1][ the_j ], p_data, receiver);
	    zmq_recv(receiver, &(p_like->weights[the_j]), sizeof (double), 0);
	    //printf("part  %d received\n", the_j);
	}

        if( p_data->data_ind[n]->n_nonan) {

	    if (weight(p_like, n)) {
		systematic_sampling(p_like, calc[0], n);
	    }
	    if (print_opt & PLOM_PRINT_HAT) {
		compute_hat(D_J_p_X[np1], p_par, p_data, calc, D_p_hat[np1], p_like->weights, n, p_data->times[np1]);
	    }
	    resample_X(p_like->select[n], &(D_J_p_X[np1]), &(D_J_p_X_tmp[np1]), p_data);

	} else if (print_opt & PLOM_PRINT_HAT) {

	    compute_hat_nn(D_J_p_X[np1], &p_par, p_data, calc, D_p_hat[np1], 1, n, p_data->times[np1]);

	}

    }
}


void run_SMC_zmq_inproc(struct s_X ***D_J_p_X, struct s_X ***D_J_p_X_tmp, struct s_par *p_par, struct s_hat **D_p_hat, struct s_likelihood *p_like, struct s_data *p_data, struct s_calc **calc, plom_f_pred_t f_pred, int option_filter, FILE *p_file_X, FILE *p_file_hat, FILE *p_file_pred_res, const enum plom_print print_opt, void *sender, void *receiver, void *controller)
{
    int j, n, np1, t1;
    int nt;
    int the_nt;
    int size_proj = N_PAR_SV*N_CAC + p_data->p_it_only_drift->nbtot + N_TS_INC_UNIQUE;

    p_like->Llike_best = 0.0; p_like->n_all_fail = 0;

    //print t0 (n=-1)
    if ((print_opt & PLOM_PRINT_X) || (print_opt & PLOM_PRINT_HAT)) {	
	for(j=0;j<J;j++) {
	    reset_inc(D_J_p_X[0][j], p_data);
	    proj2obs(D_J_p_X[0][j], p_data);
	}       
	if (print_opt & PLOM_PRINT_HAT) {
	    compute_hat_nn(D_J_p_X[0], &p_par, p_data, calc, D_p_hat[0], 1, -1, 0.0);
	    print_p_hat(p_file_hat, NULL, D_p_hat[0], p_data, 0.0);
	}

	if (print_opt & PLOM_PRINT_X) {
	    print_X(p_file_X, &p_par, D_J_p_X[0], p_data, calc[0], 1, 0, 0, -1, 0.0);
	}
    }


    for (n=0; n < p_data->nb_obs; n++) {

#if FLAG_JSON //for the webApp, we block at every iterations to prevent the client to be saturated with msg
        if (print_opt & PLOM_PRINT_HAT) {
            if(n % 10 == 0){
                block();
            }
        }
#endif

	np1 = n+1;
	t1 = p_data->times[np1];

	//we are going to overwrite the content of the [np1] pointer: initialise it with values from [n]
	for(j=0;j<J;j++) {
	    memcpy(D_J_p_X[np1][j]->proj, D_J_p_X[n][j]->proj, size_proj * sizeof(double));
	    D_J_p_X[np1][j]->dt = D_J_p_X[n][j]->dt;
	}

	//send work           
	for (nt=0; nt<calc[0]->n_threads; nt++) {
	    zmq_send(sender, &nt, sizeof (int), ZMQ_SNDMORE);
	    zmq_send(sender, &n, sizeof (int), 0);
	}

	//get results from the workers
	for (nt=0; nt<calc[0]->n_threads; nt++) {
	    zmq_recv(receiver, &the_nt, sizeof (int), 0);	       
	}


        if(option_filter && p_data->data_ind[n]->n_nonan) {

            if(weight(p_like, n)) {
                systematic_sampling(p_like, calc[0], n);
            }

            //!! time indexes: D_J_p_X is [N_DATA+1], *D_p_hat->... are in [N_DATA] so we have to be carrefull!
	    if (print_opt & PLOM_PRINT_HAT) {
		compute_hat(D_J_p_X[np1], p_par, p_data, calc, D_p_hat[np1], p_like->weights, n, t1);
	    }

	    resample_X(p_like->select[n], &(D_J_p_X[np1]), &(D_J_p_X_tmp[np1]), p_data);

            if (print_opt & PLOM_PRINT_PRED_RES) {
                print_prediction_residuals(p_file_pred_res, &p_par, p_data, calc[0], D_J_p_X[np1], p_like->Llike_best_n, p_like->ess_n, 1, n, t1);
            }

        } else if (print_opt & PLOM_PRINT_HAT) { 
            //we do not filter or all data ara NaN (no info).
	    compute_hat_nn(D_J_p_X[np1], &p_par, p_data, calc, D_p_hat[np1], 1, n, t1);
        }


        if (print_opt & PLOM_PRINT_HAT) {
            print_p_hat(p_file_hat, NULL, D_p_hat[np1], p_data, t1);
        }

	if (print_opt & PLOM_PRINT_X) {
	    print_X(p_file_X, &p_par, D_J_p_X[np1], p_data, calc[0], 1, 0, 0, n, t1);
	}
    } 
}
