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
 * For eval_jac
 */
double drift_derivative(double jac_tpl, double jac_der, struct s_router *r, int cac)
{

    /*
      Rational basis
      we have an equation (for instance dI/dt) named eq and let's say that we are interested in its derivative against v (we assume that v follows a diffusion)'
      The template gives us d eq/d v (jac_tpl)
      However, as v can be transform (let's say log here) we want d eq / d log(v)
      The chain rule gives us:
      d eq/ dv = d eq / d log(v) * d log(v)/dv = jac_tpl
      so
      d eq / d log(v) = ( d eq / dv ) / ( d log(v) / dv)

      so in term of C:
      d eq / d log(v) = jac_tpl / r->f_derivative(v, ..)
      jac_der is the C term of v, provided by the template

      As v (jac_der) is in the scale of s_par, in case of logit_ab
      transfo, we need to provide a and b in the scale of s_par, that
      is router min_z and max_z
     */


    if(jac_tpl){
        int g = r->map[cac];
        return jac_tpl / r->f_derivative[g](jac_der, r->min_z[g], r->max_z[g]);
    }

    return 0.0;
}



/**
 * transform p_X to x_k: concatenation of non-overlapping components of:
 * - X->proj
 * - X->obs
 * - X->drift
 * @param xk the concatenated vector
 * @param p_X the X vector
 */
void X2xk(gsl_vector *xk, struct s_X *p_X, struct s_data *p_data)
{
    struct s_iterator *p_it = p_data->p_it_only_drift;
    int i;

    //proj (without drift)
    for(i=0; i<N_PAR_SV*N_CAC; i++)
        gsl_vector_set(xk, i, p_X->proj[i]);

    //obs
    for(i=0; i<N_TS; i++) {
        gsl_vector_set(xk, N_PAR_SV*N_CAC + i, p_X->obs[i]);
    }

    //drift
    for(i=0; i<p_it->nbtot; i++) {
	gsl_vector_set(xk, N_PAR_SV*N_CAC + N_TS + i, p_X->proj[N_PAR_SV*N_CAC+i]);
    }
}

/**
 * reconstruct s_X from xk and set to 0.0 terms of proj that became negative due to Kalman.
 * Note that only N_PAR_SV*N_CAC of proj and drift are necessary. All the observed variable are useless here
 * @param p_X the X vector to be constructed
 * @param xk the concatenated vector
 */
void xk2X(struct s_X *p_X, gsl_vector *xk, struct s_data *p_data, struct s_calc *p_calc, const double t)
{
    struct s_iterator *p_it = p_data->p_it_only_drift;
    int i, cac;
    double sumsv;

    // sanitising first: recover negative state variables, and rescale if total population is reached
    for(i=0; i<N_PAR_SV*N_CAC; i++) {
	gsl_vector_set(xk, i, (gsl_vector_get(xk, i) > 0.0) ? gsl_vector_get(xk, i) : 0.0) ;
    }
    if(!POP_SIZE_EQ_SUM_SV){
	for(cac=0; cac<N_CAC; cac++){
	    sumsv = 0.0;
	    for(i=0; i<N_PAR_SV; i++) {
		sumsv += gsl_vector_get(xk, i*N_CAC +cac);
	    }
	    if( (gsl_spline_eval(p_calc->spline[0][cac],t,p_calc->acc[0][cac]) - sumsv) < 0 ){
		for(i=0; i<N_PAR_SV; i++) {
		    gsl_vector_set(xk, i*N_CAC +cac, gsl_vector_get(xk, i*N_CAC +cac)/sumsv * gsl_spline_eval(p_calc->spline[0][cac],t,p_calc->acc[0][cac]));
		}
	    }
	}
    }

    //proj (without drift)
    for(i=0; i<N_PAR_SV*N_CAC; i++) {
        p_X->proj[i] = (gsl_vector_get(xk, i) > 0.0) ? gsl_vector_get(xk, i) : 0.0 ;
	gsl_vector_set(xk, i, (gsl_vector_get(xk, i) > 0.0) ? gsl_vector_get(xk, i) : 0.0) ;
    }

    //drift
    for(i=0; i<p_it->nbtot; i++) {
	p_X->proj[N_PAR_SV*N_CAC + i] = gsl_vector_get(xk, N_PAR_SV*N_CAC + N_TS + i) ;
    }
}

/**
 * test wether a state variable (including remainder), has negative value
 * @param a list beginning by the state variables (X->proj, xk...)
 * @return a plom error code
 */
plom_err_code test_all_sv_pos(gsl_vector *xk, struct s_data *p_data, struct s_calc *p_calc, const double t)
{
    int is_err = 0;
    int i, cac;
    double sumsv;

    for(i=0; i<N_PAR_SV*N_CAC; i++) {
	if(gsl_vector_get(xk, i) < 0.0){
	    is_err = 0;
	};
    }

    for(cac=0; cac<N_CAC; cac++){
	sumsv = 0.0;
	for(i=0; i<N_PAR_SV; i++) {
	    sumsv += gsl_vector_get(xk, i*N_CAC +cac);
	}
	if( (gsl_spline_eval(p_calc->spline[0][cac],t,p_calc->acc[0][cac]) - sumsv) < 0 ){
	    is_err = 0;
	}
    }
    return (is_err) ? PLOM_ERR_LIKE: PLOM_SUCCESS;
}

/**
 * get total population in SV
 * @param a list beginning by the state variables (X->proj, xk...)
 * @return the total population
 */
double get_total_pop(double *X)
{
    int cac;
    double t_p = 0.0;
    for(cac=0;cac<N_CAC;cac++) {
        t_p += sum_SV(X, cac);
    }
    return t_p;
}


double log_transf_correc(gsl_vector *mean, gsl_matrix *var, struct s_router **routers)
{
    char str[STR_BUFFSIZE];
    int i, k;

    double p_tmp, Lp;
    p_tmp=0.0, Lp=0.0;

    int offset = 0;

    for(i=0; i<(N_PAR_SV+N_PAR_PROC+N_PAR_OBS); i++) {
        for(k=0; k<routers[i]->n_gp; k++) {
            if(gsl_matrix_get(var, offset, offset) >0.0) {

                p_tmp = 1./((*(routers[i]->f_inv_derivative[k]))(gsl_vector_get(mean, offset), routers[i]->min[k], routers[i]->max[k]));

                //check for numerical issues
                if((isinf(p_tmp)==1) || (isnan(p_tmp)==1)) {
#if FLAG_VERBOSE
                    snprintf(str, STR_BUFFSIZE, "error prob_prior computation, p=%g\n", p_tmp);
                    print_err(str);
#endif
                    p_tmp=LIKE_MIN;
                } else if(p_tmp <= LIKE_MIN) {
                    p_tmp = LIKE_MIN ;
                }
                Lp += log(p_tmp);
            }
        }
        offset++;
    }

    return(Lp);
}



/**
 *  reset incidence-related rows and columns to 0
 */
void reset_inc_cov(gsl_matrix *Ct)
{
    int oi,oii;
    for (oi=N_PAR_SV*N_CAC; oi<N_PAR_SV*N_CAC+N_TS_INC; oi++) {
        for (oii=0; oii<Ct->size1; oii++) {
            gsl_matrix_set(Ct,oi,oii,0.0);	// set row to 0
            gsl_matrix_set(Ct,oii,oi,0.0);	// set column to 0
        }
    }
}


/**
 *   run an extended Kalman filter and returns the log likelihood
 */
double run_kalman(struct s_X *p_X, struct s_best *p_best, struct s_par *p_par, struct s_kalman_update *p_kalman_update, struct s_data *p_data, struct s_calc **calc, plom_f_pred_t f_pred, int m, FILE *p_file_X, FILE *p_file_hat, FILE *p_file_pred_res, const enum plom_print print_opt)
{

    int n, np1, t0, t1;
    int ts;

    double like, log_lik, log_lik_temp;

    struct s_data_ind **data_ind = p_data->data_ind;

    t0=0;
    log_lik = 0.0;
    p_kalman_update->sc_st = 0.0;
    p_kalman_update->sc_pred_error = 0.0;

    gsl_matrix_view Ct = gsl_matrix_view_array(&(p_X->proj[N_PAR_SV*N_CAC + p_data->p_it_only_drift->nbtot + N_TS_INC_UNIQUE]), N_KAL, N_KAL);
    gsl_matrix_set_zero(&Ct.matrix);


    //debug Q:
//    struct s_kalman_specific_data *p_kalman_specific_data = (struct s_kalman_specific_data *) calc[0]->method_specific_thread_safe_data;
//    p_kalman_specific_data->eval_Q(p_kalman_specific_data->Q, p_X->proj, p_par, p_data, calc[0], p_kalman_specific_data, 1.0);
//    for(n=0; n<p_kalman_specific_data->Q->size1; n++){
//	for(nn=0; nn<=n; nn++){
//	    printf("Q[%d][%d] %g\t %s\n", n, nn, gsl_matrix_get(p_kalman_specific_data->Q, n, nn),
//		   (gsl_matrix_get(p_kalman_specific_data->Q, n, nn) == gsl_matrix_get(p_kalman_specific_data->Q, nn, n)) ? "ok" : "NO!");
//	}
//    }


    //print t0
    if ((print_opt & PLOM_PRINT_X) || (print_opt & PLOM_PRINT_HAT)) {
	reset_inc(p_X, p_data);	
	proj2obs(p_X, p_data);

	if (print_opt & PLOM_PRINT_X) {
	    printf("e\n");
	    print_X(p_file_X, &p_par, &p_X, p_data, calc[0], 1, 1, m, -1, 0.0);
	}

	if (print_opt & PLOM_PRINT_HAT) {
	    X2xk(p_kalman_update->xk, p_X, p_data);
	    print_p_hat_ekf(p_file_hat, p_data, p_par, calc[0], p_kalman_update, &Ct.matrix, -1, 0.0);
	}
    }


    //////////////////
    // for all data //
    //////////////////
    for(n=0; n < p_data->nb_obs; n++) {

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


	reset_inc(p_X, p_data);	
	reset_inc_cov(&Ct.matrix);

	// propagate X->proj (containing the covariance Ct) if populations not exploding
	if (get_total_pop(p_X->proj)<WORLD_POP) {
	    f_pred(p_X, t0, t1, p_par, p_data, calc[0]);
	} else {
	    print_err("total_pop(X->proj)>=WORLD_POP");
	}

	proj2obs(p_X, p_data);

	if (print_opt & PLOM_PRINT_X) {
	    print_X(p_file_X, &p_par, &p_X, p_data, calc[0], 1, 1, m, n, t1);
	}

	if(p_data->data_ind[n]->n_nonan){
	    X2xk(p_kalman_update->xk, p_X, p_data);

	    log_lik_temp = 0.0;

	    //Observations are assimilated one by one, which does not
	    //change the outcome but makes the code more robust to varying
	    //numbers of observations.
	    if (print_opt & PLOM_PRINT_PRED_RES) {
		print_prediction_residuals_ekf(p_file_pred_res, p_par, p_data, calc[0], p_X, p_kalman_update, &Ct.matrix, n, t1);
	    }
    
	    for(ts=0; ts< data_ind[n]->n_nonan; ts++) {
		int ts_nonan = data_ind[n]->ind_nonan[ts];
		double xk_t_ts = gsl_vector_get(p_kalman_update->xk, N_PAR_SV*N_CAC + ts_nonan);

		p_kalman_update->sc_rt = obs_var(xk_t_ts, p_par, p_data, calc[0], ts_nonan, m, t1);
	    
		//Observations are assimilated one by one so we reset p_kalman_update->ht to 0.0 at each iteration on ts
		gsl_vector_set_zero(p_kalman_update->ht); 
		eval_ht(p_kalman_update, xk_t_ts, p_par, p_data, calc[0], ts_nonan, n, t1);

		// compute gain
		ekf_gain_computation(p_kalman_update,
				     obs_mean(xk_t_ts, p_par, p_data, calc[0], ts_nonan, n, t1),
				     p_data->data[n][ts_nonan],
				     &Ct.matrix); 

		like = ekf_update(p_kalman_update, &Ct.matrix);
		log_lik_temp += log(like);	  
	    }

	    //echo back the change on xk to p_X->proj
	    xk2X(p_X, p_kalman_update->xk, p_data, calc[0], t1);
	    
	    
	    plom_err_code rc = test_all_sv_pos(p_kalman_update->xk, p_data, calc[0], t1);
	    if(rc != PLOM_SUCCESS){
		print_err("error negative compartment sizes");
		log_lik += LOG_LIKE_MIN;
	    } else {
		log_lik += log_lik_temp;
	    }
	}

	if (print_opt & PLOM_PRINT_HAT) {
	    print_p_hat_ekf(p_file_hat, p_data, p_par, calc[0], p_kalman_update, &Ct.matrix, n, t1);
	}


    } // end of for loop on n

    if (OPTION_PRIOR) {
	double log_prob_prior_value;
	plom_err_code rc = log_prob_prior(&log_prob_prior_value, p_best, p_best->mean, p_best->var, p_data);
#if FLAG_VERBOSE
	if(rc != PLOM_SUCCESS){
	    print_err("error log_prob_prior computation");
	}
#endif
        log_lik += log_prob_prior_value;
    }


    if (OPTION_TRANSF) {
        log_lik += log_transf_correc(p_best->mean, p_best->var, p_data->routers);
    }

    return log_lik;
}
