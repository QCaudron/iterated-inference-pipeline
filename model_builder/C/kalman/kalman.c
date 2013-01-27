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
        return jac_tpl / r->f_derivative(jac_der, r->min_z[g], r->max_z[g]);
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
    struct s_router **routers = p_data->routers;
    struct s_iterator *p_it = p_data->p_it_only_drift;

    // X->proj
    int i, k, offset;
    for(i=0; i<N_PAR_SV*N_CAC; i++)
        gsl_vector_set(xk, i, p_X->proj[i]);

    // X->obs
    for(i=0; i<N_TS; i++)
        gsl_vector_set(xk, N_PAR_SV*N_CAC + i, p_X->obs[i]);

    // X->drift
    offset = 0;
    for(i=0; i<p_it->length; i++) {
        for(k=0; k< routers[ p_it->ind[i] ]->n_gp; k++) {
            gsl_vector_set(xk, N_PAR_SV*N_CAC + N_TS + offset, p_X->drift[i][k]);
            offset++;
        }
    }
}

/**
 * reconstruct s_X from xk and set to 0.0 terms of proj that became negative due to Kalman.
 * Note that only N_PAR_SV*N_CAC of proj and drift are necessary. All the observed variable are useless here
 * @param p_X the X vector to be constructed
 * @param xk the concatenated vector
 */
void xk2X(struct s_X *p_X, gsl_vector *xk, struct s_data *p_data)
{
    struct s_router **routers = p_data->routers;
    struct s_iterator *p_it = p_data->p_it_only_drift;

    int i, k, offset;

    // X->proj
    for(i=0; i<N_PAR_SV*N_CAC; i++)
        p_X->proj[i] = (gsl_vector_get(xk, i) > 0.0) ? gsl_vector_get(xk, i) : 0.0 ;

    // X->drift
    offset = 0;
    for(i=0; i<p_it->length; i++) {
        for(k=0; k< routers[ p_it->ind[i] ]->n_gp; k++) {
            p_X->drift[i][k] = gsl_vector_get(xk, N_PAR_SV*N_CAC + N_TS + offset) ;
            offset++;
        }
    }
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

                p_tmp = 1./((*(routers[i]->f_inv_derivative))(gsl_vector_get(mean, offset), routers[i]->min[k], routers[i]->max[k]));

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
 *   run an extended Kalman filter and returns the log likelihood
 */
double run_kalman(struct s_X *p_X, struct s_best *p_best, struct s_par *p_par, struct s_kal *p_kal, struct s_data *p_data, struct s_calc **calc, FILE *p_file_X, int m)
{
    // loops indices
    int n, nn;		// data and nonan data indices
    double t0, t1;	// first and last times
    int ts, ts_nonan;	// time series indices

    // likelihoods
    double like, log_lik, log_lik_temp;

    struct s_data_ind **data_ind = p_data->data_ind;

    t0=0;
    log_lik = 0.0;
    p_kal->sc_st = 0.0;
    p_kal->sc_pred_error = 0.0;

    gsl_matrix_view Ct = gsl_matrix_view_array(&(p_X->proj[PLOM_SIZE_PROJ]), N_KAL, N_KAL);
    gsl_matrix_set_zero(&Ct.matrix);

    //////////////////
    // for all data //
    //////////////////
    for(n=0; n<N_DATA_NONAN; n++) {

#if FLAG_JSON //for the webApp, we block at every iterations to prevent the client to be saturated with msg
        if (OPTION_TRAJ) {
            if(n % 10 == 0){
                block();
            }
        }
#endif

        t1=p_data->times[n];

        //drifted parameters are stored in p_calc->natural_drifted_safe. We initialize the values from p_X->drift[i][k]
        drift_par(calc[0], p_data, p_X, 0, p_data->p_it_only_drift->length);

        /////////////////////////
        // for every time unit //
        /////////////////////////
        /*
         * we have to use this subloop to mimate equaly spaced time step
         * and hence set the incidence to 0 every time unit...
         */
        for(nn=t0; nn<t1; nn++) {
            store_state_current_n_nn(calc, n, nn);

            reset_inc(p_X);	// reset incidence to 0
            reset_inc_cov(&Ct.matrix);	// reset incidence covariance to 0

            // propagate X->proj (containing the covariance Ct) if populations not exploding
            if (get_total_pop(p_X->proj)<WORLD_POP) {
                f_prediction_ode_rk(p_X->proj, nn, nn+1, p_par, calc[0]);
            } else {
                print_err("total_pop(X->proj)>=WORLD_POP");
            }

            proj2obs(p_X, p_data);

            if (OPTION_TRAJ) {
                print_X(p_file_X, &p_par, &p_X, p_data, calc[0], nn+1, 1, 1, m);
            }
        } // end of for loop on nn

        // transform p_X to x_k
        X2xk(p_kal->xk, p_X, p_data);

        // from here we work only with xk

        p_kal->sc_rt = 0.0;
        log_lik_temp = 0.0;
        for(ts=0; ts< data_ind[n]->n_nonan; ts++) {

            ts_nonan = data_ind[n]->ind_nonan[ts];
            p_kal->sc_rt = obs_var(gsl_vector_get(p_kal->xk, N_PAR_SV*N_CAC + ts_nonan), p_par, p_data, calc[0], ts_nonan);

            eval_ht(p_kal->ht, p_kal->xk, p_par, p_data, calc[0], ts_nonan);

            // compute gain
            ekf_gain_computation(obs_mean(gsl_vector_get(p_kal->xk, N_PAR_SV*N_CAC +ts_nonan), p_par, p_data, calc[0], ts_nonan),
                                 p_data->data[calc[0]->current_nn][ts_nonan],
                                 &Ct.matrix, p_kal->ht, p_kal->kt, p_kal->sc_rt,
                                 &(p_kal->sc_st), &(p_kal->sc_pred_error)); //scalar sc_st and sc_pred_error will be modified so we pass their address

            like = ekf_update(p_kal->xk, &Ct.matrix, p_kal->ht, p_kal->kt, p_kal->sc_st, p_kal->sc_pred_error);
            log_lik_temp += log(like);
        }

        //echo back the change on xk to p_X->proj and p_X->drift
        xk2X(p_X, p_kal->xk, p_data);

        log_lik += log_lik_temp;
        t0=t1;

    } // end of for loop on n

    if (OPTION_PRIOR) {
        log_lik += log_prob_prior(p_best, p_best->mean, p_data);
    }

    if (OPTION_TRANSF) {
        log_lik += log_transf_correc(p_best->mean, p_best->var, p_data->routers);
    }

    return log_lik;
}
