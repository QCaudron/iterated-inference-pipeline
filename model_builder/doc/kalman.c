// # Extended Kalman Filter
// In C with adaptive-step integration.


struct s_kal
{
    gsl_vector *xk;  // concatenation of non-overlapping components of X->proj, X->obs and X->drift
    gsl_vector *kt;  // Kalman Gain vector
    gsl_vector *ht;  // Gradient of the observation function

    double sc_st;         // Innovation or residual covariance
    double sc_pred_error; // Innovation or measurement residual
    double sc_rt,         //observation process variance
};



double run_kalman(struct s_X *p_X, struct s_best *p_best, struct s_par *p_par, struct s_kal *p_kal, struct s_data *p_data, struct s_calc **calc)
{
    // loops indices
    int n, nn;		// data and nonan data indices
    double t0, t1;	// first and last times
    int ts, ts_nonan;	// time series indices

    // likelihoods
    double like, log_lik, log_lik_temp;

    // (re)set to 0
    t0=0;
    log_lik = 0.0;
    p_kal->sc_st = 0.0;
    p_kal->sc_pred_error = 0.0;

    //the Covariance matrix Ct
    gsl_matrix_view Ct = gsl_matrix_view_array(&(p_X->proj[PLOM_SIZE_PROJ]), N_KAL, N_KAL);
    gsl_matrix_set_zero(&Ct.matrix);

    //For every time step when the data are not all NaN (which is how missing values are represented)
    for(n=0; n<N_DATA_NONAN; n++) {

        //##Projection
        //We project p\_X to the next data point

        //The time of the next data point
        t1=p_data->times[n];

        //To be hidden... drifted parameters are stored in p\_calc->natural\_drifted\_safe. We initialize the values from p_X->drift[i][k]
        drift_par(calc[0], p_data, p_X, 0, p_data->p_it_only_drift->length);

        //Data points might not be equally spaced... We use this sub loop (on _nn_) to simulate equally spaced time step and ensure that incidences are regularly reset to 0.
        for(nn=t0; nn<t1; nn++) {
            store_state_current_n_nn(calc, n, nn);


            // At the end of each observation period, the incidence is brought back to 0 as well as the corresponding uncertainty reflected in the covariance matrix Ct
            reset_inc(p_X);
            reset_inc_cov(&Ct.matrix);

            //Integration of the process model from nn to nn+1 with the following formulae:
            //
            //dp_X += f(p\_X)dt
            //
            //dCt  += [Ft(p\_X)*Ct + Ct*t(Ft(p\_X)) + Q]dt
            //
            // with _f_ being the deterministic skeleton of the model, _Ft_ its Jacobian, and _Q_ the dispersion matrix

            f_prediction_ode_rk(p_X->proj, nn, nn+1, p_par, calc[0]);
            proj2obs(p_X, p_data);

        }

        //We transform _p\_X_, which was a compact representation of the state, into _x\_k_ that accounts for all the observed quantities
        X2xk(p_kal->xk, p_X, p_data);


        //##Update

        log_lik_temp = 0.0;

        //Observations are assimilated one by one, which does not change the outcome but makes the code more robust to varying numbers of observations.
        for(ts=0; ts< p_data->data_ind[n]->n_nonan; ts++) {

            ts_nonan = p_data->data_ind[n]->ind_nonan[ts];

            // ###Computation of _sc\_rt_
            // variance of the observation process:
            p_kal->sc_rt = obs_var(gsl_vector_get(p_kal->xk, N_PAR_SV*N_CAC + ts_nonan), p_par, p_data, calc[0], ts_nonan);

            // ###Evaluation of ht
            // gradient of the observation process:
            eval_ht(p_kal->ht, p_kal->xk, p_par, p_data, calc[0], ts_nonan);

            // ### Gain Computation
            //
            // sc\_pred\_error = data\_t\_ts - xk\_t\_ts (estimate at time _t_ of the observed quantity)
            //
            // sc\_st = ht' * Ct * ht + sc\_rt
            //
            // kt = Ct * ht' * sc\_st^-1
            ekf_gain_computation(obs_mean(gsl_vector_get(p_kal->xk, N_PAR_SV*N_CAC +ts_nonan), p_par, p_data, calc[0], ts_nonan), //xk_t_ts
                                 p_data->data[calc[0]->current_nn][ts_nonan], // data_t_ts
                                 &Ct.matrix, p_kal->ht, p_kal->kt, p_kal->sc_rt,
                                 &(p_kal->sc_st), &(p_kal->sc_pred_error)); // scalar sc_st and sc_pred_error will be modified so we pass their address


            // ### Update of xk, Ct, and the likelihood
            //
            // xk = xk + kt*(pred_error)
            //
            // Ct = (Id - kt*ht)*Ct
            //
            // lik *= N(sc\_pred\_error, sqrt(sc\_st))
            like = ekf_update(p_kal->xk, &Ct.matrix, p_kal->ht, p_kal->kt, p_kal->sc_st, p_kal->sc_pred_error);
            log_lik_temp += log(like);
        }

        //Echo back the change on _xk_ to _p\_X->proj_ and _p\_X->drift_
        xk2X(p_X, p_kal->xk, p_data);

        log_lik += log_lik_temp;
        t0=t1;

    } // end of for loop on n

    return log_lik;
}
