// # Extended Kalman Filter
// In C with adaptive-step integration.
// ### Compilation
// We use libplom and libplomtpl (the templates of your specific model)
// ```gcc -fopenmp -std=gnu99 -Wall -O3 -I C/core/ -L C/lib kalman.c -o smc -lplom -lplomtpl -lm -lgsl -lgslcblas -ljansson -lzmq```
// ###Usage
// ```./kalman < theta.json```
// Parameters are read from _stdin_, the log likelihood is returned


#include "plom.h"


double run_kalman(struct s_X *p_X, struct s_best *p_best, struct s_par *p_par, struct s_kal *p_kal, struct s_common *p_common, struct s_data *p_data, struct s_calc **calc, FILE *p_file_X, int m)
{
    // ##Global Variables
    // All global variables are defined in plom.h
    //Typically these values will be prompted from _stdin_ but here we assign them
    J = 100; //Number of particles
    GENERAL_ID = 0; //unique ID to generate a seed and sign the output files
    LIKE_MIN = 1e-17; LOG_LIKE_MIN = log(LIKE_MIN); //minimum value of the likelihood

    
    
    // ##Inputs (JSON)
    // All PLoM inputs are in [JSON](http://www.json.org/). JSON is handled on the C side with the lib [Jansson](http://www.digip.org/jansson/)
    // ###settings
    // settings comes from a JSON **file** (settings.json) containing:
    //
    //  - the data
    //  - all the constant representing the various dimension of your model. Most of this constant are global variables (for instance _N\_PAR\_SV_ is the number of state variables...).
    json_t *settings = load_settings(PATH_SETTINGS); //PATH_SETTINGS is defined in plom.h (defaults to settings/settings.json)
    
    // ###theta
    // theta comes from a JSON document (theta.json) **read from stdin** containing the parameter values and some information on their **transformation** and **grouping**.
    // _load\_json_  will block waiting from input from _stdin_.
    json_t *theta = load_json();
    
    
    // ##Memory allocation
    // We create the main data structure used by the various inference methods
    //
    // - _p\_data_: the **immutable** data. The data itself and all the navigation data required to iterate through the parameters.
    // - _p\_X_: an object representing the state variables, the observed variable and the drifted variables.
    // - _p\_best_: The parameters (in the **transformed** scale) and everything needed to mutate the parameters. _p\_best_ is typically used for inference methods (MIF, pMCMC, simplex...)
    // - _p\_par_: The parameters in the **natural** scale
    // - _p\_kal_: a structure containing the objects needed for the EKF: the Covariance matrix Ct, the Process Jacobian Ft, the gain kt, the prediction error, and other temporary variables that are needed in the update process
    // - _p\_common_: a structure containing the diffusion matrix Q and the jacobian of the observation process, ht.
    // - _calc_: an array containing one _p\_calc_ object per thread. _p\_calc_ objects are to support the calculations in multi-threaded environments. _p\_calc_ objects store everything related to random number generators, numerical integration methods and **mutable** states that need to be computed **in parallel**.
    
    
    struct s_data *p_data = build_data(settings, theta, 0);
    struct s_X *p_X = build_X(p_data);
    struct s_best *p_best = build_best(p_data, theta, 0);
    struct s_par *p_par = build_par(p_data);
    struct s_kal *p_kal = build_kal();
    struct s_common *p_common = build_common();
    struct s_calc **calc = build_calc(GENERAL_ID, N_PAR_SV*N_CAC +N_TS_INC_UNIQUE, func, p_data);
    
    //We are done with the json objects, we free them.
    json_decref(settings);
    json_decref(theta);
    
    
    
    //# The Extended Kalman Filter
    
    // loops indices
    int n, nn;		// data and nonan data indices
    double t0, t1;	// first and last times
    int ts, ts_nonan;	// time series indices
    
    // likelihoods
    double like;
    double log_lik;
    double log_lik_temp;
    
    struct s_data_ind **data_ind = p_data->data_ind;
    
    t0=0;
    log_lik = 0.0;
    
    //For every time step when the data are not all NaN (which is how missing values are represented)
    for(n=0; n<N_DATA_NONAN; n++) {
        
        //##Projection
        //We project \p\_X to the next data point
        
        //The time of the next data point
        t1=p_data->times[n];
        
        
        //Data points might not be equally spaced... We use this sub loop (on _nn_) to simulate equally spaced time step and ensure that incidences are regularly reset to 0.
        for(nn=t0; nn<t1; nn++) {
            store_state_current_n_nn(calc, n, nn);
            
            // At the end of each observation period, the incidence is brought back to 0 as well as the corresponding uncertainty reflected in the covariance matrix Ct
            reset_inc(p_X);	
            reset_inc_Cov(p_kal->Ct); 
            
            // The mean state and covariance are concatenated in a list to be handed to the integrating tool of the gsl library. 
            double xc[N_PAR_SV*N_CAC+N_TS_INC_UNIQUE + N_KAL*(N_KAL+1)/2];
            X2xc(xc, p_X, p_kal->Ct);
            
            //Integration of the process model from nn to nnp1 with the following formulae:
            //
            //dp_X += f(p_X)dt
            //
            //dCt  += [Ft(p_X)*Ct + Ct*t(Ft(p_X)) + Q]dt
            //
            // with f being the deterministic skeleton of the model, Ft its Jacobian, and Q the dispersion matrix
        
            f_prediction_ode_rk(xc, nn, nn+1, p_par, calc[0]);
            
            
            // Recovery of the mean state and covariance from the concatenated list xc
            int i;
            for (i=0; i<N_PAR_SV*N_CAC+N_TS_INC_UNIQUE; i++) {
                p_X->proj[i]  = xc[i];
            }
            list2sym_matrix(p_kal->Ct, xc, i);
            
        } 
        
        //We transform p_X, which was a compact representation of the state, into x_k that accounts for all the observed quantities
        X2xk(p_kal->xk, p_X, p_data);
 
        
        //##Update
        
        p_kal->sc_rt = 0.0;
        log_lik_temp = 0.0;
        //Observations are assimilated one by one, which does not change the outcome but makes the code more robust to varying numbers of observations.
        
        for(ts=0; ts< data_ind[n]->n_nonan; ts++) {
            
            // ###Computation of rt
            // variance of the observation process:
            p_kal->sc_rt = obs_var(gsl_vector_get(p_kal->xk, N_PAR_SV*N_CAC + ts_nonan), p_par, p_data, calc[0], ts_nonan);
            
            // ###Evaluation of ht
            // gradient of sthe observation process:
            eval_ht(p_common->ht, p_kal->xk, p_par, p_data, calc[0], ts_nonan);
            
            // ### Gain Computation
            //
            //   st = ht*Ct*ht + rt
            //
            //   kt = Ct*ht/st
            ekf_gain_computation(obs_mean(gsl_vector_get(p_kal->xk, N_PAR_SV*N_CAC +ts_nonan), p_par, p_data, calc[0], ts_nonan),
                                 p_data->data[calc[0]->current_nn][ts_nonan],
                                 p_kal->Ct, p_common->ht, p_kal->kt, p_kal->sc_rt,
                                 &(p_kal->sc_st), &(p_kal->sc_pred_error));
            
            // ### Update of xk, Ct, and the likelihood
            //
            // xk = xk + kt*(pred_error)
            //
            // Ct = (Id - kt*ht)*Ct
            //
            // lik *= N(pred_error,st)
            like = ekf_update(p_kal->xk, p_kal->Ct, p_common->ht, p_kal->kt, p_kal->sc_st, p_kal->sc_pred_error);
            log_lik_temp += log(like);
        }
        
        //Echo back the change on xk to p_X->proj and p_X->drift
        xk2X(p_X, p_kal->xk, p_data);
        
        log_lik += log_lik_temp;
        t0=t1;
        
    } // end of for loop on n
    
    
    //##Result (output)
    return log_lik;
}

