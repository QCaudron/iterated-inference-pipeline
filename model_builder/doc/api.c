// # Sequential Monte Carlo
// In C and multi-threaded.
// ### Compilation
// We use libplom and libplomtpl (the templates of your specific model)
// ```gcc -fopenmp -std=gnu99 -Wall -O3 -I C/core/ -L C/lib api.c -o smc -lplom -lplomtpl -lm -lgsl -lgslcblas -ljansson -lzmq```
// ###Usage
// ```./smc < theta.json```
// Parameters are read from _stdin_, the log likelihood is printed on _stdout_


#include "plom.h"

int main(void)
{
    // ##Global Variables
    // All global variables are defined in plom.h
    //Typically these values will be prompted from _stdin_ but here we assign them
    J = 100; //Number of particles
    GENERAL_ID = 0; //unique ID to generate a seed and sign the output files
    LIKE_MIN = 1e-17; LOG_LIKE_MIN = log(LIKE_MIN); //minimum value of the likelihood
    N_THREADS = omp_get_max_threads(); omp_set_num_threads(N_THREADS); //the number of thread (set to the number of core of the computer


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
    // - _calc_: an array containing one _p\_calc_ object per thread. _p\_calc_ objects are to support the calculations in multi-threaded environments. _p\_calc_ objects store everything related to random number generators, numerical integration methods and **mutable** states that need to be computed **in parallel**.
    // - _p\_par_: The parameters in the **natural** scale
    // - _J\_p\_X_: an Array of J  _p_X_ objects representing the state variables, the observed variable and the drifted variables.
    // - _J\_p\_X\_tmp_: A replication of _J\_p\_X_ that will be used for the resampling step of the particle filter
    // - _p\_best_: The parameters (in the **transformed** scale) and everything needed to mutate the parameters. _p\_best_ is typically used for inference methods (MIF, pMCMC, simplex...)
    // - _p\_like_: Everything related to the likelihood (weights, ESS...)

    struct s_data *p_data = build_data(settings, theta, 0);
    struct s_calc **calc = build_calc(GENERAL_ID, N_PAR_SV*N_CAC +N_TS_INC_UNIQUE, func, p_data);
    struct s_par *p_par = build_par(p_data);
    struct s_X **J_p_X = build_J_p_X(p_data);
    struct s_X **J_p_X_tmp = build_J_p_X(p_data);
    struct s_best *p_best = build_best(p_data, theta, 0);
    struct s_likelihood *p_like = build_likelihood();

    //We are done with the json objects, we free them.
    json_decref(settings);
    json_decref(theta);

    // ##Transformation and initialization
    //
    // - transform the parameters (log, logit...) as specified in the properties "_transformation_" of _theta.json_ and scale them to the time unit of the data
    transform_theta(p_best, NULL, NULL, p_data, 1, 1);
    // - load p\_par with the untransformed parameters (untransforming p\_best->mean). The back transformation also ensure that the parameter are in the same unit as the data.
    back_transform_theta2par(p_par, p_best->mean, p_data->p_it_all, p_data);

    // - load the states: p_par contain the initial conditions (as proportions), we copy them to the first element of J\_p\_X (J\_p\_X[0])
    linearize_and_repeat(J_p_X[0], p_par, p_data, p_data->p_it_par_sv);
    // - initial condition are proportions, we convert proportion into numbers
    prop2Xpop_size(J_p_X[0], p_data, 1);
    // - Some parameter might follow a diffusion. Diffusions add states to p\_X, that need to be initialized.
    theta_driftIC2Xdrift(J_p_X[0], p_best->mean, p_data);
    //
    //We are done with the first particle, we replicate it for the J-1 other particles.
    replicate_J_p_X_0(J_p_X, p_data);


    //# The particle filter

    int j, n, nn, nnp1, t0, t1, thread_id;
    t0=0;
    //For every time step when the data are not all NaN (which is how missing values are represented)
    for(n=0; n<N_DATA_NONAN; n++) {
        //##Projection
        //We project J\_p\_X to the next data point

        //The time of the next data point
        t1=p_data->times[n];

        //Data points might not be equally spaced... We use this sub loop (on _nn_) to simulate equally spaced time step and ensure that incidences are regularly reset to 0.
        for(nn=t0 ; nn<t1 ; nn++) {
            store_state_current_n_nn(calc, n, nn);
            nnp1 = nn+1;

            //###Numerical integration
            //In this case, we use openMP to run the integration for the J particles in parallel.
#pragma omp parallel for private(thread_id)
            for(j=0; j<J ; j++) {
                //_thread\_id_ will ensure that particles processed  concurrently (_N\_TREADS_) use different _calc_ elements
                thread_id = omp_get_thread_num();
                //reset the incidence to 0
                reset_inc(J_p_X[j]);
                //integration of the process model from nn to nnp1.
                f_prediction_with_drift_sto(J_p_X[j], nn, nnp1, p_par, p_data, calc[thread_id]);

                //###Likelihood computation
                //When we have a data point:
                if(nnp1 == t1) {

                    //- If the observation process parameter contain diffusion,
                    //  we drift the state variable associated to the diffusion (Euler Maruyama)
                    //  and we apply the drift (computed on a state variable) to the corresponding observation process parameter
                    if(N_DRIFT_PAR_OBS) {
                        compute_drift(J_p_X[j], p_par, p_data, calc[thread_id], N_DRIFT_PAR_PROC, N_DRIFT_PAR_PROC + N_DRIFT_PAR_OBS, t1-t0);
                        drift_par(calc[thread_id], p_par, p_data, J_p_X[j], N_DRIFT_PAR_PROC, N_DRIFT_PAR_PROC + N_DRIFT_PAR_OBS);
                    }

                    // - compute the observed  variables
                    proj2obs(J_p_X[j], p_data);

                    // - compute the likelihood values and store it in what will become the weight of the particle
                    p_like->weights[j] = exp(get_log_likelihood(J_p_X[j], p_par, p_data, calc[thread_id]));
                }
            }
        }

        //##Filtering
        //if all the particles have a likelihood < LIKE_MIN, we keep them all and skip the resampling sampling step
        if (weight(p_like, n)) {
            // We use systematic sampling as it is fast and leads to fewer Monte Carlo variability than multinomial sampling.
            systematic_sampling(p_like, calc[0], n);
            //Note J\_p\_X and J\_p\_X_tmp will be swapped (to avoid unnecessary copy) so we pass their **address**
            resample_X(p_like->select[n], &J_p_X, &J_p_X_tmp, p_data);
        }

        t0=t1;
    }


    //##Result (output)
    printf("log likelihood: %g\n", p_like->Llike_best);

    //#End
    //We free all the allocated objects.
    clean_calc(calc);
    clean_J_p_X(J_p_X);
    clean_J_p_X(J_p_X_tmp);
    clean_best(p_best);
    clean_par(p_par);
    clean_likelihood(p_like);
    clean_data(p_data);

    return 0;
}
