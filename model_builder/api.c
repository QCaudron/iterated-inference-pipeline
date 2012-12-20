#include "plom.h"

int main(int argc, char *argv[])
{
    J = 1; N_THREADS =2; GENERAL_ID = 0;
    omp_set_num_threads(N_THREADS);

    print_log("memory allocation and inputs loading...");
    json_t *settings = load_settings(PATH_SETTINGS);
    json_t *theta = load_json();

    struct s_data *p_data = build_data(settings, theta, 0);
    struct s_calc **calc = build_calc(GENERAL_ID, N_PAR_SV*N_CAC +N_TS_INC_UNIQUE, func, p_data);
    struct s_par *p_par = build_par(p_data);
    struct s_X **J_p_X = build_J_p_X(p_data);
    struct s_X **J_p_X_tmp = build_J_p_X(p_data);
    struct s_best *p_best = build_best(p_data, theta, 0);
    struct s_likelihood *p_like = build_likelihood();

    json_decref(settings);
    json_decref(theta);

    print_log("starting computations...");

    transform_theta(p_best, NULL, NULL, p_data, 1, 1);
    back_transform_theta2par(p_par, p_best->mean, p_data->p_it_all, p_data);
    linearize_and_repeat(J_p_X[0], p_par, p_data, p_data->p_it_par_sv);
    prop2Xpop_size(J_p_X[0], p_data, 1);
    theta_driftIC2Xdrift(J_p_X[0], p_best->mean, p_data);

    //load X_0 for the J-1 other particles
    replicate_J_p_X_0(J_p_X, p_data);

    int j, n, nn, nnp1, t0, t1, thread_id;
    t0=0;
    for(n=0; n<N_DATA_NONAN; n++) {
        t1=p_data->times[n];

        /*we have to use this subloop to mimate equaly spaced time step and hence set the incidence to 0 every time unit...*/
        for(nn=t0 ; nn<t1 ; nn++) {
            store_state_current_n_nn(calc, n, nn);
            nnp1 = nn+1;

#pragma omp parallel for private(thread_id)
            for(j=0; j<J ; j++) {
                thread_id = omp_get_thread_num();
                reset_inc(J_p_X[j]);
                f_prediction_with_drift_sto(J_p_X[j], nn, nnp1, p_par, p_data, calc[thread_id]);

                if(N_DRIFT_PAR_OBS) {
                    compute_drift(J_p_X[j], p_par, p_data, calc[thread_id], N_DRIFT_PAR_PROC, N_DRIFT_PAR_PROC + N_DRIFT_PAR_OBS, 1.0); //1.0 is delta_t (nnp1-nn)
                }

                if(nnp1 == t1) {
                    proj2obs(J_p_X[j], p_data);
                    drift_par(calc[thread_id], p_par, p_data, J_p_X[j], N_DRIFT_PAR_PROC, N_DRIFT_PAR_PROC + N_DRIFT_PAR_OBS);
                    p_like->weights[j] = exp(get_log_likelihood(J_p_X[j], p_par, p_data, calc[thread_id]));
                }
            }

        } /* end for on nn */

        if (weight(p_like, n)) {
            systematic_sampling(p_like, calc[0], n);
        }
        resample_X(p_like->select[n], &J_p_X, &J_p_X_tmp, p_data);
        t0=t1;
    } /*end of for loop on the time (n)*/

    print_log("clean up...");
    clean_calc(calc);
    clean_J_p_X(J_p_X);
    clean_J_p_X(J_p_X_tmp);
    clean_best(p_best);
    clean_par(p_par);
    clean_likelihood(p_like);
    clean_data(p_data);

    return 0;
}
