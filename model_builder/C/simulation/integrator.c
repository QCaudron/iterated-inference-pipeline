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

#include "simulation.h"

int has_failed(double *y)
{
    int i;
    for (i=0; i<(N_PAR_SV*N_CAC); i++) {
        if ( isnan(y[i]) || isinf(y[i]) ) {
            return 1;
        }
    }
    return 0;
}

int integrator(struct s_X *p_X, double *y0, double t0, double t_end, struct s_par *p_par, double abs_tol, double rel_tol, struct s_calc *p_calc, struct s_data *p_data)
{
    /*numerical integration from t0 to t_end, return 0 if success*/

    int status;
    int error=0;
    double t = t0;
    double h = p_X->dt; //h is the initial integration step size
    int i;
    double *y = p_X->proj;

    /* reset abs_tol and rel_to; */
    gsl_odeiv2_control_free(p_calc->control);
    p_calc->control = gsl_odeiv2_control_y_new(abs_tol, rel_tol); /*abs and rel error (eps_abs et eps_rel) */
    p_calc->p_par = p_par; //pass the ref to p_par so that it is available wihtin the function to integrate

    /* initialize with initial conditions and reset incidence */
    for (i=0; i< (N_PAR_SV*N_CAC) ; i++) {
        y[i]=y0[i];
    }

    /* tentative of numerical integration (with a precision of abs_tol et rel_tol) */
    while (t < t_end) {
        reset_inc(p_X, p_data);
        status = gsl_odeiv2_evolve_apply(p_calc->evolve, p_calc->control, p_calc->step, &(p_calc->sys), &t, t_end, &h, y);

        if ( (status != GSL_SUCCESS) ) { //more stringent: || has_failed(y)
            char str[STR_BUFFSIZE];
            sprintf(str, "integration failure (status: %d has_failed: %d)", status, has_failed(y));
            print_err(str);
            error=1;
            break;
        }
    }

    return error;
}


int integrate(struct s_X *p_X, double *y0, double t0, double t_end, struct s_par *p_par, double *abs_tol, double *rel_tol, struct s_calc *p_calc, struct s_data *p_data)
{
    /* recursive function that decreases abs_tol and rel_tol until integration success */

    int integration_error=1;

    while( integration_error && (*abs_tol>ABS_TOL_MIN) && (*rel_tol>REL_TOL_MIN)  ) {
        integration_error = integrator(p_X, y0, t0, t_end, p_par,  *abs_tol, *rel_tol, p_calc, p_data);
        if(integration_error) {
            *abs_tol /= 10.0;
            *rel_tol /= 10.0;
        }
    }

    return integration_error;
}


/**
 * Used for bifurcation analysis ONLY.
 */
double **get_traj_obs(struct s_X *p_X, double *y0, double t0, double t_end, double t_transiant, struct s_par *p_par, struct s_data *p_data, struct s_calc *p_calc, plom_f_pred_t f_pred)
{
    int i, ts, k;
    double **traj_obs = init2d_set0(N_TS, (int) (t_end-t0));

    FILE *p_file_X = NULL;
    if (OPTION_TRAJ) {
        p_file_X = sfr_fopen(SFR_PATH, GENERAL_ID, "X", "w", header_X, p_data);
    }

    /* initialize with initial conditions and reset incidence */
    for (i=0; i< (N_PAR_SV*N_CAC) ; i++){
        p_X->proj[i]=y0[i];
    }

    for (k= (int) t0 ; k< (int) t_end ; k++) {

        if ( t_transiant <= N_DATA ) {
            p_calc->current_nn= (k < N_DATA_PAR_FIXED) ? k : N_DATA_PAR_FIXED-1;
        }

        reset_inc(p_X, p_data);
        f_pred(p_X, k, k+1, p_par, p_data, p_calc);
        proj2obs(p_X, p_data);

        for (ts=0; ts<N_TS; ts++) {
            traj_obs[ts][k- ((int) t0)] = p_X->obs[ts];
        }

        if (OPTION_TRAJ) {
            print_X(p_file_X, &p_par, &p_X, p_data, p_calc, k+1, 1, 0, 0);
        }
    }

    if (OPTION_TRAJ) {
        sfr_fclose(p_file_X);
    }

    return traj_obs;
}



void traj(struct s_X **J_p_X, double t0, double t_end, double t_transiant, struct s_par **J_p_par, struct s_data *p_data, struct s_calc **calc, plom_f_pred_t f_pred)
{
    int j, k, nn;
    int thread_id;

    FILE *p_file_X = NULL;
    if (OPTION_TRAJ) {
        p_file_X = sfr_fopen(SFR_PATH, GENERAL_ID, "X", "w", header_X, p_data);
    }
    FILE *p_file_hat = sfr_fopen(SFR_PATH, GENERAL_ID, "hat", "w", header_hat, p_data);

    struct s_hat *p_hat = build_hat(p_data);

    //if ODE, only the first particle was used to skip the transiant
    if ( (p_data->implementation == PLOM_ODE) && (t_transiant > 0.0) ) {
        replicate_J_p_X_0(J_p_X, p_data);
    }

    for (k= (int) t0 ; k< (int) t_end ; k++) {

#if FLAG_JSON //for the webApp, we block at every iterations to prevent the client to be saturated with msg
        if (OPTION_TRAJ) {
            if(k % 10 == 0){
                block();
            }
        }
#endif

        if ( t_transiant <= N_DATA ) {
            nn= (k < N_DATA_PAR_FIXED) ? k : N_DATA_PAR_FIXED-1;
            store_state_current_n_nn(calc, nn, nn);
        }

#pragma omp parallel for private(thread_id)
        for(j=0;j<J;j++) {
            thread_id = omp_get_thread_num();
            reset_inc(J_p_X[j], p_data);

            f_pred(J_p_X[j], k, k+1, J_p_par[j], p_data, calc[thread_id]);
            proj2obs(J_p_X[j], p_data);
        }

        compute_hat_nn(J_p_X, J_p_par, p_data, calc, p_hat, 0);
        print_p_hat(p_file_hat, NULL, p_hat, p_data, k);

        if (OPTION_TRAJ && FLAG_JSON==0) {
            print_X(p_file_X, J_p_par, J_p_X, p_data, calc[0], k+1, 0, 0, 0);
        }
    }

    clean_hat(p_hat, p_data);
    if (OPTION_TRAJ) {
        sfr_fclose(p_file_X);
    }
    sfr_fclose(p_file_hat);
}
