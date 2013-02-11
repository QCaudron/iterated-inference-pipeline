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

void clean_kalman_specific_data(struct s_calc *p_calc, struct s_data *p_data)
{
    int i;

    struct s_kalman_specific_data *p =  (struct s_kalman_specific_data *) p_calc->method_specific_thread_safe_data;

    //compo_groups_drift_par_proc
    for(i=0; i<N_DRIFT; i++) {
        clean_groups_compo(p->compo_groups_drift_par_proc[i], p_data->routers[ p_data->drift[i]->ind_par_Xdrift_applied ]->n_gp );
    }
    FREE(p->compo_groups_drift_par_proc);

    gsl_matrix_free(p->FtCt);
    gsl_matrix_free(p->Ft);
    gsl_matrix_free(p->Q);

    FREE(p);
}

struct s_kalman_update *build_kalman_update(int n_kal)
{
    struct s_kalman_update *p;
    p = malloc(sizeof(struct s_kalman_update));
    if(p==NULL) {
        char str[STR_BUFFSIZE];
        sprintf(str, "Allocation impossible in file :%s line : %d",__FILE__,__LINE__);
        print_err(str);
        exit(EXIT_FAILURE);
    }

    p->xk = gsl_vector_calloc(n_kal);
    p->kt = gsl_vector_calloc(n_kal);
    p->ht = gsl_vector_calloc(n_kal);

    p->sc_st = 0.0;
    p->sc_pred_error = 0.0;
    p->sc_rt = 0.0;

    return p;
}

void clean_kalman_update(struct s_kalman_update *p)
{
    gsl_vector_free(p->xk);
    gsl_vector_free(p->kt);
    gsl_vector_free(p->ht);

    FREE(p);
}


struct s_kalman *build_kalman(json_t *settings, enum plom_implementations implementation, enum plom_noises_off noises_off, int is_bayesian, int update_covariance)
{
    char str[STR_BUFFSIZE];

    struct s_kalman *p_kalman;
    p_kalman = malloc(sizeof(struct s_kalman));
    if(p_kalman==NULL) {
        sprintf(str, "Allocation impossible in file :%s line : %d",__FILE__,__LINE__);
        print_err(str);
        exit(EXIT_FAILURE);
    }
    json_t *theta = load_json();

    p_kalman->p_data = build_data(settings, theta, implementation, noises_off, is_bayesian);

    N_KAL = N_PAR_SV*N_CAC + N_TS + p_kalman->p_data->p_it_only_drift->nbtot;
    int size_proj = N_PAR_SV*N_CAC + p_kalman->p_data->p_it_only_drift->nbtot + N_TS_INC_UNIQUE + (N_KAL*N_KAL);

    p_kalman->p_X = build_X(size_proj, N_TS, p_kalman->p_data); //proj contains Ct.
    p_kalman->p_best = build_best(p_kalman->p_data, theta, update_covariance);
    json_decref(theta);

    int n_threads =1;
    p_kalman->calc = build_calc(&n_threads, GENERAL_ID, 0, 1, size_proj, step_ode_ekf, p_kalman->p_data);
    p_kalman->p_par = build_par(p_kalman->p_data);

    p_kalman->smallest_log_like = get_smallest_log_likelihood(p_kalman->p_data->data_ind);
    p_kalman->p_kalman_update = build_kalman_update(N_KAL);

    p_kalman->calc[0]->method_specific_thread_safe_data = build_kalman_specific_data(p_kalman->calc[0], p_kalman->p_data);

    return p_kalman;
}


void clean_kalman(struct s_kalman *p_kalman)
{
    clean_X(p_kalman->p_X);
    clean_best(p_kalman->p_best);
    clean_par(p_kalman->p_par);
    clean_kalman_update(p_kalman->p_kalman_update);

    clean_kalman_specific_data(p_kalman->calc[0], p_kalman->p_data);

    clean_calc(p_kalman->calc, p_kalman->p_data);
    clean_data(p_kalman->p_data);

    FREE(p_kalman);
}
