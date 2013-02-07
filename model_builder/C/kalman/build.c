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

struct s_kalman_specific_data *build_kalman_specific_data(struct s_data *p_data, enum plom_implementations implementation,  enum plom_noises_off noises_off)
{
    int i;

    struct s_kalman_specific_data *p;
    p = malloc(sizeof(struct s_kalman_specific_data));
    if(p==NULL){
        char str[STR_BUFFSIZE];
        snprintf(str, STR_BUFFSIZE, "Allocation impossible in file :%s line : %d",__FILE__,__LINE__);
        print_err(str);
        exit(EXIT_FAILURE);
    }

    //compo_groups_drift_par_proc
    p->compo_groups_drift_par_proc = malloc(N_DRIFT* sizeof (struct s_group **));
    if(p->compo_groups_drift_par_proc==NULL) {
        char str[STR_BUFFSIZE];
        snprintf(str, STR_BUFFSIZE, "Allocation impossible in file :%s line : %d",__FILE__,__LINE__);
        print_err(str);
        exit(EXIT_FAILURE);
    }
    for(i=0; i<N_DRIFT; i++) {
        p->compo_groups_drift_par_proc[i] = get_groups_compo(p_data->routers[ p_data->drift[i]->ind_par_Xdrift_applied ]);
    }

    p->Ft = gsl_matrix_calloc(N_KAL, N_KAL);
    p->FtCt = gsl_matrix_alloc(N_KAL, N_KAL);
    p->Q = gsl_matrix_calloc(N_KAL, N_KAL);

    if ( (noises_off & (PLOM_NO_DEM_STO)) && (noises_off & (PLOM_NO_ENV_STO)) )  {
	p->eval_Q = &eval_Q_no_dem_sto_no_env_sto;
    } else if ((noises_off & PLOM_NO_DEM_STO) && !(noises_off & PLOM_NO_ENV_STO)) {
	p->eval_Q = &eval_Q_no_dem_sto;	
    } else if (!(noises_off & PLOM_NO_DEM_STO) && (noises_off & PLOM_NO_ENV_STO)) {
	p->eval_Q = &eval_Q_no_env_sto;	
    } else {
	p->eval_Q = &eval_Q;	
    }

    return p;
}


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

struct s_kal *build_kal(int n_kal)
{
    struct s_kal *p_kal;
    p_kal = malloc(sizeof(struct s_kal));
    if(p_kal==NULL) {
        char str[STR_BUFFSIZE];
        sprintf(str, "Allocation impossible in file :%s line : %d",__FILE__,__LINE__);
        print_err(str);
        exit(EXIT_FAILURE);
    }

    p_kal->xk = gsl_vector_calloc(n_kal);
    p_kal->kt = gsl_vector_calloc(n_kal);
    p_kal->ht = gsl_vector_calloc(n_kal);

    p_kal->sc_st = 0.0;
    p_kal->sc_pred_error = 0.0;
    p_kal->sc_rt = 0.0;

    return p_kal;
}

void clean_kal(struct s_kal *p_kal)
{
    gsl_vector_free(p_kal->xk);
    gsl_vector_free(p_kal->kt);
    gsl_vector_free(p_kal->ht);

    FREE(p_kal);
}


struct s_kalman *build_kalman(json_t *settings, enum plom_implementations implementation,  enum plom_noises_off noises_off, int *n_threads, int is_bayesian, int update_covariance)
{
    char str[STR_BUFFSIZE];
    int nt;

    struct s_kalman *p_kalman;
    p_kalman = malloc(sizeof(struct s_kalman));
    if(p_kalman==NULL) {
        sprintf(str, "Allocation impossible in file :%s line : %d",__FILE__,__LINE__);
        print_err(str);
        exit(EXIT_FAILURE);
    }
    json_t *theta = load_json();

    p_kalman->p_data = build_data(settings, theta, is_bayesian);

    N_KAL = N_PAR_SV*N_CAC + N_TS + p_kalman->p_data->p_it_only_drift->nbtot;
    int size_proj = N_PAR_SV*N_CAC + p_kalman->p_data->p_it_only_drift->nbtot + N_TS_INC_UNIQUE + (N_KAL*N_KAL);

    p_kalman->p_X = build_X(size_proj, N_TS, p_kalman->p_data); //proj contains Ct.
    p_kalman->p_best = build_best(p_kalman->p_data, theta, noises_off, update_covariance);
    json_decref(theta);

    p_kalman->calc = build_calc(n_threads, GENERAL_ID, implementation, 1, size_proj, func_kal, p_kalman->p_data);
    p_kalman->p_par = build_par(p_kalman->p_data);

    p_kalman->smallest_log_like = get_smallest_log_likelihood(p_kalman->p_data->data_ind);
    p_kalman->p_kal = build_kal(N_KAL);

    for(nt=0; nt< *n_threads; nt++) {
        p_kalman->calc[nt]->method_specific_thread_safe_data = build_kalman_specific_data(p_kalman->p_data, implementation, noises_off);
    }

    return p_kalman;
}


void clean_kalman(struct s_kalman *p_kalman, enum plom_implementations implementation, int n_threads)
{
    int nt;

    clean_X(p_kalman->p_X);
    clean_best(p_kalman->p_best);
    clean_par(p_kalman->p_par);
    clean_kal(p_kalman->p_kal);

    for(nt=0; nt< n_threads; nt++) {
        clean_kalman_specific_data(p_kalman->calc[nt], p_kalman->p_data);
    }

    clean_calc(p_kalman->calc, implementation);
    clean_data(p_kalman->p_data);

    FREE(p_kalman);
}
