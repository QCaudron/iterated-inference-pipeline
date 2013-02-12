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

#include "mif.h"

struct s_mif *build_mif(enum plom_implementations implementation,  enum plom_noises_off noises_off, double dt, double eps_abs, double eps_rel, double prop_L_option, int J, int *n_threads)
{
    char str[STR_BUFFSIZE];

    json_t *settings = load_settings(PATH_SETTINGS);
    json_t *theta = load_json();

    L = (int) floor(prop_L_option*N_DATA);

    struct s_mif *p_mif;
    p_mif = malloc(sizeof(struct s_mif));
    if(p_mif==NULL) {
        sprintf(str, "Allocation impossible in file :%s line : %d",__FILE__,__LINE__);
        print_err(str);
        exit(EXIT_FAILURE);
    }

    p_mif->p_data = build_data(settings, theta, implementation, noises_off, OPTION_PRIOR); //also build obs2ts
    json_decref(settings);
    int size_proj = N_PAR_SV*N_CAC + p_mif->p_data->p_it_only_drift->nbtot + N_TS_INC_UNIQUE;

    //N_DATA_NONAN is set in build_data
    if (L>N_DATA_NONAN) {
        sprintf(str, "L > N_DATA_NONAN (%d > %d). Please choose a L <= %d", L, N_DATA_NONAN, N_DATA_NONAN);
        print_err(str);
        exit(EXIT_FAILURE);
    }

    p_mif->p_best = build_best(p_mif->p_data, theta, 0);
    json_decref(theta);

    p_mif->J_p_X = build_J_p_X(size_proj, N_TS, p_mif->p_data);
    p_mif->J_p_X_tmp = build_J_p_X(size_proj, N_TS, p_mif->p_data);
    p_mif->J_p_par = build_J_p_par(p_mif->p_data);
    p_mif->p_like = build_likelihood();

    p_mif->calc = build_calc(n_threads, GENERAL_ID, dt, eps_abs, eps_rel, J, size_proj, step_ode, p_mif->p_data);

    /*MIF specific*/

    /* MIF global variables: */
    N_THETA_MIF = p_mif->p_data->p_it_par_proc_par_obs_no_drift->nbtot;
    N_IC_MIF = p_mif->p_data->p_it_par_sv_and_drift->nbtot;

    /* MIF computation variables */
    p_mif->var_theta = gsl_matrix_calloc(N_THETA_MIF, N_THETA_MIF);

    p_mif->J_theta = init2_gsl_vector_d_set0(J, N_THETA_MIF);
    p_mif->J_theta_tmp = init2_gsl_vector_d_set0(J, N_THETA_MIF);

    p_mif->J_IC_grouped = init2d_set0(J, N_IC_MIF);
    p_mif->J_IC_grouped_tmp = init2d_set0(J, N_IC_MIF);

    p_mif->D_theta_bart = init2d_set0(N_DATA_NONAN+1, N_THETA_MIF);
    p_mif->D_theta_Vt = init2d_set0(N_DATA_NONAN+1, N_THETA_MIF);

    return p_mif;
}


void clean_mif(struct s_mif *p_mif)
{
    clean_calc(p_mif->calc, p_mif->p_data);
    clean_best(p_mif->p_best);

    clean_J_p_par(p_mif->J_p_par);
    clean_J_p_X(p_mif->J_p_X);
    clean_J_p_X(p_mif->J_p_X_tmp);
    clean_likelihood(p_mif->p_like);

    clean_data(p_mif->p_data);

    gsl_matrix_free(p_mif->var_theta);

    clean2_gsl_vector_d(p_mif->J_theta, J);
    clean2_gsl_vector_d(p_mif->J_theta_tmp, J);

    clean2d(p_mif->J_IC_grouped, J);
    clean2d(p_mif->J_IC_grouped_tmp, J);

    clean2d(p_mif->D_theta_bart, N_DATA_NONAN+1);
    clean2d(p_mif->D_theta_Vt, N_DATA_NONAN+1);

    FREE(p_mif);
}
