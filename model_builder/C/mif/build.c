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

struct s_mif *build_mif(json_t *theta, enum plom_implementations implementation,  enum plom_noises_off noises_off, double dt, double eps_abs, double eps_rel, double prop_L_option, int J, int *n_threads)
{
    char str[STR_BUFFSIZE];

    json_t *settings = load_settings(PATH_SETTINGS);

    L = (int) floor(prop_L_option*N_DATA);

    struct s_mif *p;
    p = malloc(sizeof(struct s_mif));
    if(p==NULL) {
        sprintf(str, "Allocation impossible in file :%s line : %d",__FILE__,__LINE__);
        print_err(str);
        exit(EXIT_FAILURE);
    }

    p->p_data = build_data(settings, theta, implementation, noises_off, OPTION_PRIOR, -1); //also build obs2ts
    int size_proj = N_PAR_SV*N_CAC + p->p_data->p_it_only_drift->nbtot + N_TS_INC_UNIQUE;

    //N_DATA is set in build_data
    if (L>N_DATA) {
        sprintf(str, "L > N_DATA (%d > %d). Please choose a L <= %d", L, N_DATA, N_DATA);
        print_err(str);
        exit(EXIT_FAILURE);
    }

    p->p_best = build_best(p->p_data, theta);

    p->J_p_X = build_J_p_X(size_proj, N_TS, p->p_data, dt);
    p->J_p_X_tmp = build_J_p_X(size_proj, N_TS, p->p_data, dt);
    p->J_p_par = build_J_p_par(p->p_data);
    p->p_like = build_likelihood();

    p->calc = build_calc(n_threads, GENERAL_ID, eps_abs, eps_rel, J, size_proj, step_ode, p->p_data, settings);

    json_decref(settings);

    //read only
    p->calc[0]->method_specific_shared_data = gsl_matrix_calloc(p->p_best->n_to_be_estimated, p->p_best->n_to_be_estimated); //used to store the cholesky decomposition

    int nt;
    for(nt=0; nt< *n_threads; nt++){
	p->calc[nt]->method_specific_thread_safe_data = gsl_vector_calloc(p->p_best->n_to_be_estimated);  //used to store temporary vector for MVN
	if(nt>0){
	    p->calc[nt]->method_specific_shared_data = p->calc[0]->method_specific_shared_data;
	}
    }

    /*MIF specific*/

    /* MIF computation variables */
    p->J_theta = init2_gsl_vector_d_set0(J, p->p_data->p_it_all->nbtot);
    p->J_theta_tmp = init2_gsl_vector_d_set0(J, p->p_data->p_it_all->nbtot);

    p->D_theta_bart = init2d_set0(N_DATA+1, p->p_data->p_it_all->nbtot);
    p->D_theta_Vt = init2d_set0(N_DATA+1, p->p_data->p_it_all->nbtot);

    return p;
}


void clean_mif(struct s_mif *p)
{
    int nt;
    for(nt=0; nt< p->calc[0]->n_threads; nt++){
	gsl_vector_free(p->calc[nt]->method_specific_thread_safe_data);
    }
    gsl_matrix_free(p->calc[0]->method_specific_shared_data);

    clean_calc(p->calc, p->p_data);
    clean_best(p->p_best);

    clean_J_p_par(p->J_p_par);
    clean_J_p_X(p->J_p_X);
    clean_J_p_X(p->J_p_X_tmp);
    clean_likelihood(p->p_like);

    clean_data(p->p_data);

    clean2_gsl_vector_d(p->J_theta, J);
    clean2_gsl_vector_d(p->J_theta_tmp, J);

    clean2d(p->D_theta_bart, N_DATA+1);
    clean2d(p->D_theta_Vt, N_DATA+1);

    FREE(p);
}
