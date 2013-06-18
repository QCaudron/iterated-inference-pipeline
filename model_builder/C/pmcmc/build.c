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

#include "pmcmc.h"

struct s_pmcmc *build_pmcmc(json_t *theta, enum plom_implementations implementation, enum plom_noises_off noises_off, json_t *settings, double dt, double eps_abs, double eps_rel, const double freeze_forcing, double a, int m_switch, int m_epsilon, double epsilon_max, int is_smooth, double alpha, int J, int *n_threads, int nb_obs)
{
    char str[STR_BUFFSIZE];

    int nt;

    if (OPTION_PIPELINE) {
        //be sure that J is a multiple of JCHUNK
        int newJ = (int) ceil(((double) J)/ ((double) JCHUNK))*JCHUNK;
        if(newJ != J) {
            snprintf(str, STR_BUFFSIZE, "J (%d) has been set to (%d) to be a multiple of Jchunck (%d)", J, newJ, JCHUNK );
            print_log(str);
            J = newJ;
        }
    }

    struct s_pmcmc *p;
    p = malloc(sizeof(struct s_pmcmc));
    if(p==NULL) {
        sprintf(str, "Allocation impossible in file :%s line : %d",__FILE__,__LINE__);
        print_err(str);
        exit(EXIT_FAILURE);
    }


    p->p_data = build_data(settings, theta, implementation, noises_off, 1, nb_obs, "D");
    p->p_best = build_best(p->p_data, theta);

    int size_proj = N_PAR_SV*N_CAC + p->p_data->p_it_only_drift->nbtot + N_TS_INC_UNIQUE;

    p->D_J_p_X = build_D_J_p_X(size_proj, N_TS, p->p_data, dt);
    p->D_J_p_X_tmp = build_D_J_p_X(size_proj, N_TS, p->p_data, dt);
    p->p_par = build_par(p->p_data);
    p->D_p_hat_new = build_D_p_hat(p->p_data);
    p->D_p_hat_prev = build_D_p_hat(p->p_data);
    p->D_p_hat_best = build_D_p_hat(p->p_data);

    p->p_like = build_likelihood();

    p->calc = build_calc(n_threads, GENERAL_ID, eps_abs, eps_rel, J, size_proj, step_ode, freeze_forcing, -1, p->p_data, settings);

    struct s_mcmc_calc_data *p_mcmc_calc_data = build_mcmc_calc_data(p->p_best, a, m_switch, m_epsilon, epsilon_max, is_smooth, alpha);
    //store the ref for each element of calc
    for (nt=0; nt < *n_threads; nt++) {
        p->calc[nt]->method_specific_shared_data = p_mcmc_calc_data;
    }

    sprintf(str, "Starting Simforence-pmcmc with the following options: i = %d, J = %d, LIKE_MIN = %g, M = %d, N_THREADS = %d SWITCH = %d a = %g", GENERAL_ID, J, LIKE_MIN, M, *n_threads, p_mcmc_calc_data->m_switch, p_mcmc_calc_data->a);
    print_log(str);

    return p;
}


void clean_pmcmc(struct s_pmcmc *p)
{
    clean_best(p->p_best);

    clean_par(p->p_par);

    clean_D_J_p_X(p->D_J_p_X);
    clean_D_J_p_X(p->D_J_p_X_tmp);

    clean_D_p_hat(p->D_p_hat_new, p->p_data);
    clean_D_p_hat(p->D_p_hat_prev, p->p_data);
    clean_D_p_hat(p->D_p_hat_best, p->p_data);

    clean_likelihood(p->p_like);

    clean_mcmc_calc_data((struct s_mcmc_calc_data *) p->calc[0]->method_specific_shared_data);
    clean_calc(p->calc, p->p_data);
    clean_data(p->p_data);

    FREE(p);
}
