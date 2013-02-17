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

struct s_mcmc_calc_data *build_mcmc_calc_data(struct s_best *p_best, const double a, const int m_switch, const int m_epsilon, const double epsilon_max, const int is_smoothed_tunning, const double alpha)
{
    /* which parameters have to be estimated (parameters with jump_size > 0.0) */
    char str[STR_BUFFSIZE];
    int k;

    struct s_mcmc_calc_data *p;
    p = malloc(sizeof(struct s_mcmc_calc_data));
    if(p==NULL) {
	snprintf(str, STR_BUFFSIZE, "Allocation impossible in file :%s line : %d",__FILE__,__LINE__);
        print_err(str);
        exit(EXIT_FAILURE);
    }

    p->n_acceptance_rates = p_best->length;
    p->acceptance_rates = init1d_set0(p->n_acceptance_rates);
    p->smoothed_acceptance_rates = init1d_set0(p->n_acceptance_rates);

    //start with 1.0 (for non fitted parameters it will stay at 1.0 (in a way, non fitted parameters are always accepted)
    for(k=0; k< (p->n_acceptance_rates); k++) {
        p->acceptance_rates[k] = 1.0;
        p->smoothed_acceptance_rates[k] = 1.0;
    }

    p->global_acceptance_rate = 1.0;
    p->smoothed_global_acceptance_rate = 1.0;

    p->has_cycled = 1;
    p->m_full_iteration = 0;
    p->cycle_id = p_best->n_to_be_estimated -1;

    p->epsilon = 1.0;
    p->epsilon_max = epsilon_max;
    p->a = a;

    p->alpha = alpha;
    p->is_smoothed_tunning = is_smoothed_tunning;

    // iteration to swith between initial and empirical covariances
    p->m_switch = m_switch;

    int min_switch = 5*p_best->n_to_be_estimated*p_best->n_to_be_estimated;
    if (m_switch < 0) {
        p->m_switch = min_switch;
    } else if (p->m_switch < min_switch) {
	snprintf(str, STR_BUFFSIZE, "attention: covariance switching iteration (%i) is smaller than proposed one (%i)\n", m_switch, min_switch);
        print_warning(str);
    }

    p->m_epsilon = m_epsilon;

    return (p);
}


void clean_mcmc_calc_data(struct s_mcmc_calc_data *p)
{
    FREE(p->acceptance_rates);
    FREE(p->smoothed_acceptance_rates);
    FREE(p);
}

struct s_pmcmc *build_pmcmc(enum plom_implementations implementation, enum plom_noises_off noises_off, json_t *settings, double dt, double eps_abs, double eps_rel, double a, int m_switch, int m_epsilon, double epsilon_max, int is_smooth, double alpha, int update_covariance, int J, int *n_threads)
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

    json_t *theta = load_json();
    p->p_data = build_data(settings, theta, implementation, noises_off, 1); //also build obs2ts
    p->p_best = build_best(p->p_data, theta, update_covariance);
    json_decref(theta);

    int size_proj = N_PAR_SV*N_CAC + p->p_data->p_it_only_drift->nbtot + N_TS_INC_UNIQUE;

    p->D_J_p_X = build_D_J_p_X(size_proj, N_TS, p->p_data, dt);
    p->D_J_p_X_tmp = build_D_J_p_X(size_proj, N_TS, p->p_data, dt);
    p->p_par = build_par(p->p_data);
    p->D_p_hat_new = build_D_p_hat(p->p_data);
    p->D_p_hat_prev = build_D_p_hat(p->p_data);
    p->D_p_hat_best = build_D_p_hat(p->p_data);

    p->p_like = build_likelihood();

    p->calc = build_calc(n_threads, GENERAL_ID, eps_abs, eps_rel, J, size_proj, step_ode, p->p_data);

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
