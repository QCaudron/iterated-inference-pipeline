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

#include "plom.h"
#include "mcmc_util.h"

//#define FLAG_TRAJ 1 set at compilation (if needed)

int M; /*Number of pMCMC iterations*/
int JCHUNK;
int OPTION_PIPELINE; //0: false, 1: true


struct s_pmcmc
{
    /* from simforence core */
    struct s_data *p_data;
    struct s_calc **calc;
    struct s_best *p_best;

    struct s_X ***D_J_p_X; /* [N_DATA+1][J] +1 is for initial condition (one time step before first data) */
    struct s_X ***D_J_p_X_tmp; /* [N_DATA+1][J] */
    struct s_par *p_par;

    struct s_hat **D_p_hat_prev;
    struct s_hat **D_p_hat_new;
    struct s_hat **D_p_hat_best;

    struct s_likelihood *p_like;
};


/* pmcmc.c */
void run_propag(struct s_X ***D_J_p_X, struct s_X ***D_J_p_X_tmp, struct s_par *p_par, struct s_hat ***D_p_hat_new,
                struct s_likelihood *p_like, struct s_data *p_data, struct s_calc **calc, plom_f_pred_t f_pred,
                void *sender, void *receiver, void *controller);

void pmcmc(struct s_best *p_best, struct s_X ***D_J_p_X, struct s_X ***D_J_p_X_tmp, struct s_par *p_par, struct s_hat ***D_p_hat_prev, struct s_hat ***D_p_hat_new, struct s_hat **D_p_hat_best, struct s_likelihood *p_like, struct s_data *p_data, struct s_calc **calc, plom_f_pred_t f_pred, int OPTION_ACC);

/* methods.c */
void compute_best_traj(struct s_hat **D_p_hat_best, struct s_hat **D_p_hat_prev, struct s_hat **D_p_hat_new, struct s_data *p_data, double alpha, double m);

/* build.c */
struct s_pmcmc *build_pmcmc(enum plom_implementations implementation, enum plom_noises_off noises_off, json_t *settings, double dt, double eps_abs, double eps_rel, double a, int m_switch, int m_epsilon, double epsilon_max, int is_smooth, double alpha, int update_covariance, int J, int *n_threads);
void clean_pmcmc(struct s_pmcmc *p_pmcmc);
