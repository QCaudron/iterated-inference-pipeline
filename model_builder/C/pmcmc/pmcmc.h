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

//#define FLAG_TRAJ 1 set at compilation (if needed)

int M; /*Number of pMCMC iterations*/
int JCHUNK;
int OPTION_PIPELINE; //0: false, 1: true
int OPTION_FULL_UPDATE;

struct s_mcmc_calc_data
{
    double epsilon;  /**< epsilon factor */
    double epsilon_max; /**< max value for epsilon */
    double a;        /**< cooling factor */
    int m_switch;    /**< number of iterations using empirical covariance */
    int m_epsilon;       /**< number of iterations before tuning epsilon */


    double global_acceptance_rate; /**< the global acceptance rate */
    double smoothed_global_acceptance_rate; /**< as computed with exponential smoothing (http://en.wikipedia.org/wiki/Exponential_smoothing) */

    double alpha; /**< smoothing factor (The term smoothing factor is
		     something of a misnomer, as larger values of
		     alpha actually reduce the level of smoothing, and
		     in the limiting case with alpha = 1 the output
		     series is just the same as the original series
		     (with lag of one time unit). */

    int is_smoothed_tunning; /**<boolean: do we tune epsilon with the
				value of the acceptance rate obtained
				with exponential smoothing ? (1 yes 0
				no) */

    int n_acceptance_rates; /**< s_best->length */
    double *acceptance_rates; /**< [ self.n_acceptance_rates ] parameter specific acceptance rates */
    double *smoothed_acceptance_rates; /**< [ self.n_acceptance_rates ] parameter specific acceptance rates computed with exponential smoothing */

    //counters: note in the absence of the webApp we could avoid this and use the modulo operator. However when users change in real time self.n_to_be_estimated we need to resort on those counters.
    int has_cycled; /**< boolean (have we iterated on all the component self.n_to_be_estimated component of theta) */
    int m_full_iteration; /**< number of full iterations (one full iteration every self.n_to_be_estimated sub-iterations) */
    int cycle_id; /**< position in the sub loop */
};

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
gsl_matrix * propose_new_theta_and_load_X0(double *sd_fac, struct s_best *p_best, struct s_X *p_X, struct s_par *p_par, struct s_data *p_data, struct s_mcmc_calc_data *p_mcmc_calc_data, struct s_calc *p_calc, int m);

void increment_iteration_counters(struct s_mcmc_calc_data *p_mcmc_calc_data, struct s_best *p_best, const int OPTION_FULL_UPDATE);

void ran_proposal_sequential(gsl_vector *proposed, struct s_best *p_best, gsl_matrix *var, double sd_fac, struct s_calc *p_calc);

void compute_best_traj(struct s_hat **D_p_hat_best, struct s_hat **D_p_hat_prev, struct s_hat **D_p_hat_new, struct s_data *p_data, double alpha, double m);

void header_acceptance_rates(FILE *p_file, struct s_data *p_data);
void print_acceptance_rates(FILE *p_file, struct s_mcmc_calc_data *p, int m_full_iteration);

void compute_acceptance_rates(struct s_best *p_best, struct s_mcmc_calc_data *p, double is_accepted, int m);

void print_covariance(FILE *p_file_cov, gsl_matrix *covariance);

/* build.c */
struct s_mcmc_calc_data *build_mcmc_calc_data(struct s_best *p_best, const double a, const int m_switch, const int m_epsilon, const double epsilon_max, const int is_smoothed_tunning, const double alpha);
void clean_mcmc_calc_data(struct s_mcmc_calc_data *p_mcmc_calc_data);

struct s_pmcmc *build_pmcmc(enum plom_implementations implementation, enum plom_noises_off noises_off, json_t *settings, double dt, double eps_abs, double eps_rel, double a, int m_switch, int m_epsilon, double epsilon_max, int is_smooth, double alpha, int update_covariance, int J, int *n_threads);
void clean_pmcmc(struct s_pmcmc *p_pmcmc);
