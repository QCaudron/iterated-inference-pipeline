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

#define WORLD_POP 7.e9 // approx world population, used to detect divergences

/* global variables (used in case of simplex kalman) */
int M;
double CONVERGENCE_STOP_SIMPLEX;

// sizes
int N_KAL;  // N_PAR_SV*N_CAC + N_TS + iterator_only_drift.nbtot

// options
int OPTION_TRANSF;  // add log_transf_correc to log_lik


/*
   To make things easier, vectors are in lower case letters, matrix in
   capital letters and scalar prefixed by sc_
*/


// method_specific_data (accessible via p_calc->method_specific_thread_safe_data) (after casting)
struct s_kalman_specific_data
{
    gsl_matrix *Ft; /**< Jacobian matrix of the drift f(X_t,\theta) of the SDE approximated with the EKF: dX_t = f(X_t,\theta)dt + L b_t */

    struct s_group ***compo_groups_drift_par_proc;

    gsl_matrix *Q;  /**< result of L Qc L' : Dispersion matrix of the
		       SDE approximated with by the EKF: dX_t =
		       f(X_t,\theta)dt + L b_t (b_t is a Browninan
		       motion with diffusion matrix Qc). L is the
		       dispersion matrix as found in Sarkka phD thesis
		       (2006) : as many rows as state variables
		       (including observed variables), as many columns
		       as independent noises (Brownian motion) */

    gsl_matrix *FtCt;	/**< for Ft*Ct (product of the jacobian matrix and the covariance matrix) */

    void (*eval_Q) (gsl_matrix *Q, const double *X, struct s_par *p_par, struct s_data *p_data, struct s_calc *p_calc, struct s_kalman_specific_data *p_kalman_specific_data, double t);
};




struct s_kalman
{
  /* from plom core */
  struct s_data *p_data;
  struct s_calc **calc;
  struct s_best *p_best;
  struct s_par *p_par;
  struct s_X *p_X;

  double smallest_log_like; /* used in simplex kalman (see simplex.h) */

  /* kalman specific */
    struct s_kalman_update *p_kalman_update;
};

/* print.c */
void header_prediction_residuals_ekf(FILE *p_file, struct s_data *p_data);
void print_p_hat_ekf(FILE *p_file, struct s_data *p_data, struct s_kalman_update *p, gsl_matrix *Ct, int n);
void print_prediction_residuals_ekf(FILE *p_file_pred_res, struct s_par *p_par, struct s_data *p_data, struct s_calc *p_calc, struct s_X *p_X, struct s_kalman_update *p, gsl_matrix *Ct, int time);


/* build.c */
struct s_kalman_specific_data *build_kalman_specific_data(struct s_calc *p_calc, struct s_data *p_data);
void clean_kalman_specific_data(struct s_calc *p_calc, struct s_data *p_data);
struct s_kalman_update *build_kalman_update(int n_kalman_update);
void clean_kalman_update(struct s_kalman_update *p_kalman_update);

struct s_kalman *build_kalman(json_t *settings, enum plom_implementations implementation, enum plom_noises_off noises_off, int is_bayesian, int update_covariance, double dt, double eps_abs, double eps_rel);
void clean_kalman(struct s_kalman *p_kalman);

/* kalman.c */
double drift_derivative(double jac_tpl, double jac_der, struct s_router *r, int cac);
double f_simplex_kalman(const gsl_vector *x, void *params);
void xk2X(struct s_X *p_X, gsl_vector *xk, struct s_data *p_data);
void X2xk(gsl_vector *xk, struct s_X *p_X, struct s_data *p_data);
double get_total_pop(double *X);
double log_transf_correc(gsl_vector *mean, gsl_matrix *var, struct s_router **routers);
void reset_inc_cov(gsl_matrix *Ct);
double run_kalman(struct s_X *p_X, struct s_best *p_best, struct s_par *p_par, struct s_kalman_update *p_kalman_update, struct s_data *p_data, struct s_calc **calc, plom_f_pred_t f_pred, FILE *p_file_X, int m, FILE *p_file_pred_res);

/* ekf.c */
void check_and_correct_Ct(gsl_matrix *Ct, struct s_kalman_update *p);
void ekf_propag_cov(double *proj, gsl_matrix *Ft, gsl_matrix *Ct, gsl_matrix *Q, struct s_par *p_par, struct s_group ***compo_groups_drift_par_proc, double t);

void ekf_gain_computation(struct s_kalman_update *p, double xk_t_ts, double data_t_ts, gsl_matrix *Ct);
double ekf_update(struct s_kalman_update *p, gsl_matrix *Ct);

/* kalman_template.c */
int step_ode_ekf(double t, const double X[], double f[], void *params);

int cac_drift_in_cac_ts(int cac_drift, int o, int ts_unique, struct s_obs2ts **obs2ts);
void eval_jac(gsl_matrix *jac, const double *X, struct s_par *p_par, struct s_data *p_data, struct s_calc *p_calc, struct s_group ***compo_groups_drift_par_proc, double t);
void eval_ht(struct s_kalman_update * p, double x, struct s_par *p_par, struct s_data *p_data, struct s_calc *p_calc, int ts);


void eval_Q_full(gsl_matrix *Q, const double *X, struct s_par *p_par, struct s_data *p_data, struct s_calc *p_calc, struct s_kalman_specific_data *p_kalman_specific_data, double t);
void eval_Q_no_dem_sto(gsl_matrix *Q, const double *X, struct s_par *p_par, struct s_data *p_data, struct s_calc *p_calc, struct s_kalman_specific_data *p_kalman_specific_data, double t);
void eval_Q_no_env_sto(gsl_matrix *Q, const double *X, struct s_par *p_par, struct s_data *p_data, struct s_calc *p_calc, struct s_kalman_specific_data *p_kalman_specific_data, double t);
void eval_Q_no_dem_sto_no_env_sto(gsl_matrix *Q, const double *X, struct s_par *p_par, struct s_data *p_data, struct s_calc *p_calc, struct s_kalman_specific_data *p_kalman_specific_data, double t);


/* kmcmc.c */
void kmcmc(struct s_kalman *p_kalman, struct s_likelihood *p_like, struct s_mcmc_calc_data *p_mcmc_calc_data, plom_f_pred_t f_pred,  int OPTION_ACC);
