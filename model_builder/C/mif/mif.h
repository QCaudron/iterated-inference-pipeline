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

#define FREEZE (pow(MIF_a, m-1))

/*-------global variables--------*/
double MIF_a; /*cooling parameter of ionides 2006*/
double MIF_b; /*Ionides uses sqrt(20.0)*/
int L;  /*fixed lag for fixed lag smoothing (needed to infer the initial conditions)*/
int M; /*Number of MIF iterations*/
int SWITCH; /*number of the MIF iteration when we start to apply Ionides (2006) MIF update formulae instead of the simple mean accross particles */

int OPTION_IC_ONLY; /**< only fixed lag smoothing */

struct s_mif
{
    struct s_data *p_data;
    struct s_calc **calc;
    struct s_best *p_best;

    struct s_X **J_p_X; /* [J] */
    struct s_X **J_p_X_tmp; /* [J] */
    struct s_par **J_p_par; /* [J][N_G] N_G is for par_sv, par_proc, par_obs */

    struct s_likelihood *p_like;

    /*MIF specific*/
    gsl_vector **J_theta; /* [J][s_best.length] */
    gsl_vector **J_theta_tmp; /* [J][s_best.length] */

    double **D_theta_bart; /* [N_DATA+1][s_best.length] mean of theta at each time step, (N_DATA+1) because we keep values for every data point + initial condition */
    double **D_theta_Vt; /* [N_DATA+1][s_best.length] variance of theta at each time step */
};


/*function prototypes*/

/* methods.c */
void rescale_covariance_mif(struct s_best *p_best, struct s_data *p_data);
void ran_proposal_chol(theta_t *proposed, struct s_best *p_best, gsl_matrix *var, double sd_fac, struct s_calc *p_calc);
void fill_theta_bart_and_Vt_mif(double **D_theta_bart, double **D_theta_Vt, struct s_best *p_best, struct s_data *p_data, int m);
void mean_var_theta_theoretical_mif(double *theta_bart_n, double *theta_Vt_n, gsl_vector **J_theta, struct s_likelihood *p_like, struct s_data *p_data, struct s_best *p_best, int m, double var_fac, const enum plom_print print_opt);
void print_mean_var_theta_theoretical_mif(FILE *p_file, double *theta_bart_n, double *theta_Vt_n, struct s_likelihood *p_like, struct s_data *p_data, int m, int time);
void header_mean_var_theoretical_mif(FILE *p_file, struct s_data *p_data);
void resample_and_mut_theta_mif(unsigned int *select, gsl_vector **J_theta, gsl_vector **J_theta_tmp, struct s_calc **calc, struct s_data *p_data, struct s_best *p_best, double sd_fac, gsl_matrix *var_fitted, int is_mvn);
void update_fixed_lag_smoothing(struct s_best *p_best, struct s_likelihood *p_like, gsl_vector **J_theta, struct s_data *p_data);
void update_theta_best_stable_mif(struct s_best *p_best, double **D_theta_bart, struct s_data *p_data);
void update_theta_best_king_mif(struct s_best *p_best, double **D_theta_bart, double **D_theta_Vt, struct s_data *p_data, int m);
void patch_likelihood_prior(struct s_likelihood *p_like, struct s_best *p_best, gsl_vector **J_theta, struct s_data *p_data, int n, const int lag);

/* build.c */
struct s_mif *build_mif(json_t *theta, enum plom_implementations implementation,  enum plom_noises_off noises_off, double dt, double eps_abs, double eps_rel, const double freeze_forcing, double prop_L_option, int J, int *n_threads);
void clean_mif(struct s_mif *p_mif);

/* mif.c */
void mif(struct s_calc **calc, struct s_data *p_data, struct s_best *p_best, struct s_X ***J_p_X, struct s_X ***J_p_X_tmp, struct s_par **J_p_par, struct s_likelihood *p_like, gsl_vector **J_theta, gsl_vector **J_theta_tmp, double **D_theta_bart, double **D_theta_Vt, plom_f_pred_t f_pred, int is_mvn, const enum plom_print print_opt);
