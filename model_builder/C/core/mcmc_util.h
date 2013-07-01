#ifndef MCMC_UTIL_H
#define MCMC_UTIL_H

int OPTION_FULL_UPDATE; //TODO remove


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


struct s_mcmc_calc_data *build_mcmc_calc_data(struct s_best *p_best, const double a, const int m_switch, const int m_epsilon, const double epsilon_max, const int is_smoothed_tunning, const double alpha);
void clean_mcmc_calc_data(struct s_mcmc_calc_data *p_mcmc_calc_data);

gsl_matrix * get_var_and_sd_fac(double *sd_fac, struct s_best *p_best, struct s_mcmc_calc_data *p, struct s_calc *p_calc, int m);

void increment_iteration_counters(struct s_mcmc_calc_data *p_mcmc_calc_data, struct s_best *p_best, const int OPTION_FULL_UPDATE);

void ran_proposal_sequential(gsl_vector *proposed, struct s_best *p_best, gsl_matrix *var, double sd_fac, struct s_calc *p_calc);

void compute_best_traj(struct s_hat **D_p_hat_best, struct s_hat **D_p_hat_prev, struct s_hat **D_p_hat_new, struct s_data *p_data, double alpha, double m);

void header_acceptance_rates(FILE *p_file, struct s_data *p_data);
void print_acceptance_rates(FILE *p_file, struct s_mcmc_calc_data *p, int m_full_iteration);

void compute_acceptance_rates(struct s_best *p_best, struct s_mcmc_calc_data *p, double is_accepted, int m);

void header_covariance(FILE *p_file, struct s_data *p_data);
void print_covariance(FILE *p_file_cov, gsl_matrix *covariance);


#endif


