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

#ifndef PLOM_H
#define PLOM_H

#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#include <sys/time.h>

#include <getopt.h>
#include <unistd.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_multimin.h>

#include <gsl/gsl_spline.h>

#include <gsl/gsl_sort.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>

#include <jansson.h> //json

//parallel computing ability
#include <zmq.h>
#include <pthread.h>

#define FREE(ppp) do {   \
        free( ppp );	\
        ppp = NULL;	\
    }while(0)

enum plom_implementations {PLOM_ODE, PLOM_SDE, PLOM_PSR};
enum plom_noises_off {PLOM_NO_DEM_STO = 1 << 0, PLOM_NO_ENV_STO = 1 << 1, PLOM_NO_DRIFT = 1 << 2 }; //several noises can be turned off

enum plom_print {PLOM_PRINT_BEST = 1 << 0, PLOM_PRINT_X = 1 << 1, PLOM_PRINT_HAT = 1 << 2, PLOM_PRINT_PRED_RES = 1 << 3, PLOM_PRINT_X_SMOOTH = 1 << 4, PLOM_PRINT_ACC = 1 << 5, PLOM_PIPE = 1 << 6, PLOM_QUIET = 1 << 7 };

typedef enum {PLOM_SUCCESS=0, PLOM_ERR_LIKE=-1} plom_err_code;


#define BUFFER_SIZE (5000 * 1024)  /**< 5000 KB buffer size for settings.json inputs */
#define STR_BUFFSIZE 255 /**< buffer for log and error strings */
#define DEFAULT_PATH "./" /**< default path for non JSON output (has to be slash appended) */
#define PATH_SETTINGS "./.settings.json"

#define FLAG_PIPELINE 0 /**< zmq pipeline */

#define FLAG_DEBUG 0
#define FLAG_VERBOSE 1
#define FLAG_WARNING 0
#define FLAG_JSON 0 /**< webApp */
#define FLAG_OMP 0

#if FLAG_OMP
#include <omp.h>
#endif

#define PLOM_EPS_ABS 1e-6 /**< absolute error control for ODEs*/
#define PLOM_EPS_REL 1e-3 /**< relative error control for ODEs*/

#define ZERO_LOG 1e-17 /**< smallest value that can be log transformed without being replaced by @c ZERO_LOG */
#define ONE_LOGIT 0.999999999 /**< largest value that can be logit transformed without being replaced by @c ONE_LOGIT */


/*-------global variables--------*/

double ONE_YEAR; /**< one year in data unit */

/* algo parameters (read from the command line) */
int J; /**< number of particles */
double LIKE_MIN; /**< minimum value of likelihood (smaller value are considered 0.0)*/
double LOG_LIKE_MIN; /**< log of @c LIKE_MIN  @see LIKE_MIN */

int POP_SIZE_EQ_SUM_SV; /**< flag to specify that the state variable sum to the population size */

/* dimensions and is_ parameters (read from JSON) */
int N_C;               /**< number of cities*/
int N_AC;              /**< number of age classes*/
int N_CAC;             /**< @c N_C*N_AC */
int N_PAR_PROC;        /**< number of parameters for process model (for a given city and a given age class) */
int N_PAR_OBS;         /**< number of parameters for the observation model (for a given time series) */
int N_PAR_SV;          /**< number of state variables of the preocess model in a given city and a given age class. NOTE: N_PAR_SV do **not** includes the observed variable (incidence or prevalence) */
int N_PAR_FIXED;       /**< number of fixed paramaters of the process model */
int N_TS;              /**< total number of time series */
int N_TS_INC;          /**< size of subset of @c N_TS containing only incidence ts */
int N_TS_INC_UNIQUE;   /**< size of subset of @c N_TS containing only non-repeated (i.e no repetition per data stream) incidece data  */
int N_DATA;            /**< length of the data set (including @c NaN) */
int N_OBS_ALL;         /**< number of type of observed variables (e.g incidence_strain1, incidence_strain_2, prevalence S + prevalence I, ...) */
int N_OBS_INC;         /**< number of incidences */
int N_OBS_PREV;        /**< number of prevalences */

int N_DRIFT;           /**< number of parameter following a diffusion */

int IS_SCHOOL_TERMS;   /**< do we need school terms data */

char SFR_PATH[STR_BUFFSIZE]; /**< path where the output files will be written */
int GENERAL_ID;              /**< general identifiant to make the output files unique */

/* option and commands */
int OPTION_PRIOR;   /**< print add the logprior to the loglik outputs */


/*-------plom core generic structures. These structures are used as nucleus for all plom population based C programs--------*/

/**
 * vector of parameters (including initial condition)
 */
typedef gsl_vector theta_t;


/**
 * aggregation of parameters (e.g initial conditions, process model parameters...)
 */

struct s_iterator
{
    int length;           /**< number of parameters*/
    unsigned int *ind;    /**< [self.length] indexes of the parameters contained in the iterator */
    unsigned int *offset; /**< [self.length] an element i of @c ind starts at offset[i] in the component @c mean or @c proposed of s_best */
    int nbtot;            /**< total number of values for all the parameters present in the iterator (including all the groups) */
};


/**
 *  map an observed state variables (of the observation model) to its
 *  corresponding time series
 */
struct s_obs2ts /* [N_OBS_ALL] */
{
    int n_ts_unique;         /**< number of unique time series: unique in the sense of the process: get rid of @c STREAM repetition */
    unsigned int *n_stream;  /**< [self.n_ts_unique] number of data stream */
    unsigned int *n_cac;     /**< [self.n_ts_unique] number of @c cac */
    unsigned int ***cac;     /**< [self.n_ts_unique][self.n_cac][2] list of @c c and @c ac (hence the 2) aggregated into the @c n_ts_unique time serie */
    int offset;  /**< index of first ts (relative to N_TS) of the observed variable */
};


/**
 * For a given time, contains the number and the index of the time
 * series containing no missing values.
 */
struct s_data_ind /*[N_DATA]*/
{
    int n_nonan; /**< number of time series without NaN at that time (n_nonan<=N_TS) */
    unsigned int *ind_nonan; /**< [self.n_nonan] index of time series without NaN*/
};



/**
 * grouping structure of a given parameter and transformation
 * functions
 */
struct s_router /* [ N_PAR_SV + N_PAR_PROC + N_PAR_OBS ] */
{
    /* grouping */
    int p;             /**< number of element that can be grouped: (@c N_CAC or @c N_TS) */
    int n_gp;          /**< how many groups */
    unsigned int *map; /**< [p] map cac or ts to group index (0 -> n_gp) */

    char *name; /**< parameter name */
    char **group_name; /**< [self.n_gp] name of the groups */

    /* transformations */
    double (**f) (double, double, double); /**< [self.n_gp] transformation (log, logit...) */
    double (**f_inv) (double, double, double); /**< [self.n_gp] inverse of f (f*f_inv=identity) */

    double (**f_scale) (double); /**< [self.n_gp] scale value (10^x or 10^-x has an effect on s_par only.) */

    double (**f_derivative) (double, double, double); /**< [self.n_gp] derivative of f */
    double (**f_inv_derivative) (double, double, double); /**< [self.n_gp] derivative of f_inv */

    double multiplier; /**< multiplier to go from the intuitive scale to the unit of data *before*  duration as been converted to rates (if relevant) */
    int is_duration; /**< boolean specifying if a duration has to be converted into a rate */

    //used for logit_ab transfo
    double *min;  /**< [self.n_gp] @c a of logit_ab (minimum) in the user unit (and discarding rate_as_duration or pow10 scaling) */
    double *max;  /**< [self.n_gp] @c b of logit_ab (maximum) same scale as min  */

    double *min_z;  /**< [self.n_gp] @c a of logit_ab in the data unit and the prediction function scale (s_par) (i.e rate_as_duration has been applied...) */
    double *max_z;  /**< [self.n_gp] @c b of logit_ab same scale as min_z */
};


/**
 * contains the untransformed version of the parameters scaled in the
 * unit of the data
 **/
struct s_par /*optional [J] */
{
    int size_natural; /**< N_PAR_SV + N_PAR_PROC + N_PAR_OBS */
    double **natural; /**< [self.size_natural][self->p_data->router[i]->n_gp] (untransformed version of theta) */
};


/**
 * map the drift state variable to the parameters where the diffusion
 * has to be applied and link to the associated volatilities
 */
struct s_drift //[N_DRIFT]
{
    /* diffusion */
    int ind_par_Xdrift_applied; /**< to which parameters the drift state variable is applied */
    int ind_volatility_Xdrift;  /**< index of the volatility of the drift state variable */
    int offset; /**< offset of the drift variable in s_X.proj */
};


/**
 * Data, covariates and navigation
 *
 * Every component of s_data are constant and read only.
 */
struct s_data{

    char **ts_name; /**< [N_TS] name of the time series */
    char **cac_name; /**< [N_CAC] name of the populations */
    char *remainder_name; /**< the remainder name (if any) */

    char u_data[2]; /**< frequency of the data (D, W, B, M, Y) */

    enum plom_implementations implementation;
    enum plom_noises_off noises_off;

    /*non fitted parameters*/
    double **data;           /**< [N_DATA][N_TS] the data */
    unsigned int *indn_data_nonan;  /**< [N_DATA_NONAN] index of data where there is at least on ts !=NaN */
    unsigned int *times;     /**< [N_DATA+1] [0] + [times in days where the data were collected] */

    struct s_data_ind **data_ind; /**< [N_DATA] an array of pointers to s_data_ind*/
    struct s_obs2ts **obs2ts;     /**< [N_OBS_ALL] an array of pointers to s_obs2ts*/

    struct s_router **routers;    /**< [ N_PAR_SV + N_PAR_PROC + N_PAR_OBS ] an array of pointers to s_router (one for each parameter) */

    struct s_drift **drift;      /**< reference to s_drift */

    struct s_iterator *p_it_all;                         /**< to iterate on every parameters */
    struct s_iterator *p_it_only_drift;                  /**< to iterate on parameters following a diffusion *only* */
    struct s_iterator *p_it_par_sv;                      /**< to iterate on the initial conditions of the state variables */
    struct s_iterator *p_it_all_no_drift;                /**< to iterate on every parameters *not* following a diffusion  */
    struct s_iterator *p_it_par_proc_par_obs_no_drift;   /**< to iterate on the parameter of the process and observation models *not* following a diffusion */
    struct s_iterator *p_it_par_sv_and_drift;            /**< to iterate on the initial conditions of the state variable *and* the parameters following a diffusion */
    struct s_iterator *p_it_noise;                       /**< to iterate on environmental stochasticity noises *only* */
    struct s_iterator *p_it_all_no_theta_remainder;      /**< to iterate on every parameters but the theta remainder */
    struct s_iterator *p_it_theta_remainder;             /**< to iterate on the theta remainder */

    int nb_obs; /**< length of the data to be used for inference */
    int nb_obs_nonan; /**< length of the data to be used for inference discarding lines where all ts are NaN */

    /* if school terms */
    unsigned int *n_terms;   /**< [N_CAC] number of terms for cac */
    double *prop_school;     /**< proportion of the year taken by school */
    double ***school_terms;  /**< [N_CAC][n_terms][2] 2 for begin and end. School_terms are in **years())* */

    double **waifw;          /**< [N_AC][N_AC] Who Acquire Infection From Whom matrix */

    /**
     *  pointer to a structure of additional data (or variable/objects)
     *  needed for specific methods.
     *
     *  Note: for thread safe data see s_calc
     */
    void *method_specific_data;
};


/**
 * Everything needed to perform computations (possibly in parallel)
 * and store transiant states in a thread-safe way
 */

struct s_calc /*[N_THREADS] : for parallel computing we need N_THREADS = omp_get_max_threads() replication of the structure...*/
{
    int n_threads; /**< the total number of threads */
    int thread_id; /**< the id of the thread where the computation are being run */
    
    gsl_rng *randgsl; /**< random number generator */

    /////////////////
    //implementations
    /////////////////

    /* Euler multinomial */
    double **prob; /*[N_PAR_SV][number of output from the compartment]*/
    unsigned int ***inc; /* [N_PAR_SV][N_CAC][number of destinations] increments vector, we keep N_CAC to be able to compute incidences matching time series that can be a summation of different cities and age classes */

    /* Gillespie */
    //  double **reaction; /*reaction matrix*/

    /* ODE*/
    const gsl_odeiv2_step_type * T;
    gsl_odeiv2_control * control;
    gsl_odeiv2_step * step;
    gsl_odeiv2_evolve * evolve;
    gsl_odeiv2_system sys;
    double *yerr;

    /* SDE */
    double *y_pred; /**< used to store y predicted for Euler Maruyama */

    //  double *gravity; /*[NUMC] additional term specific to measles model*/

    //multi-threaded sorting
    double *to_be_sorted;  /**< [J] array of the J particle to be sorted*/
    size_t *index_sorted;  /**< [J] index of the sorted weights used to compute 95% confidence interval */

    //interpolators for covariates
    gsl_interp_accel ***acc; /**< [N_PAR_FIXED][N_CAC || N_TS] an array of array of pointer to gsl_interp_accel */
    gsl_spline ***spline;    /**< [N_PAR_FIXED][N_CAC || N_TS] an array of array of pointer to gsl_spline */   
    unsigned int *n_spline;  /**< N_CAC or N_TS for each spline */

    double *pop_size_t0; /**< [N_CAC] the population size at time 0 (will be calculated with interpolated values if !pop_size_eq_sum_sv )*/

    /**
     *  Reference to s_par used to pass s_par to some GSL function
     *  that only accept an *void. Such function only received s_calc.
     *
     *  This reference should not be used outside from the integration
     *  functions f_prediction_ode_rk, f_prediction_with_drift_deter
     *  and f_prediction_with_drift_sto). Outside these functions,
     *  s_par is not guaranted to be defined.
     */
    struct s_par *p_par;

    struct s_data *p_data; /**< ref to s_data (same reason as the ref to s_par) */

    /** method specific *thread-safe* data */
    void *method_specific_thread_safe_data;

    /** this is *not* thread safe!  */
    void *method_specific_shared_data;
};

struct s_thread_smc
{
    void *context; ///< zmq context
    int thread_id;    
    int J_chunk;
    int J;
    struct s_data *p_data;
    struct s_par *p_par;
    struct s_X ***D_J_p_X;
    struct s_calc *p_calc;
    struct s_likelihood *p_like;
};

struct s_thread_mif
{
    void *context; ///< zmq context
    int thread_id;    
    int J_chunk;
    int J;
    struct s_data *p_data;
    struct s_par ** J_p_par;
    struct s_X ***J_p_X;
    struct s_calc *p_calc;
    struct s_likelihood *p_like;
};

struct s_thread_predict
{
    void *context; ///< zmq context
    int thread_id;    
    int J_chunk;
    int J;
    struct s_data *p_data;
    struct s_par ** J_p_par;
    struct s_X **J_p_X;
    struct s_calc *p_calc;
};



/**
 * Measuring duration of events
 */
struct s_duration
{
    unsigned int d; ///< days
    unsigned int h; ///< hours
    unsigned int m; ///< minutes
    double s;       ///< seconds
};


/**
 * best estimate of the parameters and everything needed to generate
 * new proposed values of the parameters (including the priors)
 */
struct s_best {

    int length;                    /**< total number of parameters including the initial condtions (sum of n_gp for all the routers) */
    theta_t *mean;                 /**< [self.length] current best estimate of the parameters (including initial conditions) */
    theta_t *proposed;             /**< [self.length] used in MCMC algorithm to store the proposed value */

    gsl_matrix *var;               /**< [self.length][self.length] The variance-covariance matrix */

    int n_to_be_estimated;         /**< nb of parameters that have to be estimated (parameters with jump_size > 0.0) */
    unsigned int *to_be_estimated; /**< [self.length] index of self.mean component that have to be estimated. Note: [self.length] and not [self.n_to_be_estimated] because in the webApp, user can uleash jump_sizes set to 0.0 > 0.0  */

    unsigned int *is_estimated; /**< [self.length] boolean: 1 = is estimated, 0 = is not  */

    /* "follow" property of theta.json
       E.g.
       value: {
         I0: {guess: 0.1},
         I1: {follow: "I0"},
       }

       A follower will inherit from all the property of the parameter
        it follows (prior, transformation, **grouping**, guess, min,
        max) but its sd_transf will be set to 0.0.  However a
        parameter can only follow a parameter of the same category
        (categories are process model parameters (including initial
        conditions) and observation process parameters) the reason is
        that the groups are different between process model parameters
        and observation process parameters.
    */
    int n_follow;
    unsigned int *follower;        /**< [self.n_follow] index of follower */
    unsigned int *follow;          /**< [self.n_follow] index of parameter being followed by the follower */

    unsigned int *is_follower;     /**< [length] boolean: index of mean, 1 = follow, 0 = doesn't follow */

    /* used to store states necessary for computation of the sampling covariance in MCMC algo */
    double *mean_sampling;         /**< [self.length] Em(X) 1st order mean needed to compute the sampling covariance */
    gsl_matrix *var_sampling;      /**< [self.length][self.length] Sampling covariance */

    /*priors*/
    double **par_prior; /**< parameters of the priors [self.length][2] 2: 2 parameters (a, b). Meaning of a and b depend on the prior:
                           - double gsl_ran_flat_pdf (double x, double a, double b)
                           - double gsl_ran_gamma_pdf (double x, double a, double b)
                           - double gsl_ran_gaussian_pdf (double x, double sigma) <-wraped to be  (double x, double mu, double sigma)
                           - double gsl_ran_beta_pdf (double x, double a, double b)
                           - double gsl_ran_lognormal_pdf (double x, double zeta, double sigma) */

    double (**prior) (double x, double a, double b); /**< [self.length] an array of function pointer pointing to the priors of each parameter */
};


/**
 * The state variables
 */
struct s_X /* optionaly [N_DATA+1][J] for MIF and pMCMC "+1" is for initial condition (one time step before first data)  */
{
    double *proj;    /**< [self.size_proj] x integrated (projected) (ODE, MARKOV...) */
    double *obs;     /**< [self.size_obs] x observed  matching the data (N_TS time series) */

    double dt;        /**< the integration time step (for ODE solved with adaptive time step solvers) */
    double dt0;       /**< the integration time step initially picked by the user */
};


/**
 * The best estimates of the states and the observed variables
 * (weighted average of projected values, weighted by the likelihood)
 */

struct s_hat /* ([N_DATA+1]) */
{
    double *state;         /**< [N_PAR_SV*N_CAC] best estimates of the states variables */
    double **state_95;     /**< [N_PAR_SV*N_CAC][2] 2.5% and 97.5% quantile of the estimates of the states variables*/

    double *remainder;         /**< [N_CAC] best estimates of the remainder */
    double **remainder_95;     /**< [N_CAC][2] 2.5% and 97.5% quantile of the estimates of the remainder*/

    double *obs;           /**< [N_TS] best estimates of the observed states variables */
    double **obs_95;       /**< [N_TS][2] 2.5% and 97.5% quantile of the estimates of the observed states variables*/

    double *drift;         /**< [p_data->p_it_only_drift->nbtot] best estimates of the diffusion */
    double **drift_95;     /**< [p_data->p_it_only_drift->nbtot][2] 2.5% and 97.5% quantile of the estimates of the diffusion */
};

struct s_kalman_update
{
    gsl_vector *xk;  /**< [N_KAL] concatenation of non-overlapping components of X->proj, X->obs and X->drift */
    gsl_vector *kt;  /**< [N_KAL] Kalman Gain vector */
    gsl_vector *ht;  /**< [N_KAL] Gradient of the observation function */

    double sc_st;         /**< Innovation or residual covariance */
    double sc_pred_error; /**< Innovation or measurement residual */
    double sc_rt;         /**< observation process variance */

    //temporary variables
    gsl_vector *v_n_kal;  /**< temporary vector of size N_KAL */
    gsl_matrix *M_symm_n_kal;  /**< temporary symmetric matrix of size N_KAL */
    gsl_matrix *M_symm_n_kal2;  /**< another temporary symmetric matrix of size N_KAL */

    gsl_eigen_symmv_workspace *w_eigen_vv_nkal;  /**< workspace to compute eigen values and eigen vector for symmetric matrix of size N_KAL */
    gsl_vector *eval_nkal; /**< eigen values of symmetric matrix of size N_KAL  */
    gsl_matrix *evec_nkal; /**<eigen vector of of symmetric matrix of size N_KAL */
};




/**
 * likelihood values and associated quantities
 */
struct s_likelihood
{
    double ess_n;               /**< effective sample size at n (sum(weight))^2 / sum(weight^2)*/
    double Llike_best_n ;       /**< log likelihood for the best parameter at n*/
    double Llike_best;          /**< log likelihood for the best parameter*/
    double *weights;            /**< [J] the weights */

    unsigned int **select;      /**< N_DATA_NONAN][J] select is a vector with the indexes of the resampled particles. Note that we keep @c N_DATA_NONAN values to keep genealogies */

    int n_all_fail;             /**< number of times when every particles had like < LIKE_MIN within one iteration */

    /* for bayesian methods */
    double Llike_prev;
    double Llike_new;
};


/**
 * composition of a group in term of @c cac or @c ts
 */
struct s_group
{
    int length; /**< nb of element of the group */
    unsigned int *elements; /**< element id */
};



/**
 * prediction function
 */
typedef void (*plom_f_pred_t) (struct s_X *, double, double, struct s_par *, struct s_data *, struct s_calc *);


/*-------plom core functions--------*/

/*init_d.c*/
double *init1d_set0(int n);
double **init2d_set0(int n, int p);
double **init2d_var_set0(int n, unsigned int *p);
double ***init3d_set0(int n, int p1, int p2);
double ***init3d_var_set0(int n, unsigned int *p1, unsigned int **p2);
double ***init3d_varp1_set0(int n, unsigned int *p1, int p2);
double ***init3d_varp2_set0(int n, unsigned int p1, unsigned int *p2);
double ****init4d_set0(int n, int p1, int p2, int p3);

void clean2d(double **tab, int n);
void clean3d(double ***tab, int n, int p1);
void clean3d_var(double ***tab, int n, unsigned int *p1);
void clean4d(double ****tab, int n, int p1, int p2);

/*init_u.c*/
unsigned int *init1u_set0(int n);
unsigned int **init2u_set0(int n, int p);
unsigned int **init2u_var_set0(int n, unsigned int *p);
unsigned int ***init3u_set0(int n, int p1, int p2);
unsigned int ***init3u_var_set0(int n, unsigned int *p1, unsigned int **p2);
unsigned int ***init3u_varp1_set0(int n, unsigned int *p1, int p2);
unsigned int ***init3u_varp2_set0(int n, unsigned int p1, unsigned int *p2);
unsigned int ****init4u_set0(int n, int p1, int p2, int p3);

void clean2u(unsigned int **tab, int n);
void clean3u(unsigned int ***tab, int n, int p1);
void clean3u_var(unsigned int ***tab, int n, unsigned int *p1);
void clean4u(unsigned int ****tab, int n, int p1, int p2);

/* init_gsl_vec.c */
gsl_vector **init2_gsl_vector_d_set0(int n, int p);
void clean2_gsl_vector_d(gsl_vector **tab, int n);

/* init_st.c */
size_t *init1st_set0(int n);
size_t **init2st_set0(int n, int p);
void clean2st(size_t **tab, int n);

/* init_c.c */
char *init1c(int n);
char **init2c(int n, int p);
void clean2c(char **tab, int n);

/*load.c*/
unsigned int n2d(char *filename, int Np);
unsigned int n2d_nan(char *filename, int Np);
void load1u(unsigned int *tab, char *filename);
void load1d(double *tab, char *filename);
void load2u(unsigned int **tab, char *filename, int Np);
void load2d(double **tab, char *filename, int Np);
void load2d_nan(double **tab, char *filename, int Np);

void load_best(struct s_best *p_best, struct s_data *p_data, json_t *theta, int update_guess);
void load_covariance(gsl_matrix *covariance, json_t *covjson, struct s_data *p_data);
json_t *load_settings(const char *path);

/* build.c */

struct s_iterator *plom_iterator_new();
void plom_add_offset_iterator(struct s_iterator *p_it, unsigned int *all_offset, struct s_router **routers);

struct s_iterator *plom_iterator_all_new(unsigned int *all_offset, struct s_router **routers);
struct s_iterator *plom_iterator_all_no_theta_remainder_new(unsigned int *all_offset, struct s_router **routers, int ind_theta_remainder);
struct s_iterator *plom_iterator_theta_remainder_new(unsigned int *all_offset, struct s_router **routers, int ind_theta_remainder);
struct s_iterator *plom_iterator_only_drift_new(unsigned int *all_offset, struct s_router **routers, struct s_drift **drift, const enum plom_noises_off noises_off);
struct s_iterator *plom_iterator_par_sv_new(unsigned int *all_offset, struct s_router **routers);
struct s_iterator *plom_iterator_all_no_drift_new(unsigned int *all_offset, struct s_router **routers, struct s_drift **drift, const enum plom_noises_off noises_off);
struct s_iterator *plom_iterator_par_proc_par_obs_no_drift_new(unsigned int *all_offset, struct s_router **routers, struct s_drift **drift, const enum plom_noises_off noises_off);
struct s_iterator *plom_iterator_par_sv_and_drift_new(unsigned int *all_offset, struct s_router **routers, struct s_drift **drift, const enum plom_noises_off noises_off);
struct s_iterator *plom_iterator_noise_new(unsigned int *all_offset, struct s_router **routers, json_t *settings);

void clean_iterator(struct s_iterator *p_it);
struct s_router *build_router(const json_t *par, const char *par_key, const json_t *partition, const json_t *order, const char *link_key, const char *u_data, int is_bayesian);
void clean_router(struct s_router *p_router);
struct s_router **build_routers(int *ind_theta_remainder, json_t *settings, json_t *theta, const char *u_data, int is_bayesian);
void clean_routers(struct s_router **routers);
int index_of_json_array(const json_t *array, const char *element);
struct s_par *build_par(struct s_data *p_data);
void clean_par(struct s_par *p_par);
struct s_par **build_J_p_par(struct s_data *p_data);
void clean_J_p_par(struct s_par **J_p_par);
struct s_obs2ts **build_obs2ts(json_t *json_obs2ts);
void clean_obs2ts(struct s_obs2ts **obs2ts);
struct s_drift **build_drift(json_t *json_drift, struct s_router **routers);
void clean_drift(struct s_drift **drift);
struct s_data *build_data(json_t *settings, json_t *theta, enum plom_implementations implementation, enum plom_noises_off noises_off, int is_bayesian, int nb_obs, const char *u_data);
void clean_data(struct s_data *p_data);

struct s_calc **build_calc(int *n_threads, int general_id, double eps_abs, double eps_rel, int J, int dim_ode, int (*func_step_ode) (double, const double *, double *, void *), const double freeze_forcing, const int t_max, struct s_data *p_data, json_t *settings);
struct s_calc *build_p_calc(int n_threads, int thread_id, int seed, double eps_abs, double eps_rel, int dim_ode, int (*func_step_ode) (double, const double *, double *, void *), const double freeze_forcing, const int t_max, struct s_data *p_data, json_t *settings);

void clean_p_calc(struct s_calc *p_calc, struct s_data *p_data);
void clean_calc(struct s_calc **calc, struct s_data *p_data);

struct s_X *build_X(int size_proj, int size_obs, struct s_data *p_data, double dt);
struct s_X **build_J_p_X(int size_proj, int size_obs, struct s_data *p_data, double dt);
struct s_X ***build_D_J_p_X(int size_proj, int size_obs, struct s_data *p_data, double dt);
void clean_X(struct s_X *p_X);
void clean_J_p_X(struct s_X **J_p_X);
void clean_D_J_p_X(struct s_X ***D_J_p_X);
struct s_hat **build_D_p_hat(struct s_data *p_data);
void clean_D_p_hat(struct s_hat **D_p_hat, struct s_data *p_data);
struct s_hat *build_hat(struct s_data *p_data);
void clean_hat(struct s_hat *p_hat, struct s_data *p_data);
struct s_likelihood *build_likelihood(void);
void clean_likelihood(struct s_likelihood *p_like);
struct s_best *build_best(struct s_data *p_data, json_t *theta);
void clean_best(struct s_best *p_best);


/*webio.c*/
void ask_update();
void block();
void update_walk_rates(struct s_best *p_best, struct s_data *p_data);
/*print.c*/
int plom_unlink_done(const char* path, const int general_id);
void plom_print_done(json_t *theta, struct s_data *p_data, struct s_best *p_best, const char* path, const int general_id, const enum plom_print print_opt);
FILE *plom_fopen(const char* path, const int general_id, const char* file_name, const char *mode, void (*header)(FILE*, struct s_data *), struct s_data *p_data);
void plom_fclose(FILE *p_file);
void print_warning(char *msg);
void print_log(char *msg);
void print_err(char *msg);

void print_p_X(FILE *p_file, json_t *json_print, struct s_X *p_X, struct s_par *p_par, struct s_data *p_data, struct s_calc *p_calc, int j_or_m, const int n, const double t);
void print_trace(FILE *p_file_trace, int m, struct s_best *p_best, struct s_data *p_data, double log_like);
void print_p_hat(FILE *p_file, json_t *json_print, struct s_hat *p_hat, struct s_data *p_data, const double t);

void print_hat(FILE *p_file, struct s_hat **D_p_hat, struct s_data *p_data);
void print_par(struct s_par *p_par, struct s_data *p_data);
void print_prediction_residuals(FILE *p_file_pred_res, struct s_par **J_p_par, struct s_data *p_data, struct s_calc *p_calc, struct s_X **J_p_X, double llike_t, double ess_t, int is_p_par_cst, const int n, const double t);


void sample_traj_and_print(FILE *p_file, struct s_X ***D_J_p_X, struct s_par *p_par, struct s_data *p_data, struct s_likelihood *p_like, struct s_calc *p_calc, int m);
void print_X(FILE *p_file_X, struct s_par **J_p_par, struct s_X **J_p_X, struct s_data *p_data, struct s_calc *p_calc, int is_p_par_cst, int is_m, int m, const int n, const double t);

void header_X(FILE *p_file, struct s_data *p_data);
void header_prediction_residuals(FILE *p_file, struct s_data *p_data);
void header_hat(FILE *p_file, struct s_data *p_data);
void header_trace(FILE *p_file, struct s_data *p_data);

/*prediction_util.c*/
void reset_inc(struct s_X *p_X, struct s_data *p_data);
void round_inc(struct s_X *p_X, struct s_data *p_data);

double sum_SV(const double *X_proj, int cac);
double correct_rate(double rate, double dt);

void linearize_and_repeat(struct s_X *p_X, struct s_par *p_par, struct s_data *p_data, const struct s_iterator *p_it);
    void prop2Xpop_size(struct s_X *p_X, struct s_data *p_data, struct s_calc *p_calc);
void theta_driftIC2Xdrift(struct s_X *p_X, const theta_t *best_mean, struct s_data *p_data);

plom_f_pred_t get_f_pred(enum plom_implementations implementation, enum plom_noises_off noises_off);

void f_prediction_ode(struct s_X *p_X, double t0, double t1, struct s_par *p_par, struct s_data *p_data, struct s_calc *p_calc);

void f_prediction_sde_no_dem_sto_no_env_sto(struct s_X *p_X, double t0, double t1, struct s_par *p_par, struct s_data *p_data, struct s_calc *p_calc);
void f_prediction_sde_no_dem_sto_no_drift(struct s_X *p_X, double t0, double t1, struct s_par *p_par, struct s_data *p_data, struct s_calc *p_calc);
void f_prediction_sde_no_env_sto_no_drift(struct s_X *p_X, double t0, double t1, struct s_par *p_par, struct s_data *p_data, struct s_calc *p_calc);
void f_prediction_sde_no_dem_sto(struct s_X *p_X, double t0, double t1, struct s_par *p_par, struct s_data *p_data, struct s_calc *p_calc);
void f_prediction_sde_no_env_sto(struct s_X *p_X, double t0, double t1, struct s_par *p_par, struct s_data *p_data, struct s_calc *p_calc);
void f_prediction_sde_no_drift(struct s_X *p_X, double t0, double t1, struct s_par *p_par, struct s_data *p_data, struct s_calc *p_calc);
void f_prediction_sde_full(struct s_X *p_X, double t0, double t1, struct s_par *p_par, struct s_data *p_data, struct s_calc *p_calc);

void f_prediction_psr(struct s_X *p_X, double t0, double t1, struct s_par *p_par, struct s_data *p_data, struct s_calc *p_calc);
void f_prediction_psr_no_drift(struct s_X *p_X, double t0, double t1, struct s_par *p_par, struct s_data *p_data, struct s_calc *p_calc);


void plom_ran_multinomial(const gsl_rng * r, const size_t K, unsigned int N, const double p[], unsigned int n[]);


void *jac;

double foi(double *X_c, double *waifw_ac);


/* special_functions.c */
double terms_forcing(double amplitude, double time, struct s_data *p_data, int cac);
double heaviside(double x);
double ramp(double x);


/* func.c */
int64_t s_clock (void);
struct s_duration time_exec(int64_t s1, int64_t s2);
unsigned int sum1u(unsigned int *tab, int length_tab);
void online_mean_var(double *x, int N_x, double *mean, double *var);
int get_thread_id(void);
int get_min_u(unsigned int *tab, int length_tab);
int get_max_u(unsigned int *tab, int length_tab);
void update_to_be_estimated(struct s_best *p_best);
int sanitize_n_threads(int n_threads, int J);
int plom_sanitize_nb_obs(int nb_obs, int n_data);
//void store_state_current_m(struct s_calc **calc, int m);
int in_u(int i, unsigned int *tab, int length);
int in_drift(int i, struct s_drift **drift);

/* transform.c */
double f_id(double x, double a, double b);
double f_log(double x, double a, double b);
double f_inv_log(double x, double a, double b);
double f_logit(double x, double a, double b);
double f_inv_logit(double x, double a, double b);
double f_logit_ab(double x, double a, double b);
double f_inv_logit_ab(double x, double a, double b);

double f_scale_pow10(double x);
double f_scale_pow10_neg(double x);
double f_scale_id(double x);

double f_der_log(double x, double a, double b);
double f_der_inv_log(double x, double a, double b);
double f_der_logit(double x, double a, double b);
double f_der_inv_logit(double x, double a, double b);
double f_der_logit_ab(double x, double a, double b);
double f_der_inv_logit_ab(double x, double a, double b);


double u_duration_par2u_data(const char *u_par, const char *u_data);
double get_multiplier(const char *u_data, const json_t *par, int is_print);
int is_duration(const json_t *par);
void set_f_trans(struct s_router *p_router, const char *transf, const char *prior_type, const char *u_data, int g, int is_bayesian);
void set_ab_z(struct s_router *r, int g);

void assign_f_transfo(double (**f_transfo) (double x, double mul, double a, double b), const char *f_transfo_name);
void assign_f_derivative(double (**f_derivative) (double x, double mul, double a, double b), const char *f_transfo_name);
void back_transform_theta2par(struct s_par *p_par, const theta_t *theta, const struct s_iterator *p_it, struct s_data *p_data);
double back_transform_x(double x, int g, struct s_router *r);
void transform_theta(struct s_best *p_best, struct s_data *p_data, int opt_square_diag_sd);
void square_diag_sd(struct s_best *p_best, struct s_data *p_data);

/* likelihood.c */
double get_smallest_log_likelihood(struct s_data *p_data);
double get_log_likelihood(struct s_X *p_X, struct s_par *p_par, struct s_data *p_data, struct s_calc *p_calc, const int n, const double t);
double get_sum_square(struct s_X *p_X, struct s_par *p_par, struct s_data *p_data, struct s_calc *p_calc, const int n, const double t);
double sanitize_likelihood(double like);


/* swap.c */
void swap_1d(double **A, double **tmp_A);
void swap_2d(double ***A, double ***tmp_A);
void swap_gsl_vector(gsl_vector **A, gsl_vector **tmp_A);
void swap_X(struct s_X ***X, struct s_X ***tmp_X);
void swap_D_p_hat(struct s_hat ***D_phat, struct s_hat ***tmp_D_p_hat);

/* smc.c */
int weight(struct s_likelihood *p_like, int n);
void systematic_sampling(struct s_likelihood *p_like, struct s_calc *p_calc, int n);
void multinomial_sampling(struct s_likelihood *p_like, struct s_calc *p_calc, int n);
void resample_X(unsigned int *select, struct s_X ***J_p_X, struct s_X ***J_p_X_tmp, struct s_data *p_data);
void replicate_J_p_X_0(struct s_X **J_p_X, struct s_data *p_data);

void run_SMC(struct s_X ***D_J_p_X, struct s_X ***D_J_p_X_tmp, struct s_par *p_par, struct s_hat **D_p_hat, struct s_likelihood *p_like, struct s_data *p_data, struct s_calc **calc, plom_f_pred_t f_pred, int option_filter, FILE *p_file_X, FILE *p_file_hat, FILE *p_file_pred_res, const enum plom_print print_opt);

void run_SMC_zmq(struct s_X ***D_J_p_X, struct s_X ***D_J_p_X_tmp, struct s_par *p_par, struct s_hat **D_p_hat, struct s_likelihood *p_like, struct s_data *p_data, struct s_calc **calc, plom_f_pred_t f_pred, int Jchunk, const enum plom_print print_opt, void *sender, void *receiver, void *controller);

void run_SMC_zmq_inproc(struct s_X ***D_J_p_X, struct s_X ***D_J_p_X_tmp, struct s_par *p_par, struct s_hat **D_p_hat, struct s_likelihood *p_like, struct s_data *p_data, struct s_calc **calc, plom_f_pred_t f_pred, int option_filter, FILE *p_file_X, FILE *p_file_hat, FILE *p_file_pred_res, const enum plom_print print_opt, void *sender, void *receiver, void *controller);

/* metropolis_hastings_prior.c */
int metropolis_hastings(struct s_best *p_best, struct s_likelihood *p_like, double *alpha, struct s_data *p_data, struct s_calc *p_calc, gsl_matrix *var, double sd_fac, int is_mvn);

plom_err_code log_prob_prior(double *log_like, struct s_best *p_best, gsl_vector *mean, gsl_matrix *var, struct s_data *p_data);
double normal_prior(double x, double min, double max);
double pseudo_unif_prior(double x, double min, double max);

/* proposal.c */
void propose_safe_theta_and_load_X0(theta_t *proposed, struct s_best *p_best, gsl_matrix *var, double sd_fac, struct s_par *p_par, struct s_X *p_X, struct s_data *p_data, struct s_calc *p_calc, void (*ran_proposal) (theta_t *proposed, struct s_best *p_best, gsl_matrix *var, double sd_fac, struct s_calc *p_calc));

void ran_proposal(theta_t *proposed, struct s_best *p_best, gsl_matrix *var, double sd_fac, struct s_calc *p_calc);

int plom_check_IC_assign_theta_remainder(theta_t *proposed, struct s_data *p_data);
plom_err_code log_prob_proposal(double *log_like, struct s_best *p_best, theta_t *proposed, theta_t *mean, gsl_matrix *var, double sd_fac, struct s_data *p_data, int is_mvn);
void apply_following_constraints(theta_t *proposed, struct s_best *p_best, struct s_data *p_data);

/* hat.c */
void get_CI95(double *hat_95, const double *to_be_sorted, size_t *index_sorted, double *weights);
void compute_hat(struct s_X **J_p_X, struct s_par *p_par, struct s_data *p_data, struct s_calc **calc, struct s_hat *p_hat, double *weights, const int n, const double t);
void compute_hat_nn(struct s_X **J_p_X, struct s_par **J_p_par, struct s_data *p_data, struct s_calc **calc, struct s_hat *p_hat, int is_p_par_cst, const int n, const double t);
void plom_plug_hat(struct s_best *p_best, struct s_hat *p_hat, struct s_data *p_data);


/* json.c */
json_t *fast_get_json_object(const json_t *container, const char *obj_name);
double fast_get_json_real_from_object(json_t *object, const char *key);
const char *fast_get_json_string_from_object(const json_t *object, const char *key);
json_t *fast_get_json_array(const json_t *container, const char *array_name);
double fast_get_json_real_from_array(json_t *array, int i, const char *array_name);
const char *fast_get_json_string_from_array(json_t *array, int i, const char *array_name);
int fast_get_json_boolean(json_t *container, char *obj_name);
int fast_get_json_integer(json_t *container, char *obj_name);
double fast_get_json_real(json_t *container, char *obj_name);
unsigned int *fast_load_fill_json_1u(json_t *array, char *array_name);
unsigned int **fast_load_fill_json_2u(json_t *array, char *array_name);
unsigned int ***fast_load_fill_json_3u(json_t *array, char *array_name);
double *fast_load_fill_json_1d(json_t *array, char *array_name);
double **fast_load_fill_json_2d(json_t *array, char *array_name);
double ***fast_load_fill_json_3d(json_t *array, char *array_name);
json_t *load_json(void);
json_t *plom_theta_remainder_new(json_t* theta);

/* drift.c */
void compute_drift(struct s_X *p_X, struct s_par *p_par, struct s_data *p_data, struct s_calc *p_calc);

/* group.c */
void get_c_ac(int cac, int *c, int *ac);
struct s_group **get_groups_compo(struct s_router *p_router);
void clean_groups_compo(struct s_group **compo, int n_gp);


/* simplex.c */
void transfer_estimated(struct s_best *p_best, const gsl_vector *x, struct s_data *p_data);
void simplex(struct s_best *p_best, struct s_data *p_data, void *p_params_simplex, double (*f_simplex)(const gsl_vector *, void *), double CONVERGENCE_STOP_SIMPLEX, int M, enum plom_print print_opt);

/* zhelpers.c */
void send_par(void *socket, const struct s_par *p_par, struct s_data *p_data, int zmq_options);
void recv_par(struct s_par *p_par, struct s_data *p_data, void *socket);
void send_X(void *socket, const struct s_X *p_X, struct s_data *p_data, int zmq_options);
void recv_X(struct s_X *p_X, struct s_data *p_data, void *socket);


/* mvn.c */
int rmvnorm(const gsl_rng *r, const int n, const gsl_vector *mean, const gsl_matrix *var, gsl_vector *result);
double dmvnorm(const int n, const gsl_vector *x, const gsl_vector *mean, const gsl_matrix *var);
void plom_rmvnorm(gsl_vector *proposed, struct s_best *p_best, gsl_matrix *var, double sd_fac, struct s_calc *p_calc);
double plom_dmvnorm(struct s_best *p_best, theta_t *proposed, gsl_vector *mean, gsl_matrix *var, double sd_fac);
void eval_var_emp(struct s_best *p_best, double m);

/* templated */
void build_psr(struct s_calc *p);

double likelihood(double x, struct s_par *p_par, struct s_data *p_data, struct s_calc *p_calc, const int ts, const int n, const double t);
double obs_mean(double x, struct s_par *p_par, struct s_data *p_data, struct s_calc *p_calc, const int ts, const int n, const double t);
double obs_var(double x, struct s_par *p_par, struct s_data *p_data, struct s_calc *p_calc, const int ts, const int n, const double t);
double observation(double x, struct s_par *p_par, struct s_data *p_data, struct s_calc *p_calc, const int ts, const int n, const double t);

void proj2obs(struct s_X *p_X, struct s_data *p_data);

void step_psr(struct s_X *p_X, double t, struct s_par *p_par, struct s_data *p_data, struct s_calc *p_calc);

int step_ode(double t, const double X[], double f[], void *params);

void step_sde_full(struct s_X *p_X, double t, struct s_par *p_par, struct s_data *p_data, struct s_calc *p_calc);
void step_sde_no_dem_sto(struct s_X *p_X, double t, struct s_par *p_par, struct s_data *p_data, struct s_calc *p_calc);
void step_sde_no_env_sto(struct s_X *p_X, double t, struct s_par *p_par, struct s_data *p_data, struct s_calc *p_calc);
void step_sde_no_dem_sto_no_env_sto(struct s_X *p_X, double t, struct s_par *p_par, struct s_data *p_data, struct s_calc *p_calc);


/* worker_inproc.c */
void *worker_routine_smc_inproc(void *params);
void *worker_routine_mif_inproc(void *params);
void *worker_routine_predict_inproc(void *params);

#endif
