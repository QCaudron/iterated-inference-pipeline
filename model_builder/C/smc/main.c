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

int main(int argc, char *argv[])
{
    int ch;
    char str[STR_BUFFSIZE];

    /* set default values for the options */

    char sfr_help_string[] =
        "PloM Sequential Monte Carlo\n"
        "usage:\n"
        "smc [implementation] [--no_dem_sto] [--no_env_sto] [--no_drift]\n"
        "                     [--traj] [-p, --path <path>] [-i, --id <integer>] [-P, --N_THREAD <integer>]\n"
        "                     [-t, --no_filter] [-b, --no_best] [-h, --no_hat]\n"
        "                     [-s, --DT <float>] [--eps_abs <float>] [--eps_rel <float>]\n"
        "                     [-l, --LIKE_MIN <float>] [-J <integer>] [--prior]\n"
        "                     [--help]\n"
        "where implementation is 'ode', 'sde' or 'psr' (default)\n"
        "options:\n"
	"\n"
        "--no_dem_sto       turn off demographic stochasticity (if possible)\n"
        "--no_env_sto       turn off environmental stochasticity (if any)\n"
        "--no_drift         turn off drift (if any)\n"
	"\n"
        "-s, --DT           integration time step\n"
	"--eps_abs          Absolute error for adaptive step-size contro\n"
	"--eps_rel          Relative error for adaptive step-size contro\n"
	"\n"
        "-i, --id           general id (unique integer identifier that will be appended to the output files)\n"
        "-p, --path         path where the outputs will be stored\n"
        "-P, --N_THREAD     number of threads to be used (default to the number of cores)\n"
	"\n"
        "-l, --LIKE_MIN     particles with likelihood smaller that LIKE_MIN are considered lost\n"
        "-J                 number of particles\n"
	"\n"
        "--traj             print the trajectories\n"
        "-t, --no_filter    do not filter\n"
        "-b, --no_best      do not write best_<general_id>.output file\n"
        "-h, --no_hat       do not write hat_<general_id>.output file\n"
        "-r, --no_pred_res  do not write pred_res_<general_id>.output file (prediction residuals)\n"
	"\n"
        "--prior            add log(prior) to the estimated log likelihood\n"
	"\n"
        "--help             print the usage on stdout\n";

    int filter = 1;
    int output_best = 1, output_hat = 1, output_pred_res =1;
    double dt = 0.0, eps_abs = PLOM_EPS_ABS, eps_rel = PLOM_EPS_REL;

    enum plom_implementations implementation;
    enum plom_noises_off noises_off = 0;

    
    OPTION_TRAJ = 0;
    OPTION_PRIOR = 0;
    GENERAL_ID =0;
    snprintf(SFR_PATH, STR_BUFFSIZE, "%s", DEFAULT_PATH);
    J=1;
    LIKE_MIN = 1e-17;
    LOG_LIKE_MIN = log(1e-17);
    int n_threads=omp_get_max_threads();
    

    while (1) {
        static struct option long_options[] =
            {
                /* These options set a flag. */
                {"traj", no_argument,       &OPTION_TRAJ, 1},
                /* These options don't set a flag We distinguish them by their indices (that are also the short option names). */
                {"no_dem_sto", no_argument,       0, 'x'},
                {"no_env_sto", no_argument,       0, 'y'},
                {"no_drift",   no_argument,       0, 'z'},

                {"DT",         required_argument, 0, 's'},
                {"eps_abs",    required_argument, 0, 'v'},
                {"eps_rel",    required_argument, 0, 'w'},

                {"help",       no_argument,       0, 'e'},
                {"path",       required_argument, 0, 'p'},
                {"prior",  no_argument, &OPTION_PRIOR, 1},
                {"id",         required_argument, 0, 'i'},
                {"N_THREAD",   required_argument, 0, 'P'},
                {"no_filter",  no_argument,       0, 't'},
                {"no_best",    no_argument,       0, 'b'},
                {"no_hat",     no_argument,       0, 'h'},
                {"no_pred_res",no_argument,       0, 'r'},

                {"LIKE_MIN",   required_argument, 0, 'l'},

                {0, 0, 0, 0}
            };
        /* getopt_long stores the option index here. */
        int option_index = 0;

        ch = getopt_long (argc, argv, "xyzs:v:w:p:i:J:l:tbhrP:", long_options, &option_index);

        /* Detect the end of the options. */
        if (ch == -1)
            break;

        switch (ch) {
        case 0:
            /* If this option set a flag, do nothing else now. */
            if (long_options[option_index].flag != 0) {
                break;
            }
            break;

        case 'x':
            noises_off = noises_off | PLOM_NO_DEM_STO;
            break;
        case 'y':
            noises_off = noises_off | PLOM_NO_ENV_STO;
            break;
        case 'z':
            noises_off = noises_off | PLOM_NO_DRIFT;
            break;

        case 's':
            dt = atof(optarg);
            break;
        case 'v':
            eps_abs = atof(optarg);
            break;
        case 'w':
            eps_rel = atof(optarg);
            break;


        case 'e':
            print_log(sfr_help_string);
            return 1;

        case 'p':
            snprintf(SFR_PATH, STR_BUFFSIZE, "%s", optarg);
            break;
        case 'i':
            GENERAL_ID = atoi(optarg);
            break;
        case 't':
            filter = 0;
            break;
        case 'J':
            J = atoi(optarg);
            break;
        case 'P':
            n_threads = atoi(optarg);
            break;
        case 'l':
            LIKE_MIN = atof(optarg);
            LOG_LIKE_MIN = log(LIKE_MIN);
            break;
        case 'b':
            output_best = 0;
            break;
        case 'h':
            output_hat = 0;
            break;
        case 'r':
            output_pred_res = 0;
            break;

        case '?':
            /* getopt_long already printed an error message. */
            return 1;

        default:
            snprintf(str, STR_BUFFSIZE, "Unknown option '-%c'\n", optopt);
            print_err(str);
            return 1;
        }
    }
    argc -= optind;
    argv += optind;

    if(argc == 0) {
	implementation = PLOM_PSR;
    } else {
        if (!strcmp(argv[0], "ode")) {
            implementation = PLOM_ODE;
            noises_off = noises_off | PLOM_NO_DEM_STO| PLOM_NO_ENV_STO | PLOM_NO_DRIFT;
        } else if (!strcmp(argv[0], "sde")) {
            implementation = PLOM_SDE;
        } else if (!strcmp(argv[0], "psr")) {
            implementation = PLOM_PSR;
        } else {
            print_log(sfr_help_string);
            return 1;
        }
    }

    json_t *settings = load_settings(PATH_SETTINGS);

    json_t *theta = load_json();
    struct s_data *p_data = build_data(settings, theta, implementation, noises_off, OPTION_PRIOR);
    json_decref(settings);

    int size_proj = N_PAR_SV*N_CAC + p_data->p_it_only_drift->nbtot + N_TS_INC_UNIQUE;

    struct s_par *p_par = build_par(p_data);
    struct s_hat **D_p_hat = build_D_p_hat(p_data);
    struct s_X ***D_J_p_X = build_D_J_p_X(size_proj, N_TS, p_data, dt);
    struct s_X ***D_J_p_X_tmp = build_D_J_p_X(size_proj, N_TS, p_data, dt);
    struct s_best *p_best = build_best(p_data, theta, 0);
    json_decref(theta);
    struct s_likelihood *p_like = build_likelihood();

    struct s_calc **calc = build_calc(&n_threads, GENERAL_ID, eps_abs, eps_rel, J, size_proj, step_ode, p_data);

    FILE *p_file_X = (OPTION_TRAJ==1) ? sfr_fopen(SFR_PATH, GENERAL_ID, "X", "w", header_X, p_data): NULL;
    FILE *p_file_pred_res = (output_pred_res==1) ? sfr_fopen(SFR_PATH, GENERAL_ID, "pred_res", "w", header_prediction_residuals, p_data): NULL;


#if FLAG_VERBOSE
    snprintf(str, STR_BUFFSIZE, "Starting Plom-smc with the following options: i = %d, J = %d, LIKE_MIN = %g, N_THREADS = %d", GENERAL_ID, J, LIKE_MIN, n_threads);
    print_log(str);

    int64_t time_begin, time_end;
    time_begin = s_clock();
#endif

    transform_theta(p_best, NULL, NULL, p_data, 1, 1);

    back_transform_theta2par(p_par, p_best->mean, p_data->p_it_all, p_data);
    linearize_and_repeat(D_J_p_X[0][0], p_par, p_data, p_data->p_it_par_sv);
    prop2Xpop_size(D_J_p_X[0][0], p_data);
    theta_driftIC2Xdrift(D_J_p_X[0][0], p_best->mean, p_data);

    replicate_J_p_X_0(D_J_p_X[0], p_data);

    run_SMC(D_J_p_X, D_J_p_X_tmp, p_par, D_p_hat, p_like, p_data, calc, get_f_pred(implementation, noises_off), filter, p_file_X, p_file_pred_res);


    if (OPTION_PRIOR) {
        p_like->Llike_best += log_prob_prior(p_best, p_best->mean, p_best->var, p_data);
    }


#if FLAG_VERBOSE
    time_end = s_clock();
    struct s_duration t_exec = time_exec(time_begin, time_end);
    if(filter){
        snprintf(str, STR_BUFFSIZE, "logV: %g\tcomputed with %d particles in:= %dd %dh %dm %gs n_all_fail: %d", p_like->Llike_best, (int) J, t_exec.d, t_exec.h, t_exec.m, t_exec.s, p_like->n_all_fail);
    }  else{
        snprintf(str, STR_BUFFSIZE, "Done with %d particles in:= %dd %dh %dm %gs", (int) J, t_exec.d, t_exec.h, t_exec.m, t_exec.s);
    }

    print_log(str);
#endif

    if (p_file_X) {
        sfr_fclose(p_file_X);
    }

    if (p_file_pred_res) {
        sfr_fclose(p_file_pred_res);
    }

    if (output_hat) {
        FILE *p_file_hat = sfr_fopen(SFR_PATH, GENERAL_ID, "hat", "w", header_hat, p_data);
        print_hat(p_file_hat, D_p_hat, p_data);
        sfr_fclose(p_file_hat);
    }

    if (output_best) {
        FILE *p_file_best = sfr_fopen(SFR_PATH, GENERAL_ID, "best", "w", header_best, p_data);
        print_best(p_file_best, 0, p_best, p_data, p_like->Llike_best);
        sfr_fclose(p_file_best);
    }

#if FLAG_VERBOSE
    print_log("clean up...");
#endif

    clean_calc(calc, p_data);
    clean_D_J_p_X(D_J_p_X);
    clean_D_J_p_X(D_J_p_X_tmp);
    clean_D_p_hat(D_p_hat, p_data);
    clean_best(p_best);
    clean_par(p_par);
    clean_likelihood(p_like);
    clean_data(p_data);

    return 0;
}
