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

int main(int argc, char *argv[])
{
    char ch;
    char str[STR_BUFFSIZE];

    /* set default values for the options */
    char sfr_help_string[] =
        "PLoM pMCMC\n"
        "usage:\n"
        "pmcmc [implementation] [--no_dem_sto] [--no_env_sto] [--no_drift]\n"
        "                [-s, --DT <float || 0.25 day>] [--eps_abs <float || 1e-6>] [--eps_rel <float || 1e-3>]\n"
        "                [--full] [--traj] [-k, --n_traj <int || 1000>] [--acc] [-p, --path <path>] [-i, --id <integer || 0>] [-P, --N_THREAD <integer || N_CPUs>]\n"
        "                [-l, --LIKE_MIN <float || 1e-17>] [-J <integer || 1>] [-M, --iter <integer || 10>]\n"
        "                [-C --cov] [-a --cooling <float || 0.999>] [-S --switch <int || 5*n_par_fitted^2 >] "
        "                [-E --epsilon <int || 50>] [--epsilon_max <float || 2.0>] [--smooth] [--alpha <float || 0.02>]"
        "                [-Z, --zmq] [-c, --chunk <integer>]\n"
        "                [--help]\n"
        "where implementation is 'ode', 'sde' or 'psr' (default)\n"
        "options:\n"
        "\n"
        "--no_dem_sto       turn off demographic stochasticity (if possible)\n"
        "--no_env_sto       turn off environmental stochasticity (if any)\n"
        "--no_drift         turn off drift (if any)\n"
        "\n"
        "-s, --DT           Initial integration time step\n"
        "--eps_abs          Absolute error for adaptive step-size contro\n"
        "--eps_rel          Relative error for adaptive step-size contro\n"
        "\n"
        "--full             full update MVN mode\n"
        "-a, --cooling      cooling factor for sampling covariance live tuning\n"
        "-S, --switch       select switching iteration from initial covariance to empirical one\n"
        "-E, --epsilon      select number of burnin iterations before tuning epsilon\n"
        "--epsilon_max      maximum value allowed for epsilon\n"
        "--smooth           tune epsilon with the value of the acceptance rate obtained with exponential smoothing\n"
        "--alpha            smoothing factor of exponential smoothing used to compute the smoothed acceptance rate (low values increase degree of smoothing)\n"
        "\n"
        "--traj             print the smoothed trajectories\n"
        "-n, --n_traj       number of trajectories stored\n"
        "--acc              print the acceptance rate\n"
        "-C, --cov          load an initial covariance from the settings\n"
        "-p, --path         path where the outputs will be stored\n"
        "-i, --id           general id (unique integer identifier that will be appended to the output files)\n"
        "-P, --N_THREAD     number of threads to be used (default to the number of cores)\n"
        "-s, --DT           integration time step\n"
        "-l, --LIKE_MIN     particles with likelihood smaller that LIKE_MIN are considered lost\n"
        "-J                 number of particles\n"
        "-M, --iter         number of pMCMC iterations\n"
        "-Z, --zmq          dispatch particles across machines using a zeromq pipeline\n"
        "-c, --chunk        number of particles send to each machine\n"
	"-o, --nb_obs         number of observations to be fitted (for tempering)"
        "--help             print the usage on stdout\n";

    double dt = 0.0, eps_abs = PLOM_EPS_ABS, eps_rel = PLOM_EPS_REL;
    int load_cov = 0;
    int m_switch = -1;
    int m_eps = 50;
    double a = 0.999;
    double epsilon_max = 2.0;
    double alpha = 0.02;
    int n_traj = 1000;

    static int is_smooth = 0;

    JCHUNK=1;
    GENERAL_ID =0;
    snprintf(SFR_PATH, STR_BUFFSIZE, "%s", DEFAULT_PATH);
    J=1;
    LIKE_MIN = 1e-17;
    LOG_LIKE_MIN = log(1e-17);
    M = 10;
    int n_threads = omp_get_max_threads();
    OPTION_PIPELINE = 0;
    OPTION_FULL_UPDATE = 0;
    N_DATA_FORCED = -1;

    enum plom_implementations implementation;
    enum plom_noises_off noises_off = 0;
    enum plom_print print_opt = 0;

    while (1) {
        static struct option long_options[] =
            {
                /* These options set a flag. */
                {"full", no_argument, &OPTION_FULL_UPDATE, 1},

                /* These options don't set a flag We distinguish them by their indices (that are also the short option names). */
                {"traj",       no_argument,       0, 'j'},
                {"acc",        no_argument,       0, 'r'},
                {"no_dem_sto", no_argument,       0, 'x'},
                {"no_env_sto", no_argument,       0, 'y'},
                {"no_drift",   no_argument,       0, 'z'},

                {"DT",         required_argument, 0, 's'},
                {"eps_abs",    required_argument, 0, 'v'},
                {"eps_rel",    required_argument, 0, 'w'},

                {"n_traj",     required_argument, 0, 'n'},

                {"switch",     required_argument,   0, 'S'},
                {"epsilon",     required_argument,   0, 'E'},
                {"cooling",     required_argument,   0, 'a'},
                {"smooth",    no_argument, &is_smooth, 1},
                {"epsilon_max", required_argument, 0, 'f'},
                {"alpha",    required_argument, 0, 'g'},

                {"cov", no_argument, 0, 'C'},

                {"help", no_argument,  0, 'e'},
                {"path",    required_argument, 0, 'p'},
                {"id",    required_argument, 0, 'i'},
                {"N_THREAD",  required_argument,       0, 'P'},

                {"LIKE_MIN",     required_argument,   0, 'l'},
                {"iter",     required_argument,   0, 'M'},
                {"zmq",     no_argument,   0, 'Z'},
                {"chunk",     required_argument,   0, 'c'},
		{"nb_obs", required_argument,  0, 'o'},

                {0, 0, 0, 0}
            };
        /* getopt_long stores the option index here. */
        int option_index = 0;

        ch = getopt_long (argc, argv, "rjxyzs:v:w:Ci:J:l:M:p:c:P:ZS:E:a:f:g:n:o:", long_options, &option_index);

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
	case 'o':
	    N_DATA_FORCED = atoi(optarg);
            break;
        case 'a':
            a = atof(optarg);
            break;
        case 'S':
            m_switch = atoi(optarg);
            break;
        case 'E':
            m_eps = atoi(optarg);
            break;
        case 'f':
            epsilon_max = atof(optarg);
            break;
        case 'g':
            alpha = atof(optarg);
            break;
        case 'e':
            print_log(sfr_help_string);
            return 1;
        case 'C':
            load_cov = 1;
            break;
        case 'Z':
            OPTION_PIPELINE = 1;
            break;
        case 'p':
            snprintf(SFR_PATH, STR_BUFFSIZE, "%s", optarg);
            break;
        case 'P':
            n_threads = atoi(optarg);
            break;
        case 'i':
            GENERAL_ID = atoi(optarg);
            break;
        case 'J':
            J = atoi(optarg);
            break;
        case 'c':
            JCHUNK = atoi(optarg);
            break;
        case 'l':
            LIKE_MIN = atof(optarg);
            LOG_LIKE_MIN = log(LIKE_MIN);
            break;
        case 'M':
            M = atoi(optarg);
            break;
        case 'j':
	    print_opt |= PLOM_PRINT_X_SMOOTH;
            break;
        case 'r':
	    print_opt |= PLOM_PRINT_ACC;
            break;
        case 'n':
            n_traj = atoi(optarg);
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

    int update_covariance = ( (load_cov == 1) && (OPTION_FULL_UPDATE == 1)); //do we load the covariance ?
    struct s_pmcmc *p_pmcmc = build_pmcmc(implementation, noises_off, settings,
                                          dt, eps_abs, eps_rel,
                                          a, m_switch, m_eps, epsilon_max, is_smooth, alpha,
                                          update_covariance, J, &n_threads);
    json_decref(settings);

    transform_theta(p_pmcmc->p_best, p_pmcmc->p_data, !update_covariance);
    gsl_vector_memcpy(p_pmcmc->p_best->proposed, p_pmcmc->p_best->mean);


    n_traj = GSL_MIN(M, n_traj);
    int thin_traj = (int) ( (double) M / (double) n_traj); //the thinning interval

    pmcmc(p_pmcmc->p_best, p_pmcmc->D_J_p_X, p_pmcmc->D_J_p_X_tmp, p_pmcmc->p_par, &(p_pmcmc->D_p_hat_prev), &(p_pmcmc->D_p_hat_new), p_pmcmc->D_p_hat_best, p_pmcmc->p_like, p_pmcmc->p_data, p_pmcmc->calc, get_f_pred(implementation, noises_off), print_opt, thin_traj);

    FILE *p_file_hat = sfr_fopen(SFR_PATH, GENERAL_ID, "hat", "w", header_hat, p_pmcmc->p_data);
    print_hat(p_file_hat, p_pmcmc->D_p_hat_best, p_pmcmc->p_data);
    sfr_fclose(p_file_hat);

    // print empirical covariance
    FILE *p_file_cov = sfr_fopen(SFR_PATH, GENERAL_ID, "covariance", "w", NULL, NULL);
    print_covariance(p_file_cov, p_pmcmc->p_best->var_sampling);
    sfr_fclose(p_file_cov);

#if FLAG_VERBOSE
    print_log("clean up...\n");
#endif

    clean_pmcmc(p_pmcmc);

    return 0;
}
