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
    char plom_help_string[] =
        "PLOM pMCMC\n"
        "usage:\n"
        "pmcmc [implementation] [--no_dem_sto] [--no_white_noise] [--no_diff]\n"
        "                [-s, --DT <float || 0.25 day>] [--eps_abs <float || 1e-6>] [--eps_rel <float || 1e-3>]\n"
        "                [-g, --freeze_forcing <float>]\n"
        "                [--full] [-n, --n_traj <int || 1000>] [--acc] [-p, --path <path>] [-i, --id <integer || 0>] [-N, --n_thread <integer || 1>]\n"
        "                [-l, --LIKE_MIN <float || 1e-17>] [-J <integer || 1>] [-M, --iter <integer || 10>]\n"
        "                [-a --cooling <float || 0.999>] [-S --switch <int || 5*n_par_fitted^2 >] "
        "                [-E --epsilon <int || 50>] [--epsilon_max <float || 50.0>] [--smooth] [--alpha <float || 0.02>]"
        "                [-Z, --zmq] [-c, --chunk <integer>]\n"
	"                [-q, --quiet] [-P, --pipe]"
        "                [-h, --help]\n"
        "where implementation is 'ode', 'sde' or 'psr' (default)\n"
        "options:\n"
        "-q, --quiet        no verbosity\n"
        "-P, --pipe         pipe mode (echo theta.json on stdout)\n"
        "\n"
        "--no_dem_sto       turn off demographic stochasticity (if possible)\n"
        "--no_white_noise       turn off environmental stochasticity (if any)\n"
        "--no_diff         turn off drift (if any)\n"
        "\n"
        "-s, --DT           Initial integration time step\n"
        "--eps_abs          Absolute error for adaptive step-size contro\n"
        "--eps_rel          Relative error for adaptive step-size contro\n"
        "-g, --freeze_forcing freeze the metadata to their value at the specified time\n"
        "\n"
        "--full             full update MVN mode\n"
        "-a, --cooling      cooling factor for sampling covariance live tuning\n"
        "-S, --switch       select switching iteration from initial covariance to empirical one\n"
        "-E, --epsilon      select number of burnin iterations before tuning epsilon\n"
        "--epsilon_max      maximum value allowed for epsilon\n"
        "--smooth           tune epsilon with the value of the acceptance rate obtained with exponential smoothing\n"
        "--alpha            smoothing factor of exponential smoothing used to compute the smoothed acceptance rate (low values increase degree of smoothing)\n"
        "\n"
        "-n, --n_traj       number of trajectories stored\n"
        "--acc              print the acceptance rate\n"
        "-p, --path         path where the outputs will be stored\n"
        "-i, --id           general id (unique integer identifier that will be appended to the output files)\n"
        "-N, --n_thread     number of threads to be used\n"
        "-s, --DT           integration time step\n"
        "-l, --LIKE_MIN     particles with likelihood smaller that LIKE_MIN are considered lost\n"
        "-J                 number of particles\n"
        "-M, --iter         number of pMCMC iterations\n"
        "-Z, --zmq          dispatch particles across machines using a zeromq pipeline\n"
        "-c, --chunk        number of particles send to each machine\n"
	"-o, --nb_obs       number of observations to be fitted (for tempering)"
        "-h, --help         print the usage on stdout\n";

    double dt = 0.0, eps_abs = PLOM_EPS_ABS, eps_rel = PLOM_EPS_REL;
    int m_switch = -1;
    int m_eps = 50;
    double a = 0.999;
    double epsilon_max = 50.0;
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
    
    int n_threads = 1;

    OPTION_PIPELINE = 0;
    OPTION_FULL_UPDATE = 0;
    int nb_obs = -1;
    double freeze_forcing = -1.0;

    enum plom_implementations implementation;
    enum plom_noises_off noises_off = 0;
    enum plom_print print_opt = 0;

    while (1) {
        static struct option long_options[] =
            {
                /* These options set a flag. */
                {"full", no_argument, &OPTION_FULL_UPDATE, 1},

                /* These options don't set a flag We distinguish them by their indices (that are also the short option names). */
                {"acc",        no_argument,       0, 'r'},
                {"no_dem_sto", no_argument,       0, 'x'},
                {"no_white_noise", no_argument,       0, 'y'},
                {"no_diff",   no_argument,       0, 'z'},

                {"DT",         required_argument, 0, 's'},
                {"eps_abs",    required_argument, 0, 'v'},
                {"eps_rel",    required_argument, 0, 'w'},

                {"freeze_forcing", required_argument, 0, 'g'},
                {"n_traj",     required_argument, 0, 'n'},

                {"switch",     required_argument,   0, 'S'},
                {"epsilon",     required_argument,   0, 'E'},
                {"cooling",     required_argument,   0, 'a'},
                {"smooth",    no_argument, &is_smooth, 1},
                {"epsilon_max", required_argument, 0, 'f'},
                {"alpha",    required_argument, 0, 'G'},

                {"help", no_argument,  0, 'h'},
                {"path",    required_argument, 0, 'p'},
                {"id",    required_argument, 0, 'i'},
                {"n_thread",  required_argument,       0, 'N'},

                {"LIKE_MIN",     required_argument,   0, 'l'},
                {"iter",     required_argument,   0, 'M'},
                {"zmq",     no_argument,   0, 'Z'},
                {"chunk",     required_argument,   0, 'c'},
		{"nb_obs", required_argument,  0, 'o'},

                {"quiet",  no_argument,       0, 'q'},
                {"pipe",  no_argument,       0, 'P'},

                {0, 0, 0, 0}
            };
        /* getopt_long stores the option index here. */
        int option_index = 0;

        ch = getopt_long (argc, argv, "qPhrxyzs:v:w:i:J:l:M:p:c:N:ZS:E:a:f:G:n:o:g:", long_options, &option_index);

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
        case 'g':
            freeze_forcing = atof(optarg);
            break;
	case 'o':
	    nb_obs = atoi(optarg);
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
        case 'G':
            alpha = atof(optarg);
            break;
        case 'h':
            print_log(plom_help_string);
            return 1;
        case 'Z':
            OPTION_PIPELINE = 1;
            break;
        case 'p':
            snprintf(SFR_PATH, STR_BUFFSIZE, "%s", optarg);
            break;
        case 'N':
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
        case 'r':
	    print_opt |= PLOM_PRINT_ACC;
            break;
        case 'n':
            n_traj = atoi(optarg);
            break;

        case 'q':
	    print_opt |= PLOM_QUIET;
            break;
        case 'P':
	    print_opt |= PLOM_PIPE | PLOM_QUIET;
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
            print_log(plom_help_string);
            return 1;
        }
    }

    plom_unlink_done(SFR_PATH, GENERAL_ID);
    json_t *settings = load_settings(PATH_SETTINGS);

    json_t *theta = load_json();
    struct s_pmcmc *p_pmcmc = build_pmcmc(theta, implementation, noises_off, settings,
                                          dt, eps_abs, eps_rel, freeze_forcing,
                                          a, m_switch, m_eps, epsilon_max, is_smooth, alpha,
                                          J, &n_threads, nb_obs);
    json_decref(settings);

    gsl_vector_memcpy(p_pmcmc->p_best->proposed, p_pmcmc->p_best->mean);

    n_traj = GSL_MIN(M, n_traj);
    int thin_traj = (int) ( (double) M / (double) n_traj); //the thinning interval

    if(n_traj>0){       
	print_opt |= PLOM_PRINT_X_SMOOTH;
    }
    
    int64_t time_begin, time_end;
    if (!(print_opt & PLOM_QUIET)) {
	sprintf(str, "Starting plom-pmcmc with the following options: i = %d, J = %d, LIKE_MIN = %g, M = %d, n_threads = %d", GENERAL_ID, J, LIKE_MIN, M, n_threads);
	print_log(str);
	
	time_begin = s_clock();
    }

    pmcmc(p_pmcmc->p_best, p_pmcmc->D_J_p_X, p_pmcmc->D_J_p_X_tmp, p_pmcmc->p_par, &(p_pmcmc->D_p_hat_prev), &(p_pmcmc->D_p_hat_new), p_pmcmc->D_p_hat_best, p_pmcmc->p_like, p_pmcmc->p_data, p_pmcmc->calc, get_f_pred(implementation, noises_off), print_opt, thin_traj);


    //TODO compute quantile online (https://github.com/plom-io/plom-sfi/issues/9)
    //FILE *p_file_hat = plom_fopen(SFR_PATH, GENERAL_ID, "hat", "w", header_hat, p_pmcmc->p_data);
    //print_hat(p_file_hat, p_pmcmc->D_p_hat_best, p_pmcmc->p_data);
    //plom_fclose(p_file_hat);

    // print empirical covariance
    FILE *p_file_cov = plom_fopen(SFR_PATH, GENERAL_ID, "covariance", "w", header_covariance, p_pmcmc->p_data);
    print_covariance(p_file_cov, p_pmcmc->p_best->var_sampling, p_pmcmc->p_data);
    plom_fclose(p_file_cov);

    plom_print_done(theta, p_pmcmc->p_data, p_pmcmc->p_best, SFR_PATH, GENERAL_ID, print_opt);

    if (!(print_opt & PLOM_QUIET)) {
	time_end = s_clock();
	struct s_duration t_exec = time_exec(time_begin, time_end);
	snprintf(str, STR_BUFFSIZE, "Done in:= %dd %dh %dm %gs", t_exec.d, t_exec.h, t_exec.m, t_exec.s);
	print_log(str);
    }


    if (!(print_opt & PLOM_QUIET)) {
	print_log("clean up...");
    }

    json_decref(theta);

    clean_pmcmc(p_pmcmc);

    return 0;
}
