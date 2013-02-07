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

#include "kalman.h"

int main(int argc, char *argv[])
{
    char ch;
    char str[STR_BUFFSIZE];

    /* set default values for the options */
    char sfr_help_string[] =
        "PloM kmcmc\n"
        "usage:\n"
        "kmcmc [implementation] [--no_dem_sto] [--no_env_sto] [--no_drift]\n"
        "                       [--full] [--traj] [-p, --path <path>] [-i, --id <integer>] [-P, --N_THREAD <integer>]\n"
        "                       [-s, --DT <float>] [-l, --LIKE_MIN <float>] [-M, --iter <integer>]\n"
        "                       [-c --cov] [-a --cooling <float>] [-S --switch <int>]\n"
        "                       [--help]\n"
        "where implementation is 'sde' (default)\n"
        "options:\n"
        "--no_dem_sto       turn off demographic stochasticity (if possible)\n"
        "--no_env_sto       turn off environmental stochasticity (if any)\n"
        "--no_drift         turn off drift (if any)\n"
        "--full             full update MVN mode\n"
        "--traj             print the trajectories\n"
        "-c, --cov          load an initial covariance from the settings\n"
        "-p, --path         path where the outputs will be stored\n"
        "-i, --id           general id (unique integer identifier that will be appended to the output files)\n"
        "-P, --N_THREAD     number of threads to be used (default to the number of cores)\n"
        "-s, --DT           integration time step\n"
        "-l, --LIKE_MIN     particles with likelihood smaller that LIKE_MIN are considered lost\n"
        "-M, --iter         number of pMCMC iterations\n"
        "-a, --cooling      cooling rate for sampling covariance live tuning\n"
        "-S, --switch       select switching iteration from initial covariance to empirical one\n"
        "--help             print the usage on stdout\n";

    int has_dt_be_specified = 0;
    double dt_option;
    int load_cov = 0;
    int m_switch = -1;
    double a = 0.999;


    GENERAL_ID =0;
    snprintf(SFR_PATH, STR_BUFFSIZE, "%s", DEFAULT_PATH);
    LIKE_MIN = 1e-17;
    LOG_LIKE_MIN = log(1e-17);
    M = 10;
    OPTION_FULL_UPDATE = 0;
    OPTION_PRIOR = 0;
    OPTION_TRANSF = 0;
    OPTION_TRAJ = 0;

    int n_threads = 1;
    J = 1; //not an option, needed for print_X


    enum plom_implementations implementation;
    enum plom_noises_off noises_off = 0;


    while (1) {
        static struct option long_options[] =
            {
                /* These options set a flag. */
                {"traj", no_argument,       &OPTION_TRAJ, 1},
                {"full", no_argument, &OPTION_FULL_UPDATE, 1},
                {"prior",       no_argument, &OPTION_PRIOR,       1},
                {"transf",      no_argument, &OPTION_TRANSF,      1},

                /* These options don't set a flag We distinguish them by their indices (that are also the short option names). */
		{"no_dem_sto", no_argument,       0, 'x'},
		{"no_env_sto", no_argument,       0, 'y'},
		{"no_drift",   no_argument,       0, 'z'},

                {"cov",         no_argument, 0, 'c'},

                {"help",        no_argument,        0, 'e'},
                {"path",        required_argument,  0, 'p'},
                {"id",          required_argument,  0, 'i'},
                {"N_THREAD",    required_argument,  0, 'P'},

                {"DT",          required_argument,   0, 's'},
                {"LIKE_MIN",    required_argument,   0, 'l'},
                {"iter",        required_argument,   0, 'M'},
                {"switch",      required_argument,   0, 'S'},
                {"cooling",     required_argument,   0, 'a'},

                {0, 0, 0, 0}
            };
        /* getopt_long stores the option index here. */
        int option_index = 0;

        ch = getopt_long (argc, argv, "ci:l:M:p:s:P:S:a:", long_options, &option_index);

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

        case 'e':
            print_log(sfr_help_string);
            return 1;

        case 'c':
            load_cov = 1;
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
        case 'l':
            LIKE_MIN = atof(optarg);
            LOG_LIKE_MIN = log(LIKE_MIN);
            break;
        case 'M':
            M = atoi(optarg);
            break;
        case 's':
            dt_option = atof(optarg);
            has_dt_be_specified =1;
            break;
        case 'a':
            a = atof(optarg);
            break;
        case 'S':
            m_switch = atoi(optarg);
            break;

        case '?':
            /* getopt_long already printed an error message. */
            break;

        default:
            snprintf(str, STR_BUFFSIZE, "Unknown option '-%c'\n", optopt);
            print_err(str);
            return 1;
        }
    }
    argc -= optind;
    argv += optind;


    if(argc == 0) {
	implementation = PLOM_ODE; //with Kalman the SDE uses f_pred of PLOM_ODE (OK will do better)...
    } else {
        if (!strcmp(argv[0], "sde")) {
            implementation = PLOM_ODE;
        } else {
            print_log(sfr_help_string);
            return 1;
        }
    }


    json_t *settings = load_settings(PATH_SETTINGS);

    if (has_dt_be_specified) {
        DT = dt_option;
    }

    //IMPORTANT: update DELTA_STO so that DT = 1.0/DELTA_STO
    DELTA_STO = round(1.0/DT);
    DT = 1.0/ ((double) DELTA_STO);

    int update_covariance = ( (load_cov == 1) && (OPTION_FULL_UPDATE == 1)); //do we load the covariance ?
    struct s_kalman *p_kalman = build_kalman(settings, implementation,  noises_off, &n_threads, 1, update_covariance);
    json_decref(settings);

    sanitize_best_to_prior(p_kalman->p_best, p_kalman->p_data);

    transform_theta(p_kalman->p_best, NULL, NULL, p_kalman->p_data, 1, !update_covariance);
    gsl_vector_memcpy(p_kalman->p_best->proposed, p_kalman->p_best->mean);

    struct s_likelihood *p_like = build_likelihood();
    struct s_pmcmc_calc_data *p_mcmc_calc_data = build_pmcmc_calc_data(p_kalman->p_best, a, m_switch,1);

    p_kalman->calc[0]->method_specific_shared_data = p_mcmc_calc_data; //needed to use ran_proposal_sequential in pmcmc

    kmcmc(p_kalman, p_like, p_mcmc_calc_data, get_f_pred(implementation, noises_off), implementation);

    // print empirical covariance
    if (OPTION_FULL_UPDATE) {
        FILE *p_file_cov = sfr_fopen(SFR_PATH, GENERAL_ID, "covariance", "w", NULL, NULL);
        print_covariance(p_file_cov, p_kalman->p_best->var_sampling);
        sfr_fclose(p_file_cov);
    }

#if FLAG_VERBOSE
    print_log("clean up...\n");
#endif

    clean_pmcmc_calc_data(p_mcmc_calc_data);
    clean_likelihood(p_like);
    clean_kalman(p_kalman, implementation, n_threads);

    return 0;
}
