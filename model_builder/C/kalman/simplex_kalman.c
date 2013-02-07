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

double f_simplex_kalman(const gsl_vector *x, void *params)
{
    struct s_kalman *p = (struct s_kalman *) params;

    transfer_estimated(p->p_best, x, p->p_data);

    back_transform_theta2par(p->p_par, p->p_best->mean, p->p_data->p_it_all, p->p_data);
    linearize_and_repeat(p->p_X, p->p_par, p->p_data, p->p_data->p_it_par_sv);
    prop2Xpop_size(p->p_X, p->p_data, PLOM_ODE);
    theta_driftIC2Xdrift(p->p_X, p->p_best->mean, p->p_data);

    double log_like = run_kalman(p->p_X, p->p_best, p->p_par, p->p_kal, p->p_data, p->calc, f_prediction_ode,  NULL, 0);

    // "-": simplex minimizes hence the "-"
    return - log_like;
}



int main(int argc, char *argv[])
{

    char ch;
    char str[STR_BUFFSIZE];

    char sfr_help_string[] =
        "PloM ksimplex\n"
        "usage:\n"
        "ksimplex [implementation] [--no_dem_sto] [--no_env_sto] [--no_drift]\n"
        "                          [-p, --path <path>] [-i, --id <integer>]\n"
        "                          [-s, --DT <float>] [--prior] [--transf]\n"
        "                          [-l, --LIKE_MIN <float>] [-S, --size <float>] [-M, --iter <integer>]\n"
        "                          [--help]\n"
        "where implementation is 'sde' (default)\n"
        "options:\n"
        "--no_dem_sto       turn off demographic stochasticity (if possible)\n"
        "--no_env_sto       turn off environmental stochasticity (if any)\n"
        "--no_drift         turn off drift (if any)\n"
        "--prior            to maximize posterior density in natural space\n"
        "--transf           to maximize posterior density in transformed space (if combined with --prior)\n"
        "-p, --path         path where the outputs will be stored\n"
        "-i, --id           general id (unique integer identifier that will be appended to the output files)\n"
        "-s, --DT           integration time step\n"
        "-l, --LIKE_MIN     particles with likelihood smaller that LIKE_MIN are considered lost\n"
        "-M, --iter         maximum number of iterations\n"
        "-S, --size         simplex size used as a stopping criteria\n"
        "-b, --no_traces    do not write the traces\n"
        "--help             print the usage on stdout\n";


    // simplex options
    M = 10;
    CONVERGENCE_STOP_SIMPLEX = 1e-6;

    // general options
    int has_dt_be_specified = 0;
    double dt_option;

    GENERAL_ID =0;
    snprintf(SFR_PATH, STR_BUFFSIZE, "%s", DEFAULT_PATH);
    J=1;
    LIKE_MIN = 1e-17;
    LOG_LIKE_MIN = log(1e-17);

    // options
    OPTION_TRAJ = 0;
    OPTION_PRIOR = 0;
    OPTION_TRANSF = 0;

    int option_no_trace = 0;
    int n_threads = 1;

    enum plom_implementations implementation;
    enum plom_noises_off noises_off = 0;


    static struct option long_options[] = {
        {"help",       no_argument,       0, 'e'},
        {"no_trace",   no_argument,  0, 'b'},

	{"no_dem_sto", no_argument,       0, 'x'},
	{"no_env_sto", no_argument,       0, 'y'},
	{"no_drift",   no_argument,       0, 'z'},

        {"path",       required_argument, 0, 'p'},
        {"id",         required_argument, 0, 'i'},

        {"prior",  no_argument, &OPTION_PRIOR,  1},
        {"transf", no_argument, &OPTION_TRANSF, 1},

        {"DT",       required_argument, 0, 's'},
        {"LIKE_MIN", required_argument, 0, 'l'},
        {"iter",     required_argument,   0, 'M'},
        {"size",     required_argument,   0, 'S'},

        {0, 0, 0, 0}
    };

    int option_index = 0;
    while ((ch = getopt_long (argc, argv, "i:l:s:p:S:M:b", long_options, &option_index)) != -1) {
        switch (ch) {
        case 0:
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

        case 'b':
            option_no_trace = 1;
            break;

        case 'p':
            snprintf(SFR_PATH, STR_BUFFSIZE, "%s", optarg);
            break;
        case 'i':
            GENERAL_ID = atoi(optarg);
            break;
        case 's':
            dt_option = atof(optarg);
            has_dt_be_specified =1;
            break;

        case 'l':
            LIKE_MIN = atof(optarg);
            LOG_LIKE_MIN = log(LIKE_MIN);
            break;
        case 'M':
            M = atoi(optarg);
            break;
        case 'S':
            CONVERGENCE_STOP_SIMPLEX = atof(optarg);
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

#if FLAG_VERBOSE
    snprintf(str, STR_BUFFSIZE, "Starting PloM ksimplex with the following options: i = %d, LIKE_MIN = %g DT = %g DELTA_STO = %g", GENERAL_ID, LIKE_MIN, DT, DELTA_STO );
    print_log(str);
#endif


    struct s_kalman *p_kalman = build_kalman(settings, implementation,  noises_off, &n_threads, OPTION_PRIOR, 0);
    json_decref(settings);

    if (OPTION_PRIOR) {
        sanitize_best_to_prior(p_kalman->p_best, p_kalman->p_data);
    }

    transform_theta(p_kalman->p_best, NULL, NULL, p_kalman->p_data, 1, 1);

    simplex(p_kalman->p_best, p_kalman->p_data, p_kalman, f_simplex_kalman, CONVERGENCE_STOP_SIMPLEX, M, option_no_trace);

#if FLAG_VERBOSE
    print_log("clean up...\n");
#endif

    clean_kalman(p_kalman, implementation, n_threads);

    return 0;
}
