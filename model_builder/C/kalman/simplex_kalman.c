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
    double log_like;

    transfer_estimated(p->p_best, x, p->p_data);

    if(plom_check_IC_assign_theta_remainder(p->p_best->mean, p->p_data)){
#if FLAG_WARNING
        print_warning("IC constraint has not been respected: pop_IC>pop_size at t=0 minimal likelihood value has been assigned");
#endif
	log_like = p->smallest_log_like;
    } else {
	back_transform_theta2par(p->p_par, p->p_best->mean, p->p_data->p_it_all, p->p_data);
	linearize_and_repeat(p->p_X, p->p_par, p->p_data, p->p_data->p_it_par_sv);
	prop2Xpop_size(p->p_X, p->p_data, p->calc[0]);
	theta_driftIC2Xdrift(p->p_X, p->p_best->mean, p->p_data);

	//reset dt
	p->p_X->dt = p->p_X->dt0;

	log_like = run_kalman(p->p_X, p->p_best, p->p_par, p->p_kalman_update, p->p_data, p->calc, f_prediction_ode,  0, NULL, NULL, NULL, 0);    
    }

    // "-": simplex minimizes hence the "-"
    return - log_like;
}


int main(int argc, char *argv[])
{
    char ch;
    char str[STR_BUFFSIZE];

    char plom_help_string[] =
        "PLOM ksimplex\n"
        "usage:\n"
        "ksimplex [implementation] [--no_dem_sto] [--no_white_noise] [--no_diff]\n"
        "                          [-s, --DT <float>] [--eps_abs <float>] [--eps_rel <float>]\n"
        "                          [-p, --path <path>] [-i, --id <integer>]\n"
        "                          [-g, --freeze_forcing <float>]\n"
        "                          [--prior] [--transf]\n"
        "                          [-l, --LIKE_MIN <float>] [-S, --size <float>] [-M, --iter <integer>]\n"
	"                          [-q, --quiet] [-P, --pipe]"
        "                          [-h, --help]\n"
        "where implementation is 'sde' (default)\n"
        "options:\n"
	"\n"
        "-q, --quiet          no verbosity\n"
        "-P, --pipe           pipe mode (echo theta.json on stdout)\n"
	"\n"
        "--no_dem_sto       turn off demographic stochasticity (if possible)\n"
        "--no_white_noise       turn off environmental stochasticity (if any)\n"
        "--no_diff         turn off drift (if any)\n"
	"\n"
        "-s, --DT           Initial integration time step\n"
	"--eps_abs          Absolute error for adaptive step-size control\n"
	"--eps_rel          Relative error for adaptive step-size control\n"
        "-g, --freeze_forcing  freeze the metadata to their value at the specified time\n"
	"\n"
        "--prior            to maximize posterior density in natural space\n"
        "--transf           to maximize posterior density in transformed space (if combined with --prior)\n"
        "-p, --path         path where the outputs will be stored\n"
        "-i, --id           general id (unique integer identifier that will be appended to the output files)\n"
        "-l, --LIKE_MIN     particles with likelihood smaller that LIKE_MIN are considered lost\n"
        "-M, --iter         maximum number of iterations\n"
        "-S, --size         simplex size used as a stopping criteria\n"
        "-b, --no_traces    do not write the traces\n"
	"-o, --nb_obs       number of observations to be fitted (for tempering)"
        "--help             print the usage on stdout\n";

    // simplex options
    M = 10;
    CONVERGENCE_STOP_SIMPLEX = 1e-6;

    // general options
    GENERAL_ID =0;
    snprintf(SFR_PATH, STR_BUFFSIZE, "%s", DEFAULT_PATH);
    J=1;
    LIKE_MIN = 1e-17;
    LOG_LIKE_MIN = log(1e-17);
    int nb_obs = -1;
    double freeze_forcing = -1.0;

    // options
    OPTION_PRIOR = 0;
    OPTION_TRANSF = 0;

    double dt = 0.0, eps_abs = PLOM_EPS_ABS, eps_rel = PLOM_EPS_REL;

    enum plom_print print_opt = PLOM_PRINT_BEST;

    enum plom_implementations implementation;
    enum plom_noises_off noises_off = 0;


    static struct option long_options[] = {
        {"help",       no_argument,       0, 'h'},
        {"no_trace",   no_argument,  0, 'b'},

	{"no_dem_sto", no_argument,       0, 'x'},
	{"no_white_noise", no_argument,       0, 'y'},
	{"no_diff",   no_argument,       0, 'z'},


	{"DT",         required_argument, 0, 's'},
	{"eps_abs",    required_argument, 0, 'v'},
	{"eps_rel",    required_argument, 0, 'w'},

	{"freeze_forcing", required_argument, 0, 'g'},

        {"path",       required_argument, 0, 'p'},
        {"id",         required_argument, 0, 'i'},

        {"prior",  no_argument, &OPTION_PRIOR,  1},
        {"transf", no_argument, &OPTION_TRANSF, 1},

        {"LIKE_MIN", required_argument, 0, 'l'},
        {"iter",     required_argument,   0, 'M'},
        {"size",     required_argument,   0, 'S'},
	{"nb_obs", required_argument,  0, 'o'},

	{"quiet",  no_argument,       0, 'q'},
	{"pipe",  no_argument,       0, 'P'},

        {0, 0, 0, 0}
    };

    int option_index = 0;
    while ((ch = getopt_long (argc, argv, "qPhxyzs:v:w:i:l:p:S:M:o:bg:", long_options, &option_index)) != -1) {
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
        case 's':
            dt = atof(optarg);
            break;
        case 'v':
            eps_abs = atof(optarg);
            break;
        case 'w':
            eps_rel = atof(optarg);
            break;
        case 'h':
            print_log(plom_help_string);
            return 1;
        case 'b':
            print_opt &= ~PLOM_PRINT_BEST;
            break;
        case 'p':
            snprintf(SFR_PATH, STR_BUFFSIZE, "%s", optarg);
            break;
        case 'i':
            GENERAL_ID = atoi(optarg);
            break;
        case 'g':
            freeze_forcing = atof(optarg);
            break;
	case 'o':
	    nb_obs = atoi(optarg);
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
	implementation = PLOM_ODE; //with Kalman the SDE uses f_pred of PLOM_ODE (OK will do better)...
    } else {
        if (!strcmp(argv[0], "sde")) {
            implementation = PLOM_ODE;
        } else {
            print_log(plom_help_string);
            return 1;
        }
    }
    plom_unlink_done(SFR_PATH, GENERAL_ID);
    json_t *settings = load_settings(PATH_SETTINGS);

    int64_t time_begin, time_end;
    if (!(print_opt & PLOM_QUIET)) {
	snprintf(str, STR_BUFFSIZE, "Starting plom-ksimplex with the following options: i = %d, LIKE_MIN = %g", GENERAL_ID, LIKE_MIN);
	print_log(str);
	time_begin = s_clock();
    }

    json_t *theta = load_json();
    struct s_kalman *p_kalman = build_kalman(theta, settings, implementation, noises_off, OPTION_PRIOR, dt, eps_abs, eps_rel, freeze_forcing, nb_obs);
    json_decref(settings);

    simplex(p_kalman->p_best, p_kalman->p_data, p_kalman, f_simplex_kalman, CONVERGENCE_STOP_SIMPLEX, M, print_opt);

    if (!(print_opt & PLOM_QUIET)) {
	time_end = s_clock();
	struct s_duration t_exec = time_exec(time_begin, time_end);
	snprintf(str, STR_BUFFSIZE, "Done in:= %dd %dh %dm %gs", t_exec.d, t_exec.h, t_exec.m, t_exec.s);
	print_log(str);
    }

    plom_print_done(theta, p_kalman->p_data, p_kalman->p_best, SFR_PATH, GENERAL_ID, print_opt);

    if (!(print_opt & PLOM_QUIET)) {
	print_log("clean up...");
    }

    json_decref(theta);

    clean_kalman(p_kalman);

    return 0;
}
