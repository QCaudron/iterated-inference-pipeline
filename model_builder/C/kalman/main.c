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

    char plom_help_string[] =
        "PLOM Kalman\n"
        "usage:\n"
        "kalman [implementation] [--no_dem_sto] [--no_white_noise] [--no_diff]\n"
        "                        [-s, --DT <float>] [--eps_abs <float>] [--eps_rel <float>]\n"
        "                        [-g, --freeze_forcing <float>]\n"
        "                        [-r, --traj] [-p, --path <path>] [-i, --id <integer>]\n"
        "                        [-b, --no_trace] [-e, --no_hat] [--prior] [--transf]\n"
	"                        [-q, --quiet] [-P, --pipe]"
        "                        [-h, --help]\n"
        "where implementation is 'sde' (default)\n"
        "options:\n"
	"\n"
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
        "-g, --freeze_forcing  freeze the metadata to their value at the specified time\n"
	"\n"
        "-r, --traj         print the trajectories\n"
        "--prior            add log(prior) to the estimated loglik\n"
        "--transf           add log(JacobianDeterminant(transf)) to the estimated loglik. (combined to --prior, gives posterior density in transformed space)\n"
        "-p, --path         path where the outputs will be stored\n"
        "-b, --no_trace     do not write trace_<id>.output file\n"
        "-h, --no_hat       do not write hat_<general_id>.output file\n"
        "-d, --no_pred_res  do not write pred_res_<general_id>.output file (prediction residuals)\n"
        "-i, --id           general id (unique integer identifier that will be appended to the output files)\n"
        "-l, --LIKE_MIN     particles with likelihood smaller that LIKE_MIN are considered lost\n"
	"-o, --nb_obs       number of observations to be fitted (for tempering)"
        "-h, --help         print the usage on stdout\n";


    // general options
    GENERAL_ID =0;
    snprintf(SFR_PATH, STR_BUFFSIZE, "%s", DEFAULT_PATH);
    LIKE_MIN = 1e-17;
    LOG_LIKE_MIN = log(LIKE_MIN);
    enum plom_print print_opt = PLOM_PRINT_BEST | PLOM_PRINT_HAT | PLOM_PRINT_PRED_RES;

    // options
    OPTION_PRIOR = 0;
    OPTION_TRANSF = 0;
    int nb_obs = -1;
    double freeze_forcing = -1.0;

    double dt = 0.0, eps_abs = PLOM_EPS_ABS, eps_rel = PLOM_EPS_REL;

    J = 1; //not an option, needed for print_X



    enum plom_implementations implementation;
    enum plom_noises_off noises_off = 0;


    static struct option long_options[] = {
	{"traj",       no_argument,       0, 'r'},
	{"no_dem_sto", no_argument,       0, 'x'},
	{"no_white_noise", no_argument,       0, 'y'},
	{"no_diff",   no_argument,       0, 'z'},

	{"DT",         required_argument, 0, 's'},
	{"eps_abs",    required_argument, 0, 'v'},
	{"eps_rel",    required_argument, 0, 'w'},
	{"freeze_forcing", required_argument, 0, 'g'},
        {"help",       no_argument,       0, 'h'},
        {"path",       required_argument, 0, 'p'},
        {"id",         required_argument, 0, 'i'},
	{"no_trace",   no_argument,       0, 'b'},
	{"no_hat",     no_argument,       0, 'e'},
	{"no_pred_res",no_argument,       0, 'd'},

        {"prior", no_argument, &OPTION_PRIOR, 1},
        {"transf", no_argument, &OPTION_TRANSF, 1},
	{"nb_obs", required_argument,  0, 'o'},

	{"quiet",  no_argument,       0, 'q'},
	{"pipe",  no_argument,       0, 'P'},

        {"LIKE_MIN",   required_argument, 0, 'l'},

        {0, 0, 0, 0}
    };

    int option_index = 0;
    while ((ch = getopt_long (argc, argv, "qPxyzs:v:w:i:l:p:rbhedo:g:", long_options, &option_index)) != -1) {
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
	case 'o':
	    nb_obs = atoi(optarg);
            break;

        case 'g':
            freeze_forcing = atof(optarg);
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

        case 'p':
            snprintf(SFR_PATH, STR_BUFFSIZE, "%s", optarg);
            break;
        case 'i':
            GENERAL_ID = atoi(optarg);
            break;
        case 'l':
            LIKE_MIN = atof(optarg);
            LOG_LIKE_MIN = log(LIKE_MIN);
            break;
        case 'r':
	    print_opt |= PLOM_PRINT_X;
            break;
        case 'b':
	    print_opt &= ~PLOM_PRINT_BEST;
            break;
        case 'e':
	    print_opt &= ~PLOM_PRINT_HAT;
            break;
        case 'd':
	    print_opt &= ~PLOM_PRINT_PRED_RES;
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

    if (!(print_opt & PLOM_QUIET)) {
	snprintf(str, STR_BUFFSIZE, "Starting plom-Kalman with the following options: i = %d, LIKE_MIN = %g", GENERAL_ID, LIKE_MIN);
	print_log(str);
    }

    json_t *theta = load_json();
    struct s_kalman *p_kalman = build_kalman(theta, settings, implementation, noises_off, OPTION_PRIOR, dt, eps_abs, eps_rel, freeze_forcing, nb_obs);
    json_decref(settings);

    int64_t time_begin, time_end;
    if (!(print_opt & PLOM_QUIET)) {
	time_begin = s_clock();
    }

    back_transform_theta2par(p_kalman->p_par, p_kalman->p_best->mean, p_kalman->p_data->p_it_all, p_kalman->p_data);
    linearize_and_repeat(p_kalman->p_X, p_kalman->p_par, p_kalman->p_data, p_kalman->p_data->p_it_par_sv);
    prop2Xpop_size(p_kalman->p_X, p_kalman->p_data, p_kalman->calc[0]);
    theta_driftIC2Xdrift(p_kalman->p_X, p_kalman->p_best->mean, p_kalman->p_data);

    FILE *p_file_X = (print_opt & PLOM_PRINT_X) ? plom_fopen(SFR_PATH, GENERAL_ID, "X", "w", header_X, p_kalman->p_data): NULL;
    FILE *p_file_hat = (print_opt & PLOM_PRINT_HAT) ? plom_fopen(SFR_PATH, GENERAL_ID, "hat", "w", header_hat, p_kalman->p_data): NULL;
    FILE *p_file_pred_res = (print_opt & PLOM_PRINT_PRED_RES) ? plom_fopen(SFR_PATH, GENERAL_ID, "pred_res", "w", header_prediction_residuals_ekf, p_kalman->p_data): NULL;

    double log_like = run_kalman(p_kalman->p_X, p_kalman->p_best, p_kalman->p_par, p_kalman->p_kalman_update, p_kalman->p_data, p_kalman->calc, f_prediction_ode, 0, p_file_X, p_file_hat, p_file_pred_res, print_opt);


    if (print_opt & PLOM_PRINT_X) {
        plom_fclose(p_file_X);
    }

    if (print_opt & PLOM_PRINT_HAT) {
        plom_fclose(p_file_hat);
    }

    if (print_opt & PLOM_PRINT_PRED_RES) {
        plom_fclose(p_file_pred_res);
    }

    if (!(print_opt & PLOM_QUIET)) {
	time_end = s_clock();
	struct s_duration t_exec = time_exec(time_begin, time_end);
	sprintf(str, "logV: %g", log_like);
	print_log(str);

	sprintf(str, "Done in:= %dd %dh %dm %gs", t_exec.d, t_exec.h, t_exec.m, t_exec.s);
	print_log(str);
    }

    if (print_opt & PLOM_PRINT_BEST) {
        FILE *p_file_trace = plom_fopen(SFR_PATH, GENERAL_ID, "trace", "w", header_trace, p_kalman->p_data);
        print_trace(p_file_trace, 0, p_kalman->p_best, p_kalman->p_data, log_like);
        plom_fclose(p_file_trace);
    }

    plom_print_done(theta, p_kalman->p_data, p_kalman->p_best, SFR_PATH, GENERAL_ID, print_opt);

    if (!(print_opt & PLOM_QUIET)) {
	print_log("clean up...");
    }

    json_decref(theta);
    clean_kalman(p_kalman);

    return 0;
}
