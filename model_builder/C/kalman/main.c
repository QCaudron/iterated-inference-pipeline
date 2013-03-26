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

    char sfr_help_string[] =
        "PloM Kalman\n"
        "usage:\n"
        "kalman [implementation] [--no_dem_sto] [--no_env_sto] [--no_drift]\n"
        "                        [-s, --DT <float>] [--eps_abs <float>] [--eps_rel <float>]\n"
        "                        [--traj] [-p, --path <path>] [-i, --id <integer>]\n"
        "                        [-b, --no_best] [--prior] [--transf]\n"
        "                        [--help]\n"
        "where implementation is 'sde' (default)\n"
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
        "--traj             print the trajectories\n"
        "--prior            add log(prior) to the estimated loglik\n"
        "--transf           add log(JacobianDeterminant(transf)) to the estimated loglik. (combined to --prior, gives posterior density in transformed space)\n"
        "-p, --path         path where the outputs will be stored\n"
        "-b, --no_best      do not write best_<general_id>.output file\n"
        "-h, --no_hat       do not write hat_<general_id>.output file\n"
        "-r, --no_pred_res  do not write pred_res_<general_id>.output file (prediction residuals)\n"
        "-i, --id           general id (unique integer identifier that will be appended to the output files)\n"
        "-l, --LIKE_MIN     particles with likelihood smaller that LIKE_MIN are considered lost\n"
        "--help             print the usage on stdout\n";


    // general options
    GENERAL_ID =0;
    snprintf(SFR_PATH, STR_BUFFSIZE, "%s", DEFAULT_PATH);
    LIKE_MIN = 1e-17;
    LOG_LIKE_MIN = log(LIKE_MIN);
    int output_best = 1, output_hat = 1, output_pred_res =1;

    // options
    OPTION_TRAJ = 0;
    OPTION_PRIOR = 0;
    OPTION_TRANSF = 0;

    double dt = 0.0, eps_abs = PLOM_EPS_ABS, eps_rel = PLOM_EPS_REL;

    J = 1; //not an option, needed for print_X



    enum plom_implementations implementation;
    enum plom_noises_off noises_off = 0;


    static struct option long_options[] = {
        {"traj", no_argument,       &OPTION_TRAJ, 1},
	{"no_dem_sto", no_argument,       0, 'x'},
	{"no_env_sto", no_argument,       0, 'y'},
	{"no_drift",   no_argument,       0, 'z'},

	{"DT",         required_argument, 0, 's'},
	{"eps_abs",    required_argument, 0, 'v'},
	{"eps_rel",    required_argument, 0, 'w'},

        {"help",       no_argument,       0, 'e'},
        {"path",       required_argument, 0, 'p'},
        {"id",         required_argument, 0, 'i'},
	{"no_best",    no_argument,       0, 'b'},
	{"no_hat",     no_argument,       0, 'h'},
	{"no_pred_res",no_argument,       0, 'r'},


        {"traj", no_argument, &OPTION_TRAJ, 1},
        {"prior", no_argument, &OPTION_PRIOR, 1},
        {"transf", no_argument, &OPTION_TRANSF, 1},

        {"LIKE_MIN",   required_argument, 0, 'l'},

        {0, 0, 0, 0}
    };

    int option_index = 0;
    while ((ch = getopt_long (argc, argv, "xyzs:v:w:i:l:p:b", long_options, &option_index)) != -1) {
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

        case 'e':
            print_log(sfr_help_string);
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

#if FLAG_VERBOSE
    snprintf(str, STR_BUFFSIZE, "Starting Plom Kalman with the following options: i = %d, LIKE_MIN = %g", GENERAL_ID, LIKE_MIN );
    print_log(str);
#endif

    struct s_kalman *p_kalman = build_kalman(settings, implementation, noises_off, OPTION_PRIOR, 0, dt, eps_abs, eps_rel);
    json_decref(settings);

    transform_theta(p_kalman->p_best, p_kalman->p_data, 1);

#if FLAG_VERBOSE
    int64_t time_begin, time_end;
    time_begin = s_clock();
#endif

    back_transform_theta2par(p_kalman->p_par, p_kalman->p_best->mean, p_kalman->p_data->p_it_all, p_kalman->p_data);
    linearize_and_repeat(p_kalman->p_X, p_kalman->p_par, p_kalman->p_data, p_kalman->p_data->p_it_par_sv);
    prop2Xpop_size(p_kalman->p_X, p_kalman->p_data);
    theta_driftIC2Xdrift(p_kalman->p_X, p_kalman->p_best->mean, p_kalman->p_data);

    FILE *p_file_X = (OPTION_TRAJ==1) ? sfr_fopen(SFR_PATH, GENERAL_ID, "X", "w", header_X, p_kalman->p_data): NULL;
    FILE *p_file_pred_res = (output_pred_res==1) ? sfr_fopen(SFR_PATH, GENERAL_ID, "pred_res", "w", header_prediction_residuals_ekf, p_kalman->p_data): NULL;
    FILE *p_file_hat = (output_hat==1) ? sfr_fopen(SFR_PATH, GENERAL_ID, "hat", "w", header_hat, p_kalman->p_data): NULL;

    double log_like = run_kalman(p_kalman->p_X, p_kalman->p_best, p_kalman->p_par, p_kalman->p_kalman_update, p_kalman->p_data, p_kalman->calc, f_prediction_ode, p_file_X, 0, p_file_pred_res);

    if (p_file_X) {
        sfr_fclose(p_file_X);
    }

    if (p_file_pred_res) {
        sfr_fclose(p_file_pred_res);
    }

    if (p_file_hat) {
        sfr_fclose(p_file_hat);
    }


#if FLAG_VERBOSE
    time_end = s_clock();
    struct s_duration t_exec = time_exec(time_begin, time_end);
    sprintf(str, "logV: %g\t computed in:= %dd %dh %dm %gs\n", log_like, t_exec.d, t_exec.h, t_exec.m, t_exec.s);
    print_log(str);
#endif

    if(output_best) {
        FILE *p_file_best = sfr_fopen(SFR_PATH, GENERAL_ID, "best", "w", header_best, p_kalman->p_data);
        print_best(p_file_best, 0, p_kalman->p_best, p_kalman->p_data, log_like);
        sfr_fclose(p_file_best);
    }

#if FLAG_VERBOSE
    print_log("clean up...\n");
#endif

    clean_kalman(p_kalman);

    return 0;
}
