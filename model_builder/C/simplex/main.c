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

#include "simplex.h"

int main(int argc, char *argv[])
{
    char ch;
    char str[STR_BUFFSIZE];

    /* set default values for the options */

    char plom_help_string[] =
        "PLOM Simplex\n"
        "usage:\n"
        "simplex [implementation] [-p, --path <path>] [-i, --id <integer>] [-q, --least_square]\n"
        "                         [-s, --DT <float>] [--eps_abs <float>] [--eps_rel <float>]\n"
        "                         [-l, --LIKE_MIN <float>] [-S, --size <float>] [-M, --iter <integer>]  [-o, --nb_obs <integer>] [--prior]\n"
	"                         [-b, --no_trace]\n"
        "                         [-g, --freeze_forcing <float>]\n"
	"                         [-q, --quiet] [-P, --pipe]"
        "                         [-h, --help]\n"
        "where implementation is 'ode' (default)\n"
        "options:\n"
	"\n"
        "-q, --quiet          no verbosity\n"
        "-P, --pipe           pipe mode (echo theta.json on stdout)\n"
	"\n"
        "-s, --DT             Initial integration time step\n"
	"--eps_abs            Absolute error for adaptive step-size contro\n"
	"--eps_rel            Relative error for adaptive step-size contro\n"
        "-g, --freeze_forcing freeze the metadata to their value at the specified time\n"
	"\n" 
	"-q, --least_square   optimize the sum of square instead of the likelihood\n"
        "-p, --path           path where the outputs will be stored\n"
        "-i, --id             general id (unique integer identifier that will be appended to the output files)\n"
        "-l, --LIKE_MIN       likelihood smaller that LIKE_MIN are considered 0.0\n"
        "-M, --iter           maximum number of iterations\n"
        "-S, --size           simplex size used as a stopping criteria\n"
	"-b, --no_trace       do not write trace_<id>.output file\n"
        "-o, --nb_obs         number of observations to be fitted (for tempering)"
        "--prior              add log(prior) to the estimated log likelihood\n"
        "h, --help            print the usage on stdout\n";


    int M = 10;
    double CONVERGENCE_STOP_SIMPLEX = 1e-6;

    enum plom_implementations implementation;
    enum plom_noises_off noises_off = PLOM_NO_DEM_STO | PLOM_NO_ENV_STO | PLOM_NO_DRIFT;

    enum plom_print print_opt = PLOM_PRINT_BEST;

    double dt = 0.0, eps_abs = PLOM_EPS_ABS, eps_rel = PLOM_EPS_REL;

    GENERAL_ID =0;
    snprintf(SFR_PATH, STR_BUFFSIZE, "%s", DEFAULT_PATH);
    J=1;
    LIKE_MIN = 1e-17;
    LOG_LIKE_MIN = log(1e-17);
    OPTION_LEAST_SQUARE = 0;
    OPTION_PRIOR = 0;
    int nb_obs = -1;
    double freeze_forcing = -1.0;

    while (1) {
        static struct option long_options[] =
            {
                /* These options don't set a flag We distinguish them by their indices (that are also the short option names). */

		{"DT",         required_argument, 0, 's'},
		{"eps_abs",    required_argument, 0, 'v'},
		{"eps_rel",    required_argument, 0, 'w'},

                {"freeze_forcing", required_argument, 0, 'g'},
                {"help", no_argument,  0, 'h'},
                {"least_square", no_argument,  0, 'Q'},
                {"no_trace", no_argument,  0, 'b'},
                {"prior", no_argument, &OPTION_PRIOR, 1},
                {"path",    required_argument, 0, 'p'},
                {"id",    required_argument, 0, 'i'},
                {"LIKE_MIN",     required_argument,   0, 'l'},
                {"iter",     required_argument,   0, 'M'},
                {"size",     required_argument,   0, 'S'},
		{"nb_obs", required_argument,  0, 'o'},

                {"quiet",  no_argument,       0, 'q'},
                {"pipe",  no_argument,       0, 'P'},

                {0, 0, 0, 0}
            };
        /* getopt_long stores the option index here. */
        int option_index = 0;

        ch = getopt_long (argc, argv, "qPhs:v:w:Qi:l:M:S:p:o:bg:", long_options, &option_index);
	
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
        case 'b':
            print_opt &= ~PLOM_PRINT_BEST;
            break;
	case 'o':
	    nb_obs = atoi(optarg);
            break;
        case 'Q':
	    OPTION_LEAST_SQUARE = 1;
            break;
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
	implementation = PLOM_ODE;
	noises_off = noises_off | PLOM_NO_DEM_STO| PLOM_NO_ENV_STO | PLOM_NO_DRIFT;
    } else {
        if (!strcmp(argv[0], "ode")) {
            implementation = PLOM_ODE;
	    noises_off = noises_off | PLOM_NO_DEM_STO| PLOM_NO_ENV_STO | PLOM_NO_DRIFT;
	} else {
            print_log(plom_help_string);
            return 1;
        }
    }

    int64_t time_begin, time_end;
    if (!(print_opt & PLOM_QUIET)) {
	sprintf(str, "Starting plom-simplex with the following options: i = %d, LIKE_MIN = %g, M = %d, CONVERGENCE_STOP_SIMPLEX = %g", GENERAL_ID, LIKE_MIN, M, CONVERGENCE_STOP_SIMPLEX);
	print_log(str);
	time_begin = s_clock();
    }

    plom_unlink_done(SFR_PATH, GENERAL_ID);
    json_t *theta = load_json();
    struct s_simplex *p_simplex = build_simplex(theta, implementation, noises_off, GENERAL_ID, OPTION_PRIOR, dt, eps_abs, eps_rel, freeze_forcing, nb_obs);

    if (M == 0) {
        //simply return the sum of square or the log likelihood (can be used to do slices especially with least square where smc can't be used'...)

        int k;
        p_simplex->p_best->n_to_be_estimated = p_simplex->p_best->length;
        for(k=0; k< p_simplex->p_best->n_to_be_estimated; k++){
            p_simplex->p_best->to_be_estimated[k] = k;
        }

        FILE *p_file_trace = plom_fopen(SFR_PATH, GENERAL_ID, "trace", "w", header_trace, p_simplex->p_data);
        print_trace(p_file_trace, 0, p_simplex->p_best, p_simplex->p_data, f_simplex(p_simplex->p_best->mean, p_simplex));
        plom_fclose(p_file_trace);
    } else {
        //run the simplex algo
        simplex(p_simplex->p_best, p_simplex->p_data, p_simplex, f_simplex, CONVERGENCE_STOP_SIMPLEX, M, print_opt);
    }

    if (!(print_opt & PLOM_QUIET)) {
	time_end = s_clock();
	struct s_duration t_exec = time_exec(time_begin, time_end);
	snprintf(str, STR_BUFFSIZE, "Done in:= %dd %dh %dm %gs", t_exec.d, t_exec.h, t_exec.m, t_exec.s);
	print_log(str);
    }

    plom_print_done(theta, p_simplex->p_data, p_simplex->p_best, SFR_PATH, GENERAL_ID, print_opt);

    if (!(print_opt & PLOM_QUIET)) {
	print_log("clean up...");
    }

    json_decref(theta);

    clean_simplex(p_simplex);

    return 0;
}
