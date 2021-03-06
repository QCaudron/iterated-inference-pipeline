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

    char plom_help_string[] =
        "PLOM Sequential Monte Carlo\n"
        "usage:\n"
        "smc [implementation] [--no_dem_sto] [--no_white_noise] [--no_diff]\n"
        "                     [-r, --traj] [-p, --path <path>] [-i, --id <integer>] [-N, --n_thread <integer>]\n"
        "                     [-g, --freeze_forcing <float>]\n"
        "                     [-t, --no_filter] [-b, --no_trace] [-e, --no_hat]\n"
        "                     [-s, --DT <float>] [--eps_abs <float>] [--eps_rel <float>]\n"
        "                     [-l, --LIKE_MIN <float>] [-J <integer>] [--prior]\n"
	"                     [-q, --quiet] [-P, --pipe]"
        "                     [-h, --help]\n"
        "where implementation is 'ode', 'sde' or 'psr' (default)\n"
        "options:\n"
	"\n"
        "-q, --quiet        no verbosity\n"
        "-P, --pipe         pipe mode (echo theta.json on stdout)\n"
	"\n"
        "--no_dem_sto       turn off demographic stochasticity (if possible)\n"
        "--no_white_noise   turn off environmental stochasticity (if any)\n"
        "--no_diff          turn off drift (if any)\n"
	"\n"
        "-s, --DT           integration time step\n"
	"--eps_abs          Absolute error for adaptive step-size control\n"
	"--eps_rel          Relative error for adaptive step-size control\n"
        "-g, --freeze_forcing  freeze the metadata to their value at the specified time\n"
	"\n"
        "-i, --id           general id (unique integer identifier that will be appended to the output files)\n"
        "-p, --path         path where the outputs will be stored\n"
        "-N, --n_thread     number of threads to be used\n"
	"\n"
        "-l, --LIKE_MIN     particles with likelihood smaller that LIKE_MIN are considered lost\n"
        "-J                 number of particles\n"
	"-o, --nb_obs       number of observations to be fitted (for tempering)"
	"\n"
        "-r, --traj             print the trajectories\n"
        "-t, --no_filter    do not filter\n"
        "-b, --no_trace     do not write trace_<id>.output file\n"
        "-e, --no_hat       do not write hat_<id>.output file\n"
        "-d, --no_pred_res  do not write pred_res_<general_id>.output file (prediction residuals)\n"
	"\n"
        "--prior            add log(prior) to the estimated log likelihood\n"
	"\n"
        "-h, --help         print the usage on stdout\n";

    int filter = 1;
    double dt = 0.0, eps_abs = PLOM_EPS_ABS, eps_rel = PLOM_EPS_REL;

    enum plom_implementations implementation;
    enum plom_noises_off noises_off = 0;
    enum plom_print print_opt = PLOM_PRINT_BEST | PLOM_PRINT_HAT | PLOM_PRINT_PRED_RES;
   
    OPTION_PRIOR = 0;
    GENERAL_ID =0;
    snprintf(SFR_PATH, STR_BUFFSIZE, "%s", DEFAULT_PATH);
    J=1;
    LIKE_MIN = 1e-17;
    LOG_LIKE_MIN = log(1e-17);
    int nb_obs = -1;
    double freeze_forcing = -1.0;

    int n_threads = 1;

    while (1) {
        static struct option long_options[] =
            {
                /* These options don't set a flag We distinguish them by their indices (that are also the short option names). */
                {"traj",       no_argument,       0, 'r'},
                {"no_dem_sto", no_argument,       0, 'x'},
                {"no_white_noise", no_argument,       0, 'y'},
                {"no_diff",   no_argument,       0, 'z'},

                {"DT",         required_argument, 0, 's'},
                {"eps_abs",    required_argument, 0, 'v'},
                {"eps_rel",    required_argument, 0, 'w'},
		{"nb_obs", required_argument,  0, 'o'},

                {"freeze_forcing", required_argument, 0, 'g'},

                {"help",       no_argument,       0, 'h'},
                {"path",       required_argument, 0, 'p'},
                {"prior",  no_argument, &OPTION_PRIOR, 1},
                {"id",         required_argument, 0, 'i'},
                {"n_thread",   required_argument, 0, 'N'},
                {"no_filter",  no_argument,       0, 't'},
                {"no_trace",    no_argument,       0, 'b'},
                {"no_hat",     no_argument,       0, 'e'},
                {"no_pred_res",no_argument,       0, 'd'},

                {"quiet",  no_argument,       0, 'q'},
                {"pipe",  no_argument,       0, 'P'},

                {"LIKE_MIN",   required_argument, 0, 'l'},

                {0, 0, 0, 0}
            };
        /* getopt_long stores the option index here. */
        int option_index = 0;

        ch = getopt_long (argc, argv, "qPxyzs:v:w:p:i:J:l:trbehrN:o:g:", long_options, &option_index);

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
            noises_off |= PLOM_NO_DEM_STO;
            break;
        case 'y':
            noises_off |= PLOM_NO_ENV_STO;
            break;
        case 'z':
            noises_off |= PLOM_NO_DRIFT;
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
	    nb_obs = atoi(optarg);
            break;

        case 'g':
            freeze_forcing = atof(optarg);
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
        case 't':
            filter = 0;
            break;
        case 'J':
            J = atoi(optarg);
            break;
        case 'N':
            n_threads = atoi(optarg);
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
    struct s_data *p_data = build_data(settings, theta, implementation, noises_off, OPTION_PRIOR, nb_obs, "D");


    int size_proj = N_PAR_SV*N_CAC + p_data->p_it_only_drift->nbtot + N_TS_INC_UNIQUE;

    struct s_par *p_par = build_par(p_data);
    struct s_hat **D_p_hat = build_D_p_hat(p_data);
    struct s_X ***D_J_p_X = build_D_J_p_X(size_proj, N_TS, p_data, dt);
    struct s_X ***D_J_p_X_tmp = build_D_J_p_X(size_proj, N_TS, p_data, dt);
    struct s_best *p_best = build_best(p_data, theta);
    struct s_likelihood *p_like = build_likelihood();

    struct s_calc **calc = build_calc(&n_threads, GENERAL_ID, eps_abs, eps_rel, J, size_proj, step_ode, freeze_forcing, -1, p_data, settings);
    json_decref(settings);

    FILE *p_file_X = (print_opt & PLOM_PRINT_X) ? plom_fopen(SFR_PATH, GENERAL_ID, "X", "w", header_X, p_data): NULL;
    FILE *p_file_hat = (print_opt & PLOM_PRINT_HAT) ? plom_fopen(SFR_PATH, GENERAL_ID, "hat", "w", header_hat, p_data): NULL;
    FILE *p_file_pred_res = (print_opt & PLOM_PRINT_PRED_RES) ? plom_fopen(SFR_PATH, GENERAL_ID, "pred_res", "w", header_prediction_residuals, p_data): NULL;

    int64_t time_begin, time_end;
    if (!(print_opt & PLOM_QUIET)) {
	snprintf(str, STR_BUFFSIZE, "Starting plom-smc with the following options: i = %d, J = %d, LIKE_MIN = %g, n_threads = %d", GENERAL_ID, J, LIKE_MIN, n_threads);
	print_log(str);

	time_begin = s_clock();
    }

    back_transform_theta2par(p_par, p_best->mean, p_data->p_it_all, p_data);
    linearize_and_repeat(D_J_p_X[0][0], p_par, p_data, p_data->p_it_par_sv);
    prop2Xpop_size(D_J_p_X[0][0], p_data, calc[0]);
    theta_driftIC2Xdrift(D_J_p_X[0][0], p_best->mean, p_data);

    replicate_J_p_X_0(D_J_p_X[0], p_data);


#if FLAG_OMP
    run_SMC(D_J_p_X, D_J_p_X_tmp, p_par, D_p_hat, p_like, p_data, calc, get_f_pred(implementation, noises_off), filter, p_file_X, p_file_hat, p_file_pred_res, print_opt);
#else
    void *context = zmq_ctx_new ();

    void *sender = zmq_socket (context, ZMQ_PUSH);
    zmq_bind (sender, "inproc://server_sender");

    void *receiver = zmq_socket (context, ZMQ_PULL);
    zmq_bind (receiver, "inproc://server_receiver");

    void *controller = zmq_socket (context, ZMQ_PUB);
    zmq_bind (controller, "inproc://server_controller");

    struct s_thread_smc *p_thread_smc = malloc(n_threads*sizeof(struct s_thread_smc));
    pthread_t *worker = malloc(n_threads*sizeof(pthread_t));
    int nt, id;
    int J_chunk = J/n_threads;
	
    for (nt = 0; nt < n_threads; nt++) {
	p_thread_smc[nt].thread_id = nt;       	    	    
	p_thread_smc[nt].J_chunk = J_chunk;
	p_thread_smc[nt].J = J;
	p_thread_smc[nt].p_data = p_data;
	p_thread_smc[nt].p_par = p_par;
	p_thread_smc[nt].D_J_p_X = D_J_p_X;
	p_thread_smc[nt].p_calc = calc[nt];	
	p_thread_smc[nt].p_like = p_like;
	p_thread_smc[nt].context = context;
	pthread_create (&worker[nt], NULL, worker_routine_smc_inproc, (void*) &p_thread_smc[nt]);	
    }

    //wait that all worker are connected
    for (nt = 0; nt < n_threads; nt++) {
	zmq_recv(receiver, &id, sizeof (int), 0);
    }

    run_SMC_zmq_inproc(D_J_p_X, D_J_p_X_tmp, p_par, D_p_hat, p_like, p_data, calc, get_f_pred(implementation, noises_off), filter, p_file_X, p_file_hat, p_file_pred_res, print_opt, sender, receiver, controller);

    zmq_send (controller, "KILL", 5, 0);        
    zmq_close (sender);
    zmq_close (receiver);
    zmq_close (controller);

    for(nt = 0; nt < n_threads; nt++){
	pthread_join(worker[nt], NULL);
    }

    free(worker);
    free(p_thread_smc);

    zmq_ctx_destroy (context);
#endif

    if (OPTION_PRIOR) {
	double log_prob_prior_value;
	plom_err_code rc = log_prob_prior(&log_prob_prior_value, p_best, p_best->mean, p_best->var, p_data);
	if(rc != PLOM_SUCCESS){
#if FLAG_WARNING
	    print_err("error log_prob_prior computation");
#endif
	}

        p_like->Llike_best += log_prob_prior_value;
    }

    if (!(print_opt & PLOM_QUIET)) {
	time_end = s_clock();
	struct s_duration t_exec = time_exec(time_begin, time_end);
	if(filter){
	    snprintf(str, STR_BUFFSIZE, "logV: %g, J: %d,  n_all_fail: %d", p_like->Llike_best, (int) J, p_like->n_all_fail);
	}
	print_log(str);	
	snprintf(str, STR_BUFFSIZE, "Done in:= %dd %dh %dm %gs", t_exec.d, t_exec.h, t_exec.m, t_exec.s);
	print_log(str);
    }

    if (print_opt & PLOM_PRINT_X) {
        plom_fclose(p_file_X);
    }

    if (print_opt & PLOM_PRINT_HAT) {
        plom_fclose(p_file_hat);
    }

    if (print_opt & PLOM_PRINT_PRED_RES) {
        plom_fclose(p_file_pred_res);
    }

    if (print_opt & PLOM_PRINT_BEST) {
        FILE *p_file_trace = plom_fopen(SFR_PATH, GENERAL_ID, "trace", "w", header_trace, p_data);
        print_trace(p_file_trace, 0, p_best, p_data, p_like->Llike_best);
        plom_fclose(p_file_trace);
    }

    plom_print_done(theta, p_data, p_best, SFR_PATH, GENERAL_ID, print_opt);   

    if (!(print_opt & PLOM_QUIET)) {
	print_log("clean up...");
    }

    json_decref(theta);

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
