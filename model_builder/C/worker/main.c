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

struct s_thread_params
{
    int thread_id;
    int n_threads;
    struct s_data *p_data;
    struct s_calc *p_calc;
    double dt;
    char *host;
    void *context;
};

void *worker_routine (void *params) {
    char str[STR_BUFFSIZE];
    int j, jrcv, n, np1, t0, t1;

    struct s_thread_params *p = (struct s_thread_params *) params;

    // Socket to server controller
    void *server_controller = zmq_socket (p->context, ZMQ_SUB);
    snprintf(str, STR_BUFFSIZE, "tcp://%s:%d", p->host, 5559);
    zmq_connect (server_controller, str);
    zmq_setsockopt (server_controller, ZMQ_SUBSCRIBE, "", 0);

    //  Socket to receive messages (particles) from the server
    void *server_receiver = zmq_socket (p->context, ZMQ_PULL);
    snprintf(str, STR_BUFFSIZE, "tcp://%s:%d", p->host, 5557);
    zmq_connect (server_receiver, str);

    //  Socket to send messages (results) to the server
    void *server_sender = zmq_socket (p->context, ZMQ_PUSH);
    snprintf(str, STR_BUFFSIZE, "tcp://%s:%d", p->host, 5558);
    zmq_connect (server_sender, str);
   
    struct s_data *p_data = p->p_data;
    struct s_par *p_par = build_par(p_data);
    int size_proj = N_PAR_SV*N_CAC + p_data->p_it_only_drift->nbtot + N_TS_INC_UNIQUE;
    struct s_X *p_X = build_X(size_proj, N_TS, p_data, p->dt);
    struct s_calc *p_calc = p->p_calc;

    double like = 0.0;

    plom_f_pred_t f_pred = get_f_pred(p_data->implementation, p_data->noises_off);

    zmq_pollitem_t items [] = {
        { server_receiver, 0, ZMQ_POLLIN, 0 },
        { server_controller, 0, ZMQ_POLLIN, 0 }
    };


    while (1) {
        zmq_poll (items, 2, -1);
        if (items [0].revents & ZMQ_POLLIN) {
	    
            //get a particle from the server
	    zmq_recv(server_receiver, &n, sizeof (int), 0);
            np1 = n+1;
            t0 = p_data->times[n];
            t1 = p_data->times[np1];

	    recv_par(p_par, p_data, server_receiver);

            for(j=0; j<J; j++) {
		zmq_recv(server_receiver, &jrcv, sizeof (int), 0);
                //printf("j: %d jrcv %d\n", j, jrcv);
                recv_X(p_X, p_data, server_receiver);

                //do the computations..
                reset_inc(p_X, p_data);
		f_pred(p_X, t0, t1, p_par, p_data, p_calc);

		proj2obs(p_X, p_data);

		if(p_data->data_ind[n]->n_nonan) {
		    like = exp(get_log_likelihood(p_X, p_par, p_data, p_calc, n, t1));
		}

                //send results
		zmq_send(server_sender, &jrcv, sizeof (int), ZMQ_SNDMORE);    
		send_X(server_sender, p_X, p_data, ZMQ_SNDMORE);
		zmq_send(server_sender, &like, sizeof (double), 0);
                //printf("j: %d jrcv %d sent back\n", j, jrcv);
            }
        }

        //controller commands:
        if (items [1].revents & ZMQ_POLLIN) {
	    char buf [256];
	    zmq_recv(server_controller, buf, 256, 0);           

            if(strcmp(buf, "KILL") == 0) {
                break;  //  Exit loop
            }
        }
    }

    zmq_close (server_receiver);
    zmq_close (server_sender);
    zmq_close (server_controller);

    clean_par(p_par);
    clean_X(p_X);
    clean_p_calc(p_calc, p_data);

    return NULL;
}

int main(int argc, char *argv[])
{
    char ch;
    char str[STR_BUFFSIZE];
    char host[STR_BUFFSIZE] = "127.0.0.1";
    int nt;

    /* set default values for the options */

    char plom_help_string[] =
        "PLOM Worker\n"
        "usage:\n"
        "worker [implementation] [--no_dem_sto] [--no_white_noise] [--no_diff]\n"
        "                        [-s, --DT <float>] [--eps_abs <float>] [--eps_rel <float>]\n"
        "                        [-i, --id <integer>] [-h, --host <hostname>] [-N, --n_thread <integer>]\n"
	"                        [-g, --freeze_forcing <float>]\n"
        "                        [-l, --LIKE_MIN <float>] [-J <integer>]\n"
	"                        [-q, --quiet]"
        "                        [-h, --help]\n"
        "where implementation is 'ode', 'sde' or 'psr' (default)\n"
        "options:\n"
	"\n"
        "-q, --quiet        no verbosity\n"
	"\n"
        "--no_dem_sto       turn off demographic stochasticity (if possible)\n"
        "--no_white_noise   turn off environmental stochasticity (if any)\n"
        "--no_diff          turn off drift (if any)\n"
	"\n"
        "-s, --DT           Initial integration time step\n"
	"--eps_abs          Absolute error for adaptive step-size control\n"
	"--eps_rel          Relative error for adaptive step-size control\n"
        "-g, --freeze_forcing freeze the metadata to their value at the specified time\n"
	"\n"
        "-i, --id           general id (unique integer identifier that will be appended to the output files)\n"
        "-h, --host         domain name or IP address of the particule server (defaults to 127.0.0.1)\n"
        "-N, --n_thread     number of threads to be used (defaults to the number of cores)\n"
        "-l, --LIKE_MIN     particles with likelihood smaller that LIKE_MIN are considered lost\n"
        "-c  -Jchunk        size of the chunk of particles\n"
	"-o, --nb_obs       number of observations to be fitted (for tempering)"
        "--help             print the usage on stdout\n";


    double dt = 0.0, eps_abs = PLOM_EPS_ABS, eps_rel = PLOM_EPS_REL;
    GENERAL_ID =0;
    J = 1; //here J is actualy Jchunk!

    enum plom_print print_opt = 0;
    
#if FLAG_OMP
    int n_threads = omp_get_max_threads();       
#else
    int n_threads = 1;
#endif

    LIKE_MIN = 1e-17;
    LOG_LIKE_MIN = log(1e-17);
    int nb_obs = -1;
    double freeze_forcing = -1.0;

    enum plom_implementations implementation;
    enum plom_noises_off noises_off = 0;

    while (1) {
        static struct option long_options[] =
            {
                {"help",       no_argument,       0, 'h'},

                {"no_dem_sto", no_argument,       0, 'x'},
                {"no_white_noise", no_argument,       0, 'y'},
                {"no_diff",   no_argument,       0, 'z'},

		{"DT",         required_argument, 0, 's'},
		{"eps_abs",    required_argument, 0, 'v'},
		{"eps_rel",    required_argument, 0, 'w'},

                {"freeze_forcing", required_argument, 0, 'g'},

                {"id",         required_argument, 0, 'i'},
                {"N_THREAD",   required_argument, 0, 'P'},
                {"host",       required_argument, 0, 'h'},
                {"Jchunk",     required_argument, 0, 'c'},
		{"nb_obs", required_argument,  0, 'o'},

                {"LIKE_MIN",   required_argument, 0, 'l'},

                {"quiet",  no_argument,       0, 'q'},

                {0, 0, 0, 0}
            };
        /* getopt_long stores the option index here. */
        int option_index = 0;

        ch = getopt_long (argc, argv, "hqxyzs:v:w:i:N:c:h:l:o:g:", long_options, &option_index);

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
            print_log(plom_help_string);
            return 1;

        case 'i':
            GENERAL_ID = atoi(optarg);
            break;

        case 'c':
            J = atoi(optarg);
            break;

        case 'P':
            n_threads = atoi(optarg);
            break;

        case 'h':
            snprintf(host, STR_BUFFSIZE, "%s", optarg);
            break;

        case 'l':
            LIKE_MIN = atof(optarg);
            LOG_LIKE_MIN = log(LIKE_MIN);
            break;


        case 'g':
            freeze_forcing = atof(optarg);
            break;

	case 'o':
	    nb_obs = atoi(optarg);
            break;

        case 'q':
	    print_opt |= PLOM_QUIET;
            break;

        case '?':
            /* getopt_long already printed an error message. */
	    exit(0);
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

    json_t *settings = load_settings(PATH_SETTINGS);
    json_t *theta = load_json();

    n_threads = sanitize_n_threads(n_threads, J);
    n_threads = 1;

    if (!(print_opt & PLOM_QUIET)) {
	snprintf(str, STR_BUFFSIZE, "Starting plom-worker with the following options: i = %d, LIKE_MIN = %g, N_THREADS = %d", GENERAL_ID, LIKE_MIN, n_threads);
	print_log(str);
    }

    struct s_data *p_data = build_data(settings, theta, implementation, noises_off, 1, nb_obs, "D");
    json_decref(theta);

    if (!(print_opt & PLOM_QUIET)) {
	print_log("setting up zmq context...");
    }

    void *context = zmq_ctx_new();

    if (!(print_opt & PLOM_QUIET)) {
	print_log("starting the threads...");
    }

    int size_proj = N_PAR_SV*N_CAC + p_data->p_it_only_drift->nbtot + N_TS_INC_UNIQUE;

    struct s_thread_params *p_thread_params = malloc(n_threads*sizeof(struct s_thread_params));
    pthread_t *worker = malloc(n_threads*sizeof(pthread_t));
    for (nt = 0; nt < n_threads; nt++) {
        p_thread_params[nt].thread_id = nt;
        p_thread_params[nt].n_threads = n_threads;
        p_thread_params[nt].dt = dt;
	p_thread_params[nt].p_calc = build_p_calc(n_threads, nt, GENERAL_ID, eps_abs, eps_rel, size_proj, step_ode, freeze_forcing, -1, p_data, settings);
        p_thread_params[nt].host = host;
        p_thread_params[nt].p_data = p_data;
        p_thread_params[nt].context = context;
        pthread_create (&worker[nt], NULL, worker_routine, (void*) &p_thread_params[nt]);
    }

    json_decref(settings);

    for(nt = 0; nt < n_threads; nt++){
        pthread_join(worker[nt], NULL);
    }

    if (!(print_opt & PLOM_QUIET)) {
	print_log("clean up...");
    }

    free(worker);
    free(p_thread_params);

    clean_data(p_data);

    if (!(print_opt & PLOM_QUIET)) {
	print_log("closing zmq sockets...");
    }

    zmq_ctx_destroy(context);

    return 0;
}
