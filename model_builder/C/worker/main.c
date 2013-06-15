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
    double dt;
    double eps_abs;
    double eps_rel;
    char *host;
    void *context;
};

void *worker_routine (void *params) {
    char str[STR_BUFFSIZE];
    int j, jrcv, n, nn, nnp1, t1;

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
    struct s_calc *p_calc = build_p_calc(p->n_threads, p->thread_id, GENERAL_ID, p->eps_abs, p->eps_rel, size_proj, step_ode, p_data);

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
	    zmq_recv(server_receiver, &nn, sizeof (int), 0);
            nnp1 = nn+1;
            t1 = p_data->times[n];

            p_calc->current_n = n;
            p_calc->current_nn = nn;

	    recv_par(p_par, p_data, server_receiver);

            for(j=0; j<J; j++) {
		zmq_recv(server_receiver, &jrcv, sizeof (int), 0);
                //printf("j: %d jrcv %d\n", j, jrcv);
                recv_X(p_X, p_data, server_receiver);

                //do the computations..
                reset_inc(p_X, p_data);
				f_pred(p_X, nn, nnp1, p_par, p_data, p_calc);
                proj2obs(p_X, p_data);

                if(nnp1 == t1) {
                    like = exp(get_log_likelihood(p_X, p_par, p_data, p_calc));
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

	    snprintf(str, STR_BUFFSIZE, "worker %d: controller sent: %s", p->thread_id, buf);
	    print_log(str);

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

    snprintf(str, STR_BUFFSIZE, "thread %d done", p->thread_id);
    print_log(str);

    return NULL;
}

int main(int argc, char *argv[])
{
    char ch;
    char str[STR_BUFFSIZE];
    char host[STR_BUFFSIZE] = "127.0.0.1";
    int nt;

    /* set default values for the options */

    char sfr_help_string[] =
        "PLOM Worker\n"
        "usage:\n"
        "worker [implementation] [--no_dem_sto] [--no_white_noise] [--no_diff]\n"
        "                        [-s, --DT <float>] [--eps_abs <float>] [--eps_rel <float>]\n"
        "                        [-i, --id <integer>] [-h, --host <hostname>] [-P, --N_THREAD <integer>]\n"
        "                        [-l, --LIKE_MIN <float>] [-J <integer>]\n"
        "                        [--help]\n"
        "where implementation is 'ode', 'sde' or 'psr' (default)\n"
        "options:\n"
	"\n"
        "--no_dem_sto       turn off demographic stochasticity (if possible)\n"
        "--no_white_noise   turn off environmental stochasticity (if any)\n"
        "--no_diff          turn off drift (if any)\n"
	"\n"
        "-s, --DT           Initial integration time step\n"
	"--eps_abs          Absolute error for adaptive step-size contro\n"
	"--eps_rel          Relative error for adaptive step-size contro\n"
	"\n"
        "-i, --id           general id (unique integer identifier that will be appended to the output files)\n"
        "-h, --host         domain name or IP address of the particule server (defaults to 127.0.0.1)\n"
        "-P, --N_THREAD     number of threads to be used (defaults to the number of cores)\n"
        "-l, --LIKE_MIN     particles with likelihood smaller that LIKE_MIN are considered lost\n"
        "-c  -Jchunk        size of the chunk of particles\n"
	"-o, --nb_obs       number of observations to be fitted (for tempering)"
        "--help             print the usage on stdout\n";


    double dt = 0.0, eps_abs = PLOM_EPS_ABS, eps_rel = PLOM_EPS_REL;
    GENERAL_ID =0;
    J = 1; //here J is actualy Jchunk!
    
#if FLAG_OMP
    int n_threads = omp_get_max_threads();       
#else
    int n_threads = 1;
#endif

    LIKE_MIN = 1e-17;
    LOG_LIKE_MIN = log(1e-17);
    OPTION_TRAJ = 0;
    int nb_obs = -1;

    enum plom_implementations implementation;
    enum plom_noises_off noises_off = 0;

    while (1) {
        static struct option long_options[] =
            {
                {"help",       no_argument,       0, 'e'},

                {"no_dem_sto", no_argument,       0, 'x'},
                {"no_white_noise", no_argument,       0, 'y'},
                {"no_diff",   no_argument,       0, 'z'},

		{"DT",         required_argument, 0, 's'},
		{"eps_abs",    required_argument, 0, 'v'},
		{"eps_rel",    required_argument, 0, 'w'},

                {"id",         required_argument, 0, 'i'},
                {"N_THREAD",   required_argument, 0, 'P'},
                {"host",       required_argument, 0, 'h'},
                {"Jchunk",     required_argument, 0, 'c'},
		{"nb_obs", required_argument,  0, 'o'},

                {"LIKE_MIN",   required_argument, 0, 'l'},

                {0, 0, 0, 0}
            };
        /* getopt_long stores the option index here. */
        int option_index = 0;

        ch = getopt_long (argc, argv, "xyzs:v:w:i:P:c:h:l:o:", long_options, &option_index);

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
            print_log(sfr_help_string);
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

	case 'o':
	    nb_obs = atoi(optarg);
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
            print_log(sfr_help_string);
            return 1;
        }
    }

    json_t *settings = load_settings(PATH_SETTINGS);
    json_t *theta = load_json();

    n_threads = sanitize_n_threads(n_threads, J);
    n_threads = 1;

#if FLAG_VERBOSE
    snprintf(str, STR_BUFFSIZE, "Starting Plom-worker with the following options: i = %d, LIKE_MIN = %g, N_THREADS = %d", GENERAL_ID, LIKE_MIN, n_threads);
    print_log(str);
#endif
    struct s_data *p_data = build_data(settings, theta, implementation, noises_off, 1, nb_obs);
    json_decref(settings);
    json_decref(theta);

#if FLAG_VERBOSE
    print_log("setting up zmq context...");
#endif

    void *context = zmq_ctx_new();

#if FLAG_VERBOSE
    print_log("starting the threads...");
#endif

    struct s_thread_params *p_thread_params = malloc(n_threads*sizeof(struct s_thread_params));
    pthread_t *worker = malloc(n_threads*sizeof(pthread_t));
    for (nt = 0; nt < n_threads; nt++) {
        p_thread_params[nt].thread_id = nt;
        p_thread_params[nt].n_threads = n_threads;
        p_thread_params[nt].dt = dt;
        p_thread_params[nt].eps_abs = eps_abs;
        p_thread_params[nt].eps_rel = eps_rel;
        p_thread_params[nt].host = host;
        p_thread_params[nt].p_data = p_data;
        p_thread_params[nt].context = context;
        pthread_create (&worker[nt], NULL, worker_routine, (void*) &p_thread_params[nt]);
	snprintf(str, STR_BUFFSIZE, "worker %d started", nt);
	print_log(str);        
    }

    for(nt = 0; nt < n_threads; nt++){
        pthread_join(worker[nt], NULL);
    }

#if FLAG_VERBOSE
    print_log("clean up...");
#endif
    free(worker);
    free(p_thread_params);

    clean_data(p_data);

#if FLAG_VERBOSE
    print_log("closing zmq sockets...");
#endif

    zmq_ctx_destroy(context);

    return 0;
}
