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
#include <pthread.h>

struct s_thread_params
{
    int thread_id;
    int n_threads;
    struct s_data *p_data;
    enum plom_implementations implementation;
    enum plom_noises_off noises_off;
    char *IPv4;
    void *context;
};

void *worker_routine (void *params) {
    char str[STR_BUFFSIZE];
    int j, jrcv, rc, n, nn, nnp1, t1;

    struct s_thread_params *p = (struct s_thread_params *) params;

    // Socket to server controller
    void *server_controller = zmq_socket (p->context, ZMQ_SUB);
    snprintf(str, STR_BUFFSIZE, "tcp://%s:%d", p->IPv4, 5559);
    zmq_connect (server_controller, str);
    zmq_setsockopt (server_controller, ZMQ_SUBSCRIBE, "", 0);

    //  Socket to receive messages (particles) from the server
    void *server_receiver = zmq_socket (p->context, ZMQ_PULL);
    snprintf(str, STR_BUFFSIZE, "tcp://%s:%d", p->IPv4, 5557);
    zmq_connect (server_receiver, str);

    //  Socket to send messages (results) to the server
    void *server_sender = zmq_socket (p->context, ZMQ_PUSH);
    snprintf(str, STR_BUFFSIZE, "tcp://%s:%d", p->IPv4, 5558);
    zmq_connect (server_sender, str);



    struct s_data *p_data = p->p_data;
    struct s_par *p_par = build_par(p_data);
    int size_proj = N_PAR_SV*N_CAC + p_data->p_it_only_drift->nbtot + N_TS_INC_UNIQUE;
    struct s_X *p_X = build_X(size_proj, N_TS, p_data);
    struct s_calc *p_calc = build_p_calc(p->n_threads, p->thread_id, GENERAL_ID, p->implementation, size_proj, func, p_data);

    double like = 0.0;

    plom_f_pred_t f_pred = get_f_pred(p->implementation, p->noises_off);

    zmq_pollitem_t items [] = {
        { server_receiver, 0, ZMQ_POLLIN, 0 },
        { server_controller, 0, ZMQ_POLLIN, 0 }
    };

    while (1) {
        zmq_poll (items, 2, -1);
        if (items [0].revents & ZMQ_POLLIN) {

            //get a particle from the server
            n = recv_int(server_receiver);
            nn = recv_int(server_receiver);
            nnp1 = nn+1;
            t1 = p_data->times[n];

            p_calc->current_n = n;
            p_calc->current_nn = nn;

            recv_par(p_par, p_data, server_receiver);

            for(j=0; j<J; j++) {
                jrcv = recv_int(server_receiver);
                //printf("j: %d jrcv %d from %d n:%d nn%d\n", j, jrcv, p->thread_id, p_calc->current_n, p_calc->current_nn);
                recv_X(p_X, p_data, server_receiver);

                //do the computations..
                reset_inc(p_X, p_data);
		f_pred(p_X, nn, nnp1, p_par, p_data, p_calc);
                proj2obs(p_X, p_data);

                if(nnp1 == t1) {
                    like = exp(get_log_likelihood(p_X, p_par, p_data, p_calc));
                }

                //send results
                rc = send_int(server_sender, jrcv, ZMQ_SNDMORE);
                rc = send_X(server_sender, p_X, p_data, ZMQ_SNDMORE);
                rc = send_double(server_sender, like, 0);
            }
        }

        //controller commands:
        if (items [1].revents & ZMQ_POLLIN) {
            char *c_str = recv_str(server_controller);

            if(strcmp(c_str, "KILL") == 0) {
                printf("worker %d: controller sent: %s\n", p->thread_id, c_str);
                free(c_str);
                break;  //  Exit loop
            }
        }
    }

    zmq_close (server_receiver);
    zmq_close (server_sender);
    zmq_close (server_controller);

    clean_par(p_par);
    clean_X(p_X);
    clean_p_calc(p_calc, p->implementation);

    printf("thread %d done\n", p->thread_id);

    return NULL;
}

int main(int argc, char *argv[])
{
    char ch;
    char str[STR_BUFFSIZE];
    char IPv4[STR_BUFFSIZE] = "127.0.0.1";
    int nt;

    /* set default values for the options */

    char sfr_help_string[] =
        "PloM Worker\n"
        "usage:\n"
        "worker [implementation] [--no_dem_sto] [--no_env_sto] [--no_drift]\n"
        "                        [-i, --id <integer>] [-I, --IPv4 <ip address or DNS>] [-P, --N_THREAD <integer>]\n"
        "                        [-s, --DT <float>] [-l, --LIKE_MIN <float>] [-J <integer>]\n"
        "                        [--help]\n"
        "where implementation is 'ode', 'sde' or 'psr' (default)\n"
        "options:\n"
        "--no_dem_sto       turn off demographic stochasticity (if possible)\n"
        "--no_env_sto       turn off environmental stochasticity (if any)\n"
        "--no_drift         turn off drift (if any)\n"
        "-i, --id           general id (unique integer identifier that will be appended to the output files)\n"
        "-I, --IPv4         ip address or DNS of the particle server\n"
        "-P, --N_THREAD     number of threads to be used (defaults to the number of cores)\n"
        "-s, --DT           integration time step\n"
        "-l, --LIKE_MIN     particles with likelihood smaller that LIKE_MIN are considered lost\n"
        "-J  -Jchunk        size of the chunk of particles\n"
        "--help             print the usage on stdout\n";


    int has_dt_be_specified = 0;
    double dt_option;
    GENERAL_ID =0;
    J = 1; //here J is actualy Jchunk!
    int n_threads = omp_get_max_threads();
    LIKE_MIN = 1e-17;
    LOG_LIKE_MIN = log(1e-17);
    OPTION_TRAJ = 0;


    enum plom_implementations implementation;
    enum plom_noises_off noises_off = 0;


    while (1) {
        static struct option long_options[] =
            {
                {"help",       no_argument,       0, 'e'},
                {"no_dem_sto", no_argument,       0, 'x'},
                {"no_env_sto", no_argument,       0, 'y'},
                {"no_drift",   no_argument,       0, 'z'},
                {"id",         required_argument, 0, 'i'},
                {"N_THREAD",   required_argument, 0, 'P'},
                {"IPv4",       required_argument, 0, 'I'},
                {"Jchunk",     required_argument, 0, 'J'},
                {"DT",         required_argument, 0, 's'},
                {"LIKE_MIN",   required_argument, 0, 'l'},

                {0, 0, 0, 0}
            };
        /* getopt_long stores the option index here. */
        int option_index = 0;

        ch = getopt_long (argc, argv, "i:P:J:s:I:l:", long_options, &option_index);

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

        case 'i':
            GENERAL_ID = atoi(optarg);
            break;

        case 'J':
            J = atoi(optarg);
            break;

        case 'P':
            n_threads = atoi(optarg);
            break;

        case 'I':
            snprintf(IPv4, STR_BUFFSIZE, "%s", optarg);
            break;

        case 'l':
            LIKE_MIN = atof(optarg);
            LOG_LIKE_MIN = log(LIKE_MIN);
            break;

        case 's':
            dt_option = atof(optarg);
            has_dt_be_specified =1;
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

    if (has_dt_be_specified) {
        DT = dt_option;
    }
    DELTA_STO = round(1.0/DT);
    DT = 1.0/ ((double) DELTA_STO);

    n_threads = sanitize_n_threads(n_threads, J);

#if FLAG_VERBOSE
    snprintf(str, STR_BUFFSIZE, "Starting Plom-worker with the following options: i = %d, LIKE_MIN = %g, DT = %g, DELTA_STO = %g N_THREADS = %d", GENERAL_ID, LIKE_MIN, DT, DELTA_STO, n_threads);
    print_log(str);
#endif

    struct s_data *p_data = build_data(settings, theta, 1);
    json_decref(settings);
    json_decref(theta);

#if FLAG_VERBOSE
    print_log("setting up zmq context...");
#endif

    void *context = zmq_init (1);

#if FLAG_VERBOSE
    print_log("starting the threads...");
#endif

    struct s_thread_params *p_thread_params = malloc(n_threads*sizeof(struct s_thread_params));
    pthread_t *worker = malloc(n_threads*sizeof(pthread_t));
    for (nt = 0; nt < n_threads; nt++) {
        p_thread_params[nt].thread_id = nt;
        p_thread_params[nt].n_threads = n_threads;
        p_thread_params[nt].implementation = implementation;
        p_thread_params[nt].noises_off = noises_off;
        p_thread_params[nt].IPv4 = IPv4;
        p_thread_params[nt].p_data = p_data;
        p_thread_params[nt].context = context;
        pthread_create (&worker[nt], NULL, worker_routine, (void*) &p_thread_params[nt]);
        printf("thread %d started\n", nt);
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

    zmq_term (context);

    return 0;
}
