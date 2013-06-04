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

#include "simulation.h"

int main(int argc, char *argv[])
{
    char ch;
    char str[255];
    int i, j;


    /* set default values for the options */
    char sfr_help_string[] =
        "PLOM Simulation\n"
        "usage:\n"
        "simul [implementation] [--no_dem_sto] [--no_white_noise] [--no_diff]\n"
        "                       [-s, --DT <float>] [--eps_abs <float>] [--eps_rel <float>]\n"
        "                       [--traj] [--predict] [-p, --path <path>] [-i, --id <integer>] [-P, --N_THREAD <integer>]\n"
        "                       [-b, --bif] [--continue] [-l, --lyap]\n"
        "                       [-o, --t0 <integer>] [-D, --tend <integer>] [-T --transiant <integer>]\n"
        "                       [-B, --block <integer>] [-x, --precision <float>] [-J <integer>]\n"
        "                       [--help]\n"
        "where implementation is 'ode' (default), 'sde' or 'psr'\n"
        "options:\n"
        "\n"
        "--no_dem_sto       turn off demographic stochasticity (if possible)\n"
        "--no_white_noise       turn off environmental stochasticity (if any)\n"
        "--no_diff         turn off drift (if any)\n"
        "\n"
        "-s, --DT           Initial integration time step\n"
        "--eps_abs          Absolute error for adaptive step-size contro\n"
        "--eps_rel          Relative error for adaptive step-size contro\n"
        "\n"
        "--traj             print the trajectories\n"
        "--continue         print the final states in a bifurcation analysis to allow continuation\n"
        "--predict          load an array of theta for Bayesian prediction\n"
        "-p, --path         path where the outputs will be stored\n"
        "-i, --id           general id (unique integer identifier that will be appended to the output files)\n"
        "-P, --N_THREAD     number of threads to be used (default to the number of cores)\n"
        "-b, --bif          run a bifurcation analysis\n"
        "-l, --lyap         compute lyapunov exponents\n"
        "-d, --period_dyn   compute period (dynamical system def)\n"
        "-f, --fft          compute period (FFT)\n"
        "-o, --t0           time when the simulation starts\n"
        "-D, --tend         time when the simulation ends\n"
        "-T, --transiant    skip a transiant of the specified duration\n"
        "-B, --block        tuning parameter for max and min detection (has to be an odd number)\n"
        "-x, --precision    smallest significant difference to detect variation for min and max detections\n"
        "-J                 number of realisations\n"
        "--help             print the usage on stdout\n";


    enum plom_implementations implementation;
    enum plom_noises_off noises_off = 0;

    double t0 = 0.0, t_end = 0.0, t_transiant = 0.0;
    int nn0 = 0; //for PAR_FIXED: t can be > N_DATA_PAR_FIXED: For transiant and lyap, we use p_calc->current_nn = t0 if t0 < N_DATA_PAR_FIXED. For traj_obs, we let p_calc->current_nn vary starting from nn0 and up to N_DATA_PAR_FIXED. After N_DATA_PAR_FIXED, the last value is recycled

    double dt = 0.0, eps_abs = PLOM_EPS_ABS, eps_rel = PLOM_EPS_REL;

    OPTION_TRAJ = 0;
    int OPTION_LYAP = 0;
    int OPTION_BIF = 0;
    static int OPTION_CONTINUE = 0;
    static int OPTION_PREDICT = 0;
    int OPTION_PERIOD_DYNAMICAL_SYTEM = 0;
    int OPTION_FFT = 0;

    PRECISION = 1.0e-2;
    N_BLOC = 5;
    GENERAL_ID =0;
    snprintf(SFR_PATH, STR_BUFFSIZE, "%s", DEFAULT_PATH);
    J=1;

    int n_threads=omp_get_max_threads();

    while (1) {
        static struct option long_options[] =
            {
                /* These options set a flag. */
                {"traj", no_argument,       &OPTION_TRAJ, 1},
                {"continue", no_argument,   &OPTION_CONTINUE, 1},
                {"predict", no_argument,   &OPTION_PREDICT, 1},
                /* These options don't set a flag We distinguish them by their indices (that are also the short option names). */
                {"no_dem_sto", no_argument,       0, 'x'},
                {"no_white_noise", no_argument,       0, 'y'},
                {"no_diff",   no_argument,       0, 'z'},

                {"DT",         required_argument, 0, 's'},
                {"eps_abs",    required_argument, 0, 'v'},
                {"eps_rel",    required_argument, 0, 'w'},

                {"help", no_argument,  0, 'e'},
                {"path",    required_argument, 0, 'p'},
                {"id",    required_argument, 0, 'i'},
                {"N_THREAD",   required_argument, 0, 'P'},

                {"bif",    no_argument, 0, 'b'},
                {"lyap",    no_argument, 0, 'l'},
                {"period_dyn",    no_argument, 0, 'd'},
                {"fft",    no_argument, 0, 'f'},

                {"t0",    required_argument, 0, 'o'},
                {"tend",    required_argument, 0, 'D'},
                {"transiant",    required_argument, 0, 'T'},
                {"block",    required_argument, 0, 'B'},
                {"precision",    required_argument, 0, 'r'},

                {0, 0, 0, 0}
            };
        /* getopt_long stores the option index here. */
        int option_index = 0;

        ch = getopt_long (argc, argv, "xyzs:v:w:B:r:i:J:s:D:T:bldfp:o:P:", long_options, &option_index);

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

        case 'p':
            snprintf(SFR_PATH, STR_BUFFSIZE, "%s", optarg);
            break;
        case 'i':
            GENERAL_ID = atoi(optarg);
            break;
        case 'o':
            t0 = floor(atof(optarg));
            break;
        case 'D':
            t_end = ceil(atof(optarg));
            break;
        case 'r':
            PRECISION = atof(optarg);
            break;
        case 'T':
            t_transiant = atof(optarg);
            break;
        case 'b':
            OPTION_BIF = 1;
            break;
        case 'l':
            OPTION_LYAP = 1;
            break;
        case 'd':
            OPTION_PERIOD_DYNAMICAL_SYTEM = 1;
            break;
        case 'f':
            OPTION_FFT = 1;
            break;
        case 'B':
            N_BLOC = atoi(optarg);
            break;
        case 'J':
            J = atoi(optarg);
            break;
        case 'P':
            n_threads = atoi(optarg);
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
        implementation = PLOM_ODE;
        noises_off = noises_off | PLOM_NO_DEM_STO| PLOM_NO_ENV_STO | PLOM_NO_DRIFT;
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

        if(OPTION_LYAP && (implementation != PLOM_ODE)){
            print_err("Lyapunov exponents can only be computed with ODE implementation");
            exit(EXIT_FAILURE);
        }
    }

    if (t0>t_end) {
        snprintf(str, STR_BUFFSIZE,  "t0 = %g > t_end = %g, now quiting", t0, t_end);
        print_err(str);
        exit(EXIT_FAILURE);
    }

    json_t *settings = load_settings(PATH_SETTINGS);
    json_t *thetas = load_json();
    json_t *theta = (OPTION_PREDICT) ? json_array_get(thetas, 0): thetas;


    if((OPTION_BIF || OPTION_LYAP) && (J>1)) {
        J=1;
        print_log("for bifurcation analysis and Lyapunov exp. comupations, J must be 1 !!");
    }

    if( OPTION_PREDICT && (J != json_array_size(thetas)) ){
        J=json_array_size(thetas);
        print_log("for Bayesian prediction J is set to the length of the list of thetas");
    }


    struct s_data *p_data = build_data(settings, theta, implementation, noises_off, 0, -1);
    json_decref(settings);

    int size_proj = N_PAR_SV*N_CAC + p_data->p_it_only_drift->nbtot + N_TS_INC_UNIQUE;

    struct s_par **J_p_par = build_J_p_par(p_data);
    struct s_X **J_p_X = build_J_p_X(size_proj, N_TS, p_data, dt);
    struct s_best *p_best = build_best(p_data, theta);

    struct s_calc **calc = build_calc(&n_threads, GENERAL_ID, eps_abs, eps_rel, J, size_proj, step_ode, p_data);

    double *y0 = init1d_set0(N_PAR_SV*N_CAC + N_TS_INC_UNIQUE);
    double abs_tol = eps_abs, rel_tol = eps_rel;


    //if t_transiant > N_DATA: we ensure constant pop size by settings mu_d = mu_b in case of variable birth and death reates
    if (t_transiant > N_DATA) {
        ensure_cst_pop_size(p_data);
    }

    //setting p_calc->current_nn
    if ( t0 < N_DATA_PAR_FIXED ) {
        nn0 = t0;
    } else {
        nn0 = N_DATA_PAR_FIXED -1;
    }
    store_state_current_n_nn(calc, 0, nn0);


#if FLAG_VERBOSE
    snprintf(str, STR_BUFFSIZE, "Starting Plom with the following options: i = %d, t0 = %g, t_end = %g, t_transiant = %g, N_BLOC = %d, PRECISION = %g, N_THREADS = %d, J = %d", GENERAL_ID, t0, t_end, t_transiant, N_BLOC, PRECISION, n_threads, J);
    print_log(str);
#endif


    if (OPTION_FFT) {
#if FLAG_VERBOSE
        print_warning("FFT requested, to compute FFT we change the trajectory length to its nearest upper power of 2 number...");
#endif
        snprintf(str, STR_BUFFSIZE, "t_end-t0: %g -> %g, the number of stored value being %g", t_end-t0, nextpow2((t_end-t0)),  nextpow2((t_end-t0)));
        print_warning(str);

        t_end = t0 + nextpow2((t_end-t0));
    }

    if(OPTION_PREDICT){
        for(j=0; j<J; j++) {
            load_best(p_best, p_data, json_array_get(thetas, j), 1);
            transform_theta(p_best, p_data, 1);
            back_transform_theta2par(J_p_par[j], p_best->mean, p_data->p_it_all, p_data);
            linearize_and_repeat(J_p_X[j], J_p_par[j], p_data, p_data->p_it_par_sv);
            prop2Xpop_size(J_p_X[j], p_data);
            theta_driftIC2Xdrift(J_p_X[j], p_best->mean, p_data);
        }
    } else{
        transform_theta(p_best, p_data, 1);
        back_transform_theta2par(J_p_par[0], p_best->mean, p_data->p_it_all, p_data);
        linearize_and_repeat(J_p_X[0], J_p_par[0], p_data, p_data->p_it_par_sv);
        prop2Xpop_size(J_p_X[0], p_data);
        theta_driftIC2Xdrift(J_p_X[0], p_best->mean, p_data);

        replicate_J_p_X_0(J_p_X, p_data);
        for(j=0; j<J; j++) {
            back_transform_theta2par(J_p_par[j], p_best->mean, p_data->p_it_all, p_data);
        }
    }

    json_decref(thetas);

    for (i=0; i<(N_PAR_SV*N_CAC); i++){
        y0[i] = J_p_X[0]->proj[i];
    }

    plom_f_pred_t f_pred = get_f_pred(p_data->implementation, p_data->noises_off);

    /**************************************************/
    /* SKIP TRANSIANT (not possible with prediction)  */
    /**************************************************/
    if ( (!OPTION_PREDICT) && (t_transiant > 0.0) ) {
        t_transiant = floor(t_transiant/ONE_YEAR)*ONE_YEAR + t0;
#if FLAG_VERBOSE
        snprintf(str, STR_BUFFSIZE,  "skipping transiant... (Note that t_transiant has been ajusted to %g to respect seasonality and t0)", t_transiant );
        print_warning(str);
#endif

        if (implementation == PLOM_ODE) {
            //only the first particle is used to skip transiant
            if ( integrate(J_p_X[0], y0, t0, t_transiant, J_p_par[0], &abs_tol, &rel_tol, calc[0], p_data) ) {
                print_err("integration error, the program will now quit");
                exit(EXIT_FAILURE);
            }
        } else {
            int thread_id, ip1;
            for(i=t0 ; i<t_transiant ; i++) {
                ip1 = i+1;
#pragma omp parallel for private(thread_id)
                for(j=0; j<J; j++) {
                    thread_id = omp_get_thread_num();
                    reset_inc(J_p_X[j], p_data);
                    f_pred(J_p_X[j], i, ip1, J_p_par[0], p_data, calc[thread_id]);
                }
            }
        }
    }


    if(!(OPTION_BIF || OPTION_LYAP) && (t_end==0)) { //ONLY INITIAL CONDITIONS

        //we need to integrate for at least 1 time step so that
        //incidence is reset to 0 after transiant (transiant did not
        //reset incidences every week in case of PLOM_ODE)
        traj(J_p_X, t0, t0+1, t_transiant, J_p_par, p_data, calc, f_pred);

    } else if(!(OPTION_BIF || OPTION_LYAP) && OPTION_TRAJ && (t_end>t0)) { //ONLY TRAJ

        traj(J_p_X, t0, t_end, t_transiant, J_p_par, p_data, calc, f_pred);

    } else if ( (!OPTION_PREDICT) && (OPTION_BIF || OPTION_LYAP) ) {

        /*store (potentialy new, if t_transiant > 0.0) initial conditions*/
        for (i=0; i<(N_PAR_SV*N_CAC); i++){
            y0[i] = J_p_X[0]->proj[i];
        }
        reset_inc(J_p_X[0], p_data);

        /****************************************/
        /*************BIF************************/
        /****************************************/

        if (OPTION_BIF || OPTION_PERIOD_DYNAMICAL_SYTEM || OPTION_FFT) {
            int ts;

#if FLAG_VERBOSE
            print_log("bifurcation analysis");
#endif

            double **traj_obs = get_traj_obs(J_p_X[0], y0, t0, t_end, t_transiant, J_p_par[0], p_data, calc[0], f_pred); //[N_TS][(t_end-t0)]

            for (ts=0; ts< N_TS; ts++) {

                if (OPTION_PERIOD_DYNAMICAL_SYTEM){
                    period_dynamical_system(traj_obs[ts], (int) (t_end-t0), ts);
                }

                if (OPTION_BIF){
                    max_min(traj_obs[ts], J_p_par[0], p_data, calc[0], t0, (int) (t_end-t0), ts);
                }

                if (OPTION_FFT){
                    //fourrier has to be last as it modified traj_obs in place!!
                    fourrier_power_spectrum(traj_obs[ts], (int) (t_end-t0), ts);
                }
            }

            clean2d(traj_obs, N_TS);

            //print hat for continuation
            if(OPTION_CONTINUE) {
                struct s_hat *p_hat = build_hat(p_data);
                compute_hat_nn(J_p_X, J_p_par, p_data, calc, p_hat, 1);
                FILE *p_file_hat = sfr_fopen(SFR_PATH, GENERAL_ID, "hat", "w", header_hat, p_data);
                print_p_hat(p_file_hat, NULL, p_hat, p_data, 0);
                sfr_fclose(p_file_hat);
                clean_hat(p_hat, p_data);
            }
        }

        /****************************************/
        /*************LYAP***********************/
        /****************************************/

        if (OPTION_LYAP) {

#if FLAG_VERBOSE
            print_log("Lyapunov exponents computation...");
#endif
            store_state_current_n_nn(calc, 0, nn0);
            lyapunov(calc[0], J_p_par[0], y0, t0, t_end, abs_tol, rel_tol, J_p_X[0]->dt);
        }

    } //end OPTION_BIF or OPTION_LYAP


#if FLAG_VERBOSE
    print_log("clean up...");
#endif

    FREE(y0);

    clean_J_p_X(J_p_X);
    clean_best(p_best);
    clean_J_p_par(J_p_par);
    clean_calc(calc, p_data);
    clean_data(p_data);

    return 0;
}
