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

/* automatically generated code: order of the parameters */
{% for o in order.var %}
#define ORDER_{{ o|safe }} {{ forloop.counter0 }}{% endfor %}

{% for o in order.universe %}
#define ORDER_{{ o.name|safe }} {{ o.order }}{% endfor %}

{% for o in order.drift %}
#define ORDER_{{ o|safe }} {{ forloop.counter0 }}{% endfor %}

{% for o in order.data %}
#define ORDER_{{ o|safe }} {{ forloop.counter0 }}{% endfor %}


struct s_kalman_specific_data *build_kalman_specific_data(struct s_data *p_data, enum plom_noises_off noises_off)
{
    int i;

    struct s_kalman_specific_data *p;
    p = malloc(sizeof(struct s_kalman_specific_data));
    if(p==NULL){
        char str[STR_BUFFSIZE];
        snprintf(str, STR_BUFFSIZE, "Allocation impossible in file :%s line : %d",__FILE__,__LINE__);
        print_err(str);
        exit(EXIT_FAILURE);
    }

    //compo_groups_drift_par_proc
    p->compo_groups_drift_par_proc = malloc(N_DRIFT_PAR_PROC* sizeof (struct s_group **));
    if(p->compo_groups_drift_par_proc==NULL) {
        char str[STR_BUFFSIZE];
        snprintf(str, STR_BUFFSIZE, "Allocation impossible in file :%s line : %d",__FILE__,__LINE__);
        print_err(str);
        exit(EXIT_FAILURE);
    }
    for(i=0; i<N_DRIFT_PAR_PROC; i++) {
        p->compo_groups_drift_par_proc[i] = get_groups_compo(p_data->routers[ p_data->p_drift->ind_par_Xdrift_applied[i] ]);
    }

    p->FtCt = gsl_matrix_alloc(N_KAL, N_KAL);
    p->Ft = gsl_matrix_calloc(N_KAL, N_KAL);


    int n_noise;
    if(noises_off & PLOM_NO_DEM_STO) { //demographic stochasticity is **not** taken into account
        n_noise = {{ Q.deter.s }}*N_CAC;
    } else {
        n_noise = {{ Q.sto.s }}*N_CAC;
    }
    n_noise += p_data->p_it_only_drift->nbtot;

    if(n_noise == 0){
        print_err("kalman methods must be used with at least one brownian motion, try running with the sto command or add noise into your process model");
        exit(EXIT_FAILURE);
    }

    p->diag_Qc = init1d_set0(n_noise);
    p->L = gsl_matrix_calloc(N_KAL, n_noise);
    eval_L(p->L, p_data, noises_off);
    p->LQc = gsl_matrix_calloc(N_KAL, n_noise);

    return p;
}



/**
 * the function used by f_prediction_ode_rk:
 * dX/dt = f(t, X, params)
 *
 * it is an expanded copy of core/prediction.c/func
 */
int func_kal(double t, const double X[], double f[], void *params)
{
    struct s_calc *p_calc = (struct s_calc *) params;
    struct s_par *p_par = p_calc->p_par;  /* syntaxic shortcut */
    struct s_data *p_data = p_calc->p_data;
    struct s_obs2ts **obs2ts = p_data->obs2ts;		/* syntaxic shortcut */
    struct s_router **routers = p_data->routers;	/* syntaxic shortcut */

    struct s_kalman_specific_data *p_kalman_specific_data = (struct s_kalman_specific_data *) p_calc->method_specific_thread_safe_data;

    gsl_matrix *Ft = p_kalman_specific_data->Ft;
    gsl_matrix *Q = p_calc->Q;
    struct s_group ***compo_groups_drift_par_proc = p_kalman_specific_data->compo_groups_drift_par_proc;
    gsl_matrix *FtCt = p_kalman_specific_data->FtCt;
    gsl_matrix_const_view Ct   = gsl_matrix_const_view_array(&X[N_PAR_SV*N_CAC+N_TS_INC_UNIQUE],N_KAL,N_KAL);
    gsl_matrix_view res2 = gsl_matrix_view_array(&f[N_PAR_SV*N_CAC+N_TS_INC_UNIQUE],N_KAL,N_KAL);

    int c, ac, cac, n_cac, ts, o;
    double sum_inc = 0.0;
    int offset;

    const int nn = p_calc->current_nn;
    double **par = p_par->natural;
    double ***covar = p_data->par_fixed;

    {% if current_p %}
    for (c=0;c<N_C;c++) {
        for (ac=0;ac<N_AC;ac++) {
            p_par->current_p[c][ac]=get_current_pop_size(X,c,ac);
        }
    }
    {% endif %}

    double _r[N_CAC][{{print_ode.caches|length}}];
    {% if print_ode.sf %}
    double _sf[N_CAC][{{print_ode.sf|length}}];{% endif %}
    for(cac=0;cac<N_CAC;cac++) {
        {% for sf in print_ode.sf %}
        _sf[cac][{{ forloop.counter0 }}] = {{ sf|safe }};{% endfor %}

        {% for cache in print_ode.caches %}
        _r[cac][{{ forloop.counter0 }}] = {{ cache|safe }};{% endfor %}
    }


    for (c=0;c<N_C;c++) {
        for(ac=0; ac<N_AC; ac++) {
            cac = c*N_AC+ac;

            {{ print_ode.sys|safe }}
        }
    }

    /*automaticaly generated code:*/
    /*compute incidence:integral between t and t+1*/

    offset=0;
    {% for eq in print_ode.obs %}
    o = {{ eq.true_ind_obs|safe }};

    for (ts=0; ts<obs2ts[o]->n_ts_unique; ts++) {
        sum_inc = 0.0;
        for(n_cac=0; n_cac<obs2ts[o]->n_cac[ts]; n_cac++) {
            c = obs2ts[o]->cac[ts][n_cac][0];
            ac = obs2ts[o]->cac[ts][n_cac][1];
            cac = c*N_AC+ac;

            sum_inc += {{ eq.right_hand_side|safe }};
        }

        f[N_PAR_SV*N_CAC +offset] = sum_inc;
        offset++;
    }
    {% endfor %}

    ////////////////
    // covariance //
    ////////////////

    // evaluate Q and jacobian
    eval_Q(Q, X, p_par, p_data, p_calc, p_kalman_specific_data, t);
    eval_jac(Ft, X, p_par, p_data, p_calc, compo_groups_drift_par_proc, t);

    // compute Ft*Ct+Ct*Ft'+Q
    gsl_matrix_set_zero(&res2.matrix);
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, Ft, &Ct.matrix, 0.0, FtCt);
    gsl_matrix_add(&res2.matrix,FtCt);
    gsl_matrix_transpose (FtCt);
    gsl_matrix_add(&res2.matrix,FtCt);
    gsl_matrix_add(&res2.matrix,Q);

    return GSL_SUCCESS;
}



int cac_drift_in_cac_ts(int cac_drift, int o, int ts_unique, struct s_obs2ts **obs2ts)
{
    /* helper function for fourth part of eval_jac computation:
       return 1 if cac_drift is in cac of the considered time serie (o, ts_unique) */

    int n_cac_ts, c_ts, ac_ts, cac_ts;

    for(n_cac_ts=0; n_cac_ts< obs2ts[o]->n_cac[ts_unique]; n_cac_ts++) {
        c_ts = obs2ts[o]->cac[ts_unique][n_cac_ts][0];
        ac_ts = obs2ts[o]->cac[ts_unique][n_cac_ts][1];
        cac_ts = c_ts*N_AC+ac_ts;

        if(cac_ts == cac_drift) {
            return 1;
        }
    }

    return 0;
}

void eval_jac(gsl_matrix *Ft, const double *X, struct s_par *p_par, struct s_data *p_data, struct s_calc *p_calc, struct s_group ***compo_groups_drift_par_proc, double t)
{
    //X is p_X->proj

    int c, ac, cac;
    int ts, ts_unique, stream, n_cac;

    /*t is the current time in the unit of the data */

    //syntaxic shortcut
    struct s_obs2ts **obs2ts = p_data->obs2ts;
    struct s_router **routers = p_data->routers;  /* syntaxic shortcut */

    {% if is_drift %}
    int d, g;
    struct s_drift *p_drift =  p_data->p_drift;
    {% endif %}

    //the automaticaly generated code may need these variables
    const int nn = p_calc->current_nn;

    double **par = p_par->natural;
    double ***covar = p_data->par_fixed;

    //some terms are always 0: derivative of the ODE (excluding the observed variable) against the observed variable, derivative of the dynamic of the observed variable against the observed variables, derivative of the drift eq.
    gsl_matrix_set_zero(Ft);


    double _rj[N_CAC][{{ jacobian.caches|length }}];
    {% if jacobian.sf %}
    double _sf[N_CAC][{{ jacobian.sf|length }}];{% endif %}
    for(cac=0; cac<N_CAC; cac++){
        {% for sf in jacobian.sf %}
        _sf[cac][{{ forloop.counter0 }}] = {{ sf|safe }};{% endfor %}

        {% for cache in jacobian.caches %}
        _rj[cac][{{ forloop.counter0 }}] = {{ cache|safe }};{% endfor %}
    }

    //first non null part of the jacobian matrix: derivative of the ODE (excluding the observed variable) against the state variable only ( automaticaly generated code )
    for(c=0; c<N_C; c++) {
        for(ac=0; ac<N_AC; ac++) {
            cac = c*N_AC+ac;
            {% for jac_i in jacobian.jac %}
            {% for jac_ii in jac_i %}
            gsl_matrix_set(Ft, {{ forloop.parentloop.counter0 }}*N_CAC+cac, {{ forloop.counter0 }}*N_CAC+cac, _rj[cac][{{ jac_ii|safe }}]);
            {% endfor %}
            {% endfor %}
        }
    }

    //second non null part of the jacobian matrix: derivative of the dynamic of the observed variable against the state variable only ( automaticaly generated code )
    ts = 0;
    {% for jac_i in jacobian.jac_obs %}
    for(ts_unique=0; ts_unique < obs2ts[{{ forloop.counter0 }}]->n_ts_unique; ts_unique++) {
        for(stream=0; stream < obs2ts[{{ forloop.counter0 }}]->n_stream[ts_unique]; stream++) {

            for(n_cac=0; n_cac< obs2ts[{{ forloop.counter0 }}]->n_cac[ts_unique]; n_cac++) {
                c = obs2ts[{{ forloop.counter0 }}]->cac[ts_unique][n_cac][0];
                ac = obs2ts[{{ forloop.counter0 }}]->cac[ts_unique][n_cac][1];
                cac = c*N_AC+ac;

                {% for jac_ii in jac_i %}
                gsl_matrix_set(Ft, N_PAR_SV*N_CAC+ts, {{ forloop.counter0 }}*N_CAC+cac, _rj[cac][{{ jac_ii|safe }}]);
                {% endfor %}

            }

            ts++;
        }
    }
    {% endfor %}

    //third non null part of the jacobian matrix: derivative of the ODE (excluding the observed variable) against the drift variable (automaticaly generated code)
    //non null only for 'cac' present in group of Xdrift: compo_groups_drift_par_proc gives these 'cac'
    {% for jac_i in jacobian.jac_drift %}
    d = 0;
    {% for jac_ii in jac_i %}
    for(g=0; g< routers[ p_drift->ind_par_Xdrift_applied[{{ forloop.counter0 }}] ]->n_gp; g++) {
        for(n_cac=0; n_cac< compo_groups_drift_par_proc[{{ forloop.counter0 }}][g]->size; n_cac++) {
            cac = compo_groups_drift_par_proc[{{ forloop.counter0 }}][g]->elements[n_cac];
            get_c_ac(cac, &c, &ac);
            gsl_matrix_set(Ft,
                           {{ forloop.parentloop.counter0 }}*N_CAC+cac,
                           N_PAR_SV*N_CAC + N_TS + d,
                           drift_derivative(_rj[cac][{{ jac_ii.value|safe }}], {{ jac_ii.der|safe }}, routers[ORDER_{{ jac_ii.name|safe }}], cac));
        }
        d++;
    }
    {% endfor %}
    {% endfor %}


    //fourth non null part of the jacobian matrix: derivative of N_TS agains drift (automaticaly generated code)
    //this is not optimal at all and will be slow as hell: we should cache the intersection of cac contained in an observation and cac contained in a group of X_drift.
    //anyone tempted ?
    ts = 0;
    {% for jac_i in jacobian.jac_obs_drift %}
    for(ts_unique=0; ts_unique < obs2ts[{{ forloop.counter0 }}]->n_ts_unique; ts_unique++) {
        for(stream=0; stream < obs2ts[{{ forloop.counter0 }}]->n_stream[ts_unique]; stream++) {
            d = 0;
            {% for jac_ii in jac_i %}
            for(g=0; g< routers[ p_drift->ind_par_Xdrift_applied[{{ forloop.counter0 }}] ]->n_gp; g++) {
                double sum_tmp = 0.0;
                for(n_cac=0; n_cac< compo_groups_drift_par_proc[{{ forloop.counter0 }}][g]->size; n_cac++) {
                    cac = compo_groups_drift_par_proc[{{ forloop.counter0 }}][g]->elements[n_cac];
                    if(cac_drift_in_cac_ts(cac, {{ forloop.parentloop.counter0 }}, ts_unique, obs2ts)) {
                        get_c_ac(cac, &c, &ac);
                        sum_tmp += drift_derivative(_rj[cac][{{ jac_ii.value|safe }}], {{ jac_ii.der|safe }}, routers[ORDER_{{ jac_ii.name|safe }}], cac);
                    }
                }
                gsl_matrix_set(Ft, N_PAR_SV*N_CAC + ts, N_PAR_SV*N_CAC + N_TS + d, sum_tmp);
                d++;
            }
            {% endfor %}
            ts++;
        }
    }
    {% endfor %}

}

/**
 * derivative of the mean of the observation process against state
 * variables and observed variables
 */
void eval_ht(gsl_vector *ht, gsl_vector *xk, struct s_par *p_par, struct s_data *p_data, struct s_calc *p_calc, int ts)
{
    struct s_router **routers = p_data->routers;  /* syntaxic shortcut */

    //the automaticaly generated code may need these variables
    int n, nn;
    n = p_calc->current_n;
    nn = p_calc->current_nn;
    double t;
    t = (double) p_data->times[n];

    double **par = p_par->natural;
    double ***covar = p_data->par_fixed;

    gsl_vector_set_zero(ht);
    //derivative against state variable are always nul so we focus on the derivative against the observed variable

    double x = gsl_vector_get(xk, N_PAR_SV*N_CAC +ts); //the derivative are templated (automaticaly generated code) and are a function of "x". we link "x" to the right observed variable.
    gsl_vector_set(ht, N_PAR_SV*N_CAC +ts, {{ jac_proc_obs|safe }});
}



/**
 * evaluate dispersion matrix L (as described in Sarkka phD (2006))
 */
void eval_L(gsl_matrix *L, struct s_data *p_data, enum plom_noises_off noises_off)
{
    int c, ac, cac;
    int ts, ts_unique, stream, n_cac;
    struct s_obs2ts **obs2ts = p_data->obs2ts;

    {% for command, Ls in Q.items %}

    {% if command == 'sto' %}
    if(!(noises_off & PLOM_NO_DEM_STO)) { //demographic stochasticity
    {% else %}
    } else { //no demographic stochasticity (only env sto and drift)
    {% endif %}

        //////////////////
        // process part //
        //////////////////

        {% for Ls_i in Ls.Ls_proc %}
        for(cac=0; cac<N_CAC; cac++) {
            {% for Ls_ii in Ls_i %}
            gsl_matrix_set(L, {{ forloop.parentloop.counter0 }}*N_CAC+cac, {{ forloop.counter0 }}*N_CAC+cac, {{ Ls_ii|safe }});
            {% endfor %}
        }
        {% endfor %}

        //////////////////////
        // observation part //
        //////////////////////

        ts = 0;
        {% for Ls_i in Ls.Ls_obs %}
        for(ts_unique=0; ts_unique < obs2ts[{{ forloop.counter0 }}]->n_ts_unique; ts_unique++) {
            for(stream=0; stream < obs2ts[{{ forloop.counter0 }}]->n_stream[ts_unique]; stream++) {

                for(n_cac=0; n_cac< obs2ts[{{ forloop.counter0 }}]->n_cac[ts_unique]; n_cac++) {
                    c = obs2ts[{{ forloop.counter0 }}]->cac[ts_unique][n_cac][0];
                    ac = obs2ts[{{ forloop.counter0 }}]->cac[ts_unique][n_cac][1];
                    cac = c*N_AC+ac;

                    {% for Ls_ii in Ls_i %}
                    gsl_matrix_set(L, N_PAR_SV*N_CAC+ts, {{ forloop.counter0 }}*N_CAC+cac, {{ Ls_ii|safe }});
                    {% endfor %}
                }

                ts++;
            }
        }
        {% endfor %}

        ////////////////
        // drift part //
        ////////////////

        struct s_iterator *p_it = p_data->p_it_only_drift;
        int i;
        for(i=0; i<p_it->nbtot; i++) {
            gsl_matrix_set(L, N_PAR_SV*N_CAC+N_TS + i, {{ Ls.Ls_proc.0|length }}*N_CAC + i, 1.0);
        }

    {% if command == 'deter' %}
    }
    {% endif %}

    {% endfor %}

}


/**
 * evaluate the diagonal of Qc
 */
void eval_diag_Qc(double *diag_Qc, const double *X, struct s_par *p_par, struct s_data *p_data, struct s_calc *p_calc, double t, enum plom_noises_off noises_off)
{
    // X is p_X->proj
    int cac;

    // syntaxic shortcuts
    struct s_router **routers = p_data->routers;

    //the automaticaly generated code may need these variables
    const int nn = p_calc->current_nn;
    double **par = p_par->natural;
    double ***covar = p_data->par_fixed;

    {% for command, Ls in Q.items %}

    {% if command == 'sto' %}
    if(!(noises_off & PLOM_NO_DEM_STO)) { //demographic stochasticity
    {% else %}
    } else { //no demographic stochasticity (only env sto and drift)
    {% endif %}

        for(cac=0; cac<N_CAC; cac++) {
            {% if Ls.sf %}
            double _sf[{{ Ls.sf|length }}];
            {% endif %}

            {% for sf in Ls.sf %}
            _sf[{{ forloop.counter0 }}] = {{ sf|safe }};{% endfor %}

            {% for term in Ls.diag_Qc %}
            diag_Qc[{{ forloop.counter0 }}*N_CAC+cac] = {{ term|safe }};{% endfor %}
        }

        //////////////////////////////
        // drift term (volatility^2)//
        //////////////////////////////
        struct s_iterator *p_it = p_data->p_it_only_drift;
        int i, k;
        int offset = 0;
        for(i=0; i<p_it->length; i++) {
            for(k=0; k< routers[ p_it->ind[i] ]->n_gp; k++) {
                diag_Qc[{{ Ls.diag_Qc|length }}*N_CAC + offset] = pow(par[ p_data->p_drift->ind_volatility_Xdrift[i] ][k],2);
                offset++;
            }
        }

    {% if command == 'deter' %}
    }
    {% endif %}

    {% endfor %}
}


void eval_Q(gsl_matrix *Q, const double *X, struct s_par *p_par, struct s_data *p_data, struct s_calc *p_calc, struct s_kalman_specific_data *p_kalman_specific_data, double t)
{
    // reset Q
    gsl_matrix_set_zero(Q);

    // matrices initializations
    double *diag_Qc = p_kalman_specific_data->diag_Qc;
    gsl_matrix *L = p_kalman_specific_data->L;
    gsl_matrix *LQc = p_kalman_specific_data->LQc;


    eval_diag_Qc(diag_Qc, X, p_par, p_data, p_calc, t, 0); //TODO last 0 is noises_off (from enum plom_noises)

    // L*Qc
    int row, col;
    for(row=0; row<L->size1; row++) {
        for(col=0; col<L->size2; col++) {
            gsl_matrix_set(LQc, row, col, gsl_matrix_get(L, row, col)*diag_Qc[col]);
        }
    }

    // Q = L Qc L'
    int status = gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, LQc, L, 0.0, Q);
#if FLAG_VERBOSE
    if(status) {
        fprintf(stderr, "error: %s\n", gsl_strerror (status));
    }
#endif

}
