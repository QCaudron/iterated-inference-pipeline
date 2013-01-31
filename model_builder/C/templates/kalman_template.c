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
    gsl_matrix *Q = p_kalman_specific_data->Q;
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
    Gsl_matrix_set_zero(&res2.matrix);
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


void eval_ht(gsl_vector *ht, gsl_vector *xk, struct s_par *p_par, struct s_data *p_data, struct s_calc *p_calc, int ts)
{
    /* derivative of the mean of the observation process against state variables and observed variables */

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
 * get the total number of reactions
 * ie. size of F
 */
int init_REAC(void)
{
    return N_CAC*{{ stoichiometric.rnb }};
}


/**
 * evaluate stoichiometric matrix
 */
void eval_S(gsl_matrix *S, struct s_obs2ts **obs2ts)
{
    int c, ac, cac;
    int ts, ts_unique, stream, n_cac;


    //////////////////
    // dynamic part //
    //////////////////

    {% for S_i in stoichiometric.S_sv %}
    for(c=0; c<N_C; c++) {
        for(ac=0; ac<N_AC; ac++) {
            cac = c*N_AC+ac;

            {% for S_ii in S_i %}
            gsl_matrix_set(S, {{ forloop.parentloop.counter0 }}*N_CAC+cac, {{ forloop.counter0 }}*N_CAC+cac, {{ S_ii|safe }});
            {% endfor %}
        }
    }
    {% endfor %}

    //////////////////////
    // observation part //
    //////////////////////

    ts = 0;
    {% for S_i in stoichiometric.S_ov %}
    for(ts_unique=0; ts_unique < obs2ts[{{ forloop.counter0 }}]->n_ts_unique; ts_unique++) {
        for(stream=0; stream < obs2ts[{{ forloop.counter0 }}]->n_stream[ts_unique]; stream++) {

            for(n_cac=0; n_cac< obs2ts[{{ forloop.counter0 }}]->n_cac[ts_unique]; n_cac++) {
                c = obs2ts[{{ forloop.counter0 }}]->cac[ts_unique][n_cac][0];
                ac = obs2ts[{{ forloop.counter0 }}]->cac[ts_unique][n_cac][1];
                cac = c*N_AC+ac;

                {% for S_ii in S_i %}
                gsl_matrix_set(S, N_PAR_SV*N_CAC+ts, {{ forloop.counter0 }}*N_CAC+cac, {{ S_ii|safe }});
                {% endfor %}

            }

            ts++;
        }
    }
    {% endfor %}
}


/**
 * evaluate normalized force of infection matrix
 */
void eval_F(double *F, const double *X, struct s_par *p_par, struct s_data *p_data, struct s_calc *p_calc, struct s_group ***compo_groups_drift_par_proc, double t)
{
    // X is p_X->proj
    // t is the current time in the unit of the data


    int c, ac, cac;

    // syntaxic shortcuts
    struct s_router **routers = p_data->routers;

    //the automaticaly generated code may need these variables
    const int nn = p_calc->current_nn;

    double **par = p_par->natural;
    double ***covar = p_data->par_fixed;


    // fill F
    for(c=0; c<N_C; c++) {
        for(ac=0; ac<N_AC; ac++) {
            cac = c*N_AC+ac;

            {% for F_i in stoichiometric.F %}
            F[{{ forloop.counter0 }}*N_CAC+cac] = ({{ F_i|safe }});
            {% endfor %}
        }
    }
}


/**
 * evaluate demographic stochasticity bloc of Q (so-called G)
 */
int eval_G(gsl_matrix *G, const double *X, struct s_par *p_par, struct s_data *p_data, struct s_calc *p_calc, struct s_kalman_specific_data *p_kalman_specific_data, double t)
{

    // matrices initializations
    double *F = p_kalman_specific_data->F;
    gsl_matrix *S = p_kalman_specific_data->S;
    gsl_matrix *SF = p_kalman_specific_data->SF;

    ////////////
    // eval F //
    ////////////

    struct s_group ***compo_groups_drift_par_proc = p_kalman_specific_data->compo_groups_drift_par_proc;
    eval_F(F, X, p_par, p_data, p_calc, compo_groups_drift_par_proc, t);

    //////////////
    // eval S*F //
    //////////////

    int row, col; // cell indices
    for(row=0; row<N_KAL; row++) {
        for(col=0; col<p_kalman_specific_data->N_REAC; col++) {
            gsl_matrix_set(SF, row, col, gsl_matrix_get(S, row, col)*F[col]);
        }
    }

    ///////////////////////
    // eval G = S*F*t(S) //
    ///////////////////////

    int status = gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, SF, S, 0.0, G);
#if FLAG_VERBOSE
    if(status) {
        fprintf(stderr, "error: %s\n", gsl_strerror (status));
    }
#endif

    return status;
}

void eval_Q(gsl_matrix *Q, const double *X, struct s_par *p_par, struct s_data *p_data, struct s_calc *p_calc, struct s_kalman_specific_data *p_kalman_specific_data, double t)
{

    struct s_router **routers = p_data->routers;
    struct s_iterator *p_it = p_data->p_it_only_drift;

    int i, k, cac, offset;
    const int nn = p_calc->current_nn;

    double **par = p_par->natural;
    double ***covar = p_data->par_fixed;

    // reset Q
    gsl_matrix_set_zero(Q);

    ////////////////
    // drift term //
    ////////////////

    offset = 0;
    for(i=0; i<p_it->length; i++) {
        for(k=0; k< routers[ p_it->ind[i] ]->n_gp; k++) {
            // set volatility^2 on diagonal
            gsl_matrix_set(Q, N_PAR_SV*N_CAC + N_TS + offset, N_PAR_SV*N_CAC + N_TS + offset, pow(par[ p_data->p_drift->ind_volatility_Xdrift[i] ][k],2) );
            offset++;
        }
    }

    ///////////////////////////////
    // non-correlated noise term //
    ///////////////////////////////



    {% if noise_Q.Q_proc or noise_Q.Q_obs %}
    int j;
    double term;
    {% endif %}

    {% if noise_Q.Q_proc %}

    /*
      Q_proc contains only term involving state variables. We just
      replicate those term for every cac
     */

    {% if noise_Q.sf %}
    double _sf[N_CAC][{{ noise_Q.sf|length }}];
    {% endif %}

    for (cac=0; cac<N_CAC; cac++) {
        {% for sf in noise_Q.sf %}
        _sf[cac][{{ forloop.counter0 }}] = {{ sf|safe }};{% endfor %}

        {% for x in noise_Q.Q_proc %}
        i = {{ x.i }} * N_CAC + cac;
        j = {{ x.j }} * N_CAC + cac;
        term = {{ x.sign}} pow({{ x.rate|safe }}, 2);
        gsl_matrix_set(Q, i, j, term + gsl_matrix_get(Q, i, j));
        {% endfor %}
    }
    {% endif %}



    {% if noise_Q.Q_obs %}
    /*
      Q_obs contains only term involving at least one observed
      variable. The expansion is difficult: Q is expressed in terms of
      time series and not observed variable we only know the
      relationship between state variable and observed variable
      (that's what Q_obs gives us). We need to recreate time series
      from index of observed variable (x.to.ind and x.from.ind).

      for one given observed variable, s_obs2ts gives us the time
      series involved (and the cac involved)
    */

    int ts_unique, ts_unique_j;
    int stream, stream_j;
    int ts, ts_j;
    int cac_j;
    int n_cac, n_cac_j;

    struct s_obs2ts **obs2ts = p_data->obs2ts;
    struct s_obs2ts *p_obs2ts, *p_obs2ts_j;

    {% for x in noise_Q.Q_obs %}

    {% if x.i.is_obs and not x.j.is_obs or x.j.is_obs and not x.i.is_obs   %}

    /*
      x contains a state variable and an observed variable the
      observed variable correspond to different ts (given by
      obs2ts). For each ts, obs2ts also gives us the cac observed.
    */

    p_obs2ts = obs2ts[{% if x.i.is_obs %}{{ x.i.ind }}{% else %}{{ x.j.ind }}{% endif %}];
    ts=0;
    for(ts_unique=0; ts_unique < p_obs2ts->n_ts_unique; ts_unique++) {
        for(stream=0; stream < p_obs2ts->n_stream[ts_unique]; stream++) {
            for(n_cac=0; n_cac< p_obs2ts->n_cac[ts_unique]; n_cac++) {
                cac = p_obs2ts->cac[ts_unique][n_cac][0]*N_AC + p_obs2ts->cac[ts_unique][n_cac][1];

                term = {{ x.sign}} pow({{ x.rate|safe }}, 2);

                i = {% if not x.i.is_obs %}{{ x.i.ind }} * N_CAC + cac{% else %}N_PAR_SV*N_CAC + p_obs2ts->offset + ts{% endif %};
                j = {% if not x.j.is_obs %}{{ x.j.ind }} * N_CAC + cac{% else %}N_PAR_SV*N_CAC + p_obs2ts->offset + ts{% endif %};

                gsl_matrix_set(Q, i, j, term + gsl_matrix_get(Q, i, j));
            }
            ts++;
        }
    }

    {% else %}

    /*
      x contains 2 observed variables: complicated case, we have to
      find the common cac in between x.i.ind and x.j.ind
    */

    ts=0;
    p_obs2ts = obs2ts[{{ x.i.ind }}];
    p_obs2ts_j = obs2ts[{{ x.j.ind }}];

    for(ts_unique=0; ts_unique < p_obs2ts->n_ts_unique; ts_unique++) {
        for(stream=0; stream < p_obs2ts->n_stream[ts_unique]; stream++) {

            ts_j=0;
            for(ts_unique_j=0; ts_unique_j < p_obs2ts_j->n_ts_unique; ts_unique_j++) {
                for(stream_j=0; stream_j < p_obs2ts_j->n_stream[ts_unique_j]; stream_j++) {

                    i = N_PAR_SV*N_CAC + p_obs2ts->offset + ts;
                    j = N_PAR_SV*N_CAC + p_obs2ts_j->offset + ts_j;

                    //we determine the cac in common between the 2 ts (indexed i and j in Q)
                    for(n_cac=0; n_cac< p_obs2ts->n_cac[ts_unique]; n_cac++) {
                        cac = p_obs2ts->cac[ts_unique][n_cac][0]*N_AC + p_obs2ts->cac[ts_unique][n_cac][1];

                        for(n_cac_j=0; n_cac_j< p_obs2ts_j->n_cac[ts_unique_j]; n_cac_j++) {
                            cac_j = p_obs2ts_j->cac[ts_unique_j][n_cac_j][0]*N_AC + p_obs2ts_j->cac[ts_unique_j][n_cac_j][1];

                            if(cac == cac_j){
                                //x.rate is a function of cac
                                term = {{ x.sign}} pow({{ x.rate|safe }}, 2);
                                gsl_matrix_set(Q, i, j, term + gsl_matrix_get(Q, i, j));
                            }
                        }
                    }

                    ts_j++;
                }
            }

            ts++;
        }
    }

    {% endif %}

    {% endfor %}
    {% endif %}


    ////////////////////////////////////
    // demographic stochasticity term //
    ////////////////////////////////////

    // if demographic stochasticity is taken into account
    if(!COMMAND_DETER) {
        gsl_matrix *G = p_kalman_specific_data->G;
        eval_G(G, X, p_par, p_data, p_calc, p_kalman_specific_data, t);
        gsl_matrix_add(Q, G); // Q <- Q+G
    }
}
