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
int step_ode_ekf(double t, const double X[], double f[], void *params)
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

    gsl_matrix_const_view Ct   = gsl_matrix_const_view_array(&X[N_PAR_SV*N_CAC + p_data->p_it_only_drift->nbtot + N_TS_INC_UNIQUE], N_KAL, N_KAL);
    gsl_matrix_view ff = gsl_matrix_view_array(&f[N_PAR_SV*N_CAC + p_data->p_it_only_drift->nbtot + N_TS_INC_UNIQUE], N_KAL, N_KAL);

    int i, c, ac, cac, n_cac, ts, o;
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

    double _r[N_CAC][{{step_ode_sde.caches|length}}];

    {% if step_ode_sde.sf %}
    double _sf[N_CAC][{{step_ode_sde.sf|length}}];{% endif %}

    {% if is_drift  %}
    double drifted[N_DRIFT][N_CAC];
    int is_drift = ! (p_data->noises_off & PLOM_NO_DRIFT);
    {% endif %}

    for(cac=0;cac<N_CAC;cac++) {
	{% if is_drift %}
	for(i=0; i<N_DRIFT; i++){
	    int ind_drift = p_data->drift[i]->ind_par_Xdrift_applied;
	    if(is_drift){
		int g = routers[ind_drift]->map[cac];
		drifted[i][cac] = back_transform_x(X[p_data->drift[i]->offset + g], g, routers[ind_drift]);
	    } else {
		drifted[i][cac] = par[ind_drift][routers[ind_drift]->map[cac]];
	    }
	}
	{% endif %}

        {% for sf in step_ode_sde.sf %}
        _sf[cac][{{ forloop.counter0 }}] = {{ sf|safe }};{% endfor %}

        {% for cache in step_ode_sde.caches %}
        _r[cac][{{ forloop.counter0 }}] = {{ cache|safe }};{% endfor %}
    }

    for (c=0;c<N_C;c++) {
        for(ac=0; ac<N_AC; ac++) {
            cac = c*N_AC+ac;

	    {% for eq in step_ode_sde.func.ode.proc.system %}
	    f[{{eq.index}}*N_CAC+cac] = {{ eq.eq|safe }};{% endfor %}
        }
    }

    {% if is_drift  %}
    //TODO: drift of the diffusion
    for(i=N_PAR_SV*N_CAC; i<(N_PAR_SV*N_CAC + p_data->p_it_only_drift->nbtot); i++){
        f[i] = 0.0;
    }
    {% endif %}

    /*automaticaly generated code:*/
    /*compute incidence:integral between t and t+1*/

    offset = N_PAR_SV*N_CAC + p_data->p_it_only_drift->nbtot;
    {% for eq in step_ode_sde.func.ode.obs %}
    o = {{ eq.index|safe }};

    for (ts=0; ts<obs2ts[o]->n_ts_unique; ts++) {
        sum_inc = 0.0;
        for(n_cac=0; n_cac<obs2ts[o]->n_cac[ts]; n_cac++) {
            c = obs2ts[o]->cac[ts][n_cac][0];
            ac = obs2ts[o]->cac[ts][n_cac][1];
            cac = c*N_AC+ac;

            sum_inc += {{ eq.eq|safe }};
        }

        f[offset] = sum_inc;
        offset++;
    }
    {% endfor %}

    ////////////////
    // covariance //
    ////////////////

    // evaluate Q and jacobian
    p_kalman_specific_data->eval_Q(Q, X, p_par, p_data, p_calc, p_kalman_specific_data, t);
    eval_jac(Ft, X, p_par, p_data, p_calc, compo_groups_drift_par_proc, t);

    // compute Ft*Ct+Ct*Ft'+Q 
    //here Ct is symmetrical and transpose(FtCt) == transpose(Ct)transpose(Ft) == Ct transpose(Ft)
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, Ft, &Ct.matrix, 0.0, FtCt);
    for(i=0; i< ff.matrix.size1; i++){
	for(c=0; c< ff.matrix.size2; c++){
	    gsl_matrix_set(&ff.matrix, 
			   i,
			   c, 
			   gsl_matrix_get(FtCt, i, c) + gsl_matrix_get(FtCt, c, i) + gsl_matrix_get(Q, i, c));

	}
    }

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
    int i, d, g;
    struct s_drift **drift =  p_data->drift;
    double drifted[N_DRIFT][N_CAC];
    int is_drift = ! (p_data->noises_off & PLOM_NO_DRIFT);
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
	{% if is_drift %}
	for(i=0; i<N_DRIFT; i++){
	    int ind_drift = drift[i]->ind_par_Xdrift_applied;
	    if(is_drift){
		g = routers[ind_drift]->map[cac];
		drifted[i][cac] = back_transform_x(X[drift[i]->offset + g], g, routers[ind_drift]);
	    } else {
		drifted[i][cac] = par[ind_drift][routers[ind_drift]->map[cac]];
	    }
	}
	{% endif %}


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


    {% if is_drift %}
    if(is_drift){

	//third non null part of the jacobian matrix: derivative of the ODE (excluding the observed variable) against the drift variable (automaticaly generated code)
	//non null only for 'cac' present in group of Xdrift: compo_groups_drift_par_proc gives these 'cac'
	{% for jac_i in jacobian.jac_drift %}
	d = 0;
	{% for jac_ii in jac_i %}
	for(g=0; g< routers[ drift[{{ forloop.counter0 }}]->ind_par_Xdrift_applied ]->n_gp; g++) {
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
		for(g=0; g< routers[ drift[{{ forloop.counter0 }}]->ind_par_Xdrift_applied ]->n_gp; g++) {
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
    {% endif %}

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



{% for noises_off, tpl in calc_Q.items %}
/**
 * evaluate Q ({{ noises_off }})
 */
void eval_Q_{{ noises_off }}(gsl_matrix *Q, const double *X, struct s_par *p_par, struct s_data *p_data, struct s_calc *p_calc, struct s_kalman_specific_data *p_kalman_specific_data, double t)
{

    struct s_router **routers = p_data->routers;
    //the automaticaly generated code may need these variables
    const int nn = p_calc->current_nn;
    double **par = p_par->natural;
    double ***covar = p_data->par_fixed;
    int i, k, cac, offset;

    {% if is_drift  %}
    int is_drift = ! (p_data->noises_off & PLOM_NO_DRIFT);
    {% endif %}

    //////////////////////////////////////////////////////////////
    // demographic stochasticity and white noise terms (if any) //
    //////////////////////////////////////////////////////////////

    {% if tpl.Q_proc or tpl.Q_obs %}
    int j;
    double term;

    {% if is_drift  %}
    double drifted[N_DRIFT][N_CAC];
    {% endif %}

    {% if tpl.sf %}
    double _sf[N_CAC][{{ tpl.sf|length }}];
    {% endif %}

    /*
      Q_proc contains only term involving state variables. We just
      replicate those term for every cac
     */

    for (cac=0; cac<N_CAC; cac++) {

	{% if is_drift %}
	for(i=0; i<N_DRIFT; i++){
	    int ind_drift = p_data->drift[i]->ind_par_Xdrift_applied;
	    if(is_drift){
		int g = routers[ind_drift]->map[cac];
		drifted[i][cac] = back_transform_x(X[p_data->drift[i]->offset + g], g, routers[ind_drift]);
	    } else {
		drifted[i][cac] = par[ind_drift][routers[ind_drift]->map[cac]];
	    }
	}
	{% endif %}

        {% for sf in tpl.sf %}
        _sf[cac][{{ forloop.counter0 }}] = {{ sf|safe }};{% endfor %}

        {% for x in tpl.Q_proc %}
        i = {{ x.i }} * N_CAC + cac;
        j = {{ x.j }} * N_CAC + cac;
        term = {{ x.rate|safe }};

        gsl_matrix_set(Q, i, j, term);

	{% if x.i != x.j %}
        gsl_matrix_set(Q, j, i, term);
	{% endif %}

        {% endfor %}
    }

    {% endif %}



    {% if tpl.Q_obs %}
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

    {% for x in tpl.Q_obs %}

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

		term = {{ x.rate|safe }};

                gsl_matrix_set(Q, i, j, term);

		{% if x.i.ind != x.j.ind %}
		gsl_matrix_set(Q, j, i, term);
		{% endif %}

            }
            ts++;
        }
    }

    {% else %}

    /*
      x contains 2 observed variables: complicated case, we have to
      find the common cac in between x.i.ind and x.j.ind
    */


    {% if x.i.ind != x.j.ind %}

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
				term = {{ x.rate|safe }};
                                gsl_matrix_set(Q, i, j, term);

				{% if x.i.ind != x.j.ind %}
				gsl_matrix_set(Q, j, i, term);
				{% endif %}
                            }
                        }
                    }

                    ts_j++;
                }
            }
            ts++;
        }
    }


    {% else %}

    ts=0;
    p_obs2ts = obs2ts[{{ x.i.ind }}];

    for(ts_unique=0; ts_unique < p_obs2ts->n_ts_unique; ts_unique++) {
        for(stream=0; stream < p_obs2ts->n_stream[ts_unique]; stream++) {

	    i = N_PAR_SV*N_CAC + p_obs2ts->offset + ts;

	    for(n_cac=0; n_cac< p_obs2ts->n_cac[ts_unique]; n_cac++) {
		cac = p_obs2ts->cac[ts_unique][n_cac][0]*N_AC + p_obs2ts->cac[ts_unique][n_cac][1];
		gsl_matrix_set(Q, i, i, {{ x.rate|safe }});		
	    }   
	    ts++;
	}
    }

    {% endif %}

    {% endif %}

    {% endfor %}
    {% endif %}


    {% if is_drift  %}
    ///////////////////////////////////////////
    // drift term (volatility^2 on diagonal) //
    ///////////////////////////////////////////
    if(is_drift){
	struct s_iterator *p_it = p_data->p_it_only_drift;
	offset = N_PAR_SV*N_CAC + N_TS;
	for(i=0; i<p_it->length; i++) {
	    for(k=0; k< routers[ p_it->ind[i] ]->n_gp; k++) {
		gsl_matrix_set(Q, offset, offset, pow(par[ p_data->drift[i]->ind_volatility_Xdrift ][k], 2) );
		offset++;
	    }
	}
    }
    {% endif %}

}

{% endfor %}
