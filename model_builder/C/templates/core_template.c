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
 * Alloc memory for the psr implementation
 */
void build_psr(struct s_calc *p)
{
    unsigned int tab[N_PAR_SV+2]; //+2 for U and DU of the universes

    /*automaticaly generated code: dimension of prob and inc*/
    {% for x in psr %}
    tab[ORDER_{{x.state|safe}}] = {{x.nb_reaction|safe}};{% endfor %}

    p->prob = init2d_var_set0(N_PAR_SV+2, tab);
    p->inc = init3u_varp2_set0(N_PAR_SV+2, N_CAC, tab);

    //  p->gravity = init1d_set0(N_C);
}


void proj2obs(struct s_X *p_X, struct s_data *p_data)
{
    int o, ind_obs, ind_proj_inc;
    int n_ts_unique_o, n_stream_o_ts;

    struct s_obs2ts **obs2ts = p_data->obs2ts;

    {%if list_obs_prev %}
    int c, ac, cac, n_cac_o_ts;
    double sum_prev;
    double *X = p_X->proj;
    {% endif %}

    ind_obs = 0;
    ind_proj_inc = N_PAR_SV*N_CAC + p_data->p_it_only_drift->nbtot;

    /* extend incidence: duplicate ts with multiple data streams */
    for(o=0; o<N_OBS_INC; o++) {
        for(n_ts_unique_o=0; n_ts_unique_o< (obs2ts[o])->n_ts_unique; n_ts_unique_o++) {
            for(n_stream_o_ts=0; n_stream_o_ts< (obs2ts[o])->n_stream[n_ts_unique_o]; n_stream_o_ts++) {
                p_X->obs[ind_obs++] = p_X->proj[ind_proj_inc];
            }
            ind_proj_inc++;
        }
    }

    /* add prevalence: aggregate across c and ac to match ts and repeat to tacle multiple data streams */
    {% for prev in list_obs_prev %}

    o = N_OBS_INC + {{ forloop.counter0}};

    for(n_ts_unique_o=0; n_ts_unique_o< (obs2ts[o])->n_ts_unique; n_ts_unique_o++) {

        /* compute potentialy aggregated prevalence for the citie and age classes for the time serie */
        sum_prev = 0.0;
        for(n_cac_o_ts=0; n_cac_o_ts< (obs2ts[o])->n_cac[n_ts_unique_o]; n_cac_o_ts++) { //how many cities and age classes in this time serie
            c = (obs2ts[o])->cac[n_ts_unique_o][n_cac_o_ts][0];
            ac = (obs2ts[o])->cac[n_ts_unique_o][n_cac_o_ts][1];
	    cac = c*N_AC+ac;

            sum_prev += {{ prev|safe }};
        }

        /* repeat as many times as data stream */
        for(n_stream_o_ts=0; n_stream_o_ts< (obs2ts[o])->n_stream[n_ts_unique_o]; n_stream_o_ts++) {
            p_X->obs[ind_obs++] = sum_prev;
        }
    }

    {% endfor %}
}

//stepping functions for Poisson System with stochastic rates (psr)
void step_psr(struct s_X *p_X, double t, struct s_par *p_par, struct s_data *p_data, struct s_calc *p_calc)
{
    /* t is the time in unit of the data */

    struct s_obs2ts **obs2ts = p_data->obs2ts;  /* syntaxic shortcut */
    struct s_router **routers = p_data->routers;   /* syntaxic shortcut */

    int c, ac, cac, n_cac, ts, o;
    double sum_inc = 0.0;
    int offset;

    double sum, one_minus_exp_sum;

    double **par = p_par->natural;   

    double *X = p_X->proj;
    double dt = p_X->dt;


    /*automaticaly generated code:*/
    /*0-declaration of noise terms (if any)*/
    {% for n in white_noise %}
    double {{ n.name|safe }};{% endfor %}

    double _r[{{ step_psr.caches|length }}];
    {% if step_psr.sf %}
    double _sf[{{ step_psr.sf|length }}];{% endif %}

    {% if is_drift  %}
    int i;
    double drifted[N_DRIFT][N_CAC];
    int is_drift = ! (p_data->noises_off & PLOM_NO_DRIFT);
    {% endif %}


    for(c=0;c<N_C;c++) {
        for(ac=0;ac<N_AC;ac++) {
            cac = c*N_AC+ac;

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

            /*1-generate noise increments (if any) (automaticaly generated code)*/
            {% if white_noise %}
            if(p_data->noises_off & PLOM_NO_ENV_STO){
                {% for n in white_noise %}
                {{ n.name|safe }} = 1.0;{% endfor %}
            } else {
                {% for n in white_noise %}
                {{ n.name|safe }} = gsl_ran_gamma(p_calc->randgsl, (dt)/ pow(par[ORDER_{{ n.sd|safe }}][routers[ORDER_{{ n.sd|safe }}]->map[cac]], 2), pow(par[ORDER_{{ n.sd|safe }}][routers[ORDER_{{ n.sd|safe }}]->map[cac]], 2))/dt;{% endfor %}
            }
            {% endif %}

            /*2-generate process increments (automaticaly generated code)*/
            {% for sf in step_psr.sf %}
            _sf[{{ forloop.counter0 }}] = {{ sf|safe }};{% endfor %}

            {% for cache in step_psr.caches %}
            _r[{{ forloop.counter0 }}] = {{ cache|safe }};{% endfor %}

            {{ step_psr.code|safe }}

            /*3-multinomial drawn (automaticaly generated code)*/
            {% for draw in psr_multinomial %}
            plom_ran_multinomial(p_calc->randgsl, {{ draw.nb_exit|safe }}, (unsigned int) X[ORDER_{{ draw.state|safe }}*N_CAC+cac], p_calc->prob[ORDER_{{ draw.state|safe }}], p_calc->inc[ORDER_{{ draw.state|safe }}][cac]);{% endfor %}

            /*4-update state variables (automaticaly generated code)*/
            //use inc to cache the Poisson draw as thew might be re-used for the incidence computation
            {% for draw in step_psr.poisson %}
            {{ draw|safe }};{% endfor %}

            {{ step_psr.update|safe }}

        }/*end for on ac*/
    } /*end for on c*/

    /*compute incidence:integral between t and t+1 (automaticaly generated code)*/

    offset = N_PAR_SV*N_CAC + p_data->p_it_only_drift->nbtot;
    {% for eq in obs_inc_step_psr %}
    o = {{ eq.true_ind_obs|safe }};

    for(ts=0; ts<obs2ts[o]->n_ts_unique; ts++) {
        sum_inc = 0.0;
        for(n_cac=0; n_cac<obs2ts[o]->n_cac[ts]; n_cac++) {
            c = obs2ts[o]->cac[ts][n_cac][0];
            ac = obs2ts[o]->cac[ts][n_cac][1];
            cac = c*N_AC+ac;

            sum_inc += {{ eq.right_hand_side|safe }};
        }
        X[offset] += sum_inc;
        offset++;
    }
    {% endfor %}
}



//stepping functions for ODE and SDEs

{% for noises_off, func in step_ode_sde.func.items %}
{% if noises_off == 'ode'%}
int step_ode(double t, const double X[], double f[], void *params)
{% else %}
void step_sde_{{ noises_off }}(struct s_X *p_X, double t, struct s_par *p_par, struct s_data *p_data, struct s_calc *p_calc)
{% endif %}
{

    {% if noises_off == 'ode'%}
    struct s_calc *p_calc = (struct s_calc *) params;
    struct s_data *p_data = p_calc->p_data;
    struct s_par *p_par = p_calc->p_par;
    {% else %}
    double *X = p_X->proj;
    double dt = p_X->dt;
    double *f = p_calc->y_pred;
    {% endif %}

    struct s_obs2ts **obs2ts = p_data->obs2ts;
    struct s_router **routers = p_data->routers;

    int i, c, ac, cac, n_cac, ts, o;
    double sum_inc = 0.0;
    int offset;
    
    double **par = p_par->natural;

    double _r[N_CAC][{{ step_ode_sde.caches|length }}];

    {% if step_ode_sde.sf %}
    double _sf[N_CAC][{{ step_ode_sde.sf|length }}];{% endif %}

    {% for noise in func.proc.noises %}
    double {{ noise|safe }}[N_CAC];{% endfor %}

    {% if is_drift  %}
    double drifted[N_DRIFT][N_CAC];
    {% if noises_off != 'ode'%}
    int is_drift = ! (p_data->noises_off & PLOM_NO_DRIFT);
    {% endif %}
    {% endif %}


    for(cac=0;cac<N_CAC;cac++){
        {% if is_drift %}
        for(i=0; i<N_DRIFT; i++){
            int ind_drift = p_data->drift[i]->ind_par_Xdrift_applied;
            {% if noises_off != 'ode'%}
            if(is_drift){
                int g = routers[ind_drift]->map[cac];
                drifted[i][cac] = back_transform_x(X[p_data->drift[i]->offset + g], g, routers[ind_drift]);
            } else {
                drifted[i][cac] = par[ind_drift][routers[ind_drift]->map[cac]];
            }
            {% else %}
            drifted[i][cac] = par[ind_drift][routers[ind_drift]->map[cac]];
            {% endif %}
        }
        {% endif %}


        {% for sf in step_ode_sde.sf %}
        _sf[cac][{{ forloop.counter0 }}] = {{ sf|safe }};{% endfor %}

        {% for cache in step_ode_sde.caches %}
        _r[cac][{{ forloop.counter0 }}] = {{ cache|safe }};{% endfor %}

        {% for noise in func.proc.noises %}
        {{ noise|safe }}[cac] = sqrt(dt)*gsl_ran_ugaussian(p_calc->randgsl);{% endfor %}
    }

    for(c=0;c<N_C;c++) {
        for(ac=0; ac<N_AC; ac++) {
            cac = c*N_AC+ac;

            /*automaticaly generated code:*/
            /*ODE system*/
            {% for eq in func.proc.system %}
            f[{{eq.index}}*N_CAC+cac] {% if noises_off == 'ode'%}={% else %}= X[{{eq.index}}*N_CAC+cac] + {% endif %} {{ eq.eq|safe }};{% endfor %}
        }
    }


    //TODO: drift of the diffusion
    //for(i=N_PAR_SV*N_CAC; i<(N_PAR_SV*N_CAC + p_data->p_it_only_drift->nbtot); i++){
    //    f[i] = 0.0;
    //}

    /*automaticaly generated code:*/
    /*compute incidence:integral between t and t+1*/
    offset = N_PAR_SV*N_CAC + p_data->p_it_only_drift->nbtot;

    {% for eq in func.obs %}
    o = {{ eq.index|safe }};

    for (ts=0; ts<obs2ts[o]->n_ts_unique; ts++) {
        sum_inc = 0.0;
        for (n_cac=0; n_cac<obs2ts[o]->n_cac[ts]; n_cac++) {
            c = obs2ts[o]->cac[ts][n_cac][0];
            ac = obs2ts[o]->cac[ts][n_cac][1];
            cac = c*N_AC+ac;

            sum_inc += {{ eq.eq|safe }};
        }

        f[offset] {% if noises_off == 'ode'%}={% else %}= X[offset] + {% endif %} sum_inc;
        offset++;
    }
    {% endfor %}

    {% if noises_off == 'ode'%}
    return GSL_SUCCESS;
    {% else %}
    //y_pred (f) -> X (and we ensure that X is > 0.0)
    for(i=0; i<N_PAR_SV*N_CAC; i++){
        X[i] =  (f[i] < 0.0) ? 0.0 : f[i];
    }

    for(i=N_PAR_SV*N_CAC + p_data->p_it_only_drift->nbtot; i<N_PAR_SV*N_CAC + p_data->p_it_only_drift->nbtot +N_TS_INC_UNIQUE; i++){
        X[i] = (f[i] < 0.0) ? 0.0 : f[i];
    }
    {% endif %}

}
{% endfor %}


double likelihood(double x, struct s_par *p_par, struct s_data *p_data, struct s_calc *p_calc, const int ts, const int n, const double t)
{
    /*x is the predicted value from the model that we contrast with a time serie ts.
      Note: user should not use this function but get_log_likelihood
    */

    struct s_router **routers = p_data->routers;

    double like; /* likelihood value */

    double y = p_data->data[n][ts];

    double **par = p_par->natural;

    /*automaticaly generated code*/
    double gsl_mu = {{ proc_obs.mean|safe }};
    double gsl_sd = sqrt( {{ proc_obs.var|safe }} );

    if (y > 0.0) {
        like=gsl_cdf_gaussian_P(y+0.5-gsl_mu, gsl_sd)-gsl_cdf_gaussian_P(y-0.5-gsl_mu, gsl_sd);
    } else {
        like=gsl_cdf_gaussian_P(y+0.5-gsl_mu, gsl_sd);
    }

    return sanitize_likelihood(like);
}

double obs_mean(double x, struct s_par *p_par, struct s_data *p_data, struct s_calc *p_calc, const int ts, const int n, const double t)
{
  /*x is the predicted value from the model that we contrast with a time serie ts*/
  struct s_router **routers = p_data->routers;

  double **par = p_par->natural;
  
  /*automaticaly generated code*/
  double mu = {{ proc_obs.mean|safe }};

  return mu;
}

double obs_var(double x, struct s_par *p_par, struct s_data *p_data, struct s_calc *p_calc, const int ts, const int n, const double t)
{
  /*x is the predicted value from the model that we contrast with a time serie ts*/
  struct s_router **routers = p_data->routers;

  double **par = p_par->natural;

  /*automaticaly generated code*/
  double var = {{ proc_obs.var|safe }};

  return var;
}


double observation(double x, struct s_par *p_par, struct s_data *p_data, struct s_calc *p_calc, const int ts, const int n, const double t)
{
  /*x is the predicted value from the model that we contrast with a time serie ts*/
  struct s_router **routers = p_data->routers;

  double **par = p_par->natural;  

  /*return an observation of the process model*/

  /*automaticaly generated code*/
  double gsl_mu= {{ proc_obs.mean|safe }};
  double gsl_sd=sqrt({{ proc_obs.var|safe }});

  double yobs= gsl_mu+gsl_ran_gaussian(p_calc->randgsl, gsl_sd);

  return (yobs >0) ? yobs : 0.0;
}
