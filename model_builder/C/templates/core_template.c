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


void build_markov(struct s_calc *p)
{
unsigned int tab[N_PAR_SV+2]; //+2 for U and DU of the universes

/*automaticaly generated code: dimension of prob and inc*/
{% for x in buildmarkov %}
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
    int c, ac, n_cac_o_ts;
    double sum_prev;
    {% endif %}

    ind_obs = 0;
    ind_proj_inc = N_PAR_SV*N_CAC;

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

            sum_prev += {{ prev|safe }};
        }

        /* repeat as many times as data stream */
        for(n_stream_o_ts=0; n_stream_o_ts< (obs2ts[o])->n_stream[n_ts_unique_o]; n_stream_o_ts++) {
            p_X->obs[ind_obs++] = sum_prev;
        }
    }

    {% endfor %}
}

void step_euler_multinomial(double *X, double t, struct s_par *p_par, struct s_data *p_data, struct s_calc *p_calc)
{
    /* t is the time in unit of the data */

    struct s_obs2ts **obs2ts = p_data->obs2ts;  /* syntaxic shortcut */
    struct s_router **routers = p_data->routers;   /* syntaxic shortcut */

    int c, ac, cac, n_cac, ts, o;
    double sum_inc = 0.0;
    int offset;

    const int nn = p_calc->current_nn;

    double sum, one_minus_exp_sum;

    double **par = p_par->natural;
    double ***covar = p_data->par_fixed;


    /*automaticaly generated code:*/
    /*0-declaration of noise terms*/
    {% for n in gamma_noise %}
    double {{ n.0|safe }};{% endfor %}

    //ac=0;

    //      /*compute sum i !=j N_j^nu/d_ij^gamma*I^j */
    //      for(c=0;c<N_C;c++)
    //	{
    //    p_calc->gravity[c]=0.0;
    //    for(cc=0;cc<N_C;cc++)
    //      {
    //        if(cc !=c)
    //		{
    //            p_calc->gravity[c] += (pow(X[ORDER_I*N_CAC+cc*N_AC+ac], p_par->proc[ORDER_g_nu][cc][ac])/pow(p_data->mat_d[c][cc], p_par->proc[ORDER_gamma][cc][ac]));
    //		}
    //      }
    //    p_calc->gravity[c] *= p_par->proc[ORDER_iota][c][ac]*pow(p_data->p_t[n][c], p_par->proc[ORDER_g_mu][c][ac]);
    //	}

    for(c=0;c<N_C;c++) {
        for(ac=0;ac<N_AC;ac++) {
            cac = c*N_AC+ac;

            /*1-generate noise increments (automaticaly generated code)*/
            {% for n in gamma_noise %}
            {{ n.0|safe }} = gsl_ran_gamma(p_calc->randgsl, (DT)/ pow(par[ORDER_{{ n.1|safe }}][routers[ORDER_{{ n.1|safe }}]->map[cac]], 2), pow(par[ORDER_{{ n.1|safe }}][routers[ORDER_{{ n.1|safe }}]->map[cac]], 2))/DT;{% endfor %}

            /*2-generate process increments (automaticaly generated code)*/
            {{ print_prob|safe }}

            /*3-multinomial drawn (automaticaly generated code)*/
            {{ print_multinomial|safe }}

            /*4-update state variables (automaticaly generated code)*/
            //use inc to cache the Poisson draw as thew might be re-used for the incidence computation
            {% for draw in print_update.poisson %}
            {{ draw|safe }};{% endfor %}

            {{ print_update.Cstring|safe }}

        }/*end for on ac*/
    } /*end for on c*/

    /*compute incidence:integral between t and t+1 (automaticaly generated code)*/

    offset = 0;
    {% for eq in eq_obs_inc_markov %}
    o = {{ eq.true_ind_obs|safe }};

    for(ts=0; ts<obs2ts[o]->n_ts_unique; ts++) {
        sum_inc = 0.0;
        for(n_cac=0; n_cac<obs2ts[o]->n_cac[ts]; n_cac++) {
            c = obs2ts[o]->cac[ts][n_cac][0];
            ac = obs2ts[o]->cac[ts][n_cac][1];
            cac = c*N_AC+ac;

            sum_inc += {{ eq.right_hand_side|safe }};
        }
        X[N_PAR_SV*N_CAC +offset] += sum_inc;
        offset++;
    }

    {% endfor %}
}



int func(double t, const double X[], double f[], void *params)
{

    struct s_calc *p_calc = (struct s_calc *) params;
    struct s_par *p_par = p_calc->p_par;
    struct s_data *p_data = p_calc->p_data;
    struct s_obs2ts **obs2ts = p_data->obs2ts;
    struct s_router **routers = p_data->routers;

    int c, ac, cac, n_cac, ts, o;
    double sum_inc = 0.0;
    int offset;

    const int nn = p_calc->current_nn;

    double **par = p_par->natural;
    double ***covar = p_data->par_fixed;

    //  /*gravity model parameters*/
    //  double gravity;
    //  int cc;

    for(c=0;c<N_C;c++) {
        for(ac=0; ac<N_AC; ac++) {
            cac = c*N_AC+ac;

            //      /*compute sum i !=j N_j^nu/d_ij^gamma*I^j */
            //      gravity=0.0;
            //      for(cc=0;cc<N_C;cc++)
            //	{
            //    if(cc !=c)
            //      {
            //        gravity += (pow(X[ORDER_I*N_CAC+cc*N_AC+ac], p_par->proc[ORDER_g_nu][cc][ac])/pow(p_data->mat_d[c][cc], p_par->proc[ORDER_gamma][cc][ac]));
            //      }
            //	}
            //      gravity *= iota*pow(p_data->p_t[n][c], p_par->proc[ORDER_g_mu][c][ac]);

            /*automaticaly generated code:*/
            /*ODE system*/

            {{ print_ode|safe }}
        }
    }

    /*automaticaly generated code:*/
    /*compute incidence:integral between t and t+1*/

    offset=0;
    {% for eq in eq_obs_inc_ode %}
    o = {{ eq.true_ind_obs|safe }};

    for (ts=0; ts<obs2ts[o]->n_ts_unique; ts++) {
        sum_inc = 0.0;
        for (n_cac=0; n_cac<obs2ts[o]->n_cac[ts]; n_cac++) {
            c = obs2ts[o]->cac[ts][n_cac][0];
            ac = obs2ts[o]->cac[ts][n_cac][1];
            cac = c*N_AC+ac;

            sum_inc += {{ eq.right_hand_side|safe }};
        }

        f[N_PAR_SV*N_CAC +offset] = sum_inc;
        offset++;
    }
    {% endfor %}

    return GSL_SUCCESS;
}


double likelihood(double x, struct s_par *p_par, struct s_data *p_data, struct s_calc *p_calc, int ts)
{
    /*x is the predicted value from the model that we contrast with a time serie ts.
      Note: user should not use this function but get_log_likelihood
    */

    int n, nn;
    double t;
    n = p_calc->current_n;
    nn = p_calc->current_nn;
    t = (double) p_data->times[n];

    struct s_router **routers = p_data->routers;

    double like; /* likelihood value */

    double y = p_data->data[nn][ts];

    double **par = p_par->natural;
    double ***covar = p_data->par_fixed;

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


double obs_mean(double x, struct s_par *p_par, struct s_data *p_data, struct s_calc *p_calc, int ts)
{
  /*x is the predicted value from the model that we contrast with a time serie ts*/
  struct s_router **routers = p_data->routers;

  double mu;
  int n, nn;
  double t;
  n = p_calc->current_n;
  nn = p_calc->current_nn;
  t = (double) (nn+1);

  double **par = p_par->natural;
  double ***covar = p_data->par_fixed;

  /*automaticaly generated code*/
  mu = {{ proc_obs.mean|safe }};

  return mu;
}

double obs_var(double x, struct s_par *p_par, struct s_data *p_data, struct s_calc *p_calc, int ts)
{
  /*x is the predicted value from the model that we contrast with a time serie ts*/
  struct s_router **routers = p_data->routers;

  double var;
  int n, nn;
  double t;
  n = p_calc->current_n;
  nn = p_calc->current_nn;
  t = (double) (nn+1);

  double **par = p_par->natural;
  double ***covar = p_data->par_fixed;

  /*automaticaly generated code*/
  var = {{ proc_obs.var|safe }};

  return var;
}


double observation(double x, struct s_par *p_par, struct s_data *p_data, struct s_calc *p_calc, int ts)
{
  /*x is the predicted value from the model that we contrast with a time serie ts*/
  struct s_router **routers = p_data->routers;

  int n, nn;
  double t;
  n = p_calc->current_n;
  nn = p_calc->current_nn;
  t = (double) (nn+1);

  double **par = p_par->natural;
  double ***covar = p_data->par_fixed;

  /*return an observation of the process model*/

  /*automaticaly generated code*/
  double gsl_mu= {{ proc_obs.mean|safe }};
  double gsl_sd=sqrt({{ proc_obs.var|safe }});

  double yobs= gsl_mu+gsl_ran_gaussian(p_calc->randgsl, gsl_sd);

  return (yobs >0) ? yobs : 0.0;
}
