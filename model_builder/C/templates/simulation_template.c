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

/* automatically generated code: order of the parameters */
{% for o in order.var %}
#define ORDER_{{ o|safe }} {{ forloop.counter0 }}{% endfor %}

{% for o in order.universe %}
#define ORDER_{{ o.name|safe }} {{ o.order }}{% endfor %}

{% for o in order.drift %}
#define ORDER_{{ o|safe }} {{ forloop.counter0 }}{% endfor %}

{% for o in order.data %}
#define ORDER_{{ o|safe }} {{ forloop.counter0 }}{% endfor %}

int step_lyap (double t, const double X[], double f[], void *params)
{
    struct s_calc *p_calc = (struct s_calc *) params;
    struct s_par *p_par = p_calc->p_par;  /* syntaxic shortcut */
    struct s_data *p_data = p_calc->p_data;

    int i;
    int c, ac, cac;    
    double **par = p_par->natural;

    struct s_router **routers = p_data->routers;  /* syntaxic shortcut */

    /* non linear system (automatically generated code)*/
    double _r[N_CAC][{{ step_ode_sde.caches|length }}];

    {% if step_ode_sde.sf %}
    double _sf[N_CAC][{{ step_ode_sde.sf|length }}];{% endif %}

    {% if is_drift  %}
    double drifted[N_DRIFT][N_CAC];
    {% endif %}

    for(cac=0;cac<N_CAC;cac++) {
	{% if is_drift %}
	for(i=0; i<N_DRIFT; i++){
	    int ind_drift = p_data->drift[i]->ind_par_Xdrift_applied;
	    drifted[i][cac] = par[ind_drift][routers[ind_drift]->map[cac]];
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

    /* linear system: product of jacobian matrix (DIM*DIM) per

       | y[1*DIM+0]       y[1*DIM+1] ...     y[1*DIM+(DIM-1)]   |
       | y[2*DIM+0]       y[2*DIM+1] ...     y[2*DIM+(DIM-1)]   |
       | ...                                                    |
       | y[DIM*DIM+0]     y[DIM*DIM+1] ...  y[DIM*DIM+(DIM-1)]  |

       (automaticaly generated code)
    */

    double _rj[N_CAC][{{jacobian.caches_jac_only|length}}];
    for(cac=0; cac<N_CAC; cac++){
        {% for cache in jacobian.caches_jac_only %}
        _rj[cac][{{ forloop.counter0 }}] = {{ cache|safe }};{% endfor %}
    }


    for(c=0; c<N_C; c++) {
        for(ac=0; ac<N_AC; ac++) {
            cac = c*N_AC+ac;
            for(i=0; i<(N_PAR_SV*N_CAC); i++) {
                {% for jac_n in jacobian.jac_only %}
                //printf("%d %d %d %d\n", N_PAR_SV*N_CAC, ({{ forloop.counter0 }}*N_CAC+ cac)*N_PAR_SV*N_CAC, i,  N_PAR_SV*N_CAC+ ({{ forloop.counter0 }}*N_CAC+ cac)*N_PAR_SV*N_CAC +i);
                f[N_PAR_SV*N_CAC+ ({{ forloop.counter0 }}*N_CAC+ cac)*N_PAR_SV*N_CAC +i] = {% for jac_np in jac_n %}+(_rj[cac][{{ jac_np|safe }}])*X[N_PAR_SV*N_CAC+ ({{ forloop.counter0 }}*N_CAC+ cac)*N_PAR_SV*N_CAC +i]{% endfor %};
                {% endfor %}
            }
        }
    }

    return GSL_SUCCESS;
}
