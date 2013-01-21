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

void ensure_cst_pop_size(struct s_data *p_data)
{
    int nn, cac;

    {% if 'mu_b' in order.data and 'mu_d' in order.data %}
    print_warning("variable birth and death rate (mu_b and mu_d) detected in covariates. mu_d have been set to mu_b to ensure a constant population size to analyze the attractor");

    for (nn=0; nn < N_DATA_PAR_FIXED; nn++) {
        for (cac=0; cac < N_CAC; cac++) {
            p_data->par_fixed[ORDER_mu_d][nn][cac] = p_data->par_fixed[ORDER_mu_b][nn][cac];
        }
    }
    {% endif %}
}


int func_lyap (double t, const double X[], double f[], void *params)
{

    struct s_calc *p_calc = (struct s_calc *) params;
    struct s_par *p_par = p_calc->p_par;  /* syntaxic shortcut */
    struct s_data *p_data = p_calc->p_data;

    int i;
    int c, ac, cac;
    const int nn = p_calc->current_nn;
    double **par = p_par->natural;
    double ***covar = p_data->par_fixed;

    struct s_router **routers = p_data->routers;  /* syntaxic shortcut */

    {% if current_p %}
    for(c=0;c<N_C;c++){
        for(ac=0;ac<N_AC;ac++){
            p_par->current_p[c][ac]=get_current_pop_size(X,c,ac);
        }
    }
    {% endif %}

    /* non linear system (automaticaly generated code)*/


    double _r[N_CAC][{{print_ode.caches|length}}];
    for(cac=0;cac<N_CAC;cac++) {
        {% for cache in print_ode.caches %}
        _r[cac][{{ forloop.counter0 }}] = {{ cache|safe }};{% endfor %}
    }

    for(c=0;c<N_C;c++){
        for(ac=0; ac<N_AC; ac++){
            cac = c*N_AC+ac;
            {{ print_ode.sys|safe }}
        }
    }


    /* linear system: product of jacobian matrix (DIM*DIM) per

       | y[1*DIM+0]       y[1*DIM+1] ...     y[1*DIM+(DIM-1)]   |
       | y[2*DIM+0]       y[2*DIM+1] ...     y[2*DIM+(DIM-1)]   |
       | ...                                                    |
       | y[DIM*DIM+0]     y[DIM*DIM+1] ...  y[DIM*DIM+(DIM-1)]  |

       (automaticaly generated code)
    */


    {% for jac_n in jacobian.jac %}
    for(c=0; c<N_C; c++) {
        for(ac=0; ac<N_AC; ac++) {
            cac = c*N_AC+ac;
            for(i=0; i<(N_PAR_SV*N_CAC); i++) {
                //printf("%d %d %d %d\n", N_PAR_SV*N_CAC, ({{ forloop.counter0 }}*N_CAC+ cac)*N_PAR_SV*N_CAC, i,  N_PAR_SV*N_CAC+ ({{ forloop.counter0 }}*N_CAC+ cac)*N_PAR_SV*N_CAC +i);
                f[N_PAR_SV*N_CAC+ ({{ forloop.counter0 }}*N_CAC+ cac)*N_PAR_SV*N_CAC +i] = {% for jac_np in jac_n %}+({{ jac_np|safe }})*X[N_PAR_SV*N_CAC+ ({{ forloop.counter0 }}*N_CAC+ cac)*N_PAR_SV*N_CAC +i]{% endfor %};
            }
        }
    }
    {% endfor %}

    return GSL_SUCCESS;
}
