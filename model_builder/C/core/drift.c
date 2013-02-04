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

/**
 *   Euler Maruyama
 *   ind_drift_start and ind_drift_end are used to separate (if needed) in between drift on par_proc and drift on par_obs
 */

void compute_drift(struct s_X *p_X, struct s_par *p_par, struct s_data *p_data, struct s_calc *p_calc, double delta_t)
{

    int i, k;

    struct s_router **routers = p_data->routers;

    for(i=0; i<p_data->p_it_only_drift->length; i++) {
        struct s_drift *p_drift = p_data->drift[i];
        int ind_par_Xdrift_applied = p_drift->ind_par_Xdrift_applied;
        int ind_volatility_Xdrift = p_drift->ind_volatility_Xdrift;
        for(k=0; k< routers[ind_par_Xdrift_applied]->n_gp; k++) {
            p_X->proj[p_drift->offset +k] += p_par->natural[ ind_volatility_Xdrift ][k]*sqrt(delta_t)*gsl_ran_ugaussian(p_calc->randgsl);
        }
    }
}
