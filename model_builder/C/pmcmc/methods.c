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

#include "pmcmc.h"

/**
 * recursive expression for the average so that it can be used in
 * real time by print_hat (zmq and co...)
 */
void compute_best_traj(struct s_hat **D_p_hat_best, struct s_hat **D_p_hat_prev, struct s_hat **D_p_hat_new, struct s_data *p_data, double alpha, double m)
{
    /* recursive expression for the average so that it can be used in
       real time by print_hat (zmq and co...) */

    int n, i, ts;

    for(n=0; n<N_DATA; n++) {
        //sv
        for(i=0 ; i<N_PAR_SV*N_CAC ; i++) {
            D_p_hat_best[n]->state[i] = ((m-1.0)/m)*D_p_hat_best[n]->state[i] + (1.0/m)*(alpha*D_p_hat_new[n]->state[i] + (1.0-alpha)*D_p_hat_prev[n]->state[i]);
            D_p_hat_best[n]->state_95[i][0] = ((m-1.0)/m)*D_p_hat_best[n]->state_95[i][0] + (1.0/m)*(alpha*D_p_hat_new[n]->state_95[i][0] + (1.0-alpha)*D_p_hat_prev[n]->state_95[i][0]);
            D_p_hat_best[n]->state_95[i][1] = ((m-1.0)/m)*D_p_hat_best[n]->state_95[i][1] + (1.0/m)*(alpha*D_p_hat_new[n]->state_95[i][1] + (1.0-alpha)*D_p_hat_prev[n]->state_95[i][1]);
        }

        //ts
        for(ts=0; ts< N_TS; ts++) {
            D_p_hat_best[n]->obs[ts] = ((m-1.0)/m)*D_p_hat_best[n]->obs[ts] + (1.0/m)*(alpha*D_p_hat_new[n]->obs[ts] + (1.0-alpha)*D_p_hat_prev[n]->obs[ts]);
            D_p_hat_best[n]->obs_95[ts][0] = ((m-1.0)/m)*D_p_hat_best[n]->obs_95[ts][0] + (1.0/m)*(alpha*D_p_hat_new[n]->obs_95[ts][0] + (1.0-alpha)*D_p_hat_prev[n]->obs_95[ts][0]);
            D_p_hat_best[n]->obs_95[ts][1] = ((m-1.0)/m)*D_p_hat_best[n]->obs_95[ts][1] + (1.0/m)*(alpha*D_p_hat_new[n]->obs_95[ts][1] + (1.0-alpha)*D_p_hat_prev[n]->obs_95[ts][1]);
        }

        //drift
        for(i=0; i< p_data->p_it_only_drift->nbtot; i++) {
            D_p_hat_best[n]->drift[i] = ((m-1.0)/m)*D_p_hat_best[n]->drift[i] + (1.0/m)*(alpha*D_p_hat_new[n]->drift[i] + (1.0-alpha)*D_p_hat_prev[n]->drift[i]);
            D_p_hat_best[n]->drift_95[i][0] = ((m-1.0)/m)*D_p_hat_best[n]->drift_95[i][0] + (1.0/m)*(alpha*D_p_hat_new[n]->drift_95[i][0] + (1.0-alpha)*D_p_hat_prev[n]->drift_95[i][0]);
            D_p_hat_best[n]->drift_95[i][1] = ((m-1.0)/m)*D_p_hat_best[n]->drift_95[i][1] + (1.0/m)*(alpha*D_p_hat_new[n]->drift_95[i][1] + (1.0-alpha)*D_p_hat_prev[n]->drift_95[i][1]);
        }
    }

}


