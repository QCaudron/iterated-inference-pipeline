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

#include "simplex.h"

double f_simplex(const gsl_vector *x, void *params)
{
    /* function to **minimize** */
    int n, t0, t1;
    double fitness;

    /* syntax shortcuts */
    struct s_simplex *p_params_simplex = (struct s_simplex *) params;
    struct s_data *p_data = p_params_simplex->p_data;
    struct s_best *p_best = p_params_simplex->p_best;
    struct s_calc **calc = p_params_simplex->calc;
    struct s_par *p_par = p_params_simplex->p_par;
    struct s_X *p_X = p_params_simplex->p_X;

    transfer_estimated(p_best, x, p_data);

    back_transform_theta2par(p_par, p_best->mean, p_data->p_it_all, p_data);
    linearize_and_repeat(p_X, p_par, p_data, p_data->p_it_par_sv);
    prop2Xpop_size(p_X, p_data, calc[0]);
    theta_driftIC2Xdrift(p_X, p_best->mean, p_data);

    //reset dt
    p_X->dt = p_X->dt0;

    /* if the initial conditions do not respect the constraint we set
       the log likelihood to the smallest possible value:
       smallest_log_like */


    fitness=0.0;

    if (check_IC(p_X, p_data, calc[0]) == 0) {
        for(n=0; n< p_data->nb_obs; n++) {

            t0=p_data->times[n];
            t1=p_data->times[n+1];

	    store_state_current_n(calc, n);
	    reset_inc(p_X, p_data); //reset incidence to 0
	    f_prediction_ode(p_X, t0, t1, p_par, p_data, calc[0]);

	    if(p_data->data_ind[n]->n_nonan){
		proj2obs(p_X, p_data);

		if (OPTION_LEAST_SQUARE) {
		    fitness += get_sum_square(p_X, p_par, p_data, calc[0]);
		} else {
		    fitness += get_log_likelihood(p_X, p_par, p_data, calc[0]);
		}
	    }
        } /*end of for loop on n*/
	
	if (OPTION_PRIOR) {
	    double log_prob_prior_value;
	    plom_err_code rc = log_prob_prior(&log_prob_prior_value, p_best, p_best->mean, p_best->var, p_data);
#if FLAG_VERBOSE
	    if(rc != PLOM_SUCCESS){
		print_err("error log_prob_prior computation");
	    }
#endif

	    fitness += log_prob_prior_value;
	}


    } else { //new IC do not respect the constraint:
#if FLAG_VERBOSE
        print_err("IC constraint has not been respected: pop_IC>pop_size at t=0 minimal likelihood value has been assigned");
#endif
//	sanitize_IC(p_X->proj, p_data->pop_size_t0);
        if(OPTION_LEAST_SQUARE) {
            fitness = BIG_NUMBER*p_data->nb_obs;
        } else {
            fitness = p_params_simplex->smallest_log_like;
        }
    }


    if(!OPTION_LEAST_SQUARE) {
        fitness = -fitness; //GSL simplex algo minimizes so we multiply by -1 in case of log likelihood
    }

    return fitness;
}
