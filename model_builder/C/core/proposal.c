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

void apply_following_constraints(theta_t *proposed, struct s_best *p_best, struct s_data *p_data)
{
    int i, k;
    unsigned int *offset = p_data->p_it_all->offset;

    for(i=0; i<p_best->n_follow; i++){
        for(k=0; k<p_data->routers[p_best->follower[i]]->n_gp; k++){
            gsl_vector_set(proposed,
                           offset[p_best->follower[i]] + k,
                           gsl_vector_get(proposed, offset[p_best->follow[i]] + k)
                           );
        }
    }
}

/**
 * generate a new value of theta that respects the constraints on
 * the initial conditions.
 */
void propose_safe_theta_and_load_X0(theta_t *proposed, struct s_best *p_best, gsl_matrix *var, double sd_fac, struct s_par *p_par, struct s_X *p_X, struct s_data *p_data, struct s_calc *p_calc, void (*ran_proposal) (theta_t *proposed, struct s_best *p_best, gsl_matrix *var, double sd_fac, struct s_calc *p_calc))
{

    do
	{
            (*ran_proposal)(proposed, p_best, var, sd_fac, p_calc);
        }
    while (plom_check_IC_assign_theta_remainder(proposed, p_data) > 0);

    //take into account following constraints
    apply_following_constraints(proposed, p_best, p_data);

    //load_X0 (p_X->proj)
    back_transform_theta2par(p_par, proposed, p_data->p_it_par_sv, p_data);
    linearize_and_repeat(p_X, p_par, p_data, p_data->p_it_par_sv);
    prop2Xpop_size(p_X, p_data, p_calc);

    //reset dt to dt0
    p_X->dt = p_X->dt0;

    //load X drift
    theta_driftIC2Xdrift(p_X, proposed, p_data);
}



/**
 * generate a randam vector proposed. Proposed is drawn from a MVN law
 * with mean "p_best->mean" and **diagonal** matrix var "var"
 */
void ran_proposal(theta_t *proposed, struct s_best *p_best, gsl_matrix *var, double sd_fac, struct s_calc *p_calc)
{
    //here we assume a diagonal sigma_proposal so we draw iid gaussians. mvn case is a straight extension (see mvn.c)
    int i, k;

    for (i=0; i< p_best->n_to_be_estimated ; i++) {
	k = p_best->to_be_estimated[i];
	gsl_vector_set(proposed, k, gsl_vector_get(p_best->mean, k) + gsl_ran_gaussian(p_calc->randgsl, sd_fac*sqrt(gsl_matrix_get(var, k, k))));
    }
}


/**
 * return the number of errors (=number of cac where sum of IC in prop
 * of it_all_no_theta_remainder > 1.0) and assign theta remainder (if
 * it exists)
 */
int plom_check_IC_assign_theta_remainder(theta_t *proposed, struct s_data *p_data)
{
    int i, cac;
    int cnt_error = 0;

    struct s_iterator *p_it = p_data->p_it_par_sv;
    
    int ind_theta_remainder = (p_data->p_it_theta_remainder->length) ? p_data->p_it_theta_remainder->ind[0]: -1;
   
    for (cac=0; cac < N_CAC; cac++) {
	double prop_cac = 0.0;
	for (i=0; i< p_it->length ; i++) {
	    struct s_router *p_r = p_data->routers[p_it->ind[i]];
	    if(i != ind_theta_remainder){
		prop_cac += back_transform_x(gsl_vector_get(proposed, p_it->offset[i] + p_r->map[cac]), p_r->map[cac], p_r);
	    }
        }

        if(prop_cac > 1.0) {
            cnt_error++;
        } else if(ind_theta_remainder != -1) {
	    struct s_router *p_r = p_data->routers[ind_theta_remainder];
	    gsl_vector_set(proposed, p_data->p_it_theta_remainder->offset[0] +cac, (*(p_r->f))(1.0-prop_cac, p_r->min[cac], p_r->max[cac] )); //theta_remainder has a grouping of variable_population
	}
    }

    return cnt_error;
}


/**
 * compute the log of the prob of the proposal
 */

plom_err_code log_prob_proposal(double *log_like, struct s_best *p_best, theta_t *proposed, theta_t *mean, gsl_matrix *var, double sd_fac, struct s_data *p_data, int is_mvn)
{

    int i, k, offset;

    double p_tmp, Lp;
    p_tmp=0.0, Lp=0.0;

    if (is_mvn) {
        p_tmp = plom_dmvnorm(p_best, proposed, mean, var, sd_fac);
    }

    for(i=0; i<p_data->p_it_all->length; i++) {
	struct s_router *p_router = p_data->routers[p_data->p_it_all->ind[i]];
        for(k=0; k<p_router->n_gp; k++) {
	    offset = p_data->p_it_all->offset[i]+k;

	    if(p_best->is_estimated[offset]) {
                if (!is_mvn) {
                    p_tmp = gsl_ran_gaussian_pdf((gsl_vector_get(proposed, offset)-gsl_vector_get(mean, offset)), sd_fac*sqrt(gsl_matrix_get(var, offset, offset)));
                }

                /*
                  Change of variable formula:
                  Y = r(X) (r is assumed to be increasing) with inverse r_inv
                  X has f for density
                  Y has g for density

                  g(y) = f (r_inv (y)) * d r_inv(y)/dy
                  In our case, the proposal takes place on a transformed scale router.f. We know g(y) but we want f(r_inv(y))

                  we therefore divide g(y) by d r_inv(y)/dy.

                  Note that in the multivariate case, we need to
                  divide by the determinant of the Jacobian
                  matrix. However in our case the Jacobian is
                  diagonal so the determinant is the product of the
                  diagonal terms so everything generalizes nicely
                */

                p_tmp /= (*(p_router->f_inv_derivative))(gsl_vector_get(proposed, offset), p_router->min[k], p_router->max[k]);

                //check for numerical issues
                if( (isnan(p_tmp)==1) || (isinf(p_tmp)==1) || (p_tmp<0.0) ) {
		    return PLOM_ERR_LIKE;
                }

		if (!is_mvn) {
		    Lp += log(p_tmp);
		}
            }
        }
    }

    if (is_mvn) {
	Lp = log(p_tmp);
    }

    //check AGAIN for numerical issues (taking log could have created issues)
    if((isinf(Lp)==1) || (isnan(Lp)==1)) {
	return PLOM_ERR_LIKE;
    }
    

    *log_like = Lp;

    return PLOM_SUCCESS;
}
