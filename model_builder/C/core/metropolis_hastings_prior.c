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
 * return accepted (1) or rejected (0)
 */
int metropolis_hastings(struct s_best *p_best, struct s_likelihood *p_like, double *alpha, struct s_data *p_data, struct s_calc *p_calc, gsl_matrix *var, double sd_fac, int is_mvn)
{
    double ran;

    double Lproposal_new, Lproposal_prev, Lprior_new, Lprior_prev;
    plom_err_code rc_proposal_new = log_prob_proposal(&Lproposal_new, p_best, p_best->proposed, p_best->mean, var, sd_fac, p_data,  is_mvn); /* q{ theta* | theta(i-1) }*/
    plom_err_code rc_proposal_prev = log_prob_proposal(&Lproposal_prev, p_best, p_best->mean, p_best->proposed, var, sd_fac, p_data, is_mvn); /* q{ theta(i-1) | theta* }*/
    plom_err_code rc_prior_new = log_prob_prior(&Lprior_new, p_best, p_best->proposed, var, p_data); /* p{theta*} */
    plom_err_code rc_prior_prev = log_prob_prior(&Lprior_prev, p_best, p_best->mean, var, p_data);  /* p{theta(i-1)} */

    if( (rc_proposal_new == PLOM_SUCCESS) && (rc_proposal_prev == PLOM_SUCCESS) && (rc_prior_new == PLOM_SUCCESS) && (rc_prior_prev == PLOM_SUCCESS) ) {

        // ( p{theta*}(y)  p{theta*} ) / ( p{theta(i-1)}(y) p{theta(i-1)} )  *  q{ theta(i-1) | theta* } / q{ theta* | theta(i-1) }
        *alpha = exp( (p_like->Llike_new - p_like->Llike_prev + Lproposal_prev - Lproposal_new + Lprior_new - Lprior_prev) );

        ran = gsl_ran_flat(p_calc->randgsl, 0.0, 1.0);

        if(ran < *alpha) {
            return 1; //accepted
        }
    } else {
        *alpha = 0.0;
        return 0; //rejected
    }

    return 0;
}


/**
 * log_prob_prior always compute a logged value (in the same way as
 * sanitize likelihood would have done) even if it fails (by that I
 * means it doesn't immediatly return on failure). This is usefull for
 * the --prior option.
 */
plom_err_code log_prob_prior(double *log_like, struct s_best *p_best, gsl_vector *mean, gsl_matrix *var, struct s_data *p_data)
{
    int i, k, offset;

    int is_err = 0;

    double p_tmp, Lp;
    p_tmp=0.0, Lp=0.0;
    double back_transformed; //prior are on the natural scale, so we transform the parameter

    for(i=0; i<p_data->p_it_all->length; i++) {
	struct s_router *p_router = p_data->routers[p_data->p_it_all->ind[i]];
        for(k=0; k<p_router->n_gp; k++) {
	    offset = p_data->p_it_all->offset[i] + k;
	    if(p_best->is_estimated[offset]) {
                back_transformed = (*(p_router->f_inv[k]))(gsl_vector_get(mean, offset), p_router->min[k], p_router->max[k]);
                p_tmp = (*(p_best->prior[offset]))(back_transformed, p_best->par_prior[offset][0], p_best->par_prior[offset][1]);

                //check for numerical issues
                if( (isnan(p_tmp)==1) || (isinf(p_tmp)==1) || (p_tmp<0.0) ) {
		    is_err =1;    
		    p_tmp = LIKE_MIN;
		} else {
		    p_tmp = (p_tmp <= LIKE_MIN) ? LIKE_MIN : p_tmp;
		}

                Lp += log(p_tmp); 
            }
        }
    }

    *log_like = Lp;
    return (is_err) ? PLOM_ERR_LIKE: PLOM_SUCCESS;
}


double normal_prior(double x, double min, double max)
{
    double mean = (max+min)/2.0;
    double sd = (max-min)/4.0;
    return gsl_ran_gaussian_pdf((x-mean), sd);
}


/**
 * this function is not exactly a flat function between min an max,
 * but a twice-differentiable function, to prevent optimization
 * methods from bumping into boundaries and creating numerical
 * problems
 *
 * TODO: give users control on delta
 */
double pseudo_unif_prior(double x, double min, double max)
{
  double delta = 0.01*(max-min);
  double alpha = 1.0/(max-min-2.0/3.0*delta);

  double a = -alpha/(pow(delta,2.0));
  double b = -2.0*(min+delta)*a;
  double c = alpha + a*pow((min+delta),2.0);

  if (x<=min || x>=max){
    return 0;
  }
  else if (x < min + delta) {
    return a*pow(x,2.0) + b*x + c;
  }
  else if (x < max - delta) {
    return alpha;
  }
  else {
    return a*pow(x-(max-min-2.0*delta),2.0) + b*(x - (max-min-2.0*delta)) + c;
  }

}
