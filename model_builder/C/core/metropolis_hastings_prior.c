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

int metropolis_hastings(struct s_best *p_best, struct s_likelihood *p_like, double *alpha, struct s_data *p_data, struct s_calc *p_calc, gsl_matrix *var, double sd_fac, int is_mvn)
{
    double ran;

    //log_prob_proposal return -1 if numerical error during computation
    double Lproposal_new = log_prob_proposal(p_best, p_best->proposed, p_best->mean, var, sd_fac, p_data,  is_mvn); /* q{ theta* | theta(i-1) }*/
    double Lproposal_prev = log_prob_proposal(p_best, p_best->mean, p_best->proposed, var, sd_fac, p_data, is_mvn); /* q{ theta(i-1) | theta* }*/


    if( (Lproposal_new >0.0) && (Lproposal_prev >0.0) ) {
	//    if(check_prior(p_best, p_best->mean, var, p_data) && (Lproposal_new >0.0) && (Lproposal_prev >0.0) ) {
        double Lprior_new = log_prob_prior(p_best, p_best->proposed, var, p_data); /* p{theta*} */
        double Lprior_prev = log_prob_prior(p_best, p_best->mean, var, p_data);  /* p{theta(i-1)} */



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

int check_prior(struct s_best *p_best, gsl_vector *mean, gsl_matrix *var, struct s_data *p_data)
{
  /*check if mean is compatible with the prior (i.e if prob_prior >0)*/

  int i, k;
  double p=0.0;
  double back_transformed; //prior are on the natural scale, so we transform the parameter

  struct s_router **routers = p_data->routers;
  int offset = 0;

  for(i=0; i<p_data->p_it_all->length; i++) {
      for(k=0; k<routers[i]->n_gp; k++) {
          if(p_best->is_estimated[offset]) {
              back_transformed = (*(routers[i]->f_inv[k]))(gsl_vector_get(mean, offset), routers[i]->min[k], routers[i]->max[k]);

              p = (*(p_best->prior[offset]))(back_transformed, p_best->par_prior[offset][0], p_best->par_prior[offset][1]);
              if(p == 0.0) {
                  return 0; //FALSE
              }
          }
          offset++;
      }
  }

  return 1; //TRUE
}


double log_prob_prior(struct s_best *p_best, gsl_vector *mean, gsl_matrix *var, struct s_data *p_data)
{
    int i, k;

    double p_tmp, Lp;
    p_tmp=0.0, Lp=0.0;
    double back_transformed; //prior are on the natural scale, so we transform the parameter

    struct s_router **routers = p_data->routers;
    int offset = 0;

    for(i=0; i<p_data->p_it_all->length; i++) {
        for(k=0; k<routers[i]->n_gp; k++) {
	    if(p_best->is_estimated[offset]) {
                back_transformed = (*(routers[i]->f_inv[k]))(gsl_vector_get(mean, offset), routers[i]->min[k], routers[i]->max[k]);
                p_tmp = (*(p_best->prior[offset]))(back_transformed, p_best->par_prior[offset][0], p_best->par_prior[offset][1]);
                Lp += log(sanitize_likelihood(p_tmp));
            }
            offset++;
        }
    }

    return(Lp);
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
