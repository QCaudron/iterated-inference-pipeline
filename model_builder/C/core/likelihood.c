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


double get_smallest_log_likelihood(struct s_data *p_data)
{
  int n;
  double smallest_log_like = 0.0;

  for(n=0 ; n<p_data->nb_obs_nonan; n++) {
      smallest_log_like += p_data->data_ind[n]->n_nonan;
  }

  return smallest_log_like * LOG_LIKE_MIN;
}

/**
 *   Return sum of square (used for least square). The sum is computed **only** on ts != NaN
 */
double get_sum_square(struct s_X *p_X, struct s_par *p_par, struct s_data *p_data, struct s_calc *p_calc, const int n, const double t)
{
    int ts, ts_nonan;
    double ss = 0.0;

    /* syntax shortcuts */
    struct s_data_ind * p_data_ind_n = p_data->data_ind[n];

    for(ts=0; ts< p_data_ind_n->n_nonan; ts++) {
        ts_nonan = p_data_ind_n->ind_nonan[ts];
        ss += pow( p_data->data[n][ts_nonan] - obs_mean(p_X->obs[ts_nonan], p_par, p_data, p_calc, ts_nonan, n, t), 2);
    }
    
    return ss;
}

/**
 * Return sum log likelihood. The sum is computed **only** on ts != NaN
 */
double get_log_likelihood(struct s_X *p_X, struct s_par *p_par, struct s_data *p_data, struct s_calc *p_calc, const int n, const double t)
{
  int ts, ts_nonan;
  double loglike = 0.0;

  /* syntax shortcuts */
  struct s_data_ind * p_data_ind_n = p_data->data_ind[n];

  for(ts=0; ts< p_data_ind_n->n_nonan; ts++) {
      ts_nonan = p_data_ind_n->ind_nonan[ts];      
      loglike += log(likelihood(p_X->obs[ts_nonan], p_par, p_data, p_calc, ts_nonan, n, t));
  }

  return loglike;
}



/**
 *  checks for numerical issues...
 */
double sanitize_likelihood(double like)
{
    if ((isinf(like)==1) || (isnan(like)==1) || (like<0.0) ) { //error
#if FLAG_WARNING
        char str[STR_BUFFSIZE];
        sprintf(str, "error likelihood computation, like=%g", like);
        print_warning(str);
#endif
        return LIKE_MIN;
    } else {
        /*we avoid 0.0 to avoid nan when taking log and we ensure a uniform likelihood scale by making sure that everything below LIKE_MIN is LIKE_MIN */

        return (like <= LIKE_MIN) ? LIKE_MIN : like;
    }

}
