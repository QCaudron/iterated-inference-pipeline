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

struct s_simplex *build_simplex(json_t *theta, enum plom_implementations implementation,  enum plom_noises_off noises_off, int general_id, int is_bayesian, double dt, double eps_abs, double eps_rel, int nb_obs)
{
  struct s_simplex *p_simplex;
  p_simplex = malloc(sizeof(struct s_simplex));
  if(p_simplex==NULL) {
        char str[STR_BUFFSIZE];
        sprintf(str, "Allocation impossible in file :%s line : %d",__FILE__,__LINE__);
        print_err(str);
        exit(EXIT_FAILURE);
  }

  json_t *settings = load_settings(PATH_SETTINGS);
  p_simplex->p_data = build_data(settings, theta, implementation, noises_off, is_bayesian, nb_obs); //also build obs2ts
  json_decref(settings);

  int size_proj = N_PAR_SV*N_CAC + p_simplex->p_data->p_it_only_drift->nbtot + N_TS_INC_UNIQUE;

  p_simplex->p_par = build_par(p_simplex->p_data);
  p_simplex->p_X = build_X(size_proj, N_TS, p_simplex->p_data, dt);
  p_simplex->p_best = build_best(p_simplex->p_data, theta);

#if FLAG_OMP
  int n_threads = omp_get_max_threads();       
#else
  int n_threads = 1;
#endif

  p_simplex->calc = build_calc(&n_threads, general_id, eps_abs, eps_rel, 1, size_proj, step_ode, p_simplex->p_data);
  p_simplex->smallest_log_like = get_smallest_log_likelihood(p_simplex->p_data->data_ind);

  return p_simplex;
}


void clean_simplex(struct s_simplex *p_simplex)
{
    clean_calc(p_simplex->calc, p_simplex->p_data);
    clean_X(p_simplex->p_X);
    clean_best(p_simplex->p_best);
    clean_par(p_simplex->p_par);
    clean_data(p_simplex->p_data);

    FREE(p_simplex);
}
