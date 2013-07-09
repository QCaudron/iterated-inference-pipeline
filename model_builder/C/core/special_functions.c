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

double terms_forcing(double amplitude, double time, struct s_data *p_data, int cac)
{
  /*time is in the same unit as data, !!! school terms are in year !!! */

  double s_t;
  int k;

  int isholliday=0;
  double time_in_year = time / ONE_YEAR;
  double part_year = time_in_year - floor(time_in_year);

  for(k=0; k< p_data->n_terms[cac] ; k++) {
      if( ((part_year >= p_data->school_terms[cac][k][0]) && (part_year <= p_data->school_terms[cac][k][1])) ) {
	  isholliday=1;
	  break;
      }
  }

  if (isholliday) { //holidays
      s_t = amplitude/(1.0+amplitude*(1.0-p_data->prop_school[cac]));
  } else { //schools are in
      s_t = 1.0/(1.0+amplitude*(1.0-p_data->prop_school[cac]));
  }
  
  return s_t;
}

/**
 * The Heaviside step function or unit step function
 * x is typicaly t-t_intervention
 */
double heaviside(double x)
{    
    return (x < 0.0) ? 0.0 : 1.0;
}

/**
 * The ramp function 
 */
double ramp(double x)
{
    return (x >= 0) ? x : 0.0;
}
