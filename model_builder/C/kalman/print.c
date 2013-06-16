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

#include "kalman.h"

void header_prediction_residuals_ekf(FILE *p_file, struct s_data *p_data)
{
    int ts;

    fprintf(p_file, "time,");

    for(ts=0; ts<N_TS; ts++) {
	fprintf(p_file, "mean:%s,res:%s%s", p_data->ts_name[ts], p_data->ts_name[ts], (ts < (N_TS-1)) ? ",": "");
    }

    fprintf(p_file, "\n");
}

void print_p_hat_ekf(FILE *p_file, struct s_data *p_data, struct s_par *p_par, struct s_calc *p_calc, struct s_kalman_update *p, gsl_matrix *Ct, int t)
{
    int i, ts;
    double x;
    double xobs, sdobs;

#if FLAG_JSON
    json_t *json_print_n = json_array();
    json_array_append_new(json_print_n, json_integer(t));
#else
    fprintf(p_file, "%d,", t);
#endif

    /* par_sv */
    for(i=0; i< N_PAR_SV*N_CAC; i++) {
	x = gsl_vector_get(p->xk,i) - 1.96*sqrt(gsl_matrix_get(Ct, i, i));
#if FLAG_JSON
        json_array_append_new(json_print_n, json_real(x));
#else
        fprintf(p_file,"%g,", x);
#endif

        x = gsl_vector_get(p->xk,i);
#if FLAG_JSON
        json_array_append_new(json_print_n, json_real(x));
#else
        fprintf(p_file,"%g,", x);
#endif

	x = gsl_vector_get(p->xk,i) + 1.96*sqrt(gsl_matrix_get(Ct, i, i));
#if FLAG_JSON
        json_array_append_new(json_print_n, json_real(x));
#else
        fprintf(p_file,"%g,", x);
#endif
    }

    /* ts */
    for(ts= 0; ts<N_TS; ts++) {
	i = N_PAR_SV*N_CAC+ts;
	xobs = obs_mean(gsl_vector_get(p->xk, i), p_par, p_data, p_calc, ts);
	sdobs = sqrt(var_f_x(gsl_matrix_get(Ct, i, i), gsl_vector_get(p->xk, i), p_par, p_data, p_calc, ts));
	
	x = xobs - 1.96*sdobs;
#if FLAG_JSON
        json_array_append_new(json_print_n, json_real(x));
#else
        fprintf(p_file,"%g,", x);
#endif

#if FLAG_JSON
        json_array_append_new(json_print_n, json_real(xobs));
#else
        fprintf(p_file,"%g,", xobs);
#endif

	x = xobs + 1.96*sdobs;
#if FLAG_JSON
        json_array_append_new(json_print_n, json_real(x));
#else
        fprintf(p_file,"%g%s", x, (ts< (N_TS-1)) ? ",": "");
#endif
    }

    /* drift */
    for(i=N_PAR_SV*N_CAC + N_TS; i<N_PAR_SV*N_CAC+ N_TS + p_data->p_it_only_drift->nbtot; i++) {
        
	x = gsl_vector_get(p->xk,i) - 1.96*sqrt(gsl_matrix_get(Ct, i, i));
#if FLAG_JSON
        json_array_append_new(json_print_n, json_real(x));
#else
        fprintf(p_file,",%g,", x);
#endif

	x = gsl_vector_get(p->xk,i);
#if FLAG_JSON
        json_array_append_new(json_print_n, json_real(x));
#else
        fprintf(p_file,"%g,", x);
#endif

	x = gsl_vector_get(p->xk,i) + 1.96*sqrt(gsl_matrix_get(Ct, i, i));
#if FLAG_JSON
        json_array_append_new(json_print_n, json_real(x));
#else
        fprintf(p_file,"%g", x);
#endif
    }

#if FLAG_JSON
    json_t *root = json_pack("{s,s,s,o}", "flag", "hat", "msg", json_print_n);
    json_dumpf(root, stdout, JSON_COMPACT); printf("\n");
    fflush(stdout);
    json_decref(root);
#else
    fprintf(p_file, "\n");
#endif

}


void print_prediction_residuals_ekf(FILE *p_file_pred_res, struct s_par *p_par, struct s_data *p_data, struct s_calc *p_calc, struct s_X *p_X, struct s_kalman_update *p, gsl_matrix *Ct, int n, int time)
{
    int ts;
    double var_obs, pred, y, res, fCt;

#if FLAG_JSON
    json_t *root;
    json_t *json_print = json_array();
#endif

#if FLAG_JSON
    json_array_append_new(json_print, json_integer(time));
#else
    fprintf(p_file_pred_res,"%d,", time);
#endif

    for(ts=0; ts<N_TS; ts++) {
	var_obs = obs_var(p_X->obs[ts], p_par, p_data, p_calc, ts);        
	y = p_data->data[n][ts];
	pred = obs_mean(gsl_vector_get(p->xk, N_PAR_SV*N_CAC+ts), p_par, p_data, p_calc, ts);
	fCt = var_f_x(gsl_matrix_get(Ct, N_PAR_SV*N_CAC + ts, N_PAR_SV*N_CAC + ts), gsl_vector_get(p->xk, N_PAR_SV*N_CAC+ts), p_par, p_data, p_calc, ts);
	res = (y - pred)/sqrt(fCt + var_obs);

#if FLAG_JSON
	json_array_append_new(json_print, json_real(pred));
	json_array_append_new(json_print, (isnan(res)==1)? json_null() : json_real(res));
#else
	fprintf(p_file_pred_res, "%g,%g%s", pred, res, (ts < (N_TS-1)) ? ",": "\n");	
#endif
    }
    
#if FLAG_JSON
    root = json_pack("{s,s,s,o}", "flag", "pred_res", "msg", json_print);
    json_dumpf(root, stdout, JSON_COMPACT); printf("\n");
    json_decref(root);
#endif

}
