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



/**
 * get variance of the remainder for cac
 * Var(sum(Xi)) = sum(var(Xi)) + 2*sum_{i<j} (cov(Xi, Xj))
 */
double get_var_remainder(struct s_kalman_update *p, gsl_matrix *Ct, int cac)
{

    int i, j, icac, jcac;
    double sum_V = 0.0;
    double sum_C = 0.0;

    for(j=0; j< N_PAR_SV; j++) {
	jcac = j*N_CAC + cac;
	sum_V += gsl_matrix_get(Ct, jcac, jcac);
	for(i=0; i< j; i++) {
	    icac = i*N_CAC + cac;	    
	    sum_V += gsl_matrix_get(Ct, icac, jcac);
	}
    }

    return sum_V + 2.0* sum_C;
}


void print_p_hat_ekf(FILE *p_file, struct s_data *p_data, struct s_par *p_par, struct s_calc *p_calc, struct s_kalman_update *p, gsl_matrix *Ct, const int n, const double t)
{
    int i, ts;
    double x;
    double xobs, sdobs;

#if FLAG_JSON
    json_t *json_print_n = json_array();
    json_array_append_new(json_print_n, json_real(t));
#else
    fprintf(p_file, "%g,", t);
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

    /* remainder (if any) */
    if(!POP_SIZE_EQ_SUM_SV){
	int cac;
	double sumsv, rem, sd;

	for(cac=0; cac<N_CAC; cac++){
	    sumsv = 0.0;
	    for(i=0; i<N_PAR_SV; i++) {
		sumsv += gsl_vector_get(p->xk, i*N_CAC +cac);
	    }
	    	    
	    rem = gsl_spline_eval(p_calc->spline[0][cac],t,p_calc->acc[0][cac]) - sumsv;
	    sd = sqrt(get_var_remainder(p, Ct, cac));

	    x = rem - 1.96*sd;
#if FLAG_JSON
	    json_array_append_new(json_print_n, json_real(x));
#else
	    fprintf(p_file,"%g,", x);
#endif

#if FLAG_JSON
	    json_array_append_new(json_print_n, json_real(rem));
#else
	    fprintf(p_file,"%g,", rem);
#endif

	    x = rem + 1.96*sd;
#if FLAG_JSON
	    json_array_append_new(json_print_n, json_real(x));
#else
	    fprintf(p_file,"%g,", x);
#endif
	}
    }


    /* ts */
    for(ts= 0; ts<N_TS; ts++) {
	i = N_PAR_SV*N_CAC+ts;
	xobs = obs_mean(gsl_vector_get(p->xk, i), p_par, p_data, p_calc, ts, n, t);
	sdobs = sqrt(var_f_x(gsl_matrix_get(Ct, i, i), gsl_vector_get(p->xk, i), p_par, p_data, p_calc, ts, n, t));
	
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


void print_prediction_residuals_ekf(FILE *p_file_pred_res, struct s_par *p_par, struct s_data *p_data, struct s_calc *p_calc, struct s_X *p_X, struct s_kalman_update *p, gsl_matrix *Ct, const int n, const double t)
{
    int ts;
    double var_obs, pred, y, res, fCt;

#if FLAG_JSON
    json_t *root;
    json_t *json_print = json_array();
#endif

#if FLAG_JSON
    json_array_append_new(json_print, json_real(t));
#else
    fprintf(p_file_pred_res,"%g,", t);
#endif

    for(ts=0; ts<N_TS; ts++) {
	var_obs = obs_var(p_X->obs[ts], p_par, p_data, p_calc, ts, n, t);        
	y = p_data->data[n][ts];
	pred = obs_mean(gsl_vector_get(p->xk, N_PAR_SV*N_CAC+ts), p_par, p_data, p_calc, ts, n, t);
	fCt = var_f_x(gsl_matrix_get(Ct, N_PAR_SV*N_CAC + ts, N_PAR_SV*N_CAC + ts), gsl_vector_get(p->xk, N_PAR_SV*N_CAC+ts), p_par, p_data, p_calc, ts, n, t);
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
