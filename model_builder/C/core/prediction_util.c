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
 * reset incidences to 0
 */
void reset_inc(struct s_X *p_X, struct s_data *p_data)
{
    double *X = p_X->proj;
    int offset = N_PAR_SV*N_CAC + p_data->p_it_only_drift->nbtot;
    int oi;
    for (oi=0; oi<N_TS_INC_UNIQUE; oi++) {
        X[offset + oi]=0.0; /*incidence*/
    }
}

/**
 * round incidences
 */
void round_inc(struct s_X *p_X, struct s_data *p_data)
{
    double *X = p_X->proj;
    int offset = N_PAR_SV*N_CAC + p_data->p_it_only_drift->nbtot;
    int oi;
    for (oi=0; oi<N_TS_INC_UNIQUE; oi++) {
        X[offset + oi]=round(X[offset + oi]); /*incidence*/
    }
}


/**
 * return the sum of the state variable of the composant @c proj of @c p_X for @c cac
 */

double sum_SV(const double *X_proj, int cac)
{
    int i;
    double current_p=0.0;

    for(i=0; i<N_PAR_SV; i++) {
        current_p += X_proj[i*N_CAC +cac];
    }

    return current_p;
}


/**
   used for euler multinomial integrarion. When duration of
   infection is close to the time step duration, the method becomes
   inacurate (the waiting time is geometric instead of
   exponential. So we ensure that the rate has the correct magnitude
   by correcting it
*/
double correct_rate(double rate, double dt)
{
    return -log(1.0-rate*dt)/dt;
}


/**
 *  take **grouped** initial condition contained in @c p_par and ungroup them in @c p_X
 */
void linearize_and_repeat(struct s_X *p_X, struct s_par *p_par, struct s_data *p_data, const struct s_iterator *p_it)
{
    double *expanded = p_X->proj;
    int i, k;
    struct s_router **routers = p_data->routers;
    int offset = 0;

    for(i=0; i< p_it->length; i++) {
        for(k=0; k< routers[ p_it->ind[i] ]->p; k++) {
            expanded[offset++] = p_par->natural[ p_it->ind[i] ][ routers[ p_it->ind[i] ]->map[k] ];
        }
    }
}

/**
 *   From proportion of initial conditons to population size.  If
 *   @c POP_SIZE_EQ_SUM_SV the last state is replaced by
 *   pop_size - sum_every_state_except_the_last.
 */
void prop2Xpop_size(struct s_X *p_X, struct s_data *p_data)
{

    double *Xpop_size = p_X->proj;
    double *pop_size_t0 = p_data->pop_size_t0;
    int i, cac;

    for (i=0; i< N_PAR_SV ; i++) {
        for (cac=0; cac<N_CAC; cac++) {
            Xpop_size[i*N_CAC+cac] = Xpop_size[i*N_CAC+cac] * pop_size_t0[cac];
            if(p_data->implementation == PLOM_PSR){ //rounding for exact methods
                Xpop_size[i*N_CAC+cac] = round(Xpop_size[i*N_CAC+cac]);
            }
        }
    }

    if (POP_SIZE_EQ_SUM_SV) {
        for (cac=0; cac<N_CAC; cac++) {
            Xpop_size[ (N_PAR_SV-1)*N_CAC + cac ] = pop_size_t0[cac] - (sum_SV(p_X->proj, cac)  - Xpop_size[ (N_PAR_SV-1)*N_CAC + cac ]);
        }
    }
}


/* load p_X->drift from drift IC (contained in p_best) */
void theta_driftIC2Xdrift(struct s_X *p_X, const theta_t *best_mean, struct s_data *p_data)
{
    int i, k;
    struct s_iterator *p_it = p_data->p_it_only_drift;
    struct s_router **routers = p_data->routers;
    struct s_drift **drift = p_data->drift;

    for(i=0; i< p_it->length; i++) {
        for(k=0; k< routers[ p_it->ind[i] ]->n_gp; k++) {
            p_X->proj[drift[i]->offset + k] = gsl_vector_get(best_mean, p_it->offset[i] +k );
        }
    }
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//prediction functions of prototype void f_prediction_ode(struct s_X *p_X, double t0, double t1, struct s_par *p_par, struct s_data *p_data, struct s_calc *p_calc)
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

plom_f_pred_t get_f_pred(enum plom_implementations implementation, enum plom_noises_off noises_off)
{ 
    if (implementation == PLOM_ODE) {
	return &f_prediction_ode;

    } else if (implementation == PLOM_SDE){

	if (noises_off == (PLOM_NO_DEM_STO | PLOM_NO_ENV_STO | PLOM_NO_DRIFT) ) {
	    return &f_prediction_ode;
	} else if (noises_off == (PLOM_NO_DEM_STO | PLOM_NO_ENV_STO) ) {
	    return &f_prediction_sde_no_dem_sto_no_env_sto;
	} else if (noises_off == (PLOM_NO_DEM_STO | PLOM_NO_DRIFT) ) {
	    return &f_prediction_sde_no_dem_sto_no_drift;
	} else if (noises_off == (PLOM_NO_ENV_STO | PLOM_NO_DRIFT) ) {
	    return &f_prediction_sde_no_env_sto_no_drift;
	} else if (noises_off == PLOM_NO_DEM_STO ) {
	    return &f_prediction_sde_no_dem_sto;
	} else if (noises_off == PLOM_NO_ENV_STO ) {
	    return &f_prediction_sde_no_env_sto;	   
	} else if (noises_off == PLOM_NO_DRIFT ) {
	    return &f_prediction_sde_no_drift;
	} else {
	    return &f_prediction_sde_full;
	}

    } else if (implementation == PLOM_PSR){
	//no_sto_env is handled within the step funciton
	if(noises_off & PLOM_NO_DRIFT){
	    return &f_prediction_psr_no_drift;
	} else {
	    return &f_prediction_psr;
	}
    }

    return NULL;
}


void f_prediction_ode(struct s_X *p_X, double t0, double t1, struct s_par *p_par, struct s_data *p_data, struct s_calc *p_calc)
{
    double t=t0;
    double h = p_X->dt; //h is the initial integration step size
    p_calc->p_par = p_par; //pass the ref to p_par so that it is available wihtin the function to integrate

    double *y = p_X->proj;

    while (t < t1) {
        int status = gsl_odeiv2_evolve_apply (p_calc->evolve, p_calc->control, p_calc->step, &(p_calc->sys), &t, t1, &h, y);
        if (status != GSL_SUCCESS) {
            char str[STR_BUFFSIZE];
            sprintf(str, "error (%d) integration time step is too large to match desired precision", status);
            print_err(str);
            exit(EXIT_FAILURE);
        }
    }
}



void f_prediction_sde_no_dem_sto_no_env_sto(struct s_X *p_X, double t0, double t1, struct s_par *p_par, struct s_data *p_data, struct s_calc *p_calc)
{
    double t = t0;

    while (t < t1) {
	step_sde_no_dem_sto_no_env_sto(p_X, t, p_par, p_data, p_calc);
        if (N_DRIFT) {
            compute_drift(p_X, p_par, p_data, p_calc);
        }

	t += p_X->dt;
    }
}

void f_prediction_sde_no_dem_sto_no_drift(struct s_X *p_X, double t0, double t1, struct s_par *p_par, struct s_data *p_data, struct s_calc *p_calc)
{
    double t = t0;

    while (t < t1) {
	step_sde_no_dem_sto(p_X, t, p_par, p_data, p_calc);

	t += p_X->dt;
    }
}


void f_prediction_sde_no_env_sto_no_drift(struct s_X *p_X, double t0, double t1, struct s_par *p_par, struct s_data *p_data, struct s_calc *p_calc)
{
    double t = t0;

    while (t < t1) {
	step_sde_no_env_sto(p_X, t, p_par, p_data, p_calc);

	t += p_X->dt;
    }
}


void f_prediction_sde_no_dem_sto(struct s_X *p_X, double t0, double t1, struct s_par *p_par, struct s_data *p_data, struct s_calc *p_calc)
{
    double t = t0;

    while (t < t1) {
	step_sde_no_dem_sto(p_X, t, p_par, p_data, p_calc);
        if (N_DRIFT) {
            compute_drift(p_X, p_par, p_data, p_calc);
        }

	t += p_X->dt;
    }
}


void f_prediction_sde_no_env_sto(struct s_X *p_X, double t0, double t1, struct s_par *p_par, struct s_data *p_data, struct s_calc *p_calc)
{
    double t = t0;

    while (t < t1) {
	step_sde_no_env_sto(p_X, t, p_par, p_data, p_calc);
        if (N_DRIFT) {
            compute_drift(p_X, p_par, p_data, p_calc);
        }

	t += p_X->dt;
    }
}

void f_prediction_sde_no_drift(struct s_X *p_X, double t0, double t1, struct s_par *p_par, struct s_data *p_data, struct s_calc *p_calc)
{
    double t = t0;

    while (t < t1) {
	step_sde_full(p_X, t, p_par, p_data, p_calc);

	t += p_X->dt;
    }
}

void f_prediction_sde_full(struct s_X *p_X, double t0, double t1, struct s_par *p_par, struct s_data *p_data, struct s_calc *p_calc)
{
    double t = t0;

    while (t < t1) {
	step_sde_full(p_X, t, p_par, p_data, p_calc);
        if (N_DRIFT) {
            compute_drift(p_X, p_par, p_data, p_calc);
        }

	t += p_X->dt;
    }
}



void f_prediction_psr(struct s_X *p_X, double t0, double t1, struct s_par *p_par, struct s_data *p_data, struct s_calc *p_calc)
{
    double t = t0;

    while (t < t1) {
        step_psr(p_X, t, p_par, p_data, p_calc);
        if (N_DRIFT) {
            compute_drift(p_X, p_par, p_data, p_calc);
        }
	t += p_X->dt;
    }

}



void f_prediction_psr_no_drift(struct s_X *p_X, double t0, double t1, struct s_par *p_par, struct s_data *p_data, struct s_calc *p_calc)
{
    double t = t0;

    while (t < t1) {
        step_psr(p_X, t, p_par, p_data, p_calc);
	t += p_X->dt;
    }
}





/**
 * Modified version of gsl_ran_multinomial to avoid a loop. We avoid
 * to recompute the total sum of p (called norm in GSL) as it will
 * always be 1.0 with simforence (no rounding error by construction)
 * @see step_euler_multinomial.
 */
void plom_ran_multinomial (const gsl_rng * r, const size_t K, unsigned int N, const double p[], unsigned int n[])
{
    size_t k;
    double sum_p = 0.0;

    unsigned int sum_n = 0;

    for (k = 0; k < K; k++) {
        if (p[k] > 0.0) {
            n[k] = gsl_ran_binomial (r, p[k] / (1.0 - sum_p), N - sum_n);
        }
        else {
            n[k] = 0;
        }

        sum_p += p[k];
        sum_n += n[k];
    }
}



/**
   computes sum mij*Ij for model with age classes.
   X_c is X[i*N_CAC + c*N_AC]
   waifw_ac is p_data->waifw_ac[ac]
*/
double foi(double *X_c, double *waifw_ac)
{
    int ac;
    double sum = 0.0;

    for (ac=0; ac< N_AC; ac++) {
        sum += waifw_ac[ac] * X_c[ac];
    }

    return sum;
}
