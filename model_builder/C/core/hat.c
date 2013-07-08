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
 *  fill hat_95[2] with the 95% confidence interval (lower value in
 *  hat_95[0] and upper one in hat_95[1]). to_be_sorted is an array
 *  of the J particle to be sorted and weights their weights.

 * NOTE: if weights is NULL, 1.0/J will be assumed as a weight
 */

void get_CI95(double *hat_95, const double *to_be_sorted, size_t *index_sorted, double *weights)
{
    int k;
    double weight_cum;

    double invJ = 1.0/ ((double) J);

    //get the index of to_be_sorted.
    gsl_sort_index(index_sorted, to_be_sorted, 1, J); //to_be_sorted is not modified (i.e. not sorted in place), the index of the sorting are put in index_sorted.

    //cumulate sorted weight until we reach 2.5 % and take the corresponding value in to_be_sorted
    k=0;
    weight_cum = 0.0;
    while(weight_cum < 0.025) {
        weight_cum += (weights) ? weights[index_sorted[k]]: invJ;
        k++;
    }
    hat_95[0] = to_be_sorted[index_sorted[((k-1) <0) ? 0 : k-1]];

    //cumulate sorted weight until we reach 97.5 % and take the corresponding value in to_be_sorted
    k=0;
    weight_cum = 0.0;
    while(weight_cum < 0.975) {
        weight_cum += (weights) ? weights[index_sorted[k]]: invJ;
        k++;
    }
    hat_95[1] = to_be_sorted[index_sorted[((k-1) <0) ? 0 : k-1]];
}


/**
 *  compute estimation of the state variable (p_hat->state) and the
 *  observation (D_p_hat->obs). The 95% confidence interval (CI) are
 *  also computed and stored in D_p_hat->state_95 and
 *  D_p_hat->obs_95. Note that the estimations are computed by a
 *  weighted average (each value is weighted by it's likelihood
 *  value).
 */

void compute_hat(struct s_X **J_p_X, struct s_par *p_par, struct s_data *p_data, struct s_calc **calc, struct s_hat *p_hat, double *weights, const int n, const double t)
{
    //TODO weights = 1/J when no information

    int j, i, k, ts;
    int thread_id;

    struct s_router **routers = p_data->routers;

    /* par_sv */
#if FLAG_OMP
#pragma omp parallel for private(thread_id, j)
#endif
    for(i=0; i<(N_PAR_SV*N_CAC); i++) {
#if FLAG_OMP
        thread_id = omp_get_thread_num();
#else
	thread_id = 0;
#endif
        /* J_p_X contains the particles at t1 projected from
           t0. At t1 we have data so we know the weights hence we
           compute a weighted average */

        p_hat->state[i] = 0.0;
        for(j=0;j<J;j++) {
            calc[thread_id]->to_be_sorted[j] = J_p_X[j]->proj[i]; //This sucks... gsl_sort_index requires an array to be sorted and our particles are in J_p_X[t1][j]->proj[i] so we use an helper array (calc[thread_id]->to_be_sorted)
            p_hat->state[i] += calc[thread_id]->to_be_sorted[j]*weights[j];
        }

        get_CI95(p_hat->state_95[i], calc[thread_id]->to_be_sorted, calc[thread_id]->index_sorted, weights);

    }

    /* remainder */
    if(!POP_SIZE_EQ_SUM_SV){
	double pop_size;
	for(i=0; i<N_CAC; i++) {	    
	    pop_size = gsl_spline_eval(calc[0]->spline[0][i], t, calc[0]->acc[0][i]);
	    p_hat->remainder[i] = 0.0;
	    for(j=0;j<J;j++) {
		calc[0]->to_be_sorted[j] = pop_size - sum_SV(J_p_X[j]->proj, i);
		p_hat->remainder[i] += calc[0]->to_be_sorted[j]*weights[j];
	    }

	    get_CI95(p_hat->remainder_95[i], calc[0]->to_be_sorted, calc[0]->index_sorted, weights);
	}
    }


    /* obs [N_TS] same thing as for state except that we use obs_mean()
       on p_X->obs */
#if FLAG_OMP
#pragma omp parallel for private(thread_id, j)
#endif
    for(ts=0; ts<N_TS; ts++) {
#if FLAG_OMP
        thread_id = omp_get_thread_num();
#else
	thread_id = 0;
#endif


        /* weighted average */
        p_hat->obs[ts] = 0.0;
        for(j=0;j<J;j++) {
            calc[thread_id]->to_be_sorted[j] = obs_mean(J_p_X[j]->obs[ts], p_par, p_data, calc[thread_id], ts, n, t);
            p_hat->obs[ts] += calc[thread_id]->to_be_sorted[j]*weights[j];
        }

        get_CI95(p_hat->obs_95[ts], calc[thread_id]->to_be_sorted, calc[thread_id]->index_sorted, weights);

    } /* end for on ts */

    /* drift */
    struct s_drift **drift = p_data->drift;
    int offset = 0;
    for (i=0; i< p_data->p_it_only_drift->length ; i++) {
        int ind_par_Xdrift_applied = drift[i]->ind_par_Xdrift_applied;
#if FLAG_OMP
#pragma omp parallel for private(thread_id, j) //we parallelize k and not i as in most cases there are only one single diffusion
#endif
        for (k=0; k< routers[ ind_par_Xdrift_applied ]->n_gp; k++) {
#if FLAG_OMP
	    thread_id = omp_get_thread_num();
#else
	    thread_id = 0;
#endif


            p_hat->drift[offset+k] = 0.0;

            for (j=0; j<J; j++) {
                calc[thread_id]->to_be_sorted[j] = (*(routers[ ind_par_Xdrift_applied ]->f_inv))( J_p_X[j]->proj[drift[i]->offset + k], routers[ind_par_Xdrift_applied]->min[k], routers[ind_par_Xdrift_applied]->max[k]);
                p_hat->drift[offset+k] += calc[thread_id]->to_be_sorted[j]*weights[j];
            }
            get_CI95(p_hat->drift_95[offset+k], calc[thread_id]->to_be_sorted, calc[thread_id]->index_sorted, weights);

        }

        offset += routers[ ind_par_Xdrift_applied ]->n_gp;
    } /* end for on i */
}



/**
 * Compute hat at nn
 * 2 versions are possible: all the particles have the same
 * parameters (J of J_p_par=1, is_p_par_cst=1) or the particles have
 * different parameters values (J of J_p_par =J, is_p_par_cst=0).
 */
void compute_hat_nn(struct s_X **J_p_X, struct s_par **J_p_par, struct s_data *p_data, struct s_calc **calc, struct s_hat *p_hat, int is_p_par_cst, const int n, const double t)
{
    int j, i, k, ts;
    int thread_id;
    
    struct s_router **routers = p_data->routers;

    /* par_sv */
#if FLAG_OMP
#pragma omp parallel for private(thread_id, j)
#endif
    for(i=0; i<(N_PAR_SV*N_CAC); i++) {
#if FLAG_OMP
        thread_id = omp_get_thread_num();
#else
	thread_id = 0;
#endif

        p_hat->state[i] = 0.0;
        for(j=0;j<J;j++) {
            calc[thread_id]->to_be_sorted[j] = J_p_X[j]->proj[i]; //This sucks... gsl_sort_index requires an array to be sorted and our particles are in J_p_X[t1][j]->proj[i] so we use an helper array (calc[thread_id]->to_be_sorted)
            p_hat->state[i] += calc[thread_id]->to_be_sorted[j];
        }
        p_hat->state[i] /= ((double) J);

        get_CI95(p_hat->state_95[i], calc[thread_id]->to_be_sorted, calc[thread_id]->index_sorted, NULL);
    }

    /* remainder */
    if(!POP_SIZE_EQ_SUM_SV){
	double pop_size;;
	for(i=0; i<N_CAC; i++) {	    
	    pop_size = gsl_spline_eval(calc[0]->spline[0][i],t,calc[0]->acc[0][i]);
	    p_hat->remainder[i] = 0.0;
	    for(j=0;j<J;j++) {
		calc[0]->to_be_sorted[j] = pop_size - sum_SV(J_p_X[j]->proj, i);
		p_hat->remainder[i] += calc[0]->to_be_sorted[j];
	    }
	    p_hat->remainder[i] /= ((double) J);

	    get_CI95(p_hat->remainder_95[i], calc[0]->to_be_sorted, calc[0]->index_sorted, NULL);
	}
    }

    /* obs [N_TS] same thing as for state except that we use obs_mean()
       on p_X->obs */

    if(is_p_par_cst){ //J_p_par -> J_p_par[0]

#if FLAG_OMP
#pragma omp parallel for private(thread_id, j)
#endif
	for(ts=0; ts<N_TS; ts++) {
#if FLAG_OMP
	    thread_id = omp_get_thread_num();
#else
	    thread_id = 0;
#endif


            /* empirical average */
            p_hat->obs[ts] = 0.0;
            for(j=0;j<J;j++) {
                calc[thread_id]->to_be_sorted[j] = obs_mean(J_p_X[j]->obs[ts], J_p_par[0], p_data, calc[thread_id], ts, n, t);
                p_hat->obs[ts] += calc[thread_id]->to_be_sorted[j];
            }
            p_hat->obs[ts] /= ((double) J);

            get_CI95(p_hat->obs_95[ts], calc[thread_id]->to_be_sorted, calc[thread_id]->index_sorted, NULL);

        } /* end for on ts */
    } else { //same as above, the only difference is that J_p_par -> J_p_par[j]

#if FLAG_OMP
#pragma omp parallel for private(thread_id, j)
#endif
        for(ts=0; ts<N_TS; ts++) {
#if FLAG_OMP
	    thread_id = omp_get_thread_num();
#else
	    thread_id = 0;
#endif

            /* empirical average */
            p_hat->obs[ts] = 0.0;
            for(j=0;j<J;j++) {
                calc[thread_id]->to_be_sorted[j] = obs_mean(J_p_X[j]->obs[ts], J_p_par[j], p_data, calc[thread_id], ts, n, t);
                p_hat->obs[ts] += calc[thread_id]->to_be_sorted[j];
            }
            p_hat->obs[ts] /= ((double) J);

            get_CI95(p_hat->obs_95[ts], calc[thread_id]->to_be_sorted, calc[thread_id]->index_sorted, NULL);

        } /* end for on ts */

    }


    /* drift */
    struct s_drift **drift = p_data->drift;
    int offset = 0;
    for (i=0; i< p_data->p_it_only_drift->length ; i++) {
        int ind_par_Xdrift_applied = drift[i]->ind_par_Xdrift_applied;

#if FLAG_OMP
#pragma omp parallel for private(thread_id, j) //we parallelize k and not i as in most cases there are only one single diffusion
#endif        
	for (k=0; k< routers[ ind_par_Xdrift_applied ]->n_gp; k++) {
#if FLAG_OMP
	    thread_id = omp_get_thread_num();
#else
	    thread_id = 0;
#endif


            p_hat->drift[offset+k] = 0.0;
            for(j=0; j<J; j++) {
                calc[thread_id]->to_be_sorted[j] = (*(routers[ ind_par_Xdrift_applied ]->f_inv))(J_p_X[j]->proj[drift[i]->offset + k], routers[ind_par_Xdrift_applied]->min[k], routers[ind_par_Xdrift_applied]->max[k]);
                p_hat->drift[offset+k] += calc[thread_id]->to_be_sorted[j];
            }
            p_hat->drift[offset+k] /= ((double) J);

            get_CI95(p_hat->drift_95[offset+k], calc[thread_id]->to_be_sorted, calc[thread_id]->index_sorted, NULL);
        }

        offset += routers[ ind_par_Xdrift_applied ]->n_gp;
    } /* end for on i */
}




/**
 * Plug hat into best->mean (initial conditions and drift)
 */
void plom_plug_hat(struct s_best *p_best, struct s_hat *p_hat, struct s_data *p_data){
    int i, ii,  g, p, offset;
    double nb_pop_in_group;
    double pop_size;
    double x;

    struct s_iterator *p_it_sv = p_data->p_it_par_sv;

    for(i=0; i<p_it_sv->length; i++){
	struct s_router *p_router = p_data->routers[p_it_sv->ind[i]];
	for(g=0; g< p_router->n_gp; g++){
	    if(p_best->is_estimated[p_it_sv->offset[i]+g]){
		nb_pop_in_group = 0.0;
		x = 0.0;	    	    
		for(p=0; p< p_router->p; p++){ //for all the population		
		    if(p_router->map[p] == g){ //k is in an element of the goup g
			//get pop size
			pop_size = 0.0;
			for(ii=0; ii<p_it_sv->length; ii++){		
			    pop_size += p_hat->state[ii*N_CAC + p];
			}
			if(!POP_SIZE_EQ_SUM_SV){
			    pop_size += p_hat->remainder[p];
			}

			x += (p_hat->state[i*N_CAC + p] / pop_size);
			nb_pop_in_group += 1.0;
		    }
		}
		x /= nb_pop_in_group;
		printf("%s %s %g %d\n", p_router->name, p_router->group_name[g], x, p_it_sv->offset[i]+g);
		//set the value in p_best
		gsl_vector_set(p_best->mean, p_it_sv->offset[i]+g, (*(p_router->f))(x, p_router->min[g], p_router->max[g]) );
	    }
	}
    }

    //drift
    struct s_iterator *p_it_drift = p_data->p_it_only_drift;
    offset = 0;
    for(i=0; i<p_it_drift->length; i++){
	struct s_router *p_router = p_data->routers[p_it_drift->ind[i]];
	for(g=0; g< p_router->n_gp; g++){
	    if(p_best->is_estimated[p_it_drift->offset[i]+g]){
		gsl_vector_set(p_best->mean, p_it_drift->offset[i]+g, (*(p_router->f))(p_hat->drift[offset], p_router->min[g], p_router->max[g]));
	    }
	    offset++;
	}
    }
}
