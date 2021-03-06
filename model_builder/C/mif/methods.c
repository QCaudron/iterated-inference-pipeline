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

#include "mif.h"


/**
 * if x is a prameter and y an initial condition and n the number of data points:
 *  rescale cov(x,y) by 1/sqrt(n)
 *  rescale cov(x,x) by (1/sqrt(n))*(1/sqrt(n))
 */

void rescale_covariance_mif(struct s_best *p_best, struct s_data *p_data)
{
    int i, k, ii, kk;
    struct s_router **routers = p_data->routers;
    struct s_iterator *p_it_mif = p_data->p_it_par_proc_par_obs_no_drift; //parameters fitted with MIF (as opposed to fixed lag smoothing)
    struct s_iterator *p_it_fls = p_data->p_it_par_sv_and_drift; //parameters fitted with fixed lag smoothing (fls)

    int row, col;
    double inv_n_data = 1.0/((double) p_data->nb_obs);
    double sqrt_inv_n_data = 1.0/sqrt((double) p_data->nb_obs);

    //mif, mif terms: rescale by 1/p_data->nb_obs

    for(i=0; i<p_it_mif->length; i++) {
        for(k=0; k< routers[p_it_mif->ind[i]]->n_gp; k++) {
            row = p_it_mif->offset[i]+k;

            for(ii=0; ii<p_it_mif->length; ii++) {
                for(kk=0; kk< routers[p_it_mif->ind[ii]]->n_gp; kk++) {
                    col = p_it_mif->offset[ii]+kk;

                    gsl_matrix_set(p_best->var, row, col, gsl_matrix_get(p_best->var, row, col)*inv_n_data);
                }
            }
        }
    }

    //mif, fls and fls, mif terms: rescale by 1/sqrt(p_data->nb_obs)

    for(i=0; i<p_it_mif->length; i++) {
        for(k=0; k< routers[p_it_mif->ind[i]]->n_gp; k++) {
            row = p_it_mif->offset[i]+k;

            for(ii=0; ii<p_it_fls->length; ii++) {
                for(kk=0; kk< routers[p_it_fls->ind[ii]]->n_gp; kk++) {
                    col = p_it_fls->offset[ii]+kk;

                    gsl_matrix_set(p_best->var, row, col, gsl_matrix_get(p_best->var, row, col)*sqrt_inv_n_data);
                    gsl_matrix_set(p_best->var, col, row, gsl_matrix_get(p_best->var, col, row)*sqrt_inv_n_data);

                }
            }
        }
    }
}




void ran_proposal_chol(theta_t *proposed, struct s_best *p_best, gsl_matrix *var, double sd_fac, struct s_calc *p_calc)
{
    int k;
    gsl_vector *ugaussian = (gsl_vector *) p_calc->method_specific_thread_safe_data;
    gsl_matrix *chol = (gsl_matrix *) p_calc->method_specific_shared_data;

    for(k=0; k<p_best->n_to_be_estimated; k++) {
        gsl_vector_set(ugaussian, k, gsl_ran_ugaussian(p_calc->randgsl));
    }

    gsl_blas_dtrmv(CblasLower, CblasNoTrans, CblasNonUnit, chol, ugaussian);

    gsl_vector_memcpy(proposed, p_best->mean);
    for(k=0; k<p_best->n_to_be_estimated; k++) {
        gsl_vector_set(proposed,
                       p_best->to_be_estimated[k],
                       gsl_vector_get(p_best->mean, p_best->to_be_estimated[k]) + gsl_vector_get(ugaussian, k));
    }
}

/**
 *  load theta_bart and theta_Vt for t0 and all subsequent time,
 *  we do this for time>t0 so that parameters that do not vary are
 *  printed correctly...
 */

void fill_theta_bart_and_Vt_mif(double **D_theta_bart, double **D_theta_Vt, struct s_best *p_best, struct s_data *p_data, int m)
{
    int n, i, k;
    int offset;

    struct s_iterator *p_it = p_data->p_it_all;
    struct s_router **routers = p_data->routers;

    for(n=0; n<=p_data->nb_obs; n++) {
        for(i=0; i<p_it->length; i++) {
            for(k=0; k< routers[ p_it->ind[i] ]->n_gp; k++) {
                offset = p_it->offset[i]+k;
                D_theta_bart[n][offset] = gsl_vector_get(p_best->mean, offset);
                D_theta_Vt[n][offset] = 0.0;
                if(n == 0) {
                    D_theta_Vt[n][offset] += ( pow(MIF_b*FREEZE, 2)*gsl_matrix_get(p_best->var, offset, offset) );
                }
            }
        }
    }
}


/**
 * Compute filtered mean and prediction var of particles at time
 * n. We take weighted averages with "weights" for the filtered mean
 * (in order to reduce monte-carlo variability) and use a numericaly
 * stable online algo for the variance.
 */

void mean_var_theta_theoretical_mif(double *theta_bart_n, double *theta_Vt_n, gsl_vector **J_theta, struct s_likelihood *p_like, struct s_data *p_data, struct s_best *p_best, int m, double var_fac, const enum plom_print print_opt)
{
    int i, j, k;

    struct s_iterator *p_it;
    if (print_opt & PLOM_PRINT_PRED_RES){
        p_it= p_data->p_it_all;
    } else {
        p_it= p_data->p_it_par_proc_par_obs_no_drift; //only this one is truely needed
    }

    struct s_router **routers = p_data->routers;

    double kn, M2, avg, delta; //for variance computations

    int offset;

    for(i=0; i<p_it->length; i++) {
        for(k=0; k< routers[ p_it->ind[i] ]->n_gp; k++) {
            offset = p_it->offset[i]+k;

            if(p_best->is_estimated[offset]) {

                theta_bart_n[offset]=0.0;

                kn=0.0;
                avg=0.0;
                M2=0.0;

                for(j=0 ; j<J ; j++) {

                    //variance computation
                    kn += 1.0;
                    delta = gsl_vector_get(J_theta[j], offset) - avg;
                    avg += (delta / kn);
                    M2 += delta*(gsl_vector_get(J_theta[j], offset) - avg);

                    //weighted average for filtered mean
                    theta_bart_n[offset] += p_like->weights[j]*gsl_vector_get(J_theta[j], offset);
                }
                theta_Vt_n[offset] = M2/(kn -1.0);

                if( (theta_Vt_n[offset]<0.0) || (isinf(theta_Vt_n[offset])==1) || (isnan(theta_Vt_n[offset])==1)) {
                    theta_Vt_n[offset]=0.0;
#if FLAG_WARNING
                    print_warning("error in variance computation");
#endif
                }

                /*we add theoretical variance corresponding to
                  mutation of theta to reduce Monte Carlo variability
                  AND ensure that Vt_n is > 0.0 (so that the crappy
                  Ionides formulae don't crash (even if it will even
                  with that)')*/
                //TODO check that this is valid in MVN case
                theta_Vt_n[offset] += var_fac*gsl_matrix_get(p_best->var, offset, offset);
            }
        }
    }
}


void print_mean_var_theta_theoretical_mif(FILE *p_file, double *theta_bart_n, double *theta_Vt_n, struct s_likelihood *p_like, struct s_data *p_data, int m, int time)
{
    int i, k, offset;

    struct s_iterator *p_it = p_data->p_it_all_no_theta_remainder;
    struct s_router **routers = p_data->routers;

#if FLAG_JSON
    json_t *root;
    json_t *json_print = json_array();
#endif

#if FLAG_JSON
    json_array_append_new(json_print, json_integer(m));
    json_array_append_new(json_print, json_integer(time));
#else
    fprintf(p_file,"%d,%d,", m, time);
#endif

    for(i=0; i<p_it->length; i++) {
        for(k=0; k < routers[p_it->ind[i]]->n_gp; k++) {
	    offset = p_it->offset[i]+k;
#if FLAG_JSON
	    json_array_append_new(json_print, json_real(theta_bart_n[offset]));
	    json_array_append_new(json_print, json_real(theta_Vt_n[offset]));
#else
	    fprintf(p_file, "%g,%g,", theta_bart_n[offset], theta_Vt_n[offset]);
#endif
	}
    }

    /* ess */
#if FLAG_JSON
    json_array_append_new(json_print, isnan(p_like->ess_n) ? json_null() : json_real(p_like->ess_n));
#else
    fprintf(p_file, "%g\n", p_like->ess_n);
#endif

#if FLAG_JSON
    root = json_pack("{s,s,s,o}", "flag", "mif", "msg", json_print);
    json_dumpf(root, stdout, JSON_COMPACT); printf("\n");
    fflush(stdout);
    json_decref(root);
#endif

}


void header_mean_var_theoretical_mif(FILE *p_file, struct s_data *p_data)
{
    int i, g;
    struct s_iterator *p_it = p_data->p_it_all_no_theta_remainder;
    struct s_router **routers = p_data->routers;

    fprintf(p_file,"index,time,");
    for(i=0; i<p_it->length; i++) {
        const char *name = routers[p_it->ind[i]]->name;
        for(g=0; g< routers[p_it->ind[i]]->n_gp; g++) {
            const char *group = routers[p_it->ind[i]]->group_name[g];
            fprintf(p_file, "mean:%s:%s,var:%s:%s,", name, group, name, group);
        }
    }
    fprintf(p_file,"ess\n");
}


void resample_and_mut_theta_mif(unsigned int *select, gsl_vector **J_theta, gsl_vector **J_theta_tmp, struct s_calc **calc, struct s_data *p_data, struct s_best *p_best, double sd_fac, gsl_matrix *var_fitted, int is_mvn)
{
    int i, j, k, offset;
    int thread_id;

    struct s_iterator *p_it_mif = p_data->p_it_par_proc_par_obs_no_drift; //parameters fitted with MIF (as opposed to fixed lag smoothing)
    struct s_iterator *p_it_fls = p_data->p_it_par_sv_and_drift; //parameters fitted with fixed lag smoothing (fls)
    struct s_router **routers = p_data->routers;

    //resample
    for(j=0; j<J; j++) {
        gsl_vector_memcpy(J_theta_tmp[j], J_theta[select[j]]);
    }

    if(is_mvn){

        gsl_matrix *chol = (gsl_matrix *) calc[0]->method_specific_shared_data;

        gsl_matrix_memcpy(chol, var_fitted);
        gsl_matrix_scale(chol, sd_fac*sd_fac);

        // eval decomposition
        int status = gsl_linalg_cholesky_decomp(chol);
        if(status == GSL_EDOM) {
            // error: matrix not positive
            print_err("Covariance matrix is not positive definite");
        }

        //resample and (possibly) mutate
#if FLAG_OMP
#pragma omp parallel for private(thread_id, i, k, offset)
#endif        
	for(j=0; j<J; j++) {
#if FLAG_OMP
	    thread_id = omp_get_thread_num();
#else
	    thread_id = 0;
#endif


            for(k=0; k<p_best->n_to_be_estimated; k++) {
                gsl_vector_set((gsl_vector *) calc[thread_id]->method_specific_thread_safe_data, k, gsl_ran_ugaussian(calc[thread_id]->randgsl));
            }
            gsl_blas_dtrmv(CblasLower, CblasNoTrans, CblasNonUnit, chol, (gsl_vector *) calc[thread_id]->method_specific_thread_safe_data);

            //resample and mutate
            for(k=0; k<p_best->n_to_be_estimated; k++) {
                gsl_vector_set(J_theta[j],
                               p_best->to_be_estimated[k],
                               gsl_vector_get(J_theta_tmp[j], p_best->to_be_estimated[k]) + gsl_vector_get((gsl_vector *) calc[thread_id]->method_specific_thread_safe_data, k));
            }

            //resample non fitted component of p_it_mif //TODO: optimize
            for(i=0; i<p_it_mif->length; i++) {
                for(k=0; k< routers[p_it_mif->ind[i]]->n_gp; k++) {
                    offset = p_it_mif->offset[i]+k;
                    if(!p_best->is_estimated[offset]){
                        gsl_vector_set(J_theta[j], offset, gsl_vector_get(J_theta_tmp[j], offset)); //resample only
                    }
                }
            }

            //for initial values and drift par we only resample so we cancel the effect of the mutation
            for(i=0; i<p_it_fls->length; i++) {
                for(k=0; k< routers[p_it_fls->ind[i]]->n_gp; k++) {
                    offset = p_it_fls->offset[i]+k;
                    gsl_vector_set(J_theta[j], offset, gsl_vector_get(J_theta_tmp[j], offset)); //resample only
                }
            }
        }

    } else { //diagonal sigma

#if FLAG_OMP
#pragma omp parallel for private(thread_id, i, k, offset)
#endif
        for(j=0; j<J; j++) {
#if FLAG_OMP
	    thread_id = omp_get_thread_num();
#else
	    thread_id = 0;
#endif

            //resample and mutate
            for(i=0; i<p_it_mif->length; i++) {
                for(k=0; k< routers[p_it_mif->ind[i]]->n_gp; k++) {
                    offset = p_it_mif->offset[i]+k;

                    if(p_best->is_estimated[offset]) {
                        gsl_vector_set(J_theta[j], offset, gsl_vector_get(J_theta_tmp[j], offset) + gsl_ran_gaussian(calc[thread_id]->randgsl, sd_fac*sqrt(gsl_matrix_get(p_best->var, offset, offset)))); //resample and mut
                    } else {
                        gsl_vector_set(J_theta[j], offset, gsl_vector_get(J_theta_tmp[j], offset)); //resample only
                    }
                }
            }

            //resample only
            for(i=0; i<p_it_fls->length; i++) {
                for(k=0; k< routers[p_it_fls->ind[i]]->n_gp; k++) {
                    offset = p_it_fls->offset[i]+k;
                    gsl_vector_set(J_theta[j], offset, gsl_vector_get(J_theta_tmp[j], offset)); //resample only
                }
            }
        }
    }

}


void update_fixed_lag_smoothing(struct s_best *p_best, struct s_likelihood *p_like, gsl_vector **J_theta, struct s_data *p_data)
{
    int i, j, k;
    struct s_router **routers = p_data->routers;
    double *weights = p_like->weights;
    struct s_iterator *p_it = p_data->p_it_par_sv_and_drift;

    for(i=0; i<p_it->length; i++) {
        for(k=0; k< routers[ p_it->ind[i] ]->n_gp; k++) {
            int offset = p_it->offset[i]+k;
            gsl_vector_set(p_best->mean, offset, 0.0);
            for(j=0; j<J; j++) {
                gsl_vector_set(p_best->mean, offset, gsl_vector_get(p_best->mean, offset) + gsl_vector_get(J_theta[j], offset)*weights[j]);
            }
        }
    }
}



/**
 * update theta_best in a numericaly stable way for the first
 * iterations (as suggested in Ionides et al. 2006)
 */
void update_theta_best_stable_mif(struct s_best *p_best, double **D_theta_bart, struct s_data *p_data)
{
    int i, k;
    struct s_iterator *p_it = p_data->p_it_par_proc_par_obs_no_drift;
    struct s_router **routers = p_data->routers;

    int n;
    double tmp;

    int offset;

    for(i=0; i<p_it->length; i++) {
        for(k=0; k< routers[ p_it->ind[i] ]->n_gp; k++) {
            offset = p_it->offset[i]+k;
            if(p_best->is_estimated[offset]) {
                tmp = 0.0;
                for(n=0; n<p_data->nb_obs_nonan; n++){		    
		    tmp += D_theta_bart[p_data->indn_data_nonan[n] + 1][offset];
                }

                gsl_vector_set(p_best->mean, offset, tmp / ((double) p_data->nb_obs_nonan) );
            }
        }
    }
}


/**
 * The MIF update formulae Ionides et al 2006 PNAS (doesn't work although the authors said it should)
 * theta_{m+1} = theta_m + V(t1) * sum_n{1...n} [ 1/V(t_n) (theta(t_n) - theta(t_{n-1})) ]
 * NOTE: D_theta_bart is in [N_DATA+1] so everything has to be shiftted by 1 
 */
void update_theta_best_king_mif(struct s_best *p_best, double **D_theta_bart, double **D_theta_Vt, struct s_data *p_data, int m)
{
    int i, k;
    struct s_iterator *p_it = p_data->p_it_par_proc_par_obs_no_drift;
    struct s_router **routers = p_data->routers;

    int n, nn, nnp1;
    double tmp;

    int offset;

    for(i=0; i<p_it->length; i++) {
        for(k=0; k< routers[ p_it->ind[i] ]->n_gp; k++) {
            offset = p_it->offset[i]+k;
            if(p_best->is_estimated[offset]) {
		//from initial condition (before first data point) to first data point
		nnp1 = p_data->indn_data_nonan[0]; //first entry with data
                tmp = ( (D_theta_bart[nnp1 + 1][offset]-D_theta_bart[0][offset]) / D_theta_Vt[nnp1+1][offset] ); //+1 as D_theta_bart is in N_DATA+1

		//from data point to data point
                for(n=1; n< p_data->nb_obs_nonan; n++){		    
		    nn = p_data->indn_data_nonan[n-1]; //previous entry with data
		    nnp1 = p_data->indn_data_nonan[n]; //current entry with data
		    tmp += ( (D_theta_bart[nnp1+1][offset]-D_theta_bart[nn+1][offset]) / D_theta_Vt[nnp1+1][offset] ); //+1 as D_theta_bart is in N_DATA+1
                }

                gsl_vector_set(p_best->mean,
                               offset,
                               (gsl_vector_get(p_best->mean, offset) + ((p_data->times[0] + MIF_b*MIF_b)*FREEZE*FREEZE*gsl_matrix_get(p_best->var, offset, offset)*tmp) ) );  //Cf Ionides 2006 for ((p_data->times[0] + MIF_b*MIF_b)*FREEZE*FREEZE*gsl_matrix_get(p_best->var, offset, offset)
            }
        }
    }
}


/**
 * The MIF can be used to find modes of the posterior density before
 * launching a pmcmc or kmcmc. This is controled with the --prior
 * option.
 *
 * This function multiply the loglikelihood weight of particle j at
 * observation i by 1/NbObs * logprior(theta_i^j)
 */
void patch_likelihood_prior(struct s_likelihood *p_like, struct s_best *p_best, gsl_vector **J_theta, struct s_data *p_data, int n, const int lag)
{
    int i, j, k;
    struct s_router **routers = p_data->routers;
    struct s_router *p_router;
    struct s_iterator *p_it_mif = p_data->p_it_par_proc_par_obs_no_drift; //parameters fitted with MIF (as opposed to fixed lag smoothing)
    struct s_iterator *p_it_fls = p_data->p_it_par_sv_and_drift; //parameters fitted with fixed lag smoothing (fls)

    double back_transformed; //prior are on the natural scale , so we transform the parameter into this scale...
    double p_tmp;

    int offset;

    // weights are multiplied by prior(theta_j)^(1/p_data->nb_obs) for proc parameters fitted with MIF (as opposed to fixed lag smoothing)
    for(i=0; i<p_it_mif->length; i++) {
        p_router = routers[p_it_mif->ind[i]];
        for(k=0; k< p_router->n_gp; k++) {
            offset = p_it_mif->offset[i]+k;
            if(p_best->is_estimated[offset]) {
                for(j=0; j<J; j++) {
                    back_transformed = (*(p_router->f_inv))(gsl_vector_get( J_theta[j], offset), p_router->min[k], p_router->max[k]);
                    p_tmp = (*(p_best->prior[offset]))(back_transformed, p_best->par_prior[offset][0], p_best->par_prior[offset][1]);
                    p_tmp = sanitize_likelihood(p_tmp);

                    p_like->weights[j] *= pow(p_tmp, 1.0/ ((double) p_data->nb_obs));
                }
            }
        }
    }

    // weights are multiplied by prior(theta_j)^(1/lag) for parameters fitted with fixed lag smoothing
    if(n<lag){
        for(i=0; i<p_it_fls->length; i++) {
            p_router = routers[p_it_fls->ind[i]];
            for(k=0; k< p_router->n_gp; k++) {
                offset = p_it_fls->offset[i]+k;
                if(p_best->is_estimated[offset]) {
                    for(j=0; j<J; j++) {
                        back_transformed = (*(p_router->f_inv))(gsl_vector_get( J_theta[j], offset), p_router->min[k], p_router->max[k]);
                        p_tmp = (*(p_best->prior[offset]))(back_transformed, p_best->par_prior[offset][0], p_best->par_prior[offset][1]);
                        p_tmp = sanitize_likelihood(p_tmp);

                        p_like->weights[j] *= pow(p_tmp, 1.0/lag);
                    }
                }
            }
        }
    }

}
