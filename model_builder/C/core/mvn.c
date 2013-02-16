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
 * Adapted from: Multivariate Normal density function and random
 * number generator Using GSL from Ralph dos Santos Silva
 * Copyright (C) 2006

 * multivariate normal distribution random number generator
 *
 * @param n      dimension of the random vetor
 * @param mean   vector of means of size n
 * @param var    variance matrix of dimension n x n
 * @param result output variable with a sigle random vector normal distribution generation
 */
int rmvnorm(const gsl_rng *r, const int n, const gsl_vector *mean, const gsl_matrix *var, gsl_vector *result)
{

    int k;
    gsl_matrix *work = gsl_matrix_alloc(n,n);

    gsl_matrix_memcpy(work,var);
    gsl_linalg_cholesky_decomp(work);

    for(k=0; k<n; k++)
        gsl_vector_set( result, k, gsl_ran_ugaussian(r) );

    gsl_blas_dtrmv(CblasLower, CblasNoTrans, CblasNonUnit, work, result);
    gsl_vector_add(result,mean);

    gsl_matrix_free(work);

    return 0;
}


/**
 * Adapted from: Multivariate Normal density function and random
 * number generator Using GSL from Ralph dos Santos Silva
 * Copyright (C) 2006
 *
 * multivariate normal density function
 *
 * @param n	dimension of the random vetor
 * @param mean	vector of means of size n
 * @param var	variance matrix of dimension n x n
 */

double dmvnorm(const int n, const gsl_vector *x, const gsl_vector *mean, const gsl_matrix *var)
{
    int s;
    double ax,ay;
    gsl_vector *ym, *xm;
    gsl_matrix *work = gsl_matrix_alloc(n,n),
        *winv = gsl_matrix_alloc(n,n);
    gsl_permutation *p = gsl_permutation_alloc(n);

    gsl_matrix_memcpy( work, var );
    gsl_linalg_LU_decomp( work, p, &s );
    gsl_linalg_LU_invert( work, p, winv );
    ax = gsl_linalg_LU_det( work, s );
    gsl_matrix_free( work );
    gsl_permutation_free( p );

    xm = gsl_vector_alloc(n);
    gsl_vector_memcpy( xm, x);
    gsl_vector_sub( xm, mean );
    ym = gsl_vector_alloc(n);
    gsl_blas_dsymv(CblasUpper,1.0,winv,xm,0.0,ym);
    gsl_matrix_free( winv );
    gsl_blas_ddot( xm, ym, &ay);
    gsl_vector_free(xm);
    gsl_vector_free(ym);
    ay = exp(-0.5*ay)/sqrt( pow((2*M_PI),n)*ax );

    return ay;
}


void plom_rmvnorm(gsl_vector *proposed, struct s_best *p_best, gsl_matrix *var, double sd_fac, struct s_calc *p_calc)
{
    int n = p_best->n_to_be_estimated;
    int i,k;
    gsl_matrix *work = gsl_matrix_alloc(n,n);
    double value;
    //fill work with elements of var with jump_sizes 0.0 and scale var with sd_fac^2
    for(i=0; i<n; i++) {
        for(k=0; k<n; k++) {
            value = gsl_matrix_get(var, p_best->to_be_estimated[i], p_best->to_be_estimated[k]);
            value *= sd_fac*sd_fac;
            gsl_matrix_set(work, i, k, value);
        }
    }

    // eval decomposition
    int status = gsl_linalg_cholesky_decomp(work);
    if(status == GSL_EDOM) {
        // error: matrix not positive
        print_err("Covariance matrix is not positive definite");
    }

    //working vector of size n (needed because proposed also contains
    //values whose jump_size >0.0)
    gsl_vector *xwork = gsl_vector_alloc(n);
    for(k=0; k<n; k++){
        gsl_vector_set(xwork, k, gsl_ran_ugaussian(p_calc->randgsl) );
    }

    gsl_blas_dtrmv(CblasLower, CblasNoTrans, CblasNonUnit, work, xwork);


    gsl_vector_memcpy(proposed, p_best->mean);
    //fill component of proposed with jump_size >0.0 with mean + xwork
    for (k=0; k<n; k++) {
        value = gsl_vector_get(p_best->mean, p_best->to_be_estimated[k]);
        value += gsl_vector_get(xwork, k);
        gsl_vector_set(proposed, p_best->to_be_estimated[k], value);
    }

    gsl_matrix_free(work);
    gsl_vector_free(xwork);

}


double plom_dmvnorm(struct s_best *p_best, theta_t *proposed, gsl_vector *mean, gsl_matrix *var, double sd_fac)
{
    int i, k;
    int n = p_best->n_to_be_estimated;

    int s;
    double ax,ay;
    gsl_vector *ym, *xm;
    gsl_matrix *work = gsl_matrix_alloc(n,n),
        *winv = gsl_matrix_alloc(n,n);
    gsl_permutation *p = gsl_permutation_alloc(n);

    //fill work with elements of var with jump_sizes 0.0 and scale var with sd_fac^2
    for(i=0; i<n; i++) {
        for(k=0; k<n; k++) {
            double value = gsl_matrix_get(var, p_best->to_be_estimated[i], p_best->to_be_estimated[k]);
            value *= sd_fac*sd_fac;
            gsl_matrix_set(work, i, k, value);
        }
    }

    gsl_linalg_LU_decomp( work, p, &s );
    gsl_linalg_LU_invert( work, p, winv );
    ax = gsl_linalg_LU_det( work, s );
    gsl_matrix_free( work );
    gsl_permutation_free( p );


    xm = gsl_vector_alloc(n);
    //fill xm with component of proposed - mean whose jump_size are >0.0
    for(k=0; k<n; k++) {
        double value = gsl_vector_get(proposed, p_best->to_be_estimated[k]) ;
        value -= gsl_vector_get(proposed, p_best->to_be_estimated[k]);
        gsl_vector_set(xm, k, value);
    }

    ym = gsl_vector_alloc(n);
    gsl_blas_dsymv(CblasUpper,1.0,winv,xm,0.0,ym);
    gsl_matrix_free( winv );
    gsl_blas_ddot( xm, ym, &ay);
    gsl_vector_free(xm);
    gsl_vector_free(ym);

    ay = exp(-0.5*ay)/sqrt( pow((2*M_PI),n)*ax );

    return ay;
}

/**
 * evaluate the empirical covariance matrix using a stable one-pass
 * algorithm: see
 * http://en.wikipedia.org/wiki/Algorithms%5Ffor%5Fcalculating%5Fvariance#Covariance
 */
void eval_var_emp(struct s_best *p_best, double m)
{
    int i, k;

    double *x_bar = p_best->mean_sampling;
    gsl_vector *x = p_best->mean;
    gsl_matrix *cov = p_best->var_sampling;
    double val;

    /* lower triangle and diagonal */
    for (i=0; i < p_best->length; i++) {
        for (k=0; k <= i; k++) {
	    //Cov_n = C_n/n
	    //C_n = sum_{i=1,n} (x_i-x_bar_n)(y_i-y_bar_n)
	    //C_n = C_{n-1} + (n-1)/n(x_n -x_bar_{n-1})(y_n -y_bar_{n-1})
	    //val = (n-1)/n(x_n -x_bar_{n-1})(y_n -y_bar_{n-1})
            val = ((m - 1.0) / m) * (gsl_vector_get(x, i) - x_bar[i]) * (gsl_vector_get(x, k) - x_bar[k]);

	    //C_n = C_{n-1} + val
	    //NOTE: we track the covariance and not C_n. However C_n =  covariance * n  so C_{n-1} = cov_{n-1} * (n-1) hence the fomula below
            gsl_matrix_set(cov, i, k, (((m-1.0)*gsl_matrix_get(cov, i, k)) + val ) / m);
        }
    }
    
    /* fill upper triangle */
    for (i=0; i < p_best->length; i++) {
        for (k=0; k < i; k++) {
            val = gsl_matrix_get(cov, i, k);
            gsl_matrix_set(cov, k, i, val);
        }
    }

    //update x_bar
    for (i=0; i < p_best->length; i++) {
        x_bar[i] += (gsl_vector_get(x, i) - x_bar[i]) / m;
    }

}

