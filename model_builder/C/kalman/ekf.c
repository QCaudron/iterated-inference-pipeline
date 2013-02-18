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
 *    MERCHANT ABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public
 *    License along with plom.  If not, see
 *    <http://www.gnu.org/licenses/>.
 *************************************************************************/

#include "kalman.h"

/**
 *  reset incidence-related rows and columns to 0
 */
void reset_inc_cov(gsl_matrix *Ct)
{
    int oi,oii;
    for (oi=N_PAR_SV*N_CAC; oi<N_PAR_SV*N_CAC+N_TS_INC; oi++) {
        for (oii=0; oii<N_KAL; oii++) {
            gsl_matrix_set(Ct,oi,oii,0.0);	// set row to 0
            gsl_matrix_set(Ct,oii,oi,0.0);	// set column to 0
        }
    }
}

/**
 * Brings Ct back to being symetric and semi-definite positive,
 * in case numerical instabilities made it lose these properties.
 *
 * In theory the EKF shouldn't need this,
 * and there is no right way to bring back a wrong covariance matrix to being right,
 * but this is the most natural way to go as far as I know.
 * We shouldn't need this anymore with the Square-Root Unscented Kalman Filter.
 */
void check_and_correct_Ct(gsl_matrix *Ct)
{
    char str[STR_BUFFSIZE];

    int i,j;
    int status;

    //////////////////////
    // STEP 1: SYMMETRY //
    //////////////////////
    // Ct = (Ct + Ct')/2


    gsl_matrix *temp = gsl_matrix_alloc (N_KAL, N_KAL);        // temporary matrix
    gsl_matrix_memcpy(temp, Ct);	// temp = Ct

    for(i=0; i< Ct->size1; i++){
	for(j=0; j< Ct->size2; j++){
	    gsl_matrix_set(Ct, i, j, ((gsl_matrix_get(temp, i, j) + gsl_matrix_get(temp, j, i)) / 2.0) ); 	   
	}
    }
    
    ////////////////////////
    // STEP 2: POSITIVITY //
    ////////////////////////
    // Bringing negative eigen values of Ct back to zero

    // compute the eigen values and vectors of Ct
    gsl_vector *eval = gsl_vector_alloc (N_KAL);	       // to store the eigen values    
    gsl_matrix *evec = gsl_matrix_alloc (N_KAL, N_KAL);        // to store the eigen vectors
    gsl_eigen_symmv_workspace *w = gsl_eigen_symmv_alloc (N_KAL);

    //IMPORTANT: The diagonal and lower triangular part of Ct are destroyed during the computation so we do the computation on temp
    status = gsl_eigen_symmv(temp, eval, evec, w);
    if (status) {
	sprintf(str, "error: %s\n", gsl_strerror (status));
	print_err(str);
    }
    gsl_eigen_symmv_free (w);

    gsl_matrix_set_zero(temp);

    int change_basis = 0;
    // eval = max(eval,0) and diag(eval)
    for (i=0; i<N_KAL; i++) {
        if (gsl_vector_get(eval,i) < 0.0) {
            gsl_vector_set(eval, i, 0.0); // keeps bringing negative eigen values to 0
	    change_basis = 1;
        }
        gsl_matrix_set(temp, i, i, gsl_vector_get(eval, i)); // puts the eigen values in the diagonal of temp
    }

    if(change_basis){
	gsl_matrix *temp2 = gsl_matrix_alloc (N_KAL, N_KAL);       // temporary matrix

	//////////////////
	// basis change //
	//////////////////

	// Ct = evec * temp * evec'

	// temp2 = 1.0*temp*t(evec) + 0.0*temp2
	status = gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, temp, evec, 0.0, temp2);
	if (status) {
	    sprintf(str, "error: %s\n", gsl_strerror (status));
	    print_err(str);
	}

	// Ct = 1.0*evec*temp2 + 0.0*temp2;
	status = gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, evec, temp2, 0.0, Ct);
	if (status) {
	    sprintf(str, "error: %s\n", gsl_strerror (status));
	    print_err(str);
	}

	gsl_matrix_free(temp2);
    }

    gsl_vector_free(eval);
    gsl_matrix_free(evec);
    gsl_matrix_free(temp);
}

/**
 * Computation of the EKF gain kt for observation data_t_ts and obs
 * jacobian ht, given estimate xk_t_ts and current covariance Ct
 */
void ekf_gain_computation(struct s_kalman_update *p, double xk_t_ts, double data_t_ts, gsl_matrix *Ct)
{
    char str[STR_BUFFSIZE];

    int status;
    gsl_vector *workn = gsl_vector_calloc(N_KAL); // allocating space for a temporary work vector of size N_KAL, initialized to zero

    // pred_error = double data_t_ts - xk_t_ts
    p->sc_pred_error = data_t_ts - xk_t_ts;

    // positivity and symetry could have been lost when propagating Ct
    check_and_correct_Ct(Ct);

    // sc_st = ht' * Ct * ht + sc_rt
    /*
     * here ht is a column vector to fit gsl standards,
     * rather than a row vector as it should be in theory,
     * which explains the slightly different formula we use
     */

    // workn = Ct*ht
    status = gsl_blas_dgemv(CblasNoTrans, 1.0, Ct, p->ht, 0.0, workn);
    if (status) {
        sprintf(str, "error: %s\n", gsl_strerror (status));
        print_err(str);
    }

    // sc_st = ht' * workn;
    status = gsl_blas_ddot(p->ht, workn, &(p->sc_st));
    if (status) {
        sprintf(str, "error: %s\n", gsl_strerror (status));
        print_err(str);
    }
    // sc_st = sc_st + sc_rt ;
    p->sc_st += p->sc_rt;

    // kt = Ct * ht' * sc_st^-1
    double sc_stm1 = 1.0/ (p->sc_st);
    status = gsl_blas_dgemv(CblasNoTrans, sc_stm1, Ct, p->ht, 0.0, p->kt);
    if (status) {
        sprintf(str, "error: %s\n", gsl_strerror (status));
        print_err(str);
    }

    // clear working variables:
    gsl_vector_free(workn);
}


/**
 * Update of the state vector xt and covariance matrix Ct, for a given
 * gain s_kalman_update.kt
 */
double ekf_update(struct s_kalman_update *p, gsl_matrix *Ct)
{
    char str[STR_BUFFSIZE];
    int status;

    gsl_vector *workn = gsl_vector_calloc(N_KAL); // allocating space for a temporary work vector of size N_KAL (dimension of the state vector), initialized to zero

    double like;

    //////////////////
    // state update //
    //////////////////
    // xk += kt * pred_error

    status = gsl_blas_daxpy(p->sc_pred_error, p->kt, p->xk);
    if (status) {
        sprintf(str, "error: %s\n", gsl_strerror (status));
        print_err(str);
    }

    ///////////////////////
    // covariance update //
    ///////////////////////
    // Ct = Ct - kt * ht * Ct

    // workn = Ct' * ht
    status = gsl_blas_dgemv(CblasTrans, 1.0, Ct, p->ht, 0.0, workn);
    if (status) {
        sprintf(str, "error: %s\n", gsl_strerror (status));
        print_err(str);
    }

    gsl_matrix_view Workn1 = gsl_matrix_view_vector (p->kt, p->kt->size, 1);
    gsl_matrix_view Workn2 = gsl_matrix_view_vector (workn, workn->size, 1);

    status = gsl_blas_dgemm(CblasNoTrans, CblasTrans, -1.0, &Workn2.matrix, &Workn1.matrix, 1.0, Ct); // Ct = Ct - Workn2*Workn1';
    if (status) {
        sprintf(str, "error: %s\n", gsl_strerror (status));
        print_err(str);
    }

    // positivity and symmetry could have been lost when updating Ct
    check_and_correct_Ct(Ct);

    // clear working variables:
    gsl_vector_free(workn);

    like = gsl_ran_gaussian_pdf(p->sc_pred_error, sqrt(p->sc_st));

    return sanitize_likelihood(like);
}
