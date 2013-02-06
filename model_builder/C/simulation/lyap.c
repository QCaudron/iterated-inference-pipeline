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

#include "simulation.h"

void lyapunov(struct s_calc *p_calc, struct s_par *p_par, double *y0, double t0, double t_end, double abs_tol, double rel_tol)
{
    double t = t0, t1 = t_end;

    p_calc->p_par = p_par; //pass the ref to p_par so that it is available wihtin the function to integrate

    const gsl_odeiv2_step_type * T = gsl_odeiv2_step_rkf45;
    gsl_odeiv2_control * control_lyap = gsl_odeiv2_control_y_new (abs_tol, rel_tol); /*abs and rel error (eps_abs et eps_rel) */
    gsl_odeiv2_step * step_lyap = gsl_odeiv2_step_alloc (T, N_PAR_SV*N_CAC*(N_PAR_SV*N_CAC+1));
    gsl_odeiv2_evolve * evolve_lyap = gsl_odeiv2_evolve_alloc (N_PAR_SV*N_CAC*(N_PAR_SV*N_CAC+1));
    double h = DT;

    gsl_odeiv2_system sys_lyap = {func_lyap, jac_lyap, N_PAR_SV*N_CAC*(N_PAR_SV*N_CAC+1), p_calc};

    int i, k;

    double lyap[(N_PAR_SV*N_CAC)];
    for(i=0; i<(N_PAR_SV*N_CAC); i++)
        lyap[i] = 0.0;

    double y[N_PAR_SV*N_CAC*(N_PAR_SV*N_CAC+1)];

    //we start from y0 (on the attractor)
    for(i=0; i<(N_PAR_SV*N_CAC); i++){
        y[i]=y0[i];
    }

    //linear system is initialized by identity matrix
    for(i=(N_PAR_SV*N_CAC); i<(N_PAR_SV*N_CAC*(N_PAR_SV*N_CAC+1)); i++){
        y[i]=0.0;
    }
    for(i=1;i<=(N_PAR_SV*N_CAC);i++){
        y[i*((N_PAR_SV*N_CAC)+1)-1]=1.0;
    }


    while (t < t1){
        gsl_odeiv2_evolve_apply (evolve_lyap, control_lyap, step_lyap, &sys_lyap, &t, t1, &h, y);
        gram_schmidt_normalize(y, lyap);
    }

    if(fabs(lyap[0]/t) < SEUIL_LYAP){
        printf("we keep going until: %g\n", 2*t1);
        t1 *= 2;

        while (t < t1){
            gsl_odeiv2_evolve_apply (evolve_lyap, control_lyap, step_lyap, &sys_lyap, &t, t1, &h, y);
            gram_schmidt_normalize(y, lyap);
        }
    }

    //write results
    FILE *p_lyap = sfr_fopen(SFR_PATH, GENERAL_ID, "lyap", "w", NULL, NULL);

    fprintf(p_lyap, "%g,", lyap[0]/t);
    for(k=1; k<((N_PAR_SV*N_CAC)-1); k++)
        fprintf(p_lyap, "%g,", lyap[k]/t );
    fprintf(p_lyap, "%g\n", lyap[(N_PAR_SV*N_CAC)-1]/t );

    sfr_fclose(p_lyap);

    //clean
    gsl_odeiv2_control_free(control_lyap);
    gsl_odeiv2_evolve_free(evolve_lyap);
    gsl_odeiv2_step_free(step_lyap);
}

void gram_schmidt_normalize(double *y, double *lyap)
{
    /*
      renormalization with GRAM-SCHMIDT method and Lyap exp. computation
      Lyapunov exp. are computed in log base 2 (ln(x)/ln(2)) since it then indicates the rate of loss of bits of information
    */

    int j,k,l;
    double znorm[(N_PAR_SV*N_CAC)];
    double coef[(N_PAR_SV*N_CAC)];

    /*normalize first vector*/
    znorm[0]=0.0;
    for(j=1; j<=(N_PAR_SV*N_CAC); j++){
        znorm[0] += y[(N_PAR_SV*N_CAC)*j]*y[(N_PAR_SV*N_CAC)*j]; /*column vector in the matrix*/
    }
    znorm[0]=sqrt(znorm[0]);

    for(j=1; j<=(N_PAR_SV*N_CAC); j++){
        y[(N_PAR_SV*N_CAC)*j] = y[(N_PAR_SV*N_CAC)*j]/znorm[0];
    }

    /*we generate the other orthonornal vectors (from the second to (N_PAR_SV*N_CAC) th)*/
    for(j=1; j<(N_PAR_SV*N_CAC); j++) {
            /*we generate the j-1 coefficients resulting from the inner product*/
        for(k=0; k<=(j-1); k++) { /*loop on the coefficients*/

            coef[k]=0.0;
            for(l=1; l<=(N_PAR_SV*N_CAC); l++) {
                coef[k] += y[(N_PAR_SV*N_CAC)*l+j]*y[(N_PAR_SV*N_CAC)*l+k]; /*inner product <v2,v1> */
            }
        }

        /*we generate the new vector with the j-1 coefficients */
        for(k=1; k<=(N_PAR_SV*N_CAC); k++){ /*loop on the coefficients*/
            for(l=0; l<=j-1; l++) {
                y[(N_PAR_SV*N_CAC)*k+j] = y[(N_PAR_SV*N_CAC)*k+j]-coef[l]*y[(N_PAR_SV*N_CAC)*k+l];
            }
        }
        /*compute its norm*/
        znorm[j] = 0.0;
        for(k=1; k<=(N_PAR_SV*N_CAC); k++) {
            znorm[j] += y[(N_PAR_SV*N_CAC)*k+j]*y[(N_PAR_SV*N_CAC)*k+j]; /*column vector in the matrix*/
        }
        znorm[j] = sqrt(znorm[j]);

        /*normalize*/
        for(k=1; k<=(N_PAR_SV*N_CAC); k++) {
            y[(N_PAR_SV*N_CAC)*k+j] = y[(N_PAR_SV*N_CAC)*k+j]/znorm[j]; /*column vector in the matrix*/
        }

    } /*end of for on j (construction of the orthonormal basis) */


    /*update Lyap exp. spectrum (Lyap exp are computed in log base 2)*/
    for(k=0; k<(N_PAR_SV*N_CAC); k++) {
        lyap[k] += log(znorm[k])/log(2.0);
    }
}
