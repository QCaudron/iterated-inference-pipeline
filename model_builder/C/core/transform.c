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



double f_id(double x, double a, double b)
{
    return x;
}


double f_log(double x, double a, double b)
{
    double safe = ( x > ZERO_LOG ) ? x : ZERO_LOG;
    return log(safe);
}

double f_inv_log(double x, double a, double b)
{
    return exp(x);
}

double f_logit(double x, double a, double b)
{
    //sanatize
    double safe = ( x > ZERO_LOG ) ? x : ZERO_LOG;
    safe = (safe < ONE_LOGIT ) ? safe : ONE_LOGIT;

    return log(safe/(1.0-safe));
}


double f_inv_logit(double x, double a, double b)
{
    return (1.0/(1.0+exp(-x)));
}


double f_logit_ab(double x, double a, double b)
{
    if (a == b)
        return x; // nothing will happen in the transformed space for x, so no need to transform it
    else{
	double ratio = (x-a)/(b-x);
	if(ratio < ZERO_LOG){
	    ratio = ZERO_LOG;
	} else if(ratio > (1.0/ZERO_LOG)) {
	    ratio = 1.0/ZERO_LOG;
	}
	return log(ratio);
    }
}


double f_inv_logit_ab(double x, double a, double b)
{
    if (a == b) {
        return x ;
    } else {
        return (b*exp(x)+a)/(1.0+exp(x));
    }
}

double f_scale_id(double x)
{
    return x;
}

double f_scale_pow10(double x)
{
    return pow(10.0, x);
}


/**
 * 10^-x with x positive (thanks to a log transfo)
 */
double f_scale_pow10_pos(double x)
{
    return pow(10.0, -x);
}


/**
 * derivative of f_log
 */
double f_der_log(double x, double a, double b)
{
    return 1.0/x;
}



/**
 * derivative of f_inv_log
 */
double f_der_inv_log(double x, double a, double b)
{
    return exp(x);
}


/**
 * derivative of f_logit
 */
double f_der_logit(double x, double a, double b)
{
    return 1.0/(x-x*x);
}

/**
 * derivative of f_inv_logit
 */
double f_der_inv_logit(double x, double a, double b)
{
    return exp(-x)/pow(1.0 + exp(-x), 2.0);
}

/**
 * derivative of f_logit_ab
 */
double f_der_logit_ab(double x, double a, double b)
{
    if (a == b) {
        return x ;
    } else {
        return (b-a)/((x-a)*(b-x));
    }
}

/**
 * derivative of f_inv_logit_ab
 */
double f_der_inv_logit_ab(double x, double a, double b)
{
    if (a == b) {
        return x ;
    } else {
        return b*exp(x)/(exp(x) + 1.0) - (a + b*exp(x))*exp(x)/pow(exp(x) + 1.0, 2.0);
    }
}




/**
 *returns a multiplier that converts *duration* from u_par to u_data
 *using the following *multipliers
 */

double u_duration_par2u_data(const char *u_par, const char *u_data)
{
    if (strcmp(u_par, "D") == 0) {

        if(strcmp(u_data, "D") == 0)
            return 1.0;
        else if (strcmp(u_data, "W") == 0)
            return 1.0/7.0;
        else if (strcmp(u_data, "M") == 0)
            return 12.0/365.0;
        else if (strcmp(u_data, "Y") == 0)
            return 1.0/365.0;

    } else if (strcmp(u_par, "W") == 0) {

        if(strcmp(u_data, "D") == 0)
            return 7.0;
        else if (strcmp(u_data, "W") == 0)
            return 1.0;
        else if (strcmp(u_data, "M") == 0)
            return 84.0/365.0;
        else if (strcmp(u_data, "Y") == 0)
            return 7.0/365.0;

    } else if (strcmp(u_par, "M") == 0) {

        if(strcmp(u_data, "D") == 0)
            return 365.0/12.0;
        else if (strcmp(u_data, "W") == 0)
            return 365.0/84.0;
        else if (strcmp(u_data, "M") == 0)
            return 1.0;
        else if (strcmp(u_data, "Y") == 0)
            return 1.0/12.0;

    } else if (strcmp(u_par, "Y") == 0) {

        if(strcmp(u_data, "D") == 0)
            return 365.0;
        else if (strcmp(u_data, "W") == 0)
            return 365.0/7.0;
        else if (strcmp(u_data, "M") == 0)
            return 12.0;
        else if (strcmp(u_data, "Y") == 0)
            return 1.0;
    } else {
        print_err("invalid units");
        exit(EXIT_FAILURE);
    }

    return -1;
}


/**
 * Check if the parameter is a rate or a duration, if so return the
 * corresponding multiplier necessary for unit conversion, if not
 * return 1.0
 * @param direction: if 0 convert from user unit to data unit, if 1 convert from data unit to user unit
 */

double get_multiplier(const char *u_data, const json_t *par, int direction)
{
    json_t *unit = json_object_get(par, "unit");
    if (unit) {
        const char *u_par = json_string_value(unit);

        if ( (strcmp(u_par, "D") == 0) || (strcmp(u_par, "W") == 0) || (strcmp(u_par, "M") == 0) || (strcmp(u_par, "Y") == 0) ) {
            json_t *type = json_object_get(par, "type");
            if (type) {
                const char *mytype = json_string_value(type);

                if ( strcmp(mytype, "rate_as_duration") == 0) { // => duration

                    return (direction) ? u_duration_par2u_data(u_data, u_par) : u_duration_par2u_data(u_par, u_data);

                } else {
                    print_err("not a valid type");
                    exit(EXIT_FAILURE);
                }
            } else { //no type but a unit => rate


                return (direction) ? u_duration_par2u_data(u_par, u_data) : u_duration_par2u_data(u_data, u_par);

            }

        } else {
            print_err("error not a valid unit");
            exit(EXIT_FAILURE);
        }

    } else { //no unit => neither rate  nor duration
        return 1.0;
    }
}


/**
 * return 1 if the parameter is a duration, otherwise 0
 */
int is_duration(const json_t *par)
{
    json_t *unit = json_object_get(par, "unit");

    if (unit) {
        json_t *type = json_object_get(par, "type");
        if (type) {
            const char *mytype = json_string_value(type);
            if ( strcmp(mytype, "rate_as_duration") == 0) { // => duration
                return 1;
            }
        }
    }

    return 0;
}


/**
 * Set the transformation function of the router p_router
 * corresponding to the parameter par
 */

void set_f_trans(struct s_router *p_router, const json_t *par, const char *u_data, int is_bayesian)
{
    const char *prior_type = fast_get_json_string_from_object(par, "prior");
    json_t *transf = json_object_get(par, "transformation");
    const char *mytransf = json_string_value(transf);



    p_router->is_duration = is_duration(par);
    p_router->multiplier = get_multiplier(u_data, par, 0);

    if ( is_bayesian && (strcmp(prior_type, "uniform") == 0)) {

        p_router->f = &f_logit_ab;
        p_router->f_inv = &f_inv_logit_ab;

        if (strcmp(mytransf, "scale_pow10") == 0) {
            p_router->f_scale = &f_scale_pow10;
        } else if (strcmp(mytransf, "scale_pow10_pos") == 0) {
            p_router->f_scale = &f_scale_pow10_pos;
        } else {
            p_router->f_scale = &f_scale_id;
        }

        p_router->f_derivative = &f_der_logit_ab;
        p_router->f_inv_derivative = &f_der_inv_logit_ab;

    } else {

        if (strcmp(mytransf, "log")==0) {

            p_router->f = &f_log;
            p_router->f_inv = &f_inv_log;
            p_router->f_derivative = &f_der_log;
            p_router->f_inv_derivative = &f_der_inv_log;
            p_router->f_scale = &f_scale_id;

        } else if (strcmp(mytransf, "logit")==0) {

            p_router->f =  &f_logit;
            p_router->f_inv = &f_inv_logit;
            p_router->f_derivative = &f_der_logit;
            p_router->f_inv_derivative = &f_der_inv_logit;
            p_router->f_scale = &f_scale_id;

        } else if (strcmp(mytransf, "logit_ab")==0) {

            p_router->f =  &f_logit_ab;
            p_router->f_inv = &f_inv_logit_ab;
            p_router->f_derivative = &f_der_logit_ab;
            p_router->f_inv_derivative = &f_der_inv_logit_ab;
            p_router->f_scale = &f_scale_id;

        } else if (strcmp(mytransf, "identity")==0) {

            p_router->f = &f_id;
            p_router->f_inv = &f_id;
            p_router->f_derivative = &f_id;
            p_router->f_inv_derivative = &f_id;
            p_router->f_scale = &f_scale_id;

        } else if (strcmp(mytransf, "scale_pow10")==0) {

            p_router->f =  &f_id;
            p_router->f_inv = &f_id;
            p_router->f_derivative = &f_id;
            p_router->f_inv_derivative = &f_id;
            p_router->f_scale = &f_scale_pow10;

        } else if (strcmp(mytransf, "scale_pow10_pos")==0) {

            p_router->f =  &f_log;
            p_router->f_inv = &f_inv_log;
            p_router->f_derivative = &f_der_log;
            p_router->f_inv_derivative = &f_der_inv_log;
            p_router->f_scale = &f_scale_pow10_pos;

        } else {

            print_err("error transf != log, logit, logit_ab, scale_pow10");
            exit(EXIT_FAILURE);

        }
    }

}

void set_ab_z(struct s_router *r)
{
    double a, b;
    int g;
    //convert into data unit;
    for(g=0; g<r->n_gp; g++){

        //scale
        a = r->f_scale(r->min[g]);
        b = r->f_scale(r->max[g]);

        //to data unit
        a *= r->multiplier;
        b *= r->multiplier;

        //make rate (if needed)
        if(r->is_duration){
            a = 1.0/a;
            b = 1.0/b;
        }

        //be sure that a < b
        r->min_z[g] = GSL_MIN(a, b);
        r->max_z[g] = GSL_MAX(a, b);
    }
}


/**
 *   back transform selection of theta
 */
void back_transform_theta2par(struct s_par *p_par, const theta_t *theta, const struct s_iterator *p_it, struct s_data *p_data)
{

    int i, k;
    struct s_router **routers = p_data->routers;

    for(i=0; i<p_it->length; i++){
        struct s_router *r = routers[p_it->ind[i]];
        for(k=0; k< r->n_gp; k++) {
            p_par->natural[ p_it->ind[i] ][k] = back_transform_x(gsl_vector_get(theta, p_it->offset[i]+k), k, r);
        }
    }
}


double back_transform_x(double x, int g, struct s_router *r)
{
    double trans;
    //back transform
    trans= (*(r->f_inv))(x, r->min[g], r->max[g]);

    //scale
    trans= (*(r->f_scale))(trans);

    //convert unit
    trans *= r->multiplier;

    //convert to rate (if relevant) Note the 0.00001 to remain on the safe side if duration -> 0.0
    if(r->is_duration){
        trans = 1.0/(0.00001 + trans);
    }

    return trans;
}



/**
   Transform parameters entered by the user in an intuitive scale in
   the transformed scale (converting time unit to the unit of the
   data). 

   @param square_diag_sd boolean, squared the term on the diagonal of p_best->var

   @param[in, out] p_best the *untransformed* s_best structure that
   will be transformed
*/

void transform_theta(struct s_best *p_best, struct s_data *p_data, int square_diag_sd)
{
    int i, k;

    struct s_router **routers = p_data->routers;
    gsl_vector *x = p_best->mean;
    gsl_matrix *var = p_best->var;

    struct s_iterator *p_it = p_data->p_it_all;
    for (i=0; i<p_it->length; i++) {
        struct s_router *r = routers[ p_it->ind[i] ];
        for (k=0; k< r->n_gp; k++) {
            int offset = p_it->offset[i]+k;
	    gsl_vector_set(x, offset, (*(r->f))( gsl_vector_get(x, offset), r->min[k], r->max[k] ));

	    if(square_diag_sd && (gsl_matrix_get(var, offset, offset) > 0.0) ) {
		gsl_matrix_set(var, offset, offset, pow(gsl_matrix_get(var, offset, offset), 2));
	    }
	}
    }
}


void square_diag_sd(struct s_best *p_best, struct s_data *p_data)
{
    int i;
    
    gsl_matrix *var = p_best->var;
    struct s_iterator *p_it = p_data->p_it_all;
    for (i=0; i<p_it->nbtot; i++) {	    
	if (gsl_matrix_get(var, i, i) > 0.0) {
	    gsl_matrix_set(var, i, i, pow(gsl_matrix_get(var, i, i), 2));
	}
    }
}

