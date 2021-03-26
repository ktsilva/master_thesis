/*
 *  fit.h
 *
 *  Created by Catarina Silva Fernandes 2010.
 *  All rights reserved.
 *
 */
#include <time.h>
#include <fstream>
#include <sstream>
#include <cmath>
#include <ctime>
#include <locale>
#include <iostream>
#include <stdio.h>
#include <string.h>
#include "gsl/gsl_rng.h"
#include "gsl/gsl_fit.h"
#include "gsl/gsl_multifit_nlin.h"
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include "gsl/gsl_multifit.h"


void print_state (size_t iter, gsl_multifit_fdfsolver * s) {
	printf ("iter: %lu x = %g %g %g  |f(x)| = %g\n",
			iter,
			gsl_vector_get (s->x, 0), 
			gsl_vector_get (s->x, 1),
			gsl_vector_get (s->x, 2), 
			gsl_blas_dnrm2 (s->f));
}
void print_doublePL_state (size_t iter, gsl_multifit_fdfsolver * s) {
	printf ("iter: %lu rho0= %g r_s= %g alpha= %g beta= %g"
			"|f(x)| = %g\n",
			iter,
			gsl_vector_get (s->x, 0), 
			gsl_vector_get (s->x, 1),
			gsl_vector_get (s->x, 2), 
			gsl_vector_get (s->x, 3),
			gsl_blas_dnrm2 (s->f));
	
} // SETS FUNCTION FOR DOUBLE POWER LAW

struct mini_cluster {
	int n;
	double * r;
	double * density;
	double * sigma;
};



struct data_pos {
	size_t n;
	double * r;
	double * y;
	double * sigma;
	
	/*
	 when y depends on position: vector,r
	 */
 };
struct data_int {
	size_t n;
	double * y;
	double * sigma;
	
	/*
	 when y depends on an integer: 0,1,2,3...,n-1
	 */
};

int expb_f (const gsl_vector * x, void *data,  gsl_vector * f) {
	size_t n = ((struct data_int *)data)->n;
	double *y = ((struct data_int *)data)->y;
	double *sigma = ((struct data_int*) data)->sigma;
	
	double A = gsl_vector_get (x, 0);
	double lambda = gsl_vector_get (x, 1);
	double b = gsl_vector_get (x, 2);
	
	size_t i;
	
	for (i = 0; i < n; i++)
	{
		/* Model Yi = A * exp(-lambda * i) + b */
		double t = i;
		double Yi = A * exp (-lambda * t) + b;
		gsl_vector_set (f, i, (Yi - y[i])/sigma[i]);
	}
	
	return GSL_SUCCESS;
}
int expb_df (const gsl_vector * x, void *data, gsl_matrix * J) {
	size_t n = ((struct data_int *)data)->n;
	double *sigma = ((struct data_int *) data)->sigma;
	
	double A = gsl_vector_get (x, 0);
	double lambda = gsl_vector_get (x, 1);
	
	size_t i;
	
	for (i = 0; i < n; i++)
	{
		/* Jacobian matrix J(i,j) = dfi / dxj, */
		/* where fi = (Yi - yi)/sigma[i],      */
		/*       Yi = A * exp(-lambda * i) + b  */
		/* and the xj are the parameters (A,lambda,b) */
		double t = i;
		double s = sigma[i];
		double e = exp(-lambda * t);
		gsl_matrix_set (J, i, 0, e/s); 
		gsl_matrix_set (J, i, 1, -t * A * e/s);
		gsl_matrix_set (J, i, 2, 1./s);
	}
	return GSL_SUCCESS;
}
int expb_fdf (const gsl_vector * x, void *data, gsl_vector * f, gsl_matrix * J) {
	expb_f (x, data, f);
	expb_df (x, data, J);
	return GSL_SUCCESS;
}
double * cat_fit_exp(int N_unkown_parameters, double N_data_points, double * sigma, double * y){
	const gsl_multifit_fdfsolver_type *T;
	gsl_multifit_fdfsolver *s;
	int status;
	unsigned int i, iter = 0;
	const size_t n = N_data_points;
	const size_t p = N_unkown_parameters;
	
	gsl_matrix *covar = gsl_matrix_alloc (p, p);
	
	struct data_int d = { n, y, sigma};
	double x_init[3] = { 1.0, 0.0, 0.0 };
	//gsl_vector * x = gsl_vector_alloc (p);
	gsl_vector_view x = gsl_vector_view_array (x_init, p);
	
	
	gsl_multifit_function_fdf f;
	
	f.f = &expb_f;
	f.df = &expb_df;
	f.fdf = &expb_fdf;
	f.n = n;
	f.p = p;
	f.params = &d;
	
	T = gsl_multifit_fdfsolver_lmsder;
	s = gsl_multifit_fdfsolver_alloc (T, n, p);
	
	gsl_multifit_fdfsolver_set (s, &f, &x.vector);
	
	//print_state (iter, s);
	
	do{
		iter++;
		status = gsl_multifit_fdfsolver_iterate (s);
		
		//	printf ("status = %s\n", gsl_strerror (status));
		
		//	print_state (iter, s);
		
		if (status) break;
		status = gsl_multifit_test_delta (s->dx, s->x, 1e-4, 1e-4);
		//printf ("status = %s\n", gsl_strerror (status));
		
	} while (status == GSL_CONTINUE && iter < 500);
	gsl_multifit_covar (s->J, 0.0, covar);
	
#define FIT(i) gsl_vector_get(s->x, i)
#define ERR(i) sqrt(gsl_matrix_get(covar,i,i))
	
	
	double chi = gsl_blas_dnrm2(s->f);
	double dof = n - p;
	double c = GSL_MAX_DBL(1, chi / sqrt(dof)); 
	
	printf("chisq/dof = %g, chisq= %g\n",  pow(chi, 2.0) / dof, pow(chi,2.));
	
	printf ("A      = %g +/- %g\n", FIT(0), c*ERR(0));
	printf ("lambda = %g +/- %g\n", FIT(1), c*ERR(1));
	printf ("b      = %g +/- %g\n", FIT(2), c*ERR(2));
	
	
	printf ("status = %s\n", gsl_strerror (status));
	
	double * par;
	par = new double [6];
	
	par[0] = FIT(0);
	par[1] = FIT(1);
	par[2] = FIT(2);
	par[3] = c*ERR(0);
	par[4] = c*ERR(1);
	par[5] = c*ERR(2);	
	
	gsl_multifit_fdfsolver_free (s);
	gsl_matrix_free (covar);
	
	return par;
	
}



int doublePL_f(const gsl_vector * x, void *data, gsl_vector * f) {
	size_t n = ((struct data_pos *)data)->n;
	double *r = ((struct data_pos *)data)->r; 
	double *y = ((struct data_pos *)data)->y;
	double *sigma = ((struct data_pos *)data)->sigma;
	
	double rho0 = gsl_vector_get (x, 0);
	double r_s = gsl_vector_get (x, 1);
	double alpha = gsl_vector_get (x, 2);
	double beta = gsl_vector_get (x, 3);
	
	size_t i;
	
	for (i = 0; i < n; i++)
	{
		/* Model Yi = rho0 /( pow(r/r_s,alpha) * pow(1+r/r_s,beta-alpha)) */
		double Yi = rho0 * pow((r_s + r[i])/r[i],alpha) * pow(r_s/(r_s+r[i]),beta);
		gsl_vector_set (f, i, (Yi - y[i])/sigma[i]);
	}
	return GSL_SUCCESS;
}
int doublePL_df (const gsl_vector * x, void *data,  gsl_matrix * J){
	size_t n = ((struct data_pos *)data)->n;
	double *r = ((struct data_pos *)data)->r; 
	double *sigma = ((struct data_pos *)data)->sigma;
	
	double rho0 = gsl_vector_get (x, 0);
	double r_s = gsl_vector_get (x, 1);
	double alpha = gsl_vector_get (x, 2);
	double beta = gsl_vector_get (x, 3);
	
	double F [n];
	
	size_t i;
	for (i = 0; i < n; i++)
	{
		/* Model Yi = rho0 /( pow(r/r_s,alpha) * pow(1+r/r_s,beta-alpha)) */
		F[i] = rho0 * pow((r_s + r[i])/r[i],alpha) * pow(r_s/(r_s+r[i]),beta);
	}
	
	for (i = 0; i < n; i++) {
		/* Jacobian matrix J(i,j) = dfi / dxj, */
		/* where fi = (Yi - yi)/sigma[i],      */
		/* and the xj are the parameters (rho0,r_s,alpha,beta) */
		
		gsl_matrix_set(J,i,0, (F[i]/sigma[i]) * (1./rho0)); //rho0
		gsl_matrix_set(J,i,1, (F[i]/sigma[i]) * (1./(r_s + r[i]))* (alpha + ((beta * r[i])/r_s))); //r_s
		gsl_matrix_set(J,i,2, (F[i]/sigma[i]) * log((r[i]+r_s)/r[i])); //alpha
		gsl_matrix_set(J,i,3, (F[i]/sigma[i]) * log(r_s/(r[i]+r_s))); //beta
		
	}
	return GSL_SUCCESS;
}
int doublePL_fdf (const gsl_vector * x, void *data, gsl_vector * f, gsl_matrix * J) {
	doublePL_f (x, data, f);
	doublePL_df (x, data, J);
	return GSL_SUCCESS;
}
double * cat_fit_doublePL(double * x_init, void *dados) {
	
	int status;
	int N_data_points = ((struct mini_cluster *)dados)->n;
	double *r = ((struct mini_cluster *)dados)->r; 
	double *y = ((struct mini_cluster *)dados)->density;
	double *sigma = ((struct mini_cluster *)dados)->sigma;
	
	unsigned int i, iter = 0;
	const size_t n = N_data_points;
	const size_t p = 4;
	
	struct data_pos d = { n, r, y, sigma}; /* d will have the real data */
	gsl_vector_view x = gsl_vector_view_array (x_init, p); /*initial guess for parameters x_init={rho0, r_s, alpha, beta} */
	
	gsl_multifit_function_fdf f;/* defines a general system of functions with parameters and the corresponding Jacobian matrix of derivatives */	
	
	f.f = &doublePL_f; /*stores the vector result f(x,params) in f for argument x and parameters params*/
	f.df = &doublePL_df; /*stores the n-by-p matrix result J_ij in J for argument x and parameters params */
	f.fdf = &doublePL_fdf;/*make this thing work fast: compute the function and its derivative at the same time!
						   -sets the values of the f and J as above, for arguments x and parameters params.*/
	f.n = n; /*number of components of the vector f*/
	f.p = p; /*number of components of the vectors x - independent variables*/
	f.params = &d; /*a pointer to the parameters of the function*/
	
	const gsl_multifit_fdfsolver_type *T = gsl_multifit_fdfsolver_lmsder; /*Levenberg-Marquardt algorithm (ver:lmder routine in minpack)*/
	gsl_multifit_fdfsolver *s = gsl_multifit_fdfsolver_alloc (T, n, p); /*returns a pointer to a derivative solver of type T for n observations and p parameters*/
	gsl_multifit_fdfsolver_set (s, &f, &x.vector); /*initialises solver s to use the function and derivative fdf and the initial guess x.*/
	
	//print_doublePL_state (iter, s);
	
	int stop1, stop2; 
	gsl_vector *gradt;
	gradt = gsl_vector_alloc(p);
	
	do{
		iter++;
		status = gsl_multifit_fdfsolver_iterate (s);/*perform a single iteration of the solver s*/
		/* to track the progress of the solution:
		 x-current position, f-function value at the current position,
		 dx-difference between the current position and the previous position, J-J at the current postion
		 */
		
		//printf ("status = %s\n", gsl_strerror (status));
		//print_doublePL_state (iter, s);
		//printf ("status = %s \n", gsl_strerror (status));
		
		
		if (status) break;
		
		/* Deciding when to stop */
		stop1 = gsl_multifit_test_delta(s->dx, s->x, 1e-16, 1e-16);
		gsl_multifit_gradient(s->J, s->f, gradt);/*computes the gradient g of \Phi(x) = (1/2) ||F(x)||^2 from J and f (g=J^T f)*/
		stop2 = gsl_multifit_test_gradient(gradt, 1e-16);
		
	} while (stop1 == GSL_CONTINUE || (stop2==GSL_CONTINUE) && iter < 500);
	
	gsl_matrix *covar = gsl_matrix_alloc (p, p);
	gsl_multifit_covar (s->J, 0.0, covar); /*uses J to compute COVAR of the best-fit parameters.COVAR = (J^T J)^{-1}*/
	
	/* NOTA:
	 If f_i = (Y(x, t_i) - y_i) / \sigma_i -> COVAR gives the statistical error on the best-fit parameters x
	 resulting from the gaussian errors \sigma_i on data y_i.
	 This can be verified from the relation \delta f = J \delta c and the fact that the fluctuations in f from the data y_i are normalised by \sigma_i and so satisfy <\delta f \delta f^T> = I.
	 
	 If f_i = (Y(x, t_i) - y_i) -> COVAR should be multiplied by the variance of the residuals about the best-fit
	 \sigma^2 = \sum (y_i - Y(x,t_i))^2 / (n-p) to give the variance-covariance matrix \sigma^2 C.
	 This estimates the statistical error on the best-fit parameters from the scatter of the underlying data.
	 */
	
#define FIT(i) gsl_vector_get(s->x, i) /*saving params in FIT()*/
#define ERR(i) sqrt(gsl_matrix_get(covar,i,i)) /*saving COVAR in ERR() */
	
	
	double chi = gsl_blas_dnrm2(s->f); /*computes: chi=||f|| = sqrt[sum_i[f_i^2]] */
	double dof = n - p; /*degrees of freedom*/
	double c = GSL_MAX_DBL(1, chi / sqrt(dof));/*returns the MAX of the DOUBLE precision (a,b)*/ 
	printf("chisq/dof = %g, chisq= %g\n",  pow(chi, 2.0) / dof, pow(chi,2.));
	
	printf ("rho0      = %g +/- %g ; our guess: %g \n", FIT(0), c*ERR(0), gsl_vector_get(&x.vector,0));
	printf ("r_s = %g +/- %g ; our guess: %g \n", FIT(1), c*ERR(1), gsl_vector_get(&x.vector,1));
	printf ("alpha     = %g +/- %g ; our guess: %g \n", FIT(2), c*ERR(2), gsl_vector_get(&x.vector,2));
	printf ("beta     = %g +/- %g ; our guess: %g \n", FIT(3), c*ERR(3), gsl_vector_get(&x.vector,3));
	
	
	printf ("status = %s\n", gsl_strerror (status));
	
	
	/*Variance is a measure of the variability or spread in a set of data*/
	
	
	double * par;
	par = new double [8];
	
	par[0] = FIT(0);
	par[1] = FIT(1);
	par[2] = FIT(2);
	par[3] = FIT(3);
	par[4] = c*ERR(0);
	par[5] = c*ERR(1);
	par[6] = c*ERR(2);	
	par[7] = c*ERR(3);
	
	gsl_multifit_fdfsolver_free (s);
	gsl_matrix_free (covar);
	
	return par;	
	
}

int par_temp_f(const gsl_vector * x, void *data, gsl_vector * f) {
	size_t n = ((struct data_pos *)data)->n;
	double *r = ((struct data_pos *)data)->r; 
	double *y = ((struct data_pos *)data)->y;
	double *sigma = ((struct data_pos *)data)->sigma;
	
	double T_c = y[0];
	double mu = 2.;
	double T_h = gsl_vector_get (x, 0);
	double r_s = gsl_vector_get (x, 1);
	
	size_t i;
	
	for (i = 0; i < n; i++)
	{
		/* Model Yi = rho0 /( pow(r/r_s,alpha) * pow(1+r/r_s,beta-alpha)) */
		double Yi = T_c + (T_h - T_c) * pow(r[i]/r_s, mu) / (1.+ pow(r[i]/r_s, mu));
		gsl_vector_set (f, i, (Yi - y[i])/sigma[i]);
	}
	return GSL_SUCCESS;
}
int par_temp_df (const gsl_vector * x, void *data,  gsl_matrix * J){
	size_t n = ((struct data_pos *)data)->n;
	double *r = ((struct data_pos *)data)->r; 
	double *y = ((struct data_pos *)data)->y;
	double *sigma = ((struct data_pos *)data)->sigma;
	
	double T_c = y[0];
	double mu = 2.;
	double T_h = gsl_vector_get (x, 0);
	double r_s = gsl_vector_get (x, 1);
	
	double F [n];
	
	size_t i;
	for (i = 0; i < n; i++)
	{
		/* Model Yi = rho0 /( pow(r/r_s,alpha) * pow(1+r/r_s,beta-alpha)) */
		F[i] = T_c + (T_h - T_c) * pow(r[i]/r_s, mu) / (1.+ pow(r[i]/r_s, mu));
	}
	
	for (i = 0; i < n; i++) {
		/* Jacobian matrix J(i,j) = dfi / dxj, */
		/* where fi = (Yi - yi)/sigma[i],      */
		/* and the xj are the parameters (rho0,r_s,alpha,beta) */
		
		gsl_matrix_set(J,i,0, ((pow(r[i]/r_s, mu) / (1.+ pow(r[i]/r_s, mu))) /sigma[i])); //T_h
		gsl_matrix_set(J,i,1, (1./sigma[i] * (1./pow(1.+ pow(r[i]/r_s, mu),2.)) *(T_h - T_c)* log(r[i]/r_s)*pow(r[i]/r_s, mu)));//R_c
			
	}
	return GSL_SUCCESS;
}
int par_temp_fdf (const gsl_vector * x, void *data, gsl_vector * f, gsl_matrix * J) {
	par_temp_f (x, data, f);
	par_temp_df (x, data, J);
	return GSL_SUCCESS;
}


double * cat_fit_parametrisation_temp(double * x_init, void *dados) {
	
	int status;
	int N_data_points = ((struct mini_cluster *)dados)->n;
	double *r = ((struct mini_cluster *)dados)->r; 
	double *y = ((struct mini_cluster *)dados)->density;
	double *sigma = ((struct mini_cluster *)dados)->sigma;
	
	unsigned int i, iter = 0;
	const size_t n = N_data_points;
	const size_t p = 4; // we say that beta is also an unknown parameter???
	
	struct data_pos d = { n, r, y, sigma}; /* d will have the real data */
	gsl_vector_view x = gsl_vector_view_array (x_init, p); /*initial guess for parameters x_init={rho0, r_s, alpha, beta} */
	
	gsl_multifit_function_fdf f;/* defines a general system of functions with parameters and the corresponding Jacobian matrix of derivatives */	
	
	f.f = &par_temp_f; /*stores the vector result f(x,params) in f for argument x and parameters params*/
	f.df = &par_temp_df; /*stores the n-by-p matrix result J_ij in J for argument x and parameters params */
	f.fdf = &par_temp_fdf;/*make this thing work fast: compute the function and its derivative at the same time!
									-sets the values of the f and J as above, for arguments x and parameters params.*/
	f.n = n; /*number of components of the vector f*/
	f.p = p; /*number of components of the vectors x - independent variables*/
	f.params = &d; /*a pointer to the parameters of the function*/
	
	const gsl_multifit_fdfsolver_type *T = gsl_multifit_fdfsolver_lmsder; /*Levenberg-Marquardt algorithm (ver:lmder routine in minpack)*/
	gsl_multifit_fdfsolver *s = gsl_multifit_fdfsolver_alloc (T, n, p); /*returns a pointer to a derivative solver of type T for n observations and p parameters*/
	gsl_multifit_fdfsolver_set (s, &f, &x.vector); /*initialises solver s to use the function and derivative fdf and the initial guess x.*/
	
	//print_doublePL_state (iter, s);
	
	int stop1, stop2; 
	gsl_vector *gradt;
	gradt = gsl_vector_alloc(p);
	
	do{
		iter++;
		status = gsl_multifit_fdfsolver_iterate (s);/*perform a single iteration of the solver s*/
		/* to track the progress of the solution:
		 x-current position, f-function value at the current position,
		 dx-difference between the current position and the previous position, J-J at the current postion
		 */
		
		//printf ("status = %s\n", gsl_strerror (status));
		//print_doublePL_state (iter, s);
		//printf ("status = %s \n", gsl_strerror (status));
		
		
		if (status) break;
		
		/* Deciding when to stop */
		stop1 = gsl_multifit_test_delta(s->dx, s->x, 1e-16, 1e-16);
		gsl_multifit_gradient(s->J, s->f, gradt);/*computes the gradient g of \Phi(x) = (1/2) ||F(x)||^2 from J and f (g=J^T f)*/
		stop2 = gsl_multifit_test_gradient(gradt, 1e-16);
		
	} while (stop1 == GSL_CONTINUE || (stop2==GSL_CONTINUE) && iter < 500);
	
	gsl_matrix *covar = gsl_matrix_alloc (p, p);
	gsl_multifit_covar (s->J, 0.0, covar); /*uses J to compute COVAR of the best-fit parameters.COVAR = (J^T J)^{-1}*/
	
	/* NOTA:
	 If f_i = (Y(x, t_i) - y_i) / \sigma_i -> COVAR gives the statistical error on the best-fit parameters x
	 resulting from the gaussian errors \sigma_i on data y_i.
	 This can be verified from the relation \delta f = J \delta c and the fact that the fluctuations in f from the data y_i are normalised by \sigma_i and so satisfy <\delta f \delta f^T> = I.
	 
	 If f_i = (Y(x, t_i) - y_i) -> COVAR should be multiplied by the variance of the residuals about the best-fit
	 \sigma^2 = \sum (y_i - Y(x,t_i))^2 / (n-p) to give the variance-covariance matrix \sigma^2 C.
	 This estimates the statistical error on the best-fit parameters from the scatter of the underlying data.
	 */
	
#define FIT(i) gsl_vector_get(s->x, i) /*saving params in FIT()*/
#define ERR(i) sqrt(gsl_matrix_get(covar,i,i)) /*saving COVAR in ERR() */
	
	
	double chi = gsl_blas_dnrm2(s->f); /*computes: chi=||f|| = sqrt[sum_i[f_i^2]] */
	double dof = n - p; /*degrees of freedom*/
	double c = GSL_MAX_DBL(1, chi / sqrt(dof));/*returns the MAX of the DOUBLE precision (a,b)*/ 
	printf("chisq/dof = %g, chisq= %g\n",  pow(chi, 2.0) / dof, pow(chi,2.));
	
	printf ("T_h      = %g +/- %g ; our guess: %g \n", FIT(0), c*ERR(0), gsl_vector_get(&x.vector,0));
	printf ("r_c = %g +/- %g ; our guess: %g \n", FIT(1), c*ERR(1), gsl_vector_get(&x.vector,1));
	//printf ("beta     = %g +/- %g ; our guess: %g \n", FIT(2), c*ERR(2), gsl_vector_get(&x.vector,2));
	//printf ("beta     = %g +/- %g ; our guess: %g \n", FIT(3), c*ERR(3), gsl_vector_get(&x.vector,3));
	
	
	printf ("status = %s\n", gsl_strerror (status));
	
	
	/*Variance is a measure of the variability or spread in a set of data*/
	
	
	double * par;
	par = new double [4];
	
	par[0] = FIT(0);
	par[1] = FIT(1);
	
	par[2] = c*ERR(0);
	par[3] = c*ERR(1);
	
	gsl_multifit_fdfsolver_free (s);
	gsl_matrix_free (covar);
	
	return par;	
	
}

int doublePL_fix_beta_f(const gsl_vector * x, void *data, gsl_vector * f) {
	size_t n = ((struct data_pos *)data)->n;
	double *r = ((struct data_pos *)data)->r; 
	double *y = ((struct data_pos *)data)->y;
	double *sigma = ((struct data_pos *)data)->sigma;
	
	double rho0 = gsl_vector_get (x, 0);
	double r_s = gsl_vector_get (x, 1);
	double alpha = gsl_vector_get (x, 2);
	double beta = gsl_vector_get (x, 3);
	size_t i;
	
	for (i = 0; i < n; i++)
	{
		/* Model Yi = rho0 /( pow(r/r_s,alpha) * pow(1+r/r_s,beta-alpha)) */
		double Yi = rho0 * pow((r_s + r[i])/r[i],alpha) * pow(r_s/(r_s+r[i]),beta);
		gsl_vector_set (f, i, (Yi - y[i])/sigma[i]);
	}
	return GSL_SUCCESS;
}
int doublePL_fix_beta_df (const gsl_vector * x, void *data,  gsl_matrix * J){
	size_t n = ((struct data_pos *)data)->n;
	double *r = ((struct data_pos *)data)->r; 
	double *sigma = ((struct data_pos *)data)->sigma;
	
	double rho0 = gsl_vector_get (x, 0);
	double r_s = gsl_vector_get (x, 1);
	double alpha = gsl_vector_get (x, 2);
	double beta = gsl_vector_get (x, 3);

	double F [n];
	
	size_t i;
	for (i = 0; i < n; i++)
	{
		/* Model Yi = rho0 /( pow(r/r_s,alpha) * pow(1+r/r_s,beta-alpha)) */
		F[i] = rho0 * pow((r_s + r[i])/r[i],alpha) * pow(r_s/(r_s+r[i]),beta);
	}
	
	for (i = 0; i < n; i++) {
		/* Jacobian matrix J(i,j) = dfi / dxj, */
		/* where fi = (Yi - yi)/sigma[i],      */
		/* and the xj are the parameters (rho0,r_s,alpha,beta) */
		
		gsl_matrix_set(J,i,0, (F[i]/sigma[i]) * (1./rho0)); //rho0
		gsl_matrix_set(J,i,1, (F[i]/sigma[i]) * (1./(r_s + r[i]))* (alpha + ((beta * r[i])/r_s))); //r_s
		gsl_matrix_set(J,i,2, (F[i]/sigma[i]) * log((r[i]+r_s)/r[i])); //alpha
	//	gsl_matrix_set(J,i,3, (F[i]/sigma[i]) * log(r_s/(r[i]+r_s))); //beta
		
	}
	return GSL_SUCCESS;
}
int doublePL_fix_beta_fdf (const gsl_vector * x, void *data, gsl_vector * f, gsl_matrix * J) {
	doublePL_fix_beta_f (x, data, f);
	doublePL_fix_beta_df (x, data, J);
	return GSL_SUCCESS;
}
double * cat_fit_doublePL_fix_beta(double * x_init, void *dados) {
	
	int status;
	int N_data_points = ((struct mini_cluster *)dados)->n;
	double *r = ((struct mini_cluster *)dados)->r; 
	double *y = ((struct mini_cluster *)dados)->density;
	double *sigma = ((struct mini_cluster *)dados)->sigma;
	
	unsigned int i, iter = 0;
	const size_t n = N_data_points;
	const size_t p = 4; // we say that beta is also an unknown parameter???
	
	struct data_pos d = { n, r, y, sigma}; /* d will have the real data */
	gsl_vector_view x = gsl_vector_view_array (x_init, p); /*initial guess for parameters x_init={rho0, r_s, alpha, beta} */
	
	gsl_multifit_function_fdf f;/* defines a general system of functions with parameters and the corresponding Jacobian matrix of derivatives */	
	
	f.f = &doublePL_fix_beta_f; /*stores the vector result f(x,params) in f for argument x and parameters params*/
	f.df = &doublePL_fix_beta_df; /*stores the n-by-p matrix result J_ij in J for argument x and parameters params */
	f.fdf = &doublePL_fix_beta_fdf;/*make this thing work fast: compute the function and its derivative at the same time!
									-sets the values of the f and J as above, for arguments x and parameters params.*/
	f.n = n; /*number of components of the vector f*/
	f.p = p; /*number of components of the vectors x - independent variables*/
	f.params = &d; /*a pointer to the parameters of the function*/
	
	const gsl_multifit_fdfsolver_type *T = gsl_multifit_fdfsolver_lmsder; /*Levenberg-Marquardt algorithm (ver:lmder routine in minpack)*/
	gsl_multifit_fdfsolver *s = gsl_multifit_fdfsolver_alloc (T, n, p); /*returns a pointer to a derivative solver of type T for n observations and p parameters*/
	gsl_multifit_fdfsolver_set (s, &f, &x.vector); /*initialises solver s to use the function and derivative fdf and the initial guess x.*/
	
	//print_doublePL_state (iter, s);
	
	int stop1, stop2; 
	gsl_vector *gradt;
	gradt = gsl_vector_alloc(p);
	
	do{
		iter++;
		status = gsl_multifit_fdfsolver_iterate (s);/*perform a single iteration of the solver s*/
		/* to track the progress of the solution:
		 x-current position, f-function value at the current position,
		 dx-difference between the current position and the previous position, J-J at the current postion
		 */
		
		//printf ("status = %s\n", gsl_strerror (status));
		//print_doublePL_state (iter, s);
		//printf ("status = %s \n", gsl_strerror (status));
		
		
		if (status) break;
		
		/* Deciding when to stop */
		stop1 = gsl_multifit_test_delta(s->dx, s->x, 1e-16, 1e-16);
		gsl_multifit_gradient(s->J, s->f, gradt);/*computes the gradient g of \Phi(x) = (1/2) ||F(x)||^2 from J and f (g=J^T f)*/
		stop2 = gsl_multifit_test_gradient(gradt, 1e-16);
		
	} while (stop1 == GSL_CONTINUE || (stop2==GSL_CONTINUE) && iter < 500);
	
	gsl_matrix *covar = gsl_matrix_alloc (p, p);
	gsl_multifit_covar (s->J, 0.0, covar); /*uses J to compute COVAR of the best-fit parameters.COVAR = (J^T J)^{-1}*/
	
	/* NOTA:
	 If f_i = (Y(x, t_i) - y_i) / \sigma_i -> COVAR gives the statistical error on the best-fit parameters x
	 resulting from the gaussian errors \sigma_i on data y_i.
	 This can be verified from the relation \delta f = J \delta c and the fact that the fluctuations in f from the data y_i are normalised by \sigma_i and so satisfy <\delta f \delta f^T> = I.
	 
	 If f_i = (Y(x, t_i) - y_i) -> COVAR should be multiplied by the variance of the residuals about the best-fit
	 \sigma^2 = \sum (y_i - Y(x,t_i))^2 / (n-p) to give the variance-covariance matrix \sigma^2 C.
	 This estimates the statistical error on the best-fit parameters from the scatter of the underlying data.
	 */
	
#define FIT(i) gsl_vector_get(s->x, i) /*saving params in FIT()*/
#define ERR(i) sqrt(gsl_matrix_get(covar,i,i)) /*saving COVAR in ERR() */
	
	
	double chi = gsl_blas_dnrm2(s->f); /*computes: chi=||f|| = sqrt[sum_i[f_i^2]] */
	double dof = n - p; /*degrees of freedom*/
	double c = GSL_MAX_DBL(1, chi / sqrt(dof));/*returns the MAX of the DOUBLE precision (a,b)*/ 
	printf("chisq/dof = %g, chisq= %g\n",  pow(chi, 2.0) / dof, pow(chi,2.));
	
	printf ("rho0      = %g +/- %g ; our guess: %g \n", FIT(0), c*ERR(0), gsl_vector_get(&x.vector,0));
	printf ("r_s = %g +/- %g ; our guess: %g \n", FIT(1), c*ERR(1), gsl_vector_get(&x.vector,1));
	printf ("beta     = %g +/- %g ; our guess: %g \n", FIT(2), c*ERR(2), gsl_vector_get(&x.vector,2));
	//printf ("beta     = %g +/- %g ; our guess: %g \n", FIT(3), c*ERR(3), gsl_vector_get(&x.vector,3));
	
	
	printf ("status = %s\n", gsl_strerror (status));
	
	
	/*Variance is a measure of the variability or spread in a set of data*/
	
	
	double * par;
	par = new double [8];
	
	par[0] = FIT(0);
	par[1] = FIT(1);
	par[2] = FIT(2);
	par[3] = x_init[3];
	par[4] = c*ERR(0);
	par[5] = c*ERR(1);
	par[6] = c*ERR(2);
	par[7] = 0.;
	gsl_multifit_fdfsolver_free (s);
	gsl_matrix_free (covar);
	
	return par;	
	
}

int doublePL_fix_alpha_beta_f(const gsl_vector * x, void *data, gsl_vector * f) {
	size_t n = ((struct data_pos *)data)->n;
	double *r = ((struct data_pos *)data)->r; 
	double *y = ((struct data_pos *)data)->y;
	double *sigma = ((struct data_pos *)data)->sigma;
	
	double rho0 = gsl_vector_get (x, 0);
	double r_s = gsl_vector_get (x, 1);
	double alpha = gsl_vector_get (x, 2);
	double beta = gsl_vector_get (x, 3);
	size_t i;
	
	for (i = 0; i < n; i++)
	{
		/* Model Yi = rho0 /( pow(r/r_s,alpha) * pow(1+r/r_s,beta-alpha)) */
		double Yi = rho0 * pow((r_s + r[i])/r[i],alpha) * pow(r_s/(r_s+r[i]),beta);
		gsl_vector_set (f, i, (Yi - y[i])/sigma[i]);
	}
	return GSL_SUCCESS;
}
int doublePL_fix_alpha_beta_df (const gsl_vector * x, void *data,  gsl_matrix * J){
	size_t n = ((struct data_pos *)data)->n;
	double *r = ((struct data_pos *)data)->r; 
	double *sigma = ((struct data_pos *)data)->sigma;
	
	double rho0 = gsl_vector_get (x, 0);
	double r_s = gsl_vector_get (x, 1);
	double alpha = gsl_vector_get (x, 2);
	double beta = gsl_vector_get (x, 3);
	
	double F [n];
	
	size_t i;
	for (i = 0; i < n; i++)
	{
		/* Model Yi = rho0 /( pow(r/r_s,alpha) * pow(1+r/r_s,beta-alpha)) */
		F[i] = rho0 * pow((r_s + r[i])/r[i],alpha) * pow(r_s/(r_s+r[i]),beta);
	}
	
	for (i = 0; i < n; i++) {
		/* Jacobian matrix J(i,j) = dfi / dxj, */
		/* where fi = (Yi - yi)/sigma[i],      */
		/* and the xj are the parameters (rho0,r_s,alpha,beta) */
		
		gsl_matrix_set(J,i,0, (F[i]/sigma[i]) * (1./rho0)); //rho0
		gsl_matrix_set(J,i,1, (F[i]/sigma[i]) * (1./(r_s + r[i]))* (alpha + ((beta * r[i])/r_s))); //r_s
		//gsl_matrix_set(J,i,2, (F[i]/sigma[i]) * log((r[i]+r_s)/r[i])); //alpha
		//	gsl_matrix_set(J,i,3, (F[i]/sigma[i]) * log(r_s/(r[i]+r_s))); //beta
		
	}
	return GSL_SUCCESS;
}
int doublePL_fix_alpha_beta_fdf (const gsl_vector * x, void *data, gsl_vector * f, gsl_matrix * J) {
	doublePL_fix_alpha_beta_f (x, data, f);
	doublePL_fix_alpha_beta_df (x, data, J);
	return GSL_SUCCESS;
}
double * cat_fit_doublePL_fix_alpha_beta(double * x_init, void *dados) {
	
	int status;
	int N_data_points = ((struct mini_cluster *)dados)->n;
	double *r = ((struct mini_cluster *)dados)->r; 
	double *y = ((struct mini_cluster *)dados)->density;
	double *sigma = ((struct mini_cluster *)dados)->sigma;
	
	unsigned int i, iter = 0;
	const size_t n = N_data_points;
	const size_t p = 4; // we say that beta is also an unknown parameter???
	
	struct data_pos d = { n, r, y, sigma}; /* d will have the real data */
	gsl_vector_view x = gsl_vector_view_array (x_init, p); /*initial guess for parameters x_init={rho0, r_s, alpha, beta} */
	
	gsl_multifit_function_fdf f;/* defines a general system of functions with parameters and the corresponding Jacobian matrix of derivatives */	
	
	f.f = &doublePL_fix_alpha_beta_f; /*stores the vector result f(x,params) in f for argument x and parameters params*/
	f.df = &doublePL_fix_alpha_beta_df; /*stores the n-by-p matrix result J_ij in J for argument x and parameters params */
	f.fdf = &doublePL_fix_alpha_beta_fdf;/*make this thing work fast: compute the function and its derivative at the same time!
									-sets the values of the f and J as above, for arguments x and parameters params.*/
	f.n = n; /*number of components of the vector f*/
	f.p = p; /*number of components of the vectors x - independent variables*/
	f.params = &d; /*a pointer to the parameters of the function*/
	
	const gsl_multifit_fdfsolver_type *T = gsl_multifit_fdfsolver_lmsder; /*Levenberg-Marquardt algorithm (ver:lmder routine in minpack)*/
	gsl_multifit_fdfsolver *s = gsl_multifit_fdfsolver_alloc (T, n, p); /*returns a pointer to a derivative solver of type T for n observations and p parameters*/
	gsl_multifit_fdfsolver_set (s, &f, &x.vector); /*initialises solver s to use the function and derivative fdf and the initial guess x.*/
	
	//print_doublePL_state (iter, s);
	
	int stop1, stop2; 
	gsl_vector *gradt;
	gradt = gsl_vector_alloc(p);
	
	do{
		iter++;
		status = gsl_multifit_fdfsolver_iterate (s);/*perform a single iteration of the solver s*/
		/* to track the progress of the solution:
		 x-current position, f-function value at the current position,
		 dx-difference between the current position and the previous position, J-J at the current postion
		 */
		
		//printf ("status = %s\n", gsl_strerror (status));
		//print_doublePL_state (iter, s);
		//printf ("status = %s \n", gsl_strerror (status));
		
		
		if (status) break;
		
		/* Deciding when to stop */
		stop1 = gsl_multifit_test_delta(s->dx, s->x, 1e-16, 1e-16);
		gsl_multifit_gradient(s->J, s->f, gradt);/*computes the gradient g of \Phi(x) = (1/2) ||F(x)||^2 from J and f (g=J^T f)*/
		stop2 = gsl_multifit_test_gradient(gradt, 1e-16);
		
	} while (stop1 == GSL_CONTINUE || (stop2==GSL_CONTINUE) && iter < 500);
	
	gsl_matrix *covar = gsl_matrix_alloc (p, p);
	gsl_multifit_covar (s->J, 0.0, covar); /*uses J to compute COVAR of the best-fit parameters.COVAR = (J^T J)^{-1}*/
	
	/* NOTA:
	 If f_i = (Y(x, t_i) - y_i) / \sigma_i -> COVAR gives the statistical error on the best-fit parameters x
	 resulting from the gaussian errors \sigma_i on data y_i.
	 This can be verified from the relation \delta f = J \delta c and the fact that the fluctuations in f from the data y_i are normalised by \sigma_i and so satisfy <\delta f \delta f^T> = I.
	 
	 If f_i = (Y(x, t_i) - y_i) -> COVAR should be multiplied by the variance of the residuals about the best-fit
	 \sigma^2 = \sum (y_i - Y(x,t_i))^2 / (n-p) to give the variance-covariance matrix \sigma^2 C.
	 This estimates the statistical error on the best-fit parameters from the scatter of the underlying data.
	 */
	
#define FIT(i) gsl_vector_get(s->x, i) /*saving params in FIT()*/
#define ERR(i) sqrt(gsl_matrix_get(covar,i,i)) /*saving COVAR in ERR() */
	
	
	double chi = gsl_blas_dnrm2(s->f); /*computes: chi=||f|| = sqrt[sum_i[f_i^2]] */
	double dof = n - p; /*degrees of freedom*/
	double c = GSL_MAX_DBL(1, chi / sqrt(dof));/*returns the MAX of the DOUBLE precision (a,b)*/ 
	printf("chisq/dof = %g, chisq= %g\n",  pow(chi, 2.0) / dof, pow(chi,2.));
	
	printf ("rho0      = %g +/- %g ; our guess: %g \n", FIT(0), c*ERR(0), gsl_vector_get(&x.vector,0));
	printf ("r_s = %g +/- %g ; our guess: %g \n", FIT(1), c*ERR(1), gsl_vector_get(&x.vector,1));
	printf ("alpha     = %g +/- %g ; our guess: %g \n",x_init[2],0., gsl_vector_get(&x.vector,2));
	printf ("beta     = %g +/- %g ; our guess: %g \n", x_init[3], 0., gsl_vector_get(&x.vector,3));
	
	
	printf ("status = %s\n", gsl_strerror (status));
	
	
	/*Variance is a measure of the variability or spread in a set of data*/
	
	
	double * par;
	par = new double [8];
	
	par[0] = FIT(0);
	par[1] = FIT(1);
	par[2] = x_init[2];
	par[3] = x_init[3];
	par[4] = c*ERR(0);
	par[5] = c*ERR(1);
	par[6] = 0;
	par[7] = 0.;
	gsl_multifit_fdfsolver_free (s);
	gsl_matrix_free (covar);
	
	return par;	
	
}


int doublePL_fix_alpha_f(const gsl_vector * x, void *data, gsl_vector * f) {
	size_t n = ((struct data_pos *)data)->n;
	double *r = ((struct data_pos *)data)->r; 
	double *y = ((struct data_pos *)data)->y;
	double *sigma = ((struct data_pos *)data)->sigma;
	
	double rho0 = gsl_vector_get (x, 0);
	double r_s = gsl_vector_get (x, 1);
	double alpha = gsl_vector_get (x, 2);
	double beta = gsl_vector_get (x, 3);
	size_t i;
	
	for (i = 0; i < n; i++)
	{
		/* Model Yi = rho0 /( pow(r/r_s,alpha) * pow(1+r/r_s,beta-alpha)) */
		double Yi = rho0 * pow((r_s + r[i])/r[i],alpha) * pow(r_s/(r_s+r[i]),beta);
		gsl_vector_set (f, i, (Yi - y[i])/sigma[i]);
	}
	return GSL_SUCCESS;
} // tenho de ver isto
int doublePL_fix_alpha_df (const gsl_vector * x, void *data,  gsl_matrix * J){
	size_t n = ((struct data_pos *)data)->n;
	double *r = ((struct data_pos *)data)->r; 
	double *sigma = ((struct data_pos *)data)->sigma;
	
	double rho0 = gsl_vector_get (x, 0);
	double r_s = gsl_vector_get (x, 1);
	double alpha = gsl_vector_get (x, 2);
	double beta = gsl_vector_get (x, 3);
	
	double F [n];
	
	size_t i;
	for (i = 0; i < n; i++)
	{
		/* Model Yi = rho0 /( pow(r/r_s,alpha) * pow(1+r/r_s,beta-alpha)) */
		F[i] = rho0 * pow((r_s + r[i])/r[i],alpha) * pow(r_s/(r_s+r[i]),beta);
	}
	
	for (i = 0; i < n; i++) {
		/* Jacobian matrix J(i,j) = dfi / dxj, */
		/* where fi = (Yi - yi)/sigma[i],      */
		/* and the xj are the parameters (rho0,r_s,alpha,beta) */
		
		gsl_matrix_set(J,i,0, (F[i]/sigma[i]) * (1./rho0)); //rho0
		gsl_matrix_set(J,i,1, (F[i]/sigma[i]) * (1./(r_s + r[i]))* (alpha + ((beta * r[i])/r_s))); //r_s
		//gsl_matrix_set(J,i,2, (F[i]/sigma[i]) * log((r[i]+r_s)/r[i])); //alpha
		gsl_matrix_set(J,i,2, (F[i]/sigma[i]) * log(r_s/(r[i]+r_s))); //beta
		
	}
	return GSL_SUCCESS;
}
int doublePL_fix_alpha_fdf (const gsl_vector * x, void *data, gsl_vector * f, gsl_matrix * J) {
	doublePL_fix_alpha_f (x, data, f);
	doublePL_fix_alpha_df (x, data, J);
	return GSL_SUCCESS;
}
double * cat_fit_doublePL_fix_alpha(double * x_init, void *dados) {
	
	int status;
	int N_data_points = ((struct mini_cluster *)dados)->n;
	double *r = ((struct mini_cluster *)dados)->r; 
	double *y = ((struct mini_cluster *)dados)->density;
	double *sigma = ((struct mini_cluster *)dados)->sigma;
	
	unsigned int i, iter = 0;
	const size_t n = N_data_points;
	const size_t p = 4;
	
	struct data_pos d = { n, r, y, sigma}; /* d will have the real data */
	gsl_vector_view x = gsl_vector_view_array (x_init, p); /*initial guess for parameters x_init={rho0, r_s, alpha, beta} */
	
	gsl_multifit_function_fdf f;/* defines a general system of functions with parameters and the corresponding Jacobian matrix of derivatives */	
	
	f.f = &doublePL_fix_alpha_f; /*stores the vector result f(x,params) in f for argument x and parameters params*/
	f.df = &doublePL_fix_alpha_df; /*stores the n-by-p matrix result J_ij in J for argument x and parameters params */
	f.fdf = &doublePL_fix_alpha_fdf;/*make this thing work fast: compute the function and its derivative at the same time!
									 -sets the values of the f and J as above, for arguments x and parameters params.*/
	f.n = n; /*number of components of the vector f*/
	f.p = p; /*number of components of the vectors x - independent variables*/
	f.params = &d; /*a pointer to the parameters of the function*/
	
	const gsl_multifit_fdfsolver_type *T = gsl_multifit_fdfsolver_lmsder; /*Levenberg-Marquardt algorithm (ver:lmder routine in minpack)*/
	gsl_multifit_fdfsolver *s = gsl_multifit_fdfsolver_alloc (T, n, p); /*returns a pointer to a derivative solver of type T for n observations and p parameters*/
	gsl_multifit_fdfsolver_set (s, &f, &x.vector); /*initialises solver s to use the function and derivative fdf and the initial guess x.*/
	
	//print_doublePL_state (iter, s);
	
	int stop1, stop2; 
	gsl_vector *gradt;
	gradt = gsl_vector_alloc(p);
	
	do{
		iter++;
		status = gsl_multifit_fdfsolver_iterate (s);/*perform a single iteration of the solver s*/
		/* to track the progress of the solution:
		 x-current position, f-function value at the current position,
		 dx-difference between the current position and the previous position, J-J at the current postion
		 */
		
		//printf ("status = %s\n", gsl_strerror (status));
		//print_doublePL_state (iter, s);
		//printf ("status = %s \n", gsl_strerror (status));
		
		
		if (status) break;
		
		/* Deciding when to stop */
		stop1 = gsl_multifit_test_delta(s->dx, s->x, 1e-16, 1e-16);
		gsl_multifit_gradient(s->J, s->f, gradt);/*computes the gradient g of \Phi(x) = (1/2) ||F(x)||^2 from J and f (g=J^T f)*/
		stop2 = gsl_multifit_test_gradient(gradt, 1e-16);
		
	} while (stop1 == GSL_CONTINUE || (stop2==GSL_CONTINUE) && iter < 500);
	
	gsl_matrix *covar = gsl_matrix_alloc (p, p);
	gsl_multifit_covar (s->J, 0.0, covar); /*uses J to compute COVAR of the best-fit parameters.COVAR = (J^T J)^{-1}*/
	
	/* NOTA:
	 If f_i = (Y(x, t_i) - y_i) / \sigma_i -> COVAR gives the statistical error on the best-fit parameters x
	 resulting from the gaussian errors \sigma_i on data y_i.
	 This can be verified from the relation \delta f = J \delta c and the fact that the fluctuations in f from the data y_i are normalised by \sigma_i and so satisfy <\delta f \delta f^T> = I.
	 
	 If f_i = (Y(x, t_i) - y_i) -> COVAR should be multiplied by the variance of the residuals about the best-fit
	 \sigma^2 = \sum (y_i - Y(x,t_i))^2 / (n-p) to give the variance-covariance matrix \sigma^2 C.
	 This estimates the statistical error on the best-fit parameters from the scatter of the underlying data.
	 */
	
#define FIT(i) gsl_vector_get(s->x, i) /*saving params in FIT()*/
#define ERR(i) sqrt(gsl_matrix_get(covar,i,i)) /*saving COVAR in ERR() */
	
	
	double chi = gsl_blas_dnrm2(s->f); /*computes: chi=||f|| = sqrt[sum_i[f_i^2]] */
	double dof = n - p; /*degrees of freedom*/
	double c = GSL_MAX_DBL(1, chi / sqrt(dof));/*returns the MAX of the DOUBLE precision (a,b)*/ 
	printf("chisq/dof = %g, chisq= %g\n",  pow(chi, 2.0) / dof, pow(chi,2.));
	
	printf ("rho0      = %g +/- %g ; our guess: %g \n", FIT(0), c*ERR(0), gsl_vector_get(&x.vector,0));
	printf ("r_s = %g +/- %g ; our guess: %g \n", FIT(1), c*ERR(1), gsl_vector_get(&x.vector,1));
	printf ("alpha     = %g +/- %g ; our guess: %g \n",  x_init[2], 0., gsl_vector_get(&x.vector,2));
	printf ("beta     = %g +/- %g ; our guess: %g \n", FIT(2), c*ERR(2), gsl_vector_get(&x.vector,3));
	
	
	printf ("status = %s\n", gsl_strerror (status));
	
	
	/*Variance is a measure of the variability or spread in a set of data*/
	
	
	double * par;
	par = new double [8];
	
	par[0] = FIT(0);
	par[1] = FIT(1);
	par[2] = x_init[2];
	par[3] = FIT(2);
	par[4] = c*ERR(0);
	par[5] = c*ERR(1);
	par[6] = 0.;
	par[7] = c*ERR(2);
	
	gsl_multifit_fdfsolver_free (s);
	gsl_matrix_free (covar);
	
	return par;	
	
}

int beta_profile_f(const gsl_vector * x, void *data, gsl_vector * f) { //alterar
	size_t n = ((struct data_pos *)data)->n;
	double *r = ((struct data_pos *)data)->r; 
	double *y = ((struct data_pos *)data)->y;
	double *sigma = ((struct data_pos *)data)->sigma;
	
	double rho0 = gsl_vector_get (x, 0);
	double r_c = gsl_vector_get (x, 1);
	double beta = gsl_vector_get (x, 2);
	
	size_t i;
	
	for (i = 0; i < n; i++)
	{
		/* Model Yi = rho0 /( pow(r/r_s,alpha) * pow(1+r/r_s,beta-alpha)) */
		double Yi = rho0 * pow(1. + pow((r[i]/r_c),2.),-(3./2.)*beta);
		gsl_vector_set (f, i, (Yi - y[i])/sigma[i]);
	}
	return GSL_SUCCESS;
}
int beta_profile_df (const gsl_vector * x, void *data,  gsl_matrix * J){ //ALTERAR
	size_t n = ((struct data_pos *)data)->n;
	double *r = ((struct data_pos *)data)->r; 
	double *sigma = ((struct data_pos *)data)->sigma;
	
	double rho0 = gsl_vector_get (x, 0);
	double r_c = gsl_vector_get (x, 1);
	double beta = gsl_vector_get (x, 2);
	
	double F [n];
	
	size_t i;
	for (i = 0; i < n; i++)
	{
		/* Model Yi = rho0 * pow(1. + (r[i]/r_s),-(3./2.)*beta) */
		F[i] = rho0 * pow(1. + pow(r[i]/r_c,2.),-(3./2.)*beta);
		
	}
	
	for (i = 0; i < n; i++) {
		/* Jacobian matrix J(i,j) = dfi / dxj, */
		/* where fi = (Yi - yi)/sigma[i],      */
		/* and the xj are the parameters (rho0,r_s,alpha,beta) */
		
		gsl_matrix_set(J,i,0, (F[i]/sigma[i]) * (1./rho0)); //rho0
		//
		//    (-(3/2)*beta * (1+(r/r_s)ˆ2)ˆ(-(3/2)*beta -1) * 2*(r/r_s)*(-r/r_sˆ2) )
		//
		gsl_matrix_set(J,i,1, (F[i]/sigma[i]) * 3.*beta*(pow(r[i],2.)/pow(r_c,3.))* pow( 1.+ pow(r[i],r_c) ,-1.)); //r_s
		gsl_matrix_set(J,i,2, (F[i]/sigma[i]) * (-3./2.) * log(1. + pow(r[i]/r_c,2.))); //beta
		
	}
	return GSL_SUCCESS;
}
int beta_profile_fdf (const gsl_vector * x, void *data, gsl_vector * f, gsl_matrix * J) {
	beta_profile_f (x, data, f);
	beta_profile_df (x, data, J);
	return GSL_SUCCESS;
}
double * cat_fit_beta_profile(double *x_init, void *dados) {
	
	int status;
	int N_data_points = ((struct mini_cluster *)dados)->n;
	double *r = ((struct mini_cluster *)dados)->r; 
	double *y = ((struct mini_cluster *)dados)->density;
	double *sigma = ((struct mini_cluster *)dados)->sigma;
	
	unsigned int i, iter = 0;
	const size_t n = N_data_points;
	const size_t p = 3;
	
	struct data_pos d = { n, r, y, sigma}; //d will have the real data 
	gsl_vector_view x = gsl_vector_view_array (x_init, p); //initial guess for parameters x_init={rho0, r_s, alpha, beta} 
	
	gsl_multifit_function_fdf f;//defines a general system of functions with parameters and the corresponding Jacobian matrix of derivatives 
	
	f.f = &beta_profile_f; //stores the vector result f(x,params) in f for argument x and parameters params
	f.df = &beta_profile_df; //stores the n-by-p matrix result J_ij in J for argument x and parameters params 
	f.fdf = &beta_profile_fdf; //make this thing work fast: compute the function and its derivative at the same time!
	//-sets the values of the f and J as above, for arguments x and parameters params.
	f.n = n; //number of components of the vector f
	f.p = p; //number of components of the vectors x - independent variables
	f.params = &d; //a pointer to the parameters of the function
	
	const gsl_multifit_fdfsolver_type *T = gsl_multifit_fdfsolver_lmsder; //Levenberg-Marquardt algorithm (ver:lmder routine in minpack)
	gsl_multifit_fdfsolver *s = gsl_multifit_fdfsolver_alloc (T, n, p); //returns a pointer to a derivative solver of type T for n observations and p parameters
	
	gsl_multifit_fdfsolver_set (s, &f, &x.vector); //initialises solver s to use the function and derivative fdf and the initial guess x.
	
	//print_doublePL_state (iter, s);
	
	int stop1, stop2; 
	gsl_vector *gradt;
	gradt = gsl_vector_alloc(p);
	
	do{
		iter++;
		status = gsl_multifit_fdfsolver_iterate (s);//perform a single iteration of the solver s
		// to track the progress of the solution: x-current position, f-function value at the current position, dx-difference between the current position and the previous position, J-J at the current postion
		
		//printf ("status = %s\n", gsl_strerror (status));
		//print_doublePL_state (iter, s);
		//printf ("status = %s \n", gsl_strerror (status));
		
		
		if (status) break;
		
		// Deciding when to stop 
		stop1 = gsl_multifit_test_delta(s->dx, s->x, 1.e-4, 1.e-4);
		gsl_multifit_gradient(s->J, s->f, gradt);//computes the gradient g of \Phi(x) = (1/2) ||F(x)||^2 from J and f (g=J^T f)
		stop2 = gsl_multifit_test_gradient(gradt, 1.e-4);
		
	} while (stop1 == GSL_CONTINUE || (stop2==GSL_CONTINUE) && iter < 500);
	
	gsl_matrix *covar = gsl_matrix_alloc (p, p);
	gsl_multifit_covar (s->J, 0.0, covar); //uses J to compute COVAR of the best-fit parameters.COVAR = (J^T J)^{-1}
	
	/* NOTA:
	 If f_i = (Y(x, t_i) - y_i) / \sigma_i -> COVAR gives the statistical error on the best-fit parameters x
	 resulting from the gaussian errors \sigma_i on data y_i.
	 This can be verified from the relation \delta f = J \delta c and the fact that the fluctuations in f from the data y_i are normalised by \sigma_i and so satisfy <\delta f \delta f^T> = I.
	 
	 If f_i = (Y(x, t_i) - y_i) -> COVAR should be multiplied by the variance of the residuals about the best-fit
	 \sigma^2 = \sum (y_i - Y(x,t_i))^2 / (n-p) to give the variance-covariance matrix \sigma^2 C.
	 This estimates the statistical error on the best-fit parameters from the scatter of the underlying data.
	 */
	
#define FIT(i) gsl_vector_get(s->x, i) //saving params in FIT()
#define ERR(i) sqrt(gsl_matrix_get(covar,i,i)) //saving COVAR in ERR() 
	
	
	double chi = gsl_blas_dnrm2(s->f); //computes: chi=||f|| = sqrt[sum_i[f_i^2]] 
	double dof = n - p; //degrees of freedom
	double c = GSL_MAX_DBL(1, chi / sqrt(dof));//returns the MAX of the DOUBLE precision (a,b)*
	printf("chisq/dof = %g, chisq= %g\n",  pow(chi, 2.0) / dof, pow(chi,2.));
	
	printf ("rho0      = %g +/- %g ; our guess: %g \n", FIT(0), c*ERR(0), gsl_vector_get(&x.vector,0));
	printf ("r_s = %g +/- %g ; our guess: %g \n", FIT(1), c*ERR(1), gsl_vector_get(&x.vector,1));
	printf ("beta     = %g +/- %g ; our guess: %g \n", FIT(2), c*ERR(2), gsl_vector_get(&x.vector,2));
	
	
	printf ("status = %s\n", gsl_strerror(status));
	
	double * par;
	par = new double [6];
	
	par[0] = FIT(0);
	par[1] = FIT(1);
	par[2] = FIT(2);
	par[3] = c*ERR(0);
	par[4] = c*ERR(1);
	par[5] = c*ERR(2);	
	
	gsl_multifit_fdfsolver_free (s);
	gsl_matrix_free (covar);
	
	return par;	
	
	
}

int sersic_profile_f(const gsl_vector * x, void *data, gsl_vector * f) { //alterar - atencao p&a
	size_t n = ((struct data_pos *)data)->n;
	double *r = ((struct data_pos *)data)->r; 
	double *y = ((struct data_pos *)data)->y;
	double *sigma = ((struct data_pos *)data)->sigma;
	
	double rho0 = gsl_vector_get (x, 0);
	double r_s = gsl_vector_get (x, 1);
	double v = gsl_vector_get (x, 2);
	//double p = gsl_vector_get (x,3);
	
	double pp = .5 * (1.- 0.6097* v + 0.05563*pow(v,2.));
	double aa = r_s*pow(2,(1./v));
	size_t i;
	
	for (i = 0; i < n; i++)
	{	
		// Model Yi = rho0 *(r/a')^p' * exp(-(r/a')^v) 
		double Yi = rho0 * pow(r[i]/aa,pp) * exp (-pow(r[i]/aa,v));
		//double Yi = rho0 * pow(r[i]/a,-p) * exp (-pow(r[i]/a,v));
		gsl_vector_set (f, i, (Yi - y[i])/sigma[i]);
	}
	return GSL_SUCCESS;
}
int sersic_profile_df (const gsl_vector * x, void *data,  gsl_matrix * J){ //ALTERAR
	size_t n = ((struct data_pos *)data)->n;
	double *r = ((struct data_pos *)data)->r; 
	double *sigma = ((struct data_pos *)data)->sigma;
	
	double rho0 = gsl_vector_get (x, 0);
	double r_s = gsl_vector_get (x, 1);
	double v = gsl_vector_get (x, 2);
	//double p = gsl_vector_get (x, 3);
	double a= - 0.6097; double b = 0.05563;
	double pp = .5 * (1.+a* v + b*pow(v,2.));
	double aa = r_s*pow(2,(1./v));
	
	double F [n];
	
	size_t i;
	for (i = 0; i < n; i++)
	{
		//Model Yi = rho0 *(r/a')^p' * exp(-(r/a')^v) 
		//F[i] = rho0 * pow(r[i]/a,-p) * exp (-pow(r[i]/a,v));
		F[i] = rho0 * pow(r[i]/aa,pp) * exp (-pow(r[i]/aa,v));
	}
	
	for (i = 0; i < n; i++) {
		/* Jacobian matrix J(i,j) = dfi / dxj, */
		/* where fi = (Yi - yi)/sigma[i],      */
		/* and the xj are the parameters (rho0,r_s,alpha,beta) */
		
		gsl_matrix_set(J,i,0, (F[i]/sigma[i]) * (1./rho0)); //rho0
		gsl_matrix_set(J,i,1, (F[i]/sigma[i]) * (pp + pow(r[i]/aa,v))/r_s); //r_s
		gsl_matrix_set(J,i,2, (F[i]/sigma[i]) * .5 * (log(r[i]/r_s)*(-a-2.*b*pow(v,2.))+
													  log(2)*(((-a-2.*b*v)/v) - (pp/pow(v,2.)))+
													  pow(r[i]/r_s,v)*log(r[i]/r_s))); //v
	
	//	gsl_matrix_set(J,i,0, (F[i]/sigma[i]) * (1./rho0)); //rho0
	//	gsl_matrix_set(J,i,1, (F[i]/(sigma[i]*a))* (p -v*pow(r[i]/a,v)));
	//	gsl_matrix_set(J,i,3, (F[i]/sigma[i]) * log(a/r[i]) * pow(-r[i]/a,v));
	//	gsl_matrix_set(J,i,2, (F[i]/sigma[i]) * log(a/r[i]));
		
		
		
	}
	return GSL_SUCCESS;
}
int sersic_profile_fdf (const gsl_vector * x, void *data, gsl_vector * f, gsl_matrix * J) {
	sersic_profile_f (x, data, f);
	sersic_profile_df (x, data, J);
	return GSL_SUCCESS;
}
double * cat_fit_sersic_profile(double *x_init, void *dados) {
	
	int status;
	int N_data_points = ((struct mini_cluster *)dados)->n;
	double *r = ((struct mini_cluster *)dados)->r; 
	double *y = ((struct mini_cluster *)dados)->density;
	double *sigma = ((struct mini_cluster *)dados)->sigma;
	
	unsigned int i, iter = 0;
	const size_t n = N_data_points;
	const size_t p = 3;
	
	struct data_pos d = { n, r, y, sigma}; //d will have the real data 
	gsl_vector_view x = gsl_vector_view_array (x_init, p); //initial guess for parameters x_init={rho0, r_s, alpha, beta} 
	
	gsl_multifit_function_fdf f;//defines a general system of functions with parameters and the corresponding Jacobian matrix of derivatives 
	
	f.f = &sersic_profile_f; //stores the vector result f(x,params) in f for argument x and parameters params
	f.df = &sersic_profile_df; //stores the n-by-p matrix result J_ij in J for argument x and parameters params 
	f.fdf = &sersic_profile_fdf; //make this thing work fast: compute the function and its derivative at the same time!
	//-sets the values of the f and J as above, for arguments x and parameters params.
	f.n = n; //number of components of the vector f
	f.p = p; //number of components of the vectors x - independent variables
	f.params = &d; //a pointer to the parameters of the function
	
	const gsl_multifit_fdfsolver_type *T = gsl_multifit_fdfsolver_lmsder; //Levenberg-Marquardt algorithm (ver:lmder routine in minpack)
	gsl_multifit_fdfsolver *s = gsl_multifit_fdfsolver_alloc (T, n, p); //returns a pointer to a derivative solver of type T for n observations and p parameters
	
	gsl_multifit_fdfsolver_set (s, &f, &x.vector); //initialises solver s to use the function and derivative fdf and the initial guess x.
	
	//print_doublePL_state (iter, s);
	
	int stop1, stop2; 
	gsl_vector *gradt;
	gradt = gsl_vector_alloc(p);
	
	do{
		iter++;
		status = gsl_multifit_fdfsolver_iterate (s);//perform a single iteration of the solver s
		// to track the progress of the solution: x-current position, f-function value at the current position, dx-difference between the current position and the previous position, J-J at the current postion
		
		//printf ("status = %s\n", gsl_strerror (status));
		//print_doublePL_state (iter, s);
		//printf ("status = %s \n", gsl_strerror (status));
		
		
		if (status) break;
		
		// Deciding when to stop 
		stop1 = gsl_multifit_test_delta(s->dx, s->x, 1.e-4, 1.e-4);
		gsl_multifit_gradient(s->J, s->f, gradt);//computes the gradient g of \Phi(x) = (1/2) ||F(x)||^2 from J and f (g=J^T f)
		stop2 = gsl_multifit_test_gradient(gradt, 1.e-4);
		
	} while (stop1 == GSL_CONTINUE || (stop2==GSL_CONTINUE) && iter < 500);
	
	gsl_matrix *covar = gsl_matrix_alloc (p, p);
	gsl_multifit_covar (s->J, 0.0, covar); //uses J to compute COVAR of the best-fit parameters.COVAR = (J^T J)^{-1}
	
	/* NOTA:
	 If f_i = (Y(x, t_i) - y_i) / \sigma_i -> COVAR gives the statistical error on the best-fit parameters x
	 resulting from the gaussian errors \sigma_i on data y_i.
	 This can be verified from the relation \delta f = J \delta c and the fact that the fluctuations in f from the data y_i are normalised by \sigma_i and so satisfy <\delta f \delta f^T> = I.
	 
	 If f_i = (Y(x, t_i) - y_i) -> COVAR should be multiplied by the variance of the residuals about the best-fit
	 \sigma^2 = \sum (y_i - Y(x,t_i))^2 / (n-p) to give the variance-covariance matrix \sigma^2 C.
	 This estimates the statistical error on the best-fit parameters from the scatter of the underlying data.
	 */
	
#define FIT(i) gsl_vector_get(s->x, i) //saving params in FIT()
#define ERR(i) sqrt(gsl_matrix_get(covar,i,i)) //saving COVAR in ERR() 
	
	
	double chi = gsl_blas_dnrm2(s->f); //computes: chi=||f|| = sqrt[sum_i[f_i^2]] 
	double dof = n - p; //degrees of freedom
	double c = GSL_MAX_DBL(1, chi / sqrt(dof));//returns the MAX of the DOUBLE precision (a,b)*
	printf("chisq/dof = %g\n",  pow(chi, 2.0) / dof);
	
	printf ("rho0      = %g +/- %g ; our guess: %g \n", FIT(0), c*ERR(0), gsl_vector_get(&x.vector,0));
	printf ("a = %g +/- %g ; our guess: %g \n", FIT(1), c*ERR(1), gsl_vector_get(&x.vector,1));
	printf ("v     = %g +/- %g ; our guess: %g \n", FIT(2), c*ERR(2), gsl_vector_get(&x.vector,2));
	//printf ("p    = %g +/- %g ; our guess: %g \n", FIT(3), c*ERR(3), gsl_vector_get(&x.vector,3));
	
	
	printf ("status = %s\n", gsl_strerror(status));
	
	double * par;
	par = new double [6];
	
	par[0] = FIT(0);
	par[1] = FIT(1);
	par[2] = FIT(2);
	//par[3] = FIT(3);
	par[3] = c*ERR(0);
	par[4] = c*ERR(1);
	par[5] = c*ERR(2);	
	//par[7] = c*ERR(3);	
	
	gsl_multifit_fdfsolver_free (s);
	gsl_matrix_free (covar);
	
	return par;	
	
	
}

int sersic_profile_fix_v_f(const gsl_vector * x, void *data, gsl_vector * f) { //alterar - atencao p&a
	size_t n = ((struct data_pos *)data)->n;
	double *r = ((struct data_pos *)data)->r; 
	double *y = ((struct data_pos *)data)->y;
	double *sigma = ((struct data_pos *)data)->sigma;
	
	double rho0 = gsl_vector_get (x, 0);
	double r_s = gsl_vector_get (x, 1);
	double v = gsl_vector_get (x, 2);
	//double p = gsl_vector_get (x,3);
	
	double pp = .5 * (1.- 0.6097* v + 0.05563*pow(v,2.));
	double aa = r_s*pow(2,(1./v));
	size_t i;
	
	for (i = 0; i < n; i++)
	{	
		// Model Yi = rho0 *(r/a')^p' * exp(-(r/a')^v) 
		double Yi = rho0 * pow(r[i]/aa,pp) * exp (-pow(r[i]/aa,v));
		//double Yi = rho0 * pow(r[i]/a,-p) * exp (-pow(r[i]/a,v));
		gsl_vector_set (f, i, (Yi - y[i])/sigma[i]);
	}
	return GSL_SUCCESS;
}
int sersic_profile_fix_v_df (const gsl_vector * x, void *data,  gsl_matrix * J){ //ALTERAR
	size_t n = ((struct data_pos *)data)->n;
	double *r = ((struct data_pos *)data)->r; 
	double *sigma = ((struct data_pos *)data)->sigma;
	
	double rho0 = gsl_vector_get (x, 0);
	double r_s = gsl_vector_get (x, 1);
	double v = gsl_vector_get (x, 2);
	//double p = gsl_vector_get (x, 3);
	double a= - 0.6097; double b = 0.05563;
	double pp = .5 * (1.+a* v + b*pow(v,2.));
	double aa = r_s*pow(2,(1./v));
	
	double F [n];
	
	size_t i;
	for (i = 0; i < n; i++)
	{
		//Model Yi = rho0 *(r/a')^p' * exp(-(r/a')^v) 
		//F[i] = rho0 * pow(r[i]/a,-p) * exp (-pow(r[i]/a,v));
		F[i] = rho0 * pow(r[i]/aa,pp) * exp (-pow(r[i]/aa,v));
	}
	
	for (i = 0; i < n; i++) {
		/* Jacobian matrix J(i,j) = dfi / dxj, */
		/* where fi = (Yi - yi)/sigma[i],      */
		/* and the xj are the parameters (rho0,r_s,alpha,beta) */
		
		gsl_matrix_set(J,i,0, (F[i]/sigma[i]) * (1./rho0)); //rho0
		gsl_matrix_set(J,i,1, (F[i]/sigma[i]) * (pp + pow(r[i]/aa,v))/r_s); //r_s
	//	gsl_matrix_set(J,i,2, (F[i]/sigma[i]) * .5 * (log(r[i]/r_s)*(-a-2.*b*pow(v,2.))+
	//												  log(2)*(((-a-2.*b*v)/v) - (pp/pow(v,2.)))+
	//												  pow(r[i]/r_s,v)*log(r[i]/r_s))); //v
		
		//	gsl_matrix_set(J,i,0, (F[i]/sigma[i]) * (1./rho0)); //rho0
		//	gsl_matrix_set(J,i,1, (F[i]/(sigma[i]*a))* (p -v*pow(r[i]/a,v)));
		//	gsl_matrix_set(J,i,3, (F[i]/sigma[i]) * log(a/r[i]) * pow(-r[i]/a,v));
		//	gsl_matrix_set(J,i,2, (F[i]/sigma[i]) * log(a/r[i]));
		
		
		
	}
	return GSL_SUCCESS;
}
int sersic_profile_fix_v_fdf (const gsl_vector * x, void *data, gsl_vector * f, gsl_matrix * J) {
	sersic_profile_fix_v_f (x, data, f);
	sersic_profile_fix_v_df (x, data, J);
	return GSL_SUCCESS;
}



int optimalDoublePL_Mf(const gsl_vector * x, void *data, gsl_vector * f) {
	size_t n = ((struct data_pos *)data)->n;
	double *r = ((struct data_pos *)data)->r; 
	double *y = ((struct data_pos *)data)->y;
	double *sigma = ((struct data_pos *)data)->sigma;
	
	double rho0 = gsl_vector_get (x, 0);
	double r_s = gsl_vector_get (x, 1);
	double alpha = gsl_vector_get (x, 2);
	double beta = gsl_vector_get (x, 3);
	
	size_t i;
	double *Y; Y = new double [n];
	/* Model Yi = 4*pi*rho0 * r_s^beta * INT_0^(r_i) r^(2-alpha) * (r+r_s)^(alpha-beta) dr */
	
	for (i = 0; i < n; i++) Y[i] = 4.*pi*rho0*pow(r_s,beta)*pow(r[i],(2.-alpha))*pow((r[i]+r_s),(alpha-beta)); 
	
	Y = integralTrapezeZeroToPointP(n, r, Y);
	
	for (i = 0; i < n; i++) gsl_vector_set (f, i, (Y[i] - y[i])/sigma[i]);
	return GSL_SUCCESS;
}
int optimalDoublePL_Mdf (const gsl_vector * x, void *data,  gsl_matrix * J){
	size_t n = ((struct data_pos *)data)->n;
	double *r = ((struct data_pos *)data)->r;
	double *sigma = ((struct data_pos *)data)->sigma;
	
	double rho0 = gsl_vector_get (x, 0);
	double r_s = gsl_vector_get (x, 1);
	double alpha = gsl_vector_get (x, 2);
	double beta = gsl_vector_get (x, 3);
	
	size_t i;
	
	double *Y, *F0, *F1, *F2, *F3;
	Y = new double [n]; F0 = new double [n]; F1 = new double [n]; F2 = new double [n]; F3 = new double [n];
	
	
	for (i = 0; i < n; i++) Y[i] = 4.*pi*rho0*pow(r_s,beta)*pow(r[i],(2.-alpha))*pow((r[i]+r_s),(alpha-beta)); 
	
	for (i = 0; i < n; i++) F0[i] = Y[i] * (1./rho0); F0 = integralTrapezeZeroToPointP(n, r, F0);
	for (i = 0; i < n; i++) F1[i] = Y[i] * (1./(r[i]+r_s) * (alpha + beta * (r[i]/r_s))); F1 = integralTrapezeZeroToPointP(n, r, F1);
	for (i = 0; i < n; i++) F2[i] = Y[i] * log ((r[i]+r_s)/r[i]); F2 = integralTrapezeZeroToPointP(n, r, F2);
	for (i = 0; i < n; i++) F3[i] = Y[i] * log (r_s/(r[i]+r_s)); F3 = integralTrapezeZeroToPointP(n, r, F3);
	
	
	
	
	
	for (i = 0; i < n; i++) {
		/* Jacobian matrix J(i,j) = dfi / dxj, */
		/* where fi = (Yi - yi)/sigma[i],      */
		/*       Yi =rho0 /(pow(r[i]/r_s,alpha) * pow(1. + (r[i]/r_s),beta-alpha)) */
		/* and the xj are the parameters (rho0,r_s,alpha,beta) */
		
		gsl_matrix_set(J,i,0, F0[i]/sigma[i]); //rho0
		gsl_matrix_set(J,i,1, F1[i]/sigma[i]); //r_s
		gsl_matrix_set(J,i,2, F2[i]/sigma[i]); //alpha
		gsl_matrix_set(J,i,3, F3[i]/sigma[i]); //beta
		
	}
	return GSL_SUCCESS;
}
int optimalDoublePL_Mfdf (const gsl_vector * x, void *data, gsl_vector * f, gsl_matrix * J) {
	optimalDoublePL_Mf (x, data, f);
	optimalDoublePL_Mdf (x, data, J);
	return GSL_SUCCESS;
}

double * cat_fit_doublePL_M(double * x_init, void *dados) {
	
	int status;
	int N_data_points = ((struct mini_cluster *)dados)->n;
	double *r = ((struct mini_cluster *)dados)->r; 
	double *y = ((struct mini_cluster *)dados)->density;
	double *sigma = ((struct mini_cluster *)dados)->sigma;
	
	unsigned int i, iter = 0;
	const size_t n = N_data_points;
	const size_t p = 4;
	
	struct data_pos d = { n, r, y, sigma}; //d will have the real data 
	gsl_vector_view x = gsl_vector_view_array (x_init, p); //initial guess for parameters x_init={rho0, r_s, alpha, beta} 
	
	gsl_multifit_function_fdf f;//defines a general system of functions with parameters and the corresponding Jacobian matrix of derivatives 
	
	f.f = &optimalDoublePL_Mf; //stores the vector result f(x,params) in f for argument x and parameters params
	f.df = &optimalDoublePL_Mdf; //stores the n-by-p matrix result J_ij in J for argument x and parameters params 
	f.fdf = &optimalDoublePL_Mfdf; //make this thing work fast: compute the function and its derivative at the same time!
	//-sets the values of the f and J as above, for arguments x and parameters params.
	f.n = n; //number of components of the vector f
	f.p = p; //number of components of the vectors x - independent variables
	f.params = &d; //a pointer to the parameters of the function
	
	const gsl_multifit_fdfsolver_type *T = gsl_multifit_fdfsolver_lmsder; //Levenberg-Marquardt algorithm (ver:lmder routine in minpack)
	gsl_multifit_fdfsolver *s = gsl_multifit_fdfsolver_alloc (T, n, p); //returns a pointer to a derivative solver of type T for n observations and p parameters
	gsl_multifit_fdfsolver_set (s, &f, &x.vector); //initialises solver s to use the function and derivative fdf and the initial guess x.
	
	//print_doublePL_state (iter, s);
	
	int stop1, stop2; 
	gsl_vector *gradt;
	gradt = gsl_vector_alloc(p);
	
	do{
		iter++;
		status = gsl_multifit_fdfsolver_iterate (s);//perform a single iteration of the solver s
		// to track the progress of the solution: x-current position, f-function value at the current position, dx-difference between the current position and the previous position, J-J at the current postion
		
		//printf ("status = %s\n", gsl_strerror (status));
		//print_doublePL_state (iter, s);
		//printf ("status = %s \n", gsl_strerror (status));
		
		
		if (status) break;
		
		// Deciding when to stop 
		stop1 = gsl_multifit_test_delta(s->dx, s->x, 1e-4, 1e-4);
		gsl_multifit_gradient(s->J, s->f, gradt);//computes the gradient g of \Phi(x) = (1/2) ||F(x)||^2 from J and f (g=J^T f)
		stop2 = gsl_multifit_test_gradient(gradt, 1e-4);
		
	} while (stop1 == GSL_CONTINUE || (stop2==GSL_CONTINUE) && iter < 500);
	
	gsl_matrix *covar = gsl_matrix_alloc (p, p);
	gsl_multifit_covar (s->J, 0.0, covar); //uses J to compute COVAR of the best-fit parameters.COVAR = (J^T J)^{-1}
	
	/* NOTA:
	 If f_i = (Y(x, t_i) - y_i) / \sigma_i -> COVAR gives the statistical error on the best-fit parameters x
	 resulting from the gaussian errors \sigma_i on data y_i.
	 This can be verified from the relation \delta f = J \delta c and the fact that the fluctuations in f from the data y_i are normalised by \sigma_i and so satisfy <\delta f \delta f^T> = I.
	 
	 If f_i = (Y(x, t_i) - y_i) -> COVAR should be multiplied by the variance of the residuals about the best-fit
	 \sigma^2 = \sum (y_i - Y(x,t_i))^2 / (n-p) to give the variance-covariance matrix \sigma^2 C.
	 This estimates the statistical error on the best-fit parameters from the scatter of the underlying data.
	 */
	
#define FIT(i) gsl_vector_get(s->x, i) //saving params in FIT()
#define ERR(i) sqrt(gsl_matrix_get(covar,i,i)) //saving COVAR in ERR() 
	
	
	double chi = gsl_blas_dnrm2(s->f); //computes: chi=||f|| = sqrt[sum_i[f_i^2]] 
	double dof = n - p; //degrees of freedom
	double c = GSL_MAX_DBL(1, chi / sqrt(dof));//returns the MAX of the DOUBLE precision (a,b)*
	printf("chisq/dof = %g \n",  pow(chi, 2.0) / dof);
	
	printf ("rho0      = %g +/- %g ; our guess: %g \n", FIT(0), c*ERR(0), gsl_vector_get(&x.vector,0));
	printf ("r_s = %g +/- %g ; our guess: %g \n", FIT(1), c*ERR(1), gsl_vector_get(&x.vector,1));
	printf ("alpha     = %g +/- %g ; our guess: %g \n", FIT(2), c*ERR(2), gsl_vector_get(&x.vector,2));
	printf ("beta     = %g +/- %g ; our guess: %g \n", FIT(3), c*ERR(3), gsl_vector_get(&x.vector,3));
	
	
	printf ("status = %s\n", gsl_strerror (status));
	
	double * par;
	par = new double [8];
	
	par[0] = FIT(0);
	par[1] = FIT(1);
	par[2] = FIT(2);
	par[3] = FIT(3);
	par[4] = c*ERR(0);
	par[5] = c*ERR(1);
	par[6] = c*ERR(2);	
	par[7] = c*ERR(3);
	
	gsl_multifit_fdfsolver_free (s);
	gsl_matrix_free (covar);
	
	return par;	
	
}

int fix_alpha_DoublePL_Mf(const gsl_vector * x, void *data, gsl_vector * f) {
	size_t n = ((struct data_pos *)data)->n;
	double *r = ((struct data_pos *)data)->r; 
	double *y = ((struct data_pos *)data)->y;
	double *sigma = ((struct data_pos *)data)->sigma;
	
	double rho0 = gsl_vector_get (x, 0);
	double r_s = gsl_vector_get (x, 1);
	double alpha = gsl_vector_get (x, 3);
	double beta = gsl_vector_get (x, 2);
	
	size_t i;
	double *Y; Y = new double [n];
	/* Model Yi = 4*pi*rho0 * r_s^beta * INT_0^(r_i) r^(2-alpha) * (r+r_s)^(alpha-beta) dr */
	
	for (i = 0; i < n; i++) Y[i] = 4.*pi*rho0*pow(r_s,beta)*pow(r[i],(2.-alpha))*pow((r[i]+r_s),(alpha-beta)); 
	
	Y = integralTrapezeZeroToPointP(n, r, Y);
	
	for (i = 0; i < n; i++) gsl_vector_set (f, i, (Y[i] - y[i])/sigma[i]);
	return GSL_SUCCESS;
}
int fix_alpha_DoublePL_Mdf (const gsl_vector * x, void *data,  gsl_matrix * J){
	size_t n = ((struct data_pos *)data)->n;
	double *r = ((struct data_pos *)data)->r;
	double *sigma = ((struct data_pos *)data)->sigma;
	
	double rho0 = gsl_vector_get (x, 0);
	double r_s = gsl_vector_get (x, 1);
	double alpha = gsl_vector_get (x, 3);
	double beta = gsl_vector_get (x, 2);
	
	size_t i;
	
	//double *Y, *F0, *F1, *F2, *F3;
	//Y = new double [n]; F0 = new double [n]; F1 = new double [n]; F2 = new double [n]; F3 = new double [n];
	double *Y, *F0, *F1, *F2;
	Y = new double [n]; F0 = new double [n]; F1 = new double [n]; F2 = new double [n]; 
	
	
	for (i = 0; i < n; i++) Y[i] = 4.*pi*rho0*pow(r_s,beta)*pow(r[i],(2.-alpha))*pow((r[i]+r_s),(alpha-beta)); 
	
	for (i = 0; i < n; i++) F0[i] = Y[i] * (1./rho0); F0 = integralTrapezeZeroToPointP(n, r, F0);
	for (i = 0; i < n; i++) F1[i] = Y[i] * (1./(r[i]+r_s) * (alpha + beta * (r[i]/r_s))); F1 = integralTrapezeZeroToPointP(n, r, F1);
	//for (i = 0; i < n; i++) F2[i] = Y[i] * log ((r[i]+r_s)/r[i]); F2 = integralTrapezeZeroToPointP(n, r, F2);
	//for (i = 0; i < n; i++) F3[i] = Y[i] * log (r_s/(r[i]+r_s)); F3 = integralTrapezeZeroToPointP(n, r, F3);
	for (i = 0; i < n; i++) F2[i] = Y[i] * log (r_s/(r[i]+r_s)); F2 = integralTrapezeZeroToPointP(n, r, F2);
	
	
	
	
	
	for (i = 0; i < n; i++) {
		/* Jacobian matrix J(i,j) = dfi / dxj, */
		/* where fi = (Yi - yi)/sigma[i],      */
		/*       Yi =rho0 /(pow(r[i]/r_s,alpha) * pow(1. + (r[i]/r_s),beta-alpha)) */
		/* and the xj are the parameters (rho0,r_s,alpha,beta) */
		
		gsl_matrix_set(J,i,0, F0[i]/sigma[i]); //rho0
		gsl_matrix_set(J,i,1, F1[i]/sigma[i]); //r_s
		gsl_matrix_set(J,i,2, F2[i]/sigma[i]); //beta
		//alpha
		//gsl_matrix_set(J,i,3, F3[i]/sigma[i]); //beta
		
	}
	return GSL_SUCCESS;
}
int fix_alpha_DoublePL_Mfdf (const gsl_vector * x, void *data, gsl_vector * f, gsl_matrix * J) {
	fix_alpha_DoublePL_Mf (x, data, f);
	fix_alpha_DoublePL_Mdf (x, data, J);
	return GSL_SUCCESS;
}

double * cat_fit_doublePL_M_fix_alpha(double * x_init, void *dados) {
	
	int status;
	int N_data_points = ((struct mini_cluster *)dados)->n;
	double *r = ((struct mini_cluster *)dados)->r; 
	double *y = ((struct mini_cluster *)dados)->density;
	double *sigma = ((struct mini_cluster *)dados)->sigma;
	
	unsigned int i, iter = 0;
	const size_t n = N_data_points;
	const size_t p = 4;
	
	//change 3 with 2 so that we change alpha
	double save = x_init[3];
	x_init[3] = x_init[2];
	x_init[2] = save;
	struct data_pos d = { n, r, y, sigma}; //d will have the real data 
	gsl_vector_view x = gsl_vector_view_array (x_init, p); //initial guess for parameters x_init={rho0, r_s, alpha, beta} 
	
	gsl_multifit_function_fdf f;//defines a general system of functions with parameters and the corresponding Jacobian matrix of derivatives 
	
	f.f = &fix_alpha_DoublePL_Mf; //stores the vector result f(x,params) in f for argument x and parameters params
	f.df = &fix_alpha_DoublePL_Mdf; //stores the n-by-p matrix result J_ij in J for argument x and parameters params 
	f.fdf = &fix_alpha_DoublePL_Mfdf; //make this thing work fast: compute the function and its derivative at the same time!
	//-sets the values of the f and J as above, for arguments x and parameters params.
	f.n = n; //number of components of the vector f
	f.p = p; //number of components of the vectors x - independent variables
	f.params = &d; //a pointer to the parameters of the function
	
	const gsl_multifit_fdfsolver_type *T = gsl_multifit_fdfsolver_lmsder; //Levenberg-Marquardt algorithm (ver:lmder routine in minpack)
	gsl_multifit_fdfsolver *s = gsl_multifit_fdfsolver_alloc (T, n, p); //returns a pointer to a derivative solver of type T for n observations and p parameters
	gsl_multifit_fdfsolver_set (s, &f, &x.vector); //initialises solver s to use the function and derivative fdf and the initial guess x.
	
	//print_doublePL_state (iter, s);
	
	int stop1, stop2; 
	gsl_vector *gradt;
	gradt = gsl_vector_alloc(p);
	
	do{
		iter++;
		status = gsl_multifit_fdfsolver_iterate (s);//perform a single iteration of the solver s
		// to track the progress of the solution: x-current position, f-function value at the current position, dx-difference between the current position and the previous position, J-J at the current postion
		
		//printf ("status = %s\n", gsl_strerror (status));
		//print_doublePL_state (iter, s);
		//printf ("status = %s \n", gsl_strerror (status));
		
		
		if (status) break;
		
		// Deciding when to stop 
		stop1 = gsl_multifit_test_delta(s->dx, s->x, 1e-4, 1e-4);
		gsl_multifit_gradient(s->J, s->f, gradt);//computes the gradient g of \Phi(x) = (1/2) ||F(x)||^2 from J and f (g=J^T f)
		stop2 = gsl_multifit_test_gradient(gradt, 1e-4);
		
	} while (stop1 == GSL_CONTINUE || (stop2==GSL_CONTINUE) && iter < 500);
	
	gsl_matrix *covar = gsl_matrix_alloc (p, p);
	gsl_multifit_covar (s->J, 0.0, covar); //uses J to compute COVAR of the best-fit parameters.COVAR = (J^T J)^{-1}
	
	/* NOTA:
	 If f_i = (Y(x, t_i) - y_i) / \sigma_i -> COVAR gives the statistical error on the best-fit parameters x
	 resulting from the gaussian errors \sigma_i on data y_i.
	 This can be verified from the relation \delta f = J \delta c and the fact that the fluctuations in f from the data y_i are normalised by \sigma_i and so satisfy <\delta f \delta f^T> = I.
	 
	 If f_i = (Y(x, t_i) - y_i) -> COVAR should be multiplied by the variance of the residuals about the best-fit
	 \sigma^2 = \sum (y_i - Y(x,t_i))^2 / (n-p) to give the variance-covariance matrix \sigma^2 C.
	 This estimates the statistical error on the best-fit parameters from the scatter of the underlying data.
	 */
	
#define FIT(i) gsl_vector_get(s->x, i) //saving params in FIT()
#define ERR(i) sqrt(gsl_matrix_get(covar,i,i)) //saving COVAR in ERR() 
	
	
	double chi = gsl_blas_dnrm2(s->f); //computes: chi=||f|| = sqrt[sum_i[f_i^2]] 
	double dof = n - p; //degrees of freedom
	double c = GSL_MAX_DBL(1, chi / sqrt(dof));//returns the MAX of the DOUBLE precision (a,b)*
	printf("chisq/dof = %g \n",  pow(chi, 2.0) / dof);
	
	printf ("rho0      = %g +/- %g ; our guess: %g \n", FIT(0), c*ERR(0), gsl_vector_get(&x.vector,0));
	printf ("r_s = %g +/- %g ; our guess: %g \n", FIT(1), c*ERR(1), gsl_vector_get(&x.vector,1));
	printf ("alpha     = %g +/- %g ; our guess: %g \n",  x_init[3],0., gsl_vector_get(&x.vector,3));
	printf ("beta     = %g +/- %g ; our guess: %g \n", FIT(2), c*ERR(2), gsl_vector_get(&x.vector,2));
	
	
	printf ("status = %s\n", gsl_strerror (status));
	
	double * par;
	par = new double [8];
	
	par[0] = FIT(0);
	par[1] = FIT(1);
	par[2] = x_init[3];
	par[3] = FIT(2);
	par[4] = c*ERR(0);
	par[5] = c*ERR(1);
	par[7] = 0.;
	par[6] = c*ERR(2);
	
	
	gsl_multifit_fdfsolver_free (s);
	gsl_matrix_free (covar);
	
	return par;	
	
}

int fix_r_s_beta_DoublePL_Mf(const gsl_vector * x, void *data, gsl_vector * f) {
	size_t n = ((struct data_pos *)data)->n;
	double *r = ((struct data_pos *)data)->r; 
	double *y = ((struct data_pos *)data)->y;
	double *sigma = ((struct data_pos *)data)->sigma;
	
	double rho0 = gsl_vector_get (x, 0);
	double r_s = gsl_vector_get (x, 2);
	double alpha = gsl_vector_get (x, 1);
	double beta = gsl_vector_get (x, 3);
	
	size_t i;
	double *Y; Y = new double [n];
	/* Model Yi = 4*pi*rho0 * r_s^beta * INT_0^(r_i) r^(2-alpha) * (r+r_s)^(alpha-beta) dr */
	
	for (i = 0; i < n; i++) Y[i] = 4.*pi*rho0*pow(r_s,beta)*pow(r[i],(2.-alpha))*pow((r[i]+r_s),(alpha-beta)); 
	
	Y = integralTrapezeZeroToPointP(n, r, Y);
	
	for (i = 0; i < n; i++) gsl_vector_set (f, i, (Y[i] - y[i])/sigma[i]);
	return GSL_SUCCESS;
}
int fix_r_s_beta_DoublePL_Mdf (const gsl_vector * x, void *data,  gsl_matrix * J){
	size_t n = ((struct data_pos *)data)->n;
	double *r = ((struct data_pos *)data)->r;
	double *sigma = ((struct data_pos *)data)->sigma;
	
	double rho0 = gsl_vector_get (x, 0);
	double r_s = gsl_vector_get (x, 2);
	double alpha = gsl_vector_get (x, 1);
	double beta = gsl_vector_get (x, 3);
	
	size_t i;
	
	//double *Y, *F0, *F1, *F2, *F3;
	//Y = new double [n]; F0 = new double [n]; F1 = new double [n]; F2 = new double [n]; F3 = new double [n];
	double *Y, *F0, *F1;//, *F2;
	Y = new double [n]; F0 = new double [n]; F1 = new double [n]; //F2 = new double [n]; 
	
	
	for (i = 0; i < n; i++) Y[i] = 4.*pi*rho0*pow(r_s,beta)*pow(r[i],(2.-alpha))*pow((r[i]+r_s),(alpha-beta)); 
	
	
	for (i = 0; i < n; i++) F0[i] = Y[i] * (1./rho0); F0 = integralTrapezeZeroToPointP(n, r, F0);
	//for (i = 0; i < n; i++) F1[i] = Y[i] * (1./(r[i]+r_s) * (alpha + beta * (r[i]/r_s))); F1 = integralTrapezeZeroToPointP(n, r, F1);
	for (i = 0; i < n; i++) F1[i] = Y[i] * log ((r[i]+r_s)/r[i]); F1 = integralTrapezeZeroToPointP(n, r, F1);
	//for (i = 0; i < n; i++) F3[i] = Y[i] * log (r_s/(r[i]+r_s)); F3 = integralTrapezeZeroToPointP(n, r, F3);
	
	
	
	
	
	for (i = 0; i < n; i++) {
		/* Jacobian matrix J(i,j) = dfi / dxj, */
		/* where fi = (Yi - yi)/sigma[i],      */
		/*       Yi =rho0 /(pow(r[i]/r_s,alpha) * pow(1. + (r[i]/r_s),beta-alpha)) */
		/* and the xj are the parameters (rho0,r_s,alpha,beta) */
		
		gsl_matrix_set(J,i,0, F0[i]/sigma[i]); //rho0
		//gsl_matrix_set(J,i,1, F1[i]/sigma[i]); //r_s
		gsl_matrix_set(J,i,1, F1[i]/sigma[i]); //alpha
		//gsl_matrix_set(J,i,3, F3[i]/sigma[i]); //beta
		
	}
	return GSL_SUCCESS;
}
int fix_r_s_beta_DoublePL_Mfdf (const gsl_vector * x, void *data, gsl_vector * f, gsl_matrix * J) {
	fix_r_s_beta_DoublePL_Mf (x, data, f);
	fix_r_s_beta_DoublePL_Mdf (x, data, J);
	return GSL_SUCCESS;
}
double * cat_fit_doublePL_M_fix_r_s_beta(double * x_init, void *dados) {
	
	int status;
	int N_data_points = ((struct mini_cluster *)dados)->n;
	double *r = ((struct mini_cluster *)dados)->r; 
	double *y = ((struct mini_cluster *)dados)->density;
	double *sigma = ((struct mini_cluster *)dados)->sigma;
	
	unsigned int i, iter = 0;
	const size_t n = N_data_points;
	const size_t p = 4;
	
	//change 2 with 1 so that we change alpha and not r_s
	double save = x_init[2];
	x_init[2] = x_init[1];
	x_init[1] = save;
	struct data_pos d = { n, r, y, sigma}; //d will have the real data 
	gsl_vector_view x = gsl_vector_view_array (x_init, p); //initial guess for parameters x_init={rho0, r_s, alpha, beta} 
	
	gsl_multifit_function_fdf f;//defines a general system of functions with parameters and the corresponding Jacobian matrix of derivatives 
	
	f.f = &fix_r_s_beta_DoublePL_Mf; //stores the vector result f(x,params) in f for argument x and parameters params
	f.df = &fix_r_s_beta_DoublePL_Mdf; //stores the n-by-p matrix result J_ij in J for argument x and parameters params 
	f.fdf = &fix_r_s_beta_DoublePL_Mfdf; //make this thing work fast: compute the function and its derivative at the same time!
	//-sets the values of the f and J as above, for arguments x and parameters params.
	f.n = n; //number of components of the vector f
	f.p = p; //number of components of the vectors x - independent variables
	f.params = &d; //a pointer to the parameters of the function
	
	const gsl_multifit_fdfsolver_type *T = gsl_multifit_fdfsolver_lmsder; //Levenberg-Marquardt algorithm (ver:lmder routine in minpack)
	gsl_multifit_fdfsolver *s = gsl_multifit_fdfsolver_alloc (T, n, p); //returns a pointer to a derivative solver of type T for n observations and p parameters
	gsl_multifit_fdfsolver_set (s, &f, &x.vector); //initialises solver s to use the function and derivative fdf and the initial guess x.
	
	//print_doublePL_state (iter, s);
	
	int stop1, stop2; 
	gsl_vector *gradt;
	gradt = gsl_vector_alloc(p);
	
	do{
		iter++;
		status = gsl_multifit_fdfsolver_iterate (s);//perform a single iteration of the solver s
		// to track the progress of the solution: x-current position, f-function value at the current position, dx-difference between the current position and the previous position, J-J at the current postion
		
		//printf ("status = %s\n", gsl_strerror (status));
		//print_doublePL_state (iter, s);
		//printf ("status = %s \n", gsl_strerror (status));
		
		
		if (status) break;
		
		// Deciding when to stop 
		stop1 = gsl_multifit_test_delta(s->dx, s->x, 1e-4, 1e-4);
		gsl_multifit_gradient(s->J, s->f, gradt);//computes the gradient g of \Phi(x) = (1/2) ||F(x)||^2 from J and f (g=J^T f)
		stop2 = gsl_multifit_test_gradient(gradt, 1e-4);
		
	} while (stop1 == GSL_CONTINUE || (stop2==GSL_CONTINUE) && iter < 500);
	
	gsl_matrix *covar = gsl_matrix_alloc (p, p);
	gsl_multifit_covar (s->J, 0.0, covar); //uses J to compute COVAR of the best-fit parameters.COVAR = (J^T J)^{-1}
	
	/* NOTA:
	 If f_i = (Y(x, t_i) - y_i) / \sigma_i -> COVAR gives the statistical error on the best-fit parameters x
	 resulting from the gaussian errors \sigma_i on data y_i.
	 This can be verified from the relation \delta f = J \delta c and the fact that the fluctuations in f from the data y_i are normalised by \sigma_i and so satisfy <\delta f \delta f^T> = I.
	 
	 If f_i = (Y(x, t_i) - y_i) -> COVAR should be multiplied by the variance of the residuals about the best-fit
	 \sigma^2 = \sum (y_i - Y(x,t_i))^2 / (n-p) to give the variance-covariance matrix \sigma^2 C.
	 This estimates the statistical error on the best-fit parameters from the scatter of the underlying data.
	 */
	
#define FIT(i) gsl_vector_get(s->x, i) //saving params in FIT()
#define ERR(i) sqrt(gsl_matrix_get(covar,i,i)) //saving COVAR in ERR() 
	
	
	double chi = gsl_blas_dnrm2(s->f); //computes: chi=||f|| = sqrt[sum_i[f_i^2]] 
	double dof = n - p; //degrees of freedom
	double c = GSL_MAX_DBL(1, chi / sqrt(dof));//returns the MAX of the DOUBLE precision (a,b)*
	printf("chisq/dof = %g \n",  pow(chi, 2.0) / dof);
	
	printf ("rho0      = %g +/- %g ; our guess: %g \n", FIT(0), c*ERR(0), gsl_vector_get(&x.vector,0));
	printf ("r_s = %g +/- %g ; our guess: %g \n", x_init[2], 0., gsl_vector_get(&x.vector,1));
	printf ("alpha     = %g +/- %g ; our guess: %g \n", FIT(1),c*ERR(1), gsl_vector_get(&x.vector,3));
	printf ("beta     = %g +/- %g ; our guess: %g \n", x_init[3], 0., gsl_vector_get(&x.vector,2));
	
	
	printf ("status = %s\n", gsl_strerror (status));
	
	double * par;
	par = new double [8];
	
	par[0] = FIT(0);
	par[1] =  x_init[2];
	par[2] = FIT(1);
	par[3] = x_init[3];
	par[4] = c*ERR(0);
	par[5] = 0.;
	par[7] = c*ERR(1);
	par[6] = 0.;
	
	
	gsl_multifit_fdfsolver_free (s);
	gsl_matrix_free (covar);
	
	return par;	
	
}


int fix_alpha_beta_DoublePL_Mf(const gsl_vector * x, void *data, gsl_vector * f) {
	size_t n = ((struct data_pos *)data)->n;
	double *r = ((struct data_pos *)data)->r; 
	double *y = ((struct data_pos *)data)->y;
	double *sigma = ((struct data_pos *)data)->sigma;
	
	double rho0 = gsl_vector_get (x, 0);
	double r_s = gsl_vector_get (x, 1);
	double alpha = gsl_vector_get (x, 2);
	double beta = gsl_vector_get (x, 3);
	
	size_t i;
	double *Y; Y = new double [n];
	/* Model Yi = 4*pi*rho0 * r_s^beta * INT_0^(r_i) r^(2-alpha) * (r+r_s)^(alpha-beta) dr */
	
	for (i = 0; i < n; i++) Y[i] = 4.*pi*rho0*pow(r_s,beta)*pow(r[i],(2.-alpha))*pow((r[i]+r_s),(alpha-beta)); 
	
	Y = integralTrapezeZeroToPointP(n, r, Y);
	
	for (i = 0; i < n; i++) gsl_vector_set (f, i, (Y[i] - y[i])/sigma[i]);
	return GSL_SUCCESS;
}
int fix_alpha_beta_DoublePL_Mdf (const gsl_vector * x, void *data,  gsl_matrix * J){
	size_t n = ((struct data_pos *)data)->n;
	double *r = ((struct data_pos *)data)->r;
	double *sigma = ((struct data_pos *)data)->sigma;
	
	double rho0 = gsl_vector_get (x, 0);
	double r_s = gsl_vector_get (x, 1);
	double alpha = gsl_vector_get (x, 2);
	double beta = gsl_vector_get (x, 3);
	
	size_t i;
	
	//double *Y, *F0, *F1, *F2, *F3;
	//Y = new double [n]; F0 = new double [n]; F1 = new double [n]; F2 = new double [n]; F3 = new double [n];
	double *Y, *F0, *F1;
	Y = new double [n]; F0 = new double [n]; F1 = new double [n];
	
	
	for (i = 0; i < n; i++) Y[i] = 4.*pi*rho0*pow(r_s,beta)*pow(r[i],(2.-alpha))*pow((r[i]+r_s),(alpha-beta)); 
	
	for (i = 0; i < n; i++) F0[i] = Y[i] * (1./rho0); F0 = integralTrapezeZeroToPointP(n, r, F0);
	for (i = 0; i < n; i++) F1[i] = Y[i] * (1./(r[i]+r_s) * (alpha + beta * (r[i]/r_s))); F1 = integralTrapezeZeroToPointP(n, r, F1);
	//for (i = 0; i < n; i++) F2[i] = Y[i] * log ((r[i]+r_s)/r[i]); F2 = integralTrapezeZeroToPointP(n, r, F2);
	//for (i = 0; i < n; i++) F3[i] = Y[i] * log (r_s/(r[i]+r_s)); F3 = integralTrapezeZeroToPointP(n, r, F3);
	//for (i = 0; i < n; i++) F2[i] = Y[i] * log (r_s/(r[i]+r_s)); F2 = integralTrapezeZeroToPointP(n, r, F2);
	
	
	
	
	
	for (i = 0; i < n; i++) {
		/* Jacobian matrix J(i,j) = dfi / dxj, */
		/* where fi = (Yi - yi)/sigma[i],      */
		/*       Yi =rho0 /(pow(r[i]/r_s,alpha) * pow(1. + (r[i]/r_s),beta-alpha)) */
		/* and the xj are the parameters (rho0,r_s,alpha,beta) */
		
		gsl_matrix_set(J,i,0, F0[i]/sigma[i]); //rho0
		gsl_matrix_set(J,i,1, F1[i]/sigma[i]); //r_s
		//gsl_matrix_set(J,i,2, F2[i]/sigma[i]); //beta
		//alpha
		//gsl_matrix_set(J,i,3, F3[i]/sigma[i]); //beta
		
	}
	return GSL_SUCCESS;
}
int fix_alpha_beta_DoublePL_Mfdf (const gsl_vector * x, void *data, gsl_vector * f, gsl_matrix * J) {
	fix_alpha_beta_DoublePL_Mf (x, data, f);
	fix_alpha_beta_DoublePL_Mdf (x, data, J);
	return GSL_SUCCESS;
}

double * cat_fit_doublePL_M_fix_alpha_beta(double * x_init, void *dados) {
	
	int status;
	int N_data_points = ((struct mini_cluster *)dados)->n;
	double *r = ((struct mini_cluster *)dados)->r; 
	double *y = ((struct mini_cluster *)dados)->density;
	double *sigma = ((struct mini_cluster *)dados)->sigma;
	
	unsigned int i, iter = 0;
	const size_t n = N_data_points;
	const size_t p = 4;
	
	struct data_pos d = { n, r, y, sigma}; //d will have the real data 
	gsl_vector_view x = gsl_vector_view_array (x_init, p); //initial guess for parameters x_init={rho0, r_s, alpha, beta} 
	
	gsl_multifit_function_fdf f;//defines a general system of functions with parameters and the corresponding Jacobian matrix of derivatives 
	
	f.f = &fix_alpha_beta_DoublePL_Mf; //stores the vector result f(x,params) in f for argument x and parameters params
	f.df = &fix_alpha_beta_DoublePL_Mdf; //stores the n-by-p matrix result J_ij in J for argument x and parameters params 
	f.fdf = &fix_alpha_beta_DoublePL_Mfdf; //make this thing work fast: compute the function and its derivative at the same time!
	//-sets the values of the f and J as above, for arguments x and parameters params.
	f.n = n; //number of components of the vector f
	f.p = p; //number of components of the vectors x - independent variables
	f.params = &d; //a pointer to the parameters of the function
	
	const gsl_multifit_fdfsolver_type *T = gsl_multifit_fdfsolver_lmsder; //Levenberg-Marquardt algorithm (ver:lmder routine in minpack)
	gsl_multifit_fdfsolver *s = gsl_multifit_fdfsolver_alloc (T, n, p); //returns a pointer to a derivative solver of type T for n observations and p parameters
	gsl_multifit_fdfsolver_set (s, &f, &x.vector); //initialises solver s to use the function and derivative fdf and the initial guess x.
	
	//print_doublePL_state (iter, s);
	
	int stop1, stop2; 
	gsl_vector *gradt;
	gradt = gsl_vector_alloc(p);
	
	do{
		iter++;
		status = gsl_multifit_fdfsolver_iterate (s);//perform a single iteration of the solver s
		// to track the progress of the solution: x-current position, f-function value at the current position, dx-difference between the current position and the previous position, J-J at the current postion
		
		//printf ("status = %s\n", gsl_strerror (status));
		//print_doublePL_state (iter, s);
		//printf ("status = %s \n", gsl_strerror (status));
		
		
		if (status) break;
		
		// Deciding when to stop 
		stop1 = gsl_multifit_test_delta(s->dx, s->x, 1e-4, 1e-4);
		gsl_multifit_gradient(s->J, s->f, gradt);//computes the gradient g of \Phi(x) = (1/2) ||F(x)||^2 from J and f (g=J^T f)
		stop2 = gsl_multifit_test_gradient(gradt, 1e-4);
		
	} while (stop1 == GSL_CONTINUE || (stop2==GSL_CONTINUE) && iter < 500);
	
	gsl_matrix *covar = gsl_matrix_alloc (p, p);
	gsl_multifit_covar (s->J, 0.0, covar); //uses J to compute COVAR of the best-fit parameters.COVAR = (J^T J)^{-1}
	
	/* NOTA:
	 If f_i = (Y(x, t_i) - y_i) / \sigma_i -> COVAR gives the statistical error on the best-fit parameters x
	 resulting from the gaussian errors \sigma_i on data y_i.
	 This can be verified from the relation \delta f = J \delta c and the fact that the fluctuations in f from the data y_i are normalised by \sigma_i and so satisfy <\delta f \delta f^T> = I.
	 
	 If f_i = (Y(x, t_i) - y_i) -> COVAR should be multiplied by the variance of the residuals about the best-fit
	 \sigma^2 = \sum (y_i - Y(x,t_i))^2 / (n-p) to give the variance-covariance matrix \sigma^2 C.
	 This estimates the statistical error on the best-fit parameters from the scatter of the underlying data.
	 */
	
#define FIT(i) gsl_vector_get(s->x, i) //saving params in FIT()
#define ERR(i) sqrt(gsl_matrix_get(covar,i,i)) //saving COVAR in ERR() 
	
	
	double chi = gsl_blas_dnrm2(s->f); //computes: chi=||f|| = sqrt[sum_i[f_i^2]] 
	double dof = n - p; //degrees of freedom
	double c = GSL_MAX_DBL(1, chi / sqrt(dof));//returns the MAX of the DOUBLE precision (a,b)*
	printf("chisq/dof = %g \n",  pow(chi, 2.0) / dof);
	
	printf ("rho0      = %g +/- %g ; our guess: %g \n", FIT(0), c*ERR(0), gsl_vector_get(&x.vector,0));
	printf ("r_s = %g +/- %g ; our guess: %g \n", FIT(1), c*ERR(1), gsl_vector_get(&x.vector,1));
	printf ("alpha     = %g +/- %g ; our guess: %g \n",  x_init[2],0., gsl_vector_get(&x.vector,2));
	printf ("beta     = %g +/- %g ; our guess: %g \n", x_init[3],0., gsl_vector_get(&x.vector,3));
	
	
	printf ("status = %s\n", gsl_strerror (status));
	
	double * par;
	par = new double [8];
	
	par[0] = FIT(0);
	par[1] = FIT(1);
	par[2] = x_init[2];
	par[3] = x_init[3];
	par[4] = c*ERR(0);
	par[5] = c*ERR(1);
	par[7] = 0.;
	par[6] = 0.;
	
	
	gsl_multifit_fdfsolver_free (s);
	gsl_matrix_free (covar);
	
	return par;	
	
}


int quadratic_f(const gsl_vector * x, void *data, gsl_vector * f) {
	size_t n = ((struct data_pos *)data)->n;
	double *r = ((struct data_pos *)data)->r; 
	double *y = ((struct data_pos *)data)->y;
	double *sigma = ((struct data_pos *)data)->sigma;
	
	double a = gsl_vector_get (x, 0);
	double b = gsl_vector_get (x, 1);
	double c = gsl_vector_get (x, 2);
	
	size_t i;
	double *Y; Y = new double [n];
	/* Model Yi = a*pow(r,2.) + b*r + c */
	
	for (i = 0; i < n; i++) Y[i] = a*pow(r[i],2.) + b*r[i] + c; 
	
	for (i = 0; i < n; i++) gsl_vector_set (f, i, (Y[i] - y[i])/sigma[i]);
	return GSL_SUCCESS;
}
int quadratic_df (const gsl_vector * x, void *data,  gsl_matrix * J){
	size_t n = ((struct data_pos *)data)->n;
	double *r = ((struct data_pos *)data)->r;
	double *sigma = ((struct data_pos *)data)->sigma;
	
	double a = gsl_vector_get (x, 0);
	double b = gsl_vector_get (x, 1);
	double c = gsl_vector_get (x, 2);
	
	size_t i;
	//for (i = 0; i < n; i++)
	//{
		/* Model Yi = a*pow(r[i],2.) + b*r[i] + c */
	//	F[i] = a*pow(r[i],2.) + b*r[i] + c;
	//}
	
	for (i = 0; i < n; i++) {
		/* Jacobian matrix J(i,j) = dfi / dxj, */
		/* where fi = (Yi - yi)/sigma[i],      */
		/* and the xj are the parameters (rho0,r_s,alpha,beta) */
		
		gsl_matrix_set(J,i,0, (pow(r[i],2.)/sigma[i])); //a
		gsl_matrix_set(J,i,1, (r[i]/sigma[i])); //b
		gsl_matrix_set(J,i,2, 1./sigma[i]); //c
		
	}
	return GSL_SUCCESS;
	
}
int quadratic_fdf (const gsl_vector * x, void *data, gsl_vector * f, gsl_matrix * J) {
	quadratic_f(x, data, f);
	quadratic_df (x, data, J);
	return GSL_SUCCESS;
}
double * cat_fit_quadratic(double *x_init, void *dados) {
	
	/*	int n = ((struct mini_cluster *)dados)->n;
	 double *r = ((struct mini_cluster *)dados)->r; 
	 double *y = ((struct mini_cluster *)dados)->density;
	 double *sigma = ((struct mini_cluster *)dados)->sigma;
	 
	 //double xi, yi, ei, chisq;
	 double chisq;
	 gsl_matrix *X, *cov;
	 gsl_vector *Y, *w, *c;
	 
	 
	 X = gsl_matrix_alloc (n, 3);
	 Y = gsl_vector_alloc (n);
	 w = gsl_vector_alloc (n);
	 
	 c = gsl_vector_alloc (3);
	 cov = gsl_matrix_alloc (3, 3);
	 
	 for (int i = 0; i < n; i++)
	 {
	 //	int count = fscanf (stdin, "%lg %lg %lg",
	 //						&xi, &yi, &ei);
	 //	
	 //	if (count != 3)
	 //{
	 //	fprintf (stderr, "error reading file\n");
	 //	exit (-1);
	 //}
	 
	 //printf ("%g %g +/- %g\n", xi, yi, ei);
	 
	 gsl_matrix_set (X, i, 0, 1.0);
	 gsl_matrix_set (X, i, 1, r[i]);
	 gsl_matrix_set (X, i, 2, r[i]*r[i]);
	 
	 gsl_vector_set (Y, i, y[i]);
	 gsl_vector_set (w, i, 1.0/(sigma[i]*sigma[i]));
	 }
	 
	 {
	 gsl_multifit_linear_workspace * work 
	 = gsl_multifit_linear_alloc (n, 3);
	 gsl_multifit_wlinear (X, w, Y, c, cov,
	 &chisq, work);
	 gsl_multifit_linear_free (work);
	 }
	 
	 #define C(i) (gsl_vector_get(c,(i)))
	 #define COV(i,j) (gsl_matrix_get(cov,(i),(j)))
	 
	 {
	 printf ("# best fit: Y = %g + %g X + %g X^2\n", 
	 C(0), C(1), C(2));
	 
	 printf ("# covariance matrix:\n");
	 printf ("[ %+.5e, %+.5e, %+.5e  \n",
	 COV(0,0), COV(0,1), COV(0,2));
	 printf ("  %+.5e, %+.5e, %+.5e  \n", 
	 COV(1,0), COV(1,1), COV(1,2));
	 printf ("  %+.5e, %+.5e, %+.5e ]\n", 
	 COV(2,0), COV(2,1), COV(2,2));
	 printf ("# chisq = %g\n", chisq);
	 }
	 
	 gsl_matrix_free (X);
	 gsl_vector_free (Y);
	 gsl_vector_free (w);
	 //gsl_vector_free (c);
	 gsl_matrix_free (cov);
	 
	 
	 
	 
	 double * par;
	 par = new double [3];
	 
	 par[0] = C(0);
	 par[1] = C(1);
	 par[2] = C(2);
	 
	 return par;	
	 
	 
	 
	 
	 
	 */
	
	
	
	
	
	/////////////////
	
	
	
	int status;
	int N_data_points = ((struct mini_cluster *)dados)->n;
	double *r = ((struct mini_cluster *)dados)->r; 
	double *y = ((struct mini_cluster *)dados)->density;
	double *sigma = ((struct mini_cluster *)dados)->sigma;
	
	unsigned int i, iter = 0;
	const size_t n = N_data_points;
	const size_t p = 3;
	
	struct data_pos d = { n, r, y, sigma}; //d will have the real data 
	gsl_vector_view x = gsl_vector_view_array (x_init, p); //initial guess for parameters x_init={rho0, r_s, alpha, beta} 
	
	gsl_multifit_function_fdf f;//defines a general system of functions with parameters and the corresponding Jacobian matrix of derivatives 
	
	f.f = &quadratic_f; //stores the vector result f(x,params) in f for argument x and parameters params
	f.df = &quadratic_df; //stores the n-by-p matrix result J_ij in J for argument x and parameters params 
	f.fdf = &quadratic_fdf; //make this thing work fast: compute the function and its derivative at the same time!
	//-sets the values of the f and J as above, for arguments x and parameters params.
	f.n = n; //number of components of the vector f
	f.p = p; //number of components of the vectors x - independent variables
	f.params = &d; //a pointer to the parameters of the function
	
	const gsl_multifit_fdfsolver_type *T = gsl_multifit_fdfsolver_lmsder; //Levenberg-Marquardt algorithm (ver:lmder routine in minpack)
	gsl_multifit_fdfsolver *s = gsl_multifit_fdfsolver_alloc (T, n, p); //returns a pointer to a derivative solver of type T for n observations and p parameters
	
	for (int i =0; i < 3; i++) {
		printf("%d = %g \n", i, gsl_vector_get(&x.vector,i));
		printf("%lu\n",  n );
	}
	printf("before fdf solver \n");
	gsl_multifit_fdfsolver_set (s, &f, &x.vector); //initialises solver s to use the function and derivative fdf and the initial guess x.
	printf("after fdf solver \n");
	
	//print_doublePL_state (iter, s);
	
	int stop1, stop2; 
	gsl_vector *gradt;
	gradt = gsl_vector_alloc(p);
	
	do{
		iter++;
		status = gsl_multifit_fdfsolver_iterate (s);//perform a single iteration of the solver s
		// to track the progress of the solution: x-current position, f-function value at the current position, dx-difference between the current position and the previous position, J-J at the current postion
		
		//printf ("status = %s\n", gsl_strerror (status));
		//print_doublePL_state (iter, s);
		//printf ("status = %s \n", gsl_strerror (status));
		
		
		if (status) break;
		
		// Deciding when to stop 
		stop1 = gsl_multifit_test_delta(s->dx, s->x, 1.e-4, 1.e-4);
		gsl_multifit_gradient(s->J, s->f, gradt);//computes the gradient g of \Phi(x) = (1/2) ||F(x)||^2 from J and f (g=J^T f)
		stop2 = gsl_multifit_test_gradient(gradt, 1.e-4);
		
	} while (stop1 == GSL_CONTINUE || (stop2==GSL_CONTINUE) && iter < 500);
	
	gsl_matrix *covar = gsl_matrix_alloc (p, p);
	gsl_multifit_covar (s->J, 0.0, covar); //uses J to compute COVAR of the best-fit parameters.COVAR = (J^T J)^{-1}
	
	//NOTA:
	//	 If f_i = (Y(x, t_i) - y_i) / \sigma_i -> COVAR gives the statistical error on the best-fit parameters x
	//	 resulting from the gaussian errors \sigma_i on data y_i.
	//	 This can be verified from the relation \delta f = J \delta c and the fact that the fluctuations in f from the data y_i are normalised by \sigma_i and so satisfy <\delta f \delta f^T> = I.
	
	//	 If f_i = (Y(x, t_i) - y_i) -> COVAR should be multiplied by the variance of the residuals about the best-fit
	//	 \sigma^2 = \sum (y_i - Y(x,t_i))^2 / (n-p) to give the variance-covariance matrix \sigma^2 C.
	//	 This estimates the statistical error on the best-fit parameters from the scatter of the underlying data.
	
	
#define FIT(i) gsl_vector_get(s->x, i) //saving params in FIT()
#define ERR(i) sqrt(gsl_matrix_get(covar,i,i)) //saving COVAR in ERR() 
	
	
	double chi = gsl_blas_dnrm2(s->f); //computes: chi=||f|| = sqrt[sum_i[f_i^2]] 
	double dof = n - p; //degrees of freedom
	double c = GSL_MAX_DBL(1, chi / sqrt(dof));//returns the MAX of the DOUBLE precision (a,b)*
	printf("chisq/dof = %g, chisq= %g\n",  pow(chi, 2.0) / dof, pow(chi,2.));
	
	printf ("a      = %g +/- %g ; our guess: %g \n", FIT(0), c*ERR(0), gsl_vector_get(&x.vector,0));
	printf ("b = %g +/- %g ; our guess: %g \n", FIT(1), c*ERR(1), gsl_vector_get(&x.vector,1));
	printf ("c     = %g +/- %g ; our guess: %g \n", FIT(2), c*ERR(2), gsl_vector_get(&x.vector,2));
	
	
	printf ("status = %s\n", gsl_strerror(status));
	
	double * par;
	par = new double [6];
	
	par[0] = FIT(0);
	par[1] = FIT(1);
	par[2] = FIT(2);
	par[3] = c*ERR(0);
	par[4] = c*ERR(1);
	par[5] = c*ERR(2);	
	
	gsl_multifit_fdfsolver_free (s);
	gsl_matrix_free (covar);
	
	return par;	
}


/********* Linear regression *********/
/*perform least-squares fits to a straight line model, Y(c,x) = c_0 + c_1 x*/


double * cat_fit_linear (double * x, double * y, int n) {
	/*	computes the best-fit linear regression coefficients (c0,c1) of the model Y = c_0 + c_1 X for the datasets (x, y),
		two vectors of length n with strides xstride and ystride.
		variance-covariance matrix for the parameters (c0, c1) is estimated from the scatter of the points around the best-fit line
		returned via the parameters (cov00, cov01, cov11)
		sum of squares of the residuals from the best-fit line is returned in sumsq
		n = 4, fit[5]
		Function: int gsl_fit_linear (const double * x, const size_t xstride, const double * y, const size_t ystride, 
		size_t n, double * c0, double * c1, double * cov00, double * cov01, double * cov11, double * sumsq)
	*/
	double c0, c1, cov00, cov01, cov11, sumsq;
	
	double * fit;
	fit = new double[5];
	
	gsl_fit_linear ( x, 1, y, 1, n, 
					&c0, &c1, &cov00, &cov01, &cov11, 
					&sumsq);
	
	fit[0] = c0;
	fit[1] = c1;
	fit[2] = cov00;
	fit[3] = cov01;
	fit[4] = cov11;
	fit[5] = sumsq;
	
	
	printf("# best fit: Y = %g + %g X\n", fit[0], fit[1]);
	printf("# covariance matrix:\n");
	printf("# [ %g, %g\n#   %g, %g]\n", 
		   fit[2], fit[3], fit[3], fit[4]);
	printf("# sumsq = %g\n", fit[5]);
	
	return fit;
	
	
	
}

/* ....::.... Non-Linear regression :::....:
 f_i = ((rho0 /( pow(r/r_s,alpha) * pow(r/r_s,beta-alpha))) - y_i)/\sigma_i
 where we have chosen t_i = i. The Jacobian matrix J is the derivative of these functions with respect to the 4 parameters.
 It is given by, J_{ij} = d f_i / d x_j
 where x_0 = rho0, x_1 = r_s, x_2 = \alpha and x_3 = \beta.
 */



