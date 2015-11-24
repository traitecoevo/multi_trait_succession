#include <iostream>
#include <fstream>
#include "gsl_root_class.h"
#include <algorithm>

#define GSL_FAIL	-2
#define SINGLE_ROOT_TYPE  gsl_root_fsolver_brent
#define SINGLE_ROOT_TYPE2  gsl_root_fdfsolver_newton

gsl_root_class::gsl_root_class(void) {
}


double gsl_root_class::single_bracket(double (* function) (double X, void *pars),
	void * pars, double low, double up, double rel_tol, int max_iter, bool print) {
	if(print)
		std::cout << "starting SINGLE ROOT solve" << std::endl;
	// Initialise the solver
	int status, iter = 0;
	double r,x_low= low, x_up = up;

	const gsl_root_fsolver_type *T =SINGLE_ROOT_TYPE;
	gsl_root_fsolver *s = gsl_root_fsolver_alloc (T);
	gsl_function f= {function, pars}; // Provides the function to solve
	gsl_root_fsolver_set (s, &f, x_low, x_up);

	// Solve
	do
	{
		iter++;
		status = gsl_root_fsolver_iterate (s);
		r = gsl_root_fsolver_root(s);
		x_low = gsl_root_fsolver_x_lower(s);
		x_up = gsl_root_fsolver_x_upper(s);
		status= gsl_root_test_interval(x_low, x_up, 0,rel_tol);
		if(print)
			std::cout << "iter = " << iter << " " << x_low << " " << x_up<<std::endl;
	}
	while (status == GSL_CONTINUE && iter < max_iter);
	if(print)
		std::cout << "count = " << iter << ", status = " << gsl_strerror (status)<<std::endl;
	gsl_root_fsolver_free (s);
	return r;
}

// Find root of function using derivative method. Faster than bracketing if initial estimate of root is good
// Required function giving derivative, and fdf function, an optimised function giving f and df simultaneously
double gsl_root_class::single_deriv(double (* f) (double X, void *pars), double (* df) (double x, void * params ),
	void (* fdf) (double x, void * params,  double * f, double * df), void * pars, double xinit, double rel_tol, int max_iter, bool print) {
	if(print)
		std::cout << "starting SINGLE ROOT-deriv solve" << std::endl;
	// Initialise the solver
	int status, iter = 0;
	double r1,r2;

	const gsl_root_fdfsolver_type *T =SINGLE_ROOT_TYPE2;
	gsl_root_fdfsolver *s = gsl_root_fdfsolver_alloc (T);
	gsl_function_fdf FDF= {f, df, fdf, pars}; // Provides the function to solve
	gsl_root_fdfsolver_set (s, &FDF, xinit);

	r1 = xinit;
	if(print)
		std::cout << "iter = " << iter << " " << r1<<std::endl;
// Solve
	do
	{
		iter++;
		status = gsl_root_fdfsolver_iterate (s);
		r2 = gsl_root_fdfsolver_root(s);
		status= gsl_root_test_interval(r2, r1, 0,rel_tol);
		if(print)
			std::cout << "iter = " << iter << " " << r2 << " " << r2/r1<<std::endl;
		r1 =r2;
	}
	while (status == GSL_CONTINUE && iter < max_iter);
	if(print)
		std::cout << "count = " << iter << ", status = " << gsl_strerror (status)<<std::endl;
	gsl_root_fdfsolver_free (s);
	return r2;
}

int gsl_root_class::multi(int (* function) (const gsl_vector * x, void *params, gsl_vector * f), void *pars,
	double *X, int dim, double rel_tol, int max_iter, bool print) {
	if(print)
		std::cout << "Starting MULTI SOLVE" << std::endl;
	// Initialise the solver
	int status, i, iter = 0, n =dim;
	gsl_vector *x = gsl_vector_alloc (n);
	for(int i=0; i<n; i++) gsl_vector_set (x, i, X[i]);

	// Set solver type & intialsie
	const gsl_multiroot_fsolver_type *T = gsl_multiroot_fsolver_hybrid;
	gsl_multiroot_fsolver *s = gsl_multiroot_fsolver_alloc (T, n);
	gsl_multiroot_function f = {function, n, pars}; // Provides the function to solve
	gsl_multiroot_fsolver_set (s, &f, x);
	// Solve
	do
	{
		iter++;
		if(print){
			std::cout << "iter = " << iter;
			for(i =0; i< (int) n; i++)
				std::cout << " x[" << i << "] = " << gsl_vector_get (s->x, i) << " f(x)= " << gsl_vector_get (s->f, i);
			std::cout<<std::endl;
		}
		status = gsl_multiroot_fsolver_iterate (s);
		if (status) break; // check if solver is stuck
		status = gsl_multiroot_test_residual (s->f, rel_tol);
	}
	while (status == GSL_CONTINUE && iter < max_iter);
	if(print)
		std::cout << "GSL count = " << iter << ", status = " << gsl_strerror (status)<<std::endl;
	for(i=0; i<n; i++)  X[i] = gsl_vector_get (s->x, i);
	gsl_multiroot_fsolver_free (s);
	gsl_vector_free (x);
	return status;
}

int gsl_root_class::multi_jac(int (* func_f) (const gsl_vector * x, void *params, gsl_vector * f),
	int (*func_df) (const gsl_vector * x, void *params, gsl_matrix * J),
	int (*func_fdf) (const gsl_vector * x, void *params, gsl_vector * f, gsl_matrix * J),
	void *pars, double *X, int dim, double rel_tol, int max_iter, bool print){
	if(print)
		std::cout << "MULTIROOTSOLVE:\tstart" << std::endl;
	// Initialise the solver
	int status=-2, i, iter = 0, n =dim;
	gsl_vector *x = gsl_vector_alloc (n);
	for(int i=0; i<n; i++) gsl_vector_set (x, i, X[i]);

	// Set solver type & intialsie
	const gsl_multiroot_fdfsolver_type *T = gsl_multiroot_fdfsolver_hybridsj; // Gnewton
	gsl_multiroot_fdfsolver *s = gsl_multiroot_fdfsolver_alloc (T, n);
	gsl_multiroot_function_fdf f = {func_f, func_df, func_fdf, n, pars}; // Provides the function to solve
	gsl_multiroot_fdfsolver_set (s, &f, x);
	// Solve
	do
	{
		iter++;

		if(print){
			std::cout << "MULTIROOTSOLVE:\titer = " << iter;
			for(i =0; i< (int) n; i++)
				std::cout << " x[" << i << "] = " << gsl_vector_get (s->x, i) << " f(x)= " << gsl_vector_get (s->f, i);
			std::cout<<std::endl;
		}
		status = gsl_multiroot_fdfsolver_iterate (s);
		if (status) break; // check if solver is stuck

		// TEST IF FINISHED. DEFAULT =GSL VERSION TESTS SUM OF RESIDUALS
		// Status = gsl_multiroot_test_residual (s->f, rel_tol);
		// INSTEAD use option where conisder maximum of each function independently
		status = 0;  // Success
		for(i =0; i< (int) n; i++)
			if(fabs(gsl_vector_get (s->f, i)) > rel_tol) status = -2; // Hasn't cnverged yet
	}
	while (status == GSL_CONTINUE && iter < max_iter);
	if(print){
		std::cout << "MULTIROOTSOLVE:\titer = " << iter;
		for(i =0; i< (int) n; i++)
		std::cout << " x[" << i << "] = " << gsl_vector_get (s->x, i) << " f(x)= " << gsl_vector_get (s->f, i);
		std::cout << "\t" << gsl_strerror (status)<<std::endl;
	}
	for(i=0; i<n; i++)  X[i] = gsl_vector_get (s->x, i);
	gsl_multiroot_fdfsolver_free (s);
	gsl_vector_free(x);
	return status;
}
