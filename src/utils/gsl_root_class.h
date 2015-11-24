/////////////////////////////////////////////////
// ROOT solver class header file
//
// (C) Daniel Falster, last update 16/07/2010
//////////////////////////////////////////////////
#ifndef __GSL_ROOT_CLASS_H_
#define __GSL_ROOT_CLASS_H_

#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_matrix.h>

class gsl_root_class {
public:
	gsl_root_class();
	double single_bracket(double (* function) (double X, void *pars),
		void * pars, double low, double up, double rel_tol, int max_iter, bool print);
	double single_deriv(double (* f) (double X, void *pars), double (* df) (double x, void * params ),
		void (* fdf) (double x, void * params,  double * f, double * df), void * pars,
		double xinit, double rel_tol, int max_iter, bool print);
	int multi(int (* function) (const gsl_vector * x, void *params, gsl_vector * f), void *pars,
		double *X, int dim, double rel_tol, int max_iter, bool print);
	int multi_jac(int (* func_f) (const gsl_vector * x, void *params, gsl_vector * f),
		int (* func_df) (const gsl_vector * x, void *params, gsl_matrix * J),
		int (* func_fdf) (const gsl_vector * x, void *params, gsl_vector * f, gsl_matrix * J),
		void *pars, double *X, int dim, double rel_tol, int max_iter, bool print);
private:
};

#endif
