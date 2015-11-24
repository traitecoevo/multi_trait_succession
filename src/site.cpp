#include "Utils.h"
#include "gsl/gsl_cdf.h"
#include "gsl/gsl_interp.h"
#include "gsl/gsl_spline.h"
#include "site.h"
#include <math.h>
#include "gsl/gsl_sf_gamma.h"

// Pars for Weibull distribution
#define psi 2.0

// Pars for Exponential distribution

site::site()                // constructor
{
	// Set default as Weibull distribution
	set_type(1.0);
}

void site::set_type(int T) {
	type=T;
	solve_patchage_dist(30.0); // Default distribution has mean age =30.0 years
}

// Returns point where 0.999 of integral of p(a,t)
double site::solve_patchage_dist(double site_mean) {
	switch(type){

		default:case(0): // Single patch
			return 0.0;

		case(1): // Weibull distribution
			a_mean = site_mean;
				// Solve lam as function of site average patch age
			lam = pow(gsl_sf_gamma(1.0/psi)/psi/a_mean, psi);
				// Solve for density age zero
			p0 = psi*pow(lam, 1.0/psi)/gsl_sf_gamma(1.0/psi);
			return 2.633*a_mean /3.0*4.0;

		case(2): // Exponential distribution
			a_mean = site_mean;
					// Solve lam as function of site average patch age
			p0 = lam = 1.0/site_mean;
			return -log(0.0001)/lam;
	}
	return 0.0;
}

double site::disturbance_rate(double age) {
	switch(type){
		default:case(0): // Single patch
		return 0.0;
		case(1): // Weibull distribution
		return lam*psi*pow(age, psi-1);
		case(2): // Exponential distribution
		return lam;
	}
	return 0.0;
}


double site::patch_age_density_freq(double age) {
	switch(type){
		default:case(0):  // Single patch
		return 1.0;
		case(1):case(2): // Weibull + exponential distribution
		return p0*Pi(age);
	}
	return 0.0;

}


double site::Pi(double age) {
	switch(type){
		default:case(0): // Single patch
		return 1.0;
		case(1): // Weibull distribution
		return exp(-lam*pow(age, psi));
		case(2): // Exponential distribution
		return exp(-lam*age);
	}
	return 0.0;

}

double site::patch_weight(double time_start, double time) {
	switch(type){
		default:case(0):// Single patch
		return 1.0;
		case(1):case(2): // Weibull  / exponential distribution
		return Pi(time)/ Pi(time_start);
	}
	return 0.0;

}
