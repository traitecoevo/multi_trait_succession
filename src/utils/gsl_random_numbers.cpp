
#include "gsl_random_numbers.h"

#include <math.h>
#include <iostream>

gsl_random_numbers::gsl_random_numbers(){
	// Set seed
	gsl_rng_env_setup();
	T = gsl_rng_default; r = gsl_rng_alloc (T);
	gsl_rng_set(r, time(0));
}


// Generates uniform random number between 0-1
double gsl_random_numbers::RandomNumber_Uniform_1D(void) {
	return gsl_rng_uniform(r);
}

// Generates uniform random number in specified  range
double gsl_random_numbers::RandomNumber_Uniform_1D(double low, double high) {
	return gsl_rng_uniform(r) *(high - low) + low;
}

// Generates 2 uniform random numbers in specified range, stored these in X1 & X2
void gsl_random_numbers::RandomNumber_Uniform_2D(double &X1, double &X2, double low1, double high1, double low2, double high2) {
	X1 =RandomNumber_Uniform_1D(low1, high1);
	X2 =RandomNumber_Uniform_1D(low2, high2);
}


// Generates random number from Gaussian distribution
double gsl_random_numbers::RandomNumber_Gaussian_1D(void) {
	return gsl_ran_gaussian (r, 1);
}

// Generates random number from Gaussian distribution w variance sigma2
double gsl_random_numbers::RandomNumber_Gaussian_1D(double sigma2) {
	return gsl_ran_gaussian (r, sqrt(sigma2));
}

// Generates 2 random numbers from bivariate Gaussian distribution
void gsl_random_numbers::RandomNumber_Gaussian_2D(double &X1, double &X2, double sigma2_x, double sigma2_y, double rho) {
	gsl_ran_bivariate_gaussian (r, sqrt(sigma2_x), sqrt(sigma2_y), rho, &X1, &X2);
}

// Generates random number from Gaussian distribution
unsigned int gsl_random_numbers::RandomNumber_Poisson(double lamda) {
	return	gsl_ran_poisson(r, lamda);
}

// Returns a vector of length dim, sampled from an n-D Gaussian distribution, with no covariance. To incorporate covariance, need a more complicated routine - see here http://en.wikipedia.org/wiki/Multivariate_normal_distribution#Correlations_and_independence
// Indexing 0,...n-1.
vector<double> gsl_random_numbers::RandomNumber_Gaussian_nD_noCovar(vector<double> sigma2, int dim) {

	vector<double> x(dim, 0.0);
	for(int i= 0; i<dim; i++)
		x[i] = gsl_ran_gaussian (r, sqrt(sigma2[i]));
	return x;
}
