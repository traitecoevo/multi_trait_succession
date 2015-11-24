//*******************************************************************************************/
// Random number generation, using Gsl library
//
// (C) Daniel Falster, last update 2011.05.16
//*******************************************************************************************/
#ifndef __GSL_RAN_CLASS_H_
#define __GSL_RAN_CLASS_H_

#include <gsl/gsl_rng.h>  // For random number generation
#include <gsl/gsl_randist.h>
#include <vector>

using namespace std;

class gsl_random_numbers
{
public:
	gsl_random_numbers();
//	~gsl_random_numbers();

	double RandomNumber_Uniform_1D(void);
	double RandomNumber_Uniform_1D(double low, double high);
	void RandomNumber_Uniform_2D(double &X1, double &X2, double low1, double high1, double low2, double high2);

	double RandomNumber_Gaussian_1D(void);
	double RandomNumber_Gaussian_1D(double sigma2);
	void RandomNumber_Gaussian_2D(double &X1, double &X2, double sigma2_x, double sigma2_y, double rho);
	vector<double> RandomNumber_Gaussian_nD_noCovar(vector<double> sigma2, int dim);
	unsigned int RandomNumber_Poisson(double lamda);

private:
	gsl_rng * r;
	const gsl_rng_type * T;

};
#endif




