//*******************************************************************************************/
//  Strategy.h
//  Created by Daniel Falster, 2006-2010
//*******************************************************************************************/
#ifndef __STRATEGY_H_
#define __STRATEGY_H_

#include <vector>
#include <iostream>

#define TRAIT_DIM 4

using namespace std;

/*Class used to store parameters.*/
class sim_params {
public:
	sim_params();      // constructor
	//~sim_params();
	void reset(void);
	double log_mean_disturbance_interval;
	double Eta, Eta_c;  			// Canopy shape parameters
	double theta;							// Leaf area per sapwood area
	double a1, B1;   					// Empirical constants for scaling relationship of height - leaf area
	double a2, B2;   					// Empirical constants for scaling relationship of Leaf area - stem volume (unused)
	double a3;   							// Empirical constants: leaf area to root mass
	double a4, B4;   					// Empirical constants for scaling relationship of leaf turnover with LMA
	double b;									// Bark area per sapwood area
	double n_area;						// Leaf nitrogen per area
	double c_p1, c_p2;   			// Leaf productivity parameters  - only used when no N reallocation
	double c_Rs, c_Rb;				// Sapwood and bark respiration rate
	double c_Rr, c_Rl;  			// Root and leaf respiration rate
	double Y;									// Yield = carbon fixed in tissue per carbon assimilated;
	double c_bio;       			// Conversion factor
	double k_b, k_r;    			// Bark and root turnover rates
	double c_ext;       			// Light extinction coefficient
	double c_r1;        			// Proportion production allocated to reproduction
	double c_r2;       				// Size range across which individuals mature
	double c_acc;	    				// Accessory cost of reproduction - multiplication factor
	double Pi_0;							// Survival during dispersal
	double c_s0;							// Parameter for seedling mortality
	double c_d0;        			// Baseline structural mortality rate
	double c_d1;	    				// Coefficient for wood density in mortality function
	double c_d2;        			// Baseline for growth mortality rate
	double c_d3;							// Coefficient for dry mass production in mortality function
	double seed_mass;
	double wood_dens;

	int ENV_DIM;		// Number depths for integration of assimilation
	void print_sim_params(std::ofstream &file); // Print params file which can be loaded in matlab
	void load(string filename);
	double getValueFromFile(std::ifstream &file);
private:
};


class Strategy {
public:
	Strategy();
	Strategy(const Strategy& s);
	// Pass pointer to parameters file
	void set_sim_params(sim_params* pars);
	// Set traits of individual
	void set_traits(const double S[]);
	// Retrieve jth trait
	double get_trait(int j);
	// Calculate growth, mortality and reproduction for individual
	void Grow_Mort_Repro(vector<double>& GMR, double m, const double Env[], double t);
	// Calculate survival during germination
	double germination(double m, const double Env[], double t);
	// Total leaf area of plant mass m
	double LfAr(double m);
	// Total leaf area above focal height for plant of height and leaf mass m
	double competiton(double focal_height, double height, double m);
	// Height of plant with leaf mass m
	double Height(double m);
	// Respiration of plant with leaf mass m
	double Respiration(double m);
	// Turnover of plant with leaf mass m
	double Turnover(double m);
	// Total mass of plant with leaf mass m
	double TotalMass(double m);
	// Change in total mass of plant per unit leaf mass m
	double dTotalMass_dm(double m);
	// Return maximum possible leaf mass for individual with these traits
	double mass_max(void);
	// Return leaf mass at birth
	double mass_at_birth(void);
	// Return total offspring mass at birth
	double offspring_mass(void);

	// Gets relative height for depth slice i - used in constructing light profile
	double canopy_depth(int i);

	double Eta_c;

private:
	sim_params *p;		// Pointer to parameters

	// Relative height in canopy for sampling environment
	void set_depths(void);
	vector<double> depths;

	// Storage for traits
	vector<double> Traits;
	double lma, hmat, rho, leaf_mass_at_birth, total_mass_at_birth;

	// Constants used in calculations, dependent on traits
	void set_constants(void);
	double k_l;
	// Production functions
	double Assim(const double Env);
	double Production(const double Env[], double m);
	double mortality(double dmtdt, double LfAr);
	double r_alloc(double m);

	// Accessory functions - not actually used in solver but may be useful when relating to data
	double MassSapwood(double m);
	double MassHeartwood(double m);
	double MassRoots(double m);
	double MassStem(double m);
	double dHeight_dm(double m);
	double dMassStem_dm(double m);
};

// Wrappers for root finding of initial leaf mass
double init_mass_F(double  X, void *sim_params);
double init_mass_dF(double X, void *sim_params);
void init_mass_FDF (double X, void * sim_params,  double * f, double * df);

#endif
