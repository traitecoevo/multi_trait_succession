// --------------------------------------------------------
//  IMPLEMENTATION OF CLASSES PARAMS AND STARTEGY
//  Created by Daniel Falster on 2010-06-15.
// --------------------------------------------------------

#include "Strategy.h"
#include <math.h>
#include <iostream>
#include <fstream>
#include "gsl_root_class.h"
#include <iomanip>

sim_params::sim_params() {
	reset();
}

void sim_params::reset(void) {
	// Disturbance interval
	log_mean_disturbance_interval = log10(30.0);
	// Extinction coefficient
	c_ext = 0.5;
	// Canopy shape parameters
	Eta = 12;
	Eta_c = 1.0 - 2 / (1.0 + Eta) + 1 / (1 + 2 * Eta);
	// Ratio leaf area to sapwood area
	theta = 4669;
	// Height - leaf area scaling
	a1 = 5.44;      B1 = 0.306;
	// Leaf area - stem volume scaling
	a2 = 6.67E-5;   B2 = 1.75;
	// Root leaf scaling
	a3 = 0.07;
	// Scaling of leaf turnover with LMA
	a4 = 0.0286;	B4 = 1.71;

	// Ratio bark area : sapwood area
	b = 0.17;
	// Root and bark turnover
	k_b = 0.2; k_r = 1.0;
	// Nitrogen concentrations & leaf photosynthesis
	n_area = 1.87E-3;   // Kg/m2
	c_p1 = 150.36; c_p2 = 0.19;
	// Respiration rates
	c_Rs = 4012;		// Mol CO2 / m3 /yr
	c_Rb = 2 * c_Rs;
	c_Rr = 217;		// Mol CO2 / kg /yr
	c_Rl = 2.1E4;		// mol CO2 / kgN /yr  [=6.66e-4 * (365*24*60*60)]
	// Carbon conversion parameter
	Y =  0.7;
	c_bio = 2.45E-2;	// kg / mol [= 12E-3 / 0.49]
	// MORTALITY RATES
	Pi_0      = 0.25;	 // Survival during dispersal
	c_s0      = 0.1;	 // Parameter for seedling survival
	c_d1      = 0.0065;  // Coefficient for wood density in mortality function
	c_d0      = 0.01 / exp(-c_d1 * 608.0); // Baseline mortality rate
	c_d2      = 5.5;     // Baseline mortality rate
	c_d3      = 20.0;    // Coefficient for dry mass production in mortality function

	// REPRODUCTION
	c_r1 = 1.0;          // Maximum allocation to reproduction
	c_r2 = 50;           // Size range across which individuals mature
	c_acc = 4.0;	 // Accessory cost of reproduction - multiplication factor

	seed_mass = 3.8E-05;
	wood_dens = 608.0;

	// NUMBER OF DEPTHS SAMPLED FOR LIGHT ENVIRONMENT
	ENV_DIM = 30;
}

// Print params file which can be loaded in matlab
void sim_params::print_sim_params(ofstream &file) {
	file << "function p=params" << endl;
	file << setprecision(10) << "p.eta=" << Eta << ";" << endl;
	file << "p.theta=" << theta << ";" << endl;
	file << "p.b=" << b << ";" << endl;
	file << "p.a1=" << a1 << ";" << endl;
	file << "p.B1=" << B1 << ";" << endl;
	file << "p.a2=" << a2 << ";" << endl;
	file << "p.B2=" << B2 << ";" << endl;
	file << "p.a3=" << a3 << ";" << endl;
	file << "p.a4=" << a4 << ";" << endl;
	file << "p.B4=" << B4 << ";" << endl;
	file << "p.n_area=" << n_area << ";" << endl;
	file << "p.c_p1=" << c_p1 << ";" << endl;
	file << "p.c_p2=" << c_p2 << ";" << endl;
	file << "p.c_Rl=" << c_Rl << ";" << endl;
	file << "p.c_Rs=" << c_Rs << ";" << endl;
	file << "p.c_Rb=" << c_Rb << ";" << endl;
	file << "p.c_Rr=" << c_Rr << ";" << endl;
	file << "p.k_b=" << k_b << ";" << endl;
	file << "p.k_r=" << k_r << ";" << endl;
	file << "p.Y=" << Y << ";" << endl;
	file << "p.c_bio=" << c_bio << ";" << endl;
	file << "p.c_acc=" << c_acc << ";" << endl;
	file << "p.c_r1=" << c_r1 << ";" << endl;
	file << "p.c_r2=" << c_r2 << ";" << endl;
	file << "p.c_ext=" << c_ext << ";" << endl;
	file << "p.Pi_0=" << Pi_0 << ";" << endl;
	file << "p.c_d0=" << c_d0 << ";" << endl;
	file << "p.c_d1=" << c_d1 << ";" << endl;
	file << "p.c_d2=" << c_d2 << ";" << endl;
	file << "p.c_d3=" << c_d3 << ";" << endl;
	file << "p.seed_mass=" << seed_mass << ";" << endl;
	file << "p.wood_dens=" << wood_dens << ";" << endl;
	file << "p.c_s0=" << c_s0 << ";" << endl;
	file << "p.log_mean_disturbance_interval=" << log_mean_disturbance_interval << ";" << endl;
}

double sim_params::getValueFromFile(ifstream &file) {
	string s;
	getline(file, s,  '=' );
	getline(file, s, ';' );
	return (atof(s.c_str()));
}

void sim_params::load(string filename) {
	string word;
	ifstream file(filename.c_str()); if (!file) {cerr << "parameters file " << filename << " not opened " << endl; exit(1);}
	file >> word; file >> word;

	Eta = getValueFromFile(file);
	theta = getValueFromFile(file);
	b = getValueFromFile(file);
	a1 = getValueFromFile(file);
	B1 = getValueFromFile(file);
	a2 = getValueFromFile(file);
	B2 = getValueFromFile(file);
	a3 = getValueFromFile(file);
	a4 = getValueFromFile(file);
	B4 = getValueFromFile(file);
	n_area = getValueFromFile(file);
	c_p1 = getValueFromFile(file);
	c_p2 = getValueFromFile(file);
	c_Rl = getValueFromFile(file);
	c_Rs = getValueFromFile(file);
	c_Rb = getValueFromFile(file);
	c_Rr = getValueFromFile(file);
	k_b = getValueFromFile(file);
	k_r = getValueFromFile(file);
	Y = getValueFromFile(file);
	c_bio = getValueFromFile(file);
	c_acc = getValueFromFile(file);
	c_r1 = getValueFromFile(file);
	c_r2 = getValueFromFile(file);
	c_ext = getValueFromFile(file);
	Pi_0 = getValueFromFile(file);
	c_d0 = getValueFromFile(file);
	c_d1 = getValueFromFile(file);
	c_d2 = getValueFromFile(file);
	c_d3 = getValueFromFile(file);
	seed_mass = getValueFromFile(file);
	wood_dens = getValueFromFile(file);
	c_s0 = getValueFromFile(file);
	log_mean_disturbance_interval = getValueFromFile(file);

	file.close();
}


// Constructor
Strategy::Strategy() {
	p = NULL;
	Traits.resize(TRAIT_DIM);
}

// Copy constructor
Strategy::Strategy(const Strategy& s) {
	set_sim_params(s.p);
	Traits.resize(TRAIT_DIM);
	set_traits(&s.Traits.front());
}

// Set parameters for operations
void Strategy::set_sim_params(sim_params* pars) {
	p = pars; // Stores pointer to main params file
	set_depths();
}

void Strategy::set_depths(void) {
	if ((int) depths.size() != (*p).ENV_DIM)
	{
		// Calculate depths in canopy for environment sampling - sampling depths spaced evenly with respect to cumulative leaf area
		depths.resize((*p).ENV_DIM);
		double step = 1.0 / ((*p).ENV_DIM - 1.0);
		for (int i = 0; i < (*p).ENV_DIM; i++)
			depths[i] = pow(1.0 - sqrt(1.0 - step * i), 1.0 / (*p).Eta); // Inverse of cdf from Yokozawa
	}

}



// Set traits for individual by passing 1D vector
void Strategy::set_traits(const double S[]) {
	if (p == NULL) {
		cerr << "sim_params not set for Strategy" << endl; exit(1);
	}
	for (int i = 0; i < TRAIT_DIM; i++) Traits[i] = S[i];

	lma = pow(10, S[0]);
	hmat = pow(10, S[1]);
	rho = (*p).wood_dens;
	total_mass_at_birth = (*p).seed_mass;

	// Warning, in case someone tries to modify seed mass and wood density via trait vector
	const double eps = 1e-4;
	if (fabs(total_mass_at_birth / pow(10, S[3]) - 1) > eps) {
		cerr << "Value for seed mass is set in params object. You should use the same value in trait vector and params object."
		     << total_mass_at_birth << "\t" << pow(10, S[3]) << "\t" << fabs(total_mass_at_birth / pow(10, S[3]) - 1) << endl;
		exit(1);
	}

	if (fabs(rho / pow(10, S[2]) - 1) > eps) {
		cerr << "Value for wood density is set in params object. You should use the same value in trait vector and params object."
		     << rho << "\t" << pow(10, S[2]) << "\t" << fabs(rho / pow(10, S[2]) - 1) << endl;
		exit(1);
	}

	// Check for values outside reasonable bounds
	if (S[0] > 5.0)	  lma = pow(10, 5.0); // Limit max lma to prevent population crashes
	if (S[0] < -10.0) lma = pow(10, -10.0); // Limit min lma to prevent problems with root finding for offspring size
	// Calculate constants used for efficient calculation - must keep this
	set_constants();
	set_depths();
}

// Calculate constants, based on traits
void Strategy::set_constants(void) {
	Eta_c = 1.0 - 2 / (1.0 + (*p).Eta) + 1 / (1 + 2 * (*p).Eta);

	// Calculate trait dependent constants
	double ll = 1.0 / ( (*p).a4 * pow(lma, -(*p).B4) );

	k_l = 1.0 / ll;

	// Solve initial for leaf mass at birth
	gsl_root_class gslRootSolver; // Root solver for solving initial leaf mass
	leaf_mass_at_birth = gslRootSolver.single_deriv(&init_mass_F, &init_mass_dF, &init_mass_FDF, this, 1.0 * total_mass_at_birth, 1E-5, 80, 0);
}

// Returns relative height in canopy for depth slice i - used in constructing light profile
double Strategy::canopy_depth(int i) {
	return depths[i];
}

double Strategy::get_trait(int j) {
	if (j >= TRAIT_DIM) {
		cerr << "request for trait out of bounds in get_trait" << endl; system("pause"); exit(1);
	}
	return (Traits[j]);
}

// Calculates rates of growth, mortality and fecundity for individual with size m. Memory for GMR[5] must be allocated elsewhere
// Env[] contains estimates of light environment  at different depths in the canopy used in calculating production.
void Strategy::Grow_Mort_Repro(vector<double>& GMR, double m, const double Env[], double t) {
	GMR[3]    = Production(Env, m);		// GPP
	GMR[6]    = Respiration(m);				// Maintenance respiration
	GMR[4]    = (*p).Y * (GMR[3] - GMR[6]);	// NPP
	GMR[5]    = Turnover(m);				// Tissue turnover
	double dmdt = GMR[4] - GMR[5];			// NET PRODUCTION

	GMR[0] = (1 - r_alloc(m)) / dTotalMass_dm(m) * max(0.0, dmdt); // GROWTH   - only positive growth allowed
	GMR[1] = mortality(dmdt, m);								// MORTALITY
	GMR[2] = (*p).Pi_0 * (max(0.0, dmdt) * r_alloc(m)) / ((*p).c_acc * total_mass_at_birth); // REPRODUCTION

// Check for NaN in mortality
	if ((!(GMR[1] >= 0) && !(GMR[1] < 0)))
		cout << "Mort " << lma << "\t" << GMR[1] << "\t" << dmdt << "\t" << m << "\t" << LfAr(m) << "\t" << dmdt / LfAr(m) << "\t" << exp((*p).c_d1 * rho + (*p).c_d2 * dmdt / LfAr(m)) << endl;
}

// Calculate photosynthesis by integrating over canopy layers
double Strategy::Production(const double Env[], double m) {
	double A = 0;
	A += 0.5 * (Assim(Env[0]) + Assim(Env[(*p).ENV_DIM - 1]));
	for (int i = 1; i < (*p).ENV_DIM - 1; i++)
		A += Assim(Env[i]);
	A *= 1.0 / ((*p).ENV_DIM - 1);
	return (*p).c_bio * A * LfAr(m);
}

// Returns assimilation rate as a fraction of maximum where E is available light (0-1)
double Strategy::Assim(const double E) {
	return (*p).c_p1 * E / (E + (*p).c_p2);
}

// Returns total maintenance respiration for individual with state m
double Strategy::Respiration(double m) {
	return (*p).c_bio * (  (*p).c_Rl * (*p).n_area * LfAr(m)   +  ((*p).c_Rs + (*p).c_Rb * (*p).b) * MassSapwood(m) / rho + (*p).c_Rr * MassRoots(m));
}

// Returns loss due to tissue turnover for individual with state m
double Strategy::Turnover(double m) {
	return k_l * m  +  (*p).k_b * (*p).b * MassSapwood(m)   +   (*p).k_r * MassRoots(m);
}

// Gives fraction of production allocated to reproduction for individual with state m
double Strategy::r_alloc(double m) {
	return (*p).c_r1 / (1.0 + exp((*p).c_r2 * (1.0 - Height(m) / hmat)));
}

// Returns survival during germination.
double Strategy::germination(double m, const double Env[], double t) {
	double dmdt =  ((*p).Y * (Production(Env, m) - Respiration(m)) - Turnover(m)) / LfAr(m);
	if (dmdt > 0.0)
		return 1.0 / (pow((*p).c_s0 / dmdt, 2.0) + 1.0);
	else
		return 0.0;
}

// Returns total leaf area of individual with state m
double Strategy::LfAr(double m) {
	return m / lma;
}

// Returns height of individual with state m
double Strategy::Height(double m) {
	return  (*p).a1 * pow(LfAr(m), (*p).B1);
}

// Returns height of individual with state m
double Strategy::dHeight_dm(double m) {
	return  (*p).a1 * (*p).B1 * pow(LfAr(m), (*p).B1 - 1) / lma;
}

// Returns total leaf area above focal_height for plant of 'height' and state 'm'
double Strategy::competiton(double focal_height, double height, double m) {
	if (height == 0 || focal_height > height)
		return 0.0;
	else
		return LfAr(m) * pow(1.0 - pow(focal_height / height, (*p).Eta), 2);
}

// Returns total mass of individual with state m
double Strategy::TotalMass(double m) {
	return m + MassRoots(m) + MassStem(m);
}

// Returns stem mass of individual with state m
double Strategy::MassStem(double m) {
	return (1.0 + (*p).b) * MassSapwood(m) + MassHeartwood(m);
}

double Strategy::MassRoots(double m) {
	return (*p).a3 * LfAr(m);
}
double Strategy::MassSapwood(double m) {
	return rho / (*p).theta * (*p).a1 * (*p).Eta_c * pow(LfAr(m), 1 + (*p).B1);
}

double Strategy::MassHeartwood(double m) {
	return rho * (*p).a2 * (*p).Eta_c * pow(LfAr(m), (*p).B2);
}


// Derivative of stem mass with respect to leaf mass.
double Strategy::dMassStem_dm(double m) {
	return rho / (*p).theta * (*p).a1 * (*p).Eta_c * (1.0 + (*p).b) * (1.0 + (*p).B1) * pow(LfAr(m), (*p).B1) / lma
	       +	rho * (*p).a2 * (*p).Eta_c * (*p).B2 * pow(LfAr(m), (*p).B2 - 1) / lma;
}

// Derivative of total mass with respect to leaf mass.
double Strategy::dTotalMass_dm(double m) {
	return 1.0 + (*p).a3 / lma + dMassStem_dm(m);
}

// Returns mortality as sum of density independent, size dependent and density dependent rates
double Strategy::mortality(double DmDt, double m) {
	double M1, M2;
	// Mortality dependent on production
	if (m <= 0)
		M1 = 0;
	else {
		// Net assimilation per leaf area
		M1 = DmDt / LfAr(m);
		if (M1 > 1e5) M1 = 1e5; // Protects against very large leaf areas from low LMA
		if (M1 < -10) M1 = -10.0; // Limits maximum mortality to some very large value

		M1 = (*p).c_d2 * exp(-(*p).c_d3 * M1);
	}
	// Mortality dependent on wood density
	M2 =  (*p).c_d0 * exp(-(*p).c_d1 * rho);

	return M1 + M2;
}

// Return maximum possible leaf mass for individual with these traits
// Used by continuous solver to delete unneeded cohorts when changing trait valuess
double Strategy::mass_max() {
	double f = 0.9999; // Fraction of max size considered too big
	return lma * hmat * ((*p).c_r2 - log((1 - f) / f)) / (*p).c_r2;
}

/*functions  used to solve initial leaf mass for a given seed size. F gives difference b/w total mass
estimated from leaf mass and seed size. dF gives derivative of F, and FDF gives both */
double init_mass_F(double  X, void *sim_params) {
	Strategy * p = ( Strategy * ) sim_params;
	return p->TotalMass(X) - p->offspring_mass();
}

double init_mass_dF(double X, void *sim_params) {
	Strategy * p = ( Strategy * ) sim_params;
	return p->dTotalMass_dm(X) * 1.0;
}

void init_mass_FDF (double X, void * sim_params,  double * f, double * df) {
	Strategy * p = ( Strategy * ) sim_params;
	*f = p->TotalMass(X) - p->offspring_mass();
	*df = p->dTotalMass_dm(X) * 1.0;
}

double Strategy::mass_at_birth(void) {
	return leaf_mass_at_birth;
}

double Strategy::offspring_mass(void) {
	return total_mass_at_birth;
}
