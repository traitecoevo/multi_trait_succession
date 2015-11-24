//*******************************************************************************************/
//  EBT_Metapopulation.h
//  Created by Daniel Falster, 2006-2010
//*******************************************************************************************/

#ifndef __EBT_META_H_
#define __EBT_META_H_

#include "EBT_Base.h"

class EBT_Metapopulation : public EBT_Base
{
public:
	// Constructors
	EBT_Metapopulation();
//	~EBT_Metapopulation();

	// Functions needed by Fitness solver - override virtual functions in FitnesSolverClass (inherited by EBT_Base)
	double mutantFitness(const double Traits[], bool with_competition);
	vector<vector<double> > getPopulationProperties(void);
	vector<double> residentFitness(const double residentTraits[], const double X[], int N, int useDynamicCohorts);
	void printPopulationProfile(const double residentTraits[], const double X[], int N, string DIR, bool printRates);
	void printPopulationProfile(const double residentTraits[], const double X[], int N, string DIR);

	// Configuration
	void setup(sim_params *PhysPars);
	void set_max_cohort_fitness(double eps);
	int print;
	// Used to enforce pre-determined time steps
	void load_time_steps(string filename);

	// Interaction
	vector<double> residentFitness(const double residentTraits[], const double X[], int N, int useDynamicCohorts, bool adaptive_time_step, int print_options, string dir);
	double get_cohort_residentFitness(int spp_n);
	void print_splines(double top, double steps, bool all_RK_steps, string name);

	double maxIterations;
private:
	// Parameters
	double mean_dist, time_end;		// Length simulation
	double max_cohort_fitness;  // Limits contribution of any single cohort to fitness calculations
	ofstream debug_file;



	// Functions
	double mutant_individual_restiming(Strategy* strat, int t_start, int t_end, bool with_competition, int output, int print);
	void empty_thisPatch_store(void);

	// Data storage
	vector<Patch*> thisPatch_store;
	bool thisPatch_store_isSetup;
	vector<double> mut_fit_dy;
	vector<double> default_time_steps; // Used when enforcing pre-determined time steps

};

#endif
