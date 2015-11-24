//*******************************************************************************************/
//  EBT_BASE.h
//  Created by Daniel Falster, 2006-2010
//*******************************************************************************************/

#ifndef __EBT_BASE_H_
#define __EBT_BASE_H_

#include "FitnessSolver.h"
#include "Utils.h"
#include "Strategy.h"
#include "site.h"
#include "RungeKutta_CashKarp_45.h"

#include "gsl/gsl_interp.h"
#include "gsl/gsl_spline.h"

// Flag for debugging
#define DEBUG 0

#define EBT_MAX_N0_SPECIES 200   // Sets absolute maximum number of species. Currently only used in spp_weighting_for_splitting vector -- making this dynamic would allow you to remove limit.

struct cohort_data{double time_at_birth, mu, lam, coh_R, height, bound, bound_n0, bound_nt, bound_S, bound_R, number_of_splits;};

class Patch
{
public:
	Patch();
	~Patch();
	double age;
	vector<double> disturbance_rate, survival, Integral_of_survival, top_height, LAI;
	vector<gsl_spline*> Env_Spline;
	vector<gsl_interp_accel*> Spline_Accel;
};

class CohortRatesOfChange
{
public:
	CohortRatesOfChange();
	vector<double> d_mu, d_lambda, d_R, d_bound, d_bound_S, d_bound_nt, d_bound_R;
};

class PopulationData {
public:
	PopulationData();
	double Bottom_Coh_Mass, Seed_Rain_in_cumulative, Seed_Rain_out_cumulative;
	double H_top, H_av, Biom, LAI, GPP, NPP, Turnover, Resp, N;
};

class SpeciesData {
public:
	SpeciesData();
	Strategy strat;
	double X_start, X_end;

	vector<PopulationData> popn;
	vector<vector<cohort_data> > cohort;
	vector<CohortRatesOfChange> rates;
	vector<double> d_pi, Seed_Rain_in, Seed_Rain_out;

	stack<int> spare_cohorts;
	list<int> cohort_list;
	int no_cohorts;

	int dim_time, dim_cohort;
	void set_dimensions(int cohort_steps);

	double H_top_cumulative, H_av_cumulative, Biom_cumulative, LAI_cumulative, GPP_cumulative, NPP_cumulative, Turnover_cumulative,  Resp_cumulative, N_cumulative;

// For printing
	ofstream* OutFile;
};



class EBT_Base : public FitnessSolver {
public:
	EBT_Base();

// Configuration
	void update_sim_params(void);
	void set_eps(int which, double eps);
	void displayEps(int which);
	void set_eps_cohort_splt(double eps, double merge_times, vector<double> spp_weighting);
	void set_eps_cohort_splt(double eps,  double merge_times);


	void remove_spp(int pos);

	void set_max_no_cohorts(int howMany);
	void displayEpsCohortSplt(void);

	void set_default_cohort_introduction_times(string filename);
	void set_default_cohort_introduction_times(double end, double small_step,
		double large_step, double mulltiplier);
	void set_default_cohort_introduction_times_fixed(double end, double step);

	void print_cohort_introduction_times(int spp, string filename);

// Functions used for mutant fitness
	double mutant_individual_single_spline(Strategy* strat, int t_index, double S_end,
		bool with_competition, int output, int print);
// Printing
	void make_and_print_spline(string filename);

	site *s;						// Site parameters & functions
	vector<double> d_patch_survival, d_patch_Integral_of_survival;

protected:
// Configuration
	bool isSetup;                   // Flag to check setup before sim start
	bool continuous_seed_flow_flag;
	sim_params * p;                      // Physiological parameters
	vector<double> cohort_intro_default;
	int check_coh_error;

// Array sizes
	int no_res;                    // Number resident strategies

// Clear data
// void clear_all(void);          // Remove all data
	void set_zero(SpeciesData& data); // Clear data
	void clear(void); // Remove all species data

// Resizing arrays
	void add_spp(const double TRAITS[]);        // Allocate memory for new species

// Cohort operations during simulation
	void begin_new_cohort(int spp_no);
	void step_cohorts(double &time, double &dt, bool adaptive);

	void calculate_rates_of_change(int t_index, double time, int index);
	void update_states(int t_curr, int t_next, double delta_T, int type);
	void calculate_EBT_error_for_cohorts(double time);
	void calculate_EBT_error_for_cohorts_short(double time);
	double error_in_cohort_approximation(int n, double size, double width, double lambda, double time, bool relative_error);
	list<int>::iterator merge_cohort_up(int n, int t_index, list<int>::iterator It);
	list<int>::iterator delete_cohort(int n, int t_index, list<int>::iterator It_ext);
	list<int>::iterator collapse_cohort(int n, int t_index, list<int>::iterator It_ext);

// Environment operations
	double calculate_env_at_height(double size, int t_index);
	int make_env_spline(int t_index, int RK_index);
	double get_env_at_height_from_spline(double height, Patch* patch, int RK_index);
	void get_env_at_height_from_spline_vector(double top_height, Strategy* strat, Patch* patch,	int RK_index, vector<double>& result);
// Data storage
	Patch* thisPatch;
	Patch* nextPatch;

	vector<SpeciesData> Spp;
	SpeciesData temp_spp;
	int t_curr, t_next;

// Error controls
	double eps_ode;				// Error control on resident ode
	double eps_env;				// Used in adaptive sampling of environment
	double ODE_MINSTEPSIZE, ODE_MAXSTEPSIZE; // Maximum and minimum step size for ODE
	double eps_mutant_ode;		// Error control on ODE for mutant individual
	double eps_cohort_split;		// Threshold for cohort splitting
	double eps_cohort_merge;		// Threshold for cohort merging
	vector<double> spp_weighting_for_splitting; // Vector of multiplication factors, to apply specific split value to each spp
	int max_no_cohort_steps;
	bool USE_ENV_SPLINE;
// Error calculations
	double max_d_LAI;          // Maximum change in LAI thus far - used for error controls

// Counts number of failed steps in resident stepping & number of simulations
	int ode_fail, min_fail, num_sim;

// Temporary storage used in fitting env_sample spline
	vector<double> env_spline_data_X, env_spline_data_Y, env_dy;
	double env_dh_start;

// Temporary storage for mutant fitness
	Strategy mut;
	vector<double> Y1, Y2;
	vector<vector<double> > mut_dy;

// Temporary variables to save reallocating memory
	vector<double> GMR, GMR2, GMR3, env_sample;

// Runge Kutta class used for stepping
	RungeKutta_CashKarp_45 RK;

// For printing to file
	ofstream* OutFile_thisPatch;
	int Nfiles;

	void print_start(string dir, int option);
	void print_states(int t_index, vector<int> max_cols, int option);
	void print_end(void);
};

#endif
