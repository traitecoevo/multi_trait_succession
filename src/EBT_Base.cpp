//*******************************************************************************************/
// Implementation of class EBT_Base_base
// Author : Daniel Falster
// Date   : Feb 2007
//*******************************************************************************************/
#include "EBT_Base.h"

#include <iostream>
#include <gsl/gsl_linalg.h> // For linear algebra in cohort split


// For environment operations

// DEFAULT ERROR CONTROLS

#define MUTANT_ODE_INIT_STEP 0.01 // Initial time-step used when stepping mutant ODE

#define default_eps_ode  1e-3
#define default_eps_env  1e-5
#define default_ODE_MINSTEPSIZE  1e-4   // Minimum possible step size
#define default_ODE_MAXSTEPSIZE  2.0   // Maximum possible step size
#define default_eps_mutant_ode 1e-6
#define ENV_MINSTEP	  0.01

PopulationData::PopulationData() {
	Bottom_Coh_Mass = Seed_Rain_out_cumulative = Seed_Rain_in_cumulative = H_top = H_av = Biom = LAI = GPP = NPP = Turnover = N = 0.0;
}

Patch::Patch() {
	disturbance_rate.resize(6, 0.0);   survival.resize(6, 1.0);  Integral_of_survival.resize(6, 0.0);
	top_height.resize(6, 0.0);  LAI.resize(6, 0.0);
	Env_Spline.resize(6, '\0');
	Spline_Accel.resize(6, '\0');
}

Patch::~Patch() {
	for (int RK = 0; RK < 6; RK++)
		if (Env_Spline[RK] != '\0') {
			gsl_spline_free (Env_Spline[RK]);
			gsl_interp_accel_free (Spline_Accel[RK]);
		}
};


SpeciesData::SpeciesData() {
	cohort_list.clear();
	no_cohorts = 0;
	d_pi.resize(9);
	Seed_Rain_in.resize(9);
	Seed_Rain_out.resize(9);
	dim_time = dim_cohort = 0;
}

void SpeciesData::set_dimensions(int no_cohorts) {
	int i, time_steps = 2;
	cohort_data temp_coh;
	temp_coh.mu = temp_coh.lam = temp_coh.coh_R = temp_coh.height = temp_coh.bound = 0.0;
	temp_coh.bound_n0 = temp_coh.bound_nt = temp_coh.time_at_birth
	                                        = temp_coh.bound_S = temp_coh.bound_R = 0.0;
	temp_coh.number_of_splits = 0;
// Adds or deletes new rows to cohort as necessary
	if (dim_cohort != no_cohorts)
	{
		rates.resize(no_cohorts);
		cohort.resize(no_cohorts, vector<cohort_data>( dim_time, temp_coh));
		for (i = no_cohorts - 1; i >= dim_cohort; i--)
			spare_cohorts.push(i);
		dim_cohort = no_cohorts;
	}
// Resize time
	if (dim_time != time_steps)
	{
		popn.resize(time_steps);
		for (i = 0; i < dim_cohort; i++)
			cohort[i].resize(time_steps, temp_coh);
		dim_time   = time_steps;
	}
}

CohortRatesOfChange::CohortRatesOfChange() {
	d_mu.resize(9, 0.0);
	d_lambda.resize(9, 0.0);
	d_R.resize(9, 0.0);
	d_bound.resize(9, 0.0);
	d_bound_S.resize(9, 0.0);
	d_bound_nt.resize(9, 0.0);
	d_bound_R.resize(9, 0.0);
}

EBT_Base::EBT_Base() {
	isSetup = 0;
	GMR.resize(7);
	GMR2.resize(7);
	GMR3.resize(7);
	env_sample.resize(1);
	d_patch_survival.resize(9, 0.0);
	d_patch_Integral_of_survival.resize(9, 0.0);

	// Memory for mutant
	Y1.resize(4); Y2.resize(4);
	mut_dy.insert(mut_dy.begin(), 4, vector<double>(9, 0));

// Allocate memory for shading algorithm
	env_spline_data_X.resize(100, 0); env_spline_data_Y.resize(100, 0);
	env_dy.resize(9, 0);
	env_dh_start = 0.1;

	// ALLOCATE MEMORY FOR PATCH
	thisPatch = new Patch(); // Allocate_new_step_for_thisPatch();

	// Default error controls
	eps_ode  = default_eps_ode;
	eps_mutant_ode = default_eps_mutant_ode;
	eps_env  =  default_eps_env;
	ODE_MINSTEPSIZE = default_ODE_MINSTEPSIZE;
	ODE_MAXSTEPSIZE = default_ODE_MAXSTEPSIZE;
	set_max_no_cohorts(2000);
	USE_ENV_SPLINE = 1;

	check_coh_error = 0;
	no_res = 0;
	continuous_seed_flow_flag = 1;

	spp_weighting_for_splitting.resize(EBT_MAX_N0_SPECIES, 1.0);
}

// Function setting max number of cohorts
void EBT_Base::set_max_no_cohorts(int howMany) {
	max_no_cohort_steps = howMany;
}


// Function allowing user to vary critical controls on solver accuracy
void EBT_Base::set_eps(int which, double eps) {
	switch (which)
	{
	case (1): eps_ode  = eps; break; // Error control on resident ode
	case (2): max_no_cohort_steps = (int) eps; break; // Maximum number of steps for cohorts
	case (4): eps_mutant_ode = eps; break; // Error control on ODE for mutant individual
	case (5): eps_env  = eps; break; // Error control Used in adaptive sampling of environment
	case (6): ODE_MINSTEPSIZE  = eps; break; // Minimum step size in resident ODE
	case (7): USE_ENV_SPLINE = (bool) eps;
	default: break;
	}
	displayEps(which);
}

void EBT_Base::displayEps(int which) {
	switch (which)
	{
	case (1): cout << "EBT EPS:\t" << which << "\t" << eps_ode << endl; break; // Error control on resident ode
	case (2): cout << "EBT EPS:\t" << which << "\t" << max_no_cohort_steps << endl; break; // Maximum number of steps for cohorts
	case (4): cout << "EBT EPS:\t" << which << "\t" << eps_mutant_ode << endl; break; // Error control on ODE for mutant individual
	case (5): cout << "EBT EPS:\t" << which << "\t" << eps_env << endl; break; // Error control Used in adaptive sampling of environment
	case (6): cout << "EBT EPS:\t" << which << "\t" << ODE_MINSTEPSIZE << endl; break; // Minimum step size in resident ODE
	case (7): cout << "EBT EPS:\t" << which << "\t" << USE_ENV_SPLINE << endl; break;
	default: break;
	}

}


void EBT_Base::set_eps_cohort_splt(double eps, double merge_times) {
	vector<double> tempVec(50, 1.0);
	set_eps_cohort_splt(eps, merge_times, tempVec);
}

void EBT_Base::set_eps_cohort_splt(double eps, double merge_times, vector<double> spp_weighting) {
	for (int i = 0; i < (int) spp_weighting.size(); i++)
		spp_weighting_for_splitting[i] = spp_weighting[i];
	eps_cohort_split  = eps;
	eps_cohort_merge  = eps / merge_times;
	displayEpsCohortSplt();
}

void EBT_Base::displayEpsCohortSplt(void) {
	cout << "EBT split:\t" << eps_cohort_split << "\t" << eps_cohort_merge << "\t";
	for (int i = 0; i < 10; i++) cout << spp_weighting_for_splitting[i] << " ";
	cout << endl;
}

void EBT_Base::update_sim_params(void) {
	env_sample.resize((*p).ENV_DIM);
	mut.set_sim_params(p);
	temp_spp.strat.set_sim_params(p);
	for (int n = 0; n < no_res; n++)
		Spp[n].strat.set_sim_params(p);
}

void EBT_Base::clear(void) {
	while (no_res > 0)
		remove_spp(no_res - 1);  // Remove excess strategies
}

// Calculate times at which bottom cohorts are internalised - fixed step
void EBT_Base::set_default_cohort_introduction_times_fixed(double end, double step) {
	cohort_intro_default.clear();
	double dt = step, time = 0;
	int count = 0;
	while (time <= end)
	{
		count = 0;
		while (count < 100 && time <= end + dt)
		{
			cohort_intro_default.push_back(time);
			time += dt;
			count++;
		}
		dt *= 10;
	}
	cohort_intro_default.push_back(time + dt);
}

// Calculate times at which bottom cohorts are internalised - assuming 1:1 scaling between cohort spacing and ptach age
// Lower limit capped at small_step; upper_limit capped at large step;
// Good default values:	double small_step = 1E-5, large_step=1.0, multiplier = 0.07;
void EBT_Base::set_default_cohort_introduction_times(double end, double small_step, double large_step, double multiplier) {
	cohort_intro_default.clear();
	double dt = 0, time = 0;
	while (time <= end )
	{
		cohort_intro_default.push_back(time);
		dt = pow(2.0, floor(log2(time * multiplier)));
		time += max(min(dt, large_step), small_step);
	}
	cohort_intro_default.push_back(time + dt);
}


// Load cohort introduction times from file
void EBT_Base::set_default_cohort_introduction_times(string filename) {
	cout << "Loading cohort intro times from file " << cohort_intro_default.size() << " -> ";
	ifstream file(filename.c_str());
	if (!file) {cerr << "parameters file " << filename << " not opened " << endl; hold(); return;}
	cohort_intro_default.clear();
	double time;
	while (file >> time)
		cohort_intro_default.push_back(time);
	file.close();
	cout << cohort_intro_default.size() << "\tsuccess" << endl;
}

// Print cohort introduction times to file from file
void EBT_Base::print_cohort_introduction_times(int n, string filename) {
	ofstream file(filename.c_str());
	if (!file) {cerr << "parameters file " << filename << " not opened " << endl; hold(); return;}

	list<int>::iterator It = Spp[n].cohort_list.end();
	while (It != Spp[n].cohort_list.begin()) {
		It--;
		file << setprecision(15) << Spp[n].cohort[*It][t_curr].time_at_birth << "\t";
	}
	file.close();
}

// set all entries for population to zero
void EBT_Base::set_zero(SpeciesData& data) {	int t;
	cohort_data* Coh;
	list<int>::iterator It;
	t = 0;
	while (t <= 1) {
		// Clear population data
		data.popn[t].Seed_Rain_out_cumulative = data.popn[t].Seed_Rain_in_cumulative
		                                        = data.popn[t].Bottom_Coh_Mass =  0.0;
		for (It = data.cohort_list.begin(); It != data.cohort_list.end(); It++)
		{
			Coh = &data.cohort[*It][t];
			Coh->mu = Coh->lam = Coh->coh_R = Coh->bound = Coh->bound_S
			                                  = Coh->bound_nt = Coh->bound_R = Coh->time_at_birth = 0.0;
		}
		data.popn[t].H_top = data.popn[t].H_av = data.popn[t].Biom = data.popn[t].LAI
		                     = data.popn[t].GPP = data.popn[t].NPP = data.popn[t].Resp = data.popn[t].Turnover = data.popn[t].N = 0;
		t++;
	}
	data.cohort_list.clear();
	data.no_cohorts = 0;

	data.H_top_cumulative = data.H_av_cumulative = data.Biom_cumulative = data.LAI_cumulative
	                        = data.GPP_cumulative = data.NPP_cumulative = data.Resp_cumulative = data.Turnover_cumulative = data.N_cumulative = 0.0;

	// Make stack of available cohorts
	while (!data.spare_cohorts.empty())
		data.spare_cohorts.pop();
	for (t = data.dim_cohort - 1; t >= 0; t--)
		data.spare_cohorts.push(t);
}

void EBT_Base::add_spp(const double TRAITS[]) {
	temp_spp.strat.set_traits(TRAITS);
	Spp.insert(Spp.end(), temp_spp);  // Note temp_spp should already be set to zero
	no_res = Spp.size();
}

void EBT_Base::remove_spp(int pos) {
	int i;
	vector<SpeciesData>::iterator It = Spp.begin();
	if (pos < 0 || pos >= no_res) {cout << "Error removing spp - bad number " << pos << endl; hold();}
	else
	{
		for (i = 0; i < pos; i++) It++;
		Spp.erase(It);
	}
	no_res = Spp.size();
}


// Simulates mutant with single light environment
double EBT_Base::mutant_individual_single_spline(Strategy* strat, int t_index, double S_end, bool with_competition, int output, int print) {
#if (DEBUG >=4)
	cout << "\n\t\t\t\tENTER EBT_Base::mutant_individual" << endl;;
#endif

	double errmax, time = 0, dt = MUTANT_ODE_INIT_STEP;
	int i, success, K = 0, MAXSTEPS = 1000, rk_index;
	// Initial values
	Y1[0] = strat->mass_at_birth(); // Size
	Y1[1] = 0;				 // Log Survival
	Y1[2] = 0;				 // Seeds - net

// Survival during germination
// Set env_sample vector to maximum if no competition [MUTANT IN VIRGIN env_sample]
	if (!with_competition)
		for (i = 0; i < (*p).ENV_DIM; i++)  env_sample[i] = 1.0;
	else if (with_competition) get_env_at_height_from_spline_vector(strat->Height(Y1[0]), strat, thisPatch, 1, env_sample);
	double Pi_1 = strat->germination(Y1[0], &env_sample.front(), time);
	if (Pi_1 == 0)
		return 0.0;
	else
		Y1[1] = log(Pi_1);
	for (i = 0; i < 3; i++)  Y2[i] = Y1[i];

	// Simulate for patch length
	while (exp(Y1[1]) > S_end && K < MAXSTEPS)
	{
		K++;
		if (print)
			cout << time << " " << flush;
#if (DEBUG >=4)
		cout << "\n\t\t\t\t\t" << time << " - " << flush;
#endif
		success = 0;
		// Take step
		while (!success)
		{
#if (DEBUG >=4)
			cout << dt << " " << flush;
#endif
			for (rk_index = 1; rk_index <= 6; rk_index++) {
				// First step based on Y1, but Y2 should be equal to Y1 at start of each step
				if (with_competition)
					get_env_at_height_from_spline_vector(strat->Height(Y2[0]), strat, thisPatch, 1, env_sample);
				strat->Grow_Mort_Repro(GMR, Y2[0], &env_sample.front(), time);
				mut_dy[0][rk_index] = GMR[0];
				mut_dy[1][rk_index] = -GMR[1];
				mut_dy[2][rk_index] = exp(Y2[1]) * GMR[2];
				for (i = 0; i < 3; i++)  Y2[i] = RK.return_next_y(Y1[i], dt, mut_dy[i], rk_index);
			}
			// CALULATE ERROR
			errmax = 0;
			for (i = 0; i < 3; i++) errmax = max(errmax, mut_dy[i][0]);
			errmax /= eps_mutant_ode;

			if ((errmax < 1.0) || (dt - ODE_MINSTEPSIZE < 1e-8)) { // Successful step or reached minimal allowable step
				success = 1;
				time += dt;
				if (errmax < 1.0) // Increase step size
					dt = RK.grow_step(dt, errmax);
				else
					dt = ODE_MINSTEPSIZE;
				for (i = 0; i < 3; i++)  Y1[i] = Y2[i];
			}
			else			// Failed step - decrease step size
			{
				dt = max(ODE_MINSTEPSIZE, RK.shrink_step(dt, errmax));
				for (i = 0; i < 3; i++)  Y2[i] = Y1[i]; // Reset Y2 back to previous step
			}
		}
		if (print) cout << exp(Y1[1]) << "\t" << Y1[2] << "\t" << Y1[0] << "\t" << thisPatch->top_height[0] << endl;
	}
#if (DEBUG >=4)
	cout << "\n\t\t\t\tEXIT mutant_individual" << endl;
#endif
	if (K > 10000) {cerr << "Error, too many steps in mutant_individual" << endl; hold();}

	switch (output)
	{
	case (0): return Y1[0]; // Size
	case (1): return exp(Y1[1]); // Survival
	case (2): return Y1[2]; // Seeds - net
	}
	return 0;
}


void EBT_Base::step_cohorts(double &time, double &dt, bool adaptive) {

	list<int>::iterator It;
	int i, n, success = 0;
	double weight, errmax;
	nextPatch = new Patch();

	// Make env spline
	make_env_spline(t_curr, 1);
	// Estimate error in cohort approximations:
	// Use ver 1 - for continuous model; 2 - for metapopulation model
	switch (check_coh_error) {
case (0): default: break;
	case (1): calculate_EBT_error_for_cohorts_short(time); break;
	case (2): calculate_EBT_error_for_cohorts(time); break;
	}
	calculate_rates_of_change(t_curr, RK.return_x_after_step(time, dt, 0), 1);

	// Step population with embedded 5th order Cash-Karp Runge Kutta method
	while (!success) {
#if (DEBUG >=3)
		cout << dt << " " << flush;
#endif
		// STEP WITH RUNGE-KUTTA DRIVER
		update_states(t_curr, t_next, dt, 1);	    // RK STEP 1
		for (i = 2; i <= 6; i++) { // RK STEPS 2-6; last step calculates error estimates  saved in der[0]
			make_env_spline(t_next, i);
			calculate_rates_of_change(t_next, RK.return_x_after_step(time, dt, i - 1), i);
			update_states(t_curr, t_next, dt, i);
		}
		// CALULATE ERROR - weighted by fraction of total leaf area contributed by each cohort to total LAI
		errmax = 0.0;
		for (n = 0; n < no_res; n++) {
			for (It = Spp[n].cohort_list.begin(); It != Spp[n].cohort_list.end(); It++) {
				weight = Spp[n].strat.LfAr(Spp[n].cohort[*It][t_curr].mu)
				         * Spp[n].cohort[*It][t_curr].lam / max(thisPatch->LAI[0], 1.0E-4);

				if (It != Spp[n].cohort_list.begin())
					errmax = max(errmax, Spp[n].rates[*It].d_mu[0] * weight);
				errmax = max(errmax, Spp[n].rates[*It].d_lambda[0] * weight);
				errmax = max(errmax, Spp[n].rates[*It].d_bound[0] * weight);
				errmax = max(errmax, Spp[n].rates[*It].d_bound_nt[0] * weight);
				errmax = max(errmax, Spp[n].rates[*It].d_R[0] * weight);
			}
		}
		errmax /= eps_ode;
		if (!adaptive) {
			success = 1;
			time = time + dt;
			dt = ODE_MAXSTEPSIZE;
		}
		else if ((errmax < 1.0) || (dt - ODE_MINSTEPSIZE < 1e-8)) { // Successful step or minimal step size reached
			time = time + dt;
			if (dt - ODE_MINSTEPSIZE < 1e-8)		min_fail++;
			success = 1;
			if (errmax < 1.0) // Increase step size
				dt = min(ODE_MAXSTEPSIZE, RK.grow_step(dt, errmax));
			else
				dt = ODE_MINSTEPSIZE;
		}
		else {			// Failed step -try smaller
			ode_fail++;
			dt = max(ODE_MINSTEPSIZE, RK.shrink_step(dt, errmax));
		}
	}
	// Record maximum rate of change in LAI for use in cohort splitting
	double d_LAI = 0.0;
	for (n = 0; n < no_res; n++) d_LAI += (Spp[n].popn[t_next].LAI - Spp[n].popn[t_curr].LAI);
	max_d_LAI = max(max_d_LAI, fabs(d_LAI / (nextPatch->age - thisPatch->age)));

	// Calculate cumulative averages of aggregate properties
	double a1 = thisPatch->age, a2 =  nextPatch->age;
	double p1 = (*s).patch_age_density_freq(a1),  p2 = (*s).patch_age_density_freq(a2);
	for (n = 0; n < no_res; n++) {
		Spp[n].H_top_cumulative += 0.5 * (a2 - a1) * (p1 * Spp[n].popn[t_curr].H_top +  p2 * Spp[n].popn[t_next].H_top);
		Spp[n].H_av_cumulative += 0.5 * (a2 - a1) * (p1 * Spp[n].popn[t_curr].H_av +  p2 * Spp[n].popn[t_next].H_av);
		Spp[n].Biom_cumulative += 0.5 * (a2 - a1) * (p1 * Spp[n].popn[t_curr].Biom +  p2 * Spp[n].popn[t_next].Biom);
		Spp[n].LAI_cumulative +=  0.5 * (a2 - a1) * (p1 * Spp[n].popn[t_curr].LAI  +  p2 * Spp[n].popn[t_next].LAI);
		Spp[n].N_cumulative +=	0.5 * (a2 - a1) * (p1 * Spp[n].popn[t_curr].N    +  p2 * Spp[n].popn[t_next].N);

		// Integrated using cheap and nasty Euler integration (otherwise need rates of change at next step for NPP< GPP etc)
		Spp[n].GPP_cumulative +=  (a2 - a1) * (p1 * Spp[n].popn[t_curr].GPP);
		Spp[n].NPP_cumulative +=  (a2 - a1) * (p1 * Spp[n].popn[t_curr].NPP);
		Spp[n].Turnover_cumulative += (a2 - a1) * (p1 * Spp[n].popn[t_curr].Turnover);
		Spp[n].Resp_cumulative += (a2 - a1) * (p1 * Spp[n].popn[t_curr].Resp);
	}

	// Update indices
	if (t_curr)	{t_curr = 0; t_next = 1;}
	else		{t_curr = 1; t_next = 0;}
}

// VERSION 1: USED IN METAPOPULATION MOEL, CHECKS FOR SPLITTING and merging
// Starts at bottom and goes up list of cohorts checking for merging and splitting
void EBT_Base::calculate_EBT_error_for_cohorts(double time) {
	double bound_prev, size, lambda, bound, size_next = 0, lambda_next = 0, bound_next = 0, EBT_cohort_error;
	int reached_top;

	if (time == 0.0) return;

	// CHECK CRITERIA FOR SPLITTING & MERGING (step from bottom to top, checks if each cohort could be merged with one above)
	list<int>::iterator It;
	for (int n = 0; n < no_res; n++) {
		// Initial values for boundary cohort
		bound_prev = Spp[n].strat.mass_at_birth();
		reached_top = 0;
		for (It = Spp[n].cohort_list.begin();  It != Spp[n].cohort_list.end(); It++) {
			size = Spp[n].cohort[*It][t_curr].mu;
			lambda = Spp[n].cohort[*It][t_curr].lam;
			bound = Spp[n].cohort[*It][t_curr].bound;
			// Check merging
			It++;
			if (It != Spp[n].cohort_list.end()) {	// Retrieve values for cohort above if not top cohort
				size_next   = Spp[n].cohort[*It][t_curr].mu;
				lambda_next = Spp[n].cohort[*It][t_curr].lam;
				bound_next  = Spp[n].cohort[*It][t_curr].bound;
			}
			else
				reached_top = 1;
			It--;
			// Check for merging
			if (!reached_top &&  // Not at top
			    (Spp[n].cohort[*It][t_curr].number_of_splits == -10.0) && // Currently set to merge
			    (lambda + lambda_next > 0)) // Has some surviving individuals (otherwise error_in_cohort_approximation gives nan )
			{
				EBT_cohort_error = Spp[n].strat.LfAr(error_in_cohort_approximation(n, (lambda * size + lambda_next * size_next) / (lambda + lambda_next),  bound_next - bound_prev, lambda + lambda_next, time, 0)) / max_d_LAI; // ThisPatch->LAI[0];
				if ( EBT_cohort_error > (eps_cohort_merge * spp_weighting_for_splitting[n])) // Cohorts should not be merged
					Spp[n].cohort[*It][t_curr].number_of_splits = 0.0;
			}
			// Check cohort for splitting
			if (reached_top  // Always do for top cohort
			    || (Spp[n].cohort[*It][t_curr].number_of_splits != -10.0)) // Not flagged for merging
			{
				EBT_cohort_error = Spp[n].strat.LfAr(error_in_cohort_approximation(n, size,  bound_next - bound_prev, lambda, time, 0)) / max_d_LAI; // ThisPatch->LAI[0];
				if (EBT_cohort_error > (eps_cohort_split * spp_weighting_for_splitting[n]))
					Spp[n].cohort[*It][t_curr].number_of_splits = max(Spp[n].cohort[*It][t_curr].number_of_splits,
					    floor(pow(EBT_cohort_error / (eps_cohort_split * spp_weighting_for_splitting[n]), 0.3333)));
			}
			// Update values
			bound_prev = bound;
		}
	}
}

// VERSION 2: measures error per cohort only. Number of splits is continuous measure giving
// Starts at bottom and goes up list of cohorts checking for merging and splitting
// Only checks splitting if
void EBT_Base::calculate_EBT_error_for_cohorts_short(double time) {
	double bound_prev, size, lambda, bound, EBT_cohort_error;
	// CHECK CRITERIA FOR SPLITTING & MERGING (step from bottom to top, checks if each cohort could be merged with one above)
	list<int>::iterator It;
	for (int n = 0; n < no_res; n++) {
		// Initial values for boundary cohort
		bound_prev = Spp[n].strat.mass_at_birth();
		for (It = Spp[n].cohort_list.begin();  It != Spp[n].cohort_list.end(); It++) {
			size = Spp[n].cohort[*It][t_curr].mu;
			lambda = Spp[n].cohort[*It][t_curr].lam;
			bound = Spp[n].cohort[*It][t_curr].bound;
			if (thisPatch->LAI[0] > 0.0)
				EBT_cohort_error = error_in_cohort_approximation(n, size,  bound - bound_prev, lambda, time, 1) * Spp[n].strat.LfAr(size) / max(thisPatch->LAI[0], 1E-3);
			else
				EBT_cohort_error = 0.0;
			Spp[n].cohort[*It][t_curr].number_of_splits =  pow(EBT_cohort_error / (eps_cohort_split * spp_weighting_for_splitting[n]), 0.3333);
			// Update values
			bound_prev = bound;
		}
	}
}

// Calculates error in cohort mass accumulation
// Use flag to determine if uses relative or absolute error - caution, because relative error can become very large when x~=0.
double EBT_Base::error_in_cohort_approximation(int n, double size, double width, double lambda, double time, bool relative_error) {
	double x, x_low, x_high, ds  = size * 1e-4; // Interval for calculating second derivative

	get_env_at_height_from_spline_vector(Spp[n].strat.Height(size), &Spp[n].strat, thisPatch, 1, env_sample);
	Spp[n].strat.Grow_Mort_Repro(GMR, size, &env_sample.front(), time);
	x = GMR[0] - size * GMR[1];

	get_env_at_height_from_spline_vector(Spp[n].strat.Height(size - ds), &Spp[n].strat, thisPatch, 1, env_sample);
	Spp[n].strat.Grow_Mort_Repro(GMR, size - ds, &env_sample.front(), time);
	x_low = GMR[0] - (size - ds) * GMR[1];

	get_env_at_height_from_spline_vector(Spp[n].strat.Height(size + ds), &Spp[n].strat, thisPatch, 1, env_sample);
	Spp[n].strat.Grow_Mort_Repro(GMR, size + ds, &env_sample.front(), time);
	x_high = GMR[0] - (size + ds) * GMR[1];

	if(relative_error == 0) {
		return fabs(17.2 * width * width * df2d2x(x_high, x, x_low, ds) * lambda);
	} else 	if (relative_error == 1) {
		return fabs(17.2 * width * width * df2d2x(x_high, x, x_low, ds) / max(x, 1E-6));
	}
	return 0.0;
}


void EBT_Base::calculate_rates_of_change(int t_index, double time, int der_i) {
	double size = 0, ds;
	list<int>::iterator It;
	// thisPatch characteristics - maybe dependent on veg characteristics but not at this stage
	thisPatch->disturbance_rate[der_i - 1] = (*s).disturbance_rate(time);
	// Individual characteristics
	for (int n = 0; n < no_res; n++) {
		Spp[n].Seed_Rain_out[der_i] = 0;
		Spp[n].popn[t_index].NPP = Spp[n].popn[t_index].GPP = Spp[n].popn[t_index].Turnover = Spp[n].popn[t_index].Resp = 0.0;
		// Properties for lower boundary
		It = Spp[n].cohort_list.begin();
		if (Spp[n].cohort[*It][t_index].bound_n0 == 999) { // Set initial density for boundary
			get_env_at_height_from_spline_vector(Spp[n].strat.Height(Spp[n].strat.mass_at_birth()), &Spp[n].strat, thisPatch, der_i, env_sample);
			Spp[n].strat.Grow_Mort_Repro(GMR, Spp[n].strat.Height(Spp[n].strat.mass_at_birth()), &env_sample.front(), time);
			double S_germ = Spp[n].strat.germination(Spp[n].strat.mass_at_birth(), &env_sample.front(), time);
			Spp[n].cohort[*It][t_index].bound_S = (S_germ == 0.0) ? -250 : log(S_germ);
		}
		for (It = Spp[n].cohort_list.begin(); It != Spp[n].cohort_list.end(); It++) {
			// cohort boundaries
			size = Spp[n].cohort[*It][t_index].bound;
			get_env_at_height_from_spline_vector(Spp[n].strat.Height(size), &Spp[n].strat, thisPatch, der_i, env_sample);
			Spp[n].strat.Grow_Mort_Repro(GMR, size, &env_sample.front(), time);
			Spp[n].rates[*It].d_bound[der_i] = GMR[0];
			Spp[n].rates[*It].d_bound_S[der_i] = -GMR[1];
			Spp[n].rates[*It].d_bound_R[der_i] = GMR[2] *
			                                     exp(Spp[n].cohort[*It][t_index].bound_S) * (*s).patch_weight(Spp[n].cohort[*It][t_index].time_at_birth, time);
			ds = size * 1e-5;
			get_env_at_height_from_spline_vector(Spp[n].strat.Height(size - ds), &Spp[n].strat, thisPatch, der_i, env_sample);
			Spp[n].strat.Grow_Mort_Repro(GMR2, size - ds, &env_sample.front(), time);
			Spp[n].rates[*It].d_bound_nt[der_i] = -(GMR[1] + (GMR[0] - GMR2[0]) / ds);
			// All cohorts
			get_env_at_height_from_spline_vector(Spp[n].cohort[*It][t_index].height, &Spp[n].strat, thisPatch, der_i, env_sample);
			Spp[n].strat.Grow_Mort_Repro(GMR, Spp[n].cohort[*It][t_index].mu, &env_sample.front(), time);
			Spp[n].rates[*It].d_lambda[der_i] = -GMR[1] * Spp[n].cohort[*It][t_index].lam;
			Spp[n].rates[*It].d_mu[der_i]  = GMR[0];

			if (It == Spp[n].cohort_list.begin()) // Boundary cohort
				Spp[n].d_pi[der_i]	= GMR[0] * Spp[n].cohort[*It][t_index].lam - GMR[1] * Spp[n].popn[t_index].Bottom_Coh_Mass;

			if ((!(Spp[n].rates[*It].d_lambda[der_i] >= 0) && !(Spp[n].rates[*It].d_lambda[der_i] < 0))) {
				cout << "Err-d_lambda\t" << setprecision(5) << Spp[n].rates[*It].d_lambda[der_i] << "\t" << -GMR[1] << "\t" << GMR[0] << "\t" << GMR[2] << "\t" << Spp[n].cohort[*It][t_index].lam << "\t" << GMR[4] << "\t" << GMR[5] << "\t" << Spp[n].cohort[*It][t_index].mu << endl;
				hold();
			}
			Spp[n].rates[*It].d_R[der_i] = GMR[2] * Spp[n].cohort[*It][t_index].lam * (*s).patch_age_density_freq(time);
			// Population level (sum over all cohorts)
			Spp[n].Seed_Rain_out[der_i]		+= Spp[n].rates[*It].d_R[der_i]; // Seed rain from cohorts
			Spp[n].popn[t_index].GPP		+= GMR[3] * Spp[n].cohort[*It][t_index].lam;
			Spp[n].popn[t_index].NPP		+= GMR[4] * Spp[n].cohort[*It][t_index].lam;
			Spp[n].popn[t_index].Turnover	+= GMR[5] * Spp[n].cohort[*It][t_index].lam;
			Spp[n].popn[t_index].Resp		+= GMR[6] * Spp[n].cohort[*It][t_index].lam;
		}
		// Add seed inflow for boundary cohort
		It = Spp[n].cohort_list.begin();
		get_env_at_height_from_spline_vector(Spp[n].strat.Height(Spp[n].strat.mass_at_birth()), &Spp[n].strat, thisPatch, der_i, env_sample);
		Spp[n].strat.Grow_Mort_Repro(GMR, Spp[n].strat.mass_at_birth(), &env_sample.front(), time);

		if (!continuous_seed_flow_flag)	Spp[n].Seed_Rain_in[der_i] = Spp[n].X_start;
		else		Spp[n].Seed_Rain_in[der_i] = Spp[n].Seed_Rain_out[der_i];
		Spp[n].Seed_Rain_in[der_i] *=  Spp[n].strat.germination(Spp[n].strat.mass_at_birth(), &env_sample.front(), time);
		Spp[n].rates[*It].d_lambda[der_i] += Spp[n].Seed_Rain_in[der_i];

		if (Spp[n].cohort[*It][t_index].bound_n0 == 999) // Set initial density for boundary
			Spp[n].cohort[*It][t_index].bound_n0 = (GMR[0] == 0.0) ? 0.0 : Spp[n].Seed_Rain_in[der_i] / GMR[0];
	}
}

void EBT_Base::update_states(int t_curr, int t_next, double delta_T, int RK_index) {
	list<int>::iterator It, It2;
	int n;
	// ThisPatch state - need to implement properly
	switch (RK_index) {
	case (6): nextPatch->survival[0] = (*s).Pi(thisPatch->age + delta_T);
		nextPatch->age = RK.return_x_after_step(thisPatch->age, delta_T, RK_index);
		break;
	default: thisPatch->survival[RK_index] = (*s).Pi( RK.return_x_after_step(thisPatch->age, delta_T, RK_index)); break;
	}
	thisPatch->Integral_of_survival[RK_index - 1] = 0.0;

	// Individuals states
	for (n = 0; n < no_res; n++) {
		for (It = Spp[n].cohort_list.begin(); It != Spp[n].cohort_list.end(); It++) {
			// boundary intervals
			Spp[n].cohort[*It][t_next].bound = max(Spp[n].strat.mass_at_birth(),
			                                       RK.return_next_y(Spp[n].cohort[*It][t_curr].bound,
			                                           delta_T, Spp[n].rates[*It].d_bound, RK_index));
			Spp[n].cohort[*It][t_next].bound_S = min(Spp[n].cohort[*It][t_curr].bound_S,
			                                     RK.return_next_y(Spp[n].cohort[*It][t_curr].bound_S,
			                                         delta_T, Spp[n].rates[*It].d_bound_S, RK_index));
			Spp[n].cohort[*It][t_next].bound_nt = RK.return_next_y(Spp[n].cohort[*It][t_curr].bound_nt,
			                                      delta_T, Spp[n].rates[*It].d_bound_nt, RK_index);
			Spp[n].cohort[*It][t_next].bound_R = RK.return_next_y(Spp[n].cohort[*It][t_curr].bound_R,
			                                     delta_T, Spp[n].rates[*It].d_bound_R, RK_index);
			// Check for error in boundary R - arises from error in survival (becomes positive due on 4th step of Cash-Karp step)
			if (!(Spp[n].cohort[*It][t_next].bound_R >= 0) && !(Spp[n].cohort[*It][t_next].bound_R < 0))
				cout << "Bound_R err" << Spp[n].cohort[*It][t_curr].bound_R << "\t" << Spp[n].rates[*It].d_bound_R[RK_index] << "\t" << Spp[n].cohort[*It][t_next].bound_R << "\t" << Spp[n].cohort[*It][t_curr].bound << "\t" << endl;
			Spp[n].cohort[*It][t_next].bound_n0 = Spp[n].cohort[*It][t_curr].bound_n0;
			// Cohorts
			Spp[n].cohort[*It][t_next].lam = max(0.0,
			                                     RK.return_next_y(Spp[n].cohort[*It][t_curr].lam,
			                                         delta_T, Spp[n].rates[*It].d_lambda, RK_index));
			if ( (!(Spp[n].cohort[*It][t_next].lam >= 0) && !(Spp[n].cohort[*It][t_next].lam < 0)))
				Spp[n].cohort[*It][t_next].lam = 0;
			if (It == Spp[n].cohort_list.begin()) { // Bottom cohort
				Spp[n].popn[t_next].Bottom_Coh_Mass = max(0.0, RK.return_next_y(Spp[n].popn[t_curr].Bottom_Coh_Mass, delta_T, Spp[n].d_pi, RK_index));
				if (Spp[n].cohort[*It][t_next].lam == 0.0) {
					Spp[n].cohort[*It][t_next].mu = Spp[n].strat.mass_at_birth();
					Spp[n].cohort[*It][t_next].lam = 0; Spp[n].popn[t_next].Bottom_Coh_Mass = 0;
				}
				else {
					Spp[n].cohort[*It][t_next].mu = Spp[n].strat.mass_at_birth()
					                                + Spp[n].popn[t_next ].Bottom_Coh_Mass / Spp[n].cohort[*It][t_next].lam;
					if (Spp[n].cohort[*It][t_next].mu > Spp[n].cohort[*It][t_next].bound) {
						Spp[n].cohort[*It][t_next].mu = Spp[n].cohort[*It][t_next].bound;
						Spp[n].popn[t_next].Bottom_Coh_Mass = (Spp[n].cohort[*It][t_next].mu
						                                       - Spp[n].strat.mass_at_birth()) * Spp[n].cohort[*It][t_next].lam;
					}
				}
			}
			else {  // Internal cohort
				// Check mu > m0 & less than boundary
				Spp[n].cohort[*It][t_next].mu = min(max(Spp[n].strat.mass_at_birth(),
				                                        RK.return_next_y(Spp[n].cohort[*It][t_curr].mu,
				                                            delta_T, Spp[n].rates[*It].d_mu, RK_index)), Spp[n].cohort[*It][t_next].bound);
			}
			Spp[n].cohort[*It][t_next].coh_R	= RK.return_next_y(Spp[n].cohort[*It][t_curr].coh_R,	 delta_T, Spp[n].rates[*It].d_R, RK_index);
			Spp[n].cohort[*It][t_next].height = Spp[n].strat.Height(Spp[n].cohort[*It][t_next].mu);
			Spp[n].cohort[*It][t_next].time_at_birth  = Spp[n].cohort[*It][t_curr].time_at_birth;
			Spp[n].cohort[*It][t_next].number_of_splits = Spp[n].cohort[*It][t_curr].number_of_splits;
		}
		// Seed output / input
		Spp[n].popn[t_next].Seed_Rain_out_cumulative = RK.return_next_y(Spp[n].popn[t_curr].Seed_Rain_out_cumulative, delta_T, Spp[n].Seed_Rain_out, RK_index);
		Spp[n].popn[t_next].Seed_Rain_in_cumulative = RK.return_next_y(Spp[n].popn[t_curr].Seed_Rain_in_cumulative,  delta_T, Spp[n].Seed_Rain_in,  RK_index);

		// Total biomass & leaf area  -not this is just a sum across cohorts
		Spp[n].popn[t_next].H_av  = Spp[n].popn[t_next].Biom = Spp[n].popn[t_next].LAI = Spp[n].popn[t_next].N = 0.0;
		for (It = Spp[n].cohort_list.begin(); It != Spp[n].cohort_list.end(); It++) {
			Spp[n].popn[t_next].H_av  += Spp[n].strat.LfAr(Spp[n].cohort[*It][t_next].mu) *
			                             Spp[n].cohort[*It][t_next].lam * Spp[n].cohort[*It][t_next].height;
			Spp[n].popn[t_next].Biom += Spp[n].strat.TotalMass(Spp[n].cohort[*It][t_next].mu) * Spp[n].cohort[*It][t_next].lam;
			Spp[n].popn[t_next].LAI  += Spp[n].strat.LfAr(Spp[n].cohort[*It][t_next].mu)   * Spp[n].cohort[*It][t_next].lam;
			Spp[n].popn[t_next].N    += Spp[n].cohort[*It][t_next].lam;
		}
		Spp[n].popn[t_next].H_av  = Spp[n].strat.Eta_c * Spp[n].popn[t_next].H_av  / Spp[n].popn[t_next].LAI;
		Spp[n].popn[t_next].H_top = Spp[n].strat.Height(Spp[n].cohort[Spp[n].cohort_list.back()][t_next].bound);
	}
}

void EBT_Base::begin_new_cohort(int n) {
	int k;
	if (Spp[n].no_cohorts != 0)
	{
		k = Spp[n].cohort_list.front(); // get index for current boundary
		if (Spp[n].cohort[k][t_curr].lam > 0 &&  Spp[n].popn[t_curr].Bottom_Coh_Mass >= 0)
			Spp[n].cohort[k][t_curr].mu = Spp[n].strat.mass_at_birth() + Spp[n].popn[t_curr].Bottom_Coh_Mass / Spp[n].cohort[k][t_curr].lam;
		else
			Spp[n].cohort[k][t_curr].mu = Spp[n].strat.mass_at_birth();
		Spp[n].cohort[k][t_curr].height = Spp[n].strat.Height(Spp[n].cohort[k][t_curr].mu);
	}
	// Insert new cohort
	k = Spp[n].spare_cohorts.top();    // Get new reference value from stack
	Spp[n].spare_cohorts.pop();        // Remove from stack
	Spp[n].cohort_list.insert(Spp[n].cohort_list.begin(), k);
	// Initialise cohort
	Spp[n].cohort[k][t_curr].lam   = 0;
	Spp[n].cohort[k][t_curr].mu	  = Spp[n].strat.mass_at_birth();
	Spp[n].cohort[k][t_curr].coh_R = 0.0;
	Spp[n].cohort[k][t_curr].height = Spp[n].strat.Height(Spp[n].cohort[k][t_curr].mu);
	Spp[n].cohort[k][t_curr].bound = Spp[n].strat.mass_at_birth();

	Spp[n].cohort[k][t_curr].time_at_birth = thisPatch->age;
	Spp[n].cohort[k][t_curr].number_of_splits = -10.0;

	Spp[n].cohort[k][t_curr].bound_n0 = 999;
	Spp[n].cohort[k][t_curr].bound_S = 0.0;

	Spp[n].cohort[k][t_curr].bound_nt = 0.0;
	Spp[n].cohort[k][t_curr].bound_R = 0.0;

	Spp[n].popn[t_curr].Bottom_Coh_Mass = 0;
	Spp[n].no_cohorts++;
}

// Merges cohort with cohort above - returns pointer to cohort above
list<int>::iterator EBT_Base::merge_cohort_up(int n, int t_index, list<int>::iterator It_ext) {
	if (It_ext == Spp[n].cohort_list.end()) return It_ext;

	list<int>::iterator It = It_ext;
	double mu = Spp[n].cohort[*It][t_index].mu;
	double lam_new, lam = Spp[n].cohort[*It][t_index].lam;
	double R1 = Spp[n].cohort[*It][t_index].coh_R;

	It++;
	if (It != Spp[n].cohort_list.end())
	{
		lam_new = Spp[n].cohort[*It][t_index].lam + lam;
		if (lam_new > 0)
			Spp[n].cohort[*It][t_index].mu = (Spp[n].cohort[*It][t_index].mu * Spp[n].cohort[*It][t_index].lam + mu * lam) / lam_new;
		else
			Spp[n].cohort[*It][t_index].mu = 0.5 * (Spp[n].cohort[*It][t_index].mu + mu);
		Spp[n].cohort[*It][t_index].lam = lam_new;
		Spp[n].cohort[*It][t_index].coh_R = Spp[n].cohort[*It][t_index].coh_R + R1;
		It--;
		// Add cohort number back onto stack of spare cohorts
		Spp[n].spare_cohorts.push(*It);
		Spp[n].no_cohorts--;
		return Spp[n].cohort_list.erase(It);
	}
	else
	{
		return It;
	}
}

// Deletes cohort, adding individuals to next cohort below and ignoring conservation of mass
// After operation points to next cohort below
list<int>::iterator EBT_Base::collapse_cohort(int n, int t_index, list<int>::iterator It_ext) {
	if (It_ext != Spp[n].cohort_list.begin()) // If not bottom cohort
	{
		list<int>::iterator It = It_ext;
		It--;
		Spp[n].cohort[*It][t_index].lam +=  Spp[n].cohort[*It_ext][t_index].lam;
		// Add cohort number back onto stack of spare cohorts
		Spp[n].spare_cohorts.push(*It_ext);
		Spp[n].cohort_list.erase(It_ext);
		Spp[n].no_cohorts--;
		return It;
	}
	return It_ext;
}

// Deletes cohort and point to next cohort above - deleted individuals not added to any other cohorts
list<int>::iterator EBT_Base::delete_cohort(int n, int t_index, list<int>::iterator It_ext) {
	if (It_ext == Spp[n].cohort_list.end())	// Do nothing if already pointing at end
		return It_ext;

	list<int>::iterator It = It_ext;
	It++;
	// Add cohort number back onto stack of spare cohorts
	Spp[n].spare_cohorts.push(*It_ext);
	Spp[n].cohort_list.erase(It_ext);
	Spp[n].no_cohorts--;
	return It;
}

// Retrieve env_sample data
double EBT_Base::get_env_at_height_from_spline(double height, Patch* patch, int RK_index) {
	if (RK_index < 1 || RK_index > 6) {cout << "Bad index in get_env_spline " << RK_index << endl; abort();}
	if (height > patch->top_height[RK_index - 1])
		return 1.0;
	else
		return min(1.0, max(0.0, gsl_spline_eval(patch->Env_Spline[RK_index - 1], height, patch->Spline_Accel[RK_index - 1])));
}

void EBT_Base::get_env_at_height_from_spline_vector(double top_height, Strategy* strat, Patch* patch, int RK_index, vector<double>& result) {
	if (USE_ENV_SPLINE) {
		for (int i = 0; i < (*p).ENV_DIM; i++)
			result[i] = get_env_at_height_from_spline(top_height * ((*strat).canopy_depth(i)), patch, RK_index);
	}
	else {
		for (int i = 0; i < (*p).ENV_DIM; i++)
			result[i] = calculate_env_at_height(top_height * ((*strat).canopy_depth(i)), (RK_index == 1) ? t_curr : t_next);
	}

}

void EBT_Base::make_and_print_spline(string filename) {
	int i, tot = make_env_spline(t_curr, 1);
	cout << thisPatch->top_height[0] << endl;
	env_dh_start = 1.0;
	ofstream file(filename.c_str());
	if (!file) {cerr << " file " << filename << " not opened " << endl; hold(); return;}

	for (i = 0; i < tot; i++)
		file << setprecision(15) << env_spline_data_X[i] << "\t" << env_spline_data_Y[i] << endl;
	file.close();
}

/* Set up spline interpolation of Environment. Integrates over environment and uses step sizes determined by
adaptive algorithm to fit spline. */
int EBT_Base::make_env_spline(int t_index, int RK_index) {
	if (RK_index < 1 || RK_index > 6) {cout << "Bad index in make_env_spline " << RK_index << endl; abort();}

	int INDEX = RK_index - 1;
	// Allocate memory and set initial values. Note that env_dy is allocated by EBT_Base
	int i = 0, j, success, flag = 1;
	double dh, h = 0.0, y1 = 0.0, y2 = 0, errmax = 0;

	// Determine top height
	thisPatch->top_height[INDEX] = 1e-5;
	for (int k = 0; k < no_res; k++)
		if (Spp[k].no_cohorts > 0)
			thisPatch->top_height[INDEX] = max( thisPatch->top_height[INDEX], Spp[k].strat.Height(Spp[k].cohort[Spp[k].cohort_list.back()][t_index].bound));

	if (USE_ENV_SPLINE) // If using spline interpolation
	{
		// Use initial step size from previous spline fit
		if (thisPatch->age == 0.0 && INDEX == 0) dh = 0.1;
		else dh = max(min(env_dh_start, thisPatch->top_height[INDEX]), ENV_MINSTEP);
		if (no_res == 0)
		{
			while (i < 4) {
				env_spline_data_X[i] = h + i * dh; env_spline_data_Y[i] = 1.0;
				i++;
			}
		}
		else
		{
			while (h < thisPatch->top_height[INDEX])
			{
				success = 0;
				while (!success) // Take a single quality controlled step
				{
					// STEP WITH RUNGE-KUTTA DRIVER
					for (j = 1; j <= 6; j++) // RK STEPS -6; last step calculates error estimates  saved in der[0]
					{
						env_dy[j] = 1.0 - calculate_env_at_height(RK.return_x_after_step(h, dh, j - 1), t_index);
						y2 = RK.return_next_y(y1, dh, env_dy, j);
					}
					// CALULATE ERROR
					errmax = env_dy[0] / eps_env;
					if (errmax < 1.0 || dh == ENV_MINSTEP) // Successful step
					{
						// Check memory
						if (i + 50 > (int) env_spline_data_X.size())
						{
							env_spline_data_X.resize(i + 50, 0);
							env_spline_data_Y.resize(i + 50, 0);
						}
						success = 1;
						// Add to env_sample vector
						env_spline_data_X[i]   = RK.return_x_after_step(h, dh, 0);	env_spline_data_Y[i]   = 1.0 - env_dy[1];
						env_spline_data_X[i + 1] = RK.return_x_after_step(h, dh, 1);	env_spline_data_Y[i + 1] = 1.0 - env_dy[2];
						env_spline_data_X[i + 2] = RK.return_x_after_step(h, dh, 2);	env_spline_data_Y[i + 2] = 1.0 - env_dy[3];
						env_spline_data_X[i + 3] = RK.return_x_after_step(h, dh, 3);	env_spline_data_Y[i + 3] = 1.0 - env_dy[4];
						env_spline_data_X[i + 4] = RK.return_x_after_step(h, dh, 5);	env_spline_data_Y[i + 4] = 1.0 - env_dy[6];
						// NOTE - don't add env_dy[5] to spline because that would duplicate env_dy[1] from next step
						i += 5;   y1 = y2;   h += dh;
						if (h >= thisPatch->top_height[INDEX]) {
							env_spline_data_X[i] = h;
							env_spline_data_Y[i] = 1.0 - env_dy[5];
							i++;
						}
						// Estimate allowable step
						dh = max(ENV_MINSTEP, RK.grow_step(dh, errmax));
						// Choose last step so that includes top height + 1 extra step
						if (h + dh >= thisPatch->top_height[INDEX])
							dh = (thisPatch->top_height[INDEX] - h) / RK.return_x_after_step(0, 1.0, 5);
						// Record first successful step for resuse in next run
						if (flag)  {env_dh_start = dh;  flag = 0;}
					}
					else			// Falied step
					{
						dh = max(ENV_MINSTEP, RK.shrink_step(dh, errmax));
						// Make sure last step isn't decided here
						while (h + dh >= thisPatch->top_height[INDEX])
							dh /= 2.0;
					}

				}
			}
		}
		// Clear previous allocation if it exists
		if (thisPatch->Env_Spline[INDEX] != '\0') {
			gsl_spline_free (thisPatch->Env_Spline[INDEX]);
			gsl_interp_accel_free (thisPatch->Spline_Accel[INDEX]);
		}

		// Initiate spline
		thisPatch->Env_Spline[INDEX] = gsl_spline_alloc(gsl_interp_cspline, i);
		gsl_spline_init(thisPatch->Env_Spline[INDEX], &env_spline_data_X.front(), &env_spline_data_Y.front(), i) ;
		thisPatch->Spline_Accel[INDEX]  = gsl_interp_accel_alloc();

		// Record LAI
		if ((*p).c_ext > 0.0)  thisPatch->LAI[INDEX] = -log(gsl_spline_eval(thisPatch->Env_Spline[INDEX], 0.0, thisPatch->Spline_Accel[INDEX])) / (*p).c_ext;
		else thisPatch->LAI[INDEX] = 0.0;
	}
	else  // Not using spline interpolation
		thisPatch->LAI[INDEX] = calculate_env_at_height(0.0, t_index);
	return i;
}

// Temp gives  fractional decline in maximum available light at given size
double EBT_Base::calculate_env_at_height(double height, int t_index) {
	int n;
	list<int>::iterator It;
	double ENV_temp = 0.0; // ENV_TEMP contains cumulative leaf area above given size
	for (n = 0; n < no_res; n++)
	{
		for (It = Spp[n].cohort_list.begin(); It != Spp[n].cohort_list.end(); It++)
		{
			ENV_temp += Spp[n].strat.competiton(height, Spp[n].cohort[*It][t_index].height,
			                                    Spp[n].cohort[*It][t_index].mu) * Spp[n].cohort[*It][t_index].lam;
			if (!(exp(-(*p).c_ext * ENV_temp) >= 0) && !(exp(-(*p).c_ext * ENV_temp) < 0))
			{
				cout << "ERR\t" << t_index << setprecision(15) << " " << thisPatch->age << " " << *It
				     << " " << setprecision(4) << ENV_temp << " " << exp(-(*p).c_ext * ENV_temp)
				     << " " << !(exp(-(*p).c_ext * ENV_temp) > 0) << " " << !(exp(-(*p).c_ext * ENV_temp) < 0)
				     << "\t" << Spp[n].strat.competiton(height, Spp[n].cohort[*It][t_index].height,
				                                        Spp[n].cohort[*It][t_index].mu)
				     << " " << Spp[n].cohort[*It][t_index].lam << " "
				     << height << " " << Spp[n].cohort[*It][t_index].height << " "
				     << Spp[n].cohort[*It][t_index].mu << " "
				     << Spp[n].strat.Height(Spp[n].cohort[*It][t_index].mu) << endl;
			}
		}
	}
	if (!(exp(-(*p).c_ext * ENV_temp) >= 0) && !(exp(-(*p).c_ext * ENV_temp) < 0))
	{cerr << "Nan in calculate env_sample " << endl; abort();} // Exit(EXIT_FAILURE);}
	return exp(-(*p).c_ext * ENV_temp);	// Convert to relative light as exp(-k*L)
}


void EBT_Base::print_start(string dir, int option) {

	// Change directory
	if (!IsDirectory(dir.c_str()))     mkdir(dir.c_str());
	if (!IsDirectory(dir.c_str()))    {cerr << "Problem creating directory " << dir << endl; hold(); exit(EXIT_FAILURE);}
	chdir(dir.c_str());
// Open files
	OutFile_thisPatch = new ofstream;
	OutFile_thisPatch->open("patch_age.txt");
	*OutFile_thisPatch << "%Age\tp(a,t)\tgamma(a)\t";
	for (int i = 1; i <= 50; i++)	*OutFile_thisPatch << "H" << i << "\t";
	for (int i = 1; i <= 50; i++)	*OutFile_thisPatch << "E" << i << "\t";	*OutFile_thisPatch << endl;

	// Options for output - 0 = just print popn stats, 1 = also print states, 2 = also print rates
	if (option == 1)
		Nfiles = 1;
	else if (option == 2)
		Nfiles = 9;
	else
		Nfiles = 12;
	string names[12] = {"popn", "age", "coh_m", "coh_lam", "coh_r", "bound_m", "bound_s", "bound_n", "bound_r", "d_bound_m", "d_bound_s", "d_bound_r"};
	for (int n = 0; n < no_res; n++)
	{
		Spp[n].OutFile = new ofstream[14];
		for (int k = 0; k < Nfiles; k++)
		{
			Spp[n].OutFile[k].open((stringify(n) + "_" + names[k] + ".txt").c_str());
			if (!Spp[n].OutFile[k]) {cerr << "File " << names[k] + ".txt" << " not opened " << endl; hold(); exit(EXIT_FAILURE);}
			Spp[n].OutFile[k] << setiosflags(ios::scientific);
		}
		// Print headers
		Spp[n].OutFile[0] << "%Age\tNoCohorts\tSeed_Rain_in\tSeed_Rain_out\tSeed_Rain_in_tot\tSeed_Rain_out_tot\tPi\tN\tHeight_top\tHeight_Av\tBiom\tLAI\tGPP\tNPP\tTurnover\tRespiration" << endl;
	}
}

// PRINTES STATES & RATES OF CHANGE FOR GIVEN TIME STEP
// ASSUMES THAT FILES ARE ALREADY OPEN AND RATES OF CHANGE HAVE ALREADY BEEN CALULCATED
// NOTE USES T_NEXT, BECAUSE PRINT IS CALLED AFTER STEP_COHORTS, SO PREVIOUS STATES CONTAINED IN T_NEXT
void EBT_Base::print_states(int t_index, vector<int> max_cols, int option) {
	int count, k, n;
	double bound_low;
	list<int>::iterator It;
	// Patch age and light E
	*OutFile_thisPatch << setprecision(15) << thisPatch->age << "\t" << s->patch_age_density_freq(thisPatch->age) << "\t" << s->disturbance_rate(thisPatch->age) << "\t";
	double steps = thisPatch->top_height[0] / (50.0 - 1.0);
	for (double i = 1; i <= 50; i++)	*OutFile_thisPatch << steps*(i - 1) << "\t";
	for (double i = 1; i <= 50; i++)	*OutFile_thisPatch << get_env_at_height_from_spline(steps * (i - 1), thisPatch, 1) << "\t";
	*OutFile_thisPatch << endl;

	for (n = 0; n < no_res; n++)
	{
		Spp[n].OutFile[0] << setprecision(15) << thisPatch->age << "\t" << Spp[n].no_cohorts << "\t" << Spp[n].Seed_Rain_in[1]
		                  << "\t" << Spp[n].Seed_Rain_out[1] << "\t" << Spp[n].popn[t_index].Seed_Rain_in_cumulative << "\t"
		                  << Spp[n].popn[t_index].Seed_Rain_out_cumulative << "\t" << Spp[n].popn[t_index].Bottom_Coh_Mass
		                  << "\t" << Spp[n].popn[t_index].N << "\t" << Spp[n].popn[t_index].H_top << "\t" << Spp[n].popn[t_index].H_av
		                  << "\t" << Spp[n].popn[t_index].Biom << "\t" << Spp[n].popn[t_index].LAI << "\t" << Spp[n].popn[t_index].GPP
		                  << "\t" << Spp[n].popn[t_index].NPP << "\t" << Spp[n].popn[t_index].Turnover << "\t" << Spp[n].popn[t_index].Resp << endl;

		//%Age	NoCohorts	Seed_Rain_in	Seed_Rain_out	Seed_Rain_in_tot	Seed_Rain_out_tot	Pi	N	Height_top	Height_Av	Biom	LAI	GPP	NPP	Respiration	Turnover


		// PRINT STATES
		if (option > 1)
		{
			count = 0;
			It = Spp[n].cohort_list.end();
			while (It != Spp[n].cohort_list.begin())
			{	It--; count++;
				Spp[n].OutFile[1] << thisPatch->age - Spp[n].cohort[*It][t_index].time_at_birth << "\t";
				Spp[n].OutFile[2] << setprecision(15) << Spp[n].cohort[*It][t_index].mu << "\t";
				if (It != Spp[n].cohort_list.begin()) {
					It--; bound_low = Spp[n].cohort[*It][t_index].bound; It++;
				}
				else
					bound_low = Spp[n].strat.mass_at_birth();
				Spp[n].OutFile[3] << Spp[n].cohort[*It][t_index].lam << "\t";
				Spp[n].OutFile[4] << Spp[n].cohort[*It][t_index].coh_R << "\t";
				Spp[n].OutFile[5] << setprecision(15) << Spp[n].cohort[*It][t_index].bound << "\t";
				Spp[n].OutFile[6] << exp(Spp[n].cohort[*It][t_index].bound_S) << "\t";
				Spp[n].OutFile[7] << Spp[n].cohort[*It][t_index].bound_n0*exp(Spp[n].cohort[*It][t_index].bound_nt) << "\t";
				Spp[n].OutFile[8] << Spp[n].cohort[*It][t_index].bound_R << "\t";
			}
			// PRINT State at lower boundary
			if (Spp[n].cohort[*It][t_index].bound > Spp[n].strat.mass_at_birth())
			{
				count++;
				Spp[n].OutFile[1] << 0.0 << "\t";
				Spp[n].OutFile[2] << "NaN\t";
				Spp[n].OutFile[3] << "NaN\t";
				Spp[n].OutFile[4] << "NaN\t";
				Spp[n].OutFile[5] << Spp[n].strat.mass_at_birth() << "\t";

				// CALCULATE GROWTH & GERMINATION AT BOUNDARY
				get_env_at_height_from_spline_vector(Spp[n].strat.Height(Spp[n].strat.mass_at_birth()), &Spp[n].strat, thisPatch, 1, env_sample);
				Spp[n].OutFile[6] << Spp[n].strat.germination(Spp[n].strat.mass_at_birth(), &env_sample.front(), thisPatch->age ) << "\t";

				Spp[n].strat.Grow_Mort_Repro(GMR, Spp[n].strat.mass_at_birth(), &env_sample.front(), thisPatch->age);
				if (GMR[0] > 0.0)
					Spp[n].OutFile[7] << Spp[n].Seed_Rain_in[1] / GMR[0] << "\t";
				else
					Spp[n].OutFile[7] << 0.0 << "\t";
				Spp[n].OutFile[8] << 0.0 << "\t";
			}
			while (count < max_cols[n])
			{
				for (k = 1; k < 9; k++)
					Spp[n].OutFile[k] << "NaN" << "\t";
				count++;
			}
			for (k = 1; k < 9; k++)
				Spp[n].OutFile[k] << endl;
		}

		// PRINT RATES
		if (option > 2)
		{
			count = 0;
			It = Spp[n].cohort_list.end();
			while (It != Spp[n].cohort_list.begin())
			{	It--; count++;
				Spp[n].OutFile[9] << Spp[n].rates[*It].d_bound[1] << "\t";
				Spp[n].OutFile[10] << Spp[n].rates[*It].d_bound_S[1] << "\t";
				Spp[n].OutFile[11] << Spp[n].rates[*It].d_bound_R[1] / exp(Spp[n].cohort[*It][t_index].bound_S)*
				                   (*s).patch_weight(Spp[n].cohort[*It][t_index].time_at_birth, thisPatch->age) << "\t";
			}
			// PRINT rates at lower boundary
			if (Spp[n].cohort[*It][t_index].bound > Spp[n].strat.mass_at_birth())
			{
				count++;
				Spp[n].OutFile[9] << GMR[0] << "\t";
				Spp[n].OutFile[10] << -GMR[1] << "\t";
				Spp[n].OutFile[11] << GMR[2] << "\t";
			}
			while (count < max_cols[n])
			{
				for (k = 9; k < 12; k++)
					Spp[n].OutFile[k] << "NaN" << "\t";
				count++;
			}
			for (k = 9; k < 12; k++)
				Spp[n].OutFile[k] << endl;
		}
	}
}

void EBT_Base::print_end(void) {
	int i, j;
	// Print Strategy files
	ofstream OutFile;
	OutFile.open("Strategy.txt"); if (!OutFile) {cerr << "File profile not opened" << endl; hold(); exit(EXIT_FAILURE);}
	OutFile << "%Strat\tX_start\tX_end\t";
	for (j = 0; j < 4; j++)	OutFile << "Tr" << j << "\t";
	OutFile << endl;

	for (i = 0; i < no_res; i++)
	{
		OutFile << i << "\t" << Spp[i].X_start << "\t" << Spp[i].X_end << "\t";
		for (j = 0; j < 4; j++)
			OutFile << pow(10, Spp[i].strat.get_trait(j)) << "\t";
		OutFile << endl;
	}
	OutFile.close();
	// Print params
	OutFile.open("params.m");
	p->print_sim_params(OutFile);

	// Close files
	OutFile_thisPatch->close();
	delete OutFile_thisPatch;
	for (int n = 0; n < no_res; n++) {
		for (int k = 0; k < Nfiles; k++)
			Spp[n].OutFile[k].close();
		delete [] Spp[n].OutFile;
	}
	chdir("..");

	return;
}
