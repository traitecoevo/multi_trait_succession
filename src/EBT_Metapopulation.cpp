#include "EBT_Metapopulation.h"
#define MIN_COHORTS 30

EBT_Metapopulation::EBT_Metapopulation() {
	isSetup = 0;
	print = 0;
	continuous_seed_flow_flag = 0;
	thisPatch_store_isSetup = 0;
	mut_fit_dy.resize(9, 0);
	// Default values for error controls
	set_eps_cohort_splt(5e-1, 5);  // Without cohort weighting
	set_max_cohort_fitness(0.05);

	maxIterations = 10;
}


void EBT_Metapopulation::setup(sim_params * Pars) {
#if (DEBUG >=3)
	cout << "\t\t\tENTER EBT_Metapopulation::setup" << flush;
#endif
	// Remove all existing info
	// Clear_all();
	//---------------------------------------------------------------------------
	// Configure simulation parameters
	p = Pars;
	s = new site;
	s->set_type(1); // Weibull distribution
	set_max_cohort_fitness(0.05);
	mean_dist = 0.0; // By default set to 0.0 so that activates flag in resident fitness on first run
	//---------------------------------------------------------------------------
	// Allocate memory in Spp and thisPatch vectors, plus temp_spp
	temp_spp.set_dimensions(400);

	// Configure mutant and temp_spp
	update_sim_params();
	set_zero(temp_spp);


	// Remove existing strategies
	while (no_res > 0)	remove_spp(no_res - 1);

	isSetup = 1;
#if (DEBUG >=3)
	cout << "\t\t\tEND EBT_Metapopulation::setup" << endl;
#endif
}

void EBT_Metapopulation::set_max_cohort_fitness(double eps) {
	max_cohort_fitness = eps;
}

// Last row contains total for community
vector<vector<double> > EBT_Metapopulation::getPopulationProperties(void) {
	int i, n;
	vector<vector<double> > output(no_res + 1, vector<double>(10, 0.0));
	for (n = 0; n < no_res; n++)
	{
		output[n][0] = Spp[n].X_start;
		output[n][1] = Spp[n].N_cumulative;
		output[n][2] = Spp[n].H_top_cumulative;
		output[n][3] = Spp[n].H_av_cumulative;
		output[n][4] = Spp[n].Biom_cumulative; // Kg / m2
		output[n][5] = Spp[n].LAI_cumulative; // M2/m2
		output[n][6] = Spp[n].GPP_cumulative; // Kg /m2 /yr
		output[n][7] = Spp[n].NPP_cumulative; // Kg /m2/yr
		output[n][8] = Spp[n].Turnover_cumulative; // Kg/m2/yr
		output[n][9] = Spp[n].Resp_cumulative; // Kg/m2/yr
	}
	// Calculate whole popn values
	for (n = 0; n < no_res; n++)
	{
		for (i = 0; i < 10; i++)
			if (i == 2)		output[no_res][i] = max(output[no_res][i], output[n][i]);
			else if (i == 3)   output[no_res][i] += output[n][i] * output[n][5];
			else		output[no_res][i] += output[n][i];
	}
	output[no_res][3] = output[no_res][3] / output[no_res][5];
	return output;
}

double EBT_Metapopulation::mutantFitness(const double Traits[], bool with_competition) {
	if (!thisPatch_store_isSetup)
	{cout << "EBT:\tResident needs to be run before mutant"; return 0;}

	mut.set_traits(Traits);
	int i = 0, end_index = thisPatch_store.size() - 1;

	// For data collection
	vector<vector<double> > Data(end_index + 1, vector<double>(2, 0));
	for (i = 0; i <= end_index; i++)
	{
		Data[i][0] = thisPatch_store[i]->age;
		Data[i][1] = (*s).patch_age_density_freq(Data[i][0]) * mutant_individual_restiming(&mut, i, end_index, with_competition, 2, 0);
	}
	// Do integration
	double y = 0.0;
	for (i = 0; i < end_index - 1; i++)
		y += 0.5 * (Data[i][1] + Data[i + 1][1]) * (Data[i + 1][0] - Data[i][0]);
	return y;
}

double EBT_Metapopulation::mutant_individual_restiming(Strategy* strat, int t_start, int t_end, bool with_competition, int output, int print) {
#if (DEBUG >=4)
	cout << "\n\t\t\t\tENTER EBT_Base::mutant_individual_restiming";
#endif
	if (print) cout << "\nEBT:\tMutantIndiv: " << t_start << "\t" << t_end << "\t" << thisPatch_store[t_start]->age << endl;
	int i, rk_index, index = t_start;
	double dt, time, time_start;
	// Initial values
	time = time_start = thisPatch_store[t_start]->age;
	Y1[0] = strat->mass_at_birth(); // Size
	Y1[1] = Y1[2] = 0;			 // Log Survival, Seeds - meta

	// Set env_sample vector to maximum if no competition [MUTANT IN VIRGIN env_sample]
	if (!with_competition)
		for (i = 0; i < (*p).ENV_DIM; i++)  env_sample[i] = 1.0;
	else
		get_env_at_height_from_spline_vector(strat->Height(Y1[0]), strat, thisPatch_store[index], 1, env_sample);

	// Survival during germination
	double Pi1 = strat->germination(Y1[0], &env_sample.front(), time);
	if (Pi1 == 0)		return 0.0;
	else		Y1[1] = log(Pi1);

	for (i = 0; i < 3; i++)  Y2[i] = Y1[i];

	// Simulate for patch length
	while (index < t_end)
	{
		time = thisPatch_store[index]->age;
#if (DEBUG >=4)
		cout << "\n\t\t\t\t\t" << time << " - " << flush;
#endif

		dt = thisPatch_store[index + 1]->age - time;
		if (print) {cout << time << "\t" << dt << "\t" << Y1[0] << "\t" << Y1[1] << "\t" << Y1[2] << "\t" << endl;}
		// Take step
		for (rk_index = 1; rk_index <= 6; rk_index++)
		{
			// First step based on Y1, but Y2 should be equal to Y1 at start of each step
			if (with_competition) get_env_at_height_from_spline_vector(strat->Height(Y2[0]), strat, thisPatch_store[index], rk_index, env_sample);
			strat->Grow_Mort_Repro(GMR, Y2[0], &env_sample.front(), RK.return_x_after_step(time, dt, rk_index - 1));
			mut_dy[0][rk_index] = GMR[0];
			mut_dy[1][rk_index] = -GMR[1];
			mut_dy[2][rk_index] = exp(Y2[1]) * GMR[2] * (*s).patch_weight(time_start, RK.return_x_after_step(time, dt, rk_index - 1));
			Y2[0] = max(RK.return_next_y(Y1[0], dt, mut_dy[0], rk_index), strat->mass_at_birth());
			Y2[1] = min(RK.return_next_y(Y1[1], dt, mut_dy[1], rk_index), 0.0);
			Y2[2] = max(RK.return_next_y(Y1[2], dt, mut_dy[2], rk_index), 0.0);
			if (print)
				cout << "\t" << rk_index << "\t" << time << "\t" << dt << "\t" << Y2[0] << "\t"
				     << Y2[1] << "\t" << Y2[2] << "\t" << mut_dy[0][rk_index] << "\t"
				     << mut_dy[1][rk_index] << "\t" << mut_dy[2][rk_index] << endl;
		}
		// Update values
		for (i = 0; i < 4; i++)  Y1[i] = Y2[i];
		index++;
	}
#if (DEBUG >=4)
	cout << "\n\t\t\t\tEXIT mutant_individual_res_time" << endl;
#endif
	switch (output)
	{
	case (0): return Y1[0]; // Size
	case (1): return exp(Y1[1]); // Survival
	case (2): return Y1[2]; // Seeds
	}
	return 0;
}


// For printing of env_sample splines
void EBT_Metapopulation::print_splines(double top, double steps, bool all_RK_steps, string name) {
	double x, step = top / steps;
	ofstream OF_light(name.c_str());
	int i, r;

	// Print header
	OF_light << "%Age\tRK\tp(a)\tg(a)\tLAI\tTop\t";
	for (x = 0; x <= top; x += step)  OF_light << "H" << x << "\t";
	OF_light << endl;

	OF_light << "NaN\tNaN\tNaN\tNaN\tNaN\tNaN\t";
	for (x = 0; x <= top; x += step)  OF_light << x << "\t";
	OF_light << endl;

	for (i = 0; i < (int)thisPatch_store.size(); i++)
		for (r = 0; r <= 5 * ((all_RK_steps) ? 1 : 0); r++)
		{
			OF_light << thisPatch_store[i]->age << "\t" << r + 1 << "\t" << (*s).patch_age_density_freq(thisPatch_store[i]->age) << "\t" <<
			         thisPatch_store[i]->disturbance_rate[r] << "\t" << thisPatch_store[i]->LAI[r] << "\t" << thisPatch_store[i]->top_height[r] << "\t";
			for (x = 0; x <= top; x += step)
				OF_light << get_env_at_height_from_spline(x, thisPatch_store[i], r + 1) << "\t";
			OF_light << endl;
		}
	OF_light.close();
}

void EBT_Metapopulation::empty_thisPatch_store(void) {
	for (int i = 0; i < (int)thisPatch_store.size(); i++)
		delete thisPatch_store[i];  // Calls destructor for Patch class
	thisPatch_store.clear();
}

// Load time steps from file
void EBT_Metapopulation::load_time_steps(string filename) {
	cout << "\nLoading time step data " << default_time_steps.size() << " -> ";

	ifstream file(filename.c_str());
	if (!file) {cerr << "parameters file " << filename << " not opened " << endl; hold(); return;}
	default_time_steps.clear();
	double time;
	while (file >> time)
		default_time_steps.push_back(time);
	file.close();
	cout << default_time_steps.size() << "\tsuccess" << endl;
}


// Wrapper for resident  fitness - gives simplified interface
vector<double> EBT_Metapopulation::residentFitness(const double residentTraits[], const double X[], int N, int useDynamicCohorts) {
	return residentFitness(residentTraits, X, N, useDynamicCohorts, 1, 0, "empty");
}

void EBT_Metapopulation::printPopulationProfile(const double residentTraits[], const double X[], int N, string DIR) {
	printPopulationProfile(residentTraits, X, N, DIR, 0);
}

void EBT_Metapopulation::printPopulationProfile(const double residentTraits[], const double X[], int N, string DIR, bool printRates) {
	if (printRates)
		residentFitness(residentTraits, X, N, 2, 1, 3, DIR);
	else
		residentFitness(residentTraits, X, N, 2, 1, 2, DIR);
}


// Dynamics cohorts: 0 - use default spacing, no refinement;  1 - use default spacing with refinement;
//					2 - use existing spacing, no refinement; 3 - use existing spacing with refinement;
// Print options:    0 - no printing; 1 - print patch data only; 2 - also print states; 3 -also print rates
// Adaptive_time_step: 1 - Yes, 0 - no, use other spacing (needs to be loaded from file first)
vector<double> EBT_Metapopulation::residentFitness(const double residentTraits[], const double X[], int N, int useDynamicCohorts, bool adaptive_time_step, int print_options, string dir) {
#if (DEBUG >=3)
	cout << "\t\t\tENTER EBT_Metapopulation::residentFitness" << flush;
#endif
	int i, n, success = 0, K, K2;
	double time, dt;
	double dy, R1 = 0, R2, a1 = 0, a2, a3, R3, Rmax;
	list<int>::iterator It;

	if (useDynamicCohorts == 1 || useDynamicCohorts == 3)  check_coh_error = 2;
	else check_coh_error = 0;

	// Background checks
	if (!isSetup) {cerr << "EBT:\tError, need to run setup before simulation " << endl; hold(); exit(1);}
	for (i = 0; i < N; i++)
		if (X[i] < 0.0) {cerr << "EBT:\tNegative population size sent to EBT_resident fitness" << endl; hold(); exit(1);}

	// Check for change in disturbance interval - make sure last cohort is older than time_end
	if (mean_dist != (*p).log_mean_disturbance_interval)
	{
		time_end = (*s).solve_patchage_dist(pow(10.0, (*p).log_mean_disturbance_interval));
		mean_dist = (*p).log_mean_disturbance_interval;
		// Create default if doesn't already exist or replace existing if insufficiently long
		if ((int) cohort_intro_default.size() == 0 || cohort_intro_default[cohort_intro_default.size() - 1] < time_end - 10.0)
			set_default_cohort_introduction_times(time_end, 1E-5, 2.0, 0.2 * 10); // Approx 100 cohorts
	}
	// Check for change in simulation parameters, requiring new memory allocation
	if ((int)env_sample.size() != (*p).ENV_DIM)
	{
		update_sim_params();
	}

	// Check number species
	while (no_res > N)
		remove_spp(no_res - 1);  // Remove excess strategies
	for (i = 0; i < no_res; i++)   // Reset existing spp
		Spp[i].strat.set_traits(residentTraits + i * TRAIT_DIM );
	while (i < N) {                 // Add more strategies if needed
		add_spp(residentTraits + i * TRAIT_DIM ); i++;
	}
	vector<int> numberOfSplits(no_res, 0), numberOfMerges(no_res, 0);
	int totalNumberOfSplits = 0, totalNumberOfMerges = 0;
	// Set seed input
	for (i = 0; i < N; i++)  Spp[i].X_start = X[i];

	// Make queue of cohort introduction times - use default or existing spacing
	vector<queue<double> > cohort_intro2(no_res);
	for (n = 0; n < no_res; n++)
	{
		if (useDynamicCohorts > 1 && Spp[n].no_cohorts > 0)	// Use spacing from previous simulation
		{
			It = Spp[n].cohort_list.end();
			while (It != Spp[n].cohort_list.begin())
			{
				It--;
				cohort_intro2[n].push(Spp[n].cohort[*It][t_curr].time_at_birth);
			}
		}
		else												// Use default spacing
			for (i = 0; i < (int) cohort_intro_default.size(); i++)
				cohort_intro2[n].push(cohort_intro_default[i]);
		// Check last cohort is older than time_end
		if (cohort_intro2[n].back() < time_end)
			cohort_intro2[n].push(time_end + 1.0);
	}
	num_sim = 0;
	vector<int> print_columns(no_res);

	if (!adaptive_time_step)
		if ((int) default_time_steps.size() == 0) {
			cerr << "EBT:\tError, time steps need to be loaded first" << endl;
			exit(1);
		}

	while ( !success && num_sim < maxIterations) // Repeat until error in simulation is satisfied or reach maximum number of attempts
	{
		K = ode_fail = min_fail = 0;
		K2 = 0;			// Only used for fixed time-stepping
		num_sim++;
		t_curr = 0; t_next = 1;

		// FOR PRINTING
		if (print_options >= 1)
		{
			print_start(dir, print_options);
			for (n = 0; n < no_res; n++)
				print_columns[n] = cohort_intro2[n].size() + 1;
		}

		// Reset interpolation
		empty_thisPatch_store();
		thisPatch->age = 0;

		// Reset cohorts && Check memory available
		for (n = 0; n < no_res; n++)
		{
			if ((int)cohort_intro2[n].size()  >= Spp[n].dim_cohort)
				Spp[n].set_dimensions(cohort_intro2[n].size() + 50);
			set_zero(Spp[n]);	// Zeros and make new cohort stack
		}

		// Set initial values
		dt = 1.0;
		time = thisPatch->age = 0.0;
		max_d_LAI = 0.0;

		// Simulate until end time
		while (time < time_end) {
#if (DEBUG >=3)
			cout << "\n\t\t\t\t" << time << " - " << flush;
#endif

			// Increase minimum allowable step size
			K++;
			if (K == max_no_cohort_steps)
			{
				cerr << "EBT:\tError, too many steps in residentFitness " << time << "\t" << thisPatch_store.size()
				     << "\t" << time_end << endl;
				for (n = 0; n < no_res; n++)
					cout << setprecision(15) << residentTraits[n * TRAIT_DIM + 0] << "\t" << residentTraits[n * TRAIT_DIM + 1] << "\t" << X[n] << endl;
				exit(1);
			}

			// Internalise boundary cohort & determine time for next step

			// If using fixed stepping, draw next time slot from bank
			if (!adaptive_time_step && K2 < (int)default_time_steps.size())
			{
				while (default_time_steps[K2] <= time) K2++;
				dt = default_time_steps[K2] - time;
			}

			dt = min(time_end - time, dt);
			for (n = 0; n < no_res; n++)
			{
				if (equals(cohort_intro2[n].front(), time, 1e-8))
				{
					begin_new_cohort(n);
					cohort_intro2[n].pop();
				}
				if (!cohort_intro2[n].empty())
				{
					if (cohort_intro2[n].size() > 80000)
					{
						cout << "EBT:\t Weird queue behaviour ";
						queue<double> temp2;
						cohort_intro2[n] = temp2;
					}
					else
						dt = min(cohort_intro2[n].front() - time , dt);
				}
			}
			// Take step after which time becomes time_next,
			step_cohorts(time, dt, adaptive_time_step);
			if (print_options >= 1)
				print_states(t_next, print_columns, print_options);


			// Copy thisPatch info & allocate new thisPatch data
			thisPatch_store.push_back(thisPatch);
			thisPatch = nextPatch;
			nextPatch = '\0';
		}
		// FINISHED SIMULATION	- Take one final step

		make_env_spline(t_curr, 1);

		// PRINTING TO FILE
		if (print_options >= 1) {
			calculate_rates_of_change(t_curr, time, 1);
			print_states(t_curr, print_columns, print_options);
			print_end();
			vector<vector<double> > AGG;
			AGG = getPopulationProperties();
			ofstream OF((dir + "/stand.txt").c_str());
			OF << "%Seed_rain\tN\tH_top\tH_av\tBiom\tLAI\tGPP\tNPP\tTurn\tResp_main" << endl;

			for (n = 0; n <= no_res; n++) // Do for each species
			{
				for (int ii = 0; ii < 10; ii++)
					OF << AGG[n][ii] << "\t";
				OF << endl;
			}
			OF.close();
		}

		// Copy thisPatch info & allocate new thisPatch data
		thisPatch_store.push_back(thisPatch);
		thisPatch = new Patch();
		delete	nextPatch; nextPatch = '\0';

		// Print to screen
		if (print) {
			for (n = 0; n < no_res; n++) // Do for each species
			{
				if (num_sim == 1)
					cout << "EBT:\tSimRES\t" << n << "\t" << setprecision(4) << residentTraits[n * TRAIT_DIM + 0] << "\t"
					     << residentTraits[n * TRAIT_DIM + 1] << "\t" << residentTraits[n * TRAIT_DIM + 2] << "\t" << pow(10, (*p).log_mean_disturbance_interval) << "\t";
				else
					cout << "\t\t\t\t\t\t\t";
				cout << num_sim << "\t" << n << " " << Spp[n].no_cohorts << " " << K << "\t" << X[n] << "\t"
				     << Spp[n].popn[t_curr].Seed_Rain_out_cumulative << "\t\t" << ode_fail << " " << min_fail;
				if (n != (no_res - 1)) cout << endl;
			}

		}

		// ASSESS COHORT SPACING AND REFINE WHERE NECESSARY
		success = 1;
		if (useDynamicCohorts == 1 || useDynamicCohorts == 3)
		{

			// Check merging and splitting
			for (n = 0; n < no_res; n++)
			{
				// Rmax - sets denominator when determining contribution to fitness. Limit maximum so that
				// Residents with v high fitness don't get divided too finely.
				// I.e. If fitness > 1.0, limits contribution to max coh_fitness
				Rmax = max(Spp[n].popn[t_curr].Seed_Rain_out_cumulative / Spp[n].X_start, 1.0);
				a3 = time_end, R3 = 0.0;
				int reached_top = 0;
				for (It = Spp[n].cohort_list.begin();  It != Spp[n].cohort_list.end(); It++)
				{
					a2 = Spp[n].cohort[*It][t_curr].time_at_birth;
					R2 = (*s).patch_age_density_freq(a1) * Spp[n].cohort[*It][t_curr].bound_R;
					// Check merging
					It++;
					if (It != Spp[n].cohort_list.end())	// Retrieve values for cohort above if not top cohort
					{
						a1 = Spp[n].cohort[*It][t_curr].time_at_birth;
						R1 = (*s).patch_age_density_freq(a1) * Spp[n].cohort[*It][t_curr].bound_R;
					}
					else
						reached_top = 1;
					It--;
					// Only check if true for all previous checks - use twice actual error to be conservative about merging
					if (!reached_top && (Spp[n].cohort[*It][t_curr].number_of_splits == -10.0))
						if ( (2.0 * 0.5 * (R1 + R3) * (a3 - a1)) > max_cohort_fitness * Rmax )
							Spp[n].cohort[*It][t_curr].number_of_splits = 0.0;
					// Check all cohorts for splitting
					if ( (0.5 * (R2 + R3) * (a3 - a2)) > max_cohort_fitness * Rmax )
						Spp[n].cohort[*It][t_curr].number_of_splits = max(Spp[n].cohort[*It][t_curr].number_of_splits,
						    floor(0.5 * (R2 + R3) * (a3 - a2) / (max_cohort_fitness * Rmax)));
					// Update values
					a3 = a2; R3 = R2;
				}
			}

			// EVALUATE WHETHER TO REDO SIMULATION - more than 2 splits required
			totalNumberOfSplits = 0, totalNumberOfMerges = 0;
			for (n = 0; n < no_res; n++)
			{

				numberOfSplits[n] = numberOfMerges[n] = 0;
				for (It = Spp[n].cohort_list.begin(); It != Spp[n].cohort_list.end(); It++)
				{
					// Limit number of splits to 10
					Spp[n].cohort[*It][t_curr].number_of_splits = min(10.0, Spp[n].cohort[*It][t_curr].number_of_splits);
					numberOfSplits[n] += (int) max(0.0, Spp[n].cohort[*It][t_curr].number_of_splits);
					if ( Spp[n].cohort[*It][t_curr].number_of_splits < 0) numberOfMerges[n]++;
				}
				totalNumberOfSplits += numberOfSplits[n];
				totalNumberOfMerges += numberOfMerges[n];
			}

			if (print)
				cout << "\t" << totalNumberOfSplits << "\t" << totalNumberOfMerges << endl;

			// EXECUTE SPLITTING  / MERGING IF NEEDED
			if (totalNumberOfSplits > 2)
			{
				success = 0;
				// Make new queue of cohort introduction times
				double k, cohort_intro_time, time_next;
				int enter_next = 1;
				for (n = 0; n < no_res; n++)
				{
					double N_next = Spp[n].no_cohorts - 0.5 * numberOfMerges[n] + numberOfSplits[n];

					// Remove existing introduction times
					while (!cohort_intro2[n].empty()) cohort_intro2[n].pop();

					It = Spp[n].cohort_list.end();
					while (It != Spp[n].cohort_list.begin())
					{
						It--;
						if ( num_sim < 3  && // only allows merges in first two iterations - saves possible conflicts between merging and splitting criteria, leading to many iterations
						     !enter_next &&// only test merge if previous cohort wasn't set to merge
						     (Spp[n].cohort[*It][t_curr].number_of_splits == -10.0) && // is this cohort marked for merging?
						     (N_next  > MIN_COHORTS) && // is the estimated number of cohorts > minimum number
						     (Spp[n].cohort[*It][t_curr].time_at_birth != 0.0)) // Don't merge first cohort
							enter_next = 1; // Merge so do nothing but do enter next cohort
						else
						{
							// Add current cohort introduction time
							cohort_intro_time = Spp[n].cohort[*It][t_curr].time_at_birth;
							cohort_intro2[n].push(cohort_intro_time);
							enter_next = 0;

							// Check to see if need splitting
							if (Spp[n].cohort[*It][t_curr].number_of_splits > 0)
							{
								if (It == Spp[n].cohort_list.begin())
									time_next = time_end;
								else {
									It--;
									time_next = Spp[n].cohort[*It][t_curr].time_at_birth;
									It++;
								}

								double NS = Spp[n].cohort[*It][t_curr].number_of_splits + 1;
								// New method - maintains fewer cohorts with multiple spp

								dt = (time_next - cohort_intro_time) / NS;
								for (k = 1; k < NS; k++)
									cohort_intro2[n].push(cohort_intro_time + k * dt);
							}
						}
					}

					// Add last cohort after time end
					cohort_intro2[n].push(time_end + 1.0);
				}
			}
		}
	}

	// FINISHED - calculate fitness for output -
	// code similar to mutant fitness - to simplify, why not call mutant fitness?
	vector<double> X_out(no_res, 0);
	for (n = 0; n < no_res; n++)
	{
		mut.set_traits(residentTraits + n * TRAIT_DIM);
		int end_index = thisPatch_store.size() - 1;
		dy = 0;
		// Fitness based on cohorts
		Spp[n].X_end = Spp[n].popn[t_curr].Seed_Rain_out_cumulative;
		// Fitness from boundary points
		It = Spp[n].cohort_list.end();
		a1 = a2 = R1 = R2 = 0.0;
		i = 0;
		while (It != Spp[n].cohort_list.begin())
		{
			It--;
			while (thisPatch_store[i]->age < Spp[n].cohort[*It][t_curr].time_at_birth)
			{
				a2 = thisPatch_store[i]->age;
				R2 = (*s).patch_age_density_freq(a2) * mutant_individual_restiming(&mut, i, end_index, 1, 2, 0);
				dy += 0.5 * (R1 + R2) * (a2 - a1);
				a1 = a2; R1 = R2;
				i++;
			}
			a2 = thisPatch_store[i]->age;
			R2 = (*s).patch_age_density_freq(a2) * Spp[n].cohort[*It][t_curr].bound_R;
			dy += 0.5 * (R1 + R2) * (a2 - a1);
			a1 = a2; R1 = R2;
			i++;
		}
		X_out[n] = dy;
	}
	thisPatch_store_isSetup = 1;
#if (DEBUG >=3)
	cout << "\t\t\tEXIT EBT_Metapopulation::residentFitness" << endl;
#endif
	return X_out;
}


double EBT_Metapopulation::get_cohort_residentFitness(int spp_n) {
	return Spp[spp_n].X_end / Spp[spp_n].X_start;
}
