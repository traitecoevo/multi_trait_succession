#include "AdaptiveDynamics_nD.h"
#include "gsl_random_numbers.h"

AdaptiveDynamics_nD::AdaptiveDynamics_nD() {
	precision = 6;

	isSetup = 0;
	displayWorkingDetails = 1;
	numberOfResidents = 0;
};

void AdaptiveDynamics_nD::setDisplayWorkingDetails(int option) {
	displayWorkingDetails = option;
}

void AdaptiveDynamics_nD::setup(PopulationDynamics* myPopulationDynamicsPointer) {
	popDyn = myPopulationDynamicsPointer;
	popDyn->resetPopulationSize = 0;
	popDyn->displayWorkingDetails = 1;
	isSetup = 1;
}


// Reads traits values stored in File, as written by runStochasticModel
void AdaptiveDynamics_nD::inputTraitsFromStochFile(double TraitsToInput[], int &numRes, double &Time, string filename) {
	inputTraitsFromStochFile(TraitsToInput, numRes, Time, filename, 1E6);
}


// Reads traits values stored in File, as written by runStochasticModel
void AdaptiveDynamics_nD::inputTraitsFromStochFile(double TraitsToInput[], int &numRes, double &Time, string filename, int line) {
	cout << "Load Resident state from " << filename << endl;
	ifstream inFile(filename.c_str()); if (!inFile) {cerr << "In inputTraitsFromCEFile: file " << filename << " can't be found " << endl; exit(1);}
	int i, j;
	double *X = new double[100];
	// Get first line

	string temp;
	getline(inFile, temp);

	j = 1;
	cout << line << endl;
	while (inFile.peek() != EOF && j <= line)
	{
		inFile >> Time >> numRes;
		for (i = 0; i < numRes * TRAIT_DIM; i++) inFile >> TraitsToInput[i]; // Read in traits
		for (i = 0; i < numRes; i++)  inFile >> X[i]; // Read in population size
		getline(inFile, temp);
		j++;
	}
	// Print
	cout << "Time " << Time << endl;
	for (i = 0; i < numRes; i++)
	{
		cout << i << "\t" << X[i] << "\t";
		for (j = 0; j < TRAIT_DIM; j++)
			cout << TraitsToInput[i * TRAIT_DIM + j] << "\t";
		cout << endl;
	}

	cout << "Set cohort spacing " << numRes << endl;
	popDyn->setResidentState(TraitsToInput, X, numRes);
	cout << "Done" << endl;
	popDyn->resetEbtCohorts = 0;
	delete [] X;
}


// Runs stochastic assembly model without time-scale separation between ecological and evolutionary dynamics.
// MutationRate is the COMMUNITY LEVEL number of mutations per step of population dynamics
void AdaptiveDynamics_nD::runStochasticModelNonEquil(double startTraits[], int numRes, double mutationVariance,  double traitRange[2][3], double mutationRate, double immigrationRate, double start_time, double end_time, double printFreq, string name) {
	vector<double> covarMatrixDiagonal(4, 0);
	covarMatrixDiagonal[0] = covarMatrixDiagonal[1] = mutationVariance;
	runStochasticModelNonEquil(startTraits, numRes, covarMatrixDiagonal,  traitRange, mutationRate, immigrationRate, start_time, end_time, printFreq, name);
}

void AdaptiveDynamics_nD::runStochasticModelNonEquil(double startTraits[], int numRes, vector<double> covarMatrixDiagonal, double traitRange[2][3], double mutationRate, double immigrationRate, double start_time, double end_time, double printFreq, string name) {
	if (!isSetup)  cerr << "Stochastic Model not setup properly" << endl;

	int i, j,  max_numRes = 100;
	double *Res =  vector_NR(0, max_numRes * TRAIT_DIM);
	double *Mut =  vector_NR(0, TRAIT_DIM);
	double Rate;
	gsl_random_numbers Ran; // Random number generation
	double Time = start_time;
	vector<double> Mutation(TRAIT_DIM);

	// Vector to determine which trait mutation happens in
	vector<double> TraitEvolving(TRAIT_DIM, 0);
	for (j = 0; j < TRAIT_DIM; j++)
		if (covarMatrixDiagonal[j] > 0) TraitEvolving[j]++;

	vector<double> TraitMutationCumulativeProbability(TRAIT_DIM + 1, 0);
	for (j = 0; j < TRAIT_DIM; j++) {
		for (int k = 0; k <= j; k++ )
			TraitMutationCumulativeProbability[j] += TraitEvolving[k];
	}
	for (j = 0; j < TRAIT_DIM; j++)
		TraitMutationCumulativeProbability[j] = TraitMutationCumulativeProbability[j] / TraitMutationCumulativeProbability[TRAIT_DIM - 1];

	cout << "cumulative mutation probability: ";
	for (j = 0; j < TRAIT_DIM; j++)
		cout << TraitMutationCumulativeProbability[j] << "\t";
	cout << endl;

	// Output file
	ofstream OutFile(name.c_str()); if (!OutFile)
	{cerr << "File not opened" << endl; hold(); exit(EXIT_FAILURE);}

	// Set base traits
	for (j = 0; j < TRAIT_DIM; j++)
		Mut[j] = startTraits[j];
	for (i = 0; i < max_numRes; i++)
		for (j = 0; j < TRAIT_DIM; j++)
			Res[i * TRAIT_DIM + j] = startTraits[(i < numRes) ? (i * TRAIT_DIM + j) : j];

	// Popn dynamics
	popDyn->stopOnDemographicRoot = 0;
	popDyn->maxPopulationIterations = 1;

	cout << "start" << endl;
	double lastPrintTime = 0;
	while (Time < end_time) {

		// Step population to next time step
		popDyn->maxPopulationIterations = 1;
		popDyn->resolveEquilibrium = 1;
		popDyn->S(Mut, Res, numRes);

		// PRINT CURRENT STATE
		OutFile << scientific << Time << "\t" << numRes << "\t";
		for (i = 0; i < numRes; i++)
			for (j = 0; j < TRAIT_DIM; j++)
				OutFile << setprecision(precision) << Res[i * TRAIT_DIM + j] << "\t";
		for (i = 0; i < numRes; i++) OutFile << popDyn->getPopulationSize(i) << "\t";
		OutFile << endl;

		cout << "STEP:\t" << Time << " " << numRes << "\t";
		for (i = 0; i < numRes; i++) cout << Res[i * TRAIT_DIM + 0] << " " << Res[i * TRAIT_DIM + 1] << "\t"; cout << endl;

		// Print every population state & fitness landscape every so often, name rounded to nearest 100
		if (Time - lastPrintTime >= printFreq) {
			if (covarMatrixDiagonal[0] > 0 & covarMatrixDiagonal[1] > 0) { // print 2D landscape if both traits evolving
				outputFitnessLandscape2D(Res, numRes, traitRange, "T-" + stringify(floor(Time * 100 + 0.5) / 100), 1);
			}
			else {
				for (i = 0; i < 2; i++)
					if (covarMatrixDiagonal[i] > 0) // Only print traits with positive Variance
						outputFitnessLandscape1D(Res, numRes, Res, i, traitRange[i], "T" + stringify(i) + "-" + stringify(floor(Time * 100 + 0.5) / 100), 1);
			}
			lastPrintTime = Time;
		}

		cout << "REMOVE EXTINCT" << endl;
		// Remove any extinct populations
		while (popDyn->checkForExtinctPopns()) {
			i = 0;
			while (i < numRes) {
				if (popDyn->isExtinct(i)) {
					cout << "\t\t - extinction " << i << "\t";
					for (j = 0; j < TRAIT_DIM; j++)  cout << Res[i * TRAIT_DIM + j] << " "; cout << endl;

					popDyn->removeResident(i);
					for (int ii = i + 1; ii < numRes; ii++)
						for (j = 0; j < TRAIT_DIM; j++)
							Res[(ii - 1)*TRAIT_DIM + j] = Res[ii * TRAIT_DIM + j];
					numRes--;
				}
				else
					i++;
			}
		}

		cout << "MUTATIONS" << endl;
		// Draw mutations from existing residents, keeping community level mutation rate constant.
		// Distribute mutations across residents according to population size.

		int nMutations, nNew = 0;

		double totPopSize = 0;

		for (i = 0; i < numRes; i++)
			totPopSize += popDyn->getPopulationSize(i);
		cout << totPopSize << "\n";

		for (i = 0; i < numRes; i++) {
			// Mutation rate for species i

			Rate = popDyn->getPopulationSize(i) / totPopSize * mutationRate;

			// Draw number of events from Poison distribution
			nMutations = Ran.RandomNumber_Poisson(Rate);
			cout << "\tS" << i << "\t" << Rate << "\t" << nMutations << endl;
			for (int m = 0; m < nMutations; m++) {

				for (j = 0; j < TRAIT_DIM; j++) Mutation[j] = 0;

				// Assume only have mutation in one trait at a time. Choose which trait mutation happens in
				double pT = Ran.RandomNumber_Uniform_1D();
				j = 0;
				while (pT > TraitMutationCumulativeProbability[j]) j++;
				Mutation[j] = Ran.RandomNumber_Gaussian_1D(covarMatrixDiagonal[j]);
				for (j = 0; j < TRAIT_DIM; j++)
					Mut[j] = Res[(numRes + nNew) * TRAIT_DIM + j] = Res[i * TRAIT_DIM + j] + Mutation[j];
				cout << "\t\tMutation " + stringify(i) + "\t| ";
				for (j = 0; j < TRAIT_DIM; j++) cout << Mutation[j] << "\t"; cout << "| ";
				for (j = 0; j < TRAIT_DIM; j++) cout << Mut[j] << "\t ";	cout << "| ";
				cout << popDyn->S(Mut, Res, numRes) << endl;
				nNew++;
			}
		}
		// Immigrants
		Rate = immigrationRate;
		nMutations = Ran.RandomNumber_Poisson(Rate);
		cout << "\tI\t" << Rate << "\t" << nMutations << endl;

		for (int m = 0; m < nMutations; m++) {
			// Draw from uniform distribution,  only mutate traits with positive variance
			bool isViable = 0;
			while (!isViable) { // Mutant must have positive fitness in virgin environment
				for (i = 0; i < 2; i++)
					if (covarMatrixDiagonal[i] > 0) // Only mutate traits with positive Variance
						Mut[i] =  Ran.RandomNumber_Uniform_1D(traitRange[i][0], traitRange[i][1]);	// Draw mutation from 2-D uniform distribution
				isViable = (popDyn->S(Mut, Res, 0) > 0);
			}
			cout << "\t\tImmigrant " + stringify(m) << "\t";
			for (j = 0; j < TRAIT_DIM; j++) cout << Mut[j] << "\t"; cout << "| ";
			cout << popDyn->S(Mut, Res, numRes) << endl;

			// Copy to resident vector
			for (j = 0; j < TRAIT_DIM; j++) // Add new resident
				Res[(numRes + nNew)*TRAIT_DIM + j] = Mut[j];
			nNew++;
		}

		numRes += nNew;
		Time += 1;
	}
	OutFile.close();
	free_vector_NR(Mut, 0, TRAIT_DIM);
	free_vector_NR(Res, 0, max_numRes * TRAIT_DIM);
}

void AdaptiveDynamics_nD::outputFitnessLandscape2D(const double residents[], int numRes, double traitRange[2][3], string name, bool print_res) {
	ofstream OutFile((name + "-2D.txt").c_str()); if (!OutFile) {std::cerr << "File not opened" << std::endl; hold(); exit(EXIT_FAILURE);}
	OutFile << "%";
	for (int i = 0; i < numRes; i++)
		OutFile << residents[i * TRAIT_DIM + 0] << "\t" << residents[i * TRAIT_DIM + 1] << "\t" << popDyn->getPopulationSize(i) << "\t";
	OutFile << endl;

	double i, j;
	double* mutant =  vector_NR(0, TRAIT_DIM);
	for (int k = 0; k < TRAIT_DIM; k++) mutant[k] = residents[k];
	if (print_res)
		popDyn->printPopulationProfile(residents, numRes, (name + "-R.txt"));
	for (i = traitRange[0][0]; i <= traitRange[0][1]; i += traitRange[0][2])
	{
		mutant[0] = i;
		for (j = traitRange[1][0]; j <= traitRange[1][1]; j += traitRange[1][2])
		{
			mutant[1] = j;
			OutFile << i << "\t" << j << "\t" << popDyn->S(mutant, residents, numRes) << endl;
		}
	}
	OutFile.close();
	free_vector_NR(mutant, 0, TRAIT_DIM);
}

void AdaptiveDynamics_nD::outputFitnessLandscape1D(const double residents[], int numRes, double* startTraits, int focalTraitIndex, double traitRange[3], string name, bool print_res) {
	cout << "start outputFitnessLandscape1D" << endl;

	ofstream OutFile((name + "-1D.txt").c_str());
	if (!OutFile) {std::cerr << "File not opened" << std::endl; hold(); exit(EXIT_FAILURE);}
	OutFile << "%";
	for (int i = 0; i < numRes; i++)
		OutFile << residents[i * TRAIT_DIM + 0] << "\t" << residents[i * TRAIT_DIM + 1] << "\t" << popDyn->getPopulationSize(i) << "\t";
	OutFile << endl;
	if (print_res)
		popDyn->printPopulationProfile(residents, numRes, (name + "-R.txt"));
	vector<vector<double> > OUT = getSliceOfFitnessLandscape1D(residents, numRes, startTraits, focalTraitIndex, traitRange);
	for (int i = 0; i < (int) OUT.size(); i++)
		OutFile << OUT[i][0] << "\t" << OUT[i][1] << endl;
	OutFile.close();
	cout << "end outputFitnessLandscape1D" << endl;

}

// Returns slice of fitness along particular trait dimension
vector<vector<double> > AdaptiveDynamics_nD::getSliceOfFitnessLandscape1D(const double residents[], int numRes, const double startTraits[], int focalTraitIndex, double traitRange[3]) {
	if (!isSetup) {cerr << "fitness function not set" << endl; hold(); exit(EXIT_FAILURE);}

	// For data storage
	int i, dim = (int)ceil((traitRange[1] - traitRange[0]) / traitRange[2]) + 1;
	vector<vector<double> > OUT(dim, vector<double>(2, 0.0));

	// Copy traits
	double* mutant =  vector_NR(0, TRAIT_DIM);
	for (i = 0; i < TRAIT_DIM; i++) mutant[i] = startTraits[i];

	// Fill vector
	for (i = 0; i < dim; i++)
	{
		mutant[focalTraitIndex] = OUT[i][0] = traitRange[0] + traitRange[2] * i;
		OUT[i][1] = popDyn->S(mutant, residents, numRes);
		cout << mutant[0] << "\t" << mutant[1] << "\t" << OUT[i][1] << endl;
	}
	free_vector_NR(mutant, 0, TRAIT_DIM);
	return OUT;
}
