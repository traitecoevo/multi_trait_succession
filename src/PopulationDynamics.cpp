#include "PopulationDynamics.h"
#include "MatrixCPP.h"
#include "Collect.h"

#include <math.h>

PopulationDynamics::PopulationDynamics() {
	numberOfResidents = 0;
	residentTraits = NULL;
	populationSize = logPopulationSize = NULL;

	displayWorkingDetails = 1;
	resetPopulationSize = 1;
	resetEbtCohorts = 1;
	resolveEquilibrium = 0;
	useDynamicCohorts = 1;
	stopOnDemographicRoot = 1;

	maxPopulationSize = 2e4;
	maxLogPopulationSize = log10(maxPopulationSize);
	maxPopulationIterations = 1;

}

void PopulationDynamics::resetFitnessSolver(void) {
	resetEbtCohorts = 1;
}

void PopulationDynamics::setup(FitnessSolver* myFitnessSolverPointer) {
	myFitnessSolver = myFitnessSolverPointer;
}

void PopulationDynamics::setDisplayWorkingDetails(int option) {
	displayWorkingDetails = option;
}

void PopulationDynamics::setPopulationSize(int i, double new_N) {
	if (new_N < 0.0)
	{
		cout << "POPDYN: Error - negative pop size requested " << i << " " << new_N << endl;
		populationSize[i] = EXTINCT;
	}
	else {
		populationSize[i] = min(maxPopulationSize, new_N);
		if (maxPopulationSize <= new_N)
		{
			cout << "POPDYN: Error - max population size requested " << i << " " << new_N << endl;
//			hold();
		}
	}
	logPopulationSize[i] = log10(populationSize[i]);
}

void PopulationDynamics::setLogPopulationSize(int i, double new_logN) {
	logPopulationSize[i] = min(maxLogPopulationSize, new_logN);
	populationSize[i] = pow(10.0, logPopulationSize[i]);

	if (maxLogPopulationSize < new_logN)
	{
		cout << "POPDYN: Error - max population size requested " << i << " " << new_logN << endl;
//		hold();
	}

}


int PopulationDynamics::getnumberOfResidents(void) const
{
	return numberOfResidents;
}

double PopulationDynamics::getPopulationSize(int i) const
{
	return populationSize[i];
}


bool PopulationDynamics::isExtinct(int i) {
	return isExtinct(i, EXTINCT);
}

bool PopulationDynamics::isExtinct(int i, double threshold) {
	if (populationSize[i] < threshold)
		return 1;
	return 0;
}

bool PopulationDynamics::checkForExtinctPopns(void) {
	return checkForExtinctPopns(EXTINCT);
}

bool PopulationDynamics::checkForExtinctPopns(double threshold) {
	for (int i = 0; i < numberOfResidents; i++)
		if (isExtinct(i, threshold))
			return 1;
	return 0;
}

int PopulationDynamics::removeResident(int index) {
	if (index < 0 || index >= numberOfResidents) {
		cout << "bad index given to RemoveResident function in popdyn " << index << "\t" << numberOfResidents << endl;
		exit(1);
	}

	myFitnessSolver->remove_spp(index);      // remove from EBT

	// Copy other values across in arrays, if not the last species
	if (index < numberOfResidents - 1)
		for (int i = index + 1; i < numberOfResidents; i++) {
			populationSize[i - 1] = populationSize[i];
			logPopulationSize[i - 1] = logPopulationSize[i];
			for (int j = 0; j < TRAIT_DIM; j++)
				residentTraits[(i - 1)][j] = residentTraits[i][j];
		}
	// Decrement number of residents
	numberOfResidents--;
	return numberOfResidents;
}


void PopulationDynamics::displayResidentTraits(void) {
	int i, j;
	std::cout << "numberOfResidents\t" << numberOfResidents << std::endl;
	std::cout << "Str:\t"; for (i = 0; i < numberOfResidents; i++) std::cout << i << "\t"; std::cout << endl;
	std::cout << "Pop:\t"; for (i = 0; i < numberOfResidents; i++) std::cout << populationSize[i] << "\t"; std::cout << endl;
	std::cout << "X:\t" << endl;
	for (j = 0; j < TRAIT_DIM; j++) {
		std::cout << j << "\t";
		for (i = 0; i < numberOfResidents; i++) std::cout << residentTraits[i][j] << "\t";
		std::cout << endl;
	}
	std::cout << endl;
}

void PopulationDynamics::printPopulationProfile(const double residents[], int numRes, string name) {
	printPopulationProfile(residents, numRes, name, 0);
}

// Prints profile similar to that of ODEINTEGRATION
void PopulationDynamics::printPopulationProfile(const double residents[], int numRes, string name, bool detail) {
	solveResidentEquilibrium(residents, numRes, 0);
	ofstream OutputFile((name + ".txt").c_str()); if (!OutputFile) {std::cerr << "File not opened" << std::endl; hold(); exit(EXIT_FAILURE);}
	int i, numberOfVariables = numRes * TRAIT_DIM, precision = 6;

	// PRINT HEADER
	string white_space(precision + 3, ' ');
	OutputFile << "%T\tNumRes\t";
	for (i = 0; i < numberOfVariables; i++) OutputFile << "y"  << i << white_space << " \t";
	for (i = 0; i < numRes; i++)			OutputFile << "p"  << i << white_space << " \t";
	OutputFile << std::endl;

	OutputFile << std::scientific << "0\t" << numRes << "\t";
	for (i = 0; i < numberOfVariables; i++) OutputFile << std::setprecision(precision) << residents[i] << "\t";
	for (i = 0; i < numRes; i++) OutputFile << getPopulationSize(i) << "\t";

	OutputFile << std::endl;

	if (detail)
		myFitnessSolver->printPopulationProfile(residents, populationSize, numRes, name, 1);
}


vector<vector<double> > PopulationDynamics::getPopulationProperties(const double residents[], int numRes) {
//	solveResidentEquilibrium(residents, numRes,0);
	return myFitnessSolver->getPopulationProperties();
}

// Invasion fitness
double PopulationDynamics::S(const double mutant[], const double residents[], int numRes)      // Invasion fitness
{
	if (numRes == 0) return mutantFitness(mutant, 0);
	solveResidentEquilibrium(residents, numRes, 0);
	return mutantFitness(mutant, 1);
}


// Calculated invasion fitness, like S, but assuming residents already set
double PopulationDynamics::mutantFitness(const double mutant[], bool comp) {
	double R = myFitnessSolver->mutantFitness(mutant, comp);
	if (!(R >= 0) && !(R < 0))
	{cout << "POPDYN:\tNaN in mutant fitness" << endl; exit(1);}

	if (R == 0) return -200.00;
	return log10(R);
}

void PopulationDynamics::setResidentState(const double residents[], const double X[], int numRes) {
	int i;
	if (numRes != numberOfResidents)
	{
		residentTraits = M2Dd_resize(residentTraits, numberOfResidents, numRes, TRAIT_DIM);
		populationSize = Vd_resize(populationSize, numberOfResidents, numRes);
		logPopulationSize = Vd_resize(logPopulationSize, numberOfResidents, numRes);
		numberOfResidents = numRes;
	}
	// Copy trait values and reset popsize if necessary
	for (i = 0; i < numberOfResidents * TRAIT_DIM; i++)
		residentTraits[0][i] = residents[i];
	for (i = 0; i < numberOfResidents; i++)
	{
		populationSize[i] = X[i];
		logPopulationSize[i] = log10(X[i]);
	}
	residentFitness(3);
}

// Set trait values for resident population & find equilibrium popn sizes
// ResetResidents is a flag calling for solution to be recalculated even if trait values are the same
// B/c some other parameter may have changed
void PopulationDynamics::solveResidentEquilibrium(const double residents[], int numRes, int resetResidents) {
	int i, j, flag = 1;
	if (numRes == numberOfResidents)
	{
		for (i = 0; i < numberOfResidents * TRAIT_DIM; i++)
			if (residents[i] != residentTraits[0][i])
				flag = 0;
		if (flag && !resetResidents && !resolveEquilibrium) return;  // Residents same as previous
	}
	else    // Reallocate memory if needed
	{
		residentTraits = M2Dd_resize(residentTraits, numberOfResidents, numRes, TRAIT_DIM);
		populationSize = Vd_resize(populationSize, numberOfResidents, numRes);
		logPopulationSize = Vd_resize(logPopulationSize, numberOfResidents, numRes);

		// Set population size of new strategies to lowest value
		for (i = numberOfResidents; i < numRes; i++)
			setPopulationSize(i, EXTINCT);

		numberOfResidents = numRes;

	}

	// Copy trait values and reset popsize if necessary
	for (i = 0; i < numberOfResidents; i++) {
		flag = 0;
		// If trait values change, reset pop size and update traits
		for (j = 0; j < TRAIT_DIM; j++)
			if (residentTraits[0][i * TRAIT_DIM + j] != residents[i * TRAIT_DIM + j]) {
				flag = 1;
				residentTraits[0][i * TRAIT_DIM + j] = residents[i * TRAIT_DIM + j];
			}
		// Reset popsize if needed
		if (populationSize[i] <= EXTINCT || resetPopulationSize)
		{
			setPopulationSize(i, EXTINCT);
		}
	}

	// Solve demographic equilibrium
	int n,  count = 0, max_count = 1;
	flag = 1;

	// Check if empty environment
	if (numberOfResidents == 0) {
		residentFitnessData.resize(0);
		return;
	}

	while (flag != 0 && count < max_count)
	{
		count++;
		// Iterate residents number of times to get close to solution
		n = 0;
		while (n < maxPopulationIterations && flag)
		{
			flag = 0; n++;
			if (resetEbtCohorts) {
				cout << "resetEbtCohorts" << endl;
//				for(i=0; i< numberOfResidents; i++)  setPopulationSize(i,EXTINCT);
				residentFitnessData = residentFitness((useDynamicCohorts) ? 1 : 0);
				resetEbtCohorts = 0;
			}
			else	residentFitnessData = residentFitness((useDynamicCohorts) ? 3 : 2);
			for (i = 0; i < numberOfResidents; i++)
			{
				if (!(residentFitnessData[i] >= 0) && !(residentFitnessData[i] < 0)) {
					cout << "POPDYN:\tNaN in resident fitness" << endl; exit(1);
				}
				if (displayWorkingDetails > 1)
					cout << "POPDYN:\tSolveRes\t" << i << "\t" << setprecision(8) << populationSize[i] << "\t" << log10(residentFitnessData[i]) << endl;

				if ((populationSize[i]*residentFitnessData[i] > EXTINCT && fabs(log10(residentFitnessData[i])) > Demographic_ROOT) || !stopOnDemographicRoot)
					flag = 1;
				setPopulationSize(i, populationSize[i]*residentFitnessData[i]);
			}
		}

		// Print
		if (displayWorkingDetails > 0)
		{
			cout << "POPDYN:\tSolveRes\t" << ((flag == 1) ? "FAIL" : "SUCCESS");
			for (i = 0; i < numberOfResidents; i++)
				cout << setprecision(8) << "\t" << log10(residentFitnessData[i]);
			cout << endl;
		}

		// CHECK FOR SMALL NON-VIABLE POPN in single spp popn
		if (numberOfResidents == 1 && populationSize[0] < EXTINCT)
			flag = 0;

	}
	resolveEquilibrium = 0;
}

// Calculates fitness (as R) of N resident species in competition with each other at given population size
vector<double> PopulationDynamics::residentFitness(int useDynamicOptions) {
	vector<double> Fitness = myFitnessSolver->residentFitness(residentTraits[0], populationSize, numberOfResidents, useDynamicOptions);
	return Fitness;
}

// Wrapper for resident fitness with no arguments
vector<double> PopulationDynamics::residentFitness(void) {
	return residentFitness((useDynamicCohorts) ? 3 : 2);
}

void printResToFile(const double Res[], int numRes, vector<int> traitList, ofstream &OutFile) {
	for (int n = 0; n < numRes; n++)
	{
		OutFile << "% ";
		for (int i = 0; i < (int)traitList.size(); i++) OutFile << setprecision(6) << "Res[" << n << "][" << traitList[i] << "]=" << Res[n * TRAIT_DIM + traitList[i]] << ";\t";
		OutFile << endl;
	}
}

void printRes(const double Res[], int numRes, vector<int> traitList) {
	for (int n = 0; n < numRes; n++)
	{
		for (int i = 0; i < (int)traitList.size(); i++) cout << setprecision(6) << "Res[" << n << "][" << traitList[i] << "]=" << Res[n * TRAIT_DIM + traitList[i]] << ";\t";
		cout << endl;
	}
}
