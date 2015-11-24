#include "FitnessSolver.h"

FitnessSolver::FitnessSolver() {
	residentsAreSet = 0;
}
// Calculates fitness (as R) of N resident species in competition with each other, and with population size X
// Precondition: none
// Precondition: resident condition stored in solver, so that mutant fitness can be calculated
vector<double> FitnessSolver::residentFitness(const double residentTraits[], const double X[], int N, int useDynamicOptions) {
	residentsAreSet = 1;
	numberOfResidentSpp = N;
	vector<double> Output(numberOfResidentSpp, 1);
	return Output;
}

// Calculates mutant fitness in competition with residents. If with_competition == 0 then calculates fitness in virgin environment
// Precondition: appropriate resident condition already stored in solver, achieved by running residentFitness
// Postcondition: none
double FitnessSolver::mutantFitness(const double Traits[], bool with_competition) {
	if (!residentsAreSet) {cerr << "Mutant fitness requested without running resident s first" << endl; hold(); exit(EXIT_FAILURE);}
	return 1;
}

// Returns summary of resident population. Columns contain properties of interest, with 1 spp per row, and community total / averages on last row.
// Precondition: appropriate resident condition already stored in solver, achieved by running residentFitness
// Postcondition: none
vector<vector<double> > FitnessSolver::getPopulationProperties(void) {
	int cols = 10;
	if (!residentsAreSet) {cerr << "Population properties requested without running resident s first" << endl; hold(); exit(EXIT_FAILURE);}
	vector<vector<double> > output(numberOfResidentSpp + 1, vector<double>(cols, 0.0));
	return output;
}


// Prints details of resident population to file in directory DIR.
// Precondition: appropriate resident condition already stored in solver, achieved by running residentFitness
// Postcondition: none
void FitnessSolver::printPopulationProfile(const double residentTraits[], const double X[], int N, string DIR) {
	if (!residentsAreSet) {cerr << "Population print requested without running resident s first" << endl; hold(); exit(EXIT_FAILURE);}

	ofstream OutputFile((DIR + "/output.txt").c_str());
	OutputFile << "A little output" << endl;
}

void FitnessSolver::printPopulationProfile(const double residentTraits[], const double X[], int N, string DIR, bool detail) {
	printPopulationProfile(residentTraits, X, N, DIR);
}


void FitnessSolver::remove_spp(int index) {
	return;
}
