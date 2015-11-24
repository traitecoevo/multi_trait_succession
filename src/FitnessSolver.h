//*******************************************************************************************/
//  FitnessSolver.h
//  Created by Daniel Falster, 2006-2010
// Fitness is assumed to be Reproductive ratio, i.e. where R=1 is equilibrium
//*******************************************************************************************/

#ifndef __FITNESS_SOLVER_H_
#define __FITNESS_SOLVER_H_

#define TRAIT_DIM 4
#include <vector>
#include "Utils.h"

using namespace std;

class FitnessSolver
{
public:
	FitnessSolver();
	virtual vector<double> residentFitness(const double residentTraits[], const double X[], int N, int useDynamicOptions);
	// Calculates fitness of N resident species in competition with each other, and with population size X
	// Precondition: none
	// Precondition: resident condition stored in solver, so that mutant fitness can be calculated
	virtual double mutantFitness(const double Traits[], bool with_competition);
	// Calculates mutant fitness in competition with residents. If with_competition == 0 then calculates fitness in virgin environment
	// Precondition: appropriate resident condition already stored in solver, achieved by running residentFitness
	// Postcondition: none
	virtual vector<vector<double> > getPopulationProperties(void);
	// Returns summary of resident population. Columns contain properties of interest, with 1 spp per row, and community total / averages on last row.
	// Precondition: appropriate resident condition already stored in solver, achieved by running residentFitness
	// Postcondition: none
	virtual void printPopulationProfile(const double residentTraits[], const double X[], int N, string DIR);
	// Prints details of resident population to file in directory DIR.
	// Precondition: appropriate resident condition already stored in solver, achieved by running residentFitness
	// Postcondition: none
	virtual void printPopulationProfile(const double residentTraits[], const double X[], int N, string DIR, bool detail);
	// Prints details of resident population to file in directory DIR.
	// Precondition: appropriate resident condition already stored in solver, achieved by running residentFitness
	// Postcondition: none
	virtual void remove_spp(int index);
	// Remove species from solver
	// Precondition: 0< index < numberOfResidentSpp
	// Postcondition: numberOfResidentSpp decreased by 1
private:
	int numberOfResidentSpp;
	bool residentsAreSet;
};

#endif
