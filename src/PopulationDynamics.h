//*******************************************************************************************/
//  PopulationDynamics.h
//  Created by Daniel Falster, 2006-2010
//*******************************************************************************************/
#ifndef __POPULATION_DYNAMICS_H_
#define __POPULATION_DYNAMICS_H_

#include "FitnessSolver.h"
#include "gsl_random_numbers.h"

// Structure used by Adaptive dynamics routines - not used in PopulationDynamicsClass

// For solving population equilibrium
#define EXTINCT 1e-2                // Population size at extinction
#define LOG_EXTINCT log10(1e-2)
#define Demographic_ROOT 5e-4
#define PRINT 1

void printRes(const double mutant[], int numRes, vector<int> traitList);
void printResToFile(const double Res[], int numRes, vector<int> traitList, ofstream &OutFile);

class PopulationDynamics {
public:
	PopulationDynamics();

	// Setup
	void setup(FitnessSolver* myFitnessSolverPointer);
	void setDisplayWorkingDetails(int option);						// 1= on. 0 = off

	// Calculations of fitness and derivatives of fitness
	double S(const double mutant[], const double residents[], int numRes);                  // Invasion fitness
	vector<double> residentFitness(void);						// Wrapper called by gsl root solving for population dynamics
	vector<double> residentFitness(int useDynamicOptions);

	// Interaction
	void setPopulationSize(int i, double new_value);
	void setLogPopulationSize(int i, double new_logN);
	void resetFitnessSolver(void);
	void setResidentState(const double residents[], const double X[], int numRes);

	// Queries
	int getnumberOfResidents(void) const;
	double getPopulationSize(int i) const;
	vector<vector<double> > getPopulationProperties(const double residents[], int numRes);
	void displayResidentTraits(void);
	void printPopulationProfile(const double residents[], int numRes, string name);
	void printPopulationProfile(const double residents[], int numRes, string name, bool detail);

	bool isExtinct(int i);
	bool isExtinct(int i, double threshold);
	bool checkForExtinctPopns(void);
	bool checkForExtinctPopns(double threshold);
	int removeResident(int index);

	// Flags to control behaviour - can be accessed and modified directly
	bool resetPopulationSize;
	bool useDynamicCohorts;
	bool resetEbtCohorts;
	bool resolveEquilibrium;
	int displayWorkingDetails;
	int maxPopulationIterations;
	int stopOnDemographicRoot;

	double maxPopulationSize; // Upper limit for population size - stops solver asking for excessively large values
	double maxLogPopulationSize;


private:
	double mutantFitness(const double mutant[], bool comp);    // Calculated invasion fitness, like S, but assuming residents already set
	void solveResidentEquilibrium(const double residents[], int numRes, int resetResidents);

	// Data Storage
	int numberOfResidents;
	double *populationSize, *logPopulationSize;
	vector<double> residentFitnessData;
	FitnessSolver* myFitnessSolver;
	double** residentTraits;

};

#endif
