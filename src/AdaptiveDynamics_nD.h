//*******************************************************************************************
//  AdaptiveDynamics_nD.h
//  Created by Daniel Falster, 2006-2010

//*******************************************************************************************

#ifndef __ADAPT_ND_H_
#define __ADAPT_ND_H_

#include "PopulationDynamics.h"
#include "MatrixCPP.h"

class AdaptiveDynamics_nD {
public:
	AdaptiveDynamics_nD();
	void setup(PopulationDynamics* myPopulationDynamicsPointer);    // Pass pointer to population dynamics solver
	void setDisplayWorkingDetails(int option);						// 1= on. 0 = off

	void runStochasticModelNonEquil(double startTraits[], int numRes, double mutationVariance, double traitRange[2][3], double mutationRate, double probLongDistanceDispersal, double start_time, double end_time, double printFreq, string name);
	void runStochasticModelNonEquil(double startTraits[], int numRes, vector<double> mutationVariance, double traitRange[2][3], double mutationRate, double probLongDistanceDispersal, double start_time, double end_time, double printFreq, string name);

	void inputTraitsFromStochFile(double TraitsToInput[], int &numRes,  double &Time, string filename);
	void inputTraitsFromStochFile(double TraitsToInput[], int &numRes, double &Time, string filename, int line);

	// Fitness landscape functions copied from AD2
	void outputFitnessLandscape1D(const double residents[], int numRes, double* startTraits, int focalTraitIndex, double traitRange[3], string name, bool print_res);
	vector<vector<double> > getSliceOfFitnessLandscape1D(const double residents[], int numRes, const double startTraits[], int focalTraitIndex, double traitRange[3]);
	void outputFitnessLandscape2D(const double residents[], int numRes, double traitRange[2][3], string name, bool print_res);

private:
	PopulationDynamics* popDyn;
	bool isSetup;
	int numberOfResidents;

	int displayWorkingDetails;
	int precision;

};

#endif
