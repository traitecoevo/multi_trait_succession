//*******************************************************************************************/
//  AdaptiveDynamics_2D.h
//  Created by Daniel Falster, 2006-2010
//*******************************************************************************************/

#ifndef __ADAPT_2D_H_
#define __ADAPT_2D_H_

#include "PopulationDynamics.h"
#define VIABLE 0.01	       // Minimum fitness in virgin environment to consider viable

class AdaptiveDynamics_2D
{
public:
	AdaptiveDynamics_2D();
	~AdaptiveDynamics_2D();
	// Configuration *********************************************************************************************
	void setup(PopulationDynamics* myPopulationDynamicsPointer);
	void setDisplayWorkingDetails(int option);						// 1= on. 0 = off
	// Interaction *********************************************************************************************
	bool isViableTraitCombination(const double mutant[]) const;
	// Tests if trait combination has positive fitness in virgin environment
	double findTraitViabilityLimit(double* startTraits, double lowerBound, double upperBound);
	// Returns value of trait with index focalTraitIndex that gives zero fitness
	// POSTCONDITION: startTraits stored in residentTraits
	void followTraitViabilityContour2D(double* startTraits, double traitRange[2][3], string name);
	bool displayWorkingDetails;
	PopulationDynamics* popDyn;
	double* residentTraits;
	int indexTrait1, indexTrait2;  // Focal traits for analysis
	bool isSetup;

private:

	// Functions used for following contours *********************************************************************************************
	void followContourBothDirections(int function_option, double contour_height,	double X_init, double Y_init, double dir_init,
	                                 double dir_multiplier, double step_size, double X_lim[2], double Y_lim[2], string name);
	void followContourOneDirection(int function_option, double contour_height, double X_init, double Y_init, double& dir_init,
	                               double dir_multiplier, double xy_step, double X_lim[2], double Y_lim[2], string name);
	double evaluateContourHeight(int option, double X, double Y);

};

#endif
