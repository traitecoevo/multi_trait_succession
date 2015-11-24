#include "AdaptiveDynamics_2D.h"
#include "Collect.h"
#include "MatrixCPP.h"
#include "gsl_root_class.h"

#define PI 3.141592654

// For ROOT SOLVING ALGORITHM
#define ITMAX 100	// Maximum allowed number of iterations.
#define EPS 3.0e-8	// Machine floating-point precision.
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

double ViabilityWrapper_f(double X, void *params);

AdaptiveDynamics_2D::AdaptiveDynamics_2D() {
	isSetup = 0;
	residentTraits = Vd_alloc(TRAIT_DIM);
	indexTrait1 = 0; indexTrait2 = 1;
	displayWorkingDetails = 1;
};

AdaptiveDynamics_2D::~AdaptiveDynamics_2D() {
	delete[] residentTraits;
}

void AdaptiveDynamics_2D::setup(PopulationDynamics* myPopulationDynamicsPointer) {
	popDyn = myPopulationDynamicsPointer;
	popDyn->resetPopulationSize = 0;
	popDyn->displayWorkingDetails = 1;
	isSetup = 1;
}

void AdaptiveDynamics_2D::setDisplayWorkingDetails(int option) {
	displayWorkingDetails = option;
}


bool AdaptiveDynamics_2D::isViableTraitCombination(const double mutant[]) const
{
	return (popDyn->S(mutant, mutant, 0) > VIABLE);
}

double AdaptiveDynamics_2D::findTraitViabilityLimit(double* startTraits, double lowerBound, double upperBound) {
	if (displayWorkingDetails)   cout << "AD2::\tIn Solve Viability (" << lowerBound << " " << upperBound << ")" << flush;
	for (int i = 0; i < TRAIT_DIM; i++)  residentTraits[i] = startTraits[i];

	gsl_root_class gslRootSolver;
	double y = gslRootSolver.single_bracket(&ViabilityWrapper_f, this, lowerBound, upperBound, 1E-5, 50, 0);

	if (displayWorkingDetails)   cout << "\t" << y << endl;
	return y;
}

// Wrapper function used in root solving viability
double ViabilityWrapper_f(double X, void *params) {
	AdaptiveDynamics_2D * AD2 = (AdaptiveDynamics_2D *) params;

	AD2->residentTraits[AD2->indexTrait1] = X;
	return AD2->popDyn->S(AD2->residentTraits, AD2->residentTraits, 0) - VIABLE;
}

void AdaptiveDynamics_2D::followTraitViabilityContour2D(double* startTraits, double traitRange[2][3], string name) {
	if (!isSetup) {cerr << "fitness function not set" << endl; hold(); exit(EXIT_FAILURE);}
	cout << "creating TVP with contour plot" << endl;

	for (int i = 0; i < TRAIT_DIM; i++)  residentTraits[i] = startTraits[i];
	popDyn->S(startTraits, startTraits, 1); // call first to make sure PopulationDynamics is set
// Find starting location along bottom of x axis
	residentTraits[indexTrait2] = traitRange[1][0];
	residentTraits[indexTrait1] = traitRange[0][0];
	int b1 = isViableTraitCombination(residentTraits);   // Check viability at outside limit
	cout << residentTraits[indexTrait1] << "\t" << residentTraits[indexTrait2] << "\t" << b1 << endl;

	residentTraits[indexTrait1] = 0.5 * (traitRange[0][0] + traitRange[0][1]);
	int b2 = isViableTraitCombination(residentTraits);   // Check viability at centre
	cout << residentTraits[indexTrait1] << "\t" << residentTraits[indexTrait2] << "\t" << b2 << endl;
	if (b2 & !b1)
	{
		// Find starting point
		residentTraits[indexTrait1] = findTraitViabilityLimit(residentTraits, traitRange[0][0], 0.5 * (traitRange[0][0] + traitRange[0][1]));
		double dir = 0.25;
		// Determine step size
		double step_size = min(fabs(traitRange[0][0] - traitRange[0][1]), fabs(traitRange[1][0] - traitRange[1][1])) / 30.0;
		// Follow contour
		followContourOneDirection(0, 0, residentTraits[indexTrait1], residentTraits[indexTrait2], dir, 20, step_size, traitRange[0], traitRange[1], name);
	}
	else
		cout << "bad starting region" << endl;
}

// Wrapper for follow contour function - follows contours in both directions, joins output into single file and deletes extra files
void AdaptiveDynamics_2D::followContourBothDirections(int function_option, double contour_height, double X_init, double Y_init, double dir_init, double dir_multiplier, double step_size, double X_lim[2], double Y_lim[2], string name) {
	popDyn->resetPopulationSize = 1;
	// Follow contour in one direction
	double dir = dir_init;
	followContourOneDirection(function_option, contour_height, X_init, Y_init, dir, dir_multiplier, step_size, X_lim, Y_lim, name + "-contour1");
	// Follow contour in other direction
	dir += 0.5;
	followContourOneDirection(function_option, contour_height, X_init, Y_init, dir, 1, step_size, X_lim, Y_lim, name + "-contour2");

	// Append files output together
	vector<string> Data; string temp;
	ifstream IF((name + "-contour1.txt").c_str());
	while (getline (IF, temp))	Data.insert (Data.begin(), temp); // Store data from first file in reverse order
	IF.close(); IF.open((name + "-contour2.txt").c_str());
	getline (IF, temp); // Discard first line which is same as previous file
	while (getline (IF, temp))	Data.insert (Data.end(), temp); // Store data from 2nd file in order
	IF.close();
	ofstream OF((name + "-contour.txt").c_str());
	for (int i = 0; i < (int)Data.size(); i++) OF << Data[i] << endl; // Print all to new file
	OF.close();
	remove ((name + "-contour1.txt").c_str()); remove((name + "-contour2.txt").c_str()); // Delete files

	// Append work files together
	OF.open((name + "-contour1-work.txt").c_str(), ios::app);
	IF.open((name + "-contour2-work.txt").c_str());
	while (getline(IF, temp))	OF << temp << endl;
	OF.close(); IF.close();
	remove((name + "-contour2-work.txt").c_str());
}

// Function follows contour using modified version of Ulf's searchlight algorithm
// Assumes that initial points lie on on contour
// Finds angle of contour  (saves this in dir_init), then narrows angle down to desired accuracy from
// initial search window which is dir_multiplier wider (allows function to find initial angle quicker)
// Direction: 0.0 = right, 0.25 = up, 0.5 = left, 0.75 = down, 1= right
void AdaptiveDynamics_2D::followContourOneDirection(int function_option, double contour_height,
    double X_init, double Y_init, double& dir_init, double dir_multiplier,
    double xy_step, double X_lim[2], double Y_lim[2], string name) {
	int points, max_points = 300, result = 0;
	double X_old, Y_old, dir_old, X_new, Y_new, dir_new, dev_new, X_left, Y_left, dev_left,
	       X_right, Y_right, dev_right;
	bool is_valid = 1, first_step_flag = 1;
	popDyn->resetPopulationSize = 0;
	// Angular step - start much larger than desired angle, find direction then angle decreases
	double dir_step_final = 0.005; //0.005 --> window = 2 deg
	double dir_step = min(0.1, dir_step_final * max(dir_multiplier, 1.0));

	cout << "Follow " << xy_step << "\t" << dir_step << endl;
	// Open files
	ofstream OutFile1((name + "-work.txt").c_str());
	OutFile1 << "%points\t result\t X_new\t Y_new\t dir_step\t dir_new\t dev_new\t dev_left\t dev_right\textras" << endl;
	ofstream OutFile2((name + ".txt").c_str());
	if (! (OutFile1 && OutFile2) ) {std::cerr << "File not opened" << std::endl; hold(); exit(EXIT_FAILURE);}

	X_old = X_new = X_init; Y_old = Y_new = Y_init; dir_old = dir_new = dir_init; points = 0;

	// Print first point
	OutFile2 << points << "\t" << X_old << "\t" << Y_old << "\t"
	         << fabs(evaluateContourHeight(function_option, X_old, Y_old) - contour_height)
	         << "\t" << "" << endl;

	// Find other points
	while (is_valid && points < max_points) {
		dir_old = dir_new;
		if (result != 1) // Don't step when finding angle
		{X_old = X_new; Y_old = Y_new;}

		X_left = X_old + xy_step * cos(2.0 * PI * (dir_old + dir_step));
		Y_left = Y_old + xy_step * sin(2.0 * PI * (dir_old + dir_step));
		dev_left = fabs(evaluateContourHeight(function_option, X_left, Y_left) - contour_height);

		X_new = X_old + xy_step * cos(2.0 * PI * dir_old);
		Y_new = Y_old + xy_step * sin(2.0 * PI * dir_old);
		dev_new = fabs(evaluateContourHeight(function_option, X_new, Y_new) - contour_height);

		X_right = X_old + xy_step * cos(2.0 * PI * (dir_old - dir_step));
		Y_right = Y_old + xy_step * sin(2.0 * PI * (dir_old - dir_step));
		dev_right = fabs(evaluateContourHeight(function_option, X_right, Y_right) - contour_height);

		result = 2;
		while (is_valid && result == 2) {
			// Check still in valid region
			is_valid = ((X_new  > X_lim[0] && X_new  < X_lim[1])
			            && (Y_new  > Y_lim[0] && Y_new  < Y_lim[1]))
			           && ((X_right > X_lim[0] && X_right < X_lim[1]) && (Y_right > Y_lim[0] && Y_right < Y_lim[1]))
			           && ((X_left > X_lim[0] && X_left < X_lim[1]) && (Y_left > Y_lim[0] && Y_left < Y_lim[1]))
			           && (dev_new != -999 && dev_left != -999  && dev_right != -999);

			if (dev_left <= dev_new && dev_left < dev_right) {
				dir_new = dir_new + dir_step;
				X_right = X_new; Y_right = Y_new; dev_right = dev_new;
				X_new = X_left; Y_new = Y_left; dev_new = dev_left;
				X_left = X_old + xy_step * cos(2.0 * PI * (dir_new + dir_step));
				Y_left = Y_old + xy_step * sin(2.0 * PI * (dir_new + dir_step));
				dev_left = fabs(evaluateContourHeight(function_option, X_left, Y_left) - contour_height);
			}
			else if (dev_right <= dev_new && dev_right < dev_left) {
				dir_new = dir_new - dir_step;
				X_left = X_new; Y_left = Y_new; dev_left = dev_new;
				X_new = X_right; Y_new = Y_right; dev_new = dev_right;
				X_right = X_old + xy_step * cos(2.0 * PI * (dir_new - dir_step));
				Y_right = Y_old + xy_step * sin(2.0 * PI * (dir_new - dir_step));
				dev_right = fabs(evaluateContourHeight(function_option, X_right, Y_right) - contour_height);
			}
			else if (dev_new < dev_left && dev_new < dev_right) {
				if (dir_step > dir_step_final) {	// Success, decrease angle to desired level
					dir_step = max(dir_step * 0.5, dir_step_final);
					result = 1;
				}
				else {  // Angle is fine enough - print point
					result = 0;
					if (first_step_flag) { // Saves direction of first successful step
						dir_init = dir_new;
						first_step_flag = 0;
					}
				}
			}
			else
				result = 3;
			OutFile1 << points << "\t" << result << "\t" << X_new << "\t"	<< Y_new << "\t" << dir_step << "\t" << dir_new << "\t" << dev_new << "\t" << dev_left << "\t" << dev_right << endl;
		}
		if (result == 0) {		// PRINT successful points
			points++;
			OutFile2 << points << "\t" << X_new << "\t" << Y_new << "\t" << dev_new << "\t"
			         << "" << endl;
		}
	}
	OutFile1.close();
	OutFile2.close();
}

double AdaptiveDynamics_2D::evaluateContourHeight(int option, double X, double Y) {
	switch (option) {
	case (0): // Viability of Strategy (fitness in virgin env)
		residentTraits[indexTrait1] = X; residentTraits[indexTrait2] = Y;
		return popDyn->S(residentTraits, residentTraits, 0) - VIABLE;
	default:	return 0.0;
	}
	return 0.0;
}

