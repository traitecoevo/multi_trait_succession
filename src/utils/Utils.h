//*******************************************************************************************/
//  Utils.h
//  Created by Daniel Falster, 2006-2010
//*******************************************************************************************/
#ifndef __UTILS_H_
#define __UTILS_H_

/*Standard libraries used throughout*/

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip> // for formatting
#include <math.h>
#include <cstdlib>

#include<list>
#include <stack>
#include <queue>

// For string manipulations
#include <sstream>
#include <stdexcept>

/*Libraries used for directory manipulations*/
#include <sys/stat.h>
#include <unistd.h>
#ifndef _WIN32
#define mkdir(a) mkdir(a, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH)
#endif

// Functions used in directory operations
bool FileExist(const char* FileName);
bool IsDirectory(const char* FileName);
void go_to_dir(std::string dir);

// Maths
double dfdx(double f_xhigh, double f_xlow, double h);
double df2d2x(double f_xhigh, double f_x, double f_xlow, double h);
bool equals(double d1, double d2, double precision);
double sign(double temp);


// Formatting
double set_decimal(double T, int precis);

void print_file(std::string name, double** out, int rows, int cols);
double** vec2D_to_matrix(std::vector<std::vector<double> >& Vec, int & rows, int &cols);
void hold(void);


class BadConversion : public std::runtime_error
{
	public: BadConversion(const std::string& s): std::runtime_error(s)
		{ }
};

inline std::string stringify(double x) {
	std::ostringstream o;
	if (!(o << x)) throw BadConversion("stringify(double)");
	return o.str();
}

#endif
