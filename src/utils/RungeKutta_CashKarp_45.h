//*******************************************************************************************/
// RungeKutta_CashKarp_45
//
// (C) Daniel Falster, last update 16/07/2010
//*******************************************************************************************/
#ifndef __RK45CK_H_
#define __RK45CK_H_

#include <math.h>
#include <vector>

// Constants for Runge Kutta integration steps
class RungeKutta_CashKarp_45{
public:
	RungeKutta_CashKarp_45();
	double return_next_y(const double y, double h, std::vector<double>& f, int step_type);
	double return_x_after_step(const double x, double h, int step_type);
	double grow_step(double dt, double errmax);
	double shrink_step(double dt, double errmax);

	int error_type;
	double a1, a2, a3, a4, a5;

private:
	double dy, yerr, yscal;
	double b21, b31, b32, b41, b42, b43, b51, b52, b53, b54, b61, b62, b63, b64, b65;
	double c1, c3, c4, c6;
	double dc1, dc3, dc4, dc5, dc6;
};

#endif
