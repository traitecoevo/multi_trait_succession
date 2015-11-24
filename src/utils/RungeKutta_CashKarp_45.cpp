
#include "RungeKutta_CashKarp_45.h"

// For RK integration
#define TINY 1e-30
#define SAFETY 0.9
#define PGROW -0.2
#define PSHRINK -0.25
#define ERRCON 1.89e-4

RungeKutta_CashKarp_45::RungeKutta_CashKarp_45(){
	// Error_type
	error_type = 0;

	// Constants determined by Cash & Karp
	a1 =0.2;			a2=0.3;				a3=0.6;			a4=1.0;					a5= 7.0/8.0;
    b21=0.2;
	b31=3.0/40.0;		b32=9.0/40.0;
	b41=0.3;			b42=-0.9;			b43=1.2;
	b51=-11.0/54.0;	b52=2.5;			b53=-70.0/27.0;	b54=35.0/27.0;
	b61=1631.0/55296.0; b62=175.0/512.0;	b63=575.0/13824.0;	b64=44275.0/110592.0;  b65=253.0/4096.0;
	c1=37.0/378.0;							c3=250.0/621.0;	c4=125.0/594.0;	   c6=512.0/1771.0;
	dc1=37.0/378.0-2825.0/27648.0;			dc3=250.0/621.0-18575.0/48384.0; dc4=125.0/594.0-13525.0/55296.0; dc5=-277.0/14336.0;  dc6=512.0/1771.0-0.25;
    }

double RungeKutta_CashKarp_45::return_x_after_step(const double x, double h, int step_type)
	 {
	 switch(step_type){
		  case(0): return x;
		  // Rk_level 1
		  case(1): return (x + a1*h);
		  // Rk_level 2
		  case(2): return (x + a2*h);
		  // Rk_level 3
		  case(3): return (x + a3*h);
		  // Rk_level 4
		  case(4): return (x + a4*h);
		  // Rk_level 5
		  case(5): return (x + a5*h);
		  // Rk_level 6 - final step
		  case(6): return (x + h);
		  }
	return 0.0;
	}


double RungeKutta_CashKarp_45::return_next_y(const double y, double h, std::vector<double>& f, int step_type)
	 {
	 switch(step_type){
		  // Euler step
		  case(0): return (y + h*f[0]);
		  // Rk_level 1
		  case(1): return (y + h*(b21*f[1]));
		  // Rk_level 2
		  case(2): return (y + h*(b31*f[1]+b32*f[2]));
		  // Rk_level 3
		  case(3): return (y + h*(b41*f[1]+b42*f[2]+b43*f[3]));
		  // Rk_level 4
		  case(4): return (y + h*(b51*f[1]+b52*f[2]+b53*f[3]+b54*f[4]));
		  // Rk_level 5
		  case(5): return (y + h*(b61*f[1]+b62*f[2]+b63*f[3]+b64*f[4]+b65*f[5]));
		  // Rk_level 6 - final step
		  case(6): {
				   // Difference for final step based on weight rates of change
				   dy = h*( c1*f[1]+  c3*f[3]+  c4*f[4]+ c6*f[6]);
				   // Difference between 4th orf & 5th embedded steps
				   yerr = h*(dc1*f[1]+ dc3*f[3]+ dc4*f[4]+ dc5*f[5]+dc6*f[6]);
				   // Scaling of error
				   switch(error_type){
						  default:case(0): yscal = fabs(y) + fabs(dy) + TINY; break;
						  }
				   // Store estimate of error in f[0] and yerr
				   f[0] =  fabs(yerr/yscal);
				   f[7] =  yerr;
				   f[8] = dy/h;
		       return (y + dy);
				   }
		  }
	return 0.0;
	}

 double RungeKutta_CashKarp_45::grow_step(double dt, double errmax)
	{
	if(errmax > ERRCON)
		return SAFETY*dt*pow(errmax, PGROW);
	else
		return 5.0*dt;
     }

 double RungeKutta_CashKarp_45::shrink_step(double dt, double errmax)
	{
	return std::max(SAFETY*dt*pow(errmax, PSHRINK), 0.1*dt);
	}
