//*******************************************************************************************/
//  Site.h
//  Created by Daniel Falster, 2006-2010
//*******************************************************************************************/
#ifndef __SITE_H_
#define __SITE_H_

/*Class used to store parameters. Introduced as global variable in main*/
class site {
public:
  site();                // constructor
  void set_type(int T);
  // Functions used for patch age_distribution

  double solve_patchage_dist(double site_mean);
  double disturbance_rate(double age);
  double patch_age_density_freq(double age);

  double Pi(double age);
  double patch_weight(double time_start, double time);
  double p0, lam, a_mean;

private:
  int type;  //0 = no meta population; 1 = Weibull distribution; 2= Exponential distribution


};
#endif
