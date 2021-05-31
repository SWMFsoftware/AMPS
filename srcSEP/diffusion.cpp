//SEP diffusion models



#include "sep.h"


//========= Roux2004AJ (LeRoux-2004-AJ) =============================
double SEP::Diffusion::Roux2004AJ::D_mu_mu(double mu) {
  static double c=7.0*Pi/8.0*2.0E3/(0.01*_AU_)*(0.2*0.2);

  return c*(1.0-mu*mu);
} 

double SEP::Diffusion::Roux2004AJ::dD_mu_mu_mu(double mu) {
  static double c=7.0*Pi/8.0*2.0E3/(0.01*_AU_)*(0.2*0.2);

  return -2.0*c*mu;
} 

