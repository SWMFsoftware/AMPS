//$Id: ChargeExchange_LindsayStebbings.cpp,v 1.1 2015/03/31 21:38:21 dborovik Exp $
#include "ChargeExchange.h"

double ChargeExchange::LindsayStebbings::LifeTime(int      spec,
                                                  double* vParticle,
                                                  double* vPlasma,
                                                  double   PlasmaTemperature,
                                                  double   PlasmaNumberDensity) {
  // for notation see Lindsay & Stebbings, JGR, vol. 110, A12213, 2005
  double a1 = 4.15, a2 = 0.531, a3 = 67.3, v_th, v_rel,v_rel2, energy, omega=0.0, sigma;

  // this model is for atomic hydrogen and three charge exchange species only
  if(spec != _H_SPEC_ && spec != _H_ENA_V1_SPEC_ && spec != _H_ENA_V2_SPEC_ && spec != _H_ENA_V3_SPEC_) return 1E+100;


  // calculating V_rel according to eq. 8 of Heerkhuisen et al. 2006
  v_th = sqrt(2.0 * Kbol * PlasmaTemperature / _MASS_(_H_));
  for(int idim = 0; idim < 3; idim++ )
    omega += (vParticle[idim]-vPlasma[idim])*(vParticle[idim]-vPlasma[idim]);
  omega = sqrt(omega) / v_th;
  if(omega > 1E-2)
    v_rel = v_th* ( exp(-omega*omega)/sqrtPi + (omega + 0.5/omega)*erf(omega));
  else
    // Taylor expansion in vicinity of omega = 0 to avoid divison by zero;
    // relative error < 1E-4
    v_rel = v_th* 2.0 * (1.0 + omega*omega/3.0) / sqrtPi;

  v_rel2=v_rel*v_rel;

  if(v_rel2 < 1E-8) return 1E+20; // to avoid division by zero
  energy = 0.5 * _MASS_(_H_) * v_rel2 /eV2J / 1000; // energy in keV
  double c;
  if (energy > 4.0)
    c = pow(1 - exp( -a3/energy ), 4.5);
  else
    c = 1.0; // error < 1E-6
  if(energy > 1E-10)
    sigma = (a1 - a2*log(energy)) * (a1 - a2*log(energy)) * c * 1E-20; // m^2
  else
    // to avoid computing log() in vicinity of energy = 0:
    // substitute log(energy) with (-10*log(10))
    sigma = (a1 + a2*23.026) * (a1 + a2*23.026) * 1E-20; // m^2
  return 1.0 / (PlasmaNumberDensity * v_rel * sigma);
}
