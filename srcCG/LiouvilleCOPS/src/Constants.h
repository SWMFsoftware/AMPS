//physical constants (SI units)
#ifndef CONSTANTS_H
#define CONSTANTS_H

namespace Const {
  const double Pi=3.14159265358979323846;
  const double SqrtPi=1.77245385090551602730;

  const double kBoltzmann=1.380649E-23;     //J/K (CODATA 2018, exact)
  const double AMU=1.66053906660E-27;       //kg  (CODATA 2018)
  const double AU=149597870700.0;           //m   (IAU 2012, exact)

  //species considered by the srcCG source model
  const int nSpecies=2;
  const int H2O=0,CO2=1;

  const char SpeciesName[nSpecies][8]={"H2O","CO2"};
  const double SpeciesMass[nSpecies]={18.01528*AMU,44.0095*AMU}; //molecular masses
}

#endif
