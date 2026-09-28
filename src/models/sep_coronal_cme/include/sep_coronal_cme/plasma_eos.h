#ifndef SEP_CORONAL_CME_PLASMA_EOS_H
#define SEP_CORONAL_CME_PLASMA_EOS_H

#include "sep_status.h"

#include <string>
#include <vector>

namespace SEP { namespace CoronalCME {

struct IonSpecies {
  std::string id;
  double abundancePerProton = 0.0;
  double chargeNumber = 0.0;
  double massKg = 0.0;
  double temperatureK = 0.0;
};

struct PlasmaState {
  double protonNumberDensityM3 = 0.0;
  double electronNumberDensityM3 = 0.0;
  double massDensityKgM3 = 0.0;
  double pressurePa = 0.0;
  double soundSpeedMPerS = 0.0;
  double alfvenSpeedMPerS = 0.0;
  double fastSpeedMPerS = 0.0;
};

// Converts measured electron density using explicit ion abundance and charge
// state. No implicit mean molecular weight or alpha abundance is permitted.
Core::Result<PlasmaState> EvaluatePlasmaFromElectronDensity(
    double electronDensityM3, double electronTemperatureK,
    const std::vector<IonSpecies>& ions, bool includeElectronMass,
    double gammaAdiabatic, double magneticFieldT,
    double propagationCosine = 1.0);

} }  // namespace SEP::CoronalCME
#endif
