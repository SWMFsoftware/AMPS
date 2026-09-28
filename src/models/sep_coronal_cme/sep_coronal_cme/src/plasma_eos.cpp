#include "sep_coronal_cme/plasma_eos.h"
#include "sep_coronal_cme/constants.h"

#include <algorithm>
#include <cmath>

namespace SEP { namespace CoronalCME {
Core::Result<PlasmaState> EvaluatePlasmaFromElectronDensity(double ne,double te,
    const std::vector<IonSpecies>& ions,bool electronMass,double gamma,double b,
    double cosine) {
  if(!(std::isfinite(ne)&&ne>0.0&&std::isfinite(te)&&te>0.0&&gamma>1.0&&
       std::isfinite(b)&&std::isfinite(cosine)&&std::abs(cosine)<=1.0)||ions.empty())
    return Core::Result<PlasmaState>::Failure(Core::StatusCode::InvalidConfiguration,
        "EOS requires positive density/temperature, gamma>1, and explicit ions");
  double charge=0.0,mass=0.0,ionPressureFactor=0.0;
  for(const auto& ion:ions) {
    if(!(ion.abundancePerProton>=0.0&&ion.chargeNumber>0.0&&ion.massKg>0.0&&ion.temperatureK>0.0))
      return Core::Result<PlasmaState>::Failure(Core::StatusCode::InvalidConfiguration,
          "ion abundance, charge, mass, and temperature are explicit positive values");
    charge+=ion.chargeNumber*ion.abundancePerProton;
    mass+=ion.massKg*ion.abundancePerProton;
    ionPressureFactor+=ion.abundancePerProton*ion.temperatureK;
  }
  if(!(charge>0.0)) return Core::Result<PlasmaState>::Failure(
      Core::StatusCode::InvalidConfiguration,"composition has no positive charge density");
  PlasmaState result; result.electronNumberDensityM3=ne;
  result.protonNumberDensityM3=ne/charge;
  result.massDensityKgM3=result.protonNumberDensityM3*mass+
      (electronMass?Constants::kElectronMassKg*ne:0.0);
  result.pressurePa=Constants::kBoltzmannJPerK*(
      result.protonNumberDensityM3*ionPressureFactor+ne*te);
  result.soundSpeedMPerS=std::sqrt(gamma*result.pressurePa/result.massDensityKgM3);
  result.alfvenSpeedMPerS=std::abs(b)/std::sqrt(
      Constants::kVacuumPermeabilityHPerM*result.massDensityKgM3);
  const double va2=result.alfvenSpeedMPerS*result.alfvenSpeedMPerS;
  const double cs2=result.soundSpeedMPerS*result.soundSpeedMPerS;
  result.fastSpeedMPerS=std::sqrt(0.5*(va2+cs2+
      std::sqrt(std::max(0.0,(va2+cs2)*(va2+cs2)-4*va2*cs2*cosine*cosine))));
  return Core::Result<PlasmaState>::Success(result);
}
} }
