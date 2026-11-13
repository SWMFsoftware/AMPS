// pic.h must precede Moon/Exosphere headers: it provides the generated AMPS
// feature and species macros (_NA_SPEC_, AVX mode, source switches) that those
// headers intentionally require to be defined before they are parsed.
#include "pic.h"

#include "SodiumCase1.h"

#include "Moon.h"
#include "MoonInput.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace Moon {
namespace Species {
namespace SodiumCase1 {
namespace {

const Runtime::Configuration* gConfiguration=nullptr;

}  // namespace

bool ConfigureProduction(const Runtime::Configuration& configuration,
                         cInternalSphericalData* logicalSurface,
                         std::string* error) {
  if (error==nullptr) return false;
  error->clear();
  if (!configuration.sodiumExosphereEnabled) {
    gConfiguration=nullptr;
    return true;
  }
  if (logicalSurface==nullptr) {
    *error="Tenishev-2013 Case 1 requires an allocated logical lunar surface";
    return false;
  }
#if _EXOSPHERE_SOURCE__IMPACT_VAPORIZATION_ != _EXOSPHERE_SOURCE__ON_
  *error="Tenishev-2013 Case 1 requires Source:ImpactVaporization=on in "
         "input/moon.input and a clean reconfiguration/rebuild";
  return false;
#endif
#if _EXOSPHERE_SOURCE__PHOTON_STIMULATED_DESPRPTION_ != _EXOSPHERE_SOURCE__ON_
  *error="Tenishev-2013 Case 1 requires Source:PhotonStimulatedDesorption=on "
         "in input/moon.input and a clean reconfiguration/rebuild";
  return false;
#endif
#if _EXOSPHERE_SOURCE__THERMAL_DESORPTION_ != _EXOSPHERE_SOURCE__ON_
  *error="Tenishev-2013 Case 1 requires Source:ThermalDesorption=on in "
         "input/moon.input and a clean reconfiguration/rebuild";
  return false;
#endif
#if _EXOSPHERE_SOURCE__SOLAR_WIND_SPUTTERING_ == _EXOSPHERE_SOURCE__ON_
  *error="Tenishev-2013 Case 1 excludes direct solar-wind sputtering; set "
         "Source:SolarWindSputtering=off and cleanly rebuild";
  return false;
#endif
#if _PIC_PHOTOLYTIC_REACTIONS_MODE_ != _PIC_PHOTOLYTIC_REACTIONS_MODE_ON_
  *error="Tenishev-2013 Case 1 requires PhotolyticReactions=on in "
         "input/moon.input and a clean reconfiguration/rebuild";
  return false;
#endif
#if _NA_SPEC_ < 0
  *error="Tenishev-2013 Case 1 requires neutral Na in SpeciesList";
  return false;
#endif

  // Source switches alone are insufficient: PSD and TD still obtain their
  // coefficients from the generated AMPS source arrays.  Check the exact
  // Case-1 build contract before emitting a particle, so a later edit to
  // input/moon.input cannot leave a misleading paper label in the run receipt.
  const auto same=[](double value,double expected) {
    return std::abs(value-expected)<=
        32.0*std::numeric_limits<double>::epsilon()*
        std::max(std::abs(value),std::abs(expected));
  };
  if (!same(Exosphere::SourceProcesses::PhotonStimulatedDesorption::
                PhotonStimulatedDesorption_PhotonFlux_1AU,2.0e18) ||
      !same(Exosphere::SourceProcesses::PhotonStimulatedDesorption::
                PhotonStimulatedDesorption_CrossSection[_NA_SPEC_],3.0e-25) ||
      !same(Exosphere::SourceProcesses::ThermalDesorption::
                ThermalDesorption_uThermal[_NA_SPEC_],1.85*eV2J) ||
      !same(Exosphere::SourceProcesses::ThermalDesorption::
                ThermalDesorption_VibrationalFrequency[_NA_SPEC_],1.0e13) ||
      !same(Exosphere::SourceProcesses::ImpactVaporization::
                ImpactVaporization_HeliocentricDistance,_AU_) ||
      Exosphere::SourceProcesses::ImpactVaporization::
          ImpactVaporization_SourceRatePowerIndex!=0.0) {
    *error="the generated PSD/TD/MIV coefficients do not match Tenishev-2013 "
           "Case 1; restore the documented input/moon.input values, remove "
           "build, rerun Config.pl, and rebuild";
    return false;
  }

  const int count=PIC::BC::InternalBoundary::Sphere::TotalSurfaceElementNumber;
  if (count<=0 || logicalSurface->SurfaceElementPopulation==nullptr ||
      logicalSurface->SurfaceElementArea==nullptr) {
    *error="Tenishev-2013 Case 1 found an unallocated surface inventory";
    return false;
  }
  for (int species=0;species<PIC::nTotalSpecies;++species)
    for (int element=0;element<count;++element) {
      const double density=species==_NA_SPEC_ ?
          configuration.sodiumSurfaceDensityM2 : 0.0;
      const double area=logicalSurface->SurfaceElementArea[element];
      if (!std::isfinite(area) || area<=0.0) {
        *error="Tenishev-2013 Case 1 found a non-positive logical-surface area";
        return false;
      }
      // The stored inventory has units atoms, whereas the paper and input
      // receipt prescribe a uniform areal abundance in atoms m^-2.
      logicalSurface->SurfaceElementPopulation[species][element]=density*area;
    }
  gConfiguration=&configuration;
  return true;
}

double ImpactVaporizationRate(int species) {
  if (!Runtime::HasConfiguration() ||
      !Runtime::GetConfiguration().sodiumExosphereEnabled ||
      species!=_NA_SPEC_) return 0.0;
  return Runtime::GetConfiguration().sodiumImpactVaporizationRateS;
}

double ImpactVaporizationTemperatureK(int species) {
  if (!Runtime::HasConfiguration() ||
      !Runtime::GetConfiguration().sodiumExosphereEnabled ||
      species!=_NA_SPEC_) return 0.0;
  return Runtime::GetConfiguration().sodiumImpactVaporizationTemperatureK;
}

bool PrescribedSurfaceDensityM2(int species,int surfaceElement,
                                double* numberDensityM2) {
  if (gConfiguration==nullptr || numberDensityM2==nullptr ||
      surfaceElement<0) return false;
  // Return an active prescription for every compiled species.  Case 1 is a
  // neutral-Na model, so any additional compiled species has zero reservoir;
  // allowing it to fall through to a generated legacy abundance would invent
  // an undeclared source.
  *numberDensityM2=species==_NA_SPEC_ ?
      gConfiguration->sodiumSurfaceDensityM2 : 0.0;
  return true;
}

double PhotoionizationLifetimeS() {
  if (Runtime::HasConfiguration() &&
      Runtime::GetConfiguration().sodiumExosphereEnabled)
    return Runtime::GetConfiguration().sodiumPhotoionizationLifetimeS;
  return 3600.0*5.8/std::pow(0.4,2);
}

bool ProductionEnabled() { return gConfiguration!=nullptr; }

}  // namespace SodiumCase1
}  // namespace Species
}  // namespace Moon
