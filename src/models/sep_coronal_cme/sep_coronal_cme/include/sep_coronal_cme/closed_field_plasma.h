#ifndef SEP_CORONAL_CME_CLOSED_FIELD_PLASMA_H
#define SEP_CORONAL_CME_CLOSED_FIELD_PLASMA_H

#include "sep_coronal_cme/vector_math.h"
#include "sep_status.h"

namespace SEP { namespace CoronalCME {

struct HydrostaticState { double densityKgM3=0.0,pressurePa=0.0,temperatureScale=0.0; };
double RotatingEffectivePotential(Vec3 positionM,Vec3 omegaRadPerS);
Core::Result<HydrostaticState> IsothermalHydrostatic(double baseDensityKgM3,
    double basePressurePa,double potentialDifferenceM2S2);
Core::Result<HydrostaticState> PolytropicHydrostatic(double baseDensityKgM3,
    double basePressurePa,double gammaClosed,double potentialDifferenceM2S2);
Core::Status CheckFootpointCompatibility(double pressureFromFootpointA,
    double pressureFromFootpointB,double relativeTolerance);

} }  // namespace SEP::CoronalCME
#endif
