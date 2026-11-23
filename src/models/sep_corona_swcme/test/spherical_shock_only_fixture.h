#ifndef SEP_CORONA_SWCME_SPHERICAL_SHOCK_ONLY_VERIFICATION_FIXTURE_H
#define SEP_CORONA_SWCME_SPHERICAL_SHOCK_ONLY_VERIFICATION_FIXTURE_H

#include "sep_corona_swcme/spherical_shock_only.h"
#include "sep_coronal_cme/constants.h"

namespace SEP {
namespace CoronaSwcme {
namespace SphericalShockOnly {
namespace Verification {
// TEST INPUTS ONLY: inherited Stage14 analytic dipole/isothermal controls,
// explicitly identified as manufactured. These helpers are shared so a thin
// bridge tests the SAME physical realization, never to supply production
// defaults or manufacture magnetogram/solar-wind observational validation.
namespace C = CoronalCME;
constexpr double Rs = C::Constants::kSolarRadiusM, Au = 149597870700.;
inline Definition Configuration(int polarity = 1) {
  Definition d;
  // Exact inherited Stage14 analytic-verification anchors, not measurements
  // or event predictions. The fixture's explicit10-AU provider guard permits
  // independent8/10-AU sensitivity checks without changing physical5-AU loss.
  d.ambient.composition = {5. / 3., 1.5e6, 1.5e6, 6e6, .04, false};
  auto& a = d.ambient.ambient;
  a.sourceSurfaceRadiusM = 2.5 * Rs;
  a.referenceRadiusM = Au;
  a.electronDensityAtReferenceM3 = 5e6;
  a.closedBaseElectronDensityM3 = 1e14;
  a.rotationRateRadPerS = 2.86533e-6;
  a.minimumMagneticFieldT = 1e-14;
  a.traceStepM = 1e6;
  a.gradientRelativeStep = 1e-5;
  a.traceMaximumSteps = 200000;
  a.windTablePoints = 1025;
  a.sphericalOpenMassFlux = false;
  a.harmonics = polarity == 1 ? "1,0,0.0001,0" : "1,0,-0.0001,0";
  d.ambient.support = {0, 1200, Rs, Rs, 10 * Au, 1e-6, 2000};
  d.ambient.coordinateFrame = "HCI-like-inertial";
  d.ambient.physicsFingerprint = "S13-explicit-manufactured-dipole-isothermal-v1";
  d.physicalInnerRadiusM = 1.05 * Rs;
  d.physicalOuterRadiusM = 5 * Au;
  d.inputIdentity = "S13-independent-core-verification-v1";
  d.front.mode = swcme::kinematics::Mode::Ballistic;
  d.front.r0_m = 4 * Rs;
  d.front.V0_m_s = 0;
  d.front.Vsw_m_s = 400000;
  d.front.Gamma_m_inv = 0;
  return d;
}
inline Numerics IntegrationControls() {
  return {0.1, 1e-9, 20, 1e8, 65536};
}
inline OpenTubeIsoPolyFamilyInput ThermalFamily() {
  OpenTubeIsoPolyFamilyInput family;
  family.lowerSourceColatitudeRad=.24;
  family.upperSourceColatitudeRad=.68;
  family.inputIdentity="S13-manufactured-thermal-band";
  auto& p=family.population;
  p.massLoadingKgPerSWb=2e-6;
  p.referenceDensityKgM3=1e-20;
  p.electron={1.15,10*Rs,18*Rs};
  p.proton={1.25,12*Rs,22*Rs};
  p.alpha={1.20,12*Rs,25*Rs};
  p.angularTraceToleranceRad=1e-12;
  return family;
}
}  // namespace Verification
}  // namespace SphericalShockOnly
}  // namespace CoronaSwcme
}  // namespace SEP
#endif
