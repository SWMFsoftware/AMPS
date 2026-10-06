#ifndef SEP_CORONA_SWCME_SHEATH_DIAGNOSTICS_H
#define SEP_CORONA_SWCME_SHEATH_DIAGNOSTICS_H

#include "sep_corona_swcme/sheath_model.h"

#include <vector>

namespace SEP { namespace CoronaSwcme {

// Numerical controls for differentiating a smooth part of the material map.
// The three coordinates are (theta,phi,tau); the time derivative is taken at
// fixed labels and is therefore a material derivative.  A caller must keep
// every stencil on one smooth, already-admitted branch.  Shocks and contacts
// are checked with jump conditions, never by differentiating through them.
struct SheathDiagnosticOptions {
  double angularStepRad = 0.0;
  double admissionTimeStepS = 0.0;
  double evolutionTimeStepS = 0.0;
  // SI heliocentric gravitational parameter.  Set to zero only in a declared
  // force-free manufactured fixture; the production event uses solar gravity.
  double gravitationalParameterM3S2 = 1.32712440018e20;
};

struct SheathResidualSample {
  SheathMaterialLabel label;
  double epochS = 0.0;
  CoronalCME::Vec3 positionM;
  CoronalCME::Vec3 inertiaPaPerM;
  CoronalCME::Vec3 pressureGradientPaPerM;
  CoronalCME::Vec3 lorentzPaPerM;
  CoronalCME::Vec3 gravityPaPerM;
  CoronalCME::Vec3 momentumResidualPaPerM;
  double momentumScalePaPerM = 0.0;
  double localMomentumRatio = 0.0;
  // Conservative ideal-MHD total-energy residual in W m^-3:
  //   D E/Dt + E div(U) + div[(p+B^2/2mu0)U-(U.B)B/mu0] - rho g.U.
  // E excludes gravitational potential, whose work is the last term.
  double energyResidualWPerM3 = 0.0;
  double energyScaleWPerM3 = 0.0;
  double inertiaWorkWPerM3 = 0.0;
  double pressureGradientWorkWPerM3 = 0.0;
  double lorentzWorkWPerM3 = 0.0;
  double gravityWorkWPerM3 = 0.0;
  double residualForceWorkWPerM3 = 0.0;
};

struct WeightedSheathDiagnosticPoint {
  SheathMaterialLabel label;
  // Physical current volume represented by this point.  This must come from
  // the curved map metric, not area times an unweighted radial thickness.
  double volumeM3 = 0.0;
};

struct SheathResidualReport {
  std::vector<SheathResidualSample> samples;
  std::size_t unsupportedSamples = 0;
  double sampledVolumeM3 = 0.0;
  double integratedAbsoluteMomentumResidualN = 0.0;
  double integratedMomentumScaleN = 0.0;
  double integratedMomentumRatio = 0.0;
  double volumeWeightedLocalMomentumP99 = 0.0;
  double signedEnergyResidualW = 0.0;
  double absoluteEnergyResidualW = 0.0;
  double integratedEnergyScaleW = 0.0;
  double signedResidualForceWorkW = 0.0;
  double absoluteResidualForceWorkW = 0.0;
  double signedInertiaWorkW = 0.0;
  double signedPressureGradientWorkW = 0.0;
  double signedLorentzWorkW = 0.0;
  double signedGravityWorkW = 0.0;
  double forceWorkRatio = 0.0;
};

// These diagnostics do not modify or qualify the prescribed diagnostic map.
// They are intentionally tied to a smooth material deformation and are not the
// conservation algorithm for the selected per-ray finite-volume piston model.
// That model requires face-flux/source/piston-work ledgers, shock-zone
// convergence, well-balance tests and ray-assembly/divergence diagnostics from
// docs/BG3D4_PISTON_CLOSURE.md.  Reusing a PASS here as piston evidence would
// mix two different governing systems.
//
// For the current map, these routines form
// every differential term from independent finite differences of the public
// Evaluate interface, so exact Cauchy identities inside the implementation
// cannot make the momentum or energy result pass by construction.
Core::Result<SheathResidualSample> EvaluateSheathResidual(
    const ShockFedSheathModel&,const SheathMaterialLabel&,double epochS,
    double gammaAdiabatic,const SheathDiagnosticOptions&);

Core::Result<SheathResidualReport> IntegrateSheathResiduals(
    const ShockFedSheathModel&,const std::vector<WeightedSheathDiagnosticPoint>&,
    double epochS,double gammaAdiabatic,const SheathDiagnosticOptions&);

} } // namespace SEP::CoronaSwcme

#endif
