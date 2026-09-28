#include "sep_coronal_cme/model_configuration.h"

namespace SEP {
namespace CoronalCME {

const char* Name(RunIntent value) noexcept {
  switch (value) {
    case RunIntent::ProductionShockInjection: return "production-shock-injection";
    case RunIntent::AnalyticVerification: return "analytic-verification";
  }
  return "invalid";
}

const char* Name(TransportModel value) noexcept {
  switch (value) {
    case TransportModel::BallisticVerification: return "ballistic-verification";
    case TransportModel::Parker: return "parker";
    case TransportModel::FocusedPitchAngleDiffusion:
      return "focused-pitch-angle-diffusion";
    case TransportModel::FocusedDiscreteScattering:
      return "focused-discrete-scattering";
  }
  return "invalid";
}

const char* Name(TransportFrame value) noexcept {
  switch (value) {
    case TransportFrame::Inertial: return "inertial";
    case TransportFrame::RigidCorotating: return "rigid-corotating";
  }
  return "invalid";
}

const char* Name(SolarRotationModel value) noexcept {
  switch (value) {
    case SolarRotationModel::Rigid: return "rigid";
    case SolarRotationModel::LatitudeDependentVerification:
      return "latitude-dependent-verification";
  }
  return "invalid";
}

const char* Name(WindModel value) noexcept {
  switch (value) {
    case WindModel::FluxTubePolytropic: return "flux-tube-polytropic";
    case WindModel::EmpiricalKinematic: return "empirical-kinematic";
  }
  return "invalid";
}

const char* Name(WindEnergyClosure value) noexcept {
  switch (value) {
    case WindEnergyClosure::Isothermal: return "isothermal";
    case WindEnergyClosure::Polytropic: return "polytropic";
    case WindEnergyClosure::EmpiricalProfile: return "empirical-profile";
  }
  return "invalid";
}

const char* Name(ClosedFieldModel value) noexcept {
  switch (value) {
    case ClosedFieldModel::IsothermalHydrostatic:
      return "isothermal-hydrostatic";
    case ClosedFieldModel::PolytropicHydrostatic:
      return "polytropic-hydrostatic";
  }
  return "invalid";
}

const char* Name(InterfaceRepresentation value) noexcept {
  switch (value) {
    case InterfaceRepresentation::SharpOneSided: return "sharp-one-sided";
    case InterfaceRepresentation::FiniteWidthVolume:
      return "finite-width-volume";
  }
  return "invalid";
}

const char* Name(InterfacePolicy value) noexcept {
  switch (value) {
    case InterfacePolicy::DiagnosticKinematic: return "diagnostic-kinematic";
    case InterfacePolicy::BoundedApproximation: return "bounded-approximation";
    case InterfacePolicy::StationaryTangentialDiscontinuity:
      return "stationary-td-equilibrium";
  }
  return "invalid";
}

}  // namespace CoronalCME
}  // namespace SEP
