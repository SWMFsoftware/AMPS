// ============================================================================
// srcSEP3D/main_lib.cpp
//
// Phase R2 AMPS application boundary.
//
// The production boundary now owns the typed Runtime introduced in R2.  Both
// standalone and coupled hosts install a validated immutable configuration and
// use the same Runtime transitions; this file does not parse process arguments
// or parameter files.  Mesh construction and transport physics belong to
// later phases, so AMPS callbacks that would cross those boundaries still stop
// explicitly rather than executing the retired prototype defaults.
// ============================================================================

#include "SEP3D.h"

#include <cstdlib>
#include <iostream>

namespace {

[[noreturn]] void StopAtUnimplementedPhase(const char* entryPoint,
                                            const char* requiredPhase) {
  std::cerr
      << "[srcSEP3D:R2] " << entryPoint << " reached the implemented Runtime "
      << "boundary, but " << requiredPhase << " is not implemented. The "
      << "retired wedge, Maxwellian source, and legacy sampler remain removed; "
      << "no placeholder physics was executed.\n";
  std::abort();
}

} // namespace

// ---------------------------------------------------------------------------
// Required Exosphere hooks
//
// AMPS links these symbols for applications based on the Exosphere module.
// They are inert compatibility hooks, not srcSEP3D physics.  Named casts make
// the intentionally unused inputs explicit and keep strict-warning builds
// clean.
// ---------------------------------------------------------------------------
double Exosphere::OrbitalMotion::GetTAA(SpiceDouble et) {
  (void)et;
  return 0.0;
}

int Exosphere::ColumnIntegral::GetVariableList(char* variableList) {
  (void)variableList;
  return 0;
}

void Exosphere::ColumnIntegral::ProcessColumnIntegrationVector(
    double* result, int resultLength) {
  (void)result;
  (void)resultLength;
}

double Exosphere::GetSurfaceTemperature(double cosSubsolarAngle,
                                         double* position) {
  (void)cosSubsolarAngle;
  (void)position;
  return 0.0;
}

char Exosphere::SO_FRAME[_MAX_STRING_LENGTH_PIC_] = "HCI_like_inertial";
char Exosphere::ObjectName[_MAX_STRING_LENGTH_PIC_] = "Sun";

void Exosphere::ColumnIntegral::CoulumnDensityIntegrant(
    double* result, int resultLength, double* position,
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node) {
  (void)result;
  (void)resultLength;
  (void)position;
  (void)node;
}

double Exosphere::SurfaceInteraction::StickingProbability(
    int spec, double& reemissionParticleFraction, double temperature) {
  (void)spec;
  (void)temperature;
  reemissionParticleFraction = 0.0;
  return 0.0;
}

SEP3D::RuntimeModel::Runtime& SEP3D::ApplicationRuntime() {
  // AMPS exposes process-level application callbacks, so one process-owned
  // Runtime is the matching ownership scope.  The object's state and counters
  // are encapsulated and change only through typed, transactional methods.
  static RuntimeModel::Runtime runtime;
  return runtime;
}

SEP3D::Core::Status SEP3D::ConfigureApplication(
    const std::shared_ptr<const RuntimeModel::RunConfiguration3D>& configuration) {
  return ApplicationRuntime().Configure(configuration);
}

void SEP3D::Init_BeforeParser() {
  // Deliberately empty.  Coupled entry points must not inspect argc/argv or
  // AMPS_PARAM.in; the host supplies a resolved RunConfiguration3D through
  // ConfigureApplication before invoking mesh setup.
}

double localResolution(double* position) {
  (void)position;
  StopAtUnimplementedPhase("localResolution", "Phase M mesh resolution");
}

double InitLoadMeasure(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node) {
  (void)node;
  // The callback is retained only to satisfy the AMPS application interface.
  // Returning a uniform weight is safe for a link check; amps_init_mesh stops
  // before the callback can be installed in an executable run.
  return 1.0;
}

bool TrajectoryTrackingCondition(double* position, double* velocity, int spec,
                                 void* particleData) {
  (void)position;
  (void)velocity;
  (void)spec;
  (void)particleData;
  return false;
}

void amps_init_mesh() {
  if (SEP3D::ApplicationRuntime().state() ==
      SEP3D::RuntimeModel::LifecycleState::Created) {
    std::cerr
        << "[srcSEP3D:R2] amps_init_mesh requires the host to install an "
        << "immutable RunConfiguration3D before AMPS allocates mesh storage.\n";
    std::abort();
  }
  StopAtUnimplementedPhase("amps_init_mesh", "Phase M mesh construction");
}

void amps_init() {
  StopAtUnimplementedPhase("amps_init", "background and transport phases");
}

int amps_time_step() {
  StopAtUnimplementedPhase("amps_time_step", "transport mover phases");
}
