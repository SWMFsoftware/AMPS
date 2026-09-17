// ============================================================================
// srcSEP3D/main_lib.cpp
//
// Phase R0 AMPS application boundary.
//
// R0 has one narrow purpose: establish a compile-clean, truthful production
// tree before a new three-dimensional Runtime or transport mover is added.
// The former file silently retained a narrow axisymmetric wedge, uniform
// resolution, Maxwellian prepopulation, and placeholder sampling.  Running
// those defaults could look like a valid srcSEP3D calculation even though they
// belonged to the retired prototype.
//
// Therefore the AMPS-required symbols are retained for production linking,
// but every entry point that would start a simulation terminates with an
// explicit R0 diagnostic.  Phase R2 will replace this guard with the typed
// Runtime lifecycle; Phase M will supply the real mesh and refinement law.
// No placeholder physics is executed in the interim.
// ============================================================================

#include "SEP3D.h"

#include <cstdlib>
#include <iostream>

namespace {

[[noreturn]] void StopBeforePrototypePhysics(const char* entryPoint) {
  std::cerr
      << "[srcSEP3D:R0] " << entryPoint << " was called, but Phase R0 is a "
      << "compile/link baseline only. The retired wedge, Maxwellian source, "
      << "and legacy sampler have been removed. Complete the Runtime and mesh "
      << "phases before executing a physical simulation.\n";
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

void SEP3D::Init_BeforeParser() {
  // Deliberately empty in R0.  Configuration ownership is introduced with
  // the Runtime lifecycle rather than through file-scope mutable globals.
}

double localResolution(double* position) {
  (void)position;
  StopBeforePrototypePhysics("localResolution");
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
  StopBeforePrototypePhysics("amps_init_mesh");
}

void amps_init() {
  StopBeforePrototypePhysics("amps_init");
}

int amps_time_step() {
  StopBeforePrototypePhysics("amps_time_step");
}
