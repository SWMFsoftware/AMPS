// ============================================================================
// srcSEP3D section of the shared post-compile input file
//
// This layer deliberately knows nothing about AMPS or MPI.  It expands the
// shared file syntax (#include, ! comments, and continued logical lines), then
// interprets only ``#section begin: sep3d``.  Keeping text handling here lets
// the native boundary run the parser after Init_BeforeParser while unit tests
// exercise exactly the same implementation without allocating a mesh.
// ============================================================================

#ifndef SEP3D_RUNTIME_APPLICATION_INPUT_H
#define SEP3D_RUNTIME_APPLICATION_INPUT_H

#include "../core/sep3d_types.h"

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace SEP3D {
namespace RuntimeModel {

struct Sep3dApplicationInput {
  // A shared-section run must have a finite, explicit iteration horizon.
  // RunConfiguration3D's large library default is useful to coupled hosts,
  // but silently inheriting it in a standalone particle-producing job would
  // make a typo look like an effectively unbounded source calculation.
  std::uint64_t maximumTimeSteps = 0;

  // This maps to SourceOptions::samplesPerStep.  The existing srcSEP3D source
  // contract interprets the value per compiled species; the parser must not
  // silently change that established particle-accounting rule.
  std::uint64_t particlesPerIteration = 0;

  // These selectors are mandatory in shared-section mode.  Keeping the
  // authority names in the parsed record prevents an inline reduced-front
  // block from becoming active merely because it happened to be present.
  std::string shockModel;
  std::string backgroundPlasmaModel;
  std::string sourceModel;

  // One global AMPS step is derived only after the distributed mesh exists:
  // dt=f*h_min/v_max.  The margin is dimensionless and lies in (0,1]; speed
  // and the explicitly chosen normalization radius use SI units.
  double maximumParticleSpeedMPerS = 0.0;
  double timeStepMarginFactor = 0.0;
  double sourceNormalizationRadiusM = 0.0;

  // The particle-source subsection is explicit even though only the first
  // statistical representation is currently executable.  Energies are total
  // kinetic energy per particle [J].  ``compression-ratio`` derives the
  // isotropic DSA phase-space exponent q=3X/(X-1) independently on every
  // accepted triangle; ``constant`` uses fixedPhaseSpacePowerIndex.
  std::string particleWeightingModel;
  std::string momentumPowerLawModel;
  double minimumInjectionEnergyJ = 0.0;
  double maximumInjectionEnergyJ = 0.0;
  double fixedPhaseSpacePowerIndex = 0.0;
  std::uint64_t maximumInjectionEventsPerSpeciesPerStep = 0;

  // The reduced model's complete canonical v1.1 assignment layer is carried
  // in memory. Relative magnetic/PFSS assets resolve from the physical file
  // containing the subsection begin directive, so moving the process working
  // directory cannot change which bytes initialize the field.
  std::string reducedShockConfiguration;
  std::string reducedShockAssetDirectory;

  // Resolved provenance is diagnostic rather than an alternative physics
  // identity.  The numeric value enters RunConfiguration3D's fingerprint;
  // these fields make the summary and parser failures traceable to text.
  std::string rootFile;
  std::string valueFile;
  std::size_t valueLine = 0;
  std::vector<std::string> expandedFiles;
};

// Parse one shared input tree transactionally.  ``result`` is unchanged on
// failure.  Include paths are relative to the file containing the directive;
// recursive cycles and a depth greater than 64 are rejected with provenance.
Core::Status ParseSep3dApplicationInput(
    const std::string& path, Sep3dApplicationInput* result);

// Human-readable rank-zero receipt emitted only after parsing and immutable
// configuration replacement both succeed.
std::string Sep3dApplicationInputSummary(
    const Sep3dApplicationInput& input);

}  // namespace RuntimeModel
}  // namespace SEP3D

#endif  // SEP3D_RUNTIME_APPLICATION_INPUT_H
