#ifndef SRCSEP_UTIL_SEP_INITIALIZATION_H
#define SRCSEP_UTIL_SEP_INITIALIZATION_H

#include "sep_common_header_path.h"
#include SRCSEP_SEP_COMMON_HEADER(sep_transport_common.h)

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace SEP {
namespace Initialization {

// Small AMPS-independent vector used while parsing and validating startup
// input.  Keeping PIC and MPI out of this layer lets the complete input and
// mesh-resolution laws run before AMPS allocates global state.
struct Vec3 {
  double x = 0.0;
  double y = 0.0;
  double z = 0.0;
  Vec3() {}
  Vec3(double xValue, double yValue, double zValue)
      : x(xValue), y(yValue), z(zValue) {}
};

// Raw provider assignment transported to the canonical SWCME1D resolver.
// This startup layer deliberately does not interpret model keys or units: the
// model-owned resolver remains the single scientific authority.  Keeping line
// provenance here lets a malformed physical value identify its input record.
struct SwcmeAssignment {
  std::string key;
  std::string value;
  std::size_t line = 0;
};

enum class RefinementProfile { Linear, PowerLaw, Smoothstep };
enum class TubeRadiusMode { PhysicalConstant, ConstantAngularWidth };

// SI-only schema consumed by srcSEP's standalone runtime initialization.
// The spatial mesh supports the embedded field line; particles themselves
// remain on the one-dimensional PIC field-line representation.
struct Configuration {
  unsigned schemaVersion = 1;

  // Explicit numerical initialization.  Version 2 never derives either value
  // from mesh resolution or an unrelated boundary source: those choices are
  // campaign physics and must be supplied by the operator.
  double timeStepS = 0.0;
  std::uint64_t macroparticlesPerStep = 0;
  // Common base statistical weight installed for every species in AMPS'
  // compiled SpeciesList.  Species identity/mass/charge never come from this
  // post-compile file; injection applies any species-specific correction via
  // the already validated source table.
  double particleWeight = 0.0;

  // srcSEP samples a one-dimensional field line at a heliocentric radius.
  // The radius is SI and is installed into the retained field-line sampler.
  double observerHeliocentricRadiusM = 0.0;

  // Initialization products are written only after the corresponding AMPS
  // mesh/field line exists.  Explicit paths avoid an undocumented working-
  // directory convention and are part of the startup fingerprint.
  std::string meshTecplotFile;
  std::string fieldLineTecplotFile;

  Vec3 parkerOriginM;
  Vec3 parkerInitialPointM;
  double parkerLengthM = 0.0;
  std::uint64_t parkerPointCount = 0;
  double solarWindSpeedMPerS = 0.0;
  double solarRotationRateRadPerS = 0.0;

  double innerRadiusM = 0.0;
  double outerRadiusM = 0.0;
  double globalCellSizeM = 0.0;
  double minimumCellSizeM = 0.0;
  unsigned maximumMeshLevel = 0;

  bool solarRefinementEnabled = true;
  double solarSurfaceCellSizeM = 0.0;
  double solarTransitionOuterRadiusM = 0.0;
  RefinementProfile solarProfile = RefinementProfile::Smoothstep;
  double solarExponent = 1.0;

  bool tubeRefinementEnabled = true;
  double tubeReferenceRadiusM = 0.0;
  double tubeRadiusAtReferenceM = 0.0;
  TubeRadiusMode tubeRadiusMode = TubeRadiusMode::ConstantAngularWidth;
  double tubeCenterCellSizeM = 0.0;
  RefinementProfile tubeProfile = RefinementProfile::Smoothstep;
  double tubeExponent = 1.0;

  // A schema-v2 file supplies every active canonical SWCME setting.  The
  // application forwards this layer unchanged to SW1DAdapter::Configure().
  std::vector<SwcmeAssignment> swcmeAssignments;
};

// Parse a complete version-1 or version-2 INI document.  Version 2 adds the
// explicit numerical/source/observer/output/SWCME initialization contract.
// Unknown or duplicate sections and keys are errors; no schema-required value
// is silently inherited from a production default.
Transport::Status ParseText(const std::string& text, Configuration* result);
Transport::Status LoadFile(const std::string& path, Configuration* result);
Transport::Status Validate(const Configuration& configuration);

// Retarget the two initialization products to one command-line-selected
// directory while preserving the leaf names declared in [output].  This is a
// pure configuration transform; the MPI-aware runtime creates parent
// directories immediately before collective output begins.
Transport::Status ApplyOutputDirectoryOverride(
    const std::string& directory, Configuration* configuration);

// Install exactly once before amps_init_mesh().  Absence means legacy mode and
// intentionally preserves every historical hard-coded srcSEP mesh choice.
Transport::Status Install(const Configuration& configuration);
bool HasActive();
const Configuration& Active();

// Stable configuration identity printed before AMPS initialization.
std::string Fingerprint(const Configuration& configuration);

// Generate an outward Parker curve at uniform arc-length spacing. pointCount
// includes both endpoints. The midpoint tangent method is second-order in the
// geometry while preserving each requested segment length to roundoff.
Transport::Status BuildParkerLine(const Configuration& configuration,
                                  std::vector<Vec3>* points);

// Write the finite one-dimensional transport mesh as a Tecplot ordered POINT
// zone.  The writer validates/builds the complete curve before opening the
// destination, so invalid physics never leaves a partial visualization file.
Transport::Status WriteParkerLineTecplot(
    const Configuration& configuration, const std::string& path);

double TubeRadiusM(double heliocentricRadiusM,
                   const Configuration& configuration);
double TubeDistanceM(const Vec3& positionM,
                     const Configuration& configuration);
double RequestedCellSizeM(const Vec3& positionM,
                          const Configuration& configuration);

const char* Name(RefinementProfile profile);
const char* Name(TubeRadiusMode mode);

}  // namespace Initialization
}  // namespace SEP

#endif  // SRCSEP_UTIL_SEP_INITIALIZATION_H
