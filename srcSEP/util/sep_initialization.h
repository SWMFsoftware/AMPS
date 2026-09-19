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

enum class RefinementProfile { Linear, PowerLaw, Smoothstep };
enum class TubeRadiusMode { PhysicalConstant, ConstantAngularWidth };

// SI-only schema consumed by srcSEP's standalone runtime initialization.
// The spatial mesh supports the embedded field line; particles themselves
// remain on the one-dimensional PIC field-line representation.
struct Configuration {
  unsigned schemaVersion = 1;

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
};

// Parse a complete version-1 INI document. Unknown/duplicate sections and
// keys are errors; no value is silently inherited from a production default.
Transport::Status ParseText(const std::string& text, Configuration* result);
Transport::Status LoadFile(const std::string& path, Configuration* result);
Transport::Status Validate(const Configuration& configuration);

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
