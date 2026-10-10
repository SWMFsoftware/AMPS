#ifndef SRCMOON_LUNAR_SURFACE_H
#define SRCMOON_LUNAR_SURFACE_H

#include "MoonInput.h"

#include <array>
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace Moon {
namespace Surface {

// Subset of the detached PDS label that controls decoding and coordinates.
// Text fields are retained verbatim (after removing PDS quoting/units) so the
// run log can report the archive-defined frame instead of inventing one.
struct LolaMetadata {
  std::string productId;
  std::string dataSetId;
  std::string projection;
  std::string coordinateSystemType;
  std::string coordinateSystemName;
  std::string positiveLongitudeDirection;
  std::string sampleType;
  std::string unit;
  std::size_t lines = 0;
  std::size_t lineSamples = 0;
  unsigned int sampleBits = 0;
  double scalingFactor = 0.0;
  double referenceRadiusM = 0.0;
  double mapResolutionPixelsPerDegree = 0.0;
  double minimumLatitudeDeg = 0.0;
  double maximumLatitudeDeg = 0.0;
  double westernmostLongitudeDeg = 0.0;
  double easternmostLongitudeDeg = 0.0;
};

// Read-only, in-memory view of a supported LOLA raster.  The stored vector is
// the native signed integer DN array; ScalingFactor is applied during lookup.
// OFFSET is interpreted as the reference radius, as declared by LDEM_4.LBL,
// and is added only when constructing a physical vertex radius.
class LolaDem {
 public:
  // Load and validate the configured detached label/image pair.  On failure,
  // the object is not replaced with a partially decoded product.
  bool Load(const Runtime::Configuration& configuration, std::string* error);
  // Bilinear elevation in metres at east-positive longitude and
  // planetocentric latitude in degrees.  Longitude is periodic; latitude is
  // clamped to the first/last pixel centres at the poles.
  double ElevationM(double eastLongitudeDeg,
                    double planetocentricLatitudeDeg) const;
  const LolaMetadata& metadata() const { return metadata_; }

 private:
  LolaMetadata metadata_;
  std::vector<std::int16_t> elevationDn_;
};

// Geometry remains in the label's body-fixed axes.  verticesM are Cartesian
// metres from the lunar centre and faces are zero-based vertex indices.
struct Triangulation {
  std::vector<std::array<double, 3>> verticesM;
  std::vector<std::array<std::size_t, 3>> faces;
  unsigned int subdivisionLevel = 0;
  // Great-circle edge length measured before radial topography displacement.
  double maximumEdgeLengthM = 0.0;
  double minimumRadiusM = 0.0;
  double maximumRadiusM = 0.0;
};

// Construct an icosphere whose maximum great-circle edge on the reference
// sphere is no larger than requestedResolutionM.  This avoids the polar
// convergence of latitude/longitude triangulations.
bool BuildIcosphere(double referenceRadiusM, double requestedResolutionM,
                    Triangulation* result, std::string* error);

// Shift each vertex radially to referenceRadiusM + interpolated elevation.
// The native DEM is neither smoothed nor pre-resampled; this function does not
// change vertex directions or connectivity.
bool ApplyLolaTopography(const LolaDem& dem, Triangulation* mesh,
                         std::string* error);

// CEA uses one-based connectivity because the AMPS reader subtracts one.
// Tecplot is an inspection artifact in SI units and is not read by AMPS.
bool WriteCeaSurface(const Triangulation& mesh, const std::string& path,
                     std::string* error);
bool WriteTecplotSurface(const Triangulation& mesh, const std::string& path,
                         std::string* error);

// Production bridge.  In LOLA mode rank zero generates both files and every
// rank reloads the CEA representation through the normal AMPS surface loader
// before the AMR tree is initialized.  Sphere mode leaves the legacy analytic
// boundary unchanged.
bool PrepareProductionSurface(const Runtime::Configuration& configuration,
                              std::string* error);
// These accessors describe process initialization state, not a new runtime
// geometry switch.  Geometry cannot be changed after AMR construction.
bool RealisticSurfaceActive();
double MinimumLoadedSurfaceRadiusM();

}  // namespace Surface
}  // namespace Moon

#endif  // SRCMOON_LUNAR_SURFACE_H
