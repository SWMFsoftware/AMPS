#ifndef SRCMOON_MOON_INPUT_H
#define SRCMOON_MOON_INPUT_H

#include <string>

namespace Moon {
namespace Runtime {

// Selects the geometric boundary used by the particle mover.  Sphere keeps
// the historical cInternalSphericalData boundary; Lola replaces only the
// collision geometry with a generated triangulation.  A logical spherical
// grid is retained in both modes for the existing source and inventory code.
enum class SurfaceGeometry {
  Sphere,
  Lola
};

// All lengths are SI.  Paths are absolute after parsing; relative paths are
// resolved against the directory containing the .in file, never against the
// process working directory.
struct Configuration {
  SurfaceGeometry surfaceGeometry = SurfaceGeometry::Sphere;
  std::string inputFile;             // Absolute path of the parsed receipt.
  std::string spiceRoot;             // Root containing cspice/ and Kernels/.
  std::string spiceToolkitDirectory; // <spiceRoot>/cspice.
  std::string spiceKernelDirectory;  // <spiceRoot>/Kernels.
  std::string lolaProductId;         // Must equal PRODUCT_ID in the PDS label.
  std::string lolaImageFile;         // Detached binary PDS image.
  std::string lolaLabelFile;         // Authoritative detached PDS label.
  std::string surfaceCeaFile;        // Generated AMPS CEA-long mesh.
  std::string surfaceTecplotFile;    // Generated Tecplot FETRIANGLE mesh.
  double surfaceMeshResolutionM = 0.0; // Maximum pre-topography arc length.
};

// Parse only the unique ``moon`` section from a shared AMPS .in file.  The
// operation is transactional: result is unchanged when false is returned.
bool ParseApplicationInput(const std::string& path, Configuration* result,
                           std::string* error);

// Install the immutable process configuration before amps_init().  Installation
// is intentionally one-shot: silently replacing paths or geometry after MPI
// and mesh initialization would make the run receipt disagree with the model.
bool InstallConfiguration(const Configuration& configuration,
                          std::string* error);
bool HasConfiguration();
// When no application input was provided, return a read-only default object.
// Callers must test HasConfiguration() before interpreting its empty paths.
const Configuration& GetConfiguration();

const char* SurfaceGeometryName(SurfaceGeometry geometry);

}  // namespace Runtime
}  // namespace Moon

#endif  // SRCMOON_MOON_INPUT_H
