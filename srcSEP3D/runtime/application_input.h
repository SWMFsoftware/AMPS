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
  // This maps to SourceOptions::samplesPerStep.  The existing srcSEP3D source
  // contract interprets the value per compiled species; the parser must not
  // silently change that established particle-accounting rule.
  std::uint64_t particlesPerIteration = 0;

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
