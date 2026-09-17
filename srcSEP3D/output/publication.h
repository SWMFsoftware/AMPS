// ============================================================================
// Phase-O transactional, self-describing output publication.
//
// A complete sequence is written beneath a new staging directory. Every data
// stream is closed and hashed, then a manifest is written last. One directory
// rename publishes the bundle atomically, so readers see either the previous
// complete sequence or the new complete sequence, never a mixed set of files.
// ============================================================================

#ifndef SEP3D_OUTPUT_PUBLICATION_H
#define SEP3D_OUTPUT_PUBLICATION_H

#include "sampling.h"

#include <cstdint>
#include <map>
#include <string>

namespace SEP3D {
namespace Output {

struct PublicationMetadata {
  std::uint64_t sequence = 0;
  double simulationTimeS = 0.0;
  std::uint64_t snapshotGeneration = 0;
  std::string configurationFingerprint;
  std::string codeIdentity;
  std::string snapshotFingerprint;
};

struct PublicationResult {
  Core::Status status;
  std::string directory;
  std::map<std::string, std::string> artifactHashes;
};

struct ParsedPublication {
  PublicationMetadata metadata;
  std::map<std::string, std::string> artifactHashes;
};

PublicationResult Publish(const std::string& outputDirectory,
                          const std::string& prefix,
                          const PublicationMetadata& metadata,
                          const SamplingSnapshot& snapshot);

// Parse and verify a publication independently of the writer. This checks the
// manifest schema, unit-bearing CSV headers, and every declared artifact hash.
// output is replaced only after the entire directory validates.
Core::Status ParseAndVerifyPublication(const std::string& directory,
                                       ParsedPublication* output);

std::string HashFileFNV1a64(const std::string& path,
                            Core::Status* status = nullptr);

}  // namespace Output
}  // namespace SEP3D

#endif  // SEP3D_OUTPUT_PUBLICATION_H
