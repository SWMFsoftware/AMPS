#ifndef SEP_UTIL_SEP_TRANSACTIONAL_OUTPUT_H
#define SEP_UTIL_SEP_TRANSACTIONAL_OUTPUT_H

#include "sep_common_header_path.h"
#include SRCSEP_SEP_COMMON_HEADER(sep_transport_common.h)

#include <cstdint>
#include <string>

namespace SEP {
namespace Output {

struct ArtifactMetadata {
  std::string schemaVersion;
  std::uint64_t recordCount = 0;
  std::string configurationFingerprint;
};

struct WriteResult {
  Transport::Status status;
  std::string finalPath;
  std::string checksum;
};

// Reject absolute/traversing path components before any filesystem mutation.
// The root itself is supplied by the trusted driver configuration.
Transport::Status ValidateRelativePath(const std::string& relativePath);
Transport::Status EnsureOutputDirectory(const std::string& trustedRoot,
                                        const std::string& relativeDirectory);

// Write payload and an authenticated completion manifest transactionally.  The
// payload is flushed and fsync'd before the temporary file is atomically renamed;
// readers can reject artifacts with a missing/mismatched .manifest sidecar.
WriteResult WriteTransactional(const std::string& trustedRoot,
                               const std::string& relativePath,
                               const std::string& payload,
                               const ArtifactMetadata& metadata);
Transport::Status ValidateCompletedArtifact(const std::string& finalPath,
                                            std::string* payload);
std::string Checksum64(const std::string& payload);

}  // namespace Output
}  // namespace SEP

#endif  // SEP_UTIL_SEP_TRANSACTIONAL_OUTPUT_H
