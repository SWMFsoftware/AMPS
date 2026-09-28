#ifndef SEP_CORONAL_CME_SHA256_H
#define SEP_CORONAL_CME_SHA256_H

#include <string>

namespace SEP {
namespace CoronalCME {
namespace Internal {

// Dependency-free SHA-256 used for configuration/restart identities.  Asset
// *contents* are expected to be checksummed by their ingest layer; this helper
// hashes the canonical resolved manifest bytes that bind those checksums to
// the selected equations and numerical controls.
std::string Sha256Hex(const std::string& bytes);

}  // namespace Internal
}  // namespace CoronalCME
}  // namespace SEP

#endif  // SEP_CORONAL_CME_SHA256_H
