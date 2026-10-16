#ifndef SHIELDSIM_CONTENT_HASH_HH
#define SHIELDSIM_CONTENT_HASH_HH

/* Stable SHA-256 content identifiers used in archived Stage-1 provenance. */

#include <string>

namespace shield_sim {

// Returns the lower-case identifier "sha256:<64 hexadecimal digits>".
std::string Sha256(const std::string& bytes);

} // namespace shield_sim

#endif
