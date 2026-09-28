#ifndef SEP_CORONAL_CME_CONFIGURATION_PARSER_H
#define SEP_CORONAL_CME_CONFIGURATION_PARSER_H

#include "sep_coronal_cme/model_configuration.h"

#include <string>

namespace SEP {
namespace CoronalCME {

// Parse a complete application input using the mandatory two-pass dispatch.
// Pass one recognizes only section/assignment syntax and locates the unique
// run.schema_version.  Pass two then applies exactly that version's grammar.
// Schema 1--4 content is returned byte-for-byte to the frozen application
// parser; this shared library never reinterprets a legacy deck.
Core::Result<VersionedConfiguration> ParseConfiguration(
    const std::string& inputBytes);

// Validate a configuration that was populated by a typed coupled host.  This
// follows the same physical rules as the text parser and performs no I/O.
Core::Status ValidateConfiguration(const ModelConfiguration& configuration);

// Canonical schema-5 identity.  The serializer is sorted by dotted key and
// length-prefixes names/values, avoiding locale, whitespace, and delimiter
// ambiguities.  Provenance paths and comments are intentionally excluded;
// their independently verified content checksums remain included.
std::string ComputePhysicsFingerprint(
    const std::map<std::string, std::string>& normalizedAssignments);

// Stable semantic comparison used by restart gates.  A mismatch returns a
// DataIntegrityFailure before a provider, AMPS node, or particle is allocated.
Core::Status CheckRestartIdentity(const std::string& expectedFingerprint,
                                  const std::string& actualFingerprint,
                                  const std::string& identityCategory);

// Reports whether the selected provider is implemented in the requested
// development stage. Parsing and capability activation are intentionally
// separate: Stage 0 freezes later selectors, while this call prevents them
// from falling through to a currently available but physically different
// provider.
Core::Status CheckCapabilityAvailability(
    const ModelConfiguration& configuration, int completedStage);

}  // namespace CoronalCME
}  // namespace SEP

#endif  // SEP_CORONAL_CME_CONFIGURATION_PARSER_H
