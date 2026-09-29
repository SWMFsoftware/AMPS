#ifndef SEP_COMMON_SEP_FIELD_LINE_BUNDLE_IO_H
#define SEP_COMMON_SEP_FIELD_LINE_BUNDLE_IO_H

#include "sep_field_line_exchange.h"

#include <string>

namespace SEP { namespace FieldLine {

// The target is a directory containing one canonical JSON manifest and one
// checksummed tabular member per line.  Publication renames a complete sibling
// temporary directory, so a failed member write cannot expose a partial bundle.
Core::Status WriteBundleTransactional(const FieldLineSet& set,
                                      const std::string& targetDirectory);
Core::Result<FieldLineSet> ReadBundle(const std::string& targetDirectory);
std::string BundleIdentity(const FieldLineSet& set);

} }  // namespace SEP::FieldLine

#endif  // SEP_COMMON_SEP_FIELD_LINE_BUNDLE_IO_H
