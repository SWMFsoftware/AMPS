#ifndef SRCSEP_FIELD_LINE_BUNDLE_ADAPTER_H
#define SRCSEP_FIELD_LINE_BUNDLE_ADAPTER_H

#include "sep_field_line_bundle_io.h"
#include "sep_status.h"

#include <functional>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace SEP1D { namespace Adapters {

// Application code provides the publisher that converts a validated neutral
// NodeState into the existing SEP::Background::SnapshotStore representation.
// Keeping that callback here prevents this adapter from duplicating or
// shadowing sep_common headers in an application-local util directory.
using SnapshotPublisher = std::function<SEP::Core::Status(
    const SEP::FieldLine::FieldLineSet&)>;

// Typed representation of one repeated [line_mesh.ID] section.  The stable
// mesh ID is configuration identity, whereas stableLineId is the immutable
// physical identity stored in the imported bundle.  Keeping those two names
// separate prevents vector position or request order from becoming physical
// identity during restart and ensemble comparisons.
struct LineMeshRequest {
  std::string stableMeshId;
  std::string stableLineId;
  double startArcLengthM = 0.0;
  double endArcLengthM = 0.0;
  std::size_t pointCount = 0;
};

// Convert the normalized application-parser assignments into typed repeated
// line meshes.  Keys outside line_mesh.* are ignored so this routine can be
// called on the complete normalized srcSEP input.  Each concrete record must
// contain exactly the five schema fields documented in the Stage-10 README;
// missing fields, unknown fields, duplicate physical lines, unit-bearing
// numbers, non-finite values, reversed intervals, and fewer than two points
// are rejected before a background snapshot or particle is allocated.
SEP::Core::Result<std::vector<LineMeshRequest>> ParseLineMeshRequests(
    const std::map<std::string, std::string>& normalizedAssignments);

// Resolve every request by exact stable line ID and conservatively resample it
// with sep_common's bounded, positive, one-sided interpolation.  There is no
// nearest-line fallback and no extrapolation beyond the exported interval.
SEP::Core::Result<std::vector<SEP::FieldLine::LineRecord>> BuildLineMeshes(
    const SEP::FieldLine::FieldLineSet& bundle,
    const std::vector<LineMeshRequest>& requests);

class FieldLineBundleAdapter {
 public:
  SEP::Core::Status LoadAndPublish(const std::string& directory,
                                   const SnapshotPublisher& publisher);
  SEP::Core::Result<SEP::FieldLine::NodeState> Evaluate(
      const std::string& stableLineId, double arcLengthM) const;
  const std::string& BundleIdentity() const noexcept { return identity_; }

 private:
  std::shared_ptr<const SEP::FieldLine::FieldLineSet> bundle_;
  std::string identity_;
};

} }  // namespace SEP1D::Adapters

#endif  // SRCSEP_FIELD_LINE_BUNDLE_ADAPTER_H
