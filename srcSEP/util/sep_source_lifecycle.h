#ifndef SRCSEP_UTIL_SEP_SOURCE_LIFECYCLE_H
#define SRCSEP_UTIL_SEP_SOURCE_LIFECYCLE_H

#include "sep_common_header_path.h"
#include SRCSEP_SEP_COMMON_HEADER(sep_transport_common.h)

#include <cstdint>
#include <cstdio>
#include <string>
#include <vector>

namespace SEP { namespace SourceLifecycle {

// Each counter belongs to a physical source channel, independent of its MPI
// owner. The fixed-line bridge verifies the complete ordered initial geometry;
// remapping physical lines requires an explicit stable-line mapping instead.
enum class Kind : std::uint64_t {
  FixedBeginning = 1, Shock = 2, InnerBoundary = 3,
  SolarWindPopulation = 4, SolarWindBeginning = 5, MaxwellianBeginning = 6,
  PopulationControl = 7
};
struct Channel {
  Kind kind = Kind::FixedBeginning;
  std::uint64_t line = 0, species = 0;
};
struct Event {
  Channel channel;
  std::uint64_t campaign = 0, index = 0;
};
struct Counter { Channel channel; std::uint64_t next = 0; };
struct State {
  std::uint64_t schema = 1, campaign = 0;
  std::string fixedLineLayout;
  std::vector<Counter> counters;
};
class Scheduler {
 public:
  Transport::Status Initialize(std::uint64_t campaign,
                               const std::string& fixedLineLayout);
  Transport::Status Begin(const Channel& channel, Event* event);
  State Capture() const;
  Transport::Status Restore(const State& state, std::uint64_t campaign,
                            const std::string& fixedLineLayout);
  void Swap(Scheduler* other);
 private:
  bool initialized_ = false;
  State state_;
};

// Identity contains no sampled momentum, coordinates, rank, or wall time.
// Purpose-tagged source streams use this same immutable birth tuple.
std::uint64_t ParticleIdentity(const Event& event, std::uint64_t macro);
Transport::KeyedRandomStream Random(const Event& event, std::uint64_t macro,
                                    std::uint64_t purpose);
Transport::Status WriteState(FILE* file, const State& state);
Transport::Status ReadState(FILE* file, State* state);

// Native application bridge (implemented in field_line.cpp). Fixed standalone
// geometry uses its exact coordinate namespace; the supported SWMF origin-grid
// producer must supply explicit global line ids and immutable grid metadata.
// Moving vertices alone do not replace the latter's source identity. Missing
// mappings or unsupported remaps reject before reserving a new source event.
// All ranks must call BeginNativeEvent in the same physical dispatch, before
// an owner-rank return. These hooks are ready for a combined field-line
// checkpoint; registering them alone cannot save field-line particle lists.
Event BeginNativeEvent(Kind kind, int species, int line);
Transport::Status InitializeNativeFresh();
State CaptureNativeState();
Transport::Status RestoreNativeState(const State& state);
Transport::Status ValidateNativeState(const State& state);
Transport::Status PrepareNativeRestore(const State& state, Scheduler* staged);
void CommitNativeRestore(Scheduler* staged);

} }  // namespace SEP::SourceLifecycle
#endif
