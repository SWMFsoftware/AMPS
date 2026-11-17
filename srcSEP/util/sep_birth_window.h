#ifndef SRCSEP_UTIL_SEP_BIRTH_WINDOW_H
#define SRCSEP_UTIL_SEP_BIRTH_WINDOW_H

#include "sep_common_header_path.h"
#include SRCSEP_SEP_COMMON_HEADER(sep_transport_common.h)
#include <cstdint>
#include <cstddef>
#include <type_traits>

namespace SEP { namespace Transport { namespace Birth {

// Numerical PIC clock endpoints and an explicit physical duration [s]. Freeze
// this before the host timestep changes its clock; a newborn may not infer the
// departure epoch from whichever counter mode happens to be compiled.
struct DriverWindow {
  std::uint64_t schema=0;
  double startEpochS=0.0,durationS=0.0,endEpochS=0.0;
};
enum class Origin : std::uint64_t {
  InitialCondition=1,StepStartSource=2,UniformInStepSource=3
};
// A zero-filled/legacy record has no declared birth history. Only explicit
// fresh initial-condition/source APIs may create schema one. The original
// source window remains after its pending first movement has been consumed.
struct Record {
  std::uint64_t schema=0;
  Origin origin=Origin::InitialCondition;
  std::uint64_t pendingFirstMove=0;
  DriverWindow parent;
  double firstMoveDurationS=0.0;
};
static_assert(std::is_trivially_copyable<Record>::value,
    "native birth history must remain a pointer-free value record");
// Schema-one native storage deliberately contains only eight-byte scalar
// slots. Prove that memcpy cannot copy uninitialized padding into particle
// records, reproducibility hashes or checkpoints. A different ABI must define
// a new fieldwise representation rather than silently persisting padding.
static_assert(sizeof(std::uint64_t)==8&&sizeof(double)==8&&sizeof(Origin)==8&&
    sizeof(DriverWindow)==32&&offsetof(DriverWindow,startEpochS)==8&&
    offsetof(DriverWindow,durationS)==16&&offsetof(DriverWindow,endEpochS)==24&&
    sizeof(Record)==64&&offsetof(Record,origin)==8&&
    offsetof(Record,pendingFirstMove)==16&&offsetof(Record,parent)==24&&
    offsetof(Record,firstMoveDurationS)==56,
    "birth schema one requires a padding-free eight-byte scalar layout");
// Time is represented relative to the declared parent end: departure is
// endEpochS-firstMoveDurationS. Keep that duration separate even when their
// floating-point subtraction rounds to endEpochS at a late absolute epoch.
struct FirstMovePlan {
  DriverWindow parent;
  double durationS=0.0;
};
enum class HostStepMode {SingleGlobal,SpeciesGlobal};
struct HostClockPolicy {
  HostStepMode mode=HostStepMode::SingleGlobal;
  std::uint64_t speciesCount=1;
  bool counterAtEnd=true,externallyOwned=false;
};
// Mirrors the audited PIC::TimeStep/SimulationTime::Update capability, not an
// inferred synchronization law. Single-global advances at start(counter OFF)
// or end(ON); species-global advances only for one species with counter ON.
// An externally owned SWMF epoch needs a distinct explicit host-end contract.
Status ValidateHostClockPolicy(const HostClockPolicy& policy);

Status FreezeWindow(double startEpochS,double durationS,DriverWindow* output);
Status ValidateWindow(const DriverWindow& window);
Status InitialCondition(double epochS,Record* output);
Status StepStartBirth(const DriverWindow& window,Record* output);
Status UniformBirth(const DriverWindow& window,double phaseOpen01,Record* output);
Status ValidateRecord(const Record& record);
Status PlanFirstMove(const Record& record,const DriverWindow& frozenWindow,
                     double requestedHostDurationS,FirstMovePlan* output);
// Returns a staged value only. The native adapter must publish it together
// with the successfully advanced phase space and complete RNG state.
Status StageConsumedFirstMove(const Record& record,const FirstMovePlan& plan,
                              Record* output);
// One immutable particle-phase snapshot must cover the entire frozen parent
// window. A shorter interval is an explicit incompatibility, not permission to
// truncate/resample the birth phase or use a second provider realization.
ScalarResult SnapshotRemaining(const FirstMovePlan& plan,double elapsedS,
                               double snapshotFromS,double snapshotUntilS);
// A prescribed late-birth remainder can be shorter than the adaptive
// refinement floor. Accept that whole terminal interval, while retaining the
// declared floor for every stability limit that actually subdivides it.
ScalarResult SelectFirstMoveSubstep(const FirstMovePlan& plan,double elapsedS,
    const std::vector<StepLimit>& limits,double minimumAdaptiveStepS,
    StepDiagnostics* diagnostics=nullptr);

} } }
#endif
