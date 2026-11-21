#include "sep_birth_window.h"
#include <algorithm>
#include <cmath>

namespace SEP { namespace Transport { namespace Birth { namespace {
Status Invalid(const char* message) {
  return Status::Error(StatusCode::InvalidParticleState,message);
}
bool SameWindow(const DriverWindow& a,const DriverWindow& b) {
  return a.schema==b.schema&&a.startEpochS==b.startEpochS&&
      a.durationS==b.durationS&&a.endEpochS==b.endEpochS;
}
Status NewScheduled(const DriverWindow& window,Origin origin,double duration,
                    Record* output) {
  if(!output)return Status::Error(StatusCode::InvalidArgument,"birth record output is required");
  const auto valid=ValidateWindow(window);if(!valid.ok())return valid;
  if(!std::isfinite(duration)||duration<=0.0||duration>window.durationS)
    return Invalid("source birth has no representable positive remaining duration; no redraw or zero-duration substitute is permitted");
  Record candidate;candidate.schema=1;candidate.origin=origin;
  candidate.pendingFirstMove=1;candidate.parent=window;
  candidate.firstMoveDurationS=duration;
  const auto recordValid=ValidateRecord(candidate);if(!recordValid.ok())return recordValid;
  *output=candidate;return Status::Ok();
}
} // namespace

Status ValidateHostClockPolicy(const HostClockPolicy& policy) {
  if(policy.speciesCount==0)return Status::Error(StatusCode::InvalidArgument,"source clock policy has no species");
  if(policy.externallyOwned)return Status::Error(StatusCode::UnsupportedConfiguration,
      "externally owned coupling clock requires an explicit source-window host-end contract; PIC advancement cannot be assumed");
  if(policy.mode==HostStepMode::SingleGlobal)return Status::Ok();
  if(policy.mode!=HostStepMode::SpeciesGlobal)return Status::Error(StatusCode::InvalidArgument,"unknown source clock timestep mode");
  if(policy.speciesCount!=1||!policy.counterAtEnd)
    return Status::Error(StatusCode::UnsupportedConfiguration,
        "compiled PIC species-global clock advances only for one species with the global counter enabled; equal dt alone does not implement synchronization");
  return Status::Ok();
}

Status ValidateWindow(const DriverWindow& window) {
  if(window.schema!=1||!std::isfinite(window.startEpochS)||
      !std::isfinite(window.durationS)||window.durationS<=0.0||
      !std::isfinite(window.endEpochS)||window.endEpochS<=window.startEpochS||
      window.endEpochS!=window.startEpochS+window.durationS)
    return Invalid("source driver window is missing, incompatible, or cannot advance the authoritative numerical clock");
  return Status::Ok();
}
Status ValidateExternalHostWindow(const HostClockPolicy& policy,
    const DriverWindow& window,double observedDepartureS,double nominalDurationS) {
  // This overload proves a concrete host interval instead of weakening the
  // unbound externallyOwned rejection above. Only single-global counter-OFF
  // SWMF has an implemented external phase transition; other policies could
  // double-increment or omit physical time and must remain unsupported.
  if(!policy.externallyOwned||policy.speciesCount==0||
     policy.mode!=HostStepMode::SingleGlobal||policy.counterAtEnd)
    return Status::Error(StatusCode::UnsupportedConfiguration,
        "explicit external source window requires the audited single-global SWMF clock with PIC counter advancement disabled");
  const auto status=ValidateWindow(window);if(!status.ok())return status;
  // Both quantities are supplied in SI seconds by the same initialized model
  // and coupler ledger. Exact comparisons detect stale windows and attempts
  // to substitute a shorter host remainder; accepting an epsilon here would
  // hide a physical duration/sampling-normalization change.
  if(window.startEpochS!=observedDepartureS||window.durationS!=nominalDurationS)
    return Invalid("external source window disagrees with the physical departure or unchanged nominal AMPS timestep");
  return Status::Ok();
}
Status FreezeWindow(double startEpochS,double durationS,DriverWindow* output) {
  if(!output)return Status::Error(StatusCode::InvalidArgument,"driver window output is required");
  DriverWindow candidate;candidate.schema=1;candidate.startEpochS=startEpochS;
  candidate.durationS=durationS;candidate.endEpochS=startEpochS+durationS;
  const auto valid=ValidateWindow(candidate);if(!valid.ok())return valid;
  *output=candidate;return Status::Ok();
}
Status InitialCondition(double epochS,Record* output) {
  if(!output||!std::isfinite(epochS))return Status::Error(StatusCode::InvalidArgument,
      "explicit initial-condition birth requires a finite epoch and output");
  Record candidate;candidate.schema=1;candidate.origin=Origin::InitialCondition;
  candidate.parent.schema=1;candidate.parent.startEpochS=epochS;
  candidate.parent.endEpochS=epochS;
  *output=candidate;return Status::Ok();
}
Status StepStartBirth(const DriverWindow& window,Record* output) {
  return NewScheduled(window,Origin::StepStartSource,window.durationS,output);
}
Status UniformBirth(const DriverWindow& window,double phaseOpen01,Record* output) {
  if(!std::isfinite(phaseOpen01)||phaseOpen01<=0.0||phaseOpen01>=1.0)
    return Status::Error(StatusCode::InvalidArgument,"uniform birth phase must lie strictly inside (0,1)");
  // Multiplication by the remaining fraction is authoritative. Computing an
  // absolute rounded birth epoch first and subtracting it can erase a real
  // late birth's entire first interval, particularly at epoch 1e12 seconds.
  return NewScheduled(window,Origin::UniformInStepSource,
      (1.0-phaseOpen01)*window.durationS,output);
}
Status ValidateRecord(const Record& record) {
  if(record.schema!=1||record.pendingFirstMove>1)
    return Invalid("missing or incompatible birth-history schema; a previous realization cannot be initialized as a fresh birth");
  if(record.origin==Origin::InitialCondition) {
    if(record.pendingFirstMove||record.parent.schema!=1||
        !std::isfinite(record.parent.startEpochS)||
        record.parent.endEpochS!=record.parent.startEpochS||
        record.parent.durationS!=0.0||record.firstMoveDurationS!=0.0)
      return Invalid("initial-condition birth history contains a source window or pending source movement");
    return Status::Ok();
  }
  if(record.origin!=Origin::StepStartSource&&record.origin!=Origin::UniformInStepSource)
    return Invalid("unknown birth origin");
  const auto valid=ValidateWindow(record.parent);if(!valid.ok())return valid;
  if(!std::isfinite(record.firstMoveDurationS)||record.firstMoveDurationS<=0.0||
      record.firstMoveDurationS>record.parent.durationS||
      (record.origin==Origin::StepStartSource&&record.firstMoveDurationS!=record.parent.durationS)||
      (record.origin==Origin::UniformInStepSource&&record.firstMoveDurationS>=record.parent.durationS))
    return Invalid("scheduled source birth has invalid remaining first-move duration");
  return Status::Ok();
}
Status PlanFirstMove(const Record& record,const DriverWindow& frozenWindow,
                     double requestedHostDurationS,FirstMovePlan* output) {
  if(!output)return Status::Error(StatusCode::InvalidArgument,"first-move plan output is required");
  const auto valid=ValidateRecord(record);if(!valid.ok())return valid;
  const auto windowValid=ValidateWindow(frozenWindow);if(!windowValid.ok())return windowValid;
  if(!record.pendingFirstMove||record.origin==Origin::InitialCondition)
    return Invalid("particle has no pending scheduled first move");
  if(!SameWindow(record.parent,frozenWindow)||requestedHostDurationS!=frozenWindow.durationS)
    return Status::Error(StatusCode::UnsupportedConfiguration,
        "pending birth requires its original explicit driver window and common host timestep; unequal species timestep synchronization is unsupported");
  FirstMovePlan candidate;candidate.parent=record.parent;
  candidate.durationS=record.firstMoveDurationS;
  *output=candidate;return Status::Ok();
}
Status StageConsumedFirstMove(const Record& record,const FirstMovePlan& plan,Record* output) {
  if(!output)return Status::Error(StatusCode::InvalidArgument,"consumed birth record output is required");
  FirstMovePlan expected;const auto valid=PlanFirstMove(record,plan.parent,
      plan.parent.durationS,&expected);if(!valid.ok())return valid;
  if(plan.durationS!=expected.durationS)
    return Invalid("first-move commit does not match the planned birth duration");
  Record candidate=record;candidate.pendingFirstMove=0;
  *output=candidate;return Status::Ok();
}
ScalarResult SnapshotRemaining(const FirstMovePlan& plan,double elapsedS,
                               double snapshotFromS,double snapshotUntilS) {
  ScalarResult result;const auto valid=ValidateWindow(plan.parent);
  if(!valid.ok()){result.status=valid;return result;}
  if(!std::isfinite(plan.durationS)||plan.durationS<=0.0||
      plan.durationS>plan.parent.durationS||!std::isfinite(elapsedS)||
      elapsedS<0.0||elapsedS>plan.durationS||!std::isfinite(snapshotFromS)||
      std::isnan(snapshotUntilS)||snapshotFromS>plan.parent.startEpochS||
      snapshotUntilS<plan.parent.endEpochS) {
    result.status=Invalid("snapshot must cover the complete frozen birth window and elapsed time must remain within the planned first move");
    return result;
  }
  result.value=(snapshotUntilS-plan.parent.endEpochS)+(plan.durationS-elapsedS);
  if(std::isnan(result.value)||result.value<0.0||
      (std::isfinite(snapshotUntilS)&&!std::isfinite(result.value))) {
    result.status=Invalid("birth-relative snapshot interval is invalid");return result;
  }
  result.status=Status::Ok();return result;
}
ScalarResult SelectFirstMoveSubstep(const FirstMovePlan& plan,double elapsedS,
    const std::vector<StepLimit>& limits,double minimumAdaptiveStepS,
    StepDiagnostics* diagnostics) {
  // This also validates the plan, elapsed interval, and explicit parent clock
  // window without consulting or mutating the global PIC simulation clock.
  const auto remaining=SnapshotRemaining(plan,elapsedS,plan.parent.startEpochS,
      plan.parent.endEpochS);
  if(!remaining.status.ok())return remaining;
  if(!std::isfinite(minimumAdaptiveStepS)||minimumAdaptiveStepS<0.0) {
    ScalarResult invalid;invalid.status=Status::Error(StatusCode::InvalidArgument,
        "minimum adaptive substep must be finite and nonnegative");return invalid;
  }
  // The ordinary selector still validates every requested physical/stability
  // limit. For a remaining interval below the refinement floor, its only
  // admissible result is that complete interval: any smaller limit remains
  // below the adjusted threshold and raises StepUnderflow as before.
  return SelectSubstep(remaining.value,limits,
      std::min(minimumAdaptiveStepS,remaining.value),diagnostics);
}
} } }
