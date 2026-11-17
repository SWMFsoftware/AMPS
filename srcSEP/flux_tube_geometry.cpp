#include "sep.h"
#include "flux_tube_checkpoint_state.h"
#include "transport_common.h"
#include "util/sep_background_runtime.h"
#include "util/sep_run_configuration.h"

#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>
#include <limits>
#include <utility>

namespace {

// The historical implementation implicitly assigned a one-metre radius to the
// first point.  Preserve that normalization explicitly as an SI reference area
// while replacing the dimensionally incorrect "radius grows as r squared"
// rule.  Applications with a measured tube area should set it before startup.
double g_reference_area_m2 = Pi;
SEP::FieldLine::FluxTubeGeometry::ExplicitAreaProfile g_explicit_area = NULL;
std::vector<SEP::FieldLine::FluxTubeGeometryCore::MagneticFluxRecord>
    g_field_line_flux;
// This serial invalidates outstanding software restore plans; it is never a
// physical key, source event, random seed or age and therefore is not persisted.
std::uint64_t g_flux_revision=0;
void RequireFluxRevision() {
  if(g_flux_revision==std::numeric_limits<std::uint64_t>::max())
    throw std::runtime_error("flux-owner mutation serial exhausted before state change");
}

double MagneticMagnitudeT(PIC::FieldLine::cFieldLineVertex* vertex) {
  if (vertex == NULL) return 0.0;
  const double* B_T = vertex->GetMagneticField();
  if (B_T == NULL) return 0.0;
  const double magnitude_T =
      std::sqrt(B_T[0]*B_T[0] + B_T[1]*B_T[1] + B_T[2]*B_T[2]);
  return std::isfinite(magnitude_T) ? magnitude_T : 0.0;
}

double ExplicitAreaM2(const double* x_m, int iFieldLine) {
  if (g_explicit_area == NULL || x_m == NULL) return 0.0;
  const double area_m2 = g_explicit_area(x_m, iFieldLine);
  if (!std::isfinite(area_m2) || area_m2 <= 0.0) {
    throw std::runtime_error(
        "FluxTubeGeometry explicit area profile returned an invalid area");
  }
  return area_m2;
}

}  // namespace

void SEP::FieldLine::FluxTubeGeometry::SetMagneticFluxWb(
    int iFieldLine,double magnetic_flux_Wb,const char* provenance,
    unsigned long long generation) {
  if (iFieldLine<0 || !std::isfinite(magnetic_flux_Wb) ||
      magnetic_flux_Wb<=0.0 || provenance==NULL || provenance[0]=='\0')
    throw std::invalid_argument("FluxTubeGeometry magnetic-flux record is invalid");
  FluxTubeGeometryCore::MagneticFluxRecord record;
  record.magneticFluxWb=magnetic_flux_Wb;
  record.generation=static_cast<std::uint64_t>(generation);
  record.provenance=provenance;
  const SEP::Transport::Status status=
      FluxTubeGeometryCore::ValidateMagneticFluxRecord(record);
  if (!status.ok()) throw std::invalid_argument(status.message);
  RequireFluxRevision();
  if(g_field_line_flux.size()<=static_cast<std::size_t>(iFieldLine)) {
    // Growth and provenance allocation are staged privately. An allocation
    // failure cannot leave new holes or partially replaced flux in the owner.
    auto staged=g_field_line_flux;staged.resize(static_cast<std::size_t>(iFieldLine)+1);
    staged[static_cast<std::size_t>(iFieldLine)]=std::move(record);
    g_field_line_flux.swap(staged);
  } else std::swap(g_field_line_flux[static_cast<std::size_t>(iFieldLine)],record);
  ++g_flux_revision;
}

double SEP::FieldLine::FluxTubeGeometry::GetMagneticFluxWb(int iFieldLine) {
  if (iFieldLine<0 || iFieldLine>=PIC::FieldLine::nFieldLine)
    throw std::invalid_argument("FluxTubeGeometry field-line id is invalid");
  if (g_field_line_flux.size()<=static_cast<std::size_t>(iFieldLine) ||
      !FluxTubeGeometryCore::ValidateMagneticFluxRecord(
          g_field_line_flux[static_cast<std::size_t>(iFieldLine)]).ok()) {
    PIC::FieldLine::cFieldLineSegment* first=
        PIC::FieldLine::FieldLinesAll[iFieldLine].GetFirstSegment();
    const double seed_B=first ? MagneticMagnitudeT(first->GetBegin()) : 0.0;
    if (!(seed_B>0.0))
      throw std::runtime_error("FluxTubeGeometry cannot derive Phi_i without seed |B|");
    // Compatibility input is converted exactly once to an owned flux record.
    // Every downstream area/volume consumer then uses Phi_i/|B| and therefore
    // cannot choose a different reference vertex or normalization.
    SetMagneticFluxWb(iFieldLine,g_reference_area_m2*seed_B,
        "configured seed area times first-vertex magnetic field",0);
  }
  return g_field_line_flux[static_cast<std::size_t>(iFieldLine)].magneticFluxWb;
}

void SEP::FieldLine::FluxTubeGeometry::SetReferenceAreaM2(
    double reference_area_m2) {
  if (!std::isfinite(reference_area_m2) || reference_area_m2 <= 0.0) {
    throw std::invalid_argument(
        "FluxTubeGeometry reference area must be finite and positive");
  }
  if(g_reference_area_m2!=reference_area_m2) {
    RequireFluxRevision();g_reference_area_m2=reference_area_m2;++g_flux_revision;
  }
}

void SEP::FieldLine::FluxTubeGeometry::SetExplicitAreaProfile(
    ExplicitAreaProfile profile) {
  if (profile == NULL) {
    throw std::invalid_argument(
        "FluxTubeGeometry explicit area profile cannot be null");
  }
  if(g_explicit_area!=profile){RequireFluxRevision();g_explicit_area=profile;++g_flux_revision;}
}

void SEP::FieldLine::FluxTubeGeometry::ClearExplicitAreaProfile() {
  if(g_explicit_area){RequireFluxRevision();g_explicit_area=NULL;++g_flux_revision;}
}

namespace SEP { namespace FieldLine { namespace FluxTubeGeometry { namespace {
Transport::Status CheckpointError(const std::string& text) {
  return Transport::Status::Error(Transport::StatusCode::InvalidParticleState,text);
}
Transport::Status CheckpointQuiescent() {
  if(g_explicit_area)return Transport::Status::Error(Transport::StatusCode::UnsupportedConfiguration,
      "flux checkpoint cannot serialize an opaque area callback; explicit typed area ownership is required");
  if(Transport::PICAdapter::NativeDriverWindowActive())
    return CheckpointError("flux checkpoint requires a completed native driver window");
  const auto snapshot=Background::SnapshotStore::Instance().Current();
  if(snapshot)Background::SnapshotStore::Instance().AssertProviderMayWrite(snapshot->provider());
  return Transport::Status::Ok();
}
} // namespace
struct CheckpointRestorePlan::Implementation {
  FluxTubeState::State state;
  std::uint64_t expectedRevision=0;
  bool ready=false;
};
CheckpointRestorePlan::CheckpointRestorePlan()=default;
CheckpointRestorePlan::~CheckpointRestorePlan()=default;
CheckpointRestorePlan::CheckpointRestorePlan(CheckpointRestorePlan&&)noexcept=default;
CheckpointRestorePlan& CheckpointRestorePlan::operator=(CheckpointRestorePlan&&)noexcept=default;
bool HasOpaqueAreaProfile(){return g_explicit_area!=NULL;}
Transport::Status ValidateCheckpointState(const FluxTubeState::State& state) {
  try {
    const auto quiet=CheckpointQuiescent();if(!quiet.ok())return quiet;
    return FluxTubeState::Validate(state,static_cast<std::size_t>(PIC::FieldLine::nFieldLine),
        Run::Active().get().fieldLineSeedAreaM2.value);
  }catch(const std::exception& error){return CheckpointError(error.what());}
}
Transport::Status CaptureCheckpointState(FluxTubeState::State* output) {
  if(!output)return CheckpointError("flux checkpoint output is absent");
  try {
    const auto quiet=CheckpointQuiescent();if(!quiet.ok())return quiet;
    FluxTubeState::State candidate;candidate.referenceAreaM2=g_reference_area_m2;
    candidate.records=g_field_line_flux;
    const auto valid=ValidateCheckpointState(candidate);if(!valid.ok())return valid;
    *output=std::move(candidate);return Transport::Status::Ok();
  }catch(const std::exception& error){return CheckpointError(error.what());}
}
Transport::Status StageCheckpointRestore(const FluxTubeState::State& state,CheckpointRestorePlan* plan) {
  if(!plan)return CheckpointError("flux restore-plan output is absent");
  plan->implementation_.reset();
  const auto valid=ValidateCheckpointState(state);if(!valid.ok())return valid;
  try {
    RequireFluxRevision();
    auto staged=std::unique_ptr<CheckpointRestorePlan::Implementation>(new CheckpointRestorePlan::Implementation);
    staged->state=state;staged->expectedRevision=g_flux_revision;
    plan->implementation_=std::move(staged);return Transport::Status::Ok();
  }catch(const std::exception& error){return CheckpointError(error.what());}
}
Transport::Status ValidateCheckpointReady(const CheckpointRestorePlan& plan) {
  if(!plan.implementation_)return CheckpointError("flux restore plan is absent");
  auto& staged=*plan.implementation_;staged.ready=false;
  const auto valid=ValidateCheckpointState(staged.state);if(!valid.ok())return valid;
  if(staged.expectedRevision!=g_flux_revision||g_flux_revision==std::numeric_limits<std::uint64_t>::max())
    return CheckpointError("flux owner changed after restore staging");
  staged.ready=true;return Transport::Status::Ok();
}
bool CommitCheckpointRestore(CheckpointRestorePlan* plan)noexcept {
  if(!plan||!plan->implementation_)return false;
  auto& staged=*plan->implementation_;
  if(!staged.ready||staged.expectedRevision!=g_flux_revision||g_explicit_area||
      g_flux_revision==std::numeric_limits<std::uint64_t>::max())return false;
  try {const auto quiet=CheckpointQuiescent();if(!quiet.ok())return false;}
  catch(...){return false;}
  try {
    if(PIC::FieldLine::nFieldLine<0||staged.state.records.size()>static_cast<std::size_t>(PIC::FieldLine::nFieldLine)||
        staged.state.referenceAreaM2!=Run::Active().get().fieldLineSeedAreaM2.value)return false;
  }catch(...){return false;}
  // The enclosing all-rank transaction has validated shape/configuration and
  // joined all writers. Commit swaps already allocated values and consumes the
  // plan; deleting the old vector allocates nothing and preserves no codeptr.
  std::swap(g_reference_area_m2,staged.state.referenceAreaM2);
  g_field_line_flux.swap(staged.state.records);++g_flux_revision;
  plan->implementation_.reset();return true;
}
} } }

double SEP::FieldLine::FluxTubeGeometry::AreaAtVertexM2(
    PIC::FieldLine::cFieldLineVertex* vertex,
    int iFieldLine) {
  if (vertex == NULL || iFieldLine < 0 ||
      iFieldLine >= PIC::FieldLine::nFieldLine) {
    throw std::invalid_argument("FluxTubeGeometry vertex or field-line id is invalid");
  }

  const double local_abs_B_T = MagneticMagnitudeT(vertex);
  if (local_abs_B_T > 0.0) {
    FluxTubeGeometryCore::MagneticFluxRecord record;
    record.magneticFluxWb=GetMagneticFluxWb(iFieldLine);
    record.provenance="production field-line flux registry";
    return FluxTubeGeometryCore::AreaFromMagneticFlux(
        record,Units::MagneticFieldT(local_abs_B_T)).Value();
  }

  // A missing magnetic magnitude does not justify silently fabricating an r^2
  // profile.  Requiring an explicit area callback makes the physical closure a
  // visible application configuration decision.
  const double explicit_area_m2 = ExplicitAreaM2(vertex->GetX(), iFieldLine);
  if (explicit_area_m2 > 0.0) return explicit_area_m2;
  throw std::runtime_error(
      "FluxTubeGeometry requires finite |B| or an explicit area profile");
}

double SEP::FieldLine::FluxTubeGeometry::AreaAtSegmentFractionM2(
    PIC::FieldLine::cFieldLineSegment* segment,
    int iFieldLine,
    double fraction) {
  if (segment == NULL || !std::isfinite(fraction) || fraction < 0.0 ||
      fraction > 1.0 || iFieldLine < 0 ||
      iFieldLine >= PIC::FieldLine::nFieldLine) {
    throw std::invalid_argument("FluxTubeGeometry segment fraction is invalid");
  }
  double local_B_T[3] = {0.0,0.0,0.0};
  segment->GetMagneticField(fraction, local_B_T);
  const double local_abs_B_T = std::sqrt(
      local_B_T[0]*local_B_T[0] + local_B_T[1]*local_B_T[1] +
      local_B_T[2]*local_B_T[2]);

  if (std::isfinite(local_abs_B_T) && local_abs_B_T > 0.0) {
    FluxTubeGeometryCore::MagneticFluxRecord record;
    record.magneticFluxWb=GetMagneticFluxWb(iFieldLine);
    record.provenance="production field-line flux registry";
    return FluxTubeGeometryCore::AreaFromMagneticFlux(
        record,Units::MagneticFieldT(local_abs_B_T)).Value();
  }

  double x_m[3] = {0.0,0.0,0.0};
  segment->GetCartesian(x_m, fraction);
  const double explicit_area_m2 = ExplicitAreaM2(x_m, iFieldLine);
  if (explicit_area_m2 > 0.0) return explicit_area_m2;
  throw std::runtime_error(
      "FluxTubeGeometry requires finite |B| or an explicit area profile");
}

double SEP::FieldLine::FluxTubeGeometry::PartialSegmentVolumeM3(
    PIC::FieldLine::cFieldLineSegment* segment,
    int iFieldLine,
    double begin_fraction,
    double end_fraction) {
  if (segment == NULL) {
    throw std::invalid_argument("FluxTubeGeometry segment cannot be null");
  }
  if (!std::isfinite(begin_fraction) || !std::isfinite(end_fraction) ||
      begin_fraction < 0.0 || end_fraction > 1.0 ||
      end_fraction <= begin_fraction) {
    throw std::invalid_argument(
        "FluxTubeGeometry fractions must satisfy 0 <= begin < end <= 1");
  }

  // Simpson integration samples the actual flux-conserving area at both ends
  // and the midpoint.  It therefore converges for curved B(s) without replacing
  // A|B|=constant by an unrelated conical-radius assumption.
  const double middle_fraction = 0.5*(begin_fraction+end_fraction);
  const double begin_area_m2 =
      AreaAtSegmentFractionM2(segment,iFieldLine,begin_fraction);
  const double middle_area_m2 =
      AreaAtSegmentFractionM2(segment,iFieldLine,middle_fraction);
  const double end_area_m2 =
      AreaAtSegmentFractionM2(segment,iFieldLine,end_fraction);
  const double partial_length_m =
      (end_fraction-begin_fraction)*segment->GetLength();
  return FluxTubeGeometryCore::IntegrateSimpsonArea(
      Units::AreaM2(begin_area_m2), Units::AreaM2(middle_area_m2),
      Units::AreaM2(end_area_m2),
      Units::LengthM(partial_length_m)).Value();
}

double SEP::FieldLine::FluxTubeGeometry::SegmentVolumeM3(
    PIC::FieldLine::cFieldLineSegment* segment,
    int iFieldLine) {
  return PartialSegmentVolumeM3(segment, iFieldLine, 0.0, 1.0);
}

double SEP::FieldLine::FluxTubeGeometry::SweptVolumeM3(
    PIC::FieldLine::cFieldLineSegment* segment,
    int iFieldLine,
    double fraction,
    double normal_speed_m_s,
    double time_s) {
  return FluxTubeGeometryCore::SweptVolume(
      Units::AreaM2(AreaAtSegmentFractionM2(segment, iFieldLine, fraction)),
      Units::SpeedMPerS(normal_speed_m_s), Units::TimeS(time_s)).Value();
}
