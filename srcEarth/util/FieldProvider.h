#ifndef _SRC_EARTH_UTIL_FIELD_PROVIDER_H_
#define _SRC_EARTH_UTIL_FIELD_PROVIDER_H_

//======================================================================================
// FieldProvider.h
//======================================================================================
// Backend-neutral contract for every magnetic/electric field used by the Earth
// cutoff, access, flux, and spectrum calculations.
//
// Roadmap Step 3 requires the particle tracer to be independent of the origin of the
// field.  A direct DIPOLE/IGRF/Tsyganenko evaluator and a cell-centred PIC::CPLR/SWMF
// mesh therefore expose the same two concepts:
//
//   IFieldProvider  - a mutable source that can freeze a requested epoch/driver state;
//   IFieldSnapshot  - the immutable, read-only object sampled by one trajectory batch.
//
// The common contract is intentionally independent of PIC, MPI, Geopack, and the
// Tsyganenko interfaces.  It can be compiled and unit-tested as a small C++ library.
// Backend adapters may contain those dependencies, but trajectory/product code must
// see only the metadata and sample status declared here.
//
// Phase-1 time semantics
// ----------------------
// A snapshot is static or quasi-static: B (and optional E) are held fixed from the
// first trajectory in a batch through the last.  A series of such snapshots is NOT a
// time-dependent characteristic.  The immutable flag and snapshot identifier make
// this limitation visible in output provenance and allow cutoff and flux products to
// verify that they used exactly the same field state.
//
// Unit/frame contract
// -------------------
// Positions are metres, B is tesla, and E is volt/metre.  The coordinate frame is
// explicit and is currently GSM for production Earth calculations.  No adapter is
// allowed to return nT, Re, km, or an implicit frame while claiming this contract.
//======================================================================================

#include <cmath>
#include <cstdint>
#include <iomanip>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>

namespace Earth {
namespace Field {

enum class CoordinateFrame {
  Unknown,
  GSM,
  GEO,
  GSE,
  SM
};

enum class InterpolationMode {
  Unknown,
  DirectAnalytic,
  DirectEmpirical,
  CellCenteredLinear,
  CellCenteredLinearDerivedElectric
};

// FieldSampleStatus deliberately separates physical domain exclusion from numerical
// or configuration failures.  A caller may classify OutsideDomain as an outer-boundary
// trajectory event; it must never silently reinterpret StaleEpoch, SourceUnavailable,
// InterpolationFailure, or NonFiniteValue as physical magnetic shielding.
enum class FieldSampleStatus {
  Valid,
  InvalidRequest,
  OutsideDomain,
  StaleEpoch,
  SourceUnavailable,
  InterpolationFailure,
  NonFiniteValue
};

inline const char* CoordinateFrameName(CoordinateFrame frame) {
  switch (frame) {
    case CoordinateFrame::GSM: return "GSM";
    case CoordinateFrame::GEO: return "GEO";
    case CoordinateFrame::GSE: return "GSE";
    case CoordinateFrame::SM:  return "SM";
    default:                   return "UNKNOWN";
  }
}

inline const char* InterpolationModeName(InterpolationMode mode) {
  switch (mode) {
    case InterpolationMode::DirectAnalytic: return "DIRECT_ANALYTIC";
    case InterpolationMode::DirectEmpirical: return "DIRECT_EMPIRICAL";
    case InterpolationMode::CellCenteredLinear: return "CELL_CENTERED_LINEAR";
    case InterpolationMode::CellCenteredLinearDerivedElectric:
      return "CELL_CENTERED_LINEAR_DERIVED_E";
    default: return "UNKNOWN";
  }
}

inline const char* FieldSampleStatusName(FieldSampleStatus status) {
  switch (status) {
    case FieldSampleStatus::Valid:                return "VALID";
    case FieldSampleStatus::InvalidRequest:       return "INVALID_REQUEST";
    case FieldSampleStatus::OutsideDomain:        return "OUTSIDE_DOMAIN";
    case FieldSampleStatus::StaleEpoch:           return "STALE_EPOCH";
    case FieldSampleStatus::SourceUnavailable:    return "SOURCE_UNAVAILABLE";
    case FieldSampleStatus::InterpolationFailure: return "INTERPOLATION_FAILURE";
    case FieldSampleStatus::NonFiniteValue:       return "NONFINITE_VALUE";
    default:                                      return "UNKNOWN";
  }
}

struct AxisAlignedDomainSI {
  bool enabled{false};
  double minimum_m[3]{0.0,0.0,0.0};
  double maximum_m[3]{0.0,0.0,0.0};

  bool Contains(const double position_m[3]) const {
    if (!enabled) return true;
    if (position_m==nullptr) return false;
    for (int d=0;d<3;++d) {
      if (!std::isfinite(position_m[d]) ||
          position_m[d]<minimum_m[d] || position_m[d]>maximum_m[d]) return false;
    }
    return true;
  }
};

struct SnapshotMetadata {
  int schemaVersion{1};
  std::string snapshotId;
  std::string sourceId;
  std::string modelName;
  std::string epochUTC;
  CoordinateFrame frame{CoordinateFrame::Unknown};
  InterpolationMode interpolation{InterpolationMode::Unknown};
  AxisAlignedDomainSI domain;

  // These strings are intentionally fixed by the contract.  Keeping them in the
  // metadata makes serialized provenance self-describing and catches accidental
  // nT/km adapters during validation.
  std::string positionUnit{"m"};
  std::string magneticFieldUnit{"T"};
  std::string electricFieldUnit{"V/m"};

  bool magneticFieldAvailable{false};
  bool electricFieldAvailable{false};
  bool immutableDuringBatch{true};
  bool valid{false};
  std::string validityMessage;
};

struct FieldQuery {
  double position_m[3]{0.0,0.0,0.0};

  // Empty means "use the frozen snapshot epoch".  Supplying an epoch is a strong
  // synchronization check: a different value returns StaleEpoch instead of quietly
  // sampling the wrong field state.
  std::string epochUTC;
  bool requireElectricField{false};
  bool enforceDomain{true};
};

struct FieldSample {
  double magneticField_T[3]{0.0,0.0,0.0};
  double electricField_V_m[3]{0.0,0.0,0.0};
  FieldSampleStatus status{FieldSampleStatus::SourceUnavailable};
  InterpolationMode interpolation{InterpolationMode::Unknown};
  std::string snapshotId;
  std::string message;

  bool ok() const { return status==FieldSampleStatus::Valid; }
};

inline bool FiniteVector3(const double value[3]) {
  return value!=nullptr && std::isfinite(value[0]) &&
         std::isfinite(value[1]) && std::isfinite(value[2]);
}

// Deterministic FNV-1a identifier.  It is a provenance key, not a cryptographic data
// hash.  Providers must include every field-defining driver in the canonicalState
// string so two physically different snapshots cannot share an identifier.
inline std::string MakeSnapshotId(const std::string& sourceId,
                                  const std::string& epochUTC,
                                  const std::string& canonicalState) {
  const std::string input=sourceId+"\n"+epochUTC+"\n"+canonicalState;
  std::uint64_t hash=14695981039346656037ULL;
  for (std::string::const_iterator it=input.begin();it!=input.end();++it) {
    hash^=static_cast<unsigned char>(*it);
    hash*=1099511628211ULL;
  }
  std::ostringstream out;
  out << "field-v1-" << std::hex << std::setw(16) << std::setfill('0') << hash;
  return out.str();
}

inline FieldSampleStatus ValidateMetadata(const SnapshotMetadata& metadata,
                                          std::string* message=nullptr) {
  const auto fail=[&](FieldSampleStatus status,const char* text) {
    if (message!=nullptr) *message=text;
    return status;
  };

  if (!metadata.valid)
    return fail(FieldSampleStatus::SourceUnavailable,
                metadata.validityMessage.empty() ? "field snapshot is not valid" :
                metadata.validityMessage.c_str());
  if (metadata.snapshotId.empty() || metadata.sourceId.empty() ||
      metadata.modelName.empty() || metadata.epochUTC.empty())
    return fail(FieldSampleStatus::InvalidRequest,
                "field metadata is missing snapshot/source/model/epoch identity");
  if (metadata.frame==CoordinateFrame::Unknown)
    return fail(FieldSampleStatus::InvalidRequest,
                "field metadata has an unknown coordinate frame");
  if (metadata.positionUnit!="m" || metadata.magneticFieldUnit!="T" ||
      metadata.electricFieldUnit!="V/m")
    return fail(FieldSampleStatus::InvalidRequest,
                "field metadata does not satisfy the SI unit contract");
  if (!metadata.magneticFieldAvailable)
    return fail(FieldSampleStatus::SourceUnavailable,
                "magnetic field is unavailable in this snapshot");
  if (!metadata.immutableDuringBatch)
    return fail(FieldSampleStatus::InvalidRequest,
                "field snapshot is not immutable for the trajectory batch");
  if (metadata.domain.enabled) {
    for (int d=0;d<3;++d) {
      if (!std::isfinite(metadata.domain.minimum_m[d]) ||
          !std::isfinite(metadata.domain.maximum_m[d]) ||
          metadata.domain.minimum_m[d]>metadata.domain.maximum_m[d])
        return fail(FieldSampleStatus::InvalidRequest,
                    "field metadata has an invalid domain bound");
    }
  }

  if (message!=nullptr) message->clear();
  return FieldSampleStatus::Valid;
}

inline FieldSampleStatus ValidateQuery(const SnapshotMetadata& metadata,
                                       const FieldQuery& query,
                                       std::string* message=nullptr) {
  FieldSampleStatus status=ValidateMetadata(metadata,message);
  if (status!=FieldSampleStatus::Valid) return status;
  if (!FiniteVector3(query.position_m)) {
    if (message!=nullptr) *message="field query position is non-finite";
    return FieldSampleStatus::InvalidRequest;
  }
  if (!query.epochUTC.empty() && query.epochUTC!=metadata.epochUTC) {
    if (message!=nullptr) *message="field query epoch differs from frozen snapshot epoch";
    return FieldSampleStatus::StaleEpoch;
  }
  if (query.requireElectricField && !metadata.electricFieldAvailable) {
    if (message!=nullptr) *message="electric field was requested but is unavailable";
    return FieldSampleStatus::SourceUnavailable;
  }
  if (query.enforceDomain && !metadata.domain.Contains(query.position_m)) {
    if (message!=nullptr) *message="field query position is outside snapshot validity domain";
    return FieldSampleStatus::OutsideDomain;
  }
  if (message!=nullptr) message->clear();
  return FieldSampleStatus::Valid;
}

class IFieldSnapshot {
public:
  virtual ~IFieldSnapshot() = default;
  virtual const SnapshotMetadata& Metadata() const = 0;
  virtual FieldSample Sample(const FieldQuery& query) const = 0;
};

struct SnapshotRequest {
  std::string epochUTC;
  std::string requestId;
};

class IFieldProvider {
public:
  virtual ~IFieldProvider() = default;
  virtual std::string SourceId() const = 0;
  virtual std::shared_ptr<const IFieldSnapshot>
      CreateSnapshot(const SnapshotRequest& request) = 0;
};

inline void RequireSameSnapshot(const SnapshotMetadata& expected,
                                const SnapshotMetadata& actual,
                                const char* context) {
  if (expected.snapshotId.empty() || actual.snapshotId.empty() ||
      expected.snapshotId!=actual.snapshotId) {
    std::ostringstream msg;
    msg << (context!=nullptr ? context : "field product")
        << ": cutoff/access and flux/spectrum did not use the same field snapshot"
        << " (expected='" << expected.snapshotId
        << "', actual='" << actual.snapshotId << "').";
    throw std::runtime_error(msg.str());
  }
}

} // namespace Field
} // namespace Earth

#endif // _SRC_EARTH_UTIL_FIELD_PROVIDER_H_
