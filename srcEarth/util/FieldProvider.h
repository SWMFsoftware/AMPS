#ifndef _SRC_EARTH_UTIL_FIELD_PROVIDER_H_
#define _SRC_EARTH_UTIL_FIELD_PROVIDER_H_

//======================================================================================
// FieldProvider.h
//======================================================================================
// Backend-neutral, dependency-free contract for magnetic/electric fields used by
// cutoff, access, flux, and energy-spectrum calculations.
//
// Step-3 design rule
// ------------------
// A provider may own mutable setup state, but CreateSnapshot() must return a frozen,
// read-only field state.  One trajectory batch samples exactly one snapshot.  Direct
// analytic/empirical fields and compact Mode3D/SWMF fields expose the same metadata,
// status, unit, and identity contract without changing the particle mover interface.
//
// This header intentionally has no MPI, PIC, SPICE, Geopack, or SWMF dependency.  The
// production adapters live beside their backends; the contract and its failure modes
// can therefore be validated by a strict standalone C++11 unit test.
//
// Physical contract
// -----------------
//   position : metres
//   B        : tesla
//   E        : volt/metre
//   frame    : explicit (GSM in current production Earth paths)
//
// A field error is not a geomagnetic-access decision.  OUTSIDE_DOMAIN may be consumed
// by a boundary locator, but stale state, unavailable sources, interpolation failures,
// and non-finite values must never be silently reclassified as forbidden access.
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

  // Fixed strings make serialized provenance self-describing and allow validation to
  // reject a backend that accidentally exposes nT, km, Re, or an implicit E unit.
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

  // Empty means "use the snapshot epoch".  A nonempty value is a synchronization
  // assertion; it never asks an immutable snapshot to advance to another epoch.
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

// Stable FNV-1a provenance identifier.  It is not a cryptographic checksum.  The
// canonicalState argument must contain every physical driver and discretization
// choice that can change sampled values.  Output names and request labels must not be
// included: two requests for the same frozen physical field must receive the same ID.
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
  const auto fail=[&](FieldSampleStatus status,const std::string& text) {
    if (message!=nullptr) *message=text;
    return status;
  };

  if (!metadata.valid)
    return fail(FieldSampleStatus::SourceUnavailable,
                metadata.validityMessage.empty() ?
                "field snapshot is not valid" : metadata.validityMessage);
  if (metadata.schemaVersion!=1)
    return fail(FieldSampleStatus::InvalidRequest,
                "unsupported field-snapshot metadata schema");
  if (metadata.snapshotId.empty() || metadata.sourceId.empty() ||
      metadata.modelName.empty() || metadata.epochUTC.empty())
    return fail(FieldSampleStatus::InvalidRequest,
                "field metadata is missing snapshot/source/model/epoch identity");
  if (metadata.frame==CoordinateFrame::Unknown)
    return fail(FieldSampleStatus::InvalidRequest,
                "field metadata has an unknown coordinate frame");
  if (metadata.interpolation==InterpolationMode::Unknown)
    return fail(FieldSampleStatus::InvalidRequest,
                "field metadata has an unknown interpolation mode");
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
  // Diagnostic label only.  A provider may log it, but it must not make physically
  // identical snapshots acquire different snapshot IDs.
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
  std::string expectedError,actualError;
  const FieldSampleStatus expectedStatus=ValidateMetadata(expected,&expectedError);
  const FieldSampleStatus actualStatus=ValidateMetadata(actual,&actualError);
  bool sameDomain=expected.domain.enabled==actual.domain.enabled;
  if (sameDomain && expected.domain.enabled) {
    for (int d=0;d<3;++d) {
      sameDomain=sameDomain &&
          expected.domain.minimum_m[d]==actual.domain.minimum_m[d] &&
          expected.domain.maximum_m[d]==actual.domain.maximum_m[d];
    }
  }

  // The ID is the primary physical-state key, but compare the complete public
  // contract as a second line of defence.  If an adapter accidentally reuses an ID
  // after changing frame, interpolation, units, capabilities, or domain, cutoff and
  // flux must fail synchronization instead of silently combining unlike fields.
  const bool sameContract=
      expected.schemaVersion==actual.schemaVersion &&
      expected.sourceId==actual.sourceId &&
      expected.modelName==actual.modelName &&
      expected.epochUTC==actual.epochUTC &&
      expected.frame==actual.frame &&
      expected.interpolation==actual.interpolation &&
      expected.positionUnit==actual.positionUnit &&
      expected.magneticFieldUnit==actual.magneticFieldUnit &&
      expected.electricFieldUnit==actual.electricFieldUnit &&
      expected.magneticFieldAvailable==actual.magneticFieldAvailable &&
      expected.electricFieldAvailable==actual.electricFieldAvailable &&
      expected.immutableDuringBatch==actual.immutableDuringBatch &&
      sameDomain;
  if (expectedStatus!=FieldSampleStatus::Valid ||
      actualStatus!=FieldSampleStatus::Valid ||
      expected.snapshotId!=actual.snapshotId ||
      !sameContract) {
    std::ostringstream msg;
    msg << (context!=nullptr ? context : "field product")
        << ": field snapshot changed or is invalid"
        << " (expectedId='" << expected.snapshotId
        << "', actualId='" << actual.snapshotId
        << "', expectedStatus=" << FieldSampleStatusName(expectedStatus)
        << ", actualStatus=" << FieldSampleStatusName(actualStatus)
        << ", sameContract=" << (sameContract ? "true" : "false") << ").";
    throw std::runtime_error(msg.str());
  }
}

} // namespace Field
} // namespace Earth

#endif // _SRC_EARTH_UTIL_FIELD_PROVIDER_H_
