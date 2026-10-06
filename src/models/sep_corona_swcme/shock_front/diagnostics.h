#ifndef SEP_CORONA_SWCME_SHOCK_FRONT_DIAGNOSTICS_H
#define SEP_CORONA_SWCME_SHOCK_FRONT_DIAGNOSTICS_H

#include "provider.h"

#include <cstdint>
#include <memory>
#include <string>
#include <vector>

namespace SEP { namespace CoronaSwcme { namespace ShockFront {

struct FrontIntersection {
  std::uint64_t branchId = 0;
  std::size_t segmentIndex = 0;
  double segmentFraction = 0.0;
  double arcLengthM = 0.0;
  CoronalCME::Vec3 positionM;
  bool tangent = false;
  int signedPolarityAlongTrace = 0;
  ShockRecord front;
};

// Intersect a caller-supplied instantaneous field-line polyline with the
// analytical finite leading SSE surface.  Segment/sphere roots are solved
// analytically, including a double tangent root; rear-sphere roots and roots
// outside finite support are discarded by the production geometry query.
Core::Result<std::vector<FrontIntersection>> IntersectPolylineWithFront(
    const Provider& provider,const std::vector<CoronalCME::Vec3>& pointsM,
    double epochS,double distanceToleranceM);

enum class ObserverPassageKind { Transverse, Grazing };
struct ObserverPassage {
  double timeS = 0.0;
  CoronalCME::Vec3 positionM;
  ObserverPassageKind kind = ObserverPassageKind::Transverse;
  double relativeNormalSpeedMPerS = 0.0;
  ShockRecord front;
};

// Find all passages of a constant-velocity HCI observer over a closed time
// interval.  Sign-changing roots are bracketed; local minima of |g| are also
// minimized so a geometric graze is not lost merely because g keeps one sign.
Core::Result<std::vector<ObserverPassage>> FindObserverPassages(
    const Provider& provider,CoronalCME::Vec3 positionAtReferenceM,
    CoronalCME::Vec3 velocityMPerS,double referenceTimeS,double beginTimeS,
    double endTimeS,double scanStepS,double timeToleranceS,
    double distanceToleranceM);

struct ReducedRestartState {
  std::string eventIdentity;
  double epochS = 0.0;
  std::uint64_t generation = 0;
  Phase phase = Phase::CoronalHistory;
  double handoffTimeS = 0.0;
  double handoffRadiusM = 0.0;
  double handoffSpeedMPerS = 0.0;
};

Core::Result<ReducedRestartState> MakeRestartState(const Provider& provider);
Core::Result<std::shared_ptr<const Epoch>> RestoreRestartState(
    Provider* provider,const ReducedRestartState& state);

// Machine-readable text helpers do no filesystem I/O.  The application owns
// output paths and collective file policy; shared physics owns names, units,
// validity and the rule that an absent downstream is written as JSON null or
// an empty CSV field rather than a physical-looking zero.
std::string SerializeEpochJson(const Epoch& epoch);
std::string SerializeSurfaceCsv(const Epoch& epoch);

struct ConnectedBranchSample {
  std::uint64_t branchId = 0;
  double timeS = 0.0;
  CoronalCME::Vec3 positionM;
  double obliquityRad = 0.0;
  double fastMach = 0.0;
  bool valid = false;
};
struct ConnectedBranchDerivative {
  bool valid = false;
  CoronalCME::Vec3 intersectionVelocityMPerS;
  double obliquityRateRadPerS = 0.0;
  double fastMachRatePerS = 0.0;
  std::string reason;
};
ConnectedBranchDerivative DifferentiateConnectedBranch(
    const ConnectedBranchSample& previous,const ConnectedBranchSample& current,
    const ConnectedBranchSample& next);

enum class EnsembleOutcome { KnownAcceptedShock, KnownNoShock, NumericalUnknown };
struct EnsembleMember {
  std::string commonLabel;
  double weight = 0.0;
  EnsembleOutcome outcome = EnsembleOutcome::NumericalUnknown;
};
struct EnsembleLedger {
  double acceptedWeight = 0.0;
  double noShockWeight = 0.0;
  double unknownWeight = 0.0;
  double acceptedProbabilityLower = 0.0;
  double acceptedProbabilityUpper = 0.0;
};
Core::Result<EnsembleLedger> SummarizeEnsemble(
    const std::vector<EnsembleMember>& members,double normalizationTolerance);

const char* Name(ObserverPassageKind) noexcept;
const char* Name(EnsembleOutcome) noexcept;

} } } // namespace SEP::CoronaSwcme::ShockFront

#endif
