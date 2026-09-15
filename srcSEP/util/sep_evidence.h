#ifndef SEP_UTIL_SEP_EVIDENCE_H
#define SEP_UTIL_SEP_EVIDENCE_H

#include "sep_configuration_matrix.h"
#include "sep_transport_common.h"

#include <cstdint>
#include <string>
#include <vector>

namespace SEP {
namespace Evidence {

// Evidence levels are ordered by what actually executed, not by how realistic
// a test description sounds.  SourceIntegration deliberately sits below a
// linked native shell, while SWMF replay and held-out observations remain
// independent external gates rather than automatic promotions of native tests.
enum class Level {
  AnalyticalCore,
  SourceIntegration,
  NativeAmps,
  SwmfReplay,
  ObservationalValidation
};
enum class GateStatus { Pass, Fail, Blocked, Incomplete };

const char* LevelName(Level level);
const char* GateStatusName(GateStatus status);

struct Artifact {
  std::string command;
  std::string path;
  std::string checksum;
};

struct Claim {
  std::string id;
  std::string statement;
  Level requiredLevel = Level::AnalyticalCore;
  Level observedLevel = Level::AnalyticalCore;
  GateStatus status = GateStatus::Incomplete;
  std::vector<Artifact> artifacts;
};

Transport::Status ValidateClaim(const Claim& claim);
std::string RenderClaimTable(const std::vector<Claim>& claims);

struct NativeHarnessRequest {
  ConfigurationMatrix::Combination combination;
  std::string expectedBinary;
  std::string configurationFingerprint;
  std::string sourceGeneration;
  int threadCount = 1;
  int rankCount = 1;
  std::string boundaryCase;
};

struct NativeObservation {
  bool productionMoverEntered = false;
  bool coefficientAdapterEntered = false;
  bool turbulenceDriverEntered = false;
  std::uint64_t queueFlushOwners = 0;
  std::string executedBinary;
  std::string configurationFingerprint;
  std::string sourceGeneration;
  int threadCount = 0;
  int rankCount = 0;
};

class NativeAdapter {
 public:
  virtual ~NativeAdapter() {}
  // Implementations must create the smallest native field line and particle
  // through host APIs, advance the requested registered mover, and fill only
  // observations actually emitted by production hooks.
  virtual Transport::Status Execute(const NativeHarnessRequest& request,
                                    NativeObservation* observation) = 0;
};

struct NativeHarnessResult {
  Transport::Status status;
  NativeObservation observation;
  Level evidenceLevel = Level::SourceIntegration;
};

NativeHarnessResult RunNativeHarness(const NativeHarnessRequest& request,
                                     NativeAdapter* adapter,
                                     Level executingLevel);

// ValidateNativeObservation is the WP34 trust boundary.  A dependency-light
// double can test orchestration but cannot set NativeAmps evidence; promotion
// requires observations emitted by the linked production binary and exact
// identity/layout agreement with the request.
Transport::Status ValidateNativeObservation(const NativeHarnessRequest& request,
                                             const NativeObservation& observed,
                                             Level executingLevel);

}  // namespace Evidence
}  // namespace SEP

#endif  // SEP_UTIL_SEP_EVIDENCE_H
