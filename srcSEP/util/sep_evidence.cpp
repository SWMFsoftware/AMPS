#include "sep_evidence.h"

#include <sstream>

namespace SEP {
namespace Evidence {

const char* LevelName(Level level) {
  switch (level) {
    case Level::AnalyticalCore: return "analytical-core";
    case Level::SourceIntegration: return "source-integration";
    case Level::NativeAmps: return "native-amps";
    case Level::SwmfReplay: return "swmf-replay";
    case Level::ObservationalValidation: return "observational-validation";
  }
  return "unknown";
}

const char* GateStatusName(GateStatus status) {
  switch (status) {
    case GateStatus::Pass: return "PASS";
    case GateStatus::Fail: return "FAIL";
    case GateStatus::Blocked: return "BLOCKED";
    case GateStatus::Incomplete: return "INCOMPLETE";
  }
  return "ERROR";
}

Transport::Status ValidateClaim(const Claim& claim) {
  if (claim.id.empty() || claim.statement.empty())
    return Transport::Status::Error(Transport::StatusCode::InvalidArgument,
                                    "evidence claim requires an ID and statement");
  if (claim.status == GateStatus::Pass &&
      static_cast<int>(claim.observedLevel) < static_cast<int>(claim.requiredLevel))
    return Transport::Status::Error(Transport::StatusCode::UnsupportedConfiguration,
        "claim cannot PASS below its required evidence level");
  if (claim.status == GateStatus::Pass && claim.artifacts.empty())
    return Transport::Status::Error(Transport::StatusCode::InvalidArgument,
                                    "passing claim requires a reproducible artifact");
  for (std::size_t i = 0; i < claim.artifacts.size(); ++i)
    if (claim.artifacts[i].command.empty() || claim.artifacts[i].path.empty())
      return Transport::Status::Error(Transport::StatusCode::InvalidArgument,
                                      "evidence artifact lacks command or path");
  return Transport::Status::Ok();
}

std::string RenderClaimTable(const std::vector<Claim>& claims) {
  std::ostringstream out;
  out << "| Claim | Required evidence | Observed evidence | Status |\n"
      << "| --- | --- | --- | --- |\n";
  for (std::size_t i = 0; i < claims.size(); ++i)
    out << "| " << claims[i].id << " | " << LevelName(claims[i].requiredLevel)
        << " | " << LevelName(claims[i].observedLevel) << " | "
        << GateStatusName(claims[i].status) << " |\n";
  return out.str();
}

Transport::Status ValidateNativeObservation(const NativeHarnessRequest& request,
                                             const NativeObservation& observed,
                                             Level executingLevel) {
  if (executingLevel != Level::NativeAmps)
    return Transport::Status::Error(Transport::StatusCode::UnsupportedConfiguration,
        "native-adapter evidence requires execution by the linked AMPS binary");
  const Transport::Status compatibility =
      ConfigurationMatrix::Preflight(request.combination, false);
  if (!compatibility.ok()) return compatibility;
  if (request.expectedBinary.empty() || request.configurationFingerprint.empty() ||
      request.sourceGeneration.empty() || request.threadCount <= 0 ||
      request.rankCount <= 0 || request.boundaryCase.empty())
    return Transport::Status::Error(Transport::StatusCode::InvalidArgument,
                                    "native harness request is incomplete");
  if (!observed.productionMoverEntered || !observed.coefficientAdapterEntered ||
      !observed.turbulenceDriverEntered)
    return Transport::Status::Error(Transport::StatusCode::OutOfDomain,
        "native harness did not traverse every required production adapter");
  if (observed.queueFlushOwners != 1)
    return Transport::Status::Error(Transport::StatusCode::OutOfDomain,
        "exactly one production phase must own coupling-queue flush");
  if (observed.executedBinary != request.expectedBinary ||
      observed.configurationFingerprint != request.configurationFingerprint ||
      observed.sourceGeneration != request.sourceGeneration ||
      observed.threadCount != request.threadCount ||
      observed.rankCount != request.rankCount)
    return Transport::Status::Error(Transport::StatusCode::OutOfDomain,
        "native observation identity or decomposition differs from request");
  return Transport::Status::Ok();
}

NativeHarnessResult RunNativeHarness(const NativeHarnessRequest& request,
                                     NativeAdapter* adapter,
                                     Level executingLevel) {
  NativeHarnessResult result;
  result.evidenceLevel = executingLevel;
  if (!adapter) {
    result.status = Transport::Status::Error(
        Transport::StatusCode::InvalidArgument,
        "native harness requires an adapter implementation");
    return result;
  }
  result.status = adapter->Execute(request, &result.observation);
  if (!result.status.ok()) return result;
  result.status = ValidateNativeObservation(request, result.observation,
                                            executingLevel);
  return result;
}

}  // namespace Evidence
}  // namespace SEP
