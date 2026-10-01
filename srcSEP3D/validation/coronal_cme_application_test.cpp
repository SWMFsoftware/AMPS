#include "coronal_cme_application_test.h"

#include "sep_coronal_cme/particle_source.h"
#include "sep_coronal_cme/runtime_integration.h"

#include <algorithm>
#include <cctype>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <map>
#include <sstream>

namespace SEP3D { namespace Validation { namespace {

namespace fs = std::filesystem;
namespace SCCM = SEP::CoronalCME;

std::string Upper(std::string value) {
  std::transform(value.begin(), value.end(), value.begin(),
      [](unsigned char c) { return static_cast<char>(std::toupper(c)); });
  return value;
}

NativeTestResult Result(const NativeTestDescriptor& descriptor,
                        NativeTestStatus status, std::string message) {
  NativeTestResult result;
  result.id = descriptor.id;
  result.name = descriptor.name;
  result.status = status;
  result.message = std::move(message);
  return result;
}

std::string JsonString(const std::string& input) {
  std::ostringstream output;
  output << '"';
  for (unsigned char c : input) {
    switch (c) {
      case '"': output << "\\\""; break;
      case '\\': output << "\\\\"; break;
      case '\n': output << "\\n"; break;
      case '\r': output << "\\r"; break;
      case '\t': output << "\\t"; break;
      default:
        if (c < 0x20) {
          output << "\\u" << std::hex << std::setw(4) << std::setfill('0')
                 << static_cast<unsigned>(c) << std::dec;
        } else {
          output << static_cast<char>(c);
        }
    }
  }
  output << '"';
  return output.str();
}

NativeTestResult EvaluateOne(const NativeTestDescriptor& descriptor,
                             const NativeApplicationState& state) {
  const std::string& id = descriptor.id;
  if (id == "NAT3D01") {
    NativeTestResult result = Result(descriptor,
        state.globalAllocatedBlocks > 0 && state.globalPhysicalCells > 0 &&
                state.solarBoundaryRegistered
            ? NativeTestStatus::Pass : NativeTestStatus::Fail,
        "real AMPS mesh allocation and solar-boundary registration");
    result.metrics = {{"mpi_ranks", state.mpiRankCount},
                      {"allocated_blocks", state.globalAllocatedBlocks},
                      {"physical_cells", state.globalPhysicalCells}};
    return result;
  }
  if (id == "NAT3D02")
    return Result(descriptor,
        state.backgroundReady && state.turbulenceReady &&
                state.finiteBackgroundAndTurbulence
            ? NativeTestStatus::Pass : NativeTestStatus::Fail,
        "background and directional turbulence are finite on physical cells");
  if (id == "NAT3D03")
    return Result(descriptor, state.finiteBackgroundDerivatives
        ? NativeTestStatus::Pass : NativeTestStatus::Fail,
        "AMR background derivatives and focusing state are finite");
  if (id == "NAT3D09") {
    NativeTestResult result = Result(descriptor,
        state.globalAllocatedBlocks > 0 ? NativeTestStatus::Pass
                                        : NativeTestStatus::Fail,
        "distributed AMPS block ownership is nonempty and globally reduced");
    result.metrics = {{"allocated_blocks", state.globalAllocatedBlocks},
                      {"completed_steps", state.completedSteps}};
    return result;
  }
  if (id == "NAT3D10")
    return Result(descriptor,
        state.plannedActiveLeaves > 0 && state.activeMaskInstalled
            ? NativeTestStatus::Pass : NativeTestStatus::Fail,
        "preflight resource plan and active-use mask reached allocation");
  if (id == "NAT3D11") {
    if (!state.shockRequired)
      return Result(descriptor, NativeTestStatus::Skip,
                    "input selects no live shock authority");
    return Result(descriptor,
        state.shockReady && state.shockBackgroundGeneration ==
                                state.backgroundGeneration
            ? NativeTestStatus::Pass : NativeTestStatus::Fail,
        "shock and background generations are coherent at a joined boundary");
  }
  if (id == "NAT3D12")
    return Result(descriptor,
        state.initializationProductsExist && state.initializationProductsFinite
            ? NativeTestStatus::Pass : NativeTestStatus::Fail,
        "initialization products exist and contain no non-finite numeric token");
  if (id == "MPI3D01") {
    if (state.mpiRankCount < 2)
      return Result(descriptor, NativeTestStatus::Skip,
                    "requires at least two MPI ranks");
    return Result(descriptor, state.mpiFingerprintConsistent
        ? NativeTestStatus::Pass : NativeTestStatus::Fail,
        "all ranks report the same configuration/background identity");
  }
  if (id == "MPI3D02") {
    if (!state.restartConfigured)
      return Result(descriptor, NativeTestStatus::Skip,
                    "requires a validated --restart checkpoint");
    return Result(descriptor, state.completedSteps > 0 &&
                                  state.mpiFingerprintConsistent
        ? NativeTestStatus::Pass : NativeTestStatus::Fail,
        "restarted MPI state advanced with a consistent identity");
  }

  if (id == "SCCM3D01") {
    SCCM::InitializationLedger ledger;
    ledger.completedMask = state.initializationMask;
    ledger.backgroundGeneration = state.backgroundGeneration;
    ledger.shockBackgroundGeneration = state.shockBackgroundGeneration;
    const SEP::Core::Status checked =
        SCCM::ValidateInitializationForOutput(ledger);
    return Result(descriptor, checked.ok() ? NativeTestStatus::Pass
                                           : NativeTestStatus::Fail,
                  checked.ok() ? "all ten SCCM initialization conditions close"
                               : checked.message);
  }
  if (id == "SCCM3D02") {
    std::vector<SCCM::SpeciesWeight> weights;
    for (const auto& species : state.species)
      weights.push_back({species.chemicalSymbol + "#" +
                             std::to_string(species.compiledSlot),
                         species.compiledSlot, species.particleWeight});
    const SEP::Core::Status checked = SCCM::ValidateAllSpeciesWeights(
        weights, static_cast<int>(state.species.size()));
    bool positiveSteps = true;
    for (const auto& species : state.species)
      positiveSteps = positiveSteps && species.timeStepS > 0.0;
    return Result(descriptor, checked.ok() && positiveSteps
        ? NativeTestStatus::Pass : NativeTestStatus::Fail,
        checked.ok() && positiveSteps
            ? "every compiled AMPS species has positive weight and time step"
            : (checked.ok() ? "a compiled species has no positive time step"
                            : checked.message));
  }
  if (id == "SCCM3D03") {
    std::vector<SCCM::CompiledSpeciesIdentity> species;
    for (const auto& item : state.species) {
      SCCM::CompiledSpeciesIdentity identity;
      identity.stableId = item.chemicalSymbol + "#" +
                          std::to_string(item.compiledSlot);
      identity.compiledSlot = item.compiledSlot;
      identity.chemicalSymbol = item.chemicalSymbol;
      identity.massKg = item.massKg;
      identity.chargeC = item.chargeC;
      // Zero means not applicable/unknown (for example an electron).  The
      // source binding validates compiled identity, not an energy-per-nucleon
      // observer conversion, and therefore must not invent A=1 for leptons.
      identity.nucleonCount = 0;
      identity.sourceEnabled = state.sourceEnabled;
      species.push_back(identity);
    }
    const SEP::Core::Status checked = SCCM::ValidateSourceSpecies(
        species, static_cast<int>(state.species.size()));
    return Result(descriptor, checked.ok() ? NativeTestStatus::Pass
                                           : NativeTestStatus::Fail,
                  checked.ok() ? "all compiled species bind to the SCCM source"
                               : checked.message);
  }
  if (id == "SCCM3D04") {
    // Both the identity full-domain plan and a pruned corridor are valid.
    // Require installation/allocation evidence independently of whether any
    // leaves were removed; a wide corridor may legitimately retain all leaves.
    std::string failure;
    if (!state.solarBoundaryRegistered)
      failure = "solar absorption boundary is not registered";
    else if (!state.activeMaskInstalled)
      failure = "active-region plan was not installed and verified";
    else if (!state.activeRegionAllocationVerified)
      failure = "active-region block allocation was not verified";
    else if (state.activeRegionMode != "full-domain" &&
             state.activeRegionMode != "parker-tube")
      failure = "active-region mode is missing or unsupported";
    else if (state.plannedActiveLeaves == 0)
      failure = "active-region plan contains no physical active leaves";
    else if (state.globalAllocatedBlocks != state.plannedActiveLeaves)
      failure = "allocated owner-block count differs from planned active leaves";
    else if (state.activeRegionPruningApplied !=
             (state.plannedInactiveLeaves > 0))
      failure = "pruning flag disagrees with the planned inactive-leaf count";
    else if (state.plannedSolarInteriorLeaves > state.plannedInactiveLeaves)
      failure = "solar-interior leaves are not all excluded by the mask";
    else if (state.activeRegionMode == "full-domain" &&
             state.plannedInactiveLeaves != state.plannedSolarInteriorLeaves)
      failure = "full-domain mode removed leaves outside the solid Sun";

    NativeTestResult result = Result(descriptor,
        failure.empty() ? NativeTestStatus::Pass : NativeTestStatus::Fail,
        failure.empty()
            ? "solar boundary and " + state.activeRegionMode +
                  " plan/allocation are verified (pruning=" +
                  (state.activeRegionPruningApplied ? "yes)" : "no)")
            : failure);
    result.metrics = {
        {"plan_installed", state.activeMaskInstalled ? 1.0 : 0.0},
        {"pruning_applied", state.activeRegionPruningApplied ? 1.0 : 0.0},
        {"allocation_verified", state.activeRegionAllocationVerified ? 1.0 : 0.0},
        {"planned_active_leaves", state.plannedActiveLeaves},
        {"planned_inactive_leaves", state.plannedInactiveLeaves},
        {"planned_solar_interior_leaves", state.plannedSolarInteriorLeaves},
        {"allocated_blocks", state.globalAllocatedBlocks}};
    return result;
  }
  if (id == "SCCM3D05")
    return Result(descriptor,
        state.backgroundGeneration > 0 &&
                state.shockBackgroundGeneration == state.backgroundGeneration
            ? NativeTestStatus::Pass : NativeTestStatus::Fail,
        "SCCM background/shock generation identity is current");
  if (id == "SCCM3D06")
    return Result(descriptor,
        state.initializationProductsExist && state.initializationProductsFinite
            ? NativeTestStatus::Pass : NativeTestStatus::Fail,
        "SCCM initialization output is finite and physically populated");
  if (id == "SCCM3D07") {
    if (state.expectedMpiRanks > 0 &&
        state.expectedMpiRanks != state.mpiRankCount)
      return Result(descriptor, NativeTestStatus::Fail,
                    "MPI rank count differs from --expect-mpi-ranks");
    return Result(descriptor, state.mpiFingerprintConsistent
        ? NativeTestStatus::Pass : NativeTestStatus::Fail,
        "SCCM identity and collective state agree on every MPI rank");
  }
  return Result(descriptor, NativeTestStatus::Error,
                "native test has no evaluation callback");
}

}  // namespace

const std::vector<NativeTestDescriptor>& CoronalCmeNativeTests() {
  static const std::vector<NativeTestDescriptor> tests = {
      {"NAT3D01", "AMPS mesh integration", "real AMR allocation and solar sphere"},
      {"NAT3D02", "Cell background integration", "finite plasma/IMF/turbulence"},
      {"NAT3D03", "AMR gradient integration", "finite derivatives and focusing"},
      {"NAT3D09", "Particle load balance", "collective block/runtime evidence"},
      {"NAT3D10", "Production resource budgets", "preflight and active mask"},
      {"NAT3D11", "Coupled snapshot/shock schedule", "joined generations"},
      {"NAT3D12", "Production product grammar", "finite initialization products"},
      {"MPI3D01", "Multi-rank state reproducibility", "rank identity agreement"},
      {"MPI3D02", "Multi-rank restart continuation", "restart advancement"},
      {"SCCM3D01", "SCCM initialization ledger", "all ten initialization gates", "sep-corona"},
      {"SCCM3D02", "SCCM species numerics", "all-species weights and steps", "sep-corona"},
      {"SCCM3D03", "SCCM source species binding", "compiled AMPS table coverage", "sep-corona"},
      {"SCCM3D04", "SCCM mesh boundary", "sphere and active-region contract", "sep-corona"},
      {"SCCM3D05", "SCCM provider generations", "background/shock coherence", "sep-corona"},
      {"SCCM3D06", "SCCM initialization output", "finite populated products", "sep-corona"},
      {"SCCM3D07", "SCCM MPI identity", "collective fingerprint agreement", "sep-corona"},
  };
  return tests;
}

Core::Status SelectCoronalCmeNativeTests(
    bool allTests, const std::vector<std::string>& requested,
    std::vector<NativeTestDescriptor>* selected,
    const std::string& suite) {
  if (selected == nullptr)
    return Core::Status(Core::StatusCode::InvalidInput,
                        "native test selection output is null");
  if (!suite.empty() && (allTests || !requested.empty()))
    return Core::Status(Core::StatusCode::InvalidInput,
                        "native suite and individual/all selectors conflict");
  if (!suite.empty() && suite != "sep-corona")
    return Core::Status(Core::StatusCode::InvalidInput,
                        "unknown native test suite '" + suite + "'");
  std::map<std::string, NativeTestDescriptor> known;
  for (const auto& descriptor : CoronalCmeNativeTests())
    known.emplace(descriptor.id, descriptor);
  std::vector<NativeTestDescriptor> candidate;
  if (allTests) {
    candidate = CoronalCmeNativeTests();
  } else if (!suite.empty()) {
    // Iterate the authoritative registry on every invocation. No fixed ID
    // range, count, Python catalogue or input-deck list controls membership.
    for (const auto& descriptor : CoronalCmeNativeTests())
      if (descriptor.suite == suite) candidate.push_back(descriptor);
  } else {
    for (const std::string& raw : requested) {
      const auto found = known.find(Upper(raw));
      if (found == known.end())
        return Core::Status(Core::StatusCode::InvalidInput,
                            "unknown native test ID '" + raw + "'");
      if (std::none_of(candidate.begin(), candidate.end(),
              [&](const NativeTestDescriptor& item) {
                return item.id == found->first;
              })) candidate.push_back(found->second);
    }
  }
  if (candidate.empty())
    return Core::Status(Core::StatusCode::InvalidInput,
                        "native test selection is empty");
  *selected = std::move(candidate);
  return Core::Status::OK();
}

std::vector<NativeTestResult> EvaluateCoronalCmeNativeTests(
    const NativeApplicationState& state,
    const std::vector<NativeTestDescriptor>& selected) {
  std::vector<NativeTestResult> results;
  results.reserve(selected.size());
  for (const auto& descriptor : selected)
    results.push_back(EvaluateOne(descriptor, state));
  return results;
}

const char* Name(NativeTestStatus status) noexcept {
  switch (status) {
    case NativeTestStatus::Pass: return "PASS";
    case NativeTestStatus::Fail: return "FAIL";
    case NativeTestStatus::Skip: return "SKIP";
    case NativeTestStatus::Error: return "ERROR";
  }
  return "ERROR";
}

Core::Status WriteNativeTestJson(
    const std::string& path, const NativeApplicationState& state,
    const std::vector<NativeTestResult>& results) {
  if (path.empty())
    return Core::Status(Core::StatusCode::InvalidInput,
                        "native JSON path is empty");
  const fs::path target(path);
  std::error_code error;
  if (!target.parent_path().empty())
    fs::create_directories(target.parent_path(), error);
  if (error)
    return Core::Status(Core::StatusCode::Error,
                        "cannot create native-test output directory: " +
                            error.message());
  const fs::path temporary = target.string() + ".tmp";
  std::ofstream out(temporary);
  if (!out)
    return Core::Status(Core::StatusCode::Error,
                        "cannot open native-test temporary JSON");
  out << "{\n  \"schema\": \"srcsep-component-tests-v1\",\n"
      << "  \"configuration_fingerprint\": "
      << JsonString(state.configurationFingerprint) << ",\n"
      << "  \"mpi_ranks\": " << state.mpiRankCount << ",\n"
      << "  \"completed_steps\": " << state.completedSteps << ",\n"
      << "  \"providers\": {\"input_schema_version\": " << state.inputSchemaVersion
      << ", \"background\": " << JsonString(state.backgroundAuthority)
      << ", \"shock\": " << JsonString(state.shockAuthority) << "},\n"
      // Additive evidence preserves the existing report schema/runner ABI.
      << "  \"active_region\": {\"mode\": "
      << JsonString(state.activeRegionMode)
      << ", \"plan_installed\": " << (state.activeMaskInstalled ? "true" : "false")
      << ", \"pruning_applied\": "
      << (state.activeRegionPruningApplied ? "true" : "false")
      << ", \"allocation_verified\": "
      << (state.activeRegionAllocationVerified ? "true" : "false")
      << ", \"planned_active_leaves\": " << state.plannedActiveLeaves
      << ", \"planned_inactive_leaves\": " << state.plannedInactiveLeaves
      << ", \"planned_solar_interior_leaves\": " << state.plannedSolarInteriorLeaves
      << ", \"allocated_blocks\": " << state.globalAllocatedBlocks << "},\n"
      << "  \"results\": [\n";
  for (std::size_t i = 0; i < results.size(); ++i) {
    const auto& result = results[i];
    out << "    {\"id\": " << JsonString(result.id)
        << ", \"name\": " << JsonString(result.name)
        << ", \"status\": " << JsonString(Name(result.status))
        << ", \"message\": " << JsonString(result.message)
        << ", \"metrics\": [";
    for (std::size_t j = 0; j < result.metrics.size(); ++j) {
      if (j) out << ',';
      out << "{\"name\":" << JsonString(result.metrics[j].first)
          << ",\"value\":" << std::setprecision(17)
          << result.metrics[j].second << '}';
    }
    out << "], \"artifacts\": [";
    for (std::size_t j = 0; j < result.artifacts.size(); ++j) {
      if (j) out << ',';
      out << JsonString(result.artifacts[j]);
    }
    out << "]}" << (i + 1 == results.size() ? "\n" : ",\n");
  }
  out << "  ]\n}\n";
  out.flush();
  if (!out.good())
    return Core::Status(Core::StatusCode::Error,
                        "failed while writing native-test JSON");
  out.close();
  fs::rename(temporary, target, error);
  if (error) {
    fs::remove(temporary);
    return Core::Status(Core::StatusCode::Error,
                        "cannot publish native-test JSON: " + error.message());
  }
  return Core::Status::OK();
}

int NativeTestExitCode(const std::vector<NativeTestResult>& results) {
  bool failed = false;
  for (const auto& result : results) {
    if (result.status == NativeTestStatus::Error) return 2;
    if (result.status == NativeTestStatus::Fail) failed = true;
  }
  return failed ? 1 : 0;
}

} }  // namespace SEP3D::Validation
