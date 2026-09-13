#include "tests.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace {

SEP::Testing::Result BooleanResult(bool passed, const char* passMessage,
                                   const char* failMessage) {
  SEP::Testing::Result result;
  result.status = passed ? SEP::Testing::Status::Pass
                         : SEP::Testing::Status::Fail;
  result.message = passed ? passMessage : failMessage;
  result.metrics.push_back(
      {"assertion_failures", passed ? 0.0 : 1.0, 0.0, "<=", "count"});
  return result;
}

SEP::Testing::Result RunDxx() {
  SEP::Testing::Result result = BooleanResult(DxxTest(),
      "analytical and quadrature spatial-diffusion checks passed",
      "spatial-diffusion analytical or quadrature check exceeded tolerance");
  result.configuration.push_back("speed_rounds_m_per_s=1e5,1e6,1e7,1e8");
  result.configuration.push_back("quadrature_panels=1000000");
  result.configuration.push_back("relative_tolerance=1e-5");
  return result;
}

SEP::Testing::Result RunFteConvection() {
  // Registry adapters use fixed seeds because their process exits after the
  // selected test set.  Re-seeding before each stochastic test makes results
  // independent of earlier selector overlap/order without perturbing any
  // subsequent production simulation (there is none in test-only mode).
  const std::uint64_t seed = 1002;
  rnd_seed(static_cast<int>(seed));
  SEP::Testing::Result result = BooleanResult(FTE_Convectoin(),
      "focused-transport convection preserved the expected displacement and velocity",
      "focused-transport convection exceeded its existing numerical tolerance");
  result.hasSeed = true;
  result.seed = seed;
  result.configuration.push_back("trials=10000");
  result.configuration.push_back("plasma_velocity_m_per_s=0,0,0");
  result.configuration.push_back("dt_s=1");
  return result;
}

SEP::Testing::Result RunParkerConvection() {
  const std::uint64_t seed = 1001;
  rnd_seed(static_cast<int>(seed));
  SEP::Testing::Result result = BooleanResult(ParkerModelMoverTest_convection(),
      "Parker convection and adiabatic momentum assertions passed",
      "Parker convection or adiabatic momentum assertion failed");
  result.hasSeed = true;
  result.seed = seed;
  result.configuration.push_back("trials=10000");
  result.configuration.push_back("density_jump=1_to_4");
  result.configuration.push_back("dt_s=1");
  return result;
}

SEP::Testing::Result RunParkerDistribution() {
  const std::uint64_t seed = 1003;
  rnd_seed(static_cast<int>(seed));
  SEP::Testing::Result result = BooleanResult(
      ParkerModelMoverTest_const_plasma_field(),
      "stochastic Parker campaign produced and wrote an in-range histogram",
      "stochastic Parker campaign produced no in-range histogram or output failed");
  result.hasSeed = true;
  result.seed = seed;
  result.configuration.push_back("trials=4000000");
  result.configuration.push_back("dt_s=3");
  result.configuration.push_back("histogram_bins=100");
  result.artifacts.push_back("dxParker.dat");
  return result;
}

SEP::Testing::Result RunScatteringBeyondOneAu() {
  // The historical routine performs 400,000 stochastic particle histories and
  // writes multiple diagnostics, but it defines no numerical acceptance rule
  // and includes an E=0 case whose mean-free-path expression is singular.  A
  // truthful registry must not infer PASS merely because such a routine
  // returned.  It remains discoverable/individually selectable and returns SKIP
  // until a later physics-validation step supplies an approved reference and
  // finite-energy domain; --run-test-manager retains the historical execution.
  SEP::Testing::Result result;
  result.status = SEP::Testing::Status::Skip;
  result.message =
      "legacy diagnostic has no authoritative acceptance criterion; use "
      "--run-test-manager for historical output execution";
  result.hasSeed = false;
  result.configuration.push_back("legacy_energy_grid_MeV=0,50,100,150");
  result.configuration.push_back("registry_execution=disabled_pending_reference");
  result.artifacts.push_back("rmax-E=<energy>MeV.dat");
  result.artifacts.push_back("time-E=<energy>MeV.dat");
  return result;
}

SEP::Testing::Result RunTurbulenceEnergyClosure() {
  // This assertion uses the public production helper and an independently
  // evaluated magnetic-energy expression.  It is intentionally initialization
  // free, deterministic, and non-mutating, so the turbulence group has one real
  // bounded component assertion without pretending to replace the later
  // transport/cascade/reflection validation campaign.
  const double magneticFieldTesla = 5.0e-9;
  const double fractionalFluctuation = 0.2;
  const double expected =
      std::pow(magneticFieldTesla * fractionalFluctuation, 2) /
      (2.0 * SEP::AlfvenTurbulence_Kolmogorov::WaveEnergyConstants::MU0);
  const double actual =
      SEP::AlfvenTurbulence_Kolmogorov::CalculateTypicalWaveEnergyDensity1AU(
          magneticFieldTesla, fractionalFluctuation);
  const double relativeError = std::fabs(actual - expected) / expected;
  const double tolerance = 32.0 * std::numeric_limits<double>::epsilon();

  SEP::Testing::Result result = BooleanResult(
      std::isfinite(actual) && relativeError <= tolerance,
      "wave-energy density matches the independent magnetic-pressure closure",
      "wave-energy density violates the magnetic-pressure closure");
  result.metrics.push_back(
      {"relative_error", relativeError, tolerance, "<=", "dimensionless"});
  result.configuration.push_back("B0_T=5e-9");
  result.configuration.push_back("deltaB_over_B0=0.2");
  return result;
}

SEP::Testing::Descriptor MakeDescriptor(
    const char* id, const char* name, const char* group,
    const char* description, SEP::Testing::InitializationLevel initialization,
    SEP::Testing::RuntimeClass runtime, const char* seedPolicy,
    const char* stateIsolation, SEP::Testing::TestCallback callback) {
  SEP::Testing::Descriptor descriptor;
  descriptor.id = id;
  descriptor.name = name;
  descriptor.group = group;
  descriptor.description = description;
  descriptor.initialization = initialization;
  descriptor.supportedBuildModes = "serial and MPI standalone executable";
  descriptor.runtime = runtime;
  descriptor.seedPolicy = seedPolicy;
  descriptor.stateIsolation = stateIsolation;
  descriptor.callback = callback;
  return descriptor;
}

int StatusSeverity(SEP::Testing::Status status) {
  switch (status) {
    case SEP::Testing::Status::Pass: return 0;
    case SEP::Testing::Status::Skip: return 1;
    case SEP::Testing::Status::Fail: return 2;
    case SEP::Testing::Status::Error: return 3;
  }
  return 3;
}

SEP::Testing::Status StatusFromSeverity(int severity) {
  if (severity == 0) return SEP::Testing::Status::Pass;
  if (severity == 1) return SEP::Testing::Status::Skip;
  if (severity == 2) return SEP::Testing::Status::Fail;
  return SEP::Testing::Status::Error;
}

}  // namespace

const SEP::Testing::Registry& ComponentTestRegistry() {
  // Function-local static construction avoids cross-translation-unit
  // initialization ordering.  Registry construction also validates every
  // descriptor and deterministically sorts this deliberately unsorted source
  // list before it becomes observable through --list-tests.
  static const SEP::Testing::Registry registry({
      MakeDescriptor("TURB01", "Alfven wave-energy closure", "turbulence",
          "Compare the production 1-AU wave-energy helper with an independent magnetic-pressure expression.",
          SEP::Testing::InitializationLevel::None,
          SEP::Testing::RuntimeClass::Routine, "deterministic; no RNG",
          "read-only pure calculation", RunTurbulenceEnergyClosure),
      MakeDescriptor("PARKER02", "Parker stochastic displacement histogram", "parker",
          "Run the legacy constant-plasma Parker campaign and validate histogram/output completion.",
          SEP::Testing::InitializationLevel::FieldLineModel,
          SEP::Testing::RuntimeClass::Extended, "fixed registry seed 1003",
          "adapter reseeds; particle and segment lists are cleaned; process exits after tests",
          RunParkerDistribution),
      MakeDescriptor("DXX01", "Spatial diffusion coefficient", "diffusion",
          "Compare GetDxx with its constant-coefficient solution and independent numerical quadrature.",
          SEP::Testing::InitializationLevel::FieldLineModel,
          SEP::Testing::RuntimeClass::Routine, "deterministic; no RNG",
          "pitch-angle coefficient function pointer is restored before return", RunDxx),
      MakeDescriptor("SCAT01", "Scattering beyond 1 AU", "scattering",
          "Inventory the legacy return-probability diagnostic and its reference-acceptance limitation.",
          SEP::Testing::InitializationLevel::FieldLineModel,
          SEP::Testing::RuntimeClass::Extended, "uses the configured AMPS random stream; seed not exposed",
          "registry returns SKIP without mutation; legacy TestManager still executes it",
          RunScatteringBeyondOneAu),
      MakeDescriptor("FTE01", "Focused-transport convection", "transport",
          "Check field-line streaming displacement and velocity invariance in a static plasma fixture.",
          SEP::Testing::InitializationLevel::FieldLineModel,
          SEP::Testing::RuntimeClass::Routine, "fixed registry seed 1002",
          "adapter reseeds; plasma/coefficient/particle/list state is restored; process exits after tests", RunFteConvection),
      MakeDescriptor("PARKER01", "Parker convection and adiabatic momentum", "parker",
          "Check stationary line coordinates and the density-driven analytical momentum update.",
          SEP::Testing::InitializationLevel::FieldLineModel,
          SEP::Testing::RuntimeClass::Routine, "fixed registry seed 1001",
          "adapter reseeds; plasma/coefficient/particle/list state is restored; process exits after tests", RunParkerConvection),
  });
  return registry;
}

SEP::Testing::InitializationLevel RequiredInitializationLevel(
    const std::vector<const SEP::Testing::Descriptor*>& selected) {
  SEP::Testing::InitializationLevel required =
      SEP::Testing::InitializationLevel::None;
  for (const SEP::Testing::Descriptor* descriptor : selected) {
    if (descriptor &&
        descriptor->initialization ==
            SEP::Testing::InitializationLevel::FieldLineModel) {
      required = SEP::Testing::InitializationLevel::FieldLineModel;
    }
  }
  return required;
}

int RunSelectedComponentTests(
    const std::vector<const SEP::Testing::Descriptor*>& selected,
    std::ostream& out) {
  SEP::Testing::Summary summary;
  const SEP::Testing::Registry& registry = ComponentTestRegistry();

  for (const SEP::Testing::Descriptor* descriptor : selected) {
    SEP::Testing::Result result = registry.RunOne(*descriptor);

    // Every required rank executes the callback.  Reduce the categorical result
    // by severity and duration by maximum so rank-local failure cannot be hidden
    // by a root PASS and the reported time reflects the slowest participant.
    int initialized = 0;
    MPI_Initialized(&initialized);
    if (initialized) {
      const int localSeverity = StatusSeverity(result.status);
      int globalSeverity = localSeverity;
      double globalDuration = result.elapsedSeconds;
      MPI_Allreduce(&localSeverity, &globalSeverity, 1, MPI_INT, MPI_MAX,
                    MPI_GLOBAL_COMMUNICATOR);
      MPI_Allreduce(&result.elapsedSeconds, &globalDuration, 1, MPI_DOUBLE,
                    MPI_MAX, MPI_GLOBAL_COMMUNICATOR);
      if (globalSeverity != localSeverity) {
        result.message = "one or more MPI ranks reported a more severe outcome";
      }
      result.status = StatusFromSeverity(globalSeverity);
      result.elapsedSeconds = globalDuration;
    }

    if (PIC::ThisThread == 0) SEP::Testing::PrintResult(result, out);
    summary.results.push_back(result);
    if (result.status == SEP::Testing::Status::Pass) summary.passed++;
    else if (result.status == SEP::Testing::Status::Fail) summary.failed++;
    else if (result.status == SEP::Testing::Status::Skip) summary.skipped++;
    else summary.errors++;
  }

  if (PIC::ThisThread == 0) {
    out << "Component-test summary: PASS=" << summary.passed
        << " FAIL=" << summary.failed << " SKIP=" << summary.skipped
        << " ERROR=" << summary.errors << '\n';
  }
  return summary.ExitCode();
}
