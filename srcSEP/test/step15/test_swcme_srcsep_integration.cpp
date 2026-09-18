// The runner supplies AMPS/src/models/swcme as an include root. A bare public
// include proves this validation consumes the same canonical provider as the
// two production SEP applications rather than an application-local copy.
#include "swcme_sep_interface.hpp"
#include "../../util/sep_background_snapshot.h"
#include "../../util/sep_parker_core.h"
#include "../../util/sep_test_registry.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

namespace {

// VAL04 exercises the production Parker stepping function with a background
// obtained from the real SWCME 1-D model and its public SEP adapter.  Setting
// kappa to zero intentionally isolates the coupling contract: any difference
// from the independently accumulated characteristic must come from unit,
// sign, epoch, ownership, or field-selection errors rather than Monte-Carlo
// sampling noise.
class ZeroSpatialDiffusion final
    : public SEP::Transport::SpatialDiffusionProvider {
 public:
  SEP::Transport::SpatialDiffusionSample Evaluate(double, double) const override {
    SEP::Transport::SpatialDiffusionSample sample;
    sample.status = SEP::Transport::Status::Ok();
    sample.kappaParallelM2PerS = 0.0;
    sample.dKappaParallelDsMPerS = 0.0;
    sample.provenance = "validation:zero-spatial-diffusion-v1";
    return sample;
  }
};

double RelativeError(double actual, double expected) {
  return std::fabs(actual - expected) /
      std::max(std::fabs(expected), std::numeric_limits<double>::min());
}

SEP::Testing::Result RunVal04SwcmeReplay() {
  const std::uint64_t seed = 1504001;
  const double first_epoch_s = 2.0 * 3600.0;
  const double cadence_s = 300.0;
  const std::size_t cadence_count = 24;
  const double reference_radius_m = 0.40 * swcme::constants::AU_M;
  const double initial_momentum_kg_m_per_s = 4.0e-19;

  // Every physical control used by the coupled case is explicit here and is
  // repeated in the structured Result below.  ShockOnly + Source is SWCME's
  // controlled SEP configuration: the shock is an injection surface and is
  // not simultaneously represented as a compression in the transport field.
  swcme1d::Params parameters;
  parameters.V_sw_kms = 410.0;
  parameters.n1AU_cm3 = 5.0;
  parameters.B1AU_nT = 5.0;
  parameters.T_K = 1.0e5;
  parameters.gamma_ad = 5.0 / 3.0;
  parameters.sin_theta = 1.0;
  parameters.kinematics_mode = swcme::kinematics::Mode::DBM;
  parameters.r0_Rs = 20.0;
  parameters.V0_sh_kms = 1450.0;
  parameters.Gamma_kmInv = 1.0e-7;
  parameters.region_mode = swcme::regions::Mode::ShockOnly;
  parameters.shock_acceleration_mode = swcme::acceleration::Mode::Source;
  parameters.relative_source_weight_per_area = 1.0;

  swcme::sep::SpectrumConfig spectrum;
  spectrum.kinetic_energy_min_MeV = 1.0;
  spectrum.kinetic_energy_max_MeV = 1000.0;
  spectrum.reference_energy_MeV = 10.0;
  spectrum.normalization = swcme::sep::NormalizationMode::RelativeOnly;
  const swcme::sep::Interface1D swcme_adapter(parameters, spectrum);
  const std::string fingerprint = SEP::Background::FingerprintConfiguration(
      swcme_adapter.resolved_manifest());

  SEP::Background::SnapshotStore& snapshots =
      SEP::Background::SnapshotStore::Instance();
  snapshots.ResetForTests();
  ZeroSpatialDiffusion diffusion;
  SEP::Transport::ParkerState production_state(
      0.0, initial_momentum_kg_m_per_s);
  double reference_arc_length_m = 0.0;
  double reference_momentum_kg_m_per_s = initial_momentum_kg_m_per_s;
  double maximum_arc_relative_error = 0.0;
  double maximum_momentum_relative_error = 0.0;
  double minimum_density_m3 = std::numeric_limits<double>::infinity();
  double maximum_density_m3 = 0.0;
  bool status_ok = true;

  for (std::size_t cadence = 0; cadence < cadence_count; ++cadence) {
    const double epoch_s = first_epoch_s + cadence * cadence_s;
    const double query_radius_m = reference_radius_m + reference_arc_length_m;
    const swcme::sep::Interface1D::PreparedStep prepared =
        swcme_adapter.prepare(epoch_s);
    swcme::sep::BackgroundState background;
    const swcme::ModelStatus background_status =
        swcme_adapter.evaluate_background(prepared, query_radius_m, background);
    if (!background_status.ok()) {
      status_ok = false;
      break;
    }

    minimum_density_m3 = std::min(minimum_density_m3, background.density_m3);
    maximum_density_m3 = std::max(maximum_density_m3, background.density_m3);

    // A new immutable snapshot is published at every actual SWCME cadence.
    // The mover consumes it under ParticleReadPhase, which proves that the
    // provider/ownership/fingerprint/generation contract used by production
    // background adapters is exercised rather than bypassed by the test.
    const SEP::Background::BackgroundSnapshot snapshot(
        SEP::Background::Provider::Swcme,
        SEP::Background::Ownership::ModelOwned,
        epoch_s, epoch_s, epoch_s + cadence_s,
        static_cast<std::uint64_t>(cadence + 1), fingerprint,
        "SWCME Interface1D resolved configuration and prepared state");
    snapshots.Publish(snapshot);
    {
      SEP::Background::ParticleReadPhase read_phase =
          snapshots.BeginParticleRead(epoch_s);
      const SEP::Background::BackgroundSnapshot& mover_snapshot =
          snapshots.AcquireForMover();
      status_ok = status_ok && mover_snapshot.provider() ==
          SEP::Background::Provider::Swcme;
      status_ok = status_ok && mover_snapshot.configuration_fingerprint() ==
          fingerprint && mover_snapshot.Covers(epoch_s);

      // A keyed stream is still supplied because it is part of the production
      // mover API.  With kappa=0 no random value contributes to the result,
      // making this coupled replay deterministic and bitwise reproducible.
      SEP::Transport::KeyedRandomStream random(seed, 0, cadence, 15041);
      const SEP::Transport::ParkerIncrement increment =
          SEP::Transport::AdvanceParker(
              production_state,
              SEP::Transport::ParkerBackground(
                  background.velocity_m_s[0],
                  background.div_velocity_s_inv),
              1.0e7, cadence_s, diffusion, random);
      if (!increment.status.ok()) {
        status_ok = false;
        break;
      }
      production_state = increment.state;
    }

    // This characteristic is deliberately accumulated outside AdvanceParker.
    // It uses only SI fields emitted by SWCME and the mathematical definitions
    // ds/dt=U and dp/dt=-(div U)p/3, so the comparison catches adapter unit or
    // sign errors without sharing the production mover implementation.
    reference_arc_length_m += background.velocity_m_s[0] * cadence_s;
    reference_momentum_kg_m_per_s *=
        std::exp(-background.div_velocity_s_inv * cadence_s / 3.0);
    maximum_arc_relative_error = std::max(
        maximum_arc_relative_error,
        RelativeError(production_state.arcLengthM, reference_arc_length_m));
    maximum_momentum_relative_error = std::max(
        maximum_momentum_relative_error,
        RelativeError(production_state.momentumKgMPerS,
                      reference_momentum_kg_m_per_s));
  }

  // The same prepared production model must also yield a physical SEP source;
  // otherwise this would only be a background getter smoke test.  The source
  // check demonstrates the real SWCME shock/source contract at the coupling
  // boundary, while the characteristic above demonstrates field consumption.
  const swcme::sep::Interface1D::PreparedStep source_step =
      swcme_adapter.prepare(first_epoch_s);
  swcme::sep::SEPSourceState source;
  const swcme::ModelStatus source_status =
      swcme_adapter.source_at_shock(source_step, source);
  const bool source_ok = source_status.ok() && source.active &&
      source.acceleration_mode == swcme::acceleration::Mode::Source &&
      source.compression > 1.0 && source.fast_mach > 1.0 &&
      std::isfinite(source.q_phase_space);
  snapshots.ResetForTests();

  const double floating_tolerance =
      128.0 * std::numeric_limits<double>::epsilon();
  const bool pass = status_ok && source_ok &&
      std::isfinite(minimum_density_m3) && minimum_density_m3 > 0.0 &&
      maximum_density_m3 >= minimum_density_m3 &&
      maximum_arc_relative_error <= floating_tolerance &&
      maximum_momentum_relative_error <= floating_tolerance;

  SEP::Testing::Result result;
  result.status = pass ? SEP::Testing::Status::Pass : SEP::Testing::Status::Fail;
  result.message = pass
      ? "actual SWCME states drive srcSEP Parker transport and source boundary"
      : "SWCME-to-srcSEP replay violated a coupling or characteristic invariant";
  result.hasSeed = true;
  result.seed = seed;
  result.configuration.push_back("producer=swcme::sep::Interface1D");
  result.configuration.push_back("consumer=SEP::Transport::AdvanceParker");
  result.configuration.push_back("kinematics=DBM;r0_Rs=20;V0_sh_kms=1450;Gamma_kmInv=1e-7");
  result.configuration.push_back("region=SHOCK_ONLY;acceleration=SOURCE");
  result.configuration.push_back("first_epoch_s=7200;cadence_s=300;cadence_count=24");
  result.configuration.push_back("reference_radius_AU=0.40;kappa_parallel_m2_per_s=0");
  result.configuration.push_back("configuration_fingerprint=" + fingerprint);
  result.metrics.push_back({"arc_length_relative_error",
      maximum_arc_relative_error, floating_tolerance, "<=", "dimensionless"});
  result.metrics.push_back({"momentum_relative_error",
      maximum_momentum_relative_error, floating_tolerance, "<=", "dimensionless"});
  result.metrics.push_back({"minimum_background_density",
      minimum_density_m3, 0.0, ">", "m^-3"});
  result.metrics.push_back({"active_source_records",
      source_ok ? 1.0 : 0.0, 1.0, "==", "count"});
  result.metrics.push_back({"assertion_failures",
      pass ? 0.0 : 1.0, 0.0, "<=", "count"});
  return result;
}

SEP::Testing::Descriptor Val04Descriptor() {
  SEP::Testing::Descriptor descriptor;
  descriptor.id = "VAL04";
  descriptor.name = "Real SWCME-to-srcSEP background replay";
  descriptor.group = "validation";
  descriptor.description =
      "Replay actual SWCME SI background/source states through srcSEP Parker transport.";
  descriptor.initialization = SEP::Testing::InitializationLevel::None;
  descriptor.supportedBuildModes = "source-only C++17 sanitizer runner";
  descriptor.runtime = SEP::Testing::RuntimeClass::Routine;
  descriptor.seedPolicy = "fixed API-completeness seed 1504001; stochastic term disabled";
  descriptor.stateIsolation =
      "test-owned SWCME model and reset background snapshot store";
  descriptor.callback = RunVal04SwcmeReplay;
  return descriptor;
}

}  // namespace

int main() {
  const SEP::Testing::Registry registry({Val04Descriptor()});
  const std::vector<const SEP::Testing::Descriptor*> selected =
      registry.Select({"VAL04"}, {}, false);
  const SEP::Testing::Summary summary = registry.Run(selected, std::cout);
  std::string error;
  if (!SEP::Testing::WriteJsonSummary(
          summary, "step15-swcme-results.json", &error) ||
      !SEP::Testing::WriteJUnitSummary(
          summary, "step15-swcme-results.xml", &error)) {
    std::cerr << "Report error: " << error << '\n';
    return 2;
  }
  return summary.ExitCode();
}
