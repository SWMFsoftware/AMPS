#include "sep_background_runtime.h"
#include "sep_provider_radius_frame.h"

#include "../sep.h"
#include "../adapters/reduced_shock_background_adapter.h"
#include "../adapters/swcme1d_adapter.h"
#include "sep_run_configuration.h"
#include "sep_analytic_radial_background.h"
#include "sep_standard_batsrus_waves.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>
#include <stdexcept>

#if _PIC_COUPLER_MODE_ == _PIC_COUPLER_MODE__SWMF_
#include "amps2swmf.h"
#endif

namespace SEP {
namespace Background {

namespace {

std::string RunIdentityExcludingMover() {
  // The full run/checkpoint fingerprint intentionally includes transport
  // choice and provenance. Background identity excludes precisely those two
  // coordinates: canonicalize a PRIVATE COPY, never install or execute this
  // sentinel mover and never modify the authoritative frozen configuration.
  // Every other existing run coordinate retains its original identity effect.
  auto projected=Run::Active().get();
  projected.mover.value=Mover::ProductionMover::Parker;
  projected.mover.source=Run::ValueSource::Default;
  return Run::Fingerprint(projected);
}

double SwmfEpochSeconds() {
#if _PIC_COUPLER_MODE_ == _PIC_COUPLER_MODE__SWMF_
  return AMPS2SWMF::MagneticFieldLineUpdate::LastCouplingTime;
#else
  throw std::logic_error(
      "SWMF background requested by a build without the SWMF coupler");
#endif
}

double PreviousSwmfEpochSeconds(double currentEpochS) {
#if _PIC_COUPLER_MODE_ == _PIC_COUPLER_MODE__SWMF_
  return AMPS2SWMF::MagneticFieldLineUpdate::SecondCouplingFlag
      ? AMPS2SWMF::MagneticFieldLineUpdate::LastLastCouplingTime
      : currentEpochS;
#else
  return currentEpochS;
#endif
}

}  // namespace

Provider ConfiguredProvider() {
#if _PIC_COUPLER_MODE_ == _PIC_COUPLER_MODE__SWMF_
  if (PIC::CPLR::SWMF::BlCouplingFlag) return Provider::Swmf;
#endif

  // The reduced event is an explicit application selection and therefore
  // precedes the legacy particle-source shock switch.  Keeping these controls
  // orthogonal is essential: this background-only stage must not reinterpret
  // a baseline SWCME/analytic particle source as a reduced-front source.
  if (SEP::ReducedShock::Enabled()) return Provider::ReducedShock;

  return SEP::ShockModelType == SEP::cShockModelType::SwCme1d
             ? Provider::Swcme
             : Provider::Analytic;
}

double SimulationTimeSeconds() {
  // PIC advances standalone time; the SWMF boundary supplies its explicitly
  // bound fixed-step clock phase. This is the sole read-only srcSEP clock API:
  // no Set(), Advance(), or independently accumulated elapsed-time variable.
  return PIC::SimulationTime::Get();
}

std::string CurrentConfigurationFingerprint() {
  std::ostringstream canonical;

  // Only settings that define the background realization belong here.  The
  // particle mover and its scattering parameterization are intentionally not
  // included, preserving identical provider state in cross-mover comparisons.
  canonical << "schema=srcsep-background-v3"
            << ";provider=" << ProviderName(ConfiguredProvider())
            << ";particle-coordinate-frame=" << AMPS::FieldLine::CoordinateFrameName(ParticleCoordinateFrame())
            << ";coefficient-radius-frame=" << RadiusFrame::ContractIdentity(ConfiguredProvider())
            << ";reduced-event="
            << (ReducedShock::Enabled() ? ReducedShock::EventIdentity() : "none")
            << ";domain=" << SEP::DomainType
            << ";imf=" << SEP::ModeIMF
            << ";shock="
            << (SEP::ShockModelType == SEP::cShockModelType::SwCme1d
                    ? "swcme-1d"
                    : "analytic-1d")
            // The complete frozen run fingerprint includes the D01 policy and
            // fallback SI values. Naming the policy again keeps human-readable
            // provenance intelligible without decoding the hash.
            << ";run-configuration-excluding-mover=" << RunIdentityExcludingMover()
            << ";swcme-failure-policy="
            << SEP::SW1DAdapter::FailurePolicyName(
                   SEP::SW1DAdapter::GetFailurePolicy())
            << ";turbulence-active="
            << (SEP::AlfvenTurbulence_Kolmogorov::ActiveFlag ? 1 : 0)
            << ";turbulence-representation="
            << (SEP::AlfvenTurbulence_Kolmogorov::WaveNumberResolved::IsActive()
                    ? "wave-number-resolved"
                    : "integrated")
            << ";field-line-mode=" << _PIC_FIELD_LINE_MODE_
            << ";coupler-mode=" << _PIC_COUPLER_MODE_;

  if(ConfiguredProvider()==Provider::Analytic) {
    // Only the physical origin/speed enters from the initialization owner.
    // Its complete fingerprint also contains the mover and output choices,
    // which must not redefine this background when transport is switched.
    canonical<<";analytic-radial="<<AnalyticRadial::CanonicalDeclaration(
        AnalyticRadial::DeclaredNativeWind());
  }
  if(ConfiguredProvider()==Provider::ReducedShock) {
    // The legacy adapter's ambiguous numberDensityM3 remains an electron
    // compatibility field, whereas native transport now stores explicit np.
    // Prevent an older ne-based checkpoint from sharing this interpretation.
    canonical<<";reduced-native-density=proton-v1;reduced-query=interval-v1";
  }
  if(ConfiguredProvider()==Provider::Swmf) {
    canonical<<";imported-wave-contract="<<Waves::StandardBatsrus::ContractIdentity();
#if _PIC_COUPLER_MODE_ == _PIC_COUPLER_MODE__SWMF_
    // Include the actual provider normalization/frame and the reconstruction
    // law in immutable provenance. A same-wire-layout checkpoint cannot
    // silently resume the previous unresolved-origin/density-proxy contract.
    canonical<<";swmf-origin-grid="<<AMPS2SWMF::FieldLineData::SourceGridDeclaration()
      <<";swmf-gradient=cartesian-inverse-distance-affine-pivoted-QR-v1"
      <<";swmf-gradient-neighbors=32..512;swmf-gradient-min-lines=4;swmf-gradient-cond-limit=1e6";
#endif
  }

  return FingerprintConfiguration(canonical.str());
}

AMPS::FieldLine::CoordinateFrame ParticleCoordinateFrame() {
  const auto provider=ConfiguredProvider();
  const auto frame=DeclaredParticleCoordinateFrame(provider);
  if(provider==Provider::Swmf) {
#if _PIC_COUPLER_MODE_ == _PIC_COUPLER_MODE__SWMF_
    // Geometry provenance is validated, not used to choose a different
    // particle equation. A fixed imported grid cannot silently restore U*dt
    // in the srcSEP field-line movers; its owner must supply material geometry.
    if(AMPS2SWMF::FieldLineData::ParticleCoordinateFrame()!=frame)
      throw std::logic_error("srcSEP requires plasma-comoving SWMF field-line geometry");
#else
    throw std::logic_error("SWMF coordinate-frame declaration is unavailable in this build");
#endif
  }
  return frame;
}

void PublishModelOwnedSnapshot(Provider provider, double epoch_seconds,
                               double valid_until_seconds,
                               const std::string& provenance,
                               std::uint64_t explicit_generation) {
  if (provider != Provider::Analytic && provider != Provider::Swcme &&
      provider != Provider::ReducedShock) {
    throw std::invalid_argument(
        "model-owned publication accepts only analytic, SWCME, or reduced-shock providers");
  }

  if (ConfiguredProvider() != provider) {
    throw std::logic_error(
        "model-owned publication does not match the configured provider");
  }

  SnapshotStore& store = SnapshotStore::Instance();
  const std::shared_ptr<const BackgroundSnapshot> current = store.Current();

  // Analytic and SWCME profile updates do not rebuild field-line geometry, so
  // they retain its generation.  The first imported/constructed line set is
  // generation 1; zero remains reserved as an invalid/uninitialized sentinel.
  // Analytic/SWCME updates historically retain one geometry generation.  The
  // reduced provider has an event-defined physical generation even though its
  // line topology stays fixed, so its caller supplies that nonzero identity.
  std::uint64_t generation=explicit_generation!=0
      ? explicit_generation
      : (current?current->field_line_generation():UINT64_C(1));
  // A material line's geometry generation changes when it actually moves or
  // changes topology, independently of a provider cache/output cadence. The
  // initial registered geometry is generation one in snapshot numbering.
  for(int line=0;line<PIC::FieldLine::nFieldLine;++line) {
    const auto& geometry=PIC::FieldLine::FieldLinesAll[line];
    if(!geometry.ExternallyManagedGeometry)continue;
    if(geometry.GeometryGeneration==std::numeric_limits<std::uint64_t>::max())
      throw std::logic_error("material snapshot geometry generation overflow");
    generation=std::max(generation,geometry.GeometryGeneration+1);
  }

  BackgroundSnapshot snapshot(
      provider, Ownership::ModelOwned, epoch_seconds, epoch_seconds,
      valid_until_seconds, generation, CurrentConfigurationFingerprint(),
      provenance,
      current ? current->current_physical_epoch_seconds() : epoch_seconds,
      epoch_seconds);
  store.Publish(snapshot);
}

void PrepareSnapshotForParticleStep() {
  SnapshotStore& store = SnapshotStore::Instance();
  const Provider configured_provider = ConfiguredProvider();
  const double simulation_time = SimulationTimeSeconds();

  if (!std::isfinite(simulation_time)) {
    throw std::logic_error("PIC simulation clock returned a non-finite value");
  }
  // Zero bulk particle advection requires an actual material geometry owner.
  // Provider names and nonzero endpoint U do not establish that owner. Reject
  // an unregistered analytical/imported line before a production particle can
  // silently lose its bulk motion; tests must satisfy the same invariant.
  for(int line=0;line<PIC::FieldLine::nFieldLine;++line)
    if(!PIC::FieldLine::FieldLinesAll[line].ExternallyManagedGeometry)
      throw std::logic_error("srcSEP plasma-comoving particles require an explicitly registered material geometry owner");

  std::shared_ptr<const BackgroundSnapshot> current = store.Current();

  if (!current) {
    if (configured_provider == Provider::Swmf) {
      const double coupling_epoch = SwmfEpochSeconds();
      if (!std::isfinite(coupling_epoch) || coupling_epoch > simulation_time) {
        throw std::logic_error(
            "no usable SWMF coupling epoch is available for the particle step");
      }

      // The newest imported coupling state remains authoritative until SWMF
      // supplies another generation.  It is explicitly read-only and valid from
      // its coupling epoch through that as-yet-unknown replacement time.
      store.Publish(BackgroundSnapshot(
          Provider::Swmf, Ownership::ImportedReadOnly, coupling_epoch,
          coupling_epoch, std::numeric_limits<double>::infinity(), UINT64_C(1),
          CurrentConfigurationFingerprint(), "SWMF magnetic-field-line import",
          PreviousSwmfEpochSeconds(coupling_epoch), coupling_epoch));
    } else {
      // Library users that do not pass through standalone main.cpp may still
      // run a static analytic background.  SWCME, however, must publish its
      // prepared StepState explicitly so metadata and physical state share the
      // same epoch.
      if (configured_provider == Provider::Swcme ||
          configured_provider == Provider::ReducedShock) {
        throw std::logic_error(
            "time-dependent model particle step requested before its prepared state was published");
      }

      store.Publish(BackgroundSnapshot(
          Provider::Analytic, Ownership::ModelOwned, simulation_time,
          simulation_time, std::numeric_limits<double>::infinity(), UINT64_C(1),
          CurrentConfigurationFingerprint(), "standalone analytic background",
          simulation_time, simulation_time));
    }

    current = store.Current();
  }

  if (current->provider() != configured_provider) {
    throw std::logic_error(
        "configured provider differs from the published background; use an explicit handoff");
  }

  if (configured_provider == Provider::Swmf) {
    const double coupling_epoch = SwmfEpochSeconds();
    if (!std::isfinite(coupling_epoch) || coupling_epoch > simulation_time) {
      throw std::logic_error(
          "SWMF coupling epoch is invalid for the authoritative simulation clock");
    }

    if (coupling_epoch != current->epoch_seconds()) {
      // A changed SWMF coupling epoch denotes newly imported authoritative
      // field-line/plasma storage.  Incrementing the generation lets movers and
      // diagnostics identify that replacement even if topology is unchanged.
      store.Publish(BackgroundSnapshot(
          Provider::Swmf, Ownership::ImportedReadOnly, coupling_epoch,
          coupling_epoch, std::numeric_limits<double>::infinity(),
          current->field_line_generation() + UINT64_C(1),
          CurrentConfigurationFingerprint(), "SWMF magnetic-field-line import",
          PreviousSwmfEpochSeconds(coupling_epoch), coupling_epoch));
      current = store.Current();
    }
  }

  if (!current->Covers(simulation_time)) {
    throw std::logic_error(
        "published background validity interval does not cover PIC simulation time");
  }
}

void PublishLocalEvolutionHandoff(double epoch_seconds,
                                  double valid_until_seconds,
                                  const std::string& provenance) {
  SnapshotStore& store = SnapshotStore::Instance();
  const std::shared_ptr<const BackgroundSnapshot> current = store.Current();

  if (!current || current->provider() != Provider::Swmf ||
      current->ownership() != Ownership::ImportedReadOnly) {
    throw std::logic_error(
        "local evolution handoff requires a current read-only SWMF snapshot");
  }

  // The caller must first copy the imported arrays into private storage.  This
  // metadata handoff then makes the ownership transition explicit and records
  // the source generation in provenance supplied by that caller.
  BackgroundSnapshot handoff(
      Provider::LocalEvolution, Ownership::HandoffCopy, epoch_seconds,
      epoch_seconds, valid_until_seconds,
      current->field_line_generation() + UINT64_C(1),
      current->configuration_fingerprint(), provenance,
      current->current_physical_epoch_seconds(), epoch_seconds);
  store.PublishHandoff(handoff, Provider::Swmf);
}

}  // namespace Background
}  // namespace SEP
