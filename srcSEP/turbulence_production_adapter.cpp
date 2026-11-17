#include "turbulence_production_adapter.h"
#include "polarized_wave_native_adapter.h"

#include "sep.h"
#include "transport_common.h"
#include "util/sep_background_runtime.h"
#include "util/sep_turbulence_core.h"
#include "util/sep_background_derivatives.h"
#include "util/sep_resonant_collision.h"
#include "adapters/swcme1d_adapter.h"
#include "adapters/reduced_shock_background_adapter.h"

#include <cmath>
#include <memory>
#include <sstream>
#include <vector>

namespace SEP {
namespace Turbulence {
namespace PICAdapter {
namespace {

ProductionLedger gLastStepLedger;
std::vector<ShockContribution> gShockContributions;
ShockSourceDiagnostics gShockDiagnostics;

double TotalPicWaveEnergyJ() {
  long double total = 0.0L;
  for (int fieldLine = 0; fieldLine < PIC::FieldLine::nFieldLine; ++fieldLine) {
    PIC::FieldLine::cFieldLine* line =
        &PIC::FieldLine::FieldLinesAll[fieldLine];
    for (int i = 0; i < line->GetTotalSegmentNumber(); ++i) {
      double* energy = line->GetSegment(i)->GetDatum_ptr(
          AlfvenTurbulence_Kolmogorov::CellIntegratedWaveEnergy);
      if (energy) total += energy[0] + energy[1];
    }
  }
  return static_cast<double>(total);
}

Transport::Status ImportFieldLine(int fieldLineId, State* state) {
  const auto authority=Polarized::PICAdapter::LegacyMutationPreflight("legacy two-branch core import");
  if(!authority.ok())return authority;
  namespace FL = PIC::FieldLine;
  FL::cFieldLine* line = &FL::FieldLinesAll[fieldLineId];
  const int segmentCount = line->GetTotalSegmentNumber();
  if (segmentCount <= 0)
    return Transport::Status::Error(Transport::StatusCode::OutOfDomain,
                                    "turbulence field line has no segments");

  state->configuration = ActiveConfiguration();
  state->cells.assign(static_cast<std::size_t>(segmentCount), CellState());
  for (int i = 0; i < segmentCount; ++i) {
    FL::cFieldLineSegment* segment = line->GetSegment(i);
    if (!segment)
      return Transport::Status::Error(Transport::StatusCode::OutOfDomain,
                                      "turbulence segment is missing");
    CellState& cell = state->cells[static_cast<std::size_t>(i)];
    cell.lengthM = segment->GetLength();
    cell.volumeM3 = SEP::FieldLine::FluxTubeGeometry::SegmentVolumeM3(
        segment, fieldLineId);

    double* energy = segment->GetDatum_ptr(
        AlfvenTurbulence_Kolmogorov::CellIntegratedWaveEnergy);
    if (!energy)
      return Transport::Status::Error(Transport::StatusCode::InvalidParticleState,
                                      "integrated turbulence datum is absent");
    cell.ePlusJ = energy[0];
    cell.eMinusJ = energy[1];

    // Apply queued shock records to pending source slots, never directly to the
    // host datum.  Advance() consumes these slots once and records the applied
    // amount in the same limiter/conservation ledger as every other operator.
    for (std::size_t source=0;source<gShockContributions.size();++source) {
      const ShockContribution& contribution=gShockContributions[source];
      if (contribution.fieldLine==fieldLineId && contribution.segment==i) {
        cell.pendingShockPlusJ+=contribution.plusJ;
        cell.pendingShockMinusJ+=contribution.minusJ;
      }
    }

    double magneticField[3] = {0.0, 0.0, 0.0};
    line->GetMagneticField(magneticField, static_cast<double>(i) + 0.5);
    cell.magneticFieldT = Vector3D::Length(magneticField);
    double numberDensity = 0.0;
    segment->GetPlasmaDensity(0.5, numberDensity);
    cell.massDensityKgPerM3 =
        numberDensity * PIC::CPLR::SWMF::MeanPlasmaAtomicMass;
    // The same provider composition used by particle scattering determines
    // wave propagation. Electron number density is not generally rho/m_p.
    const auto snapshot=Background::SnapshotStore::Instance().Current();
    if(!snapshot)return Transport::Status::Error(
        Transport::StatusCode::InvalidParticleState,
        "turbulence mapping requires a published background generation");
    double position[3];segment->GetCartesian(position,0.5);
    if(snapshot->provider()==Background::Provider::Swcme) {
      const auto query=SW1DAdapter::QueryAtRadius(Vector3D::Length(position),false);
      if(query.code!=SW1DAdapter::StatusCode::Success||
          query.epoch_seconds!=snapshot->epoch_seconds()||
          !query.sample.magnetic_and_composition_available())
        return Transport::Status::Error(Transport::StatusCode::InvalidParticleState,
            "turbulence SWCME composition differs from frozen background");
      cell.massDensityKgPerM3=numberDensity*query.sample.mass_density_kg_m3()/
          query.sample.number_density_m3();
    } else if(snapshot->provider()==Background::Provider::ReducedShock) {
      ReducedShock::AmbientSample sample;std::string error;
      if(ReducedShock::Metadata().epochS!=snapshot->epoch_seconds()||
          ReducedShock::Metadata().generation!=snapshot->field_line_generation()||
          !ReducedShock::EvaluateAmbient({{position[0],position[1],position[2]}},&sample,&error))
        return Transport::Status::Error(Transport::StatusCode::InvalidParticleState,
            "turbulence reduced composition differs from frozen background: "+error);
      cell.massDensityKgPerM3=numberDensity*sample.massDensityKgM3/sample.numberDensityM3;
    }
    if (!(cell.magneticFieldT > 0.0) ||
        !(cell.massDensityKgPerM3 > 0.0))
      return Transport::Status::Error(Transport::StatusCode::InvalidParticleState,
          "turbulence mapping requires positive B and mass density");
    cell.alfvenSpeedMPerS = cell.magneticFieldT /
        std::sqrt(VacuumPermeability * cell.massDensityKgPerM3);

    double velocity0[3] = {0.0, 0.0, 0.0};
    double velocity1[3] = {0.0, 0.0, 0.0};
    segment->GetBegin()->GetPlasmaVelocity(velocity0);
    segment->GetEnd()->GetPlasmaVelocity(velocity1);
    double direction[3] = {0.0, 0.0, 0.0};
    segment->GetDir(direction);
    const auto polarity=Transport::MagneticPolarityAlongFieldLine(
        magneticField,direction);
    if(!polarity.status.ok())return polarity.status;
    cell.magneticPolarityAlongLine=polarity.value>0?1:-1;
    cell.plasmaSpeedMPerS = 0.5 *
        (Vector3D::DotProduct(velocity0, direction) +
         Vector3D::DotProduct(velocity1, direction));

    if (state->configuration.representation == Representation::Spectral) {
      double* spectrum = segment->GetDatum_ptr(
          AlfvenTurbulence_Kolmogorov::WaveNumberResolved::SpectralWaveEnergy);
      if (!spectrum)
        return Transport::Status::Error(
            Transport::StatusCode::InvalidParticleState,
            "spectral turbulence authority is absent");
      const auto copied=Polarized::PICAdapter::CopyLegacyNativeSpectrum(
          state->configuration,spectrum,
          static_cast<std::size_t>(AlfvenTurbulence_Kolmogorov::WaveNumberResolved::SpectralWaveEnergy.length),
          &cell.spectralEnergyJ);
      if(!copied.ok())return copied;
    }
  }

  // Reflection needs d ln(V_A)/ds at cell centers.  One-sided differences are
  // used only at physical boundaries; all interior cells use a centered metric.
  for (int i = 0; i < segmentCount; ++i) {
    const int left = i == 0 ? i : i - 1;
    const int right = i + 1 == segmentCount ? i : i + 1;
    if (left == right) continue;
    double separationM = 0.0;
    for (int j = left; j < right; ++j)
      separationM += 0.5 * (state->cells[j].lengthM +
                            state->cells[j + 1].lengthM);
    state->cells[i].dLnAlfvenSpeeddsPerM =
        (std::log(state->cells[right].alfvenSpeedMPerS) -
         std::log(state->cells[left].alfvenSpeedMPerS)) / separationM;
  }
  const std::shared_ptr<const Background::BackgroundSnapshot> snapshot =
      Background::SnapshotStore::Instance().Current();
  if (!snapshot)
    return Transport::Status::Error(Transport::StatusCode::InvalidParticleState,
                                    "turbulence advance has no background snapshot");
  state->fieldLineGeneration = snapshot->field_line_generation();
  state->epochS = Background::SimulationTimeSeconds();
  std::ostringstream provenance;
  provenance << "PIC field-line adapter;configuration="
             << Background::CurrentConfigurationFingerprint();
  state->provenance = provenance.str();
  if (state->configuration.source == Source::SwmfInitialThenEvolveLocal)
    return HandoffImportedState(
        state, state->epochS, snapshot->configuration_fingerprint());
  return InitializeState(state);
}

Transport::Status ExportFieldLine(int fieldLineId, const State& state) {
  const auto authority=Polarized::PICAdapter::LegacyMutationPreflight("legacy two-branch core export");
  if(!authority.ok())return authority;
  if(state.configuration.representation==Representation::Spectral) {
    const auto grid=Polarized::PICAdapter::LegacyFixedSpectralPreflight(state.configuration,"legacy spectral core export");
    if(!grid.ok())return grid;
  }
  namespace FL = PIC::FieldLine;
  FL::cFieldLine* line = &FL::FieldLinesAll[fieldLineId];
  for (std::size_t i = 0; i < state.cells.size(); ++i) {
    FL::cFieldLineSegment* segment = line->GetSegment(static_cast<int>(i));
    double* energy = segment->GetDatum_ptr(
        AlfvenTurbulence_Kolmogorov::CellIntegratedWaveEnergy);
    double* derived = segment->GetDatum_ptr(
        AlfvenTurbulence_Kolmogorov::WaveEnergyDensity);
    if (!energy || !derived)
      return Transport::Status::Error(Transport::StatusCode::InvalidParticleState,
                                      "turbulence export datum is absent");
    const DerivedCell view = DeriveCell(state.cells[i]);
    if (!view.status.ok()) return view.status;
    energy[0] = state.cells[i].ePlusJ;
    energy[1] = state.cells[i].eMinusJ;
    derived[0] = view.wPlusJPerM3;
    derived[1] = view.wMinusJPerM3;
    derived[2] = view.crossHelicity;
    if (state.configuration.representation == Representation::Spectral) {
      double* spectrum = segment->GetDatum_ptr(
          AlfvenTurbulence_Kolmogorov::WaveNumberResolved::SpectralWaveEnergy);
      for (std::size_t k = 0; k < state.cells[i].spectralEnergyJ.size(); ++k)
        spectrum[k] = state.cells[i].spectralEnergyJ[k];
    }
  }
  return Transport::Status::Ok();
}

}  // namespace

CheckpointCache CaptureCheckpointCache() {
  CheckpointCache cache;cache.ledger=gLastStepLedger;
  cache.shockContributions=gShockContributions;cache.shockDiagnostics=gShockDiagnostics;
  return cache;
}
Transport::Status ValidateCheckpointCache(const CheckpointCache& cache) {
  const double entries[]={cache.ledger.shockSourceJ,cache.ledger.particleExchangeJ,
    cache.ledger.boundaryExchangeJ,cache.ledger.physicalDissipationJ,
    cache.ledger.limiterCorrectionJ,cache.ledger.closureResidualJ};
  for(double value:entries)if(!std::isfinite(value))return Transport::Status::Error(
      Transport::StatusCode::InvalidParticleState,"nonfinite turbulence checkpoint ledger");
  for(const auto& source:cache.shockContributions)
    if(source.fieldLine<0||source.fieldLine>=PIC::FieldLine::nFieldLine||
        source.segment<0||source.segment>=PIC::FieldLine::FieldLinesAll[source.fieldLine].GetTotalSegmentNumber()||
        !std::isfinite(source.plusJ)||!std::isfinite(source.minusJ)||
        source.plusJ<0.0||source.minusJ<0.0||source.provenance.empty())
      return Transport::Status::Error(Transport::StatusCode::InvalidParticleState,
          "invalid turbulence source checkpoint history");
  return Transport::Status::Ok();
}
void RestoreCheckpointCache(CheckpointCache* cache) {
  gLastStepLedger=cache->ledger;gShockDiagnostics=cache->shockDiagnostics;
  gShockContributions.swap(cache->shockContributions);
}

Transport::Status QueueShockContribution(const ShockContribution& c) {
  if (c.fieldLine<0 || c.segment<0 || !std::isfinite(c.plusJ) ||
      !std::isfinite(c.minusJ) || c.plusJ<0.0 || c.minusJ<0.0 ||
      c.provenance.empty()) {
    ++gShockDiagnostics.invalidPhysics;
    return Transport::Status::Error(Transport::StatusCode::InvalidArgument,
                                    "shock contribution record is invalid");
  }
  gShockContributions.push_back(c);
  ++gShockDiagnostics.accepted;
  return Transport::Status::Ok();
}

void ClearShockContributions() {
  gShockContributions.clear();
  gShockDiagnostics=ShockSourceDiagnostics();
}

const ShockSourceDiagnostics& LastShockSourceDiagnostics() {
  return gShockDiagnostics;
}

void RecordShockSourceRejection(bool geometryFailure) {
  if (geometryFailure) ++gShockDiagnostics.invalidGeometry;
  else ++gShockDiagnostics.noIntersection;
}

Transport::Status Advance(double dtS, double shockRadiusBeforeM,
                          double shockRadiusAfterM) {
  if (!(dtS > 0.0) || !std::isfinite(dtS))
    return Transport::Status::Error(Transport::StatusCode::InvalidArgument,
                                    "production turbulence dt must be positive");

  // Capability checks precede queue reduction, shock sources and every native
  // wave write. The approved polarized closure is currently a neutral tested
  // kernel; the legacy two-array updater cannot silently stand in for it.
  const auto capability=Transport::PICAdapter::ProductionClosurePreflight();
  if(!capability.ok())return capability;

  // This is the transaction boundary for all particle records.  No mover holds
  // a live pointer in the queue, and no turbulence operator can observe a
  // partially reduced worker set.
  Transport::PICAdapter::FlushWaveContributions();
  gLastStepLedger = ProductionLedger();
  ClearShockContributions();
  if (!AlfvenTurbulence_Kolmogorov::ActiveFlag)
    return Transport::Status::Ok();
  if (!EvolvesLocally(ActiveConfiguration().source))
    return Transport::Status::Ok();

  // The existing shock and resonant-particle managers are retained as narrow
  // source adapters.  They run inside this transaction before import; every
  // spatial transport/reflection/cascade mutation below is owned by the common
  // turbulence core and recorded in its signed energy ledger.
  if (ActiveConfiguration().shockInjectionEnabled &&
      std::isfinite(shockRadiusBeforeM) &&
      std::isfinite(shockRadiusAfterM) &&
      shockRadiusBeforeM != shockRadiusAfterM) {
    SEP::ParticleSource::ShockWave::ShockTurbulenceEnergyInjection(
        shockRadiusBeforeM, shockRadiusAfterM, dtS);
    if (AlfvenTurbulence_Kolmogorov::WaveNumberResolved::IsActive())
      AlfvenTurbulence_Kolmogorov::WaveNumberResolved::
          ProjectIntegratedEnergyToSpectrum(
              AlfvenTurbulence_Kolmogorov::CellIntegratedWaveEnergy);
  }
  const double afterShockJ = TotalPicWaveEnergyJ();

  if (AlfvenTurbulence_Kolmogorov::ParticleCouplingMode &&
      !SEP::Mover::CurrentCapabilities().evolvesWaveStateDirectly) {
    PIC::FieldLine::Parallel::MPIAllReduceDatumStoredAtEdge(
        AlfvenTurbulence_Kolmogorov::G_plus_streaming);
    PIC::FieldLine::Parallel::MPIAllReduceDatumStoredAtEdge(
        AlfvenTurbulence_Kolmogorov::G_minus_streaming);
    if (AlfvenTurbulence_Kolmogorov::WaveNumberResolved::IsActive()) {
      AlfvenTurbulence_Kolmogorov::WaveNumberResolved::
          WaveParticleCouplingManager(
              AlfvenTurbulence_Kolmogorov::CellIntegratedWaveEnergy, dtS);
      AlfvenTurbulence_Kolmogorov::WaveNumberResolved::
          UpdateIntegratedEnergyFromSpectrum(
              AlfvenTurbulence_Kolmogorov::CellIntegratedWaveEnergy);
    }
    else {
      AlfvenTurbulence_Kolmogorov::IsotropicSEP::WaveParticleCouplingManager(
          AlfvenTurbulence_Kolmogorov::CellIntegratedWaveEnergy, dtS);
    }
  }
  gLastStepLedger.particleExchangeJ =
      TotalPicWaveEnergyJ() - afterShockJ;

  for (int fieldLine = 0; fieldLine < PIC::FieldLine::nFieldLine; ++fieldLine) {
    State state;
    Transport::Status status = ImportFieldLine(fieldLine, &state);
    if (!status.ok()) return status;
    // Sources were imported after the legacy source adapters, so the common
    // core must not apply a second particle/shock increment.
    state.configuration.coupling = CouplingPolicy::Disabled;
    state.configuration.shockInjectionEnabled = false;
    const StepResult result = Turbulence::Advance(&state, dtS);
    if (!result.status.ok()) return result.status;
    gLastStepLedger.boundaryExchangeJ +=
        result.ledger.innerBoundaryJ + result.ledger.outerBoundaryJ;
    gLastStepLedger.shockSourceJ += result.ledger.shockSourceJ;
    gLastStepLedger.physicalDissipationJ +=
        result.ledger.physicalDissipationJ;
    gLastStepLedger.limiterCorrectionJ +=
        result.ledger.limiterCorrectionJ;
    gLastStepLedger.closureResidualJ += result.ledger.closureResidualJ;
    gLastStepLedger.coreSubcycles += result.diagnostics.subcycles;
    // Export deterministic work and correction counts with the physical
    // ledger.  Wall-clock timing belongs in the harness, whereas these counters
    // identify which operator caused a cost or robustness regression.
    gLastStepLedger.limiterActivations +=
        result.diagnostics.limiterActivations;
    gLastStepLedger.rejectedUpdates += result.diagnostics.rejectedUpdates;
    gLastStepLedger.advectionFaceUpdates +=
        result.diagnostics.advectionFaceUpdates;
    gLastStepLedger.reflectionCellUpdates +=
        result.diagnostics.reflectionCellUpdates;
    gLastStepLedger.cascadeCellUpdates +=
        result.diagnostics.cascadeCellUpdates;
    gLastStepLedger.spectralBinUpdates +=
        result.diagnostics.spectralBinUpdates;
    status = ExportFieldLine(fieldLine, state);
    if (!status.ok()) return status;
  }

  if (AlfvenTurbulence_Kolmogorov::WaveNumberResolved::IsActive())
    PIC::FieldLine::Parallel::MPIAllGatherDatumStoredAtEdge(
        AlfvenTurbulence_Kolmogorov::WaveNumberResolved::SpectralWaveEnergy);
  PIC::FieldLine::Parallel::MPIAllGatherDatumStoredAtEdge(
      AlfvenTurbulence_Kolmogorov::CellIntegratedWaveEnergy);
  return Transport::Status::Ok();
}

const ProductionLedger& LastStepLedger() { return gLastStepLedger; }

}  // namespace PICAdapter
}  // namespace Turbulence
}  // namespace SEP
