#include "turbulence_production_adapter.h"

#include "sep.h"
#include "transport_common.h"
#include "util/sep_background_runtime.h"
#include "util/sep_turbulence_core.h"

#include <cmath>
#include <memory>
#include <sstream>
#include <vector>

namespace SEP {
namespace Turbulence {
namespace PICAdapter {
namespace {

ProductionLedger gLastStepLedger;

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

    double magneticField[3] = {0.0, 0.0, 0.0};
    line->GetMagneticField(magneticField, static_cast<double>(i) + 0.5);
    cell.magneticFieldT = Vector3D::Length(magneticField);
    double numberDensity = 0.0;
    segment->GetPlasmaDensity(0.5, numberDensity);
    cell.massDensityKgPerM3 =
        numberDensity * PIC::CPLR::SWMF::MeanPlasmaAtomicMass;
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
      const std::size_t bins = state->configuration.spectralBins;
      cell.spectralEnergyJ.assign(spectrum, spectrum + 2 * bins);
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
    energy[0] = state.cells[i].ePlusJ;
    energy[1] = state.cells[i].eMinusJ;
    const DerivedCell view = DeriveCell(state.cells[i]);
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

Transport::Status Advance(double dtS, double shockRadiusBeforeM,
                          double shockRadiusAfterM) {
  if (!(dtS > 0.0) || !std::isfinite(dtS))
    return Transport::Status::Error(Transport::StatusCode::InvalidArgument,
                                    "production turbulence dt must be positive");

  // This is the transaction boundary for all particle records.  No mover holds
  // a live pointer in the queue, and no turbulence operator can observe a
  // partially reduced worker set.
  Transport::PICAdapter::FlushWaveContributions();
  gLastStepLedger = ProductionLedger();
  if (!AlfvenTurbulence_Kolmogorov::ActiveFlag)
    return Transport::Status::Ok();
  if (!EvolvesLocally(ActiveConfiguration().source))
    return Transport::Status::Ok();

  // The existing shock and resonant-particle managers are retained as narrow
  // source adapters.  They run inside this transaction before import; every
  // spatial transport/reflection/cascade mutation below is owned by the common
  // turbulence core and recorded in its signed energy ledger.
  const double beforeSourcesJ = TotalPicWaveEnergyJ();
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
  gLastStepLedger.shockSourceJ = afterShockJ - beforeSourcesJ;

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
    gLastStepLedger.physicalDissipationJ +=
        result.ledger.physicalDissipationJ;
    gLastStepLedger.limiterCorrectionJ +=
        result.ledger.limiterCorrectionJ;
    gLastStepLedger.closureResidualJ += result.ledger.closureResidualJ;
    gLastStepLedger.coreSubcycles += result.diagnostics.subcycles;
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
