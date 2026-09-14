#include "turbulence_production_adapter.h"

#include "sep.h"
#include "transport_common.h"
#include "util/sep_background_runtime.h"
#include "util/sep_runtime_contracts.h"
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
std::vector<ShockContribution> gShockContributions;
ShockSourceDiagnostics gShockDiagnostics;
RuntimeContracts::TurbulenceRuntimeStore gRuntimeStore;

Transport::Status ReadFieldLineView(int fieldLineId, bool importWaveAuthority,
                                    State* state) {
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
    // Host arrays are imported only when the runtime owner is first created.
    // Later calls sample geometry/background fields but leave the authoritative
    // E+/E-, spectra, pending transactions, and ledgers in gRuntimeStore.
    cell.ePlusJ = importWaveAuthority ? energy[0] : 0.0;
    cell.eMinusJ = importWaveAuthority ? energy[1] : 0.0;

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
      if (importWaveAuthority)
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
  if (!importWaveAuthority) return Transport::Status::Ok();
  if (state->configuration.source == Source::SwmfInitialThenEvolveLocal)
    return HandoffImportedState(
        state, state->epochS, snapshot->configuration_fingerprint());
  Transport::Status status = InitializeState(state);
  if (status.ok() && state->sourceChecksum.empty())
    state->sourceChecksum = snapshot->configuration_fingerprint();
  return status;
}

Transport::Status EnsureRuntimeLine(int fieldLineId) {
  const std::shared_ptr<const Background::BackgroundSnapshot> snapshot =
      Background::SnapshotStore::Instance().Current();
  if (!snapshot)
    return Transport::Status::Error(Transport::StatusCode::InvalidParticleState,
                                    "turbulence owner has no background snapshot");
  State sampled;
  const bool first = !gRuntimeStore.Contains(fieldLineId);
  Transport::Status status = ReadFieldLineView(fieldLineId, first, &sampled);
  if (!status.ok()) return status;
  if (first) {
    RuntimeContracts::TurbulenceLineRecord record;
    record.fieldLineId = fieldLineId;
    record.ownerRank = PIC::ThisThread;
    record.geometryGeneration = snapshot->field_line_generation();
    record.configurationFingerprint = snapshot->configuration_fingerprint();
    record.state = sampled;
    return gRuntimeStore.Install(record);
  }
  Turbulence::EnergyLedger remapLedger;
  status = gRuntimeStore.RefreshBackground(
      fieldLineId, snapshot->field_line_generation(),
      snapshot->configuration_fingerprint(), sampled.cells,
      RuntimeContracts::PendingRemapPolicy::RejectPending, &remapLedger);
  if (!status.ok()) return status;
  RuntimeContracts::TurbulenceLineRecord* record =
      gRuntimeStore.FindMutable(fieldLineId);
  record->state.epochS = Background::SimulationTimeSeconds();
  return Transport::Status::Ok();
}

double RelativisticKineticEnergyJ(double momentumKgMPerS, double massKg) {
  const double c = 299792458.0;
  const double rest = massKg * c * c;
  return std::hypot(momentumKgMPerS * c, rest) - rest;
}

std::uint64_t CouplingTransactionId(
    const Transport::PICAdapter::CouplingRecord& record) {
  // FNV-1a hashes only stable physical keys.  Rank, worker, queue address, and
  // particle-buffer handle are deliberately absent, making the transaction
  // ordering invariant under MPI/OpenMP decomposition.
  std::uint64_t hash = UINT64_C(1469598103934665603);
  const std::uint64_t values[] = {
      record.stableParticleId, record.snapshotGeneration, record.eventIndex,
      record.intervalIndex, static_cast<std::uint64_t>(record.fieldLineId + 1),
      static_cast<std::uint64_t>(record.species + 1),
      static_cast<std::uint64_t>(record.eventType)};
  for (std::size_t i = 0; i < sizeof(values) / sizeof(values[0]); ++i)
    for (unsigned byte = 0; byte < 8; ++byte) {
      hash ^= static_cast<unsigned char>((values[i] >> (8 * byte)) & 0xffU);
      hash *= UINT64_C(1099511628211);
    }
  return hash == 0 ? UINT64_C(1) : hash;
}

Transport::Status ApplyQueuedParticleTransactions(
    const std::vector<Transport::PICAdapter::CouplingRecord>& records) {
  std::vector<RuntimeContracts::CouplingTransaction> transactions;
  transactions.reserve(records.size());
  for (std::size_t i = 0; i < records.size(); ++i) {
    const Transport::PICAdapter::CouplingRecord& source = records[i];
    // Parker and deterministic focused intervals have no branch-resolved
    // wave-frame event. Their adiabatic momentum change belongs to background
    // work and must not be assigned to an arbitrary resonant wave branch.
    if (!source.pitchAngleResolved ||
        source.eventType !=
            Transport::PICAdapter::CouplingEventType::ScatteringEvent ||
        (source.resonantBranch != -1 && source.resonantBranch != 1)) continue;
    RuntimeContracts::TurbulenceLineRecord* line =
        gRuntimeStore.FindMutable(source.fieldLineId);
    if (!line)
      return Transport::Status::Error(Transport::StatusCode::OutOfDomain,
                                      "coupling record names no runtime wave owner");
    RuntimeContracts::CouplingTransaction target;
    target.transactionId = CouplingTransactionId(source);
    target.particleId = source.stableParticleId;
    target.generation = source.snapshotGeneration;
    target.event = source.eventIndex;
    target.interval = source.intervalIndex;
    target.fieldLineId = source.fieldLineId;
    const double midpoint = 0.5 * (source.startCoordinate +
                                   source.finishCoordinate);
    const double bounded = std::max(0.0, std::min(
        midpoint, static_cast<double>(line->state.cells.size()) - 1.0e-12));
    target.cell = static_cast<std::size_t>(std::floor(bounded));
    target.branch = source.resonantBranch;
    target.spectralBin = 0;
    const double massKg = PIC::MolecularData::GetMass(source.species);
    target.particleEnergyChangeJ = source.statisticalWeight *
        (RelativisticKineticEnergyJ(source.postMomentumKgMPerS, massKg) -
         RelativisticKineticEnergyJ(source.preMomentumKgMPerS, massKg));
    target.turbulenceIdentity = source.turbulenceStateIdentity;
    transactions.push_back(target);
  }
  const RuntimeContracts::CouplingBatchResult result =
      RuntimeContracts::ApplyCouplingTransactions(
          transactions, RuntimeContracts::InvalidTransactionPolicy::RejectBatch,
          &gRuntimeStore);
  return result.status;
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

  gLastStepLedger = ProductionLedger();
  if (!AlfvenTurbulence_Kolmogorov::ActiveFlag)
    return Transport::Status::Ok();
  if (!EvolvesLocally(ActiveConfiguration().source))
    return Transport::Status::Ok();

  // WP42 establishes every persistent owner before the particle transaction is
  // drained.  The current background generation is therefore the same one
  // recorded by mover-produced CouplingRecord objects.  Same-generation calls
  // refresh only geometry coefficients; wave energy is never reimported.
  for (int fieldLine = 0; fieldLine < PIC::FieldLine::nFieldLine; ++fieldLine) {
    const Transport::Status status = EnsureRuntimeLine(fieldLine);
    if (!status.ok()) return status;
  }

  // WP43 is an all-or-nothing commit.  Worker queues are deterministically
  // reduced without touching legacy G arrays, converted to branch-resolved
  // physical energy transactions, and staged in the authoritative store.
  std::vector<Transport::PICAdapter::CouplingRecord> particleRecords;
  Transport::Status status =
      Transport::PICAdapter::DrainWaveContributions(&particleRecords);
  if (!status.ok()) return status;
  status = ApplyQueuedParticleTransactions(particleRecords);
  if (!status.ok()) return status;

  // Shock geometry is evaluated after the persistent owners are synchronized.
  // The legacy geometry adapter now emits typed source records only; it never
  // mutates wave arrays.  Records are validated before any pending slot changes.
  ClearShockContributions();
  if (ActiveConfiguration().shockInjectionEnabled &&
      std::isfinite(shockRadiusBeforeM) &&
      std::isfinite(shockRadiusAfterM) &&
      shockRadiusBeforeM != shockRadiusAfterM) {
    SEP::ParticleSource::ShockWave::ShockTurbulenceEnergyInjection(
        shockRadiusBeforeM, shockRadiusAfterM, dtS);
  }
  for (std::size_t i = 0; i < gShockContributions.size(); ++i) {
    const ShockContribution& source = gShockContributions[i];
    RuntimeContracts::TurbulenceLineRecord* line =
        gRuntimeStore.FindMutable(source.fieldLine);
    if (!line || source.segment < 0 ||
        static_cast<std::size_t>(source.segment) >= line->state.cells.size())
      return Transport::Status::Error(
          Transport::StatusCode::OutOfDomain,
          "shock source record does not match the runtime turbulence owner");
  }
  for (std::size_t i = 0; i < gShockContributions.size(); ++i) {
    const ShockContribution& source = gShockContributions[i];
    Turbulence::CellState& cell = gRuntimeStore.FindMutable(
        source.fieldLine)->state.cells[static_cast<std::size_t>(source.segment)];
    cell.pendingShockPlusJ += source.plusJ;
    cell.pendingShockMinusJ += source.minusJ;
  }

  for (int fieldLine = 0; fieldLine < PIC::FieldLine::nFieldLine; ++fieldLine) {
    RuntimeContracts::TurbulenceLineRecord* owner =
        gRuntimeStore.FindMutable(fieldLine);
    if (!owner)
      return Transport::Status::Error(Transport::StatusCode::OutOfDomain,
                                      "runtime turbulence owner disappeared");
    const StepResult result = Turbulence::Advance(&owner->state, dtS);
    if (!result.status.ok()) return result.status;
    gLastStepLedger.boundaryExchangeJ +=
        result.ledger.innerBoundaryJ + result.ledger.outerBoundaryJ;
    gLastStepLedger.shockSourceJ += result.ledger.shockSourceJ;
    gLastStepLedger.particleExchangeJ += result.ledger.particleExchangeJ;
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
    status = ExportFieldLine(fieldLine, owner->state);
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

Transport::Status SerializeRuntimeStore(std::string* text) {
  return gRuntimeStore.Serialize(text);
}

Transport::Status RestoreRuntimeStore(const std::string& text) {
  return gRuntimeStore.Deserialize(text);
}

std::size_t RuntimeOwnerCount() { return gRuntimeStore.size(); }

void ResetRuntimeStoreForTests() { gRuntimeStore.Clear(); }

}  // namespace PICAdapter
}  // namespace Turbulence
}  // namespace SEP
