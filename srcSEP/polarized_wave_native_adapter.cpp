#include "polarized_wave_native_adapter.h"
#include "sep.h"
#include "util/sep_background_runtime.h"
#include "util/sep_common_header_path.h"
#include SRCSEP_SEP_COMMON_HEADER(sep_coefficient_physics.h)
#include "adapters/swcme1d_adapter.h"
#include "adapters/reduced_shock_background_adapter.h"

#include <algorithm>
#include <atomic>
#include <climits>
#include <cmath>
#include <cstring>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <tuple>
#include <utility>

namespace SEP { namespace Turbulence { namespace Polarized { namespace PICAdapter { namespace {
namespace W=AlfvenTurbulence_Kolmogorov;
namespace FL=PIC::FieldLine;
Transport::Status Error(const std::string& message) {
  return Transport::Status::Error(Transport::StatusCode::InvalidParticleState,message);
}
Transport::Status Unsupported(const std::string& message) {
  return Transport::Status::Error(Transport::StatusCode::UnsupportedConfiguration,message);
}
std::atomic<std::uint64_t> gNextOwnerId{1};
NativeStateOwner gOwner;
// Configuration identity binds both physical grid and all retained operator
// settings even while production evolution is disabled. A restore cannot
// silently acquire the same arrays under a different coupling/boundary policy.
std::string ConfigurationIdentity(const Configuration& config) {
  std::ostringstream text;text<<std::setprecision(17)<<"srcsep-polarized-native-configuration-v3;N="<<config.spectralBins
      <<";first-center="<<config.spectralKMinPerM<<";last-center="<<config.spectralKMaxPerM
      <<";initial-Lplus=0.5;initial-Lminus=0.5;source="<<SourceName(config.source)
      <<';'<<static_cast<int>(config.representation)<<';'<<static_cast<int>(config.coupling)
      <<';'<<static_cast<int>(config.innerBoundary.policy)<<';'<<config.innerBoundary.value
      <<';'<<static_cast<int>(config.outerBoundary.policy)<<';'<<config.outerBoundary.value
      <<';'<<config.advectionEnabled<<';'<<config.reflectionEnabled<<';'<<config.cascadeEnabled
      <<';'<<config.shockInjectionEnabled<<';'<<config.periodicBoundaries<<';'<<config.cflSafety
      <<';'<<static_cast<int>(config.operatorSplitting)<<';'<<config.operatorAccuracySafety
      <<';'<<config.maximumSourceFraction<<';'<<config.maximumCascadeFraction<<';'<<config.minimumSubstepS
      <<';'<<config.maximumSubsteps<<';'<<config.reflectionCoefficient<<';'<<config.cascadeCoefficient
      <<';'<<config.perpendicularCorrelationLengthM<<';'<<config.electronHeatingFraction
      <<';'<<config.diagnosticCadence<<';'<<config.conservationRelativeTolerance;
  return Background::FingerprintConfiguration(text.str());
}
bool SameConfiguration(const Configuration& a,const Configuration& b)noexcept {
  const auto values=[](const Configuration& c) {
    return std::tie(c.source,c.representation,c.coupling,c.innerBoundary.policy,c.innerBoundary.value,
        c.outerBoundary.policy,c.outerBoundary.value,c.advectionEnabled,c.reflectionEnabled,c.cascadeEnabled,
        c.shockInjectionEnabled,c.periodicBoundaries,c.cflSafety,c.operatorSplitting,c.operatorAccuracySafety,
        c.maximumSourceFraction,c.maximumCascadeFraction,c.minimumSubstepS,c.maximumSubsteps,
        c.reflectionCoefficient,c.cascadeCoefficient,c.perpendicularCorrelationLengthM,c.electronHeatingFraction,
        c.spectralKMinPerM,c.spectralKMaxPerM,c.spectralBins,c.diagnosticCadence,c.conservationRelativeTolerance);
  };
  return values(a)==values(b);
}
Transport::Status Collective(const Transport::Status& local) {
  // Each rank completes its local, allocation-capable preparation first.
  // No rank may enter the publication phase if any peer rejected its image.
  int success=local.ok()?1:0,all=0;
  MPI_Allreduce(&success,&all,1,MPI_INT,MPI_MIN,MPI_GLOBAL_COMMUNICATOR);
  return all?Transport::Status::Ok():local.ok()?Error("another MPI rank rejected polarized initialization before publication"):local;
}
bool SameIdentity(const std::string& local) {
  char reference[96]={0};
  int fits=local.size()<sizeof(reference)?1:0,all=0;
  MPI_Allreduce(&fits,&all,1,MPI_INT,MPI_MIN,MPI_GLOBAL_COMMUNICATOR);if(!all)return false;
  if(PIC::ThisThread==0)std::memcpy(reference,local.data(),local.size());
  MPI_Bcast(reference,sizeof(reference),MPI_CHAR,0,MPI_GLOBAL_COMMUNICATOR);
  int equal=std::strlen(reference)==local.size()&&std::memcmp(reference,local.data(),local.size())==0?1:0;
  MPI_Allreduce(&equal,&all,1,MPI_INT,MPI_MIN,MPI_GLOBAL_COMMUNICATOR);return all!=0;
}
struct DiagnosticRow {
  double* integrated=nullptr;double* density=nullptr;double* spectrum=nullptr;
  std::uint64_t line=0,segment=0;
  std::vector<double> before,after;
};
Transport::Status StageDiagnostics(const Snapshot& state,std::vector<DiagnosticRow>* rows) {
  // Native E+/E-, W+/W-/cross-helicity and two-branch spectral buffers become
  // DERIVED readback only. Keep their old bytes for rollback; four independent
  // immutable arrays in Snapshot are the sole physical energy authority.
  // Dynamic2*N dimensions are checked before touching a host buffer, so a
  // legacy NK mismatch is an explicit capability failure rather than a copy.
  const auto bins=state.grid.centersPerM.size();
  if(W::CellIntegratedWaveEnergy.length!=2||W::WaveEnergyDensity.length!=3||
      W::WaveNumberResolved::SpectralWaveEnergy.length!=static_cast<int>(2*bins)||
      !W::WaveNumberResolved::SpectralWaveEnergy.is_active())
    return Unsupported("native polarized diagnostic layout differs from the frozen physical grid; no resize/remap allowed");
  std::vector<DiagnosticRow> staged;staged.reserve(state.cells.size());
  for(const auto& cell:state.cells) {
    if(cell.fieldLine>=static_cast<std::uint64_t>(FL::nFieldLine)||
        cell.segment>=static_cast<std::uint64_t>(FL::FieldLinesAll[cell.fieldLine].GetTotalSegmentNumber()))
      return Error("polarized native cell namespace differs from fixed field-line geometry");
    auto* segment=FL::FieldLinesAll[cell.fieldLine].GetSegment(static_cast<int>(cell.segment));
    if(!segment)return Error("polarized native segment is absent");
    DiagnosticRow row;row.line=cell.fieldLine;row.segment=cell.segment;
    row.integrated=segment->GetDatum_ptr(W::CellIntegratedWaveEnergy);
    row.density=segment->GetDatum_ptr(W::WaveEnergyDensity);
    row.spectrum=segment->GetDatum_ptr(W::WaveNumberResolved::SpectralWaveEnergy);
    if(!row.integrated||!row.density||!row.spectrum)return Error("polarized native diagnostic buffer is absent");
    const auto total=DeriveBranches(cell);if(!total.status.ok())return total.status;
    row.after={total.ePlusJ,total.eMinusJ,total.ePlusJ/cell.volumeM3,total.eMinusJ/cell.volumeM3,
        total.ePlusJ+total.eMinusJ==0?0:(total.ePlusJ-total.eMinusJ)/(total.ePlusJ+total.eMinusJ)};
    for(int branch=0;branch<2;++branch)for(std::size_t j=0;j<bins;++j)
      row.after.push_back(cell.magneticEnergyJ[2*branch][j]+cell.magneticEnergyJ[2*branch+1][j]);
    row.before.insert(row.before.end(),row.integrated,row.integrated+2);
    row.before.insert(row.before.end(),row.density,row.density+3);
    row.before.insert(row.before.end(),row.spectrum,row.spectrum+2*bins);
    for(double value:row.after)if(!std::isfinite(value))return Error("derived polarized diagnostic overflowed");
    staged.push_back(std::move(row));
  }
  rows->swap(staged);return Transport::Status::Ok();
}
void WriteDiagnostics(const std::vector<DiagnosticRow>& rows,bool restore)noexcept {
  for(const auto& row:rows) {
    const auto& values=restore?row.before:row.after;
    std::copy(values.begin(),values.begin()+2,row.integrated);
    std::copy(values.begin()+2,values.begin()+5,row.density);
    std::copy(values.begin()+5,values.end(),row.spectrum);
  }
}
bool Readback(const std::vector<DiagnosticRow>& rows) {
  for(const auto& row:rows)
    if(!std::equal(row.after.begin(),row.after.begin()+2,row.integrated)||
       !std::equal(row.after.begin()+2,row.after.begin()+5,row.density)||
       !std::equal(row.after.begin()+5,row.after.end(),row.spectrum))return false;
  return true;
}
NativeDiagnostics DeriveDiagnostics(const Snapshot& snapshot) {
  NativeDiagnostics result;result.cells=snapshot.cells.size();result.bins=snapshot.grid.centersPerM.size();
  for(const auto& cell:snapshot.cells) {
    const auto branches=DeriveBranches(cell);
    if(branches.ePlusJ>0)result.maximumBranchRelativeRoundoff=std::max(result.maximumBranchRelativeRoundoff,
        std::fabs(cell.initialRoundoffPlusJ)/branches.ePlusJ);
    if(branches.eMinusJ>0)result.maximumBranchRelativeRoundoff=std::max(result.maximumBranchRelativeRoundoff,
        std::fabs(cell.initialRoundoffMinusJ)/branches.eMinusJ);
    result.maximumAbsoluteInitialRoundoffJ=std::max(result.maximumAbsoluteInitialRoundoffJ,
        std::max(std::fabs(cell.initialRoundoffPlusJ),std::fabs(cell.initialRoundoffMinusJ)));
  }
  return result;
}
bool PointersStillCurrent(const std::vector<DiagnosticRow>& rows,std::size_t bins) {
  // Plans retain raw datum pointers only inside the enclosing quiescent
  // transaction. Recheck both shape and storage identity to reject a relocated
  // or reconfigured field-line buffer before any commit writes.
  if(W::CellIntegratedWaveEnergy.length!=2||W::WaveEnergyDensity.length!=3||
      W::WaveNumberResolved::SpectralWaveEnergy.length!=static_cast<int>(2*bins))return false;
  for(const auto& row:rows) {
    if(row.line>=static_cast<std::uint64_t>(FL::nFieldLine)||
        row.segment>=static_cast<std::uint64_t>(FL::FieldLinesAll[row.line].GetTotalSegmentNumber()))return false;
    auto* segment=FL::FieldLinesAll[row.line].GetSegment(static_cast<int>(row.segment));
    if(!segment||segment->GetDatum_ptr(W::CellIntegratedWaveEnergy)!=row.integrated||
        segment->GetDatum_ptr(W::WaveEnergyDensity)!=row.density||
        segment->GetDatum_ptr(W::WaveNumberResolved::SpectralWaveEnergy)!=row.spectrum)return false;
  }
  return true;
}
} // namespace

NativeStateOwner::NativeStateOwner():lifetimeId_(gNextOwnerId.fetch_add(1)) {}

std::string NativeBackgroundIdentity(const Background::BackgroundSnapshot& bg) {
  std::ostringstream text;text<<std::setprecision(17)<<Background::ProviderName(bg.provider())
      <<":ownership="<<Background::OwnershipName(bg.ownership())
      <<":generation="<<bg.field_line_generation()<<":epoch_s="<<bg.epoch_seconds()
      <<":valid_from="<<bg.valid_from_seconds()<<":valid_until="<<bg.valid_until_seconds()
      <<":previous_physical_epoch="<<bg.previous_physical_epoch_seconds()
      <<":current_physical_epoch="<<bg.current_physical_epoch_seconds()
      <<":configuration="<<bg.configuration_fingerprint();return text.str();
}
std::string NativePhysicalImageIdentity(const std::vector<NativePhysicalCell>& cells) {
  // Bind full endpoint primitives as well as midpoint samples. Opposite
  // endpoint perturbations can cancel at the midpoint yet change the mover's
  // interpolated background, so a midpoint-only hash would be insufficient.
  std::ostringstream text;text<<std::setprecision(17)<<"srcsep-native-wave-input-image-v2;cells="<<cells.size();
  for(const auto& cell:cells) {
    text<<';'<<cell.fieldLine<<','<<cell.segment<<','<<cell.volumeM3<<','<<cell.meanMagneticFieldT<<','<<cell.massDensityKgM3;
    for(const auto* vector:{&cell.positionM,&cell.tangent,&cell.magneticFieldT,&cell.plasmaVelocityMPerS})
      for(double value:*vector)text<<','<<value;
    for(const auto* endpoint:{&cell.begin,&cell.end}) {
      text<<','<<endpoint->numberDensityM3;
      for(const auto* vector:{&endpoint->positionM,&endpoint->magneticFieldT,&endpoint->plasmaVelocityMPerS})
        for(double value:*vector)text<<','<<value;
    }
  }
  return Background::FingerprintConfiguration(text.str());
}
Transport::Status ReadNativePhysicalCells(const Background::BackgroundSnapshot& bg,
    std::vector<NativePhysicalCell>* output) {
  // Composition is provider-owned: rho=n*meanMass only for the explicit
  // analytic interface. SWCME/reduced providers supply their own rho/n at the
  // exact published epoch. No guessed thermal/composition closure is added.
  if(!output)return Error("native physical image destination is absent");
  if(bg.provider()==Background::Provider::Swmf||bg.provider()==Background::Provider::LocalEvolution)
    return Unsupported("four-spectrum initialization has no authoritative SWMF branch/polarization handoff contract");
  std::vector<NativePhysicalCell> cells;
  for(int line=0;line<FL::nFieldLine;++line)for(int i=0;i<FL::FieldLinesAll[line].GetTotalSegmentNumber();++i) {
    auto* segment=FL::FieldLinesAll[line].GetSegment(i);if(!segment)return Error("native physical segment is absent");
    double field[3],position[3],density=0;
    FL::FieldLinesAll[line].GetMagneticField(field,i+0.5);
    segment->GetCartesian(position,0.5);segment->GetPlasmaDensity(0.5,density);
    NativePhysicalCell cell;cell.fieldLine=line;cell.segment=i;
    cell.volumeM3=FieldLine::FluxTubeGeometry::SegmentVolumeM3(segment,line);
    cell.meanMagneticFieldT=Vector3D::Length(field);
    double tangent[3],velocity0[3],velocity1[3];segment->GetDir(tangent);
    segment->GetBegin()->GetPlasmaVelocity(velocity0);segment->GetEnd()->GetPlasmaVelocity(velocity1);
    for(int d=0;d<3;++d) {
      cell.positionM[d]=position[d];cell.tangent[d]=tangent[d];cell.magneticFieldT[d]=field[d];
      cell.plasmaVelocityMPerS[d]=0.5*(velocity0[d]+velocity1[d]);
    }
    auto endpoint=[](FL::cFieldLineVertex* vertex,NativeEndpoint* state) {
      vertex->GetX(state->positionM.data());vertex->GetPlasmaVelocity(state->plasmaVelocityMPerS.data());
      const double* magnetic=vertex->GetDatum_ptr(FL::DatumAtVertexMagneticField);
      std::copy(magnetic,magnetic+3,state->magneticFieldT.begin());vertex->GetPlasmaDensity(state->numberDensityM3);
    };
    endpoint(segment->GetBegin(),&cell.begin);endpoint(segment->GetEnd(),&cell.end);
    cell.massDensityKgM3=density*PIC::CPLR::SWMF::MeanPlasmaAtomicMass;
    if(bg.provider()==Background::Provider::Swcme) {
      const auto query=SW1DAdapter::QueryAtRadius(Vector3D::Length(position),false);
      if(query.code!=SW1DAdapter::StatusCode::Success||query.epoch_seconds!=bg.epoch_seconds()||
          !query.sample.magnetic_and_composition_available())return Error("native four-state SWCME composition/epoch is unavailable");
      cell.massDensityKgM3=density*query.sample.mass_density_kg_m3()/query.sample.number_density_m3();
    } else if(bg.provider()==Background::Provider::ReducedShock) {
      ReducedShock::AmbientSample sample;std::string error;
      if(ReducedShock::Metadata().epochS!=bg.epoch_seconds()||
          ReducedShock::Metadata().generation!=bg.field_line_generation()||
          !ReducedShock::EvaluateAmbient({{position[0],position[1],position[2]}},&sample,&error))
        return Error("native four-state reduced-background composition/epoch is unavailable: "+error);
      cell.massDensityKgM3=density*sample.massDensityKgM3/sample.numberDensityM3;
    }
    if(!std::isfinite(cell.volumeM3)||cell.volumeM3<=0||!std::isfinite(cell.meanMagneticFieldT)||
        cell.meanMagneticFieldT<=0||!std::isfinite(cell.massDensityKgM3)||cell.massDensityKgM3<=0)
      return Error("native four-state background requires positive finite V/B/rho");
    cells.push_back(cell);
  }
  if(cells.empty())return Error("native four-state physical image has no segments");
  output->swap(cells);return Transport::Status::Ok();
}
Transport::Status ValidatePhysicalImage(const Snapshot& waves,const std::vector<NativePhysicalCell>& cells) {
  if(waves.cells.size()!=cells.size())return Error("four-wave/native physical image cell counts differ");
  for(std::size_t i=0;i<cells.size();++i) {
    const auto& a=waves.cells[i];const auto& b=cells[i];
    if(a.fieldLine!=b.fieldLine||a.segment!=b.segment||a.volumeM3!=b.volumeM3||
        a.meanMagneticFieldT!=b.meanMagneticFieldT||a.massDensityKgM3!=b.massDensityKgM3)
      return Error("four-wave V/B/rho or cell identity differs from independently read native/background image");
    for(const auto* vector:{&b.positionM,&b.tangent,&b.magneticFieldT,&b.plasmaVelocityMPerS})
      for(double value:*vector)if(!std::isfinite(value))return Error("native/background physical vector is nonfinite");
    const double norm=std::sqrt(b.magneticFieldT[0]*b.magneticFieldT[0]+
        b.magneticFieldT[1]*b.magneticFieldT[1]+b.magneticFieldT[2]*b.magneticFieldT[2]);
    if(norm!=b.meanMagneticFieldT)return Error("native/background magnetic vector does not match its declared magnitude");
    for(const auto* endpoint:{&b.begin,&b.end}) {
      if(!std::isfinite(endpoint->numberDensityM3)||endpoint->numberDensityM3<=0)
        return Error("native/background endpoint number density is invalid");
      for(const auto* vector:{&endpoint->positionM,&endpoint->magneticFieldT,&endpoint->plasmaVelocityMPerS})
        for(double value:*vector)if(!std::isfinite(value))return Error("native/background endpoint vector is nonfinite");
    }
  }
  return Transport::Status::Ok();
}

Transport::Status ConfigureNativeDiagnosticsBeforeFieldLineInit() {
  // Datum offsets become fixed when FieldLine::Init allocates native storage.
  // Resolve the configured physical N before then; never resize a live datum
  // or reinterpret its existing2*NK bytes as a new grid.
  const auto& config=ActiveConfiguration();
  const auto grid=ResolveCenteredLogGrid(config.spectralBins,config.spectralKMinPerM,config.spectralKMaxPerM);
  if(!grid.status.ok())return grid.status;
  auto& spectral=W::WaveNumberResolved::SpectralWaveEnergy;
  auto& rates=W::WaveNumberResolved::SpectralWaveEnergyExchangeRate;
  const std::size_t columns=2*config.spectralBins,rateColumns=columns*W::WaveNumberResolved::nSpectralEnergyExchangeRateComponents;
  if(columns>INT_MAX||rateColumns>INT_MAX)return Error("native diagnostic spectral allocation exceeds integer buffer dimensions");
  if(spectral.is_active()||rates.is_active()) {
    if(spectral.length!=static_cast<int>(columns)||rates.length!=static_cast<int>(rateColumns))
      return Unsupported("live native diagnostic spectrum cannot be resized or remapped");
    return Transport::Status::Ok();
  }
  spectral.length=static_cast<int>(columns);rates.length=static_cast<int>(rateColumns);
  return Transport::Status::Ok();
}

Transport::Status NativeStateOwner::InitializeFromIntegrated(const Configuration& config) {
  std::lock_guard<std::mutex> lock(mutex_);
  Transport::Status local=ValidateConfiguration(config);SnapshotResult prepared;
  std::vector<DiagnosticRow> rows;std::vector<NativePhysicalCell> physical;
  std::shared_ptr<const Background::BackgroundSnapshot> background;
  NativeDiagnostics diagnostic;
  try {
    if(current_)local=Unsupported("native four-state authority already exists; reinitialization/remap requires a physical transaction");
    if(local.ok()&&(config.representation!=Representation::Spectral||!EvolvesLocally(config.source)))
      local=Unsupported("native four-state initialization requires a declared local spectral source");
    background=Background::SnapshotStore::Instance().Current();
    if(local.ok()&&!background)local=Error("native four-state initialization requires a published background realization");
    if(local.ok()&&!background->Covers(Background::SimulationTimeSeconds()))
      local=Error("native four-state initialization background does not cover the authoritative simulation clock");
    if(local.ok()) {
      Background::SnapshotStore::Instance().AssertProviderMayWrite(background->provider());
      local=ReadNativePhysicalCells(*background,&physical);
    }
    const auto grid=ResolveCenteredLogGrid(config.spectralBins,config.spectralKMinPerM,config.spectralKMaxPerM);
    if(local.ok()&&!grid.status.ok())local=grid.status;
    if(local.ok()) {
      Initialization input;input.grid.centersPerM=grid.grid.centersPerM;input.grid.logWidth=grid.grid.logWidth;
      input.plus.known=true;input.minus.known=true;input.plus.leftFraction=0.5;input.minus.leftFraction=0.5;
      input.plus.provenance=input.minus.provenance="model-owner decision: equal INITIAL L/R in each branch; subsequent four spectra independent";
      input.powerLawIndex=5.0/3.0;input.vacuumPermeabilityHPerM=Transport::CoefficientPhysics::MagneticEnergyPermeabilitySI;
      input.epochS=background->epoch_seconds();input.validUntilS=background->valid_until_seconds();
      input.geometryGeneration=background->field_line_generation();
      input.backgroundIdentity=NativeBackgroundIdentity(*background)+":native-physical="+NativePhysicalImageIdentity(physical);
      input.configurationFingerprint=ConfigurationIdentity(config);
      input.initializerProvenance="existing InitializeWaveEnergyFromPhysicalParameters amplitude/profile/Eplus-Eminus; Kolmogorov bin integral; grid-v2 centers";
      for(const auto& source:physical) {
        InitialCell cell;cell.fieldLine=source.fieldLine;cell.segment=source.segment;cell.volumeM3=source.volumeM3;
        cell.meanMagneticFieldT=source.meanMagneticFieldT;cell.massDensityKgM3=source.massDensityKgM3;
        const auto* energy=FL::FieldLinesAll[cell.fieldLine].GetSegment(static_cast<int>(cell.segment))->GetDatum_ptr(W::CellIntegratedWaveEnergy);
        if(!energy){local=Error("existing physical integrated initializer output is absent");break;}
        cell.ePlusJ=energy[0];cell.eMinusJ=energy[1];input.cells.push_back(cell);
      }
      if(local.ok())prepared=PrepareInitialization(input);
      if(local.ok()&&!prepared.status.ok())local=prepared.status;
      if(local.ok())local=ValidatePhysicalImage(*prepared.snapshot,physical);
      if(local.ok())local=StageDiagnostics(*prepared.snapshot,&rows);
      if(local.ok())diagnostic=DeriveDiagnostics(*prepared.snapshot);
    }
  } catch(const std::exception& error) {local=Error(std::string("native four-state preparation failed: ")+error.what());}
  auto agreed=Collective(local);if(!agreed.ok())return agreed;
  if(!SameIdentity(prepared.snapshot->identity))return Error("native four-state physical cells/bin energies differ across MPI ranks; nothing published");
  diagnostic.mpiIdentityAgrees=true;
  // All allocations, primitive checks, shape checks and collective identity
  // comparison precede the first diagnostic write. Roll back on failed readback.
  WriteDiagnostics(rows,false);
  local=Readback(rows)?Transport::Status::Ok():Error("native four-state diagnostic readback differs from staged bin moments");
  agreed=Collective(local);
  if(!agreed.ok()){WriteDiagnostics(rows,true);return agreed;}
  diagnostic.readbackAgrees=true;
  current_=std::move(prepared.snapshot);diagnostics_=diagnostic;
  MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);return Transport::Status::Ok();
}

struct NativeRestorePlan::Implementation {
  // This is an internal transaction component, not a second mutable owner.
  // Shared pointers pin both old and candidate lifetimes; lifetimeId also
  // rejects address reuse if an old owner is destroyed and reconstructed.
  const NativeStateOwner* owner=nullptr;
  std::uint64_t ownerLifetime=0;
  std::shared_ptr<const Snapshot> expectedCurrent,candidate;
  std::shared_ptr<const Background::BackgroundSnapshot> expectedLiveBackground;
  Configuration configuration;
  std::string configurationIdentity,livePhysicalIdentity;
  std::vector<DiagnosticRow> diagnostics;
  NativeDiagnostics derivedDiagnostics;
  bool locallyReady=false,collectivelyReady=false;
};
NativeRestorePlan::NativeRestorePlan()=default;NativeRestorePlan::~NativeRestorePlan()=default;
NativeRestorePlan::NativeRestorePlan(NativeRestorePlan&&)noexcept=default;
NativeRestorePlan& NativeRestorePlan::operator=(NativeRestorePlan&&)noexcept=default;
Transport::Status NativeStateOwner::StageRestore(const Snapshot& decoded,
    const Background::BackgroundSnapshot& bg,const std::vector<NativePhysicalCell>& physical,
    NativeRestorePlan* plan)const {
  if(!plan)return Error("native four-state restore-plan destination is absent");
  plan->implementation_.reset();
  try {
  // Pure snapshots permit any explicitly declared positive permeability for
  // controlled physics tests. The native SI model has one canonical constant;
  // internally consistent but noncanonical decoded values must still reject.
  if(decoded.vacuumPermeabilityHPerM!=Transport::CoefficientPhysics::MagneticEnergyPermeabilitySI)
    return Unsupported("restored native wave magnetic permeability differs from the canonical SI energy contract");
  const auto prepared=PrepareSnapshotForRestore(decoded);if(!prepared.status.ok())return prepared.status;
  if(decoded.backgroundIdentity!=NativeBackgroundIdentity(bg)+":native-physical="+NativePhysicalImageIdentity(physical)||
      decoded.geometryGeneration!=bg.field_line_generation()||
      decoded.epochS!=bg.epoch_seconds()||decoded.validUntilS!=bg.valid_until_seconds())
    return Error("restored four-wave epoch/provider/geometry differs from decoded background image");
  auto status=ValidatePhysicalImage(decoded,physical);if(!status.ok())return status;
  if(decoded.initialPlus.leftFraction!=0.5||decoded.initialMinus.leftFraction!=0.5)
    return Unsupported("restored native initial-polarization policy differs from explicit model-owner equal INITIAL L/R decision");
  const auto& config=ActiveConfiguration();
  const auto grid=ResolveCenteredLogGrid(config.spectralBins,config.spectralKMinPerM,config.spectralKMaxPerM);
  if(!grid.status.ok())return grid.status;
  if(decoded.grid.centersPerM!=grid.grid.centersPerM||decoded.grid.edgesPerM!=grid.grid.edgesPerM||
      decoded.configurationFingerprint!=ConfigurationIdentity(config))
    return Unsupported("restored four-state grid/configuration differs from configured native authority; no remap performed");
  auto staged=std::unique_ptr<NativeRestorePlan::Implementation>(new NativeRestorePlan::Implementation);
  staged->owner=this;staged->ownerLifetime=lifetimeId_;staged->candidate=prepared.snapshot;
  staged->configuration=config;
  staged->configurationIdentity=ConfigurationIdentity(config);
  staged->derivedDiagnostics=DeriveDiagnostics(*prepared.snapshot);
  staged->expectedLiveBackground=Background::SnapshotStore::Instance().Current();
  if(!staged->expectedLiveBackground)return Error("native four-state restore requires a live background realization");
  Background::SnapshotStore::Instance().AssertProviderMayWrite(staged->expectedLiveBackground->provider());
  std::vector<NativePhysicalCell> livePhysical;
  status=ReadNativePhysicalCells(*staged->expectedLiveBackground,&livePhysical);if(!status.ok())return status;
  staged->livePhysicalIdentity=NativePhysicalImageIdentity(livePhysical);
  status=StageDiagnostics(decoded,&staged->diagnostics);if(!status.ok())return status;
  {std::lock_guard<std::mutex> lock(mutex_);staged->expectedCurrent=current_;}
  plan->implementation_=std::move(staged);return Transport::Status::Ok();
  } catch(const std::exception& error) {
    return Error(std::string("native four-state restore preparation failed before publication: ")+error.what());
  }
}
Transport::Status NativeStateOwner::ValidateReadyForCommit(const NativeRestorePlan& plan)const {
  // Every failed readiness check revokes prior approval. Full endpoint images
  // and provider read-phase state are rechecked while allocations/queries are
  // still allowed, before the outer collective and allocation-free commit.
  if(!plan.implementation_)return Error("native four-state restore plan is absent");
  auto& staged=*plan.implementation_;staged.locallyReady=false;staged.collectivelyReady=false;
  try {
  std::lock_guard<std::mutex> lock(mutex_);
  if(staged.owner!=this||staged.ownerLifetime!=lifetimeId_||staged.expectedCurrent!=current_)
    return Error("native four-state restore plan is stale or belongs to another owner lifetime");
  if(staged.configurationIdentity!=ConfigurationIdentity(ActiveConfiguration()))
    return Unsupported("native four-state configuration changed after restore staging");
  const auto live=Background::SnapshotStore::Instance().Current();
  if(live!=staged.expectedLiveBackground)return Error("native background identity changed after wave restore staging");
  Background::SnapshotStore::Instance().AssertProviderMayWrite(live->provider());
  if(!PointersStillCurrent(staged.diagnostics,staged.candidate->grid.centersPerM.size()))
    return Error("native diagnostic layout or storage changed after wave restore staging");
  std::vector<NativePhysicalCell> physical;
  const auto status=ReadNativePhysicalCells(*live,&physical);if(!status.ok())return status;
  if(NativePhysicalImageIdentity(physical)!=staged.livePhysicalIdentity)
    return Error("native wave-input endpoint primitives changed after restore staging");
  staged.locallyReady=true;return Transport::Status::Ok();
  } catch(const std::exception& error) {
    return Error(std::string("native four-state restore readiness rejected before publication: ")+error.what());
  }
}
Transport::Status NativeStateOwner::ConfirmCollectiveReady(const NativeRestorePlan& plan)const {
  const auto agreed=Collective(ValidateReadyForCommit(plan));if(!agreed.ok())return agreed;
  if(!SameIdentity(plan.implementation_->candidate->identity))
    return Error("restored four-state identities differ across MPI ranks; nothing published");
  plan.implementation_->derivedDiagnostics.mpiIdentityAgrees=true;
  plan.implementation_->collectivelyReady=true;return Transport::Status::Ok();
}
bool NativeStateOwner::CommitRestore(NativeRestorePlan* plan)noexcept {
  if(!plan||!plan->implementation_)return false;
  std::lock_guard<std::mutex> lock(mutex_);
  auto& prepared=*plan->implementation_;
  if(prepared.owner!=this||prepared.ownerLifetime!=lifetimeId_||prepared.expectedCurrent!=current_||
      !prepared.locallyReady||!prepared.collectivelyReady||
      !SameConfiguration(prepared.configuration,ActiveConfiguration())||
      !PointersStillCurrent(prepared.diagnostics,prepared.candidate->grid.centersPerM.size()))return false;
  try {
    const auto live=Background::SnapshotStore::Instance().Current();
    if(live!=prepared.expectedLiveBackground)return false;
    Background::SnapshotStore::Instance().AssertProviderMayWrite(live->provider());
  } catch(...) {return false;}
  // The enclosing transaction holds the documented quiescent interval after
  // ConfirmCollectiveReady. No allocation, physics query or collective occurs
  // here; all such checks belong to the preceding all-rank readiness gate.
  WriteDiagnostics(prepared.diagnostics,false);
  if(!Readback(prepared.diagnostics)){WriteDiagnostics(prepared.diagnostics,true);return false;}
  prepared.derivedDiagnostics.readbackAgrees=true;
  current_=std::move(prepared.candidate);diagnostics_=prepared.derivedDiagnostics;
  plan->implementation_.reset();return true;
}
std::shared_ptr<const Snapshot> NativeStateOwner::Acquire()const {std::lock_guard<std::mutex> lock(mutex_);return current_;}
NativeDiagnostics NativeStateOwner::Diagnostics()const {std::lock_guard<std::mutex> lock(mutex_);return diagnostics_;}
Transport::Status InitializeNativePolarizedState(){return gOwner.InitializeFromIntegrated(ActiveConfiguration());}
bool HasNativePolarizedState(){return static_cast<bool>(gOwner.Acquire());}
std::shared_ptr<const Snapshot> AcquireNativePolarizedState(){return gOwner.Acquire();}
NativeStateOwner& NativePolarizedStateOwner(){return gOwner;}
Transport::Status LegacyMutationPreflight(const char* operation) {
  return LegacyMutationPreflight(operation,gOwner);
}
Transport::Status LegacyMutationPreflight(const char* operation,const NativeStateOwner& authority) {
  if(!operation||!*operation)return Error("legacy mutation capability requires an operation context");
  return authority.Acquire()?Unsupported(std::string("legacy mutation cannot overwrite four-spectrum authority: ")+operation):
      Transport::Status::Ok();
}
void RequireLegacyMutationAllowed(const char* operation) {
  const auto status=LegacyMutationPreflight(operation);
  if(!status.ok())throw std::logic_error(status.message);
}
Transport::Status CopyLegacyNativeSpectrum(const Configuration& config,const double* source,
    std::size_t length,std::vector<double>* output) {
  const auto capability=LegacyFixedSpectralPreflight(config,"two-branch native import");if(!capability.ok())return capability;
  const auto grid=ResolveCenteredLogGrid(config.spectralBins,config.spectralKMinPerM,config.spectralKMaxPerM);
  if(!grid.status.ok())return grid.status;
  CenteredLogGrid requested;requested.centersPerM=grid.grid.centersPerM;requested.logWidth=grid.grid.logWidth;
  CenteredLogGrid native;native.logWidth=W::IsotropicSEP::DLNK;
  for(int j=0;j<W::IsotropicSEP::NK;++j)native.centersPerM.push_back(W::IsotropicSEP::K_MIN*std::exp(j*native.logWidth));
  return CopyLegacyBranchSpectrum(requested,native,source,length,output);
}
Transport::Status LegacyFixedSpectralPreflight(const Configuration& config,const char* operation) {
  const auto authority=LegacyMutationPreflight(operation);if(!authority.ok())return authority;
  if(config.spectralBins!=static_cast<std::size_t>(W::IsotropicSEP::NK)||
      config.spectralKMinPerM!=W::IsotropicSEP::K_MIN||config.spectralKMaxPerM!=W::IsotropicSEP::K_MAX||
      (W::WaveNumberResolved::SpectralWaveEnergy.is_active()&&
       W::WaveNumberResolved::SpectralWaveEnergy.length!=2*W::IsotropicSEP::NK))
    return Unsupported(std::string("legacy fixed-NK spectral operation cannot copy/remap the configured physical grid: ")+operation);
  return Transport::Status::Ok();
}
} } } }
