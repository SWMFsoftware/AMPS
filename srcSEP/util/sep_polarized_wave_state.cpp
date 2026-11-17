#include "sep_polarized_wave_state.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>
#include <set>
#include <sstream>
#include <utility>

namespace SEP { namespace Turbulence { namespace Polarized { namespace {
Transport::Status Error(const std::string& message) {
  return Transport::Status::Error(Transport::StatusCode::InvalidArgument,message);
}
bool Positive(double x){return std::isfinite(x)&&x>0;}
bool Nonnegative(double x){return std::isfinite(x)&&x>=0;}
class Fingerprint {
 public:
  void Add(std::uint64_t value) {
    for(int byte=0;byte<8;++byte){hash_^=value&255;hash_*=UINT64_C(1099511628211);value>>=8;}
  }
  void Add(double value) {
    static_assert(sizeof(double)==sizeof(std::uint64_t),"wave identity requires binary64 doubles");
    static_assert(std::numeric_limits<double>::is_iec559,"wave identity requires IEEE floating point");
    std::uint64_t bits=0;std::memcpy(&bits,&value,sizeof(bits));Add(bits);
  }
  void Add(const std::string& text) {
    Add(static_cast<std::uint64_t>(text.size()));
    for(unsigned char byte:text){hash_^=byte;hash_*=UINT64_C(1099511628211);}
  }
  std::string Text()const {std::ostringstream out;out<<"polarized-v2:"<<std::hex<<hash_;return out.str();}
 private:std::uint64_t hash_=UINT64_C(1469598103934665603);
};
long double Sum(const std::vector<double>& values) {
  long double result=0,correction=0;
  for(double value:values) {
    const long double adjusted=static_cast<long double>(value)-correction;
    const long double updated=result+adjusted;
    correction=(updated-result)-adjusted;result=updated;
  }
  return result;
}
Transport::Status CheckPolarization(const BranchPolarization& p,const char* branch) {
  if(!p.known)return Transport::Status::Error(Transport::StatusCode::UnsupportedConfiguration,
      std::string("existing initializer supplies only E+/E-; explicit left/right polarization is missing for ")+branch);
  if(!std::isfinite(p.leftFraction)||p.leftFraction<0||p.leftFraction>1||p.provenance.empty())
    return Error(std::string("polarization requires an explicit fraction in [0,1] and provenance for ")+branch);
  return Transport::Status::Ok();
}
Transport::Status CheckDomainMetadata(const Transport::ResonantCollision::DomainContract& domain) {
  const double numbers[]={domain.maximumDeltaBOverB,domain.maximumKThermalIonGyroradius,
      domain.thermalIonGyroradiusM,domain.maximumAlfvenSpeedOverC,domain.ionInertialLengthM,
      domain.maximumKIonInertialLength,domain.ionCyclotronFrequencyRadPerS,
      domain.maximumWaveFrequencyOverIonCyclotron};
  for(double number:numbers)if(!Nonnegative(number))
    return Error("polarized domain metadata must be finite and nonnegative; zero denotes unavailable");
  using Authority=Transport::ResonantCollision::DomainAuthority;
  return domain.authority==Authority::Unspecified||domain.authority==Authority::ExplicitDeclaredLimits||
      domain.authority==Authority::ValidationStudy?Transport::Status::Ok():Error("polarized domain authority is unknown");
}
std::string SnapshotIdentity(const Snapshot& state) {
  Fingerprint f;f.Add(state.schema);f.Add(state.revision);f.Add(state.geometryGeneration);
  f.Add(state.epochS);f.Add(state.validUntilS);f.Add(state.vacuumPermeabilityHPerM);
  f.Add(state.initialPowerLawIndex);f.Add(state.initialPlus.leftFraction);f.Add(state.initialMinus.leftFraction);
  f.Add(state.initialPlus.provenance);f.Add(state.initialMinus.provenance);
  f.Add(state.backgroundIdentity);f.Add(state.configurationFingerprint);f.Add(state.initializerProvenance);
  f.Add(static_cast<std::uint64_t>(state.grid.edgesPerM.size()));f.Add(state.grid.logWidth);
  for(double center:state.grid.centersPerM)f.Add(center);
  for(double edge:state.grid.edgesPerM)f.Add(edge);
  f.Add(static_cast<std::uint64_t>(state.cells.size()));
  for(const auto& cell:state.cells) {
    f.Add(cell.fieldLine);f.Add(cell.segment);f.Add(cell.volumeM3);f.Add(cell.meanMagneticFieldT);
    f.Add(cell.massDensityKgM3);f.Add(cell.alfvenSpeedMPerS);
    f.Add(cell.initialRoundoffPlusJ);f.Add(cell.initialRoundoffMinusJ);
    for(const auto& component:cell.magneticEnergyJ) {
      f.Add(static_cast<std::uint64_t>(component.size()));for(double value:component)f.Add(value);
    }
    const auto& d=cell.domain;f.Add(static_cast<std::uint64_t>(d.authority));
    f.Add(d.provenance);f.Add(d.plasmaCompositionIdentity);f.Add(d.thermalStateIdentity);
    f.Add(d.maximumDeltaBOverB);f.Add(d.maximumKThermalIonGyroradius);f.Add(d.thermalIonGyroradiusM);
    f.Add(d.maximumAlfvenSpeedOverC);f.Add(d.ionInertialLengthM);f.Add(d.maximumKIonInertialLength);
    f.Add(d.ionCyclotronFrequencyRadPerS);f.Add(d.maximumWaveFrequencyOverIonCyclotron);
  }
  return f.Text();
}
Transport::Status Project(double total,const std::vector<double>& weights,
                           std::vector<double>* bins) {
  bins->resize(weights.size());
  for(std::size_t j=0;j<weights.size();++j) {
    (*bins)[j]=total*weights[j];
    if(!Nonnegative((*bins)[j])||(total>0&&weights[j]>0&&(*bins)[j]==0))
      return Transport::Status::Error(Transport::StatusCode::UnresolvedCoefficient,
          "positive polarized bin energy underflowed or overflowed");
  }
  return Transport::Status::Ok();
}
} // namespace

GridResult ResolveCenteredLogGrid(const CenteredLogGrid& input) {
  GridResult result;
  if(input.centersPerM.size()<2||!Positive(input.logWidth)) {
    result.status=Error("centered log grid needs at least two centers and explicit positive Delta ln k");
    return result;
  }
  const double tolerance=128*std::numeric_limits<double>::epsilon();
  for(std::size_t j=0;j<input.centersPerM.size();++j) {
    if(!Positive(input.centersPerM[j])||(j&&input.centersPerM[j]<=input.centersPerM[j-1])) {
      result.status=Error("wave-number centers must be finite, positive and ordered");return result;
    }
    if(j) {
      const double expected=input.centersPerM.front()*std::exp(j*input.logWidth);
      if(!Positive(expected)||std::fabs(input.centersPerM[j]-expected)>
          tolerance*std::max(input.centersPerM[j],expected)) {
        result.status=Error("wave-number centers disagree with their declared uniform log width");return result;
      }
    }
  }
  Grid candidate;candidate.centersPerM=input.centersPerM;candidate.logWidth=input.logWidth;
  candidate.edgesPerM.resize(input.centersPerM.size()+1);
  const double logarithm=std::log(input.centersPerM.front());
  for(std::size_t j=0;j<candidate.edgesPerM.size();++j) {
    candidate.edgesPerM[j]=std::exp(logarithm+(static_cast<double>(j)-0.5)*input.logWidth);
    if(!Positive(candidate.edgesPerM[j])||(j&&candidate.edgesPerM[j]<=candidate.edgesPerM[j-1])) {
      result.status=Error("centered-grid edge conversion overflowed or lost resolution");return result;
    }
  }
  for(std::size_t j=0;j<input.centersPerM.size();++j)
    if(!(candidate.edgesPerM[j]<input.centersPerM[j]&&
         input.centersPerM[j]<candidate.edgesPerM[j+1])) {
      result.status=Error("center does not lie strictly inside its derived physical bin");return result;
    }
  result.grid=std::move(candidate);result.status=Transport::Status::Ok();return result;
}

GridResult ResolveCenteredLogGrid(std::size_t bins,double first,double last) {
  GridResult result;
  // This is an explicit allocation bound, never a physical spectral cutoff.
  if(bins<2||bins>10000000||!Positive(first)||!Positive(last)||last<=first) {
    result.status=Error("centered spectral grid requires 2..10000000 bins and finite positive ordered centers");return result;
  }
  CenteredLogGrid input;input.logWidth=std::log(last/first)/(bins-1);
  if(!Positive(input.logWidth)) {
    result.status=Error("centered spectral logarithmic width is not representable");return result;
  }
  input.centersPerM.resize(bins);
  for(std::size_t j=0;j<bins;++j)input.centersPerM[j]=first*std::exp(j*input.logWidth);
  return ResolveCenteredLogGrid(input);
}

Transport::Status CopyLegacyBranchSpectrum(const CenteredLogGrid& requested,
    const CenteredLogGrid& stored,const double* source,std::size_t sourceLength,
    std::vector<double>* destination) {
  if(!destination)return Error("legacy spectral destination is absent");
  const auto a=ResolveCenteredLogGrid(requested),b=ResolveCenteredLogGrid(stored);
  if(!a.status.ok())return a.status;
  if(!b.status.ok())return b.status;
  if(a.grid.centersPerM!=b.grid.centersPerM||a.grid.logWidth!=b.grid.logWidth||
      sourceLength!=2*b.grid.centersPerM.size())
    return Transport::Status::Error(Transport::StatusCode::UnsupportedConfiguration,
        "legacy native spectrum count/physical centered grid differs from configuration; no copy or remap performed");
  if(!source)return Error("legacy spectral source is absent");
  std::vector<double> candidate;candidate.reserve(sourceLength);
  for(std::size_t j=0;j<sourceLength;++j) {
    if(!Nonnegative(source[j]))return Error("legacy magnetic spectral energy is nonphysical");
    candidate.push_back(source[j]);
  }
  destination->swap(candidate);return Transport::Status::Ok();
}

BranchDiagnostics DeriveBranches(const Cell& cell) {
  BranchDiagnostics result;
  const auto bins=cell.magneticEnergyJ[0].size();
  if(bins==0){result.status=Error("polarized branch diagnostic has no authoritative bins");return result;}
  for(const auto& component:cell.magneticEnergyJ) {
    if(component.size()!=bins){result.status=Error("polarized diagnostic array lengths differ");return result;}
    for(double value:component)if(!Nonnegative(value)) {
      result.status=Error("polarized diagnostic energy is nonphysical");return result;
    }
  }
  result.ePlusJ=static_cast<double>(Sum(cell.magneticEnergyJ[0])+Sum(cell.magneticEnergyJ[1]));
  result.eMinusJ=static_cast<double>(Sum(cell.magneticEnergyJ[2])+Sum(cell.magneticEnergyJ[3]));
  result.status=Nonnegative(result.ePlusJ)&&Nonnegative(result.eMinusJ)&&Nonnegative(result.ePlusJ+result.eMinusJ)?Transport::Status::Ok():
      Error("polarized branch sum overflowed");
  return result;
}

SnapshotResult PrepareInitialization(const Initialization& input) {
  SnapshotResult result;
  result.status=CheckPolarization(input.plus,"E+");if(!result.status.ok())return result;
  result.status=CheckPolarization(input.minus,"E-");if(!result.status.ok())return result;
  const auto resolved=ResolveCenteredLogGrid(input.grid);
  if(!resolved.status.ok()){result.status=resolved.status;return result;}
  if(!std::isfinite(input.powerLawIndex)||!Positive(input.vacuumPermeabilityHPerM)||
      !std::isfinite(input.epochS)||std::isnan(input.validUntilS)||input.validUntilS<input.epochS||
      input.geometryGeneration==0||input.backgroundIdentity.empty()||
      input.configurationFingerprint.empty()||input.initializerProvenance.empty()||input.cells.empty()) {
    result.status=Error("polarized initialization lacks physical grid, epoch, ownership or initializer provenance");return result;
  }
  const std::size_t bins=input.grid.centersPerM.size();
  std::vector<double> weights(bins),logWeights(bins);
  for(std::size_t j=0;j<bins;++j)logWeights[j]=(1-input.powerLawIndex)*std::log(input.grid.centersPerM[j]);
  const double maximum=*std::max_element(logWeights.begin(),logWeights.end());
  if(!std::isfinite(maximum)){result.status=Error("initial power-law exponent overflowed");return result;}
  for(std::size_t j=0;j<bins;++j) {
    weights[j]=std::exp(logWeights[j]-maximum);
    if(!Positive(weights[j])) {
      result.status=Transport::Status::Error(Transport::StatusCode::UnresolvedCoefficient,
          "positive initial spectral weight is outside representable dynamic range");return result;
    }
  }
  const long double normalization=Sum(weights);
  for(double& value:weights) {
    value=static_cast<double>(value/normalization);
    if(!Positive(value)) {
      result.status=Transport::Status::Error(Transport::StatusCode::UnresolvedCoefficient,
          "positive normalized spectral weight underflowed");return result;
    }
  }
  auto candidate=std::make_shared<Snapshot>();
  candidate->grid=resolved.grid;candidate->epochS=input.epochS;candidate->validUntilS=input.validUntilS;
  candidate->geometryGeneration=input.geometryGeneration;
  candidate->vacuumPermeabilityHPerM=input.vacuumPermeabilityHPerM;
  candidate->backgroundIdentity=input.backgroundIdentity;
  candidate->configurationFingerprint=input.configurationFingerprint;
  candidate->initializerProvenance=input.initializerProvenance;
  candidate->initialPlus=input.plus;candidate->initialMinus=input.minus;
  candidate->initialPowerLawIndex=input.powerLawIndex;
  std::set<std::pair<std::uint64_t,std::uint64_t>> identities;
  for(const auto& initial:input.cells) {
    if(!identities.insert({initial.fieldLine,initial.segment}).second) {
      result.status=Error("duplicate field-line/segment identity in polarized initialization");return result;
    }
    if(!Positive(initial.volumeM3)||!Positive(initial.meanMagneticFieldT)||
        !Positive(initial.massDensityKgM3)||!Nonnegative(initial.ePlusJ)||!Nonnegative(initial.eMinusJ)) {
      result.status=Error("polarized cell requires positive volume/B/rho and nonnegative magnetic branch energy");return result;
    }
    // An unset (zero) domain number is explicit unavailable metadata, not a
    // numerical physical limit. Malformed numbers cannot enter a supposedly
    // structurally valid snapshot even though qualification is deferred.
    result.status=CheckDomainMetadata(initial.domain);if(!result.status.ok())return result;
    Cell cell;cell.fieldLine=initial.fieldLine;cell.segment=initial.segment;
    cell.volumeM3=initial.volumeM3;cell.meanMagneticFieldT=initial.meanMagneticFieldT;
    cell.massDensityKgM3=initial.massDensityKgM3;
    cell.alfvenSpeedMPerS=cell.meanMagneticFieldT/
        std::sqrt(input.vacuumPermeabilityHPerM*cell.massDensityKgM3);
    cell.domain=initial.domain;
    if(!Positive(cell.alfvenSpeedMPerS)||cell.alfvenSpeedMPerS>=299792458.0) {
      result.status=Error("polarized cell Alfven speed is outside the declared MHD/plasma frame");return result;
    }
    const double totals[4]={initial.ePlusJ*input.plus.leftFraction,
        initial.ePlusJ*(1-input.plus.leftFraction),
        initial.eMinusJ*input.minus.leftFraction,initial.eMinusJ*(1-input.minus.leftFraction)};
    for(int component=0;component<4;++component) {
      const double branch=component<2?initial.ePlusJ:initial.eMinusJ;
      const double fraction=component==0?input.plus.leftFraction:component==1?
          1-input.plus.leftFraction:component==2?input.minus.leftFraction:1-input.minus.leftFraction;
      if(branch>0&&fraction>0&&totals[component]==0) {
        result.status=Transport::Status::Error(Transport::StatusCode::UnresolvedCoefficient,
            "nonzero requested helicity underflowed during initial partition");return result;
      }
      result.status=Project(totals[component],weights,&cell.magneticEnergyJ[component]);
      if(!result.status.ok())return result;
    }
    const auto diagnostic=DeriveBranches(cell);
    if(!diagnostic.status.ok()){result.status=diagnostic.status;return result;}
    cell.initialRoundoffPlusJ=diagnostic.ePlusJ-initial.ePlusJ;
    cell.initialRoundoffMinusJ=diagnostic.eMinusJ-initial.eMinusJ;
    const double tolerance=128*std::numeric_limits<double>::epsilon();
    if(std::fabs(cell.initialRoundoffPlusJ)>tolerance*initial.ePlusJ||
        std::fabs(cell.initialRoundoffMinusJ)>tolerance*initial.eMinusJ) {
      result.status=Transport::Status::Error(Transport::StatusCode::UnresolvedCoefficient,
          "polarized spectral projection exceeds numerical roundoff bounds");return result;
    }
    candidate->cells.push_back(std::move(cell));
  }
  std::sort(candidate->cells.begin(),candidate->cells.end(),[](const Cell& a,const Cell& b) {
    return std::make_pair(a.fieldLine,a.segment)<std::make_pair(b.fieldLine,b.segment);
  });
  candidate->identity=SnapshotIdentity(*candidate);
  candidate->structurallyValid=true;candidate->qualifiedForProduction=false;
  result.snapshot=std::move(candidate);result.status=Transport::Status::Ok();return result;
}

Transport::Status ValidateSnapshot(const Snapshot& state) {
  if(state.schema!=StateSchema||state.revision==0||!state.structurallyValid||state.qualifiedForProduction||
      state.geometryGeneration==0||!std::isfinite(state.epochS)||std::isnan(state.validUntilS)||
      state.validUntilS<state.epochS||!Positive(state.vacuumPermeabilityHPerM)||
      !std::isfinite(state.initialPowerLawIndex)||state.backgroundIdentity.empty()||
      state.configurationFingerprint.empty()||state.initializerProvenance.empty()||state.cells.empty())
    return Error("polarized snapshot schema, metadata or implementation capability is invalid");
  auto status=CheckPolarization(state.initialPlus,"E+");if(!status.ok())return status;
  status=CheckPolarization(state.initialMinus,"E-");if(!status.ok())return status;
  CenteredLogGrid centered;centered.centersPerM=state.grid.centersPerM;centered.logWidth=state.grid.logWidth;
  const auto grid=ResolveCenteredLogGrid(centered);
  if(!grid.status.ok())return grid.status;
  if(grid.grid.edgesPerM!=state.grid.edgesPerM)return Error("polarized snapshot physical edges disagree with declared centers/width");
  std::pair<std::uint64_t,std::uint64_t> previous;bool first=true;
  for(const auto& cell:state.cells) {
    const auto id=std::make_pair(cell.fieldLine,cell.segment);
    if((!first&&id<=previous)||!Positive(cell.volumeM3)||!Positive(cell.meanMagneticFieldT)||
        !Positive(cell.massDensityKgM3)||!Positive(cell.alfvenSpeedMPerS)||cell.alfvenSpeedMPerS>=299792458.0||
        !std::isfinite(cell.initialRoundoffPlusJ)||!std::isfinite(cell.initialRoundoffMinusJ))
      return Error("polarized snapshot physical cell or canonical identity ordering is invalid");
    previous=id;first=false;
    const double va=cell.meanMagneticFieldT/std::sqrt(state.vacuumPermeabilityHPerM*cell.massDensityKgM3);
    if(!Positive(va)||std::fabs(va-cell.alfvenSpeedMPerS)>
        128*std::numeric_limits<double>::epsilon()*std::max(va,cell.alfvenSpeedMPerS))
      return Error("polarized snapshot Alfven speed disagrees with its B/rho/permeability");
    status=CheckDomainMetadata(cell.domain);if(!status.ok())return status;
    for(const auto& component:cell.magneticEnergyJ)if(component.size()!=state.grid.centersPerM.size())
      return Error("polarized snapshot four-spectrum length differs from its physical grid");
    const auto branches=DeriveBranches(cell);if(!branches.status.ok())return branches.status;
  }
  return state.identity==SnapshotIdentity(state)?Transport::Status::Ok():
      Error("polarized snapshot identity does not bind its four arrays, physical metadata and signed ledgers");
}

SnapshotResult PrepareSnapshotForRestore(const Snapshot& decoded) {
  SnapshotResult result;result.status=ValidateSnapshot(decoded);
  if(result.status.ok())result.snapshot=std::make_shared<const Snapshot>(decoded);
  return result;
}

EvolutionCandidateResult PrepareEvolutionCandidate(const Snapshot& before,
    const std::vector<Cell>& cells,const EvolutionMetadata& metadata) {
  EvolutionCandidateResult result;
  // Overflow is checked before allocation or increment. Even a malformed
  // decoded header cannot wrap revision back to a seemingly fresh state.
  if(before.revision==std::numeric_limits<std::uint64_t>::max()) {
    result.status=Error("polarized evolution revision is exhausted");return result;
  }
  result.status=ValidateSnapshot(before);if(!result.status.ok())return result;
  if(!std::isfinite(metadata.epochS)||std::isnan(metadata.validUntilS)||
      metadata.epochS<before.epochS||metadata.validUntilS<metadata.epochS||
      metadata.geometryGeneration==0||metadata.backgroundIdentity.empty()||metadata.operationProvenance.empty()||
      cells.size()!=before.cells.size()) {
    result.status=Error("polarized evolution needs explicit forward epoch/background, operation provenance and fixed cell count");return result;
  }
  for(std::size_t i=0;i<cells.size();++i) {
    const auto& old=before.cells[i];const auto& next=cells[i];
    if(next.fieldLine!=old.fieldLine||next.segment!=old.segment||
        std::memcmp(&next.initialRoundoffPlusJ,&old.initialRoundoffPlusJ,sizeof(double))||
        std::memcmp(&next.initialRoundoffMinusJ,&old.initialRoundoffMinusJ,sizeof(double))) {
      result.status=Error("polarized evolution cannot reorder cells or alter the signed initial projection ledger");return result;
    }
    for(const auto& component:next.magneticEnergyJ)if(component.size()!=before.grid.centersPerM.size()) {
      result.status=Error("polarized evolution cannot resize/remap the exact spectral grid");return result;
    }
  }
  try {
    auto candidate=std::make_shared<Snapshot>(before);
    candidate->cells=cells;candidate->revision=before.revision+1;
    candidate->epochS=metadata.epochS;candidate->validUntilS=metadata.validUntilS;
    candidate->geometryGeneration=metadata.geometryGeneration;candidate->backgroundIdentity=metadata.backgroundIdentity;
    candidate->qualifiedForProduction=false;
    // This is the single canonical builder used by initialization and restore
    // validation; callers never duplicate its bit-level physical identity.
    candidate->identity=SnapshotIdentity(*candidate);
    result.status=ValidateSnapshot(*candidate);if(!result.status.ok())return result;
    result.operation.beforeIdentity=before.identity;result.operation.afterIdentity=candidate->identity;
    result.operation.operationProvenance=metadata.operationProvenance;
    result.snapshot=std::move(candidate);return result;
  } catch(const std::exception& error) {
    result.status=Error(std::string("polarized evolution staging failed before publication: ")+error.what());return result;
  }
}

CollisionViewResult CollisionView(const Snapshot& state,std::size_t cellIndex) {
  CollisionViewResult result;
  if(!state.structurallyValid||state.schema!=StateSchema||cellIndex>=state.cells.size()) {
    result.status=Error("collision view requires a valid polarized snapshot and cell index");return result;
  }
  const auto& cell=state.cells[cellIndex];auto& spectrum=result.spectrum;
  spectrum.kEdgesPerM=state.grid.edgesPerM;spectrum.magneticEnergyJ=cell.magneticEnergyJ;
  spectrum.polarizationResolved=true;spectrum.volumeM3=cell.volumeM3;
  spectrum.meanMagneticFieldT=cell.meanMagneticFieldT;spectrum.alfvenSpeedMPerS=cell.alfvenSpeedMPerS;
  spectrum.vacuumPermeabilityHPerM=state.vacuumPermeabilityHPerM;
  spectrum.epochS=state.epochS;spectrum.generation=state.geometryGeneration;
  std::ostringstream identity;identity<<state.identity<<":r"<<state.revision
      <<":line"<<cell.fieldLine<<":segment"<<cell.segment;
  spectrum.identity=identity.str();spectrum.domain=cell.domain;
  result.status=Transport::ResonantCollision::ValidateSpectrum(spectrum);return result;
}

Transport::Status StateOwner::PublishInitialization(const Initialization& initialization) {
  const auto candidate=PrepareInitialization(initialization);
  if(!candidate.status.ok())return candidate.status;
  std::lock_guard<std::mutex> lock(mutex_);
  if(current_)return Transport::Status::Error(Transport::StatusCode::UnsupportedConfiguration,
      "polarized state already initialized; replacement requires a physical update/restart transaction");
  current_=candidate.snapshot;return Transport::Status::Ok();
}
std::shared_ptr<const Snapshot> StateOwner::Acquire()const {
  std::lock_guard<std::mutex> lock(mutex_);return current_;
}

} } }
