#include "sep_polarized_wave_state_validation.h"
#include "sep_polarized_wave_state.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>

namespace SEP { namespace Testing { namespace {
namespace PW=Turbulence::Polarized;
namespace RC=Transport::ResonantCollision;
bool Near(double a,double b,double relative=1e-12) {
  return std::isfinite(a)&&std::isfinite(b)&&
      std::fabs(a-b)<=relative*std::max(std::fabs(b),1e-300);
}
Result Finish(bool pass,const char* message) {
  Result result;result.status=pass?Status::Pass:Status::Fail;result.message=message;
  result.configuration.push_back("capability=non-PIC-polarized-state-owner;production_qualification=false");
  result.configuration.push_back("initialization=existing-Eplus-Eminus-profile-plus-explicit-polarization");
  return result;
}
PW::Initialization Initialization() {
  PW::Initialization input;
  input.grid.logWidth=std::log(1e-2/1e-8)/127;
  for(int j=0;j<128;++j)input.grid.centersPerM.push_back(1e-8*std::exp(j*input.grid.logWidth));
  input.plus.known=true;input.plus.leftFraction=0.2;input.plus.provenance="explicit manufactured + polarization";
  input.minus.known=true;input.minus.leftFraction=0.7;input.minus.provenance="explicit manufactured - polarization";
  input.vacuumPermeabilityHPerM=1.256637062e-6;input.epochS=0;input.validUntilS=1;
  input.geometryGeneration=1;input.backgroundIdentity="controlled-analytic:g1:t0";
  input.configurationFingerprint="controlled-initializer-config";
  input.initializerProvenance="InitializeWaveEnergyFromPhysicalParameters: manufactured integrated outputs";
  PW::InitialCell cell;cell.fieldLine=1;cell.segment=3;cell.volumeM3=1e6;
  cell.meanMagneticFieldT=5e-9;cell.massDensityKgM3=5e6*1.67262192369e-27;
  cell.ePlusJ=4e-7;cell.eMinusJ=2e-7;
  cell.domain.authority=RC::DomainAuthority::ValidationStudy;
  cell.domain.provenance="initialization mathematical verification; not MHD qualification";
  input.cells.push_back(cell);cell.segment=2;cell.ePlusJ*=1.5;cell.eMinusJ*=1.5;
  input.cells.push_back(cell);return input;
}
Result UnknownPolarization() {
  PW::StateOwner owner;auto input=Initialization();input.plus.known=false;
  const auto missing=owner.PublishInitialization(input);
  bool pass=!missing.ok()&&missing.code==Transport::StatusCode::UnsupportedConfiguration&&!owner.Acquire();
  input.plus.known=true;input.minus.known=false;
  pass=pass&&!owner.PublishInitialization(input).ok()&&!owner.Acquire();
  input.minus.known=true;
  const auto valid=owner.PublishInitialization(input);const auto state=owner.Acquire();
  pass=pass&&valid.ok()&&state&&state->structurallyValid&&!state->qualifiedForProduction;
  if(state) {
    const auto& cell=state->cells[0];const auto total=PW::DeriveBranches(cell);
    const double leftPlus=std::accumulate(cell.magneticEnergyJ[0].begin(),cell.magneticEnergyJ[0].end(),0.0);
    const double leftMinus=std::accumulate(cell.magneticEnergyJ[2].begin(),cell.magneticEnergyJ[2].end(),0.0);
    pass=pass&&total.status.ok()&&Near(leftPlus/total.ePlusJ,0.2)&&Near(leftMinus/total.eMinusJ,0.7);
  }
  return Finish(pass,"unknown L/R fractions publish no state; explicit unequal fractions are preserved independently in both propagation branches");
}
Result CenteredGridContract() {
  const auto input=Initialization();const auto grid=PW::ResolveCenteredLogGrid(input.grid);
  if(!grid.status.ok())return Finish(false,"declared centered log grid failed to resolve");
  bool pass=grid.grid.edgesPerM.size()==129;
  const auto& edges=grid.grid.edgesPerM;
  pass=pass&&Near(edges.front(),1e-8*std::exp(-input.grid.logWidth/2))&&
      Near(edges.back(),1e-2*std::exp(input.grid.logWidth/2));
  double maximumError=0;
  for(std::size_t j=1;j<128;++j) {
    const double expected=std::sqrt(input.grid.centersPerM[j-1]*input.grid.centersPerM[j]);
    maximumError=std::max(maximumError,std::fabs(edges[j]/expected-1));
  }
  pass=pass&&maximumError<1e-13;
  auto malformed=input.grid;malformed.logWidth=0;
  pass=pass&&!PW::ResolveCenteredLogGrid(malformed).status.ok();
  malformed=input.grid;malformed.centersPerM[17]*=1.01;
  pass=pass&&!PW::ResolveCenteredLogGrid(malformed).status.ok();
  auto result=Finish(pass,"legacy logarithmic centers and explicit width define geometric interior edges and outer half bins without truncation");
  result.metrics.push_back({"maximum_edge_relative_error",maximumError,1e-13,"<","dimensionless"});
  return result;
}
Result ExactBinIntegral() {
  const auto input=Initialization();const auto prepared=PW::PrepareInitialization(input);
  if(!prepared.status.ok())return Finish(false,"four-spectrum initialization failed");
  const auto& snapshot=*prepared.snapshot;const auto& cell=snapshot.cells[0];
  const auto& edges=snapshot.grid.edgesPerM;const double q=input.powerLawIndex;
  const double entire=(std::pow(edges.front(),1-q)-std::pow(edges.back(),1-q))/(q-1);
  const auto original=std::find_if(input.cells.begin(),input.cells.end(),[&](const PW::InitialCell& c){return c.segment==cell.segment;});
  const double totals[4]={original->ePlusJ*input.plus.leftFraction,original->ePlusJ*(1-input.plus.leftFraction),
      original->eMinusJ*input.minus.leftFraction,original->eMinusJ*(1-input.minus.leftFraction)};
  double maximumError=0;
  for(int component=0;component<4;++component)for(std::size_t j=0;j+1<edges.size();++j) {
    const double integral=(std::pow(edges[j],1-q)-std::pow(edges[j+1],1-q))/(q-1);
    const double expected=totals[component]*integral/entire;
    maximumError=std::max(maximumError,std::fabs(cell.magneticEnergyJ[component][j]/expected-1));
  }
  const auto derived=PW::DeriveBranches(cell);
  const bool pass=maximumError<1e-12&&derived.status.ok()&&Near(derived.ePlusJ,original->ePlusJ)&&
      Near(derived.eMinusJ,original->eMinusJ)&&
      Near(derived.ePlusJ-original->ePlusJ,cell.initialRoundoffPlusJ)&&
      Near(derived.eMinusJ-original->eMinusJ,cell.initialRoundoffMinusJ);
  auto result=Finish(pass,"existing Kolmogorov center weights equal exact finite log-bin integrals and retain integrated branch energy with explicit rounding residuals");
  result.metrics.push_back({"maximum_bin_integral_relative_error",maximumError,1e-12,"<","dimensionless"});
  return result;
}
Result ResonantAuthority() {
  const auto input=Initialization();const auto prepared=PW::PrepareInitialization(input);
  if(!prepared.status.ok())return Finish(false,"owner initialization failed before collision-view test");
  const auto view=PW::CollisionView(*prepared.snapshot,0);
  if(!view.status.ok())return Finish(false,"authoritative collision view was unavailable");
  RC::Particle particle;particle.massKg=1.67262192369e-27;particle.chargeC=1.602176634e-19;
  const double speed=2e7,c=299792458.0;
  const double momentum=particle.massKg*speed/std::sqrt((1-speed/c)*(1+speed/c));
  particle.parallelMomentumKgMPerS=0.4*momentum;
  particle.perpendicularMomentumKgMPerS=std::sqrt(0.84)*momentum;
  const auto orbit=RC::PrepareOrbit(view.spectrum,RC::Branch::AlongB,particle);
  if(!orbit.status.ok())return Finish(false,"validation-study orbit failed");
  const auto baseline=RC::EvaluateResonance(view.spectrum,orbit.orbit,0.4);
  if(!baseline.status.ok()||!baseline.active)return Finish(false,"manufactured resonant bin was inactive");
  auto selected=view.spectrum;selected.identity="controlled-resonant-bin-perturbation";
  selected.magneticEnergyJ[1][baseline.bin]*=4;
  const auto changed=RC::EvaluateResonance(selected,orbit.orbit,0.4);
  auto unselected=view.spectrum;unselected.identity="controlled-unselected-helicity-perturbation";
  unselected.magneticEnergyJ[0][baseline.bin]*=4;
  const auto unchanged=RC::EvaluateResonance(unselected,orbit.orbit,0.4);
  const bool pass=changed.status.ok()&&unchanged.status.ok()&&
      Near(changed.frequencyPerS/baseline.frequencyPerS,4)&&Near(unchanged.frequencyPerS,baseline.frequencyPerS)&&
      baseline.helicity==RC::Helicity::Right&&!prepared.snapshot->qualifiedForProduction;
  return Finish(pass,"collision kernels consume actual authoritative polarized bins, responding only to the selected helicity/bin rather than an aggregate assumed spectrum");
}
Result PublicationTransaction() {
  auto input=Initialization();PW::StateOwner owner;
  auto invalid=input;invalid.cells.back().volumeM3=0;
  bool pass=!owner.PublishInitialization(invalid).ok()&&!owner.Acquire();
  invalid=input;invalid.cells.back().domain.maximumAlfvenSpeedOverC=std::numeric_limits<double>::quiet_NaN();
  pass=pass&&!owner.PublishInitialization(invalid).ok()&&!owner.Acquire();
  invalid=input;invalid.cells.push_back(invalid.cells.front());
  pass=pass&&!owner.PublishInitialization(invalid).ok()&&!owner.Acquire();
  pass=pass&&owner.PublishInitialization(input).ok();const auto published=owner.Acquire();
  if(!published)return Finish(false,"valid candidate failed to publish");
  const auto identity=published->identity;
  input.cells.front().ePlusJ*=9;
  pass=pass&&!owner.PublishInitialization(input).ok()&&owner.Acquire()==published&&
      published->identity==identity;
  auto reordered=Initialization();std::reverse(reordered.cells.begin(),reordered.cells.end());
  const auto canonical=PW::PrepareInitialization(reordered);
  pass=pass&&canonical.status.ok()&&canonical.snapshot->identity==identity;
  auto unknownDomain=Initialization();for(auto& cell:unknownDomain.cells)cell.domain={};
  const auto structural=PW::PrepareInitialization(unknownDomain);
  pass=pass&&structural.status.ok()&&structural.snapshot->structurallyValid&&!structural.snapshot->qualifiedForProduction;
  if(structural.status.ok()) {
    const auto collision=PW::CollisionView(*structural.snapshot,0);
    pass=pass&&collision.status.ok()&&!RC::ValidateDomain(collision.spectrum).ok();
  }
  return Finish(pass,"failed or repeated initializations never replace immutable authority; canonical identity is input-order independent and structural validity does not imply physical qualification");
}
Result NativeLayoutGuard() {
  const auto baseline=Initialization().grid;
  std::vector<double> source(256,2),destination={3,5,7};const auto saved=destination;
  bool pass=PW::CopyLegacyBranchSpectrum(baseline,baseline,source.data(),source.size(),&destination).ok()&&
      destination==source;
  destination=saved;
  for(std::size_t count:{127u,129u,257u}) {
    const auto changed=PW::ResolveCenteredLogGrid(count,1e-8,1e-2);
    PW::CenteredLogGrid requested;requested.centersPerM=changed.grid.centersPerM;
    requested.logWidth=changed.grid.logWidth;
    // nullptr proves rejection precedes even the first native buffer read.
    pass=pass&&changed.status.ok()&&!PW::CopyLegacyBranchSpectrum(requested,baseline,
        nullptr,source.size(),&destination).ok()&&destination==saved;
  }
  auto wrongK=baseline;for(double& value:wrongK.centersPerM)value*=0.01;
  pass=pass&&!PW::CopyLegacyBranchSpectrum(wrongK,baseline,nullptr,source.size(),&destination).ok()&&destination==saved&&
      !PW::CopyLegacyBranchSpectrum(baseline,baseline,nullptr,258,&destination).ok()&&destination==saved&&
      !PW::ResolveCenteredLogGrid(1,1e-8,1e-2).status.ok();
  source.back()=-1;
  pass=pass&&!PW::CopyLegacyBranchSpectrum(baseline,baseline,source.data(),source.size(),&destination).ok()&&destination==saved;
  return Finish(pass,"legacy native count/physical-grid mismatches and malformed energies reject before source reads or destination mutation, including N>NK");
}
Result RestoreValueContract() {
  const auto initialized=PW::PrepareInitialization(Initialization());
  if(!initialized.status.ok())return Finish(false,"initialization prerequisite for immutable restore failed");
  auto decoded=*initialized.snapshot;
  const auto restored=PW::PrepareSnapshotForRestore(decoded);
  bool pass=restored.status.ok()&&restored.snapshot.get()!=initialized.snapshot.get()&&
      restored.snapshot->identity==initialized.snapshot->identity&&
      restored.snapshot->cells[0].magneticEnergyJ==initialized.snapshot->cells[0].magneticEnergyJ;
  // Decoder storage remains mutable; a successful stage must own a deep copy.
  decoded.cells[0].magneticEnergyJ[0][0]*=2;
  pass=pass&&restored.snapshot->cells[0].magneticEnergyJ[0][0]==initialized.snapshot->cells[0].magneticEnergyJ[0][0]&&
      !PW::PrepareSnapshotForRestore(decoded).status.ok();
  decoded=*initialized.snapshot;decoded.cells[0].initialRoundoffPlusJ=1;
  pass=pass&&!PW::PrepareSnapshotForRestore(decoded).status.ok();
  decoded=*initialized.snapshot;decoded.grid.edgesPerM[0]*=0.99;
  pass=pass&&!PW::PrepareSnapshotForRestore(decoded).status.ok();
  decoded=*initialized.snapshot;decoded.revision++;
  pass=pass&&!PW::PrepareSnapshotForRestore(decoded).status.ok();
  return Finish(pass,"immutable restore stages the exact four arrays and signed ledgers; changed bins, edges, revision or ledger reject rather than reconstructing a power law or equal helicities");
}
Result EvolutionCandidateContract() {
  const auto initialized=PW::PrepareInitialization(Initialization());
  if(!initialized.status.ok())return Finish(false,"initialization prerequisite for canonical evolution failed");
  const auto before=initialized.snapshot;auto cells=before->cells;
  PW::EvolutionMetadata target;target.epochS=0.5;target.validUntilS=1;
  target.geometryGeneration=before->geometryGeneration;target.backgroundIdentity="controlled-analytic:g1:t0.5";
  target.operationProvenance="manufactured single-helicity edit; bookkeeping verification, not production qualification";
  const double old=cells[0].magneticEnergyJ[1][7];cells[0].magneticEnergyJ[1][7]*=1.25;
  const auto evolved=PW::PrepareEvolutionCandidate(*before,cells,target);
  bool pass=evolved.status.ok()&&evolved.snapshot&&evolved.snapshot->revision==before->revision+1&&
      evolved.snapshot->identity!=before->identity&&PW::ValidateSnapshot(*evolved.snapshot).ok()&&
      evolved.snapshot->grid.edgesPerM==before->grid.edgesPerM&&
      evolved.snapshot->initializerProvenance==before->initializerProvenance&&
      evolved.snapshot->configurationFingerprint==before->configurationFingerprint&&!evolved.snapshot->qualifiedForProduction&&
      evolved.operation.beforeIdentity==before->identity&&evolved.operation.afterIdentity==evolved.snapshot->identity&&
      evolved.operation.operationProvenance==target.operationProvenance;
  // A staged update owns its own arrays and retains unequal polarizations;
  // changing caller scratch storage or another component cannot alter it.
  if(evolved.snapshot) {
    cells[0].magneticEnergyJ[1][7]*=2;
    pass=pass&&Near(evolved.snapshot->cells[0].magneticEnergyJ[1][7],old*1.25)&&
        evolved.snapshot->cells[0].magneticEnergyJ[0]==before->cells[0].magneticEnergyJ[0]&&
        before->cells[0].magneticEnergyJ[1][7]==old;
  }
  cells=before->cells;std::swap(cells[0],cells[1]);
  pass=pass&&!PW::PrepareEvolutionCandidate(*before,cells,target).status.ok();
  cells=before->cells;cells[0].initialRoundoffPlusJ+=1e-20;
  pass=pass&&!PW::PrepareEvolutionCandidate(*before,cells,target).status.ok();
  cells=before->cells;cells[0].magneticEnergyJ[3].pop_back();
  pass=pass&&!PW::PrepareEvolutionCandidate(*before,cells,target).status.ok();
  cells=before->cells;cells[0].magneticEnergyJ[0][0]=-1;
  pass=pass&&!PW::PrepareEvolutionCandidate(*before,cells,target).status.ok();
  cells=before->cells;auto invalid=target;invalid.epochS=-1;
  pass=pass&&!PW::PrepareEvolutionCandidate(*before,cells,invalid).status.ok();
  invalid=target;invalid.operationProvenance.clear();
  pass=pass&&!PW::PrepareEvolutionCandidate(*before,cells,invalid).status.ok();
  auto exhausted=*before;exhausted.revision=std::numeric_limits<std::uint64_t>::max();
  const auto overflow=PW::PrepareEvolutionCandidate(exhausted,cells,target);
  pass=pass&&!overflow.status.ok()&&overflow.status.message.find("exhausted")!=std::string::npos;
  return Finish(pass,"canonical staged evolution binds all four arrays with one identity, preserves initialization/configuration, rejects invalid edits and revision wrap, and publishes no owner state");
}
Descriptor Make(const char* id,const char* name,TestCallback callback) {
  Descriptor descriptor;descriptor.id=id;descriptor.name=name;descriptor.group="polarized-wave-state";
  descriptor.description=name;descriptor.initialization=InitializationLevel::None;
  descriptor.supportedBuildModes="native registry and identical dependency-light sanitizer callback";
  descriptor.runtime=RuntimeClass::Routine;descriptor.seedCapability=SeedCapability::Deterministic;
  descriptor.repeatSupported=true;descriptor.seedPolicy="deterministic manufactured initial conditions";
  descriptor.stateIsolation="private immutable four-spectrum snapshots and separate owner per callback";
  descriptor.callback=callback;return descriptor;
}
} // namespace
std::vector<Descriptor> PolarizedWaveStateDescriptors() {
  return {Make("PWS01","Explicit initial wave polarization capability",UnknownPolarization),
          Make("PWS02","Declared logarithmic center-to-edge contract",CenteredGridContract),
          Make("PWS03","Exact initial polarized bin-energy integrals",ExactBinIntegral),
          Make("PWS04","Authoritative resonant polarized-bin consumption",ResonantAuthority),
          Make("PWS05","Transactional immutable wave initialization",PublicationTransaction),
          Make("PWS06","Legacy native spectral layout transaction",NativeLayoutGuard),
          Make("PWS07","Exact immutable four-spectrum restore value",RestoreValueContract),
          Make("PWS08","Canonical staged four-spectrum evolution candidate",EvolutionCandidateContract)};
}
} }
