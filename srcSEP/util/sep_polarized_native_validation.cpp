#include "sep_polarized_native_validation.h"
#include "../polarized_wave_native_adapter.h"
#include "../sep.h"
#include "sep_background_runtime.h"
#include "sep_common_header_path.h"
#include SRCSEP_SEP_COMMON_HEADER(sep_coefficient_registry.h)

#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>
#include <new>
#include <sstream>
#include <type_traits>

namespace SEP { namespace Testing { namespace {
namespace PW=Turbulence::Polarized;
namespace A=PW::PICAdapter;
namespace W=AlfvenTurbulence_Kolmogorov;
namespace FL=PIC::FieldLine;
struct Checks {
  bool pass=true;std::string message;
  void Require(bool ok,const std::string& text) {
    pass=pass&&ok;if(!ok){if(!message.empty())message+="; ";message+=text;}
  }
  bool Status(const Transport::Status& status){Require(status.ok(),status.message);return status.ok();}
};
struct Row {double* data=nullptr;std::vector<double> original;};
class Fixture {
 public:
  Checks& check;std::vector<Row> rows;
  Turbulence::Configuration previous=Turbulence::ActiveConfiguration(),configuration=previous;
  bool ready=false;
  explicit Fixture(Checks& checks):check(checks) {
    configuration.source=Turbulence::Source::SelfConsistentSpectral;
    configuration.representation=Turbulence::Representation::Spectral;
    if(!check.Status(Turbulence::SetActiveConfiguration(configuration)))return;
    for(int line=0;line<FL::nFieldLine;++line)for(int i=0;i<FL::FieldLinesAll[line].GetTotalSegmentNumber();++i) {
      auto* segment=FL::FieldLinesAll[line].GetSegment(i);
      for(auto* datum:{&W::CellIntegratedWaveEnergy,&W::WaveEnergyDensity,&W::WaveNumberResolved::SpectralWaveEnergy}) {
        Row row;row.data=segment->GetDatum_ptr(*datum);
        if(!row.data){check.Require(false,"native diagnostic buffer absent");return;}
        row.original.assign(row.data,row.data+datum->length);rows.push_back(std::move(row));
      }
    }
    ready=true;
  }
  ~Fixture() {
    for(const auto& row:rows)std::copy(row.original.begin(),row.original.end(),row.data);
    check.Require(Turbulence::SetActiveConfiguration(previous).ok(),"cannot restore native fixture turbulence configuration");
  }
  void InitializePhysicalProfile() {
    const double amplitude=Transport::Coefficient::ActiveConfiguration().prescribedDeltaBOverB;
    W::InitializeWaveEnergyFromPhysicalParameters(W::CellIntegratedWaveEnergy,
        W::WaveEnergyConstants::TYPICAL_B0_1AU,amplitude,amplitude,-2.0,false,true);
    FL::Parallel::MPIGatherDatumStoredAtEdge(W::CellIntegratedWaveEnergy,0);
    FL::Parallel::MPIBcastDatumStoredAtEdge(W::CellIntegratedWaveEnergy,0);
  }
  std::vector<double> Bytes()const {
    std::vector<double> result;for(const auto& row:rows)
      result.insert(result.end(),row.data,row.data+row.original.size());return result;
  }
};
bool SameBytes(const std::vector<double>& a,const std::vector<double>& b) {
  return a.size()==b.size()&&(a.empty()||std::memcmp(a.data(),b.data(),a.size()*sizeof(double))==0);
}
Result Finish(const Checks& check,const std::string& message) {
  Result result;result.status=check.pass?Status::Pass:Status::Fail;
  result.message=check.pass?message:check.message;
  result.configuration.push_back("authority=four-independent-magnetic-J-spectra;initial-LR=explicit-0.5-per-branch");
  result.configuration.push_back("grid=resolved-log-centers-plus-outer-half-width;production_qualification=false");
  result.configuration.push_back("physical_evolution=not-enabled;test=collective-native-initialization-or-restore-value-contract");
  return result;
}
bool Prerequisite(Result* result) {
  if(W::WaveNumberResolved::IsActive()&&W::WaveNumberResolved::SpectralWaveEnergy.is_active())return true;
  result->status=Status::Skip;
  result->message="requires native --turbulence-model wave-number-resolved so resolved 2*N diagnostic storage is registered before FieldLine::Init";
  return false;
}
Result NativeMoments() {
  Result missing;if(!Prerequisite(&missing))return missing;
  Checks check;double maximumBinError=0,maximumMomentError=0;std::uint64_t cells=0,bins=0;
  double maximumInitialMomentError=0;
  {
    Fixture fixture(check);
    if(fixture.ready) {
      fixture.InitializePhysicalProfile();A::NativeStateOwner owner;
      // Capture the physical initializer before the owner derives/replaces
      // native diagnostics. A common rescaling of all four arrays must fail.
      std::vector<std::array<double,2>> initialEnergy;
      for(int line=0;line<FL::nFieldLine;++line)for(int i=0;i<FL::FieldLinesAll[line].GetTotalSegmentNumber();++i) {
        const auto* energy=FL::FieldLinesAll[line].GetSegment(i)->GetDatum_ptr(W::CellIntegratedWaveEnergy);
        initialEnergy.push_back({{energy[0],energy[1]}});
      }
      if(check.Status(owner.InitializeFromIntegrated(fixture.configuration))) {
        const auto state=owner.Acquire();check.Require(state&&!state->qualifiedForProduction,"native initialization falsely claims physical qualification");
        if(state) {
          check.Status(PW::ValidateSnapshot(*state));cells=state->cells.size();bins=state->grid.centersPerM.size();
          const double q=5.0/3.0;
          const double normalization=(std::pow(state->grid.edgesPerM.front(),1-q)-
              std::pow(state->grid.edgesPerM.back(),1-q))/(q-1);
          std::size_t cellIndex=0;
          for(const auto& cell:state->cells) {
            const auto branch=PW::DeriveBranches(cell);check.Status(branch.status);
            const double* native=FL::FieldLinesAll[cell.fieldLine].GetSegment(static_cast<int>(cell.segment))->GetDatum_ptr(W::CellIntegratedWaveEnergy);
            for(int sigma=0;sigma<2;++sigma) {
              const double total=sigma?branch.eMinusJ:branch.ePlusJ;
              const double original=initialEnergy[cellIndex][sigma];
              const double rounding=sigma?cell.initialRoundoffMinusJ:cell.initialRoundoffPlusJ;
              if(original>0)maximumInitialMomentError=std::max(maximumInitialMomentError,
                  std::fabs(static_cast<double>((static_cast<long double>(total)-rounding-original)/original)));
              else check.Require(total==0&&rounding==0,"zero physical amplitude acquired four-state energy");
              if(total>0)maximumMomentError=std::max(maximumMomentError,std::fabs(native[sigma]/total-1));
              for(std::size_t j=0;j<bins;++j) {
                const double integral=(std::pow(state->grid.edgesPerM[j],1-q)-
                    std::pow(state->grid.edgesPerM[j+1],1-q))/(q-1);
                const double expected=0.5*original*integral/normalization;
                for(int helicity=0;helicity<2;++helicity)if(expected>0)
                  maximumBinError=std::max(maximumBinError,std::fabs(cell.magneticEnergyJ[2*sigma+helicity][j]/expected-1));
              }
            }
            for(int component=1;component<4;++component)
              check.Require(cell.magneticEnergyJ[component].data()!=cell.magneticEnergyJ[0].data(),"helicity arrays share mutable storage");
            ++cellIndex;
          }
          const auto diagnostics=owner.Diagnostics();
          check.Require(diagnostics.mpiIdentityAgrees&&diagnostics.readbackAgrees,"MPI physical/bin identity or native readback failed");
          check.Require(maximumBinError<1e-12&&maximumMomentError==0,"native bin-integral or branch-moment readback differs from independent reference");
          check.Require(maximumInitialMomentError<=2*std::numeric_limits<double>::epsilon(),
              "four-state branch energy minus signed rounding ledger differs from pre-initialization physical amplitude");
        }
      }
    }
  }
  auto result=Finish(check,"existing physical-profile initializer feeds four independent initial spectra; exact bin integrals and native branch moments agree on every MPI rank");
  result.metrics.push_back({"physical_cells",static_cast<double>(cells),1,">=","count"});
  result.metrics.push_back({"resolved_bins",static_cast<double>(bins),2,">=","count"});
  result.metrics.push_back({"maximum_independent_bin_integral_error",maximumBinError,1e-12,"<","dimensionless"});
  result.metrics.push_back({"maximum_native_moment_readback_error",maximumMomentError,0,"==","dimensionless"});
  result.metrics.push_back({"maximum_preinitialization_amplitude_error",maximumInitialMomentError,
      2*std::numeric_limits<double>::epsilon(),"<=","dimensionless"});return result;
}
Result RejectedNativeInitialization() {
  Result missing;if(!Prerequisite(&missing))return missing;
  Checks check;
  {
    Fixture fixture(check);
    if(fixture.ready) {
      fixture.InitializePhysicalProfile();const auto before=fixture.Bytes();
      auto mismatched=fixture.configuration;mismatched.spectralBins++;
      A::NativeStateOwner rejected;
      check.Require(!rejected.InitializeFromIntegrated(mismatched).ok()&&!rejected.Acquire()&&SameBytes(before,fixture.Bytes()),
          "N mismatch touched native host arrays or published four-state authority");
      auto* energy=FL::FieldLinesAll[0].GetSegment(0)->GetDatum_ptr(W::CellIntegratedWaveEnergy);
      const double old=energy[0];if(PIC::ThisThread==0)energy[0]=-1;
      const auto invalid=rejected.InitializeFromIntegrated(fixture.configuration);energy[0]=old;
      check.Require(!invalid.ok()&&!rejected.Acquire()&&SameBytes(before,fixture.Bytes()),
          "one-rank invalid input escaped collective rejection or changed diagnostics");
      A::NativeStateOwner owner;
      check.Require(A::LegacyMutationPreflight("empty native authority",owner).ok(),
          "legacy capability rejected an absent four-spectrum authority");
      if(check.Status(owner.InitializeFromIntegrated(fixture.configuration))) {
        const auto committed=owner.Acquire();const auto image=fixture.Bytes();
        const auto capability=A::LegacyMutationPreflight("private native authority",owner);
        check.Require(!capability.ok()&&capability.code==Transport::StatusCode::UnsupportedConfiguration&&
            owner.Acquire()==committed&&SameBytes(image,fixture.Bytes()),
            "live four-spectrum authority allowed legacy mutation or capability query changed state");
        mismatched=fixture.configuration;mismatched.spectralKMinPerM*=0.01;
        check.Require(!owner.InitializeFromIntegrated(mismatched).ok()&&owner.Acquire()==committed&&SameBytes(image,fixture.Bytes()),
            "live four-state grid was silently remapped or diagnostics changed on failed reinitialization");
      }
      std::vector<double> destination={3,5,7};const auto sentinel=destination;
      mismatched=fixture.configuration;mismatched.spectralBins=W::IsotropicSEP::NK+1;
      check.Require(!A::CopyLegacyNativeSpectrum(mismatched,nullptr,2*W::IsotropicSEP::NK,&destination).ok()&&destination==sentinel,
          "legacy N>NK guard read null source or changed destination");
      mismatched=fixture.configuration;mismatched.spectralKMinPerM*=0.01;
      check.Require(!A::CopyLegacyNativeSpectrum(mismatched,nullptr,2*W::IsotropicSEP::NK,&destination).ok()&&destination==sentinel,
          "legacy physical k mismatch escaped pre-read rejection");
    }
  }
  return Finish(check,"malformed counts, one-rank invalid energy and live grid changes reject collectively before immutable-owner publication or native buffer mutation");
}
Result NativeRestore() {
  Result missing;if(!Prerequisite(&missing))return missing;
  Checks check;
  {
    Fixture fixture(check);
    if(fixture.ready) {
      fixture.InitializePhysicalProfile();A::NativeStateOwner owner;
      if(check.Status(owner.InitializeFromIntegrated(fixture.configuration))) {
        const auto captured=owner.Acquire();const auto background=Background::SnapshotStore::Instance().Current();
        std::vector<A::NativePhysicalCell> physical;
        if(check.Status(A::ReadNativePhysicalCells(*background,&physical))) {
          A::NativeRestorePlan plan;
          check.Status(owner.StageRestore(*captured,*background,physical,&plan));
          if(check.Status(owner.ConfirmCollectiveReady(plan))) {
            check.Require(owner.CommitRestore(&plan),"collectively ready native restore failed to commit");
            check.Require(!owner.CommitRestore(&plan),"consumed native restore plan committed twice");
            const auto restored=owner.Acquire();
            check.Require(restored.get()!=captured.get()&&restored->identity==captured->identity&&
                restored->cells[0].magneticEnergyJ==captured->cells[0].magneticEnergyJ,
                "restore reconstructed or aliased four independently owned arrays");
            const auto before=fixture.Bytes();auto malformed=*captured;
            malformed.cells[0].magneticEnergyJ[1][0]*=2;
            check.Require(!owner.StageRestore(malformed,*background,physical,&plan).ok()&&SameBytes(before,fixture.Bytes()),
                "changed four-array value was accepted under stale identity or touched diagnostics");
            auto wrongPhysical=physical;for(double& value:wrongPhysical[0].magneticFieldT)value*=-1;
            check.Require(!owner.StageRestore(*captured,*background,wrongPhysical,&plan).ok()&&SameBytes(before,fixture.Bytes()),
                "changed physical field orientation was accepted or changed native arrays");
            wrongPhysical=physical;
            const double shift=std::min(wrongPhysical[0].begin.numberDensityM3,wrongPhysical[0].end.numberDensityM3)*0.125;
            wrongPhysical[0].begin.numberDensityM3+=shift;wrongPhysical[0].end.numberDensityM3-=shift;
            check.Require(!owner.StageRestore(*captured,*background,wrongPhysical,&plan).ok()&&SameBytes(before,fixture.Bytes()),
                "equal-and-opposite endpoint density changes escaped the wave-input physical-image identity");
            PW::Initialization alternate;
            alternate.grid.centersPerM=captured->grid.centersPerM;alternate.grid.logWidth=captured->grid.logWidth;
            alternate.plus=captured->initialPlus;alternate.minus=captured->initialMinus;
            alternate.powerLawIndex=captured->initialPowerLawIndex;
            alternate.vacuumPermeabilityHPerM=2*captured->vacuumPermeabilityHPerM;
            alternate.epochS=captured->epochS;alternate.validUntilS=captured->validUntilS;
            alternate.geometryGeneration=captured->geometryGeneration;
            alternate.backgroundIdentity=captured->backgroundIdentity;
            alternate.configurationFingerprint=captured->configurationFingerprint;
            alternate.initializerProvenance=captured->initializerProvenance;
            for(const auto& cell:captured->cells) {
              PW::InitialCell input;input.fieldLine=cell.fieldLine;input.segment=cell.segment;
              input.volumeM3=cell.volumeM3;input.meanMagneticFieldT=cell.meanMagneticFieldT;
              input.massDensityKgM3=cell.massDensityKgM3;input.domain=cell.domain;
              const auto moments=PW::DeriveBranches(cell);input.ePlusJ=moments.ePlusJ;input.eMinusJ=moments.eMinusJ;
              alternate.cells.push_back(input);
            }
            const auto noncanonical=PW::PrepareInitialization(alternate);
            if(check.Status(noncanonical.status)) {
              check.Status(PW::ValidateSnapshot(*noncanonical.snapshot));
              check.Require(!owner.StageRestore(*noncanonical.snapshot,*background,physical,&plan).ok()&&SameBytes(before,fixture.Bytes()),
                  "internally valid but noncanonical native magnetic permeability was accepted");
            }

            check.Status(owner.StageRestore(*captured,*background,physical,&plan));
            auto changed=fixture.configuration;changed.cflSafety*=0.9;
            check.Status(Turbulence::SetActiveConfiguration(changed));
            const auto configReady=owner.ValidateReadyForCommit(plan);const bool configCommitted=owner.CommitRestore(&plan);
            check.Status(Turbulence::SetActiveConfiguration(fixture.configuration));
            check.Require(!configReady.ok()&&!configCommitted&&SameBytes(before,fixture.Bytes()),
                "configuration changed after staging but restore readiness/commit accepted it");

            check.Status(owner.StageRestore(*captured,*background,physical,&plan));
            auto* segment=FL::FieldLinesAll[0].GetSegment(0);
            char* original=segment->GetAssociatedDataBufferPointer();
            std::vector<char> relocated(original,original+segment->AssociatedDataLength());
            segment->SetAssociatedDataBufferPointer(relocated.data());
            const auto pointerReady=owner.ValidateReadyForCommit(plan);const bool pointerCommitted=owner.CommitRestore(&plan);
            segment->SetAssociatedDataBufferPointer(original);
            check.Require(!pointerReady.ok()&&!pointerCommitted&&SameBytes(before,fixture.Bytes()),
                "relocated native diagnostic buffer escaped staged-pointer readiness checks");

            check.Status(owner.StageRestore(*captured,*background,physical,&plan));
            {
              auto read=Background::SnapshotStore::Instance().BeginParticleRead(Background::SimulationTimeSeconds());
              check.Require(!owner.ValidateReadyForCommit(plan).ok()&&!owner.CommitRestore(&plan)&&SameBytes(before,fixture.Bytes()),
                  "native wave restore was accepted during an active particle read phase");
            }
            // Qualify Commit's own defenses independently: first obtain a
            // collectively ready plan, then mutate one prerequisite and call
            // Commit directly, without ValidateReady clearing readiness flags.
            check.Status(owner.StageRestore(*captured,*background,physical,&plan));
            if(check.Status(owner.ConfirmCollectiveReady(plan))) {
              check.Status(Turbulence::SetActiveConfiguration(changed));
              const bool committed=owner.CommitRestore(&plan);
              check.Status(Turbulence::SetActiveConfiguration(fixture.configuration));
              check.Require(!committed&&SameBytes(before,fixture.Bytes()),"ready-plan commit ignored changed configuration");
            }
            check.Status(owner.StageRestore(*captured,*background,physical,&plan));
            if(check.Status(owner.ConfirmCollectiveReady(plan))) {
              segment->SetAssociatedDataBufferPointer(relocated.data());
              const bool committed=owner.CommitRestore(&plan);segment->SetAssociatedDataBufferPointer(original);
              check.Require(!committed&&SameBytes(before,fixture.Bytes()),"ready-plan commit ignored relocated native buffers");
            }
            check.Status(owner.StageRestore(*captured,*background,physical,&plan));
            if(check.Status(owner.ConfirmCollectiveReady(plan))) {
              auto read=Background::SnapshotStore::Instance().BeginParticleRead(Background::SimulationTimeSeconds());
              check.Require(!owner.CommitRestore(&plan)&&SameBytes(before,fixture.Bytes()),"ready-plan commit wrote during particle read phase");
            }
            // An owner address may be reused. Its independent lifetime identity
            // must still reject the previous owner's staged restore plan.
            typename std::aligned_storage<sizeof(A::NativeStateOwner),alignof(A::NativeStateOwner)>::type storage;
            auto* first=new(&storage) A::NativeStateOwner;
            A::NativeRestorePlan previousLifetime;
            check.Status(first->StageRestore(*captured,*background,physical,&previousLifetime));
            check.Status(first->ConfirmCollectiveReady(previousLifetime));
            first->~NativeStateOwner();auto* second=new(&storage) A::NativeStateOwner;
            const bool lifetimeCommitted=second->CommitRestore(&previousLifetime);
            check.Require(!lifetimeCommitted&&!second->ValidateReadyForCommit(previousLifetime).ok(),
                "restore plan survived destruction and same-address reconstruction of its owner");
            second->~NativeStateOwner();
            check.Require(SameBytes(before,fixture.Bytes()),"rejected readiness fixtures changed native buffers");
            const auto diagnostics=owner.Diagnostics();
            check.Require(diagnostics.mpiIdentityAgrees&&diagnostics.readbackAgrees,
                "restore reported unmeasured MPI identity or native readback status");
          }
        }
      }
    }
  }
  return Finish(check,"native staging preserves exact four arrays; collective readiness measures MPI identity, rejects stale configuration/storage/owner/read-phase and binds endpoint primitives before allocation-free publication");
}
Result PhysicalInitializerRejection() {
  Result missing;if(!Prerequisite(&missing))return missing;
  Checks check;
  {
    Fixture fixture(check);
    if(fixture.ready) {
      fixture.InitializePhysicalProfile();const auto before=fixture.Bytes();
      auto initialize=[&](double first,double last,double power) {
        return W::InitializeWaveEnergyFromPhysicalParametersChecked(W::CellIntegratedWaveEnergy,
            W::WaveEnergyConstants::TYPICAL_B0_1AU,first,last,power,false,false);
      };
      check.Require(!initialize(-1,0.4,-2).ok()&&SameBytes(before,fixture.Bytes()),
          "negative initial amplitude was clipped or changed host values");
      const int lastIndex=FL::FieldLinesAll[0].GetTotalSegmentNumber()-1;
      auto* last=FL::FieldLinesAll[0].GetSegment(lastIndex)->GetEnd();
      double* field=last->GetDatum_ptr(FL::DatumAtVertexMagneticField);
      const double oldField[3]={field[0],field[1],field[2]};
      if(PIC::ThisThread==0)std::fill(field,field+3,0);
      const auto invalidB=initialize(0.4,0.4,-2);std::copy(oldField,oldField+3,field);
      check.Require(!invalidB.ok()&&SameBytes(before,fixture.Bytes()),
          "invalid last-cell B selected a fallback or left earlier initialized cells changed");
      double oldX[3],beginX[3];last->GetX(oldX);FL::FieldLinesAll[0].GetFirstVertex()->GetX(beginX);
      if(PIC::ThisThread==0)last->SetX(beginX);
      const auto invalidRadius=initialize(0.4,0.4,-2);last->SetX(oldX);
      check.Require(!invalidRadius.ok()&&SameBytes(before,fixture.Bytes()),
          "degenerate radial endpoints acquired a fabricated offset or changed host values");
      check.Require(!initialize(0.4,0.4,std::numeric_limits<double>::quiet_NaN()).ok()&&SameBytes(before,fixture.Bytes()),
          "nonfinite profile exponent was replaced by a fallback");
      check.Status(initialize(0,0,-2));
      for(int line=0;line<FL::nFieldLine;++line)for(int i=0;i<FL::FieldLinesAll[line].GetTotalSegmentNumber();++i) {
        const auto* energy=FL::FieldLinesAll[line].GetSegment(i)->GetDatum_ptr(W::CellIntegratedWaveEnergy);
        check.Require(energy[0]==0&&energy[1]==0,"explicit zero physical amplitude did not remain zero");
      }
    }
  }
  return Finish(check,"physical initializer rejects malformed amplitude, terminal B and radial geometry collectively before any host write, while explicit zero amplitude remains physical zero");
}
Descriptor Make(const char* id,const char* name,TestCallback callback) {
  Descriptor descriptor;descriptor.id=id;descriptor.name=name;descriptor.description=name;
  descriptor.group="polarized-wave-native";descriptor.initialization=InitializationLevel::FieldLineModel;
  descriptor.runtime=RuntimeClass::Routine;descriptor.seedCapability=SeedCapability::Deterministic;
  descriptor.repeatSupported=true;descriptor.seedPolicy="deterministic exact moments and explicit physical profile";
  descriptor.supportedBuildModes="native spectral model with resolved diagnostic layout; run with --turbulence-model wave-number-resolved";
  descriptor.stateIsolation="private immutable owners; all touched native diagnostic arrays/configuration restored";
  descriptor.callback=callback;return descriptor;
}
} // namespace
std::vector<Descriptor> PolarizedNativeDescriptors() {
  return {Make("PNW01","Native four-spectrum moment/bin/MPI readback",NativeMoments),
          Make("PNW02","Collective native four-state prepublication rejection",RejectedNativeInitialization),
          Make("PNW03","Exact native four-array staged restore contract",NativeRestore),
          Make("PNW04","Physical wave initializer transactional rejection",PhysicalInitializerRejection)};
}
} }
