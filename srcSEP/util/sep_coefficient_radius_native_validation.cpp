#include "sep_coefficient_radius_native_validation.h"
#include "../coefficient_providers.h"
#include "sep_initialization.h"
#include "sep_background_runtime.h"
#include "sep_provider_radius_frame.h"
#include "sep_native_startup_input.h"
#include "../adapters/parallel_diffusion_adapter.h"
#include "sep_common_header_path.h"
#include SRCSEP_SEP_COMMON_HEADER(parallel_diffusion/parallel_diffusion.h)
#include <array>
#include <cmath>
#include <exception>
#include <limits>
#include <memory>

namespace SEP { namespace Testing { namespace {
namespace T=Transport;
namespace A=T::PICAdapter;
namespace C=T::Coefficient;
namespace FL=PIC::FieldLine;
struct Checks {
  unsigned failures=0;std::string message;
  void Require(bool condition,const std::string& reason) {
    if(condition)return;
    ++failures;if(!message.empty())message+="; ";message+=reason;
  }
  void Near(double actual,double expected,double tolerance,const char* reason) {
    Require(std::isfinite(actual)&&std::isfinite(expected)&&expected>0&&std::abs(actual/expected-1.0)<=tolerance,reason);
  }
};
Result NativeRadius() {
  Result unavailable;
  if(!Initialization::HasActive()) {
    unavailable.status=Status::Skip;unavailable.message="RADIUSN01 requires opted-in schema-6 --test-input startup; run test/run_native_coefficient_radius_tests.py";return unavailable;
  }
  Checks check;const auto& initialization=Initialization::Active();
  const auto immutableFingerprint=Initialization::Fingerprint(initialization);
  const auto background=Background::SnapshotStore::Instance().Current();
  int species=-1;for(int i=0;i<PIC::nTotalSpecies;++i)
    if(PIC::MolecularData::GetMass(i)>0&&PIC::MolecularData::GetElectricCharge(i)!=0){species=i;break;}
  if(Background::ConfiguredProvider()!=Background::Provider::Analytic||!background||
      !FL::FieldLinesAll||FL::nFieldLine<1||!FL::FieldLinesAll[0].GetFirstSegment()||species<0) {
    unavailable.status=Status::Skip;unavailable.message="requires an initialized analytic native line and a charged species";return unavailable;
  }
  const double origin[3]={initialization.parkerOriginM.x,initialization.parkerOriginM.y,initialization.parkerOriginM.z};
  check.Require(std::hypot(origin[0],std::hypot(origin[1],origin[2]))>0&&initialization.solarWindSpeedMPerS!=400000.0,
      "isolated radius fixture must declare translated origin and nondefault wind in its immutable input");
  const long int particleCount=PIC::ParticleBuffer::GetAllPartNum();
  const auto diagnostics=A::GetCoefficientProviderDiagnostics();
  double measured=0,expected=0,oldGlobal=0,radius=0,globalRadius=0,divergence=0;
  {
    struct Restore {
      Checks& check;C::Configuration original=C::ActiveConfiguration();
      explicit Restore(Checks& c):check(c){}
      ~Restore(){const auto status=C::SetActiveConfiguration(original);check.Require(status.ok(),"cannot restore original coefficient configuration");}
    } restore(check);
    // Explicit controlled coefficient law; this is not a mutation of the
    // immutable startup or a new plasma/source assumption. QLT1's closed
    // correlation-length law isolates r from numerical pitch quadrature.
    auto configuration=restore.original;configuration.source=C::SourceMode::Prescribed;
    configuration.meanFreePath=C::MeanFreePathKind::Qlt1;
    configuration.prescribedDeltaBOverB=0.05;configuration.correlationLengthAt1AuM=1e9;
    configuration.invalidPolicy=C::InvalidPolicy::Fail;
    const auto configured=C::SetActiveConfiguration(configuration);
    check.Require(configured.ok(),configured.message);
    if(configured.ok()) {
      A::ParticleContext context;context.segment=FL::FieldLinesAll[0].GetFirstSegment();
      context.state.fieldLineId=0;context.state.coordinate=0.5;context.state.species=species;
      context.state.massKg=PIC::MolecularData::GetMass(species);
      context.snapshotGeneration=background->field_line_generation();context.snapshotFromS=background->valid_from_seconds();
      context.snapshotUntilS=background->valid_until_seconds();context.particleStepEpochS=background->epoch_seconds();
      context.stepPrepared=true;context.moveDurationS=0.0;
      double position[3];FL::FieldLinesAll[0].GetCartesian(position,context.state.coordinate);
      radius=std::hypot(position[0]-origin[0],std::hypot(position[1]-origin[1],position[2]-origin[2]));
      globalRadius=std::hypot(position[0],std::hypot(position[1],position[2]));
      const double speed=1e6,c=299792458.0;
      const double momentum=context.state.massKg*speed/std::sqrt((1.0-speed/c)*(1.0+speed/c));
      A::LocalBackgroundView sampled;
      {
        auto read=Background::SnapshotStore::Instance().BeginParticleRead(background->epoch_seconds());
        const auto local=A::EvaluateLocalBackground(context,&sampled);check.Require(local.ok(),local.message);
        if(local.ok()) {
          check.Require(sampled.heliocentricRadiusKnown&&
              sampled.radiusFrameIdentity==Background::RadiusFrame::ContractIdentity(Background::Provider::Analytic),
              "native background did not declare its analytic coefficient radial frame");
          check.Near(sampled.heliocentricRadiusM,radius,2e-14,"native background radius differs from independently translated native position");
          for(unsigned d=0;d<3;++d)check.Require(std::isfinite(sampled.heliocentricPositionM[d])&&
              std::abs(sampled.heliocentricPositionM[d]-(position[d]-origin[d]))<=2e-14*radius,
              "native background solar-relative vector differs from independently translated native position");
          divergence=sampled.velocityDivergencePerS;
          check.Near(divergence,2.0*initialization.solarWindSpeedMPerS/radius,2e-12,"native gradient and configured-origin radius do not describe one declared radial flow");
          A::PICMeanFreePathProvider provider(context,sampled.identity);
          const auto sample=provider.Evaluate(0.0,momentum,0.5);check.Require(sample.status.ok(),sample.status.message);
          check.Require(sample.provenance.find("radius-frame="+sampled.radiusFrameIdentity)!=std::string::npos&&
              sample.provenance.find("heliocentric-radius-m=")!=std::string::npos,
              "coefficient provenance omits the actual provider radius/frame contract");
          measured=sample.lambdaParallelM;
          // Independent declared QLT1 law: Lc scales linearly with the
          // PROVIDER radius and lambda=Lc/r^2_amplitude*(rL/Lc)^(1/3).
          // Neither production coefficient kernel nor native radius helper
          // supplies these expected values. Charge/mass come from PIC's actual
          // compiled species table; no inferred isotope metadata is used.
          const double charge=std::abs(PIC::MolecularData::GetElectricCharge(species));
          const double larmor=std::sqrt(2.0/3.0)*momentum/(charge*sampled.magneticFieldMagnitudeT);
          auto law=[&](double r) {const double correlation=configuration.correlationLengthAt1AuM*r/configuration.spectrum.referenceRadiusM;
            return correlation/(configuration.prescribedDeltaBOverB*configuration.prescribedDeltaBOverB)*std::pow(larmor/correlation,1.0/3.0);};
          expected=law(radius);oldGlobal=law(globalRadius);
          check.Near(measured,expected,2e-12,"actual QLT1 provider used the global-zero radius instead of the immutable analytic origin");
          check.Require(std::abs(oldGlobal/expected-1.0)>0.1,"global-zero negative control did not materially distinguish the requested radius law");
        }
      }
      // A controlled tensor replaces only grad(U). The radial frame must be
      // established before that branch, without imposing fresh endpoint-U
      // validation or substituting a different origin for existing fixtures.
      const std::array<double,9> zero{{0,0,0,0,0,0,0,0,0}};
      {
        A::ScopedVelocityGradientFixture tensor(zero);
        auto read=Background::SnapshotStore::Instance().BeginParticleRead(background->epoch_seconds());
        A::LocalBackgroundView controlled;const auto status=A::EvaluateLocalBackground(context,&controlled);
        check.Require(status.ok()&&controlled.heliocentricRadiusKnown,"explicit derivative fixture lost its provider-defined radius");
        if(status.ok())check.Near(controlled.heliocentricRadiusM,radius,2e-14,"derivative fixture changed the provider radial origin");
      }
      // Raw imported coordinates plus an explicit tensor STILL do not locate
      // their physical Sun. This synthetic metadata-only case must reject the
      // radius-dependent coefficient rather than inherit the configured
      // analytic origin from the surrounding process. All native arrays and
      // the original immutable background pointer are restored unchanged.
      {
        struct SnapshotRestore {
          Checks& check;std::shared_ptr<const Background::BackgroundSnapshot> original;
          SnapshotRestore(Checks& c):check(c),original(Background::SnapshotStore::Instance().Current()){}
          ~SnapshotRestore() {
            try {auto& store=Background::SnapshotStore::Instance();store.ResetForTests();
              store.Publish(*original);store.RestoreCheckpoint(original);}
            catch(const std::exception&){check.Require(false,"cannot restore background after unknown imported-radius control");}
          }
        } restoreSnapshot(check);
        auto& store=Background::SnapshotStore::Instance();store.ResetForTests();
        store.Publish(Background::BackgroundSnapshot(Background::Provider::Swmf,Background::Ownership::ImportedReadOnly,
            background->epoch_seconds(),background->epoch_seconds(),std::numeric_limits<double>::infinity(),1,
            "controlled-imported-unmapped-radius","explicit tensor with unknown imported radial origin"));
        auto unknown=context;unknown.snapshotGeneration=1;
        A::ScopedVelocityGradientFixture tensor(zero);
        {auto read=store.BeginParticleRead(background->epoch_seconds());
         A::LocalBackgroundView imported;const auto status=A::EvaluateLocalBackground(unknown,&imported);
         check.Require(status.ok()&&!imported.heliocentricRadiusKnown,"imported coordinate array acquired an inferred coefficient origin");
         A::PICMeanFreePathProvider provider(unknown,"controlled-imported-unmapped-radius");
         const auto sample=provider.Evaluate(0.0,momentum,0.5);
         check.Require(sample.status.code==T::StatusCode::UnsupportedConfiguration&&sample.status.message.find("provider radial frame")!=std::string::npos,
             "unknown imported radius did not reject the actual coefficient consumer explicitly");}
      }
    }
  }
  const auto after=A::GetCoefficientProviderDiagnostics();
  check.Require(Initialization::Fingerprint(Initialization::Active())==immutableFingerprint&&
      Background::SnapshotStore::Instance().Current()==background&&PIC::ParticleBuffer::GetAllPartNum()==particleCount,
      "read-only coefficient fixture changed startup/provider/particle ownership");
  check.Require(after.invalidSamples==diagnostics.invalidSamples&&after.ballisticSubstitutions==diagnostics.ballisticSubstitutions&&
      after.amplitudeRegularizations==diagnostics.amplitudeRegularizations&&after.quadratureEvaluations==diagnostics.quadratureEvaluations,
      "closed QLT1 radius fixture unexpectedly changed coefficient diagnostics");
  unsigned failures=0;MPI_Allreduce(&check.failures,&failures,1,MPI_UNSIGNED,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
  Result result;result.status=failures?Status::Fail:Status::Pass;result.message=failures?check.message:
      "actual immutable translated-origin/nondefault-wind startup drives the independently normalized radius-dependent QLT1 provider";
  result.configuration={"evidence_level=actual-native-runtime-input-and-coefficient-provider",
      "initialization_override=false;initialization_source=opted-in-schema6--test-input;source_and_host_steps=none",
      "expected=explicit-QLT1-correlation-law;old-global-zero-radius=independent-negative-control",
      "full-production-material-grid-imported-provider-and-source-qualified=false"};
  result.metrics={{"assertion_failures",static_cast<double>(failures),0,"<=","count"},
      {"provider_radius_expected",radius,0,">","m"},{"old_global_radius",globalRadius,0,">","m"},
      {"lambda_actual",measured,expected,"reference","m"},{"lambda_expected",expected,0,">","m"},
      {"lambda_old_global_negative",oldGlobal,expected,"different","m"},
      {"declared_wind_speed",initialization.solarWindSpeedMPerS,400000,"different","m/s"},
      {"native_divergence",divergence,0,">","1/s"}};
  return result;
}
// This separate process fixture exercises the optional library through the
// actual Parker-facing provider. Its expected radial law is evaluated here
// independently; no library evaluator, background radius helper or mutable
// initialization override supplies the reference. The full schema-6 input
// explicitly installs this model before native allocation.
Result NativeParallelRadius() {
  namespace PD=SEP::ParallelDiffusion;
  Result unavailable;
  if(!Initialization::HasActive()||!ParallelDiffusionBinding::IsInstalled()||
      !PD::HasActiveConfiguration()) {
    unavailable.status=Status::Skip;unavailable.message="RADIUSN02 requires opted-in schema-6 power_law_lambda startup; run radius runner --case RADIUSN02";return unavailable;
  }
  Checks check;const auto& initialization=Initialization::Active();
  const auto immutableFingerprint=Initialization::Fingerprint(initialization);
  const auto background=Background::SnapshotStore::Instance().Current();
  int species=-1;for(int i=0;i<PIC::nTotalSpecies;++i)
    if(PIC::MolecularData::GetMass(i)>0&&PIC::MolecularData::GetElectricCharge(i)!=0){species=i;break;}
  if(Background::ConfiguredProvider()!=Background::Provider::Analytic||!background||
      !FL::FieldLinesAll||FL::nFieldLine<1||!FL::FieldLinesAll[0].GetFirstSegment()||species<0) {
    unavailable.status=Status::Skip;unavailable.message="requires actual analytic native line and charged species";return unavailable;
  }
  const auto library=PD::GetActiveConfiguration();
  const auto libraryFingerprint=PD::ConfigurationFingerprint(library);
  const auto& power=library.powerLawLambda;
  check.Require(library.model==PD::ModelId::PowerLawLambda&&power.independentVariable==PD::IndependentVariable::Speed&&
      power.lambda0M==3e9&&power.independentReferenceSI==2e6&&power.independentExponent==0.75&&
      power.useRadialFactor&&power.radius0M==1.495978707e11&&power.radialExponent==1.5&&!power.useFieldFactor,
      "immutable library fixture did not install its explicitly declared speed/radius power law");
  check.Require(C::ActiveConfiguration().spatial==C::SpatialKind::ParallelDiffusionLibrary&&
      C::ActiveConfiguration().source==C::SourceMode::Prescribed,
      "actual Parker coefficient selection differs from the explicit library/prescribed fixture");
  const double origin[3]={initialization.parkerOriginM.x,initialization.parkerOriginM.y,initialization.parkerOriginM.z};
  check.Require(std::hypot(origin[0],std::hypot(origin[1],origin[2]))>0&&initialization.solarWindSpeedMPerS==725000,
      "isolated library fixture lost its translated origin/nondefault wind");
  const auto particleCount=PIC::ParticleBuffer::GetAllPartNum();
  const auto diagnostics=A::GetCoefficientProviderDiagnostics();
  A::ParticleContext context;context.segment=FL::FieldLinesAll[0].GetFirstSegment();
  context.state.fieldLineId=0;context.state.coordinate=0.5;context.state.species=species;
  context.state.massKg=PIC::MolecularData::GetMass(species);
  context.snapshotGeneration=background->field_line_generation();context.snapshotFromS=background->valid_from_seconds();
  context.snapshotUntilS=background->valid_until_seconds();context.particleStepEpochS=background->epoch_seconds();
  context.stepPrepared=true;context.moveDurationS=0;
  double position[3];FL::FieldLinesAll[0].GetCartesian(position,context.state.coordinate);
  const double radius=std::hypot(position[0]-origin[0],std::hypot(position[1]-origin[1],position[2]-origin[2]));
  const double globalRadius=std::hypot(position[0],std::hypot(position[1],position[2]));
  const double speed=1e6;
  // The two physical factors and kappa=v*lambda/3 are the declared model,
  // independently applied to the native position in the input's frame.
  const double lambda=3e9*std::pow(speed/2e6,0.75)*std::pow(radius/1.495978707e11,1.5);
  const double expected=speed*lambda/3.0;
  const double oldGlobal=expected*std::pow(globalRadius/radius,1.5);
  double measured=0,derivative=0;
  {
    auto read=Background::SnapshotStore::Instance().BeginParticleRead(background->epoch_seconds());
    A::PICSpatialDiffusionProvider provider(context,"explicit-runtime-power-law-radius");
    const auto sample=provider.Evaluate(0,speed);check.Require(sample.status.ok(),sample.status.message);
    if(sample.status.ok()) {
      measured=sample.kappaParallelM2PerS;derivative=sample.dKappaParallelDsMPerS;
      check.Near(measured,expected,2e-12,"actual parallel-diffusion library received global coordinates instead of solar-relative provider coordinates");
      check.Require(std::isfinite(derivative),"actual library's refined along-line derivative is nonfinite");
    }
  }
  check.Require(std::abs(oldGlobal/expected-1.0)>0.1,"library global-zero negative control did not distinguish the translated radial factor");
  // An imported Cartesian record and explicit tensor still do not define a
  // solar origin. Even this installed radius-dependent library must reject
  // before evaluating its center/stencil in an invented coordinate frame.
  {
    struct Restore {
      Checks& check;std::shared_ptr<const Background::BackgroundSnapshot> original;
      explicit Restore(Checks& c):check(c),original(Background::SnapshotStore::Instance().Current()){}
      ~Restore(){try{auto& store=Background::SnapshotStore::Instance();store.ResetForTests();
        store.Publish(*original);store.RestoreCheckpoint(original);}
        catch(const std::exception&){check.Require(false,"cannot restore background after library unknown-origin control");}}
    } restore(check);
    auto& store=Background::SnapshotStore::Instance();store.ResetForTests();
    store.Publish(Background::BackgroundSnapshot(Background::Provider::Swmf,Background::Ownership::ImportedReadOnly,
        background->epoch_seconds(),background->epoch_seconds(),std::numeric_limits<double>::infinity(),1,
        "controlled-imported-unmapped-library-radius","explicit tensor with unknown imported solar origin"));
    auto unknown=context;unknown.snapshotGeneration=1;
    const std::array<double,9> zero{{0,0,0,0,0,0,0,0,0}};A::ScopedVelocityGradientFixture tensor(zero);
    {auto read=store.BeginParticleRead(background->epoch_seconds());
     A::PICSpatialDiffusionProvider provider(unknown,"controlled-unknown-library-radius");
     const auto sample=provider.Evaluate(0,speed);
     check.Require(sample.status.code==T::StatusCode::UnsupportedConfiguration&&sample.status.message.find("provider radial frame")!=std::string::npos,
         "unknown imported radius did not reject the actual library consumer with typed frame context");}
  }
  const auto after=A::GetCoefficientProviderDiagnostics();
  check.Require(Initialization::Fingerprint(Initialization::Active())==immutableFingerprint&&
      PD::ConfigurationFingerprint(PD::GetActiveConfiguration())==libraryFingerprint&&
      Background::SnapshotStore::Instance().Current()==background&&PIC::ParticleBuffer::GetAllPartNum()==particleCount,
      "read-only library radius fixture changed startup/library/background/particle ownership");
  check.Require(after.invalidSamples==diagnostics.invalidSamples&&after.ballisticSubstitutions==diagnostics.ballisticSubstitutions&&
      after.amplitudeRegularizations==diagnostics.amplitudeRegularizations&&after.quadratureEvaluations==diagnostics.quadratureEvaluations,
      "library radius fixture unexpectedly changed coefficient diagnostics");
  unsigned failures=0;MPI_Allreduce(&check.failures,&failures,1,MPI_UNSIGNED,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
  Result result;result.status=failures?Status::Fail:Status::Pass;result.message=failures?check.message:
      "actual immutable power-law library consumes translated native provider coordinates and rejects unknown imported origins";
  result.configuration={"evidence_level=actual-native-runtime-input-and-Parker-library-provider",
      "initialization_override=false;initialization_source=opted-in-schema6--test-input;source_and_host_steps=none",
      "expected=independent-speed-and-radius-power-law;old-global-zero=independent-negative-control",
      "full-production-material-grid-imported-provider-and-source-qualified=false"};
  result.metrics={{"assertion_failures",static_cast<double>(failures),0,"<=","count"},
      {"provider_radius_expected",radius,0,">","m"},{"old_global_radius",globalRadius,0,">","m"},
      {"kappa_actual",measured,expected,"reference","m^2/s"},{"kappa_expected",expected,0,">","m^2/s"},
      {"kappa_old_global_negative",oldGlobal,expected,"different","m^2/s"},
      {"native_refined_kappa_derivative",derivative,0,"finite","m/s"}};
  return result;
}
} // namespace
std::vector<Descriptor> CoefficientRadiusNativeDescriptors() {
  Descriptor descriptor;descriptor.id="RADIUSN01";descriptor.name="Immutable translated-origin native coefficient radius";
  descriptor.description=descriptor.name;descriptor.group="coefficient-radius-native";
  descriptor.initialization=InitializationLevel::FieldLineModel;descriptor.runtime=RuntimeClass::Extended;
  descriptor.seedCapability=SeedCapability::Deterministic;descriptor.repeatSupported=true;descriptor.seedPolicy="actual declared runtime geometry and explicit QLT1 parameters; no draws";
  descriptor.stateIsolation="separate opted-in schema-6 --test-input process; read-only geometry/background; typed coefficient configuration restored";
  descriptor.supportedBuildModes="initialized analytic model with nonzero runtime origin/nondefault wind; no mutable initialization test override";
  descriptor.tags.push_back(NativeStartupInput::DescriptorTag());
  descriptor.callback=NativeRadius;
  auto parallel=descriptor;parallel.id="RADIUSN02";parallel.name="Immutable translated-origin native parallel-library radius";
  parallel.description=parallel.name;parallel.callback=NativeParallelRadius;
  parallel.seedPolicy="actual installed speed/radius power-law and declared native geometry; no draws";
  parallel.stateIsolation="separate tagged schema6 startup; immutable library/native geometry; metadata-only control restored";
  return {descriptor,parallel};
}
} }
