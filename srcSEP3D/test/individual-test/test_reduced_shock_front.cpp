// Reduced shock-front application-boundary checks.
//
// These tests deliberately instantiate the production srcSEP3D factory and
// adapter rather than calling only the shared analytical provider.  They are
// still portable prerequisites: native owner/received-ghost and MPI evidence
// is collected separately by RSH25--RSH28 in the linked AMPS executable.

#include "sep3d_test_registry.h"

#include "background_factory.h"
#include "background_snapshot.h"
#include "configuration_io.h"
#include "run_configuration.h"
#include "shock_front_background_adapter.h"
#include "provider.h"
#include "reduced_front_output.h"

#include <cmath>
#include <filesystem>
#include <memory>
#include <string>
#include <vector>

namespace {

namespace A=SEP3D::Adapters;
namespace B=SEP3D::Background;
namespace C=SEP3D::Core;
namespace R=SEP3D::RuntimeModel;
namespace SF=SEP::CoronaSwcme::ShockFront;
using Result=SEP3D::Testing::Result;

Result Finish(bool good,const std::string& message) {
  Result result;
  result.status=good?SEP3D::Testing::Status::Pass:SEP3D::Testing::Status::Fail;
  result.message=message;return result;
}

C::Status Fixture(std::shared_ptr<const R::RunConfiguration3D>* configuration,
    std::shared_ptr<B::BackgroundProvider>* provider) {
  R::RunConfiguration3DOptions options;
  options.inputSchemaVersion=4;
  options.intent=R::RunIntent::TransportOnly;
  options.background=R::BackgroundAuthority::RuntimeModel;
  options.backgroundModelId="sep-corona-swcme-shock-front-v1";
  std::filesystem::path event;
  for(const auto& base:{std::filesystem::path("."),std::filesystem::path("srcSEP3D"),
      std::filesystem::path("../srcSEP3D")}) {
    const auto candidate=base/"examples/shock-front/handoff_smoke.event";
    if(std::filesystem::exists(candidate)) {event=std::filesystem::canonical(candidate);break;}
  }
  if(event.empty())return C::Status::Error(
      "cannot locate the reduced handoff event from this working directory");
  options.backgroundModelAssetPath=event.string();
  options.shock=R::ShockAuthority::None;
  options.source.enabled=false;
  options.memoryModel.particlesPerCell=0.0;
  options.populationControl=R::PopulationControlMode::Off;
  auto status=R::RunConfiguration3D::Create(options,configuration);
  if(status.ok())status=R::CreateBackgroundProvider(**configuration,provider);
  return status;
}

Result FactoryAndAmbient() {
  std::shared_ptr<const R::RunConfiguration3D> configuration;
  std::shared_ptr<B::BackgroundProvider> provider;
  auto status=Fixture(&configuration,&provider);
  if(!status.ok())return Finish(false,status.message);
  status=provider->Prepare(0.0);
  if(!status.ok())return Finish(false,status.message);

  // Query two widely separated physical regions.  The reduced contract says
  // both are canonical ambient references; an interior point must never be
  // relabelled as a sheath/ejecta value merely because it lies behind the
  // prescribed front.
  // The synthetic map is a +Z dipole; sample its open polar tube.  Its
  // equatorial exterior is a genuine magnetic null and is tested as an
  // explicit failure by the shared suite rather than floored here.
  const auto low=provider->Evaluate({0,0,2.0*C::Const::R_sun});
  const auto oneAu=provider->Evaluate({0,0,C::Const::AU});
  const auto* metadata=provider->PreparedMetadata();
  const auto reduced=std::dynamic_pointer_cast<A::ShockFrontBackgroundAdapter>(provider);
  const auto epoch=reduced?reduced->FrontEpoch():nullptr;
  const bool good=configuration->options().source.enabled==false&&
      configuration->options().shock==R::ShockAuthority::None&&low.valid&&
      oneAu.valid&&low.numberDensityM3>0&&oneAu.numberDensityM3>0&&
      low.B.Norm()>0&&oneAu.B.Norm()>0&&metadata&&metadata->generation==1&&
      metadata->provider==B::ProviderKind::RuntimeModel&&epoch&&
      epoch->generation==1&&!provider->ResolvedManifest().empty();
  return Finish(good,"production factory publishes one zero-particle ambient-reference epoch without a downstream-volume claim");
}

Result EpochAndRollback() {
  // A native particle CFL can be far shorter than run.background_dt_s.  Each
  // forward surface state must nevertheless receive a fresh generation; the
  // nominal cadence is only an absolute lower bound on that identity.
  {
    std::shared_ptr<const R::RunConfiguration3D> fineConfiguration;
    std::shared_ptr<B::BackgroundProvider> fineProvider;
    auto fineStatus=Fixture(&fineConfiguration,&fineProvider);
    if(!fineStatus.ok()||!fineProvider->Prepare(0.0).ok()||
        fineProvider->PreparedMetadata()->generation!=1||
        !fineProvider->Prepare(2.0).ok()||
        fineProvider->PreparedMetadata()->generation!=2)
      return Finish(false,"sub-cadence forward epochs reused one generation");
  }
  std::shared_ptr<const R::RunConfiguration3D> configuration;
  std::shared_ptr<B::BackgroundProvider> provider;
  auto status=Fixture(&configuration,&provider);
  if(!status.ok())return Finish(false,status.message);
  status=provider->Prepare(60.0);if(!status.ok())return Finish(false,status.message);
  const auto reduced=std::dynamic_pointer_cast<A::ShockFrontBackgroundAdapter>(provider);
  if(!reduced)return Finish(false,"factory did not return the reduced adapter");
  const auto before=reduced->FrontEpoch();
  if(!before)return Finish(false,"front epoch was not committed");

  B::BackgroundSnapshotBuilder builder;
  std::shared_ptr<const B::BackgroundSnapshot> first,second;
  status=builder.Build(*provider,{{0,0,2*C::Const::R_sun},{0,0,C::Const::AU}},&first);
  if(status.ok())status=provider->Prepare(120.0);
  if(status.ok())status=builder.Build(*provider,{{0,0,2*C::Const::R_sun},{0,0,C::Const::AU}},&second);
  if(!status.ok())return Finish(false,status.message);
  const auto after=reduced->FrontEpoch();
  const auto committed=after;
  const auto rejected=provider->Prepare(119.0); // same rounded generation, older epoch

  const bool good=before->trajectory.phase==SF::Phase::CoronalHistory&&after&&
      after->trajectory.phase==SF::Phase::SwcmeOuter&&
      before->eventIdentity==after->eventIdentity&&
      before->generation==2&&after->generation==3&&
      first->metadata().generation==2&&second->metadata().generation==3&&
      !rejected.ok()&&reduced->FrontEpoch()==committed&&
      reduced->FrontEpoch()->generation==3&&
      first->samples()[0].generation==2&&second->samples()[0].generation==3;
  return Finish(good,"pre/post-handoff epochs retain one authority and a rejected candidate preserves the committed front and ambient snapshots");
}

Result DeckPathResolution() {
  // Resolve the same translated deck once from the AMPS root and once from
  // srcSEP3D.  background.model_asset is relative to the deck, never the
  // process working directory; both parses must therefore freeze the exact
  // same canonical event path and physics identity.
  const auto original=std::filesystem::current_path();
  struct Restore { std::filesystem::path path; ~Restore(){
    std::error_code ignored;std::filesystem::current_path(path,ignored);} } restore{original};
  std::filesystem::path root=original;
  if(!std::filesystem::exists(root/"CODEX_REDUCED_SHOCK_TASK.txt"))root=original.parent_path();
  if(!std::filesystem::exists(root/"CODEX_REDUCED_SHOCK_TASK.txt"))
    return Finish(false,"cannot locate AMPS root for deck-path test");
  R::RunConfiguration3DOptions fromRoot,fromApplication;
  R::RunConfiguration3DOptions longFromRoot,longFromApplication;
  std::filesystem::current_path(root);
  auto status=R::LoadConfigurationFile(
      "srcSEP3D/examples/shock-front/handoff_smoke.in",&fromRoot);
  if(!status.ok())return Finish(false,"root parse: "+status.message);
  status=R::LoadConfigurationFile(
      "srcSEP3D/examples/shock-front/corona_to_1au.in",&longFromRoot);
  if(!status.ok())return Finish(false,"long root parse: "+status.message);
  std::filesystem::current_path(root/"srcSEP3D");
  status=R::LoadConfigurationFile("examples/shock-front/handoff_smoke.in",
      &fromApplication);
  if(!status.ok())return Finish(false,"application parse: "+status.message);
  status=R::LoadConfigurationFile("examples/shock-front/corona_to_1au.in",
      &longFromApplication);
  if(!status.ok())return Finish(false,"long application parse: "+status.message);
  std::shared_ptr<const R::RunConfiguration3D> a,b,longA,longB;
  status=R::RunConfiguration3D::Create(fromRoot,&a);
  if(status.ok())status=R::RunConfiguration3D::Create(fromApplication,&b);
  if(status.ok())status=R::RunConfiguration3D::Create(longFromRoot,&longA);
  if(status.ok())status=R::RunConfiguration3D::Create(longFromApplication,&longB);
  const bool good=status.ok()&&fromRoot.backgroundModelAssetPath==
      fromApplication.backgroundModelAssetPath&&a->physics_fingerprint()==
      b->physics_fingerprint()&&longFromRoot.backgroundModelAssetPath==
      longFromApplication.backgroundModelAssetPath&&longA->physics_fingerprint()==
      longB->physics_fingerprint()&&a->physics_fingerprint()!=longA->physics_fingerprint()&&
      fromRoot.memoryModel.particlesPerCell==0&&
      longFromRoot.memoryModel.particlesPerCell==0&&!fromRoot.source.enabled&&
      !longFromRoot.source.enabled&&fromRoot.shock==R::ShockAuthority::None&&
      longFromRoot.shock==R::ShockAuthority::None&&
      fromRoot.coordinateFrame=="HCI"&&longFromRoot.coordinateFrame=="HCI"&&
      longFromRoot.requestedTimeStepS==600&&longFromRoot.maximumTimeSteps==340;
  return Finish(good,status.ok()?"translated smoke and separate 1-AU decks resolve identical checksummed events from AMPS-root and srcSEP3D working directories":status.message);
}

Result PositiveProductionExample() {
  std::filesystem::path deck;
  for(const auto& base:{std::filesystem::path("."),std::filesystem::path("srcSEP3D"),
      std::filesystem::path("../srcSEP3D")}) {
    const auto candidate=base/"examples/shock-front/positive_1au.in";
    if(std::filesystem::exists(candidate)) {deck=std::filesystem::canonical(candidate);break;}
  }
  if(deck.empty())return Finish(false,"cannot locate positive 1-AU example");
  R::RunConfiguration3DOptions options;
  auto status=R::LoadConfigurationFile(deck.string(),&options);
  std::shared_ptr<const R::RunConfiguration3D> configuration;
  if(status.ok())status=R::RunConfiguration3D::Create(options,&configuration);
  std::shared_ptr<B::BackgroundProvider> background;
  if(status.ok())status=R::CreateBackgroundProvider(*configuration,&background);
  const auto adapter=std::dynamic_pointer_cast<A::ShockFrontBackgroundAdapter>(background);
  const auto shared=adapter?adapter->SharedProvider():nullptr;
  if(!status.ok()||!shared)return Finish(false,status.ok()?
      "positive example did not select the shared reduced adapter":status.message);
  const auto handoff=shared->HandoffTimeS();
  const auto endpoint=shared->EndpointTimeS();
  if(!handoff.ok()||!endpoint.ok())return Finish(false,
      "positive example has no exact handoff/endpoint roots");
  status=background->Prepare(endpoint.value);
  const auto epoch=adapter->FrontEpoch();
  const auto observer=shared->EvaluateFrontPoint(
      shared->Event().observerPositionM,endpoint.value,UINT64_C(1));
  const std::string surface=epoch?SEP3D::Output::SerializeReducedFrontTecplot(
      *epoch,shared->Event()):std::string{};
  const bool good=status.ok()&&configuration->options().intent==
      R::RunIntent::ShockPropagation&&configuration->options().shock==
      R::ShockAuthority::None&&!configuration->options().source.enabled&&
      configuration->options().memoryModel.particlesPerCell==0&&
      configuration->options().activeRegion==R::ActiveRegionMode::FullDomain&&
      configuration->options().requestedTimeStepS==600&&
      configuration->options().maximumTimeSteps==208&&
      std::fabs(handoff.value-10200.0)<1e-8&&
      std::fabs(endpoint.value-124200.0)<1e-6&&epoch&&
      epoch->trajectory.phase==SF::Phase::SwcmeOuter&&observer.ok()&&
      epoch->vertices.size()==static_cast<std::size_t>(
          shared->Event().polarCells*shared->Event().azimuthCells+1)&&
      epoch->triangles.size()==static_cast<std::size_t>(
          (2*shared->Event().polarCells-1)*shared->Event().azimuthCells)&&
      epoch->records.size()==epoch->triangles.size()&&
      observer.value.status==SF::FrontStatus::SolvedFastShock&&
      observer.value.fastMach>1&&observer.value.downstreamValid&&
      surface.find("ZONETYPE=FETRIANGLE")!=std::string::npos&&
      surface.find("DATAPACKING=BLOCK")!=std::string::npos&&
      surface.find("VARLOCATION=([4-39]=CELLCENTERED)")!=std::string::npos&&
      surface.find("\"triangle_stable_id\"")!=std::string::npos&&
      surface.find("\"planar_chord_area_m2\"")!=std::string::npos&&
      surface.find("\"theta_Bn_rad\"")!=std::string::npos&&
      surface.find("\"theta_Bn_valid\"")!=std::string::npos&&
      surface.find("\"magnetic_compression_valid\"")!=std::string::npos&&
      surface.find("\"rho2_kg_m3\"")!=std::string::npos&&
      surface.find("no_shock_fill=\"upstream-ambient-visualization-placeholder\"")!=
          std::string::npos&&surface.find("nan")==std::string::npos&&
      surface.find("volume_role=\"ambient-reference-only\"")!=std::string::npos;
  return Finish(good,good?
      "positive full-domain zero-particle deck reaches an independently evaluated accepted 1-AU shock and exports validity-aware RH surface limits":
      "positive example trajectory, acceptance, or surface contract is inconsistent");
}

} // namespace

std::vector<SEP3D::Testing::Descriptor> RegisterReducedShockFrontTests() {
  using D=SEP3D::Testing::Descriptor;
  auto make=[](const char* id,const char* name,SEP3D::Testing::TestCallback callback) {
    D d;d.id=id;d.name=name;d.group="RSHAPP";
    d.description="portable srcSEP3D boundary prerequisite for the reduced shock-front provider";
    d.initialization=SEP3D::Testing::InitializationLevel::None;
    d.supportedBuildModes="standalone-no-AMPS";
    d.runtime=SEP3D::Testing::RuntimeClass::Routine;d.seedPolicy="deterministic";
    d.stateIsolation="fresh immutable configuration/provider per callback";
    d.callback=std::move(callback);return d;
  };
  return {make("RSHAPP01","Reduced factory and ambient publication",FactoryAndAmbient),
      make("RSHAPP02","Reduced epoch handoff and rollback",EpochAndRollback),
      make("RSHAPP03","Reduced deck path resolution",DeckPathResolution),
      make("RSHAPP04","Positive production-style 1-AU example",PositiveProductionExample)};
}
