// Propagation acceptance tests use the production parser/provider and publication
// boundary. No PIC/MPI mock is accepted as proof of a native coupled campaign.
#include "sep3d_test_registry.h"
#include "configuration_io.h"
#include "source_runtime.h"
#include "shock_history.h"
#include <cmath>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <sstream>

namespace {
namespace R=SEP3D::RuntimeModel;
namespace A=SEP3D::Adapters;
namespace O=SEP3D::Output;
using Result=SEP3D::Testing::Result;
Result Finish(bool good,const std::string& message) {
  Result result; result.status=good?SEP3D::Testing::Status::Pass:SEP3D::Testing::Status::Fail;
  result.message=message; return result;
}
Result RunCME3D03() {
  R::RunConfiguration3DOptions options;
  auto status=R::LoadConfigurationFile("examples/sep3d_swcme_20rs_1au.in",&options);
  if (!status.ok()) return Finish(false,"propagation example parse failed: "+status.message);
  std::shared_ptr<const R::RunConfiguration3D> config;
  status=R::RunConfiguration3D::Create(options,&config);
  if (!status.ok()) return Finish(false,"propagation freeze failed: "+status.message);
  std::string summary;
  status=R::BuildDryRunSummary(*config,&summary);
  if (!status.ok() || summary.find("run_intent=shock-propagation")==std::string::npos ||
      summary.find("source_enabled=false")==std::string::npos)
    return Finish(false,"native preflight could not plan the submitted propagation geometry: "+status.message);
  if (options.intent!=R::RunIntent::ShockPropagation || options.source.enabled || !options.observers.empty() ||
      options.stopShockRadiusM<=SEP3D::Core::Const::AU || options.maximumTimeSteps!=7200)
    return Finish(false,"propagation example did not preserve its source-off/stop contract");
  // These conflicts must be caught before AMPS allocates a mesh. The optional
  // radius target cannot silently replace the complete finite step budget.
  auto invalid=[&](R::RunConfiguration3DOptions bad) {
    std::shared_ptr<const R::RunConfiguration3D> unused;
    return !R::RunConfiguration3D::Create(bad,&unused).ok();
  };
  auto bad=options; bad.source.enabled=true; if(!invalid(bad)) return Finish(false,"enabled source accepted");
  bad=options;bad.inputSchemaVersion=3;if(!invalid(bad))return Finish(false,"schema-3 propagation accepted");
  bad=options;bad.restartInputPath="old.restart";if(!invalid(bad))return Finish(false,"propagation restart accepted");
  bad=options;bad.stopShockRadiusM=2*SEP3D::Core::Const::AU;if(!invalid(bad))return Finish(false,"stop beyond outer boundary accepted");
  bad=options;bad.stopShockRadiusM=20*SEP3D::Core::Const::R_sun;if(!invalid(bad))return Finish(false,"stop at launch accepted");
  std::ifstream example("examples/sep3d_swcme_20rs_1au.in");
  const std::string deck((std::istreambuf_iterator<char>(example)),std::istreambuf_iterator<char>());
  auto rejects=[&](const std::string& from,const std::string& to) {
    std::string changed=deck;const auto position=changed.find(from);
    if(position==std::string::npos) return false;
    changed.replace(position,from.size(),to);
    R::RunConfiguration3DOptions unused;return !R::ParseConfigurationText(changed,&unused).ok();
  };
  if(!rejects("event.valid_until = 432000 s","event.valid_until = 86400 s") ||
     !rejects("event.launch_epoch = 0 s","event.launch_epoch = 3600 s"))
    return Finish(false,"propagation accepted incomplete launch-to-stop temporal coverage");
  // Check the installed canonical provider against the independently integrated
  // signed drag law r=r0+w*t+log(1+gamma*(v0-w)*t)/gamma for v0>w.
  // This calculation is a portable oracle, never exported as native evidence.
  std::shared_ptr<A::ShockProvider> provider;
  status=A::CreateStandaloneSwcmeShockProvider(*config,&provider);
  if(!status.ok())return Finish(false,"source-free provider failed: "+status.message);
  double previous=0; bool crossed=false;
  const double gamma=1e-11,wind=400000,excess=600000;
  for(std::uint64_t tick=0;tick<=options.maximumTimeSteps;++tick) {
    const double time=tick*options.requestedTimeStepS;
    const auto state=provider->Evaluate(time);
    const double expected=20*SEP3D::Core::Const::R_sun+wind*time+std::log1p(gamma*excess*time)/gamma;
    const double expectedSpeed=wind+excess/(1+gamma*excess*time);
    if(!state.status.ok() || !state.active || !state.patches.empty() || state.generation!=tick+1 ||
       std::fabs(state.radiusM-expected)>1e-9*expected ||
       std::fabs(state.radialSpeedMPerS-expectedSpeed)>1e-9*expectedSpeed || (tick && state.radiusM<=previous))
      return Finish(false,"canonical source-off front disagrees with the independent drag oracle or creates patches");
    previous=state.radiusM;
    if(state.radiusM>=options.stopShockRadiusM){crossed=true;break;}
  }
  return Finish(crossed,"schema-4 propagation freezes a finite source-off run; canonical DBM reaches the declared stop without injection patches");
}
Result RunCME3D04() {
  namespace fs=std::filesystem;
  const fs::path path=fs::temp_directory_path()/("sep3d-propagation-writer-"+
      std::to_string(std::chrono::steady_clock::now().time_since_epoch().count())+".csv");
  std::error_code error;fs::remove(path,error);
  O::ShockHistorySample sample;
  sample.radiusM=20*SEP3D::Core::Const::R_sun;sample.speedMPerS=700000;
  sample.active=true;sample.generation=1;sample.providerIdentity="test-oracle";sample.configurationFingerprint="frozen";
  O::ShockHistoryWriter writer;
  if(!writer.Open(path,60).ok()||!writer.Append(sample).ok())return Finish(false,"tick-zero publication failed");
  auto bad=sample;bad.tick=2;bad.timeS=120;bad.generation=3;bad.radiusM+=84000000;
  if(writer.Append(bad).ok())return Finish(false,"history gap accepted");
  sample.tick=1;sample.timeS=60;sample.generation=2;sample.radiusM+=42000000;
  bad=sample;bad.injections=1;if(writer.Append(bad).ok())return Finish(false,"nonzero cumulative injection accepted");
  bad=sample;bad.mpiRadiusSpreadM=2;if(writer.Append(bad).ok())return Finish(false,"MPI geometry disagreement accepted");
  bad=sample;bad.configurationFingerprint="changed";if(writer.Append(bad).ok())return Finish(false,"changed provider identity accepted");
  bad=sample;bad.timeS=61;if(writer.Append(bad).ok())return Finish(false,"non-tick clock accepted");
  if(!writer.Append(sample).ok()||!writer.Close().ok())return Finish(false,"completed history close failed");
  O::ShockHistoryWriter overwrite;if(overwrite.Open(path,60).ok())return Finish(false,"existing history overwritten");
  std::ifstream file(path);std::string line;unsigned lines=0;while(std::getline(file,line))++lines;
  fs::remove(path,error);
  return Finish(lines==3,"history closes tick-zero/next-tick rows and rejects clock gaps, injection, MPI disagreement, changed identity and overwrites");
}
}
std::vector<SEP3D::Testing::Descriptor> RegisterShockPropagationTests() {
  using D=SEP3D::Testing::Descriptor;
  auto make=[](const char* id,const char* name,SEP3D::Testing::TestCallback callback){
    D d;d.id=id;d.name=name;d.group="CME3D";d.description="Source-free installed SWCME propagation prerequisite";
    d.initialization=SEP3D::Testing::InitializationLevel::None;d.supportedBuildModes="standalone-no-AMPS";
    d.runtime=SEP3D::Testing::RuntimeClass::Routine;d.seedPolicy="deterministic, no random draws";
    d.stateIsolation="fresh provider and temporary telemetry per callback";d.callback=std::move(callback);return d;
  };
  return {make("CME3D03","Propagation configuration/canonical provider",RunCME3D03),
          make("CME3D04","Native telemetry publication contract",RunCME3D04)};
}
