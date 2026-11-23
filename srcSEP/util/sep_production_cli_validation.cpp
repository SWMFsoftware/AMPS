#include "sep_production_cli_validation.h"
#include "sep_cli.h"
#include "../sep.h"
#include "../turbulence_wave_number_resolved/turbulence_wave_number_resolved.h"

#include <algorithm>
#include <cstdio>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace SEP { namespace Testing { namespace {
namespace C=Util::CLI;
namespace W=AlfvenTurbulence_Kolmogorov;
namespace K=Transport::Coefficient;

struct Checks {
  unsigned count=0,failures=0;
  std::string message;
  void Require(bool condition,const std::string& text) {
    ++count;if(condition)return;
    ++failures;if(message.size()<2000){if(!message.empty())message+="; ";message+=text;}
  }
  Result Finish(const char* success) const {
    Result result;result.status=failures?Status::Fail:Status::Pass;
    result.message=failures?message:success;
    result.configuration.push_back("evidence_level=native-linked-production-cli-software-contract");
    result.configuration.push_back("particle_evolution=not-executed;physics-qualification=not-claimed");
    result.metrics.push_back({"assertion_failures",static_cast<double>(failures),0.0,"==","count"});
    result.metrics.push_back({"independent_assertions",static_cast<double>(count),0.0,">","count"});
    return result;
  }
};

bool Parse(std::vector<std::string> arguments,C::Options& options,std::string* error=nullptr) {
  arguments.insert(arguments.begin(),"srcsep-production-cli-fixture");
  std::vector<char*> argv;
  for(auto& argument:arguments)argv.push_back(&argument[0]);
  std::ostringstream out,err;
  const bool valid=C::ParseCommandLine(static_cast<int>(argv.size()),argv.data(),options,out,err);
  if(error)*error=err.str();
  return valid;
}
C::Options Parsed(const std::vector<std::string>& arguments) {
  C::Options options;std::string error;
  if(!Parse(arguments,options,&error))throw std::runtime_error("valid PCLI fixture was rejected: "+error);
  return options;
}

// Test callbacks run after worker joins. Save all state actually mutated by
// ApplyTurbulenceOptions, including the legacy compatibility callback/scalar.
// No new mover implementation pointer is exposed for tests to bypass checks.
class State {
 public:
  bool coupling=W::ParticleCouplingMode,cascade=W::Cascade::active,reflection=W::Reflection::active;
  int injection=FieldLine::InjectionParameters::nParticlesPerIteration;
  Mover::ProductionMover mover=Mover::CurrentProductionMover();
  decltype(W::WaveNumberResolved::TurbulenceModelMode) representation=W::WaveNumberResolved::TurbulenceModelMode;
  K::Configuration coefficients=K::ActiveConfiguration();
  Transport::NumericalTolerances tolerances=Transport::ActiveNumericalTolerances();
  Turbulence::Configuration turbulence=Turbulence::ActiveConfiguration();
  Diffusion::fGetPitchAngleDiffusionCoefficient callback=Diffusion::GetPitchAngleDiffusionCoefficient;
  double constant=Diffusion::ConstPitchAngleDiffusionValue;
  bool restored=false;
  State() {
    // Selecting a mover resets its diagnostic count by public contract. This
    // fixture cannot truthfully restore a nonzero count through that API, so
    // reject before any mutation instead of discarding existing evidence.
    if(Mover::CompletedDispatchCount()!=0)
      throw std::runtime_error("PCLI state fixture requires a quiescent zero mover-dispatch count");
  }
  bool Restore() {
    if(restored)return true;
    W::ParticleCouplingMode=coupling;W::Cascade::active=cascade;W::Reflection::active=reflection;
    FieldLine::InjectionParameters::nParticlesPerIteration=injection;
    Mover::SelectProductionMover(mover);W::WaveNumberResolved::TurbulenceModelMode=representation;
    Diffusion::GetPitchAngleDiffusionCoefficient=callback;Diffusion::ConstPitchAngleDiffusionValue=constant;
    const bool a=K::SetActiveConfiguration(coefficients).ok();
    const bool b=Transport::SetActiveNumericalTolerances(tolerances).ok();
    const bool c=Turbulence::SetActiveConfiguration(turbulence).ok();
    restored=a&&b&&c;return restored;
  }
  ~State(){if(!restored)Restore();}
};

std::string ActualState() {
  std::ostringstream out;
  out.precision(17);
  out<<W::ParticleCouplingMode<<':'<<W::Cascade::active<<':'<<W::Reflection::active<<':'
      <<FieldLine::InjectionParameters::nParticlesPerIteration<<':'
      <<Mover::Describe(Mover::CurrentProductionMover()).canonicalName<<':'
      <<static_cast<int>(W::WaveNumberResolved::TurbulenceModelMode)<<':'
      <<K::ConfigurationFingerprint(K::ActiveConfiguration())<<':'
      <<Diffusion::ConstPitchAngleDiffusionValue<<':'
      <<Transport::ActiveNumericalTolerances().geometryFraction<<':'
      <<Turbulence::ActiveConfiguration().diagnosticCadence;
  return out.str();
}

void Applied(const C::Options& options,Checks& checks) {
  C::ApplyTurbulenceOptions(options);
  checks.Require(W::ParticleCouplingMode==options.particleCouplingMode,"coupling did not reach production global");
  checks.Require(W::Cascade::active==options.cascadeActive,"cascade did not reach production global");
  checks.Require(W::Reflection::active==options.reflectionActive,"reflection did not reach production global");
  checks.Require(FieldLine::InjectionParameters::nParticlesPerIteration==options.injectionParticlesPerIteration,
      "particle count did not reach native injection global");
  checks.Require(Mover::CurrentProductionMover()==options.particleMover,"wrong native mover selected");
  checks.Require(W::WaveNumberResolved::TurbulenceModelMode==
      (options.turbulenceModel==C::Options::TurbulenceModel::Integrated?
       W::WaveNumberResolved::ModelMode::Integrated:W::WaveNumberResolved::ModelMode::WaveNumberResolved),
      "wrong native turbulence representation selected");
  checks.Require(K::ConfigurationFingerprint(K::ActiveConfiguration())==K::ConfigurationFingerprint(options.coefficients),
      "actual coefficient configuration differs from parsed controls");
  checks.Require(Turbulence::ActiveConfiguration().diagnosticCadence==static_cast<std::size_t>(options.spectralOutputInterval),
      "diagnostic cadence did not reach the authoritative turbulence configuration");
}

void Artifact(Result& result,const std::string& leaf,const std::string& text) {
  const auto& directory=GetExecutionContext().artifactDirectory;
  if(directory.empty())throw std::runtime_error("PCLI native evidence directory is missing");
  const auto path=directory+"/"+leaf,temporary=path+".tmp";
  std::ofstream output(temporary.c_str(),std::ios::binary);output<<text;output.flush();output.close();
  if(output.fail()||std::rename(temporary.c_str(),path.c_str())!=0)
    throw std::runtime_error("cannot commit native production CLI evidence: "+path);
  result.artifacts.push_back(path);
}

Result Defaults() {
  Checks checks;const auto options=Parsed({});
  checks.Require(options.particleCouplingMode&&options.cascadeActive&&options.reflectionActive,"original Boolean defaults changed");
  checks.Require(options.particleMover==Mover::ProductionMover::FocusedTransportDiffusion&&
      options.turbulenceModel==C::Options::TurbulenceModel::Integrated,"original FTE/integrated default mapping changed");
  checks.Require(options.injectionParticlesPerIteration==300&&options.spectralOutputInterval==100&&!options.runTestManager,
      "original injection/cadence/TestManager defaults changed");
  checks.Require(!C::IsComponentTestExecutionRequested(options)&&!options.printHelp&&!options.listMovers&&!options.listTests,
      "no-argument production selection was redirected into development mode");
  checks.Require(options.inputPath.empty(),"empty private fixture unexpectedly inherited amps.in");
  Applied(options,checks);
  {std::ofstream input("amps.in");input<<"deliberately invalid: parser must only select this filename\n";}
  const auto fileDefault=Parsed({});std::remove("amps.in");
  checks.Require(fileDefault.inputPath=="amps.in"&&fileDefault.inputPathFromDefault,
      "documented existing ./amps.in production default changed");
  return checks.Finish("production defaults and their actual native application are preserved");
}

Result Help() {
  Checks checks;const auto before=ActualState();
  for(const auto& option:{"-h","--help"}) {
    C::Options options;
    checks.Require(Parse({option,"--unknown-after-help"},options)&&options.printHelp,
        "help did not exit parser before later invalid tokens");
  }
  std::ostringstream help;C::PrintHelp("amps",help);
  for(const auto& required:{"--coupling","--cascade","--cascase","--reflection","--particle-mover",
      "--turbulence-model","--particles-per-iteration","--spectrum-output-interval","--test-manager",
      "--testmanager","--run-test-manager","--no-test-manager","--list-tests","--tests",
      "--test-seed","--test-repeat","--test-output","--test-format","--fail-fast","--keep-test-artifacts"})
    checks.Require(help.str().find(required)!=std::string::npos,std::string("help omits ")+required);
  checks.Require(help.str().find("particles-per-iteration=300")!=std::string::npos&&
      help.str().find("spectrum-output-interval=100")!=std::string::npos,"help omits production defaults");
  checks.Require(ActualState()==before,"help changed actual production state");
  auto result=checks.Finish("complete CLI help and parser early-return leave native production state unchanged");
  Artifact(result,"PCLI02_help.txt",help.str());
  result.configuration.push_back("executable_pre_MPI_exit=separate-fresh-process-runner-gate");return result;
}

Result Booleans() {
  Checks checks;
  const std::vector<std::pair<std::string,bool>> words={{"on",true},{"TRUE",true},{"Yes",true},{"1",true},
      {"ENABLE",true},{"enabled",true},{"off",false},{"FALSE",false},{"No",false},{"0",false},{"DISABLE",false},{"disabled",false}};
  for(const auto& spelling:words)for(const auto& flag:{"--coupling","--coupling-mode","--particle-coupling",
      "--cascade","--cascase","--reflection","--test-manager","--testmanager"}) {
    const auto options=Parsed({flag,spelling.first});Applied(options,checks);
    const std::string key=flag;
    const bool actual=key.find("coupling")!=std::string::npos?W::ParticleCouplingMode:
        (key=="--cascade"||key=="--cascase")?W::Cascade::active:
        key=="--reflection"?W::Reflection::active:options.runTestManager;
    checks.Require(actual==spelling.second,"Boolean spelling/alias selected the wrong actual state");
  }
  for(const auto& flag:{"--no-coupling","--no-particle-coupling","--no-cascade","--no-cascase","--no-reflection"}) {
    const auto options=Parsed({flag});Applied(options,checks);
    const std::string key=flag;
    checks.Require(key.find("coupling")!=std::string::npos?!W::ParticleCouplingMode:
        key.find("casc")!=std::string::npos?!W::Cascade::active:!W::Reflection::active,
        "negative shortcut did not reach actual global");
  }
  return checks.Finish("all frozen Boolean vocabularies, aliases and negative shortcuts reach native state");
}

Result Models() {
  Checks checks;
  const std::vector<std::pair<std::string,C::Options::TurbulenceModel>> representations={
      {"integrated",C::Options::TurbulenceModel::Integrated},{"LEGACY",C::Options::TurbulenceModel::Integrated},
      {"old",C::Options::TurbulenceModel::Integrated},{"branch-integrated",C::Options::TurbulenceModel::Integrated},
      {"wave-number-resolved",C::Options::TurbulenceModel::WaveNumberResolved},{"wavenumber-resolved",C::Options::TurbulenceModel::WaveNumberResolved},
      {"k-resolved",C::Options::TurbulenceModel::WaveNumberResolved},{"K",C::Options::TurbulenceModel::WaveNumberResolved},
      {"spectral",C::Options::TurbulenceModel::WaveNumberResolved},{"new",C::Options::TurbulenceModel::WaveNumberResolved}};
  for(const auto& item:representations)for(bool equal:{false,true}) {
    const auto options=Parsed(equal?std::vector<std::string>{"--turbulence-model="+item.first}:
        std::vector<std::string>{"--turbulence-model",item.first});
    checks.Require(options.turbulenceModel==item.second,"turbulence alias selected the wrong enum");Applied(options,checks);
  }
  const std::vector<std::pair<std::string,Mover::ProductionMover>> movers={{"PARKER",Mover::ProductionMover::Parker},
      {"FTE-DMUMU",Mover::ProductionMover::FocusedTransportDiffusion},{"fte-mfp",Mover::ProductionMover::FocusedTransportMeanFreePath}};
  for(const auto& item:movers)for(bool equal:{false,true}) {
    const auto options=Parsed(equal?std::vector<std::string>{"--particle-mover="+item.first}:
        std::vector<std::string>{"--particle-mover",item.first});
    checks.Require(options.particleMover==item.second,"canonical mover selected wrong enum");Applied(options,checks);
  }
  for(const auto& retired:{"FTE","Parker3D","boris","drift"}) {
    C::Options options;checks.Require(!Parse({"--particle-mover",retired},options),"retired mover name accepted");
  }
  auto result=checks.Finish("model value syntax/case/aliases and the three canonical native mover selections are exact");
  result.configuration.push_back("historical_mover_aliases=retired-in-Step14;see-MIGRATION_MANIFEST");return result;
}

Result Numbers() {
  Checks checks;
  for(const auto& alias:{"--particles-per-iteration","--n-particles-per-iteration","--injection-particles-per-iteration","--field-line-injection-particles"}) {
    const auto options=Parsed({std::string(alias)+"=17"});Applied(options,checks);
    checks.Require(FieldLine::InjectionParameters::nParticlesPerIteration==17,"injection alias did not apply17");
  }
  for(const auto& alias:{"--spectrum-output-interval","--wave-number-output-interval","--turbulence-spectrum-output-interval","--spectral-output-interval"}) {
    const auto options=Parsed({alias,"0"});Applied(options,checks);
    checks.Require(Turbulence::ActiveConfiguration().diagnosticCadence==0,"cadence zero did not disable authoritative schedule");
  }
  const auto maximum=std::to_string(std::numeric_limits<int>::max());
  for(const auto& value:{maximum,std::string("+17"),std::string(" 17")}) {
    C::Options options;checks.Require(Parse({"--particles-per-iteration",value},options),"frozen signed/leading-whitespace positive grammar changed");
  }
  for(const auto& value:{"0","-0","-1","17x","17 ","1e2","2147483648","999999999999999999999"}) {
    C::Options options;checks.Require(!Parse({"--particles-per-iteration",value},options),"invalid positive injection integer accepted");
  }
  for(const auto& value:{"-1","3.0","9x","2147483648","999999999999999999999"}) {
    C::Options options;checks.Require(!Parse({"--spectrum-output-interval",value},options),"invalid nonnegative cadence integer accepted");
  }
  return checks.Finish("numeric aliases, zero/sign/whitespace policy and integer overflow preserve the current grammar");
}

Result Invalid() {
  Checks checks;const auto before=ActualState();
  const std::vector<std::vector<std::string>> invalid={{"--unregistered"},{"--coupling"},{"--coupling","perhaps"},
      {"--particle-mover"},{"--turbulence-model","unknown"},{"--no-cascade=on"},{"--no-coupling=off"},
      {"--run-test-manager=off"},{"positional-input"},{"--reflection","off","--unregistered"},
      {"--test-repeat","0"},{"--test-seed","-1"},{"--test-format","csv"}};
  for(const auto& arguments:invalid) {
    C::Options options;std::string error;
    checks.Require(!Parse(arguments,options,&error)&&error.find("ERROR:")!=std::string::npos,
        "malformed command was accepted or lacked contextual error");
    checks.Require(ActualState()==before,"failed partial parse mutated actual production globals");
  }
  return checks.Finish("unknown, missing, positional and flag-value failures leave applied production state untouched");
}

Result Application() {
  Checks checks;
  const auto options=Parsed({"--no-coupling","--no-cascade","--no-reflection","--particle-mover","fte-dmumu",
      "--turbulence-model","spectral","--particles-per-iteration","23","--spectrum-output-interval","0",
      "--pitch-angle-diffusion-provider","constant","--constant-dmumu","0.25","--transport-geometry-fraction","0.125"});
  Applied(options,checks);
  checks.Require(Diffusion::ConstPitchAngleDiffusionValue==0.25&&
      Diffusion::GetPitchAngleDiffusionCoefficient==Diffusion::Constant::GetPitchAngleDiffusionCoefficient,
      "constant SI coefficient or legacy compatibility callback not applied");
  checks.Require(Transport::ActiveNumericalTolerances().geometryFraction==0.125,"actual numerical controls not applied");
  auto result=checks.Finish("actual production globals, registry mover, SI provider and cadence match parsed settings");
  result.metrics.push_back({"applied_injection_count",static_cast<double>(FieldLine::InjectionParameters::nParticlesPerIteration),23.0,"==","count"});
  result.metrics.push_back({"applied_constant_Dmumu",Diffusion::ConstPitchAngleDiffusionValue,0.25,"==","1/s"});
  result.configuration.push_back("actual_state="+ActualState());return result;
}

Result Precedence() {
  Checks checks;C::Options options;
  Turbulence::Configuration file;
  file.source=Turbulence::Source::SelfConsistentSpectral;file.representation=Turbulence::Representation::Spectral;
  file.coupling=Turbulence::CouplingPolicy::Disabled;file.cascadeEnabled=false;file.reflectionEnabled=false;
  file.reflectionCoefficient=0.375;
  C::ApplyInputFileTurbulence(file,0.03,0.7,7,options);
  checks.Require(Parse({"--cascade","on","--cascade","off","--reflection=on","--turbulence-model","integrated",
      "--spectrum-output-interval=11"},options),"CLI reparse on actual input layer failed");
  Applied(options,checks);
  checks.Require(!W::Cascade::active&&W::Reflection::active&&!W::ParticleCouplingMode,
      "CLI last occurrence or input inheritance wrong in native globals");
  checks.Require(Turbulence::ActiveConfiguration().reflectionCoefficient==0.375&&
      Turbulence::ActiveConfiguration().diagnosticCadence==11&&options.turbulenceFromInputFile,
      "unchanged input value or explicit CLI cadence precedence lost");
  return checks.Finish("actual runtime-input layering, CLI reparse and last-occurrence precedence reach native state");
}

Result Reporting() {
  Checks checks;
  auto options=Parsed({"--particle-mover","parker","--turbulence-model","spectral","--coupling","on",
      "--cascade","off","--reflection","off","--particles-per-iteration","31","--spectrum-output-interval","0"});
  Applied(options,checks);std::ostringstream report;C::PrintTurbulenceOptions(options,report);
  for(const auto& expected:{"coupling: on","cascade:            off","reflection:               off","wave-number-resolved",
      "particle mover:                parker","particles/iteration:  31","0 iteration(s) (disabled)"})
    checks.Require(report.str().find(expected)!=std::string::npos,std::string("effective configuration omitted ")+expected);
  // The preserved three-mover registry currently advertises streaming for
  // every entry. That makes the historical missing-capability warning branch
  // unreachable through a valid CLI selection. Check the declared software
  // contract explicitly; this is not physical four-spectrum qualification.
  for(const auto& mover:{"parker","fte-dmumu","fte-mfp"}) {
    const auto supported=Parsed({"--particle-mover",mover,"--turbulence-model","spectral","--coupling","on"});
    checks.Require(Mover::Describe(supported.particleMover).capabilities.accumulatesWaveStreaming,
        "preserved streaming capability table changed and requires warning-gate review");
    std::ostringstream advertised;C::PrintTurbulenceOptions(supported,advertised);
    checks.Require(advertised.str().find("WARNING:")==std::string::npos,
        "report emitted a missing-capability warning for an advertised streaming mover");
  }
  options=Parsed({"--no-coupling","--turbulence-model","spectral","--particle-mover","parker"});
  Applied(options,checks);std::ostringstream quiet;C::PrintTurbulenceOptions(options,quiet);
  checks.Require(quiet.str().find("WARNING:")==std::string::npos,"disabled coupling emitted irrelevant streaming-source warning");
  auto result=checks.Finish("configuration output reports actual controls, disabled cadence and warnings consistent with advertised capabilities");
  result.configuration.push_back("legacy_missing_streaming_warning=unreachable-with-current-three-mover-table;branch-unexercised");
  result.configuration.push_back("advertised_streaming_capability=not-proof-of-physically-qualified-four-spectrum-streaming");
  Artifact(result,"PCLI09_effective_configuration.txt",report.str());return result;
}

Result Legacy() {
  Checks checks;
  checks.Require(!Parsed({}).runTestManager,"legacy TestManager default changed");
  for(const auto& arguments:std::vector<std::vector<std::string>>{{"--test-manager","on"},{"--testmanager=TRUE"},{"--run-test-manager"}}) {
    const auto options=Parsed(arguments);
    checks.Require(options.runTestManager&&!C::IsComponentTestExecutionRequested(options),"legacy diagnostic enable was redirected to selectable test-only mode");
  }
  for(const auto& arguments:std::vector<std::vector<std::string>>{{"--test-manager","off"},{"--testmanager=0"},
      {"--run-test-manager","--no-test-manager"}})
    checks.Require(!Parsed(arguments).runTestManager,"legacy diagnostic disable/precedence changed");
  const auto before=ActualState();
  for(const auto& selector:std::vector<std::vector<std::string>>{{"--test","PCLI01"},
      {"--test-group","production-cli"},{"--tests","PCLI01"},{"--all-tests"}}) {
    auto arguments=std::vector<std::string>{"--run-test-manager"};
    arguments.insert(arguments.end(),selector.begin(),selector.end());
    C::Options mixed;std::string error;
    checks.Require(!Parse(arguments,mixed,&error)&&error.find("TestManager")!=std::string::npos,
        "enabled legacy TestManager was ignored by component-only dispatch");
    checks.Require(ActualState()==before,"mixed-mode parse error mutated actual production globals");
  }
  const auto disabled=Parsed({"--run-test-manager","--no-test-manager","--test","PCLI01"});
  checks.Require(!disabled.runTestManager&&C::IsComponentTestExecutionRequested(disabled),
      "explicit final legacy disable must preserve component mode");
  auto result=checks.Finish("legacy diagnostic switches retain default-off, aliases and standalone production-dispatch semantics");
  result.configuration.push_back("heavy_TestManager_execution=not-invoked;separate-legacy-diagnostic-gate");return result;
}

Result Matrix() {
  Checks checks;unsigned cases=0;
  for(const auto& mover:{"parker","fte-dmumu","fte-mfp"})for(unsigned bits=0;bits<64;++bits) {
    const auto options=Parsed({"--particle-mover",mover,"--turbulence-model",bits&1?"spectral":"integrated",
        "--coupling",bits&2?"on":"off","--cascade",bits&4?"on":"off","--reflection",bits&8?"on":"off",
        "--particles-per-iteration",bits&16?"17":"41","--spectrum-output-interval",bits&32?"0":"7"});
    Applied(options,checks);++cases;
  }
  auto result=checks.Finish("192 mover/representation/operator/count/cadence combinations apply coherently to actual native model state");
  result.metrics.push_back({"applied_configuration_cases",static_cast<double>(cases),192.0,"==","count"});
  result.configuration.push_back("startup_after_allocation=separate-fresh-process-runner-gate");
  result.configuration.push_back("production_particle_steps=unqualified-polarized-closure;not-executed");return result;
}

Result Isolation() {
  Checks checks;const auto before=ActualState();
  // Persistence controls remain operational parser state. Both value forms
  // and last-occurrence precedence must work without changing any model or
  // source global before a complete material-run preflight.
  const auto checkpoint=Parsed({"--input","declared.in","--reduced-shock-event","declared.event",
      "--checkpoint-save","first.bin","--checkpoint-save=final.bin","--checkpoint-at-step=4",
      "--resume-checkpoint","previous.bin"});
  checks.Require(checkpoint.checkpointSavePath=="final.bin"&&checkpoint.checkpointAtStep==4&&
      checkpoint.checkpointResumePath=="previous.bin","checkpoint parser lost value or precedence");
  for(const auto& arguments:std::vector<std::vector<std::string>>{
      {"--checkpoint-save","missing-step.bin"},{"--checkpoint-at-step=4"},
      {"--input","declared.in","--reduced-shock-event","declared.event","--checkpoint-at-step=0","--checkpoint-save=x"},
      {"--resume-checkpoint=x"},
      {"--input","declared.in","--reduced-shock-event","declared.event","--resume-checkpoint=x","--test=PCLI01"}}) {
    C::Options rejected;checks.Require(!Parse(arguments,rejected),"invalid checkpoint controls entered a run");
  }
  for(unsigned invocation=0;invocation<32;++invocation) {
    const auto options=Parsed({"--coupling","off","--particles-per-iteration","19","--test-manager","on"});
    checks.Require(options.injectionParticlesPerIteration==19&&options.runTestManager,"private parser state changed under repetition");
    const auto fresh=Parsed({});
    checks.Require(fresh.injectionParticlesPerIteration==300&&fresh.particleCouplingMode&&!fresh.runTestManager,
        "one Options instance leaked into a fresh parser invocation");
    checks.Require(ActualState()==before,"parser repetition changed applied production state");
  }
  auto result=checks.Finish("repeated private parser invocations leave native globals and fresh production defaults unchanged");
  result.configuration.push_back("coupled_CLI_ownership=separate-source-callsite-audit;no-SWMF-execution-claim");return result;
}

Descriptor Make(const char* id,const char* name,TestCallback callback) {
  Descriptor descriptor;descriptor.id=id;descriptor.name=name;descriptor.group="production-cli";
  descriptor.description=name;descriptor.initialization=InitializationLevel::None;
  descriptor.supportedBuildModes="native linked srcSEP, workers quiescent, zero prior mover dispatches";
  descriptor.seedPolicy="deterministic software inputs; no random draws";
  descriptor.seedCapability=SeedCapability::Deterministic;descriptor.repeatSupported=true;
  descriptor.stateIsolation="complete scope for every ApplyTurbulenceOptions global; fresh process supplies startup/coupled-ownership evidence";
  descriptor.callback=[callback]() {
    State scope;Result result;
    try{result=callback();}
    catch(const std::exception& error){result.status=Status::Error;result.message=error.what();}
    if(!scope.Restore()){result.status=Status::Error;result.message="cannot restore original production CLI configuration";}
    return result;
  };
  return descriptor;
}
} // namespace

std::vector<Descriptor> ProductionCliDescriptors() {
  return {Make("PCLI01","Production defaults and native application",Defaults),
      Make("PCLI02","Complete help and side-effect-free parser exit",Help),
      Make("PCLI03","Boolean vocabulary and compatibility aliases",Booleans),
      Make("PCLI04","Turbulence and canonical mover value syntax",Models),
      Make("PCLI05","Numeric control boundaries and aliases",Numbers),
      Make("PCLI06","Unknown and malformed command rejection",Invalid),
      Make("PCLI07","Actual production model control application",Application),
      Make("PCLI08","Input-layer and CLI ordering precedence",Precedence),
      Make("PCLI09","Effective configuration and coupling warning",Reporting),
      Make("PCLI10","Legacy TestManager control mapping",Legacy),
      Make("PCLI11","Native applied production option matrix",Matrix),
      Make("PCLI12","Parser repetition and state isolation",Isolation)};
}
} }
