#include "sep_birth_window_validation.h"
#include "sep_birth_window.h"
#include "amps_fixed_step.h"
#include <cmath>
#include <limits>

namespace SEP { namespace Testing { namespace {
namespace B=Transport::Birth;
struct Checks {
  unsigned failures=0;
  void Require(bool condition){if(!condition)++failures;}
  Result Finish(const char* success) const {
    Result result;result.status=failures?Status::Fail:Status::Pass;
    result.message=failures?"birth-window contract assertion failed":success;
    result.configuration={"evidence_level=dependency-light-time-lifecycle-core",
        "native_source_adapter_integration=not-qualified-by-this-test",
        "time_units=s;birth_phase=open-(0,1);no-clipping-or-redraw"};
    result.metrics.push_back({"assertion_failures",static_cast<double>(failures),0,"<=","count"});
    return result;
  }
};
bool Same(const B::Record& a,const B::Record& b) {
  return a.schema==b.schema&&a.origin==b.origin&&a.pendingFirstMove==b.pendingFirstMove&&
      a.parent.schema==b.parent.schema&&a.parent.startEpochS==b.parent.startEpochS&&
      a.parent.durationS==b.parent.durationS&&a.parent.endEpochS==b.parent.endEpochS&&
      a.firstMoveDurationS==b.firstMoveDurationS;
}
Result DriverWindow() {
  Checks checks;B::DriverWindow window;
  checks.Require(B::FreezeWindow(100.0,0.5,&window).ok());
  checks.Require(window.schema==1&&window.startEpochS==100.0&&window.durationS==0.5&&window.endEpochS==100.5);
  const auto original=window;
  for(const double duration:{0.0,-1.0,std::numeric_limits<double>::infinity(),std::numeric_limits<double>::quiet_NaN()}) {
    checks.Require(!B::FreezeWindow(100.0,duration,&window).ok());
    checks.Require(window.startEpochS==original.startEpochS&&window.endEpochS==original.endEpochS&&window.durationS==original.durationS);
  }
  checks.Require(!B::FreezeWindow(1.0e12,1.0e-8,&window).ok());
  checks.Require(!B::ValidateWindow(B::DriverWindow()).ok());
  return checks.Finish("driver timing is explicit and cannot infer or round away a complete host step");
}
Result LateBirth() {
  Checks checks;B::DriverWindow window;B::Record record;B::FirstMovePlan plan;
  checks.Require(B::FreezeWindow(1.0e12,0.01,&window).ok());
  const double phase=std::nextafter(1.0,0.0);
  const double expected=(1.0-phase)*0.01;
  checks.Require(B::UniformBirth(window,phase,&record).ok());
  checks.Require(record.firstMoveDurationS==expected&&expected>0.0);
  // This is the failure mechanism of the old absolute-epoch design. The
  // rounded display epoch loses all of the physically prescribed remainder.
  checks.Require(window.endEpochS-record.firstMoveDurationS==window.endEpochS);
  checks.Require(B::PlanFirstMove(record,window,0.01,&plan).ok());
  const auto remaining=B::SnapshotRemaining(plan,0.0,window.startEpochS,window.endEpochS);
  checks.Require(remaining.status.ok()&&remaining.value==expected&&plan.durationS==expected);
  auto result=checks.Finish("late birth keeps a positive first movement even when its absolute epoch rounds to the step end");
  result.metrics.push_back({"preserved_first_move_duration",plan.durationS,0,">","s"});
  return result;
}
Result FirstMoveCommit() {
  Checks checks;B::DriverWindow window;B::Record record,consumed;B::FirstMovePlan plan;
  checks.Require(B::FreezeWindow(10.0,4.0,&window).ok());
  checks.Require(B::UniformBirth(window,0.25,&record).ok());
  checks.Require(B::PlanFirstMove(record,window,4.0,&plan).ok());
  checks.Require(plan.durationS==3.0&&record.pendingFirstMove==1);
  checks.Require(B::StageConsumedFirstMove(record,plan,&consumed).ok());
  checks.Require(record.pendingFirstMove==1&&consumed.pendingFirstMove==0&&
      consumed.firstMoveDurationS==3.0&&B::ValidateRecord(consumed).ok());
  checks.Require(!B::StageConsumedFirstMove(consumed,plan,&record).ok());
  checks.Require(B::StepStartBirth(window,&record).ok());
  checks.Require(B::PlanFirstMove(record,window,4.0,&plan).ok()&&plan.durationS==4.0);
  return checks.Finish("a source gets exactly its remaining duration and consumption is staged without changing its original history");
}
Result Rejection() {
  Checks checks;B::DriverWindow window,other;B::Record record,sentinel;
  B::FirstMovePlan plan;
  checks.Require(B::FreezeWindow(10.0,4.0,&window).ok());
  checks.Require(B::FreezeWindow(14.0,4.0,&other).ok());
  checks.Require(B::UniformBirth(window,0.5,&record).ok());
  checks.Require(B::InitialCondition(99.0,&sentinel).ok());
  const auto saved=sentinel;
  for(const double phase:{0.0,1.0,-0.1,std::numeric_limits<double>::quiet_NaN(),
      std::numeric_limits<double>::denorm_min()}) {
    checks.Require(!B::UniformBirth(window,phase,&sentinel).ok());checks.Require(Same(sentinel,saved));
  }
  plan.durationS=123.0;
  checks.Require(!B::PlanFirstMove(record,other,4.0,&plan).ok()&&plan.durationS==123.0);
  checks.Require(!B::PlanFirstMove(record,window,2.0,&plan).ok()&&plan.durationS==123.0);
  checks.Require(!B::ValidateRecord(B::Record()).ok());
  auto invalid=record;invalid.schema=2;
  checks.Require(!B::StageConsumedFirstMove(invalid,plan,&sentinel).ok()&&Same(sentinel,saved));
  return checks.Finish("unknown history, incompatible windows, invalid phases, and unequal host timesteps reject without output mutation");
}
Result SnapshotInterval() {
  Checks checks;B::DriverWindow window;B::Record record;B::FirstMovePlan plan;
  checks.Require(B::FreezeWindow(10.0,4.0,&window).ok());
  checks.Require(B::UniformBirth(window,0.25,&record).ok());
  checks.Require(B::PlanFirstMove(record,window,4.0,&plan).ok());
  const auto initial=B::SnapshotRemaining(plan,0.0,10.0,14.0);
  const auto middle=B::SnapshotRemaining(plan,1.5,10.0,14.0);
  const auto final=B::SnapshotRemaining(plan,3.0,10.0,14.0);
  checks.Require(initial.status.ok()&&initial.value==3.0&&middle.status.ok()&&
      middle.value==1.5&&final.status.ok()&&final.value==0.0);
  checks.Require(!B::SnapshotRemaining(plan,0.0,11.0,14.0).status.ok());
  checks.Require(!B::SnapshotRemaining(plan,0.0,10.0,std::nextafter(14.0,0.0)).status.ok());
  checks.Require(!B::SnapshotRemaining(plan,3.1,10.0,14.0).status.ok());
  const auto open=B::SnapshotRemaining(plan,0.0,10.0,std::numeric_limits<double>::infinity());
  checks.Require(open.status.ok()&&std::isinf(open.value));
  return checks.Finish("immutable snapshot coverage is checked against the complete parent window and the birth-relative interval");
}
Result InitialHistory() {
  Checks checks;B::Record record;B::DriverWindow window;B::FirstMovePlan plan;
  checks.Require(B::InitialCondition(1.0e12,&record).ok());
  checks.Require(B::ValidateRecord(record).ok()&&!record.pendingFirstMove&&
      record.origin==B::Origin::InitialCondition&&record.parent.startEpochS==1.0e12);
  checks.Require(B::FreezeWindow(1.0e12,0.01,&window).ok());
  checks.Require(!B::PlanFirstMove(record,window,0.01,&plan).ok());
  auto invalid=record;invalid.pendingFirstMove=1;
  checks.Require(!B::ValidateRecord(invalid).ok());
  invalid=record;invalid.origin=static_cast<B::Origin>(999);
  checks.Require(!B::ValidateRecord(invalid).ok());
  return checks.Finish("explicit initial conditions carry history without inventing a scheduled physical birth or previous age");
}
Result TerminalInterval() {
  Checks checks;B::DriverWindow window;B::Record record;B::FirstMovePlan plan;
  checks.Require(B::FreezeWindow(1.0e12,0.01,&window).ok());
  checks.Require(B::UniformBirth(window,std::nextafter(1.0,0.0),&record).ok());
  checks.Require(B::PlanFirstMove(record,window,0.01,&plan).ok());
  const double floorS=1.0e-12;
  const auto old=Transport::SelectSubstep(plan.durationS,{},floorS,nullptr);
  checks.Require(old.status.code==Transport::StatusCode::StepUnderflow);
  Transport::StepDiagnostics diagnostics;
  const auto whole=B::SelectFirstMoveSubstep(plan,0.0,{},floorS,&diagnostics);
  const auto bounded=B::SelectFirstMoveSubstep(plan,0.0,{{"snapshot",plan.durationS}},floorS,&diagnostics);
  const auto refined=B::SelectFirstMoveSubstep(plan,0.0,
      {{"diffusion",0.5*plan.durationS}},floorS,&diagnostics);
  checks.Require(whole.status.ok()&&whole.value==plan.durationS&&
      bounded.status.ok()&&bounded.value==plan.durationS&&
      refined.status.code==Transport::StatusCode::StepUnderflow&&
      diagnostics.selections==3&&diagnostics.underflowRejects==1);
  checks.Require(!B::SelectFirstMoveSubstep(plan,0.0,{},-1.0).status.ok());
  checks.Require(!B::SelectFirstMoveSubstep(plan,0.0,{},
      std::numeric_limits<double>::quiet_NaN()).status.ok());
  checks.Require(!B::SelectFirstMoveSubstep(plan,0.0,{},
      std::numeric_limits<double>::infinity()).status.ok());
  return checks.Finish("a prescribed positive birth remainder is completed without weakening adaptive refinement rejection");
}
Result HostClockCapability() {
  Checks checks;B::HostClockPolicy policy;
  for(const bool counter:{false,true}) {
    policy.counterAtEnd=counter;policy.mode=B::HostStepMode::SingleGlobal;policy.speciesCount=3;
    checks.Require(B::ValidateHostClockPolicy(policy).ok());
    policy.mode=B::HostStepMode::SpeciesGlobal;
    checks.Require(!B::ValidateHostClockPolicy(policy).ok());
    policy.speciesCount=1;
    checks.Require(B::ValidateHostClockPolicy(policy).ok()==counter);
  }
  policy.mode=B::HostStepMode::SingleGlobal;policy.externallyOwned=true;
  checks.Require(!B::ValidateHostClockPolicy(policy).ok());
  policy.externallyOwned=false;policy.speciesCount=0;
  checks.Require(!B::ValidateHostClockPolicy(policy).ok());
  policy.speciesCount=1;policy.mode=static_cast<B::HostStepMode>(999);
  checks.Require(!B::ValidateHostClockPolicy(policy).ok());
  return checks.Finish("audited PIC clock policy rejects unimplemented species synchronization and externally owned host-end assumptions before source mutation");
}
Result BoundExternalClock() {
  // Independent exact-second oracle for the observed failing 20 -> 21 s
  // interval. Binding that full interval must add support while leaving an
  // unbound external clock, wrong departure/dt and unaudited policies rejected.
  Checks checks;B::DriverWindow window;B::HostClockPolicy policy;
  policy.externallyOwned=true;policy.counterAtEnd=false;
  checks.Require(B::FreezeWindow(20.0,1.0,&window).ok());
  checks.Require(!B::ValidateHostClockPolicy(policy).ok());
  checks.Require(B::ValidateExternalHostWindow(policy,window,20.0,1.0).ok());
  checks.Require(!B::ValidateExternalHostWindow(policy,window,21.0,1.0).ok());
  checks.Require(!B::ValidateExternalHostWindow(policy,window,20.0,0.5).ok());
  policy.counterAtEnd=true;
  checks.Require(!B::ValidateExternalHostWindow(policy,window,20.0,1.0).ok());
  policy.counterAtEnd=false;policy.mode=B::HostStepMode::SpeciesGlobal;
  checks.Require(!B::ValidateExternalHostWindow(policy,window,20.0,1.0).ok());
  return checks.Finish("external clock requires an explicit full nominal window and audited nonadvancing PIC policy");
}
Result CoupledWholeSteps() {
  // Configure a separate private fixture once with a binary-exact 3/4 s dt;
  // never modify a running AMPS timestep. The listed endpoints/counts follow
  // directly from whole multiples of 3/4, independently of Plan's arithmetic.
  Checks checks;AMPS::Coupling::FixedStepContract contract;
  contract.Configure(0.75);double completed=20.0,pic=20.0,sampled=0.0;
  const double limits[]={20.5,21.0,21.5,22.0,23.0};
  const unsigned expected[]={0,1,1,0,2};
  const double endpoints[]={20.0,20.75,21.5,21.5,23.0};
  for(unsigned i=0;i<5;++i) {
    unsigned steps=0;
    while(contract.Plan(completed,limits[i],0.75)) {
      checks.Require(contract.Current().startEpochS==completed&&contract.Current().durationS==0.75);
      pic=completed;contract.EnterParticlePhase(pic,0.75);
      // This controlled exposure counter represents one actual sampling call
      // per completed physical step; host acknowledgements never enter it.
      sampled+=0.75;contract.CompleteParticlePhase(pic,0.75);
      completed=contract.Commit(pic,0.75);++steps;
    }
    checks.Require(steps==expected[i]&&completed==endpoints[i]&&pic==completed);
  }
  checks.Require(sampled==3.0&&contract.NominalDurationS()==0.75&&!contract.Active());
  return checks.Finish("nondivisible host limits retain only complete fixed-duration steps and no exposure on skips");
}
Result CoupledClockRejection() {
  // Exercise failure atomicity and phase ownership with private values only.
  // Failed dt requests must not alter the nominal duration. The adjacent
  // representable limit below 21 s cannot admit [20,21]; the exact limit can.
  // Repeated entry/completion shortcuts must fail even if an endpoint matches.
  Checks checks;AMPS::Coupling::FixedStepContract contract;contract.Configure(1.0);
  auto rejects=[&](const std::function<void()>& operation) {
    bool rejected=false;try{operation();}catch(const std::exception&){rejected=true;}
    checks.Require(rejected);
  };
  rejects([&]{contract.Configure(0.5);});
  rejects([&]{contract.Plan(20.0,21.0,0.5);});
  rejects([&]{contract.Plan(20.0,std::numeric_limits<double>::quiet_NaN(),1.0);});
  checks.Require(!contract.Plan(20.0,std::nextafter(21.0,20.0),1.0));
  checks.Require(contract.Plan(20.0,21.0,1.0));
  double pic=20.0;
  rejects([&]{contract.Commit(pic,1.0);});
  rejects([&]{contract.Plan(20.0,21.0,1.0);});
  contract.EnterParticlePhase(pic,1.0);checks.Require(pic==21.0);
  rejects([&]{contract.EnterParticlePhase(pic,1.0);});
  pic=20.0;rejects([&]{contract.CompleteParticlePhase(pic,1.0);});
  pic=21.0;contract.CompleteParticlePhase(pic,1.0);
  checks.Require(contract.Commit(pic,1.0)==21.0&&!contract.Active());
  rejects([&]{contract.Plan(1.0e20,1.0e20,1.0);});
  checks.Require(contract.NominalDurationS()==1.0);
  return checks.Finish("clock phases, representability and constant timestep reject before invalid publication");
}
} // namespace
std::vector<Descriptor> BirthWindowDescriptors() {
  const char* names[]={"Explicit frozen source driver window","Late birth retains exact positive remaining duration",
      "Staged first movement consumption","Atomic birth history rejection","Birth-relative immutable snapshot interval",
      "Initial-condition and scheduled-source distinction","Prescribed terminal interval and adaptive refinement floor",
      "Actual host clock advancement capability","Explicit external full-step binding",
      "Coupled fixed-step scheduling and exposure","Coupled clock phase and timestep rejection"};
  const std::function<Result()> callbacks[]={DriverWindow,LateBirth,FirstMoveCommit,Rejection,SnapshotInterval,InitialHistory,TerminalInterval,HostClockCapability,BoundExternalClock,CoupledWholeSteps,CoupledClockRejection};
  std::vector<Descriptor> result;
  for(unsigned i=0;i<sizeof(callbacks)/sizeof(callbacks[0]);++i) {
    // Preserve BIRTH01--08 identifiers and append the coupled value contracts;
    // two-digit formatting must not rename case 10 to the unrelated BIRTH010.
    Descriptor descriptor;descriptor.id=std::string("BIRTH")+(i<9?"0":"")+std::to_string(i+1);descriptor.name=names[i];
    descriptor.description=names[i];descriptor.group="birth-window";
    descriptor.initialization=InitializationLevel::None;descriptor.runtime=RuntimeClass::Routine;
    descriptor.supportedBuildModes="dependency-light and native registry; time lifecycle core only";
    descriptor.seedCapability=SeedCapability::Deterministic;descriptor.repeatSupported=true;
    descriptor.seedPolicy="deterministic bit-extreme and exact interval cases; no random draws";
    descriptor.stateIsolation="private value records; no native globals or clock writes";
    descriptor.callback=callbacks[i];result.push_back(descriptor);
  }
  return result;
}
} }
