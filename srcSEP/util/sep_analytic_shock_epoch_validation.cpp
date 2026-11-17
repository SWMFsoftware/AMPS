#include "sep_analytic_shock_epoch_validation.h"
#include "sep_analytic_shock_epoch.h"
#include <cmath>
#include <limits>

namespace SEP { namespace Testing { namespace {
namespace S=Shock;
Result Finish(bool pass,const std::string& message) {
  Result result;result.status=pass?Status::Pass:Status::Fail;result.message=message;
  result.metrics.push_back({"assertion_failures",pass?0.0:1.0,0,"<=","count"});
  result.configuration={"scope=immutable-explicit-analytic-trajectory-values;no-PIC-host-or-source-frame-qualification"};return result;
}
S::TrajectoryConfiguration Constant() {
  S::TrajectoryConfiguration result;result.launchEpochS=2;result.launchRadiusM=1e10;
  S::SpeedKnot a;a.radiusM=1e10;a.speedMPerS=1e6;
  S::SpeedKnot b;b.radiusM=3e10;b.speedMPerS=1e6;result.knots={a,b};return result;
}
Result ExplicitDeparture() {
  const auto configuration=Constant();S::SourceEpochView earlier,later,repeated;
  bool pass=S::EvaluateSourceAtEpoch(configuration,2.5,&earlier).ok()&&
      S::EvaluateSourceAtEpoch(configuration,100,&later).ok()&&
      S::EvaluateSourceAtEpoch(configuration,2.5,&repeated).ok();
  pass=pass&&earlier.trajectory.epochS==2.5&&earlier.trajectory.radiusM==1e10+5e5&&
      earlier.speedMPerS==1e6&&later.trajectory.radiusM==1e10+98e6&&
      repeated.trajectory.radiusM==earlier.trajectory.radiusM&&configuration.launchRadiusM==1e10;
  return Finish(pass,"constant-speed departure views match independent r0+v(t-t0) and ignore unrelated query order without a mutable cache/clock");
}
Result LinearSpeed() {
  auto configuration=Constant();configuration.launchEpochS=0;configuration.launchRadiusM=1e9;
  configuration.knots[0].radiusM=1e9;configuration.knots[0].speedMPerS=2e5;
  configuration.knots[1].radiusM=2e9;configuration.knots[1].speedMPerS=4e5;
  S::SourceEpochView view;bool pass=S::EvaluateSourceAtEpoch(configuration,100,&view).ok();
  // v(r)=2e-4*r has independent exponential r(t)=r0*exp(2e-4*t).
  // This interval remains below the next knot, so no endpoint extension is
  // part of the reference or an unapproved source extrapolation choice.
  const double expectedRadius=1e9*std::exp(0.02),expectedSpeed=2e5*std::exp(0.02);
  pass=pass&&std::fabs(view.trajectory.radiusM/expectedRadius-1)<=2e-14&&
      std::fabs(view.speedMPerS/expectedSpeed-1)<=2e-14&&S::ValidateSourceEpochView(view).ok();
  return Finish(pass,"explicit linear-speed trajectory view agrees with independent exponential radius and speed at the same epoch");
}
Result AtomicInvalidEpoch() {
  const auto configuration=Constant();S::SourceEpochView view;
  bool pass=S::EvaluateSourceAtEpoch(configuration,3,&view).ok();
  const auto original=view;
  for(double epoch:{1.0,std::numeric_limits<double>::quiet_NaN()}) {
    const auto rejected=S::EvaluateSourceAtEpoch(configuration,epoch,&view);
    pass=pass&&!rejected.ok()&&view.trajectory.radiusM==original.trajectory.radiusM&&
        view.trajectory.epochS==original.trajectory.epochS&&view.speedMPerS==original.speedMPerS;
  }
  auto invalid=configuration;invalid.knots[1].radiusM=invalid.knots[0].radiusM;
  pass=pass&&!S::EvaluateSourceAtEpoch(invalid,3,&view).ok()&&view.trajectory.radiusM==original.trajectory.radiusM;
  return Finish(pass,"before-launch/nonfinite epochs and invalid governing-law inputs reject without replacing an existing source view");
}
Result MixedEpochRejection() {
  S::SourceEpochView view;bool pass=S::EvaluateSourceAtEpoch(Constant(),3,&view).ok();
  auto wrong=view;wrong.trajectory.epochS=4;pass=pass&&!S::ValidateSourceEpochView(wrong).ok();
  wrong=view;wrong.speedMPerS*=1.01;pass=pass&&!S::ValidateSourceEpochView(wrong).ok();
  wrong=view;wrong.trajectory.radiusM+=1;pass=pass&&!S::ValidateSourceEpochView(wrong).ok();
  return Finish(pass,"mixed epoch/radius/speed source views reject instead of silently substituting a different trajectory state");
}
Result PreparedOwnerBinding() {
  const auto configuration=Constant();S::SourceEpochView owned;
  bool pass=S::EvaluateSourceAtEpoch(configuration,2,&owned).ok();
  const auto owner=owned.trajectory;
  std::string before,after;
  pass=pass&&S::SerializeTrajectory(owner,&before).ok()&&
      S::ValidateSourceDepartureView(owned,&owner,2).ok()&&
      !S::ValidateSourceDepartureView(owned,nullptr,2).ok();
  // Different future speed knots can share the exact same departure radius
  // and speed. This distinguishes full governing-law ownership from a weak
  // point-only comparison while leaving the general math evaluator legitimate.
  auto foreignConfiguration=configuration;foreignConfiguration.knots[1].speedMPerS*=1.5;
  S::SourceEpochView foreign,later;
  pass=pass&&S::EvaluateSourceAtEpoch(foreignConfiguration,2,&foreign).ok()&&
      S::ValidateSourceEpochView(foreign).ok()&&foreign.trajectory.radiusM==owned.trajectory.radiusM&&
      foreign.speedMPerS==owned.speedMPerS&&!S::ValidateSourceDepartureView(foreign,&owner,2).ok();
  pass=pass&&S::EvaluateSourceAtEpoch(configuration,3,&later).ok()&&
      !S::ValidateSourceDepartureView(later,&owner,2).ok()&&
      !S::ValidateSourceDepartureView(owned,&later.trajectory,3).ok()&&
      !S::ValidateSourceDepartureView(owned,&owner,3).ok()&&
      !S::ValidateSourceDepartureView(owned,&owner,std::numeric_limits<double>::quiet_NaN()).ok();
  auto corruptedOwner=owner;corruptedOwner.radiusM+=1;
  pass=pass&&!S::ValidateSourceDepartureView(owned,&corruptedOwner,2).ok()&&
      S::SerializeTrajectory(owner,&after).ok()&&before==after&&owned.trajectory.radiusM==owner.radiusM;
  return Finish(pass,"prepared source authority rejects null/foreign/same-point-different-law/stale-epoch owners without altering inputs; arbitrary-time pure evaluation remains valid");
}
} // namespace
std::vector<Descriptor> AnalyticShockEpochDescriptors() {
  const char* names[]={"Explicit source-departure trajectory values","Independent linear-speed source epoch","Atomic invalid source-epoch rejection","Mixed source-view provenance rejection","Prepared source-owner/departure provenance"};
  const TestCallback callbacks[]={ExplicitDeparture,LinearSpeed,AtomicInvalidEpoch,MixedEpochRejection,PreparedOwnerBinding};
  std::vector<Descriptor> result;
  for(unsigned index=0;index<5;++index) {
    Descriptor d;d.id="ASE0"+std::to_string(index+1);d.name=names[index];d.description=names[index];
    d.group="analytic-shock-epoch";d.initialization=InitializationLevel::None;d.runtime=RuntimeClass::Routine;
    d.seedCapability=SeedCapability::Deterministic;d.repeatSupported=true;d.seedPolicy="independent closed-form references; no random draws";
    d.stateIsolation="private immutable trajectory configuration/views; no native cache or PIC clock";
    d.supportedBuildModes="dependency-light and native identical callbacks";d.callback=callbacks[index];result.push_back(d);
  }
  return result;
}
} }
