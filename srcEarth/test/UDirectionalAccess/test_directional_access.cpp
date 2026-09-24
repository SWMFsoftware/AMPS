#include "../../util/AdaptiveDirectAccess.h"
#include "../../util/DirectionalAccess.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

int failures=0;

void Check(bool condition,const std::string& message) {
  if (!condition) {
    ++failures;
    std::cerr << "FAIL: " << message << "\n";
  }
}

bool Near(double actual,double expected,double relative=1.0e-12,
          double absolute=1.0e-12) {
  return std::fabs(actual-expected)<=absolute+
      relative*std::max(std::fabs(actual),std::fabs(expected));
}

std::vector<std::size_t> Realized(
    const EarthUtil::AdaptiveDirectAccessGrid& grid,
    const std::vector<int>& states) {
  std::vector<std::size_t> result;
  for (std::size_t i=0;i<grid.candidate_GV.size();++i)
    if (states[i]>=0) result.push_back(i);
  return result;
}

Earth::Trajectory::ExitState MakeExitState(double rigidity) {
  Earth::Trajectory::ExitState state;
  state.x_exit_m[0]=12.0;
  state.x_exit_m[1]=-3.0;
  state.x_exit_m[2]=8.0;
  state.p_exit_SI[0]=2.0e-19;
  state.p_exit_SI[1]=0.0;
  state.p_exit_SI[2]=0.0;
  state.v_exit_unit[0]=1.0;
  state.v_exit_unit[1]=0.0;
  state.v_exit_unit[2]=0.0;
  state.cosAlpha=0.25;
  state.traceTimeAtExit_s=2.5;
  state.rigidityAtExit_GV=rigidity;
  state.valid=true;
  return state;
}

Earth::DirectionalAccess::Sample MakeSample(
    double rigidity,EarthUtil::CutoffSampleState state) {
  Earth::DirectionalAccess::Sample sample;
  sample.rigidity_GV=rigidity;
  sample.energy_MeV=rigidity*100.0;
  sample.directionWeight_sr=0.25;
  sample.responseWeight=1.0;
  sample.state=state;
  if (state==EarthUtil::CutoffSampleState::Allowed) {
    sample.termination=
        Earth::GridlessMode::TrajectoryTermination::OuterBoundaryAllowed;
    sample.exitState=MakeExitState(rigidity);
  }
  else if (state==EarthUtil::CutoffSampleState::PhysicalForbidden) {
    sample.termination=
        Earth::GridlessMode::TrajectoryTermination::InnerBoundaryForbidden;
  }
  else {
    sample.termination=Earth::GridlessMode::TrajectoryTermination::TimeLimit;
  }
  return sample;
}

double FullSphereWeight(double resolution_deg) {
  const int nLon=static_cast<int>(std::floor(360.0/resolution_deg+0.5));
  const int nLat=static_cast<int>(std::floor(180.0/resolution_deg+0.5))+1;
  double sum=0.0;
  for (int ilat=0;ilat<nLat;++ilat)
    for (int ilon=0;ilon<nLon;++ilon)
      sum+=Earth::DirectionalAccess::RegularLonLatCellWeightSr(
          resolution_deg,-90.0+ilat*resolution_deg,resolution_deg);
  return sum;
}

} // namespace

int main() {
  using EarthUtil::AdaptiveDirectAccessControls;
  using EarthUtil::AdaptiveDirectAccessReport;
  using EarthUtil::CutoffSampleState;

  const double fourPi=4.0*std::acos(-1.0);

  // U-F14a: the exact spherical-cell reference must close to 4*pi at two distinct
  // resolutions, including half-height polar caps.  This is an analytical reference,
  // not a self-comparison with another numerical quadrature.
  Check(Near(FullSphereWeight(10.0),fourPi,2.0e-14),
        "10-degree regular sky weights close to analytic 4*pi");
  Check(Near(FullSphereWeight(30.0),fourPi,2.0e-14),
        "30-degree regular sky weights close to analytic 4*pi");
  bool invalidCellRejected=false;
  try {
    (void)Earth::DirectionalAccess::RegularLonLatCellWeightSr(
        361.0,0.0,10.0);
  }
  catch (const std::invalid_argument&) {
    invalidCellRejected=true;
  }
  Check(invalidCellRejected,"invalid angular cells are rejected");

  // U-F14b: compare the production adaptive sampler with a closed-form Heaviside
  // access function A(R)=H(R-4 GV).  The resolved bracket must contain the exact
  // cutoff and meet the requested absolute target without evaluating the full tree.
  const EarthUtil::AdaptiveDirectAccessGrid stepGrid=
      EarthUtil::BuildAdaptiveDirectAccessGrid(std::vector<double>{1.0,8.0},12);
  std::vector<int> stepStates(stepGrid.candidate_GV.size(),-1);
  AdaptiveDirectAccessControls stepControls;
  stepControls.guardDepth=1;
  stepControls.absoluteTolerance_GV=0.01;
  const AdaptiveDirectAccessReport stepReport=
      EarthUtil::EvaluateAdaptiveDirectAccessDirectionDetailed(
          stepGrid,stepControls,stepStates,0,
          [](double r,std::size_t) { return r>=4.0 ? 1 : 0; });
  const std::vector<std::size_t> realized=Realized(stepGrid,stepStates);
  double lo=0.0;
  double hi=0.0;
  for (std::size_t k=0;k+1<realized.size();++k) {
    const std::size_t ia=realized[k];
    const std::size_t ib=realized[k+1];
    if (stepStates[ia]==0 && stepStates[ib]==1) {
      lo=stepGrid.candidate_GV[ia];
      hi=stepGrid.candidate_GV[ib];
    }
  }
  Check(lo<4.0 && hi>=4.0 && hi-lo<=0.01,
        "adaptive bracket contains the analytic 4-GV cutoff at target width");
  Check(stepReport.targetReached &&
        stepReport.evaluations<static_cast<int>(stepGrid.candidate_GV.size()),
        "error-controlled refinement converges without evaluating the full tree");
  Check(EarthUtil::CountAdaptiveDirectAccessSamples(
            stepStates,0,stepStates.size())==realized.size(),
        "reported sparse sample count equals the realized candidate count");

  // Relative-tolerance reference at a different rigidity scale.  This catches an
  // implementation that applies the absolute tolerance in all intervals or evaluates
  // relative tolerance against the wrong unit/coordinate.
  const EarthUtil::AdaptiveDirectAccessGrid relativeGrid=
      EarthUtil::BuildAdaptiveDirectAccessGrid(std::vector<double>{10.0,100.0},16);
  std::vector<int> relativeStates(relativeGrid.candidate_GV.size(),-1);
  AdaptiveDirectAccessControls relativeControls;
  relativeControls.guardDepth=0;
  relativeControls.relativeTolerance=1.0e-3;
  const AdaptiveDirectAccessReport relativeReport=
      EarthUtil::EvaluateAdaptiveDirectAccessDirectionDetailed(
          relativeGrid,relativeControls,relativeStates,0,
          [](double r,std::size_t) { return r>=40.0 ? 1 : 0; });
  Check(relativeReport.targetReached &&
        relativeReport.maxAmbiguousWidth_GV<=0.1,
        "relative tolerance controls the analytic transition bracket");

  // U-F14c: equal seed endpoints do not establish monotonicity.  A guard midpoint
  // exposes this hidden forbidden island; both transitions must then be refined.
  const EarthUtil::AdaptiveDirectAccessGrid pocketGrid=
      EarthUtil::BuildAdaptiveDirectAccessGrid(std::vector<double>{1.0,16.0},11);
  std::vector<int> pocketStates(pocketGrid.candidate_GV.size(),-1);
  AdaptiveDirectAccessControls pocketControls;
  pocketControls.guardDepth=1;
  pocketControls.absoluteTolerance_GV=0.02;
  const AdaptiveDirectAccessReport pocketReport=
      EarthUtil::EvaluateAdaptiveDirectAccessDirectionDetailed(
          pocketGrid,pocketControls,pocketStates,0,
          [](double r,std::size_t) {
            return (r>=3.5 && r<=4.5) ? 0 : 1;
          });
  int forbiddenSamples=0;
  for (int state:pocketStates) if (state==0) ++forbiddenSamples;
  Check(forbiddenSamples>0 && pocketReport.refinedIntervals>1 &&
        pocketReport.targetReached,
        "guard refinement finds and converges a non-monotone forbidden pocket");

  // A genuinely oscillatory analytic penumbra exercises multiple access islands.  No
  // single-cutoff or monotonic shortcut can satisfy the exact five-transition check.
  const EarthUtil::AdaptiveDirectAccessGrid oscillatoryGrid=
      EarthUtil::BuildAdaptiveDirectAccessGrid(
          std::vector<double>{1.0,2.0,4.0,8.0,16.0},10);
  std::vector<int> oscillatoryStates(oscillatoryGrid.candidate_GV.size(),-1);
  AdaptiveDirectAccessControls oscillatoryControls;
  oscillatoryControls.guardDepth=1;
  oscillatoryControls.absoluteTolerance_GV=0.02;
  const AdaptiveDirectAccessReport oscillatoryReport=
      EarthUtil::EvaluateAdaptiveDirectAccessDirectionDetailed(
          oscillatoryGrid,oscillatoryControls,oscillatoryStates,0,
          [](double r,std::size_t) {
            return ((r>=1.3 && r<=1.7) ||
                    (r>=2.7 && r<=3.2) || r>=10.0) ? 1 : 0;
          });
  const std::vector<std::size_t> oscillatoryRealized=
      Realized(oscillatoryGrid,oscillatoryStates);
  int observedTransitions=0;
  for (std::size_t k=0;k+1<oscillatoryRealized.size();++k)
    if (oscillatoryStates[oscillatoryRealized[k]]!=
        oscillatoryStates[oscillatoryRealized[k+1]]) ++observedTransitions;
  Check(observedTransitions==5 && oscillatoryReport.refinedIntervals>5 &&
        oscillatoryReport.targetReached,
        "adaptive sampling preserves all five analytic oscillatory transitions");

  // U-F14d: unresolved support is response weighted, while hard work exhaustion is
  // an explicit failure.  Neither case may be relabelled physical forbidden merely to
  // make an access/convergence gate pass.
  std::vector<int> unresolvedStates(stepGrid.candidate_GV.size(),-1);
  AdaptiveDirectAccessControls unresolvedControls;
  unresolvedControls.guardDepth=1;
  unresolvedControls.absoluteTolerance_GV=0.02;
  unresolvedControls.responseWeight=[](double r) { return r*r; };
  const AdaptiveDirectAccessReport unresolvedReport=
      EarthUtil::EvaluateAdaptiveDirectAccessDirectionDetailed(
          stepGrid,unresolvedControls,unresolvedStates,0,
          [](double r,std::size_t) {
            if (r>=3.0 && r<=5.0) return 2;
            return r>5.0 ? 1 : 0;
          });
  Check(unresolvedReport.responseWeightedUnresolvedSupport>0.0 &&
        unresolvedReport.responseWeightedUnresolvedSupport<=1.0,
        "unresolved support is response weighted and bounded");

  std::vector<int> limitedStates(stepGrid.candidate_GV.size(),-1);
  AdaptiveDirectAccessControls limitedControls;
  limitedControls.guardDepth=1;
  limitedControls.absoluteTolerance_GV=0.02;
  limitedControls.maximumSamples=2;
  const AdaptiveDirectAccessReport limitedReport=
      EarthUtil::EvaluateAdaptiveDirectAccessDirectionDetailed(
          stepGrid,limitedControls,limitedStates,0,
          [](double r,std::size_t) { return r>=4.0 ? 1 : 0; });
  Check(limitedReport.maximumSamplesReached && !limitedReport.targetReached &&
        limitedReport.evaluations==2,
        "maximum-sample exhaustion is an explicit non-converged result");

  const EarthUtil::AdaptiveDirectAccessGrid shallowGrid=
      EarthUtil::BuildAdaptiveDirectAccessGrid(std::vector<double>{1.0,8.0},2);
  std::vector<int> shallowStates(shallowGrid.candidate_GV.size(),-1);
  AdaptiveDirectAccessControls shallowControls;
  shallowControls.absoluteTolerance_GV=1.0e-8;
  const AdaptiveDirectAccessReport shallowReport=
      EarthUtil::EvaluateAdaptiveDirectAccessDirectionDetailed(
          shallowGrid,shallowControls,shallowStates,0,
          [](double r,std::size_t) { return r>=4.0 ? 1 : 0; });
  Check(shallowReport.maximumDepthReached && !shallowReport.targetReached,
        "maximum-depth exhaustion is an explicit non-converged result");

  bool badResponseRejected=false;
  try {
    std::vector<int> badStates(stepGrid.candidate_GV.size(),-1);
    AdaptiveDirectAccessControls badControls;
    badControls.guardDepth=1;
    badControls.responseWeight=[](double) { return -1.0; };
    (void)EarthUtil::EvaluateAdaptiveDirectAccessDirectionDetailed(
        stepGrid,badControls,badStates,0,
        [](double r,std::size_t) { return r>=4.0 ? 1 : 0; });
  }
  catch (const std::runtime_error&) {
    badResponseRejected=true;
  }
  Check(badResponseRejected,"negative detector-response weight is rejected");

  // Determinism is an MPI-storage requirement: identical input must realize the same
  // candidate slots and report independent of prior calls.
  std::vector<int> deterministicA(stepGrid.candidate_GV.size(),-1);
  std::vector<int> deterministicB(stepGrid.candidate_GV.size(),-1);
  const AdaptiveDirectAccessReport deterministicReportA=
      EarthUtil::EvaluateAdaptiveDirectAccessDirectionDetailed(
          stepGrid,stepControls,deterministicA,0,
          [](double r,std::size_t) { return r>=4.0 ? 1 : 0; });
  const AdaptiveDirectAccessReport deterministicReportB=
      EarthUtil::EvaluateAdaptiveDirectAccessDirectionDetailed(
          stepGrid,stepControls,deterministicB,0,
          [](double r,std::size_t) { return r>=4.0 ? 1 : 0; });
  Check(deterministicA==deterministicB &&
        deterministicReportA.evaluations==deterministicReportB.evaluations &&
        Near(deterministicReportA.estimatedError_GV,
             deterministicReportB.estimatedError_GV,0.0,0.0),
        "adaptive candidate realization and report are deterministic");

  // U-F15: reproduce lower/effective/upper cutoff and penumbra from saved rows alone.
  // F,A,F,A,A on R=1..5 gives lower=2, upper=4, one GV of linearly reconstructed
  // allowed support inside [2,4], and therefore effective=3 GV.
  std::vector<Earth::DirectionalAccess::Sample> samples;
  samples.push_back(MakeSample(1.0,CutoffSampleState::PhysicalForbidden));
  samples.push_back(MakeSample(2.0,CutoffSampleState::Allowed));
  samples.push_back(MakeSample(3.0,CutoffSampleState::PhysicalForbidden));
  samples.push_back(MakeSample(4.0,CutoffSampleState::Allowed));
  samples.push_back(MakeSample(5.0,CutoffSampleState::Allowed));
  const Earth::DirectionalAccess::CutoffDiagnostics reconstructed=
      Earth::DirectionalAccess::ReconstructCutoff(samples);
  Check(reconstructed.reconstructable && Near(reconstructed.lower_GV,2.0) &&
        Near(reconstructed.effective_GV,3.0) &&
        Near(reconstructed.effectiveLower_GV,3.0) &&
        Near(reconstructed.effectiveUpper_GV,3.0) &&
        Near(reconstructed.upper_GV,4.0) &&
        Near(reconstructed.penumbraWidth_GV,2.0) &&
        reconstructed.transitions==3 && reconstructed.allowedIntervals==2,
        "saved oscillatory A(E,Omega) reconstructs analytic cutoff diagnostics");

  // An independently calculated unresolved-interval reference checks conservative
  // bounds rather than only field presence.  F,A,U,F,A,A on R=1..6 yields nominal
  // Rc_eff=4.5 GV, lower/upper uncertainty bounds 2.5/4.5 GV, and 2/5 response support.
  std::vector<Earth::DirectionalAccess::Sample> bounded;
  bounded.push_back(MakeSample(1.0,CutoffSampleState::PhysicalForbidden));
  bounded.push_back(MakeSample(2.0,CutoffSampleState::Allowed));
  bounded.push_back(MakeSample(3.0,CutoffSampleState::Unresolved));
  bounded.push_back(MakeSample(4.0,CutoffSampleState::PhysicalForbidden));
  bounded.push_back(MakeSample(5.0,CutoffSampleState::Allowed));
  bounded.push_back(MakeSample(6.0,CutoffSampleState::Allowed));
  const Earth::DirectionalAccess::CutoffDiagnostics boundedResult=
      Earth::DirectionalAccess::ReconstructCutoff(bounded);
  Check(boundedResult.reconstructable &&
        Near(boundedResult.lower_GV,2.0) &&
        Near(boundedResult.upper_GV,5.0) &&
        Near(boundedResult.effective_GV,4.5) &&
        Near(boundedResult.effectiveLower_GV,2.5) &&
        Near(boundedResult.effectiveUpper_GV,4.5) &&
        Near(boundedResult.responseWeightedUnresolvedSupport,0.4),
        "unresolved intervals produce independently derived conservative bounds");

  // Negative contract tests prevent a producer from satisfying Step 5 with only an
  // access bit or with inconsistent phase-space fields.
  bool missingExitRejected=false;
  std::vector<Earth::DirectionalAccess::Sample> invalid=samples;
  invalid[1].exitState.valid=false;
  try {
    (void)Earth::DirectionalAccess::ReconstructCutoff(invalid);
  }
  catch (const std::invalid_argument&) {
    missingExitRejected=true;
  }
  Check(missingExitRejected,
        "allowed sample without its complete boundary state is rejected");

  bool inconsistentExitRejected=false;
  invalid=samples;
  invalid[1].exitState.v_exit_unit[0]=-1.0;
  try {
    (void)Earth::DirectionalAccess::ReconstructCutoff(invalid);
  }
  catch (const std::invalid_argument&) {
    inconsistentExitRejected=true;
  }
  Check(inconsistentExitRejected,
        "oppositely directed exit momentum and velocity are rejected");

  bool terminationMismatchRejected=false;
  invalid=samples;
  invalid[0].termination=
      Earth::GridlessMode::TrajectoryTermination::TimeLimit;
  try {
    (void)Earth::DirectionalAccess::ReconstructCutoff(invalid);
  }
  catch (const std::invalid_argument&) {
    terminationMismatchRejected=true;
  }
  Check(terminationMismatchRejected,
        "physical-forbidden state with unresolved termination is rejected");

  bool orderingRejected=false;
  invalid=samples;
  invalid[2].rigidity_GV=invalid[1].rigidity_GV;
  try {
    (void)Earth::DirectionalAccess::ReconstructCutoff(invalid);
  }
  catch (const std::invalid_argument&) {
    orderingRejected=true;
  }
  Check(orderingRejected,"duplicate/non-increasing rigidity rows are rejected");

  if (failures!=0) {
    std::cerr << "UDirectionalAccess: " << failures << " failure(s)\n";
    return EXIT_FAILURE;
  }
  std::cout << "UDirectionalAccess: PASS (U-F14, U-F15)\n";
  return EXIT_SUCCESS;
}
