#include "../../util/AdaptiveDirectAccess.h"
#include "../../util/DirectionalAccess.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
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

std::vector<std::size_t> Realized(const EarthUtil::AdaptiveDirectAccessGrid& grid,
                                  const std::vector<int>& states) {
  std::vector<std::size_t> result;
  for (std::size_t i=0;i<grid.candidate_GV.size();++i)
    if (states[i]>=0) result.push_back(i);
  return result;
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
    sample.termination=Earth::GridlessMode::TrajectoryTermination::OuterBoundaryAllowed;
    sample.exitState.valid=true;
    sample.exitState.rigidityAtExit_GV=rigidity;
  }
  else if (state==EarthUtil::CutoffSampleState::PhysicalForbidden) {
    sample.termination=Earth::GridlessMode::TrajectoryTermination::InnerBoundaryForbidden;
  }
  else {
    sample.termination=Earth::GridlessMode::TrajectoryTermination::TimeLimit;
  }
  return sample;
}

} // namespace

int main() {
  using EarthUtil::AdaptiveDirectAccessControls;
  using EarthUtil::AdaptiveDirectAccessReport;
  using EarthUtil::CutoffSampleState;

  // The exact spherical-cell formula must close to 4*pi on a complete regular sky,
  // including the half-height polar caps. This is the angular normalization used by
  // both production DIRECT_ACCESS writers.
  double fullSphereWeight=0.0;
  const double resolution_deg=10.0;
  for (int ilat=0;ilat<=18;++ilat)
    for (int ilon=0;ilon<36;++ilon)
      fullSphereWeight+=Earth::DirectionalAccess::RegularLonLatCellWeightSr(
          resolution_deg,-90.0+ilat*resolution_deg,resolution_deg);
  Check(Near(fullSphereWeight,4.0*std::acos(-1.0),2.0e-14),
        "regular lon/lat direction weights close to 4*pi");

  // U-F13a: an analytic Heaviside access function is the reference solution.  The
  // final realized forbidden/allowed bracket must contain Rc=4 GV and satisfy the
  // requested absolute error target without evaluating the complete candidate tree.
  const EarthUtil::AdaptiveDirectAccessGrid stepGrid=
      EarthUtil::BuildAdaptiveDirectAccessGrid(std::vector<double>{1.0,8.0},12);
  std::vector<int> stepStates(stepGrid.candidate_GV.size(),-1);
  AdaptiveDirectAccessControls controls;
  controls.guardDepth=1;
  controls.absoluteTolerance_GV=0.01;
  const AdaptiveDirectAccessReport stepReport=
      EarthUtil::EvaluateAdaptiveDirectAccessDirectionDetailed(
          stepGrid,controls,stepStates,0,
          [](double r,std::size_t) { return r>=4.0 ? 1 : 0; });
  const std::vector<std::size_t> realized=Realized(stepGrid,stepStates);
  double lo=0.0,hi=0.0;
  for (std::size_t k=0;k+1<realized.size();++k) {
    const std::size_t ia=realized[k],ib=realized[k+1];
    if (stepStates[ia]==0 && stepStates[ib]==1) {
      lo=stepGrid.candidate_GV[ia];
      hi=stepGrid.candidate_GV[ib];
    }
  }
  Check(lo<4.0 && hi>=4.0 && hi-lo<=0.01,
        "adaptive step-function bracket contains the analytic cutoff at target width");
  Check(stepReport.targetReached &&
        stepReport.evaluations<static_cast<int>(stepGrid.candidate_GV.size()),
        "error-driven refinement converges without evaluating the full tree");

  // U-F13b: equal seed endpoints do not imply monotonicity.  A guard midpoint exposes
  // this narrow forbidden island, after which both transitions are refined.
  const EarthUtil::AdaptiveDirectAccessGrid pocketGrid=
      EarthUtil::BuildAdaptiveDirectAccessGrid(std::vector<double>{1.0,16.0},11);
  std::vector<int> pocketStates(pocketGrid.candidate_GV.size(),-1);
  controls.absoluteTolerance_GV=0.02;
  const AdaptiveDirectAccessReport pocketReport=
      EarthUtil::EvaluateAdaptiveDirectAccessDirectionDetailed(
          pocketGrid,controls,pocketStates,0,
          [](double r,std::size_t) { return (r>=3.5 && r<=4.5) ? 0 : 1; });
  int forbiddenSamples=0;
  for (int state:pocketStates) if (state==0) ++forbiddenSamples;
  Check(forbiddenSamples>0 && pocketReport.refinedIntervals>1,
        "guard refinement finds a non-monotone pocket with equal endpoint states");

  // U-F13c: an analytic multi-band curve exercises a genuinely oscillatory penumbra.
  // Each low-energy allowed island is hidden between forbidden seed endpoints and is
  // exposed by its guard midpoint; the final high-energy branch is endpoint-visible.
  const EarthUtil::AdaptiveDirectAccessGrid oscillatoryGrid=
      EarthUtil::BuildAdaptiveDirectAccessGrid(
          std::vector<double>{1.0,2.0,4.0,8.0,16.0},10);
  std::vector<int> oscillatoryStates(oscillatoryGrid.candidate_GV.size(),-1);
  controls.absoluteTolerance_GV=0.02;
  const AdaptiveDirectAccessReport oscillatoryReport=
      EarthUtil::EvaluateAdaptiveDirectAccessDirectionDetailed(
          oscillatoryGrid,controls,oscillatoryStates,0,
          [](double r,std::size_t) {
            return ((r>=1.3 && r<=1.7) || (r>=2.7 && r<=3.2) || r>=10.0)
                ? 1 : 0;
          });
  const std::vector<std::size_t> oscillatoryRealized=
      Realized(oscillatoryGrid,oscillatoryStates);
  int observedTransitions=0;
  for (std::size_t k=0;k+1<oscillatoryRealized.size();++k)
    if (oscillatoryStates[oscillatoryRealized[k]]!=
        oscillatoryStates[oscillatoryRealized[k+1]]) ++observedTransitions;
  Check(observedTransitions==5 && oscillatoryReport.refinedIntervals>5,
        "adaptive sampling preserves all five analytic oscillatory transitions");

  // U-F13d: unresolved intervals contribute an explicit response-weighted support;
  // a hard trajectory budget reports failure instead of masquerading as convergence.
  std::vector<int> unresolvedStates(stepGrid.candidate_GV.size(),-1);
  controls.absoluteTolerance_GV=0.02;
  controls.responseWeight=[](double r) { return r*r; };
  const AdaptiveDirectAccessReport unresolvedReport=
      EarthUtil::EvaluateAdaptiveDirectAccessDirectionDetailed(
          stepGrid,controls,unresolvedStates,0,
          [](double r,std::size_t) {
            if (r>=3.0 && r<=5.0) return 2;
            return r>5.0 ? 1 : 0;
          });
  Check(unresolvedReport.responseWeightedUnresolvedSupport>0.0 &&
        unresolvedReport.responseWeightedUnresolvedSupport<=1.0,
        "unresolved support is response weighted and bounded");

  std::vector<int> limitedStates(stepGrid.candidate_GV.size(),-1);
  controls.responseWeight=std::function<double(double)>();
  controls.maximumSamples=2;
  const AdaptiveDirectAccessReport limitedReport=
      EarthUtil::EvaluateAdaptiveDirectAccessDirectionDetailed(
          stepGrid,controls,limitedStates,0,
          [](double r,std::size_t) { return r>=4.0 ? 1 : 0; });
  Check(limitedReport.maximumSamplesReached && !limitedReport.targetReached &&
        limitedReport.evaluations==2,
        "maximum-sample exhaustion is an explicit non-converged result");

  // U-F13e: reproduce lower/effective/upper cutoff and penumbra from the saved
  // directional-access rows alone.  The curve contains two allowed islands:
  // F,A,F,A,A on R=1..5, giving lower=2, upper=4, one GV of allowed support in the
  // [2,4] penumbra, and therefore effective=3 GV.
  std::vector<Earth::DirectionalAccess::Sample> samples;
  samples.push_back(MakeSample(1.0,CutoffSampleState::PhysicalForbidden));
  samples.push_back(MakeSample(2.0,CutoffSampleState::Allowed));
  samples.push_back(MakeSample(3.0,CutoffSampleState::PhysicalForbidden));
  samples.push_back(MakeSample(4.0,CutoffSampleState::Allowed));
  samples.push_back(MakeSample(5.0,CutoffSampleState::Allowed));
  const Earth::DirectionalAccess::CutoffDiagnostics reconstructed=
      Earth::DirectionalAccess::ReconstructCutoff(samples);
  Check(reconstructed.reconstructable && Near(reconstructed.lower_GV,2.0) &&
        Near(reconstructed.effective_GV,3.0) && Near(reconstructed.upper_GV,4.0) &&
        Near(reconstructed.penumbraWidth_GV,2.0),
        "saved A(E,Omega) reconstructs lower/effective/upper cutoff and penumbra");

  bool missingExitRejected=false;
  samples[1].exitState.valid=false;
  try {
    (void)Earth::DirectionalAccess::ReconstructCutoff(samples);
  }
  catch (const std::invalid_argument&) {
    missingExitRejected=true;
  }
  Check(missingExitRejected,
        "an allowed saved sample without its boundary exit state is rejected");

  if (failures!=0) {
    std::cerr << "UDirectionalAccess: " << failures << " failure(s)\n";
    return EXIT_FAILURE;
  }
  std::cout << "UDirectionalAccess: PASS (U-F13)\n";
  return EXIT_SUCCESS;
}
