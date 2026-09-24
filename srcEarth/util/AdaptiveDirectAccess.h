#ifndef _EARTH_ADAPTIVE_DIRECT_ACCESS_H_
#define _EARTH_ADAPTIVE_DIRECT_ACCESS_H_

//======================================================================================
// AdaptiveDirectAccess.h
//======================================================================================
//
// Shared adaptive rigidity sampler for the directional A(R,Omega) validation product.
//
// Motivation
// ----------
// A dense DIRECT_ACCESS run evaluates the same fixed rigidity list at every sky
// direction.  That is robust but wasteful: most directions are smoothly ALLOWED or
// PHYSICAL_FORBIDDEN across broad portions of the detector-response range, while the
// scientifically interesting structure is concentrated near access transitions and
// unresolved trace-limit brackets.
//
// This helper keeps the validation semantics explicit while reducing trajectory count:
//
//   1. Every user-supplied seed rigidity is always evaluated.  The caller therefore
//      controls the coarse global coverage and must include the detector-response
//      endpoints.
//
//   2. A configurable guard depth probes geometric midpoints even when the two current
//      endpoint states agree.  The default C19 guard depth of one samples the midpoint
//      of every seed interval and substantially reduces the risk of missing a narrow
//      allowed/forbidden pocket whose two coarse endpoints happen to have the same
//      classification.
//
//   3. After the guard probes, refinement continues only when an interval is visibly
//      ambiguous, i.e. its endpoint states differ.  This includes resolved transitions
//      and resolved<->UNRESOLVED boundaries.  An interval whose two endpoints are both
//      UNRESOLVED is *not* recursively exploded to the full tree: repeated rigidity
//      bisection cannot cure a trajectory time/path/step cap, and the Python fold already
//      carries the whole unresolved interval as an uncertainty bound.
//
//   4. The full candidate tree is deterministic and identical on every MPI rank.  A
//      direction evaluates only the nodes it needs; untouched candidate slots remain
//      -1.  This is important because Mode3D and GRIDLESS can MPI_MAX-reduce fixed-size
//      sentinel arrays even though each direction has a different adaptive sample set.
//
// The sampler deliberately does NOT assume monotonic access.  Multiple transitions are
// retained whenever the guard/refinement probes expose them.  The Python C19 fold still
// treats every sampled ALLOWED<->FORBIDDEN interval as an uncertainty bracket instead of
// inventing a fractional transmission ramp.  Thus adaptive sampling reduces work while
// preserving the existing response-weighted convergence gate.
//======================================================================================

#include <algorithm>
#include <cstdint>
#include <type_traits>
#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <stdexcept>
#include <vector>

namespace EarthUtil {

struct AdaptiveDirectAccessGrid {
  std::vector<double> seed_GV;
  std::vector<double> candidate_GV;
  std::vector<std::size_t> seedCandidateIndex;
  int maxDepth{0};
};

inline double AdaptiveDirectAccessMidpointGV(double a,double b) {
  // Geometric midpoint is natural for positive rigidity and preserves comparable
  // fractional resolution across the response span.
  return std::sqrt(a*b);
}

inline void AppendAdaptiveCandidates_(double a,double b,int depth,int maxDepth,
                                      std::vector<double>& values) {
  if (depth>=maxDepth) return;
  const double m=AdaptiveDirectAccessMidpointGV(a,b);
  values.push_back(m);
  AppendAdaptiveCandidates_(a,m,depth+1,maxDepth,values);
  AppendAdaptiveCandidates_(m,b,depth+1,maxDepth,values);
}

inline std::size_t FindAdaptiveDirectAccessNode(const std::vector<double>& grid,
                                                 double value) {
  const auto it=std::lower_bound(grid.begin(),grid.end(),value);
  auto acceptable=[&](std::vector<double>::const_iterator p) {
    if (p==grid.end()) return false;
    const double scale=std::max(1.0,std::fabs(value));
    return std::fabs(*p-value)<=2.0e-13*scale;
  };
  if (acceptable(it)) return static_cast<std::size_t>(it-grid.begin());
  if (it!=grid.begin()) {
    const auto prev=it-1;
    if (acceptable(prev)) return static_cast<std::size_t>(prev-grid.begin());
  }
  throw std::runtime_error("adaptive direct-access candidate node lookup failed");
}

inline AdaptiveDirectAccessGrid BuildAdaptiveDirectAccessGrid(
    const std::vector<double>& seed_GV,int maxDepth) {
  if (seed_GV.size()<2)
    throw std::runtime_error("adaptive direct access requires at least two seed rigidities");
  if (maxDepth<0 || maxDepth>20)
    throw std::runtime_error("adaptive direct-access max depth must be in [0,20]");

  for (std::size_t i=0;i<seed_GV.size();++i) {
    if (!(seed_GV[i]>0.0) || !std::isfinite(seed_GV[i]))
      throw std::runtime_error("adaptive direct-access seed rigidities must be finite and positive");
    if (i>0 && !(seed_GV[i]>seed_GV[i-1]))
      throw std::runtime_error("adaptive direct-access seed rigidities must be strictly increasing");
  }

  AdaptiveDirectAccessGrid result;
  result.seed_GV=seed_GV;
  result.maxDepth=maxDepth;
  result.candidate_GV=seed_GV;
  for (std::size_t i=0;i+1<seed_GV.size();++i)
    AppendAdaptiveCandidates_(seed_GV[i],seed_GV[i+1],0,maxDepth,result.candidate_GV);

  std::sort(result.candidate_GV.begin(),result.candidate_GV.end());
  std::vector<double> unique;
  unique.reserve(result.candidate_GV.size());
  for (double value:result.candidate_GV) {
    if (unique.empty()) {
      unique.push_back(value);
      continue;
    }
    const double scale=std::max(1.0,std::fabs(value));
    if (std::fabs(value-unique.back())>2.0e-13*scale) unique.push_back(value);
  }
  result.candidate_GV.swap(unique);

  result.seedCandidateIndex.reserve(seed_GV.size());
  for (double value:seed_GV)
    result.seedCandidateIndex.push_back(
        FindAdaptiveDirectAccessNode(result.candidate_GV,value));
  return result;
}

// Per-trajectory diagnostics for the sparse DIRECT_ACCESS product.
//
// ``slot`` is the global flattened [location][sky-cell][candidate-rigidity] index.
// Only actually evaluated adaptive nodes create records, so the memory footprint scales
// with the number of trajectories rather than with the full depth-6 candidate tree.
// The struct intentionally contains only POD fields: Mode3D and GRIDLESS gather it as
// raw MPI_BYTE records after all worker threads have joined.
struct DirectAccessSampleDiagnostic {
  std::uint64_t slot{0};
  int terminationCode{-1};
  double traceTime_s{0.0};
  double traceDistance_Re{0.0};
  int steps{0};
  int retryCount{0};

  // Trace-budget convergence provenance.  These fields let the C19 postprocessor
  // distinguish the normal 300-s classification from an unresolved-only extended
  // result without rerunning or guessing from the final trace time.  A value of
  // traceExtensionCount==0 means the normal primary budget was sufficient.
  int primaryTerminationCode{-1};
  double primaryTraceTime_s{0.0};
  int traceExtensionCount{0};
  double initialTraceLimit_s{0.0};
  double finalTraceLimit_s{0.0};

  int mirrorPoints{0};
  int bounceCycles{0};
  int driftRevolutions{0};
  double driftAngle_deg{0.0};
  double driftMeanRadiusChange_Re{0.0};
  int trapMechanism{0};       // 0=None, 1=Bounce, 2=Drift
  double momentumRelativeSpread{0.0};

  // Complete outer-boundary state for an ALLOWED characteristic.  These members are
  // deliberately stored in the sparse per-sample record rather than in a parallel
  // array: MPI gathers DirectAccessSampleDiagnostic as raw bytes, and keeping the
  // state beside its deterministic slot prevents an adaptive row from being matched
  // to the wrong trajectory after rank reduction.  Non-allowed samples retain the
  // zero initialization and exitStateValid==0.
  int exitStateValid{0};
  double xExit_m[3]{0.0,0.0,0.0};
  double pExit_SI[3]{0.0,0.0,0.0};
  double vExitUnit[3]{0.0,0.0,0.0};
  double cosAlphaExit{0.0};
  double traceTimeAtExit_s{0.0};
  double rigidityAtExit_GV{0.0};

  // Convergence metadata is repeated on every realized row in one direction.  This
  // makes a saved A(E,Omega) curve self-describing even when rows are split, filtered,
  // or consumed without the original input deck.  Dense runs use the defaults below.
  int adaptiveRefinedIntervals{0};
  double adaptiveEstimatedError_GV{0.0};
  double adaptiveMaxAmbiguousWidth_GV{0.0};
  int adaptiveTargetReached{1};
  int adaptiveMaxSamplesReached{0};
  double responseWeightedUnresolvedSupport{0.0};
};

// Diagnostics are gathered with MPI_BYTE rather than a custom MPI datatype.  Keep
// this compile-time guard next to the record definition so adding a non-POD member
// cannot silently make the byte-wise gather invalid.
static_assert(std::is_trivially_copyable<DirectAccessSampleDiagnostic>::value,
              "DirectAccessSampleDiagnostic must remain trivially copyable");

// Error-control policy for one sky direction.  Depth and sample count are hard work
// bounds; the absolute/relative tolerances are scientific convergence targets.  A run
// that exhausts a hard bound before reaching its target is reported as non-converged --
// it is never silently promoted to a resolved access curve.
struct AdaptiveDirectAccessControls {
  int guardDepth{1};
  double absoluteTolerance_GV{0.0};
  double relativeTolerance{0.0};
  int maximumSamples{0}; // 0 = no cap beyond the deterministic candidate tree

  // Optional detector/spectrum weight for unresolved-support accounting.  It does not
  // affect which trajectory states are computed and therefore cannot bias refinement
  // toward a desired validation result.  Unity is used when no response is supplied.
  std::function<double(double)> responseWeight;
};

struct AdaptiveDirectAccessReport {
  int evaluations{0};
  int refinedIntervals{0};
  double estimatedError_GV{0.0};
  double maxAmbiguousWidth_GV{0.0};
  double responseWeightedUnresolvedSupport{0.0};
  bool targetReached{true};
  bool maximumSamplesReached{false};
  bool maximumDepthReached{false};
};

// Detailed Step-5 sampler.  The state vector is a fixed candidate-tree slice whose
// unevaluated entries remain -1.  The classifier receives both rigidity and candidate
// index so production callers can attach diagnostics to the exact global slot without
// a floating-point lookup.
template<class Classifier>
inline AdaptiveDirectAccessReport EvaluateAdaptiveDirectAccessDirectionDetailed(
    const AdaptiveDirectAccessGrid& grid,
    const AdaptiveDirectAccessControls& controls,
    std::vector<int>& states,
    std::size_t base,
    Classifier classify,
    int unresolvedState=2) {
  if (controls.guardDepth<0 || controls.guardDepth>grid.maxDepth)
    throw std::runtime_error("adaptive direct-access guard depth must be in [0,maxDepth]");
  if (controls.maximumSamples>0 &&
      controls.maximumSamples<static_cast<int>(grid.seedCandidateIndex.size()))
    throw std::runtime_error(
        "adaptive direct-access maximum samples cannot be smaller than the seed count");
  if (controls.absoluteTolerance_GV<0.0 || controls.relativeTolerance<0.0 ||
      !std::isfinite(controls.absoluteTolerance_GV) ||
      !std::isfinite(controls.relativeTolerance))
    throw std::runtime_error(
        "adaptive direct-access tolerances must be finite and non-negative");
  if (base+grid.candidate_GV.size()>states.size())
    throw std::runtime_error("adaptive direct-access state slice exceeds output array");

  AdaptiveDirectAccessReport report;
  const auto canEvaluate=[&]() {
    return controls.maximumSamples<=0 || report.evaluations<controls.maximumSamples;
  };
  const auto stateAt=[&](std::size_t idx) -> int {
    int& slot=states[base+idx];
    if (slot<0) {
      if (!canEvaluate()) {
        report.maximumSamplesReached=true;
        // Do not write a synthetic state into the product.  The missing -1 slot is
        // intentionally omitted by the sparse writer, while targetReached=false on
        // every realized row records that the requested curve was not converged.
        return unresolvedState;
      }
      slot=classify(grid.candidate_GV[idx],idx);
      ++report.evaluations;
    }
    return slot;
  };

  // Coarse seeds establish common response support for every direction.  The parser
  // prevents maximumSamples from being smaller than this mandatory set.
  for (std::size_t idx:grid.seedCandidateIndex) (void)stateAt(idx);

  std::function<void(double,std::size_t,double,std::size_t,int)> refine;
  refine=[&](double a,std::size_t ia,double b,std::size_t ib,int depth) {
    const int sa=stateAt(ia);
    const int sb=stateAt(ib);
    // Keep the explicit name `visibleAmbiguity`: C19's unchanged architecture gate
    // verifies that refinement is driven only by a sampled state change and never by
    // an assumed monotonic cutoff hidden from the saved access curve.
    const bool visibleAmbiguity=(sa!=sb);
    const bool guardProbe=(depth<controls.guardDepth);
    if (!visibleAmbiguity && !guardProbe) return;

    const double width=b-a;
    const double tolerance=std::max(
        controls.absoluteTolerance_GV,
        controls.relativeTolerance*std::max(std::fabs(a),std::fabs(b)));
    if (tolerance>0.0 && width<=tolerance) return;
    if (depth>=grid.maxDepth) {
      if (visibleAmbiguity) report.maximumDepthReached=true;
      return;
    }
    if (!canEvaluate()) {
      report.maximumSamplesReached=true;
      return;
    }

    const double m=AdaptiveDirectAccessMidpointGV(a,b);
    const std::size_t im=FindAdaptiveDirectAccessNode(grid.candidate_GV,m);
    (void)stateAt(im);
    ++report.refinedIntervals;

    // Each child makes its own guard/ambiguity decision.  Consequently multiple
    // allowed islands and non-monotone penumbrae survive; the algorithm never reduces
    // the curve to a single assumed cutoff during trajectory generation.
    refine(a,ia,m,im,depth+1);
    refine(m,im,b,ib,depth+1);
  };

  for (std::size_t i=0;i+1<grid.seed_GV.size();++i) {
    refine(grid.seed_GV[i],grid.seedCandidateIndex[i],
           grid.seed_GV[i+1],grid.seedCandidateIndex[i+1],0);
  }

  // Construct conservative a-posteriori indicators from realized adjacent samples.
  // Every unequal-state interval contributes its entire width to the access-error
  // support.  Intervals touching UNRESOLVED additionally contribute their weighted
  // width to the observable-specific unresolved fraction.
  std::vector<std::size_t> realized;
  realized.reserve(static_cast<std::size_t>(report.evaluations));
  for (std::size_t i=0;i<grid.candidate_GV.size();++i)
    if (states[base+i]>=0) realized.push_back(i);

  double totalWeightedWidth=0.0;
  double unresolvedWeightedWidth=0.0;
  for (std::size_t k=0;k+1<realized.size();++k) {
    const std::size_t ia=realized[k];
    const std::size_t ib=realized[k+1];
    const double a=grid.candidate_GV[ia];
    const double b=grid.candidate_GV[ib];
    const double width=b-a;
    const int sa=states[base+ia];
    const int sb=states[base+ib];

    double weight=1.0;
    if (controls.responseWeight) {
      weight=controls.responseWeight(AdaptiveDirectAccessMidpointGV(a,b));
      if (!(weight>=0.0) || !std::isfinite(weight))
        throw std::runtime_error(
            "adaptive direct-access response weight must be finite and non-negative");
    }
    totalWeightedWidth+=weight*width;
    if (sa!=sb) {
      report.estimatedError_GV+=width;
      report.maxAmbiguousWidth_GV=std::max(report.maxAmbiguousWidth_GV,width);
    }
    if (sa==unresolvedState || sb==unresolvedState)
      unresolvedWeightedWidth+=weight*width;
  }
  if (totalWeightedWidth>0.0)
    report.responseWeightedUnresolvedSupport=
        unresolvedWeightedWidth/totalWeightedWidth;

  // The recursion itself uses a local relative tolerance.  Comparing the final width
  // with the largest requested tolerance is safe here because any premature stop also
  // sets one of the explicit hard-limit flags below.
  const double globalTolerance=std::max(
      controls.absoluteTolerance_GV,
      controls.relativeTolerance*std::max(
          std::fabs(grid.candidate_GV.front()),
          std::fabs(grid.candidate_GV.back())));
  report.targetReached=!report.maximumSamplesReached &&
      !report.maximumDepthReached &&
      (globalTolerance<=0.0 || report.maxAmbiguousWidth_GV<=globalTolerance);
  return report;
}

template<class Classifier>
inline int EvaluateAdaptiveDirectAccessDirection(
    const AdaptiveDirectAccessGrid& grid,
    int guardDepth,
    std::vector<int>& states,
    std::size_t base,
    Classifier classify,
    int unresolvedState=2) {
  // Compatibility wrapper: zero tolerances preserve the historical depth/guard-only
  // evaluation pattern.  Existing callers and C19 convergence gates therefore keep
  // their exact sampling semantics until they opt into the Step-5 controls.
  AdaptiveDirectAccessControls controls;
  controls.guardDepth=guardDepth;
  return EvaluateAdaptiveDirectAccessDirectionDetailed(
      grid,controls,states,base,classify,unresolvedState).evaluations;
}

inline std::size_t CountAdaptiveDirectAccessSamples(const std::vector<int>& states,
                                                     std::size_t base,
                                                     std::size_t count) {
  if (base+count>states.size())
    throw std::runtime_error("adaptive direct-access sample-count slice exceeds array");
  std::size_t n=0;
  for (std::size_t i=0;i<count;++i) if (states[base+i]>=0) ++n;
  return n;
}

} // namespace EarthUtil

#endif
