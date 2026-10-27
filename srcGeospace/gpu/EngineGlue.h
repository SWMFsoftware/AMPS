#ifndef _SRC_GEOSPACE_GPU_ENGINE_GLUE_H_
#define _SRC_GEOSPACE_GPU_ENGINE_GLUE_H_

//======================================================================================
// EngineGlue.h -- conversions between the legacy trajectory contract and the
// engine POD types, plus the optional bitwise verification against the legacy tracer.
//
// AMPS_EARTH_ENGINE_VERIFY=1 makes every engine trajectory also run through the original
// legacy tracer on the calling thread and compares termination, step count, time,
// path length, exit state and trap diagnostics bit for bit.  Mismatches are printed
// (first 20 per rank) and summarized at exit.
//======================================================================================

#include "BatchTracer.h"
#include "../util/TrajectoryContract.h"
#include "../util/amps_param_parser.h"

#include <cctype>
#include <cstdlib>
#include <cstring>
#include <string>

namespace Geospace {
namespace GPU {

// GeospaceUtil::ToUpper(s)==upper without allocating: the same std::toupper per byte
// (ASCII "C" locale) followed by an exact comparison.  Used on per-trajectory paths.
inline bool EqualsUpperAscii(const std::string& s,const char* upper) {
  std::size_t i=0;
  for (;i<s.size();++i) {
    if (upper[i]=='\0') return false;
    if (static_cast<char>(std::toupper(static_cast<unsigned char>(s[i])))!=upper[i]) return false;
  }
  return upper[i]=='\0';
}

inline EgTrapConfig MakeTrapConfig(const GeospaceUtil::AmpsParam& prm,double earthRadius_m) {
  EgTrapConfig t;
  std::memset(&t,0,sizeof(t));
  t.enabled=prm.numerics.trapDetection ? 1 : 0;
  t.minMirrorPoints=prm.numerics.trapMinMirrorPoints;
  t.minBounceCycles=prm.numerics.trapMinBounceCycles;
  t.outerMargin_m=prm.numerics.trapOuterMargin_Re*earthRadius_m;
  t.radialEnvelopeTolerance_m=prm.numerics.trapRadialGrowthTolerance_Re*earthRadius_m;
  t.energyRelativeTolerance=prm.numerics.trapEnergyRelativeTolerance;
  t.parallelDeadband=prm.numerics.trapParallelDeadband;
  t.driftEnabled=prm.numerics.trapDriftDetection ? 1 : 0;
  t.minDriftRevolutions=prm.numerics.trapMinDriftRevolutions;
  t.driftRadialAbsoluteTolerance_m=prm.numerics.trapDriftRadialGrowthTolerance_Re*earthRadius_m;
  t.driftRadialRelativeTolerance=prm.numerics.trapDriftRadialRelativeTolerance;
  t.driftLatitudeTolerance=prm.numerics.trapDriftLatitudeTolerance;
  t.driftPitchCos2Tolerance=prm.numerics.trapDriftPitchCos2Tolerance;
  t.driftMaxMeanRadiusChange_m=prm.numerics.trapDriftMaxMeanRadiusChange_Re*earthRadius_m;
  t.driftProfileBins=prm.numerics.trapDriftProfileBins;
  t.driftMinProfileCoverage=prm.numerics.trapDriftMinProfileCoverage;
  t.driftMinMatchedBinFraction=prm.numerics.trapDriftMinMatchedBinFraction;
  t.driftWindowRevolutions=prm.numerics.trapDriftWindowRevolutions;
  t.driftMinStableWindowComparisons=prm.numerics.trapDriftMinStableWindowComparisons;
  t.driftMaxRecurrenceLagRevolutions=prm.numerics.trapDriftMaxRecurrenceLagRevolutions;
  // The legacy tracer reads AMPS_TRAP_TURN_DIAGNOSTICS per trajectory; the environment does
  // not change during a run, so it is read once.
  static const int recordTurns=[]() {
    const char* diag=std::getenv("AMPS_TRAP_TURN_DIAGNOSTICS");
    return (diag!=nullptr && std::string(diag)=="1") ? 1 : 0;
  }();
  t.recordDriftTurnStates=recordTurns;
  return t;
}

// Copy the engine result into a legacy result that already carries the provenance
// fields (mover, backward-time mode, snapshot fingerprint) set by the caller.
inline void FillTrajectoryResult(const EgResult& r,const std::vector<EgTurnState>& turns,
                                 Geospace::Trajectory::Result& out) {
  out.termination=static_cast<Geospace::GridlessMode::TrajectoryTermination>(r.termination);
  out.traceTime_s=r.traceTime_s;
  out.traceDistance_m=r.traceDistance_m;
  out.steps=r.steps;
  out.mirrorPoints=r.mirrorPoints;
  out.bounceCycles=r.bounceCycles;
  out.driftRevolutions=r.driftRevolutions;
  out.driftAngle_rad=r.driftAngle_rad;
  out.driftMeanRadiusChange_m=r.driftMeanRadiusChange_m;
  out.driftMeanRadiusSpan_m=r.driftMeanRadiusSpan_m;
  out.driftStableComparisons=r.driftStableComparisons;
  out.driftRecurrenceLagRevolutions=r.driftRecurrenceLagRevolutions;
  out.minimumOuterMargin_m=r.minimumOuterMargin_m;
  out.lastBounceEnvelopeChange_m=r.lastBounceEnvelopeChange_m;
  out.driftComparisonAvailable=r.driftComparisonAvailable!=0;
  out.driftPreviousCoveragePassed=r.driftPreviousCoveragePassed!=0;
  out.driftCurrentCoveragePassed=r.driftCurrentCoveragePassed!=0;
  out.driftSecularDriftPassed=r.driftSecularDriftPassed!=0;
  out.driftProfilesRecurred=r.driftProfilesRecurred!=0;
  out.driftPreviousCoverageFraction=r.driftPreviousCoverageFraction;
  out.driftCurrentCoverageFraction=r.driftCurrentCoverageFraction;
  out.driftCommonBins=r.driftCommonBins;
  out.driftMatchedBins=r.driftMatchedBins;
  out.driftRadialMismatchBins=r.driftRadialMismatchBins;
  out.driftLatitudeMismatchBins=r.driftLatitudeMismatchBins;
  out.driftPitchMismatchBins=r.driftPitchMismatchBins;
  out.driftMatchedBinFraction=r.driftMatchedBinFraction;
  out.driftMaximumNormalizedRadialResidual=r.driftMaximumNormalizedRadialResidual;
  out.driftMaximumNormalizedLatitudeResidual=r.driftMaximumNormalizedLatitudeResidual;
  out.driftMaximumNormalizedPitchResidual=r.driftMaximumNormalizedPitchResidual;
  out.trapMechanism=r.trapMechanism;
  out.momentumRelativeSpread=r.momentumRelativeSpread;
  for (int d=0;d<3;++d) {
    out.exitState.x_exit_m[d]=r.exit.x_exit_m[d];
    out.exitState.p_exit_SI[d]=r.exit.p_exit_SI[d];
    out.exitState.v_exit_unit[d]=r.exit.v_exit_unit[d];
  }
  out.exitState.cosAlpha=r.exit.cosAlpha;
  out.exitState.traceTimeAtExit_s=r.exit.traceTimeAtExit_s;
  out.exitState.rigidityAtExit_GV=r.exit.rigidityAtExit_GV;
  out.exitState.valid=r.exit.valid!=0;
  out.driftTurnStates.clear();
  out.driftTurnStates.reserve(turns.size());
  for (const EgTurnState& s:turns) {
    Geospace::Trajectory::Result::DriftTurnState t;
    t.revolution=s.revolution;
    t.traceTime_s=s.traceTime_s;
    t.accumulatedAngle_rad=s.accumulatedAngle_rad;
    t.radius_m=s.radius_m;
    t.latitude=s.latitude;
    t.pitchCos2=s.pitchCos2;
    t.outerMargin_m=s.outerMargin_m;
    t.momentumRelativeSpread=s.momentumRelativeSpread;
    for (int d=0;d<3;++d) {
      t.x_m[d]=s.x_m[d]; t.p_unit[d]=s.p_unit[d]; t.b_unit[d]=s.b_unit[d];
    }
    out.driftTurnStates.push_back(t);
  }
}

bool VerifyEnabled();

// Bitwise comparison of an engine result with the legacy result.
void VerifyCompare(const char* where,double R_GV,const Geospace::Trajectory::Result& engine,
                   const Geospace::Trajectory::Result& legacy);

// Turn-state capacity overflow and other non-fatal engine notices (printed once).
void WarnOnce(const char* key,const std::string& message);

} // namespace GPU
} // namespace Geospace

#endif
