//======================================================================================
// test_trap_equivalence.cpp -- bitwise equivalence of the engine trap detector
// (gpu/EngineTrap.h, fixed-capacity arena) with the legacy
// detector (util/TrajectoryTrapDetector.h, std::vector state).
//
// Synthetic orbit samples (gyration + bounce + azimuthal drift with slow radial
// breathing and noise, plus escaping/irregular cases) are fed to both detectors for a
// grid of configurations: bounce-only, drift with window 1..3, recurrence lag 1..4,
// 8..48 profile bins, secular-drift veto on/off and drift-turn recording.  After every
// sample the trapped flag must agree, and at the end every diagnostic field and every
// recorded turn state must be bitwise identical.
//
// Build/run:  g++ -std=c++17 -O2 -ffp-contract=off test_trap_equivalence.cpp -o t && ./t
//======================================================================================

#include "util/TrajectoryTrapDetector.h"
#include "gpu/EngineTrap.h"

#include <cmath>
#include <cstdio>
#include <cstring>
#include <random>
#include <vector>

using namespace Geospace;

static bool Same(double a,double b) { return std::memcmp(&a,&b,sizeof(double))==0; }

static int gFailures=0;
static int gCovDrift3=0,gCovCompared=0,gCovBounceTrap=0,gCovDriftTrap=0,gCovTurns=0;
#define CHECK_D(name,a,b) do { if (!Same((a),(b))) { if (gFailures<20) std::printf("  mismatch %s: %.17g vs %.17g\n",name,(double)(a),(double)(b)); ++gFailures; } } while(0)
#define CHECK_I(name,a,b) do { if ((long long)(a)!=(long long)(b)) { if (gFailures<20) std::printf("  mismatch %s: %lld vs %lld\n",name,(long long)(a),(long long)(b)); ++gFailures; } } while(0)

struct Case {
  TrajectoryTrap::Config cfg;
  int scenario;
  unsigned seed;
};

static void Sample(int scenario,int k,std::mt19937& rng,double x[3],double p[3],double B[3]) {
  std::normal_distribution<double> n01(0.0,1.0);
  const double Re=6371.0e3;
  const double t=k*0.37;
  double L=6.6+0.05*std::sin(t*0.01);
  double driftRate=(scenario==2) ? 0.0003 : 0.002;   // scenario 2: slow drift
  if (scenario==3 || scenario==7) L+=0.0004*k;          // secular leak
  const double phi=driftRate*k+0.3*std::sin(t*3.1);   // gyro-like azimuth jitter
  const double lat=((scenario>=6) ? 0.01 : 0.2)*std::sin(t*0.21);  // bounce in latitude
  const double r=L*Re*std::cos(lat)*std::cos(lat)*(1.0+0.002*n01(rng));
  x[0]=r*std::cos(lat)*std::cos(phi);
  x[1]=r*std::cos(lat)*std::sin(phi);
  x[2]=r*std::sin(lat);
  const double pm=1.0e-19*(1.0+((scenario==4) ? 1e-3*n01(rng) : 1e-7*n01(rng)));
  // Scenarios >=6: equatorially mirroring particles (pitch cosine inside the parallel
  // deadband) so only the drift-recurrence branch can classify the orbit.
  const double pitch=(scenario>=6) ? 1.0e-8*std::cos(t*0.21) : std::cos(t*0.21+0.5);
  B[0]=-3.0*x[0]*x[2]/(r*r)*1e-7; B[1]=-3.0*x[1]*x[2]/(r*r)*1e-7;
  B[2]=(1.0-3.0*x[2]*x[2]/(r*r))*1e-7*(-1.0);
  const double bm=std::sqrt(B[0]*B[0]+B[1]*B[1]+B[2]*B[2]);
  const double e1[3]={-std::sin(phi),std::cos(phi),0.0};
  for (int d=0;d<3;++d) p[d]=pm*(pitch*B[d]/bm+std::sqrt(std::max(0.0,1.0-pitch*pitch))*e1[d]);
  if (scenario==5 && k>5000) { x[0]*=1.0+1e-4*(k-5000); }  // escaping drift
}

static void RunCase(const Case& c,int nSamples) {
  TrajectoryBoundary::Box box{{-60*6371.0e3,-30*6371.0e3,-30*6371.0e3},
                              {20*6371.0e3,30*6371.0e3,30*6371.0e3},6371.0e3};
  TrajectoryTrap::Detector ref(c.cfg,box);

  GPU::EgTrapConfig ec;
  std::memset(&ec,0,sizeof(ec));
  ec.enabled=c.cfg.enabled; ec.minMirrorPoints=c.cfg.minMirrorPoints;
  ec.minBounceCycles=c.cfg.minBounceCycles; ec.outerMargin_m=c.cfg.outerMargin_m;
  ec.radialEnvelopeTolerance_m=c.cfg.radialEnvelopeTolerance_m;
  ec.energyRelativeTolerance=c.cfg.energyRelativeTolerance;
  ec.parallelDeadband=c.cfg.parallelDeadband; ec.driftEnabled=c.cfg.driftEnabled;
  ec.minDriftRevolutions=c.cfg.minDriftRevolutions;
  ec.driftRadialAbsoluteTolerance_m=c.cfg.driftRadialAbsoluteTolerance_m;
  ec.driftRadialRelativeTolerance=c.cfg.driftRadialRelativeTolerance;
  ec.driftLatitudeTolerance=c.cfg.driftLatitudeTolerance;
  ec.driftPitchCos2Tolerance=c.cfg.driftPitchCos2Tolerance;
  ec.driftMaxMeanRadiusChange_m=c.cfg.driftMaxMeanRadiusChange_m;
  ec.driftProfileBins=c.cfg.driftProfileBins;
  ec.driftMinProfileCoverage=c.cfg.driftMinProfileCoverage;
  ec.driftMinMatchedBinFraction=c.cfg.driftMinMatchedBinFraction;
  ec.driftWindowRevolutions=c.cfg.driftWindowRevolutions;
  ec.driftMinStableWindowComparisons=c.cfg.driftMinStableWindowComparisons;
  ec.driftMaxRecurrenceLagRevolutions=c.cfg.driftMaxRecurrenceLagRevolutions;
  ec.recordDriftTurnStates=c.cfg.recordDriftTurnStates;
  const GPU::EgTrapLayout layout=GPU::MakeTrapLayout(ec,4096);
  std::vector<double> dArena(static_cast<size_t>(layout.doublesPerSlot)+1);
  std::vector<int> iArena(static_cast<size_t>(layout.intsPerSlot)+1);
  std::vector<GPU::EgTurnState> turns(static_cast<size_t>(layout.turnCapacity)+1);
  GPU::EgTrapState state;
  GPU::EgTrapDetector eng(ec,box.min,box.max,layout,state,
                          layout.doublesPerSlot>0 ? dArena.data() : nullptr,
                          layout.intsPerSlot>0 ? iArena.data() : nullptr,
                          layout.turnCapacity>0 ? turns.data() : nullptr);
  eng.Reset();

  std::mt19937 rng(c.seed);
  for (int k=0;k<nSamples;++k) {
    double x[3],p[3],B[3];
    Sample(c.scenario,k,rng,x,p,B);
    const bool a=ref.Update(x,p,B,k*0.01);
    const bool b=eng.Update(x,p,B,k*0.01);
    if (a!=b) { CHECK_I("trapped flag",a,b); return; }
    if (a) break;
  }
  if (ref.driftRevolutions()>=3) ++gCovDrift3;
  if (ref.driftComparisonDiagnostics().available) ++gCovCompared;
  if (ref.trapped() && ref.mechanism()==TrajectoryTrap::Mechanism::Bounce) ++gCovBounceTrap;
  if (ref.trapped() && ref.mechanism()==TrajectoryTrap::Mechanism::Drift) ++gCovDriftTrap;
  if (!ref.driftTurnStates().empty()) ++gCovTurns;
  CHECK_I("mechanism",static_cast<int>(ref.mechanism()),state.mechanism);
  CHECK_I("mirrorPoints",ref.mirrorPoints(),state.mirrorPoints);
  CHECK_I("bounceCycles",ref.bounceCycles(),state.bounceCycles);
  CHECK_I("driftRevolutions",ref.driftRevolutions(),state.driftRevolutions);
  CHECK_I("stableDriftComparisons",ref.stableDriftComparisons(),state.stableDriftComparisons);
  CHECK_I("lag",ref.driftRecurrenceLagRevolutions(),state.driftRecurrenceLagRevolutions);
  CHECK_D("driftAngle",ref.driftAngleRadians(),eng.driftAngleRadians());
  CHECK_D("meanRadiusChange",ref.driftMeanRadiusChangeMeters(),state.lastDriftMeanRadiusChange_m);
  CHECK_D("meanRadiusSpan",ref.driftMeanRadiusSpanMeters(),state.lastDriftMeanRadiusSpan_m);
  CHECK_D("minimumOuterMargin",ref.minimumOuterMarginMeters(),state.minimumOuterMargin);
  CHECK_D("lastBounceEnvelope",ref.lastBounceEnvelopeChangeMeters(),state.lastBounceEnvelopeChange_m);
  CHECK_D("momentumSpread",ref.momentumRelativeSpread(),eng.momentumRelativeSpread());
  const auto& d=ref.driftComparisonDiagnostics();
  const auto& e=state.lastDriftComparison;
  CHECK_I("diag.available",d.available,e.available);
  CHECK_I("diag.prevCov",d.previousCoveragePassed,e.previousCoveragePassed);
  CHECK_I("diag.curCov",d.currentCoveragePassed,e.currentCoveragePassed);
  CHECK_I("diag.secular",d.secularDriftPassed,e.secularDriftPassed);
  CHECK_I("diag.recurred",d.profilesRecurred,e.profilesRecurred);
  CHECK_D("diag.prevCovFrac",d.previousCoverageFraction,e.previousCoverageFraction);
  CHECK_D("diag.curCovFrac",d.currentCoverageFraction,e.currentCoverageFraction);
  CHECK_I("diag.common",d.commonBins,e.commonBins);
  CHECK_I("diag.matched",d.matchedBins,e.matchedBins);
  CHECK_I("diag.radialMis",d.radialMismatchBins,e.radialMismatchBins);
  CHECK_I("diag.latMis",d.latitudeMismatchBins,e.latitudeMismatchBins);
  CHECK_I("diag.pitchMis",d.pitchMismatchBins,e.pitchMismatchBins);
  CHECK_D("diag.matchedFrac",d.matchedBinFraction,e.matchedBinFraction);
  CHECK_D("diag.maxR",d.maximumNormalizedRadialResidual,e.maximumNormalizedRadialResidual);
  CHECK_D("diag.maxLat",d.maximumNormalizedLatitudeResidual,e.maximumNormalizedLatitudeResidual);
  CHECK_D("diag.maxPitch",d.maximumNormalizedPitchResidual,e.maximumNormalizedPitchResidual);
  const auto& rt=ref.driftTurnStates();
  CHECK_I("turnStates",static_cast<long long>(rt.size()),state.turnCount);
  for (size_t i=0;i<rt.size() && static_cast<int>(i)<state.turnCount;++i) {
    CHECK_I("turn.rev",rt[i].revolution,turns[i].revolution);
    CHECK_D("turn.t",rt[i].traceTime_s,turns[i].traceTime_s);
    CHECK_D("turn.angle",rt[i].accumulatedAngle_rad,turns[i].accumulatedAngle_rad);
    CHECK_D("turn.r",rt[i].radius_m,turns[i].radius_m);
    CHECK_D("turn.margin",rt[i].outerMargin_m,turns[i].outerMargin_m);
    for (int dd=0;dd<3;++dd) {
      CHECK_D("turn.x",rt[i].x_m[dd],turns[i].x_m[dd]);
      CHECK_D("turn.p",rt[i].p_unit[dd],turns[i].p_unit[dd]);
      CHECK_D("turn.b",rt[i].b_unit[dd],turns[i].b_unit[dd]);
    }
  }
}

int main() {
  int nCases=0;
  const int windows[]={1,2,3};
  const int lags[]={1,2,4};
  const int bins[]={8,24,48};
  for (int scenario=0;scenario<8;++scenario)
  for (int drift=0;drift<2;++drift)
  for (int w:windows) for (int lag:lags) for (int nb:bins)
  for (int veto=0;veto<2;++veto) {
    if (!drift && (w!=1 || lag!=1 || nb!=24 || veto)) continue;
    Case c;
    c.cfg.enabled=true;
    c.cfg.outerMargin_m=1.0*6371.0e3;
    c.cfg.radialEnvelopeTolerance_m=0.05*6371.0e3;
    c.cfg.driftEnabled=(drift!=0);
    c.cfg.minDriftRevolutions=3;
    c.cfg.driftRadialAbsoluteTolerance_m=0.2*6371.0e3;
    c.cfg.driftProfileBins=nb;
    c.cfg.driftWindowRevolutions=w;
    c.cfg.driftMaxRecurrenceLagRevolutions=lag;
    c.cfg.driftMaxMeanRadiusChange_m=veto ? 0.05*6371.0e3 : 0.0;
    c.cfg.recordDriftTurnStates=(nb==24);
    c.scenario=scenario;
    c.seed=1234u+static_cast<unsigned>(nCases);
    RunCase(c,60000);
    ++nCases;
  }
  std::printf("trap detector equivalence: %d configurations, %d mismatches\n",nCases,gFailures);
  std::printf("coverage: >=3 drift revolutions %d, drift comparisons %d, bounce-trapped %d, "
              "drift-trapped %d, with turn states %d\n",gCovDrift3,gCovCompared,
              gCovBounceTrap,gCovDriftTrap,gCovTurns);
  return gFailures==0 ? 0 : 1;
}
