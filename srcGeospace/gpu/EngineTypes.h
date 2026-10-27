#ifndef _SRC_GEOSPACE_GPU_ENGINE_TYPES_H_
#define _SRC_GEOSPACE_GPU_ENGINE_TYPES_H_

//======================================================================================
// EngineTypes.h -- plain-old-data types exchanged between the application drivers (host)
// and the batched trajectory engine (host or device).
//
// Nothing in this header owns memory or uses the C++ standard library containers; the
// types cross the boundary between the g++-compiled AMPS objects and an optional
// nvc++/OpenACC-compiled engine library.
//======================================================================================

#include "EngineConfig.h"

namespace Geospace {
namespace GPU {

// Termination codes: identical numbering to Geospace::GridlessMode::TrajectoryTermination.
enum EgTermination : int {
  kTermOuterBoundaryAllowed=0,
  kTermInnerBoundaryForbidden=1,
  kTermMagneticallyTrappedForbidden=2,
  kTermTimeLimit=3,
  kTermStepLimit=4,
  kTermDistanceLimit=5,
  kTermInvalidTimeStep=6,
  kTermInvalidField=7,
  kTermNumericalFailure=8,
  kTermDriftTrappedForbidden=9,
  kTermCancelled=100              // engine-internal: request cancelled (speculation)
};

// Mover codes: identical order to Geospace::Trajectory::Mover.
enum EgMover : int {
  kMoverBoris=0, kMoverHC4=1, kMoverRK2=2, kMoverRK4=3, kMoverRK6=4,
  kMoverGC2=5, kMoverGC4=6, kMoverGC6=7, kMoverHybrid=8
};

enum EgPolicy : int { kPolicyStructured=0, kPolicyLegacy=1 };

// The gridless and Mode3D legacy tracers differ in (a) the legacy time-step formula
// (Re-based vs metre-based distances), (b) field-availability checks and (c) the
// optional Shue outer boundary.  The variant selects the exact legacy semantics.
enum EgVariant : int { kVariantGridless=0, kVariantMode3D=1 };

enum EgFieldKind : int {
  kFieldNone=0,      // FIELD_MODEL=NONE (zero field)
  kFieldDipole=1,    // analytic centered dipole (device)
  kFieldMesh=2,      // Mode3D compact global array + AMPS row stencil (device)
  kFieldHost=3       // Geopack/Tsyganenko and other host-only evaluators
};

enum EgError : int {
  kErrNone=0,
  kErrMeshConsistency=1,     // legacy: fatal exit() in CheckedGlobalCellIndex_
  kErrRowStencilOverflow=2,  // AMPS: fatal exit() in cRowStencil::AddCell
  kErrNegativeWeight=3,      // AMPS: fatal exit() in cRowStencil::AddCell
  kErrDipoleField=4,         // legacy Mode3D: fatal exit() for invalid analytic dipole
  kErrHostField=5,           // exception thrown by a host field evaluator
  kErrTurnStateOverflow=6    // diagnostic drift-turn buffer capacity exceeded
};

struct EgTrapConfig {
  int enabled;
  int minMirrorPoints;
  int minBounceCycles;
  double outerMargin_m;
  double radialEnvelopeTolerance_m;
  double energyRelativeTolerance;
  double parallelDeadband;
  int driftEnabled;
  int minDriftRevolutions;
  double driftRadialAbsoluteTolerance_m;
  double driftRadialRelativeTolerance;
  double driftLatitudeTolerance;
  double driftPitchCos2Tolerance;
  double driftMaxMeanRadiusChange_m;
  int driftProfileBins;
  double driftMinProfileCoverage;
  double driftMinMatchedBinFraction;
  int driftWindowRevolutions;
  int driftMinStableWindowComparisons;
  int driftMaxRecurrenceLagRevolutions;
  int recordDriftTurnStates;
};

// Shue magnetopause outer boundary (Mode3D BOUNDARY_TYPE=SHUE), SI units.
struct EgShue {
  int enabled;
  double boxMin_m[3];
  double boxMax_m[3];
  double r0_Re;
  double alpha;
  double earthRadius_m;
  double tailCapX_m;
};

// Run-level description shared by many requests.
struct EgContext {
  int variant;
  int adaptiveDt;
  double boundaryEventTolerance_m;

  // Physical stopping surfaces in SI, exactly as built by the legacy tracer.
  double boxMin_m[3];
  double boxMax_m[3];
  double rInner_m;
  // Gridless legacy time-step selector works in Re (DomainBoxRe).
  double boxRe[7];          // xMin,xMax,yMin,yMax,zMin,zMax,rInner in Re
  // Mode3D: configured AMPS domain in metres (1000*prm.domain.*); used to separate an
  // expected out-of-box mesh sample from an AMR coverage hole.
  double fieldDomainMin_m[3];
  double fieldDomainMax_m[3];

  EgShue shue;
  EgTrapConfig trap;

  int fieldKind;
  // Dipole moment [A m^2] prepared exactly as the legacy evaluator of the variant
  // (gridless: M*m_hat of Dipole::MakeParams; Mode3D: moment*sin/cos(tilt)).
  double dipoleMx;
  double dipoleMy;
  double dipoleMz;
  int dipoleStrict;          // Mode3D analytic dipole: invalid sample is fatal
  int sampleDipoleError;     // Mode3D mesh DIPOLE run: accumulate |Bmesh-Bdip|/|Bdip|
};

struct EgRequest {
  double x0_m[3];
  double v0_unit[3];
  double R_GV;
  double qTrace_C;          // signed charge used by the mover (reversed for antiparticle)
  double qabs_C;            // |physical charge|
  double m0_kg;
  double maxTime_s;         // effective time budget of this attempt
  double maxDist_m;         // <=0: no path-length limit
  double dtTrace_s;         // DT_TRACE of this attempt (halved by numerical retry)
  int maxSteps;
  int mover;
  int policy;
  int captureExit;
  int context;              // index into the context table
  int hostField;            // host evaluator index (kFieldHost only)
  long long tag;            // caller cookie, returned unchanged
};

struct EgExitState {
  double x_exit_m[3];
  double p_exit_SI[3];
  double v_exit_unit[3];
  double cosAlpha;
  double traceTimeAtExit_s;
  double rigidityAtExit_GV;
  int valid;
};

struct EgResult {
  long long tag;
  int termination;
  int steps;
  double traceTime_s;
  double traceDistance_m;

  int mirrorPoints;
  int bounceCycles;
  int driftRevolutions;
  double driftAngle_rad;
  double driftMeanRadiusChange_m;
  double driftMeanRadiusSpan_m;
  int driftStableComparisons;
  int driftRecurrenceLagRevolutions;
  double minimumOuterMargin_m;
  double lastBounceEnvelopeChange_m;

  int driftComparisonAvailable;
  int driftPreviousCoveragePassed;
  int driftCurrentCoveragePassed;
  int driftSecularDriftPassed;
  int driftProfilesRecurred;
  double driftPreviousCoverageFraction;
  double driftCurrentCoverageFraction;
  int driftCommonBins;
  int driftMatchedBins;
  int driftRadialMismatchBins;
  int driftLatitudeMismatchBins;
  int driftPitchMismatchBins;
  double driftMatchedBinFraction;
  double driftMaximumNormalizedRadialResidual;
  double driftMaximumNormalizedLatitudeResidual;
  double driftMaximumNormalizedPitchResidual;
  int trapMechanism;
  double momentumRelativeSpread;

  EgExitState exit;

  int error;                // EgError
  int turnStateCount;

  // Diagnostics: which INVALID_FIELD check fired (0 none; 1 trap initialisation,
  // 2 Shue position invalid/unavailable, 3 exit population at a Shue start-of-step
  // escape, 4 exit population at a box start-of-step escape, 5 HYBRID branch
  // selection, 6 time-step selection, 7 mover step, 8 exit population at a Shue chord
  // escape, 9 Shue chord mesh-unavailable/invalid, 10 exit population at a box chord
  // escape, 11 trap-detector field sample) and the position at termination [m].
  int invalidReason;
  double finalX_m[3];
  int unavailableKind;          // first unavailable field sample: 1 no leaf, 2 empty row
  double unavailableX_m[3];

  // Mode3D DIPOLE mesh-interpolation error samples taken by this trajectory.
  unsigned long long dipoleErrorCount;
  double dipoleErrorSum;
  double dipoleErrorMax;
  double dipoleErrorMaxLocation_m[3];
};

// One recorded drift turn (AMPS_TRAP_TURN_DIAGNOSTICS=1), SI units.
struct EgTurnState {
  int revolution;
  double traceTime_s;
  double accumulatedAngle_rad;
  double x_m[3];
  double p_unit[3];
  double b_unit[3];
  double radius_m;
  double latitude;
  double pitchCos2;
  double outerMargin_m;
  double momentumRelativeSpread;
};

} // namespace GPU
} // namespace Geospace

#endif
