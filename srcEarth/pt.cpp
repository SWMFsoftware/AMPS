//functions used fot particle tracking 

#include "pic.h"
#include "Earth.h"

namespace Earth {
namespace ParticleTracker {

// Runtime controls for particle-trajectory initialization.  They are deliberately
// stored in the Earth model layer rather than in AMPS compile-time macros so that
// an individual run can enable/disable trajectory output from the input file or
// CLI without rebuilding AMPS.
//
// Defaults preserve the historical behavior of this file: if the AMPS tracker is
// compiled in and no model-level configuration is applied, the first ~40000
// particles satisfying the legacy near-Earth condition can be tracked.
static bool RuntimeTrajectoryTrackingEnabled = true;
static int  MaxRuntimeTrajectories           = 40000;
static bool ForceTrackInjectedParticles      = false;
static int  nRuntimeTrajectoriesInitialized  = 0;

void ConfigureRuntimeTrajectoryTracking(bool enabled,
                                        int maxTrajectories,
                                        bool forceInjectedParticles) {
  RuntimeTrajectoryTrackingEnabled = enabled;
  MaxRuntimeTrajectories           = maxTrajectories;
  ForceTrackInjectedParticles      = forceInjectedParticles;

  // A new model run should start counting from zero.  The function is called
  // during 3d_forward initialization before any boundary particles are injected.
  nRuntimeTrajectoriesInitialized  = 0;
}

bool TrajectoryTrackingCondition(double *x,double *v,int spec,void *ParticleData) {
  PIC::ParticleTracker::cParticleData *ParticleTrajectoryRecord;

  ParticleTrajectoryRecord=(PIC::ParticleTracker::cParticleData*)(PIC::ParticleTracker::ParticleDataRecordOffset+((PIC::ParticleBuffer::byte*)ParticleData));

  // If the particle already carries a trajectory record, keep tracking it.
  // This preserves AMPS' expected semantics when a tracked particle is moved or
  // re-evaluated after initialization.
  if (ParticleTrajectoryRecord->TrajectoryTrackingFlag==true) {
    return true;
  }

  // Model-level gate.  This is the requested runtime flag layered on top of
  // _PIC_PARTICLE_TRACKER_MODE_: the compile-time AMPS flag must enable the
  // tracker, and this runtime flag must also allow initialization.
  if (RuntimeTrajectoryTrackingEnabled==false) return false;
  if (MaxRuntimeTrajectories<=0) return false;
  if (nRuntimeTrajectoriesInitialized>=MaxRuntimeTrajectories) return false;

  // Legacy condition: track only particles born inside 5 Re.  The 3d_forward
  // boundary source can override this with ForceTrackInjectedParticles because
  // injected particles start on the outer boundary, not near Earth.
  const bool legacyNearEarthCondition = (Vector3D::Length(x)<5*_EARTH__RADIUS_);

  if (ForceTrackInjectedParticles || legacyNearEarthCondition) {
    nRuntimeTrajectoriesInitialized++;
    return true;
  }

  return false;
}

} // namespace ParticleTracker
} // namespace Earth
