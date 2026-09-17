// ============================================================================
// srcSEP3D/SEP3D.h
//
// Production umbrella header through Phases M/B/T/P/A/O. The AMPS boundary
// owns one typed Runtime whose configuration is supplied by the standalone or
// SWMF host.  Background snapshots and turbulence providers are installed as
// immutable/typed objects; process arguments and parameter files never leak
// into the physics layers. The one production mover declaration is provided by
// amps/amps_particle_adapter.h; numerical transport remains AMPS-independent.
//
// LAYER: L3 (application).  L3 may include AMPS.  Lower layers under core/
// and background/ must remain independent of pic.h and mpi.h.
// ============================================================================

#ifndef SEP3D_H_
#define SEP3D_H_

#include "pic.h"
#include "Exosphere.h"
#include "constants.h"
#include "SpiceEmptyDefinitions.h"

#include "core/sep3d_types.h"
#include "background/bg_provider.h"
#include "background/background_snapshot.h"
#include "runtime/runtime.h"
#include "turbulence/turbulence_provider.h"
#include "amps/amps_mover_status.h"
#include "amps/amps_particle_adapter.h"

namespace SEP3D {

// Return the process-owned lifecycle object used by both standalone and SWMF
// hosts.  All restart/output counters live inside this object; the production
// boundary maintains no independent cadence statics that could diverge after
// restart.  A host must install one immutable configuration before mesh setup.
RuntimeModel::Runtime& ApplicationRuntime();
Core::Status ConfigureApplication(
    const std::shared_ptr<const RuntimeModel::RunConfiguration3D>& configuration);

// Coupled hosts install imported data before amps_init().  A standalone Parker
// run may omit both calls: amps_init() builds a frozen analytic snapshot and a
// prescribed Kolmogorov provider from the immutable RunConfiguration3D.
Core::Status InstallBackgroundSnapshot(
    const std::shared_ptr<const Background::BackgroundSnapshot>& snapshot);
Core::Status InstallTurbulenceProvider(
    const std::shared_ptr<Turbulence::TurbulenceProvider>& provider);

// AMPS calls this before its legacy parser.  srcSEP3D intentionally performs
// no argument or AMPS_PARAM.in parsing here: a standalone driver or the SWMF
// coupler resolves input and calls ConfigureApplication explicitly.
void Init_BeforeParser();

} // namespace SEP3D

// AMPS application hooks.  Their definitions live in main_lib.cpp, the only
// retained production application source besides main.cpp.
void amps_init_mesh();
void amps_init();
int amps_time_step();
double localResolution(double* x);
double InitLoadMeasure(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node);
bool TrajectoryTrackingCondition(double* x, double* v, int spec,
                                 void* particleData);

#endif // SEP3D_H_
