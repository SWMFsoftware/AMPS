// ============================================================================
// srcSEP3D/SEP3D.h
//
// Production umbrella header through Phases M/B/T/P/A/O. The AMPS boundary
// owns one typed Runtime whose configuration is supplied by the standalone or
// SWMF host.  Background snapshots and turbulence providers are installed as
// immutable/typed objects; process arguments and parameter files never leak
// into the physics layers. The generated R01 hook supplies the one mover
// declaration to pic_mover.cpp; numerical transport remains AMPS-independent.
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

#include "runtime/runtime.h"

#include <memory>

// Keep this AMPS-facing umbrella dependency-light.  pic.h includes SEP3D.h in
// every AMPS translation unit, including generic interface and mesh sources
// whose compiler command does not contain src/models/sep_common or SWCME
// include directories.  Concrete provider/source/restart headers therefore
// belong in main_lib.cpp or main.cpp, never in this transitive public header.
namespace SEP3D {
namespace Background { class BackgroundSnapshot; }
namespace Turbulence { class TurbulenceProvider; }
namespace Adapters { class ShockProvider; }
namespace Output { struct RestartState; }
}

namespace SEP3D {

// Return the process-owned lifecycle object used by both standalone and SWMF
// hosts.  All restart/output counters live inside this object; the production
// boundary maintains no independent cadence statics that could diverge after
// restart.  A host must install one immutable configuration before mesh setup.
RuntimeModel::Runtime& ApplicationRuntime();
Core::Status ConfigureApplication(
    const std::shared_ptr<const RuntimeModel::RunConfiguration3D>& configuration);

// Coupled hosts install initial imported data before amps_init(). At a later
// joined SnapshotReady boundary the same calls stage the next candidate; the
// R03 coordinator validates background and turbulence together and swaps them
// only after collective readiness. A standalone Parker run may omit both
// calls: amps_init() builds the corresponding immutable analytic providers.
Core::Status InstallBackgroundSnapshot(
    const std::shared_ptr<const Background::BackgroundSnapshot>& snapshot);
Core::Status InstallTurbulenceProvider(
    const std::shared_ptr<Turbulence::TurbulenceProvider>& provider);
Core::Status InstallShockProvider(
    const std::shared_ptr<Adapters::ShockProvider>& provider);
// Install a completely validated R07 candidate before mesh construction.
// The standalone driver calls RestoreRestartBeforeMesh first; this function
// retains the particle/provider/observer payload until amps_init() can restore
// AMPS ownership without changing any stochastic identity.
Core::Status InstallRestartState(const Output::RestartState& state);

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
