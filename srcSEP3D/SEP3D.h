// ============================================================================
// srcSEP3D/SEP3D.h
//
// Production umbrella header through Phase R2.  The AMPS boundary owns one
// typed Runtime whose configuration is supplied by the standalone or SWMF
// host.  It also includes the adapter whose compile-time assertions protect
// the mover return-code ABI. Prototype mover and sampler declarations remain
// absent until their later implementation phases.
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
#include "runtime/runtime.h"
#include "amps/amps_mover_status.h"

namespace SEP3D {

// Return the process-owned lifecycle object used by both standalone and SWMF
// hosts.  All restart/output counters live inside this object; the production
// boundary maintains no independent cadence statics that could diverge after
// restart.  A host must install one immutable configuration before mesh setup.
RuntimeModel::Runtime& ApplicationRuntime();
Core::Status ConfigureApplication(
    const std::shared_ptr<const RuntimeModel::RunConfiguration3D>& configuration);

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
