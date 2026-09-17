// ============================================================================
// srcSEP3D/SEP3D.h
//
// Production umbrella header after Phase R0 (production-tree rebaseline).
// It declares only the AMPS application entry points and includes the one
// adapter whose compile-time assertions protect the mover return-code ABI.
// Physics/runtime objects are added by later phases; prototype mover and
// sampler declarations are intentionally absent.
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
#include "amps/amps_mover_status.h"

namespace SEP3D {

// Phase R0 intentionally provides no runnable transport Runtime.  The
// application entry points exist so the production objects and archives can
// be compiled and linked by AMPS, but they stop with an explicit diagnostic if
// executed.  Phase R2 replaces that stop with the typed Runtime lifecycle.
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
