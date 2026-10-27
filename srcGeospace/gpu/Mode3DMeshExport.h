#ifndef _SRC_GEOSPACE_GPU_MODE3D_MESH_EXPORT_H_
#define _SRC_GEOSPACE_GPU_MODE3D_MESH_EXPORT_H_

//======================================================================================
// Mode3DMeshExport.h -- publish the Mode3D AMR tree + compact magnetic field to the
// batched trajectory engine (kFieldMesh requests).
//
// WHAT IS EXPORTED
// ----------------
// The legacy Mode3D code evaluates B in cMode3DMeshFieldEval::GetB_T() with
//
//   node = PIC::Mesh::mesh->findTreeNode(x,lastNode)                  (pointer tree)
//   GlobalMagneticField::InterpolateMagneticField(x,node,B)
//     -> PIC::InterpolationRoutines::CellCentered::Linear::InitStencil(x,node,row)
//     -> B = sum_s w_s * GlobalMagneticField_[3*(Temp_ID*Nc + i + Nx*(j+Ny*k))]
//
// The engine runs the device port of exactly these routines (gpu/EngineMesh.h) on a
// flattened copy of the tree.  ExportAndUpload() builds that copy (EgMeshHost):
//
//   * every tree node reachable from rootTree through downNode[] (depth-first, child
//     order 0..7, root = index 0) with xmin/xmax, xMinGlobalIndex,
//     NodeGeometricSizeIndex, upNode/downNode[8] as indices (-1 for NULL),
//     RefinmentLevel, min/maxNeibRefinmentLevel, IsUsedInCalculationFlag and Temp_ID;
//   * mesh->xGlobalMin/xGlobalMax/dx_max_refinment and _BLOCK_CELLS_X/Y/Z_;
//   * GlobalMagneticField_ (B, not copied here: the pool copies it on upload),
//     GlobalCellPresence_ converted to 0/1 (exactly the `!=1` test of
//     CheckedGlobalCellIndex_), GlobalUsedLeafBlocks_, GlobalInteriorCellCount_.
//
// All exported values are copied bit for bit, so the engine search/stencil/weights see
// the same doubles as the AMPS routines.
//
// WHEN IT IS UPLOADED
// -------------------
// BatchTracer::SetMesh() may only run while no batched trajectory is in flight.  The
// Mode3D drivers therefore call EnsureUploaded() on the rank's main thread before any
// fiber starts (Mode3D::Run after the field batch is frozen, and at the start of
// RunCutoffRigidity).  The resident copy is keyed by the GlobalMagneticField snapshot
// generation (bumped by every Assemble/Import/Redefine/Clear) plus a checksum of the
// exported tree topology/Temp_ID values, so a re-published field or a modified tree is
// always re-uploaded.  Per-trajectory code checks IsCurrent() (one atomic load); a stale
// mesh inside a fiber is a fatal error (it would mean a field publication while
// trajectories are in flight, which the application forbids through BeginFrozenFieldBatch).
//
// SELF-CHECK
// ----------
// With AMPS_EARTH_ENGINE_VERIFY=1 (or AMPS_EARTH_MESH_VERIFY=1) every upload is
// followed by a bitwise comparison of EgMesh::Find against findTreeNode (with and
// without a start-node hint) and of EgMesh::Interpolate against
// GlobalMagneticField::InterpolateMagneticField at AMPS_EARTH_MESH_VERIFY_SAMPLES
// (default 60000) deterministic sample points: uniform points over the global domain
// and points placed at/near block faces, edges, corners and cell-centre lines of
// randomly chosen leaves (these exercise the trilinear, multiblock and blended fine/
// coarse stencil branches).  The summary line is
//
//   [earth-engine][mesh-verify] generation=G nodes=N samples=S find-mismatch=a
//        hint-mismatch=b interp-mismatch=c
//
// printed by world rank 0 and by any rank that found a mismatch.
//======================================================================================

namespace Geospace {
namespace GPU {
namespace Mode3DMesh {

// Upload the currently published compact generation if the engine does not already
// hold it.  Main thread, outside fibers, no batched trajectory in flight.  Returns true
// when the current generation is resident after the call, false when no compact field
// generation is published (or the engine is disabled: AMPS_EARTH_ENGINE=LEGACY).
// Thread-safe (serialized by a mutex) for the lazy synchronous path.
bool EnsureUploaded(const char* label);

// True when the resident engine mesh belongs to the currently published generation.
// One atomic load plus one integer read; callable from any thread/fiber.
bool IsCurrent();

// Generation held by the engine (0 = none).
unsigned long long ResidentGeneration();

} // namespace Mode3DMesh
} // namespace GPU
} // namespace Geospace

#endif
