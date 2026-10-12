// ============================================================================
// AMPS binding of srcSEP3D point sampling ('#begin sampling: points').
//
// The sampler is registered with the AMPS user-defined sampling hook
//   PIC::Sampling::ExternalSamplingLocalVariables::RegisterSamplingRoutine
// (src/pic/pic.h), so it runs inside the AMPS sampling cycle:
//
//   * Sample() is called by PIC::Sampling::Sampling() on every rank at every
//     iteration while AMPS sampling is enabled (src/pic/pic.cpp, the
//     external-sampler loop).  This is at the start of PIC::TimeStep(), so it
//     sees the particles after the previous step's move and srcSEP3D
//     injection, at simulation time PIC::SimulationTime::Get().
//   * Print(N) is called on every rank when AMPS writes its output N, i.e.
//     when PIC::CollectingSampleCounter reached PIC::RequiredSampleLength (or
//     the wall-time alarm fired).  Its file therefore carries the same index
//     `.out=N` as the AMPS mesh files pic.<species>.s=<s>.out=N.dat and is
//     normalized by the same sampling length N = PIC::LastSampleLength.
//
// Interval bookkeeping.  AMPS resets PIC::CollectingSampleCounter to 0 before
// calling the printers, and does not call them at all when an output is
// suppressed (SupressOutputFlag, SkipOutputStep, minIterationNumberFor-
// DataOutput) although the counter still restarts.  Sample() therefore opens
// a new interval itself whenever it sees CollectingSampleCounter == 0; an
// interval that ended without a Print() still has its fluence added to the
// cumulative fluence (collective reduction at that moment), so no particle
// time is lost from the fluence.  Print() checks that the number of
// Sample() calls equals PIC::LastSampleLength and stops otherwise.
//
// Ownership and MPI.  Each point is bound once to its AMR tree node.  At each
// call only the rank that currently owns the node (node->Thread) samples it,
// because AMPS may rebalance blocks after every output.  Sums are kept per
// point, so moving a block between ranks inside an interval is harmless:
// Print() combines all ranks with MPI_Reduce(MPI_SUM) on rank 0, which
// normalizes, adds the interval fluence to the cumulative fluence and writes
// the Tecplot file (temporary file + rename).  The summation order across
// ranks can change the last bits of a value between rank counts (accepted by
// the user, 2026-10-11).
//
// Typed failure.  Start-up stops when a species symbol is not compiled or a
// point lies outside the allocated AMPS mesh or in an inactive region (user
// decision: fatal at start-up).  During the run, an unreadable particle, a
// non-positive cell volume or an invalid cell background stops the run with
// the point and the reason; nothing is skipped silently.
// ============================================================================

#ifndef SEP3D_AMPS_POINT_SAMPLER_H
#define SEP3D_AMPS_POINT_SAMPLER_H

#include "pic.h"

#include "../runtime/run_configuration.h"

#include <string>
#include <vector>

namespace SEP3D {
namespace AMPS {
namespace PointSampler {

// Background of the cell containing 'positionM' in the pinned snapshot:
// bulk velocity U [m/s] and unit magnetic direction b.  Supplied by
// main_lib.cpp (the stored-background reader used by the movers), so the
// inertial-frame transform sees exactly the background of the transport.
using CellBackgroundResolver = Core::Status (*)(
    const Core::Vec3& positionM, cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node,
    Core::Vec3* bulkVelocityMPerS, Core::Vec3* magneticDirection);

// Bind the configured species and points and register the AMPS sampling
// pair.  Call once on every rank after the mesh, the cell volumes, the
// compiled species table and the time steps exist and before the first
// PIC::TimeStep().  A disabled configuration registers nothing.  The output
// file is <PIC::OutputDataFileDirectory>/<filePrefix>.point-sampling.out=N.dat.
Core::Status Install(const RuntimeModel::PointSamplingOptions& options,
                     const std::vector<RuntimeModel::CompiledSpeciesRecord>& compiled,
                     const std::string& filePrefix,
                     CellBackgroundResolver resolver);

bool Installed();

// The callbacks registered with AMPS (exposed for diagnostics and tests).
void Sample();
void Print(int dataOutputFileNumber);

}  // namespace PointSampler
}  // namespace AMPS
}  // namespace SEP3D

#endif  // SEP3D_AMPS_POINT_SAMPLER_H
