// ============================================================================
// Phase M AMPS-independent mesh and storage model
//
// The production AMPS adapter and the standalone verifier must make the same
// refinement decisions.  This file therefore owns the domain presets,
// resolution law, deterministic leaf ordering, memory estimate, and the byte
// buffer used by component tests.  No AMPS tree or MPI type appears here.
// ============================================================================

#ifndef SEP3D_MESH_MODEL_H
#define SEP3D_MESH_MODEL_H

#include "../core/sep3d_types.h"
#include "../runtime/run_configuration.h"

#include <array>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <vector>

namespace SEP3D {
namespace Mesh {

// The AMPS mesh is Cartesian, while the physical domain is a heliocentric
// shell.  The cube below exactly encloses the requested outer sphere; the
// inner sphere is retained as physical boundary metadata and is not silently
// removed from the Cartesian allocation.
struct DomainBounds {
  Core::Vec3 minimumM;
  Core::Vec3 maximumM;
  double innerRadiusM = 0.0;
  double outerRadiusM = 0.0;
};

struct ResolutionConfiguration {
  double innerRadiusM = 20.0 * Core::Const::R_sun;
  double outerRadiusM = Core::Const::AU;
  double minimumCellSizeM = 0.01 * Core::Const::AU;
  double backgroundCellSizeM = 0.25 * Core::Const::AU;

  // The radial law is h(r)=h_min*r/r_inner, clipped to [h_min,h_bg].  It
  // therefore doubles exactly when radius doubles and joins the background
  // value continuously at r_inner*h_bg/h_min.
  bool enableRadialRefinement = true;

  // Optional Parker-spiral tube.  tubeLongitudeRad is the centreline
  // longitude at innerRadiusM; tubeColatitudeRad is measured from +Z.
  bool enableTubeRefinement = false;
  double tubeLongitudeRad = 0.0;
  double tubeColatitudeRad = 0.5 * Core::Const::kPi;
  int tubePolarity = 1;
  double tubeCoreRadiusM = 0.01 * Core::Const::AU;
  double tubeShoulderRadiusM = 0.03 * Core::Const::AU;
  double tubeCellSizeM = 0.01 * Core::Const::AU;
  double solarWindSpeedMPerS = Core::Const::V_sw_default;
  double solarRotationRateRadPerS = Core::Const::Omega_sun;

  unsigned cellsPerBlockEdge = 4;
  unsigned maximumLevel = 5;
  std::size_t blockOverheadBytes = 1024;
  std::size_t memoryBudgetBytes = std::size_t{4} * 1024 * 1024 * 1024;
};

Core::Status Validate(const ResolutionConfiguration& configuration);
DomainBounds MakeDomain(RuntimeModel::DomainPreset preset,
                        double innerRadiusM, double requestedOuterRadiusM);

// Analytic tube geometry used by both refinement and tests.  Simultaneously
// rotating a point and tubeLongitudeRad leaves TubeDistanceM invariant.
Core::Vec3 ParkerTubeDirection(double radiusM,
                               const ResolutionConfiguration& configuration);
double TubeDistanceM(const Core::Vec3& positionM,
                     const ResolutionConfiguration& configuration);
double RequestedCellSizeM(const Core::Vec3& positionM,
                          const ResolutionConfiguration& configuration);

struct LeafBlock {
  Core::Vec3 minimumM;
  Core::Vec3 maximumM;
  unsigned level = 0;
  std::uint64_t path = 0;       // three-bit child labels from the root
  std::uint64_t globalLeaf = 0; // assigned after deterministic spatial sort
  std::uint64_t firstCell = 0;
  int ownerRank = 0;
};

struct MeshSummary {
  std::uint64_t leafCount = 0;
  std::uint64_t cellCount = 0;
  std::vector<std::uint64_t> leavesByLevel;
  std::size_t estimatedBytes = 0;
};

class StandaloneOctree final {
 public:
  Core::Status Build(const DomainBounds& domain,
                     const ResolutionConfiguration& resolution,
                     const RuntimeModel::StorageLayout& storage,
                     int ownerRanks);

  const std::vector<LeafBlock>& leaves() const { return leaves_; }
  const MeshSummary& summary() const { return summary_; }
  bool IsBalanced() const;

 private:
  std::vector<LeafBlock> leaves_;
  MeshSummary summary_;
};

// Public negative-control hook used by MSH3D05 and by an AMPS leaf-manifest
// audit.  It accepts a caller-supplied leaf list so a deliberately unbalanced
// pair can prove that the checker itself is live.
bool AreLeavesBalanced(const std::vector<LeafBlock>& leaves);

std::size_t EstimateMemoryBytes(const MeshSummary& topology,
                                const ResolutionConfiguration& resolution,
                                const RuntimeModel::StorageLayout& storage);

// A small standalone storage image used to prove two Phase-M invariants:
// offsets are frozen before allocation, and only the deterministic owner may
// fill a cell.  Production AMPS storage uses the same StorageLayout but owns
// its memory itself.
class CellStorage final {
 public:
  Core::Status Allocate(const StandaloneOctree& mesh,
                        const RuntimeModel::StorageLayout& layout);
  Core::Status Write(std::uint64_t globalCell, int callerRank,
                     std::size_t relativeOffset, const void* source,
                     std::size_t bytes);
  Core::Status Read(std::uint64_t globalCell, std::size_t relativeOffset,
                    void* destination, std::size_t bytes) const;
  int Owner(std::uint64_t globalCell) const;
  std::size_t bytes_per_cell() const { return bytesPerCell_; }

 private:
  std::vector<unsigned char> bytes_;
  std::vector<int> owners_;
  std::size_t bytesPerCell_ = 0;
};

// Least-squares gradients accept neighbours at mixed distances, which is the
// relevant coarse/fine AMR-boundary case.  The solve fails explicitly for a
// rank-deficient stencil rather than returning a plausible zero gradient.
Core::Status ReconstructScalarGradient(
    const Core::Vec3& center, double centerValue,
    const std::vector<Core::Vec3>& neighbourPositions,
    const std::vector<double>& neighbourValues,
    Core::Vec3* gradient);
Core::Status ReconstructVectorGradient(
    const Core::Vec3& center, const Core::Vec3& centerValue,
    const std::vector<Core::Vec3>& neighbourPositions,
    const std::vector<Core::Vec3>& neighbourValues,
    Core::Tensor3* gradient);

}  // namespace Mesh
}  // namespace SEP3D

#endif  // SEP3D_MESH_MODEL_H
