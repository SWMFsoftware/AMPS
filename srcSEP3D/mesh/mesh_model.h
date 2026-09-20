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

#include "../core/parker_geometry.h"
#include "../core/sep3d_types.h"
#include "../runtime/run_configuration.h"

#include <array>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <string>
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
  Core::Vec3 originM;
  RuntimeModel::InnerBoundaryMode innerBoundary =
      RuntimeModel::InnerBoundaryMode::Absorb;
  RuntimeModel::OuterBoundaryMode outerBoundary =
      RuntimeModel::OuterBoundaryMode::Escape;
};

struct ResolutionConfiguration {
  // All radius and Parker-geometry calculations are relative to this origin.
  // Keeping it in the AMPS-independent record makes translated-domain tests
  // possible even though the current analytic/SWMF background contract still
  // requires the production origin to be heliocentric zero.
  Core::Vec3 originM = {0.0, 0.0, 0.0};
  double innerRadiusM = 20.0 * Core::Const::R_sun;
  double outerRadiusM = Core::Const::AU;
  double minimumCellSizeM = 0.01 * Core::Const::AU;
  double backgroundCellSizeM = 0.25 * Core::Const::AU;

  // The surface target and transition radius define an explicit degradation
  // region.  The selected named profile is monotone on [inner,transition] and
  // equals the global size outside it.
  bool enableRadialRefinement = true;
  double solarSurfaceCellSizeM = 0.01 * Core::Const::AU;
  double solarRefinementOuterRadiusM = 0.25 * Core::Const::AU;
  RuntimeModel::RefinementProfile solarRefinementProfile =
      RuntimeModel::RefinementProfile::Smoothstep;
  double solarRefinementExponent = 1.0;

  // Optional Parker-spiral tube.  tubeLongitudeRad is the centreline
  // longitude at innerRadiusM; tubeColatitudeRad is measured from +Z.
  bool enableTubeRefinement = false;
  double tubeLongitudeRad = 0.0;
  double tubeColatitudeRad = 0.5 * Core::Const::kPi;
  double tubeReferenceRadiusM = Core::Const::AU;
  double tubeRadiusAtReferenceM = 0.03 * Core::Const::AU;
  RuntimeModel::TubeRadiusMode tubeRadiusMode =
      RuntimeModel::TubeRadiusMode::ConstantAngularWidth;
  double tubeCellSizeM = 0.01 * Core::Const::AU;
  RuntimeModel::RefinementProfile tubeTransverseProfile =
      RuntimeModel::RefinementProfile::Smoothstep;
  double tubeTransverseExponent = 1.0;
  double solarWindSpeedMPerS = Core::Const::V_sw_default;
  double solarRotationRateRadPerS = Core::Const::Omega_sun;
  Core::Vec3 rotationAxis = {0.0, 0.0, 1.0};

  // Explicit finite line requested by schema version 2.  The fast AMR law
  // evaluates the corresponding analytic curve, so changing pointCount does
  // not change the physical mesh.  BuildParkerCenterline materializes these
  // samples for initialization output and verification.
  Core::Vec3 parkerInitialPointM = {20.0 * Core::Const::R_sun, 0.0, 0.0};
  double parkerLengthM = Core::Const::AU;
  std::uint64_t parkerPointCount = 4001;

  unsigned cellsPerBlockEdge = 4;
  unsigned maximumLevel = 5;
  std::size_t blockOverheadBytes = 1024;
  std::size_t memoryBudgetBytes = std::size_t{4} * 1024 * 1024 * 1024;
  RuntimeModel::MemoryModelOptions memoryModel;
};

Core::Status Validate(const ResolutionConfiguration& configuration);
DomainBounds MakeDomain(
    const RuntimeModel::RunConfiguration3DOptions& configuration);

// Return the physical radius of the refined tube at a given heliocentric
// radius.  ConstantAngularWidth scales linearly with radius; PhysicalConstant
// retains the reference cross section everywhere.
double TubeRadiusM(double radiusM,
                   const ResolutionConfiguration& configuration);

// Analytic tube geometry used by both refinement and tests.  Simultaneously
// rotating a point and tubeLongitudeRad leaves TubeDistanceM invariant.
Core::Vec3 ParkerTubeDirection(double radiusM,
                               const ResolutionConfiguration& configuration);
Core::Vec3 ParkerTubeTangent(double radiusM,
                             const ResolutionConfiguration& configuration);
double TubeDistanceM(const Core::Vec3& positionM,
                     const ResolutionConfiguration& configuration);
double RequestedCellSizeM(const Core::Vec3& positionM,
                          const ResolutionConfiguration& configuration);

// Generate pointCount points separated by uniform requested arc length.  A
// midpoint tangent step is used instead of a first-order Euler step so a
// coarse diagnostic line remains faithful to the analytic Parker curve.
Core::Status BuildParkerCenterline(
    const ResolutionConfiguration& configuration,
    std::vector<Core::Vec3>* points);

// Serialize the exact finite centreline used by initialization.  Geometry is
// completed before the destination is opened, so a bad configuration cannot
// leave a plausible-looking partial Tecplot product.
Core::Status WriteParkerCenterlineTecplot(
    const ResolutionConfiguration& configuration, const std::string& path);

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

struct MemoryEstimate {
  std::size_t baseCellBytes = 0;
  std::size_t associatedDataBytes = 0;
  std::size_t samplingBytes = 0;
  std::size_t nodeBytes = 0;
  std::size_t blockBytes = 0;
  std::size_t particleBytes = 0;
  std::size_t communicationAndHaloBytes = 0;
  std::size_t subtotalBytes = 0;
  std::size_t safetyMarginBytes = 0;
  std::size_t totalBytes = 0;
};

struct RefinementPreflight {
  double minimumRequestedCellM = 0.0;
  double maximumRequestedCellM = 0.0;
  Core::Vec3 minimumLocationM;
  Core::Vec3 maximumLocationM;
  double tubeRadiusAtReferenceM = 0.0;
  std::vector<std::uint64_t> estimatedBlocksByLevel;
  MemoryEstimate memory;
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
MemoryEstimate EstimateWholeRunMemory(
    const MeshSummary& topology,
    const ResolutionConfiguration& resolution,
    const RuntimeModel::StorageLayout& storage);
Core::Status BuildRefinementPreflight(
    const DomainBounds& domain,
    const ResolutionConfiguration& resolution,
    const RuntimeModel::StorageLayout& storage,
    RefinementPreflight* report);

// Classify one segment crossing.  Direction matters: moving from outside into
// the shell is not an escape, while crossing inward through the solar sphere
// is absorption.  ImportedCoverage uses the same geometric crossing but keeps
// its distinct status message at the host boundary.
Core::Status ClassifyBoundaryCrossing(const Core::Vec3& previousM,
                                      const Core::Vec3& currentM,
                                      const DomainBounds& domain);

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
