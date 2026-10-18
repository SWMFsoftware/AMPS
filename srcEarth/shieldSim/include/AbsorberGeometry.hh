#ifndef SHIELDSIM_ABSORBER_GEOMETRY_HH
#define SHIELDSIM_ABSORBER_GEOMETRY_HH

/*
 * Geometry-neutral absorber contract.
 *
 * Scorers consume stable region IDs and logical-volume pointers returned here;
 * they do not depend on incidental copy numbers. The bounded analytic-CSG and
 * convex tessellated-mesh providers implement this interface; future general
 * CAD or voxel providers can do the same without changing scoring/provenance
 * code. See README.md for each provider's explicit validation boundary.
 */

#include "CaseDefinition.hh"
#include "StrictJson.hh"

#include <string>
#include <vector>

class G4LogicalVolume;
class G4StepPoint;
class G4VSolid;

namespace shield_sim {

class MaterialFactory;

struct SensitiveRegion {
  std::string stableId;
  G4LogicalVolume* logicalVolume = nullptr;
  double massKg = 0.0;
  std::string materialHash;
};

struct GeometryValidationReport {
  bool valid = false;
  std::vector<std::string> errors;
};

struct RegionMassProperties {
  std::string stableId;
  double volumeCm3 = 0.0;
  double massKg = 0.0;
  Vector3d centroidWorldMm;
};

struct GeometryBounds {
  Vector3d localFullDimensionsMm;
  RigidTransform localToWorld;
  // Explicit extrema are required because a provider origin need not coincide
  // with the center of its conservative bound (slab.v1, for example, places
  // its origin on the upstream face).
  Vector3d localMinimumMm;
  Vector3d localMaximumMm;
};

struct AbsorberBuildResult {
  std::string provider = "slab.v1";
  std::string geometryHash;
  std::vector<SensitiveRegion> sensitiveRegions;
};

class AbsorberGeometry {
public:
  virtual ~AbsorberGeometry() = default;
  virtual const std::string& Id() const = 0;
  virtual const std::string& Hash() const = 0;
  virtual GeometryValidationReport Validate() const = 0;
  virtual std::vector<std::string> SensitiveRegionIds() const = 0;
  virtual RegionMassProperties MassProperties(const std::string& regionId) const = 0;
  virtual GeometryBounds BoundsAndReferenceFrame() const = 0;
  // Surface-folded responses need a provider-defined entry surface and its
  // normalization area. The transport/quantity layers deliberately do not
  // inspect a G4Box, tessellated solid, or voxel representation. Every shape
  // provider must define which exterior crossings constitute entry into each
  // sensitive region and the corresponding physical area in cm2.
  virtual bool IsEntryBoundary(const std::string& regionId,
                               const G4StepPoint& postStepPoint) const = 0;
  virtual double EntryNormalizationAreaCm2(const std::string& regionId) const = 0;
  // Return the scalar-fluence contribution of one accepted crossing in cm^-2.
  // For a planar patch this is 1/(A*|n.u|), not merely 1/A; the projected-
  // cosine correction is required for oblique tracks. Zero means the crossing
  // is not on the provider's declared entry surface.
  virtual double EntryFluenceWeightPerCm2(
      const std::string& regionId,const G4StepPoint& postStepPoint) const = 0;
  virtual Json CanonicalGeometryDescription() const = 0;
  virtual AbsorberBuildResult Build(G4LogicalVolume* world,
                                    MaterialFactory& materials,
                                    double maxStepMm) const = 0;
};

class SlabAbsorberGeometry final : public AbsorberGeometry {
public:
  explicit SlabAbsorberGeometry(SlabAbsorberDefinition definition)
    : definition_(std::move(definition)) {}
  const std::string& Id() const override { return definition_.id; }
  const std::string& Hash() const override { return definition_.geometryHash; }
  GeometryValidationReport Validate() const override;
  std::vector<std::string> SensitiveRegionIds() const override;
  RegionMassProperties MassProperties(const std::string& regionId) const override;
  GeometryBounds BoundsAndReferenceFrame() const override;
  bool IsEntryBoundary(const std::string& regionId,
                       const G4StepPoint& postStepPoint) const override;
  double EntryNormalizationAreaCm2(const std::string& regionId) const override;
  double EntryFluenceWeightPerCm2(
      const std::string& regionId,const G4StepPoint& postStepPoint) const override;
  Json CanonicalGeometryDescription() const override;
  AbsorberBuildResult Build(G4LogicalVolume* world,
                            MaterialFactory& materials,
                            double maxStepMm) const override;

private:
  SlabAbsorberDefinition definition_;
};

// Analytic CSG provider introduced by Stage 10. Regions are placed as sibling
// volumes in the world and therefore must be disjoint. Boolean nodes describe
// material, including shells, rather than relying on implicit mother-volume
// material replacement. This is easier to audit and makes overlap failures
// fatal. The v1 scoring contract admits exactly one sensitive region; passive
// regions can use different explicitly resolved materials.
class AnalyticCsgAbsorberGeometry final : public AbsorberGeometry {
public:
  explicit AnalyticCsgAbsorberGeometry(SlabAbsorberDefinition definition)
    : definition_(std::move(definition)) {}
  const std::string& Id() const override { return definition_.id; }
  const std::string& Hash() const override { return definition_.geometryHash; }
  GeometryValidationReport Validate() const override;
  std::vector<std::string> SensitiveRegionIds() const override;
  RegionMassProperties MassProperties(const std::string& regionId) const override;
  GeometryBounds BoundsAndReferenceFrame() const override;
  bool IsEntryBoundary(const std::string& regionId,
                       const G4StepPoint& postStepPoint) const override;
  double EntryNormalizationAreaCm2(const std::string& regionId) const override;
  double EntryFluenceWeightPerCm2(
      const std::string& regionId,const G4StepPoint& postStepPoint) const override;
  Json CanonicalGeometryDescription() const override;
  AbsorberBuildResult Build(G4LogicalVolume* world,
                            MaterialFactory& materials,
                            double maxStepMm) const override;

private:
  const AnalyticCsgRegionDefinition& SensitiveRegionDefinition(
      const std::string& stableId) const;
  SlabAbsorberDefinition definition_;
  // Populated during Build. Geant4 owns the solid; the provider retains a
  // non-owning pointer only to evaluate the declared entry surface after
  // navigation has placed a step at the boundary.
  mutable G4VSolid* sensitiveSolid_ = nullptr;
  mutable double sensitiveSurfaceAreaCm2_ = 0.0;
};

class TessellatedMeshAbsorberGeometry final : public AbsorberGeometry {
public:
  explicit TessellatedMeshAbsorberGeometry(SlabAbsorberDefinition definition)
    : definition_(std::move(definition)) {}
  const std::string& Id() const override { return definition_.id; }
  const std::string& Hash() const override { return definition_.geometryHash; }
  GeometryValidationReport Validate() const override;
  std::vector<std::string> SensitiveRegionIds() const override;
  RegionMassProperties MassProperties(const std::string& regionId) const override;
  GeometryBounds BoundsAndReferenceFrame() const override;
  bool IsEntryBoundary(const std::string& regionId,
                       const G4StepPoint& postStepPoint) const override;
  double EntryNormalizationAreaCm2(const std::string& regionId) const override;
  double EntryFluenceWeightPerCm2(
      const std::string& regionId,const G4StepPoint& postStepPoint) const override;
  Json CanonicalGeometryDescription() const override;
  AbsorberBuildResult Build(G4LogicalVolume* world,
                            MaterialFactory& materials,
                            double maxStepMm) const override;

private:
  const TessellatedMeshRegionDefinition& SensitiveRegionDefinition(
      const std::string& stableId) const;
  SlabAbsorberDefinition definition_;
  mutable G4VSolid* sensitiveSolid_ = nullptr;
  mutable double sensitiveEntryAreaCm2_ = 0.0;
};

} // namespace shield_sim

#endif
