#ifndef SHIELDSIM_ABSORBER_GEOMETRY_HH
#define SHIELDSIM_ABSORBER_GEOMETRY_HH

/*
 * Geometry-neutral absorber contract.
 *
 * Scorers consume stable region IDs and logical-volume pointers returned here;
 * they do not depend on incidental copy numbers.  A future CAD/mesh/voxel
 * provider can implement this interface while retaining the same scoring and
 * provenance code.  See README.md for the validation obligations placed on a
 * future arbitrary-shape provider.
 */

#include "CaseDefinition.hh"

#include <string>
#include <vector>

class G4LogicalVolume;

namespace shield_sim {

class MaterialFactory;

struct SensitiveRegion {
  std::string stableId;
  G4LogicalVolume* logicalVolume = nullptr;
  double massKg = 0.0;
  std::string materialHash;
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
  AbsorberBuildResult Build(G4LogicalVolume* world,
                            MaterialFactory& materials,
                            double maxStepMm) const override;

private:
  SlabAbsorberDefinition definition_;
};

} // namespace shield_sim

#endif
