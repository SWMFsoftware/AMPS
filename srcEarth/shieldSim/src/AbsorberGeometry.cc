#include "AbsorberGeometry.hh"

#include "MaterialFactory.hh"

#include <G4Box.hh>
#include <G4LogicalVolume.hh>
#include <G4PVPlacement.hh>
#include <G4SystemOfUnits.hh>
#include <G4ThreeVector.hh>
#include <G4UserLimits.hh>

namespace shield_sim {

AbsorberBuildResult SlabAbsorberGeometry::Build(G4LogicalVolume* world,
                                                MaterialFactory& materials,
                                                double maxStepMm) const {
  const std::string stableId = "absorber/" + definition_.id + "/sensitive";
  auto* material = materials.Resolve(definition_.materialId, "absorber");
  auto* solid = new G4Box((stableId + "/solid").c_str(),
                          definition_.sizeXmm * mm / 2.0,
                          definition_.sizeYmm * mm / 2.0,
                          definition_.thicknessMm * mm / 2.0);
  auto* logical = new G4LogicalVolume(solid, material, (stableId + "/logical").c_str());
  if (maxStepMm > 0.0) logical->SetUserLimits(new G4UserLimits(maxStepMm * mm));
  const double centerZ = definition_.frontZmm + definition_.thicknessMm / 2.0;
  new G4PVPlacement(nullptr, G4ThreeVector(definition_.centerXmm * mm,
                                           definition_.centerYmm * mm, centerZ * mm),
                    logical, (stableId + "/physical").c_str(), world, false, 0, true);

  AbsorberBuildResult result;
  result.geometryHash = definition_.geometryHash;
  SensitiveRegion region;
  region.stableId = stableId;
  region.logicalVolume = logical;
  region.massKg = logical->GetMass(true, false) / kg;
  region.materialHash = definition_.materialHash;
  result.sensitiveRegions.push_back(region);
  return result;
}

} // namespace shield_sim
