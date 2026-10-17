#ifndef SHIELDSIM_DETECTOR_CONSTRUCTION_HH
#define SHIELDSIM_DETECTOR_CONSTRUCTION_HH

/* ============================================================================
 * DetectorConstruction.hh
 *
 * Geant4 geometry for a slab-shielding calculation.
 *
 * Versioned geometry layout:
 *
 *   explicit source plane -> ordered 1--5 layer shield -> one absorber slab
 *
 * Case-v1 volumes are axis-aligned; case-v2 source, stack, and absorber frames
 * carry explicit local-to-world rigid transforms. Current solids remain
 * G4Box-based planar layers/slabs with explicit dimensions. Shield materials
 * come from the immutable MaterialLibrary; absorber construction is delegated
 * through AbsorberGeometry. Main rebuilds this world separately for every
 * requested absorber. The legacy CLI branch is retained below it with its
 * historical single-shield/default-world behavior.
 * ========================================================================== */

#include "ShieldSimConfig.hh"

#include <G4VUserDetectorConstruction.hh>
#include <G4ThreeVector.hh>

#include <vector>
#include <memory>

class G4Box;
class G4LogicalVolume;
class G4VPhysicalVolume;
namespace shield_sim { class MaterialFactory; class AbsorberGeometry; }

class DetectorConstruction : public G4VUserDetectorConstruction {
public:
  explicit DetectorConstruction(const Options& opts);
  ~DetectorConstruction() override;

  // Sweep mode updates only the active shield thickness and asks Geant4 to
  // rebuild the geometry through G4RunManager::ReinitializeGeometry().
  void SetShieldThickness(G4double t);
  void SetActiveAbsorber(std::size_t index);

  const Options& GetOptions() const;
  G4LogicalVolume* GetShieldLV() const;
  const std::vector<G4LogicalVolume*>& GetShieldLVs() const;
  const std::vector<G4LogicalVolume*>& GetScoringLVs() const;
  G4Box* GetWorldBox() const;
  G4double GetSourcePlaneHalfX() const;
  G4double GetSourcePlaneHalfY() const;
  G4double GetSourcePlaneCenterX() const;
  G4double GetSourcePlaneCenterY() const;
  G4double GetSourcePlaneZ() const;
  // Case-v2 source planes may be rotated.  These accessors expose the exact
  // local-to-world frame used by geometry construction and source sampling.
  G4ThreeVector GetSourcePlaneOrigin() const;
  G4ThreeVector GetSourcePlaneAxisX() const;
  G4ThreeVector GetSourcePlaneAxisY() const;
  G4ThreeVector GetSourcePlaneNormal() const;

  G4VPhysicalVolume* Construct() override;

private:
  Options fOpts;
  G4Box* fWorldBox=nullptr;
  G4LogicalVolume* fShieldLV=nullptr;
  std::vector<G4LogicalVolume*> fShieldLVs;
  std::vector<G4LogicalVolume*> fScoringLVs;
  std::unique_ptr<shield_sim::MaterialFactory> fMaterialFactory;
  std::unique_ptr<shield_sim::AbsorberGeometry> fAbsorberGeometry;
};

#endif // SHIELDSIM_DETECTOR_CONSTRUCTION_HH
