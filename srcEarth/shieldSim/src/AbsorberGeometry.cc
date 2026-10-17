#include "AbsorberGeometry.hh"

#include "MaterialFactory.hh"

#include <G4Box.hh>
#include <G4AffineTransform.hh>
#include <G4LogicalVolume.hh>
#include <G4PVPlacement.hh>
#include <G4RotationMatrix.hh>
#include <G4SystemOfUnits.hh>
#include <G4StepPoint.hh>
#include <G4ThreeVector.hh>
#include <G4TouchableHistory.hh>
#include <G4Transform3D.hh>
#include <G4UserLimits.hh>

#include <cmath>
#include <stdexcept>

namespace shield_sim {

namespace {
std::string StableRegionId(const SlabAbsorberDefinition& definition){
  return "absorber/"+definition.id+"/sensitive";
}

G4RotationMatrix GeantObjectRotation(const Quaterniond& input){
  // Use the direct-object G4Transform3D placement API.  The older pointer
  // overload names its matrix as a frame rotation and therefore uses the
  // inverse convention, which is too easy to misuse for input quaternions.
  const Quaterniond q=input.Canonical();
  G4RotationMatrix rotation;
  rotation.set(
      G4ThreeVector(1.0-2.0*(q.y*q.y+q.z*q.z),2.0*(q.x*q.y+q.w*q.z),2.0*(q.x*q.z-q.w*q.y)),
      G4ThreeVector(2.0*(q.x*q.y-q.w*q.z),1.0-2.0*(q.x*q.x+q.z*q.z),2.0*(q.y*q.z+q.w*q.x)),
      G4ThreeVector(2.0*(q.x*q.z+q.w*q.y),2.0*(q.y*q.z-q.w*q.x),1.0-2.0*(q.x*q.x+q.y*q.y)));
  return rotation;
}
}

GeometryValidationReport SlabAbsorberGeometry::Validate() const {
  GeometryValidationReport report;
  if(definition_.id.empty()) report.errors.push_back("stable absorber ID is empty");
  if(definition_.materialId.empty() || definition_.materialHash.empty())
    report.errors.push_back("resolved material identity/hash is missing");
  if(!std::isfinite(definition_.sizeXmm) || definition_.sizeXmm<=0.0 ||
     !std::isfinite(definition_.sizeYmm) || definition_.sizeYmm<=0.0 ||
     !std::isfinite(definition_.thicknessMm) || definition_.thicknessMm<=0.0)
    report.errors.push_back("slab dimensions must be finite and positive");
  if(!std::isfinite(definition_.densityGPerCm3) || definition_.densityGPerCm3<=0.0)
    report.errors.push_back("resolved density must be finite and positive");
  try { (void)definition_.transform.rotation.Canonical(); }
  catch(const std::exception& error){ report.errors.push_back(error.what()); }
  report.valid=report.errors.empty();
  return report;
}

std::vector<std::string> SlabAbsorberGeometry::SensitiveRegionIds() const {
  return {StableRegionId(definition_)};
}

RegionMassProperties SlabAbsorberGeometry::MassProperties(const std::string& regionId) const {
  if(regionId!=StableRegionId(definition_))
    throw std::invalid_argument("unknown slab sensitive-region ID '"+regionId+"'");
  const double volumeCm3=(definition_.sizeXmm/10.0)*(definition_.sizeYmm/10.0)*
                         (definition_.thicknessMm/10.0);
  return {regionId,volumeCm3,volumeCm3*definition_.densityGPerCm3/1000.0,
          definition_.transform.TransformPoint({0.0,0.0,definition_.thicknessMm/2.0})};
}

GeometryBounds SlabAbsorberGeometry::BoundsAndReferenceFrame() const {
  return {{definition_.sizeXmm,definition_.sizeYmm,definition_.thicknessMm},
          definition_.transform};
}

bool SlabAbsorberGeometry::IsEntryBoundary(
    const std::string& regionId,const G4StepPoint& postStepPoint) const {
  return EntryFluenceWeightPerCm2(regionId,postStepPoint)>0.0;
}

double SlabAbsorberGeometry::EntryFluenceWeightPerCm2(
    const std::string& regionId,const G4StepPoint& postStepPoint) const {
  if(regionId!=StableRegionId(definition_))
    throw std::invalid_argument("unknown slab sensitive-region ID '"+regionId+"'");
  const auto touchable=postStepPoint.GetTouchableHandle();
  const auto* history=touchable->GetHistory();
  if(!history) return 0.0;
  const G4AffineTransform& globalToLocal=history->GetTopTransform();
  const G4ThreeVector localPosition=
      globalToLocal.TransformPoint(postStepPoint.GetPosition());
  const G4ThreeVector localDirection=
      globalToLocal.TransformAxis(postStepPoint.GetMomentumDirection());
  // The slab reference surface is its upstream local-z face. Side entries and
  // downstream-face backscatter do not belong to this response surface. Re-entry
  // through the upstream face is a new fluence crossing and is accepted.
  if(localDirection.z()<=0.0 ||
     std::abs(localPosition.z()+definition_.thicknessMm*mm/2.0)>1.0*um)
    return 0.0;
  // Particle fluence is defined per area normal to the particle direction.
  // A track crossing a fixed planar patch therefore represents 1/(A*cosTheta)
  // rather than 1/A. localDirection is a unit vector and +z is the inward
  // surface normal for slab.v1.
  return 1.0/(EntryNormalizationAreaCm2(regionId)*localDirection.z());
}

double SlabAbsorberGeometry::EntryNormalizationAreaCm2(
    const std::string& regionId) const {
  if(regionId!=StableRegionId(definition_))
    throw std::invalid_argument("unknown slab sensitive-region ID '"+regionId+"'");
  return (definition_.sizeXmm/10.0)*(definition_.sizeYmm/10.0);
}

Json SlabAbsorberGeometry::CanonicalGeometryDescription() const {
  const Quaterniond q=definition_.transform.rotation.Canonical();
  return Json(Json::Object{
    {"geometry_hash",definition_.geometryHash},{"material_hash",definition_.materialHash},
    {"provider","slab.v1"},{"region_id",StableRegionId(definition_)},
    {"size_x_mm",definition_.sizeXmm},{"size_y_mm",definition_.sizeYmm},
    {"thickness_mm",definition_.thicknessMm},
    {"transform",Json(Json::Object{
      {"translation_mm",Json(Json::Array{definition_.transform.translationMm.x,
                                          definition_.transform.translationMm.y,
                                          definition_.transform.translationMm.z})},
      {"unit_quaternion_wxyz",Json(Json::Array{q.w,q.x,q.y,q.z})}})}
  });
}

AbsorberBuildResult SlabAbsorberGeometry::Build(G4LogicalVolume* world,
                                                MaterialFactory& materials,
                                                double maxStepMm) const {
  const auto validation=Validate();
  if(!validation.valid)
    throw std::runtime_error("invalid slab absorber geometry: "+validation.errors.front());
  const std::string stableId = StableRegionId(definition_);
  auto* material = materials.Resolve(definition_.materialId, "absorber");
  auto* solid = new G4Box((stableId + "/solid").c_str(),
                          definition_.sizeXmm * mm / 2.0,
                          definition_.sizeYmm * mm / 2.0,
                          definition_.thicknessMm * mm / 2.0);
  auto* logical = new G4LogicalVolume(solid, material, (stableId + "/logical").c_str());
  if (maxStepMm > 0.0) logical->SetUserLimits(new G4UserLimits(maxStepMm * mm));
  const Vector3d center=definition_.transform.TransformPoint(
      {0.0,0.0,definition_.thicknessMm/2.0});
  new G4PVPlacement(G4Transform3D(GeantObjectRotation(definition_.transform.rotation),
                                 G4ThreeVector(center.x * mm,center.y * mm,center.z * mm)),
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
