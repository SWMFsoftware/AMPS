#include "AbsorberGeometry.hh"

#include "MaterialFactory.hh"

#include <G4Box.hh>
#include <G4AffineTransform.hh>
#include <G4IntersectionSolid.hh>
#include <G4LogicalVolume.hh>
#include <G4PVPlacement.hh>
#include <G4RotationMatrix.hh>
#include <G4Sphere.hh>
#include <G4SubtractionSolid.hh>
#include <G4SystemOfUnits.hh>
#include <G4StepPoint.hh>
#include <G4TessellatedSolid.hh>
#include <G4ThreeVector.hh>
#include <G4TriangularFacet.hh>
#include <G4TouchableHistory.hh>
#include <G4Transform3D.hh>
#include <G4Tubs.hh>
#include <G4UnionSolid.hh>
#include <G4UserLimits.hh>
#include <G4VSolid.hh>

#include <algorithm>
#include <cmath>
#include <limits>
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

std::string CsgStableRegionId(const SlabAbsorberDefinition& definition,
                              const std::string& regionId){
  return "absorber/"+definition.id+"/region/"+regionId;
}

G4ThreeVector ToG4(const Vector3d& value){
  return {value.x*mm,value.y*mm,value.z*mm};
}

G4Transform3D ToG4Transform(const RigidTransform& transform){
  return {GeantObjectRotation(transform.rotation),ToG4(transform.translationMm)};
}

// Build a CSG tree recursively. Geant4's solid store owns every returned
// object. Node IDs are globally unique within an absorber, so diagnostic names
// remain stable even when two regions use the same primitive type.
G4VSolid* BuildCsgSolid(const AnalyticCsgNodeDefinition& node,
                       const std::string& namePrefix){
  const std::string name=namePrefix+"/node/"+node.nodeId;
  if(node.type=="box")
    return new G4Box(name.c_str(),node.sizeXmm*mm/2.0,node.sizeYmm*mm/2.0,
                     node.sizeZmm*mm/2.0);
  if(node.type=="sphere")
    return new G4Sphere(name.c_str(),node.innerRadiusMm*mm,node.outerRadiusMm*mm,
                        0.0,CLHEP::twopi,0.0,CLHEP::pi);
  if(node.type=="cylinder")
    return new G4Tubs(name.c_str(),node.innerRadiusMm*mm,node.outerRadiusMm*mm,
                      node.lengthMm*mm/2.0,0.0,CLHEP::twopi);
  G4VSolid* left=BuildCsgSolid(*node.left,namePrefix);
  G4VSolid* right=BuildCsgSolid(*node.right,namePrefix);
  const G4Transform3D rightTransform=ToG4Transform(node.rightTransform);
  if(node.type=="union")
    return new G4UnionSolid(name.c_str(),left,right,rightTransform);
  if(node.type=="intersection")
    return new G4IntersectionSolid(name.c_str(),left,right,rightTransform);
  if(node.type=="subtraction")
    return new G4SubtractionSolid(name.c_str(),left,right,rightTransform);
  throw std::runtime_error("unsupported analytic CSG node type '"+node.type+"'");
}

double PrimitiveVolumeCm3(const AnalyticCsgNodeDefinition& node){
  constexpr double piValue=3.141592653589793238462643383279502884;
  if(node.type=="box") return node.sizeXmm*node.sizeYmm*node.sizeZmm/1000.0;
  if(node.type=="sphere")
    return (4.0*piValue/3.0)*(std::pow(node.outerRadiusMm,3)-
                              std::pow(node.innerRadiusMm,3))/1000.0;
  if(node.type=="cylinder")
    return piValue*(node.outerRadiusMm*node.outerRadiusMm-
                    node.innerRadiusMm*node.innerRadiusMm)*node.lengthMm/1000.0;
  return std::numeric_limits<double>::quiet_NaN();
}

double PrimitiveSurfaceAreaCm2(const AnalyticCsgNodeDefinition& node){
  constexpr double piValue=3.141592653589793238462643383279502884;
  double squareMm=0.0;
  if(node.type=="box")
    squareMm=2.0*(node.sizeXmm*node.sizeYmm+node.sizeXmm*node.sizeZmm+
                  node.sizeYmm*node.sizeZmm);
  else if(node.type=="sphere")
    squareMm=4.0*piValue*(node.outerRadiusMm*node.outerRadiusMm+
                          node.innerRadiusMm*node.innerRadiusMm);
  else if(node.type=="cylinder")
    squareMm=2.0*piValue*((node.outerRadiusMm+node.innerRadiusMm)*node.lengthMm+
                          node.outerRadiusMm*node.outerRadiusMm-
                          node.innerRadiusMm*node.innerRadiusMm);
  else return 0.0;
  return squareMm/100.0;
}

bool NearlyEqualRelative(double actual,double expected,double tolerance){
  return std::abs(actual-expected)<=tolerance*std::max(std::abs(expected),1.0e-30);
}

Vector3d Cross(const Vector3d& a,const Vector3d& b){
  return {a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x};
}

double Magnitude(const Vector3d& value){ return std::sqrt(Dot(value,value)); }

double TriangleAreaCm2(const Vector3d& a,const Vector3d& b,const Vector3d& c){
  return 0.5*Magnitude(Cross(b-a,c-a))/100.0;
}

// Geant4 navigation gives a point on the tessellated boundary with finite
// floating-point tolerance. Test both the supporting plane and barycentric
// coordinates so another coplanar triangle cannot masquerade as a tagged one.
bool PointOnTriangle(const Vector3d& point,const Vector3d& a,const Vector3d& b,
                     const Vector3d& c,double toleranceMm){
  const Vector3d ab=b-a,ac=c-a,ap=point-a;
  const Vector3d normal=Cross(ab,ac);
  const double normalMagnitude=Magnitude(normal);
  const double geometricTolerance=std::max(toleranceMm,1.0e-6);
  if(normalMagnitude<=0.0 ||
     std::abs(Dot(normal,ap))>geometricTolerance*normalMagnitude) return false;
  const double d00=Dot(ab,ab),d01=Dot(ab,ac),d11=Dot(ac,ac);
  const double d20=Dot(ap,ab),d21=Dot(ap,ac);
  const double denominator=d00*d11-d01*d01;
  if(denominator<=0.0) return false;
  const double v=(d11*d20-d01*d21)/denominator;
  const double w=(d00*d21-d01*d20)/denominator;
  const double u=1.0-v-w;
  const double scale=std::max({Magnitude(ab),Magnitude(ac),Magnitude(c-b),1.0e-30});
  const double barycentricTolerance=geometricTolerance/scale;
  return u>=-barycentricTolerance && v>=-barycentricTolerance &&
         w>=-barycentricTolerance;
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
          definition_.transform,
          {-definition_.sizeXmm/2.0,-definition_.sizeYmm/2.0,0.0},
          { definition_.sizeXmm/2.0, definition_.sizeYmm/2.0,
            definition_.thicknessMm}};
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

const AnalyticCsgRegionDefinition&
AnalyticCsgAbsorberGeometry::SensitiveRegionDefinition(
    const std::string& stableId) const {
  for(const auto& region:definition_.csgRegions)
    if(region.role=="sensitive" &&
       CsgStableRegionId(definition_,region.id)==stableId) return region;
  throw std::invalid_argument("unknown analytic CSG sensitive-region ID '"+stableId+"'");
}

GeometryValidationReport AnalyticCsgAbsorberGeometry::Validate() const {
  GeometryValidationReport report;
  if(definition_.provider!="analytic_csg.v1" || definition_.providerVersion!="1.0.0")
    report.errors.push_back("analytic CSG provider identity/version is inconsistent");
  if(definition_.id.empty() || definition_.geometryHash.empty() ||
     definition_.geometryContentHash.empty())
    report.errors.push_back("absorber ID or canonical content hash is missing");
  std::size_t sensitive=0;
  for(const auto& region:definition_.csgRegions){
    if(region.role=="sensitive") ++sensitive;
    if(region.id.empty() || region.materialId.empty() || region.materialHash.empty() ||
       !region.solid)
      report.errors.push_back("region identity, material, or solid is missing");
    if(!std::isfinite(region.expectedVolumeCm3) || region.expectedVolumeCm3<=0.0 ||
       !std::isfinite(region.densityGPerCm3) || region.densityGPerCm3<=0.0)
      report.errors.push_back("region volume and density must be finite and positive");
    if(region.role=="sensitive"){
      const double analytic=PrimitiveVolumeCm3(*region.solid);
      if(!std::isfinite(analytic) ||
         !NearlyEqualRelative(analytic,region.expectedVolumeCm3,
                              region.massRelativeTolerance))
        report.errors.push_back("sensitive primitive volume disagrees with declared mass properties");
      const double area=PrimitiveSurfaceAreaCm2(*region.solid);
      if(!std::isfinite(area) || area<=0.0)
        report.errors.push_back("sensitive primitive has no finite exterior entry area");
      const double centroidNorm=std::sqrt(
          region.expectedCentroidMm.x*region.expectedCentroidMm.x+
          region.expectedCentroidMm.y*region.expectedCentroidMm.y+
          region.expectedCentroidMm.z*region.expectedCentroidMm.z);
      const double scale=std::max({region.solid->boundsMaximumMm.x-
                                      region.solid->boundsMinimumMm.x,
                                  region.solid->boundsMaximumMm.y-
                                      region.solid->boundsMinimumMm.y,
                                  region.solid->boundsMaximumMm.z-
                                      region.solid->boundsMinimumMm.z,1.0});
      if(centroidNorm>region.massRelativeTolerance*scale)
        report.errors.push_back("sensitive primitive centroid disagrees with its analytic origin");
    }
    try { (void)region.transform.rotation.Canonical(); }
    catch(const std::exception& error){ report.errors.push_back(error.what()); }
  }
  if(sensitive!=1)
    report.errors.push_back("analytic_csg.v1 requires exactly one sensitive region");
  try { (void)definition_.transform.rotation.Canonical(); }
  catch(const std::exception& error){ report.errors.push_back(error.what()); }
  report.valid=report.errors.empty();
  return report;
}

std::vector<std::string> AnalyticCsgAbsorberGeometry::SensitiveRegionIds() const {
  std::vector<std::string> result;
  for(const auto& region:definition_.csgRegions)
    if(region.role=="sensitive")
      result.push_back(CsgStableRegionId(definition_,region.id));
  return result;
}

RegionMassProperties AnalyticCsgAbsorberGeometry::MassProperties(
    const std::string& regionId) const {
  const auto& region=SensitiveRegionDefinition(regionId);
  const RigidTransform world=ComposeTransforms(definition_.transform,region.transform);
  return {regionId,region.expectedVolumeCm3,
          region.expectedVolumeCm3*region.densityGPerCm3/1000.0,
          world.TransformPoint(region.expectedCentroidMm)};
}

GeometryBounds AnalyticCsgAbsorberGeometry::BoundsAndReferenceFrame() const {
  return {{definition_.sizeXmm,definition_.sizeYmm,definition_.thicknessMm},
          definition_.transform,definition_.localBoundsMinimumMm,
          definition_.localBoundsMaximumMm};
}

bool AnalyticCsgAbsorberGeometry::IsEntryBoundary(
    const std::string& regionId,const G4StepPoint& postStepPoint) const {
  return EntryFluenceWeightPerCm2(regionId,postStepPoint)>0.0;
}

double AnalyticCsgAbsorberGeometry::EntryNormalizationAreaCm2(
    const std::string& regionId) const {
  const auto& region=SensitiveRegionDefinition(regionId);
  if(region.entrySurface=="local_negative_z_face"){
    if(region.solid->type=="box")
      return region.solid->sizeXmm*region.solid->sizeYmm/100.0;
    constexpr double piValue=3.141592653589793238462643383279502884;
    return piValue*(region.solid->outerRadiusMm*region.solid->outerRadiusMm-
                    region.solid->innerRadiusMm*region.solid->innerRadiusMm)/100.0;
  }
  if(sensitiveSurfaceAreaCm2_>0.0) return sensitiveSurfaceAreaCm2_;
  return PrimitiveSurfaceAreaCm2(*region.solid);
}

double AnalyticCsgAbsorberGeometry::EntryFluenceWeightPerCm2(
    const std::string& regionId,const G4StepPoint& postStepPoint) const {
  const auto& region=SensitiveRegionDefinition(regionId);
  if(!sensitiveSolid_) return 0.0;
  const auto touchable=postStepPoint.GetTouchableHandle();
  const auto* history=touchable->GetHistory();
  if(!history) return 0.0;
  const G4AffineTransform& globalToLocal=history->GetTopTransform();
  const G4ThreeVector localPosition=globalToLocal.TransformPoint(postStepPoint.GetPosition());
  const G4ThreeVector localDirection=
      globalToLocal.TransformAxis(postStepPoint.GetMomentumDirection()).unit();
  if(region.entrySurface=="local_negative_z_face"){
    const double halfLength=region.solid->type=="box"
        ? region.solid->sizeZmm*mm/2.0 : region.solid->lengthMm*mm/2.0;
    if(localDirection.z()<=0.0 ||
       std::abs(localPosition.z()+halfLength)>1.0*um) return 0.0;
    return 1.0/(EntryNormalizationAreaCm2(regionId)*localDirection.z());
  }
  if(sensitiveSolid_->Inside(localPosition)!=kSurface) return 0.0;
  const G4ThreeVector outward=sensitiveSolid_->SurfaceNormal(localPosition).unit();
  const double inwardCosine=-outward.dot(localDirection);
  if(inwardCosine<=1.0e-12) return 0.0;
  return 1.0/(EntryNormalizationAreaCm2(regionId)*inwardCosine);
}

Json AnalyticCsgAbsorberGeometry::CanonicalGeometryDescription() const {
  return definition_.canonicalGeometry;
}

AbsorberBuildResult AnalyticCsgAbsorberGeometry::Build(
    G4LogicalVolume* world,MaterialFactory& materials,double maxStepMm) const {
  const auto validation=Validate();
  if(!validation.valid)
    throw std::runtime_error("invalid analytic CSG absorber geometry: "+
                             validation.errors.front());
  sensitiveSolid_=nullptr;
  sensitiveSurfaceAreaCm2_=0.0;
  AbsorberBuildResult result;
  result.provider=definition_.provider;
  result.geometryHash=definition_.geometryHash;
  std::vector<G4PVPlacement*> placements;
  for(std::size_t i=0;i<definition_.csgRegions.size();++i){
    const auto& region=definition_.csgRegions[i];
    const std::string stableId=CsgStableRegionId(definition_,region.id);
    G4VSolid* solid=BuildCsgSolid(*region.solid,stableId);
    auto* material=materials.Resolve(region.materialId,"absorber");
    auto* logical=new G4LogicalVolume(solid,material,(stableId+"/logical").c_str());
    if(maxStepMm>0.0) logical->SetUserLimits(new G4UserLimits(maxStepMm*mm));
    const RigidTransform placement=ComposeTransforms(definition_.transform,region.transform);
    auto* physical=new G4PVPlacement(ToG4Transform(placement),logical,
        (stableId+"/physical").c_str(),world,false,static_cast<G4int>(i),false);
    placements.push_back(physical);

    const double actualVolume=solid->GetCubicVolume()/cm3;
    if(!std::isfinite(actualVolume) || !NearlyEqualRelative(
         actualVolume,region.expectedVolumeCm3,region.massRelativeTolerance))
      throw std::runtime_error("analytic CSG region '"+region.id+
          "' Geant4 volume disagrees with declared mass_properties.volume");
    if(region.role=="sensitive"){
      sensitiveSolid_=solid;
      sensitiveSurfaceAreaCm2_=solid->GetSurfaceArea()/cm2;
      const double expectedMassKg=region.expectedVolumeCm3*region.densityGPerCm3/1000.0;
      const double actualMassKg=logical->GetMass(true,false)/kg;
      if(!std::isfinite(actualMassKg) || !NearlyEqualRelative(
           actualMassKg,expectedMassKg,region.massRelativeTolerance))
        throw std::runtime_error("analytic CSG sensitive region mass disagrees with "
                                 "declared volume and resolved density");
      result.sensitiveRegions.push_back(
          {stableId,logical,actualMassKg,region.materialHash});
    }
  }
  // A fixed, explicit overlap sampling count is a validation gate, not a
  // geometry-repair operation. Any detected sibling/world overlap is fatal.
  for(auto* placement:placements)
    if(placement->CheckOverlaps(2000,0.0,false,1))
      throw std::runtime_error("analytic CSG geometry contains an overlap for volume '"+
                               std::string(placement->GetName())+"'");
  return result;
}

const TessellatedMeshRegionDefinition&
TessellatedMeshAbsorberGeometry::SensitiveRegionDefinition(
    const std::string& stableId) const {
  for(const auto& region:definition_.meshRegions)
    if(region.role=="sensitive" &&
       CsgStableRegionId(definition_,region.id)==stableId) return region;
  throw std::invalid_argument("unknown tessellated-mesh sensitive-region ID '"+
                              stableId+"'");
}

GeometryValidationReport TessellatedMeshAbsorberGeometry::Validate() const {
  GeometryValidationReport report;
  if(definition_.provider!="tessellated_mesh.v1" || definition_.providerVersion!="1.0.0")
    report.errors.push_back("tessellated-mesh provider identity/version is inconsistent");
  if(definition_.id.empty() || definition_.geometryHash.empty() ||
     definition_.geometryContentHash.empty())
    report.errors.push_back("absorber ID or canonical content hash is missing");
  std::size_t sensitive=0;
  for(const auto& region:definition_.meshRegions){
    if(region.role=="sensitive") ++sensitive;
    if(region.id.empty() || region.materialId.empty() || region.materialHash.empty() ||
       region.verticesMm.size()<4 || region.facets.size()<4)
      report.errors.push_back("region identity, material, vertices, or facets are missing");
    if(!std::isfinite(region.exactVolumeCm3) || region.exactVolumeCm3<=0.0 ||
       !std::isfinite(region.exactSurfaceAreaCm2) || region.exactSurfaceAreaCm2<=0.0 ||
       !std::isfinite(region.densityGPerCm3) || region.densityGPerCm3<=0.0)
      report.errors.push_back("mesh volume, surface area, and density must be finite and positive");
    if(!NearlyEqualRelative(region.exactVolumeCm3,region.expectedVolumeCm3,
                            region.massRelativeTolerance))
      report.errors.push_back("exact mesh volume disagrees with declared mass properties");
    if(region.role=="sensitive" && region.entrySurfaceType!="entire_exterior" &&
       region.entrySurfaceType!="tagged_facets")
      report.errors.push_back("sensitive mesh entry-surface type is unsupported");
    if(region.role=="sensitive" && region.entrySurfaceType=="tagged_facets" &&
       region.entryFacetIds.empty())
      report.errors.push_back("tagged mesh entry surface is empty");
    try { (void)region.transform.rotation.Canonical(); }
    catch(const std::exception& error){ report.errors.push_back(error.what()); }
  }
  if(sensitive!=1)
    report.errors.push_back("tessellated_mesh.v1 requires exactly one sensitive region");
  try { (void)definition_.transform.rotation.Canonical(); }
  catch(const std::exception& error){ report.errors.push_back(error.what()); }
  report.valid=report.errors.empty();
  return report;
}

std::vector<std::string> TessellatedMeshAbsorberGeometry::SensitiveRegionIds() const {
  std::vector<std::string> result;
  for(const auto& region:definition_.meshRegions)
    if(region.role=="sensitive")
      result.push_back(CsgStableRegionId(definition_,region.id));
  return result;
}

RegionMassProperties TessellatedMeshAbsorberGeometry::MassProperties(
    const std::string& regionId) const {
  const auto& region=SensitiveRegionDefinition(regionId);
  const RigidTransform world=ComposeTransforms(definition_.transform,region.transform);
  return {regionId,region.exactVolumeCm3,
          region.exactVolumeCm3*region.densityGPerCm3/1000.0,
          world.TransformPoint(region.exactCentroidMm)};
}

GeometryBounds TessellatedMeshAbsorberGeometry::BoundsAndReferenceFrame() const {
  return {{definition_.sizeXmm,definition_.sizeYmm,definition_.thicknessMm},
          definition_.transform,definition_.localBoundsMinimumMm,
          definition_.localBoundsMaximumMm};
}

bool TessellatedMeshAbsorberGeometry::IsEntryBoundary(
    const std::string& regionId,const G4StepPoint& postStepPoint) const {
  return EntryFluenceWeightPerCm2(regionId,postStepPoint)>0.0;
}

double TessellatedMeshAbsorberGeometry::EntryNormalizationAreaCm2(
    const std::string& regionId) const {
  const auto& region=SensitiveRegionDefinition(regionId);
  if(region.entrySurfaceType=="entire_exterior") return region.exactSurfaceAreaCm2;
  double area=0.0;
  for(const auto& facet:region.facets)
    if(std::find(region.entryFacetIds.begin(),region.entryFacetIds.end(),facet.id)!=
       region.entryFacetIds.end()){
      area+=TriangleAreaCm2(region.verticesMm[facet.vertexIndices[0]],
                            region.verticesMm[facet.vertexIndices[1]],
                            region.verticesMm[facet.vertexIndices[2]]);
    }
  return area;
}

double TessellatedMeshAbsorberGeometry::EntryFluenceWeightPerCm2(
    const std::string& regionId,const G4StepPoint& postStepPoint) const {
  const auto& region=SensitiveRegionDefinition(regionId);
  if(!sensitiveSolid_) return 0.0;
  const auto touchable=postStepPoint.GetTouchableHandle();
  const auto* history=touchable->GetHistory();
  if(!history) return 0.0;
  const G4AffineTransform& globalToLocal=history->GetTopTransform();
  const G4ThreeVector g4Position=globalToLocal.TransformPoint(postStepPoint.GetPosition());
  const G4ThreeVector g4Direction=
      globalToLocal.TransformAxis(postStepPoint.GetMomentumDirection()).unit();
  const Vector3d position{g4Position.x()/mm,g4Position.y()/mm,g4Position.z()/mm};
  const Vector3d direction{g4Direction.x(),g4Direction.y(),g4Direction.z()};
  if(region.entrySurfaceType=="tagged_facets"){
    for(const auto& facet:region.facets){
      if(std::find(region.entryFacetIds.begin(),region.entryFacetIds.end(),facet.id)==
         region.entryFacetIds.end()) continue;
      const Vector3d& a=region.verticesMm[facet.vertexIndices[0]];
      const Vector3d& b=region.verticesMm[facet.vertexIndices[1]];
      const Vector3d& c=region.verticesMm[facet.vertexIndices[2]];
      if(!PointOnTriangle(position,a,b,c,region.validationToleranceMm)) continue;
      const Vector3d normal=Cross(b-a,c-a);
      const double inwardCosine=-Dot((1.0/Magnitude(normal))*normal,direction);
      if(inwardCosine<=1.0e-12) return 0.0;
      return 1.0/(EntryNormalizationAreaCm2(regionId)*inwardCosine);
    }
    return 0.0;
  }
  if(sensitiveSolid_->Inside(g4Position)!=kSurface) return 0.0;
  const G4ThreeVector outward=sensitiveSolid_->SurfaceNormal(g4Position).unit();
  const double inwardCosine=-outward.dot(g4Direction);
  if(inwardCosine<=1.0e-12) return 0.0;
  return 1.0/(EntryNormalizationAreaCm2(regionId)*inwardCosine);
}

Json TessellatedMeshAbsorberGeometry::CanonicalGeometryDescription() const {
  return definition_.canonicalGeometry;
}

AbsorberBuildResult TessellatedMeshAbsorberGeometry::Build(
    G4LogicalVolume* world,MaterialFactory& materials,double maxStepMm) const {
  const auto validation=Validate();
  if(!validation.valid)
    throw std::runtime_error("invalid tessellated-mesh absorber geometry: "+
                             validation.errors.front());
  sensitiveSolid_=nullptr;
  sensitiveEntryAreaCm2_=0.0;
  AbsorberBuildResult result;
  result.provider=definition_.provider;
  result.geometryHash=definition_.geometryHash;
  std::vector<G4PVPlacement*> placements;
  for(std::size_t i=0;i<definition_.meshRegions.size();++i){
    const auto& region=definition_.meshRegions[i];
    const std::string stableId=CsgStableRegionId(definition_,region.id);
    auto* solid=new G4TessellatedSolid((stableId+"/solid").c_str());
    for(const auto& facet:region.facets){
      const Vector3d& a=region.verticesMm[facet.vertexIndices[0]];
      const Vector3d& b=region.verticesMm[facet.vertexIndices[1]];
      const Vector3d& c=region.verticesMm[facet.vertexIndices[2]];
      if(!solid->AddFacet(new G4TriangularFacet(ToG4(a),ToG4(b),ToG4(c),ABSOLUTE)))
        throw std::runtime_error("Geant4 rejected mesh facet '"+facet.id+
                                 "' in region '"+region.id+"'");
    }
    solid->SetSolidClosed(true);
    auto* material=materials.Resolve(region.materialId,"absorber");
    auto* logical=new G4LogicalVolume(solid,material,(stableId+"/logical").c_str());
    if(maxStepMm>0.0) logical->SetUserLimits(new G4UserLimits(maxStepMm*mm));
    const RigidTransform placement=ComposeTransforms(definition_.transform,region.transform);
    auto* physical=new G4PVPlacement(ToG4Transform(placement),logical,
        (stableId+"/physical").c_str(),world,false,static_cast<G4int>(i),false);
    placements.push_back(physical);

    const double actualVolume=solid->GetCubicVolume()/cm3;
    if(!std::isfinite(actualVolume) ||
       !NearlyEqualRelative(actualVolume,region.exactVolumeCm3,
                            region.massRelativeTolerance))
      throw std::runtime_error("tessellated-mesh region '"+region.id+
          "' Geant4 volume disagrees with the exact oriented-facet volume");
    if(region.role=="sensitive"){
      sensitiveSolid_=solid;
      sensitiveEntryAreaCm2_=EntryNormalizationAreaCm2(stableId);
      const double expectedMassKg=region.exactVolumeCm3*region.densityGPerCm3/1000.0;
      const double actualMassKg=logical->GetMass(true,false)/kg;
      if(!std::isfinite(actualMassKg) ||
         !NearlyEqualRelative(actualMassKg,expectedMassKg,region.massRelativeTolerance))
        throw std::runtime_error("tessellated-mesh sensitive region mass disagrees with "
                                 "the exact mesh volume and resolved density");
      result.sensitiveRegions.push_back(
          {stableId,logical,actualMassKg,region.materialHash});
    }
  }
  for(auto* placement:placements)
    if(placement->CheckOverlaps(2000,0.0,false,1))
      throw std::runtime_error("tessellated-mesh geometry contains an overlap for volume '"+
                               std::string(placement->GetName())+"'");
  return result;
}

} // namespace shield_sim
