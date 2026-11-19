#include "pic.h"

#include "MoonSurfaceEnvironment.h"

#include "Moon.h"
#include "MoonSurfaceTemperature.h"
#include "MoonTransportClock.h"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstring>
#include <exception>
#include <numeric>

#include <memory>

#ifndef _NO_SPICE_CALLS_
#include "SpiceUsr.h"
#endif

namespace Moon {
namespace SurfaceEnvironment {
namespace {

constexpr std::size_t kLeafFaces = 8;

std::unique_ptr<TerrainRayTracer> gProductionTerrain;
const Runtime::Configuration* gProductionConfiguration=nullptr;
double gParticleStepBeginEtS=0.0;
thread_local bool gCollisionFrameValid=false;
// P08 state is installed explicitly and remains immutable during transport.
// The ledger/resolver user data are application-owned and must outlive every
// impact callback. No collector is inferred from neutral surface parameters.
std::vector<Species::IonCollector::StateDefinition> gIonCollectorStates;
IonCollectorStateResolver gIonCollectorResolver=nullptr;
void* gIonCollectorUser=nullptr;
Species::IonCollector::Ledger* gIonCollectorLedger=nullptr;
double gIonCollectorElementaryChargeC=ElectronCharge;
std::string gIonCollectorChargeSourceIdentity;
thread_local double gCollisionNormalBody[3]={0.0,0.0,0.0};
thread_local double gCollisionLsoToBody[36]={0.0};
thread_local double gCollisionBodyToLso[36]={0.0};
std::atomic<std::uint64_t> gProductionCollisionCount{0};

double Dot(const std::array<double,3>& a,const std::array<double,3>& b) {
  return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}

std::array<double,3> Cross(const std::array<double,3>& a,
                           const std::array<double,3>& b) {
  return {{a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],
           a[0]*b[1]-a[1]*b[0]}};
}

std::array<double,3> Subtract(const std::array<double,3>& a,
                              const std::array<double,3>& b) {
  return {{a[0]-b[0],a[1]-b[1],a[2]-b[2]}};
}

double Norm(const std::array<double,3>& value) {
  return std::sqrt(Dot(value,value));
}

bool ValidVector(const std::array<double,3>& value) {
  return std::isfinite(value[0]) && std::isfinite(value[1]) &&
      std::isfinite(value[2]);
}

bool PositiveFinite(double value) {
  return std::isfinite(value) && value>0.0;
}

bool BoxIntersection(const std::array<double,3>& origin,
                     const std::array<double,3>& direction,
                     const std::array<double,3>& minimum,
                     const std::array<double,3>& maximum,
                     double lowerLimit,double upperLimit) {
  for (int dimension=0;dimension<3;++dimension) {
    // Outward floating-point enclosure: a mathematically closed shared
    // vertex must not be culled by a rounded ray/slab comparison before the
    // triangle predicate sees it. This is a coordinate-scaled arithmetic
    // bound, not a terrain-resolution or fitted physical tolerance.
    const double rounding=64.0*std::numeric_limits<double>::epsilon()*
        std::max({std::abs(origin[dimension]),std::abs(minimum[dimension]),
                  std::abs(maximum[dimension])});
    const double low=minimum[dimension]-rounding;
    const double high=maximum[dimension]+rounding;
    if (direction[dimension]==0.0) {
      if (origin[dimension]<low || origin[dimension]>high) return false;
      continue;
    }
    double first=(low-origin[dimension])/direction[dimension];
    double second=(high-origin[dimension])/direction[dimension];
    if (first>second) std::swap(first,second);
    lowerLimit=std::max(lowerLimit,first);
    upperLimit=std::min(upperLimit,second);
    if (lowerLimit>upperLimit) return false;
  }
  return true;
}

bool TriangleIntersection(const Surface::Triangulation& mesh,
                          std::size_t faceIndex,
                          const std::array<double,3>& origin,
                          const std::array<double,3>& direction,
                          double minimumDistance,double maximumDistance,
                          RayHit* candidate) {
  const auto& face=mesh.faces[faceIndex];
  const auto& a=mesh.verticesM[face[0]];
  const auto& b=mesh.verticesM[face[1]];
  const auto& c=mesh.verticesM[face[2]];
  const auto edge0=Subtract(b,a);
  const auto edge1=Subtract(c,a);
  const auto p=Cross(direction,edge1);
  const double determinant=Dot(edge0,p);
  // A scale-relative parallel test avoids assigning an arbitrary length unit
  // to a dimensionless barycentric determinant.
  const double scale=Norm(edge0)*Norm(edge1)*Norm(direction);
  if (!(scale>0.0) || std::abs(determinant)<=
      64.0*std::numeric_limits<double>::epsilon()*scale) return false;
  const double inverse=1.0/determinant;
  const auto relative=Subtract(origin,a);
  const double u=Dot(relative,p)*inverse;
  // A closed triangle contains its edges and vertices. Strict u/v>=0 tests
  // created real holes when radial rays met shared LOLA vertices. Enclose
  // the barycentric arithmetic error using machine epsilon, origin/edge
  // scale and determinant conditioning. Interior/exterior points farther
  // than this arithmetic uncertainty retain the original classification.
  const double barycentricRoundoff=64.0*std::numeric_limits<double>::epsilon()*
      (1.0+Norm(relative)/std::min(Norm(edge0),Norm(edge1)))*
      scale/std::abs(determinant);
  if (u < -barycentricRoundoff || u>1.0+barycentricRoundoff) return false;
  const auto q=Cross(relative,edge0);
  const double v=Dot(direction,q)*inverse;
  if (v < -barycentricRoundoff || u+v>1.0+barycentricRoundoff) return false;
  const double distance=Dot(edge1,q)*inverse;
  if (distance<minimumDistance || distance>maximumDistance) return false;
  auto normal=Cross(edge0,edge1);
  const double length=Norm(normal);
  if (!(length>0.0)) return false;
  for (double& component:normal) component/=length;
  candidate->hit=true;
  candidate->distanceM=distance;
  candidate->face=faceIndex;
  for (int dimension=0;dimension<3;++dimension)
    candidate->positionM[dimension]=origin[dimension]+
        distance*direction[dimension];
  candidate->outwardNormal=normal;
  return true;
}

bool EvaluateSignedTerrainDistance(double segmentTimeS,
    const double* segmentPositionLsoM,const double* segmentVelocityLsoMPerS,
    const double* segmentAccelerationLsoMPerS2,
    double absoluteEpochS,double* signedDistanceM,RayHit* radialHit,
    std::array<double,3>* positionBodyM,std::string* error) {
#ifdef _NO_SPICE_CALLS_
  // A moving body-fixed facet cannot be located in the LSO frame without the
  // epoch-dependent SPICE transform.  Do not substitute a frozen or identity
  // transform: that would silently change the physical collision geometry.
  (void)segmentTimeS;
  (void)segmentPositionLsoM;
  (void)segmentVelocityLsoMPerS;
  (void)segmentAccelerationLsoMPerS2;
  (void)absoluteEpochS;
  (void)signedDistanceM;
  (void)radialHit;
  (void)positionBodyM;
  if (error!=nullptr)
    *error="moving LOLA collision geometry requires a SPICE-enabled build";
  return false;
#else
  if (gProductionTerrain==nullptr || signedDistanceM==nullptr ||
      radialHit==nullptr || positionBodyM==nullptr || error==nullptr)
    return false;
  double positionLso[3];
  for (int dimension=0;dimension<3;++dimension)
    positionLso[dimension]=segmentPositionLsoM[dimension]+
        segmentTimeS*segmentVelocityLsoMPerS[dimension]+
        0.5*segmentTimeS*segmentTimeS*segmentAccelerationLsoMPerS2[dimension];
  SpiceDouble rotation[3][3];
  pxform_c(Frames::SolarOrbital,Frames::BodyFixed,absoluteEpochS,rotation);
  for (int row=0;row<3;++row) {
    (*positionBodyM)[row]=0.0;
    for (int column=0;column<3;++column)
      (*positionBodyM)[row]+=rotation[row][column]*positionLso[column];
  }
  const double radius=Norm(*positionBodyM);
  if (!(radius>0.0)) {
    *error="moving LOLA collision query reached the lunar centre";
    return false;
  }
  std::array<double,3> direction=*positionBodyM;
  for (double& component:direction) component/=radius;
  const std::array<double,3> origin{{0.0,0.0,0.0}};
  if (!gProductionTerrain->IntersectRay(origin,direction,0.0,
          std::numeric_limits<double>::infinity(),
          std::numeric_limits<std::size_t>::max(),radialHit,error)) return false;
  if (!radialHit->hit) {
    *error="body-fixed radial ray did not intersect the production LOLA mesh";
    return false;
  }
  *signedDistanceM=radius-radialHit->distanceM;
  return true;
#endif
}

}  // namespace

bool TerrainRayTracer::Initialize(const Surface::Triangulation& triangulation,
                                  std::string* error) {
  if (error==nullptr) return false;
  error->clear();
  if (triangulation.verticesM.empty() || triangulation.faces.empty()) {
    *error="terrain ray tracer requires a nonempty triangulation";
    return false;
  }
  for (const auto& face:triangulation.faces) {
    if (face[0]>=triangulation.verticesM.size() ||
        face[1]>=triangulation.verticesM.size() ||
        face[2]>=triangulation.verticesM.size() ||
        face[0]==face[1] || face[1]==face[2] || face[0]==face[2]) {
      *error="terrain ray tracer received invalid face connectivity";
      return false;
    }
  }
  mesh_=triangulation;
  faceOrder_.resize(mesh_.faces.size());
  std::iota(faceOrder_.begin(),faceOrder_.end(),0);
  nodes_.clear();
  nodes_.reserve(2*mesh_.faces.size());
  BuildNode(0,faceOrder_.size());
  return true;
}

int TerrainRayTracer::BuildNode(std::size_t begin,std::size_t end) {
  const int nodeIndex=static_cast<int>(nodes_.size());
  nodes_.push_back(Node{});
  Node bounds;
  bounds.begin=begin;
  bounds.end=end;
  bounds.minimum={{std::numeric_limits<double>::infinity(),
      std::numeric_limits<double>::infinity(),
      std::numeric_limits<double>::infinity()}};
  bounds.maximum={{-std::numeric_limits<double>::infinity(),
      -std::numeric_limits<double>::infinity(),
      -std::numeric_limits<double>::infinity()}};
  std::array<double,3> centroidMinimum=bounds.minimum;
  std::array<double,3> centroidMaximum=bounds.maximum;
  for (std::size_t ordered=begin;ordered<end;++ordered) {
    const auto& face=mesh_.faces[faceOrder_[ordered]];
    std::array<double,3> centroid{{0.0,0.0,0.0}};
    for (std::size_t vertex:face) for (int dimension=0;dimension<3;++dimension) {
      const double value=mesh_.verticesM[vertex][dimension];
      bounds.minimum[dimension]=std::min(bounds.minimum[dimension],value);
      bounds.maximum[dimension]=std::max(bounds.maximum[dimension],value);
      centroid[dimension]+=value/3.0;
    }
    for (int dimension=0;dimension<3;++dimension) {
      centroidMinimum[dimension]=std::min(centroidMinimum[dimension],
                                          centroid[dimension]);
      centroidMaximum[dimension]=std::max(centroidMaximum[dimension],
                                          centroid[dimension]);
    }
  }
  nodes_[nodeIndex]=bounds;
  if (end-begin<=kLeafFaces) return nodeIndex;
  int axis=0;
  for (int dimension=1;dimension<3;++dimension)
    if (centroidMaximum[dimension]-centroidMinimum[dimension]>
        centroidMaximum[axis]-centroidMinimum[axis]) axis=dimension;
  const std::size_t middle=begin+(end-begin)/2;
  std::nth_element(faceOrder_.begin()+begin,faceOrder_.begin()+middle,
      faceOrder_.begin()+end,[this,axis](std::size_t left,std::size_t right) {
        double leftCentroid=0.0,rightCentroid=0.0;
        for (std::size_t vertex:mesh_.faces[left])
          leftCentroid+=mesh_.verticesM[vertex][axis];
        for (std::size_t vertex:mesh_.faces[right])
          rightCentroid+=mesh_.verticesM[vertex][axis];
        return leftCentroid<rightCentroid;
      });
  const int left=BuildNode(begin,middle);
  const int right=BuildNode(middle,end);
  nodes_[nodeIndex].left=left;
  nodes_[nodeIndex].right=right;
  return nodeIndex;
}

bool TerrainRayTracer::IntersectRay(const std::array<double,3>& originM,
    const std::array<double,3>& direction,double minimumDistanceM,
    double maximumDistanceM,std::size_t excludedFace,RayHit* hit,
    std::string* error) const {
  if (error==nullptr || hit==nullptr) return false;
  error->clear();
  *hit=RayHit{};
  const double directionLength=Norm(direction);
  if (nodes_.empty() || !ValidVector(originM) || !ValidVector(direction) ||
      !(directionLength>0.0) || !std::isfinite(minimumDistanceM) ||
      minimumDistanceM<0.0 || std::isnan(maximumDistanceM) ||
      maximumDistanceM<minimumDistanceM) {
    *error="invalid terrain ray-intersection request";
    return false;
  }
  std::array<double,3> unitDirection=direction;
  for (double& component:unitDirection) component/=directionLength;
  std::vector<int> stack(1,0);
  while (!stack.empty()) {
    const int nodeIndex=stack.back();
    stack.pop_back();
    const Node& node=nodes_[nodeIndex];
    if (!BoxIntersection(originM,unitDirection,node.minimum,node.maximum,
                         minimumDistanceM,
                         std::min(maximumDistanceM,hit->distanceM))) continue;
    if (node.left>=0) {
      stack.push_back(node.left);
      stack.push_back(node.right);
      continue;
    }
    for (std::size_t ordered=node.begin;ordered<node.end;++ordered) {
      const std::size_t face=faceOrder_[ordered];
      if (face==excludedFace) continue;
      RayHit candidate;
      if (TriangleIntersection(mesh_,face,originM,unitDirection,
              minimumDistanceM,std::min(maximumDistanceM,hit->distanceM),
              &candidate)) *hit=candidate;
    }
  }
  return true;
}

bool TerrainRayTracer::IsSunVisible(
    const std::array<double,3>& surfacePositionM,
    const std::array<double,3>& outwardNormal,
    const std::array<double,3>& unitSunDirection,double rayOriginOffsetM,
    std::size_t sourceFace,bool* visible,std::string* error) const {
  if (visible==nullptr || error==nullptr) return false;
  const double normalLength=Norm(outwardNormal);
  const double sunLength=Norm(unitSunDirection);
  if (!ValidVector(surfacePositionM) || !(normalLength>0.0) ||
      !(sunLength>0.0) || !std::isfinite(rayOriginOffsetM) ||
      rayOriginOffsetM<=0.0) {
    *error="invalid terrain-illumination request";
    return false;
  }
  std::array<double,3> normal=outwardNormal,sun=unitSunDirection,origin;
  for (int dimension=0;dimension<3;++dimension) {
    normal[dimension]/=normalLength;
    sun[dimension]/=sunLength;
    origin[dimension]=surfacePositionM[dimension]+
        rayOriginOffsetM*normal[dimension];
  }
  // A facet whose outward normal faces away from the Sun is self-shadowed
  // before any terrain ray is launched.
  if (Dot(normal,sun)<=0.0) {
    *visible=false;
    error->clear();
    return true;
  }
  RayHit obstruction;
  if (!IntersectRay(origin,sun,0.0,std::numeric_limits<double>::infinity(),
                    sourceFace,&obstruction,error)) return false;
  *visible=!obstruction.hit;
  return true;
}

bool ApparentSolarSystemGeometryAtPoint(double et,
    const std::array<double,3>& pointM,const std::string& abcorr,
    std::array<double,3>* sunRelativeM,std::array<double,3>* earthRelativeM,
    std::string* error) {
  if (error==nullptr) return false;
  error->clear();
  if (sunRelativeM==nullptr || earthRelativeM==nullptr || sunRelativeM==earthRelativeM ||
      !std::isfinite(et) || !ValidVector(pointM) ||
      (abcorr!="NONE" && abcorr!="LT" && abcorr!="LT+S" && abcorr!="CN" && abcorr!="CN+S")) {
    *error="invalid fixed-site illumination geometry: finite TDB epoch/point, distinct outputs and explicit reception correction required";
    return false;
  }
#ifdef _NO_SPICE_CALLS_
  *error="fixed-site illumination geometry requires a SPICE-enabled build; no geometric/apparent fallback is supplied";
  return false;
#else
  // NAIF spkcpo_c: obspos is km, constant in obsref relative to obsctr.
  // OBSERVER evaluates the non-inertial output axes at reception et. CSPICE
  // then computes light time AND observer-velocity stellar aberration for
  // this actual rotating site. Correcting the Moon-centred target vector
  // first and subtracting pointM is not equivalent for non-NONE corrections.
  // Primary API/semantics: https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/spkcpo_c.html
  SpiceDouble observerKm[3],sun[6],earth[6],lightTime=0;
  for (int k=0;k<3;++k) observerKm[k]=pointM[k]/1000.0;
  if (failed_c()) {
    *error="SPICE already has a pending error before the fixed-site geometry query";
    return false; // Do not erase an error owned by a preceding caller.
  }
  // This bool API must report a kernel/frame failure without publishing one
  // target and leaving the other stale. Preserve the caller's error policy;
  // only errors originating in these two calls are cleared before restoring
  // it. No kernel pool, ephemeris epoch or print setting is changed.
  SpiceChar previousAction[SPICE_ERROR_SMSGLN];
  SpiceChar returnAction[]="RETURN";
  erract_c("GET",sizeof(previousAction),previousAction);
  erract_c("SET",0,returnAction);
  spkcpo_c("SUN",et,Frames::BodyFixed,"OBSERVER",abcorr.c_str(),observerKm,
           ObjectName,Frames::BodyFixed,sun,&lightTime);
  if (!failed_c())
    spkcpo_c("EARTH",et,Frames::BodyFixed,"OBSERVER",abcorr.c_str(),observerKm,
             ObjectName,Frames::BodyFixed,earth,&lightTime);
  const bool failed=failed_c();
  SpiceChar detail[SPICE_ERROR_LMSGLN]={0};
  if (failed) { getmsg_c("LONG",sizeof(detail),detail);reset_c(); }
  erract_c("SET",0,previousAction);
  if (failed) {
    *error=std::string("SPICE fixed-site Sun/Earth geometry failed: ")+detail;
    return false;
  }
  std::array<double,3> candidateSun,candidateEarth;
  for (int k=0;k<3;++k) { candidateSun[k]=1000.0*sun[k];candidateEarth[k]=1000.0*earth[k]; }
  const double sunRangeM=std::hypot(candidateSun[0],candidateSun[1],candidateSun[2]);
  const double earthRangeM=std::hypot(candidateEarth[0],candidateEarth[1],candidateEarth[2]);
  if (!ValidVector(candidateSun) || !ValidVector(candidateEarth) ||
      !PositiveFinite(sunRangeM) || !PositiveFinite(earthRangeM)) {
    *error="SPICE fixed-site geometry returned a nonfinite or zero target distance";
    return false;
  }
  *sunRelativeM=candidateSun;*earthRelativeM=candidateEarth;
  return true;
#endif
}

bool EvaluateSolarDisc(const TerrainRayTracer* terrain,
    const SolarDiscRequest& request,SolarDiscSample* output,
    std::string* error,const std::function<double(double)>& absorptivity) {
  if (output==nullptr || error==nullptr) return false;
  error->clear();
  const double sunDistance=Norm(request.pointToSunM);
  const double normalLength=Norm(request.outwardNormal);
  const double earthDistance=Norm(request.pointToEarthM);
  if (!ValidVector(request.surfacePositionM) ||
      !ValidVector(request.outwardNormal) ||
      !ValidVector(request.pointToSunM) || !ValidVector(request.pointToEarthM) ||
      !PositiveFinite(normalLength) || !PositiveFinite(sunDistance) ||
      !PositiveFinite(request.solarRadiusM) ||
      !(sunDistance>request.solarRadiusM) ||
      !std::isfinite(request.opaqueEarthRadiusM) || request.opaqueEarthRadiusM<0 ||
      (request.opaqueEarthRadiusM>0 && (!PositiveFinite(earthDistance) ||
          !(earthDistance>request.opaqueEarthRadiusM))) ||
      // A surface self-intersection displacement must never leak into a
      // gas-particle observation operator. In volume mode the same production
      // rays start at the exact requested position, without an artificial
      // facet horizon or any ignored physical BVH triangle.
      (request.volumePointOrigin ?
          (request.applyLocalSelfHorizon || request.rayOriginOffsetM!=0.0 ||
           request.sourceFace!=std::numeric_limits<std::size_t>::max()) :
          !PositiveFinite(request.rayOriginOffsetM)) || request.radialPoints<1 ||
      request.azimuthPoints<4 || request.azimuthPoints%4!=0 ||
      static_cast<std::size_t>(request.radialPoints)>
          std::numeric_limits<std::size_t>::max()/request.azimuthPoints) {
    *error="invalid finite-disc geometry/grid: explicit positive radii, "
        "external observer and radial>=1, azimuth>=4 divisible by four required; "
        "surface mode needs positive origin offset, volume mode needs zero "
        "offset, applyLocalSelfHorizon=false and no ignored source facet";
    return false;
  }
  std::array<double,3> centre=request.pointToSunM,normal=request.outwardNormal;
  for (int k=0;k<3;++k) { centre[k]/=sunDistance; normal[k]/=normalLength; }
  // Pick the least aligned Cartesian axis; no pole/special-frame singularity.
  int axis=0;
  for (int k=1;k<3;++k) if (std::abs(centre[k])<std::abs(centre[axis])) axis=k;
  std::array<double,3> auxiliary{{0,0,0}}; auxiliary[axis]=1;
  auto tangent=Cross(centre,auxiliary);
  const double tangentLength=Norm(tangent);
  for (double& x:tangent) x/=tangentLength;
  const auto second=Cross(centre,tangent);
  const double sinRadius=request.solarRadiusM/sunDistance;
  const double pi=std::acos(-1.0);
  SolarDiscSample candidate;
  candidate.rayCount=static_cast<std::size_t>(request.radialPoints)*request.azimuthPoints;
  const double weight=1.0/static_cast<double>(candidate.rayCount);
  double fullSolidAngleWeight=0.0;
  // Use sin(alpha)^2/(1+cos(alpha)) for 1-cos(alpha), retaining precision for
  // the real Sun's small angular radius. This is geometry, not a replacement
  // for the selected D04 1-AU spectrum/rate-distance convention.
  candidate.unobscuredSolidAngleSr=2*pi*sinRadius*sinRadius/
      (1+std::sqrt(1-sinRadius*sinRadius));
  for (int ring=0;ring<request.radialPoints;++ring) {
    // Midpoint in r^2 gives equal PROJECTED disc areas. With uniform radiance,
    // dOmega=dA_projected/cos(rho); therefore irradiance needs the 1/cos(rho)
    // Jacobian below. Averaging cos(incidence) alone would be a different law.
    const double sine=sinRadius*std::sqrt((ring+0.5)/request.radialPoints);
    const double cosine=std::sqrt(1.0-sine*sine);
    const double solidAngleWeight=weight/cosine;
    for (int azimuth=0;azimuth<request.azimuthPoints;++azimuth) {
      const double phi=2*pi*(azimuth+0.5)/request.azimuthPoints;
      fullSolidAngleWeight+=solidAngleWeight;
      std::array<double,3> ray;
      for (int k=0;k<3;++k)
        ray[k]=cosine*centre[k]+sine*(std::cos(phi)*tangent[k]+std::sin(phi)*second[k]);
      const double incidence=std::max(-1.0,std::min(1.0,Dot(normal,ray)));
      bool terrainVisible=!request.applyLocalSelfHorizon || incidence>0;
      if (terrain && terrainVisible) {
        if (request.applyLocalSelfHorizon) {
          if (!terrain->IsSunVisible(request.surfacePositionM,normal,ray,
                  request.rayOriginOffsetM,request.sourceFace,&terrainVisible,error)) return false;
        } else {
          std::array<double,3> origin=request.surfacePositionM;
          if (!request.volumePointOrigin)
            for (int k=0;k<3;++k) origin[k]+=request.rayOriginOffsetM*normal[k];
          RayHit obstruction;
          if (!terrain->IntersectRay(origin,ray,0,std::numeric_limits<double>::infinity(),
                  request.sourceFace,&obstruction,error)) return false;
          terrainVisible=!obstruction.hit;
        }
      }
      bool earthVisible=true;
      if (request.opaqueEarthRadiusM>0) {
        // Intersect the SAME ray with the opaque Earth. Multiplying separate
        // terrain/occultation area fractions would be wrong for overlapping
        // masks (e.g. a ridge covers the very part eclipsed by Earth).
        const double along=Dot(request.pointToEarthM,ray);
        const double perpendicular2=std::max(0.0,earthDistance*earthDistance-along*along);
        const double discriminant=request.opaqueEarthRadiusM*request.opaqueEarthRadiusM-perpendicular2;
        if (along>0 && discriminant>=0) {
          const double nearEarth=along-std::sqrt(discriminant);
          const double alongSun=Dot(request.pointToSunM,ray);
          const double sunDiscriminant=std::max(0.0,request.solarRadiusM*request.solarRadiusM-
              std::max(0.0,sunDistance*sunDistance-alongSun*alongSun));
          const double nearSun=alongSun-std::sqrt(sunDiscriminant);
          earthVisible=!(nearEarth>=0 && nearEarth<nearSun);
        }
      }
      if (earthVisible) {
        candidate.earthVisibleAreaFraction+=weight;
        candidate.earthVisibleSolidAngleFraction+=solidAngleWeight;
      }
      if (!terrainVisible) continue;
      candidate.terrainVisibleAreaFraction+=weight;
      candidate.terrainVisibleSolidAngleFraction+=solidAngleWeight;
      candidate.sampledAnyTerrainVisible=true;
      if (!earthVisible) continue;
      candidate.jointlyVisibleAreaFraction+=weight;
      candidate.jointlyVisibleSolidAngleFraction+=solidAngleWeight;
      candidate.sampledAnyJointlyVisible=true;
      if (incidence<=0) continue; // Visible sky does not imply inward facet heating.
      const double projected=weight*incidence/cosine;
      const double absorption=absorptivity ? absorptivity(incidence):1.0;
      if (!std::isfinite(absorption) || absorption<0 || absorption>1) {
        *error="finite-disc directional absorptivity is outside [0,1]";
        return false;
      }
      candidate.projectedIncidence+=projected;
      candidate.absorbedProjectedIncidence+=projected*absorption;
    }
  }
  // Enclose accumulated unit-fraction roundoff, not a physical correction.
  candidate.terrainVisibleAreaFraction=std::min(1.0,candidate.terrainVisibleAreaFraction);
  candidate.earthVisibleAreaFraction=std::min(1.0,candidate.earthVisibleAreaFraction);
  candidate.jointlyVisibleAreaFraction=std::min(1.0,candidate.jointlyVisibleAreaFraction);
  // Normalize the SAME jointly masked rays, not independent Earth and terrain
  // fractions. Clamping encloses only accumulated floating-point roundoff;
  // no missing angular data or source intensity is filled/renormalized here.
  candidate.terrainVisibleSolidAngleFraction=std::min(1.0,
      candidate.terrainVisibleSolidAngleFraction/fullSolidAngleWeight);
  candidate.earthVisibleSolidAngleFraction=std::min(1.0,
      candidate.earthVisibleSolidAngleFraction/fullSolidAngleWeight);
  candidate.jointlyVisibleSolidAngleFraction=std::min(1.0,
      candidate.jointlyVisibleSolidAngleFraction/fullSolidAngleWeight);
  *output=candidate;
  return true;
}

bool TemporalVisibility::Add(const SolarDiscSample& sample,std::string* error) {
  if (error==nullptr) return false;
  error->clear();
  // Forward bound for positive quadrature accumulation, plus fixed per-ray
  // geometry arithmetic. This scales with the ACTUAL configured grid; a
  // legitimate fine non-power-of-two grid must not fail a fixed-epsilon test.
  const double arithmetic=(64.0+static_cast<double>(sample.rayCount))*
      std::numeric_limits<double>::epsilon();
  const double rounding=arithmetic<1 ? arithmetic/(1-arithmetic):
      std::numeric_limits<double>::infinity();
  if (!sample.rayCount || epochCount==std::numeric_limits<std::size_t>::max() ||
      !std::isfinite(sample.jointlyVisibleAreaFraction) ||
      sample.jointlyVisibleAreaFraction<0 || sample.jointlyVisibleAreaFraction>1 ||
      !std::isfinite(sample.projectedIncidence) || sample.projectedIncidence<0 ||
      arithmetic>=1 || sample.projectedIncidence>1+rounding ||
      sample.sampledAnyJointlyVisible!=(sample.jointlyVisibleAreaFraction>0)) {
    *error="invalid equal-epoch visibility sample or counter overflow";
    return false;
  }
  ++epochCount;
  if (sample.sampledAnyJointlyVisible) ++sampledAnyVisibleCount;
  meanVisibleAreaFraction+=(sample.jointlyVisibleAreaFraction-meanVisibleAreaFraction)/epochCount;
  meanProjectedIncidence+=(sample.projectedIncidence-meanProjectedIncidence)/epochCount;
  return true;
}

double TemporalVisibility::SampledAnyDiscOccurrenceFraction() const {
  return epochCount ? static_cast<double>(sampledAnyVisibleCount)/epochCount:
      std::numeric_limits<double>::quiet_NaN();
}

bool InitializeProductionTerrain(const Surface::Triangulation& triangulation,
                                 const Runtime::Configuration& configuration,
                                 std::string* error) {
  if (error==nullptr) return false;
  auto candidate=std::make_unique<TerrainRayTracer>();
  if (!candidate->Initialize(triangulation,error)) return false;
  gProductionTerrain=std::move(candidate);
  gProductionConfiguration=&configuration;
  return true;
}

bool ProductionTerrainInitialized() { return gProductionTerrain!=nullptr; }

const TerrainRayTracer* ProductionTerrain() {
  return gProductionTerrain.get();
}

bool ProductionIonCollectorInstalled() { return gIonCollectorResolver!=nullptr; }

bool InstallProductionIonCollector(
    const std::vector<Species::IonCollector::StateDefinition>& states,
    IonCollectorStateResolver resolver,void* user,
    Species::IonCollector::Ledger* ledger,std::string* error) {
  return InstallProductionIonCollector(states,resolver,user,ledger,error,
      ElectronCharge,"AMPS legacy ElectronCharge");
}
bool InstallProductionIonCollector(
    const std::vector<Species::IonCollector::StateDefinition>& states,
    IonCollectorStateResolver resolver,void* user,
    Species::IonCollector::Ledger* ledger,std::string* error,
    double elementaryChargeC,const std::string& elementaryChargeSourceIdentity) {
  if(!error)return false;
#if _COMPILATION_MODE_ != _COMPILATION_MODE__MPI_
  *error="ideal ion collector currently requires MPI-only execution (OMP=1); shared-ledger thread publication is unqualified";return false;
#else
  if(!std::isfinite(elementaryChargeC) || elementaryChargeC<=0 || elementaryChargeSourceIdentity.empty() ||
     ProductionIonCollectorInstalled() || !Runtime::HasConfiguration() ||
     Runtime::GetConfiguration().ionSurfaceBoundary!="ideal-collector" ||
     Runtime::GetConfiguration().particleTracingMover=="legacy" || !resolver || !ledger || !ledger->records.empty()) {
    *error="install ideal ion collector exactly once with explicit selector, physical-time tracer, state resolver and an empty ledger; unsupported/default mode cannot install a collector";return false;
  }
  if(!Species::IonCollector::ValidateStates(states,error))return false;
  for(const auto& state:states) {
    const auto& m=state.material;
    if(m.compiledSpecies>=PIC::nTotalSpecies ||
       PIC::MolecularData::GetMass(m.compiledSpecies)!=m.modelMassKg ||
       PIC::MolecularData::GetElectricCharge(m.compiledSpecies)!=m.elementaryCharge*elementaryChargeC) {
      *error="collector state inventory disagrees with actual compiled species/mass/charge; configure matching isotope/reference masses or supply an explicitly matching inventory, never alter ledger masses alone";return false;
    }
  }
  try {
    auto candidate=states;
    gIonCollectorStates.swap(candidate);
  } catch(const std::exception&) {
    *error="collector state-inventory allocation failed; provide more memory or reduce the explicit state catalog";return false;
  }
  // Publish only after every fallible validation/allocation. A rejected
  // installation leaves the previously unsupported physical policy intact.
  gIonCollectorElementaryChargeC=elementaryChargeC;
  gIonCollectorChargeSourceIdentity=elementaryChargeSourceIdentity;
  gIonCollectorUser=user;gIonCollectorLedger=ledger;gIonCollectorResolver=resolver;
  error->clear();return true;
#endif
}

bool HandleProductionIonImpact(long int particle,const double* position,
    const double* velocity,void* nodePointer,
    const Species::IonCollector::ImpactFrame& frame) {
  if(particle<0 || particle>=PIC::ParticleBuffer::GetMaxNPart() ||
     !PIC::ParticleBuffer::IsParticleAllocated(particle))
    exit(__LINE__,__FILE__,"ion-boundary dispatch requires an actually allocated incident PIC particle");
  const int species=PIC::ParticleBuffer::GetI(particle);
  if(species<0 || species>=PIC::nTotalSpecies)
    exit(__LINE__,__FILE__,"ion-boundary incident particle has an invalid compiled species");
  const double charge=PIC::MolecularData::GetElectricCharge(species);
  if(charge==0)return false; // Preserve neutral arithmetic/RNG/callback order.
  if(!std::isfinite(charge) || !Runtime::HasConfiguration() ||
     Runtime::GetConfiguration().ionSurfaceBoundary!="ideal-collector")
    exit(__LINE__,__FILE__,"charged lunar surface impact is unsupported by default: no neutral Na/Ar sticking, neutralization, implantation or heat law is applied. For mathematical first-impact diagnostics explicitly select ion_surface_boundary=ideal-collector and install a matching state/energy inventory; otherwise disable ion production or supply qualified ion-surface physics.");
  if(!ProductionIonCollectorInstalled() || !gIonCollectorLedger || !position || !velocity || !nodePointer)
    exit(__LINE__,__FILE__,"ideal ion collector was selected without an installed state/formation-energy resolver, ledger or actual impact node. An input switch alone cannot supply missing physical states/energies; install explicit source-bound state data before enabling ion transport.");
  auto* node=static_cast<cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>*>(nodePointer);
  if(!node->block)exit(__LINE__,__FILE__,"ideal collector impact has no allocated owning AMR block");
  std::size_t stateIndex=0;std::string error;
  if(!gIonCollectorResolver(particle,nodePointer,gIonCollectorUser,&stateIndex,&error) ||
     stateIndex>=gIonCollectorStates.size()) {
    const std::string message="ideal collector cannot identify actual incident source/product state; do not infer a ground state or formation energy from species alone: "+error;
    exit(__LINE__,__FILE__,message.c_str());
  }
  Species::IonCollector::Incident incident;
  incident.state=gIonCollectorStates[stateIndex];
  const auto& identity=incident.state.material;
  if(identity.compiledSpecies!=species || PIC::MolecularData::GetMass(species)!=identity.modelMassKg ||
     charge!=identity.elementaryCharge*gIonCollectorElementaryChargeC)
    exit(__LINE__,__FILE__,"incident collector state disagrees with actual particle species/mass/charge; source-state crosswalk is incompatible");
  incident.representedNumber=node->block->GetLocalParticleWeight(species)*
      PIC::ParticleBuffer::GetIndividualStatWeightCorrection(particle);
  incident.elementaryChargeC=gIonCollectorElementaryChargeC;
  for(int d=0;d<3;++d) { incident.positionLsoM[d]=position[d];incident.velocityLsoMS[d]=velocity[d]; }
  incident.impactEpochEtS=PIC::Mover::ForceCutSurface::GetParticleEpochS(particle);
  incident.frame=frame;
  if(!Species::IonCollector::Collect(incident,gIonCollectorLedger,&error))
    exit(__LINE__,__FILE__,error.c_str());
  // The caller returns Delete. We deliberately leave bytes/allocator/cell-list
  // ownership to the mover, after successful ledger publication exactly once.
  return true;
}

void BeginProductionParticleStep(double ephemerisTimeBeginS) {
  if (!std::isfinite(ephemerisTimeBeginS))
    exit(__LINE__,__FILE__,"invalid moving-surface particle-step epoch");
  gParticleStepBeginEtS=ephemerisTimeBeginS;
  // Synchronously prepare direct-photon factors before source callbacks can
  // run in particle threads. The callbacks then read immutable coefficients
  // at one declared SOURCE-BEGIN epoch, with no CSPICE/RNG/cache mutation.
  if (Moon::SurfaceTemperature::ProductionStateInitialized()) {
    std::string error;
    if (!Moon::SurfaceTemperature::RefreshProductionPhotonFluxCache(ephemerisTimeBeginS,&error))
      exit(__LINE__,__FILE__,error.c_str());
  }
}
double ProductionSourceEpochS() { return gParticleStepBeginEtS; }
double ProductionInjectionEpochS() {
  if(Moon::CommonTransportTimeSelected()) {
    const auto* clock=Moon::GetCommonTransportClock();TransportTime::BirthTime birth;
    std::string error;const double offset=std::isfinite(PIC::Mover::MovingSurfaceInjectionTimeOffsetS)?
      PIC::Mover::MovingSurfaceInjectionTimeOffsetS:0.0;
    if(!clock || !clock->BirthFromOffset(offset,&birth,&error)) {
      const std::string message="common-time source birth requires an active normal transport step and an actual offset within [0,dt]: "+error+
        ". Call normal Begin/Commit hooks, retain explicit birth offsets and restart from the declared initial state; do not reuse the frozen external epoch as transport time.";
      exit(__LINE__,__FILE__,message.c_str());
    }
    return birth.epochEtS;
  }
  return gParticleStepBeginEtS+(
      std::isfinite(PIC::Mover::MovingSurfaceInjectionTimeOffsetS) ?
      PIC::Mover::MovingSurfaceInjectionTimeOffsetS:0.0);
}

bool ApplyFrozenTracingBodyTransform(const Runtime::Configuration& configuration,
    double* transform,std::string* error) {
  if(!error)return false;
  error->clear();
  if(!transform){*error="null body-state transform";return false;}
  if(configuration.particleTracingMover=="legacy" ||
      configuration.simulationEnvironmentMode!="frozen")return true;
  for(int i=0;i<36;++i)if(!std::isfinite(transform[i])) {
    *error="nonfinite frozen body-state transform";return false;
  }
  for(int row=3;row<6;++row)for(int column=0;column<3;++column)
    transform[6*row+column]=0;
  return true;
}

bool ProductionInjectionBodyToLsoTransform(double* transform,std::string* error) {
  if(!error)return false;
  error->clear();
  if(!transform || !gProductionConfiguration) {
    *error="source body transform requires initialized production terrain/configuration";return false;
  }
#ifdef _NO_SPICE_CALLS_
  *error="source body transform requires SPICE";return false;
#else
  const auto& configuration=*gProductionConfiguration;
  const bool frozen=configuration.particleTracingMover!="legacy" &&
      configuration.simulationEnvironmentMode=="frozen";
  const double epoch=frozen?ProductionSourceEpochS():ProductionInjectionEpochS();
  SpiceDouble candidate[6][6];
  sxform_c(Frames::BodyFixed,Frames::SolarOrbital,epoch,candidate);
  if(!ApplyFrozenTracingBodyTransform(configuration,&candidate[0][0],error))return false;
  std::memcpy(transform,candidate,sizeof(candidate));return true;
#endif
}

bool CommitProductionMovingSurfaceCollision(long int particle,
    const double* selectedNormalLso,double* normalBody,
    double* lsoToBody36,double* bodyToLso36,std::string* error) {
  if(!error)return false;
  error->clear();
  if(!gProductionConfiguration || gProductionConfiguration->particleTracingMover=="legacy" ||
      !selectedNormalLso || !normalBody || !lsoToBody36 || !bodyToLso36) {
    *error="accepted tracer collision requires explicit production configuration and outputs";return false;
  }
  const double physicalEpoch=PIC::Mover::ForceCutSurface::GetParticleEpochS(particle);
  if(!std::isfinite(physicalEpoch)){*error="accepted collision has no physical particle epoch";return false;}
#ifdef _NO_SPICE_CALLS_
  *error="accepted collision frame requires SPICE";return false;
#else
  const auto& configuration=*gProductionConfiguration;
  const double epoch=configuration.simulationEnvironmentMode=="frozen"?
      ProductionSourceEpochS():physicalEpoch;
  SpiceDouble inverse[6][6],forward[6][6];
  sxform_c(Frames::SolarOrbital,Frames::BodyFixed,epoch,inverse);
  sxform_c(Frames::BodyFixed,Frames::SolarOrbital,epoch,forward);
  if(!ApplyFrozenTracingBodyTransform(configuration,&inverse[0][0],error) ||
      !ApplyFrozenTracingBodyTransform(configuration,&forward[0][0],error))return false;
  double normal[3]={};
  for(int row=0;row<3;++row)for(int column=0;column<3;++column)
    normal[row]+=inverse[row][column]*selectedNormalLso[column];
  const double length=std::hypot(normal[0],normal[1],normal[2]);
  if(!PositiveFinite(length)){*error="accepted collision normal is nonfinite or zero";return false;}
  for(double& value:normal)value/=length;
  // Publish only after the complete accepted frame is valid. Adaptive trial
  // queries cannot overwrite this transaction or increment its physical count.
  std::memcpy(normalBody,normal,sizeof(normal));
  std::memcpy(lsoToBody36,inverse,sizeof(inverse));std::memcpy(bodyToLso36,forward,sizeof(forward));
  gProductionCollisionCount.fetch_add(1,std::memory_order_relaxed);
  return true;
#endif
}

bool FindProductionMovingSurfaceIntersection(long int particle,
    const double* segmentPositionLsoM,const double* segmentVelocityLsoMPerS,
    const double* segmentAccelerationLsoMPerS2,double elapsedFlightTimeS,
    double maximumFlightTimeS,double& hitFlightTimeS,
    double* hitPositionLsoM,double* outwardNormalLso,
    double* surfaceVelocityLsoMPerS) {
#ifdef _NO_SPICE_CALLS_
  // Initialization rejects this production combination as well.  Retain the
  // explicit callback guard so a future caller cannot bypass initialization
  // and accidentally treat the rotating terrain as stationary.
  (void)particle;
  (void)segmentPositionLsoM;
  (void)segmentVelocityLsoMPerS;
  (void)elapsedFlightTimeS;
  (void)maximumFlightTimeS;
  (void)hitFlightTimeS;
  (void)hitPositionLsoM;
  (void)outwardNormalLso;
  (void)surfaceVelocityLsoMPerS;
  exit(__LINE__,__FILE__,
      "moving LOLA collision geometry requires a SPICE-enabled build");
  return false;
#else
  (void)segmentAccelerationLsoMPerS2;
  if (gProductionTerrain==nullptr || gProductionConfiguration==nullptr)
    return false;
  if (segmentPositionLsoM==nullptr || segmentVelocityLsoMPerS==nullptr ||
      segmentAccelerationLsoMPerS2==nullptr ||
      hitPositionLsoM==nullptr || outwardNormalLso==nullptr ||
      surfaceVelocityLsoMPerS==nullptr ||
      !std::isfinite(elapsedFlightTimeS) || elapsedFlightTimeS<0.0 ||
      !PositiveFinite(maximumFlightTimeS))
    exit(__LINE__,__FILE__,"invalid moving LOLA collision segment");

  const Runtime::Configuration& configuration=*gProductionConfiguration;
  const bool localTracing=configuration.particleTracingMover!="legacy";
  const bool frozenTracing=localTracing && configuration.simulationEnvironmentMode=="frozen";
  // The production RK2 position extension is x(t)=x0+v0*t+a0*t^2/2.
  // Searching a straight ray while the mover advances an accelerated segment
  // gives the wrong collision epoch/position. Bound its particle arc speed
  // over this segment by |v0|+|a0|*dt when choosing geometric subdivisions.
  const double speed=std::sqrt(
      segmentVelocityLsoMPerS[0]*segmentVelocityLsoMPerS[0]+
      segmentVelocityLsoMPerS[1]*segmentVelocityLsoMPerS[1]+
      segmentVelocityLsoMPerS[2]*segmentVelocityLsoMPerS[2])+
      maximumFlightTimeS*std::sqrt(
          segmentAccelerationLsoMPerS2[0]*segmentAccelerationLsoMPerS2[0]+
          segmentAccelerationLsoMPerS2[1]*segmentAccelerationLsoMPerS2[1]+
          segmentAccelerationLsoMPerS2[2]*segmentAccelerationLsoMPerS2[2]);
  if (!(speed>0.0)) return false;
  const std::size_t intervals=std::max<std::size_t>(1,
      static_cast<std::size_t>(std::ceil(speed*maximumFlightTimeS/
          configuration.surfaceCollisionMaximumSegmentLengthM)));
  double absoluteSegmentStart=ProductionInjectionEpochS()+elapsedFlightTimeS;
  if (localTracing) {
    const auto clock=PIC::Mover::ForceCutSurface::GetTimeContext(particle);
    if (!std::isfinite(clock.InitialEpochS))
      exit(__LINE__,__FILE__,"local-time LOLA tracing needs a physical particle epoch");
    absoluteSegmentStart=clock.InitialEpochS+elapsedFlightTimeS;
  }
  // The explicit frozen tracing profile holds ALL ephemeris-driven geometry
  // at its declared source epoch as well as the force fields. Physical local
  // seconds still advance in particle bytes. The legacy moving-surface path
  // retains its old epoch semantics. Future evolving tracing must supply its
  // corresponding epoch-aware force provider, not mix frozen forces with an
  // unrelated global-step pose.
  const auto geometryEpoch=[&](double localSeconds) {
    return frozenTracing ? ProductionSourceEpochS() : absoluteSegmentStart+localSeconds;
  };
  std::string error;
  double lowerTime=0.0,lowerDistance=0.0;
  RayHit lowerHit;
  std::array<double,3> lowerBody;
  if (!EvaluateSignedTerrainDistance(0.0,segmentPositionLsoM,
          segmentVelocityLsoMPerS,segmentAccelerationLsoMPerS2,geometryEpoch(0),&lowerDistance,
          &lowerHit,&lowerBody,&error))
    exit(__LINE__,__FILE__,error.c_str());
  if (lowerDistance<-
      configuration.surfaceCollisionDistanceToleranceM)
    exit(__LINE__,__FILE__,
        "particle segment begins below the moving LOLA surface; reduce the "
        "particle time step or collision maximum segment length");

  double upperTime=0.0,upperDistance=lowerDistance;
  RayHit upperHit;
  std::array<double,3> upperBody;
  bool bracketed=false;
  for (std::size_t interval=1;interval<=intervals;++interval) {
    upperTime=maximumFlightTimeS*static_cast<double>(interval)/intervals;
    if (!EvaluateSignedTerrainDistance(upperTime,segmentPositionLsoM,
            segmentVelocityLsoMPerS,segmentAccelerationLsoMPerS2,geometryEpoch(upperTime),
            &upperDistance,&upperHit,&upperBody,&error))
      exit(__LINE__,__FILE__,error.c_str());
    if (lowerDistance>
            configuration.surfaceCollisionDistanceToleranceM &&
        upperDistance<=0.0) {
      bracketed=true;
      break;
    }
    lowerTime=upperTime;
    lowerDistance=upperDistance;
    lowerHit=upperHit;
    lowerBody=upperBody;
  }
  if (!bracketed) return false;

  // Refine the first outside-to-inside bracket.  Both an SI time and an SI
  // signed-distance tolerance must be met; reaching the iteration limit is a
  // fatal resolution/configuration error rather than an approximate impact.
  for (int iteration=0;
       iteration<configuration.surfaceCollisionMaximumIterations;
       ++iteration) {
    const double middleTime=0.5*(lowerTime+upperTime);
    double middleDistance=0.0;
    RayHit middleHit;
    std::array<double,3> middleBody;
    if (!EvaluateSignedTerrainDistance(middleTime,segmentPositionLsoM,
            segmentVelocityLsoMPerS,segmentAccelerationLsoMPerS2,geometryEpoch(middleTime),
            &middleDistance,&middleHit,&middleBody,&error))
      exit(__LINE__,__FILE__,error.c_str());
    if (middleDistance>0.0) {
      lowerTime=middleTime;
      lowerDistance=middleDistance;
      lowerHit=middleHit;
      lowerBody=middleBody;
    }
    else {
      upperTime=middleTime;
      upperDistance=middleDistance;
      upperHit=middleHit;
      upperBody=middleBody;
    }
    if (upperTime-lowerTime<=
            configuration.surfaceCollisionTimeToleranceS &&
        std::min(std::abs(lowerDistance),std::abs(upperDistance))<=
            configuration.surfaceCollisionDistanceToleranceM) {
      hitFlightTimeS=std::abs(lowerDistance)<=std::abs(upperDistance) ?
          lowerTime : upperTime;
      const RayHit& hit=std::abs(lowerDistance)<=std::abs(upperDistance) ?
          lowerHit : upperHit;
      const std::array<double,3>& body=
          std::abs(lowerDistance)<=std::abs(upperDistance) ?
              lowerBody : upperBody;
      for (int dimension=0;dimension<3;++dimension)
        hitPositionLsoM[dimension]=segmentPositionLsoM[dimension]+
            hitFlightTimeS*segmentVelocityLsoMPerS[dimension]+
            0.5*hitFlightTimeS*hitFlightTimeS*segmentAccelerationLsoMPerS2[dimension];
      const double epoch=geometryEpoch(hitFlightTimeS);
      SpiceDouble bodyToLso[6][6],bodyState[6]={body[0],body[1],body[2],
                                                0.0,0.0,0.0};
      SpiceDouble lsoToBody[6][6];
      SpiceDouble lsoState[6];
      sxform_c(Frames::BodyFixed,Frames::SolarOrbital,epoch,bodyToLso);
      sxform_c(Frames::SolarOrbital,Frames::BodyFixed,epoch,lsoToBody);
      if (frozenTracing) {
        // A deliberately frozen pose has zero surface velocity relative to
        // LSO. Retain the measured rotation but remove only the derivative
        // block in BOTH state maps, keeping impact/re-emission transforms
        // consistent with the actual boundary trajectory used above.
        if(!ApplyFrozenTracingBodyTransform(configuration,&bodyToLso[0][0],&error) ||
            !ApplyFrozenTracingBodyTransform(configuration,&lsoToBody[0][0],&error))
          exit(__LINE__,__FILE__,error.c_str());
      }
      mxvg_c(bodyToLso,bodyState,6,6,lsoState);
      for (int row=0;row<3;++row) {
        outwardNormalLso[row]=0.0;
        for (int column=0;column<3;++column)
          outwardNormalLso[row]+=bodyToLso[row][column]*
              hit.outwardNormal[column];
        surfaceVelocityLsoMPerS[row]=lsoState[row+3];
      }
      if(!localTracing) {
      std::memcpy(gCollisionNormalBody,hit.outwardNormal.data(),
                  3*sizeof(double));
      std::memcpy(gCollisionLsoToBody,lsoToBody,36*sizeof(double));
      std::memcpy(gCollisionBodyToLso,bodyToLso,36*sizeof(double));
      gCollisionFrameValid=true;
      gProductionCollisionCount.fetch_add(1,std::memory_order_relaxed);
      }
      return true;
    }
  }
  exit(__LINE__,__FILE__,
      "moving LOLA collision root did not converge; reduce "
      "surface_collision_maximum_segment_length_m or increase the explicit "
      "iteration limit after a documented I14 convergence study");
  return false;
#endif
}

bool LastProductionCollisionFrame(double* outwardNormalBodyFixed,
                                  double* lsoToBodyStateTransform36,
                                  double* bodyToLsoStateTransform36) {
  if (!gCollisionFrameValid || outwardNormalBodyFixed==nullptr ||
      lsoToBodyStateTransform36==nullptr ||
      bodyToLsoStateTransform36==nullptr) return false;
  std::memcpy(outwardNormalBodyFixed,gCollisionNormalBody,3*sizeof(double));
  std::memcpy(lsoToBodyStateTransform36,gCollisionLsoToBody,
              36*sizeof(double));
  std::memcpy(bodyToLsoStateTransform36,gCollisionBodyToLso,
              36*sizeof(double));
  gCollisionFrameValid=false;
  return true;
}

std::uint64_t ProductionCollisionCount() {
  return gProductionCollisionCount.load(std::memory_order_relaxed);
}

}  // namespace SurfaceEnvironment
}  // namespace Moon
