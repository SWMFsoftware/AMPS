#include "MoonTrajectoryExposure.h"
#include "Moon.h"
#include "MoonInput.h"
#include "MoonSurfaceEnvironment.h"
#include <cmath>
#include <limits>

namespace Moon { namespace Chemistry { namespace Exposure {
namespace T=Species::TrajectoryExposure;
namespace {
bool Reject(std::string* error,const char* text){if(error)*error=text;return false;}
T::Vector Rotate(const double* rotation,const T::Vector&x){T::Vector result{};for(int r=0;r<3;++r)for(int c=0;c<3;++c)result[r]+=rotation[6*r+c]*x[c];return result;}

bool Allowed(std::string* error){
  if(!Runtime::HasConfiguration()||Runtime::GetConfiguration().particleTracingMover=="legacy"||
      Runtime::GetConfiguration().simulationEnvironmentMode!="frozen")
    return Reject(error,"trajectory photon geometry requires the declared frozen environment and physical-time tracer");
  return true;
}
bool SunValid(const T::Vector& sun,std::string* error){
  for(double x:sun)if(!std::isfinite(x))return Reject(error,"nonfinite production Sun vector");
  return std::hypot(sun[0],sun[1],sun[2])>_RADIUS_(_SUN_)||
      Reject(error,"production Sun vector is uninitialized/invalid; use explicit CreateAtEpoch before first particle step");
}
bool Surface(T::Field& f,std::array<double,36> forward,std::array<double,36> inverse,std::string* error){
  if(Runtime::GetConfiguration().surfaceGeometry==Runtime::SurfaceGeometry::Sphere)
    f.opaqueSpheres.push_back({{{0,0,0}},_RADIUS_(_TARGET_)});
  else {
    const auto* terrain=SurfaceEnvironment::ProductionTerrain();
    if(!terrain)return Reject(error,"LOLA photon exposure requires the initialized actual production terrain/BVH");
    // A frozen body rotates axes without dR/dt*x material motion. Reuse the
    // exact production policy, preserving source/collision frame conventions.
    if(!SurfaceEnvironment::ApplyFrozenTracingBodyTransform(Runtime::GetConfiguration(),forward.data(),error)||
        !SurfaceEnvironment::ApplyFrozenTracingBodyTransform(Runtime::GetConfiguration(),inverse.data(),error))return false;
    for(int r=0;r<3;++r)for(int c=0;c<3;++c){long double identity=0;for(int k=0;k<3;++k)identity+=forward[6*r+k]*inverse[6*k+c];
      if(std::abs(identity-(r==c?1:0))>128*std::numeric_limits<double>::epsilon())
        return Reject(error,"actual frozen terrain frame rotations are missing/inconsistent");}
    const auto& mesh=terrain->triangulation();f.terrainFaces.reserve(mesh.faces.size());
    for(const auto& indices:mesh.faces){T::Face face;for(int k=0;k<3;++k){const auto&v=mesh.verticesM[indices[k]];face.verticesM[k]=Rotate(forward.data(),{{v[0],v[1],v[2]}});}f.terrainFaces.push_back(face);}
    f.terrainVisibility=[terrain,inverse](const T::Vector&x,const T::Vector&sun,bool*visible,std::string*diagnostic){
      const auto origin=Rotate(inverse.data(),x),target=Rotate(inverse.data(),sun);T::Vector direction;
      for(int d=0;d<3;++d)direction[d]=target[d]-origin[d];
      const double range=std::hypot(direction[0],direction[1],direction[2]);
      if(!(range>0)||!std::isfinite(range))return Reject(diagnostic,"invalid actual volume-point to Sun ray");
      for(double&v:direction)v/=range;
      SurfaceEnvironment::RayHit hit;
      // Actual volume point: no surface-origin offset, local-horizon cosine
      // or excluded facet. Skipping a face here could hide a real gas shadow.
      if(!terrain->IntersectRay(origin,direction,0,range,std::numeric_limits<std::size_t>::max(),&hit,diagnostic))return false;
      *visible=!hit.hit;return true;
    };
  }
  return true;
}
}
bool FrozenProductionField::Create(const T::Controls& controls,bool earth,
    std::shared_ptr<const FrozenProductionField>*output,std::string*error){
  if(!output||!error||!Allowed(error))return false;
  auto candidate=std::make_shared<FrozenProductionField>();candidate->controls_=controls;
  auto& f=candidate->field_;f.frame=Frames::SolarOrbital;
  f.frozenGeometryEpochS=SurfaceEnvironment::ProductionSourceEpochS();f.referenceDistanceM=_AU_;
  for(int d=0;d<3;++d)f.sunPositionM[d]=Moon::xSun_SO[d];
  if(!SunValid(f.sunPositionM,error))return false;
  if(earth)f.opaqueSpheres.push_back({{{Moon::xEarth_SO[0],Moon::xEarth_SO[1],Moon::xEarth_SO[2]}},_RADIUS_(_EARTH_)});
  std::array<double,36> forward{},inverse{};
  for(int i=0;i<36;++i){forward[i]=Moon::OrbitalMotion::IAU_to_SO_TransformationMartix[i/6][i%6];inverse[i]=Moon::OrbitalMotion::SO_to_IAU_TransformationMartix[i/6][i%6];}
  if(!Surface(f,forward,inverse,error))return false;
  *output=std::move(candidate);error->clear();return true;
}
bool FrozenProductionField::CreateAtEpoch(const T::Controls& controls,bool earth,double et,
    std::shared_ptr<const FrozenProductionField>*output,std::string*error){
  if(!output||!error||!Allowed(error))return false;
  if(!std::isfinite(et))return Reject(error,"explicit frozen field requires finite TDB ephemeris seconds");
#if defined(_NO_SPICE_CALLS_)
  return Reject(error,"explicit normal frozen field requires furnished production SPICE kernels");
#else
  // A pre-buffer registration cannot borrow xSun_SO/xEarth_SO/step caches:
  // those may only be populated by the FIRST later amps_time_step. Query
  // actual furnished ephemerides and body transforms at ONE explicit ET;
  // do not advance global time or change a production cache to mimic startup.
  auto candidate=std::make_shared<FrozenProductionField>();candidate->controls_=controls;
  auto& f=candidate->field_;f.frame=Frames::SolarOrbital;f.frozenGeometryEpochS=et;f.referenceDistanceM=_AU_;
  SpiceDouble sun[6],earthState[6],light=0,forwardState[6][6],inverseState[6][6];
  spkezr_c("SUN",et,Frames::SolarOrbital,"NONE","MOON",sun,&light);
  if(earth)spkezr_c("EARTH",et,Frames::SolarOrbital,"NONE","MOON",earthState,&light);
  sxform_c(Frames::BodyFixed,Frames::SolarOrbital,et,forwardState);
  sxform_c(Frames::SolarOrbital,Frames::BodyFixed,et,inverseState);
  if(failed_c()){reset_c();return Reject(error,"furnished normal SPICE cannot provide explicit same-epoch SUN/EARTH/body field");}
  for(int d=0;d<3;++d)f.sunPositionM[d]=1000*sun[d];
  if(!SunValid(f.sunPositionM,error))return false;
  if(earth){T::Vector centre{{1000*earthState[0],1000*earthState[1],1000*earthState[2]}};
    for(double x:centre)if(!std::isfinite(x))return Reject(error,"nonfinite same-epoch Earth state");
    f.opaqueSpheres.push_back({centre,_RADIUS_(_EARTH_)});}
  std::array<double,36> forward{},inverse{};
  for(int i=0;i<36;++i){forward[i]=forwardState[i/6][i%6];inverse[i]=inverseState[i/6][i%6];}
  if(!Surface(f,forward,inverse,error))return false;
  *output=std::move(candidate);error->clear();return true;
#endif
}
bool FrozenProductionField::Build(const Controller::ExposureRequest& request,
    const std::vector<T::PhotonChannel>& photons,const std::vector<T::ConstantChannel>& constants,
    T::History*output,std::string*error)const{
  if(!std::isfinite(request.physicalEpochS)||!request.node||!request.node->block)
    return Reject(error,"trajectory exposure needs an actual physical particle epoch and allocated production node");
  return T::BuildPhotonHistory(request.arc,request.maximumS,field_,photons,constants,controls_,output,error);
}
}}}
