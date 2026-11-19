// Numerical fixtures in convenient SI units, never lunar observations. The
// production integrator/clock are linked unchanged; Python supplies separate
// closed-form inverse-square and shadow-duration oracles.
#include "../../species/TrajectoryExposure.h"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
namespace T=Moon::Species::TrajectoryExposure;
namespace P=Moon::Species::PiecewiseHazard;
int main(int argc,char**argv){
  const std::string mode=argc>1?argv[1]:"linear";
  T::Arc a;T::Field f;f.frame="ANALYTICAL-CARTESIAN";f.sunPositionM={{10,0,0}};f.referenceDistanceM=10;
  T::Controls c;c.rootTimeToleranceS=1e-12;c.absoluteDepthTolerance=1e-10L;c.relativeDepthTolerance=1e-8L;c.maximumSegments=100000;c.maximumGeometryIntervals=1000000;
  std::vector<T::PhotonChannel> photons{{7,2}};std::vector<T::ConstantChannel> constants{{9,.3}};
  if(mode=="linear"||mode=="resource")a.velocityMS[0]=1;
  else if(mode=="accelerated")a.accelerationMS2[0]=1;
  else if(mode=="sphere"||mode=="thin-sphere"){
    f.sunPositionM={{100,0,0}};f.referenceDistanceM=100;
    const double radius=mode=="sphere"?1:.001;
    f.opaqueSpheres.push_back({{{50,0,0}},radius});a.positionM={{0,-3,0}};a.velocityMS={{0,6,0}};
  }
  else if(mode=="terrain"||mode=="callback-error"||mode=="far-sun"){
    T::Face p{{{{5,-.1,-1}},{{5,.1,-1}},{{5,.1,1}}}},q{{{{5,-.1,-1}},{{5,.1,1}},{{5,-.1,1}}}};
    f.terrainFaces={p,q};a.positionM={{0,-1,0}};a.velocityMS={{0,2,0}};
    if(mode=="far-sun"){f.sunPositionM={{1e11,1e11,1e11}};f.referenceDistanceM=std::sqrt(3.)*1e11;
      a.positionM={{0,-6,-5}};photons[0].rateAtReferenceS=1;constants.clear();}
    f.terrainVisibility=[&](const T::Vector&x,const T::Vector&sun,bool*visible,std::string*error){
      if(mode=="callback-error"){*error="declared geometry dependency unavailable";return false;}
      const long double s=(5.L-x[0])/(static_cast<long double>(sun[0])-x[0]),y=x[1]+s*(static_cast<long double>(sun[1])-x[1]),z=x[2]+s*(static_cast<long double>(sun[2])-x[2]);
      *visible=!(s>=0&&s<=1&&std::abs(y)<=.1&&std::abs(z)<=1);return true;
    };
  }
  else if(mode=="constants"){photons.clear();}
  else if(mode=="close-roots"){
    f.sunPositionM={{100,0,0}};f.referenceDistanceM=100;f.opaqueSpheres.push_back({{{50,0,0}},4.5e-7});
    a.positionM={{0,-.003,0}};a.velocityMS={{0,.006,0}};photons[0].rateAtReferenceS=1;constants.clear();
    c.rootTimeToleranceS=1e-4;c.absoluteDepthTolerance=2.5e-4;c.relativeDepthTolerance=0;
  }
  else if(mode=="position-projection"||mode=="position-projection-budget"){
    // A sub-ULP flight at a translated observer coordinate is representable
    // as an exact quadratic arc, but all callback Vector positions round to
    // the same binary64 point. The true line reaches the screen edge at
    // t=2*halfwidth/velocity=0.2 s; its callback projection stays dark.
    const double halfwidth=1e-12;
    T::Face p{{{{5,-halfwidth,-1}},{{5,halfwidth,-1}},{{5,halfwidth,1}}}},q{{{{5,-halfwidth,-1}},{{5,halfwidth,1}},{{5,-halfwidth,1}}}};
    f.terrainFaces={p,q};f.sunPositionM={{10,-1e6,0}};f.referenceDistanceM=2e6;
    a.positionM={{0,1e6,0}};a.velocityMS={{0,1e-11,0}};photons[0].rateAtReferenceS=1;constants.clear();
    c.rootTimeToleranceS=1e-3;c.absoluteDepthTolerance=mode=="position-projection"?2:.2L;c.relativeDepthTolerance=0;
    f.terrainVisibility=[halfwidth](const T::Vector&x,const T::Vector&sun,bool*visible,std::string*){
      const long double u=(5.L-x[0])/(static_cast<long double>(sun[0])-x[0]);
      const long double y=x[1]+u*(static_cast<long double>(sun[1])-x[1]);
      const long double z=x[2]+u*(static_cast<long double>(sun[2])-x[2]);
      *visible=!(u>=0&&u<=1&&std::abs(y)<=halfwidth&&std::abs(z)<=1);return true;};
  }
  else if(mode=="negative")photons[0].rateAtReferenceS=-1;
  else if(mode=="duplicate")constants[0].reaction=7;
  else if(mode=="nan")a.positionM[0]=std::numeric_limits<double>::quiet_NaN();
  else if(mode=="no-frame")f.frame.clear();
  else if(mode=="missing-callback")f.terrainFaces.push_back({});
  else if(mode=="root-budget"){
    f.opaqueSpheres.push_back({{{5,0,0}},1});a.positionM={{0,-4,0}};a.velocityMS={{0,8,0}};
    c.rootTimeToleranceS=.01;c.absoluteDepthTolerance=1e-20L;c.relativeDepthTolerance=0;
  }
  else if(mode!="constant")return 2;
  if(mode=="resource")c.maximumSegments=1;
  T::History history;history.interpolationDepthErrorBound=-123;std::string error;
  const bool ok=T::BuildPhotonHistory(a,1,f,photons,constants,c,&history,&error);
  std::cout<<std::setprecision(21)<<"{\"evaluated\":"<<(ok?"true":"false")<<",\"unchanged\":"<<(history.interpolationDepthErrorBound==-123?"true":"false");
  if(ok){P::Exposure value;const auto status=P::Integrate(history.segments,0,1,&value,&error);if(status!=P::Status::Complete)return 2;
    std::cout<<",\"depth\":"<<value.integratedHazard<<",\"bound\":"<<history.interpolationDepthErrorBound<<",\"shadow_bound\":"<<history.shadowBoundaryDepthErrorBound<<",\"segments\":"<<history.segments.size()<<",\"roots\":"<<history.silhouetteRoots;
    if(mode=="far-sun"){
      const long double visible=1-.1L/(1-5.L/1e11L);
      auto factor=[&](long double t){long double q2=0;for(int d=0;d<3;++d){const long double r=static_cast<long double>(f.sunPositionM[d])-a.positionM[d]-t*a.velocityMS[d];q2+=r*r;}return static_cast<long double>(f.referenceDistanceM)*f.referenceDistanceM/q2;};
      std::cout<<",\"independent_lower\":"<<visible*factor(0)<<",\"independent_upper\":"<<visible*factor(1);
    }
  }
  std::cout<<"}\n";std::cerr<<error<<'\n';return 0;
}
