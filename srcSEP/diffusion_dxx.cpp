
#include <algorithm>
#include <cmath>
#include <ctgmath>
#include <limits>

#include "sep.h"
#include "util/sep_coefficient_physics.h"
#include "util/sep_coefficient_registry.h"

//static variables from c_D_x_x
//template<class T> double SEP::Diffusion::cD_x_x<T>::speed;

/*template<class T>
double SEP::Diffusion::cD_x_x<T>::p=0.0;

template<class T>
double SEP::Diffusion::cD_x_x<T>::W[2]={0.0,0.0};

template<class T>
double SEP::Diffusion::cD_x_x<T>::AbsB=0.0;

template<class T>
double SEP::Diffusion::cD_x_x<T>::xLocation[3]={0.0,0.0,0.0};

template<class T>
double SEP::Diffusion::cD_x_x<T>::vAlfven=0.0;

template<class T>
double SEP::Diffusion::cD_x_x<T>::B[3]={0.0,0.0,0.0};
    
template<class T>
PIC::FieldLine::cFieldLineSegment*  SEP::Diffusion::cD_x_x<T>::Segment=NULL;*/

namespace DxxInternalNumerics {
  SEP::Transport::CoefficientPhysics::SpatialDiffusionResult Evaluate(
      double speedMPerS, int spec, double fieldLineCoordinate,
      PIC::FieldLine::cFieldLineSegment* segment) {
    namespace CP=SEP::Transport::CoefficientPhysics;
    CP::SpatialDiffusionResult result;
    if (SEP::Diffusion::GetPitchAngleDiffusionCoefficient==NULL) {
      result.status=SEP::Transport::Status::Ok();
      result.valueState=CP::ValueState::Finite;
      result.kappaParallelM2PerS=0.0;
      return result;
    }
    CP::SpatialQuadratureConfiguration configuration=
        SEP::Transport::Coefficient::ActiveConfiguration().spatialQuadrature;
    configuration.gapPolicy=
        SEP::Transport::Coefficient::ActiveConfiguration().resonanceGapPolicy;
    return CP::IntegrateSpatialDiffusion(
        speedMPerS,
        [=](double mu) {
          CP::PitchAngleResult pitch;
          const double shape=std::max(0.0,1.0-mu*mu);
          double derivative=0.0;
          SEP::Diffusion::GetPitchAngleDiffusionCoefficient(
              pitch.dMuMuPerS,derivative,mu,speedMPerS*mu,
              speedMPerS*sqrt(shape),spec,fieldLineCoordinate,segment);
          pitch.dDmuMuDmuPerS=derivative;
          pitch.valueState=CP::ValueState::Finite;
          pitch.status=std::isfinite(pitch.dMuMuPerS) &&
              pitch.dMuMuPerS>=0.0
              ? SEP::Transport::Status::Ok()
              : SEP::Transport::Status::Error(
                    SEP::Transport::StatusCode::InvalidCoefficient,
                    "legacy Dmumu callback returned invalid quadrature input");
          return pitch;
        },configuration);
  }
} 

void SEP::Diffusion::GetDxx(double& D,double &dDxx_dx,double v,int spec,double FieldLineCoord,PIC::FieldLine::cFieldLineSegment *Segment,int iFieldLine) {
  namespace FL = PIC::FieldLine;
  namespace CP = SEP::Transport::CoefficientPhysics;
  const CP::SpatialDiffusionResult center=
      DxxInternalNumerics::Evaluate(v,spec,FieldLineCoord,Segment);
  if (!center.status.ok()) {
    D=std::numeric_limits<double>::quiet_NaN();
    dDxx_dx=std::numeric_limits<double>::quiet_NaN();
    return;
  }
  D=center.kappaParallelM2PerS;
  if (center.valueState==CP::ValueState::Ballistic) {
    // The legacy signature cannot carry the typed state.  Positive infinity is
    // retained only as a compatibility serialization; production adapters use
    // SpatialDiffusionResult and reject infinite kappa before Parker stepping.
    dDxx_dx=0.0;
    return;
  }

  double h=SEP::Transport::ActiveNumericalTolerances().geometryFraction*
      Segment->GetLength();
  bool derivativeFound=false;
  for (int refinement=0;refinement<10;refinement++) {
    PIC::FieldLine::cFieldLineSegment* plusSegment=Segment;
    PIC::FieldLine::cFieldLineSegment* minusSegment=Segment;
    const double plusCoordinate=FL::FieldLinesAll[iFieldLine].move(
        FieldLineCoord,h,plusSegment);
    const double minusCoordinate=FL::FieldLinesAll[iFieldLine].move(
        FieldLineCoord,-h,minusSegment);
    plusSegment=FL::FieldLinesAll[iFieldLine].GetSegment(plusCoordinate);
    minusSegment=FL::FieldLinesAll[iFieldLine].GetSegment(minusCoordinate);
    const CP::SpatialDiffusionResult plus=plusSegment ?
        DxxInternalNumerics::Evaluate(v,spec,plusCoordinate,plusSegment) :
        CP::SpatialDiffusionResult();
    const CP::SpatialDiffusionResult minus=minusSegment ?
        DxxInternalNumerics::Evaluate(v,spec,minusCoordinate,minusSegment) :
        CP::SpatialDiffusionResult();
    if (plus.status.ok() && minus.status.ok() &&
        plus.valueState==CP::ValueState::Finite &&
        minus.valueState==CP::ValueState::Finite) {
      dDxx_dx=(plus.kappaParallelM2PerS-minus.kappaParallelM2PerS)/(2.0*h);
      derivativeFound=std::isfinite(dDxx_dx);
      if (derivativeFound) break;
    }
    else if (plus.status.ok() && plus.valueState==CP::ValueState::Finite) {
      dDxx_dx=(plus.kappaParallelM2PerS-D)/h;
      derivativeFound=std::isfinite(dDxx_dx);
      if (derivativeFound) break;
    }
    else if (minus.status.ok() && minus.valueState==CP::ValueState::Finite) {
      dDxx_dx=(D-minus.kappaParallelM2PerS)/h;
      derivativeFound=std::isfinite(dDxx_dx);
      if (derivativeFound) break;
    }
    h*=0.5;
  }
  if (!derivativeFound) dDxx_dx=std::numeric_limits<double>::quiet_NaN();
} 

//====================================================================================================
//calcualte particle's mean free path
double SEP::Diffusion::GetMeanFreePath(double v,int spec,double FieldLineCoord,PIC::FieldLine::cFieldLineSegment *Segment,int iFieldLine) {
  double D,dDxx_dx; 

  GetDxx(D,dDxx_dx,v,spec,FieldLineCoord,Segment,iFieldLine);
  return 3.0*D/v;
}
