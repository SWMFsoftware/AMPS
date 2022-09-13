
#include <cmath>
#include <ctgmath>

#include "sep.h"
#include "quadrature.h"

namespace DxxInternalNumerics {
  thread_local double v;
  thread_local int spec;
  thread_local double FieldLineCoord;
  thread_local PIC::FieldLine::cFieldLineSegment *Segment;

  double Integrant(double *mu) {
    double D,dD_dmu,vParallel,vNorm;
    double t=1.0-mu[0]*mu[0];

    vParallel=v*mu[0];
    vNorm=v*sqrt(t); 

    SEP::Diffusion::GetPitchAngleDiffusionCoefficient(D,dD_dmu,mu[0],vParallel,vNorm,spec,FieldLineCoord,Segment);

    if (D==0.0) {
      //for debugging: catch the issue in the debugger by pacing a breat point in calculation of the D_mu_mu
      SEP::Diffusion::GetPitchAngleDiffusionCoefficient(D,dD_dmu,mu[0],vParallel,vNorm,spec,FieldLineCoord,Segment);
    }

    return t*t/D;
  }
} 

void SEP::Diffusion::GetDxx(double& D,double &dDxx_dx,double v,int spec,double FieldLineCoord,PIC::FieldLine::cFieldLineSegment *Segment,int iFieldLine) {
  namespace FL = PIC::FieldLine;
  double s,ds,D1,D0;

  DxxInternalNumerics::v=v;
  DxxInternalNumerics::spec=spec;
  DxxInternalNumerics::FieldLineCoord=FieldLineCoord;
  DxxInternalNumerics::Segment=Segment;

  double xmin[]={-1.0};
  double xmax[]={1.0};

  D=v*v/8.0*Quadrature::Gauss::Cube::GaussLegendre(1,3,DxxInternalNumerics::Integrant,xmin,xmax);  

  ds=Segment->GetLength();

  DxxInternalNumerics::FieldLineCoord=FL::FieldLinesAll[iFieldLine].move(FieldLineCoord,ds);
  DxxInternalNumerics::Segment=FL::FieldLinesAll[iFieldLine].GetSegment(DxxInternalNumerics::FieldLineCoord);

  if (DxxInternalNumerics::FieldLineCoord<0.0) {
    dDxx_dx=0.0;
    return;
  }
  else {
    D1=v*v/8.0*Quadrature::Gauss::Cube::GaussLegendre(1,3,DxxInternalNumerics::Integrant,xmin,xmax);
  }

  DxxInternalNumerics::FieldLineCoord=FL::FieldLinesAll[iFieldLine].move(FieldLineCoord,-ds);
  DxxInternalNumerics::Segment=FL::FieldLinesAll[iFieldLine].GetSegment(DxxInternalNumerics::FieldLineCoord);

  if (DxxInternalNumerics::FieldLineCoord<0.0) {
    dDxx_dx=0.0;
    return;
  }
  else {
    D0=v*v/8.0*Quadrature::Gauss::Cube::GaussLegendre(1,3,DxxInternalNumerics::Integrant,xmin,xmax);
  }

  dDxx_dx=(D1-D0)/(2.0*ds);
} 

//====================================================================================================
//calcualte particle's mean free path
double SEP::Diffusion::GetMeanFreePath(double v,int spec,double FieldLineCoord,PIC::FieldLine::cFieldLineSegment *Segment,int iFieldLine) {
  double D,dDxx_dx; 

  GetDxx(D,dDxx_dx,v,spec,FieldLineCoord,Segment,iFieldLine);
  return 3.0*D/v;
}
