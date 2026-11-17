//$Id$

/*
 * Surface_particle_interaction.cpp
 *
 *  Created on: Nov 29, 2016
 *      Author: vtenishe
 *
 * The gas-surface interaction models; their equations and sources are given in Surface.h.
 */

#include <cmath>
#include <cstdio>

#include "Surface.h"

//the parameter tables; NaN means "not set" (CheckParameters() stops the run when a selected model has an unset
//parameter). Set them at run time (not in static initializers of other files: they are initialized here).
double Surface::MaxwellReflection::AccommodationCoefficient[PIC::nTotalSpecies];
double Surface::CLL::TangentialMomentumAccommodationCoefficient[PIC::nTotalSpecies];
double Surface::CLL::NormalPartKineticEnergyAccommodationCoefficient[PIC::nTotalSpecies];
double Surface::QuasiSpecularReflection::CLL::TangentialMomentumAccommodationCoefficient[PIC::nTotalSpecies];
double Surface::QuasiSpecularReflection::CLL::NormalPartKineticEnergyAccommodationCoefficient[PIC::nTotalSpecies];
double Surface::QuasiSpecularReflection::Nocilla::DriftSpeed[PIC::nTotalSpecies];
double Surface::QuasiSpecularReflection::Nocilla::DriftAngle[PIC::nTotalSpecies];
double Surface::QuasiSpecularReflection::Nocilla::Temperature[PIC::nTotalSpecies];
double Surface::QuasiSpecularReflection::Schamberg::EnergyAccommodationCoefficient[PIC::nTotalSpecies];
double Surface::QuasiSpecularReflection::Schamberg::Nu[PIC::nTotalSpecies];
double Surface::QuasiSpecularReflection::Schamberg::BeamHalfAngle[PIC::nTotalSpecies];
int Surface::QuasiSpecularReflection::Nocilla::Mode[PIC::nTotalSpecies];   //Direct (0) unless set
double Surface::QuasiSpecularReflection::Nocilla::NormalMomentumAccommodation[PIC::nTotalSpecies];
double Surface::QuasiSpecularReflection::Nocilla::TangentialMomentumAccommodation[PIC::nTotalSpecies];
double Surface::QuasiSpecularReflection::Nocilla::EnergyAccommodation[PIC::nTotalSpecies];
long int Surface::QuasiSpecularReflection::Nocilla::nInconsistentReemissions=0;
long int Surface::QuasiSpecularReflection::Nocilla::nColdLimitReemissions=0;

namespace {
  struct cInitUnsetParameters {
    cInitUnsetParameters() {
      for (int s=0;s<PIC::nTotalSpecies;s++) {
        Surface::MaxwellReflection::AccommodationCoefficient[s]=NAN;
        Surface::CLL::TangentialMomentumAccommodationCoefficient[s]=NAN;
        Surface::CLL::NormalPartKineticEnergyAccommodationCoefficient[s]=NAN;
        Surface::QuasiSpecularReflection::CLL::TangentialMomentumAccommodationCoefficient[s]=NAN;
        Surface::QuasiSpecularReflection::CLL::NormalPartKineticEnergyAccommodationCoefficient[s]=NAN;
        Surface::QuasiSpecularReflection::Nocilla::DriftSpeed[s]=NAN;
        Surface::QuasiSpecularReflection::Nocilla::DriftAngle[s]=NAN;
        Surface::QuasiSpecularReflection::Nocilla::Temperature[s]=NAN;
        Surface::QuasiSpecularReflection::Schamberg::EnergyAccommodationCoefficient[s]=NAN;
        Surface::QuasiSpecularReflection::Schamberg::Nu[s]=NAN;
        Surface::QuasiSpecularReflection::Schamberg::BeamHalfAngle[s]=NAN;
        Surface::QuasiSpecularReflection::Nocilla::NormalMomentumAccommodation[s]=NAN;
        Surface::QuasiSpecularReflection::Nocilla::TangentialMomentumAccommodation[s]=NAN;
        Surface::QuasiSpecularReflection::Nocilla::EnergyAccommodation[s]=NAN;
      }
    }
  } InitUnsetParameters;

  //a uniform random number in (0,1) for logarithms (rnd() of the Mersenne twister can return 0)
  inline double RndOpen() {
    double r;

    do r=rnd(); while (r<=0.0);
    return r;
  }

  //stop the run when the parameters of the model of the species are not valid
  inline void RequireValidParameters(int spec,int model) {
    std::string error;

    if (Surface::CheckParameters(spec,model,&error)==false) {
      error="Error: "+error;
      exit(__LINE__,__FILE__,error.c_str());
    }
  }

  //set a parameter table for one species or for all of them (spec<0): assign(s,backup,restore) stores the new values
  //(and keeps the old ones in backup) or restores the old ones; the new values are checked as parameters of 'model'
  //(the selected model of the species is not changed)
  template <class F> bool SetForSpecies(int spec,int model,std::string *error,F assign) {
    if (spec>=PIC::nTotalSpecies) {
      if (error!=NULL) *error="species "+std::to_string(spec)+" does not exist";
      return false;
    }

    for (int s=0;s<PIC::nTotalSpecies;s++) if ((spec<0)||(s==spec)) {
      double backup[4];
      std::string reason;

      assign(s,backup,false);

      if (Surface::CheckParameters(s,model,&reason)==false) {
        assign(s,backup,true);
        if (error!=NULL) *error=reason;
        return false;
      }
    }

    return true;
  }

  inline double Dot(const double *a,const double *b) {return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];}
}


//particle reflection in the model of the specular reflection
void Surface::SpecularReflection::Reflect(double *vInit,CutCell::cTriangleFace *TriangleCutFace) {
  double c;

  c=vInit[0]*TriangleCutFace->ExternalNormal[0]+vInit[1]*TriangleCutFace->ExternalNormal[1]+vInit[2]*TriangleCutFace->ExternalNormal[2];
  vInit[0]-=2.0*c*TriangleCutFace->ExternalNormal[0];
  vInit[1]-=2.0*c*TriangleCutFace->ExternalNormal[1];
  vInit[2]-=2.0*c*TriangleCutFace->ExternalNormal[2];
}

int Surface::SpecularReflection::Processor(long int ptr,double* xInit,double* vInit,CutCell::cTriangleFace *TriangleCutFace,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* startNode) {
  Reflect(vInit,TriangleCutFace);
  return _PARTICLE_REJECTED_ON_THE_FACE_;
}

//adsorption: the particle is removed
int Surface::Adsorption::Processor(long int ptr,double* xInit,double* vInit,CutCell::cTriangleFace *TriangleCutFace,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* startNode) {
  return _PARTICLE_DELETED_ON_THE_FACE_;
}

//diffuse reflection: the flux-weighted half-range Maxwellian at the temperature 'temp' (Surface.h)
void Surface::DiffuseReflection::Reemit(double *vInit,double temp,int spec,CutCell::cTriangleFace *TriangleCutFace) {
  double *e0,*e1,*e2,ksi_n_r,ksi_t1_r,ksi_t2_r,beta,phi,c;

  if (!(temp>0.0)||std::isinf(temp)) {
    char msg[200];

    snprintf(msg,sizeof(msg),"Error: diffuse re-emission at the invalid temperature %e K",temp);
    exit(__LINE__,__FILE__,msg);
  }

  beta=sqrt(PIC::MolecularData::GetMass(spec)/(2.0*Kbol*temp));

  //get the internal frame of reference related to the surface element
  e0=TriangleCutFace->e0Orthogonal;
  e1=TriangleCutFace->e1Orthogonal;
  e2=TriangleCutFace->ExternalNormal;

  //the normal component of the reflected velocity: the density ~ u exp(-beta^2 u^2), u > 0
  ksi_n_r=sqrt(-log(rnd()))/beta;

  //tangenetial component of thereflected velocity: two normal components with the variance 1/(2 beta^2) (Box-Muller)
  phi=PiTimes2*rnd();
  c=sqrt(-log(rnd()))/beta;

  ksi_t1_r=c*cos(phi);
  ksi_t2_r=c*sin(phi);

  //combine the reflection velocity vector
  for (int idim=0;idim<3;idim++) vInit[idim]=ksi_t1_r*e0[idim]+ksi_t2_r*e1[idim]+ksi_n_r*e2[idim];
}

int Surface::DiffuseReflection::Processor(long int ptr,double* xInit,double* vInit,CutCell::cTriangleFace *TriangleCutFace,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* startNode) {
  int spec=PIC::ParticleBuffer::GetI(ptr);

  Reemit(vInit,Surface::GetSurfaceTemperature(TriangleCutFace,startNode),spec,TriangleCutFace);
  return _PARTICLE_REJECTED_ON_THE_FACE_;
}

//Maxwell model
bool Surface::MaxwellReflection::SetAccommodationCoefficient(int spec,double alpha,std::string *error) {
  return SetForSpecies(spec,Model::Maxwell,error,[&](int s,double *backup,bool restore) {
    if (restore) AccommodationCoefficient[s]=backup[0];
    else backup[0]=AccommodationCoefficient[s],AccommodationCoefficient[s]=alpha;
  });
}

int Surface::MaxwellReflection::Processor(long int ptr,double* xInit,double* vInit,CutCell::cTriangleFace *TriangleCutFace,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* startNode) {
  int res,spec;

  spec=PIC::ParticleBuffer::GetI(ptr);
  if (!((AccommodationCoefficient[spec]>=0.0)&&(AccommodationCoefficient[spec]<=1.0))) RequireValidParameters(spec,Model::Maxwell);

  if (rnd()<AccommodationCoefficient[spec]) res=Surface::DiffuseReflection::Processor(ptr,xInit,vInit,TriangleCutFace,startNode);
  else res=Surface::SpecularReflection::Processor(ptr,xInit,vInit,TriangleCutFace,startNode);

  return res;
}

//CLL model (Surface.h): Lord's sampling of the Cercignani-Lampis kernel
void Surface::CLL::Reflect(double *vInit,double alpha_n,double sigma_t,double Temp,int spec,CutCell::cTriangleFace *TriangleCutFace) {
  const double *e0=TriangleCutFace->e0Orthogonal,*e1=TriangleCutFace->e1Orthogonal,*n=TriangleCutFace->ExternalNormal;
  const double a=sqrt(2.0*Kbol*Temp/PIC::MolecularData::GetMass(spec));   //the most probable speed at the wall
  const double alpha_t=sigma_t*(2.0-sigma_t);   //the tangential energy accommodation coefficient
  double vn,vt[3],r,th,w_m,un,g0,g1;
  int idim;

  //the incident normal component (negative: toward the wall) and the incident tangential velocity vector
  vn=Dot(vInit,n);
  for (idim=0;idim<3;idim++) vt[idim]=vInit[idim]-vn*n[idim];

  //the normal component (Rice distribution): w = sqrt(r^2 + w_m^2 + 2 r w_m cos th), r = sqrt(-alpha_n ln R1),
  //th = 2 pi R2, w_m = sqrt(1-alpha_n) |v.n|/a
  w_m=sqrt(1.0-alpha_n)*fabs(vn)/a;

  if (alpha_n>0.0) {
    r=sqrt(-alpha_n*log(RndOpen()));
    th=PiTimes2*rnd();
    un=a*sqrt(r*r+w_m*w_m+2.0*r*w_m*cos(th));
  }
  else un=fabs(vn);   //alpha_n=0: the normal component is reflected unchanged

  //the tangential components: the mean (1-sigma_t) v_t, and an isotropic Gaussian in the tangent plane with the
  //variance alpha_t a^2/2 of each component: r' = sqrt(-alpha_t ln R3), th' = 2 pi R4
  if (alpha_t>0.0) {
    r=a*sqrt(-alpha_t*log(RndOpen()));
    th=PiTimes2*rnd();
    g0=r*cos(th),g1=r*sin(th);
  }
  else g0=0.0,g1=0.0;

  for (idim=0;idim<3;idim++) vInit[idim]=(1.0-sigma_t)*vt[idim]+g0*e0[idim]+g1*e1[idim]+un*n[idim];
}

bool Surface::CLL::SetCoefficients(int spec,double alpha_n,double sigma_t,std::string *error) {
  return SetForSpecies(spec,Model::CLL,error,[&](int s,double *backup,bool restore) {
    if (restore) {
      NormalPartKineticEnergyAccommodationCoefficient[s]=backup[0];
      TangentialMomentumAccommodationCoefficient[s]=backup[1];
    }
    else {
      backup[0]=NormalPartKineticEnergyAccommodationCoefficient[s],backup[1]=TangentialMomentumAccommodationCoefficient[s];
      NormalPartKineticEnergyAccommodationCoefficient[s]=alpha_n,TangentialMomentumAccommodationCoefficient[s]=sigma_t;
    }
  });
}

int Surface::CLL::Processor(long int ptr,double* xInit,double* vInit,CutCell::cTriangleFace *TriangleCutFace,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* startNode) {
  int spec=PIC::ParticleBuffer::GetI(ptr);
  double alpha_n=NormalPartKineticEnergyAccommodationCoefficient[spec],sigma_t=TangentialMomentumAccommodationCoefficient[spec];

  if (!((alpha_n>=0.0)&&(alpha_n<=1.0)&&(sigma_t>=0.0)&&(sigma_t<=2.0))) RequireValidParameters(spec,Model::CLL);

  Reflect(vInit,alpha_n,sigma_t,Surface::GetSurfaceTemperature(TriangleCutFace,startNode),spec,TriangleCutFace);
  return _PARTICLE_REJECTED_ON_THE_FACE_;
}

//quasi-specular CLL: the CLL kernel with its own coefficients
bool Surface::QuasiSpecularReflection::CLL::SetCoefficients(int spec,double alpha_n,double sigma_t,std::string *error) {
  return SetForSpecies(spec,Model::QuasiSpecularCLL,error,[&](int s,double *backup,bool restore) {
    if (restore) {
      NormalPartKineticEnergyAccommodationCoefficient[s]=backup[0];
      TangentialMomentumAccommodationCoefficient[s]=backup[1];
    }
    else {
      backup[0]=NormalPartKineticEnergyAccommodationCoefficient[s],backup[1]=TangentialMomentumAccommodationCoefficient[s];
      NormalPartKineticEnergyAccommodationCoefficient[s]=alpha_n,TangentialMomentumAccommodationCoefficient[s]=sigma_t;
    }
  });
}

int Surface::QuasiSpecularReflection::CLL::Processor(long int ptr,double* xInit,double* vInit,CutCell::cTriangleFace *TriangleCutFace,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* startNode) {
  int spec=PIC::ParticleBuffer::GetI(ptr);
  double alpha_n=NormalPartKineticEnergyAccommodationCoefficient[spec],sigma_t=TangentialMomentumAccommodationCoefficient[spec];

  if (!((alpha_n>=0.0)&&(alpha_n<=1.0)&&(sigma_t>=0.0)&&(sigma_t<=1.0))) RequireValidParameters(spec,Model::QuasiSpecularCLL);

  Surface::CLL::Reflect(vInit,alpha_n,sigma_t,Surface::GetSurfaceTemperature(TriangleCutFace,startNode),spec,TriangleCutFace);
  return _PARTICLE_REJECTED_ON_THE_FACE_;
}

int Surface::QuasiSpecularReflection::Processor(long int ptr,double* xInit,double* vInit,CutCell::cTriangleFace *TriangleCutFace,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* startNode) {
  return CLL::Processor(ptr,xInit,vInit,TriangleCutFace,startNode);
}

//Nocilla model (Surface.h)
//A sample of x > 0 with the density g(x) ~ x exp(-(x-s)^2), s >= 0 (the normal speed ratio of the molecules that cross
//a plane from a drifting Maxwellian).  With y = x - s, g = (y+s) exp(-y^2) on y > -s is the mixture of
//  P1: y exp(-y^2),        y >= 0,     weight 1/2                                  y = sqrt(-ln R)
//  P2: s exp(-y^2),        y >= 0,     weight s sqrt(pi)/2                         y = |N(0,1/2)|
//  P3: (y+s) exp(-y^2),   -s < y < 0,  weight s sqrt(pi)/2 erf(s) - (1-exp(-s^2))/2
//      (y+s = s sqrt(R), the density ~ (y+s) on (-s,0), accepted with the probability exp(-y^2))
//whose weights add up to the normalization (exp(-s^2) + s sqrt(pi) (1+erf(s)))/2 of g.
double Surface::QuasiSpecularReflection::Nocilla::SampleFluxNormalSpeedRatio(double s) {
  if (s<0.0) {
    //a drift toward the wall: the density ~ x exp(-x^2) exp(-2|s| x), x > 0, sampled by rejection, exactly:
    //|s| <= 1: x from ~ x exp(-x^2) (x = sqrt(-ln R)) accepted with the probability exp(-2|s| x);
    //|s| > 1:  x from ~ x exp(-2|s| x) (the gamma distribution of shape 2, x = -(ln R1 + ln R2)/(2|s|)) accepted with
    //          the probability exp(-x^2)
    const double a=-s;
    double x;

    if (a<=1.0) {
      do x=sqrt(-log(RndOpen())); while (rnd()>exp(-2.0*a*x));
    }
    else {
      do x=-(log(RndOpen())+log(RndOpen()))/(2.0*a); while (rnd()>exp(-x*x));
    }

    return x;
  }

  if (s==0.0) return sqrt(-log(RndOpen()));

  const double w1=0.5,w2=0.5*sqrt(Pi)*s,w3=w2*erf(s)-0.5*(1.0-exp(-s*s));
  const double r=rnd()*(w1+w2+w3);
  double y;

  if (r<w1) y=sqrt(-log(RndOpen()));
  else if (r<w1+w2) y=fabs(sqrt(-log(RndOpen()))*cos(PiTimes2*rnd()));
  else {
    do y=s*sqrt(rnd())-s; while (rnd()>exp(-y*y));
  }

  return y+s;
}

bool Surface::QuasiSpecularReflection::Nocilla::SetParameters(int spec,double U_r,double theta_r,double T_r,std::string *error) {
  return SetForSpecies(spec,Model::QuasiSpecularNocilla,error,[&](int s,double *backup,bool restore) {
    if (restore) DriftSpeed[s]=backup[0],DriftAngle[s]=backup[1],Temperature[s]=backup[2],Mode[s]=(int)backup[3];
    else {
      backup[0]=DriftSpeed[s],backup[1]=DriftAngle[s],backup[2]=Temperature[s],backup[3]=Mode[s];
      DriftSpeed[s]=U_r,DriftAngle[s]=theta_r,Temperature[s]=T_r,Mode[s]=Direct;
    }
  });
}

bool Surface::QuasiSpecularReflection::Nocilla::SetAccommodationCoefficients(int spec,double sigma_n,double sigma_t,double epsilon,std::string *error) {
  return SetForSpecies(spec,Model::QuasiSpecularNocilla,error,[&](int s,double *backup,bool restore) {
    if (restore) {
      NormalMomentumAccommodation[s]=backup[0],TangentialMomentumAccommodation[s]=backup[1],EnergyAccommodation[s]=backup[2];
      Mode[s]=(int)backup[3];
    }
    else {
      backup[0]=NormalMomentumAccommodation[s],backup[1]=TangentialMomentumAccommodation[s],backup[2]=EnergyAccommodation[s];
      backup[3]=Mode[s];
      NormalMomentumAccommodation[s]=sigma_n,TangentialMomentumAccommodation[s]=sigma_t,EnergyAccommodation[s]=epsilon;
      Mode[s]=AccommodationCoefficients;
    }
  });
}

//e^{x^2} erfc(x), x >= 0 (directly below 26; the asymptotic series above, relative error < 1e-13)
static double ScaledErfc(double x) {
  if (x<26.0) return exp(x*x)*erfc(x);

  double sum=1.0,term=1.0;

  for (int k=1;k<=5;k++) term*=-(2.0*k-1.0)/(2.0*x*x),sum+=term;
  return sum/(x*sqrt(Pi));
}

//<x> and <x^2> of the density ~ x exp(-(x-L)^2), x > 0 (Cook & Hoffbauer Eqs. 16, 17 with U sin(theta_u)=0): with
//D = e^{-L^2} + sqrt(pi) L (1+erf L), <x> = [L e^{-L^2} + sqrt(pi)(1/2+L^2)(1+erf L)]/D and
//<x^2> = [(1+L^2) e^{-L^2} + sqrt(pi)(3L/2+L^3)(1+erf L)]/D; for L < 0 all terms are multiplied by e^{L^2}
//(1+erf L = erfc(-L), sqrt(pi) e^{L^2} erfc(-L) = sqrt(pi) ScaledErfc(-L)) to avoid the underflow
void Surface::QuasiSpecularReflection::Nocilla::FluxNormalMoments(double L,double *mean,double *meanSquare) {
  double D,N1,N2;

  if (L>=0.0) {
    const double E=exp(-L*L),S=sqrt(Pi)*erfc(-L);

    D=E+L*S,N1=L*E+(0.5+L*L)*S,N2=(1.0+L*L)*E+(1.5*L+L*L*L)*S;
  }
  else {
    const double g=sqrt(Pi)*ScaledErfc(-L);

    D=1.0+L*g,N1=L+(0.5+L*L)*g,N2=(1.0+L*L)+(1.5*L+L*L*L)*g;
  }

  *mean=N1/D,*meanSquare=N2/D;
}

int Surface::QuasiSpecularReflection::Nocilla::SolveAccommodation(double vn,double vt,double sigma_n,double sigma_t,double epsilon,
    double Tw,double mass,double *Lambda,double *c_r,double *A,double *Q) {
  //vn = |v.n|, vt = |v_t| of the incident molecule
  const double vmean_w=sqrt(Pi*Kbol*Tw/(2.0*mass));   //the mean normal velocity of the diffuse re-emission
  const double V2=(1.0-epsilon)*(vn*vn+vt*vt)+epsilon*4.0*Kbol*Tw/mass;
  const double LambdaMin=-1000.0,LambdaMax=1.0E4;     //F(-1000) = 1.0e6, F(1e4) - 1 = 1.5e-8
  double mean,meanSquare,F;

  *A=(1.0-sigma_n)*vn+sigma_n*vmean_w;
  *Q=(1.0-sigma_t)*vt;

  const double B=V2-(*Q)*(*Q);

  if (!(*A>0.0)||!(B>(*A)*(*A))) {
    //no distribution with a positive temperature has these moments
    *Lambda=HUGE_VAL,*c_r=0.0;
    return 2;
  }

  const double r=B/((*A)*(*A));

  FluxNormalMoments(LambdaMax,&mean,&meanSquare);
  if (r<=(meanSquare+1.0)/(mean*mean)) {
    *Lambda=HUGE_VAL,*c_r=0.0;
    return 1;
  }

  //F(Lambda) = (<x^2>+1)/<x>^2 decreases monotonically: bisection
  double lo=LambdaMin,hi=LambdaMax;

  FluxNormalMoments(lo,&mean,&meanSquare);
  if (r>=(meanSquare+1.0)/(mean*mean)) hi=lo;   //beyond the range (V2-Q^2 > 1e6 A^2): the limit LambdaMin

  for (int i=0;(i<200)&&(hi-lo>1.0E-13*(1.0+fabs(lo)+fabs(hi)));i++) {
    double mid=0.5*(lo+hi);

    FluxNormalMoments(mid,&mean,&meanSquare);
    F=(meanSquare+1.0)/(mean*mean);

    if (F>r) lo=mid;
    else hi=mid;
  }

  *Lambda=0.5*(lo+hi);
  FluxNormalMoments(*Lambda,&mean,&meanSquare);
  *c_r=(*A)/mean;
  return 0;
}

void Surface::QuasiSpecularReflection::Nocilla::ReemitAccommodation(double *vInit,double sigma_n,double sigma_t,double epsilon,double Tw,int spec,CutCell::cTriangleFace *TriangleCutFace) {
  const double *e0=TriangleCutFace->e0Orthogonal,*e1=TriangleCutFace->e1Orthogonal,*n=TriangleCutFace->ExternalNormal;
  double vn,t[3],b[3],l,x,y1,y2,r,th,Lambda,c_r,A,Q;
  int idim;

  //the unit tangential direction t of the incident velocity (the plane of incidence) and b = n x t
  vn=Dot(vInit,n);
  for (idim=0;idim<3;idim++) t[idim]=vInit[idim]-vn*n[idim];
  l=sqrt(Dot(t,t));

  if (l>0.0) for (idim=0;idim<3;idim++) t[idim]/=l;
  else {
    double psi=PiTimes2*rnd();
    for (idim=0;idim<3;idim++) t[idim]=cos(psi)*e0[idim]+sin(psi)*e1[idim];
  }

  b[0]=n[1]*t[2]-n[2]*t[1];
  b[1]=n[2]*t[0]-n[0]*t[2];
  b[2]=n[0]*t[1]-n[1]*t[0];

  int code=SolveAccommodation(fabs(vn),l,sigma_n,sigma_t,epsilon,Tw,PIC::MolecularData::GetMass(spec),&Lambda,&c_r,&A,&Q);

  if (code!=0) {
    //the cold limit: the mean velocity
    if (code==1) {
      #if _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
      #pragma omp atomic
      #endif
      nColdLimitReemissions++;
    }
    else {
      #if _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
      #pragma omp atomic
      #endif
      nInconsistentReemissions++;
    }

    for (idim=0;idim<3;idim++) vInit[idim]=A*n[idim]+Q*t[idim];
    return;
  }

  x=SampleFluxNormalSpeedRatio(Lambda);
  r=sqrt(-log(RndOpen()));
  th=PiTimes2*rnd();
  y1=r*cos(th),y2=r*sin(th);

  for (idim=0;idim<3;idim++) vInit[idim]=c_r*(x*n[idim]+y1*t[idim]+y2*b[idim])+Q*t[idim];
}

void Surface::QuasiSpecularReflection::Nocilla::Reemit(double *vInit,double U_r,double theta_r,double T_r,int spec,CutCell::cTriangleFace *TriangleCutFace) {
  const double *e0=TriangleCutFace->e0Orthogonal,*e1=TriangleCutFace->e1Orthogonal,*n=TriangleCutFace->ExternalNormal;
  const double c_r=sqrt(2.0*Kbol*T_r/PIC::MolecularData::GetMass(spec));
  double vn,t[3],b[3],l,x,y1,y2,r,th;
  int idim;

  //the unit tangential direction t of the incident velocity (the plane of incidence) and b = n x t
  vn=Dot(vInit,n);
  for (idim=0;idim<3;idim++) t[idim]=vInit[idim]-vn*n[idim];
  l=sqrt(Dot(t,t));

  if (l>0.0) for (idim=0;idim<3;idim++) t[idim]/=l;
  else {
    //normal incidence: the plane of incidence is not defined; a random direction keeps the re-emission axisymmetric
    double psi=PiTimes2*rnd();
    for (idim=0;idim<3;idim++) t[idim]=cos(psi)*e0[idim]+sin(psi)*e1[idim];
  }

  b[0]=n[1]*t[2]-n[2]*t[1];
  b[1]=n[2]*t[0]-n[0]*t[2];
  b[2]=n[0]*t[1]-n[1]*t[0];

  //the normal speed ratio x ~ x exp(-(x-s)^2), s = U_r cos(theta_r)/c_r, and the two tangential thermal components
  x=SampleFluxNormalSpeedRatio(U_r*cos(theta_r)/c_r);
  r=sqrt(-log(RndOpen()));
  th=PiTimes2*rnd();
  y1=r*cos(th),y2=r*sin(th);

  for (idim=0;idim<3;idim++) vInit[idim]=c_r*(x*n[idim]+y1*t[idim]+y2*b[idim])+U_r*sin(theta_r)*t[idim];
}

int Surface::QuasiSpecularReflection::Nocilla::Processor(long int ptr,double* xInit,double* vInit,CutCell::cTriangleFace *TriangleCutFace,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* startNode) {
  int spec=PIC::ParticleBuffer::GetI(ptr);

  if (CheckParameters(spec,Model::QuasiSpecularNocilla)==false) RequireValidParameters(spec,Model::QuasiSpecularNocilla);

  if (Mode[spec]==AccommodationCoefficients) {
    ReemitAccommodation(vInit,NormalMomentumAccommodation[spec],TangentialMomentumAccommodation[spec],EnergyAccommodation[spec],
        Surface::GetSurfaceTemperature(TriangleCutFace,startNode),spec,TriangleCutFace);
  }
  else Reemit(vInit,DriftSpeed[spec],DriftAngle[spec],Temperature[spec],spec,TriangleCutFace);

  return _PARTICLE_REJECTED_ON_THE_FACE_;
}

//Schamberg model (Surface.h)
double Surface::QuasiSpecularReflection::Schamberg::BeamShapeFunction(double phi0) {
  //Phi(phi0) = (1-x^2)/(1-4x^2) (sin(2 phi0)/2 - x)/(sin(phi0) - x), x = 2 phi0/pi: the mean cosine of the angle to the
  //axis for the flux ~ cos(pi phi/(2 phi0)) per unit solid angle, phi <= phi0 (with k = pi/(2 phi0):
  //int_0^phi0 cos(k phi) sin(phi) cos(phi) dphi = (1 - (k/2) sin(2 phi0))/(4-k^2),
  //int_0^phi0 cos(k phi) sin(phi) dphi = (1 - k sin(phi0))/(1-k^2)).  The limits phi0 -> 0 (1) and phi0 -> pi/2 (2/3),
  //and x = 1/2, are removable singularities: there the ratio is evaluated by quadrature.
  if (phi0<=0.0) return 1.0;

  const double x=2.0*phi0/Pi;

  if ((phi0>1.0E-3)&&(fabs(x-0.5)>1.0E-3)&&(fabs(x-1.0)>1.0E-3)) {
    return (1.0-x*x)/(1.0-4.0*x*x)*(0.5*sin(2.0*phi0)-x)/(sin(phi0)-x);
  }

  //Gauss-Legendre quadrature of the two integrals (20 points on [0,phi0], the integrands are smooth)
  static const double gx[10]={0.0765265211334973,0.2277858511416451,0.3737060887154195,0.5108670019508271,0.6360536807265150,
                              0.7463319064601508,0.8391169718222188,0.9122344282513259,0.9639719272779138,0.9931285991850949};
  static const double gw[10]={0.1527533871307258,0.1491729864726037,0.1420961093183820,0.1316886384491766,0.1181945319615184,
                              0.1019301198172404,0.0832767415767048,0.0626720483341091,0.0406014298003869,0.0176140071391521};
  double num=0.0,den=0.0;

  for (int i=0;i<10;i++) for (int sgn=-1;sgn<=1;sgn+=2) {
    double phi=0.5*phi0*(1.0+sgn*gx[i]),w=cos(0.5*Pi*phi/phi0)*sin(phi)*gw[i];

    num+=w*cos(phi),den+=w;
  }

  return num/den;
}

bool Surface::QuasiSpecularReflection::Schamberg::SetParameters(int spec,double alpha,double nu,double phi0,std::string *error) {
  return SetForSpecies(spec,Model::QuasiSpecularSchamberg,error,[&](int s,double *backup,bool restore) {
    if (restore) EnergyAccommodationCoefficient[s]=backup[0],Nu[s]=backup[1],BeamHalfAngle[s]=backup[2];
    else {
      backup[0]=EnergyAccommodationCoefficient[s],backup[1]=Nu[s],backup[2]=BeamHalfAngle[s];
      EnergyAccommodationCoefficient[s]=alpha,Nu[s]=nu,BeamHalfAngle[s]=phi0;
    }
  });
}

void Surface::QuasiSpecularReflection::Schamberg::Reflect(double *vInit,double alpha,double nu,double phi0,double Temp,int spec,CutCell::cTriangleFace *TriangleCutFace) {
  const double *e0=TriangleCutFace->e0Orthogonal,*n=TriangleCutFace->ExternalNormal;
  const double mass=PIC::MolecularData::GetMass(spec);
  double vn,t[3],b[3],axis[3],b1[3],d[3],lt,Vi2,Vr,sin_i,sin_r,cos_r;
  int idim;

  //the incident speed, the angle of incidence (from the normal) and the tangential direction t (plane of incidence)
  vn=Dot(vInit,n);
  for (idim=0;idim<3;idim++) t[idim]=vInit[idim]-vn*n[idim];
  lt=sqrt(Dot(t,t));
  Vi2=Dot(vInit,vInit);

  if (lt>0.0) for (idim=0;idim<3;idim++) t[idim]/=lt;
  else for (idim=0;idim<3;idim++) t[idim]=e0[idim];   //normal incidence: the axis is the normal (sin(theta_r)=0)

  sin_i=(Vi2>0.0) ? lt/sqrt(Vi2) : 0.0;
  if (sin_i>1.0) sin_i=1.0;

  //the re-emitted speed: V_r^2 = (1-alpha) V_i^2 + alpha 4 k Tw/m (the wall temperature is used only when alpha > 0)
  Vr=sqrt((1.0-alpha)*Vi2+((alpha>0.0) ? alpha*4.0*Kbol*Temp/mass : 0.0));

  //the beam axis: sin(theta_r) = sin(theta_i)^nu, in the plane of incidence, forward
  sin_r=pow(sin_i,nu);
  cos_r=sqrt(1.0-sin_r*sin_r);

  for (idim=0;idim<3;idim++) {
    axis[idim]=cos_r*n[idim]+sin_r*t[idim];
    b1[idim]=-sin_r*n[idim]+cos_r*t[idim];  //perpendicular to the axis in the plane of incidence
  }

  b[0]=n[1]*t[2]-n[2]*t[1];
  b[1]=n[2]*t[0]-n[0]*t[2];
  b[2]=n[0]*t[1]-n[1]*t[0];

  if (phi0<=0.0) {
    for (idim=0;idim<3;idim++) d[idim]=axis[idim];
  }
  else {
    //the direction in the cone: cos(phi) uniform in [cos(phi0),1] (uniform per solid angle) accepted with the
    //probability cos(pi phi/(2 phi0)), the azimuth psi uniform; a direction below the surface is sampled again
    double cos_phi,sin_phi,phi,psi;

    do {
      do {
        cos_phi=1.0-rnd()*(1.0-cos(phi0));
        phi=acos(cos_phi);
      }
      while (rnd()>cos(0.5*Pi*phi/phi0));

      sin_phi=sqrt(1.0-cos_phi*cos_phi);
      psi=PiTimes2*rnd();

      for (idim=0;idim<3;idim++) d[idim]=cos_phi*axis[idim]+sin_phi*(cos(psi)*b1[idim]+sin(psi)*b[idim]);
    }
    while (Dot(d,n)<=0.0);
  }

  for (idim=0;idim<3;idim++) vInit[idim]=Vr*d[idim];
}

int Surface::QuasiSpecularReflection::Schamberg::Processor(long int ptr,double* xInit,double* vInit,CutCell::cTriangleFace *TriangleCutFace,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* startNode) {
  int spec=PIC::ParticleBuffer::GetI(ptr);

  if (CheckParameters(spec,Model::QuasiSpecularSchamberg)==false) RequireValidParameters(spec,Model::QuasiSpecularSchamberg);

  const double alpha=EnergyAccommodationCoefficient[spec];
  Reflect(vInit,alpha,Nu[spec],BeamHalfAngle[spec],(alpha>0.0) ? Surface::GetSurfaceTemperature(TriangleCutFace,startNode) : 0.0,spec,TriangleCutFace);
  return _PARTICLE_REJECTED_ON_THE_FACE_;
}


