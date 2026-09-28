#include "sep_coronal_cme/pfss_harmonics.h"

#include "sep_coronal_cme/constants.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <utility>

namespace SEP { namespace CoronalCME { namespace {

struct Basis {
  double y=0.0, theta=0.0, theta2=0.0, phi=0.0, thetaPhi=0.0, phi2=0.0;
};

double FactorialRatio(int lower, int upper) {
  // Computes lower!/upper! without factorial overflow.
  double result = 1.0;
  for (int value=lower+1; value<=upper; ++value) result /= value;
  return result;
}

double AssociatedLegendre(int l, int m, double x) {
  double pmm=1.0;
  if (m>0) {
    const double root=std::sqrt(std::max(0.0,1.0-x*x));
    double factor=1.0;
    for (int k=1;k<=m;++k) { pmm *= -factor*root; factor += 2.0; }
  }
  if (l==m) return pmm;
  double pmmp1=x*(2*m+1)*pmm;
  if (l==m+1) return pmmp1;
  double previous=pmm, current=pmmp1;
  for (int n=m+2;n<=l;++n) {
    const double next=((2*n-1)*x*current-(n+m-1)*previous)/(n-m);
    previous=current; current=next;
  }
  return current;
}

Basis RealBasis(int l, int m, bool sine, double theta, double phi) {
  const double x=std::cos(theta), sinTheta=std::sin(theta);
  const double p=AssociatedLegendre(l,m,x);
  const double previous=l>m ? AssociatedLegendre(l-1,m,x) : 0.0;
  const double norm=std::sqrt((2.0*l+1.0)/(4.0*Constants::kPi)*
      (m==0 ? 1.0 : 2.0)*FactorialRatio(l-m,l+m));
  // PFSS evaluation at exact poles is obtained as a limiting value. A small
  // denominator is used only in the derivative recurrence; the physical
  // caller normally avoids storing a singular spherical basis at a pole.
  const double safeSin=std::abs(sinTheta)>1.0e-12 ? sinTheta :
      (sinTheta>=0.0 ? 1.0e-12 : -1.0e-12);
  const double pTheta=(l*x*p-(l+m)*previous)/safeSin;
  const double pTheta2=-x/safeSin*pTheta-
      (l*(l+1.0)-m*m/(safeSin*safeSin))*p;
  const double angle=m*phi;
  const double trig=sine ? std::sin(angle) : std::cos(angle);
  const double trigPhi=m*(sine ? std::cos(angle) : -std::sin(angle));
  const double trigPhi2=-m*m*trig;
  return {norm*p*trig,norm*pTheta*trig,norm*pTheta2*trig,
          norm*p*trigPhi,norm*pTheta*trigPhi,norm*p*trigPhi2};
}

Vec3 SphericalToCartesian(double theta,double phi,double br,double bt,double bp) {
  const double st=std::sin(theta),ct=std::cos(theta),sp=std::sin(phi),cp=std::cos(phi);
  return {br*st*cp+bt*ct*cp-bp*sp,
          br*st*sp+bt*ct*sp+bp*cp,br*ct-bt*st};
}

bool Finite(double value) { return std::isfinite(value); }

}  // namespace

Core::Result<PfssHarmonics> PfssHarmonics::Create(double rs,double rb,
    std::vector<HarmonicCoefficient> coefficients, MonopolePolicy policy,
    double tolerance) {
  if (!(Finite(rs)&&Finite(rb)&&rs>0.0&&rb>rs&&tolerance>=0.0))
    return Core::Result<PfssHarmonics>::Failure(Core::StatusCode::InvalidConfiguration,
        "PFSS requires finite 0<R_sun<R_b and a nonnegative monopole tolerance");
  std::map<std::pair<int,int>,bool> seen;
  std::vector<HarmonicCoefficient> selected;
  for (const auto& coefficient:coefficients) {
    if (coefficient.degree<0||coefficient.order<0||coefficient.order>coefficient.degree||
        !Finite(coefficient.cosineT)||!Finite(coefficient.sineT)||
        (coefficient.order==0&&coefficient.sineT!=0.0))
      return Core::Result<PfssHarmonics>::Failure(Core::StatusCode::InvalidConfiguration,
          "invalid real spherical-harmonic coefficient");
    if (!seen.emplace(std::make_pair(coefficient.degree,coefficient.order),true).second)
      return Core::Result<PfssHarmonics>::Failure(Core::StatusCode::InvalidConfiguration,
          "duplicate harmonic degree/order");
    if (coefficient.degree==0) {
      if (std::abs(coefficient.cosineT)>tolerance&&policy==MonopolePolicy::Reject)
        return Core::Result<PfssHarmonics>::Failure(Core::StatusCode::InvalidConfiguration,
            "PFSS monopole exceeds the declared tolerance");
      continue; // removal is explicit and the solve never retains l=0
    }
    selected.push_back(coefficient);
  }
  PfssHarmonics result; result.solarRadiusM_=rs; result.sourceSurfaceRadiusM_=rb;
  result.coefficients_=std::move(selected);
  return Core::Result<PfssHarmonics>::Success(std::move(result));
}

Core::Result<SphericalField> PfssHarmonics::Evaluate(double r,double theta,double phi) const {
  if (!(Finite(r)&&Finite(theta)&&Finite(phi)&&r>=solarRadiusM_&&
        r<=sourceSurfaceRadiusM_&&theta>=0.0&&theta<=Constants::kPi))
    return Core::Result<SphericalField>::Failure(Core::StatusCode::OutOfDomain,
        "PFSS query is outside its spherical shell");
  SphericalField field;
  const double sinTheta=std::sin(theta);
  if (std::abs(sinTheta)<1.0e-10)
    return Core::Result<SphericalField>::Failure(Core::StatusCode::OutOfDomain,
        "spherical PFSS components are singular at a coordinate pole");
  for (const auto& c:coefficients_) {
    const int l=c.degree;
    const double denominator=l*std::pow(solarRadiusM_,l-1)+
      (l+1.0)*std::pow(sourceSurfaceRadiusM_,2*l+1)*
      std::pow(solarRadiusM_,-l-2);
    const double a=-1.0/denominator;
    const double rbPower=std::pow(sourceSurfaceRadiusM_,2*l+1);
    const double f=a*(std::pow(r,l)-rbPower*std::pow(r,-l-1));
    const double fp=a*(l*std::pow(r,l-1)+(l+1.0)*rbPower*std::pow(r,-l-2));
    const double fpp=a*(l*(l-1.0)*std::pow(r,l-2)-
        (l+1.0)*(l+2.0)*rbPower*std::pow(r,-l-3));
    for (int part=0;part<(c.order==0?1:2);++part) {
      const double amplitude=part==0?c.cosineT:c.sineT;
      const Basis y=RealBasis(l,c.order,part==1,theta,phi);
      field.brT += -amplitude*fp*y.y;
      field.bThetaT += -amplitude*f/r*y.theta;
      field.bPhiT += -amplitude*f/(r*sinTheta)*y.phi;
      field.derivative[0][0] += -amplitude*fpp*y.y;
      field.derivative[0][1] += -amplitude*fp*y.theta;
      field.derivative[0][2] += -amplitude*fp*y.phi;
      const double radialTangential=fp/r-f/(r*r);
      field.derivative[1][0] += -amplitude*radialTangential*y.theta;
      field.derivative[1][1] += -amplitude*f/r*y.theta2;
      field.derivative[1][2] += -amplitude*f/r*y.thetaPhi;
      field.derivative[2][0] += -amplitude*radialTangential/sinTheta*y.phi;
      field.derivative[2][1] += -amplitude*f/r*(y.thetaPhi/sinTheta-
          y.phi*std::cos(theta)/(sinTheta*sinTheta));
      field.derivative[2][2] += -amplitude*f/(r*sinTheta)*y.phi2;
    }
  }
  return Core::Result<SphericalField>::Success(field);
}

Core::Result<Vec3> PfssHarmonics::EvaluateCartesian(Vec3 x) const {
  const double r=Norm(x);
  if (!(r>0.0)) return Core::Result<Vec3>::Failure(Core::StatusCode::OutOfDomain,
      "PFSS Cartesian query cannot be at the origin");
  const double theta=std::acos(std::max(-1.0,std::min(1.0,x.z/r)));
  const double phi=std::atan2(x.y,x.x);
  // Exact coordinate poles are evaluated at a vanishingly displaced azimuthal
  // limit so the returned Cartesian vector remains finite.
  const double safeTheta=std::max(1.0e-9,std::min(Constants::kPi-1.0e-9,theta));
  const auto field=Evaluate(r,safeTheta,phi);
  if (!field.ok()) return Core::Result<Vec3>::Failure(field.status.code,field.status.message);
  return Core::Result<Vec3>::Success(SphericalToCartesian(safeTheta,phi,
      field.value.brT,field.value.bThetaT,field.value.bPhiT));
}

Core::Result<FieldLineTopology> PfssHarmonics::Classify(Vec3 start,double step,int maximum) const {
  if (!(step>0.0&&maximum>0)) return Core::Result<FieldLineTopology>::Failure(
      Core::StatusCode::InvalidConfiguration,"field-line step/count must be positive");
  const double epsilon=2.0*step;
  auto trace=[&](double sign)->Core::Result<FieldLineTopology>{
    Vec3 x=start;
    auto direction=[&](Vec3 p)->Core::Result<Vec3>{
      const auto b=EvaluateCartesian(p);
      if (!b.ok()||Norm(b.value)<1.0e-30) return Core::Result<Vec3>::Failure(
          Core::StatusCode::NumericalFailure,"field-line trace encountered a null");
      return Core::Result<Vec3>::Success(sign*Unit(b.value));
    };
    for(int n=0;n<maximum;++n) {
      const double radius=Norm(x);
      if (radius>=sourceSurfaceRadiusM_-epsilon)
        return Core::Result<FieldLineTopology>::Success(FieldLineTopology::OpenToOuterBoundary);
      if (n>0&&radius<=solarRadiusM_+epsilon)
        return Core::Result<FieldLineTopology>::Success(FieldLineTopology::ClosedBelowOuterBoundary);
      const auto k1=direction(x); if(!k1.ok()) return Core::Result<FieldLineTopology>::Success(FieldLineTopology::Separatrix);
      const auto k2=direction(x+0.5*step*k1.value); if(!k2.ok()) return Core::Result<FieldLineTopology>::Success(FieldLineTopology::Separatrix);
      const auto k3=direction(x+0.5*step*k2.value); if(!k3.ok()) return Core::Result<FieldLineTopology>::Success(FieldLineTopology::Separatrix);
      const auto k4=direction(x+step*k3.value); if(!k4.ok()) return Core::Result<FieldLineTopology>::Success(FieldLineTopology::Separatrix);
      x=x+(step/6.0)*(k1.value+2.0*k2.value+2.0*k3.value+k4.value);
    }
    return Core::Result<FieldLineTopology>::Success(FieldLineTopology::Separatrix);
  };
  const auto forward=trace(1.0), backward=trace(-1.0);
  if (!forward.ok()) return forward;
  if (!backward.ok()) return backward;
  if (forward.value==FieldLineTopology::OpenToOuterBoundary||
      backward.value==FieldLineTopology::OpenToOuterBoundary)
    return Core::Result<FieldLineTopology>::Success(FieldLineTopology::OpenToOuterBoundary);
  if (forward.value==FieldLineTopology::ClosedBelowOuterBoundary&&
      backward.value==FieldLineTopology::ClosedBelowOuterBoundary)
    return Core::Result<FieldLineTopology>::Success(FieldLineTopology::ClosedBelowOuterBoundary);
  return Core::Result<FieldLineTopology>::Success(FieldLineTopology::Separatrix);
}

PfssHarmonics PfssHarmonics::Scaled(double factor) const {
  PfssHarmonics copy=*this;
  for(auto& c:copy.coefficients_) { c.cosineT*=factor; c.sineT*=factor; }
  return copy;
}

PfssHarmonics PfssHarmonics::HeatKernelFiltered(int degree) const {
  PfssHarmonics copy=*this;
  if(degree<=0) return copy;
  for(auto& c:copy.coefficients_) {
    const double transfer=std::exp(-c.degree*(c.degree+1.0)/(degree*(degree+1.0)));
    c.cosineT*=transfer; c.sineT*=transfer;
  }
  return copy;
}

Core::Result<std::vector<HarmonicCoefficient>> PfssHarmonics::ProjectMap(
    const std::vector<MapSample>& samples,int maximum,MonopolePolicy policy,
    double tolerance,ReconstructionMetrics* metrics) {
  if(samples.empty()||maximum<1) return Core::Result<std::vector<HarmonicCoefficient>>::Failure(
      Core::StatusCode::InvalidConfiguration,"map projection requires samples and Lmax>=1");
  double area=0.0,monopoleIntegral=0.0;
  for(const auto& s:samples) {
    if(!(s.solidAngleWeightSr>0.0&&Finite(s.radialFieldT)))
      return Core::Result<std::vector<HarmonicCoefficient>>::Failure(
          Core::StatusCode::InvalidConfiguration,"map weights and field must be finite/positive");
    area+=s.solidAngleWeightSr; monopoleIntegral+=s.radialFieldT*s.solidAngleWeightSr;
  }
  const double mean=monopoleIntegral/area;
  if(std::abs(mean)>tolerance&&policy==MonopolePolicy::Reject)
    return Core::Result<std::vector<HarmonicCoefficient>>::Failure(
        Core::StatusCode::InvalidConfiguration,"map monopole exceeds tolerance");
  std::vector<HarmonicCoefficient> result;
  for(int l=1;l<=maximum;++l) for(int m=0;m<=l;++m) {
    HarmonicCoefficient c{l,m,0.0,0.0};
    for(const auto& s:samples) {
      const double value=s.radialFieldT-(policy==MonopolePolicy::Remove?mean:0.0);
      c.cosineT+=value*RealBasis(l,m,false,s.thetaRad,s.phiRad).y*s.solidAngleWeightSr;
      if(m>0)c.sineT+=value*RealBasis(l,m,true,s.thetaRad,s.phiRad).y*s.solidAngleWeightSr;
    }
    result.push_back(c);
  }
  if(metrics) { metrics->removedMonopoleT=policy==MonopolePolicy::Remove?mean:0.0; }
  return Core::Result<std::vector<HarmonicCoefficient>>::Success(std::move(result));
}

Core::Result<ReconstructionMetrics> PfssHarmonics::CompareMap(
    const std::vector<MapSample>& samples,const std::vector<HarmonicCoefficient>& coefficients) {
  if(samples.empty()) return Core::Result<ReconstructionMetrics>::Failure(
      Core::StatusCode::InvalidConfiguration,"map comparison requires samples");
  double weight=0.0,error=0.0,rawUnsigned=0.0,reconstructedUnsigned=0.0;
  for(const auto& s:samples) {
    double value=0.0;
    for(const auto& c:coefficients) {
      value+=c.cosineT*RealBasis(c.degree,c.order,false,s.thetaRad,s.phiRad).y;
      if(c.order>0)value+=c.sineT*RealBasis(c.degree,c.order,true,s.thetaRad,s.phiRad).y;
    }
    weight+=s.solidAngleWeightSr; error+=(value-s.radialFieldT)*(value-s.radialFieldT)*s.solidAngleWeightSr;
    rawUnsigned+=std::abs(s.radialFieldT)*s.solidAngleWeightSr;
    reconstructedUnsigned+=std::abs(value)*s.solidAngleWeightSr;
  }
  ReconstructionMetrics result; result.weightedRmsT=std::sqrt(error/weight);
  result.unsignedFluxChangeFraction=rawUnsigned>0.0?
      (reconstructedUnsigned-rawUnsigned)/rawUnsigned:0.0;
  return Core::Result<ReconstructionMetrics>::Success(result);
}

Core::Result<TopologyPair> RouteTopology(FieldLineTopology pfss,
    FieldLineTopology composite,std::uint64_t generation,const std::string& authority) {
  if(pfss==FieldLineTopology::Invalid||composite==FieldLineTopology::Invalid)
    return Core::Result<TopologyPair>::Failure(Core::StatusCode::InvalidState,
        "both PFSS and composite topology classifications are required");
  if(authority!="pfss-open-to-Rb"&&authority!="composite-open-to-Ri")
    return Core::Result<TopologyPair>::Failure(Core::StatusCode::InvalidConfiguration,
        "unknown target-speed topology authority");
  TopologyPair pair{pfss,composite,generation,
      composite==FieldLineTopology::OpenToOuterBoundary?"open-wind":"closed-hydrostatic",
      authority};
  return Core::Result<TopologyPair>::Success(std::move(pair));
}

} }  // namespace SEP::CoronalCME
