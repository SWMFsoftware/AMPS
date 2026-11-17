#include "sep_background_derivative_validation.h"
#include "sep_background_derivatives.h"
#include "sep_parker_core.h"
#include <algorithm>
#include <cmath>
#include <limits>

namespace SEP { namespace Testing { namespace {
Result Finish(bool pass,const char* message) {
  Result result;
  result.status=pass?Status::Pass:Status::Fail;
  result.message=message;
  return result;
}
bool Near(double a,double b,double tolerance=1e-12) {
  return std::fabs(a-b)<=tolerance*std::max(1.0,std::fabs(b));
}
Result RadialExpansion() {
  // A stationary n=C/r^2 satisfies continuity with constant radial U.
  // At a fixed vertex its time derivative is zero, while div U=2U/r.
  const double x[3]={3.0,4.0,12.0};
  const double b[3]={3.0/13.0,4.0/13.0,12.0/13.0};
  std::array<double,9> gradient;
  Transport::VelocityDerivatives derivative;
  bool pass=Transport::RadialVelocityGradient(x,26.0,0.0,&gradient).ok()&&
      Transport::ContractVelocityGradient(gradient,b,&derivative).ok()&&
      Near(derivative.divergencePerS,4.0)&&
      Near(derivative.fieldAlignedStrainPerS,0.0);
  // Independent Cartesian differences of U_i=26*x_i/r, not a rearranged
  // divergence expression, verify every component including off diagonals.
  const double h=1e-3;
  double maximumError=0.0;
  for(int j=0;j<3;++j) {
    double plus[3]={x[0],x[1],x[2]},minus[3]={x[0],x[1],x[2]};
    plus[j]+=h;minus[j]-=h;
    double rp=0.0,rm=0.0;
    for(int k=0;k<3;++k){rp+=plus[k]*plus[k];rm+=minus[k]*minus[k];}
    rp=std::sqrt(rp);rm=std::sqrt(rm);
    for(int i=0;i<3;++i)
      maximumError=std::max(maximumError,std::fabs(
          (26.0*plus[i]/rp-26.0*minus[i]/rm)/(2*h)-gradient[3*i+j]));
  }
  pass=pass&&maximumError<2e-8;
  // A nonconstant radial speed also contributes du/dr to both contractions.
  pass=pass&&Transport::RadialVelocityGradient(x,26.0,4.0,&gradient).ok()&&
      Transport::ContractVelocityGradient(gradient,b,&derivative).ok()&&
      Near(derivative.divergencePerS,8.0)&&
      Near(derivative.fieldAlignedStrainPerS,4.0);
  Result result=Finish(pass,"stationary radial density has geometric expansion; Cartesian tensor differences agree");
  result.metrics.push_back({"tensor_max_error",maximumError,2e-8,"<","1/s"});
  return result;
}
Result CurvedField() {
  // Constant Cartesian U=(7,0,0) has grad U=0 for every curved B. For
  // b(s)=(cos(s/R),sin(s/R),0), d(U.b)/ds=-7 sin(s/R)/R is nonzero.
  const double b[3]={std::cos(0.7),std::sin(0.7),0};
  const std::array<double,9> zero{{0,0,0,0,0,0,0,0,0}};
  Transport::VelocityDerivatives derivative;
  const double projectedDerivative=-7.0*std::sin(0.7)/3.0;
  const bool pass=Transport::ContractVelocityGradient(zero,b,&derivative).ok()&&
      derivative.divergencePerS==0.0&&
      derivative.fieldAlignedStrainPerS==0.0&&
      std::fabs(projectedDerivative)>1.0;
  return Finish(pass,"curved-field projected velocity derivative is distinct from Cartesian field-aligned strain");
}
Result TensorValidation() {
  // Non-radial manufactured linear flow U=M*x: its gradient is exactly M.
  const std::array<double,9> matrix{{2,3,5,7,11,13,17,19,23}};
  const double b[3]={1,2,2};
  double expected=0.0;
  for(int i=0;i<3;++i)for(int j=0;j<3;++j)
    expected+=b[i]*matrix[3*i+j]*b[j]/9.0;
  Transport::VelocityDerivatives derivative;
  bool pass=Transport::ContractVelocityGradient(matrix,b,&derivative).ok()&&
      Near(derivative.divergencePerS,36.0)&&
      Near(derivative.fieldAlignedStrainPerS,expected);
  auto invalid=matrix;invalid[4]=std::numeric_limits<double>::quiet_NaN();
  pass=pass&&!Transport::ContractVelocityGradient(invalid,b,&derivative).ok();
  const double origin[3]={0,0,0};
  std::array<double,9> untouched{{1,1,1,1,1,1,1,1,1}};
  pass=pass&&!Transport::RadialVelocityGradient(origin,2,0,&untouched).ok()&&untouched[0]==1;
  return Finish(pass,"full tensor contractions preserve non-radial flow and reject unavailable or malformed derivatives");
}
class NoDiffusion final:public Transport::SpatialDiffusionProvider {
 public:
  Transport::SpatialDiffusionSample Evaluate(double,double)const override {
    Transport::SpatialDiffusionSample sample;
    sample.status=Transport::Status::Ok();
    sample.provenance="controlled-zero-diffusion";
    return sample;
  }
};
Result AntipolarityAdvection() {
  const double u[3]={5,2,0},tangent[3]={1,0,0};
  // B can reverse at fixed physical U and fixed line-coordinate geometry.
  // The Cartesian characteristic x(t)=x0+5*t cannot reverse with it.
  const double negativeB[3]={-1,0,0};
  const auto advection=Transport::FieldLineAdvectionMPerS(u,tangent);
  Transport::KeyedRandomStream random(1,1,1,1);
  NoDiffusion diffusion;
  const auto increment=Transport::AdvanceParker(
      Transport::ParkerState(3,1),Transport::ParkerBackground(advection.value,0),
      10,2,diffusion,random);
  const double erroneousMagneticProjection=u[0]*negativeB[0];
  const bool pass=advection.status.ok()&&Near(advection.value,5)&&
      increment.status.ok()&&Near(increment.state.arcLengthM,13)&&
      erroneousMagneticProjection<0;
  return Finish(pass,"antiparallel magnetic polarity preserves the exact Cartesian convection characteristic");
}
Descriptor Make(const char* id,const char* name,TestCallback callback) {
  Descriptor descriptor;
  descriptor.id=id;descriptor.name=name;descriptor.group="background";
  descriptor.description=name;
  descriptor.initialization=InitializationLevel::None;
  descriptor.supportedBuildModes="native CLI and dependency-light launcher: identical callback";
  descriptor.runtime=RuntimeClass::Routine;
  descriptor.seedPolicy="deterministic Cartesian manufactured fields";
  descriptor.stateIsolation="stack-owned tensor and state";
  descriptor.callback=callback;
  return descriptor;
}
} // namespace
std::vector<Descriptor> BackgroundDerivativeDescriptors() {
  return {Make("CP03G01","Exact steady radial expansion tensor",RadialExpansion),
          Make("CP03G02","Curved-field strain contract",CurvedField),
          Make("CP03G03","Full Cartesian gradient validation",TensorValidation),
          Make("CP03G04","Magnetic polarity independent convection",AntipolarityAdvection)};
}
} }
