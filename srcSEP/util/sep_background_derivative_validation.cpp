#include "sep_background_derivative_validation.h"
#include "sep_background_derivatives.h"
#include "sep_parker_core.h"
#include "amps_velocity_gradient.h"
#include "amps_fieldline_frame.h"
#include "sep_focused_transport_core.h"
#include "sep_focused_transport_mfp_core.h"
#include "sep_background_runtime.h"
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
Result ImportedAffineGradient() {
  // Independent manufactured Cartesian flow U=U0+M*(x-x0). Dyadic M and
  // metre offsets make the supplied samples exact, even at a solar-scale
  // translated origin. Off-diagonal/asymmetric entries detect tensor swaps.
  using Kernel=AMPS::Coupling::VelocityGradient;
  const double matrix[9]={2,3,5,7,11,13,17,19,23};
  std::vector<Kernel::Sample> samples;
  for(int x=-1;x<=1;x+=2)for(int y=-1;y<=1;y+=2)for(int z=-1;z<=1;z+=2) {
    Kernel::Sample s;s.line=static_cast<std::uint64_t>((x+1)+(y+1)*2+1);
    const double d[3]={x*1048576.0,y*1048576.0,z*1048576.0};
    for(int i=0;i<3;++i){s.x[i]=1073741824.0+d[i];s.u[i]=131072.0;
      for(int j=0;j<3;++j)s.u[i]+=std::ldexp(matrix[3*i+j],-20)*d[j];}
    samples.push_back(s);
  }
  Kernel kernel(samples);bool pass=true;double error=0;
  for(std::size_t i=0;i<samples.size();++i) {
    const auto result=kernel.Evaluate(i);pass=pass&&result.valid;
    for(int j=0;j<9;++j)error=std::max(error,std::abs(result.gradient[j]-std::ldexp(matrix[j],-20)));
  }
  pass=pass&&error<2e-16;
  auto result=Finish(pass,"received three-dimensional affine velocities recover the full Cartesian tensor at translated solar-scale positions");
  result.metrics.push_back({"maximum_affine_gradient_error",error,2e-16,"<","1/s"});return result;
}
Result ImportedRankRejection() {
  // Four labelled lines on z=0 still provide no derivative normal to that
  // plane. A minimum-norm/radial/zero replacement must never claim success.
  using Kernel=AMPS::Coupling::VelocityGradient;std::vector<Kernel::Sample> samples;
  for(int i=0;i<4;++i)for(int j=0;j<4;++j) {
    Kernel::Sample s;s.line=i+1;s.x={{double(i),double(j),0.0}};s.u={{double(i),double(j),1.0}};samples.push_back(s);
  }
  Kernel kernel(samples);bool pass=true;
  for(std::size_t i=0;i<samples.size();++i)pass=pass&&!kernel.Evaluate(i).valid;
  return Finish(pass,"coplanar received velocities reject without invented transverse derivatives");
}
Result ImportedGradientConvergence() {
  // At the origin grad(x^3,y^3,z^3)=0. On geometrically similar 3-D stencils,
  // the leading affine-reconstruction error is O(h^2); halving h must reduce
  // it by four. This independent polynomial oracle checks spatial error,
  // separately from affine exactness and QR roundoff/conditioning controls.
  using Kernel=AMPS::Coupling::VelocityGradient;double errors[3]={};
  bool pass=true;
  for(int level=0;level<3;++level) {
    const double h=std::ldexp(1.0,-level);std::vector<Kernel::Sample> samples;
    for(int x=-1;x<=1;++x)for(int y=-1;y<=1;++y)for(int z=-1;z<=1;++z) {
      Kernel::Sample s;s.line=(x+1)*3+(y+1)+1;s.x={{x*h,y*h,z*h}};
      for(int d=0;d<3;++d)s.u[d]=s.x[d]*s.x[d]*s.x[d];
      samples.push_back(s);
    }
    Kernel kernel(samples);const auto result=kernel.Evaluate(13);pass=pass&&result.valid;
    for(double value:result.gradient)errors[level]=std::max(errors[level],std::abs(value));
  }
  pass=pass&&errors[2]>0&&Near(errors[0]/errors[1],4.0,1e-12)&&Near(errors[1]/errors[2],4.0,1e-12);
  auto result=Finish(pass,"internal velocity-gradient reconstruction converges quadratically for a symmetric cubic manufactured field");
  result.metrics.push_back({"coarse_to_medium_error_ratio",errors[0]/errors[1],4.0,"observed","dimensionless"});return result;
}
// These coefficients explicitly disable scattering. Analytic characteristics,
// not a second implementation of the mover, supply the expected positions.
struct NoPitch : Transport::PitchAngleDiffusionProvider {
  Transport::PitchAngleDiffusionSample Evaluate(double,double,double) const override {
    Transport::PitchAngleDiffusionSample s;s.status=Transport::Status::Ok();
    s.provenance="frame-test:no-scattering";s.turbulenceStateIdentity="frame-test";return s;
  }
};
struct NoEvents : Transport::MeanFreePathProvider {
  Transport::MeanFreePathSample Evaluate(double,double,double) const override {
    Transport::MeanFreePathSample s;s.status=Transport::Status::Ok();
    s.lambdaParallelM=std::numeric_limits<double>::infinity();
    s.provenance="frame-test:ballistic";s.turbulenceStateIdentity="frame-test";return s;
  }
};
Result CoordinateFrames() {
  using namespace Transport;using namespace AMPS::FieldLine;
  const double u[3]={400000,20000,0},t[3]={1,0,0};
  const double mass=1.67262192369e-27,c=299792458.0,v=2e6,mu=0.25,dt=1;
  const auto momentum=MomentumFromSpeed(v,mass,c);bool pass=momentum.status.ok();
  NoPitch pitch;NoEvents events;NoDiffusion diffusion;
  // Exercise the actual provider-to-coordinate contract, rather than assume
  // that model-generated geometry is Eulerian. Each data authority must give
  // the identical plasma-frame characteristic. The generic fixed-grid helper
  // is checked separately below; it is not a srcSEP provider selection.
  for(auto provider:{Background::Provider::Analytic,Background::Provider::Swcme,
      Background::Provider::ReducedShock,Background::Provider::Swmf,
      Background::Provider::LocalEvolution}) {
    const auto frame=Background::DeclaredParticleCoordinateFrame(provider);
    pass=pass&&frame==CoordinateFrame::PlasmaComoving;
    const double adv=CoordinateAdvection(frame,FieldLineAdvectionMPerS(u,t).value);
    const double expected=v*mu;
    KeyedRandomStream random(1,2,3,4);
    auto a=AdvanceFocusedTransportDmumu(FocusedTransportState(0,momentum.value,mu),
        FocusedTransportBackground(0,adv,0,0),mass,c,dt,pitch,random,nullptr);
    auto b=AdvanceFocusedTransportMfp(FocusedTransportMfpState(0,momentum.value,mu),
        FocusedTransportMfpBackground(0,adv,0,0,0),mass,c,dt,dt,events,random,nullptr);
    auto p=AdvanceParker(ParkerState(0,momentum.value),ParkerBackground(adv,0),v,dt,diffusion,random);
    pass=pass&&a.status.ok()&&b.status.ok()&&p.status.ok()&&Near(a.displacementM,expected,2e-14)&&
        Near(b.displacementM,expected,2e-14)&&Near(p.displacementM,0);
  }
  pass=pass&&CoordinateAdvection(CoordinateFrame::Fixed,400000)==400000;
  try{Background::DeclaredParticleCoordinateFrame(static_cast<Background::Provider>(-1));pass=false;}
  catch(const std::exception&){}
  // A plasma-rest particle retains its material coordinate. Translating the
  // endpoints by U*dt then moves Cartesian X by exactly U*dt, never 2*U*dt.
  const double adv=CoordinateAdvection(CoordinateFrame::PlasmaComoving,400000);
  KeyedRandomStream random(4,3,2,1);
  const auto stationary=AdvanceFocusedTransportMfp(FocusedTransportMfpState(0,0,0),
      FocusedTransportMfpBackground(0,adv,0,0,0),mass,c,dt,dt,events,random,nullptr);
  pass=pass&&stationary.status.ok()&&stationary.displacementM==0;
  return Finish(pass,"all srcSEP data providers use plasma-comoving characteristics in both focused movers and Parker; generic Eulerian control retains U");
}
Result MaterialLabels() {
  using namespace AMPS::FieldLine;
  // Integer node insertion changes the index, not the material fraction.
  // Stretching/translating endpoints changes X without changing that fraction.
  const double s=RemapCoordinate(3.25,10,8),fraction=s-std::floor(s);
  const double x=100+fraction*8;
  bool pass=s==5.25&&RemapCoordinate(s,8,10)==3.25&&x==102;
  try{RemapCoordinate(std::numeric_limits<double>::quiet_NaN(),10,8);pass=false;}catch(const std::exception&){}
  try{CoordinateAdvection(CoordinateFrame::Unknown,3);pass=false;}catch(const std::exception&){}
  return Finish(pass,"material-label remap preserves fractions under prefix changes and rejects unknown frames/nonfinite coordinates");
}
Result ComovingCompression() {
  using namespace Transport;using namespace AMPS::FieldLine;
  // Homogeneous isotropic expansion has div U=3a and strain=a. Removing
  // coordinate advection must retain p(t)=p0*exp(-a*t) in all three movers.
  const double a=0.01,dt=1,mass=1.67262192369e-27,c=299792458.0,p0=1e-20;
  const std::array<double,9> g{{a,0,0,0,a,0,0,0,a}};const double bdir[3]={1,0,0};
  VelocityDerivatives d;bool pass=ContractVelocityGradient(g,bdir,&d).ok();
  const double adv=CoordinateAdvection(CoordinateFrame::PlasmaComoving,400000);
  FocusedTransportBackground bg(0,adv,0,d.divergencePerS);bg.fieldAlignedStrainPerS=d.fieldAlignedStrainPerS;bg.equationMode=FocusedEquationMode::FullGyrotropic;
  FocusedTransportMfpBackground mb(0,adv,0,d.divergencePerS,0);mb.fieldAlignedStrainPerS=d.fieldAlignedStrainPerS;mb.equationMode=FocusedEquationMode::FullGyrotropic;
  KeyedRandomStream random(1,2,3,4);NoPitch pitch;NoEvents events;NoDiffusion diffusion;
  auto f=AdvanceFocusedTransportDmumu(FocusedTransportState(0,p0,0.25),bg,mass,c,dt,pitch,random,nullptr);
  auto m=AdvanceFocusedTransportMfp(FocusedTransportMfpState(0,p0,0.25),mb,mass,c,dt,dt,events,random,nullptr);
  auto p=AdvanceParker(ParkerState(0,p0),ParkerBackground(adv,d.divergencePerS),1,dt,diffusion,random);
  const double expected=p0*std::exp(-a*dt);
  pass=pass&&f.status.ok()&&m.status.ok()&&p.status.ok()&&
      std::abs(f.state.momentumKgMPerS/expected-1)<1e-13&&
      std::abs(m.state.momentumKgMPerS/expected-1)<1e-13&&std::abs(p.state.momentumKgMPerS/expected-1)<1e-13;
  return Finish(pass,"comoving coordinate correction retains independent physical compression and strain in all movers");
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
          Make("CP03G04","Magnetic polarity independent convection",AntipolarityAdvection),
          Make("CP03G05","Received Cartesian affine velocity reconstruction",ImportedAffineGradient),
          Make("CP03G06","Underdetermined imported derivative rejection",ImportedRankRejection),
          Make("CP03G07","Received velocity gradient spatial convergence",ImportedGradientConvergence),
          Make("CP03G08","Fixed and comoving transport characteristics",CoordinateFrames),
          Make("CP03G09","Material label and geometry invariants",MaterialLabels),
          Make("CP03G10","Comoving transport retains compression",ComovingCompression)};
}
} }
