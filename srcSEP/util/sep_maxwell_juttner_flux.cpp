#include "sep_maxwell_juttner_flux.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace SEP { namespace Sources { namespace MaxwellJuttner { namespace {
using Transport::Status;
using Transport::StatusCode;
using Wide=long double;
using Vector=std::array<Wide,3>;
constexpr Wide Pi=3.141592653589793238462643383279502884L;
constexpr Wide C=static_cast<Wide>(SpeedOfLightMPerS);

Status Invalid(const char* message){return Status::Error(StatusCode::InvalidArgument,message);}
Status Unresolved(const char* message){return Status::Error(StatusCode::UnresolvedCoefficient,message);}
Wide Dot(const Vector& a,const Vector& b){return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];}
Wide Norm(const Vector& a){return std::hypot(a[0],std::hypot(a[1],a[2]));}
Vector Cross(const Vector& a,const Vector& b){return {{a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]}};}
bool Representable(Wide value,double* output,bool positive=false) {
  if(!std::isfinite(value)||std::abs(value)>std::numeric_limits<double>::max())return false;
  const double narrowed=static_cast<double>(value);
  if((positive&&!(narrowed>0.0))||(value!=0.0&&narrowed==0.0))return false;
  *output=narrowed;return true;
}
struct Prepared {
  Vector beta,normal,a,axis;
  Wide bulkGamma=1.0L,theta=0.0L,aMagnitude=0.0L,massTimesC=0.0L;
  Wide arithmeticRelativeBudget=0.0L;
};
Status Prepare(const MaterialBoundaryInput& input,Prepared* output) {
  if(!output)return Invalid("material Juttner preparation requires an output");
  if(!std::isfinite(input.properNumberDensityPerM3)||input.properNumberDensityPerM3<0.0||
      !std::isfinite(input.properTemperatureK)||input.properTemperatureK<0.0||
      !std::isfinite(input.restMassKg)||!(input.restMassKg>0.0))
    return Invalid("material Juttner requires finite nonnegative proper density/temperature and positive rest mass in SI");
  Prepared prepared;Vector bulk;
  for(unsigned i=0;i<3;++i) {
    if(!std::isfinite(input.bulkVelocityMPerS[i])||!std::isfinite(input.boundaryVelocityMPerS[i])||
        !std::isfinite(input.inwardNormal[i]))return Invalid("material Juttner vectors must be finite");
    // Exact equality declares the same material velocity, rather than a
    // tolerance that quietly turns a nonmaterial plane into a material one.
    if(input.bulkVelocityMPerS[i]!=input.boundaryVelocityMPerS[i])
      return Status::Error(StatusCode::UnsupportedConfiguration,"Juttner material-plane kernel requires boundary velocity exactly equal to bulk U; nonmaterial flux law is separate");
    bulk[i]=input.bulkVelocityMPerS[i];prepared.normal[i]=input.inwardNormal[i];prepared.beta[i]=bulk[i]/C;
  }
  const Wide speed=Norm(bulk),normalNorm=Norm(prepared.normal);
  if(!(speed<C))return Invalid("Juttner bulk velocity must be strictly subluminal");
  if(!(normalNorm>0.0L)||!std::isfinite(normalNorm))return Invalid("Juttner inward normal direction must be finite and nonzero");
  for(auto& component:prepared.normal)component/=normalNorm;
  const Wide beta=speed/C;
  prepared.bulkGamma=1.0L/std::sqrt((1.0L-beta)*(1.0L+beta));
  // Forming Gamma from a finite-precision norm is conditioned by Gamma^2:
  // d(log Gamma)/d beta=beta*Gamma^2. A transverse cap does not expose that
  // error in |a|, so checking the normal invariant alone is insufficient.
  // Reserve a conservative working-precision envelope for the norm/dot,
  // square root/division, normal rotation, radial construction and boost.
  // This is a numerical capability bound, not a physical bulk-velocity cap.
  prepared.arithmeticRelativeBudget=64.0L*std::numeric_limits<Wide>::epsilon()*
      prepared.bulkGamma*prepared.bulkGamma;
  if(!std::isfinite(prepared.arithmeticRelativeBudget)||
      prepared.arithmeticRelativeBudget>CrossingRelativeResolution)
    return Unresolved("bulk Gamma formation is too ill-conditioned for the declared material-crossing relative precision");
  const Wide normalBeta=Dot(prepared.normal,prepared.beta);
  // Lorentz-transformed material-plane normal: a=n-[Gamma/(Gamma+1)]
  // (beta dot n) beta. Keeping the cancellation in wider arithmetic protects
  // a nearly parallel relativistic cap. Its invariant length is independently
  // checked; an unresolvable boost is a numerical error, not an axis fallback.
  const Wide coefficient=prepared.bulkGamma/(prepared.bulkGamma+1.0L);
  for(unsigned i=0;i<3;++i)prepared.a[i]=prepared.normal[i]-coefficient*normalBeta*prepared.beta[i];
  prepared.aMagnitude=Norm(prepared.a);
  const Wide invariantMagnitude=std::sqrt((1.0L-normalBeta)*(1.0L+normalBeta));
  if(!(prepared.aMagnitude>0.0L)||!std::isfinite(prepared.aMagnitude)||
      std::abs(prepared.aMagnitude-invariantMagnitude)>2.0e-12L*invariantMagnitude)
    return Unresolved("material normal Lorentz transformation is not resolved to relative tolerance 2e-12");
  for(unsigned i=0;i<3;++i)prepared.axis[i]=prepared.a[i]/prepared.aMagnitude;
  prepared.massTimesC=static_cast<Wide>(input.restMassKg)*C;
  prepared.theta=static_cast<Wide>(BoltzmannJPerK)*input.properTemperatureK/(prepared.massTimesC*C);
  double theta=0.0;
  if(input.properTemperatureK>0.0&&!Representable(prepared.theta,&theta,true))
    return Unresolved("proper dimensionless Juttner temperature is outside the binary64 numerical domain");
  *output=prepared;return Status::Ok();
}

Status MeanSpeed(Wide theta,Wide* fraction,NormalizationMethod* method,Wide* remainderBound) {
  if(theta==0.0L){*fraction=0.0L;*method=NormalizationMethod::ColdZero;*remainderBound=0.0L;return Status::Ok();}
  Wide beta=0.0L,bound=0.0L;
  if(theta<=1.0L/32.0L) {
    // NIST DLMF 10.40.2, 10.40(ii), nu=2: after at least two
    // terms the remainder has the first omitted term's sign and magnitude
    // no larger than that term. Terms are formed in theta, so 1/theta and
    // exp(1/theta) never overflow for a physically cold positive temperature.
    Wide sum=1.0L,term=1.0L;bool converged=false;
    for(unsigned k=1;k<=128;++k) {
      const Wide odd=2.0L*k-1.0L;
      term*=((16.0L-odd*odd)/(8.0L*k))*theta;
      if(k>=2&&std::abs(term)<=1.0e-17L*std::abs(sum)) {
        bound=std::abs(term)/(sum-std::abs(term));converged=true;break;
      }
      if(!std::isfinite(term)||!std::isfinite(sum+term))break;
      sum+=term;
    }
    if(!converged||!(sum>0.0L))return Unresolved("bounded scaled-K2 asymptotic did not reach its declared remainder tolerance");
    beta=std::sqrt(8.0L*theta/Pi)*(1.0L+theta)/sum;
    *method=NormalizationMethod::BoundedAsymptotic;
  } else {
    const double inverse=static_cast<double>(1.0L/theta);
    try {
      if(inverse>=1.0) {
        const Wide scaled=std::exp(static_cast<Wide>(inverse))*std::cyl_bessel_k(2.0,inverse);
        beta=2.0L*theta*(1.0L+theta)/scaled;
      } else {
        // K2(z)=K0(z)+2*K1(z)/z. Form z^2 exp(z)K2(z), never
        // K2 itself, which overflows for a perfectly finite hot flux. Wider
        // products keep the small z^2 K0 term without premature underflow.
        const Wide z=inverse;
        const Wide scaled=std::exp(z)*(2.0L*z*std::cyl_bessel_k(1.0,inverse)+z*z*std::cyl_bessel_k(0.0,inverse));
        beta=2.0L*(1.0L+z)/scaled;
      }
    } catch(const std::exception&) {
      return Unresolved("C++ modified-Bessel evaluation failed for the declared proper temperature");
    }
    *method=NormalizationMethod::LibraryBessel;
  }
  // Binary64 library values carry roundoff. No clipping repairs a value above
  // the physical light-speed bound: reject it unless a wider representation
  // already lies in [0,1]. A rounded value equal to one is valid for a rate;
  // sampled four-momentum, rather than a clipped velocity, remains authoritative.
  if(!std::isfinite(beta)||!(beta>0.0L)||beta>1.0L)
    return Unresolved("Juttner mean speed is nonfinite or violates its physical light-speed bound");
  *fraction=beta;*remainderBound=bound;return Status::Ok();
}

Status DrawPrepared(const MaterialBoundaryInput& input,const Prepared& prepared,
    const std::array<double,6>& uniforms,FluxSample* output) {
  if(!output)return Invalid("material Juttner sample requires an output");
  if(!(input.properNumberDensityPerM3>0.0)||!(input.properTemperatureK>0.0))
    return Status::Error(StatusCode::OutOfDomain,"zero-density or cold material plasma has zero crossing flux and no conditional birth distribution");
  for(double value:uniforms)if(!std::isfinite(value)||!(value>0.0&&value<1.0))
    return Invalid("material Juttner sampler requires six strictly open unit variates");
  const Wide theta=prepared.theta;
  // Stable exact gamma mixture for x=gamma'-1: normalized x(x+2)e^-x/theta.
  // Compare the smaller probability directly, avoiding cancellation when the
  // minority branch is near the finite random grid's resolution.
  const bool shapeThree=theta<=1.0L?uniforms[0]<theta/(1.0L+theta):
      !(uniforms[0]<(1.0L/theta)/(1.0L+1.0L/theta));
  Wide unitGamma=-std::log(static_cast<Wide>(uniforms[1]))-std::log(static_cast<Wide>(uniforms[2]));
  if(shapeThree)unitGamma-=std::log(static_cast<Wide>(uniforms[3]));
  const Wide x=theta*unitGamma;
  const Wide magnitude=std::sqrt(x)*std::sqrt(x+2.0L);
  // Cosine-weighted incoming hemisphere about aHat. Its CDF is mu^2, not
  // uniform mu. The third exponential slot is intentionally unused by shape
  // two; the keyed wrapper still consumes exactly six uniforms per birth.
  const Wide mu=std::sqrt(static_cast<Wide>(uniforms[4]));
  const Wide transverse=std::sqrt(1.0L-static_cast<Wide>(uniforms[4]));
  const Wide phi=2.0L*Pi*uniforms[5];
  unsigned least=0;for(unsigned i=1;i<3;++i)if(std::abs(prepared.axis[i])<std::abs(prepared.axis[least]))least=i;
  Vector seed{{0.0L,0.0L,0.0L}};seed[least]=1.0L;
  Vector tangent=Cross(prepared.axis,seed);const Wide tangentNorm=Norm(tangent);
  for(auto& component:tangent)component/=tangentNorm;
  const Vector other=Cross(prepared.axis,tangent);
  Vector rest,lab;
  for(unsigned i=0;i<3;++i)rest[i]=magnitude*(mu*prepared.axis[i]+transverse*(std::cos(phi)*tangent[i]+std::sin(phi)*other[i]));
  const Wide restGamma=1.0L+x,betaDotMomentum=Dot(prepared.beta,rest);
  const Wide labGamma=prepared.bulkGamma*(restGamma+betaDotMomentum);
  const Wide boostCoefficient=prepared.bulkGamma*(prepared.bulkGamma/(prepared.bulkGamma+1.0L))*betaDotMomentum+prepared.bulkGamma*restGamma;
  for(unsigned i=0;i<3;++i)lab[i]=rest[i]+boostCoefficient*prepared.beta[i];
  // This positive invariant expression avoids subtracting nearly equal lab
  // normal and material speeds. It equals (v_lab-U) dot n by the full Lorentz
  // momentum transformation. No acceptance/rejection or velocity clamp occurs.
  const Wide inwardSpeed=C*Dot(prepared.a,rest)/labGamma;
  FluxSample sample;
  if(!Representable(x,&sample.restKineticGamma,true)||!Representable(restGamma,&sample.restGamma,true)||
      !Representable(labGamma,&sample.simulationGamma,true)||!Representable(inwardSpeed,&sample.relativeInwardSpeedMPerS,true))
    return Unresolved("sampled Juttner energy or inward speed is not numerically representable; RNG/event remain uncommitted");
  for(unsigned i=0;i<3;++i)if(!Representable(prepared.massTimesC*rest[i],&sample.restMomentumKgMPerS[i])||
      !Representable(prepared.massTimesC*lab[i],&sample.simulationMomentumKgMPerS[i]))
    return Unresolved("sampled Juttner momentum is not numerically representable; tail redraw is prohibited");
  // Component finiteness is not enough: a cold thermal increment can disappear
  // beneath a boosted component's binary64 ULP and even encode an OUTGOING
  // particle beside a positive unrounded diagnostic. Bound the effect on the
  // actual encoded momentum -> velocity map, with no subtractive v_lab-U
  // evaluation. If q_d=q+delta, g=sqrt(1+q^2), g_d=sqrt(1+q_d^2), exactly
  //   delta_g=(2q.dot.delta+delta^2)/(g_d+g),
  //   delta(v.n)/c=n.dot.delta/g_d-(n.dot.q)delta_g/(g_d*g).
  // Enclose each delta by SI narrowing plus the conditioned working-precision
  // construction/boost envelope. Include the full radial magnitude so a
  // cancellation in a rotated component cannot erase its construction error.
  // The returned gamma is an additional encoded view, so its own rounding is
  // checked too. Neither a bulk-speed floor nor max(1,speed) can hide loss of
  // the small intended crossing. Failure leaves the event/RNG uncommitted.
  Vector encoded,deltaBound;
  const Wide budget=prepared.arithmeticRelativeBudget;
  Wide normalDelta=0.0L,normalMomentumBound=0.0L,squaredNormDelta=0.0L;
  for(unsigned i=0;i<3;++i) {
    encoded[i]=static_cast<Wide>(sample.simulationMomentumKgMPerS[i])/prepared.massTimesC;
    const Wide scale=magnitude+std::abs(boostCoefficient*prepared.beta[i])+std::abs(lab[i])+std::abs(encoded[i]);
    deltaBound[i]=std::abs(encoded[i]-lab[i])+budget*scale;
    normalDelta+=std::abs(prepared.normal[i])*deltaBound[i];
    normalMomentumBound+=std::abs(prepared.normal[i]*lab[i]);
    squaredNormDelta+=2.0L*std::abs(lab[i])*deltaBound[i]+deltaBound[i]*deltaBound[i];
  }
  const Wide physicalGamma=std::hypot(1.0L,Norm(lab)),encodedGamma=std::hypot(1.0L,Norm(encoded));
  const Wide gammaDelta=squaredNormDelta/(physicalGamma+encodedGamma);
  const Wide momentumSpeedError=C*(normalDelta+normalMomentumBound*gammaDelta/physicalGamma)/encodedGamma;
  const Wide reportedGamma=sample.simulationGamma;
  const Wide gammaEncodingError=std::abs(reportedGamma-labGamma)+
      budget*prepared.bulkGamma*(restGamma+std::abs(betaDotMomentum));
  const Wide reportedSpeedError=C*(normalDelta+normalMomentumBound*gammaEncodingError/physicalGamma)/reportedGamma;
  const Wide scalarArithmeticError=budget*(inwardSpeed+C*prepared.aMagnitude*magnitude/labGamma);
  const Wide crossingError=std::max(momentumSpeedError,reportedSpeedError)+scalarArithmeticError;
  if(!std::isfinite(crossingError)||crossingError>CrossingRelativeResolution*inwardSpeed)
    return Unresolved("encoded boosted Juttner phase cannot resolve its positive thermal material crossing to relative tolerance 2e-12; no momentum repair or tail redraw");
  *output=sample;return Status::Ok();
}
} // namespace

Status EvaluateMaterialFlux(const MaterialBoundaryInput& input,FluxNormalization* output) {
  if(!output)return Invalid("material Juttner flux requires an output");
  Prepared prepared;auto status=Prepare(input,&prepared);if(!status.ok())return status;
  Wide meanFraction=0.0L,bound=0.0L;FluxNormalization result;
  status=MeanSpeed(prepared.theta,&meanFraction,&result.method,&bound);if(!status.ok())return status;
  const Wide flux=static_cast<Wide>(input.properNumberDensityPerM3)*prepared.aMagnitude*C*meanFraction/4.0L;
  if(!Representable(prepared.theta,&result.dimensionlessTemperature,input.properTemperatureK>0.0)||
      !Representable(C*meanFraction,&result.meanRestSpeedMPerS,input.properTemperatureK>0.0)||
      !Representable(prepared.aMagnitude,&result.transformedNormalMagnitude,true)||
      !Representable(flux,&result.numberFluxPerM2S,input.properNumberDensityPerM3>0.0&&input.properTemperatureK>0.0))
    return Unresolved("material Juttner SI normalization is not representable; no zero/floor substitution is allowed");
  // A positive mathematical error bound may be smaller than binary64. Round
  // that AUXILIARY upper bound upward (possibly to denorm_min), rather than
  // rejecting a perfectly representable cold physical rate or claiming zero
  // error. No physical flux, temperature or particle value is floored.
  result.asymptoticRelativeRemainderBound=static_cast<double>(bound);
  if(static_cast<Wide>(result.asymptoticRelativeRemainderBound)<bound)
    result.asymptoticRelativeRemainderBound=std::nextafter(result.asymptoticRelativeRemainderBound,std::numeric_limits<double>::infinity());
  *output=result;return Status::Ok();
}
Status SampleMaterialFluxFromUniforms(const MaterialBoundaryInput& input,
    const std::array<double,6>& uniforms,FluxSample* output) {
  Prepared prepared;auto status=Prepare(input,&prepared);if(!status.ok())return status;
  return DrawPrepared(input,prepared,uniforms,output);
}
Status SampleMaterialFlux(const MaterialBoundaryInput& input,Transport::KeyedRandomStream* random,FluxSample* output) {
  if(!random||!output)return Invalid("material Juttner keyed sample requires RNG and output");
  Prepared prepared;auto status=Prepare(input,&prepared);if(!status.ok())return status;
  if(!(input.properNumberDensityPerM3>0.0)||!(input.properTemperatureK>0.0))
    return Status::Error(StatusCode::OutOfDomain,"zero material flux has no conditional particle distribution");
  auto stagedRandom=*random;std::array<double,6> uniforms;
  for(auto& value:uniforms)value=stagedRandom.UniformOpen01();
  FluxSample candidate;status=DrawPrepared(input,prepared,uniforms,&candidate);if(!status.ok())return status;
  *random=stagedRandom;*output=candidate;return Status::Ok();
}
} } }
