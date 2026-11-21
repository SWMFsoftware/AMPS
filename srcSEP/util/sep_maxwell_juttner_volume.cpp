#include "sep_maxwell_juttner_volume.h"
#include "sep_maxwell_juttner_flux.h"
#include <algorithm>
#include <cmath>
#include <limits>

namespace SEP { namespace Sources { namespace MaxwellJuttnerVolume {
namespace {
using Wide=long double;
constexpr Wide C=MaxwellJuttner::SpeedOfLightMPerS;
constexpr Wide K=MaxwellJuttner::BoltzmannJPerK;
Transport::Status Invalid(const char* message){return Transport::Status::Error(Transport::StatusCode::InvalidArgument,message);}
Transport::Status Unresolved(const char* message){return Transport::Status::Error(Transport::StatusCode::UnresolvedCoefficient,message);}
}
Transport::Status BulkGamma(const std::array<double,3>& velocity,double* output) {
  if(!output)return Invalid("thermal volume source requires a Lorentz-factor output");
  for(double value:velocity)if(!std::isfinite(value))return Invalid("thermal volume bulk velocity is nonfinite");
  const Wide speed=std::hypot(static_cast<Wide>(velocity[0]),std::hypot(static_cast<Wide>(velocity[1]),static_cast<Wide>(velocity[2])));
  if(!(speed<C))return Invalid("thermal volume bulk velocity must be subluminal");
  const Wide beta=speed/C,gamma=1/std::sqrt((1-beta)*(1+beta));
  // This bounds arithmetic conditioning, not an assumed physical speed range.
  // An unresolved gamma is rejected rather than clipped to a convenient bulk
  // speed, because gamma also determines the simulation-volume number count.
  if(64*std::numeric_limits<Wide>::epsilon()*gamma*gamma>2e-12L)
    return Unresolved("thermal volume bulk gamma is numerically unresolved");
  *output=static_cast<double>(gamma);return Transport::Status::Ok();
}
Transport::Status Draw(const Input& input,Transport::KeyedRandomStream* random,Sample* output) {
  if(!random||!output||!std::isfinite(input.properTemperatureK)||input.properTemperatureK<0||
     !std::isfinite(input.restMassKg)||input.restMassKg<=0)
    return Invalid("thermal volume sampling requires nonnegative proper T and positive rest mass");
  Sample candidate;auto status=BulkGamma(input.bulkVelocityMPerS,&candidate.bulkGamma);if(!status.ok())return status;
  if(input.properTemperatureK==0){candidate.simulationGamma=candidate.bulkGamma;*output=candidate;return Transport::Status::Ok();}
  const Wide massC=static_cast<Wide>(input.restMassKg)*C;
  const Wide theta=K*input.properTemperatureK/(massC*C);
  if(!(theta>0)||!std::isfinite(theta)||static_cast<double>(theta)==0)
    return Unresolved("thermal volume temperature is outside numerical range");
  auto staged=*random;
  // For q=gamma'-1, the rest-volume radial density is
  // (q+1)*sqrt(q*(q+2))*exp(-q/theta). The exact envelope
  // (q+1)*(q+sqrt(2*q))*exp(-q/theta) is a positive mixture of
  // Gamma(3), Gamma(2), Gamma(5/2), Gamma(3/2) laws. Work with log mixture
  // weights to preserve both cold and ultrarelativistic limits without a
  // momentum cutoff, classical switch, Bessel underflow or missing tail.
  const Wide shapes[4]={1.5L,2.0L,2.5L,3.0L};Wide weights[4],maximum=-std::numeric_limits<Wide>::infinity();
  for(int i=0;i<4;++i){weights[i]=std::lgamma(shapes[i])+shapes[i]*std::log(theta)+(i==0||i==2?0.5L*std::log(2.0L):0.0L);maximum=std::max(maximum,weights[i]);}
  Wide total=0;for(auto& weight:weights){weight=std::exp(weight-maximum);total+=weight;}
  const Wide bulkSpeed=std::hypot(static_cast<Wide>(input.bulkVelocityMPerS[0]),
      std::hypot(static_cast<Wide>(input.bulkVelocityMPerS[1]),static_cast<Wide>(input.bulkVelocityMPerS[2])));
  for(unsigned attempt=0;attempt<4096;++attempt) {
    Wide select=staged.UniformOpen01()*total;int component=3;
    for(int i=0;i<3;++i){if(select<weights[i]){component=i;break;}select-=weights[i];}
    Wide gammaVariate=0;
    const int exponentials=component==3?3:(component==0?1:2);
    for(int i=0;i<exponentials;++i)gammaVariate-=std::log(static_cast<Wide>(staged.UniformOpen01()));
    if(component==0||component==2){const Wide normal=staged.Normal01();gammaVariate+=0.5L*normal*normal;}
    const Wide q=theta*gammaVariate;
    if(!(q>0)||!std::isfinite(q)||!std::isfinite(static_cast<double>(q))||static_cast<double>(q)==0)
      return Unresolved("thermal volume sampled energy is unrepresentable; no tail redraw");
    const Wide radialAcceptance=std::sqrt(q+2)/(std::sqrt(q)+std::sqrt(2.0L));
    if(staged.UniformOpen01()>radialAcceptance)continue;
    const Wide mu=2*staged.UniformOpen01()-1,phi=2*std::acos(-1.0L)*staged.UniformOpen01();
    const Wide sine=std::sqrt((1-mu)*(1+mu)),magnitude=std::sqrt(q)*std::sqrt(q+2);
    const Wide direction[3]={mu,sine*std::cos(phi),sine*std::sin(phi)};
    Wide betaDotMomentum=0;
    for(int d=0;d<3;++d)betaDotMomentum+=input.bulkVelocityMPerS[d]/C*magnitude*direction[d];
    // d^3p = gamma_U*(1+U.v'/c^2)*d^3p'. Consequently a snapshot of a
    // simulation-frame spatial volume is not an unweighted rest-volume draw.
    // This second exact rejection supplies its Jacobian; number density is
    // gamma_U*n0. It is independent of boundary crossing/flux weighting.
    const Wide jacobian=1+betaDotMomentum/(1+q);
    if(!(jacobian>0)||!std::isfinite(jacobian))return Unresolved("thermal volume phase Jacobian is unresolved");
    if(staged.UniformOpen01()*(1+bulkSpeed/C)>jacobian)continue;
    for(int d=0;d<3;++d){const Wide momentum=massC*magnitude*direction[d];
      const double encoded=static_cast<double>(momentum);
      if(!std::isfinite(encoded)||(momentum!=0&&encoded==0))return Unresolved("thermal volume momentum encoding is unresolved");
      candidate.restMomentumKgMPerS[d]=encoded;}
    candidate.restKineticGamma=static_cast<double>(q);
    candidate.simulationGamma=static_cast<double>(candidate.bulkGamma*((1+q)+betaDotMomentum));
    if(!std::isfinite(candidate.simulationGamma)||candidate.simulationGamma<1)
      return Unresolved("thermal volume Lorentz energy is unresolved");
    // Legitimate envelope/Jacobian rejections belong to this distribution.
    // Arithmetic failures never cause a tail redraw. Commit the complete RNG
    // state, including its normal cache, only with a representable sample.
    *random=staged;*output=candidate;return Transport::Status::Ok();
  }
  return Unresolved("thermal volume exact rejection sampler exhausted its declared attempt budget");
}
} } }
