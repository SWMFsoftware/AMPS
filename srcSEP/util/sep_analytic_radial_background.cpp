#include "sep_analytic_radial_background.h"
#include "sep_background_derivatives.h"

#include <cmath>
#include <iomanip>
#include <locale>
#include <sstream>

namespace SEP { namespace Background { namespace AnalyticRadial {

Transport::Status Evaluate(const Wind& wind,const double positionM[3],Sample* output) {
  using Transport::Status;using Transport::StatusCode;
  if(!positionM||!output||!std::isfinite(wind.speedMPerS)||wind.speedMPerS<=0)
    return Status::Error(StatusCode::InvalidArgument,
        "analytic radial wind requires an explicit finite positive speed, point and output in SI");
  double relative[3];
  for(unsigned dimension=0;dimension<3;++dimension) {
    if(!std::isfinite(wind.originM[dimension])||!std::isfinite(positionM[dimension]))
      return Status::Error(StatusCode::InvalidArgument,
          "analytic radial point and declared origin must be finite in one Cartesian frame");
    relative[dimension]=positionM[dimension]-wind.originM[dimension];
    if(!std::isfinite(relative[dimension]))
      return Status::Error(StatusCode::UnresolvedCoefficient,
          "analytic radial displacement is not representable; no origin substitution is allowed");
  }
  Sample candidate;
  candidate.radiusM=std::hypot(std::hypot(relative[0],relative[1]),relative[2]);
  if(!std::isfinite(candidate.radiusM)||candidate.radiusM<=0)
    return Status::Error(StatusCode::OutOfDomain,
        "analytic radial field is undefined at its declared origin or an unrepresentable radius");
  for(unsigned dimension=0;dimension<3;++dimension) {
    // Normalize before multiplying by speed; v*x would overflow for some
    // finite SI inputs although the resulting bounded velocity is finite.
    candidate.velocityMPerS[dimension]=wind.speedMPerS*(relative[dimension]/candidate.radiusM);
    if(!std::isfinite(candidate.velocityMPerS[dimension]))
      return Status::Error(StatusCode::UnresolvedCoefficient,"analytic radial velocity is not representable");
  }
  // Reuse the common exact radial tensor with du/dr=0. The translated point
  // is passed explicitly: global zero is not the origin of this declaration.
  const auto derivatives=Transport::RadialVelocityGradient(relative,wind.speedMPerS,0,&candidate.gradientPerS);
  if(!derivatives.ok())return derivatives;
  // A nonzero expansion lost to underflow cannot be advertised as a static
  // incompressible field. This is a numerical-domain rejection, not a floor.
  const double expansionScale=wind.speedMPerS/candidate.radiusM;
  if(!std::isfinite(expansionScale)||expansionScale<=0)
    return Status::Error(StatusCode::UnresolvedCoefficient,
        "analytic radial expansion scale is not representable");
  *output=candidate;
  return Status::Ok();
}

Transport::Status ValidateVelocity(const Wind& wind,const double positionM[3],
    const double velocityMPerS[3],double relativeTolerance) {
  using Transport::Status;using Transport::StatusCode;
  if(!velocityMPerS||!std::isfinite(relativeTolerance)||relativeTolerance<0)
    return Status::Error(StatusCode::InvalidArgument,
        "analytic radial endpoint comparison requires velocity and a finite nonnegative numerical tolerance");
  Sample expected;
  const auto status=Evaluate(wind,positionM,&expected);
  if(!status.ok())return status;
  for(unsigned dimension=0;dimension<3;++dimension) {
    if(!std::isfinite(velocityMPerS[dimension]))
      return Status::Error(StatusCode::InvalidParticleState,"native analytic radial velocity is not finite");
    if(std::fabs(velocityMPerS[dimension]-expected.velocityMPerS[dimension])/wind.speedMPerS>relativeTolerance)
      return Status::Error(StatusCode::UnsupportedConfiguration,
          "native velocity differs from the configured-origin/speed analytic radial declaration; an explicit provider is required");
  }
  return Status::Ok();
}

std::string CanonicalDeclaration(const Wind& wind) {
  std::ostringstream text;text.imbue(std::locale::classic());
  text<<std::scientific<<std::setprecision(17)<<"origin_m=";
  for(unsigned dimension=0;dimension<3;++dimension) {
    if(dimension)text<<',';
    // Signed zero denotes the same Cartesian origin, so presentation-level
    // sign bits must not give that identical physical declaration another ID.
    text<<(wind.originM[dimension]==0?0.0:wind.originM[dimension]);
  }
  text<<";speed_m_per_s="<<wind.speedMPerS;
  return text.str();
}

} } }
