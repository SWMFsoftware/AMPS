#ifndef SRCSEP_UTIL_SEP_BACKGROUND_DERIVATIVES_H
#define SRCSEP_UTIL_SEP_BACKGROUND_DERIVATIVES_H

#include "sep_common_header_path.h"
#include SRCSEP_SEP_COMMON_HEADER(sep_transport_common.h)
#include <array>
#include <cmath>

namespace SEP { namespace Transport {

// G_ij = partial U_i / partial x_j [1/s], in one Cartesian basis.
// A derivative along a curved field of U dot b contains U dot db/ds;
// it is therefore not the strain b_i G_ij b_j used by focused transport.
struct VelocityDerivatives {
  std::array<double,9> gradientPerS{{0,0,0,0,0,0,0,0,0}};
  double divergencePerS = 0.0;
  double fieldAlignedStrainPerS = 0.0;
};

// The PIC field-line particle coordinate and injected v_parallel increase
// along the geometric tangent, independent of magnetic polarity.
inline ScalarResult FieldLineAdvectionMPerS(const double velocity[3],
                                            const double tangent[3]) {
  ScalarResult result;
  double norm2=0.0,projection=0.0;
  if(!velocity||!tangent) {
    result.status=Status::Error(StatusCode::InvalidArgument,
        "field-line advection requires velocity and geometric tangent");
    return result;
  }
  for(int i=0;i<3;++i) {
    if(!std::isfinite(velocity[i])||!std::isfinite(tangent[i])) {
      result.status=Status::Error(StatusCode::InvalidArgument,
          "field-line advection vectors must be finite");
      return result;
    }
    norm2+=tangent[i]*tangent[i];projection+=velocity[i]*tangent[i];
  }
  if(!(norm2>0.0)||!std::isfinite(norm2)||!std::isfinite(projection)) {
    result.status=Status::Error(StatusCode::InvalidArgument,
        "field-line advection requires a nonzero finite tangent");
    return result;
  }
  result.value=projection/std::sqrt(norm2);
  result.status=Status::Ok();
  return result;
}

inline ScalarResult MagneticPolarityAlongFieldLine(const double magnetic[3],
                                                    const double tangent[3]) {
  ScalarResult result=FieldLineAdvectionMPerS(magnetic,tangent);
  if(!result.status.ok())return result;
  if(result.value==0.0) {
    result.status=Status::Error(StatusCode::InvalidParticleState,
        "magnetic field has no orientation along the assigned field line");
    return result;
  }
  result.value=result.value>0.0?1.0:-1.0;
  return result;
}

inline Status ContractVelocityGradient(
    const std::array<double,9>& gradient, const double b[3],
    VelocityDerivatives* result) {
  if (!result || !b)
    return Status::Error(StatusCode::InvalidArgument,
                         "velocity derivative output and direction are required");
  double norm2 = 0.0;
  for (int i=0;i<3;++i) {
    if (!std::isfinite(b[i]))
      return Status::Error(StatusCode::InvalidArgument,
                           "velocity strain direction is not finite");
    norm2 += b[i]*b[i];
  }
  if (!(norm2>0.0) || !std::isfinite(norm2))
    return Status::Error(StatusCode::InvalidArgument,
                         "velocity strain direction must be nonzero");
  VelocityDerivatives candidate;
  candidate.gradientPerS=gradient;
  for(int i=0;i<3;++i)for(int j=0;j<3;++j) {
    const double value=gradient[3*i+j];
    if(!std::isfinite(value))
      return Status::Error(StatusCode::InvalidArgument,
                           "velocity gradient contains a nonfinite component");
    candidate.fieldAlignedStrainPerS+=b[i]*value*b[j]/norm2;
    if(i==j)candidate.divergencePerS+=value;
  }
  if(!std::isfinite(candidate.divergencePerS)||
      !std::isfinite(candidate.fieldAlignedStrainPerS))
    return Status::Error(StatusCode::InvalidArgument,
                         "velocity gradient contraction overflowed");
  *result=candidate;
  return Status::Ok();
}

// Exact tensor for the declared radial field U(x)=u(r) x/r.  This includes
// both transverse geometric expansion terms, even when du/dr is zero.
inline Status RadialVelocityGradient(const double positionM[3],
    double speedMPerS, double radialDerivativePerS,
    std::array<double,9>* gradient) {
  if(!positionM||!gradient||!std::isfinite(speedMPerS)||
      !std::isfinite(radialDerivativePerS))
    return Status::Error(StatusCode::InvalidArgument,
                         "radial velocity requires finite physical inputs");
  double radius2=0.0;
  for(int i=0;i<3;++i) {
    if(!std::isfinite(positionM[i]))
      return Status::Error(StatusCode::InvalidArgument,
                           "radial velocity position is not finite");
    radius2+=positionM[i]*positionM[i];
  }
  if(!(radius2>0.0)||!std::isfinite(radius2))
    return Status::Error(StatusCode::OutOfDomain,
                         "radial velocity requires positive finite radius");
  const double radius=std::sqrt(radius2);
  const double transverse=speedMPerS/radius;
  std::array<double,9> candidate;
  for(int i=0;i<3;++i)for(int j=0;j<3;++j)
    candidate[3*i+j]=(i==j?transverse:0.0)+
        (radialDerivativePerS-transverse)*
        (positionM[i]/radius)*(positionM[j]/radius);
  for(double value:candidate)if(!std::isfinite(value))
    return Status::Error(StatusCode::InvalidArgument,
                         "radial velocity gradient overflowed");
  *gradient=candidate;
  return Status::Ok();
}

} }
#endif
