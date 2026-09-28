#include "sep_coronal_cme/interface_balance.h"
#include "sep_coronal_cme/constants.h"

#include <algorithm>
#include <cmath>

namespace SEP { namespace CoronalCME { namespace {
Vec3 Traction(const PrimitiveState& s,Vec3 n,double wn) {
  // Full vector traction is required: checking only total normal pressure
  // would miss tangential momentum and magnetic-stress residuals.
  const double magnetic=Dot(s.magneticFieldT,s.magneticFieldT)/(2*Constants::kVacuumPermeabilityHPerM);
  const double bn=Dot(s.magneticFieldT,n);
  return s.densityKgM3*wn*s.velocityMPerS+(s.pressurePa+magnetic)*n-
      (bn/Constants::kVacuumPermeabilityHPerM)*s.magneticFieldT;
}
bool Valid(const PrimitiveState& s){return s.densityKgM3>0&&s.pressurePa>0&&s.fastSpeedMPerS>0;}
}
Core::Result<SharpBalance> EvaluateSharpInterface(InterfaceKind kind,
    BalancePolicy policy,StateOrigin originA,StateOrigin originB,
    const PrimitiveState& a,const PrimitiveState& b,Vec3 normal,double vi,
    const InterfaceTolerances& t,bool evidence) {
  normal=Unit(normal);if(Norm(normal)==0||!Valid(a)||!Valid(b))return Core::Result<SharpBalance>::Failure(
      Core::StatusCode::InvalidConfiguration,"sharp interface state/normal is invalid");
  if(policy==BalancePolicy::StationaryTangentialDiscontinuity&&
      (kind!=InterfaceKind::OpenClosedSeparatrix||
       originA==StateOrigin::AnalyticComposite||originB==StateOrigin::AnalyticComposite))
    return Core::Result<SharpBalance>::Failure(Core::StatusCode::InvalidConfiguration,
        "stationary TD requires solver/imported open-closed state");
  if(policy==BalancePolicy::BoundedApproximation&&!evidence)
    return Core::Result<SharpBalance>::Failure(Core::StatusCode::DataIntegrityFailure,
        "bounded approximation lacks uncertainty/three-level convergence evidence");
  const double wa=Dot(a.velocityMPerS,normal)-vi,wb=Dot(b.velocityMPerS,normal)-vi;
  const double bna=Dot(a.magneticFieldT,normal),bnb=Dot(b.magneticFieldT,normal);
  const Vec3 jump=Traction(b,normal,wb)-Traction(a,normal,wa);
  SharpBalance result;result.tractionJumpPa=jump;result.tractionNormPa=Norm(jump);
  result.signedMassFluxJumpKgM2S=b.densityKgM3*wb-a.densityKgM3*wa;
  if(kind==InterfaceKind::OpenClosedSeparatrix){result.magneticGatePassed=
      std::abs(bna)<=t.magneticAbsoluteT+t.magneticRelative*Norm(a.magneticFieldT)&&
      std::abs(bnb)<=t.magneticAbsoluteT+t.magneticRelative*Norm(b.magneticFieldT);
    result.kinematicGatePassed=std::abs(wa)<=t.speedAbsoluteMPerS+t.speedMach*a.fastSpeedMPerS&&
      std::abs(wb)<=t.speedAbsoluteMPerS+t.speedMach*b.fastSpeedMPerS;
  }else{result.magneticGatePassed=true;result.kinematicGatePassed=true;}
  const double massScale=std::max(t.massFluxFloor,std::max(std::abs(a.densityKgM3*wa),std::abs(b.densityKgM3*wb)));
  result.massFluxGatePassed=std::abs(result.signedMassFluxJumpKgM2S)<=t.massFluxAbsolute+t.massFluxRelative*massScale;
  const double tractionScale=std::max(t.tractionFloorPa,std::max(Norm(Traction(a,normal,wa)),Norm(Traction(b,normal,wb))));
  result.tractionGatePassed=result.tractionNormPa<=t.tractionAbsolutePa+t.tractionRelative*tractionScale;
  result.policyPassed=result.magneticGatePassed&&result.kinematicGatePassed&&result.massFluxGatePassed&&
      (policy==BalancePolicy::DiagnosticKinematic||result.tractionGatePassed);
  return Core::Result<SharpBalance>::Success(result);
}
Core::Result<VolumeBalance> EvaluateVolumeMomentumResidual(Vec3 time,Vec3 divergence,
    Vec3 gravity,Vec3 declared) {
  // This volume force-density is never compared to or substituted by a sharp
  // surface traction. Its dimensions are N/m^3 rather than Pa.
  const Vec3 residual=time+divergence-gravity-declared;
  if(!(std::isfinite(residual.x)&&std::isfinite(residual.y)&&std::isfinite(residual.z)))
    return Core::Result<VolumeBalance>::Failure(Core::StatusCode::NumericalFailure,
        "volume momentum residual is nonfinite");
  return Core::Result<VolumeBalance>::Success({residual,Norm(residual)});
}
} }
