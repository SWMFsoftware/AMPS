#include "sep_absorbing_bridge.h"
#include <algorithm>
#include <cmath>

namespace SEP { namespace Transport {
Status RefineBoundaryInterval(double begin,double proposed,double minimum,
    const BoundaryIntervalCheck& check,double* accepted,BoundaryIntervalControl* history) {
  if(!std::isfinite(begin)||!std::isfinite(proposed)||proposed<=0||
      !std::isfinite(minimum)||minimum<=0||!check||!accepted||
      (history&&(!std::isfinite(history->previousAcceptedS)||history->previousAcceptedS<0||
                 history->successfulSinceGrowth>8)))
    return Status::Error(StatusCode::InvalidArgument,"invalid predictable boundary-refinement control");
  double interval=proposed;
  bool growthProbe=false;
  if(history&&history->previousAcceptedS>0&&history->previousAcceptedS<proposed) {
    interval=history->previousAcceptedS;
    // A hint selects ONLY the first proposal; the callback still grades every
    // quarter/mid/end query using its original SI error bound before RNG draws.
    // Probe a factor-two increase after eight checked intervals, rather than
    // paying for the same rejected coarse interval at every particle substep.
    // Bounded growth also recovers after a short interval near a mesh vertex.
    if(history->successfulSinceGrowth==8) {
      interval=interval<proposed/2?2*interval:proposed;
      growthProbe=true;
    }
  }
  for(unsigned attempt=0;attempt<64;++attempt) {
    bool resolved=false;const auto status=check(interval,&resolved);
    if(!status.ok())return status;
    if(resolved) {
      // Publish history only after a fully checked successful preparation.
      // Failed queries leave both the caller's interval and hint unchanged.
      *accepted=interval;
      if(history){history->previousAcceptedS=interval;history->successfulSinceGrowth=attempt||growthProbe?0:
          std::min(8U,history->successfulSinceGrowth+1);}
      return Status::Ok();
    }
    interval*=0.5;
    if(interval<minimum||begin+interval==begin)
      return Status::Error(StatusCode::StepUnderflow,
          "moving shock boundary tolerance requires an unrepresentable internal interval; physical AMPS dt remains unchanged");
  }
  return Status::Error(StatusCode::StepUnderflow,
      "moving shock boundary refinement exhausted its declared 64-attempt work limit");
}
AbsorbingBridgeResult AbsorbingLinearBridge(double a,double b,double kdt,RandomStream& random) {
  AbsorbingBridgeResult out;out.status=Status::Ok();
  if(!std::isfinite(a)||!std::isfinite(b)||!std::isfinite(kdt)||kdt<0) {
    out.status=Status::Error(StatusCode::InvalidArgument,"absorbing bridge requires finite distances and nonnegative kappa*dt");return out;
  }
  if(a<=0){out.hit=true;out.fraction=0;return out;}
  if(kdt==0){out.hit=b<=0;if(out.hit)out.fraction=a/(a-b);return out;}
  if(b>0) {
    const double u=random.UniformOpen01();
    if(!(u>0&&u<1)){out.status=Status::Error(StatusCode::InvalidParticleState,"bridge random stream violated open-unit support");return out;}
    // exp(-a*b/(kappa*dt)) is the conditional hidden-crossing probability.
    // Compare logarithms: underflow of a tiny probability cannot silently
    // discard a representable rare event or require a probability floor.
    if(std::log(u)>=-(a/std::sqrt(kdt))*(b/std::sqrt(kdt)))return out;
  }
  out.hit=true;
  const double z=random.Normal01();
  if(!std::isfinite(z)){out.status=Status::Error(StatusCode::InvalidParticleState,"bridge normal draw is not finite");return out;}
  // Given a hit, u=t_hit/(dt-t_hit) is inverse Gaussian with mean a/|b|
  // and shape a^2/(2*kappa*dt). Sample its Michael-Schucany-Haas mixture in
  // scaled variables. The rationalized small root avoids catastrophic
  // cancellation; the b=0 Levy limit needs no arbitrary endpoint floor.
  const double scale=std::max(a,std::max(std::abs(b),std::sqrt(kdt)));
  const double x=a/scale,y=std::abs(b)/scale,v=(std::sqrt(kdt)/scale)*std::abs(z);
  const double d=v*v,xy=x*y;
  const double denominator=xy+d+std::sqrt(d*(2*xy+d));
  if(!std::isfinite(denominator)){out.status=Status::Error(StatusCode::InvalidParticleState,"bridge first-passage sample is not representable");return out;}
  bool small=true;
  if(y>0) {
    const double u=random.UniformOpen01();
    if(!(u>0&&u<1)){out.status=Status::Error(StatusCode::InvalidParticleState,"bridge mixture draw violated open-unit support");return out;}
    small=u<=denominator/(denominator+xy);
  }
  out.fraction=small?x*x/(x*x+denominator):denominator/(denominator+y*y);
  if(!std::isfinite(out.fraction)||out.fraction<0||out.fraction>1)
    out.status=Status::Error(StatusCode::InvalidParticleState,"invalid conditional first-passage time");
  return out;
}
} }
