#include "sep_shock_source_core.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>

namespace SEP {
namespace Shock {
namespace {

Transport::Status Error(const std::string& text) {
  return Transport::Status::Error(Transport::StatusCode::InvalidArgument, text);
}

double Dot(const Vector3& a, const Vector3& b) {
  return a.x*b.x + a.y*b.y + a.z*b.z;
}

Vector3 Subtract(const Vector3& a, const Vector3& b) {
  Vector3 value;
  value.x = a.x-b.x; value.y = a.y-b.y; value.z = a.z-b.z;
  return value;
}

Vector3 AddScaled(const Vector3& a, const Vector3& d, double fraction) {
  Vector3 value;
  value.x = a.x+fraction*d.x;
  value.y = a.y+fraction*d.y;
  value.z = a.z+fraction*d.z;
  return value;
}

double Length(const Vector3& a) { return std::sqrt(Dot(a,a)); }

// Time to cross one linear v(r)=a+b r interval.  The logarithmic formula is
// exact when b is nonzero and log1p prevents cancellation for a nearly flat
// interval.  Constant-speed intervals use the exact limiting expression.
double TravelTime(double r0, double r1, double v0, double v1) {
  if (r1 == r0) return 0.0;
  const double slope = (v1-v0)/(r1-r0);
  if (std::fabs(slope)*(r1-r0) <=
      32.0*std::numeric_limits<double>::epsilon()*std::fabs(v0))
    return (r1-r0)/v0;
  return std::log1p((v1-v0)/v0)/slope;
}

// Invert the preceding exact travel-time integral inside one interval.
double RadiusAfter(double r0, double v0, double slope, double dt) {
  if (std::fabs(slope*dt) < 32.0*std::numeric_limits<double>::epsilon())
    return r0 + v0*dt;
  return r0 + v0*std::expm1(slope*dt)/slope;
}

}  // namespace

Transport::Status ValidateTrajectoryConfiguration(
    const TrajectoryConfiguration& c) {
  if (!std::isfinite(c.launchEpochS) || !std::isfinite(c.launchRadiusM) ||
      c.launchRadiusM < 0.0 || c.knots.empty())
    return Error("shock trajectory requires a finite launch epoch/radius and knots");
  for (std::size_t i=0; i<c.knots.size(); ++i) {
    if (!std::isfinite(c.knots[i].radiusM) ||
        !std::isfinite(c.knots[i].speedMPerS) ||
        c.knots[i].radiusM < 0.0 || c.knots[i].speedMPerS <= 0.0 ||
        (i && c.knots[i].radiusM <= c.knots[i-1].radiusM))
      return Error("shock speed knots require increasing radii and positive speeds");
  }
  return Transport::Status::Ok();
}

Transport::ScalarResult SpeedAtRadius(const TrajectoryConfiguration& c,
                                      double radiusM) {
  Transport::ScalarResult result;
  result.status = ValidateTrajectoryConfiguration(c);
  if (!result.status.ok() || !std::isfinite(radiusM) || radiusM < 0.0) {
    if (result.status.ok()) result.status = Error("shock radius is invalid");
    return result;
  }
  if (radiusM <= c.knots.front().radiusM) result.value=c.knots.front().speedMPerS;
  else if (radiusM >= c.knots.back().radiusM) result.value=c.knots.back().speedMPerS;
  else {
    std::size_t upper=1;
    while (c.knots[upper].radiusM < radiusM) ++upper;
    const SpeedKnot& a=c.knots[upper-1];
    const SpeedKnot& b=c.knots[upper];
    result.value=a.speedMPerS+(b.speedMPerS-a.speedMPerS)*
        (radiusM-a.radiusM)/(b.radiusM-a.radiusM);
  }
  result.status=Transport::Status::Ok();
  return result;
}

Transport::Status StateAtEpoch(const TrajectoryConfiguration& c,
                               double epochS, TrajectoryState* state) {
  const Transport::Status valid=ValidateTrajectoryConfiguration(c);
  if (!valid.ok()) return valid;
  if (!state || !std::isfinite(epochS) || epochS < c.launchEpochS)
    return Error("shock epoch must be finite and not precede launch");
  state->configuration=c;
  state->epochS=c.launchEpochS;
  state->radiusM=c.launchRadiusM;
  return AdvanceToEpoch(state,epochS);
}

Transport::Status AdvanceToEpoch(TrajectoryState* state, double epochS) {
  if (!state) return Error("shock trajectory state is null");
  const Transport::Status valid=ValidateTrajectoryConfiguration(state->configuration);
  if (!valid.ok()) return valid;
  if (!std::isfinite(state->epochS) || !std::isfinite(state->radiusM) ||
      !std::isfinite(epochS) || epochS < state->epochS)
    return Error("shock trajectory time must be monotone");
  // Re-evaluate from the serialized launch datum rather than composing rounded
  // increments from the current radius.  The work is bounded by the small knot
  // table and guarantees bitwise timestep-partition/restart invariance.
  state->radiusM=state->configuration.launchRadiusM;
  double remaining=epochS-state->configuration.launchEpochS;
  while (remaining > 0.0) {
    const Transport::ScalarResult speed=SpeedAtRadius(state->configuration,
                                                       state->radiusM);
    if (!speed.status.ok()) return speed.status;
    const std::vector<SpeedKnot>& knots=state->configuration.knots;
    std::size_t upper=0;
    while (upper<knots.size() && knots[upper].radiusM<=state->radiusM) ++upper;
    if (upper==knots.size()) {
      state->radiusM += speed.value*remaining;
      remaining=0.0;
      continue;
    }
    const double nextRadius=knots[upper].radiusM;
    const double nextSpeed=knots[upper].speedMPerS;
    const double toKnot=TravelTime(state->radiusM,nextRadius,speed.value,nextSpeed);
    if (remaining < toKnot) {
      const double slope=(nextSpeed-speed.value)/(nextRadius-state->radiusM);
      state->radiusM=RadiusAfter(state->radiusM,speed.value,slope,remaining);
      remaining=0.0;
    }
    else {
      state->radiusM=nextRadius;
      remaining-=toKnot;
    }
  }
  state->epochS=epochS;
  return Transport::Status::Ok();
}

Transport::Status SerializeTrajectory(const TrajectoryState& state,
                                      std::string* text) {
  if (!text) return Error("shock checkpoint output is null");
  const Transport::Status valid=ValidateTrajectoryConfiguration(state.configuration);
  if (!valid.ok()) return valid;
  std::ostringstream out;
  out << "SEP_SHOCK_TRAJECTORY 1\n" << std::setprecision(17)
      << state.configuration.launchEpochS << ' '
      << state.configuration.launchRadiusM << ' '
      << state.epochS << ' ' << state.radiusM << '\n'
      << state.configuration.knots.size() << '\n';
  for (std::size_t i=0;i<state.configuration.knots.size();++i)
    out << state.configuration.knots[i].radiusM << ' '
        << state.configuration.knots[i].speedMPerS << '\n';
  *text=out.str();
  return Transport::Status::Ok();
}

Transport::Status DeserializeTrajectory(const std::string& text,
                                        TrajectoryState* state) {
  if (!state) return Error("shock checkpoint destination is null");
  std::istringstream in(text);
  std::string magic; int version=0; std::size_t count=0;
  TrajectoryState parsed;
  if (!(in>>magic>>version) || magic!="SEP_SHOCK_TRAJECTORY" || version!=1 ||
      !(in>>parsed.configuration.launchEpochS>>parsed.configuration.launchRadiusM
          >>parsed.epochS>>parsed.radiusM>>count) || count==0)
    return Error("shock checkpoint header is invalid");
  parsed.configuration.knots.resize(count);
  for (std::size_t i=0;i<count;++i)
    if (!(in>>parsed.configuration.knots[i].radiusM
            >>parsed.configuration.knots[i].speedMPerS))
      return Error("shock checkpoint knot table is truncated");
  std::string extra;
  if (in>>extra) return Error("shock checkpoint has trailing data");
  const Transport::Status valid=ValidateTrajectoryConfiguration(parsed.configuration);
  if (!valid.ok() || parsed.epochS<parsed.configuration.launchEpochS ||
      !std::isfinite(parsed.radiusM)) return Error("shock checkpoint state is invalid");
  *state=parsed;
  return Transport::Status::Ok();
}

IntersectionResult IntersectSphere(const std::vector<Vector3>& vertices,
                                   double radiusM, double toleranceM,
                                   IntersectionPolicy policy) {
  IntersectionResult result;
  if (vertices.size()<2 || !std::isfinite(radiusM) || radiusM<=0.0 ||
      !std::isfinite(toleranceM) || toleranceM<0.0) {
    result.status=IntersectionStatus::InvalidGeometry;
    result.message="sphere intersection requires a polyline, radius, and tolerance";
    return result;
  }
  double arc=0.0;
  for (std::size_t i=0;i+1<vertices.size();++i) {
    const Vector3 d=Subtract(vertices[i+1],vertices[i]);
    const double length=Length(d);
    if (!std::isfinite(length) || length<=toleranceM) {
      result.status=IntersectionStatus::InvalidGeometry;
      result.message="polyline contains a zero-length or non-finite segment";
      return result;
    }
    const double a=Dot(d,d), b=2.0*Dot(vertices[i],d);
    const double c=Dot(vertices[i],vertices[i])-radiusM*radiusM;
    double discriminant=b*b-4.0*a*c;
    const double discTolerance=8.0*a*radiusM*std::max(toleranceM,
        std::numeric_limits<double>::epsilon()*radiusM);
    if (discriminant>=-discTolerance) {
      if (discriminant<0.0) discriminant=0.0;
      const double root=std::sqrt(discriminant);
      double roots[2]={(-b-root)/(2.0*a),(-b+root)/(2.0*a)};
      const int n=discriminant==0.0 ? 1 : 2;
      for (int j=0;j<n;++j) {
        const double fractionTolerance=toleranceM/length;
        if (roots[j]<-fractionTolerance || roots[j]>1.0+fractionTolerance) continue;
        const double fraction=std::max(0.0,std::min(1.0,roots[j]));
        Intersection hit;
        hit.segment=i; hit.fraction=fraction;
        hit.positionM=AddScaled(vertices[i],d,fraction);
        hit.arcLengthM=arc+fraction*length;
        const double radialDerivative=Dot(hit.positionM,d);
        const double orientationTolerance=radiusM*length*
            16.0*std::numeric_limits<double>::epsilon();
        hit.orientation=std::fabs(radialDerivative)<=orientationTolerance
            ? CrossingOrientation::Tangent
            : (radialDerivative>0.0 ? CrossingOrientation::Outward
                                    : CrossingOrientation::Inward);
        if (!result.intersections.empty() &&
            std::fabs(result.intersections.back().arcLengthM-hit.arcLengthM)
                <=toleranceM) continue;
        result.intersections.push_back(hit);
      }
    }
    arc+=length;
  }
  if (result.intersections.empty()) {
    result.status=IntersectionStatus::NoIntersection;
    result.message="polyline does not intersect the sphere";
    return result;
  }
  if (policy==IntersectionPolicy::FirstOutward) {
    std::vector<Intersection> outward;
    for (std::size_t i=0;i<result.intersections.size();++i)
      if (result.intersections[i].orientation==CrossingOrientation::Outward) {
        outward.push_back(result.intersections[i]); break;
      }
    if (outward.empty()) {
      result.status=IntersectionStatus::Ambiguous;
      result.message="intersections exist but none is outward";
      return result;
    }
    result.intersections.swap(outward);
  }
  else if (policy==IntersectionPolicy::First) result.intersections.resize(1);
  result.status=IntersectionStatus::Ok;
  return result;
}

TurbulenceSourceEnergy ComputeTurbulenceSource(
    const TurbulenceSourceInput& in) {
  TurbulenceSourceEnergy result;
  if (!std::isfinite(in.sweptVolumeM3) || in.sweptVolumeM3<0.0 ||
      !std::isfinite(in.upstreamMassDensityKgPerM3) ||
      in.upstreamMassDensityKgPerM3<0.0 ||
      !std::isfinite(in.shockNormalSpeedMPerS) ||
      !std::isfinite(in.upstreamNormalSpeedMPerS) ||
      !std::isfinite(in.efficiency) || in.efficiency<0.0 || in.efficiency>1.0 ||
      !std::isfinite(in.plusBranchFraction) || in.plusBranchFraction<0.0 ||
      in.plusBranchFraction>1.0) {
    result.status=Error("shock turbulence source contains invalid physical inputs");
    return result;
  }
  const double relative=std::max(0.0,
      in.shockNormalSpeedMPerS-in.upstreamNormalSpeedMPerS);
  result.totalJ=0.5*in.efficiency*in.upstreamMassDensityKgPerM3*
      in.sweptVolumeM3*relative*relative;
  result.plusJ=result.totalJ*in.plusBranchFraction;
  result.minusJ=result.totalJ-result.plusJ;
  result.status=std::isfinite(result.totalJ) ? Transport::Status::Ok()
      : Error("shock turbulence source overflowed");
  return result;
}

}  // namespace Shock
}  // namespace SEP
