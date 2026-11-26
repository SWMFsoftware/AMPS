// ============================================================================
// srcSEP3D/adapters/reduced_front_shock_surface.cpp
//
// ShockSurface implementation over the shared reduced shock-front core.  See
// the header for the geometry (generating sphere), the frame conversion and
// the consumers.  Every method converts frame points to HCI with
// frame_.ToInertialPosition(x - originM_, t), evaluates the core in HCI, and
// converts results back with FromInertial*(...) (+ originM_ for points).
// ============================================================================

#include "reduced_front_shock_surface.h"

#include "provider.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace SEP3D {
namespace Adapters {
namespace {

namespace SF = SEP::CoronaSwcme::ShockFront;

// Every failure of this adapter, including a core failure passed through
// with its message, is reported as InvalidInput.
Core::Status Invalid(const std::string& message) {
  return Core::Status(Core::StatusCode::InvalidInput, message);
}
// Component copies between the core's and srcSEP3D's vector types.
Core::Vec3 FromCore(const SEP::CoronalCME::Vec3& value) {
  return {value.x, value.y, value.z};
}
SEP::CoronalCME::Vec3 ToCore(const Core::Vec3& value) {
  return {value.x, value.y, value.z};
}
// All three components finite (rejects NaN and infinite query points).
bool Finite(const Core::Vec3& v) {
  return std::isfinite(v.x) && std::isfinite(v.y) && std::isfinite(v.z);
}

}  // namespace

Core::Status ReducedFrontShockSurface::Create(
    std::shared_ptr<const SF::Provider> provider,
    const Core::ComputationFrame& frame, const Core::Vec3& originM,
    ShockSurfacePointer* output) {
  if (!output || !provider)
    return Invalid("reduced-front shock surface needs a provider and output");
  const Core::Status valid = frame.Validate();
  if (!valid.ok()) return valid;
  // The frame conversion below assumes the core's vectors are HCI.
  if (provider->Event().coordinateFrame != "HCI")
    return Invalid("reduced-front shock surface requires an HCI event");
  std::shared_ptr<ReducedFrontShockSurface> surface(new ReducedFrontShockSurface);
  surface->provider_ = std::move(provider);
  surface->frame_ = frame;
  surface->originM_ = originM;
  *output = surface;
  return Core::Status::OK();
}

bool ReducedFrontShockSurface::Spherical() const {
  return provider_->Event().geometryModel ==
         SF::GeometryModel::HeliocentricSphere;
}

std::string ReducedFrontShockSurface::GeometryName() const {
  return Spherical() ? "heliocentric-sphere" : "finite-sse-fixed-width";
}
std::string ReducedFrontShockSurface::Identity() const {
  return provider_->Event().physicsFingerprint;
}
double ReducedFrontShockSurface::SupportBeginS() const {
  return provider_->Event().ambient.support.startS;
}
double ReducedFrontShockSurface::SupportEndS() const {
  return provider_->Event().ambient.support.endS;
}

Core::Status ReducedFrontShockSurface::GeneratingSphere(
    double epochS, Core::Vec3* centreHci, double* radiusM) const {
  // Heliocentric sphere: centre 0, radius R_sh.  SSE cap: centre
  // R_apex/(1 + sin lambda) along the CME axis, radius |c| sin lambda
  // (Provider::Sphere in src/models/sep_corona_swcme/shock_front).
  const auto sphere = provider_->Sphere(epochS);
  if (!sphere.ok()) return Invalid(sphere.status.message);
  *centreHci = FromCore(sphere.value.centerM);
  *radiusM = sphere.value.radiusM;
  return Core::Status::OK();
}

Core::Status ReducedFrontShockSurface::Level(const Core::Vec3& positionM,
                                             double epochS,
                                             double* level) const {
  if (!level || !Finite(positionM)) return Invalid("invalid level query");
  Core::Vec3 centre;
  double radius = 0.0;
  const Core::Status sphere = GeneratingSphere(epochS, &centre, &radius);
  if (!sphere.ok()) return sphere;
  const Core::Vec3 x = frame_.ToInertialPosition(positionM - originM_, epochS);
  // The heliocentric sphere: the exact signed distance.  The SSE cap: the
  // generating-sphere distance, a first-order distance near the cap; points
  // inside the generating sphere are the cap's geometric interior.
  *level = (x - centre).Norm() - radius;
  return Core::Status::OK();
}

Core::Status ReducedFrontShockSurface::Intersect(
    const Core::Vec3& x0, const Core::Vec3& x1, double epochS,
    std::vector<ShockSegmentCrossing>* crossings) const {
  if (!crossings || !Finite(x0) || !Finite(x1))
    return Invalid("invalid segment query");
  crossings->clear();
  Core::Vec3 centre;
  double radius = 0.0;
  const Core::Status sphere = GeneratingSphere(epochS, &centre, &radius);
  if (!sphere.ok()) return sphere;
  const Core::Vec3 a = frame_.ToInertialPosition(x0 - originM_, epochS);
  const Core::Vec3 b = frame_.ToInertialPosition(x1 - originM_, epochS);
  const Core::Vec3 v = b - a, r = a - centre;
  const double qa = v.Dot(v), qb = r.Dot(v), qc = r.Dot(r) - radius * radius;
  if (!(qa > 0.0)) return Core::Status::OK();
  const double discriminant = qb * qb - qa * qc;
  if (discriminant <= 0.0) return Core::Status::OK();  // miss or tangent
  // The segment x(s) = a + s v, s in [0,1], meets |x - c| = R where
  //   qa s^2 + 2 qb s + qc = 0,  qa = v.v, qb = r.v, qc = r.r - R^2, r = a - c,
  // with the reduced discriminant D = qb^2 - qa qc.  The roots are
  // s = (-qb -+ sqrt(D)) / qa.
  // Rationalized roots, as in the 1-D release-surface intersection
  // (srcSEP/util/sep_reduced_source.cpp): q = -(qb + sign(qb) sqrt(D)) never
  // subtracts nearly equal numbers, and the roots are q/qa and qc/q (their
  // product is qc/qa).  q = 0 would need qb = 0 and D = 0, which the D <= 0
  // return above already excludes; the guard only keeps the division safe.
  const double q = -qb - std::copysign(std::sqrt(discriminant), qb);
  double roots[2] = {q / qa, q == 0.0 ? -qb / qa : qc / q};
  if (roots[1] < roots[0]) std::swap(roots[0], roots[1]);
  for (double s : roots) {
    if (!(s >= 0.0 && s <= 1.0)) continue;
    const Core::Vec3 pointHci = a + v * s;
    ShockSegmentCrossing crossing;
    crossing.fraction = s;
    // Moving inward through the generating sphere enters the interior: the
    // outward radial (pointHci - centre) and the segment direction v have a
    // negative dot product there.
    crossing.direction = (pointHci - centre).Dot(v) < 0.0 ? +1 : -1;
    const Core::Vec3 pointFrame =
        frame_.FromInertialPosition(pointHci, epochS) + originM_;
    const Core::Status front = FrontAt(pointFrame, epochS, &crossing.front);
    if (!front.ok()) return front;
    // The SSE generating sphere's rear and flanks are not the front.
    if (!crossing.front.supported) continue;
    crossings->push_back(crossing);
  }
  return Core::Status::OK();
}

Core::Status ReducedFrontShockSurface::FrontAt(const Core::Vec3& positionM,
                                               double epochS,
                                               ShockFrontState* state) const {
  if (!state || !Finite(positionM)) return Invalid("invalid front query");
  *state = ShockFrontState{};
  const Core::Vec3 x = frame_.ToInertialPosition(positionM - originM_, epochS);
  const double r = x.Norm();
  if (!(r > 0.0)) return Invalid("front query at the solar centre");
  // The core parametrizes its front by heliocentric rays: the front point
  // reported for positionM is the one on the ray from the Sun through it,
  // whatever its distance from the front.  A ray outside the finite cap,
  // or a front record outside the support or below the core's physical
  // inner boundary, is "no front here": OK with supported = false.
  const auto ray = provider_->EvaluateRay(ToCore(x / r), epochS);
  if (!ray.ok()) {
    if (ray.status.code == SEP::Core::StatusCode::OutOfDomain)
      return Core::Status::OK();  // outside the finite cap: no front here
    return Invalid(ray.status.message);
  }
  const auto record = provider_->EvaluateFrontPoint(ray.value.positionM, epochS,
                                                    ray.value.stableId);
  if (!record.ok()) {
    if (record.status.code == SEP::Core::StatusCode::OutOfDomain)
      return Core::Status::OK();
    return Invalid(record.status.message);
  }
  const SF::ShockRecord& value = record.value;
  if (value.status == SF::FrontStatus::OutsideFrontSupport ||
      value.status == SF::FrontStatus::BelowPhysicalInnerBoundary)
    return Core::Status::OK();
  const Core::Vec3 pointHci = FromCore(value.geometry.positionM);
  const Core::Vec3 normalHci = FromCore(value.geometry.outwardNormal);
  state->supported = true;
  // Accepted only when the core solved a fast-mode shock with a valid
  // downstream state, forward inflow and compression; every other core
  // status (sub-fast front, unresolved weak shock, wrong branch, invalid
  // jump, unavailable ambient, non-forward inflow) is supported but not
  // accepted, so neither the source nor the absorbing boundary acts there.
  state->accepted = value.status == SF::FrontStatus::SolvedFastShock &&
                    value.downstreamValid && value.inflowMPerS > 0.0 &&
                    value.jump.compressionRatio > 1.0;
  state->positionM = frame_.FromInertialPosition(pointHci, epochS) + originM_;
  state->outwardNormal = frame_.FromInertialVector(normalHci, epochS);
  // The frame moves with velocity Omega x x at the front point; the front's
  // normal speed relative to the frame loses that component.  The inflow
  // w1 = Vn - U1.n is frame invariant (both terms lose the same amount).
  // For the heliocentric sphere n = e_r is perpendicular to Omega x x, so
  // the correction vanishes; it is nonzero on the flanks of the SSE cap.
  state->normalSpeedMPerS = value.geometry.normalSpeedMPerS -
      frame_.AngularVelocity().Cross(pointHci).Dot(normalHci);
  state->inflowMPerS = value.inflowMPerS;
  state->compressionRatio =
      state->accepted ? value.jump.compressionRatio : 1.0;
  // |cos theta_Bn| folds the angle into [0, pi/2]; min(1, .) guards acos
  // against a rounding excess.
  state->thetaBnRad = std::acos(std::min(1.0,
      std::fabs(value.signedMagneticNormalCosine)));
  state->upstreamProtonDensityM3 = value.upstream.plasma.protonNumberDensityM3;
  state->upstreamMagneticFieldT =
      frame_.FromInertialVector(FromCore(value.upstream.magneticFieldT), epochS);
  state->upstreamVelocityMPerS = frame_.FromInertialVelocity(
      FromCore(value.upstream.velocityMPerS), pointHci, epochS);
  // generation is always 0 here: the record comes from direct core queries
  // at epochS (EvaluateRay / EvaluateFrontPoint), not from a prepared Epoch.
  state->generation = 0;
  state->classification = value.reason.empty()
      ? std::string(state->accepted ? "solved-fast-shock" : "not-accepted")
      : value.reason;
  return Core::Status::OK();
}

Core::Status ReducedFrontShockSurface::ReleaseLevel(const Core::Vec3& positionM,
                                                    double epochS,
                                                    double fraction,
                                                    double* level) const {
  if (!level || !Finite(positionM) || !std::isfinite(fraction) || fraction < 0)
    return Invalid("invalid release-surface query");
  Core::Vec3 centre;
  double radius = 0.0;
  const Core::Status sphere = GeneratingSphere(epochS, &centre, &radius);
  if (!sphere.ok()) return sphere;
  const Core::Vec3 x = frame_.ToInertialPosition(positionM - originM_, epochS);
  // The 1-D srcSEP release surface: the generating sphere enlarged by
  // L_ref = f R about its own centre (sep_reduced_source.cpp).
  // R is the generating-sphere radius a: R_sh for the heliocentric sphere,
  // |c| sin(lambda) for the SSE cap (the 1-D ReleaseDistanceM also scales
  // the generating-sphere radius).  The level is the signed distance to the
  // enlarged sphere [m].
  *level = (x - centre).Norm() - (1.0 + fraction) * radius;
  return Core::Status::OK();
}

// Estimate of the largest front normal speed over [t0, t1] (clipped to the
// support): the apex speed sampled at 17 equally spaced epochs.  In HCI the
// apex has the largest normal speed of the front: the generating sphere's
// centre moves at V/(1 + s) along e_cme and its radius grows at V s/(1 + s)
// (V the apex speed, s = sin(lambda)), so a point with normal n has normal
// speed V (n.e_cme + s)/(1 + s) <= V.  The value is a sampled maximum, not a
// bound between samples; epochs where the trajectory has no state are
// skipped, and 0 is returned when none has one.
double ReducedFrontShockSurface::MaximumNormalSpeedMPerS(double t0S,
                                                         double t1S) const {
  const double begin = std::max(std::min(t0S, t1S), SupportBeginS());
  const double end = std::min(std::max(t0S, t1S), SupportEndS());
  double speed = 0.0, radius = 0.0;
  const int samples = 16;
  for (int i = 0; i <= samples; ++i) {
    const double t = begin + (end - begin) * i / samples;
    const auto state = provider_->Trajectory(t);
    if (!state.ok()) continue;
    speed = std::max(speed, std::fabs(state.value.apexSpeedMPerS));
    radius = std::max(radius, state.value.apexRadiusM);
  }
  // A non-spherical front also turns against a rotating frame.  The frame
  // velocity at a front point is at most |Omega| r <= |Omega| R_apex, which
  // is added as an upper bound.  For the heliocentric sphere Omega x x is
  // tangent to the front and does not change the normal speed.
  if (!Spherical()) speed += std::fabs(frame_.RotationRateRadPerS()) * radius;
  return speed;
}

}  // namespace Adapters
}  // namespace SEP3D
