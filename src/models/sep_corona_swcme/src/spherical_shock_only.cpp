#include "sep_corona_swcme/spherical_shock_only.h"

#include <algorithm>
#include <cmath>
#include <functional>
#include <iomanip>
#include <limits>
#include <locale>
#include <sstream>

#include "sep_coronal_cme/configuration_parser.h"
#include "sep_coronal_cme/constants.h"
#include "swcme_spherical_geometry.hpp"

namespace SEP {
namespace CoronaSwcme {
namespace SphericalShockOnly {
namespace {
using Core::StatusCode;
using CoronalCME::Cross;
using CoronalCME::Dot;
using CoronalCME::Norm;
using CoronalCME::Vec3;
bool Finite(Vec3 x) {
  return std::isfinite(x.x) && std::isfinite(x.y) && std::isfinite(x.z);
}
Vec3 Multiply(const std::array<double, 9>& a, Vec3 x) {
  return {a[0] * x.x + a[1] * x.y + a[2] * x.z, a[3] * x.x + a[4] * x.y + a[5] * x.z,
          a[6] * x.x + a[7] * x.y + a[8] * x.z};
}
Core::Result<AmbientState> OpenAt(const AmbientModel& model, Vec3 x, double t) {
  auto result = model.EvaluateWithDerivatives(x, t, 1);
  if (!result.ok())
    return result;
  if (result.value.primitive.region == AmbientRegion::PfssClosed)
    return Core::Result<AmbientState>::Failure(
        StatusCode::UnsupportedCapability,
        "spherical SHOCK_ONLY material candidate is restricted to actual open field "
        "lines");
  if (!ValidAmbientRegionSector(result.value.primitive.region,
                                result.value.primitive.magneticSector))
    return Core::Result<AmbientState>::Failure(
        StatusCode::DataIntegrityFailure,
        "spherical SHOCK_ONLY open sample has no actual signed magnetic sector");
  return result;
}
bool Covers(const Definition& d, double t) {
  return std::isfinite(t) && t >= d.ambient.support.startS &&
         t <= d.ambient.support.endS;
}
PhysicalSupport SupportAtRadius(const Definition& d, double radiusM) {
  return radiusM < d.physicalInnerRadiusM   ? PhysicalSupport::InnerGuard
         : radiusM > d.physicalOuterRadiusM ? PhysicalSupport::OuterGuard
                                            : PhysicalSupport::ParticleDomain;
}
std::string AmbientIdentity(const AmbientDefinition& d) {
  // Seal every typed physical/coverage input; a caller cannot reuse an old
  // label after altering a harmonic, species temperature or field convention.
  std::ostringstream out;
  out.imbue(std::locale::classic());
  out << std::hexfloat;
  for (const auto& s : {d.physicsFingerprint, d.coordinateFrame, d.ambient.harmonics})
    out << s.size() << ':' << s;
  const auto& c = d.composition;
  const auto& a = d.ambient;
  const auto& s = d.support;
  for (double x : {c.gammaAdiabatic,
                   c.electronTemperatureK,
                   c.protonTemperatureK,
                   c.alphaTemperatureK,
                   c.alphaToProtonNumberRatio,
                   a.sourceSurfaceRadiusM,
                   a.referenceRadiusM,
                   a.electronDensityAtReferenceM3,
                   a.closedBaseElectronDensityM3,
                   a.rotationRateRadPerS,
                   a.minimumMagneticFieldT,
                   a.traceStepM,
                   a.gradientRelativeStep,
                   s.startS,
                   s.endS,
                   s.solarRadiusM,
                   s.firstValidPlasmaRadiusM,
                   s.coverageRadiusM,
                   s.rootToleranceS,
                   s.maximumApexAccelerationMPerS2})
    out << x << ';';
  out << c.includeElectronMass << ';' << a.sphericalOpenMassFlux << ';'
      << a.traceMaximumSteps << ';' << a.windTablePoints;
  return CoronalCME::ComputeContentChecksum(out.str());
}
std::string EventIdentity(const Definition& d,const OpenTubeIsoPolyFamily* resolvedWind) {
  std::ostringstream out;
  out.imbue(std::locale::classic());
  out << std::hexfloat;
  for (const auto& s : {d.inputIdentity, d.ambient.physicsFingerprint})
    out << s.size() << ':' << s;
  out << d.physicalInnerRadiusM << ';' << d.physicalOuterRadiusM << ';'
      << static_cast<unsigned>(d.front.mode) << ';'
      << static_cast<unsigned>(d.front.extrapolation);
  for (double x : {d.front.r0_m, d.front.V0_m_s, d.front.Vsw_m_s, d.front.Gamma_m_inv})
    out << ';' << x;
  for (const auto* knots : {&d.front.data_time_s, &d.front.data_radius_m}) {
    out << ';' << knots->size();
    for (double x : *knots)
      out << ';' << x;
  }
  if(resolvedWind)out << ";resolved-wind:" << resolvedWind->Identity();
  return CoronalCME::ComputeContentChecksum(out.str());
}
Core::Status Controls(const Numerics& c) {
  for (double x : {c.positionToleranceM, c.deformationTolerance, c.maximumTimeStepS,
                   c.maximumArcStepM})
    if (!std::isfinite(x) || x <= 0)
      return Core::Status::Failure(
          StatusCode::InvalidConfiguration,
          "material integration requires explicit positive SI accuracy/work controls");
  if (c.maximumAttempts == 0)
    return Core::Status::Failure(
        StatusCode::InvalidConfiguration,
        "material integration requires an explicit nonzero work bound");
  return Core::Status::Success();
}

struct LocalFastShock {
  double mach=0;
  std::optional<CoronalCME::MhdShockSolution> jump;
};
Core::Result<LocalFastShock> FastShock(const AmbientPrimitive& p,Vec3 normal,
    double speed,double gamma) {
  using Return=Core::Result<LocalFastShock>;
  const CoronalCME::MhdPrimitiveState primitive{p.plasma.massDensityKgM3,
      p.plasma.pressurePa,p.velocityMPerS,p.magneticFieldT};
  const auto characteristic=CoronalCME::EvaluateMhdCharacteristics(primitive,normal,gamma);
  if(!characteristic.ok())return Return::Failure(characteristic.status.code,
      characteristic.status.message);
  LocalFastShock result;
  // Normal inflow is measured in the shock frame. The characteristic uses
  // the declared adiabatic derivative for waves, independently of the steady
  // wind's density exponent or any implied heating rate.
  result.mach=(speed-Dot(p.velocityMPerS,normal))/characteristic.value.fastSpeedMPerS;
  if(result.mach>1) {
    const auto jump=CoronalCME::SolveObliqueFastShock(primitive,normal,speed,gamma);
    if(!jump.ok())return Return::Failure(jump.status.code,jump.status.message);
    result.jump=jump.value;
  }
  return Return::Success(std::move(result));
}
}  // namespace

Core::Result<std::shared_ptr<const Event>> Event::Create(const Definition& input) {
  using Return = Core::Result<std::shared_ptr<const Event>>;
  if (input.inputIdentity.empty() || input.ambient.physicsFingerprint.empty() ||
      input.ambient.coordinateFrame != "HCI-like-inertial" ||
      !std::isfinite(input.physicalInnerRadiusM) ||
      !std::isfinite(input.physicalOuterRadiusM) ||
      input.physicalInnerRadiusM < input.ambient.support.solarRadiusM ||
      input.physicalOuterRadiusM <= input.physicalInnerRadiusM ||
      input.ambient.support.firstValidPlasmaRadiusM > input.physicalInnerRadiusM ||
      input.ambient.support.coverageRadiusM < input.physicalOuterRadiusM ||
      input.ambient.support.startS < 0 || input.ambient.ambient.sphericalOpenMassFlux ||
      input.ambient.ambient.isothermalFluxTubeOpenWind)
    return Return::Failure(
        StatusCode::InvalidConfiguration,
        "spherical open-field candidate requires explicit identities, HCI frame, "
        "physical/guard coverage and magnetic-flux density");
  if (swcme::kinematics::validate(input.front) != swcme::kinematics::Status::Ok ||
      input.front.extrapolation != swcme::kinematics::ExtrapolationPolicy::OutsideTime)
    return Return::Failure(StatusCode::InvalidConfiguration,
                           "spherical candidate requires supported canonical front "
                           "kinematics without extrapolation");
  // The canonical front laws are nondecreasing. Both endpoints therefore
  // certify front radial coverage; local trajectory/inserted-node guard
  // coverage is independently checked at every requested RK state.
  for (double t : {input.ambient.support.startS, input.ambient.support.endS}) {
    const auto front = swcme::kinematics::evaluate(input.front, t);
    if (front.status != swcme::kinematics::Status::Ok ||
        front.radius_m < input.ambient.support.firstValidPlasmaRadiusM ||
        front.radius_m > input.ambient.support.coverageRadiusM)
      return Return::Failure(StatusCode::OutOfDomain,
                             "canonical spherical front does not cover the complete "
                             "declared event interval");
  }
  auto event = std::shared_ptr<Event>(new Event);
  event->definition_ = input;
  event->definition_.ambient.physicsFingerprint = AmbientIdentity(input.ambient);
  const auto ambient = AmbientModel::Create(event->definition_.ambient);
  if (!ambient.ok())
    return Return::Failure(ambient.status.code, ambient.status.message);
  event->ambient_ = ambient.value;
  // One immutable wind supplies all resolved spatial samples, geometry and
  // front queries. Failure to prepare never falls back to the old kinematic
  // wind, which would change upstream speed, pressure and shock eligibility.
  if(input.resolvedWind) {
    const auto wind=OpenTubeIsoPolyFamily::Prepare(event->ambient_,*input.resolvedWind);
    if(!wind.ok())return Return::Failure(wind.status.code,wind.status.message);
    event->resolvedWind_=wind.value;
  }
  event->identity_ = EventIdentity(event->definition_,event->resolvedWind_.get());
  return Return::Success(std::move(event));
}

Core::Result<AmbientState> Event::AmbientAt(Vec3 x,double t,std::vector<double>* speciesHeat) const {
  if(!resolvedWind_)return OpenAt(*ambient_,x,t);
  auto sample=resolvedWind_->EvaluateWithDerivatives(x,t,1);
  using Return=Core::Result<AmbientState>;
  if(!sample.ok())return Return::Failure(sample.status.code,sample.status.message);
  if(speciesHeat)*speciesHeat=std::move(sample.value.speciesNetHeatWPerM3);
  return Return::Success(std::move(sample.value.ambient));
}

Core::Result<Sample> Event::Evaluate(Vec3 x, double t) const {
  using Return = Core::Result<Sample>;
  if (!Covers(definition_, t))
    return Return::Failure(StatusCode::OutOfDomain,
                           "spherical candidate time is outside its actual history");
  Sample result;
  const auto state = AmbientAt(x,t,&result.speciesNetHeatWPerM3);
  if (!state.ok())
    return Return::Failure(state.status.code, state.status.message);
  const auto front = swcme::kinematics::evaluate(definition_.front, t);
  if (front.status != swcme::kinematics::Status::Ok)
    return Return::Failure(StatusCode::OutOfDomain,
                           "spherical front has no history at the requested epoch");
  result.ambient = state.value;
  const double r = Norm(x);
  result.support = SupportAtRadius(definition_, r);
  // The spherical front is positioned in HCI at the requested epoch.
  // A smooth ambient sample behind it remains a reference state. Only the
  // surface query exports the one-sided RH states, so callers cannot use
  // this branch to infer a post-shock volume or a sheath evolution law.
  result.side = r > front.radius_m    ? FrontSide::Upstream
                : r == front.radius_m ? FrontSide::OnFront
                                      : FrontSide::BehindFrontReferenceOnly;
  const auto& a = state.value;
  const auto& p = a.primitive;
  const Vec3 omega{0, 0, definition_.ambient.ambient.rotationRateRadPerS};
  const Vec3 corotation = Cross(omega, x), relative = p.velocityMPerS - corotation;
  const double divergence = a.gradientU[0] + a.gradientU[4] + a.gradientU[8];
  // The canonical ambient is stationary in its rotating magnetic frame.
  // Eulerian derivatives are scalar_t=-grad(scalar).corotation and
  // vector_t=Omega x vector-grad(vector)*corotation. These exact frame
  // identities distinguish true material evolution from a frozen snapshot.
  const double rhoRate = Dot(a.gradientDensityKgM4, relative);
  const double pressureRate = Dot(a.gradientPressurePaPerM, relative);
  const Vec3 acceleration =
      Cross(omega, p.velocityMPerS) + Multiply(a.gradientU, relative);
  const Vec3 magneticRate =
      Cross(omega, p.magneticFieldT) + Multiply(a.gradientB, relative);
  result.continuityResidualKgM3S = rhoRate + p.plasma.massDensityKgM3 * divergence;
  result.inductionResidualTPerS = magneticRate -
                                  Multiply(a.gradientU, p.magneticFieldT) +
                                  divergence * p.magneticFieldT;
  const Vec3 curl{a.gradientB[7] - a.gradientB[5], a.gradientB[2] - a.gradientB[6],
                  a.gradientB[3] - a.gradientB[1]};
  const Vec3 lorentz =
      Cross(curl, p.magneticFieldT) / CoronalCME::Constants::kVacuumPermeabilityHPerM;
  const Vec3 gravity =
      -CoronalCME::Constants::kSolarGravitationalParameterM3PerS2 / (r * r * r) * x;
  // This is the SI momentum residual a sustaining model would have to
  // provide, with pressure and Lorentz forces removed from rho*(a-g).
  // Project along the signed B and dot with U for the work diagnostic.
  // Neither diagnostic applies a force or changes the trajectory.
  result.requiredExternalForceNPerM3 =
      p.plasma.massDensityKgM3 * (acceleration - gravity) + a.gradientPressurePaPerM -
      lorentz;
  result.fieldAlignedForcePerMassMPerS2 =
      Dot(result.requiredExternalForceNPerM3, p.magneticFieldT) /
      (p.plasma.massDensityKgM3 * Norm(p.magneticFieldT));
  result.requiredMechanicalPowerWPerM3 =
      Dot(p.velocityMPerS, result.requiredExternalForceNPerM3);
  const double gamma = definition_.ambient.composition.gammaAdiabatic;
  // The full pressure equation is valid for the stated shared adiabatic
  // index, whether T varies or not. It measures the heat required by the
  // prescribed thermal law, including compressional work and advection.
  result.requiredThermalHeatWPerM3 =
      (pressureRate + gamma * p.plasma.pressurePa * divergence) / (gamma - 1);
  if(!resolvedWind_)result.requiredIsothermalHeatWPerM3=result.requiredThermalHeatWPerM3;
  if (!Finite(result.inductionResidualTPerS) ||
      !Finite(result.requiredExternalForceNPerM3) ||
      !std::isfinite(result.continuityResidualKgM3S) ||
      !std::isfinite(result.fieldAlignedForcePerMassMPerS2) ||
      !std::isfinite(result.requiredMechanicalPowerWPerM3) ||
      !std::isfinite(result.requiredThermalHeatWPerM3))
    return Return::Failure(StatusCode::NumericalFailure,
                           "material balance diagnostic is nonfinite");
  return Return::Success(std::move(result));
}

Core::Result<Motion> Event::Integrate(Vec3 initial, double begin, double duration,
                                      const Numerics& controls) const {
  using Return = Core::Result<Motion>;
  using State = std::array<double, 12>;
  const auto checked = Controls(controls);
  if (!checked.ok())
    return Return::Failure(checked.code, checked.message);
  if (!std::isfinite(duration) || !Covers(definition_, begin) ||
      !Covers(definition_, begin + duration))
    return Return::Failure(StatusCode::OutOfDomain,
                           "material trajectory interval is outside declared coverage");
  const auto initialSample = AmbientAt(initial,begin);
  if (!initialSample.ok())
    return Return::Failure(initialSample.status.code, initialSample.status.message);
  // The interface refinement advances by metres or less at heliocentric
  // radii of order 10^9 m. Evolve the spatial profile's displacements and
  // F-I so that the local error estimate sees the RK increments without
  // losing them to roundoff in the large coordinate/identity offsets.
  const Vec3 origin = resolvedWind_ ? initial : Vec3{};
  const double identityOffset = resolvedWind_ ? 1 : 0;
  const auto position=[&](const State& x) {return origin+Vec3{x[0],x[1],x[2]};};
  const auto deformation=[&](const State& x,unsigned i,unsigned j) {
    return x[3+3*i+j]+(i==j?identityOffset:0);
  };
  State start{};
  start[0] = initial.x-origin.x;
  start[1] = initial.y-origin.y;
  start[2] = initial.z-origin.z;
  for (unsigned i = 0; i < 3; ++i)
    start[3 + 3 * i + i] = 1-identityOffset;
  std::size_t attempts = 0;
  struct StageRates {
    double maximumSpeed = 0;
    double maximumDeformationRate = 0;
    bool coronal = false;
    bool exterior = false;
  };
  const auto derivative = [&](const State& state, double relativeTime, State& out,
                              StageRates& rates) {
    const auto sample = AmbientAt(position(state),begin+relativeTime);
    if (!sample.ok())
      return sample.status;
    const auto& a = sample.value;
    const Vec3 v = a.primitive.velocityMPerS;
    out[0] = v.x;
    out[1] = v.y;
    out[2] = v.z;
    rates.maximumSpeed = std::max(rates.maximumSpeed, Norm(v));
    rates.coronal |= a.primitive.region == AmbientRegion::PfssOpen;
    rates.exterior |= a.primitive.region == AmbientRegion::ParkerExterior;
    // Variational equation for THE SAME flow: Fdot=grad(U)*F. A density
    // ratio alone is never relabeled as a complete three-dimensional map.
    for (unsigned i = 0; i < 3; ++i)
      for (unsigned j = 0; j < 3; ++j) {
        out[3 + 3 * i + j] = 0;
        for (unsigned k = 0; k < 3; ++k)
          out[3 + 3 * i + j] += a.gradientU[3 * i + k] * deformation(state,k,j);
        rates.maximumDeformationRate =
            std::max(rates.maximumDeformationRate,std::abs(out[3 + 3 * i + j]));
      }
    for (double value : out)
      if (!std::isfinite(value))
        return Core::Status::Failure(StatusCode::NumericalFailure,
                                     "material variational derivative is nonfinite");
    return Core::Status::Success();
  };
  const auto rk = [&](const State& original, double t, double h, State& out,
                      StageRates& rates) {
    State k1{}, k2{}, k3{}, k4{}, temporary{};
    auto status = derivative(original, t, k1, rates);
    if (!status.ok())
      return status;
    for (unsigned i = 0; i < 12; ++i)
      temporary[i] = original[i] + h * k1[i] / 2;
    status = derivative(temporary, t + h / 2, k2, rates);
    if (!status.ok())
      return status;
    for (unsigned i = 0; i < 12; ++i)
      temporary[i] = original[i] + h * k2[i] / 2;
    status = derivative(temporary, t + h / 2, k3, rates);
    if (!status.ok())
      return status;
    for (unsigned i = 0; i < 12; ++i)
      temporary[i] = original[i] + h * k3[i];
    status = derivative(temporary, t + h, k4, rates);
    if (!status.ok())
      return status;
    for (unsigned i = 0; i < 12; ++i)
      out[i] = original[i] + h * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]) / 6;
    return Core::Status::Success();
  };
  // Reserve part of the caller's total error budget for the one-sided
  // gradient at the physical PFSS/Parker interface. U is continuous but
  // its azimuthal derivative changes there. Across the interface RK4's
  // smooth Richardson estimate is invalid: an unresolved jump in Fdot
  // gives O(h) error and an error budget proportional to h would never
  // let that crossing converge. Accumulate an absolute allowance for
  // interface panels and keep the usual shared budget for smooth panels.
  const double smoothFraction = resolvedWind_ ? 0.5 : 1;
  double crossingPositionBudget = (1-smoothFraction)*controls.positionToleranceM;
  double crossingDeformationBudget = (1-smoothFraction)*controls.deformationTolerance;
  const double sourceRadius=definition_.ambient.ambient.sourceSurfaceRadiusM;
  const auto radiusOf=[&](const State& x) {return Norm(position(x));};
  std::function<Core::Status(const State&, double, double, State&)> advance;
  advance = [&](const State& original, double low, double high, State& out) {
    if (++attempts > controls.maximumAttempts)
      return Core::Status::Failure(
          StatusCode::NumericalFailure,
          "material trajectory exhausted its declared work bound");
    const double h = high - low, middle = low + h / 2;
    if (middle == low || middle == high)
      return Core::Status::Failure(
          StatusCode::NumericalFailure,
          "material trajectory cannot represent a required numerical subdivision");
    State coarse{}, half{}, fine{};
    StageRates rates;
    bool split = std::abs(h) > controls.maximumTimeStepS;
    if (!split) {
      auto s = rk(original, low, h, coarse, rates);
      if (!s.ok())
        return s;
      s = rk(original, low, h / 2, half, rates);
      if (!s.ok())
        return s;
      s = rk(half, middle, h / 2, fine, rates);
      if (!s.ok())
        return s;
      double positionError = 0, deformationError = 0;
      for (unsigned i = 0; i < 3; ++i)
        positionError += std::pow(fine[i] - coarse[i], 2);
      positionError = std::sqrt(positionError) / 15;
      for (unsigned i = 3; i < 12; ++i)
        deformationError =
            std::max(deformationError, std::abs(fine[i] - coarse[i]) / 15);
      const bool handoff=resolvedWind_&&((rates.coronal&&rates.exterior)||
          ((radiusOf(original)<sourceRadius)!=(radiusOf(fine)<sourceRadius)));
      split = rates.maximumSpeed * std::abs(h) > controls.maximumArcStepM;
      if(handoff) {
        // Use a conservative first-order rate envelope from both quadratures
        // for the tiny unresolved interval. Twice the measured maximum
        // bounds the difference between two increments at that rate. It is
        // charged only when this panel is accepted. No change to host dt,
        // epoch, branch side or magnetic field is used to cross the interface.
        const double positionCost=2*std::abs(h)*rates.maximumSpeed;
        const double deformationCost=2*std::abs(h)*rates.maximumDeformationRate;
        split=split||positionCost>crossingPositionBudget||deformationCost>crossingDeformationBudget;
        if(!split) {
          crossingPositionBudget-=positionCost;
          crossingDeformationBudget-=deformationCost;
        }
      } else {
        const double share = smoothFraction*std::abs(h/duration);
        split=split||positionError>controls.positionToleranceM*share||
            deformationError>controls.deformationTolerance*share;
      }
    }
    if (!split) {
      out = fine;
      return Core::Status::Success();
    }
    State mid{};
    auto s = advance(original, low, middle, mid);
    if (!s.ok())
      return s;
    return advance(mid, middle, high, out);
  };
  State final = start;
  if (duration != 0) {
    const auto status = advance(start, 0, duration, final);
    if (!status.ok())
      return Return::Failure(status.code, status.message);
  }
  Motion result;
  result.positionM = position(final);
  const auto endpoint = AmbientAt(result.positionM,begin+duration);
  if (!endpoint.ok())
    return Return::Failure(endpoint.status.code, endpoint.status.message);
  result.velocityMPerS = endpoint.value.primitive.velocityMPerS;
  for (unsigned i = 0; i < 3; ++i)
    for (unsigned j = 0; j < 3; ++j)
      result.deformation[3*i+j] = deformation(final,i,j);
  const auto& f = result.deformation;
  // A positive volume ratio J=det(F) is required by the mass and Cauchy
  // magnetic updates. Compute its cancellation-prone products in extended
  // precision before checking the value exposed at the material boundary.
  const long double determinant =
      static_cast<long double>(f[0]) * (static_cast<long double>(f[4]) * f[8] -
                                        static_cast<long double>(f[5]) * f[7]) -
      static_cast<long double>(f[1]) * (static_cast<long double>(f[3]) * f[8] -
                                        static_cast<long double>(f[5]) * f[6]) +
      static_cast<long double>(f[2]) * (static_cast<long double>(f[3]) * f[7] -
                                        static_cast<long double>(f[4]) * f[6]);
  result.jacobian = static_cast<double>(determinant);
  if (!Finite(result.positionM) || !std::isfinite(result.jacobian) ||
      result.jacobian <= 0)
    return Return::Failure(StatusCode::NumericalFailure,
                           "material trajectory is nonfinite or orientation reversing");
  return Return::Success(std::move(result));
}
Core::Result<Motion> Event::Move(Vec3 x, double begin, double end,
                                 const Numerics& n) const {
  if (!Covers(definition_, end))
    return Core::Result<Motion>::Failure(
        StatusCode::OutOfDomain,
        "material trajectory endpoint is outside declared history");
  return Integrate(x, begin, end - begin, n);
}
Core::Result<Motion> Event::AdvanceFixed(Vec3 x, double begin, double end, double dt,
                                         const Numerics& n) const {
  if (!std::isfinite(dt) || dt <= 0 || end <= begin || end != begin + dt)
    return Core::Result<Motion>::Failure(
        StatusCode::InvalidConfiguration,
        "a complete unchanged physical AMPS time step is required");
  return Integrate(x, begin, dt, n);
}

Core::Result<Shock> Event::Surface(Vec3 ray, double t) const {
  using Return = Core::Result<Shock>;
  if (!Covers(definition_, t) || !Finite(ray) || !(Norm(ray) > 0))
    return Return::Failure(
        StatusCode::OutOfDomain,
        "spherical surface requires a covered epoch and finite nondegenerate ray");
  const auto front = swcme::kinematics::evaluate(definition_.front, t);
  if (front.status != swcme::kinematics::Status::Ok)
    return Return::Failure(StatusCode::OutOfDomain,
                           "spherical surface kinematics has no declared state");
  const Vec3 unit = ray / Norm(ray);
  swcme::geometry::SphereSurface geometry;
  const auto status = swcme::geometry::spherical_surface(
      front.radius_m, {{unit.x, unit.y, unit.z}}, front.speed_m_s, geometry);
  if (!status.ok())
    return Return::Failure(StatusCode::NumericalFailure, status.summary());
  Shock result;
  result.epochS = t;
  result.radiusM = geometry.radius_m;
  result.normalSpeedMPerS = geometry.normal_speed_m_s;
  result.outwardNormal = unit;
  result.positionM = {geometry.position_m[0], geometry.position_m[1],
                      geometry.position_m[2]};
  const auto sample = Evaluate(result.positionM, t);
  if (!sample.ok())
    return Return::Failure(sample.status.code, sample.status.message);
  result.ambient = sample.value;
  // Surface owns the exact canonical radius. Recomputing the norm of its
  // Cartesian projection can round by an ULP to either side, especially for
  // a non-axis ray. This is an authority tag for the already constructed
  // surface, not a widened position tolerance or a displaced plasma query.
  result.ambient.side = FrontSide::OnFront;
  // Physical-boundary membership has the same canonical-radius authority.
  // A projected point one ULP inside 1.05 Rs or outside 5 AU must not turn an
  // exactly boundary-owned surface into a guard-only surface. Arbitrary
  // Evaluate queries continue to classify their actual Cartesian radius.
  result.ambient.support = SupportAtRadius(definition_, geometry.radius_m);
  auto physical=FastShock(sample.value.ambient.primitive,unit,front.speed_m_s,
      definition_.ambient.composition.gammaAdiabatic);
  if(!physical.ok())return Return::Failure(physical.status.code,physical.status.message);
  result.fastMach=physical.value.mach;
  result.physicalFast=physical.value.jump.has_value();
  result.jump=std::move(physical.value.jump);
  return Return::Success(std::move(result));
}

Core::Result<TubeShock> Event::SurfaceOnTube(const OpenTubeIsoPolyWind& wind,
    double azimuth,double epoch) const {
  using Return=Core::Result<TubeShock>;
  if(resolvedWind_)return Return::Failure(StatusCode::UnsupportedCapability,
      "resolved spherical events use their owned spatial wind through Surface or Connect");
  if(!wind.UsesAmbient(ambient_))return Return::Failure(StatusCode::InvalidConfiguration,
      "resolved wind and spherical event do not own the same immutable ambient inputs");
  if(!Covers(definition_,epoch))return Return::Failure(StatusCode::OutOfDomain,
      "spherical tube surface epoch is outside declared coverage");
  const auto front=swcme::kinematics::evaluate(definition_.front,epoch);
  if(front.status!=swcme::kinematics::Status::Ok)return Return::Failure(StatusCode::OutOfDomain,
      "spherical tube surface has no canonical front at the requested epoch");
  auto upstream=wind.Evaluate(front.radius_m,azimuth,epoch);
  if(!upstream.ok())return Return::Failure(upstream.status.code,upstream.status.message);
  // Intersect the canonical sphere at its exact radius with the resolved
  // monotone tube. The normal comes from the actual wind-tube position in HCI.
  const Vec3 normal=upstream.value.positionM/Norm(upstream.value.positionM);
  auto physical=FastShock(upstream.value.primitive,normal,front.speed_m_s,
      definition_.ambient.composition.gammaAdiabatic);
  if(!physical.ok())return Return::Failure(physical.status.code,physical.status.message);
  return Return::Success({epoch,front.radius_m,front.speed_m_s,physical.value.mach,normal,
      std::move(upstream.value),SupportAtRadius(definition_,front.radius_m),
      std::move(physical.value.jump)});
}

Core::Result<std::vector<Connection>> Event::Connect(const std::vector<Vec3>& vertices,
                                                     double epoch,
                                                     double tolerance) const {
  using Return = Core::Result<std::vector<Connection>>;
  if (vertices.size() < 2 || !std::isfinite(tolerance) || tolerance <= 0 ||
      !Covers(definition_, epoch))
    return Return::Failure(StatusCode::InvalidConfiguration,
                           "spherical line connection requires actual same-epoch "
                           "vertices and explicit SI residual tolerance");
  const auto front = swcme::kinematics::evaluate(definition_.front, epoch);
  if (front.status != swcme::kinematics::Status::Ok)
    return Return::Failure(
        StatusCode::OutOfDomain,
        "line connection has no front history at the requested epoch");
  std::vector<Connection> result;
  double length = 0;
  for (std::size_t segment = 0; segment + 1 < vertices.size(); ++segment) {
    const Vec3 a = vertices[segment], delta = vertices[segment + 1] - a;
    if (!Finite(a) || !Finite(delta) || !(Norm(delta) > 0))
      return Return::Failure(
          StatusCode::InvalidConfiguration,
          "field-line polyline has a nonfinite or zero-length segment");
    // Long-double coefficients start from exact input doubles. A negative
    // discriminant is a genuine miss, never promoted to a tangent by the
    // caller's residual tolerance; tangent and double crossings stay distinct.
    const long double ax = a.x, ay = a.y, az = a.z, dx = delta.x, dy = delta.y,
                      dz = delta.z, r = front.radius_m;
    const long double A = dx * dx + dy * dy + dz * dz,
                      B = 2 * (ax * dx + ay * dy + az * dz),
                      C = ax * ax + ay * ay + az * az - r * r;
    const long double disc = B * B - 4 * A * C;
    if (disc >= 0) {
      std::vector<long double> roots;
      if (disc == 0)
        roots.push_back(-B / (2 * A));
      else {
        const long double q = -.5L * (B + std::copysign(std::sqrt(disc), B));
        roots.push_back(q / A);
        roots.push_back(C / q);
        std::sort(roots.begin(), roots.end());
      }
      for (long double f : roots) {
        if (f < 0 || f > 1)
          continue;
        const double fraction = static_cast<double>(f);
        // Interior shared endpoints belong to the following segment, exactly
        // as native floor(S) ownership; retain the final outer endpoint once.
        if (fraction == 1 && segment + 2 < vertices.size())
          continue;
        const Vec3 point = a + fraction * delta;
        Connection connection;
        connection.segment = segment;
        connection.fraction = fraction;
        connection.arcLengthM = length + fraction * Norm(delta);
        connection.tangent = disc == 0;
        connection.surfaceResidualM = std::abs(Norm(point) - front.radius_m);
        if (connection.surfaceResidualM > tolerance)
          return Return::Failure(
              StatusCode::NumericalFailure,
              "spherical line root failed its unchanged spatial residual tolerance");
        const auto shock = Surface(point, epoch);
        if (!shock.ok())
          return Return::Failure(shock.status.code, shock.status.message);
        connection.shock = shock.value;
        result.push_back(std::move(connection));
      }
    }
    length += Norm(delta);
  }
  return Return::Success(std::move(result));
}

}  // namespace SphericalShockOnly
}  // namespace CoronaSwcme
}  // namespace SEP
