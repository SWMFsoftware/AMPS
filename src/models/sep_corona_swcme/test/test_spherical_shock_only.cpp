#include <algorithm>
#include <cmath>
#include <functional>
#include <iomanip>
#include <iostream>
#include <map>
#include <stdexcept>

#include "sep_corona_swcme/open_tube_isopoly.h"
#include "sep_corona_swcme/spherical_shock_only.h"
#include "sep_coronal_cme/constants.h"
#include "spherical_limit_reference.h"
#include "spherical_shock_only_fixture.h"

namespace C = SEP::CoronalCME;
namespace A = SEP::CoronaSwcme;
namespace S = A::SphericalShockOnly;
namespace {
constexpr double Rs = C::Constants::kSolarRadiusM, Au = 149597870700.;
void Require(bool v, const std::string& why) {
  if (!v)
    throw std::runtime_error(why);
}
void Near(double a, double b, double rel, const std::string& why) {
  Require(std::isfinite(a) && std::isfinite(b) &&
              std::abs(a - b) <= rel * std::max(std::abs(a), std::abs(b)),
          why);
}
std::shared_ptr<const S::Event> Event(const S::Definition& d) {
  auto r = S::Event::Create(d);
  Require(r.ok(), r.status.message);
  return r.value;
}

C::Vec3 Apply(const std::array<double, 9>& f, C::Vec3 b) {
  return {f[0] * b.x + f[1] * b.y + f[2] * b.z, f[3] * b.x + f[4] * b.y + f[5] * b.z,
          f[6] * b.x + f[7] * b.y + f[8] * b.z};
}
S::Definition Fixture(int polarity = 1) {
  return S::Verification::Configuration(polarity);
}
S::Numerics Controls() {
  return S::Verification::IntegrationControls();
}
namespace R = spherical_limit_reference;
R::Profile MatchedReference(const S::Definition& input, int polarity) {
  Require(input.ambient.ambient.harmonics ==
              (polarity == 1 ? "1,0,0.0001,0" : "1,0,-0.0001,0"),
          "matched1D oracle requires its explicitly declared pure dipole input");
  const auto& a = input.ambient.ambient;
  const auto& c = input.ambient.composition;
  R::Parameters p;
  p.solarRadius = input.ambient.support.solarRadiusM;
  p.sourceRadius = a.sourceSurfaceRadiusM;
  p.referenceRadius = a.referenceRadiusM;
  p.maximumRadius = input.ambient.support.coverageRadiusM;
  p.electronDensityAtReference = a.electronDensityAtReferenceM3;
  p.electronTemperature = c.electronTemperatureK;
  p.protonTemperature = c.protonTemperatureK;
  p.alphaTemperature = c.alphaTemperatureK;
  p.alphaRatio = c.alphaToProtonNumberRatio;
  p.rotation = a.rotationRateRadPerS;
  p.gamma = c.gammaAdiabatic;
  p.dipoleCoefficient = polarity * 1e-4L;
  p.includeElectronMass = c.includeElectronMass;
  return R::Profile(p);
}
void ScaledMatch(R::Real actual, R::Real expected, R::Real physicalScale,
                 R::Real tolerance, const std::string& name) {
  Require(std::isfinite(actual) && std::isfinite(expected) && physicalScale > 0 &&
              std::isfinite(physicalScale) &&
              std::abs(actual - expected) <= tolerance * physicalScale,
          "independent matched1D/3D " + name);
}
void MatchedSphericalPrimitives() {
  R::Real maximumPrimitive = 0, maximumGradient = 0;
  for (int polarity : {1, -1}) {
    const auto config = Fixture(polarity);
    const auto e = Event(config);
    const auto reference = MatchedReference(config, polarity);
    const double source = config.ambient.ambient.sourceSurfaceRadiusM;
    for (int hemisphere : {1, -1})
      for (double radius :
           {1.05 * Rs, 1.1 * Rs, 2 * Rs, source - 1e5, source, source + 1e5,
            static_cast<double>(reference.CriticalRadius()), .1 * Au, Au, 5 * Au}) {
        const auto r = reference.At(radius, hemisphere);
        const auto sample = e->Evaluate({0, 0, hemisphere * radius}, 100);
        Require(sample.ok(), sample.status.message);
        const auto& a = sample.value.ambient;
        const auto& state = a.primitive;
        const auto scalar = [&](R::Real actual, R::Real expected, R::Real scale,
                                const std::string& name) {
          ScaledMatch(actual, expected, scale, R::PrimitiveTolerance, name);
          maximumPrimitive =
              std::max(maximumPrimitive, std::abs(actual - expected) / scale);
        };
        scalar(state.plasma.massDensityKgM3, r.density, r.density, "mass density");
        scalar(state.plasma.pressurePa, r.pressure, r.pressure,
               "total species pressure");
        scalar(state.plasma.electronNumberDensityM3, r.electronDensity,
               r.electronDensity, "electron density");
        scalar(state.plasma.protonNumberDensityM3, r.protonDensity, r.protonDensity,
               "proton density");
        const double velocity[] = {state.velocityMPerS.x, state.velocityMPerS.y,
                                   state.velocityMPerS.z};
        const double magnetic[] = {state.magneticFieldT.x, state.magneticFieldT.y,
                                   state.magneticFieldT.z};
        const double densityGradient[] = {
            a.gradientDensityKgM4.x, a.gradientDensityKgM4.y, a.gradientDensityKgM4.z};
        const double pressureGradient[] = {a.gradientPressurePaPerM.x,
                                           a.gradientPressurePaPerM.y,
                                           a.gradientPressurePaPerM.z};
        for (unsigned component = 0; component < 3; ++component) {
          scalar(velocity[component], r.velocity[component], r.speed, "HCI velocity");
          scalar(magnetic[component], r.magnetic[component], std::abs(r.magneticZ),
                 "signed HCI B");
          ScaledMatch(densityGradient[component], r.densityGradient[component],
                      r.density / r.radius, R::GradientTolerance,
                      "Cartesian density gradient");
          ScaledMatch(pressureGradient[component], r.pressureGradient[component],
                      r.pressure / r.radius, R::GradientTolerance,
                      "Cartesian pressure gradient");
        }
        for (unsigned component = 0; component < 9; ++component) {
          const R::Real uScale = r.speed / r.radius,
                        bScale = std::abs(r.magneticZ) / r.radius;
          ScaledMatch(a.gradientU[component], r.velocityGradient[component], uScale,
                      R::GradientTolerance, "full Cartesian velocity tensor");
          ScaledMatch(a.gradientB[component], r.magneticGradient[component], bScale,
                      R::GradientTolerance, "full Cartesian magnetic tensor");
          maximumGradient = std::max(
              {maximumGradient,
               std::abs(a.gradientU[component] - r.velocityGradient[component]) /
                   uScale,
               std::abs(a.gradientB[component] - r.magneticGradient[component]) /
                   bScale});
        }
      }
  }
  std::cout << std::setprecision(17) << "matched1D3D_max_primitive_scaled="
            << static_cast<double>(maximumPrimitive)
            << " max_gradient_scaled=" << static_cast<double>(maximumGradient) << '\n';
}
void MatchedSphericalFlow() {
  R::Real maximumPosition = 0, maximumDeformation = 0, maximumJacobian = 0;
  for (int polarity : {1, -1}) {
    const auto config = Fixture(polarity);
    const auto e = Event(config);
    const auto reference = MatchedReference(config, polarity);
    const double source = config.ambient.ambient.sourceSurfaceRadiusM;
    for (int hemisphere : {1, -1})
      for (const auto& role : std::vector<std::pair<double, double>>{
               {1.1 * Rs, 1}, {source - 1000, 1}, {.1 * Au, 100}}) {
        const auto expected = reference.Advance(role.first, role.second);
        const auto actual = e->AdvanceFixed({0, 0, hemisphere * role.first}, 100,
                                            100 + role.second, role.second, Controls());
        Require(actual.ok(), actual.status.message);
        const auto& map = actual.value;
        const R::Real positionError =
            std::hypot(std::hypot(static_cast<R::Real>(map.positionM.x),
                                  static_cast<R::Real>(map.positionM.y)),
                       map.positionM.z - hemisphere * expected.radius);
        Require(positionError < R::PositionToleranceM,
                "independent radial travel-time/3D fixed-step position mismatch");
        maximumPosition = std::max(maximumPosition, positionError);
        for (unsigned i = 0; i < 9; ++i) {
          ScaledMatch(map.deformation[i], expected.deformation[i], 1,
                      R::DeformationTolerance,
                      "full travel-time/tube-area/rotation deformation");
          maximumDeformation =
              std::max(maximumDeformation,
                       std::abs(map.deformation[i] - expected.deformation[i]));
        }
        ScaledMatch(map.jacobian, expected.jacobian, expected.jacobian,
                    R::DeformationTolerance, "one-dimensional tube-map Jacobian");
        maximumJacobian =
            std::max(maximumJacobian, std::abs(map.jacobian / expected.jacobian - 1));
      }
  }
  std::cout << std::setprecision(17)
            << "matched1D3D_max_position_m=" << static_cast<double>(maximumPosition)
            << " max_F=" << static_cast<double>(maximumDeformation)
            << " max_J_relative=" << static_cast<double>(maximumJacobian) << '\n';
}
void MatchedParallelShock() {
  R::Real maximumError = 0;
  for (int polarity : {1, -1}) {
    auto config = Fixture(polarity);
    config.front.V0_m_s = 1.5e6;
    const auto e = Event(config);
    const auto reference = MatchedReference(config, polarity);
    for (int hemisphere : {1, -1}) {
      const auto expected =
          reference.ParallelJump(config.front.r0_m, config.front.V0_m_s, hemisphere);
      const auto actual = e->Surface({0, 0, static_cast<double>(hemisphere)}, 0);
      Require(actual.ok() && actual.value.physicalFast && actual.value.jump,
              "matched spherical parallel front lacks an admissible actual RH jump");
      const auto& jump = *actual.value.jump;
      const auto check = [&](R::Real a, R::Real b, const std::string& name) {
        ScaledMatch(a, b, std::abs(b), R::PrimitiveTolerance, name);
        maximumError = std::max(maximumError, std::abs((a - b) / b));
      };
      check(jump.compressionRatio, expected.compression, "parallel compression");
      check(jump.downstream.massDensityKgM3, expected.density,
            "parallel downstream density");
      check(jump.downstream.pressurePa, expected.pressure,
            "parallel downstream total pressure");
      check(jump.downstream.velocityMPerS.z, expected.velocityZ,
            "parallel downstream HCI velocity");
      check(jump.downstream.magneticFieldT.z, expected.magneticZ,
            "parallel normal signed magnetic field");
      check(actual.value.fastMach, expected.fastMach,
            "parallel normal-relative fast Mach");
      Require(jump.residuals.maximum <= 1e-9,
              "matched spherical shock relaxed canonical RH residual gate");
    }
  }
  std::cout << std::setprecision(17) << "matched1D3D_max_parallel_jump_relative="
            << static_cast<double>(maximumError) << '\n';
}
void S13SC01() {
  const auto d = Fixture();
  const auto e = Event(d);
  auto invalid = d;
  invalid.ambient.ambient.sphericalOpenMassFlux = true;
  Require(!S::Event::Create(invalid).ok(),
          "radial density substituted for nonradial flux density");
  invalid = d;
  invalid.ambient.support.coverageRadiusM = 4 * Au;
  Require(!S::Event::Create(invalid).ok(),
          "provider does not cover physical particle domain");
  invalid = d;
  invalid.ambient.coordinateFrame = "relabeled-frame";
  Require(!S::Event::Create(invalid).ok(), "foreign frame silently relabeled");
  invalid = d;
  invalid.front.extrapolation = swcme::kinematics::ExtrapolationPolicy::Ballistic;
  Require(!S::Event::Create(invalid).ok(), "implicit front extrapolation accepted");
  auto changed = d;
  changed.ambient.ambient.electronDensityAtReferenceM3 *= 2;
  Require(Event(changed)->Identity() != e->Identity(),
          "typed ambient mutation reused stale physical event identity");
  changed = d;
  changed.front.V0_m_s = 100;
  Require(Event(changed)->Identity() != e->Identity(),
          "typed front mutation reused stale event identity");
}
void S13SC02() {
  MatchedSphericalPrimitives();
  const auto e = Event(Fixture());
  for (double r : {1.05 * Rs, .1 * Au, Au, 5 * Au}) {
    const C::Vec3 x{0, 0, r};
    const auto s = e->Evaluate(x, 0);
    Require(s.ok(), s.status.message);
    const auto& p = s.value.ambient.primitive;
    const auto direct = e->Ambient()->Evaluate(x, 0);
    Require(direct.ok(), direct.status.message);
    Require(p.plasma.massDensityKgM3 == direct.value.plasma.massDensityKgM3 &&
                p.magneticFieldT.z == direct.value.magneticFieldT.z,
            "public candidate differs from canonical ambient");
    if (r >= 2.5 * Rs)
      Near(p.plasma.massDensityKgM3 * p.velocityMPerS.z * r * r,
           e->Ambient()->MassFluxPerSteradianKgPerS(), 2e-12,
           "independent radial mass flux");
    Near(p.plasma.electronNumberDensityM3, p.plasma.protonNumberDensityM3 * 1.08, 2e-13,
         "charge-neutral proton/alpha export");
    // Independent Parker implicit equation, with the ACTUAL composition sound
    // speed, rather than a production-generated wind table as truth.
    const double cs2 = p.plasma.pressurePa / p.plasma.massDensityKgM3;
    const double rc = C::Constants::kSolarGravitationalParameterM3PerS2 / (2 * cs2);
    const double y = p.velocityMPerS.z * p.velocityMPerS.z / cs2;
    Near(y - std::log(y), 4 * std::log(r / rc) + 4 * rc / r - 3, 2e-8,
         "independent transonic Parker relation");
  }
}
void S13SC03() {
  MatchedSphericalFlow();
  const auto e = Event(Fixture());
  for (const auto& role :
       std::vector<std::pair<C::Vec3, double>>{{.1 * Au * C::Unit({.3, .4, .8}), 100.},
                                               {1.1 * Rs * C::Unit({.2, .1, 1}), 1.}}) {
    const C::Vec3 x = role.first;
    const double duration = role.second;
    const auto initial = e->Evaluate(x, 100);
    Require(initial.ok(), initial.status.message);
    const auto moved = e->AdvanceFixed(x, 100, 100 + duration, duration, Controls());
    Require(moved.ok(), moved.status.message);
    const auto final = e->Evaluate(moved.value.positionM, 100 + duration);
    Require(final.ok(), final.status.message);
    const auto& p0 = initial.value.ambient.primitive;
    const auto& p1 = final.value.ambient.primitive;
    Near(p1.plasma.massDensityKgM3 * moved.value.jacobian, p0.plasma.massDensityKgM3,
         2e-7, "independent material mass/J invariant");
    const auto frozen =
        Apply(moved.value.deformation, p0.magneticFieldT) / moved.value.jacobian;
    Require(C::Norm(frozen - p1.magneticFieldT) <= 2e-6 * C::Norm(p1.magneticFieldT),
            "independent full-vector Cauchy induction");
    const auto inverse =
        e->Move(moved.value.positionM, 100 + duration, 100, Controls());
    Require(inverse.ok(), inverse.status.message);
    Require(C::Norm(inverse.value.positionM - x) < .2,
            "inverse material map did not return to original point");
    const auto identity =
        Apply(inverse.value.deformation, Apply(moved.value.deformation, {.3, .2, .4}));
    Require(C::Norm(identity - C::Vec3{.3, .2, .4}) < 2e-8,
            "inverse deformation is not actual inverse flow");
  }
}
void S13SC04() {
  const auto e = Event(Fixture());
  const C::Vec3 x = 1.1 * Rs * C::Unit({.2, .1, 1});
  const auto sample = e->Evaluate(x, 100);
  Require(sample.ok(), sample.status.message);
  const auto& a = sample.value.ambient;
  const auto& p = a.primitive;
  const double scale = p.plasma.massDensityKgM3 * C::Norm(p.velocityMPerS) / C::Norm(x);
  Require(std::abs(sample.value.continuityResidualKgM3S) < 3e-5 * scale,
          "nonradial magnetic-flux mass continuity residual");
  Require(C::Norm(sample.value.inductionResidualTPerS) <
              3e-5 * C::Norm(p.magneticFieldT) * C::Norm(p.velocityMPerS) / C::Norm(x),
          "nonradial rotating-PFSS induction residual");
  // Independently difference the actual material velocity and field-aligned
  // pressure. This diagnostic must expose the finite nozzle mismatch rather
  // than advertise force-free Parker balance on a nonradial tube.
  const double dt = .1, h = 1e5;
  const C::Vec3 b = C::Unit(p.magneticFieldT);
  const auto plus = e->Ambient()->Evaluate(x + dt * p.velocityMPerS, 100 + dt);
  const auto minus = e->Ambient()->Evaluate(x - dt * p.velocityMPerS, 100 - dt);
  const auto bp = e->Ambient()->Evaluate(x + h * b, 100),
             bm = e->Ambient()->Evaluate(x - h * b, 100);
  Require(plus.ok() && minus.ok() && bp.ok() && bm.ok(),
          "independent momentum difference has missing open support");
  const C::Vec3 acceleration =
      (plus.value.velocityMPerS - minus.value.velocityMPerS) / (2 * dt);
  const C::Vec3 gravity =
      -C::Constants::kSolarGravitationalParameterM3PerS2 / std::pow(C::Norm(x), 3) * x;
  const double parallel = C::Dot(acceleration - gravity, b) +
                          (bp.value.plasma.pressurePa - bm.value.plasma.pressurePa) /
                              (2 * h * p.plasma.massDensityKgM3);
  Require(std::abs(parallel - sample.value.fieldAlignedForcePerMassMPerS2) < 1e-3,
          "independent field-aligned momentum force mismatch");
  Require(std::abs(parallel) > 1e-3,
          "nonradial kinematic Parker candidate hid its sustaining-force requirement");
  const double div = a.gradientU[0] + a.gradientU[4] + a.gradientU[8];
  Near(sample.value.requiredIsothermalHeatWPerM3, p.plasma.pressurePa * div, 1e-4,
       "isothermal heat differs from p div(U) under continuity");
  std::cout << "nonradial_parallel_force_m_s2=" << parallel << '\n';

  // Reuse the manufactured event's dipole and composition as controlled
  // geometry, with explicit thermal inputs for this distinct closure test.
  // Its comparator is the analytical dipole flux function and an independent
  // difference of the solved primitives. No observation is assigned by the test.
  A::OpenTubeIsoPolyRequest request;
  request.sourceColatitudeRad=.45;
  request.massLoadingKgPerSWb=2e-6;
  request.referenceDensityKgM3=1e-20;
  request.electron={1.15,10*Rs,18*Rs};
  request.proton={1.25,12*Rs,22*Rs};
  request.alpha={1.20,12*Rs,25*Rs};
  request.angularTraceToleranceRad=1e-12;
  const double source=e->Configuration().ambient.ambient.sourceSurfaceRadiusM;
  const double omega=e->Configuration().ambient.ambient.rotationRateRadPerS;
  double maximumForceResidual=0;
  for(int polarity:{1,-1}) {
    auto setup=Fixture(polarity);
    if(polarity<0)setup.front.V0_m_s=1.5e6;
    const auto event=Event(setup);
    const auto solved=A::OpenTubeIsoPolyWind::Prepare(event->Ambient(),request);
    Require(solved.ok(),solved.status.message);
    const auto surface=event->SurfaceOnTube(solved.value,.35,100);
    Require(surface.ok(),surface.status.message);
    Require(!e->SurfaceOnTube(solved.value,.35,100).ok(),
        "spherical front accepted a tube from another ambient owner");
    if(polarity>0)Require(surface.value.fastMach<1&&!surface.value.jump,
        "resolved wind emitted a jump for a sub-fast front");
    else {
      Require(surface.value.jump&&surface.value.fastMach>1,
          "resolved wind did not feed the physical fast shock");
      const auto& local=surface.value.upstream.primitive;
      const auto& jump=*surface.value.jump;
      const double incoming=surface.value.normalSpeedMPerS-
          C::Dot(local.velocityMPerS,surface.value.outwardNormal);
      Near(surface.value.fastMach,incoming/jump.upstreamCharacteristics.fastSpeedMPerS,
          1e-13,"resolved-wind shock used a different upstream state");
      Require(jump.residuals.maximum<=1e-9&&jump.entropyLogIncrement>0,
          "resolved-wind shock relaxed the unchanged RH gate");
    }
    Require(!solved.value.Evaluate(.99*Rs,0,100).ok()&&
        !solved.value.Evaluate(10*Au+1,0,100).ok()&&
        !solved.value.Evaluate(Au,0,1201).ok(),"wind extrapolated a radius or epoch");
    for(double radius:{1.05*Rs,1.8*Rs,3*Rs,solved.value.CriticalRadiusM(),
        15*Rs,20*Rs,Au,5*Au}) {
      const auto at=solved.value.Evaluate(radius,.35,100);
      Require(at.ok(),at.status.message);
      const auto& w=at.value;
      const auto& v=w.primitive;
      const auto radial=w.positionM/radius;
      const C::Vec3 relative=v.velocityMPerS-C::Cross({0,0,omega},w.positionM);
      const double bNorm=C::Norm(v.magneticFieldT);
      Near(v.plasma.massDensityKgM3*C::Norm(relative)/bNorm,request.massLoadingKgPerSWb,
          1e-12,"resolved wind leaks mass per magnetic flux");
      Require(C::Norm(C::Cross(relative,v.magneticFieldT))<1e-13*C::Norm(relative)*bNorm&&
          v.magneticSector==polarity,"resolved wind violates signed rotating ideal flux alignment");
      if(radius<source) {
        const double cosTheta=radial.z;
        Near((1-cosTheta*cosTheta)*(radius*radius+2*source*source*source/radius),
            3*source*source*std::pow(std::sin(request.sourceColatitudeRad),2),1e-9,
            "trace does not follow the independent pure-dipole flux function");
      }
      const double dr=radius*1e-5;
      const auto lo=solved.value.Evaluate(radius-dr,.35,100);
      const auto hi=solved.value.Evaluate(radius+dr,.35,100);
      Require(lo.ok()&&hi.ok(),"independent momentum stencil is outside covered tube");
      const auto& left=lo.value.balance.wind;
      const auto& right=hi.value.balance.wind;
      const double inertia=w.balance.wind.speedMPerS*(right.speedMPerS-left.speedMPerS)/(2*dr);
      const double pressure=(right.pressurePa-left.pressurePa)/(2*dr*v.plasma.massDensityKgM3);
      const double potential=(right.effectivePotentialM2S2-left.effectivePotentialM2S2)/(2*dr);
      const double scale=std::abs(inertia)+std::abs(pressure)+std::abs(potential);
      maximumForceResidual=std::max(maximumForceResidual,
          std::abs(inertia+pressure+potential)/scale);
      const double gamma=setup.ambient.composition.gammaAdiabatic;
      const double radialSpeed=C::Dot(v.velocityMPerS,radial);
      const double densityRate=(right.densityKgM3-left.densityKgM3)/(2*dr);
      for(std::size_t j=0;j<w.balance.netHeatWPerM3.size();++j) {
        const double dP=(hi.value.balance.partialPressurePa[j]-
                         lo.value.balance.partialPressurePa[j])/(2*dr);
        const double heat=radialSpeed/(gamma-1)*(dP-
            gamma*w.balance.partialPressurePa[j]/v.plasma.massDensityKgM3*densityRate);
        Near(w.balance.netHeatWPerM3[j],heat,4e-5,
            "species heating uses the wrong material derivative on the traced tube");
      }
    }
    const auto inner=solved.value.Evaluate(source-1,.35,100);
    const auto outer=solved.value.Evaluate(source+1,.35,100);
    Require(inner.ok()&&outer.ok(),"wind handoff does not span the source surface");
    const auto normal=outer.value.positionM/(source+1);
    const auto& a=inner.value.primitive;
    const auto& b=outer.value.primitive;
    Near(C::Dot(a.magneticFieldT,normal),C::Dot(b.magneticFieldT,normal),1e-7,
        "wind handoff loses normal magnetic flux");
    Near(a.plasma.massDensityKgM3*C::Dot(a.velocityMPerS,normal),
        b.plasma.massDensityKgM3*C::Dot(b.velocityMPerS,normal),1e-7,
        "wind handoff loses normal mass flux");
    Require(C::Norm(C::Cross(normal,C::Cross(a.velocityMPerS,a.magneticFieldT)-
        C::Cross(b.velocityMPerS,b.magneticFieldT)))<1e-7*C::Norm(b.velocityMPerS)*
        C::Norm(b.magneticFieldT),"wind handoff has a tangential electric-field jump");
  }
  Require(maximumForceResidual<3e-5,"resolved tube requires an unaccounted field-aligned force");
  std::cout << "isopoly_differenced_momentum_max_scaled=" << maximumForceResidual << '\n';
}
void S13SC05() {
  const auto d = Fixture();
  const auto e = Event(d);
  const double r = d.ambient.ambient.sourceSurfaceRadiusM;
  const auto below = e->Evaluate({0, 0, r - 100}, 10),
             above = e->Evaluate({0, 0, r + 100}, 10);
  Require(below.ok() && above.ok(),
          "source-surface one-sided primitives/tensors missing");
  Near(below.value.ambient.primitive.plasma.massDensityKgM3,
       above.value.ambient.primitive.plasma.massDensityKgM3, 2e-6,
       "source-surface mass-flux density continuity");
  const auto moved = e->Move({0, 0, r - 1000}, 10, 11, Controls());
  Require(moved.ok(), moved.status.message);
  Require(moved.value.positionM.z > r,
          "material flow never crossed actual source surface");
  const auto final = e->Evaluate(moved.value.positionM, 11);
  Require(final.ok(), final.status.message);
  Require(final.value.ambient.primitive.region == A::AmbientRegion::ParkerExterior,
          "crossing retained wrong physical region");
  const auto initial = e->Evaluate({0, 0, r - 1000}, 10);
  Require(initial.ok(), initial.status.message);
  Near(final.value.ambient.primitive.plasma.massDensityKgM3 * moved.value.jacobian,
       initial.value.ambient.primitive.plasma.massDensityKgM3, 2e-7,
       "source-surface material mass continuity");
  const auto transported =
      Apply(moved.value.deformation, initial.value.ambient.primitive.magneticFieldT) /
      moved.value.jacobian;
  Require(C::Norm(transported - final.value.ambient.primitive.magneticFieldT) <
              2e-6 * C::Norm(transported),
          "source-surface full-vector Cauchy flux continuity");
}
void S13SC06() {
  MatchedParallelShock();
  auto d = Fixture();
  d.front.r0_m = 1.05 * Rs;
  d.front.V0_m_s = 500;
  const auto launch = Event(d)->Surface({.05, .25, 1}, 0);
  Require(launch.ok() && launch.value.radiusM == 1.05 * Rs &&
              !launch.value.physicalFast && !launch.value.jump,
          "physical1.05Rs launch was forced to become a shock");
  Require(launch.value.ambient.support == S::PhysicalSupport::ParticleDomain,
          "canonical1.05Rs surface was misclassified as inner guard by projection "
          "rounding");
  d = Fixture();
  auto e = Event(d);
  auto slow = e->Surface({.4, .2, 1}, 0);
  Require(slow.ok() && !slow.value.physicalFast && !slow.value.jump,
          "sub-fast geometry invented a jump");
  d.front.V0_m_s = 1.5e6;
  e = Event(d);
  const auto fast = e->Surface({.1, .2, 1}, 0);
  Require(fast.ok(), fast.status.message);
  Require(fast.value.ambient.side == S::FrontSide::OnFront,
          "canonical non-axis surface projection lost its OnFront authority");
  Require(fast.value.physicalFast && fast.value.jump && fast.value.fastMach > 1,
          "real fast normal-relative oblique shock absent");
  const auto& jump = *fast.value.jump;
  Require(jump.compressionRatio > 1 && jump.entropyLogIncrement > 0 &&
              jump.residuals.maximum <= 1e-9,
          "unchanged canonical oblique-RH gate");
  const auto& a = fast.value.ambient.ambient.primitive;
  const double normal =
      fast.value.normalSpeedMPerS - C::Dot(a.velocityMPerS, fast.value.outwardNormal);
  Near(fast.value.fastMach, normal / jump.upstreamCharacteristics.fastSpeedMPerS, 1e-13,
       "shock used observer speed instead of local normal inflow");
}
void S13SC07() {
  const auto e = Event(Fixture());
  const double r = 4 * Rs;
  const auto roots =
      e->Connect({{-2 * r, 0, .5 * r}, {0, 0, .5 * r}, {2 * r, 0, .5 * r}}, 0, 1e-5);
  Require(roots.ok(), roots.status.message);
  Require(roots.value.size() == 2, "multiple sphere/curved-line roots lost");
  const double expected = std::sqrt(.75) * r;
  Near(std::abs(roots.value[0].shock.positionM.x), expected, 2e-14,
       "independent first quadratic root");
  Near(roots.value[1].shock.positionM.x, expected, 2e-14,
       "independent second quadratic root");
  const auto tangent = e->Connect({{-r, 0, r}, {r, 0, r}}, 0, 1e-5);
  Require(tangent.ok() && tangent.value.size() == 1 && tangent.value[0].tangent,
          "exact tangency lost or doubled");
  const auto miss = e->Connect({{-r, 0, r + 1}, {r, 0, r + 1}}, 0, 1e-5);
  Require(miss.ok() && miss.value.empty(),
          "outside line miss was promoted to tangent by residual tolerance");
  const auto shared = e->Connect({{0, 0, .5 * r}, {0, 0, r}, {0, 0, 2 * r}}, 0, 1e-5);
  Require(shared.ok() && shared.value.size() == 1 && shared.value[0].segment == 1 &&
              shared.value[0].fraction == 0,
          "shared endpoint native ownership duplicated/lost");
}
void S13SC08() {
  auto d = Fixture();
  const auto e = Event(d);
  const C::Vec3 x{0, 0, 5 * Au};
  const auto guard = e->Move(x, 0, 100, Controls());
  Require(guard.ok(), guard.status.message);
  const auto state = e->Evaluate(guard.value.positionM, 100);
  Require(state.ok() && state.value.support == S::PhysicalSupport::OuterGuard,
          "physical escape sphere was conflated with provider guard");
  Require(!e->AdvanceFixed(x, 0, 50, 100, Controls()).ok(),
          "shorter physical dt accepted");
  Require(!e->Evaluate({0, 0, 10 * Au + 1}, 0).ok(),
          "provider extrapolated past declared guard");
  Require(!e->Evaluate({0, 0, 1.1 * Rs}, -1).ok(),
          "provider extrapolated before actual event");
  Require(!e->Evaluate({1.1 * Rs, 0, 0}, 0).ok(),
          "closed coronal loop admitted as open transport");
  Require(!e->Evaluate({Au, 0, 0}, 0).ok(),
          "HCS/null admitted by replacing magnetic field");
  auto bad = Controls();
  bad.maximumAttempts = 1;
  Require(!e->Move(x, 0, 100, bad).ok(),
          "exhausted work limit silently returned a trajectory");
  auto shorterGuard = d;
  shorterGuard.ambient.support.coverageRadiusM = 8 * Au;
  const auto comparison = Event(shorterGuard)->Move(x, 0, 100, Controls());
  Require(comparison.ok(), comparison.status.message);
  Require(C::Norm(comparison.value.positionM - guard.value.positionM) < .1,
          "8/10-AU guard-only controls changed the physical5-AU boundary trajectory");
  Require(shorterGuard.physicalOuterRadiusM == 5 * Au,
          "guard sensitivity moved the physical escape boundary");
  auto boundary = d;
  boundary.front.r0_m = boundary.physicalOuterRadiusM;
  const auto outerSurface = Event(boundary)->Surface({.1, .15, 1}, 0);
  Require(
      outerSurface.ok() &&
          outerSurface.value.ambient.support == S::PhysicalSupport::ParticleDomain,
      "canonical5AU surface was misclassified as outer guard by projection rounding");
}
void S13SC09() {
  auto d = Fixture();
  d.front.mode = swcme::kinematics::Mode::DBM;
  d.front.V0_m_s = 1000000;
  d.front.Vsw_m_s = 400000;
  d.front.Gamma_m_inv = 1e-11;
  const auto e = Event(d);
  const double t = 100;
  const auto s = e->Surface({0, 0, 1}, t);
  Require(s.ok(), s.status.message);
  const double u = 600000, den = 1 + d.front.Gamma_m_inv * u * t;
  Near(s.value.radiusM, 4 * Rs + 400000 * t + std::log(den) / d.front.Gamma_m_inv,
       2e-13, "independent DBM radius");
  Near(s.value.normalSpeedMPerS, 400000 + u / den, 2e-13, "independent DBM speed");
  d = Fixture();
  d.front.mode = swcme::kinematics::Mode::DataDriven;
  d.front.data_time_s = {0, 600, 1200};
  d.front.data_radius_m = {4 * Rs, 4 * Rs + 6e7, 4 * Rs + 1.2e8};
  const auto table = Event(d);
  const auto middle = table->Surface({0, 0, 1}, 300);
  Require(middle.ok(), middle.status.message);
  Near(middle.value.radiusM, 4 * Rs + 3e7, 2e-13,
       "independent exact linear PCHIP radius");
  Near(middle.value.normalSpeedMPerS, 100000, 2e-13,
       "independent exact linear PCHIP derivative");
}
void S13SC10() {
  const auto plus = Event(Fixture(1)), minus = Event(Fixture(-1));
  for (C::Vec3 x : {1.1 * Rs * C::Unit({.1, .1, 1}), Au * C::Unit({.2, .3, 1}),
                    Au * C::Unit({.2, .3, -1})}) {
    const auto p = plus->Evaluate(x, 100), m = minus->Evaluate(x, 100);
    Require(p.ok() && m.ok(), "both polarity verification lacks actual open support");
    const auto& a = p.value.ambient.primitive;
    const auto& b = m.value.ambient.primitive;
    Near(a.plasma.massDensityKgM3, b.plasma.massDensityKgM3, 2e-12,
         "polarity changed physical mass flux");
    Require(C::Norm(a.magneticFieldT + b.magneticFieldT) <
                2e-12 * C::Norm(a.magneticFieldT),
            "signed B polarity not retained");
    Require(
        C::Norm(a.velocityMPerS - b.velocityMPerS) < 2e-12 * C::Norm(a.velocityMPerS),
        "polarity reversed outward plasma velocity");
    Require(a.magneticSector == -b.magneticSector,
            "actual signed sectors not reversed");
  }
}
}  // namespace
int main(int argc, char** argv) {
  const std::map<std::string, std::function<void()>> tests = {
      {"S13SC01", S13SC01}, {"S13SC02", S13SC02}, {"S13SC03", S13SC03},
      {"S13SC04", S13SC04}, {"S13SC05", S13SC05}, {"S13SC06", S13SC06},
      {"S13SC07", S13SC07}, {"S13SC08", S13SC08}, {"S13SC09", S13SC09},
      {"S13SC10", S13SC10}};
  if (argc == 2 && std::string(argv[1]) == "--list") {
    for (const auto& row : tests)
      std::cout << row.first << '\n';
    return 0;
  }
  const bool all = argc == 2 && std::string(argv[1]) == "--all";
  if (!all &&
      (argc != 3 || std::string(argv[1]) != "--test" || !tests.count(argv[2]))) {
    std::cerr << "--all | --list | --test S13SCxx required\n";
    return 2;
  }
  bool failure = false;
  for (const auto& row : tests)
    if (all || row.first == argv[2]) {
      try {
        row.second();
        std::cout << '[' << row.first << "] PASS\n";
      } catch (const std::exception& e) {
        failure = true;
        std::cerr << '[' << row.first << "] FAIL " << e.what() << '\n';
      }
    }
  return failure ? 1 : 0;
}
