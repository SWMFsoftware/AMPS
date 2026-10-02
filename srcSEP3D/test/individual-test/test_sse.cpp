// Finite-SSE application acceptance. The oracles are the canonical SWCME
// geometry, exact apex intersections and unperturbed ambient fields; these
// tests execute the actual parser, factories, mover and restart reader.
// They do not stand in for a native MPI halo-exchange test.
#include "sep3d_test_registry.h"
#include "configuration_io.h"
#include "parker_geometry.h"
#include "background_factory.h"
#include "bg_swcme.h"
#include "source_runtime.h"
#include "restart.h"
#include "swcme3d_input.hpp"
#include "swcme_shock.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <sstream>

namespace {
namespace A=SEP3D::Adapters; namespace B=SEP3D::Background;
namespace C=SEP3D::Core; namespace R=SEP3D::RuntimeModel;
namespace O=SEP3D::Output;
using Result=SEP3D::Testing::Result;
Result Finish(bool good,const std::string& message) {
  Result r; r.status=good?SEP3D::Testing::Status::Pass:SEP3D::Testing::Status::Fail;
  r.message=message; return r;
}
bool Near(double a,double b,double tolerance=2e-9) {
  return std::fabs(a-b)<=tolerance*std::max({1e-30,std::fabs(a),std::fabs(b)});
}
bool Replace(std::string* text,const std::string& a,const std::string& b) {
  const auto at=text->find(a); if(at==std::string::npos) return false;
  text->replace(at,a.size(),b); return true;
}
std::string Deck() {
  std::ifstream in("examples/sep3d_swcme_sse_mesh_background_20rs_1au.in");
  return std::string(std::istreambuf_iterator<char>(in),{});
}
struct Fixture {
  R::RunConfiguration3DOptions options;
  std::shared_ptr<const R::RunConfiguration3D> configuration;
  std::shared_ptr<B::BackgroundProvider> background;
  std::shared_ptr<A::ShockProvider> shock;
  swcme::input3d::ResolvedConfiguration canonical;
  C::Status Setup(const std::string& deck=Deck()) {
    auto s=R::ParseConfigurationText(deck,&options); if(!s.ok())return s;
    s=R::RunConfiguration3D::Create(options,&configuration); if(!s.ok())return s;
    std::vector<swcme::input3d::Assignment> assignments;
    for(const auto& raw:options.swcmeAssignments) {
      swcme::input3d::Assignment a; a.key=raw.key;a.value=raw.value;assignments.push_back(a);
    }
    const auto resolved=swcme::input3d::Resolve(assignments);
    if(!resolved.ok())return C::Status::Error(resolved.status.message);
    canonical=resolved.configuration;
    s=R::CreateBackgroundProvider(*configuration,&background);if(!s.ok())return s;
    s=background->Prepare(0);if(!s.ok())return s;
    return A::CreateStandaloneSwcmeShockProvider(*configuration,&shock);
  }
};

Result Configuration() {
  Fixture f;auto s=f.Setup();if(!s.ok())return Finish(false,s.message);
  std::string preflight;
  s=R::BuildDryRunSummary(*f.configuration,&preflight);
  if(!s.ok())return Finish(false,"finite-SSE mesh preflight: "+s.message);
  const auto state=f.shock->Evaluate(0);
  bool good=state.status.ok()&&state.geometry==A::ShockGeometryKind::FiniteSSE&&
      Near(state.halfWidthRad,40*C::Const::kPi/180)&&!f.options.source.enabled;
  // Malformed width/axis and unimplemented ellipsoid crossings must still be
  // rejected before mesh allocation. Acceptance is not a blanket guard removal.
  for(const auto& edit:std::vector<std::pair<std::string,std::string>>{
      {"geometry.shape = sse","geometry.shape = ellipsoid"},
      {"geometry.half_width_rad = 0.6981317007977318","geometry.half_width_rad = 0"},
      {"geometry.cme_direction_x = 1","geometry.cme_direction_x = 0"}}) {
    auto text=Deck();R::RunConfiguration3DOptions unused;
    good=good&&Replace(&text,edit.first,edit.second)&&!R::ParseConfigurationText(text,&unused).ok();
  }
  // A wide cap can have a valid apex but cross the 1.05-Rs inner handoff at
  // its tangent flank. Provider preparation must reject the whole epoch.
  auto wide=Deck();Replace(&wide,"geometry.half_width_rad = 0.6981317007977318",
      "geometry.half_width_rad = 1.5533430342749532");
  Fixture rejected;const auto wideStatus=rejected.Setup(wide);
  good=good&&!wideStatus.ok()&&rejected.background&&rejected.background->PreparedMetadata()==nullptr;
  return Finish(good,"SSE parses/freezes; bad axes/widths, ellipsoid and overlapping tangent-flank ejecta are rejected");
}

Result Geometry() {
  Fixture f;auto s=f.Setup();if(!s.ok())return Finish(false,s.message);
  // Rotate the canonical and host records together and use a translated
  // origin. Check radius, normal and normal speed independently over all rays.
  auto params=f.canonical.model;params.cme_dir[0]=0;params.cme_dir[1]=1;
  swcme3d::Model oracle(params);
  for(double time:{0.0,120.0,7200.0}) {
    const auto prepared=oracle.prepare_step(time);
    A::ExpandingShock host=f.shock->Evaluate(time).MoverGeometry();
    host.centerM={4e8,-7e8,3e8};host.cmeDirection={0,1,0};
    for(int angle=0;angle<=180;++angle) {
      const double alpha=angle*C::Const::kPi/180;
      const C::Vec3 u(std::sin(alpha),std::cos(alpha),0);
      double radius,normal[3];const bool exists=oracle.shape_radius_normal(prepared,u.x,u.y,u.z,radius,normal);
      A::ShockSurfacePoint point;s=A::EvaluateShockSurface(host,u,&point);
      if(!s.ok()||point.exists!=exists) return Finish(false,"host/canonical finite angular support mismatch");
      if(!exists)continue;
      const C::Vec3 n(normal);
      const double speed=prepared.V_sh_ms*radius/prepared.r_sh_m*u.Dot(n);
      // Tangent normal speeds can differ by sqrt(roundoff); use an absolute
      // speed bound there rather than a relative comparison to mathematical 0.
      if(!Near(point.radiusM,radius,angle==40?2e-7:2e-9)||
          (point.outwardNormal-n).Norm()>2e-7||
          std::fabs(point.normalSpeedMPerS-speed)>0.05) {
        std::ostringstream detail;detail.precision(17);
        detail<<"SSE geometry mismatch angle="<<angle<<" time="<<time
            <<" radius="<<point.radiusM<<"/"<<radius
            <<" normal_error="<<(point.outwardNormal-n).Norm()
            <<" speed="<<point.normalSpeedMPerS<<"/"<<speed;
        return Finish(false,detail.str());
      }
    }
  }
  return Finish(true,"543 rotated/translated/time-dependent directions agree with canonical SSE radius, normal and normal speed");
}

Result Fields() {
  Fixture f;auto s=f.Setup();if(!s.ok())return Finish(false,s.message);
  swcme3d::Model oracle(f.canonical.model);
  for(double time:{0.0,120.0,7200.0}) {
    s=f.background->Prepare(time);if(!s.ok())return Finish(false,s.message);
    const auto prepared=oracle.prepare_step(time);
    for(double angle:{0.0,15.0,30.0,39.0,41.0,90.0,180.0}) {
      const double alpha=angle*C::Const::kPi/180;
      const C::Vec3 u(std::cos(alpha),std::sin(alpha),0);
      double front,normal[3];const bool exists=oracle.shape_radius_normal(prepared,u.x,u.y,u.z,front,normal);
      const auto at=u*(exists?0.995*front:0.8*prepared.r_sh_m);
      const auto actual=f.background->Evaluate(at);
      double n,vx,vy,vz,bx,by,bz;
      const auto status=oracle.evaluate_cartesian_with_B_checked(prepared,&at.x,&at.y,&at.z,&n,&vx,&vy,&vz,&bx,&by,&bz,1);
      if(!status.ok()||!actual.valid||!Near(n,actual.numberDensityM3)||
          (actual.U-C::Vec3(vx,vy,vz)).Norm()>1e-6||
          (actual.B-C::Vec3(bx,by,bz)).Norm()>1e-16)
        return Finish(false,"finite-SSE mesh field differs from canonical primitive state");
      if(!exists&&(!Near(n,swcme::solarwind::density_m3(prepared.common.solar_wind,at.Norm()))||
          (actual.U-u*prepared.V_sw_ms).Norm()>1e-6))
        return Finish(false,"CME leaked into ambient outside its cap");
      // Full Cartesian derivatives must remain finite, including oblique
      // flank samples where using the apex width would over-size the stencil.
      for(int i=0;i<3;++i)for(int j=0;j<3;++j)
        if(!std::isfinite(actual.gradB(i,j))||!std::isfinite(actual.gradU(i,j)))
          return Finish(false,"SSE Cartesian derivatives are nonfinite");
    }
  }
  return Finish(true,"moving cap mesh primitives match canonical SWCME; outside-cap cells stay ambient and flank derivatives are finite");
}

Result WeakFlankRuntime() {
  // Exact epoch/Cartesian coordinates from the four-rank October 2 failure.
  // The remote point is ambient, while its ray's ACTUAL surface crosses from
  // sub-fast to weakly super-fast between ticks 2 and 3. Exercise both facts:
  // successful remote derivatives cannot hide a failed surface reconstruction.
  Fixture f;auto s=f.Setup();if(!s.ok())return Finish(false,s.message);
  swcme3d::Model model(f.canonical.model);
  for(int tick=0;tick<=10;++tick) {
    const double time=60.0*tick;
    s=f.background->Prepare(time);if(!s.ok())return Finish(false,s.message);
    const auto step=model.prepare_step(time);
    for(double sign:{-1.0,1.0}) {
      const C::Vec3 at(78429796951.338379,-31283096398.984375,
                      sign*2999748969.765625);
      const auto u=at/at.Norm();const double direction[3]={u.x,u.y,u.z};
      swcme3d::LocalShockState shock;
      const auto status=model.shock_state_direction_checked(step,direction,shock);
      if(!status.ok())return Finish(false,"reported SSE flank surface: "+status.summary());
      // Frozen independent 80-digit conserved-flux root, not a value
      // generated by the production cubic. The audit-only verifier is
      // test/reference/verify_sse_weak_flank.py; both mirrored rays agree.
      if(tick==3 && (!shock.has_shock || !shock.solver_converged ||
          !Near(shock.fast_mach,1.000000867696212,2e-12) ||
          std::fabs(shock.compression-1.0000011572267356)>5e-13 ||
          shock.energy_residual>1e-8 || shock.momentum_residual>1e-8))
        return Finish(false,"reported weak flank jump was lost or changed branch");
      // The complete adapter executes its Cartesian finite-difference
      // stencils, including the exact failing coordinates and nearby points.
      for(const auto& d:std::vector<C::Vec3>{{0,0,0},{1e5,0,0},{-1e5,0,0},
          {0,1e5,0},{0,-1e5,0},{0,0,1e5},{0,0,-1e5}}) {
        const auto position=at+d;const auto sample=f.background->Evaluate(position);
        if(!sample.valid || !Near(sample.numberDensityM3,
            swcme::solarwind::density_m3(step.common.solar_wind,position.Norm())) ||
            (sample.U-position.Normalized()*step.V_sw_ms).Norm()>1e-6)
          return Finish(false,"reported remote ambient/stencil point: "+sample.status.message);
      }
      // Query the true layer too; an ambient-only shortcut must not be the
      // sole reason this runtime regression passes. Its inward endpoint must
      // recover the canonical weak RH downstream density and heated pressure.
      const auto bounds=swcme::regions::make_boundaries(shock.Rdir_m,step.region_config);
      const auto layer=f.background->Evaluate(
          u*(shock.Rdir_m-0.5*bounds.smooth_shock_width_m));
      if(!layer.valid || (shock.has_shock &&
          (!Near(layer.numberDensityM3,shock.downstream_n_m3,2e-9) ||
           !Near(layer.pressurePa,shock.downstream.pressure_Pa,2e-9))))
        return Finish(false,"actual finite-CME weak layer: "+layer.status.message);
    }
  }
  return Finish(true,"reported mirrored SSE flank points and actual weak layer pass every epoch/stencil through ten 60-s steps");
}

Result AmbientSupport() {
  // Deliberately retain a genuinely unresolved, roundoff-scale fast shock.
  // Background interfaces must return ambient outside the geometric support,
  // but direct shock diagnostics AND queries inside the ICME must still fail.
  // This negative control prevents a broad catch-and-ambient fallback from
  // making the runtime test pass by erasing physical compression everywhere.
  Fixture f;auto s=f.Setup();if(!s.ok())return Finish(false,s.message);
  const C::Vec3 u=C::Vec3(78429796951.338379,-31283096398.984375,
                         2999748969.765625).Normalized();
  const double direction[3]={u.x,u.y,u.z};
  for(auto shape:{swcme3d::ShockShape::Sphere,swcme3d::ShockShape::SSE,
                  swcme3d::ShockShape::Ellipsoid}) {
    auto params=f.canonical.model;params.shape=shape;
    swcme3d::Model probe(params);const auto initial=probe.prepare_step(0);
    swcme3d::LocalShockState seed;
    if(!probe.shock_state_direction_checked(initial,direction,seed).ok())
      return Finish(false,"support-test setup shock failed");
    const auto normal=swcme::shock::detail::normalized(
        {{seed.normal[0],seed.normal[1],seed.normal[2]}});
    const double fast=swcme::shock::detail::fast_mode_speed(seed.upstream,normal,params.gamma_ad);
    const double factor=seed.Rdir_m/initial.r_sh_m*
        (normal[0]*u.x+normal[1]*u.y+normal[2]*u.z);
    params.V0_sh_kms=(swcme::shock::detail::dot(seed.upstream.velocity_m_s,normal)+
        fast*(1+0.25*swcme::shock::WEAK_SHOCK_MACH_RESOLUTION))/factor/1000;
    swcme3d::Model model(params);const auto step=model.prepare_step(0);
    swcme3d::LocalShockState unresolved;
    const auto rejected=model.shock_state_direction_checked(step,direction,unresolved);
    if(rejected.ok() || !unresolved.has_shock || unresolved.solver_converged ||
       rejected.summary().find("NUMERICALLY_UNRESOLVED_WEAK_SHOCK")==std::string::npos)
      return Finish(false,"roundoff-scale negative control was not an unresolved fast shock");
    const auto bounds=swcme::regions::make_boundaries(seed.Rdir_m,step.region_config);
    for(double radius:{2*seed.Rdir_m,0.5*(bounds.R_te_m-0.5*bounds.smooth_te_width_m)}) {
      const auto point=u*radius;double n=-17,vx=-17,vy=-17,vz=-17,bx=-17,by=-17,bz=-17,rho=-17,p=-17;
      const auto fastStatus=model.evaluate_cartesian_fast_checked(step,&point.x,&point.y,&point.z,&n,&vx,&vy,&vz,1);
      if(!fastStatus.ok() || !Near(n,swcme::solarwind::density_m3(step.common.solar_wind,radius)))
        return Finish(false,"ambient n/V interface solved an irrelevant shock");
      const auto vectorStatus=model.evaluate_cartesian_with_B_checked(step,&point.x,&point.y,&point.z,&n,&vx,&vy,&vz,&bx,&by,&bz,1);
      const auto primitiveStatus=model.evaluate_cartesian_primitive_checked(step,&point.x,&point.y,&point.z,&n,&vx,&vy,&vz,&bx,&by,&bz,&rho,&p,1);
      if(!vectorStatus.ok() || !primitiveStatus.ok() || !(rho>0 && p>0))
        return Finish(false,"ambient vector/thermodynamic interface solved an irrelevant shock");
    }
    const auto point=u*(seed.Rdir_m-0.25*bounds.smooth_shock_width_m);
    double n=-17,vx=-17,vy=-17,vz=-17,bx=-17,by=-17,bz=-17,rho=-17,p=-17;
    const auto failed=model.evaluate_cartesian_primitive_checked(step,&point.x,&point.y,&point.z,
        &n,&vx,&vy,&vz,&bx,&by,&bz,&rho,&p,1);
    if(failed.ok() || n!=-17 || vx!=-17 || bx!=-17 || p!=-17)
      return Finish(false,"unresolved in-CME shock was swallowed or mutated rejected outputs");
  }
  return Finish(true,"all three canonical shapes sample ambient without RH while genuine in-CME failures remain explicit and transactional");
}

Result Intersections() {
  A::ExpandingShock shock;shock.active=true;shock.generation=7;
  shock.geometry=A::ShockGeometryKind::FiniteSSE;
  shock.radiusAtStepStartM=10;shock.halfWidthRad=C::Const::kPi/6;
  auto hit=A::FirstShockIntersection({0.5,0,0},{11,0,0},1,shock,0);
  if(!hit.status.ok()||!hit.crossed||!Near(hit.stepFraction,9.5/10.5)||
      (hit.outwardNormal-C::Vec3(1,0,0)).Norm()>1e-12)
    return Finish(false,"SSE selected rear generating-sphere root instead of outward apex");
  if(A::FirstShockIntersection({2,0,0},{5,0,0},1,shock,0).crossed||
      A::FirstShockIntersection({-11,0,0},{-0.5,0,0},1,shock,0).crossed||
      A::FirstShockIntersection({0.5,0,0},{11,0,0},1,shock,7).crossed)
    return Finish(false,"SSE rear/outside-cap/duplicate-generation crossing was accepted");
  // A stationary particle at the nose is overtaken only when the translating,
  // growing sphere reaches it. R_apex=10+2t gives the independent t=0.5 oracle.
  shock.radialSpeedMPerS=2;
  hit=A::FirstShockIntersection({11,0,0},{11,0,0},1,shock,0);
  if(!hit.status.ok()||!hit.crossed||!Near(hit.stepFraction,0.5)||!Near(hit.normalSpeedMPerS,2))
    return Finish(false,"SSE moving-center intersection/normal speed is wrong");
  // Oblique stationary particle, independently placed on the t=0.5 outward
  // front. The translating sphere must overtake it at the same fraction as
  // the apex, but with a geometry-derived smaller normal speed.
  const double alpha=C::Const::kPi/12;
  const double directionalFactor=(std::cos(alpha)+std::sqrt(0.25-std::pow(std::sin(alpha),2)))/1.5;
  const C::Vec3 flank=C::Vec3(std::cos(alpha),std::sin(alpha),0)*(11*directionalFactor);
  hit=A::FirstShockIntersection(flank,flank,1,shock,0);
  if(!hit.crossed||!Near(hit.stepFraction,0.5)||hit.normalSpeedMPerS<=0||hit.normalSpeedMPerS>=2)
    return Finish(false,"moving SSE flank crossing or normal speed is wrong");
  shock.radialSpeedMPerS=0;
  const C::Vec3 tangentDirection(std::cos(C::Const::kPi/6),std::sin(C::Const::kPi/6),0);
  hit=A::FirstShockIntersection(tangentDirection,tangentDirection*12,1,shock,0);
  if(!hit.crossed||std::fabs(hit.normalSpeedMPerS)>1e-10)
    return Finish(false,"finite SSE tangent boundary is not handled");
  // Large origin and SI distances exercise the scaled quadratic rather than
  // passing only small unit-sphere examples. Also retain the spherical path.
  shock.centerM={3e11,-2e11,7e10};shock.radiusAtStepStartM=1e11;shock.radialSpeedMPerS=0;
  hit=A::FirstShockIntersection(shock.centerM+C::Vec3(0.5e11,0,0),
      shock.centerM+C::Vec3(1.5e11,0,0),60,shock,0);
  if(!hit.crossed||!Near(hit.stepFraction,0.5))return Finish(false,"SI-scale translated SSE crossing failed");
  shock.geometry=A::ShockGeometryKind::Sphere;
  hit=A::FirstShockIntersection(shock.centerM+C::Vec3(-1.5e11,0,0),
      shock.centerM+C::Vec3(-0.5e11,0,0),60,shock,0);
  return Finish(hit.crossed&&Near(hit.stepFraction,0.5),"rear and outside-cap roots are rejected; moving-center/SI-scale crossings and legacy spheres are exact");
}

Result Subcycling() {
  // The full requested-time adapter must retain shape/axis/width when it
  // advances the apex after every substep. A stationary proton is overtaken
  // at t=0.5; a perpendicular ray must never record a spherical surrogate hit.
  A::RequestedTimeAdvance request;auto& input=request.input;
  input.particle.stableId=1;input.particle.species=0;input.particle.statisticalWeight=1;
  input.particle.positionM={11,0,0};input.speciesMassKg=C::Const::m_p;
  input.requestedDtS=1;input.innerRadiusM=0.1;input.outerRadiusM=100;
  input.campaignSeed=11;input.shock.active=true;input.shock.generation=4;
  input.shock.geometry=A::ShockGeometryKind::FiniteSSE;
  input.shock.radiusAtStepStartM=10;input.shock.radialSpeedMPerS=2;
  input.shock.halfWidthRad=C::Const::kPi/6;
  A::LocalTransportRecord local;local.cellSizeM=1000;
  local.background.valid=true;local.background.status=C::Status::OK();
  local.background.B={1,0,0};local.background.bHat={1,0,0};local.background.absB=1;
  request.resolverContext=&local;
  request.resolveLocal=[](const A::ParticleRecord&,double,void* context,A::LocalTransportRecord* out) {
    *out=*static_cast<A::LocalTransportRecord*>(context);return C::Status::OK();
  };
  auto moved=A::AdvanceParticleRequestedTime(request);
  if(!moved.status.ok()||!Near(moved.consumedTimeS,1)||!moved.shockIntersection.crossed||
      moved.particle.lastShockGeneration!=4||moved.acceptedSubsteps<2)
    return Finish(false,"finite requested-time crossing/subcycling failed: "+moved.status.message);
  input.particle.positionM={0,11,0};
  moved=A::AdvanceParticleRequestedTime(request);
  return Finish(moved.status.ok()&&Near(moved.consumedTimeS,1)&&!moved.shockIntersection.crossed&&
      moved.particle.lastShockGeneration==0,"requested-time SSE retains moving geometry, crosses once, consumes the interval and leaves outside-cap particles uncrossed");
}

Result Source() {
  // Exercise actual finite source-surface preparation, not fabricated patch
  // records. FULL_ICME resolved compression remains source-free; injection
  // instead uses the existing SHOCK_ONLY/SOURCE representation.
  // Obtain the normalized finite curve from the propagation fixture BEFORE
  // enabling injection. That input legitimately has no particle observers;
  // the injection configuration below must provide a valid one of its own.
  Fixture base;auto setup=base.Setup();if(!setup.ok())return Finish(false,setup.message);
  const auto& options=base.options;
  C::ParkerSpiralGeometry geometry;
  geometry.sourceRadiusM=options.innerRadiusM;
  geometry.sourceLongitudeRad=options.tubeLongitudeRad;
  geometry.sourceColatitudeRad=options.tubeColatitudeRad;
  geometry.solarWindSpeedMPerS=options.parker.solarWindSpeedMPerS;
  geometry.solarRotationRateRadPerS=options.parker.solarRotationRateRadPerS;
  geometry.rotationAxis=options.parker.rotationAxis;
  double observerRadius=0;
  setup=C::ParkerCurveRadiusAtArcLengthM(0.5*options.parkerSpiralLengthM,
      geometry,&observerRadius);
  if(!setup.ok())return Finish(false,"source probe curve: "+setup.message);
  const auto localPoint=C::ParkerCurvePoint(observerRadius,geometry);
  const auto observerPoint=options.coordinateOriginM+localPoint;
  auto deck=Deck();
  Replace(&deck,"intent = shock-propagation","intent = shock-injection");
  Replace(&deck,"stop_shock_radius_m = 1.57077764235e11","stop_shock_radius_m = 0");
  Replace(&deck,"shock.region_mode = full_icme","shock.region_mode = shock_only");
  Replace(&deck,"shock.acceleration_mode = resolved_compression","shock.acceleration_mode = source");
  Replace(&deck,"provider = swcme","provider = analytic-parker");
  Replace(&deck,"macroparticle_weight = 1e23","macroparticle_weight = 6e24");
  // The disabled flag in this deck belongs to [source]; target its section to
  // avoid accidentally enabling population control or another optional model.
  const auto section=deck.find("[source]");const auto enabled=deck.find("enabled = false",section);
  if(enabled==std::string::npos)return Finish(false,"source example section is missing");
  deck.replace(enabled,std::string("enabled = false").size(),"enabled = true");
  // A generic fixed Earth position from another example need not intersect
  // this Parker corridor. Construct a named TEST probe halfway along the
  // actual finite curve instead, retaining a fixed Cartesian observer and all
  // production containment checks. Precision 17 preserves the calculated SI
  // point through text parsing. This has no bearing on an event's Earth position.
  std::ostringstream probe;probe.precision(17);
  probe<<"\n[observer.sse-source-probe]\n"
       <<"kind = fixed-cartesian\nnormalization = differential-intensity\n"
       <<"position_x_m = "<<observerPoint.x<<'\n'
       <<"position_y_m = "<<observerPoint.y<<'\n'
       <<"position_z_m = "<<observerPoint.z<<'\n'
       <<"follows_trajectory = false\n"
       <<"velocity_x_m_per_s = 0\nvelocity_y_m_per_s = 0\nvelocity_z_m_per_s = 0\n"
       <<"collection_radius_m = "<<0.01*observerRadius<<'\n'
       <<"shell_radius_m = "<<observerRadius<<'\n'
       <<"cadence_s = "<<options.requestedTimeStepS<<'\n'
       <<"energy_bins = 48\nenergy_spacing = logarithmic\npitch_angle_bins = 24\n"
       <<"minimum_energy_j = "<<options.source.minimumEnergyJ<<'\n'
       <<"maximum_energy_j = "<<options.source.maximumEnergyJ<<'\n'
       <<"minimum_mu = -1\nmaximum_mu = 1\nspecies = all\n"
       <<"products = flux,spectrum,anisotropy\n";
  deck+=probe.str();
  Fixture f;auto s=f.Setup(deck);if(!s.ok())return Finish(false,s.message);
  if(f.options.activeRegion==R::ActiveRegionMode::ParkerTube &&
      observerRadius>f.options.activeSolarSphereRadiusM) {
    // Negative control: an opposite-direction point at the same radius is
    // outside this narrow corridor and solar neighborhood. It must still be
    // rejected; accepting the source probe must not weaken the spatial guard.
    auto outside=f.options;
    outside.observers.front().positionM=options.coordinateOriginM-localPoint;
    std::shared_ptr<const R::RunConfiguration3D> unused;
    const auto rejected=R::RunConfiguration3D::Create(outside,&unused);
    if(rejected.ok() || rejected.message.find("does not intersect")==std::string::npos)
      return Finish(false,"outside-corridor observer validation was lost");
  }
  const auto state=f.shock->Evaluate(0);
  if(!state.status.ok()||state.patches.empty())return Finish(false,"canonical finite SSE injection has no active source patches");
  const double minimumCosine=std::cos(state.halfWidthRad);
  for(const auto& patch:state.patches) {
    const auto r=patch.positionM-state.centerM;
    if(r.Normalized().Dot(state.cmeDirection)<minimumCosine-1e-10||
       !patch.active||patch.compression<=1||patch.shockNormalSpeedMPerS<=0)
      return Finish(false,"source patch leaked outside the cap or lacks a physical fast shock");
  }
  return Finish(true,"finite source uses an in-corridor fixed probe and canonical fast-shock patches; outside-corridor observers remain rejected");
}

Result Restart() {
  // Minimal valid checkpoint isolates front geometry from particle sampling.
  // It must survive both schema-4 round trip and schema-3 spherical migration.
  O::RestartState state;state.configurationFingerprint="sse-config";
  state.resolvedConfigurationManifest="sse-manifest";state.storageLayoutFingerprint="layout";
  state.codeIdentity="sse-test";state.snapshotFingerprint="snapshot";
  state.baseTimeStepS=60;state.backgroundGeneration=1;state.campaignSeed=9;
  state.nextStableParticleId=1;state.activeSnapshot.complete=true;state.activeSnapshot.generation=1;
  state.activeSnapshot.authority=R::BackgroundAuthority::Swcme;
  state.shockState.active=true;state.shockState.generation=1;state.shockState.radiusM=20*C::Const::R_sun;
  state.shockState.radialSpeedMPerS=1e6;state.shockState.compressionRatio=2;
  state.shockState.geometry=A::ShockGeometryKind::FiniteSSE;state.shockState.cmeDirection={0,1,0};
  state.shockState.halfWidthRad=0.6981317007977318;state.shockState.providerIdentity="SSE-restart-marker";
  O::RestartLoadOptions options;options.expectedConfigurationFingerprint=state.configurationFingerprint;
  options.expectedCodeIdentity=state.codeIdentity;options.expectedSnapshotFingerprint=state.snapshotFingerprint;
  options.availableBackgroundGeneration=1;
  const auto name="sep3d-sse-"+std::to_string(std::chrono::steady_clock::now().time_since_epoch().count());
  const auto directory=std::filesystem::temp_directory_path()/name;
  std::filesystem::create_directories(directory);const auto file=directory/"state.bin";
  auto s=O::WriteRestart(file.string(),state);O::RestartState loaded;
  if(s.ok())s=O::ReadRestart(file.string(),options,&loaded);
  bool good=s.ok()&&loaded.shockState.geometry==state.shockState.geometry&&
      loaded.shockState.cmeDirection==state.shockState.cmeDirection&&
      loaded.shockState.halfWidthRad==state.shockState.halfWidthRad;
  state.shockState.geometry=A::ShockGeometryKind::Sphere;
  const auto legacy=directory/"legacy.bin";s=O::WriteRestart(legacy.string(),state);
  if(s.ok()) {
    std::ifstream in(legacy,std::ios::binary);std::string bytes(std::istreambuf_iterator<char>(in),{});in.close();
    const auto marker=bytes.find(state.shockState.providerIdentity);
    // Schema 4 added 4 enum bytes plus four IEEE-754 doubles before the
    // provider's 8-byte length. Remove that field block to make a real schema-3
    // fixture, then independently update payload length and FNV checksum.
    if(marker==std::string::npos||marker<44)good=false;
    else {
      bytes.erase(marker-8-36,36);bytes[16]=3;
      const auto put64=[&](std::size_t at,std::uint64_t value) {
        for(unsigned i=0;i<8;++i)bytes[at+i]=static_cast<char>((value>>(8*i))&255);
      };
      put64(8,bytes.size()-24);std::uint64_t hash=UINT64_C(14695981039346656037);
      for(std::size_t i=16;i<bytes.size()-8;++i){hash^=static_cast<unsigned char>(bytes[i]);hash*=UINT64_C(1099511628211);}
      put64(bytes.size()-8,hash);std::ofstream out(legacy,std::ios::binary|std::ios::trunc);
      out.write(bytes.data(),static_cast<std::streamsize>(bytes.size()));out.close();
      loaded=O::RestartState();s=O::ReadRestart(legacy.string(),options,&loaded);
      good=good&&s.ok()&&loaded.shockState.geometry==A::ShockGeometryKind::Sphere;
    }
  } else good=false;
  std::filesystem::remove_all(directory);
  return Finish(good,"schema-4 restores SSE shape/axis/width exactly; a valid schema-3 checkpoint migrates as Sphere");
}
} // namespace

std::vector<SEP3D::Testing::Descriptor> RegisterSseTests() {
  using D=SEP3D::Testing::Descriptor;
  const auto make=[](const char* id,const char* name,SEP3D::Testing::TestCallback callback) {
    D d;d.id=id;d.name=name;d.group="SSE3D";d.description="Finite SSE application acceptance";
    d.initialization=SEP3D::Testing::InitializationLevel::None;d.supportedBuildModes="standalone-no-AMPS";
    d.runtime=SEP3D::Testing::RuntimeClass::Routine;d.seedPolicy="deterministic";
    d.stateIsolation="fresh configuration/provider/checkpoint";d.callback=std::move(callback);return d;
  };
  return {make("SSE3D01","Finite input and inner handoff",Configuration),
      make("SSE3D02","Canonical geometry handoff",Geometry),
      make("SSE3D03","Finite mesh background and gradients",Fields),
      make("SSE3D04","Moving finite-cap particle crossings",Intersections),
      make("SSE3D05","Canonical finite source injection",Source),
      make("SSE3D06","Finite restart and schema-3 migration",Restart),
      make("SSE3D07","Finite requested-time mover subcycling",Subcycling),
      make("SSE3D08","Reported weak-flank runtime regression",WeakFlankRuntime),
      make("SSE3D09","Ambient support and unresolved-layer rejection",AmbientSupport)};
}
