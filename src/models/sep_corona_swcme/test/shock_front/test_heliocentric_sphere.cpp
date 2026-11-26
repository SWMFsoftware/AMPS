// Qualification of the v1.3 reduced front: closed heliocentric sphere,
// data-driven PCHIP history, isothermal Parker wind fitted to U(r_ref) and a
// single |B| magnetogram calibration, with the isothermal wind continued down
// the open PFSS tubes below the 2.5 Rs source surface.  Every oracle below is independent of
// the code path it checks: the SWCME data_driven kinematics evaluates the
// same knots, closed-form sphere geometry checks the mesher, finite
// differences check the analytical derivatives, and a second uncalibrated
// ambient checks the calibration linearity.
#include "provider.h"
#include "sha256.h"
#include "swcme_kinematics.hpp"

#include "sep_coronal_cme/constants.h"
#include "sep_coronal_cme/flux_tube_wind.h"
#include "sep_coronal_cme/pfss_harmonics.h"

#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <map>
#include <sstream>
#include <string>

namespace SF=SEP::CoronaSwcme::ShockFront;
namespace CME=SEP::CoronalCME;

namespace {

int failures=0,passes=0;
void Check(bool condition,const std::string& id,const std::string& message) {
  if(condition){++passes;std::cout<<"PASS "<<id<<" "<<message<<'\n';}
  else{++failures;std::cerr<<"FAIL "<<id<<" "<<message<<'\n';}
}
std::string Sci(double x){std::ostringstream o;o.precision(3);o<<std::scientific<<x;return o.str();}
bool Near(double value,double expected,double relative) {
  return std::abs(value-expected)<=relative*std::max(std::abs(expected),1e-300);
}

constexpr double kRs=6.957e8,kAu=1.495978707e11;
// 2013-04-11 knots (Lario et al. 2014 Table 2 + Wind shock), in SI.
const char* kTimes="0,599.999,719.999,1379.998,225998.393";
const std::string kRadii=[]{
  std::ostringstream o;o.precision(17);
  const double r[]={2.52,3.57,3.73,4.68};
  for(double x:r)o<<x*kRs<<',';
  o<<1.483836225e11;return o.str();
}();
const std::string kHarmonics=
    "degree,order,cosine_T,sine_T\n1,0,1e-4,0\n1,1,2e-5,-1e-5\n2,1,3e-5,1e-5\n";

std::string Event(const std::map<std::string,std::string>& overrides={},
    const std::string& omitted="") {
  std::map<std::string,std::string> v={
    {"schema","shock-front-ambient-v1.3"},
    {"profile","heliocentric-sphere-pfss-parker-pchip-v1"},
    {"run.reference_epoch","2013-04-11T07:26:44Z"},
    {"run.time_scale","UTC-resolved-to-elapsed-SI-seconds"},
    {"run.coordinate_frame","HCI"},{"run.start_s","0"},{"run.end_s","149040"},
    {"run.background_dt_s","600"},{"run.particle_mode","disabled"},
    {"domain.solar_radius_m","695700000"},
    {"domain.first_valid_plasma_radius_m","695700000"},
    {"domain.coverage_radius_m","299195741400"},
    {"ambient.source_surface_radius_m","1739250000"},
    {"ambient.reference_radius_m","149597870700"},
    {"ambient.electron_density_at_reference_m3","6237800"},
    {"ambient.closed_base_electron_density_m3","1e14"},
    {"ambient.rotation_rate_rad_per_s","2.86533e-6"},
    {"ambient.minimum_magnetic_field_t","1e-15"},
    {"ambient.trace_step_m","1e6"},{"ambient.gradient_relative_step","1e-5"},
    {"ambient.trace_maximum_steps","10000"},{"ambient.wind_table_points","1025"},
    {"ambient.open_wind_density_model","isothermal-flux-tube"},
    {"ambient.wind_speed_at_reference_m_s","377948"},
    {"ambient.magnetic_calibration_position_hci_m","7.479893535e10,0,1.295556e11"},
    {"ambient.magnetic_calibration_epoch_s","0"},
    {"ambient.magnetic_calibration_magnitude_t","6.978e-9"},
    {"plasma.gamma_adiabatic","1.6666666666666667"},
    {"plasma.electron_temperature_k","186677"},{"plasma.proton_temperature_k","62502"},
    {"plasma.alpha_temperature_k","606240"},
    {"plasma.alpha_to_proton_number_ratio","0.05426"},
    {"plasma.include_electron_mass","false"},
    {"geometry.model","heliocentric-sphere"},
    {"geometry.physical_inner_radius_m","765270000"},
    {"history.model","data-driven-pchip"},{"history.knot_times_s",kTimes},
    {"history.knot_radii_m",kRadii},{"history.validity_policy","report-only"},
    {"diagnostics.criticality_model","none"},{"numerics.polar_cells","12"},
    {"numerics.azimuth_cells","24"},{"numerics.surface_topology","triangular-sphere-v1"},
    {"numerics.weak_mach_tolerance","1e-6"},{"numerics.rh_residual_tolerance","1e-8"},
    {"endpoint.apex_radius_m","1e11"},
    {"endpoint.observer_position_hci_m","7.479893535e10,0,1.295556e11"},
    {"endpoint.require_fast_shock_at_observer","false"},
    {"assets.harmonics_file","h.csv"},
    {"assets.harmonics_sha256",CME::Internal::Sha256Hex(kHarmonics)}};
  for(const auto& item:overrides)v[item.first]=item.second;
  std::ostringstream text;text<<"#begin physics: shock_front\n";
  for(const auto& item:v)if(item.first!=omitted)
    text<<"  "<<item.first<<" = "<<item.second<<'\n';
  text<<"#end physics\n";return text.str();
}

SEP::Core::Result<std::shared_ptr<const SF::Configuration>> Resolve(
    const std::string& text) {
  return SF::ResolveConfiguration(text,[](const std::string& path) {
    if(path!="h.csv")return SEP::Core::Result<std::string>::Failure(
        SEP::Core::StatusCode::InvalidConfiguration,"unknown asset");
    return SEP::Core::Result<std::string>::Success(kHarmonics);
  },"hsp");
}

std::shared_ptr<SF::Provider> Make(const std::string& text) {
  const auto c=Resolve(text);
  if(!c.ok()){std::cerr<<c.status.message<<'\n';return nullptr;}
  const auto p=SF::Provider::Create(c.value);
  if(!p.ok()){std::cerr<<p.status.message<<'\n';return nullptr;}
  return p.value;
}

} // namespace

int main() {
  const auto provider=Make(Event());
  Check(provider!=nullptr,"HSP00","v1.3 event resolves and creates its provider");
  if(!provider)return 1;
  const auto& c=provider->Event();

  // HSP01: PCHIP equals the independent SWCME data_driven kinematics.
  {
    swcme::kinematics::Config k;k.mode=swcme::kinematics::Mode::DataDriven;
    k.data_time_s=c.knotTimesS;k.data_radius_m=c.knotRadiiM;
    // SWCME validates its (unused here) ballistic fields for every mode.
    k.r0_m=c.knotRadiiM.front();k.Vsw_m_s=4e5;
    double worstR=0,worstV=0;bool knots=true;
    for(std::size_t i=0;i<c.knotTimesS.size();++i) {
      const auto s=provider->Trajectory(std::min(c.knotTimesS[i],149040.0));
      if(c.knotTimesS[i]<=149040)knots=knots&&s.ok()&&s.value.apexRadiusM==c.knotRadiiM[i];
    }
    for(int i=0;i<=400;++i) {
      const double t=149040.0*i/400;
      const auto s=provider->Trajectory(t);
      const auto o=swcme::kinematics::evaluate(k,t);
      if(!s.ok()||o.status!=swcme::kinematics::Status::Ok){worstR=1;break;}
      worstR=std::max(worstR,std::abs(s.value.apexRadiusM-o.radius_m)/o.radius_m);
      worstV=std::max(worstV,std::abs(s.value.apexSpeedMPerS-o.speed_m_s)/
          std::max(1.0,o.speed_m_s));
    }
    Check(knots&&worstR<=1e-14&&worstV<=1e-12,"HSP01",
        "PCHIP passes through every knot and matches SWCME data_driven (rel R "+
        std::to_string(worstR)+", V "+std::to_string(worstV)+")");
  }
  // HSP02: analytical speed and acceleration equal finite differences.
  {
    double worstV=0,worstA=0;
    for(double t:{100.0,650.0,1000.0,5e3,5e4,1.2e5}) {
      const double h=1e-2;
      const auto m=provider->Trajectory(t-h),z=provider->Trajectory(t),p=provider->Trajectory(t+h);
      if(!m.ok()||!z.ok()||!p.ok()){worstV=1;break;}
      worstV=std::max(worstV,std::abs((p.value.apexRadiusM-m.value.apexRadiusM)/(2*h)-
          z.value.apexSpeedMPerS)/z.value.apexSpeedMPerS);
      worstA=std::max(worstA,std::abs((p.value.apexSpeedMPerS-m.value.apexSpeedMPerS)/(2*h)-
          z.value.apexAccelerationMPerS2)/std::max(1e-3,std::abs(z.value.apexAccelerationMPerS2)));
    }
    Check(worstV<1e-7&&worstA<1e-4,"HSP02","dR/dt and d2R/dt2 equal central differences");
  }
  // HSP03: no extrapolation and strict knot admission.
  {
    const bool outside=!provider->Trajectory(-1).ok()&&!provider->Trajectory(149041).ok();
    const bool late=!Resolve(Event({{"run.end_s","226000"}})).ok();
    const bool start=!Resolve(Event({{"history.knot_times_s","1,599.999,719.999,1379.998,225998.393"}})).ok();
    const bool monotone=!Resolve(Event({{"history.knot_radii_m","2e9,1e9,3e9,4e9,1.4e11"}})).ok();
    const bool unknown=!Resolve(Event({{"geometry.half_width_rad","1"}})).ok();
    const bool missing=!Resolve(Event({},"ambient.wind_speed_at_reference_m_s")).ok();
    const bool tube=!Resolve(Event({{"ambient.open_wind_density_model","field-aligned-flux-tube"}})).ok()&&
        !Resolve(Event({{"ambient.open_wind_density_model","spherical-radial-mass-flux"}})).ok();
    Check(outside&&late&&start&&monotone&&unknown&&missing&&tube,"HSP03",
        "out-of-knot epochs, uncovered runs, shifted/decreasing knots, unknown/missing keys and kinematic open-tube closures are rejected");
  }
  // HSP04: closed-sphere mesh is exact.
  {
    const double t=900;const auto epoch=provider->EvaluateEpoch(t,1);
    const auto apex=provider->Trajectory(t);
    bool ok=epoch.ok()&&apex.ok();
    if(ok) {
      const auto& e=*epoch.value;const double R=apex.value.apexRadiusM,V=apex.value.apexSpeedMPerS;
      ok=e.triangles.size()==2u*24u*11u&&e.records.size()==e.triangles.size()&&
          e.vertices.size()==11u*24u+2u;
      double worstNormal=0,worstSpeed=0,worstRadius=0;bool outward=true,edge=false;
      for(std::size_t i=0;ok&&i<e.records.size();++i) {
        const auto& g=e.records[i].geometry;
        worstNormal=std::max(worstNormal,CME::Norm(g.outwardNormal-g.positionM/CME::Norm(g.positionM)));
        worstRadius=std::max(worstRadius,std::abs(CME::Norm(g.positionM)-R)/R);
        worstSpeed=std::max(worstSpeed,std::abs(g.normalSpeedMPerS-V)/V);
        edge=edge||g.supportEdge;
        const auto& tri=e.triangles[i];
        const auto x0=e.vertices[tri.vertex[0]].positionM,x1=e.vertices[tri.vertex[1]].positionM,
            x2=e.vertices[tri.vertex[2]].positionM;
        outward=outward&&CME::Dot(CME::Cross(x1-x0,x2-x0),g.outwardNormal)>0;
      }
      for(const auto& vertex:e.vertices)
        worstRadius=std::max(worstRadius,std::abs(CME::Norm(vertex.positionM)-R)/R);
      ok=ok&&outward&&!edge&&worstNormal<1e-15&&worstSpeed==0&&worstRadius<1e-15&&
          Near(e.area.geometricSupportM2,4*CME::Constants::kPi*R*R,1e-14)&&
          e.vertices.back().apex&&
          CME::Norm(e.vertices.back().positionM-R*c.direction)<=1e-15*R;
    }
    Check(ok,"HSP04","closed sphere: 2Naz(Npol-1) outward faces, radial normals, Vn=dR/dt, area 4 pi R^2, observer-facing apex");
  }
  // HSP05: generating sphere and front-point membership.
  {
    const auto sphere=provider->Sphere(900);const auto apex=provider->Trajectory(900);
    bool ok=sphere.ok()&&apex.ok()&&CME::Norm(sphere.value.centerM)==0&&
        sphere.value.radiusM==apex.value.apexRadiusM&&sphere.value.minimumRayCosine==-1;
    if(ok) {
      const CME::Vec3 ray=CME::Unit({-0.3,0.8,-0.52});
      const auto on=provider->EvaluateFrontPoint(apex.value.apexRadiusM*ray,900);
      const auto off=provider->EvaluateFrontPoint(1.001*apex.value.apexRadiusM*ray,900);
      ok=on.ok()&&!off.ok()&&on.value.geometry.normalSpeedMPerS==apex.value.apexSpeedMPerS;
    }
    Check(ok,"HSP05","generating sphere is the front; on-surface points classify and off-surface points are rejected");
  }
  // HSP06: wind closure decoupled from the measured temperatures.
  {
    const auto& ambient=provider->Ambient();
    const auto ref=ambient.Evaluate(provider->Event().observerPositionM/
        CME::Norm(provider->Event().observerPositionM)*kAu,0);
    bool ok=ref.ok();
    double radialSpeed=0;
    if(ok){const CME::Vec3 n=CME::Unit(provider->Event().observerPositionM);
      radialSpeed=CME::Dot(ref.value.velocityMPerS,n);}
    ok=ok&&Near(radialSpeed,377948,1e-9);
    const auto hot=Make(Event({{"plasma.proton_temperature_k","1e6"}}));
    double worstU=0,worstFlux=0;bool pressureDiffers=true;
    for(double r:{1.2*kRs,3*kRs,20*kRs,0.5*kAu,1.5*kAu}) {
      const CME::Vec3 x=r*CME::Unit({0.5,0.2,0.84});
      const auto a=ambient.Evaluate(x,0);const auto b=hot?hot->Ambient().Evaluate(x,0):a;
      if(!a.ok()||!b.ok()){ok=false;break;}
      const CME::Vec3 n=x/r;
      const double ua=CME::Dot(a.value.velocityMPerS,n),ub=CME::Dot(b.value.velocityMPerS,n);
      worstU=std::max(worstU,std::abs(ua-ub)/ua);
      // rho U r^2 is the exterior invariant; the open shell conserves
      // rho u/|B| along tubes instead (HSP09).
      const bool exterior=a.value.region==SEP::CoronaSwcme::AmbientRegion::ParkerExterior;
      if(exterior)worstFlux=std::max(worstFlux,std::abs(a.value.plasma.massDensityKgM3*ua*r*r/
          ambient.MassFluxPerSteradianKgPerS()-1));
      pressureDiffers=pressureDiffers&&b.value.plasma.pressurePa>a.value.plasma.pressurePa&&
          (exterior||a.value.region==SEP::CoronaSwcme::AmbientRegion::PfssOpen)&&
          exterior==(r>c.ambient.ambient.sourceSurfaceRadiusM);
    }
    Check(ok&&hot&&worstU==0&&worstFlux<1e-12&&pressureDiffers,"HSP06",
        "U(1 AU) is the declared speed; U(r) is independent of measured T while pressure follows it; rho U r^2 is conserved");
  }
  // HSP07: one-factor |B| calibration is exact and linear.
  {
    const auto& ambient=provider->Ambient();
    const auto at=ambient.Evaluate(c.ambient.ambient.magneticCalibrationPositionM,0);
    const auto doubled=Make(Event({{"ambient.magnetic_calibration_magnitude_t","1.3956e-8"}}));
    bool ok=at.ok()&&doubled&&Near(CME::Norm(at.value.magneticFieldT),6.978e-9,1e-12)&&
        Near(doubled->Ambient().MagneticCalibrationFactor(),2*ambient.MagneticCalibrationFactor(),1e-12);
    double worst=0;
    for(double r:{1.15*kRs,4*kRs,0.3*kAu,1.2*kAu})for(double t:{0.0,7e4}) {
      const CME::Vec3 x=r*CME::Unit({-0.4,0.7,0.3});
      const auto a=ambient.Evaluate(x,t);const auto b=doubled?doubled->Ambient().Evaluate(x,t):a;
      if(!a.ok()||!b.ok()){ok=false;break;}
      worst=std::max(worst,CME::Norm(b.value.magneticFieldT-2*a.value.magneticFieldT)/
          CME::Norm(a.value.magneticFieldT));
      // In the open shell U and rho depend on B only through |B|/|B_r,ss|
      // and b, which the scaled coefficients reproduce to roundoff.
      const double tolerance=r<c.ambient.ambient.sourceSurfaceRadiusM?1e-12:0;
      ok=ok&&a.value.magneticSector==b.value.magneticSector&&
          CME::Norm(a.value.velocityMPerS-b.value.velocityMPerS)<=
              tolerance*CME::Norm(a.value.velocityMPerS)&&
          std::abs(a.value.plasma.massDensityKgM3-b.value.plasma.massDensityKgM3)<=
              tolerance*a.value.plasma.massDensityKgM3;
    }
    Check(ok&&worst<1e-13,"HSP07","calibrated |B| reproduces the declared magnitude; B scales linearly, sector/U/rho unchanged");
  }
  // HSP08: legacy profiles keep their wind closure (no reference-speed fit).
  {
    SEP::CoronaSwcme::AmbientDefinition d=provider->Ambient().Definition();
    // The measured 1-AU temperatures give r_c~46 Rs, which has no resolved
    // subsonic branch at the photosphere: that is why v1.3 decouples the wind.
    // A coronal-temperature legacy definition checks the inactive path.
    d.ambient.windSpeedAtReferenceMPerS=0;d.ambient.magneticCalibration=false;
    d.composition.electronTemperatureK=d.composition.protonTemperatureK=
        d.composition.alphaTemperatureK=1.4e6;
    const auto measured=SEP::CoronaSwcme::AmbientModel::Create(provider->Ambient().Definition());
    auto cold=provider->Ambient().Definition();cold.ambient.windSpeedAtReferenceMPerS=0;
    const bool coldRejected=!SEP::CoronaSwcme::AmbientModel::Create(cold).ok();
    const auto legacy=SEP::CoronaSwcme::AmbientModel::Create(d);
    bool ok=legacy.ok();
    if(ok) {
      const auto& comp=d.composition;
      // Composition sound speed p/rho for the measured partition.
      std::vector<CME::IonSpecies> ions={{"H+",1,1,CME::Constants::kProtonMassKg,comp.protonTemperatureK},
          {"He++",comp.alphaToProtonNumberRatio,2,CME::Constants::kAlphaMassKg,comp.alphaTemperatureK}};
      const auto plasma=CME::EvaluatePlasmaFromElectronDensity(d.ambient.electronDensityAtReferenceM3,
          comp.electronTemperatureK,ions,comp.includeElectronMass,comp.gammaAdiabatic,0);
      ok=plasma.ok()&&Near(legacy.value->WindSoundSpeedSquaredM2PerS2(),
          plasma.value.pressurePa/plasma.value.massDensityKgM3,1e-15)&&
          legacy.value->MagneticCalibrationFactor()==1.0;
    }
    Check(ok&&measured.ok()&&coldRejected,"HSP08",
        "inactive closures reproduce the legacy composition-temperature wind and unit magnetogram; the measured-T legacy wind is rejected, not repaired");
  }
  // HSP09: the open-shell wind is the steady corotating-frame isothermal
  // Parker wind of each PFSS tube. Independent oracles along one traced open
  // field line: corotating-frame flow parallel to B, rho u/|B| constant
  // (continuity), the isothermal Bernoulli integral u^2/2+c^2 ln rho-GM/r
  // constant (field-aligned momentum), subsonic flow, and continuity of B,
  // U and rho with the Parker exterior at the source surface.
  {
    const auto& ambient=provider->Ambient();
    const double omega=c.ambient.ambient.rotationRateRadPerS;
    const double rss=c.ambient.ambient.sourceSurfaceRadiusM;
    const double c2=ambient.WindSoundSpeedSquaredM2PerS2();
    const double gm=CME::Constants::kSolarGravitationalParameterM3PerS2;
    const double t=5000;
    const auto frame=[&](const CME::Vec3& x,const CME::Vec3& u) {
      return u-CME::Vec3{-omega*x.y,omega*x.x,0};
    };
    CME::Vec3 x=1.15*kRs*CME::Unit({0.12,0.07,0.99});
    auto first=ambient.Evaluate(x,t);
    bool ok=first.ok()&&first.value.region==SEP::CoronaSwcme::AmbientRegion::PfssOpen;
    double kMin=1e300,kMax=0,eMin=1e300,eMax=-1e300,worstAlign=0,worstMach=0;
    int samples=0;
    // Outward along sector*B by RK4 (step 2e5 m), sampling every 50 steps.
    for(int i=0;ok&&CME::Norm(x)<0.999*rss&&i<200000;++i) {
      const auto at=ambient.Evaluate(x,t);
      if(!at.ok()||at.value.region!=SEP::CoronaSwcme::AmbientRegion::PfssOpen){ok=false;break;}
      if(i%50==0) {
        const auto& a=at.value;const double r=CME::Norm(x);
        const CME::Vec3 w=frame(x,a.velocityMPerS);
        const double u=CME::Norm(w),b=CME::Norm(a.magneticFieldT),rho=a.plasma.massDensityKgM3;
        worstAlign=std::max(worstAlign,CME::Norm(CME::Cross(w,a.magneticFieldT))/(u*b));
        const double k=rho*u/b,e=(0.5*u*u+c2*std::log(rho)-gm/r)/c2;
        kMin=std::min(kMin,k);kMax=std::max(kMax,k);eMin=std::min(eMin,e);eMax=std::max(eMax,e);
        worstMach=std::max(worstMach,u*u/c2);++samples;
      }
      const auto direction=[&](const CME::Vec3& p) {
        const auto s=ambient.Evaluate(p,t);
        return s.ok()?(s.value.magneticSector*1.0)*CME::Unit(s.value.magneticFieldT):CME::Vec3{0,0,0};
      };
      const double h=2e5;
      const auto k1=direction(x),k2=direction(x+0.5*h*k1),k3=direction(x+0.5*h*k2),k4=direction(x+h*k3);
      x=x+(h/6)*(k1+2*k2+2*k3+k4);
    }
    ok=ok&&samples>20;
    // Source-surface continuity on the radial ray through the traced foot.
    const CME::Vec3 n=CME::Unit(x);
    const auto below=ambient.Evaluate((rss*(1-1e-9))*n,t),above=ambient.Evaluate((rss*(1+1e-9))*n,t);
    double jumpB=1,jumpU=1,jumpRho=1;
    if(below.ok()&&above.ok()) {
      jumpB=CME::Norm(below.value.magneticFieldT-above.value.magneticFieldT)/CME::Norm(above.value.magneticFieldT);
      jumpU=CME::Norm(below.value.velocityMPerS-above.value.velocityMPerS)/CME::Norm(above.value.velocityMPerS);
      jumpRho=std::abs(below.value.plasma.massDensityKgM3/above.value.plasma.massDensityKgM3-1);
      ok=ok&&below.value.region==SEP::CoronaSwcme::AmbientRegion::PfssOpen&&
          above.value.region==SEP::CoronaSwcme::AmbientRegion::ParkerExterior&&
          below.value.magneticSector==above.value.magneticSector;
    } else ok=false;
    Check(ok&&worstAlign<1e-12&&kMax/kMin-1<1e-6&&eMax-eMin<1e-6&&worstMach<1&&
        jumpB<1e-6&&jumpU<1e-6&&jumpRho<1e-6,"HSP09",
        "open tube: u'||B, rho u/|B| and isothermal Bernoulli constant, subsonic, B/U/rho continuous at Rss ("+
        std::to_string(samples)+" samples, dK "+Sci(kMax/kMin-1)+", dE/c2 "+Sci(eMax-eMin)+
        ", align "+Sci(worstAlign)+", max M^2 "+Sci(worstMach)+", jumps B/U/rho "+Sci(jumpB)+"/"+Sci(jumpU)+"/"+Sci(jumpRho)+")");
  }
  // HSP10: the tube closure rejects a wind supersonic at the source surface
  // and cannot be combined with the radial rho U r^2 open-shell law.
  {
    auto fast=provider->Ambient().Definition();fast.ambient.windSpeedAtReferenceMPerS=1.5e6;
    auto mixed=provider->Ambient().Definition();mixed.ambient.sphericalOpenMassFlux=true;
    auto legacy=provider->Ambient().Definition();legacy.ambient.isothermalFluxTubeOpenWind=false;
    const auto f=SEP::CoronaSwcme::AmbientModel::Create(fast);
    const bool fastRejected=!f.ok();
    const bool mixedRejected=!SEP::CoronaSwcme::AmbientModel::Create(mixed).ok();
    const auto l=SEP::CoronaSwcme::AmbientModel::Create(legacy);
    // The legacy kinematic branch keeps u=U(r) along the tube.
    bool legacyKinematic=l.ok();
    if(legacyKinematic) {
      const CME::Vec3 x=1.15*kRs*CME::Unit({0.12,0.07,0.99});
      const auto a=l.value->Evaluate(x,0);
      const double omega=legacy.ambient.rotationRateRadPerS;
      const auto radial=CME::SolveRadialIsothermalParker({CME::Norm(x)},
          std::sqrt(l.value->WindSoundSpeedSquaredM2PerS2()),l.value->CriticalRadiusM(),1);
      legacyKinematic=a.ok()&&radial.ok()&&Near(CME::Norm(a.value.velocityMPerS-
          CME::Vec3{-omega*x.y,omega*x.x,0}),radial.value[0].speedMPerS,1e-6);
    }
    Check(fastRejected&&mixedRejected&&legacyKinematic,"HSP10",
        "r_c<=Rss and the mixed closure are rejected; the inactive flag keeps the kinematic open-tube speed");
  }
  // HSP11: the recurrence PFSS evaluator used by v1.3 equals the historical
  // per-term evaluator to roundoff over the whole shell.
  {
    const auto& pfssRows=c.ambient.ambient.harmonics;
    std::vector<CME::HarmonicCoefficient> rows;
    {std::istringstream in(pfssRows);std::string row;
     while(std::getline(in,row,';')){if(row.empty())continue;std::replace(row.begin(),row.end(),',',' ');
       std::istringstream r(row);CME::HarmonicCoefficient h;r>>h.degree>>h.order>>h.cosineT>>h.sineT;rows.push_back(h);}}
    const auto pfss=CME::PfssHarmonics::Create(kRs,c.ambient.ambient.sourceSurfaceRadiusM,rows);
    double worst=0;bool ok=pfss.ok()&&!rows.empty();
    unsigned state=12345;auto uniform=[&]{state=1103515245u*state+12345u;return (state>>8)/16777216.0;};
    for(int i=0;ok&&i<2000;++i) {
      const double r=kRs+(c.ambient.ambient.sourceSurfaceRadiusM-kRs)*uniform();
      const double z=2*uniform()-1,phi=6.283185307179586*uniform(),q=std::sqrt(1-z*z);
      const CME::Vec3 x=r*CME::Vec3{q*std::cos(phi),q*std::sin(phi),z};
      const auto a=pfss.value.EvaluateCartesian(x),b=pfss.value.EvaluateFieldCartesian(x);
      if(!a.ok()||!b.ok()){ok=false;break;}
      worst=std::max(worst,CME::Norm(a.value-b.value)/CME::Norm(a.value));
    }
    Check(ok&&worst<1e-12,"HSP11","recurrence PFSS evaluator equals the per-term evaluator (max rel "+Sci(worst)+")");
  }
  // HSP12: the radius-parametrized open-tube foot agrees with the general
  // arc-length trace. The legacy kinematic flux-tube closure (generic Trace,
  // 1e6 m steps) and v1.3 (radius steps of 2e7 m) must give the same mass
  // flux per unit flux K=Fm/(Rss^2|B_r,ss|) at every open shell point; a
  // closed equatorial point stays closed and is rejected as a material state.
  {
    const auto coarse=Make(Event({{"ambient.trace_step_m","2e7"}}));
    auto legacyDefinition=provider->Ambient().Definition();
    legacyDefinition.ambient.isothermalFluxTubeOpenWind=false;
    const auto legacy=SEP::CoronaSwcme::AmbientModel::Create(legacyDefinition);
    bool ok=coarse&&legacy.ok();
    double worst=0;int open=0;
    const double omega=c.ambient.ambient.rotationRateRadPerS;
    for(int i=0;ok&&i<40;++i) {
      const double r=(1.05+1.4*i/40.0)*kRs;
      const CME::Vec3 x=r*CME::Unit({0.1+0.01*i,0.05,0.99});
      const auto a=coarse->Ambient().Evaluate(x,0),b=legacy.value->Evaluate(x,0);
      if(!a.ok()||!b.ok()){ok=false;break;}
      if(a.value.region!=b.value.region){ok=false;break;}
      if(a.value.region!=SEP::CoronaSwcme::AmbientRegion::PfssOpen)continue;
      ++open;
      const auto along=[&](const SEP::CoronaSwcme::AmbientPrimitive& p) {
        return CME::Norm(p.velocityMPerS-CME::Vec3{-omega*x.y,omega*x.x,0});
      };
      const double ka=a.value.plasma.massDensityKgM3*along(a.value)/CME::Norm(a.value.magneticFieldT);
      const double kb=b.value.plasma.massDensityKgM3*along(b.value)/CME::Norm(b.value.magneticFieldT);
      worst=std::max(worst,std::abs(ka/kb-1));
    }
    const CME::Vec3 equator=1.2*kRs*CME::Unit({1,0.3,0.02});
    const auto closedA=coarse?coarse->Ambient().Evaluate(equator,0):legacy.value->Evaluate(equator,0);
    const auto closedB=legacy.ok()?legacy.value->Evaluate(equator,0):closedA;
    ok=ok&&open>20&&closedA.ok()&&closedB.ok()&&
        closedA.value.region==SEP::CoronaSwcme::AmbientRegion::PfssClosed&&
        closedB.value.region==SEP::CoronaSwcme::AmbientRegion::PfssClosed;
    Check(ok&&worst<1e-6,"HSP12","radius-parametrized foot reproduces the general trace's K on "+
        std::to_string(open)+" open points (max rel "+Sci(worst)+"); closed equator agrees");
  }
  std::cout<<"heliocentric sphere: "<<passes<<" pass, "<<failures<<" fail\n";
  return failures==0?0:1;
}
