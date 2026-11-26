#include "provider.h"

#include "sha256.h"

#include "amps_input.h"
#include "swcme_kinematics.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <iomanip>
#include <limits>
#include <map>
#include <set>
#include <sstream>

namespace SEP { namespace CoronaSwcme { namespace ShockFront { namespace {

template<class T>
Core::Result<T> Bad(const std::string& message) {
  return Core::Result<T>::Failure(Core::StatusCode::InvalidConfiguration,message);
}

std::string Trim(std::string value) {
  const auto first=value.find_first_not_of(" \t\r\n");
  if(first==std::string::npos)return {};
  const auto last=value.find_last_not_of(" \t\r\n");
  return value.substr(first,last-first+1);
}

// Read the generic AMPS runtime-input text (src/INPUT_FILE.md) that defines
// one reduced shock-front event and return its closed key/value set.
//
// The event is one block, ``#begin physics: shock_front ... #end physics``,
// either as a stand-alone file (srcSEP --reduced-shock-event, viewers, the
// shared tests) or included into an application input (srcSEP3D), in which
// case the application passes the canonical text of that block.  The model
// vocabulary is dotted (``ambient.reference_radius_m``); a nested block
// ``#begin ambient`` is the same key space and is flattened to the dotted
// spelling (AMPS::Input::Flatten), so both spellings name one parameter.
//
// The event file is part of the physical problem definition, not a
// permissive user-preference file.  The generic grammar already rejects
// duplicate keys inside one block; flattening can still meet the same dotted
// key twice (``ambient.x`` and ``#begin ambient`` + ``x``), which is rejected
// here: accepting "last value wins" would let the runtime state depend on
// parser ordering while two ranks or restart readers appeared to share the
// same event.  Every diagnostic carries ``file:line`` of the user's text.
Core::Result<std::map<std::string,std::string>> Assignments(
    const std::string& bytes,const std::string& displayName) {
  using Return=Core::Result<std::map<std::string,std::string>>;
  AMPS::Input::Document document;
  const AMPS::Input::Status parsed=
      AMPS::Input::ParseText(bytes,displayName,".",&document);
  if(!parsed.ok())return Bad<std::map<std::string,std::string>>(parsed.message);
  const AMPS::Input::Status blocks=AMPS::Input::RequireOnlyBlocks(
      document.root,{"physics"},"the reduced shock-front event reader");
  if(!blocks.ok())return Bad<std::map<std::string,std::string>>(blocks.message);
  const AMPS::Input::Block* event=nullptr;
  for(const AMPS::Input::Block& block:document.root.blocks) {
    if(block.name!="shock_front")return Bad<std::map<std::string,std::string>>(
        AMPS::Input::ErrorAt(block.begin,"block '"+block.Label()+
            "' is not a reduced shock-front event; expected "
            "'#begin physics: shock_front'").message);
    event=&block;
  }
  if(event==nullptr)return Bad<std::map<std::string,std::string>>(
      displayName+": no '#begin physics: shock_front' block was found");
  std::map<std::string,std::string> values;
  for(const AMPS::Input::Assignment& item:AMPS::Input::Flatten(*event)) {
    if(!values.emplace(item.key,item.value).second)
      return Bad<std::map<std::string,std::string>>(AMPS::Input::ErrorAt(
          item.where,"duplicate shock-front event key '"+item.key+
          "' (a dotted key and a nested block name the same parameter)",
          item.text).message);
  }
  return Return::Success(std::move(values));
}

Core::Result<double> Number(const std::map<std::string,std::string>& values,
    const std::string& key) {
  const auto found=values.find(key);if(found==values.end())return Bad<double>(
      "missing shock-front event key '"+key+"'");
  std::istringstream input(found->second);double value=0;std::string trailing;
  if(!(input>>value)||input>>trailing||!std::isfinite(value))return Bad<double>(
      "shock-front event key '"+key+"' is not one finite SI number");
  return Core::Result<double>::Success(value);
}

Core::Result<int> Integer(const std::map<std::string,std::string>& values,
    const std::string& key) {
  const auto value=Number(values,key);if(!value.ok())return Bad<int>(value.status.message);
  if(value.value<std::numeric_limits<int>::min()||
      value.value>std::numeric_limits<int>::max()||
      value.value!=std::floor(value.value))return Bad<int>(
          "shock-front event key '"+key+"' is not an integer");
  return Core::Result<int>::Success(static_cast<int>(value.value));
}

Core::Result<std::string> Text(const std::map<std::string,std::string>& values,
    const std::string& key) {
  const auto found=values.find(key);if(found==values.end())return Bad<std::string>(
      "missing shock-front event key '"+key+"'");
  return Core::Result<std::string>::Success(found->second);
}

Core::Result<bool> Boolean(const std::map<std::string,std::string>& values,
    const std::string& key) {
  const auto value=Text(values,key);if(!value.ok())return Bad<bool>(value.status.message);
  if(value.value=="true")return Core::Result<bool>::Success(true);
  if(value.value=="false")return Core::Result<bool>::Success(false);
  return Bad<bool>("shock-front event key '"+key+"' must be true or false");
}

Core::Result<CoronalCME::Vec3> Vector(
    const std::map<std::string,std::string>& values,const std::string& key) {
  const auto value=Text(values,key);if(!value.ok())return Bad<CoronalCME::Vec3>(
      value.status.message);
  std::string bytes=value.value;std::replace(bytes.begin(),bytes.end(),',',' ');
  CoronalCME::Vec3 result;std::istringstream input(bytes);std::string trailing;
  if(!(input>>result.x>>result.y>>result.z)||input>>trailing||
      !std::isfinite(result.x)||!std::isfinite(result.y)||!std::isfinite(result.z))
    return Bad<CoronalCME::Vec3>("shock-front event key '"+key+
        "' must contain three finite Cartesian components");
  return Core::Result<CoronalCME::Vec3>::Success(result);
}

Core::Result<std::vector<double>> List(
    const std::map<std::string,std::string>& values,const std::string& key) {
  const auto value=Text(values,key);if(!value.ok())return Bad<std::vector<double>>(
      value.status.message);
  std::string bytes=value.value;std::replace(bytes.begin(),bytes.end(),',',' ');
  std::istringstream input(bytes);std::vector<double> result;double x=0;
  while(input>>x)result.push_back(x);
  if(!input.eof()||result.empty())return Bad<std::vector<double>>(
      "shock-front event key '"+key+"' must be a comma-separated list of finite SI numbers");
  for(double item:result)if(!std::isfinite(item))return Bad<std::vector<double>>(
      "shock-front event key '"+key+"' contains a nonfinite value");
  return Core::Result<std::vector<double>>::Success(std::move(result));
}

Core::Result<std::string> Harmonics(const std::string& bytes) {
  // Convert the checksummed CSV asset into the exact compact representation
  // consumed by the maintained PFSS kernel.  Coefficients are photospheric
  // radial-field amplitudes in tesla in that kernel's declared real-harmonic
  // convention.  The byte checksum is checked separately below; this
  // canonicalization is therefore not being used as a lossy integrity check.
  std::istringstream input(bytes);std::string line;std::ostringstream canonical;
  bool any=false;
  while(std::getline(input,line)) {
    const auto comment=line.find('#');if(comment!=std::string::npos)line.erase(comment);
    line=Trim(line);if(line.empty())continue;
    if(line=="degree,order,cosine_T,sine_T")continue;
    std::replace(line.begin(),line.end(),',',' ');
    int degree=0,order=0;double cosine=0,sine=0;std::string trailing;
    std::istringstream row(line);
    if(!(row>>degree>>order>>cosine>>sine)||row>>trailing||degree<1||
        std::abs(order)>degree||degree>32||!std::isfinite(cosine)||
        !std::isfinite(sine))return Bad<std::string>(
            "harmonic asset row is not l,m,cosine_T,sine_T");
    if(any)canonical<<';';
    canonical<<degree<<','<<order<<','<<std::setprecision(17)<<cosine<<','<<sine;
    any=true;
  }
  if(!any)return Bad<std::string>("harmonic asset contains no coefficients");
  return Core::Result<std::string>::Success(canonical.str());
}

bool StrictUtcEpoch(const std::string& value) {
  // The reduced provider currently converts an absolute event label to
  // elapsed SI seconds outside this kernel.  Accept exactly the UTC form that
  // makes that conversion unambiguous; accepting arbitrary nonempty text here
  // would let two clocks share one apparent HCI state.  Leap seconds (`:60`)
  // require an external timescale table and are intentionally unsupported by
  // this first profile rather than silently treated as a normal SI second.
  if(value.size()!=20||value[4]!='-'||value[7]!='-'||value[10]!='T'||
      value[13]!=':'||value[16]!=':'||value[19]!='Z')return false;
  for(std::size_t i=0;i<value.size();++i)if(i!=4&&i!=7&&i!=10&&i!=13&&
      i!=16&&i!=19&&!std::isdigit(static_cast<unsigned char>(value[i])))return false;
  const int year=std::stoi(value.substr(0,4));
  const int month=std::stoi(value.substr(5,2));
  const int day=std::stoi(value.substr(8,2));
  const int hour=std::stoi(value.substr(11,2));
  const int minute=std::stoi(value.substr(14,2));
  const int second=std::stoi(value.substr(17,2));
  if(year<1||month<1||month>12||hour>23||minute>59||second>59)return false;
  static const int days[]={0,31,28,31,30,31,30,31,31,30,31,30,31};
  int maximum=days[month];
  if(month==2&&(year%400==0||(year%4==0&&year%100!=0)))maximum=29;
  return day>=1&&day<=maximum;
}

} // namespace

const char* Name(HistoryModel value) noexcept {
  switch(value) {
    case HistoryModel::ConstantSpeed:return "constant-speed";
    case HistoryModel::QuinticPulseThenConstant:return "quintic-pulse-then-constant";
    case HistoryModel::FittedQuadratic:return "fitted-quadratic";
    case HistoryModel::DataDrivenPchip:return "data-driven-pchip";
  }
  return "unknown";
}
const char* Name(GeometryModel value) noexcept {
  return value==GeometryModel::HeliocentricSphere?"heliocentric-sphere":
      "finite-sse-fixed-width";
}
const char* Name(Phase value) noexcept {
  return value==Phase::CoronalHistory?"coronal-history":"swcme-outer";
}
const char* Name(FrontStatus value) noexcept {
  switch(value) {
    case FrontStatus::OutsideFrontSupport:return "outside-front-support";
    case FrontStatus::BelowPhysicalInnerBoundary:return "below-physical-inner-boundary";
    case FrontStatus::AmbientUnavailable:return "ambient-unavailable";
    case FrontStatus::NonForwardInflow:return "non-forward-inflow";
    case FrontStatus::SubfastFront:return "subfast-front";
    case FrontStatus::SolvedFastShock:return "solved-fast-shock";
    case FrontStatus::NumericallyUnresolvedWeakShock:return "numerically-unresolved-weak-shock";
    case FrontStatus::WrongBranch:return "wrong-branch";
    case FrontStatus::InvalidJump:return "invalid-jump";
  }
  return "unknown";
}
const char* Name(HtDiagnosticStatus value) noexcept {
  switch(value) {
    case HtDiagnosticStatus::Valid:return "valid";
    case HtDiagnosticStatus::UpstreamMagneticNull:return "upstream-magnetic-null";
    case HtDiagnosticStatus::PerpendicularNoFiniteBoost:
      return "perpendicular-no-finite-boost";
    case HtDiagnosticStatus::NumericalFailure:return "numerical-failure";
  }
  return "unknown";
}

namespace {

// Schema v1.3: a closed Sun-centred front through observed height-time knots
// over the canonical PFSS/Parker ambient whose wind is the isothermal Parker
// transonic solution fitted to U(r_ref) and whose magnetogram amplitude is
// fixed by one |B| measurement.  The key set is closed exactly as in v1.1/v1.2.
Core::Result<std::shared_ptr<const Configuration>> ResolveHeliocentricSphere(
    const std::map<std::string,std::string>& v,const AssetReader& reader) {
  using Return=Core::Result<std::shared_ptr<const Configuration>>;
  const std::set<std::string> required={
    "schema","profile","run.reference_epoch","run.time_scale",
    "run.coordinate_frame","run.start_s","run.end_s","run.background_dt_s",
    "run.particle_mode","domain.solar_radius_m",
    "domain.first_valid_plasma_radius_m","domain.coverage_radius_m",
    "ambient.source_surface_radius_m","ambient.reference_radius_m",
    "ambient.electron_density_at_reference_m3",
    "ambient.closed_base_electron_density_m3","ambient.rotation_rate_rad_per_s",
    "ambient.minimum_magnetic_field_t","ambient.trace_step_m",
    "ambient.gradient_relative_step","ambient.trace_maximum_steps",
    "ambient.wind_table_points","ambient.open_wind_density_model",
    "ambient.wind_speed_at_reference_m_s",
    "ambient.magnetic_calibration_position_hci_m",
    "ambient.magnetic_calibration_epoch_s",
    "ambient.magnetic_calibration_magnitude_t",
    "plasma.gamma_adiabatic",
    "plasma.electron_temperature_k","plasma.proton_temperature_k",
    "plasma.alpha_temperature_k","plasma.alpha_to_proton_number_ratio",
    "plasma.include_electron_mass","geometry.model",
    "geometry.physical_inner_radius_m","history.model",
    "history.knot_times_s","history.knot_radii_m","history.validity_policy",
    "diagnostics.criticality_model","numerics.polar_cells",
    "numerics.azimuth_cells","numerics.surface_topology",
    "numerics.weak_mach_tolerance","numerics.rh_residual_tolerance",
    "endpoint.apex_radius_m","endpoint.observer_position_hci_m",
    "endpoint.require_fast_shock_at_observer",
    "assets.harmonics_file","assets.harmonics_sha256"};
  for(const auto& item:v)if(!required.count(item.first))return Return::Failure(
      Core::StatusCode::InvalidConfiguration,
      "unknown shock-front event key '"+item.first+"'");
  for(const auto& key:required)if(!v.count(key))return Return::Failure(
      Core::StatusCode::InvalidConfiguration,
      "missing shock-front event key '"+key+"'");
  auto out=std::make_shared<Configuration>();
#define SF_TEXT(KEY,TARGET) do { const auto x=Text(v,KEY); if(!x.ok()) return Return::Failure(x.status.code,x.status.message); TARGET=x.value; } while(false)
#define SF_NUMBER(KEY,TARGET) do { const auto x=Number(v,KEY); if(!x.ok()) return Return::Failure(x.status.code,x.status.message); TARGET=x.value; } while(false)
#define SF_INTEGER(KEY,TARGET) do { const auto x=Integer(v,KEY); if(!x.ok()) return Return::Failure(x.status.code,x.status.message); TARGET=x.value; } while(false)
  SF_TEXT("schema",out->schema);SF_TEXT("profile",out->profile);
  SF_TEXT("run.reference_epoch",out->referenceEpoch);SF_TEXT("run.time_scale",out->timeScale);
  SF_TEXT("run.coordinate_frame",out->coordinateFrame);SF_TEXT("run.particle_mode",out->particleMode);
  SF_NUMBER("run.start_s",out->ambient.support.startS);SF_NUMBER("run.end_s",out->ambient.support.endS);
  SF_NUMBER("run.background_dt_s",out->backgroundDtS);
  SF_NUMBER("domain.solar_radius_m",out->ambient.support.solarRadiusM);
  SF_NUMBER("domain.first_valid_plasma_radius_m",out->ambient.support.firstValidPlasmaRadiusM);
  SF_NUMBER("domain.coverage_radius_m",out->ambient.support.coverageRadiusM);
  auto& a=out->ambient.ambient;
  SF_NUMBER("ambient.source_surface_radius_m",a.sourceSurfaceRadiusM);
  SF_NUMBER("ambient.reference_radius_m",a.referenceRadiusM);
  SF_NUMBER("ambient.electron_density_at_reference_m3",a.electronDensityAtReferenceM3);
  SF_NUMBER("ambient.closed_base_electron_density_m3",a.closedBaseElectronDensityM3);
  SF_NUMBER("ambient.rotation_rate_rad_per_s",a.rotationRateRadPerS);
  SF_NUMBER("ambient.minimum_magnetic_field_t",a.minimumMagneticFieldT);
  SF_NUMBER("ambient.trace_step_m",a.traceStepM);
  SF_NUMBER("ambient.gradient_relative_step",a.gradientRelativeStep);
  SF_INTEGER("ambient.trace_maximum_steps",a.traceMaximumSteps);
  SF_INTEGER("ambient.wind_table_points",a.windTablePoints);
  SF_NUMBER("ambient.wind_speed_at_reference_m_s",a.windSpeedAtReferenceMPerS);
  {const auto x=Vector(v,"ambient.magnetic_calibration_position_hci_m");
   if(!x.ok())return Return::Failure(x.status.code,x.status.message);
   a.magneticCalibrationPositionM=x.value;a.magneticCalibration=true;}
  SF_NUMBER("ambient.magnetic_calibration_epoch_s",a.magneticCalibrationEpochS);
  SF_NUMBER("ambient.magnetic_calibration_magnitude_t",a.magneticCalibrationMagnitudeT);
  {const auto x=Text(v,"ambient.open_wind_density_model");
   if(!x.ok())return Return::Failure(x.status.code,x.status.message);
   // Material lines run through open PFSS tubes below the source surface.
   // Only the tube continuation of the exterior isothermal Parker wind
   // (continuity per unit flux + isothermal Bernoulli along the tube) closes
   // mass and field-aligned momentum there; the radial rho U r^2 law does not.
   if(x.value!="isothermal-flux-tube")return Return::Failure(
       Core::StatusCode::UnsupportedCapability,
       "v1.3 requires ambient.open_wind_density_model=isothermal-flux-tube");
   a.sphericalOpenMassFlux=false;a.isothermalFluxTubeOpenWind=true;}
  SF_NUMBER("plasma.gamma_adiabatic",out->ambient.composition.gammaAdiabatic);
  SF_NUMBER("plasma.electron_temperature_k",out->ambient.composition.electronTemperatureK);
  SF_NUMBER("plasma.proton_temperature_k",out->ambient.composition.protonTemperatureK);
  SF_NUMBER("plasma.alpha_temperature_k",out->ambient.composition.alphaTemperatureK);
  SF_NUMBER("plasma.alpha_to_proton_number_ratio",out->ambient.composition.alphaToProtonNumberRatio);
  {const auto x=Boolean(v,"plasma.include_electron_mass");if(!x.ok())return Return::Failure(x.status.code,x.status.message);out->ambient.composition.includeElectronMass=x.value;}
  SF_NUMBER("geometry.physical_inner_radius_m",out->physicalInnerRadiusM);
  {const auto x=List(v,"history.knot_times_s");if(!x.ok())return Return::Failure(x.status.code,x.status.message);out->knotTimesS=x.value;}
  {const auto x=List(v,"history.knot_radii_m");if(!x.ok())return Return::Failure(x.status.code,x.status.message);out->knotRadiiM=x.value;}
  SF_INTEGER("numerics.polar_cells",out->polarCells);
  SF_INTEGER("numerics.azimuth_cells",out->azimuthCells);
  SF_TEXT("numerics.surface_topology",out->surfaceTopology);
  SF_NUMBER("numerics.weak_mach_tolerance",out->weakMachTolerance);
  SF_NUMBER("numerics.rh_residual_tolerance",out->rhResidualTolerance);
  SF_NUMBER("endpoint.apex_radius_m",out->endpointRadiusM);
  {const auto x=Vector(v,"endpoint.observer_position_hci_m");if(!x.ok())return Return::Failure(x.status.code,x.status.message);out->observerPositionM=x.value;}
  {const auto x=Boolean(v,"endpoint.require_fast_shock_at_observer");if(!x.ok())return Return::Failure(x.status.code,x.status.message);out->requireFastShockAtObserver=x.value;}
  SF_TEXT("assets.harmonics_file",out->harmonicAssetPath);
  SF_TEXT("assets.harmonics_sha256",out->harmonicAssetSha256);
#undef SF_TEXT
#undef SF_NUMBER
#undef SF_INTEGER
  const auto geometry=Text(v,"geometry.model");const auto history=Text(v,"history.model");
  const auto criticality=Text(v,"diagnostics.criticality_model");
  const auto policy=Text(v,"history.validity_policy");
  if(!geometry.ok()||geometry.value!="heliocentric-sphere"||!history.ok()||
      history.value!="data-driven-pchip"||!criticality.ok()||criticality.value!="none")
    return Return::Failure(Core::StatusCode::UnsupportedCapability,
        "v1.3 requires heliocentric-sphere, data-driven-pchip and criticality_model=none");
  out->geometryModel=GeometryModel::HeliocentricSphere;
  out->historyModel=HistoryModel::DataDrivenPchip;
  if(!policy.ok()||(policy.value!="report-only"&&
      policy.value!="require-declared-shock-coverage"))return Return::Failure(
          Core::StatusCode::InvalidConfiguration,"invalid trajectory validity policy");
  out->validityPolicy=policy.value=="report-only"?TrajectoryValidityPolicy::ReportOnly:
      TrajectoryValidityPolicy::RequireDeclaredShockCoverage;

  // Observed height-time knots: strictly increasing in both time and radius
  // (an outward front), covering the whole declared run without extrapolation.
  const auto& t=out->knotTimesS;const auto& r=out->knotRadiiM;
  bool knots=t.size()>=2&&t.size()==r.size();
  for(std::size_t i=1;knots&&i<t.size();++i)knots=t[i]>t[i-1]&&r[i]>r[i-1];
  if(!knots||t.front()!=out->ambient.support.startS||
      t.back()<out->ambient.support.endS||!(r.front()>0))
    return Return::Failure(Core::StatusCode::InvalidConfiguration,
        "v1.3 knots must be strictly increasing, start at run.start_s and cover run.end_s");
  out->knotSlopesMPerS=swcme::kinematics::pchip_slopes(t,r);
  for(double slope:out->knotSlopesMPerS)if(!std::isfinite(slope)||slope<0)
    return Return::Failure(Core::StatusCode::InvalidConfiguration,
        "v1.3 knots produce a non-monotone PCHIP front");

  const double observerNorm=CoronalCME::Norm(out->observerPositionM);
  if(!std::isfinite(observerNorm)||!(observerNorm>0))return Return::Failure(
      Core::StatusCode::InvalidConfiguration,"v1.3 observer position must be finite and nonzero");
  out->direction=out->observerPositionM/observerNorm;
  out->halfWidthRad=3.14159265358979323846;
  out->initialApexRadiusM=r.front();out->initialApexSpeedMPerS=out->knotSlopesMPerS.front();
  out->finalPulseSpeedMPerS=out->knotSlopesMPerS.back();out->handoffApexRadiusM=r.back();

  const auto harmonicBytes=reader(out->harmonicAssetPath);
  if(!harmonicBytes.ok())return Return::Failure(harmonicBytes.status.code,
      "harmonic asset '"+out->harmonicAssetPath+"': "+harmonicBytes.status.message);
  if(CoronalCME::Internal::Sha256Hex(harmonicBytes.value)!=out->harmonicAssetSha256)
    return Return::Failure(Core::StatusCode::DataIntegrityFailure,
        "harmonic asset checksum differs from the declared SHA-256");
  const auto harmonics=Harmonics(harmonicBytes.value);
  if(!harmonics.ok())return Return::Failure(harmonics.status.code,harmonics.status.message);
  a.harmonics=harmonics.value;
  out->ambient.coordinateFrame=out->coordinateFrame;

  const bool valid=out->profile=="heliocentric-sphere-pfss-parker-pchip-v1"&&
      StrictUtcEpoch(out->referenceEpoch)&&
      out->timeScale=="UTC-resolved-to-elapsed-SI-seconds"&&
      out->coordinateFrame=="HCI"&&out->particleMode=="disabled"&&
      out->physicalInnerRadiusM>=out->ambient.support.firstValidPlasmaRadiusM&&
      r.front()>=out->physicalInnerRadiusM&&
      out->backgroundDtS>0&&out->polarCells>=2&&out->azimuthCells>=4&&
      out->surfaceTopology=="triangular-sphere-v1"&&
      out->weakMachTolerance>0&&out->rhResidualTolerance>0&&
      out->endpointRadiusM>0&&
      out->endpointRadiusM<=out->ambient.support.coverageRadiusM;
  if(!valid)return Return::Failure(Core::StatusCode::InvalidConfiguration,
      "shock-front event violates the frozen v1.3 heliocentric-sphere domain");

  std::ostringstream manifest;manifest<<std::setprecision(17);
  for(const auto& item:v)manifest<<item.first<<'='<<item.second<<'\n';
  manifest<<"assets.harmonics_content_sha256="<<out->harmonicAssetSha256<<'\n';
  out->normalizedManifest=manifest.str();
  out->physicsFingerprint=CoronalCME::Internal::Sha256Hex(out->normalizedManifest);
  out->ambient.physicsFingerprint=out->physicsFingerprint;
  return Return::Success(std::move(out));
}

} // namespace

Core::Result<std::shared_ptr<const Configuration>> ResolveConfiguration(
    const std::string& inputBytes,const AssetReader& reader,
    const std::string& displayName) {
  using Return=Core::Result<std::shared_ptr<const Configuration>>;
  if(!reader)return Return::Failure(Core::StatusCode::InvalidConfiguration,
      "shock-front event requires an asset reader");
  const auto parsed=Assignments(inputBytes,displayName);
  if(!parsed.ok())return Return::Failure(parsed.status.code,parsed.status.message);
  const auto& v=parsed.value;
  {
    const auto schema=Text(v,"schema");
    if(schema.ok()&&schema.value=="shock-front-ambient-v1.3")
      return ResolveHeliocentricSphere(v,reader);
  }
  // Freeze the complete first-profile schema.  In particular, silently
  // ignoring a misspelled physical key is more dangerous than rejecting a
  // future key: it could change the front trajectory or ambient normalization
  // without changing the apparent input.  Later profiles need an explicit
  // schema/version branch rather than relaxing this closed set.
  const auto schema=Text(v,"schema");
  const bool fitted=schema.ok()&&schema.value=="shock-front-ambient-v1.2";
  std::set<std::string> required={
    "schema","profile","run.reference_epoch","run.time_scale",
    "run.coordinate_frame","run.start_s","run.end_s","run.background_dt_s",
    "run.particle_mode","domain.solar_radius_m",
    "domain.first_valid_plasma_radius_m","domain.coverage_radius_m",
    "ambient.source_surface_radius_m","ambient.reference_radius_m",
    "ambient.electron_density_at_reference_m3",
    "ambient.closed_base_electron_density_m3","ambient.rotation_rate_rad_per_s",
    "ambient.minimum_magnetic_field_t","ambient.trace_step_m",
    "ambient.gradient_relative_step","ambient.trace_maximum_steps",
    "ambient.wind_table_points","ambient.open_wind_density_model",
    "plasma.gamma_adiabatic",
    "plasma.electron_temperature_k","plasma.proton_temperature_k",
    "plasma.alpha_temperature_k","plasma.alpha_to_proton_number_ratio",
    "plasma.include_electron_mass","geometry.model","geometry.direction_hci",
    "geometry.half_width_rad","geometry.physical_inner_radius_m",
    "history.model","history.initial_apex_radius_m",
    "history.initial_apex_speed_m_s","history.final_pulse_speed_m_s",
    "history.acceleration_duration_s","handoff.apex_radius_m",
    "handoff.outer_model","handoff.drag_gamma_m_inv",
    "handoff.effective_wind_speed_m_s","handoff.trajectory_validity_policy",
    "diagnostics.criticality_model","numerics.polar_cells",
    "numerics.azimuth_cells","numerics.surface_topology",
    "numerics.weak_mach_tolerance",
    "numerics.rh_residual_tolerance","endpoint.apex_radius_m",
    "endpoint.observer_position_hci_m","endpoint.require_fast_shock_at_observer",
    "assets.harmonics_file","assets.harmonics_sha256"};
  // The fitted schema contains a single authoritative polynomial and
  // handoff epoch. The analytically derived handoff radius/speed must
  // not also be supplied as independently editable configuration.
  if(fitted) {
    for(const char* key:{"history.final_pulse_speed_m_s","history.acceleration_duration_s",
        "handoff.apex_radius_m"})required.erase(key);
    required.insert({"history.fit_end_s","history.apex_acceleration_m_s2",
        "history.minimum_handoff_radius_m","history.drag_horizon_s",
        "history.prepared_event_sha256"});
  }
  for(const auto& item:v)if(!required.count(item.first))return Return::Failure(
      Core::StatusCode::InvalidConfiguration,
      "unknown shock-front event key '"+item.first+"'");
  for(const auto& key:required)if(!v.count(key))return Return::Failure(
      Core::StatusCode::InvalidConfiguration,
      "missing shock-front event key '"+key+"'");

  auto out=std::make_shared<Configuration>();
#define SF_TEXT(KEY,TARGET) do { const auto x=Text(v,KEY); if(!x.ok()) return Return::Failure(x.status.code,x.status.message); TARGET=x.value; } while(false)
#define SF_NUMBER(KEY,TARGET) do { const auto x=Number(v,KEY); if(!x.ok()) return Return::Failure(x.status.code,x.status.message); TARGET=x.value; } while(false)
#define SF_INTEGER(KEY,TARGET) do { const auto x=Integer(v,KEY); if(!x.ok()) return Return::Failure(x.status.code,x.status.message); TARGET=x.value; } while(false)
  SF_TEXT("schema",out->schema);SF_TEXT("profile",out->profile);
  SF_TEXT("run.reference_epoch",out->referenceEpoch);SF_TEXT("run.time_scale",out->timeScale);
  SF_TEXT("run.coordinate_frame",out->coordinateFrame);SF_TEXT("run.particle_mode",out->particleMode);
  SF_NUMBER("run.start_s",out->ambient.support.startS);SF_NUMBER("run.end_s",out->ambient.support.endS);
  SF_NUMBER("run.background_dt_s",out->backgroundDtS);
  SF_NUMBER("domain.solar_radius_m",out->ambient.support.solarRadiusM);
  SF_NUMBER("domain.first_valid_plasma_radius_m",out->ambient.support.firstValidPlasmaRadiusM);
  SF_NUMBER("domain.coverage_radius_m",out->ambient.support.coverageRadiusM);
  SF_NUMBER("ambient.source_surface_radius_m",out->ambient.ambient.sourceSurfaceRadiusM);
  SF_NUMBER("ambient.reference_radius_m",out->ambient.ambient.referenceRadiusM);
  SF_NUMBER("ambient.electron_density_at_reference_m3",out->ambient.ambient.electronDensityAtReferenceM3);
  SF_NUMBER("ambient.closed_base_electron_density_m3",out->ambient.ambient.closedBaseElectronDensityM3);
  SF_NUMBER("ambient.rotation_rate_rad_per_s",out->ambient.ambient.rotationRateRadPerS);
  SF_NUMBER("ambient.minimum_magnetic_field_t",out->ambient.ambient.minimumMagneticFieldT);
  SF_NUMBER("ambient.trace_step_m",out->ambient.ambient.traceStepM);
  SF_NUMBER("ambient.gradient_relative_step",out->ambient.ambient.gradientRelativeStep);
  SF_INTEGER("ambient.trace_maximum_steps",out->ambient.ambient.traceMaximumSteps);
  SF_INTEGER("ambient.wind_table_points",out->ambient.ambient.windTablePoints);
  {const auto x=Text(v,"ambient.open_wind_density_model");
   if(!x.ok())return Return::Failure(x.status.code,x.status.message);
   if(x.value=="spherical-radial-mass-flux")
     out->ambient.ambient.sphericalOpenMassFlux=true;
   else if(x.value=="field-aligned-flux-tube")
     out->ambient.ambient.sphericalOpenMassFlux=false;
   else return Return::Failure(Core::StatusCode::UnsupportedCapability,
       "unsupported ambient.open_wind_density_model");}
  SF_NUMBER("plasma.gamma_adiabatic",out->ambient.composition.gammaAdiabatic);
  SF_NUMBER("plasma.electron_temperature_k",out->ambient.composition.electronTemperatureK);
  SF_NUMBER("plasma.proton_temperature_k",out->ambient.composition.protonTemperatureK);
  SF_NUMBER("plasma.alpha_temperature_k",out->ambient.composition.alphaTemperatureK);
  SF_NUMBER("plasma.alpha_to_proton_number_ratio",out->ambient.composition.alphaToProtonNumberRatio);
  {const auto x=Boolean(v,"plasma.include_electron_mass");if(!x.ok())return Return::Failure(x.status.code,x.status.message);out->ambient.composition.includeElectronMass=x.value;}
  {const auto x=Vector(v,"geometry.direction_hci");if(!x.ok())return Return::Failure(x.status.code,x.status.message);out->direction=x.value;}
  SF_NUMBER("geometry.half_width_rad",out->halfWidthRad);
  SF_NUMBER("geometry.physical_inner_radius_m",out->physicalInnerRadiusM);
  SF_NUMBER("history.initial_apex_radius_m",out->initialApexRadiusM);
  SF_NUMBER("history.initial_apex_speed_m_s",out->initialApexSpeedMPerS);
  if(fitted) {
    double minimumHandoffM;
    SF_NUMBER("history.fit_end_s",out->fittedEndS);
    SF_NUMBER("history.apex_acceleration_m_s2",out->fittedAccelerationMPerS2);
    SF_NUMBER("history.minimum_handoff_radius_m",minimumHandoffM);
    SF_NUMBER("history.drag_horizon_s",out->fittedDragHorizonS);
    const auto digest=Text(v,"history.prepared_event_sha256");
    const auto terminal=EvaluateFittedApex(out->fittedEndS,
        out->initialApexRadiusM,out->initialApexSpeedMPerS,out->fittedAccelerationMPerS2);
    // Speed is linear through the fit. Checking both endpoints prevents
    // an inward excursion before the matched quadratic-drag continuation.
    if(!digest.ok()||digest.value.size()!=64||
        digest.value.find_first_not_of("0123456789abcdef")!=std::string::npos||
        out->ambient.support.startS!=0||!(out->fittedEndS>0)||
        !(minimumHandoffM>0)||!std::isfinite(terminal.apexRadiusM)||
        !(terminal.apexRadiusM>=minimumHandoffM)||!(out->initialApexSpeedMPerS>0)||
        !std::isfinite(terminal.apexSpeedMPerS)||!(terminal.apexSpeedMPerS>0)||
        !(out->ambient.support.endS>=out->fittedEndS)||!(out->fittedDragHorizonS>=0)||
        out->ambient.support.endS-out->fittedEndS>out->fittedDragHorizonS)
      return Return::Failure(Core::StatusCode::InvalidConfiguration,
          "invalid fitted SSE clock, provenance, advancing history or drag coverage");
    out->handoffApexRadiusM=terminal.apexRadiusM;
    out->finalPulseSpeedMPerS=terminal.apexSpeedMPerS;
  } else {
    SF_NUMBER("history.final_pulse_speed_m_s",out->finalPulseSpeedMPerS);
    SF_NUMBER("history.acceleration_duration_s",out->accelerationDurationS);
    SF_NUMBER("handoff.apex_radius_m",out->handoffApexRadiusM);
  }
  SF_NUMBER("handoff.drag_gamma_m_inv",out->dragGammaPerM);
  SF_NUMBER("handoff.effective_wind_speed_m_s",out->effectiveTrajectoryWindMPerS);
  SF_INTEGER("numerics.polar_cells",out->polarCells);
  SF_INTEGER("numerics.azimuth_cells",out->azimuthCells);
  SF_TEXT("numerics.surface_topology",out->surfaceTopology);
  SF_NUMBER("numerics.weak_mach_tolerance",out->weakMachTolerance);
  SF_NUMBER("numerics.rh_residual_tolerance",out->rhResidualTolerance);
  SF_NUMBER("endpoint.apex_radius_m",out->endpointRadiusM);
  {const auto x=Vector(v,"endpoint.observer_position_hci_m");if(!x.ok())return Return::Failure(x.status.code,x.status.message);out->observerPositionM=x.value;}
  {const auto x=Boolean(v,"endpoint.require_fast_shock_at_observer");if(!x.ok())return Return::Failure(x.status.code,x.status.message);out->requireFastShockAtObserver=x.value;}
  SF_TEXT("assets.harmonics_file",out->harmonicAssetPath);
  SF_TEXT("assets.harmonics_sha256",out->harmonicAssetSha256);
#undef SF_TEXT
#undef SF_NUMBER
#undef SF_INTEGER

  const auto geometry=Text(v,"geometry.model");const auto history=Text(v,"history.model");
  const auto outer=Text(v,"handoff.outer_model");const auto policy=Text(v,"handoff.trajectory_validity_policy");
  const auto criticality=Text(v,"diagnostics.criticality_model");
  if(!geometry.ok()||geometry.value!="finite-sse-fixed-width"||!outer.ok()||
      outer.value!="quadratic-drag"||!criticality.ok()||criticality.value!="none")
    return Return::Failure(Core::StatusCode::UnsupportedCapability,
        "first reduced profile requires finite-sse-fixed-width, quadratic-drag and criticality_model=none");
  if(!history.ok()||(fitted? history.value!="fitted-quadratic":
      history.value!="constant-speed"&&history.value!="quintic-pulse-then-constant"))
    return Return::Failure(
          Core::StatusCode::UnsupportedCapability,"unsupported front history model");
  out->historyModel=fitted?HistoryModel::FittedQuadratic:
      (history.value=="constant-speed"?HistoryModel::ConstantSpeed:
       HistoryModel::QuinticPulseThenConstant);
  if(!policy.ok()||(policy.value!="report-only"&&
      policy.value!="require-declared-shock-coverage"))return Return::Failure(
          Core::StatusCode::InvalidConfiguration,"invalid trajectory validity policy");
  out->validityPolicy=policy.value=="report-only"?TrajectoryValidityPolicy::ReportOnly:
      TrajectoryValidityPolicy::RequireDeclaredShockCoverage;

  const auto harmonicBytes=reader(out->harmonicAssetPath);
  if(!harmonicBytes.ok())return Return::Failure(harmonicBytes.status.code,
      "harmonic asset '"+out->harmonicAssetPath+"': "+harmonicBytes.status.message);
  if(CoronalCME::Internal::Sha256Hex(harmonicBytes.value)!=out->harmonicAssetSha256)
    return Return::Failure(Core::StatusCode::DataIntegrityFailure,
        "harmonic asset checksum differs from the declared SHA-256");
  const auto harmonics=Harmonics(harmonicBytes.value);
  if(!harmonics.ok())return Return::Failure(harmonics.status.code,harmonics.status.message);
  out->ambient.ambient.harmonics=harmonics.value;
  out->ambient.coordinateFrame=out->coordinateFrame;

  // This block is the physical admissibility boundary, not merely syntax
  // validation.  All dimensional quantities are SI; direction is an HCI unit
  // vector; the finite SSE half width is below pi/2 so its leading ray root is
  // unique; and the endpoint must lie inside the ambient's declared radial
  // support.  No clamp or floor repairs an event outside these assumptions.
  const double directionNorm=CoronalCME::Norm(out->direction);
  const bool valid=(fitted||out->schema=="shock-front-ambient-v1.1")&&
      out->profile=="fixed-sse-pfss-parker-drag-v1"&&
      StrictUtcEpoch(out->referenceEpoch)&&
      out->timeScale=="UTC-resolved-to-elapsed-SI-seconds"&&
      out->coordinateFrame=="HCI"&&out->particleMode=="disabled"&&
      std::isfinite(directionNorm)&&directionNorm>0&&
      std::abs(directionNorm-1)<=1e-12&&out->halfWidthRad>0&&
      out->halfWidthRad<0.5*3.14159265358979323846&&
      out->physicalInnerRadiusM>=out->ambient.support.firstValidPlasmaRadiusM&&
      out->initialApexRadiusM>=out->physicalInnerRadiusM&&
      out->initialApexSpeedMPerS>=0&&out->finalPulseSpeedMPerS>0&&
      (fitted||out->accelerationDurationS>0)&&out->handoffApexRadiusM>out->initialApexRadiusM&&
      out->dragGammaPerM>=0&&out->effectiveTrajectoryWindMPerS>=0&&
      out->backgroundDtS>0&&out->polarCells>=2&&out->azimuthCells>=4&&
      out->surfaceTopology=="triangular-sse-cap-v1"&&
      out->weakMachTolerance>0&&out->rhResidualTolerance>0&&
      (fitted||out->endpointRadiusM>=out->handoffApexRadiusM)&&
      out->endpointRadiusM<=out->ambient.support.coverageRadiusM;
  if(!valid)return Return::Failure(Core::StatusCode::InvalidConfiguration,
      "shock-front event violates the frozen first-profile domain");
  if(out->historyModel==HistoryModel::ConstantSpeed&&
      std::abs(out->initialApexSpeedMPerS-out->finalPulseSpeedMPerS)>
      1e-12*std::max(1.0,out->finalPulseSpeedMPerS))return Return::Failure(
          Core::StatusCode::InvalidConfiguration,
          "constant-speed history requires identical initial/final speed");

  // The normalized manifest plus the independently verified asset digest is
  // the immutable physical identity carried by every epoch and restart.  The
  // sorted map order and 17-digit binary64 round-trip formatting prevent rank,
  // locale, and input-line ordering from changing that identity.
  std::ostringstream manifest;manifest<<std::setprecision(17);
  for(const auto& item:v)manifest<<item.first<<'='<<item.second<<'\n';
  manifest<<"assets.harmonics_content_sha256="<<out->harmonicAssetSha256<<'\n';
  out->normalizedManifest=manifest.str();
  out->physicsFingerprint=CoronalCME::Internal::Sha256Hex(out->normalizedManifest);
  out->ambient.physicsFingerprint=out->physicsFingerprint;
  return Return::Success(std::move(out));
}

} } } // namespace SEP::CoronaSwcme::ShockFront
