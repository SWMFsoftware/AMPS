#include "sep_corona_swcme/cme_event.h"

#include "sep_coronal_cme/configuration_parser.h"

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <iomanip>
#include <limits>
#include <set>
#include <sstream>
#include <utility>

namespace SEP { namespace CoronaSwcme { namespace {

using CoronalCME::EllipsoidKinematics;
using CoronalCME::HermiteKnot;
using CoronalCME::KinematicValue;
using EventResult=Core::Result<std::shared_ptr<const EventConfiguration>>;
constexpr double kAstronomicalUnitM=149597870700.0;

template<class T> Core::Result<T> Bad(const std::string& message) {
  return Core::Result<T>::Failure(Core::StatusCode::InvalidConfiguration,message);
}

std::string Trim(const std::string& input) {
  std::size_t begin=0,end=input.size();
  while(begin<end&&std::isspace(static_cast<unsigned char>(input[begin])))++begin;
  while(end>begin&&std::isspace(static_cast<unsigned char>(input[end-1])))--end;
  return input.substr(begin,end-begin);
}

bool Number(const std::string& text,double* value) {
  std::istringstream stream(text);
  std::string trailing;
  return static_cast<bool>(stream>>*value)&&!(stream>>trailing)&&std::isfinite(*value);
}

bool Integer(const std::string& text,int* value) {
  double parsed=0;
  if(!Number(text,&parsed)||parsed<0||parsed>10000000||std::floor(parsed)!=parsed)return false;
  *value=static_cast<int>(parsed);
  return true;
}

Core::Result<std::map<std::string,std::string>> ParseAssignments(
    const std::string& bytes) {
  std::map<std::string,std::string> values;
  std::istringstream input(bytes);
  std::string line;
  std::size_t lineNumber=0;
  while(std::getline(input,line)) {
    ++lineNumber;
    line=Trim(line);
    if(line.empty()||line[0]=='#')continue;
    const std::size_t equal=line.find('=');
    if(equal==std::string::npos||line.find('=',equal+1)!=std::string::npos)
      return Bad<std::map<std::string,std::string>>(
          "line "+std::to_string(lineNumber)+": expected one key=value assignment");
    const std::string key=Trim(line.substr(0,equal));
    const std::string value=Trim(line.substr(equal+1));
    if(key.empty()||value.empty())return Bad<std::map<std::string,std::string>>(
        "line "+std::to_string(lineNumber)+": empty key or value");
    if(!values.emplace(key,value).second)return Bad<std::map<std::string,std::string>>(
        "line "+std::to_string(lineNumber)+": duplicate key "+key);
  }
  return Core::Result<std::map<std::string,std::string>>::Success(std::move(values));
}

bool LowerHexSha256(const std::string& value) {
  if(value.size()!=64)return false;
  return std::all_of(value.begin(),value.end(),[](char c) {
    return (c>='0'&&c<='9')||(c>='a'&&c<='f');
  });
}

Core::Result<std::array<ComponentHistory,4>> ParseHistory(
    const std::string& bytes) {
  using Return=Core::Result<std::array<ComponentHistory,4>>;
  std::istringstream input(bytes);
  std::string line;
  bool header=false;
  std::array<ComponentHistory,4> result;
  std::size_t lineNumber=0;
  while(std::getline(input,line)) {
    ++lineNumber;
    line=Trim(line);
    if(line.empty()||line[0]=='#')continue;
    if(!header) {
      if(line!="time_s,center_m,center_rate_m_s,radial_m,radial_rate_m_s,"
               "lateral1_m,lateral1_rate_m_s,lateral2_m,lateral2_rate_m_s")
        return Return::Failure(Core::StatusCode::InvalidConfiguration,
            "history line "+std::to_string(lineNumber)+": unexpected CSV header");
      header=true;
      continue;
    }
    std::replace(line.begin(),line.end(),',',' ');
    std::istringstream row(line);
    std::array<double,9> x{};
    for(double& value:x)if(!(row>>value)||!std::isfinite(value))return Return::Failure(
        Core::StatusCode::InvalidConfiguration,
        "history line "+std::to_string(lineNumber)+": expected nine finite values");
    std::string trailing;
    if(row>>trailing)return Return::Failure(Core::StatusCode::InvalidConfiguration,
        "history line "+std::to_string(lineNumber)+": trailing field");
    for(std::size_t component=0;component<4;++component) {
      result[component].knots.push_back({x[0],x[1+2*component],x[2+2*component]});
    }
  }
  if(!header||result[0].knots.size()<2)return Return::Failure(
      Core::StatusCode::InvalidConfiguration,"history needs a header and at least two knots");
  for(std::size_t i=1;i<result[0].knots.size();++i)
    if(!(result[0].knots[i].timeS>result[0].knots[i-1].timeS))return Return::Failure(
        Core::StatusCode::InvalidConfiguration,"history times must be strictly increasing");
  return Return::Success(std::move(result));
}

Core::Result<AmbientInput> ParseAmbientAsset(const std::string& bytes) {
  const auto parsed=ParseAssignments(bytes);
  if(!parsed.ok())return Core::Result<AmbientInput>::Failure(
      parsed.status.code,"ambient asset: "+parsed.status.message);
  const auto& v=parsed.value;
  const std::set<std::string> required={"profile","source_surface_radius_m",
      "reference_radius_m","electron_density_at_reference_m3",
      "closed_base_electron_density_m3","rotation_rate_rad_per_s",
      "minimum_magnetic_field_t","trace_step_m","gradient_relative_step",
      "trace_maximum_steps","wind_table_points","harmonics"};
  for(const auto& key:required)if(!v.count(key))return Bad<AmbientInput>(
      "ambient asset missing key: "+key);
  for(const auto& item:v)if(!required.count(item.first))return Bad<AmbientInput>(
      "ambient asset unknown key: "+item.first);
  if(v.at("profile")!="pfss-parker-isothermal-v1")return Bad<AmbientInput>(
      "ambient asset profile does not match the selected provider");
  AmbientInput out;
  const std::vector<std::pair<std::string,double*>> numbers={
      {"source_surface_radius_m",&out.sourceSurfaceRadiusM},
      {"reference_radius_m",&out.referenceRadiusM},
      {"electron_density_at_reference_m3",&out.electronDensityAtReferenceM3},
      {"closed_base_electron_density_m3",&out.closedBaseElectronDensityM3},
      {"rotation_rate_rad_per_s",&out.rotationRateRadPerS},
      {"minimum_magnetic_field_t",&out.minimumMagneticFieldT},
      {"trace_step_m",&out.traceStepM},
      {"gradient_relative_step",&out.gradientRelativeStep}};
  for(const auto& item:numbers)if(!Number(v.at(item.first),item.second))
    return Bad<AmbientInput>("ambient asset invalid number: "+item.first);
  if(!Integer(v.at("trace_maximum_steps"),&out.traceMaximumSteps)||
      !Integer(v.at("wind_table_points"),&out.windTablePoints))
    return Bad<AmbientInput>("ambient asset has invalid integer controls");
  out.harmonics=v.at("harmonics");
  if(out.harmonics.empty())return Bad<AmbientInput>("ambient harmonics are empty");
  return Core::Result<AmbientInput>::Success(std::move(out));
}

Core::Result<RegionalInput> ParseRegionalAsset(const std::string& bytes) {
  const auto parsed=ParseAssignments(bytes);
  if(!parsed.ok())return Core::Result<RegionalInput>::Failure(
      parsed.status.code,"ejecta asset: "+parsed.status.message);
  const auto& v=parsed.value;
  const std::set<std::string> required={"profile","vector_potential_model",
      "sheath_admission_start_s","contact_apex_fraction",
      "ejecta_reference_density_kg_m3","ejecta_reference_pressure_pa",
      "axial_flux_wb","poloidal_flux_wb","minimum_jacobian",
      "maximum_integrated_force_ratio","maximum_local_force_ratio_p99",
      "maximum_force_work_ratio","added_heating"};
  for(const auto& key:required)if(!v.count(key))return Bad<RegionalInput>(
      "ejecta asset missing key: "+key);
  for(const auto& item:v)if(!required.count(item.first))return Bad<RegionalInput>(
      "ejecta asset unknown key: "+item.first);
  if(v.at("profile")!="vector-potential-material-map-v1")return Bad<RegionalInput>(
      "ejecta asset profile does not match the selected closure");
  RegionalInput out;
  out.vectorPotentialModel=v.at("vector_potential_model");
  out.addedHeating=v.at("added_heating");
  const std::vector<std::pair<std::string,double*>> numbers={
      {"sheath_admission_start_s",&out.sheathAdmissionStartS},
      {"contact_apex_fraction",&out.contactApexFraction},
      {"ejecta_reference_density_kg_m3",&out.ejectaReferenceDensityKgM3},
      {"ejecta_reference_pressure_pa",&out.ejectaReferencePressurePa},
      {"axial_flux_wb",&out.axialFluxWb},{"poloidal_flux_wb",&out.poloidalFluxWb},
      {"minimum_jacobian",&out.minimumJacobian},
      {"maximum_integrated_force_ratio",&out.maximumIntegratedForceRatio},
      {"maximum_local_force_ratio_p99",&out.maximumLocalForceRatioP99},
      {"maximum_force_work_ratio",&out.maximumForceWorkRatio}};
  for(const auto& item:numbers)if(!Number(v.at(item.first),item.second))
    return Bad<RegionalInput>("ejecta asset invalid number: "+item.first);
  return Core::Result<RegionalInput>::Success(std::move(out));
}

struct Cubic { double c0=0,c1=0,c2=0,c3=0,dt=0; };

Cubic Polynomial(const HermiteKnot& a,const HermiteKnot& b) {
  const double dt=b.timeS-a.timeS;
  return {a.value,dt*a.ratePerS,
      -3*a.value-2*dt*a.ratePerS+3*b.value-dt*b.ratePerS,
      2*a.value+dt*a.ratePerS-2*b.value+dt*b.ratePerS,dt};
}

double Value(const Cubic& p,double s) {
  return ((p.c3*s+p.c2)*s+p.c1)*s+p.c0;
}

double First(const Cubic& p,double s) {
  return (p.c1+2*p.c2*s+3*p.c3*s*s)/p.dt;
}

void AddRoots(double a,double b,double c,std::vector<double>* points) {
  const double scale=std::max({std::abs(a),std::abs(b),std::abs(c),1.0});
  if(std::abs(a)<=32*std::numeric_limits<double>::epsilon()*scale) {
    if(std::abs(b)>32*std::numeric_limits<double>::epsilon()*scale) {
      const double root=-c/b;
      if(root>0&&root<1)points->push_back(root);
    }
    return;
  }
  const double discriminant=b*b-4*a*c;
  if(discriminant<0)return;
  const double root=std::sqrt(std::max(0.0,discriminant));
  for(double s:{(-b-root)/(2*a),(-b+root)/(2*a)})
    if(s>0&&s<1)points->push_back(s);
}

std::pair<double,double> ValueRange(const Cubic& p) {
  std::vector<double> points={0,1};
  AddRoots(3*p.c3,2*p.c2,p.c1,&points);
  double low=std::numeric_limits<double>::infinity(),high=-low;
  for(double s:points){low=std::min(low,Value(p,s));high=std::max(high,Value(p,s));}
  return {low,high};
}

std::pair<double,double> RateRange(const Cubic& p) {
  std::vector<double> points={0,1};
  if(p.c3!=0) {
    const double s=-p.c2/(3*p.c3);
    if(s>0&&s<1)points.push_back(s);
  }
  double low=std::numeric_limits<double>::infinity(),high=-low;
  for(double s:points){low=std::min(low,First(p,s));high=std::max(high,First(p,s));}
  return {low,high};
}

KinematicValue Apex(const EllipsoidKinematics& k) {
  return {k.centerDistanceM.value+k.radialSemiAxisM.value,
      k.centerDistanceM.firstDerivative+k.radialSemiAxisM.firstDerivative,
      k.centerDistanceM.secondDerivative+k.radialSemiAxisM.secondDerivative};
}

Core::Result<EllipsoidKinematics> HistoryState(
    const EventConfiguration& c,double timeS) {
  std::array<KinematicValue,4> q{};
  for(std::size_t i=0;i<q.size();++i) {
    const auto value=c.components[i].At(timeS);
    if(!value.ok())return Core::Result<EllipsoidKinematics>::Failure(
        value.status.code,value.status.message);
    q[i]=value.value;
  }
  return Core::Result<EllipsoidKinematics>::Success({q[0],q[1],q[2],q[3]});
}

Core::Result<KinematicValue> Dbm(const HandoffLaw& law,
    const KinematicValue& crossing,double elapsedS) {
  using Return=Core::Result<KinematicValue>;
  if(!std::isfinite(elapsedS)||elapsedS<0||!(crossing.value>0)||
      !(crossing.firstDerivative>=0)||!(law.ambientSpeedMPerS>0)||
      !(law.dragCoefficientPerM>=0))return Return::Failure(
          Core::StatusCode::InvalidConfiguration,"invalid DBM crossing state");
  const double difference=crossing.firstDerivative-law.ambientSpeedMPerS;
  const double magnitude=std::abs(difference);
  if(law.dragCoefficientPerM==0||magnitude==0)return Return::Success({
      std::fma(crossing.firstDerivative,elapsedS,crossing.value),
      crossing.firstDerivative,0});
  const double x=law.dragCoefficientPerM*magnitude*elapsedS;
  if(!std::isfinite(x)||!(1+x>0))return Return::Failure(
      Core::StatusCode::OutOfDomain,"DBM continuation left its finite branch");
  // Sign-aware drag-based motion solves dv/dt=-gamma(v-w)|v-w| exactly.
  // Keeping the sign of v-w matters for a CME initially slower than the wind;
  // the common unsigned logarithmic formula would accelerate it the wrong way.
  // Values are heliocentric metres, seconds and m/s in the inertial HCI frame.
  const double velocityDifference=difference/(1+x);
  const double speed=law.ambientSpeedMPerS+velocityDifference;
  const double dragDistance=std::abs(x)<1e-8?
      std::copysign(magnitude*elapsedS*(1-.5*x+x*x/3-x*x*x/4+x*x*x*x/5),difference):
      std::copysign(std::log1p(x)/law.dragCoefficientPerM,difference);
  const double radius=crossing.value+law.ambientSpeedMPerS*elapsedS+dragDistance;
  const double acceleration=-law.dragCoefficientPerM*velocityDifference*
      std::abs(velocityDifference);
  if(!(std::isfinite(radius)&&std::isfinite(speed)&&std::isfinite(acceleration)&&
      radius>0&&speed>=0))return Return::Failure(Core::StatusCode::OutOfDomain,
          "DBM continuation produced a nonphysical state");
  return Return::Success({radius,speed,acceleration});
}

KinematicValue Scale(const KinematicValue& scale,double reference) {
  return {reference*scale.value,reference*scale.firstDerivative,
      reference*scale.secondDerivative};
}

KinematicValue Blend(const KinematicValue& a,const KinematicValue& b,
    double w,double dw,double ddw) {
  // This differentiates one position blend.  The dw*(b-a) and ddw terms are
  // physical transition velocity/acceleration contributions; blending already
  // computed velocities would omit them and create a radius/speed reset.
  return {(1-w)*a.value+w*b.value,
      (1-w)*a.firstDerivative+w*b.firstDerivative+dw*(b.value-a.value),
      (1-w)*a.secondDerivative+w*b.secondDerivative+
          2*dw*(b.firstDerivative-a.firstDerivative)+ddw*(b.value-a.value)};
}

Core::Result<AssetIdentity> Acquire(const std::string& role,
    const std::map<std::string,std::string>& values,const AssetReader& reader,
    std::string* bytes) {
  const std::string prefix="assets."+role;
  const std::string path=values.at(prefix+"_file");
  const std::string checksum=values.at(prefix+"_sha256");
  if(!LowerHexSha256(checksum))return Bad<AssetIdentity>(
      prefix+"_sha256 must be a lower-case SHA-256 digest");
  const auto acquired=reader(path);
  if(!acquired.ok())return Core::Result<AssetIdentity>::Failure(
      acquired.status.code,"cannot acquire "+role+" asset: "+acquired.status.message);
  const std::string actual=CoronalCME::ComputeContentChecksum(acquired.value);
  if(actual!=checksum)return Core::Result<AssetIdentity>::Failure(
      Core::StatusCode::DataIntegrityFailure,role+" asset SHA-256 mismatch");
  *bytes=acquired.value;
  return Core::Result<AssetIdentity>::Success({role,path,actual,bytes->size()});
}

} // namespace

Core::Result<KinematicValue> ComponentHistory::At(double timeS) const {
  return CoronalCME::EvaluateCubicHermiteHistory(knots,timeS);
}

const char* Name(EvolutionPhase phase) noexcept {
  switch(phase) {
    case EvolutionPhase::CoronalHistory:return "coronal-history";
    case EvolutionPhase::HandoffTransition:return "handoff-transition";
    case EvolutionPhase::SwcmeOuter:return "swcme-outer";
  }
  return "unknown";
}

Core::Result<RadialExtent> EvaluateRadialExtent(
    const EllipsoidKinematics& k,double solarRadiusM) {
  using Return=Core::Result<RadialExtent>;
  const double d=k.centerDistanceM.value,a=k.radialSemiAxisM.value;
  const double b=k.firstLateralSemiAxisM.value,z=k.secondLateralSemiAxisM.value;
  if(!(std::isfinite(d)&&std::isfinite(a)&&std::isfinite(b)&&std::isfinite(z)&&
      std::isfinite(solarRadiusM)&&d>0&&a>0&&b>0&&z>0&&solarRadiusM>0))
    return Return::Failure(Core::StatusCode::InvalidConfiguration,
        "radial extent requires positive finite geometry and solar radius");
  const auto extreme=[&](double lateralSquared,bool minimum) {
    double value=minimum?std::min((d-a)*(d-a),(d+a)*(d+a)):
        std::max((d-a)*(d-a),(d+a)*(d+a));
    if(a*a!=lateralSquared) {
      const double q=-d*a/(a*a-lateralSquared);
      if(q>=-1&&q<=1) {
        const double candidate=d*d+lateralSquared+2*d*a*q+
            (a*a-lateralSquared)*q*q;
        value=minimum?std::min(value,candidate):std::max(value,candidate);
      }
    }
    return std::sqrt(std::max(0.0,value));
  };
  RadialExtent out;
  out.minimumRadiusM=extreme(std::min(b*b,z*z),true);
  out.maximumRadiusM=extreme(std::max(b*b,z*z),false);
  out.intersectsSolarSurface=out.minimumRadiusM<=solarRadiusM&&
      out.maximumRadiusM>=solarRadiusM;
  return Return::Success(out);
}

Core::Result<EventKinematics> EventConfiguration::At(double timeS) const {
  using Return=Core::Result<EventKinematics>;
  if(!std::isfinite(timeS)||timeS<support.startS||timeS>support.endS)
    return Return::Failure(Core::StatusCode::OutOfDomain,
        "event query is outside declared time coverage");
  const auto reference=HistoryState(*this,handoff.transitionBeginS);
  if(!reference.ok())return Return::Failure(reference.status.code,reference.status.message);
  const auto crossing=Apex(reference.value);
  if(timeS<=handoff.transitionBeginS) {
    const auto state=HistoryState(*this,timeS);
    if(!state.ok())return Return::Failure(state.status.code,state.status.message);
    return Return::Success({timeS,0,EvolutionPhase::CoronalHistory,
        state.value,Apex(state.value)});
  }
  const auto apex=Dbm(handoff,crossing,timeS-handoff.transitionBeginS);
  if(!apex.ok())return Return::Failure(apex.status.code,apex.status.message);
  // Outer continuation is self-similar about the Sun: the same scalar maps
  // center and all three semiaxes.  This preserves the full triaxial shape and
  // its material parameterization; matching only the apex would be insufficient.
  const KinematicValue scale={apex.value.value/crossing.value,
      apex.value.firstDerivative/crossing.value,
      apex.value.secondDerivative/crossing.value};
  EllipsoidKinematics outer={Scale(scale,reference.value.centerDistanceM.value),
      Scale(scale,reference.value.radialSemiAxisM.value),
      Scale(scale,reference.value.firstLateralSemiAxisM.value),
      Scale(scale,reference.value.secondLateralSemiAxisM.value)};
  if(timeS>=handoff.transitionEndS)return Return::Success({timeS,1,
      EvolutionPhase::SwcmeOuter,outer,apex.value});
  const auto early=HistoryState(*this,timeS);
  if(!early.ok())return Return::Failure(early.status.code,early.status.message);
  const double duration=handoff.transitionEndS-handoff.transitionBeginS;
  const double u=(timeS-handoff.transitionBeginS)/duration;
  const double u2=u*u,u3=u2*u,u4=u3*u,u5=u4*u;
  const double w=10*u3-15*u4+6*u5;
  const double dw=(30*u2-60*u3+30*u4)/duration;
  const double ddw=(60*u-180*u2+120*u3)/(duration*duration);
  EllipsoidKinematics state={
      Blend(early.value.centerDistanceM,outer.centerDistanceM,w,dw,ddw),
      Blend(early.value.radialSemiAxisM,outer.radialSemiAxisM,w,dw,ddw),
      Blend(early.value.firstLateralSemiAxisM,outer.firstLateralSemiAxisM,w,dw,ddw),
      Blend(early.value.secondLateralSemiAxisM,outer.secondLateralSemiAxisM,w,dw,ddw)};
  return Return::Success({timeS,w,EvolutionPhase::HandoffTransition,state,Apex(state)});
}

Core::Status ValidateEventConfiguration(const EventConfiguration& c) {
  const auto fail=[](const std::string& message) {
    return Core::Status::Failure(Core::StatusCode::InvalidConfiguration,message);
  };
  if(c.profile!="pfss-parker-shock-fed-map-v1"||c.coordinateFrame!="inertial-hci"||
      c.ambientProfile!="pfss-parker-isothermal-v1"||
      c.sheathModel!="rh-ballistic-material-map-v1"||
      c.ejectaModel!="vector-potential-material-map-v1"||
      c.outerEvolution!="swcme-dbm-constant-wind-v1"||
      c.attachmentPolicy!="attached-then-detached")
    return fail("unsupported G1 profile or component selector");
  const auto& s=c.support;
  if(!(std::isfinite(s.startS)&&std::isfinite(s.endS)&&s.endS>s.startS&&
      s.solarRadiusM>0&&s.firstValidPlasmaRadiusM>=s.solarRadiusM&&
      s.coverageRadiusM>=1.05*kAstronomicalUnitM&&s.rootToleranceS>0&&
      s.rootToleranceS<(s.endS-s.startS)/10&&s.maximumApexAccelerationMPerS2>0))
    return fail("invalid time, radial support, root tolerance, or acceleration bound");
  if(!(c.handoff.transitionBeginS>s.startS&&
      c.handoff.transitionEndS>c.handoff.transitionBeginS&&
      c.handoff.transitionEndS<s.endS&&c.handoff.ambientSpeedMPerS>0&&
      c.handoff.dragCoefficientPerM>=0))return fail("invalid handoff interval or DBM law");
  const auto& p=c.composition;
  if(!(p.gammaAdiabatic>1&&p.electronTemperatureK>0&&p.protonTemperatureK>0&&
      p.alphaTemperatureK>0&&p.alphaToProtonNumberRatio>=0))
    return fail("invalid composition, temperatures, or adiabatic index");
  const auto& a=c.ambient;
  if(!(a.sourceSurfaceRadiusM>s.solarRadiusM&&
      a.sourceSurfaceRadiusM<a.referenceRadiusM&&
      a.referenceRadiusM<=s.coverageRadiusM&&
      a.electronDensityAtReferenceM3>0&&a.closedBaseElectronDensityM3>0&&
      a.rotationRateRadPerS>=0&&a.minimumMagneticFieldT>0&&a.traceStepM>0&&
      a.traceStepM<=0.05*s.solarRadiusM&&a.gradientRelativeStep>=1e-7&&
      a.gradientRelativeStep<=1e-3&&a.traceMaximumSteps>=64&&
      a.windTablePoints>=257&&a.harmonics.size()>0))
    return fail("ambient asset is outside the selected provider support");
  const auto& r=c.regional;
  if(!(r.sheathAdmissionStartS<=s.startS&&r.contactApexFraction>0&&
      r.contactApexFraction<1&&r.ejectaReferenceDensityKgM3>0&&
      r.ejectaReferencePressurePa>0&&std::isfinite(r.axialFluxWb)&&
      std::isfinite(r.poloidalFluxWb)&&
      (std::abs(r.axialFluxWb)+std::abs(r.poloidalFluxWb)>0)&&
      r.minimumJacobian>0&&r.maximumIntegratedForceRatio==2&&
      r.maximumLocalForceRatioP99==5&&r.maximumForceWorkRatio==2&&
      r.vectorPotentialModel=="axisymmetric-polynomial-a-v1"&&
      r.addedHeating=="zero"))return fail(
          "regional asset is incomplete or differs from the frozen closure budgets");
  for(const auto& h:c.components) {
    if(h.knots.size()!=c.components[0].knots.size()||h.knots.size()<2||
        h.knots.front().timeS>s.startS||
        h.knots.back().timeS<c.handoff.transitionEndS)
      return fail("component history does not cover start through handoff transition");
    for(std::size_t i=1;i<h.knots.size();++i) {
      if(h.knots[i].timeS!=c.components[0].knots[i].timeS||
          !(h.knots[i].timeS>h.knots[i-1].timeS))
        return fail("component histories do not share strictly increasing knots");
      const auto range=ValueRange(Polynomial(h.knots[i-1],h.knots[i]));
      if(!(range.first>0))return fail("geometry component is nonpositive between knots");
    }
  }
  // Exact polynomial extrema certify the pre-handoff apex value and rate; the
  // transition/DBM branch is analytic but non-polynomial and is independently
  // checked by nested deterministic scans below.
  for(std::size_t i=1;i<c.components[0].knots.size();++i) {
    Cubic apex=Polynomial(c.components[0].knots[i-1],c.components[0].knots[i]);
    const Cubic radial=Polynomial(c.components[1].knots[i-1],c.components[1].knots[i]);
    apex.c0+=radial.c0;apex.c1+=radial.c1;apex.c2+=radial.c2;apex.c3+=radial.c3;
    const auto value=ValueRange(apex),rate=RateRange(apex);
    if(value.first<s.firstValidPlasmaRadiusM||value.second>s.coverageRadiusM)
      return fail("coronal apex leaves declared radial plasma coverage");
    if(rate.first<0)return fail("coronal apex reverses between history knots");
  }
  double maxima[2]={0,0};
  for(int level=0;level<2;++level) {
    const int intervals=4096*(level+1);
    for(int i=0;i<=intervals;++i) {
      const double time=s.startS+(s.endS-s.startS)*i/intervals;
      const auto state=c.At(time);
      if(!state.ok())return state.status;
      const auto& k=state.value.ellipsoid;
      const double dimensions[]={k.centerDistanceM.value,k.radialSemiAxisM.value,
          k.firstLateralSemiAxisM.value,k.secondLateralSemiAxisM.value};
      for(double dimension:dimensions)if(!(std::isfinite(dimension)&&dimension>0))
        return fail("continuous handoff contains nonpositive geometry");
      if(!(state.value.apexRadiusM.value>=s.firstValidPlasmaRadiusM&&
          state.value.apexRadiusM.value<=s.coverageRadiusM&&
          state.value.apexRadiusM.firstDerivative>=0))
        return fail("continuous handoff leaves coverage or reverses");
      maxima[level]=std::max(maxima[level],
          std::abs(state.value.apexRadiusM.secondDerivative));
    }
  }
  const double unresolved=4*std::abs(maxima[1]-maxima[0])+256*
      std::numeric_limits<double>::epsilon()*std::max(1.0,maxima[1]);
  if(maxima[1]+unresolved>s.maximumApexAccelerationMPerS2*(1+1e-12))
    return fail("continuous handoff exceeds the frozen apex acceleration bound");
  if(c.assets.size()!=3)return fail("history, ambient, and ejecta assets are mandatory");
  std::set<std::string> roles;
  for(const auto& asset:c.assets)
    if(!roles.insert(asset.role).second||asset.path.empty()||asset.bytes==0||
        !LowerHexSha256(asset.sha256))return fail("invalid or duplicate asset identity");
  if(roles!=std::set<std::string>({"ambient","ejecta_reference","history"}))
    return fail("incomplete asset role set");
  if(c.physicsFingerprint.empty()||c.manifest.empty())return fail("missing event identity");
  const auto initial=c.At(s.startS),final=c.At(s.endS);
  if(!initial.ok()||!final.ok())return fail("attachment endpoints are outside event coverage");
  const auto initialExtent=EvaluateRadialExtent(initial.value.ellipsoid,s.solarRadiusM);
  const auto finalExtent=EvaluateRadialExtent(final.value.ellipsoid,s.solarRadiusM);
  if(!initialExtent.ok()||!finalExtent.ok()||
      !initialExtent.value.intersectsSolarSurface||
      finalExtent.value.minimumRadiusM<=s.solarRadiusM)
    return fail("attached-then-detached policy does not match exact radial extrema");
  return Core::Status::Success();
}

Core::Result<std::shared_ptr<const EventConfiguration>> ResolveEventConfiguration(
    const std::string& inputBytes,const AssetReader& reader) {
  if(!reader)return EventResult::Failure(Core::StatusCode::InvalidConfiguration,
      "event configuration requires an asset reader");
  const auto parsed=ParseAssignments(inputBytes);
  if(!parsed.ok())return EventResult::Failure(parsed.status.code,parsed.status.message);
  auto values=parsed.value;
  const std::set<std::string> required={
    "run.profile","run.coordinate_frame","run.start_s","run.end_s",
    "domain.solar_radius_m","domain.first_valid_plasma_radius_m","domain.coverage_radius_m",
    "geometry.latitude_rad","geometry.longitude_rad","geometry.lateral_tilt_rad",
    "geometry.attachment_policy",
    "assets.history_file","assets.history_sha256","assets.ambient_file","assets.ambient_sha256",
    "assets.ejecta_reference_file","assets.ejecta_reference_sha256",
    "ambient.profile","plasma.gamma_adiabatic","plasma.electron_temperature_k",
    "plasma.proton_temperature_k","plasma.alpha_temperature_k",
    "plasma.alpha_to_proton_number_ratio","plasma.include_electron_mass",
    "regions.sheath_model","regions.ejecta_model","handoff.evolution",
    "handoff.transition_begin_s","handoff.transition_end_s","handoff.ambient_speed_m_s",
    "handoff.drag_coefficient_per_m","numerics.root_tolerance_s",
    "numerics.maximum_apex_acceleration_m_s2"};
  for(const auto& key:required)if(!values.count(key))return EventResult::Failure(
      Core::StatusCode::InvalidConfiguration,"missing event key: "+key);
  for(const auto& item:values)if(!required.count(item.first))return EventResult::Failure(
      Core::StatusCode::InvalidConfiguration,"unknown event key: "+item.first);

  std::shared_ptr<EventConfiguration> event(new EventConfiguration);
  event->profile=values["run.profile"];
  event->coordinateFrame=values["run.coordinate_frame"];
  event->ambientProfile=values["ambient.profile"];
  event->sheathModel=values["regions.sheath_model"];
  event->ejectaModel=values["regions.ejecta_model"];
  event->outerEvolution=values["handoff.evolution"];
  event->attachmentPolicy=values["geometry.attachment_policy"];
  const std::vector<std::pair<std::string,double*>> numbers={
    {"run.start_s",&event->support.startS},{"run.end_s",&event->support.endS},
    {"domain.solar_radius_m",&event->support.solarRadiusM},
    {"domain.first_valid_plasma_radius_m",&event->support.firstValidPlasmaRadiusM},
    {"domain.coverage_radius_m",&event->support.coverageRadiusM},
    {"plasma.gamma_adiabatic",&event->composition.gammaAdiabatic},
    {"plasma.electron_temperature_k",&event->composition.electronTemperatureK},
    {"plasma.proton_temperature_k",&event->composition.protonTemperatureK},
    {"plasma.alpha_temperature_k",&event->composition.alphaTemperatureK},
    {"plasma.alpha_to_proton_number_ratio",&event->composition.alphaToProtonNumberRatio},
    {"handoff.transition_begin_s",&event->handoff.transitionBeginS},
    {"handoff.transition_end_s",&event->handoff.transitionEndS},
    {"handoff.ambient_speed_m_s",&event->handoff.ambientSpeedMPerS},
    {"handoff.drag_coefficient_per_m",&event->handoff.dragCoefficientPerM},
    {"numerics.root_tolerance_s",&event->support.rootToleranceS},
    {"numerics.maximum_apex_acceleration_m_s2",&event->support.maximumApexAccelerationMPerS2}};
  for(const auto& item:numbers)if(!Number(values[item.first],item.second))return EventResult::Failure(
      Core::StatusCode::InvalidConfiguration,"invalid finite number: "+item.first);
  double latitude=0,longitude=0,tilt=0;
  if(!Number(values["geometry.latitude_rad"],&latitude)||
      !Number(values["geometry.longitude_rad"],&longitude)||
      !Number(values["geometry.lateral_tilt_rad"],&tilt))return EventResult::Failure(
          Core::StatusCode::InvalidConfiguration,"invalid geometry orientation");
  const auto basis=CoronalCME::BuildRadialPrincipalBasis(latitude,longitude,tilt);
  if(!basis.ok())return EventResult::Failure(basis.status.code,basis.status.message);
  event->basis=basis.value;
  if(values["plasma.include_electron_mass"]!="true"&&
      values["plasma.include_electron_mass"]!="false")return EventResult::Failure(
          Core::StatusCode::InvalidConfiguration,
          "plasma.include_electron_mass requires true or false");
  event->composition.includeElectronMass=values["plasma.include_electron_mass"]=="true";

  std::string historyBytes,ambientBytes,ejectaBytes;
  for(const auto& request:std::vector<std::pair<std::string,std::string*>>{
      {"history",&historyBytes},{"ambient",&ambientBytes},
      {"ejecta_reference",&ejectaBytes}}) {
    const auto asset=Acquire(request.first,values,reader,request.second);
    if(!asset.ok())return EventResult::Failure(asset.status.code,asset.status.message);
    event->assets.push_back(asset.value);
    values["assets."+request.first+"_content_checksum"]=asset.value.sha256;
  }
  const auto history=ParseHistory(historyBytes);
  if(!history.ok())return EventResult::Failure(history.status.code,history.status.message);
  event->components=history.value;
  const auto ambient=ParseAmbientAsset(ambientBytes);
  if(!ambient.ok())return EventResult::Failure(ambient.status.code,ambient.status.message);
  event->ambient=ambient.value;
  const auto regional=ParseRegionalAsset(ejectaBytes);
  if(!regional.ok())return EventResult::Failure(regional.status.code,regional.status.message);
  event->regional=regional.value;
  event->normalizedAssignments=values;
  event->physicsFingerprint=CoronalCME::ComputePhysicsFingerprint(values);
  std::ostringstream manifest;
  manifest<<std::setprecision(17)<<"sep-corona-swcme-event-v1\n";
  for(const auto& item:values) {
    if(item.first.size()>=5&&item.first.compare(item.first.size()-5,5,"_file")==0)continue;
    manifest<<item.first<<'='<<item.second<<'\n';
  }
  event->manifest=manifest.str();
  const auto valid=ValidateEventConfiguration(*event);
  if(!valid.ok())return EventResult::Failure(valid.code,valid.message);
  return EventResult::Success(std::move(event));
}

} } // namespace SEP::CoronaSwcme
