#include "sep_corona_swcme/cme_event.h"
#include "sep_coronal_cme/configuration_parser.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>

namespace {

using SEP::Core::Result;
using SEP::CoronaSwcme::EventConfiguration;

void Require(bool condition,const std::string& message) {
  if(!condition)throw std::runtime_error(message);
}

std::string History(bool overshoot=false) {
  std::ostringstream out;
  out<<"time_s,center_m,center_rate_m_s,radial_m,radial_rate_m_s,"
        "lateral1_m,lateral1_rate_m_s,lateral2_m,lateral2_rate_m_s\n";
  if(overshoot) {
    // Positive endpoints with a deliberately negative cubic interior.
    out<<"0,800000000,700000,200000000,-10000000,300000000,200000,250000000,180000\n"
       <<"200,940000000,700000,260000000,10000000,340000000,200000,286000000,180000\n";
  } else {
    out<<"0,800000000,700000,200000000,300000,300000000,200000,250000000,180000\n"
       <<"200,940000000,700000,260000000,300000,340000000,200000,286000000,180000\n";
  }
  return out.str();
}

std::string AmbientAsset() {
  return "profile=pfss-parker-isothermal-v1\n"
      "source_surface_radius_m=1739250000\n"
      "reference_radius_m=149597870700\n"
      "electron_density_at_reference_m3=5000000\n"
      "closed_base_electron_density_m3=1e14\n"
      "rotation_rate_rad_per_s=2.86533e-6\n"
      "minimum_magnetic_field_t=1e-12\n"
      "trace_step_m=10000000\n"
      "gradient_relative_step=1e-5\n"
      "trace_maximum_steps=10000\n"
      "wind_table_points=1025\n"
      "harmonics=1,0,0.0002,0\n";
}

std::string EjectaAsset() {
  return "profile=vector-potential-material-map-v1\n"
      "vector_potential_model=axisymmetric-polynomial-a-v1\n"
      "sheath_admission_start_s=0\n"
      "contact_apex_fraction=0.82\n"
      "ejecta_reference_density_kg_m3=1e-13\n"
      "ejecta_reference_pressure_pa=1e-3\n"
      "axial_flux_wb=1e13\n"
      "poloidal_flux_wb=5e12\n"
      "minimum_jacobian=1e-6\n"
      "maximum_integrated_force_ratio=2\n"
      "maximum_local_force_ratio_p99=5\n"
      "maximum_force_work_ratio=2\n"
      "added_heating=zero\n";
}

std::string Configuration(const std::map<std::string,std::string>& assets,
    const std::string& sheath="rh-ballistic-material-map-v1") {
  const auto checksum=[&](const std::string& key) {
    return SEP::CoronalCME::ComputeContentChecksum(assets.at(key));
  };
  std::ostringstream out;
  out<<"run.profile=pfss-parker-shock-fed-map-v1\n"
     <<"run.coordinate_frame=inertial-hci\n"
     <<"run.start_s=0\nrun.end_s=200000\n"
     <<"domain.solar_radius_m=695700000\n"
     <<"domain.first_valid_plasma_radius_m=695700000\n"
     <<"domain.coverage_radius_m=299195741400\n"
     <<"geometry.latitude_rad=0.1\ngeometry.longitude_rad=-0.2\n"
     <<"geometry.lateral_tilt_rad=0.3\n"
     <<"geometry.attachment_policy=attached-then-detached\n"
     <<"assets.history_file=history.csv\nassets.history_sha256="<<checksum("history.csv")<<"\n"
     <<"assets.ambient_file=ambient.asset\nassets.ambient_sha256="<<checksum("ambient.asset")<<"\n"
     <<"assets.ejecta_reference_file=ejecta.asset\n"
     <<"assets.ejecta_reference_sha256="<<checksum("ejecta.asset")<<"\n"
     <<"ambient.profile=pfss-parker-isothermal-v1\n"
     <<"plasma.gamma_adiabatic=1.6666666666666667\n"
     <<"plasma.electron_temperature_k=1000000\n"
     <<"plasma.proton_temperature_k=1000000\n"
     <<"plasma.alpha_temperature_k=1000000\n"
     <<"plasma.alpha_to_proton_number_ratio=0.04\n"
     <<"plasma.include_electron_mass=true\n"
     <<"regions.sheath_model="<<sheath<<"\n"
     <<"regions.ejecta_model=vector-potential-material-map-v1\n"
     <<"handoff.evolution=swcme-dbm-constant-wind-v1\n"
     <<"handoff.transition_begin_s=100\n"
     <<"handoff.transition_end_s=200\n"
     <<"handoff.ambient_speed_m_s=400000\n"
     <<"handoff.drag_coefficient_per_m=1e-11\n"
     <<"numerics.root_tolerance_s=0.001\n"
     <<"numerics.maximum_apex_acceleration_m_s2=100000\n";
  return out.str();
}

Result<std::string> Read(const std::map<std::string,std::string>& assets,
    const std::string& path) {
  const auto found=assets.find(path);
  if(found==assets.end())return Result<std::string>::Failure(
      SEP::Core::StatusCode::DataIntegrityFailure,"missing fixture asset");
  return Result<std::string>::Success(found->second);
}

std::shared_ptr<const EventConfiguration> Resolve(
    const std::map<std::string,std::string>& assets,const std::string& input) {
  const auto result=SEP::CoronaSwcme::ResolveEventConfiguration(input,
      [&](const std::string& path){return Read(assets,path);});
  Require(result.ok(),result.status.message);
  return result.value;
}

void TestValidAndIndependentHandoff() {
  const std::map<std::string,std::string> assets={{"history.csv",History()},
      {"ambient.asset",AmbientAsset()},{"ejecta.asset",EjectaAsset()}};
  const auto event=Resolve(assets,Configuration(assets));
  Require(event->assets.size()==3,"complete asset identity missing");
  Require(event->physicsFingerprint.size()==64,"event fingerprint is not SHA-256");
  const auto before=event->At(99),begin=event->At(100),inside=event->At(150),end=event->At(200),late=event->At(200000);
  Require(before.ok()&&begin.ok()&&inside.ok()&&end.ok()&&late.ok(),"covered event query failed");
  Require(begin.value.phase==SEP::CoronaSwcme::EvolutionPhase::CoronalHistory,
      "transition begin must retain the exact coronal endpoint");
  Require(inside.value.phase==SEP::CoronaSwcme::EvolutionPhase::HandoffTransition&&
      inside.value.handoffWeight>0&&inside.value.handoffWeight<1,
      "inside transition phase/weight is wrong");
  Require(end.value.phase==SEP::CoronaSwcme::EvolutionPhase::SwcmeOuter,
      "transition end did not select outer state");
  Require(late.value.apexRadiusM.value>149597870700.0,
      "full event support does not actually reach 1 AU");
  Require(!event->At(-1).ok()&&!event->At(200001).ok(),
      "unsupported history extrapolation was accepted");
  const auto initialExtent=SEP::CoronaSwcme::EvaluateRadialExtent(
      begin.value.ellipsoid,event->support.solarRadiusM);
  const auto finalExtent=SEP::CoronaSwcme::EvaluateRadialExtent(
      late.value.ellipsoid,event->support.solarRadiusM);
  Require(initialExtent.ok()&&initialExtent.value.intersectsSolarSurface&&
      finalExtent.ok()&&!finalExtent.value.intersectsSolarSurface&&
      finalExtent.value.minimumRadiusM>event->support.solarRadiusM,
      "exact attachment/detachment extrema are wrong");
  const auto shape=SEP::CoronalCME::FixedOrientationEllipsoid::FromCenter(
      event->basis,inside.value.ellipsoid);
  Require(shape.ok(),"independent extent fixture shape failed");
  double sampledMinimum=1e300,sampledMaximum=0;
  for(int i=0;i<=2000;++i)for(int j=0;j<32;++j) {
    const double polar=3.14159265358979323846*i/2000;
    const double azimuth=2*3.14159265358979323846*j/32;
    const auto point=shape.value.Point(polar,azimuth);
    const double radius=std::sqrt(point.x*point.x+point.y*point.y+point.z*point.z);
    sampledMinimum=std::min(sampledMinimum,radius);
    sampledMaximum=std::max(sampledMaximum,radius);
  }
  const auto exactExtent=SEP::CoronaSwcme::EvaluateRadialExtent(
      inside.value.ellipsoid,event->support.solarRadiusM);
  Require(std::abs(exactExtent.value.minimumRadiusM-sampledMinimum)<2e-6*sampledMinimum&&
      std::abs(exactExtent.value.maximumRadiusM-sampledMaximum)<2e-6*sampledMaximum,
      "exact radial extrema disagree with independent surface sampling");

  // Independent centered difference checks the derivative of the one blended
  // position curve. Omitting the w-dot term fails this test by orders of
  // magnitude for the deliberately non-matched DBM/history accelerations.
  const double h=1e-3;
  const auto left=event->At(150-h),right=event->At(150+h);
  const double derivative=(right.value.apexRadiusM.value-left.value.apexRadiusM.value)/(2*h);
  Require(std::abs(derivative-inside.value.apexRadiusM.firstDerivative)<=
      1e-6*std::max(1.0,std::abs(derivative)),"handoff velocity is not the derivative of position");
}

void TestStrictIdentityAndNegatives() {
  std::map<std::string,std::string> assets={{"history.csv",History()},
      {"ambient.asset",AmbientAsset()},{"ejecta.asset",EjectaAsset()}};
  const std::string valid=Configuration(assets);
  const auto reader=[&](const std::string& path){return Read(assets,path);};

  std::string corrupt=valid;
  const std::size_t checksum=corrupt.find("assets.ambient_sha256=")+23;
  corrupt[checksum]=corrupt[checksum]=='0'?'1':'0';
  auto result=SEP::CoronaSwcme::ResolveEventConfiguration(corrupt,reader);
  Require(!result.ok()&&result.status.code==SEP::Core::StatusCode::DataIntegrityFailure,
      "corrupt asset checksum was not a data-integrity failure");

  result=SEP::CoronaSwcme::ResolveEventConfiguration(valid+"unknown.key=1\n",reader);
  Require(!result.ok(),"unknown event key accepted");
  result=SEP::CoronaSwcme::ResolveEventConfiguration(
      Configuration(assets,"legacy-phenomenology"),reader);
  Require(!result.ok()&&result.status.code==SEP::Core::StatusCode::InvalidConfiguration,
      "unsupported sheath silently selected another closure");

  auto malformedAmbient=assets;
  malformedAmbient["ambient.asset"]+="unreviewed_floor_t=1e-9\n";
  result=SEP::CoronaSwcme::ResolveEventConfiguration(Configuration(malformedAmbient),
      [&](const std::string& path){return Read(malformedAmbient,path);});
  Require(!result.ok(),"unknown ambient physics input survived checksum validation");

  auto overshootAssets=assets;
  overshootAssets["history.csv"]=History(true);
  result=SEP::CoronaSwcme::ResolveEventConfiguration(
      Configuration(overshootAssets),[&](const std::string& path){return Read(overshootAssets,path);});
  Require(!result.ok(),"negative between-knot geometry escaped continuous certification");

  auto missingAssets=assets;
  missingAssets.erase("ejecta.asset");
  result=SEP::CoronaSwcme::ResolveEventConfiguration(valid,
      [&](const std::string& path){return Read(missingAssets,path);});
  Require(!result.ok(),"missing mandatory ejecta asset accepted");

  auto renamedAssets=assets;
  renamedAssets["elsewhere/history.csv"]=renamedAssets["history.csv"];
  std::string renamed=valid;
  const std::string from="assets.history_file=history.csv";
  const std::string to="assets.history_file=elsewhere/history.csv";
  renamed.replace(renamed.find(from),from.size(),to);
  const auto a=Resolve(assets,valid),b=Resolve(renamedAssets,renamed);
  Require(a->physicsFingerprint==b->physicsFingerprint,
      "provenance path spelling changed physical identity");
}

std::string FileBytes(const std::string& path) {
  std::ifstream input(path,std::ios::binary);
  Require(static_cast<bool>(input),"cannot read frozen example: "+path);
  std::ostringstream bytes;
  bytes<<input.rdbuf();
  Require(input.good()||input.eof(),"failed reading frozen example: "+path);
  return bytes.str();
}

void TestFrozenExample() {
  const std::string input=FileBytes("examples/bg3d1/event.conf");
  const auto result=SEP::CoronaSwcme::ResolveEventConfiguration(input,
      [](const std::string& path) {
        try {
          return Result<std::string>::Success(FileBytes(path));
        } catch(const std::exception& error) {
          return Result<std::string>::Failure(
              SEP::Core::StatusCode::DataIntegrityFailure,error.what());
        }
      });
  Require(result.ok(),"frozen BG3D-1 example failed: "+result.status.message);
  const auto final=result.value->At(result.value->support.endS);
  Require(final.ok()&&final.value.apexRadiusM.value>149597870700.0,
      "frozen example does not cover its actual 1-AU propagation time");
  std::cout<<"[CMBGU02-EVIDENCE] fingerprint="<<result.value->physicsFingerprint
           <<" end_s="<<result.value->support.endS
           <<" final_apex_m="<<final.value.apexRadiusM.value<<'\n';
}

} // namespace

int main() {
  try {
    TestValidAndIndependentHandoff();
    TestStrictIdentityAndNegatives();
    TestFrozenExample();
    std::cout<<"[CMBGU01] PASS strict background configuration/assets/capabilities\n"
             <<"[CMBGU02] PASS continuous launch history and matched DBM handoff\n";
    return 0;
  } catch(const std::exception& error) {
    std::cerr<<"BG3D-1 FAIL: "<<error.what()<<'\n';
    return 1;
  }
}
