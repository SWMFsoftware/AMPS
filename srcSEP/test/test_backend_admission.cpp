#include "sep_backend_admission.h"
#include "swcme1d_adapter.h"
#include "field_line_bundle_adapter.h"
#include "swcme1d_input.hpp"
#include <cmath>
#include <filesystem>
#include <functional>
#include <iostream>
#include <limits>
#include <stdexcept>

namespace A=SEP::Background::Admission;
namespace FL=SEP::FieldLine;
namespace BA=SEP1D::Adapters;
namespace SW=SEP::SW1DAdapter;

void require(bool condition,const std::string& message) {
  if(!condition)throw std::runtime_error(message);
}

int main(int argc,char** argv) {
  if(argc!=2)return 2;
  unsigned failures=0;
  const auto run=[&](const char* id,const std::function<void()>& test) {
    try {test();std::cout<<'['<<id<<"] PASS\n";}
    catch(const std::exception& error) {++failures;std::cout<<'['<<id<<"] FAIL "<<error.what()<<'\n';}
  };
  run("SBACK01",[] {
    require(A::RequireMovingTransport(A::Backend::ReducedSseParker).ok(),"existing reduced profile lost backend admission");
    require(!A::RequireMovingTransport(static_cast<A::Backend>(99)).ok(),"unknown provider acquired a fallback");
  });
  run("SBACK02",[] {
    // Independent algebraic continuity test, not a reference generated from
    // the queried density routine: for n=C2/r^2+C4/r^4+C6/r^6 and constant U,
    // div(n U)=-U(2 C4/r^5+4 C6/r^7). Its nonzero value is the missing mass
    // source per particle mass, not numerical noise to floor at large radius.
    const auto resolved=swcme::input1d::Resolve(swcme::input1d::Preset::Fast,{});
    require(resolved.ok(),"canonical fast profile did not resolve");
    swcme1d::Model model(resolved.configuration.model);
    const auto state=model.prepare_step(0);
    const auto& wind=state.common.solar_wind;
    const double r=20*swcme1d::Rs;
    const double defect=-wind.V_sw_m_s*(2*wind.C4/std::pow(r,5)+4*wind.C6/std::pow(r,7));
    require(std::isfinite(defect)&&defect<0,"profile no longer exhibits the audited continuity defect; review its support contract");
    const double f1=r*r*swcme::solarwind::density_m3(wind,r)*wind.V_sw_m_s;
    const double f2=4*r*r*swcme::solarwind::density_m3(wind,2*r)*wind.V_sw_m_s;
    require(f1>f2,"material radial flux is unexpectedly constant");
    require(!A::RequireMovingTransport(A::Backend::Swcme1D).ok(),"nonconservative SWCME material mode admitted");
    std::cout<<"METRIC SBACK02 continuity_defect_per_m3_s="<<defect<<" radial_flux_ratio="<<f1/f2<<'\n';
  });
  run("SBACK03",[] {
    // Rejection of movement must preserve the maintained diagnostic adapter.
    // Compare actual canonical SI values and explicit epoch/state identity;
    // no second ambient model is installed to make this comparison pass.
    SW::ConfigurationRequest request;
    require(SW::Configure(request).ok(),"diagnostic adapter configure failed");
    require(SW::PrepareState(0).ok(),"diagnostic adapter preparation failed");
    const auto resolved=swcme::input1d::Resolve(swcme::input1d::Preset::Fast,{});
    require(resolved.ok(),"direct canonical profile failed");
    const auto summary=SW::GetConfigurationSummary();
    require(summary.fingerprint==resolved.configuration.fingerprint&&
        summary.normalized_manifest.find("fixed.frame=HCI_like_inertial")!=std::string::npos,
        "adapter lost its canonical configuration/frame identity");
    swcme1d::Model model(resolved.configuration.model);const auto state=model.prepare_step(0);
    for(double r:{20*swcme1d::Rs,swcme1d::AU}) {
      const auto q=SW::QueryAtRadius(r,false);
      // Sentinels make an unfilled canonical output fail visibly even if a
      // future evaluator incorrectly reports success without writing a field.
      double n=std::numeric_limits<double>::quiet_NaN(),u=n,br=n,bphi=n,b=n,div=n;
      const auto status=model.evaluate_radii_with_B_div_checked(state,&r,&n,&u,&br,&bphi,&b,&div,1);
      require(status.ok()&&q.ok()&&q.epoch_seconds==0&&q.source_state_id>0,"diagnostic query/epoch missing");
      require(q.sample.number_density_m3()==n&&q.sample.speed_m_s()==u&&q.sample.divergence_s_inv()==div&&
          q.sample.radial_field_t()==br&&q.sample.azimuthal_field_t()==bphi&&q.sample.thermal_state_available(),
          "adapter differs from its direct core authority");
    }
    // Exercise an explicitly multispecies authority too: ne is not np when
    // alpha/proton=0.1. This independent charge-neutrality/mass sum detects a
    // density-basis error that a pure-H comparison would conceal. These are
    // declared diagnostic inputs, not a newly admitted material profile.
    request.command_line_assignments={
      {"ambient.thermodynamic_closure","multi_species","SBACK03",0},
      {"ambient.alpha_to_proton_ratio","0.1","SBACK03",0},
      {"ambient.electron_temperature","120000 K","SBACK03",0},
      {"ambient.alpha_temperature","480000 K","SBACK03",0},
      {"ambient.proton_temperature","200000 K","SBACK03",0}};
    require(SW::Configure(request).ok()&&SW::PrepareState(0).ok(),"explicit multispecies diagnostics unavailable");
    const auto multi=SW::QueryAtRadius(swcme1d::AU,false);
    require(multi.ok()&&multi.sample.magnetic_and_composition_available(),"multispecies composition export missing");
    const double ne=multi.sample.number_density_m3(),np=ne/1.2,na=0.1*np;
    const double rho=swcme::constants::PROTON_MASS_KG*np+swcme::constants::ALPHA_PARTICLE_MASS_KG*na;
    const double pressure=swcme::constants::BOLTZMANN_J_K*(np*200000+ne*120000+na*480000);
    const double rounding=8*std::numeric_limits<double>::epsilon();
    require(std::abs(multi.sample.mass_density_kg_m3()-rho)<=rounding*rho&&
        std::abs(multi.sample.pressure_pa()-pressure)<=rounding*pressure&&multi.sample.proton_temperature_k()==200000,
        "canonical multispecies density/temperature basis changed");
    const auto before=SW::GetDiagnostics();
    require(!SW::PrepareState(-1).ok(),"out-of-coverage epoch accepted");
    require(SW::GetDiagnostics().prepared_state_id==before.prepared_state_id,"failed epoch changed committed diagnostic state");
  });
  run("SBACK04",[&] {
    // A manufactured neutral exchange fixture tests ownership, not physical
    // transport. Complete connection history must still NOT become material
    // plasma history. Verify the real codec and adapter callback transaction.
    FL::FieldLineSet set;set.frame="HCI";set.epochS=0;set.historyBeginS=0;set.historyEndS=5;
    set.rotationProvenance="explicit-software-fixture";
    FL::LineRecord line;line.stableLineId="line-A";line.open=line.singleSmoothSector=true;
    line.historyCoverage=FL::HistoryCoverage::Complete;line.sourceMeasureStatus=FL::SourceMeasureStatus::FiniteTube;
    line.unsignedMagneticFluxWb=2;line.quadratureWeight=2;line.observerCoverage={0,5,true};
    for(int i=0;i<2;++i) {
      FL::NodeState node;node.arcLengthM=i;node.positionM={1.0+i,0,0};node.outwardTangent={1,0,0};
      node.magneticFieldT={2,0,0};node.plasmaVelocityMPerS={10,0,0};node.massDensityKgM3=4;
      node.pressurePa=2;node.temperatureK=1e6;node.focusingLengthM=10;node.tubeAreaM2=1;
      node.magneticSector=1;node.mappingValid=true;node.sourceLabel="fixture";line.nodes.push_back(node);
    }
    set.lines.push_back(line);set.bundleId=FL::BundleIdentity(set);
    const std::string directory=(std::filesystem::path(argv[1])/"bundle").string();
    const auto written=FL::WriteBundleTransactional(set,directory);require(written.ok(),written.message);
    BA::FieldLineBundleAdapter adapter;unsigned calls=0;
    const auto publisher=[&](const FL::FieldLineSet&) {++calls;return SEP::Core::Status::Success();};
    require(adapter.LoadAndPublish(directory,publisher).ok()&&calls==1,"snapshot exchange compatibility failed");
    const auto identity=adapter.BundleIdentity();
    const auto rejected=adapter.LoadAndPublish(directory,publisher,BA::BundleUse::MaterialTransport);
    require(!rejected.ok()&&calls==1&&adapter.BundleIdentity()==identity,"material rejection published or replaced a snapshot");
    require(adapter.Evaluate("line-A",0.5).ok(),"rejection destroyed the old exchange view");
    require(!adapter.Evaluate("missing-line",0.5).ok(),"absent physical identity acquired nearest-line fallback");
    require(!adapter.LoadAndPublish(directory,publisher,static_cast<BA::BundleUse>(99)).ok()&&calls==1,
        "unknown bundle use invoked publication");
  });
  run("SBACK05",[] {
    const auto status=A::RequireMovingTransport(A::Backend::FullCoronaSwcme);
    require(!status.ok()&&status.message.find("BG3D-4")!=std::string::npos,
        "unqualified full regional core admitted");
  });
  return failures?1:0;
}
