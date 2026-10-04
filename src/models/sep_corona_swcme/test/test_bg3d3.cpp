#include "sep_corona_swcme/surface_shock.h"
#include "sep_coronal_cme/configuration_parser.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>

namespace {

using SEP::CoronaSwcme::AmbientModel;
using SEP::CoronaSwcme::EventConfiguration;
using SEP::CoronaSwcme::SurfaceShockModel;
using SEP::CoronalCME::FixedOrientationEllipsoid;
using SEP::CoronalCME::Norm;

void Require(bool condition,const std::string& message) {
  if(!condition)throw std::runtime_error(message);
}

std::string FileBytes(const std::string& path) {
  std::ifstream input(path,std::ios::binary);
  Require(static_cast<bool>(input),"cannot read shock fixture: "+path);
  std::ostringstream bytes;bytes<<input.rdbuf();
  Require(input.good()||input.eof(),"failed reading shock fixture: "+path);
  return bytes.str();
}

std::shared_ptr<const EventConfiguration> Resolve(
    const std::string& configuration,const std::string& history) {
  const auto result=SEP::CoronaSwcme::ResolveEventConfiguration(configuration,
      [&](const std::string& path) {
        try {
          if(path=="examples/bg3d1/history.csv")
            return SEP::Core::Result<std::string>::Success(history);
          return SEP::Core::Result<std::string>::Success(FileBytes(path));
        } catch(const std::exception& error) {
          return SEP::Core::Result<std::string>::Failure(
              SEP::Core::StatusCode::DataIntegrityFailure,error.what());
        }
      });
  Require(result.ok(),result.status.message);
  return result.value;
}

std::shared_ptr<const EventConfiguration> FrozenEvent() {
  return Resolve(FileBytes("examples/bg3d1/event.conf"),
      FileBytes("examples/bg3d1/history.csv"));
}

std::shared_ptr<const EventConfiguration> SlowEvent() {
  const std::string history=
      "time_s,center_m,center_rate_m_s,radial_m,radial_rate_m_s,"
      "lateral1_m,lateral1_rate_m_s,lateral2_m,lateral2_rate_m_s\n"
      "0,800000000,1000,200000000,500,300000000,400,250000000,300\n"
      "200,800200000,1000,200100000,500,300080000,400,250060000,300\n";
  std::string configuration=FileBytes("examples/bg3d1/event.conf");
  const std::string oldDigest="06f55a7bca32d815e12419f8ecabd0874eeb8df6bd83f16c2e1c8448340795e3";
  const std::string newDigest=SEP::CoronalCME::ComputeContentChecksum(history);
  const std::size_t position=configuration.find(oldDigest);
  Require(position!=std::string::npos,"history checksum not found in fixture");
  configuration.replace(position,oldDigest.size(),newDigest);
  return Resolve(configuration,history);
}

std::shared_ptr<SurfaceShockModel> Model(
    const std::shared_ptr<const EventConfiguration>& event) {
  const auto ambient=AmbientModel::Create(event);
  Require(ambient.ok(),ambient.status.message);
  const auto model=SurfaceShockModel::Create(event,ambient.value);
  Require(model.ok(),model.status.message);
  return model.value;
}

void TestGeometryAndNormalVelocity() {
  const auto event=FrozenEvent();
  const auto prepared=Model(event)->Prepare(10000,11,12,24);
  Require(prepared.ok(),"composite surface preparation failed: "+prepared.status.message);
  const auto& epoch=*prepared.value;
  Require(!epoch.frontPatches.empty()&&!epoch.contactPatches.empty()&&
      epoch.frontPatches.size()==epoch.shocks->patches.size(),
      "front/contact/shock geometry is incomplete");
  std::set<std::uint64_t> identities;
  for(std::size_t i=0;i<epoch.frontPatches.size();++i) {
    const auto& patch=epoch.frontPatches[i];
    Require(patch.geometry.physicalId>0&&identities.insert(
        patch.geometry.physicalId).second,"front identity is zero or duplicated");
    Require(epoch.shocks->patches[i].stableId==patch.geometry.physicalId&&
        patch.geometry.areaM2>0&&std::abs(Norm(patch.geometry.outwardNormal)-1)<1e-12,
        "shock/geometry identity, area, or normal is invalid");
  }
  for(const auto& patch:epoch.contactPatches) {
    const auto inside=epoch.front.Evaluate(patch.centerM);
    Require(inside.ok()&&inside.value.implicitValue<=1e-10,
        "contact patch crosses the shock front");
  }
  Require(epoch.contact.ApexRadiusM()<epoch.front.ApexRadiusM(),
      "contact/ejecta boundary is not distinct from the shock front");

  const double step=0.01;
  const auto beforeState=event->At(10000-step),afterState=event->At(10000+step);
  Require(beforeState.ok()&&afterState.ok(),"normal-speed oracle lacks history");
  const auto before=FixedOrientationEllipsoid::FromCenter(event->basis,
      beforeState.value.ellipsoid);
  const auto after=FixedOrientationEllipsoid::FromCenter(event->basis,
      afterState.value.ellipsoid);
  Require(before.ok()&&after.ok(),"normal-speed oracle shape failed");
  for(std::size_t i=0;i<epoch.frontPatches.size();i+=17) {
    const auto& patch=epoch.frontPatches[i];
    const auto x0=before.value.Point(patch.geometry.polarParameterRad,
        patch.geometry.azimuthParameterRad);
    const auto x1=after.value.Point(patch.geometry.polarParameterRad,
        patch.geometry.azimuthParameterRad);
    const double oracle=SEP::CoronalCME::Dot((x1-x0)/(2*step),
        patch.geometry.outwardNormal);
    Require(std::abs(oracle-patch.normalSpeedMPerS)<=2e-8*
        std::max(1.0,std::abs(oracle)),"complete point normal speed is wrong");
  }
}

double Area(const EventConfiguration& event,double time,int polar,int azimuth) {
  const auto state=event.At(time);
  Require(state.ok(),state.status.message);
  const auto shape=FixedOrientationEllipsoid::FromCenter(event.basis,
      state.value.ellipsoid,event.support.solarRadiusM);
  Require(shape.ok(),shape.status.message);
  const auto patches=shape.value.Tessellate(polar,azimuth);
  Require(patches.ok(),patches.status.message);
  double area=0;for(const auto& patch:patches.value)area+=patch.areaM2;
  return area;
}

void TestAreaConvergenceAndShockState() {
  const auto event=FrozenEvent();
  const double a0=Area(*event,10000,12,24),a1=Area(*event,10000,24,48);
  const double a2=Area(*event,10000,48,96),reference=Area(*event,10000,96,192);
  const double e0=std::abs(a0-reference),e1=std::abs(a1-reference),e2=std::abs(a2-reference);
  Require(e2<e1&&e1<e0,"surface area does not converge under refinement");

  auto model=Model(event);
  const auto result=model->Prepare(10000,21,24,48);
  Require(result.ok(),"fast shock surface failed: "+result.status.message);
  std::size_t fast=0,subfast=0;
  double maximumResidual=0,minimumCompression=1e300;
  for(const auto& patch:result.value->shocks->patches) {
    if(patch.fast) {
      ++fast;
      Require(patch.jump.upstreamFastMach>1&&patch.jump.compressionRatio>1&&
          patch.jump.upstreamCharacteristics.obliquityRad>=0&&
          patch.jump.upstreamCharacteristics.obliquityRad<=3.14159265358979323846/2,
          "fast patch lacks Mach/compression/obliquity diagnostics");
      maximumResidual=std::max(maximumResidual,patch.jump.residuals.maximum);
      minimumCompression=std::min(minimumCompression,patch.jump.compressionRatio);
      Require(!patch.sourceEligibleBeforeClearance&&!patch.sourceActive&&
          patch.physicalNumberRatePerS==0&&patch.physicalKineticEnergyRateW==0,
          "background shock acquired particle-source eligibility or rate");
    } else ++subfast;
  }
  Require(fast>0&&maximumResidual<=1e-9&&minimumCompression>1,
      "fast surface lacks admissible both-sided RH state");
  Require(result.value->shocks->measures.counterfactualAreaM2==0&&
      result.value->shocks->measures.sourceActiveAreaM2==0&&
      result.value->shocks->measures.counterfactualNumberRatePerS==0&&
      result.value->shocks->measures.counterfactualKineticEnergyRateW==0,
      "background shock acquired aggregate particle-source measures");

  const auto slowEvent=SlowEvent();
  const auto slow=Model(slowEvent)->Prepare(0,1,12,24);
  Require(slow.ok(),"sub-fast geometry preparation failed: "+slow.status.message);
  std::size_t slowFast=0;
  for(const auto& patch:slow.value->shocks->patches)if(patch.fast)++slowFast;
  Require(slowFast==0&&slow.value->shocks->measures.fastAreaM2==0,
      "sub-fast front was promoted into a shock");

  const auto committed=model->Current();
  Require(committed&&committed->shocks->generation==result.value->shocks->generation,
      "prepared shock epoch was not committed");
  const auto rejected=model->Prepare(-1,22,24,48);
  Require(!rejected.ok()&&model->Current()==committed,
      "failed candidate changed the committed surface/shock epoch");
  std::cout<<"[CMBGU04-EVIDENCE] patches="<<result.value->shocks->patches.size()
           <<" fast="<<fast<<" subfast="<<subfast
           <<" max_rh_residual="<<maximumResidual
           <<" min_compression="<<minimumCompression
           <<" area_errors="<<e0<<','<<e1<<','<<e2<<'\n';
}

} // namespace

int main() {
  try {
    TestGeometryAndNormalVelocity();
    TestAreaConvergenceAndShockState();
    std::cout<<"[CMBGU04] PASS complete geometry and local shock diagnostics\n";
    return 0;
  } catch(const std::exception& error) {
    std::cerr<<"BG3D-3 FAIL: "<<error.what()<<'\n';
    return 1;
  }
}
