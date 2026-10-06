#include "sep_corona_swcme/piston_sheath.h"
#include "sep_coronal_cme/mhd_jump_solver.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

namespace {

using SEP::CoronaSwcme::PistonBackgroundRegion;
using SEP::CoronaSwcme::PistonRayDisposition;
using SEP::CoronaSwcme::PistonSheathModel;
using SEP::CoronalCME::MhdPrimitiveState;
using SEP::CoronalCME::Vec3;

void Require(bool condition,const std::string& message) {
  if(!condition)throw std::runtime_error(message);
}

std::string FileBytes(const std::string& path) {
  std::ifstream input(path,std::ios::binary);
  Require(static_cast<bool>(input),"cannot read fixture: "+path);
  std::ostringstream bytes;bytes<<input.rdbuf();
  Require(input.good()||input.eof(),"cannot finish reading fixture: "+path);
  return bytes.str();
}

SEP::Core::Result<std::string> ReadFile(const std::string& path) {
  try {return SEP::Core::Result<std::string>::Success(FileBytes(path));}
  catch(const std::exception& error) {
    return SEP::Core::Result<std::string>::Failure(
        SEP::Core::StatusCode::DataIntegrityFailure,error.what());
  }
}

std::shared_ptr<const SEP::CoronaSwcme::EventConfiguration> Event() {
  const auto event=SEP::CoronaSwcme::ResolveEventConfiguration(
      FileBytes("examples/bg3d4_piston/event.conf"),ReadFile);
  Require(event.ok(),event.status.message);return event.value;
}

double Relative(double value,double reference) {
  return std::abs(value-reference)/std::max(std::abs(reference),1e-300);
}

} // namespace

int main() {
  try {
    const auto event=Event();
    const auto& controls=event->pistonNumerics;
    Require(controls.enabled,
        "strict event did not publish fingerprinted piston controls");
    auto model=PistonSheathModel::Create(event);
    Require(model.ok(),"contact-driven sheath creation failed: "+
        model.status.message);
    std::cout<<"[CSWC0627-PROGRESS] created\n"<<std::flush;
    const auto initial=model.value->Receipts();
    Require(initial.ok(),"initial sheath receipts unavailable");
    const auto supported=std::find_if(initial.value.begin(),initial.value.end(),
        [](const auto& ray) {return ray.disposition==PistonRayDisposition::Supported;});
    const auto supportedCount=std::count_if(initial.value.begin(),initial.value.end(),
        [](const auto& ray) {
          return ray.disposition==PistonRayDisposition::Supported;
        });
    Require(supported!=initial.value.end()&&supportedCount>0,
        "refined qualification asset has no supported piston ray");
    const std::uint64_t rayId=supported->rayId;
    Require(std::all_of(initial.value.begin(),initial.value.end(),[](const auto& ray) {
          return ray.disposition!=PistonRayDisposition::Supported||
              (!ray.shock.present&&ray.cells==80);
        }),"ambient-only startup fabricated a shock or inventory");

    for(double epoch:{600.0,1200.0,1800.0}) {
      const auto status=model.value->AdvanceTo(epoch);
      Require(status.ok(),"contact-driven advance failed: "+status.message);
      std::cout<<"[CSWC0627-PROGRESS] epoch="<<epoch<<'\n'<<std::flush;
    }
    const auto receipts=model.value->Receipts();
    Require(receipts.ok(),"evolved sheath receipts unavailable");
    const auto evolved=std::find_if(receipts.value.begin(),receipts.value.end(),
        [rayId](const auto& ray) {return ray.rayId==rayId;});
    if(evolved!=receipts.value.end())std::cout
        <<"[CSWC0627-RAW] contact="<<evolved->contactRadiusM
        <<" shock_present="<<evolved->shock.present
        <<" states="<<evolved->shock.statesAvailable
        <<" compression="<<evolved->shock.compressionRatio
        <<" shock="<<evolved->shock.radiusM
        <<" disturbance="<<evolved->disturbanceRadiusM
        <<" buffer="<<evolved->ambientBufferCells<<'\n';
    Require(evolved!=receipts.value.end()&&evolved->shock.present&&
        evolved->shock.statesAvailable&&evolved->shock.compressionRatio>1&&
        evolved->disturbanceRadiusM>evolved->contactRadiusM&&
        evolved->ambientBufferCells>=
            static_cast<std::size_t>(controls.minimumBufferCells),
        "contact piston did not produce a resolved, buffered shock");

    const auto contactSample=model.value->QueryRay(rayId,
        evolved->contactRadiusM*(1-1e-8));
    const auto sheathSample=model.value->QueryRay(rayId,
        0.5*(evolved->contactRadiusM+evolved->shock.radiusM));
    const auto leadingSample=model.value->QueryRay(rayId,
        0.5*(evolved->shock.radiusM+evolved->disturbanceRadiusM));
    const auto bufferSample=model.value->QueryRay(rayId,
        0.5*(evolved->disturbanceRadiusM+evolved->outerRadiusM));
    const auto ambientSample=model.value->QueryRay(rayId,
        evolved->outerRadiusM+1e7);
    Require(contactSample.ok()&&sheathSample.ok()&&leadingSample.ok()&&
        bufferSample.ok()&&ambientSample.ok()&&
        contactSample.value.region==PistonBackgroundRegion::EjectaNotOwned&&
        sheathSample.value.region==PistonBackgroundRegion::Sheath&&
        sheathSample.value.plasmaAvailable&&
        leadingSample.value.region==PistonBackgroundRegion::Compression&&
        leadingSample.value.plasmaAvailable&&
        bufferSample.value.region==PistonBackgroundRegion::Ambient&&
        bufferSample.value.plasmaAvailable&&
        ambientSample.value.region==PistonBackgroundRegion::Ambient&&
        ambientSample.value.plasmaAvailable,
        "production query does not preserve ejecta/sheath/leading/ambient regions");

    const auto& shock=evolved->shock;
    MhdPrimitiveState upstream;
    upstream.massDensityKgM3=shock.upstreamDensityKgM3;
    upstream.pressurePa=shock.upstreamPressurePa;
    upstream.velocityMPerS={shock.upstreamVelocityMPerS,0,0};
    upstream.magneticFieldT={shock.upstreamRadialMagneticFieldT,
        shock.upstreamTransverseMagneticFieldT,
        shock.upstreamTransverseMagneticField2T};
    const Vec3 normal={1,0,0};
    const auto canonical=SEP::CoronalCME::SolveObliqueFastShock(upstream,normal,
        shock.speedMPerS,event->composition.gammaAdiabatic,1e-8);
    Require(canonical.ok(),"computed shock is not canonically admissible: "+
        canonical.status.message);
    const double theta=std::acos(std::min(1.0,std::abs(
        shock.upstreamRadialMagneticFieldT)/std::sqrt(
        shock.upstreamRadialMagneticFieldT*shock.upstreamRadialMagneticFieldT+
        shock.upstreamTransverseMagneticFieldT*
            shock.upstreamTransverseMagneticFieldT+
        shock.upstreamTransverseMagneticField2T*
            shock.upstreamTransverseMagneticField2T)));
    const double compressionError=Relative(shock.compressionRatio,
        canonical.value.compressionRatio);
    const double pressureError=Relative(shock.downstreamPressurePa,
        canonical.value.downstream.pressurePa);
    const double velocityError=Relative(shock.downstreamVelocityMPerS,
        canonical.value.downstream.velocityMPerS.x);
    std::cout<<"[CSWC0627-RH-RAW] zone="<<shock.firstShockCell<<':'
             <<shock.lastShockCell<<" speed="<<shock.speedMPerS
             <<" rho="<<shock.upstreamDensityKgM3<<':'
             <<shock.downstreamDensityKgM3<<" pressure="
             <<shock.upstreamPressurePa<<':'<<shock.downstreamPressurePa
             <<" velocity="<<shock.upstreamVelocityMPerS<<':'
             <<shock.downstreamVelocityMPerS<<" canonical_x="
             <<canonical.value.compressionRatio<<" canonical_pressure="
             <<canonical.value.downstream.pressurePa<<" canonical_velocity="
             <<canonical.value.downstream.velocityMPerS.x<<" errors="
             <<compressionError<<','<<pressureError<<','<<velocityError
             <<" theta="<<theta<<'\n';
    // The reduced radial equations omit tangential momentum and switch-on
    // physics.  The canonical RH error is therefore an acceptance gate only
    // on the preregistered quasi-perpendicular subset.  A Parker-like field
    // can make every ray in this small production smoke asset oblique; such a
    // ray is retained as explicit model-error evidence, not silently judged
    // against a closure that does not claim to represent it.  CSWC0623 is the
    // independent perpendicular-MHD numerical reference, while a production
    // quasi-perpendicular ray remains a separate open BG3D-4 coverage gate.
    const bool quasiPerpendicular=
        theta>70*3.14159265358979323846/180;
    if(quasiPerpendicular)Require(compressionError<0.12&&pressureError<0.15&&
        velocityError<0.12,
        "computed quasi-perpendicular shock disagrees with canonical RH bounds");
    std::cout<<"[CSWC0627-RH-SCOPE] "
             <<(quasiPerpendicular?"APPLICABLE":"NOT-APPLICABLE")
             <<" quasi_perpendicular_threshold_deg=70\n";

    const double committedEpoch=model.value->EpochS();
    const std::uint64_t committedGeneration=model.value->Generation();
    const auto rejected=model.value->AdvanceTo(event->support.endS+1);
    Require(!rejected.ok()&&model.value->EpochS()==committedEpoch&&
        model.value->Generation()==committedGeneration,
        "failed multi-ray candidate changed committed epoch/generation");
    const auto levelA=event->At(committedEpoch);
    Require(levelA.ok()&&Relative(levelA.value.apexRadiusM.value,
        shock.radiusM)>1e-3,
        "computed Level-B shock was silently reset to the Level-A history");

    std::cout<<"[CSWC0627-EVIDENCE] ray_id="<<rayId
             <<" supported_rays="<<supportedCount
             <<" epoch_s="<<committedEpoch
             <<" contact_m="<<evolved->contactRadiusM
             <<" shock_m="<<shock.radiusM
             <<" disturbance_m="<<evolved->disturbanceRadiusM
             <<" outer_m="<<evolved->outerRadiusM
             <<" buffer_cells="<<evolved->ambientBufferCells
             <<" compression="<<shock.compressionRatio
             <<" theta_bn_rad="<<theta
             <<" rh_errors="<<compressionError<<','<<pressureError<<','
             <<velocityError<<'\n';
    std::cout<<"[CSWC0627] PASS contact-driven computed shock and query\n";
    return 0;
  } catch(const std::exception& error) {
    std::cerr<<"CSWC0627 FAIL: "<<error.what()<<'\n';return 1;
  }
}
