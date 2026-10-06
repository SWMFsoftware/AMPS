#include "sep_corona_swcme/piston_sheath.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

namespace {

using SEP::CoronaSwcme::EventConfiguration;
using SEP::CoronaSwcme::PistonBackgroundRegion;
using SEP::CoronaSwcme::PistonRayDisposition;
using SEP::CoronaSwcme::PistonSheathModel;
constexpr double kMu0=1.25663706212e-6;

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

void Replace(std::string* text,const std::string& oldValue,
    const std::string& newValue) {
  const std::size_t where=text->find(oldValue);
  Require(where!=std::string::npos&&
      text->find(oldValue,where+oldValue.size())==std::string::npos,
      "refinement fixture key is absent or duplicated: "+oldValue);
  text->replace(where,oldValue.size(),newValue);
}

struct Level {
  int cells=0;
  int firstShockCell=-1,lastShockCell=-1;
  int shockZones=0;
  double maximumQOverP=0,selectedMaximumQOverP=0;
  double cfl=0,shockRadiusM=0,compression=0,shockSpeedMPerS=0;
  double fastMach=0,contactTotalPressurePa=0;
  std::string fingerprint;
};

std::shared_ptr<const EventConfiguration> RefinedEvent(int cells,
    const std::string& cfl) {
  std::string input=FileBytes("examples/bg3d4_piston/event.conf");
  Replace(&input,"bg3d4.initial_cells=80",
      "bg3d4.initial_cells="+std::to_string(cells));
  Replace(&input,"bg3d4.minimum_buffer_cells=24",
      "bg3d4.minimum_buffer_cells="+std::to_string(24*cells/80));
  Replace(&input,"bg3d4.append_cells=64",
      "bg3d4.append_cells="+std::to_string(64*cells/80));
  Replace(&input,"bg3d4.cfl=0.2","bg3d4.cfl="+cfl);
  const auto event=SEP::CoronaSwcme::ResolveEventConfiguration(input,ReadFile);
  Require(event.ok(),"refined event rejected: "+event.status.message);
  return event.value;
}

Level Run(int cells,const std::string& cflText,double targetS=1800.0) {
  const auto event=RefinedEvent(cells,cflText);
  auto model=PistonSheathModel::Create(event);
  Require(model.ok(),"refined piston model rejected: "+model.status.message);
  const auto initial=model.value->Receipts();
  Require(initial.ok(),"refined initial receipts unavailable");
  const auto selected=std::find_if(initial.value.begin(),initial.value.end(),
      [](const auto& ray) {
        return ray.disposition==PistonRayDisposition::Supported;
      });
  Require(selected!=initial.value.end(),"refined event has no supported ray");
  const auto advanced=model.value->AdvanceTo(targetS);
  Require(advanced.ok(),"refined production advance failed: "+advanced.message);
  const auto receipts=model.value->Receipts();
  Require(receipts.ok(),"refined production receipts unavailable");
  const auto ray=std::find_if(receipts.value.begin(),receipts.value.end(),
      [&](const auto& candidate) {return candidate.rayId==selected->rayId;});
  Require(ray!=receipts.value.end()&&ray->shock.present&&
      ray->shock.statesAvailable,"refined production shock is unresolved");
  const auto contact=model.value->QueryRay(ray->rayId,
      ray->contactRadiusM*(1+1e-8));
  Require(contact.ok()&&contact.value.plasmaAvailable&&
      contact.value.region==PistonBackgroundRegion::Sheath,
      "refined contact-side sheath state is unavailable");
  const auto& s=ray->shock;
  const double b2=s.upstreamRadialMagneticFieldT*
      s.upstreamRadialMagneticFieldT+
      s.upstreamTransverseMagneticFieldT*
          s.upstreamTransverseMagneticFieldT+
      s.upstreamTransverseMagneticField2T*
          s.upstreamTransverseMagneticField2T;
  const double br2=s.upstreamRadialMagneticFieldT*
      s.upstreamRadialMagneticFieldT;
  const double sound2=event->composition.gammaAdiabatic*
      s.upstreamPressurePa/s.upstreamDensityKgM3;
  const double alfven2=b2/(kMu0*s.upstreamDensityKgM3);
  const double discriminant=std::max(0.0,
      (sound2+alfven2)*(sound2+alfven2)-
      4*sound2*br2/(kMu0*s.upstreamDensityKgM3));
  const double fast=std::sqrt(0.5*(sound2+alfven2+
      std::sqrt(discriminant)));
  const auto& cell=contact.value.plasma;
  const double contactBt2=cell.transverseMagneticFieldT*
      cell.transverseMagneticFieldT+
      cell.transverseMagneticField2T*cell.transverseMagneticField2T;
  Level out;
  out.cells=cells;out.cfl=event->pistonNumerics.cfl;
  out.firstShockCell=s.firstShockCell;out.lastShockCell=s.lastShockCell;
  out.shockZones=s.shockZoneCount;
  out.maximumQOverP=s.maximumArtificialPressureRatio;
  out.selectedMaximumQOverP=s.selectedZoneMaximumArtificialPressureRatio;
  out.shockRadiusM=s.radiusM;out.compression=s.compressionRatio;
  out.shockSpeedMPerS=s.speedMPerS;
  out.fastMach=(s.speedMPerS-s.upstreamVelocityMPerS)/fast;
  out.contactTotalPressurePa=cell.pressurePa+contactBt2/(2*kMu0);
  out.fingerprint=event->physicsFingerprint;
  std::cout<<"[CSWC0628-LEVEL] cells="<<out.cells<<" cfl="<<out.cfl
           <<" epoch_s="<<targetS
           <<" zone="<<out.firstShockCell<<':'<<out.lastShockCell
           <<" zones="<<out.shockZones<<" max_q_over_p="<<out.maximumQOverP
           <<" selected_max_q_over_p="<<out.selectedMaximumQOverP
           <<" shock_m="<<out.shockRadiusM
           <<" compression="<<out.compression
           <<" speed_m_s="<<out.shockSpeedMPerS
           <<" fast_mach="<<out.fastMach
           <<" contact_total_pressure_pa="<<out.contactTotalPressurePa
           <<" fingerprint="<<out.fingerprint<<'\n'<<std::flush;
  return out;
}

double RelativeDifference(double a,double b) {
  return std::abs(a-b)/std::max({std::abs(a),std::abs(b),1e-300});
}

struct Changes {
  double radius=0,compression=0,speed=0,mach=0,contactPressure=0;
};

Changes Difference(const Level& a,const Level& b) {
  return {RelativeDifference(a.shockRadiusM,b.shockRadiusM),
      RelativeDifference(a.compression,b.compression),
      RelativeDifference(a.shockSpeedMPerS,b.shockSpeedMPerS),
      RelativeDifference(a.fastMach,b.fastMach),
      RelativeDifference(a.contactTotalPressurePa,b.contactTotalPressurePa)};
}

bool ComponentwiseLess(const Changes& fine,const Changes& coarse,double slack) {
  return fine.radius<=slack*coarse.radius&&
      fine.compression<=slack*coarse.compression&&
      fine.speed<=slack*coarse.speed&&fine.mach<=slack*coarse.mach&&
      fine.contactPressure<=slack*coarse.contactPressure;
}

void Print(const char* name,const Changes& value) {
  std::cout<<name<<"="<<value.radius<<','<<value.compression<<','
           <<value.speed<<','<<value.mach<<','<<value.contactPressure<<'\n';
}

} // namespace

int main(int argc,char** argv) {
  try {
    // A single-level mode is retained for diagnosing a failed refinement
    // without silently rerunning or retuning the remaining preregistered
    // sequence.  It does not emit the suite PASS marker.
    if(argc==3||argc==4) {
      Run(std::stoi(argv[1]),argv[2],argc==4?std::stod(argv[3]):1800.0);
      return 0;
    }
    Require(argc==1,
        "usage: test_bg3d4_piston_refinement [cells cfl [target_s]]");
    // The first preregistered 1800-s attempt was rejected as a convergence
    // fixture because it contained 3/2/4 disjoint Q/P zones at 40/80/160
    // cells: different grids were not measuring the same shock branch.  The
    // replacement epoch is fixed physically at 3*T_ramp=3600 s, two complete
    // ramp durations after acceleration ends, and requires exactly one zone
    // at every level before any numerical change can be graded.  This changes
    // an incompatible initial/evaluation condition, not an accuracy bound.
    //
    // Frozen before executing the replacement: a captured VNR shock is formally
    // first order, so halving mass-cell width should reduce successive
    // changes.  Ten-percent monotonic slack permits a detector/roundoff floor,
    // not a growing sequence.  Fine spatial changes must be <3% in geometry,
    // <8% in jump diagnostics and <10% at the piston.  Fixed-grid RK4 time
    // changes must be <1%; a 25% monotonic slack allows spatial shock error to
    // dominate once temporal truncation reaches that floor.
    constexpr double evaluationS=3600.0;
    const Level mass0=Run(40,"0.2",evaluationS),
        mass1=Run(80,"0.2",evaluationS),
        mass2=Run(160,"0.2",evaluationS);
    const Changes massCoarse=Difference(mass0,mass1),
        massFine=Difference(mass1,mass2);
    Print("[CSWC0628-MASS-COARSE]",massCoarse);
    Print("[CSWC0628-MASS-FINE]",massFine);
    Require(mass0.shockZones==1&&mass1.shockZones==1&&mass2.shockZones==1,
        "production refinement does not identify one common shock branch");
    Require(ComponentwiseLess(massFine,massCoarse,1.1),
        "production mass-refinement changes do not decrease");
    Require(massFine.radius<0.03&&massFine.compression<0.08&&
        massFine.speed<0.08&&massFine.mach<0.08&&
        massFine.contactPressure<0.10,
        "production fine mass-refinement change exceeds frozen bounds");

    const Level time0=Run(40,"0.4",evaluationS),time1=mass0,
        time2=Run(40,"0.1",evaluationS);
    const Changes timeCoarse=Difference(time0,time1),
        timeFine=Difference(time1,time2);
    Print("[CSWC0628-TIME-COARSE]",timeCoarse);
    Print("[CSWC0628-TIME-FINE]",timeFine);
    Require(ComponentwiseLess(timeFine,timeCoarse,1.25),
        "production timestep-refinement changes grow above the spatial floor");
    Require(timeFine.radius<0.01&&timeFine.compression<0.01&&
        timeFine.speed<0.01&&timeFine.mach<0.01&&
        timeFine.contactPressure<0.01,
        "production fine timestep change exceeds frozen one-percent bound");
    Require(mass0.fingerprint!=mass1.fingerprint&&
        mass1.fingerprint!=mass2.fingerprint&&
        time0.fingerprint!=mass0.fingerprint&&
        time2.fingerprint!=mass0.fingerprint,
        "refined controls did not change immutable event identity");
    std::cout<<"[CSWC0628-MASS-TIME] PASS production-interface refinement\n";
    return 0;
  } catch(const std::exception& error) {
    std::cerr<<"CSWC0628-MASS-TIME FAIL: "<<error.what()<<'\n';return 1;
  }
}
