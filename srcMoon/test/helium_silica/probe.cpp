// Thin adapter to the actual production source-specific registry. No charge
// law, coefficient, geometry conversion or retention formula is duplicated.
#include "srcMoon/species/HeliumSilicaChargeSurvival.h"
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
namespace H=Moon::Species::HeliumSilicaChargeSurvival;
std::string Json(const std::string& text) {
  std::ostringstream out;out<<'"';
  for(unsigned char c:text) {
    if(c=='"'||c=='\\')out<<'\\'<<c;
    else if(c=='\n')out<<"\\n";
    else if(c=='\r')out<<"\\r";
    else if(c=='\t')out<<"\\t";
    else if(c<32)out<<"\\u"<<std::hex<<std::setw(4)<<std::setfill('0')<<static_cast<unsigned>(c)<<std::dec;
    else out<<c;
  }
  return out.str()+'"';
}
H::Target Target(const std::string& value) {
  if(value=="Si")return H::Target::Silicon;
  if(value=="O16")return H::Target::Oxygen16;
  if(value=="O18")return H::Target::Oxygen18;
  if(value=="O_mean")return H::Target::PublishedOxygenMean;
  return static_cast<H::Target>(-1);
}
const char* Name(H::Status status) {
  switch(status) {
    case H::Status::Evaluated:return "Evaluated";
    case H::Status::Unavailable:return "Unavailable";
    case H::Status::InvalidInput:return "InvalidInput";
    case H::Status::NumericalError:return "NumericalError";
  }
  return "ERROR";
}
void Fit(const H::Fit& fit) {
  std::cout<<"{\"fit_id\":"<<Json(fit.id)<<",\"target\":"<<Json(H::TargetName(fit.target))
    <<",\"mechanism\":"<<Json(H::MechanismName(fit.mechanism))
    <<",\"energy_min_eV\":"<<fit.minimumIncidentEnergyEv<<",\"energy_max_eV\":"<<fit.maximumIncidentEnergyEv
    <<",\"minimum_included\":"<<(fit.minimumIncluded?"true":"false")<<",\"maximum_included\":"<<(fit.maximumIncluded?"true":"false")
    <<",\"vc_m_s\":"<<fit.characteristicVelocityMS<<",\"dvc_m_s\":"<<fit.reportedVelocityUncertaintyMS
    <<",\"a\":"<<fit.correctedYieldIntercept<<",\"da\":"<<fit.reportedInterceptUncertainty
    <<",\"r_squared\":";
  if(fit.hasRSquared)std::cout<<fit.rSquared;else std::cout<<"null";
  std::cout<<",\"surface_ratio\":";
  if(fit.hasSurfaceOxygenSiliconRatio)std::cout<<fit.surfaceOxygenSiliconRatio;else std::cout<<"null";
  std::cout<<",\"surface_ratio_uncertainty\":";
  if(fit.hasSurfaceOxygenSiliconRatio)std::cout<<fit.reportedSurfaceRatioUncertainty;else std::cout<<"null";
  std::cout<<'}';
}
int main(int argc,char** argv) {
  std::cout<<std::setprecision(17);
  try {
    if(argc<4)throw std::runtime_error("usage: probe publication.pdf table.csv mode [args]");
    H::Registry registry;std::string error;auto source=H::FrozenV1Sources(argv[1],argv[2]);const std::string mode=argv[3];
    if(mode=="caller-hash")source.publicationSha256=std::string(64,'0');
    if(mode!="unloaded" && !registry.LoadPublishedTable(source,&error)) {
      std::cout<<"{\"status\":\"ERROR\",\"error\":"<<Json(error)<<"}\n";return 2;
    }
    if(mode=="catalog") {
      std::cout<<"{\"status\":\"Evaluated\",\"fits\":[";bool first=true;
      for(auto target:{H::Target::Silicon,H::Target::Oxygen16,H::Target::Oxygen18,H::Target::PublishedOxygenMean})
        for(auto mechanism:{H::Mechanism::CollisionInducedNeutralization,H::Mechanism::AugerNeutralization}) {
          H::Fit fit;const auto status=registry.GetFit(target,mechanism,&fit,&error);
          if(status==H::Status::Unavailable)continue;
          if(status!=H::Status::Evaluated)throw std::runtime_error(error);
          if(!first)std::cout<<',';
          first=false;Fit(fit);
        }
      std::cout<<"]}\n";return 0;
    }
    if(mode=="reload") {
      if(argc!=6)throw std::runtime_error("reload needs changed publication and table paths");
      const bool accepted=registry.LoadPublishedTable(H::FrozenV1Sources(argv[4],argv[5]),&error);
      H::Fit fit;const auto status=registry.SelectFit(H::Target::Silicon,3000,&fit,&error);
      const bool preserved=!accepted && registry.loaded() && registry.sources().publicationPath==argv[1] &&
        registry.sources().processedTablePath==argv[2] && status==H::Status::Evaluated && fit.characteristicVelocityMS==291000;
      std::cout<<"{\"status\":"<<Json(preserved?"Evaluated":"FAIL")<<",\"previous_registry_preserved\":"<<(preserved?"true":"false")<<"}\n";
      return preserved?0:1;
    }
    if(mode=="missing-fit") {
      H::Fit fit;fit.id="SENTINEL";const auto status=registry.GetFit(H::Target::Silicon,H::Mechanism::AugerNeutralization,&fit,&error);
      std::cout<<"{\"status\":"<<Json(Name(status))<<",\"unchanged\":"<<(fit.id=="SENTINEL"?"true":"false")<<",\"error\":"<<Json(error)<<"}\n";
      return status==H::Status::Unavailable?77:2;
    }
    H::Result result;result.logChargedFraction=12345;result.fit.id="SENTINEL";
    H::Status status;
    if(mode=="unloaded")status=registry.EvaluatePublishedCoordinate(H::Target::Silicon,3000,1e-5,&result,&error);
    else {
      if(argc!=7)throw std::runtime_error("evaluate/velocities requires target,energy_eV,inverse_velocity_s_per_m");
      const auto target=Target(argv[4]);const double energy=std::stod(argv[5]);
      if(mode=="velocities")status=registry.EvaluateFromVelocityVectors(target,energy,{{0,0,-100000}},{{50000,0,80000}},{{0,0,1}},&result,&error);
      else status=registry.EvaluatePublishedCoordinate(target,energy,std::stod(argv[6]),&result,&error);
    }
    std::cout<<"{\"status\":"<<Json(Name(status))<<",\"unchanged\":"<<((result.logChargedFraction==12345 && result.fit.id=="SENTINEL")?"true":"false")
      <<",\"physical_velocity_mapping_available\":"<<(result.physicalVelocityMappingAvailable?"true":"false")<<",\"error\":"<<Json(error);
    if(status==H::Status::Evaluated) {
      std::cout<<",\"fit\":";Fit(result.fit);
      std::cout<<",\"log_probability\":"<<result.logChargedFraction<<",\"probability\":"<<result.chargedFraction
        <<",\"probability_at_vc_plus_uncertainty\":"<<result.fractionAtVcPlusReportedUncertainty
        <<",\"probability_at_vc_minus_uncertainty\":"<<result.fractionAtVcMinusReportedUncertainty;
    }
    std::cout<<"}\n";
    return status==H::Status::Evaluated?0:status==H::Status::Unavailable?77:2;
  }
  catch(const std::exception& exception) {
    std::cout<<"{\"status\":\"ERROR\",\"error\":"<<Json(exception.what())<<"}\n";return 2;
  }
}
