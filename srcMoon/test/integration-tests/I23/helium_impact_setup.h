// Shared by the I23 MPI ensemble and linked PIC probes. Every physical value
// arrives from the Python runner, which reads it from authenticated sources:
// the NIST/CODATA state profile, the frozen D05 Voronov table, the approved
// admission record and the two frozen admission evidence files. This header
// only assembles the production HeliumElectronImpact request objects; it
// contains no rate, mass, threshold or tolerance of its own.
//
// Arguments are `key=value` tokens. `source=id|path|sha256|role` repeats.
#ifndef AMPS_I23_HELIUM_IMPACT_SETUP_H
#define AMPS_I23_HELIUM_IMPACT_SETUP_H
#include "srcMoon/species/HeliumElectronImpact.h"
#include <cmath>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace I23Setup {
namespace H=Moon::Species::HeliumElectronImpact;
namespace P=Moon::Species::PhysicalStateProfile;
namespace R=Moon::Species::ReactionExecution;
namespace B=::AMPS::ChemicalDatabase::ReservoirProducts;

inline void Require(bool value,const std::string& message) {
  if (!value) throw std::runtime_error(message);
}
inline void Require(const H::Status& status) {Require(status.ok(),status.message);}

struct Arguments {
  std::map<std::string,std::string> values;
  std::vector<P::SourceFile> sources;
  const std::string& operator[](const std::string& key) const {
    const auto found=values.find(key);
    Require(found!=values.end(),"missing I23 probe argument "+key);
    return found->second;
  }
  long double Real(const std::string& key) const {
    std::size_t used=0;const auto value=std::stold((*this)[key],&used);
    Require(used==(*this)[key].size() && std::isfinite(value),"nonnumeric I23 argument "+key);
    return value;
  }
};

inline Arguments Parse(int argc,char** argv,int first=1) {
  Arguments output;
  for (int i=first;i<argc;++i) {
    const std::string token=argv[i];
    const auto equal=token.find('=');
    Require(equal!=std::string::npos && equal>0,"I23 probe argument must be key=value: "+token);
    const auto key=token.substr(0,equal),value=token.substr(equal+1);
    if (key=="source") {
      std::vector<std::string> parts;std::stringstream stream(value);std::string part;
      while (std::getline(stream,part,'|')) parts.push_back(part);
      Require(parts.size()==4,"source argument must be id|path|sha256|role");
      output.sources.push_back({parts[0],parts[1],parts[2],parts[3]});
    }
    else Require(output.values.emplace(key,value).second,"duplicate I23 argument "+key);
  }
  return output;
}

inline void LoadProfile(const Arguments& a,P::Registry* profile) {
  P::Files files;
  files.isotopes={a["isotopes"],a["isotopes_sha256"]};
  files.constants={a["constants"],a["constants_sha256"]};
  files.states={a["states"],a["states_sha256"]};
  files.sources=a.sources;
  std::string error;
  Require(profile->Load(files,&error),"state profile refused: "+error);
}

// Returns Unavailable (missing source) distinctly so the caller can SKIP.
inline H::Status Prepare(const Arguments& a,const P::Registry& profile,int he,int ion,H::Provider* provider) {
  return H::Provider::Prepare(profile,a["d05"],a["d05_sha256"],he,ion,provider);
}

inline H::Bath Bath(const Arguments& a,long double density,const std::string& frame) {
  H::Bath bath;
  bath.thermalEnergyEv=static_cast<double>(a.Real("thermal_energy_ev"));
  bath.densityM3=static_cast<double>(density);
  bath.frameIdentity=frame;
  return bath; // zero drift and zero target velocity: the approved scope
}

// Assemble the reviewed admission. The approved observable scope, tolerances
// and evidence identities come from the frozen admission files via the runner.
inline H::Admission Admission(const Arguments& a,const P::Registry& profile,const std::string& frame) {
  H::Admission admission;
  auto& e=admission.evidence;
  e.modelIdentity=a["model_id"];
  e.sourceProfileIdentity=profile.identity();e.sourceProfileSha256=profile.identity();
  e.stateCrosswalkIdentity="srcMoon-I23-helium-state-crosswalk-v1";e.stateCrosswalkSha256=a["crosswalk_sha256"];
  e.domainEvidenceIdentity=a["approval_identity"];
  e.modelErrorEvidenceIdentity="srcMoon-I23-helium-model-error-v1";e.modelErrorEvidenceSha256=a["model_error_sha256"];
  e.frameIdentity=frame;e.timeScaleIdentity=a["time_scale"];
  e.energyReferenceIdentity="P02:separated-neutral-ground-isotope-atoms-and-free-electron-reference-zero";
  e.observableIdentity=a["observable"];
  e.maximumHeavyVelocityErrorMS=a.Real("max_heavy_velocity_error_m_s");
  e.maximumHeavyEnergyErrorJ=a.Real("max_heavy_energy_error_j");
  e.maximumRelativeMassError=a.Real("max_relative_mass_error");
  e.maximumClassicalHeavySpeedMS=a.Real("max_classical_heavy_speed_m_s");
  auto& n=admission.numerical;
  n.massAbsoluteKg=a.Real("reservoir_mass_absolute_kg");n.massRelative=a.Real("reservoir_mass_relative");
  n.timeProjectionAbsoluteS=a.Real("reservoir_time_projection_absolute_s");
  n.conservation={a.Real("reservoir_momentum_absolute_kg_m_s"),a.Real("reservoir_momentum_relative"),
                  a.Real("reservoir_energy_absolute_j"),a.Real("reservoir_energy_relative"),
                  a.Real("reservoir_unit_vector_tolerance")};
  auto& planner=admission.planner;
  planner.qualified=true;
  planner.identity=a["kinematic_identity"];
  planner.evidenceIdentity=a["crosswalk_sha256"]+";"+a["model_error_sha256"];
  planner.eventFrame=frame;
  planner.maximumClassicalSpeedMS=static_cast<double>(a.Real("max_classical_heavy_speed_m_s"));
  planner.bounds={{static_cast<double>(a.Real("ledger_mass_absolute")),static_cast<double>(a.Real("ledger_mass_relative"))},
                  {static_cast<double>(a.Real("ledger_momentum_absolute")),static_cast<double>(a.Real("ledger_momentum_relative"))},
                  {static_cast<double>(a.Real("ledger_energy_absolute")),static_cast<double>(a.Real("ledger_energy_relative"))}};
  admission.crosswalkPath=a["crosswalk"];admission.modelErrorPath=a["model_error"];
  admission.thresholdCrosswalkMaximumDifferenceJ=a.Real("threshold_crosswalk_max_difference_j");
  admission.maximumAbsoluteLogSlopeError=a.Real("max_abs_log_slope_error");
  return admission;
}
}
#endif
