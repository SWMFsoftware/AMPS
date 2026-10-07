#include "particle_normalization.h"

#include "sep_coronal_cme/constants.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>
#include <utility>

namespace SEP3D { namespace RuntimeModel { namespace {

Core::Status Invalid(const std::string& message) {
  return Core::Status(Core::StatusCode::InvalidInput,message);
}

bool Near(double value,double reference,double relative=1.0e-5) {
  return std::isfinite(value)&&std::isfinite(reference)&&reference!=0&&
      std::abs(value-reference)<=relative*std::abs(reference);
}

Core::Status SpeciesRate(
    const CompiledSpeciesRecord& species,
    const SEP::CoronaSwcme::ShockFront::IncidentParticleFlux& flux,
    double* rate) {
  using SEP::CoronalCME::Constants::kAlphaMassKg;
  using SEP::CoronalCME::Constants::kElectronMassKg;
  using SEP::CoronalCME::Constants::kElementaryChargeC;
  using SEP::CoronalCME::Constants::kProtonMassKg;
  if(rate==nullptr)return Invalid("species-rate output is null");

  // Mass and signed charge, rather than a site-dependent chemical spelling,
  // identify the three populations present in the selected ambient EOS.
  // Tight relative tolerances accommodate harmless constant-table rounding
  // but do not merge isotopes or charge states.
  if(Near(species.massKg,kElectronMassKg)&&
      Near(species.chargeC,-kElementaryChargeC))
    *rate=flux.electronRatePerS;
  else if(Near(species.massKg,kProtonMassKg)&&
      Near(species.chargeC,kElementaryChargeC))
    *rate=flux.protonRatePerS;
  else if(Near(species.massKg,kAlphaMassKg)&&
      Near(species.chargeC,2*kElementaryChargeC))
    *rate=flux.alphaRatePerS;
  else {
    std::ostringstream message;
    message << "compiled species '" << species.symbol << "' (mass="
            << species.massKg << " kg, charge=" << species.chargeC
            << " C) has no population in the reduced ambient source model";
    return Core::Status(Core::StatusCode::ConfigurationConflict,message.str());
  }
  if(!std::isfinite(*rate)||*rate<=0)
    return Invalid("compiled species '"+species.symbol+
        "' has a zero or invalid incident source rate");
  return Core::Status::OK();
}

} // namespace

Core::Status CalculateMeshGlobalTimeStep(double minimumCellSize,
    double maximumSpeed,double margin,double* timeStep) {
  if(timeStep==nullptr)return Invalid("global-time-step output is null");
  if(!std::isfinite(minimumCellSize)||minimumCellSize<=0||
      !std::isfinite(maximumSpeed)||maximumSpeed<=0||
      maximumSpeed>Core::Const::c||!std::isfinite(margin)||margin<=0||margin>1)
    return Invalid("global time step requires h_min>0, 0<v_max<=c, and margin in (0,1]");
  const double candidate=margin*minimumCellSize/maximumSpeed;
  if(!std::isfinite(candidate)||candidate<=0)
    return Invalid("global time-step calculation did not produce a positive finite value");
  *timeStep=candidate;
  return Core::Status::OK();
}

Core::Status AlignObserverCadencesToGlobalStep(
    double timeStep,std::vector<ObserverOptions>* observers) {
  if(observers==nullptr||!std::isfinite(timeStep)||timeStep<=0)
    return Invalid("observer-cadence alignment requires dt>0 and an output vector");
  std::vector<ObserverOptions> candidate=*observers;
  for(ObserverOptions& observer:candidate) {
    if(!std::isfinite(observer.cadenceS)||observer.cadenceS<=0)
      return Invalid("observer '"+observer.id+"' has an invalid requested cadence");
    const double ratio=observer.cadenceS/timeStep;
    if(!std::isfinite(ratio)||ratio>
        static_cast<double>(std::numeric_limits<std::uint64_t>::max()))
      return Invalid("observer '"+observer.id+"' cadence tick count overflows");
    const double nearest=std::round(ratio);
    const double tolerance=64*std::numeric_limits<double>::epsilon()*
        std::max(1.0,std::abs(ratio));
    const double ticks=std::abs(ratio-nearest)<=tolerance?
        std::max(1.0,nearest):std::max(1.0,std::ceil(ratio));
    const double aligned=ticks*timeStep;
    if(!std::isfinite(aligned)||aligned+64*
        std::numeric_limits<double>::epsilon()*
        std::max(1.0,std::abs(observer.cadenceS))<observer.cadenceS)
      return Invalid("observer '"+observer.id+
          "' cadence could not be aligned without increasing output frequency");
    observer.cadenceS=aligned;
  }
  *observers=std::move(candidate);
  return Core::Status::OK();
}

Core::Status CalculateSpeciesParticleNormalizations(
    const std::vector<CompiledSpeciesRecord>& species,
    const SEP::CoronaSwcme::ShockFront::IncidentParticleFlux& flux,
    double timeStep,std::uint64_t particlesPerIteration,
    std::vector<SpeciesParticleNormalization>* normalizations) {
  if(normalizations==nullptr)
    return Invalid("species-normalization output is null");
  if(species.empty()||!std::isfinite(timeStep)||timeStep<=0||
      particlesPerIteration==0)
    return Invalid("species normalization requires species, dt>0, and a positive model-particle count");

  std::vector<SpeciesParticleNormalization> candidate;
  candidate.reserve(species.size());
  for(const CompiledSpeciesRecord& compiled:species) {
    double rate=0;
    Core::Status status=SpeciesRate(compiled,flux,&rate);
    if(!status.ok())return status;
    const long double weight=static_cast<long double>(rate)*timeStep/
        static_cast<long double>(particlesPerIteration);
    if(!(weight>0)||!std::isfinite(static_cast<double>(weight)))
      return Invalid("particle weight is non-positive or overflows for species '"+
          compiled.symbol+"'");
    SpeciesParticleNormalization item;
    item.ampsIndex=compiled.ampsIndex;
    item.symbol=compiled.symbol;
    item.physicalSourceRatePerS=rate;
    item.macroparticleWeight=static_cast<double>(weight);
    candidate.push_back(std::move(item));
  }
  *normalizations=std::move(candidate);
  return Core::Status::OK();
}

} } // namespace SEP3D::RuntimeModel
