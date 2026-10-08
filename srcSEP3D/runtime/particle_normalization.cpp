#include "particle_normalization.h"

#include "sep_coronal_cme/constants.h"
#include "sep_injection_spectrum.h"
#include "../transport/keyed_random.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>
#include <unordered_set>
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

enum class AmbientPopulation { Electron, Proton, Alpha };

Core::Status IdentifyPopulation(const CompiledSpeciesRecord& species,
                                AmbientPopulation* population) {
  using SEP::CoronalCME::Constants::kAlphaMassKg;
  using SEP::CoronalCME::Constants::kElectronMassKg;
  using SEP::CoronalCME::Constants::kElementaryChargeC;
  using SEP::CoronalCME::Constants::kProtonMassKg;
  if(population==nullptr)return Invalid("ambient-population output is null");
  if(Near(species.massKg,kElectronMassKg)&&
      Near(species.chargeC,-kElementaryChargeC))
    *population=AmbientPopulation::Electron;
  else if(Near(species.massKg,kProtonMassKg)&&
      Near(species.chargeC,kElementaryChargeC))
    *population=AmbientPopulation::Proton;
  else if(Near(species.massKg,kAlphaMassKg)&&
      Near(species.chargeC,2*kElementaryChargeC))
    *population=AmbientPopulation::Alpha;
  else return Core::Status(Core::StatusCode::ConfigurationConflict,
      "compiled species '"+species.symbol+
      "' has no population in the reduced ambient source model");
  return Core::Status::OK();
}

double NumberDensity(const SEP::CoronaSwcme::ShockFront::Provider& provider,
    const SEP::CoronaSwcme::ShockFront::ShockRecord& record,
    AmbientPopulation population) {
  if(population==AmbientPopulation::Electron)
    return record.upstream.plasma.electronNumberDensityM3;
  if(population==AmbientPopulation::Proton)
    return record.upstream.plasma.protonNumberDensityM3;
  return provider.Event().ambient.composition.alphaToProtonNumberRatio*
      record.upstream.plasma.protonNumberDensityM3;
}

double MomentumFromKineticEnergy(double kineticEnergyJ,double massKg) {
  const long double c=static_cast<long double>(Core::Const::c);
  const long double k=static_cast<long double>(kineticEnergyJ);
  const long double m=static_cast<long double>(massKg);
  return static_cast<double>(std::sqrt(k*(k+2*m*c*c))/c);
}

Core::Vec3 AsCore(const SEP::CoronalCME::Vec3& value) {
  return Core::Vec3(value.x,value.y,value.z);
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

Core::Status BuildSurfaceParticleRateDistribution(
    const SEP::CoronaSwcme::ShockFront::Provider& provider,
    const SEP::CoronaSwcme::ShockFront::Epoch& epoch,
    const CompiledSpeciesRecord& species,
    SurfaceParticleRateDistribution* distribution) {
  namespace SF=SEP::CoronaSwcme::ShockFront;
  if(distribution==nullptr)
    return Invalid("surface particle-rate output is null");
  if(epoch.triangles.size()!=epoch.records.size())
    return Invalid("reduced-front triangles and shock records are not one-to-one");
  if(epoch.area.numericalFailureM2>0)
    return Core::Status(Core::StatusCode::BackgroundInvalid,
        "surface source contains numerically unresolved physical area");
  AmbientPopulation population=AmbientPopulation::Electron;
  Core::Status status=IdentifyPopulation(species,&population);
  if(!status.ok())return status;

  SurfaceParticleRateDistribution candidate;
  candidate.generation=epoch.generation;
  candidate.epochS=epoch.trajectory.timeS;
  long double cumulative=0;
  long double acceptedArea=0;
  for(std::size_t index=0;index<epoch.records.size();++index) {
    const SF::ShockRecord& record=epoch.records[index];
    const SF::SurfaceTriangle& triangle=epoch.triangles[index];
    if(record.geometry.stableId!=triangle.stableId)
      return Invalid("reduced-front triangle/record stable identities disagree");
    if(record.status!=SF::FrontStatus::SolvedFastShock)continue;
    const double area=triangle.curvedAreaM2;
    const double density=NumberDensity(provider,record,population);
    if(!record.downstreamValid||!std::isfinite(area)||area<=0||
        !std::isfinite(record.geometry.areaM2)||
        std::abs(record.geometry.areaM2-area)>
            128*std::numeric_limits<double>::epsilon()*area||
        !std::isfinite(record.inflowMPerS)||record.inflowMPerS<=0||
        !std::isfinite(density)||density<0||
        !std::isfinite(record.jump.compressionRatio)||
        record.jump.compressionRatio<=1)
      return Invalid("accepted reduced-front face has invalid area, inflow, "
          "population density, or compression");
    if(density==0)
      return Core::Status(Core::StatusCode::ConfigurationConflict,
          "compiled species '"+species.symbol+
          "' has zero abundance in the selected ambient composition");
    const long double rate=static_cast<long double>(density)*
        record.inflowMPerS*area;
    cumulative+=rate;
    acceptedArea+=area;
    if(!std::isfinite(static_cast<double>(cumulative))||
        !(rate>0))
      return Invalid("accepted-face particle rate overflowed or is non-positive");
    SurfaceFaceParticleRate face;
    face.triangleIndex=index;
    face.stableId=triangle.stableId;
    face.physicalRatePerS=static_cast<double>(rate);
    face.cumulativeRatePerS=static_cast<double>(cumulative);
    face.compressionRatio=record.jump.compressionRatio;
    candidate.faces.push_back(face);
  }
  candidate.acceptedAreaM2=static_cast<double>(acceptedArea);
  candidate.physicalRatePerS=static_cast<double>(cumulative);
  if(!std::isfinite(candidate.physicalRatePerS)||
      candidate.physicalRatePerS<0)
    return Invalid("current reduced-front source rate is invalid");
  // A valid prescribed front can be sub-fast during startup. That is a
  // physical zero-rate Poisson process, not a provider failure and never a
  // reason to manufacture compression or a floor.
  *distribution=std::move(candidate);
  return Core::Status::OK();
}

Core::Status GenerateConstantWeightSurfaceInjectionBatch(
    const SEP::CoronaSwcme::ShockFront::Provider& provider,
    const SEP::CoronaSwcme::ShockFront::Epoch& epoch,
    const CompiledSpeciesRecord& species,const SourceOptions& source,
    double macroparticleWeight,double intervalS,std::uint64_t campaignSeed,
    std::uint64_t step,SurfaceInjectionBatch* batch) {
  if(batch==nullptr)return Invalid("surface injection-batch output is null");
  if(source.weightingModel!=SourceWeightingModel::ConstantStatisticalWeight)
    return Core::Status::Reserved("non-constant reduced-front particle weights");
  if(!std::isfinite(macroparticleWeight)||macroparticleWeight<=0||
      !std::isfinite(intervalS)||intervalS<=0||campaignSeed==0||
      source.maximumMacroparticlesPerSpeciesPerStep==0||
      !std::isfinite(source.minimumEnergyJ)||source.minimumEnergyJ<=0||
      !std::isfinite(source.maximumEnergyJ)||
      source.maximumEnergyJ<=source.minimumEnergyJ)
    return Invalid("constant-weight surface injection controls are invalid");

  SurfaceInjectionBatch candidate;
  Core::Status status=BuildSurfaceParticleRateDistribution(
      provider,epoch,species,&candidate.distribution);
  if(!status.ok())return status;
  const double macroRate=candidate.distribution.physicalRatePerS/
      macroparticleWeight;
  if(macroRate==0) {
    *batch=std::move(candidate);
    return Core::Status::OK();
  }
  if(!std::isfinite(macroRate)||macroRate<0)
    return Invalid("constant-weight Poisson rate is invalid");

  Transport::RandomKey waitingKey;
  waitingKey.campaignSeed=campaignSeed;
  waitingKey.particleId=static_cast<std::uint64_t>(species.ampsIndex+1);
  waitingKey.step=step;
  waitingKey.substep=epoch.generation;
  waitingKey.purpose=Transport::RandomPurpose::ReducedSourceWaitingTime;
  Transport::KeyedRandomStream waiting(waitingKey);
  std::unordered_set<std::uint64_t> identities;
  double eventTime=0;
  for(std::uint64_t ordinal=0;;++ordinal) {
    eventTime+=-std::log(waiting.UniformOpen01())/macroRate;
    if(!(eventTime<intervalS))break;
    if(candidate.events.size()>=
        source.maximumMacroparticlesPerSpeciesPerStep)
      return Core::Status(Core::StatusCode::Error,
          "reduced-front Poisson event count reached the declared fatal "
          "maximum; no source truncation was applied");

    Transport::RandomKey identityKey=waitingKey;
    identityKey.substep=ordinal;
    identityKey.purpose=Transport::RandomPurpose::ReducedSourceIdentity;
    std::uint64_t stableId=Transport::KeyedRandomStream::Hash(identityKey,0);
    if(stableId==0)stableId=1;
    if(!identities.insert(stableId).second)
      return Core::Status(Core::StatusCode::Error,
          "reduced-front source stable-ID collision");

    Transport::RandomKey eventKey=waitingKey;
    eventKey.particleId=stableId;
    eventKey.substep=ordinal;
    eventKey.purpose=Transport::RandomPurpose::ReducedSourceFace;
    Transport::KeyedRandomStream faceRandom(eventKey);
    const double selected=faceRandom.UniformOpen01()*
        candidate.distribution.physicalRatePerS;
    const auto faceIt=std::lower_bound(candidate.distribution.faces.begin(),
        candidate.distribution.faces.end(),selected,
        [](const SurfaceFaceParticleRate& face,double value) {
          return face.cumulativeRatePerS<value;
        });
    if(faceIt==candidate.distribution.faces.end())
      return Invalid("surface source cumulative face distribution is incomplete");
    const auto& triangle=epoch.triangles[faceIt->triangleIndex];
    for(std::uint32_t vertex:triangle.vertex)
      if(vertex>=epoch.vertices.size())
        return Invalid("surface source triangle references a missing vertex");

    eventKey.purpose=Transport::RandomPurpose::ReducedSourceBarycentric;
    Transport::KeyedRandomStream positionRandom(eventKey);
    // sqrt(u) produces constant area density on a planar triangle.  The
    // provider's exact curved area selects the patch; the chord triangle is
    // the maintained geometric representation on which AMPS owns a point.
    const double root=std::sqrt(positionRandom.UniformOpen01());
    const double second=positionRandom.UniformOpen01();
    const double w0=1-root,w1=root*(1-second),w2=root*second;
    const Core::Vec3 v0=AsCore(epoch.vertices[triangle.vertex[0]].positionM);
    const Core::Vec3 v1=AsCore(epoch.vertices[triangle.vertex[1]].positionM);
    const Core::Vec3 v2=AsCore(epoch.vertices[triangle.vertex[2]].positionM);

    const double q=source.spectrumModel==
        SourceSpectrumModel::LocalCompressionDsa
        ? 3*faceIt->compressionRatio/(faceIt->compressionRatio-1)
        : source.fixedPhaseSpacePowerIndex;
    if(!std::isfinite(q)||q<=2)
      return Invalid("surface source phase-space power index is invalid");
    SEP::Injection::Spectrum spectrum;
    spectrum.measure=SEP::Injection::Measure::Momentum;
    const double minimumMomentum=MomentumFromKineticEnergy(
        source.minimumEnergyJ,species.massKg);
    const double maximumMomentum=MomentumFromKineticEnergy(
        source.maximumEnergyJ,species.massKg);
    // Sample the dimensionless coordinate y=p/p_min. A direct power-law CDF
    // on SI momenta (~1e-23 kg m/s for energetic electrons) can overflow the
    // intermediate p^(1-a) even though the physical inverse is finite. A
    // constant coordinate scale leaves the normalized law exactly unchanged.
    spectrum.minimum=1.0;
    spectrum.maximum=maximumMomentum/minimumMomentum;
    spectrum.powerIndex=q-2;
    eventKey.purpose=Transport::RandomPurpose::ReducedSourceSpectrum;
    Transport::KeyedRandomStream spectrumRandom(eventKey);
    const auto sampled=SEP::Injection::InverseCdf(
        spectrum,spectrumRandom.UniformOpen01());
    if(!sampled.status.ok())
      return Invalid("surface source momentum inverse CDF failed: "+
          sampled.status.message);

    SurfaceInjectionEvent event;
    event.stableId=stableId;
    event.triangleStableId=faceIt->stableId;
    event.triangleIndex=faceIt->triangleIndex;
    event.eventTimeS=eventTime;
    event.remainingStepFraction=(intervalS-eventTime)/intervalS;
    event.momentumKgMPerS=minimumMomentum*sampled.value;
    event.compressionRatio=faceIt->compressionRatio;
    event.positionM=w0*v0+w1*v1+w2*v2;
    candidate.events.push_back(event);
  }
  *batch=std::move(candidate);
  return Core::Status::OK();
}

Core::Status GenerateLogUniformMomentumImportanceBatch(
    const SEP::CoronaSwcme::ShockFront::Provider&,
    const SEP::CoronaSwcme::ShockFront::Epoch&,
    const CompiledSpeciesRecord&,const SourceOptions&,double,double,
    std::uint64_t,std::uint64_t,SurfaceInjectionBatch*) {
  return Core::Status::Reserved(
      "log-uniform momentum importance weighting for reduced-front injection");
}

} } // namespace SEP3D::RuntimeModel
