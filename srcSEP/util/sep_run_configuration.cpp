#include "sep_run_configuration.h"

#include "sep_configuration_matrix.h"

#include <cmath>
#include <iomanip>
#include <sstream>
#include <stdexcept>

namespace SEP {
namespace Run {
namespace {

FrozenConfiguration gActiveConfiguration;
bool gActiveConfigurationInstalled=false;

Transport::Status Error(const std::string& text) {
  return Transport::Status::Error(Transport::StatusCode::InvalidArgument,text);
}

template <class T>
LayeredValue<T> Select(const LayeredValue<T>& base,
                       const LayeredValue<T>& overlay) {
  return static_cast<int>(overlay.source)>=static_cast<int>(base.source)
      ? overlay : base;
}

std::uint64_t Fnv(const std::string& text) {
  std::uint64_t value=UINT64_C(14695981039346656037);
  for (std::size_t i=0;i<text.size();++i) {
    value^=static_cast<unsigned char>(text[i]);
    value*=UINT64_C(1099511628211);
  }
  return value;
}

std::string Canonical(const Configuration& c) {
  std::ostringstream out;
  out<<std::setprecision(17)
      <<static_cast<int>(c.mover.value)<<' '<<static_cast<int>(c.mover.source)<<' '
      <<static_cast<int>(c.shockModel.value)<<' '<<static_cast<int>(c.shockModel.source)<<' '
      <<static_cast<int>(c.scenario.value)<<' '<<static_cast<int>(c.scenario.source)<<' '
      <<c.totalIterations.value<<' '<<static_cast<int>(c.totalIterations.source)<<' '
      <<c.fieldLineSeedAreaM2.value<<' '<<static_cast<int>(c.fieldLineSeedAreaM2.source)<<' '
      <<c.shockTurbulenceEfficiency.value<<' '
      <<static_cast<int>(c.shockTurbulenceEfficiency.source)<<' '
      <<c.shockTurbulencePlusFraction.value<<' '
      <<static_cast<int>(c.shockTurbulencePlusFraction.source)<<' '
      <<c.mergeMinimum.value<<' '<<static_cast<int>(c.mergeMinimum.source)<<' '
      <<c.mergeMaximum.value<<' '<<static_cast<int>(c.mergeMaximum.source)<<' '
      <<c.populationControl.spatialBins<<' '
      <<c.populationControl.momentumBins<<' '
      <<c.populationControl.pitchBins<<' '
      <<c.populationControl.deterministicSelection<<' '
      <<c.populationControl.relativeInvariantTolerance<<' '
      <<c.populationControl.lineageSchema<<' '
      <<static_cast<int>(c.invalidSamplingPolicy.value)<<' '
      <<static_cast<int>(c.invalidSamplingPolicy.source)<<' '
      <<static_cast<int>(c.turbulence.source)<<' '
      <<static_cast<int>(c.turbulence.representation)<<' '
      <<static_cast<int>(c.turbulence.coupling)<<' '
      <<c.turbulence.advectionEnabled<<' '<<c.turbulence.reflectionEnabled<<' '
      <<c.turbulence.cascadeEnabled<<' '<<c.turbulence.shockInjectionEnabled<<' '
      <<c.turbulence.cflSafety<<' '
      <<static_cast<int>(c.turbulence.operatorSplitting)<<' '
      <<c.turbulence.operatorAccuracySafety<<' '
      <<c.turbulence.maximumSourceFraction<<' '
      <<c.turbulence.maximumCascadeFraction<<' '
      <<c.turbulence.minimumSubstepS<<' '
      <<c.turbulence.maximumSubsteps<<' '
      <<c.turbulence.reflectionCoefficient<<' '
      <<c.turbulence.cascadeCoefficient<<' '
      <<c.numericalTolerances.geometryFraction<<' '
      <<c.numericalTolerances.deterministicRelativeTolerance<<' '
      <<c.injection.campaignSeed<<' '<<c.injection.macroparticlesPerEvent<<' '
      <<c.injection.injectionEfficiency<<' '
      <<static_cast<int>(c.injection.spectrum.measure)<<' '
      <<c.injection.spectrum.minimum<<' '<<c.injection.spectrum.maximum<<' '
      <<c.injection.spectrum.powerIndex<<' '
      <<static_cast<int>(c.injection.angular)<<'\n';
  return out.str();
}

}  // namespace

Configuration Defaults() {
  Configuration c;
  c.mover.value=Mover::ProductionMover::FocusedTransportDiffusion;
  c.shockModel.value=ShockModel::Swcme1d;
  c.scenario.value=CmeScenario::Fast;
  c.totalIterations.value=UINT64_C(100000001);
  // This names the old pi-square-metre normalization explicitly.  Applications
  // should override it with their seed-surface partition; it is fingerprinted
  // and converted once to Phi_i=A_seed |B_seed| by WP26.
  c.fieldLineSeedAreaM2.value=3.14159265358979323846;
  c.shockTurbulenceEfficiency.value=0.02;
  c.shockTurbulencePlusFraction.value=0.5;
  c.mergeMinimum.value=600;
  c.mergeMaximum.value=1000;
  c.populationControl.minimumParticlesPerCell=c.mergeMinimum.value;
  c.populationControl.maximumParticlesPerCell=c.mergeMaximum.value;
  c.invalidSamplingPolicy.value=SamplingCore::InvalidParticlePolicy::Exclude;
  c.injection.campaignSeed=0;
  c.injection.macroparticlesPerEvent=300;
  c.injection.injectionEfficiency=1.0;
  c.injection.spectrum.measure=Injection::Measure::Momentum;
  c.injection.spectrum.minimum=1.0;
  c.injection.spectrum.maximum=10.0;
  c.injection.spectrum.powerIndex=4.0;
  return c;
}

Transport::Status Validate(const Configuration& c) {
  if (c.totalIterations.value==0 || !std::isfinite(c.fieldLineSeedAreaM2.value) ||
      c.fieldLineSeedAreaM2.value<=0.0 ||
      !std::isfinite(c.shockTurbulenceEfficiency.value) ||
      c.shockTurbulenceEfficiency.value<0.0 ||
      c.shockTurbulenceEfficiency.value>1.0 ||
      !std::isfinite(c.shockTurbulencePlusFraction.value) ||
      c.shockTurbulencePlusFraction.value<0.0 ||
      c.shockTurbulencePlusFraction.value>1.0 || c.mergeMinimum.value<0 ||
      c.mergeMaximum.value<c.mergeMinimum.value)
    return Error("run configuration has invalid iteration, flux seed, or merge limits");
  if (c.populationControl.minimumParticlesPerCell!=c.mergeMinimum.value ||
      c.populationControl.maximumParticlesPerCell!=c.mergeMaximum.value)
    return Error("population-control limits disagree with CLI provenance aliases");
  Transport::Status status=PopulationControl::ValidateConfiguration(
      c.populationControl);
  if (!status.ok()) return status;
  status=Transport::Coefficient::ValidateConfiguration(c.coefficients);
  if (!status.ok()) return status;
  status=Transport::ValidateNumericalTolerances(c.numericalTolerances);
  if (!status.ok()) return status;
  status=Turbulence::ValidateConfiguration(c.turbulence);
  if (!status.ok()) return status;
  ConfigurationMatrix::Combination combination;
  combination.mover=c.mover.value;
  combination.coefficientSource=c.coefficients.source;
  combination.turbulenceSource=c.turbulence.source;
  combination.coupling=c.turbulence.coupling;
  status=ConfigurationMatrix::Preflight(combination,false);
  if (!status.ok()) return status;
  return Injection::Validate(c.injection);
}

Configuration Merge(const Configuration& base,const Configuration& overlay) {
  Configuration c=overlay;
  c.mover=Select(base.mover,overlay.mover);
  c.shockModel=Select(base.shockModel,overlay.shockModel);
  c.scenario=Select(base.scenario,overlay.scenario);
  c.totalIterations=Select(base.totalIterations,overlay.totalIterations);
  c.fieldLineSeedAreaM2=Select(base.fieldLineSeedAreaM2,overlay.fieldLineSeedAreaM2);
  c.shockTurbulenceEfficiency=Select(base.shockTurbulenceEfficiency,
                                     overlay.shockTurbulenceEfficiency);
  c.shockTurbulencePlusFraction=Select(base.shockTurbulencePlusFraction,
                                       overlay.shockTurbulencePlusFraction);
  c.mergeMinimum=Select(base.mergeMinimum,overlay.mergeMinimum);
  c.mergeMaximum=Select(base.mergeMaximum,overlay.mergeMaximum);
  c.populationControl.minimumParticlesPerCell=c.mergeMinimum.value;
  c.populationControl.maximumParticlesPerCell=c.mergeMaximum.value;
  c.invalidSamplingPolicy=Select(base.invalidSamplingPolicy,overlay.invalidSamplingPolicy);
  return c;
}

std::string Fingerprint(const Configuration& c) {
  std::ostringstream out;
  out<<std::hex<<std::setw(16)<<std::setfill('0')<<Fnv(Canonical(c));
  return out.str();
}

Transport::Status FrozenConfiguration::Create(const Configuration& c,
                                              FrozenConfiguration* frozen) {
  if (!frozen) return Error("frozen configuration destination is null");
  const Transport::Status status=Validate(c);
  if (!status.ok()) return status;
  frozen->configuration_=c;
  frozen->fingerprint_=Fingerprint(c);
  return Transport::Status::Ok();
}

Transport::Status Serialize(const FrozenConfiguration& c,std::string* text) {
  if (!text) return Error("run configuration checkpoint output is null");
  *text="SEP_RUN_CONFIGURATION 2\n"+c.fingerprint()+"\n"+Canonical(c.get());
  return Transport::Status::Ok();
}

Transport::Status Deserialize(const std::string& text,FrozenConfiguration* result) {
  if (!result) return Error("run configuration checkpoint destination is null");
  std::istringstream in(text); std::string magic,fingerprint; int version=0;
  Configuration c=Defaults(); int mover,shock,scenario,sampling;
  int moverSource,shockSource,scenarioSource,totalSource,areaSource,
      shockEfficiencySource,shockSplitSource,minSource,
      maxSource,samplingSource,turbSource,turbRep,turbCoupling,measure,angular,
      turbulenceSplitting;
  if (!(in>>magic>>version>>fingerprint) || magic!="SEP_RUN_CONFIGURATION" || version!=2 ||
      !(in>>mover>>moverSource>>shock>>shockSource>>scenario>>scenarioSource
          >>c.totalIterations.value>>totalSource>>c.fieldLineSeedAreaM2.value>>areaSource
          >>c.shockTurbulenceEfficiency.value>>shockEfficiencySource
          >>c.shockTurbulencePlusFraction.value>>shockSplitSource
          >>c.mergeMinimum.value>>minSource>>c.mergeMaximum.value>>maxSource
          >>c.populationControl.spatialBins
          >>c.populationControl.momentumBins
          >>c.populationControl.pitchBins
          >>c.populationControl.deterministicSelection
          >>c.populationControl.relativeInvariantTolerance
          >>c.populationControl.lineageSchema
          >>sampling>>samplingSource>>turbSource>>turbRep>>turbCoupling
          >>c.turbulence.advectionEnabled>>c.turbulence.reflectionEnabled
          >>c.turbulence.cascadeEnabled>>c.turbulence.shockInjectionEnabled
          >>c.turbulence.cflSafety>>turbulenceSplitting
          >>c.turbulence.operatorAccuracySafety
          >>c.turbulence.maximumSourceFraction
          >>c.turbulence.maximumCascadeFraction
          >>c.turbulence.minimumSubstepS
          >>c.turbulence.maximumSubsteps
          >>c.turbulence.reflectionCoefficient
          >>c.turbulence.cascadeCoefficient
          >>c.numericalTolerances.geometryFraction
          >>c.numericalTolerances.deterministicRelativeTolerance
          >>c.injection.campaignSeed>>c.injection.macroparticlesPerEvent
          >>c.injection.injectionEfficiency>>measure
          >>c.injection.spectrum.minimum>>c.injection.spectrum.maximum
          >>c.injection.spectrum.powerIndex>>angular))
    return Error("run configuration checkpoint is invalid or truncated");
  c.mover.value=static_cast<Mover::ProductionMover>(mover);
  c.mover.source=static_cast<ValueSource>(moverSource);
  c.shockModel.value=static_cast<ShockModel>(shock);
  c.shockModel.source=static_cast<ValueSource>(shockSource);
  c.scenario.value=static_cast<CmeScenario>(scenario);
  c.scenario.source=static_cast<ValueSource>(scenarioSource);
  c.totalIterations.source=static_cast<ValueSource>(totalSource);
  c.fieldLineSeedAreaM2.source=static_cast<ValueSource>(areaSource);
  c.shockTurbulenceEfficiency.source=static_cast<ValueSource>(shockEfficiencySource);
  c.shockTurbulencePlusFraction.source=static_cast<ValueSource>(shockSplitSource);
  c.mergeMinimum.source=static_cast<ValueSource>(minSource);
  c.mergeMaximum.source=static_cast<ValueSource>(maxSource);
  c.populationControl.minimumParticlesPerCell=c.mergeMinimum.value;
  c.populationControl.maximumParticlesPerCell=c.mergeMaximum.value;
  c.invalidSamplingPolicy.value=static_cast<SamplingCore::InvalidParticlePolicy>(sampling);
  c.invalidSamplingPolicy.source=static_cast<ValueSource>(samplingSource);
  c.turbulence.source=static_cast<Turbulence::Source>(turbSource);
  c.turbulence.representation=static_cast<Turbulence::Representation>(turbRep);
  c.turbulence.coupling=static_cast<Turbulence::CouplingPolicy>(turbCoupling);
  c.turbulence.operatorSplitting=
      static_cast<Turbulence::OperatorSplitting>(turbulenceSplitting);
  c.injection.spectrum.measure=static_cast<Injection::Measure>(measure);
  c.injection.angular=static_cast<Injection::AngularDistribution>(angular);
  Transport::Status status=FrozenConfiguration::Create(c,result);
  if (!status.ok()) return status;
  return result->fingerprint()==fingerprint ? Transport::Status::Ok()
      : Error("run configuration fingerprint does not match checkpoint");
}

Transport::Status VerifyRestartCompatibility(const FrozenConfiguration& requested,
                                             const FrozenConfiguration& checkpoint) {
  return requested.fingerprint()==checkpoint.fingerprint()
      ? Transport::Status::Ok()
      : Transport::Status::Error(Transport::StatusCode::UnsupportedConfiguration,
                                 "restart run-configuration fingerprint mismatch");
}

Transport::Status InstallActive(const Configuration& configuration) {
  FrozenConfiguration candidate;
  const Transport::Status status=FrozenConfiguration::Create(configuration,&candidate);
  if (!status.ok()) return status;
  gActiveConfiguration=candidate;
  gActiveConfigurationInstalled=true;
  return Transport::Status::Ok();
}

const FrozenConfiguration& Active() {
  if (!gActiveConfigurationInstalled) {
    const Configuration defaults=Defaults();
    const Transport::Status status=FrozenConfiguration::Create(
        defaults,&gActiveConfiguration);
    if (!status.ok()) throw std::runtime_error(status.message);
    gActiveConfigurationInstalled=true;
  }
  return gActiveConfiguration;
}

const char* ValueSourceName(ValueSource source) {
  switch (source) {
    case ValueSource::Default: return "default";
    case ValueSource::InputFile: return "input-file";
    case ValueSource::CommandLine: return "command-line";
  }
  return "unknown";
}

}  // namespace Run
}  // namespace SEP
