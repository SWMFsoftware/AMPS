#include "sep_wave_normalization_validation.h"
#include "sep_turbulence_core.h"
#include "sep_focused_transport_mfp_core.h"
#include "sep_background_derivatives.h"
#include "sep_common_header_path.h"
#include SRCSEP_SEP_COMMON_HEADER(sep_coefficient_physics.h)

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace SEP {
namespace Testing {
namespace {

namespace CP = Transport::CoefficientPhysics;
const double kReferenceMu0 = 1.256637062e-6;  // AMPS constants.h SI input.
const double kPlusW = 4.0e-13;              // Magnetic energy [J/m^3].
const double kMinusW = 1.0e-13;
const double kRelativeTolerance = 128.0 * std::numeric_limits<double>::epsilon();

bool Relative(double value, double expected) {
  return std::isfinite(value) && (expected == 0.0 ? value == 0.0 :
      std::fabs(value - expected) <= kRelativeTolerance * std::fabs(expected));
}

Result Complete(bool ok, const char* message) {
  Result result;
  result.status = ok ? Status::Pass : Status::Fail;
  result.message = message;
  result.configuration.push_back("evidence_level=analytical-core");
  result.configuration.push_back("wave_energy_convention=magnetic-only");
  result.metrics.push_back({"assertion_failures", ok ? 0.0 : 1.0,
                            0.0, "<=", "count"});
  return result;
}

Turbulence::CellState Cell() {
  Turbulence::CellState cell;
  cell.lengthM = 10.0;
  cell.volumeM3 = 20.0;
  cell.ePlusJ = kPlusW * cell.volumeM3;
  cell.eMinusJ = kMinusW * cell.volumeM3;
  return cell;
}

Result IntegratedAmplitude() {
  const Turbulence::CellState cell = Cell();
  const Turbulence::DerivedCell derived = Turbulence::DeriveCell(cell);
  const CP::MagneticVariance coefficientInput =
      CP::MagneticVarianceFromEnergyDensity(kPlusW, kMinusW);
  const double expected = 2.0 * kReferenceMu0 * (kPlusW + kMinusW);
  bool ok = derived.status.ok() && coefficientInput.status.ok() &&
      Relative(derived.wPlusJPerM3, kPlusW) &&
      Relative(derived.wMinusJPerM3, kMinusW) &&
      Relative(derived.deltaB2T2, expected) &&
      Relative(coefficientInput.plusT2, 2.0 * kReferenceMu0 * kPlusW) &&
      Relative(coefficientInput.minusT2, 2.0 * kReferenceMu0 * kMinusW) &&
      derived.deltaB2T2 == coefficientInput.totalT2 &&
      Relative(derived.crossHelicity, 0.6);
  const double oldHalfAmplitude = kReferenceMu0 * (kPlusW + kMinusW);
  ok = ok && !Relative(oldHalfAmplitude, expected);
  Result result = Complete(ok,
      "integrated magnetic energy and coefficient view recover deltaB^2=2*mu0*W_B");
  result.metrics.push_back({"variance_relative_error",
      std::fabs(derived.deltaB2T2 / expected - 1.0),
      kRelativeTolerance, "<=", "dimensionless"});
  result.metrics.push_back({"old_conversion_relative_error",
      std::fabs(oldHalfAmplitude / expected - 1.0),
      kRelativeTolerance, ">", "dimensionless"});
  return result;
}

Result SourceCoefficientClosure() {
  CP::LocalInputView local;
  local.source = "self-consistent";
  local.representation = "integrated-magnetic-branch-energy";
  local.generation = 17;
  local.checksum = 19;
  local.heliocentricRadiusM = 1.495978707e11;
  local.magneticFieldT = 5.0e-9;
  local.alfvenSpeedMPerS = 0.0;
  const CP::MagneticVariance variance =
      CP::MagneticVarianceFromEnergyDensity(kPlusW, kMinusW);
  if (!variance.status.ok()) return Complete(false, "valid magnetic input was rejected");
  local.deltaBPlus2T2 = variance.plusT2;
  local.deltaBMinus2T2 = variance.minusT2;
  local.deltaB2T2 = variance.totalT2;
  CP::LocalInputView imported = local;
  imported.source = "unidentified-already-magnetic-view";
  imported.branchDirectionsKnown=false;
  imported.representation = "synthetic-magnetic-energy-unknown-direction";
  // This is an explicitly magnetic synthetic copy, not a raw SWMF export.
  // The standard BATSRUS total-energy half conversion has its own SBW gates;
  // the generic unknown-direction negative controls below remain unchanged.
  // Independently declared prescribed amplitude for the same physical state.
  CP::LocalInputView prescribed = local;
  prescribed.source = "prescribed";
  prescribed.deltaBPlus2T2 = 2.0 * kReferenceMu0 * kPlusW;
  prescribed.deltaBMinus2T2 = 2.0 * kReferenceMu0 * kMinusW;
  prescribed.deltaB2T2 = 2.0 * kReferenceMu0 * (kPlusW + kMinusW);
  CP::SpeciesProperties species;
  species.modelSpecies = 0;
  species.name = "proton";
  species.signedChargeC = 1.602176634e-19;
  species.restMassKg = 1.67262192369e-27;
  species.nucleonCount = 1.0;
  CP::SpectrumParameters spectrum;
  spectrum.kMinAtReferencePerM = 1.0e-8;
  spectrum.kMaxAtReferencePerM = 1.0e-3;
  spectrum.kMinRadialExponent = spectrum.kMaxRadialExponent = 0.0;
  const double speed = 1.0e5;
  const double mu = 0.5;
  const CP::PitchAngleResult a = CP::EvaluateJokipiiSlab(local, spectrum, species, speed, mu);
  const CP::PitchAngleResult b = CP::EvaluateJokipiiSlab(imported, spectrum, species, speed, mu);
  const CP::PitchAngleResult c = CP::EvaluateJokipiiSlab(prescribed, spectrum, species, speed, mu);
  // Independent finite-band slab reference with relativistic gyration.
  // This normalization test remains at low speed; the separate relativistic
  // regression resolves a finite resonance band at v=0.6c.
  const double beta=speed/299792458.0;
  const double omega = species.signedChargeC * local.magneticFieldT /
      species.restMassKg*std::sqrt((1.0-beta)*(1.0+beta));
  const double resonantK = omega / (speed * mu);
  const double index = 5.0 / 3.0;
  const double integral = (std::pow(spectrum.kMinAtReferencePerM, 1.0 - index) -
      std::pow(spectrum.kMaxAtReferencePerM, 1.0 - index)) / (index - 1.0);
  const double power = prescribed.deltaB2T2 / integral * std::pow(resonantK, -index);
  const double expected = (std::acos(-1.0) / 4.0) * omega * omega *
      (1.0 - mu * mu) * power /
      (local.magneticFieldT * local.magneticFieldT * speed * mu);
  const bool ok = a.status.ok() && b.status.ok() && c.status.ok() &&
      Relative(a.dMuMuPerS, expected) && Relative(b.dMuMuPerS, expected) &&
      Relative(c.dMuMuPerS, expected) &&
      Relative(a.dDmuMuDmuPerS, c.dDmuMuDmuPerS);
  Result result = Complete(ok,
      "local/unidentified/prescribed already-magnetic views produce the same normalized slab coefficient");
  result.configuration.push_back("speed_m_per_s=1e5;B_T=5e-9;mu=0.5;gyro_closure=relativistic");
  result.metrics.push_back({"Dmumu_relative_error", std::fabs(a.dMuMuPerS / expected - 1.0),
                            kRelativeTolerance, "<=", "dimensionless"});
  return result;
}

Result SpectralBinMeasure() {
  bool ok = true;
  const unsigned counts[] = {2, 4, 16, 128};
  const double expected = 2.0 * kReferenceMu0 * (kPlusW + kMinusW);
  double maximumError = 0.0;
  for (unsigned bins : counts) {
    Turbulence::State state;
    state.configuration.source = Turbulence::Source::SelfConsistentSpectral;
    state.configuration.representation = Turbulence::Representation::Spectral;
    state.configuration.spectralBins = bins;
    state.cells.push_back(Cell());
    Turbulence::CellState& cell = state.cells[0];
    cell.spectralEnergyJ.assign(2 * bins, 0.0);
    double remainingPlus = cell.ePlusJ, remainingMinus = cell.eMinusJ;
    for (unsigned k = 0; k < bins; ++k) {
      // These are already bin-integrated joules. No additional dln(k) factor
      // belongs in their sum, independently of how finely log-k is partitioned.
      cell.spectralEnergyJ[k] = k + 1 == bins ? remainingPlus : cell.ePlusJ / bins;
      cell.spectralEnergyJ[bins + k] = k + 1 == bins ? remainingMinus : cell.eMinusJ / bins;
      remainingPlus -= cell.spectralEnergyJ[k];
      remainingMinus -= cell.spectralEnergyJ[bins + k];
    }
    ok = ok && Turbulence::InitializeState(&state).ok();
    const Turbulence::DerivedCell derived = Turbulence::DeriveCell(state.cells[0]);
    maximumError = std::max(maximumError, std::fabs(derived.deltaB2T2 / expected - 1.0));
    ok = ok && derived.status.ok() && Relative(derived.deltaB2T2, expected);
  }
  Result result = Complete(ok,
      "spectral J-per-log-bin authority sums to the same magnetic variance at every bin count");
  result.configuration.push_back("spectral_bins=2,4,16,128;stored_measure=bin-integrated-J");
  result.metrics.push_back({"maximum_variance_relative_error", maximumError,
                            kRelativeTolerance, "<=", "dimensionless"});
  return result;
}

Result InvalidAndEmpty() {
  const CP::MagneticVariance empty = CP::MagneticVarianceFromEnergyDensity(0.0, 0.0);
  const CP::MagneticVariance imbalanced = CP::MagneticVarianceFromEnergyDensity(kPlusW, 0.0);
  const CP::MagneticVariance negative = CP::MagneticVarianceFromEnergyDensity(-kPlusW, 0.0);
  const CP::MagneticVariance nan = CP::MagneticVarianceFromEnergyDensity(
      std::numeric_limits<double>::quiet_NaN(), 0.0);
  const CP::MagneticVariance infinity = CP::MagneticVarianceFromEnergyDensity(
      std::numeric_limits<double>::infinity(), 0.0);
  const CP::MagneticVariance underflow = CP::MagneticVarianceFromEnergyDensity(
      std::numeric_limits<double>::denorm_min(), 0.0);
  Turbulence::CellState invalidCell = Cell();
  invalidCell.volumeM3 = 0.0;
  const Turbulence::DerivedCell badView = Turbulence::DeriveCell(invalidCell);
  invalidCell = Cell();
  invalidCell.ePlusJ = std::numeric_limits<double>::denorm_min();
  invalidCell.volumeM3 = 1.0e300;
  const Turbulence::DerivedCell densityUnderflow = Turbulence::DeriveCell(invalidCell);
  const bool ok = empty.status.ok() && empty.totalT2 == 0.0 &&
      imbalanced.status.ok() && imbalanced.minusT2 == 0.0 &&
      Relative(imbalanced.totalT2, 2.0 * kReferenceMu0 * kPlusW) &&
      !negative.status.ok() && !nan.status.ok() && !infinity.status.ok() &&
      underflow.status.code == Transport::StatusCode::StepUnderflow && !badView.status.ok() &&
      densityUnderflow.status.code == Transport::StatusCode::StepUnderflow;
  return Complete(ok,
      "empty/one-sided magnetic branches remain physical and invalid/underflow inputs fail explicitly");
}

Result RelativisticResonance() {
  CP::SpeciesProperties species;
  species.modelSpecies=0;species.name="proton";
  species.signedChargeC=1.602176634e-19;
  species.restMassKg=1.67262192369e-27;species.nucleonCount=1;
  const double lightSpeed=299792458.0,speed=0.6*lightSpeed,mu=0.5;
  // beta=3/5 gives exact gamma=5/4. The independent resonance oracle uses
  // physical momentum directly, not the production gyrofrequency helper.
  const double momentum=1.25*species.restMassKg*speed;
  const double magnetic=5e-9,variance=1e-18;
  const double expectedK=species.signedChargeC*magnetic/(momentum*mu);
  CP::LocalInputView input;
  input.source="controlled-physical-momentum";
  input.representation="finite-band-slab";input.generation=1;input.checksum=1;
  input.heliocentricRadiusM=1.495978707e11;input.magneticFieldT=magnetic;
  input.deltaB2T2=variance;input.deltaBPlus2T2=variance;
  CP::SpectrumParameters spectrum;
  spectrum.kMinRadialExponent=spectrum.kMaxRadialExponent=0;
  spectrum.kMinAtReferencePerM=0.9*expectedK;
  spectrum.kMaxAtReferencePerM=1.1*expectedK;
  spectrum.spectralIndex=3;
  const auto coefficient=CP::EvaluateJokipiiSlab(input,spectrum,species,speed,mu);
  // For P=C*k^-3, the magnetostatic slab expression simplifies to a
  // momentum-based coefficient. The erroneous rest-frequency resonance is
  // 1.25*k_res, wholly outside this band's upper edge (1.1*k_res).
  const double minimum=spectrum.kMinAtReferencePerM;
  const double maximum=spectrum.kMaxAtReferencePerM;
  const double normalization=2*variance/(1/(minimum*minimum)-1/(maximum*maximum));
  const double expected=0.25*std::acos(-1.0)*(1-mu*mu)*normalization*
      speed*mu*mu*momentum/(species.signedChargeC*magnetic*magnetic*magnetic);
  const auto rest=CP::GyrofrequencyRadPerS(magnetic,species);
  const auto relativistic=CP::RelativisticGyrofrequencyRadPerS(magnetic,species,speed);
  const auto low=CP::RelativisticGyrofrequencyRadPerS(magnetic,species,1e-6*lightSpeed);
  const CP::FlorinskiyParameters dynamic;
  const bool ok=coefficient.status.ok()&&coefficient.dMuMuPerS>0&&
      Relative(coefficient.dMuMuPerS,expected)&&rest.status.ok()&&
      relativistic.status.ok()&&Relative(relativistic.value/rest.value,0.8)&&
      low.status.ok()&&std::fabs(low.value/rest.value-1)<1e-12&&
      !CP::RelativisticGyrofrequencyRadPerS(magnetic,species,lightSpeed).status.ok()&&
      !CP::EvaluateJokipiiSlab(input,spectrum,species,lightSpeed,mu).status.ok()&&
      !CP::EvaluateJokipiiSlab(input,spectrum,species,1.01*lightSpeed,mu).status.ok()&&
      !CP::EvaluateFlorinskiySlab(input,dynamic,species,lightSpeed,mu).status.ok();
  Result result=Complete(ok,"relativistic p-based resonance and finite-band coefficient agree; low-speed limit and luminal rejection hold");
  result.configuration.push_back("beta=0.6;gamma=1.25;resonance_oracle=abs(q)*B/(p*abs(mu))");
  result.configuration.push_back("reference=Baring+2017-MNRAS464-4875-Eq19-doi:10.1093/mnras/stw2344");
  result.metrics.push_back({"relativistic_Dmumu_relative_error",
      std::fabs(coefficient.dMuMuPerS/expected-1),kRelativeTolerance,"<=","dimensionless"});
  return result;
}

class MidpointRandom final:public Transport::RandomStream {
 public:
  double UniformOpen01()override{return 0.5;}
  double Normal01()override{return 0.0;}
};
class OnePhysicalBranch final:public Transport::MeanFreePathProvider {
 public:
  OnePhysicalBranch(bool plus,double mass):plus_(plus),mass_(mass){}
  Transport::MeanFreePathSample Evaluate(double,double momentum,double)const override {
    Transport::MeanFreePathSample sample;
    const auto speed=Transport::SpeedFromMomentum(momentum,mass_,299792458.0);
    sample.status=speed.status;sample.lambdaParallelM=speed.value;
    sample.hasBranchResolvedRates=true;
    sample.nuPlusPerS=plus_?1:0;sample.nuMinusPerS=plus_?0:1;
    sample.provenance="one-physical-wave-direction";
    sample.turbulenceStateIdentity="controlled-local:B-directed:g1";
    return sample;
  }
 private:bool plus_;double mass_;
};
Result LocalBranchOrientation() {
  CP::SpeciesProperties species;
  species.modelSpecies=0;species.name="proton";
  species.signedChargeC=1.602176634e-19;
  species.restMassKg=1.67262192369e-27;species.nucleonCount=1;
  CP::LocalInputView positive;
  positive.source="controlled-local";positive.representation="B-directed-branches";
  positive.generation=1;positive.checksum=1;positive.heliocentricRadiusM=1.495978707e11;
  positive.magneticFieldT=5e-9;positive.deltaBPlus2T2=1e-18;
  positive.deltaBMinus2T2=1e-19;positive.deltaB2T2=1.1e-18;
  positive.alfvenSpeedMPerS=2e5;
  auto negative=positive;negative.magneticPolarityAlongLine=-1;
  // Reversing the B reference at fixed physical waves necessarily maps the
  // same along-line population from E+ to E-; its physical rates are unchanged.
  std::swap(negative.deltaBPlus2T2,negative.deltaBMinus2T2);
  const CP::FlorinskiyParameters parameters;
  const auto dPositive=CP::EvaluateFlorinskiyAlongFieldLine(positive,parameters,species,1e6,0.4);
  const auto dNegative=CP::EvaluateFlorinskiyAlongFieldLine(negative,parameters,species,1e6,0.4);
  bool coefficientPass=dPositive.status.ok()&&dNegative.status.ok()&&
      Relative(dPositive.dMuMuPerS,dNegative.dMuMuPerS)&&
      std::fabs(dPositive.dDmuMuDmuPerS-dNegative.dDmuMuDmuPerS)<
          1e-8*std::max(1e-30,std::fabs(dPositive.dDmuMuDmuPerS));
  negative.branchDirectionsKnown=false;
  coefficientPass=coefficientPass&&
      !CP::EvaluateFlorinskiyAlongFieldLine(negative,parameters,species,1e6,0.4).status.ok();
  const double tangent[3]={1,0,0},magnetic[3]={-1,0,0};
  const auto orientation=Transport::MagneticPolarityAlongFieldLine(magnetic,tangent);
  coefficientPass=coefficientPass&&orientation.status.ok()&&orientation.value==-1;

  const auto momentum=Transport::MomentumFromSpeed(1e6,species.restMassKg,299792458.0);
  Transport::FocusedTransportMfpState initial(0,momentum.value,0.4);
  initial.remainingOpticalDepth=0.1;
  Transport::FocusedTransportMfpBackground background;
  background.alfvenSpeedMPerS=2e5;
  OnePhysicalBranch plus(true,species.restMassKg),minus(false,species.restMassKg);
  MidpointRandom randomA,randomB;
  const auto a=Transport::AdvanceFocusedTransportMfp(initial,background,
      species.restMassKg,299792458.0,0.2,0.2,plus,randomA,nullptr);
  background.magneticPolarityAlongLine=-1;
  const auto b=Transport::AdvanceFocusedTransportMfp(initial,background,
      species.restMassKg,299792458.0,0.2,0.2,minus,randomB,nullptr);
  const bool eventPass=a.status.ok()&&b.status.ok()&&
      a.diagnostics.scatteringEvents==1&&b.diagnostics.scatteringEvents==1&&
      Relative(a.state.momentumKgMPerS,b.state.momentumKgMPerS)&&
      Relative(a.state.mu,b.state.mu)&&Relative(a.displacementM,b.displacementM)&&
      a.diagnostics.plusBranchEvents==1&&b.diagnostics.minusBranchEvents==1;

  bool advectionPass=true,checkpointPass=true;
  for(int spectral=0;spectral<2;++spectral) {
    Turbulence::State states[3];
    for(int scenario=0;scenario<3;++scenario) {
      auto& state=states[scenario];
      state.configuration.source=spectral?Turbulence::Source::SelfConsistentSpectral:
          Turbulence::Source::SelfConsistentIntegrated;
      state.configuration.representation=spectral?Turbulence::Representation::Spectral:
          Turbulence::Representation::Integrated;
      state.configuration.spectralBins=4;
      state.configuration.coupling=Turbulence::CouplingPolicy::Disabled;
      state.configuration.reflectionEnabled=false;state.configuration.cascadeEnabled=false;
      state.configuration.shockInjectionEnabled=false;state.configuration.periodicBoundaries=true;
      state.cells.resize(3);
      for(auto& cell:state.cells) {
        cell.lengthM=1;cell.volumeM3=1;cell.alfvenSpeedMPerS=1;
        cell.magneticPolarityAlongLine=scenario==0?1:-1;
      }
      if(scenario==1)state.cells[1].eMinusJ=1;
      else state.cells[1].ePlusJ=1;
      advectionPass=advectionPass&&Turbulence::InitializeState(&state).ok()&&
          Turbulence::Advance(&state,0.1).status.ok();
    }
    advectionPass=advectionPass&&Relative(states[0].cells[2].ePlusJ,0.1)&&
        Relative(states[1].cells[2].eMinusJ,0.1)&&Relative(states[2].cells[0].ePlusJ,0.1);
    std::string checkpoint;Turbulence::State restored;
    checkpointPass=checkpointPass&&Turbulence::SerializeCheckpoint(states[1],&checkpoint).ok()&&
        Turbulence::DeserializeCheckpoint(checkpoint,&restored).ok()&&
        restored.cells[0].magneticPolarityAlongLine==-1;
  }
  Result result=Complete(coefficientPass&&eventPass&&advectionPass&&checkpointPass,
      "B-directed local branches preserve physical QLT, wave-frame events and integrated/spectral transport under polarity reversal");
  result.metrics.push_back({"QLT_orientation_pass",coefficientPass?1.0:0.0,1,"==","boolean"});
  result.metrics.push_back({"MFP_orientation_pass",eventPass?1.0:0.0,1,"==","boolean"});
  result.metrics.push_back({"wave_advection_orientation_pass",advectionPass?1.0:0.0,1,"==","boolean"});
  result.metrics.push_back({"checkpoint_orientation_pass",checkpointPass?1.0:0.0,1,"==","boolean"});
  return result;
}

Descriptor Make(const char* id, const char* name, TestCallback callback) {
  Descriptor descriptor;
  descriptor.id = id;
  descriptor.name = name;
  descriptor.group = "wave-normalization";
  descriptor.description = name;
  descriptor.initialization = InitializationLevel::None;
  descriptor.supportedBuildModes = "linked CLI and dependency-light strict C++11 sanitizer runner";
  descriptor.runtime = RuntimeClass::Routine;
  descriptor.seedPolicy = "deterministic; no stochastic draws";
  descriptor.stateIsolation = "private core states/source views; no native PIC/SWMF import claim";
  descriptor.callback = callback;
  return descriptor;
}

}  // namespace

std::vector<Descriptor> WaveNormalizationDescriptors() {
  std::vector<Descriptor> result;
  result.push_back(Make("WN01", "Magnetic energy amplitude closure", IntegratedAmplitude));
  result.push_back(Make("WN02", "Source-equivalent slab coefficient", SourceCoefficientClosure));
  result.push_back(Make("WN03", "Spectral bin-energy normalization", SpectralBinMeasure));
  result.push_back(Make("WN04", "Magnetic normalization domain statuses", InvalidAndEmpty));
  result.push_back(Make("WN05", "Relativistic finite-band QLT resonance", RelativisticResonance));
  result.push_back(Make("WN06", "Local B-directed branch orientation", LocalBranchOrientation));
  return result;
}

}  // namespace Testing
}  // namespace SEP
