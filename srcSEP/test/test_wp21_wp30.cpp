#include "../util/sep_flux_tube_geometry_core.h"
#include "../util/sep_injection_spectrum.h"
#include "../util/sep_run_configuration.h"
#include "../util/sep_sampling_products.h"
#include "../util/sep_shock_source_core.h"
#include "../util/sep_transactional_output.h"

#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include <string>
#include <sys/stat.h>
#include <unistd.h>

namespace {

int Check(bool condition,const char* label) {
  if (condition) { std::cout<<"PASS "<<label<<'\n'; return 0; }
  std::cerr<<"FAIL "<<label<<'\n'; return 1;
}

bool Near(double a,double b,double relative=1.0e-12) {
  return std::fabs(a-b)<=relative*std::max(1.0,std::max(std::fabs(a),std::fabs(b)));
}

SEP::Shock::TrajectoryConfiguration ShockConfiguration() {
  SEP::Shock::TrajectoryConfiguration c;
  c.launchEpochS=0.0; c.launchRadiusM=0.05;
  const double radii[]={0.1,0.15,0.3,0.5,1.3};
  const double speeds[]={1800.0,1500.0,1500.0,1100.0,900.0};
  for (std::size_t i=0;i<5;++i) {
    SEP::Shock::SpeedKnot k; k.radiusM=radii[i];k.speedMPerS=speeds[i];
    c.knots.push_back(k);
  }
  return c;
}

}  // namespace

int main() {
  int failures=0;

  // WP21: the 0.3--0.5 interval must connect exactly to both stated knot
  // speeds.  Re-evaluation from launch gives bitwise equality after any split.
  const SEP::Shock::TrajectoryConfiguration trajectory=ShockConfiguration();
  const SEP::Transport::ScalarResult left=
      SEP::Shock::SpeedAtRadius(trajectory,0.3);
  const SEP::Transport::ScalarResult right=
      SEP::Shock::SpeedAtRadius(trajectory,0.5);
  SEP::Shock::TrajectoryState one,split;
  SEP::Shock::StateAtEpoch(trajectory,1.0e-4,&one);
  SEP::Shock::StateAtEpoch(trajectory,4.0e-5,&split);
  SEP::Shock::AdvanceToEpoch(&split,1.0e-4);
  std::string shockCheckpoint;
  SEP::Shock::SerializeTrajectory(split,&shockCheckpoint);
  SEP::Shock::TrajectoryState restarted;
  const SEP::Transport::Status shockRestart=
      SEP::Shock::DeserializeTrajectory(shockCheckpoint,&restarted);
  failures+=Check(left.status.ok() && right.status.ok() && left.value==1500.0 &&
      right.value==1100.0 && one.radiusM==split.radiusM && shockRestart.ok() &&
      restarted.radiusM==split.radiusM,"WP21 continuous partition-invariant shock trajectory/restart");

  // WP22: a line beginning outside, entering, and exiting the sphere returns
  // both oriented crossings; FirstOutward selects the physical downstream exit.
  std::vector<SEP::Shock::Vector3> line(3);
  line[0].x=-2.0; line[1].x=0.0; line[2].x=2.0;
  const SEP::Shock::IntersectionResult all=SEP::Shock::IntersectSphere(
      line,1.0,1.0e-12,SEP::Shock::IntersectionPolicy::All);
  const SEP::Shock::IntersectionResult outward=SEP::Shock::IntersectSphere(
      line,1.0,1.0e-12,SEP::Shock::IntersectionPolicy::FirstOutward);
  std::vector<SEP::Shock::Vector3> degenerate(2);
  degenerate[0].x=degenerate[1].x=2.0;
  const SEP::Shock::IntersectionResult invalid=SEP::Shock::IntersectSphere(
      degenerate,1.0,1.0e-12,SEP::Shock::IntersectionPolicy::All);
  failures+=Check(all.status==SEP::Shock::IntersectionStatus::Ok &&
      all.intersections.size()==2 &&
      all.intersections[0].orientation==SEP::Shock::CrossingOrientation::Inward &&
      outward.intersections.size()==1 && outward.intersections[0].positionM.x==1.0 &&
      invalid.status==SEP::Shock::IntersectionStatus::InvalidGeometry,
      "WP22 typed multi-crossing and degenerate intersection policy");

  // WP23: integrate the normalized dN/dp law on a fine log grid and exercise
  // the exact index-one logarithmic limit without log(p_min)^fractional-power.
  SEP::Injection::Spectrum spectrum;
  spectrum.measure=SEP::Injection::Measure::Momentum;
  spectrum.minimum=1.0e-20; spectrum.maximum=1.0e-17; spectrum.powerIndex=2.5;
  long double integral=0.0L; const int quadrature=200000;
  for (int i=0;i<quadrature;++i) {
    const double lo=spectrum.minimum*std::exp(
        std::log(spectrum.maximum/spectrum.minimum)*i/quadrature);
    const double hi=spectrum.minimum*std::exp(
        std::log(spectrum.maximum/spectrum.minimum)*(i+1)/quadrature);
    const SEP::Transport::ScalarResult density=
        SEP::Injection::ProbabilityDensity(spectrum,std::sqrt(lo*hi));
    integral+=density.value*(hi-lo);
  }
  spectrum.powerIndex=1.0;
  const SEP::Transport::ScalarResult median=SEP::Injection::InverseCdf(spectrum,0.5);
  failures+=Check(Near(static_cast<double>(integral),1.0,2.0e-8) &&
      median.status.ok() && Near(median.value,std::sqrt(spectrum.minimum*spectrum.maximum)),
      "WP23 normalized spectrum and limiting-index inverse CDF");

  // WP24: all semantic source fields are position/purpose sensitive, while an
  // identical tuple produces byte-identical streams independently of scheduling.
  SEP::Injection::RandomKey key;
  key.campaign=8;key.event=7;key.fieldLine=6;key.species=5;key.macroparticle=4;
  key.purpose=SEP::Injection::RandomPurpose::Spectrum;
  SEP::Injection::RandomKey swapped=key;
  swapped.fieldLine=key.species; swapped.species=key.fieldLine;
  SEP::Transport::KeyedRandomStream randomA=SEP::Injection::MakeRandomStream(key);
  SEP::Transport::KeyedRandomStream randomB=SEP::Injection::MakeRandomStream(key);
  failures+=Check(SEP::Injection::HashRandomKey(key)!=
      SEP::Injection::HashRandomKey(swapped) &&
      randomA.UniformOpen01()==randomB.UniformOpen01(),
      "WP24 tagged keyed-source reproducibility and field-order separation");

  // WP25: energy uses the upstream normal-relative speed, closes exactly across
  // branches, and vanishes for a stationary shock relative to the upstream flow.
  SEP::Shock::TurbulenceSourceInput source;
  source.sweptVolumeM3=3.0;source.upstreamMassDensityKgPerM3=2.0;
  source.shockNormalSpeedMPerS=10.0;source.upstreamNormalSpeedMPerS=4.0;
  source.efficiency=0.2;source.plusBranchFraction=0.25;
  const SEP::Shock::TurbulenceSourceEnergy wave=
      SEP::Shock::ComputeTurbulenceSource(source);
  source.upstreamNormalSpeedMPerS=10.0;
  const SEP::Shock::TurbulenceSourceEnergy stationary=
      SEP::Shock::ComputeTurbulenceSource(source);
  failures+=Check(wave.status.ok() && Near(wave.totalJ,21.6) &&
      wave.plusJ+wave.minusJ==wave.totalJ && stationary.totalJ==0.0,
      "WP25 normal-relative ledgered shock turbulence energy");

  // WP26: A|B| returns the owned Phi_i and conservative refinement assigns the
  // residual to the final child so the summed flux closes in double precision.
  SEP::FieldLine::FluxTubeGeometryCore::MagneticFluxRecord flux;
  flux.magneticFluxWb=12.0;flux.generation=3;flux.provenance="controlled seed";
  const SEP::Units::AreaM2 area=
      SEP::FieldLine::FluxTubeGeometryCore::AreaFromMagneticFlux(
          flux,SEP::Units::MagneticFieldT(3.0));
  std::vector<double> fractions;fractions.push_back(0.2);fractions.push_back(0.3);
  fractions.push_back(0.5);
  std::vector<SEP::FieldLine::FluxTubeGeometryCore::MagneticFluxRecord> children;
  const SEP::Transport::Status redistribution=
      SEP::FieldLine::FluxTubeGeometryCore::RedistributeMagneticFlux(
          flux,fractions,&children);
  std::string fluxCheckpoint;
  std::vector<SEP::FieldLine::FluxTubeGeometryCore::MagneticFluxRecord> restoredFlux;
  const SEP::Transport::Status fluxSerialized=
      SEP::FieldLine::FluxTubeGeometryCore::SerializeMagneticFluxTable(
          children,&fluxCheckpoint);
  const SEP::Transport::Status fluxRestored=
      SEP::FieldLine::FluxTubeGeometryCore::DeserializeMagneticFluxTable(
          fluxCheckpoint,&restoredFlux);
  failures+=Check(area.Value()*3.0==12.0 && redistribution.ok() &&
      children[0].magneticFluxWb+children[1].magneticFluxWb+
      children[2].magneticFluxWb==12.0 && children[0].generation==4 &&
      fluxSerialized.ok() && fluxRestored.ok() && restoredFlux.size()==3,
      "WP26 owned magnetic flux and conservative refinement");

  // WP27: compare the local-field relativistic Larmor formula directly and
  // prove that a superluminal diagnostic state is rejected, never capped.
  const double c=299792458.0,mass=1.67262192369e-27,charge=1.602176634e-19;
  const SEP::SamplingCore::Kinematics valid=SEP::SamplingCore::EvaluateKinematics(
      1.0e8,2.0e8,mass,charge,5.0e-9,c);
  const SEP::SamplingCore::Kinematics invalidParticle=
      SEP::SamplingCore::EvaluateKinematics(c,0.0,mass,charge,5.0e-9,c);
  const double gamma=1.0/std::sqrt(1.0-(5.0e16/(c*c)));
  failures+=Check(valid.status.ok() && Near(valid.larmorRadiusM,
      gamma*mass*2.0e8/(charge*5.0e-9)) && !invalidParticle.status.ok(),
      "WP27 local-field relativistic Larmor and invalid exclusion");

  // WP28: below/on/internal/on-upper/above boundaries have explicit semantics;
  // no event is folded into an edge bin and weighted uncertainty is retained.
  SEP::SamplingCore::BinGrid grid;grid.edges.push_back(1.0);grid.edges.push_back(2.0);
  grid.edges.push_back(4.0);
  SEP::SamplingCore::WeightedBins bins;
  SEP::SamplingCore::Initialize(grid,&bins);
  SEP::SamplingCore::Accumulate(&bins,0.5,1.0);
  SEP::SamplingCore::Accumulate(&bins,1.0,2.0);
  SEP::SamplingCore::Accumulate(&bins,2.0,3.0);
  SEP::SamplingCore::Accumulate(&bins,4.0,4.0);
  SEP::SamplingCore::NormalizationContext normalization;
  normalization.volumeM3=2.0;normalization.areaM2=4.0;
  normalization.durationS=5.0;normalization.representativeSpeedMPerS=8.0;
  const SEP::SamplingCore::ProductValues density=SEP::SamplingCore::BuildProduct(
      bins,SEP::SamplingCore::Product::NumberDensityPerEnergy,normalization);
  const SEP::SamplingCore::ProductValues fluxProduct=SEP::SamplingCore::BuildProduct(
      bins,SEP::SamplingCore::Product::DirectionalCrossingFlux,normalization);
  failures+=Check(bins.counters.underflow==1 && bins.counters.overflow==1 &&
      bins.sumWeight[0]==2.0 && bins.sumWeight[1]==3.0 &&
      SEP::SamplingCore::EffectiveSampleSize(bins,0)==1.0 &&
      density.status.ok() && density.value[0]==1.0 &&
      fluxProduct.status.ok() && fluxProduct.value[0]==0.1,
      "WP28 non-clamping bins and weighted effective sample size");

  // WP29: a completed payload validates through its checksum sidecar; corrupting
  // the payload makes the reader reject it.  All mutation stays in mkdtemp.
  char temporary[]="/tmp/srcsep-wp29-XXXXXX";
  const char* directory=::mkdtemp(temporary);
  SEP::Output::ArtifactMetadata metadata;
  metadata.schemaVersion="srcsep-test-v1";metadata.recordCount=2;
  metadata.configurationFingerprint="controlled";
  const SEP::Output::WriteResult written=SEP::Output::WriteTransactional(
      directory ? directory : "", "products/sample.dat", "alpha\nbeta\n",metadata);
  std::string payload;
  const SEP::Transport::Status completed=
      SEP::Output::ValidateCompletedArtifact(written.finalPath,&payload);
  const SEP::Output::WriteResult duplicate=SEP::Output::WriteTransactional(
      directory ? directory : "", "products/sample.dat", "replacement",metadata);
  std::ofstream corrupt(written.finalPath.c_str(),std::ios::app);corrupt<<"damage";corrupt.close();
  const SEP::Transport::Status rejected=
      SEP::Output::ValidateCompletedArtifact(written.finalPath,&payload);
  failures+=Check(directory && written.status.ok() && completed.ok() &&
      !duplicate.status.ok() &&
      !rejected.ok() && !SEP::Output::ValidateRelativePath("../escape").ok(),
      "WP29 transactional checksum/completion/path validation");
  if (directory) {
    std::remove((written.finalPath+".manifest").c_str());
    std::remove(written.finalPath.c_str());
    ::rmdir((std::string(directory)+"/products").c_str());
    ::rmdir(directory);
  }

  // WP30: CLI authority survives an input-file overlay, serialized configuration
  // round-trips with the same fingerprint, and restart mismatch is explicit.
  SEP::Run::Configuration input=SEP::Run::Defaults();
  input.totalIterations.value=20;input.totalIterations.source=SEP::Run::ValueSource::InputFile;
  SEP::Run::Configuration cli=SEP::Run::Defaults();
  cli.totalIterations.value=30;cli.totalIterations.source=SEP::Run::ValueSource::CommandLine;
  const SEP::Run::Configuration merged=SEP::Run::Merge(input,cli);
  SEP::Run::FrozenConfiguration frozen,roundTrip,other;
  const SEP::Transport::Status frozenStatus=SEP::Run::FrozenConfiguration::Create(merged,&frozen);
  std::string runCheckpoint;SEP::Run::Serialize(frozen,&runCheckpoint);
  const SEP::Transport::Status roundStatus=SEP::Run::Deserialize(runCheckpoint,&roundTrip);
  SEP::Run::Configuration changed=merged;changed.totalIterations.value=31;
  SEP::Run::FrozenConfiguration::Create(changed,&other);
  failures+=Check(merged.totalIterations.value==30 && frozenStatus.ok() &&
      roundStatus.ok() && frozen.fingerprint()==roundTrip.fingerprint() &&
      !SEP::Run::VerifyRestartCompatibility(frozen,other).ok(),
      "WP30 precedence, immutable roundtrip, and restart compatibility");

  return failures==0 ? 0 : 1;
}
