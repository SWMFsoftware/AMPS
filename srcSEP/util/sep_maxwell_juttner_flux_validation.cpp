#include "sep_maxwell_juttner_flux_validation.h"
#include "sep_maxwell_juttner_flux.h"
#include <algorithm>
#include <cmath>
#include <limits>

namespace SEP { namespace Testing { namespace {
namespace J=Sources::MaxwellJuttner;
namespace T=Transport;
using Wide=long double;
constexpr Wide Pi=3.141592653589793238462643383279502884L;
constexpr double Mass=1.67262192369e-27; // test species only, not a provider default
struct Checks {
  unsigned failures=0;std::vector<Metric> metrics;
  void Require(bool condition){if(!condition)++failures;}
  void Relative(double actual,double expected,double tolerance) {
    Require(std::isfinite(actual)&&std::isfinite(expected)&&expected!=0.0&&std::abs(actual/expected-1.0)<=tolerance);
  }
  Result Finish(const char* message,bool stochastic=false,double observable=0.0)const {
    Result result;result.status=failures?Status::Fail:Status::Pass;
    result.message=failures?"material-plane Maxwell-Juttner flux assertion failed":message;
    result.configuration={"evidence_level=full-particle-material-plane-core;native-guiding-center-source-unqualified",
        "measure=proper-density-isotropic-rest-f-d3p;simulation-area-time-relative-crossing-flux",
        "sampling=exact-gamma2-gamma3-mixture;cosine-rest-hemisphere;no-cutoff-no-tail-redraw",
        "normalization=NIST-bounded-cold-scaled-K2;C++-Bessel-hot-recurrence;independent-integral-tolerance=2e-12"};
    result.metrics=metrics;result.metrics.push_back({"assertion_failures",static_cast<double>(failures),0.0,"<=","count"});
    if(stochastic){result.hasSeed=true;result.seed=TestSeed(20261117);result.metrics.push_back({"actual_seeded_observable",observable,0.0,"record","dimensionless"});}
    return result;
  }
};
J::MaterialBoundaryInput Input(double theta) {
  J::MaterialBoundaryInput input;input.properNumberDensityPerM3=7.0;input.restMassKg=Mass;
  input.properTemperatureK=theta*Mass*J::SpeedOfLightMPerS*J::SpeedOfLightMPerS/J::BoltzmannJPerK;
  return input;
}
bool SameRandom(const T::KeyedRandomStreamState& a,const T::KeyedRandomStreamState& b) {
  return a.schema==b.schema&&a.keyIdentity==b.keyIdentity&&a.integerState==b.integerState&&
      a.hasSpareNormal==b.hasSpareNormal&&a.spareNormal==b.spareNormal;
}
bool SameSample(const J::FluxSample& a,const J::FluxSample& b) {
  return a.restKineticGamma==b.restKineticGamma&&a.restGamma==b.restGamma&&a.simulationGamma==b.simulationGamma&&
      a.restMomentumKgMPerS==b.restMomentumKgMPerS&&a.simulationMomentumKgMPerS==b.simulationMomentumKgMPerS&&
      a.relativeInwardSpeedMPerS==b.relativeInwardSpeedMPerS;
}
double Magnitude(const std::array<double,3>& p){return std::hypot(p[0],std::hypot(p[1],p[2]));}

Wide IndependentMeanBeta(Wide theta) {
  // Integrate the ORIGINAL volume radial measure, not the gamma-mixture CDF
  // or production Bessel formula. x=theta*t^2 makes the cold square-root
  // endpoint smooth. Scaling keeps all factors finite from cold to hot.
  // For y=t^2>=100, integrand in y is <=4*y^2*exp(-y), so the omitted
  // positive tail is <=4*exp(-100)*(100^2+200+2)<1.6e-39.
  const Wide oneWeight=1.0L/(1.0L+theta),thermalWeight=theta/(1.0L+theta);
  const Wide rootConstant=2.0L/(theta+2.0L),rootThermal=theta/(theta+2.0L);
  const unsigned intervals=16384;const Wide step=10.0L/intervals;
  auto integrand=[&](Wide t){const Wide y=t*t;return 2.0L*y*(oneWeight+thermalWeight*y)*std::sqrt(rootConstant+rootThermal*y)*std::exp(-y);};
  Wide sum=integrand(0.0L)+integrand(10.0L);
  for(unsigned i=1;i<intervals;++i)sum+=(i%2?4.0L:2.0L)*integrand(i*step);
  const Wide denominator=sum*step/3.0L;
  return 2.0L*std::sqrt(theta/(theta+2.0L))/denominator;
}
Result Normalization() {
  Checks check;
  // 90-decimal mpmath 1.3.0 reference, generated independently from the
  // invariant rest-volume integral/Bessel normalization, never by this kernel.
  // Primary formulas: Zenitani(2015) eq1,22-28; DLMF10.32.9/10.40.2.
  const double theta[]={1e-12,1e-6,0.01,0.03125,0.032,0.1,1.0,10.0,1e6};
  const double reference[]={1.595769121604334413778380534394095229076e-6,
      0.001595767725309058335284834826878717257008,0.158193622214308610804813896366870039373,
      0.2746027663975783947878916521951930229994,0.2777019017132188606308569897181352648955,
      0.4643453984070374973216110865907530071787,0.9056391780882437200917173595102237690094,
      0.9977958699969102860858721036228232508021,0.9999999999997500003333322282433203984483};
  double largest=0.0;
  for(unsigned i=0;i<9;++i) {
    const auto input=Input(theta[i]);J::FluxNormalization flux;
    const auto status=J::EvaluateMaterialFlux(input,&flux);check.Require(status.ok());
    if(!status.ok())continue;
    const double beta=flux.meanRestSpeedMPerS/J::SpeedOfLightMPerS;
    check.Relative(beta,reference[i],2e-12);check.Relative(beta,static_cast<double>(IndependentMeanBeta(theta[i])),2e-12);
    check.Relative(flux.numberFluxPerM2S,input.properNumberDensityPerM3*J::SpeedOfLightMPerS*reference[i]/4.0,2e-12);
    check.Require(flux.meanRestSpeedMPerS>0.0&&flux.meanRestSpeedMPerS<=J::SpeedOfLightMPerS);
    if(theta[i]<=0.03125)check.Require(flux.method==J::NormalizationMethod::BoundedAsymptotic&&flux.asymptoticRelativeRemainderBound<=1e-17);
    largest=std::max(largest,std::abs(beta/reference[i]-1.0));
  }
  auto input=Input(0.0);J::FluxNormalization cold;check.Require(J::EvaluateMaterialFlux(input,&cold).ok()&&cold.numberFluxPerM2S==0.0&&cold.meanRestSpeedMPerS==0.0);
  input=Input(0.1);input.properNumberDensityPerM3=0.0;check.Require(J::EvaluateMaterialFlux(input,&cold).ok()&&cold.numberFluxPerM2S==0.0);
  const auto small=Input(1e-12);check.Require(J::EvaluateMaterialFlux(small,&cold).ok());
  check.Relative(cold.numberFluxPerM2S,small.properNumberDensityPerM3*std::sqrt(J::BoltzmannJPerK*small.properTemperatureK/(2.0*static_cast<double>(Pi)*Mass)),2e-12);
  // Only the auxiliary analytic remainder underflows here; the relativistic
  // rate is representable. Its upper bound must round upward, not turn a
  // positive cold rate into an error or claim an exactly zero remainder.
  const auto extremelyCold=Input(1e-200);check.Require(J::EvaluateMaterialFlux(extremelyCold,&cold).ok());
  check.Require(cold.asymptoticRelativeRemainderBound>0.0);
  check.Relative(cold.meanRestSpeedMPerS/J::SpeedOfLightMPerS,std::sqrt(8e-200/static_cast<double>(Pi)),2e-12);
  check.metrics.push_back({"maximum_reference_relative_error",largest,2e-12,"<=","relative"});
  return check.Finish("proper-density SI flux agrees with independent high-precision values and original positive radial integrals across cold and relativistic temperatures");
}
Result InvariantMeasure() {
  Checks check;double largestIdentity=0.0,largestQuadrature=0.0;
  // Deliberately integrate in fixed simulation-oriented spherical axes. The
  // incoming half-space is discovered by the DIRECT lab velocity and moving
  // boundary, not by the production transformed-axis/hemisphere construction.
  for(const std::array<double,3> normal:{std::array<double,3>{{1,0,0}},{{0,1,0}},{{0.3,0.4,0.5}}}) {
    const Wide norm=std::hypot(static_cast<Wide>(normal[0]),std::hypot(static_cast<Wide>(normal[1]),static_cast<Wide>(normal[2])));
    const Wide n[3]={normal[0]/norm,normal[1]/norm,normal[2]/norm};
    const Wide bulkBeta=0.8L,gamma=1.0L/std::sqrt(1.0L-bulkBeta*bulkBeta);
    const Wide transformed=std::sqrt(1.0L-bulkBeta*bulkBeta*n[0]*n[0]);
    auto input=Input(0.2);input.inwardNormal=normal;
    input.bulkVelocityMPerS[0]=input.boundaryVelocityMPerS[0]=static_cast<double>(bulkBeta*J::SpeedOfLightMPerS);
    J::FluxNormalization flux;check.Require(J::EvaluateMaterialFlux(input,&flux).ok());check.Relative(flux.transformedNormalMagnitude,static_cast<double>(transformed),2e-14);
    for(const Wide kinetic:{0.01L,1.0L,10.0L}) {
      const Wide g=1.0L+kinetic,q=std::sqrt(kinetic*(kinetic+2.0L)),v=q/g;
      Wide integral=0.0L;const unsigned polar=512,azimuth=1024;
      for(unsigned i=0;i<polar;++i) {
        // Uniform theta with its explicit sin(theta) measure avoids the
        // square-root endpoint quadrature error of uniform-cosine midpoint
        // cells. The hemisphere is still discovered in the lab calculation;
        // no production aHat or transformed angular CDF enters this integral.
        const Wide angle=Pi*(i+0.5L)/polar,mu=std::cos(angle),sine=std::sin(angle);
        for(unsigned j=0;j<azimuth;++j) {
          const Wide phi=2.0L*Pi*(j+0.5L)/azimuth;
          const Wide rest[3]={q*mu,q*sine*std::cos(phi),q*sine*std::sin(phi)};
          const Wide labGamma=gamma*(g+bulkBeta*rest[0]);
          const Wide lab[3]={gamma*(rest[0]+bulkBeta*g),rest[1],rest[2]};
          const Wide direct=(lab[0]*n[0]+lab[1]*n[1]+lab[2]*n[2])/labGamma-bulkBeta*n[0];
          const Wide jacobian=labGamma/g;
          const Wide independentlySimplified=(n[0]*rest[0]/gamma+n[1]*rest[1]+n[2]*rest[2])/g;
          const Wide residual=std::abs(jacobian*direct-independentlySimplified);
          largestIdentity=std::max(largestIdentity,static_cast<double>(residual));
          integral+=sine*jacobian*std::max(direct,0.0L);
        }
      }
      integral*=Pi/(2.0L*polar*azimuth);
      const double relative=static_cast<double>(std::abs(integral/(transformed*v/4.0L)-1.0L));
      largestQuadrature=std::max(largestQuadrature,relative);check.Require(relative<=1.2e-5);
    }
  }
  check.Require(largestIdentity<=2e-17);
  check.metrics.push_back({"maximum_invariant_measure_residual",largestIdentity,2e-17,"<=","dimensionless"});
  check.metrics.push_back({"direct_moving_plane_angular_quadrature_relative_error",largestQuadrature,1.2e-5,"<=","relative"});
  return check.Finish("independent boosted moving-plane quadrature and pointwise Jacobian identity verify parallel, transverse and tilted material normals");
}
Result RadialAndAngular() {
  Checks check;const unsigned count=100000;double observable=0.0,maxZ=0.0;
  for(unsigned caseIndex=0;caseIndex<3;++caseIndex) {
    const double theta[]={1e-6,0.1,10.0};const double value=theta[caseIndex];auto input=Input(value);
    T::KeyedRandomStream random(TestSeed(20261117),caseIndex+1,UINT64_C(0x4d4a4658),0);
    Wide sumX=0.0L,sumMu=0.0L,sumMu2=0.0L,sumY=0.0L,sumZ=0.0L;unsigned below[5]={0,0,0,0,0};
    const double thresholds[]={0.5,1.0,2.0,4.0,8.0};
    for(unsigned i=0;i<count;++i) {
      J::FluxSample sample;const auto status=J::SampleMaterialFlux(input,&random,&sample);check.Require(status.ok());if(!status.ok())break;
      const double magnitude=Magnitude(sample.restMomentumKgMPerS),mu=sample.restMomentumKgMPerS[0]/magnitude;
      const double unit=sample.restKineticGamma/value;
      sumX+=unit;sumMu+=mu;sumMu2+=mu*mu;sumY+=sample.restMomentumKgMPerS[1]/magnitude;sumZ+=sample.restMomentumKgMPerS[2]/magnitude;
      for(unsigned k=0;k<5;++k)if(unit<=thresholds[k])++below[k];
      check.Require(sample.relativeInwardSpeedMPerS>0.0&&sample.restKineticGamma>0.0&&mu>0.0&&mu<=1.0);
      if(i==0)observable+=unit+mu;
    }
    const double p3=value/(1.0+value),mean=2.0+p3,variance=6.0+6.0*p3-mean*mean;
    const double z=std::abs(static_cast<double>(sumX/count)-mean)/std::sqrt(variance/count);check.Require(z<=6.5);maxZ=std::max(maxZ,z);
    check.Require(std::abs(static_cast<double>(sumMu/count)-2.0/3.0)<=6.5*std::sqrt((0.5-4.0/9.0)/count));
    check.Require(std::abs(static_cast<double>(sumMu2/count)-0.5)<=6.5/std::sqrt(12.0*count));
    check.Require(std::abs(static_cast<double>(sumY/count))<=6.5/std::sqrt(4.0*count)&&std::abs(static_cast<double>(sumZ/count))<=6.5/std::sqrt(4.0*count));
    for(unsigned k=0;k<5;++k) {
      const double t=thresholds[k];const double cdf=1.0-std::exp(-t)*(1.0+t+p3*t*t/2.0);
      const double error=std::abs(static_cast<double>(below[k])/count-cdf);
      check.Require(error<=6.5*std::sqrt(cdf*(1.0-cdf)/count)+1.0/count);
    }
    // Both one-component gamma shortcuts fail the declared mean gate where
    // they materially alter the selected radial law. This control would not
    // be visible from checking only positivity or the native injection count.
    if(caseIndex==0)check.Require(std::abs(mean-3.0)>20.0*std::sqrt(variance/count));
    if(caseIndex==2)check.Require(std::abs(mean-2.0)>20.0*std::sqrt(variance/count));
  }
  check.metrics.push_back({"maximum_radial_mean_standard_errors",maxZ,6.5,"<=","standard-errors"});
  check.metrics.push_back({"samples_per_temperature",count,100000,"=","count"});
  return check.Finish("predeclared finite ensembles verify exact mixture CDF/mean and independent cosine-weighted angular moments",true,observable);
}
Result ColdMoments() {
  Checks check;const unsigned count=160000;double observable=0.0,maxZ=0.0;
  for(unsigned caseIndex=0;caseIndex<2;++caseIndex) {
    const double theta=caseIndex?1e-5:1e-8;const auto input=Input(theta);
    T::KeyedRandomStream random(TestSeed(20261117),caseIndex+11,UINT64_C(0x4d4a4658),0);
    Wide mean=0.0L,normal2=0.0L,transverseY2=0.0L,transverseZ2=0.0L;
    const Wide scale=std::sqrt(theta)*J::SpeedOfLightMPerS;
    for(unsigned i=0;i<count;++i) {
      J::FluxSample sample;const auto status=J::SampleMaterialFlux(input,&random,&sample);check.Require(status.ok());if(!status.ok())break;
      Wide v[3];for(unsigned j=0;j<3;++j)v[j]=sample.restMomentumKgMPerS[j]/(static_cast<Wide>(Mass)*sample.restGamma*scale);
      mean+=v[0];normal2+=v[0]*v[0];transverseY2+=v[1]*v[1];transverseZ2+=v[2]*v[2];
      if(i==0)observable+=static_cast<double>(v[0]+v[1]);
    }
    // Compare the SAME relativistic sampler to its analytic Maxwellian
    // limit, allowing an explicit finite-temperature remainder, never a
    // switch of distributions. For x>=0, 0<=2x-beta^2<=3x^2; the exact
    // mixture moments give a <=6 theta relative second-moment remainder
    // here. sqrt(1+x/2)/(1+x) in [1-x,1] bounds the mean by 5 theta.
    const double expectedMean=std::sqrt(static_cast<double>(Pi)/2.0),actualMean=static_cast<double>(mean/count);
    const double meanSE=std::sqrt((2.0-static_cast<double>(Pi)/2.0)/count);
    const double z=std::max(0.0,std::abs(actualMean-expectedMean)-5.0*theta*expectedMean)/meanSE;
    check.Require(z<=6.5);maxZ=std::max(maxZ,z);
    check.Require(std::abs(static_cast<double>(normal2/count)-2.0)<=6.5*2.0/std::sqrt(static_cast<double>(count))+12.0*theta);
    check.Require(std::abs(static_cast<double>(transverseY2/count)-1.0)<=6.5*std::sqrt(2.0/count)+6.0*theta);
    check.Require(std::abs(static_cast<double>(transverseZ2/count)-1.0)<=6.5*std::sqrt(2.0/count)+6.0*theta);
  }
  check.metrics.push_back({"maximum_cold_normal_mean_standard_errors",maxZ,6.5,"<=","standard-errors"});
  check.metrics.push_back({"samples_per_temperature",count,160000,"=","count"});
  return check.Finish("unchanged Juttner flux sampler converges to flux-weighted Maxwellian normal and transverse moments with explicit finite-theta bounds",true,observable);
}
Result AtomicRejection() {
  Checks check;const auto good=Input(0.1);J::FluxSample sentinel;sentinel.restKineticGamma=19.0;sentinel.restMomentumKgMPerS={{1.0,2.0,3.0}};
  T::KeyedRandomStream random(5,6,7,8);random.Normal01();const auto before=random.SaveState();
  for(unsigned bad=0;bad<9;++bad) {
    auto input=good;
    switch(bad) {
      case 0:input.restMassKg=0.0;break;case 1:input.properTemperatureK=-1.0;break;
      case 2:input.properNumberDensityPerM3=-1.0;break;case 3:input.inwardNormal={{0,0,0}};break;
      case 4:input.bulkVelocityMPerS[0]=input.boundaryVelocityMPerS[0]=J::SpeedOfLightMPerS;break;
      case 5:input.boundaryVelocityMPerS[1]=1.0;break;
      case 6:input.inwardNormal[0]=std::numeric_limits<double>::quiet_NaN();break;
      case 7:input.properNumberDensityPerM3=0.0;break;
      case 8:input.properTemperatureK=0.0;break;
    }
    auto output=sentinel;check.Require(!J::SampleMaterialFlux(input,&random,&output).ok()&&SameRandom(before,random.SaveState())&&SameSample(output,sentinel));
  }
  for(double endpoint:{0.0,1.0,std::numeric_limits<double>::quiet_NaN()}) {
    std::array<double,6> uniforms{{0.5,0.5,0.5,0.5,0.5,0.5}};uniforms[3]=endpoint;auto output=sentinel;
    check.Require(!J::SampleMaterialFluxFromUniforms(good,uniforms,&output).ok()&&SameSample(output,sentinel));
  }
  // The physical tail is not rejected/redrawn to make SI momentum fit. This
  // extreme but finite case has representable cold rest momentum and an
  // unrepresentable boosted bulk momentum; staged RNG must roll back exactly.
  auto huge=good;huge.restMassKg=1e308;huge.properTemperatureK=1e308;
  huge.bulkVelocityMPerS[0]=huge.boundaryVelocityMPerS[0]=0.9*J::SpeedOfLightMPerS;
  auto output=sentinel;const auto status=J::SampleMaterialFlux(huge,&random,&output);
  check.Require(status.code==T::StatusCode::UnresolvedCoefficient&&SameRandom(before,random.SaveState())&&SameSample(output,sentinel));
  // Independent review counterexample: every component fits binary64, yet
  // the cold incoming thermal increment is below the boosted momentum's ULP.
  // At 100 decimal digits the encoded p_lab gives an OUTGOING relative normal
  // velocity (-1.90350872756e-10 m/s), while the unrounded sample's incoming
  // speed is +3.52943346859e-12 m/s. Component overflow checks alone cannot
  // qualify a crossing. Reject instead of changing p, redrawing or publishing
  // an inconsistent positive speed beside an outgoing particle.
  auto coldBoost=good;coldBoost.properNumberDensityPerM3=1e-100;
  coldBoost.properTemperatureK=1.0888195450003048e-27;
  coldBoost.bulkVelocityMPerS[0]=coldBoost.boundaryVelocityMPerS[0]=2997924.58;
  const std::array<double,6> coldUniforms{{0.5,0.5,0.5,0.5,0.5,0.3}};
  output=sentinel;
  check.Require(J::SampleMaterialFluxFromUniforms(coldBoost,coldUniforms,&output).code==T::StatusCode::UnresolvedCoefficient&&SameSample(output,sentinel));
  output=sentinel;
  check.Require(J::SampleMaterialFlux(coldBoost,&random,&output).code==T::StatusCode::UnresolvedCoefficient&&SameRandom(before,random.SaveState())&&SameSample(output,sentinel));
  // Sign alone is insufficient: these encoded particles are incoming, but
  // their crossings are hundreds/thousands times the intended thermal speed.
  // A relative test with a bulk-speed or max(1,speed) floor would miss them.
  for(double fraction:{0.1,0.5}) {
    coldBoost.bulkVelocityMPerS[0]=coldBoost.boundaryVelocityMPerS[0]=fraction*J::SpeedOfLightMPerS;
    output=sentinel;
    check.Require(J::SampleMaterialFluxFromUniforms(coldBoost,coldUniforms,&output).code==T::StatusCode::UnresolvedCoefficient&&SameSample(output,sentinel));
  }
  coldBoost.bulkVelocityMPerS=coldBoost.boundaryVelocityMPerS={{0.0,0.0,0.0}};
  J::FluxSample resolvedCold;check.Require(J::SampleMaterialFluxFromUniforms(coldBoost,coldUniforms,&resolvedCold).ok()&&resolvedCold.relativeInwardSpeedMPerS>0.0);
  // Exact-U high-precision review found ~5.78e-5 relative Gamma/boost error
  // for the next binary64 velocity below c. A transverse normal keeps |a|=1,
  // hiding this from the old normal-invariant check. Reject unresolved Gamma
  // formation before draws rather than claiming every subluminal drift fits
  // the selected working precision or inferring a physical velocity cutoff.
  for(double speed:{std::nextafter(J::SpeedOfLightMPerS,0.0),J::SpeedOfLightMPerS-0.01}) {
    auto illConditioned=Input(1.0);illConditioned.inwardNormal={{0,1,0}};
    illConditioned.bulkVelocityMPerS[0]=illConditioned.boundaryVelocityMPerS[0]=speed;
    output=sentinel;
    check.Require(J::SampleMaterialFlux(illConditioned,&random,&output).code==T::StatusCode::UnresolvedCoefficient&&SameRandom(before,random.SaveState())&&SameSample(output,sentinel));
  }
  J::FluxNormalization flux;flux.numberFluxPerM2S=17.0;auto overflow=good;overflow.properNumberDensityPerM3=std::numeric_limits<double>::max();
  check.Require(!J::EvaluateMaterialFlux(overflow,&flux).ok()&&flux.numberFluxPerM2S==17.0);
  return check.Finish("invalid physics, zero flux, closed variates and unrepresentable boosted tails reject atomically without RNG/event consumption or a floor");
}
Result KeyedLifecycle() {
  Checks check;auto input=Input(1.0);input.bulkVelocityMPerS=input.boundaryVelocityMPerS={{0.3*J::SpeedOfLightMPerS,-0.2*J::SpeedOfLightMPerS,0.1*J::SpeedOfLightMPerS}};
  input.inwardNormal={{1,2,3}};const auto seed=TestSeed(20261117);
  T::KeyedRandomStream a(seed,41,72,0);a.Normal01();auto manual=a;const auto before=a.SaveState();
  J::FluxSample first;check.Require(J::SampleMaterialFlux(input,&a,&first).ok());for(unsigned i=0;i<6;++i)manual.UniformOpen01();
  check.Require(SameRandom(a.SaveState(),manual.SaveState())&&a.SaveState().hasSpareNormal==before.hasSpareNormal&&a.SaveState().spareNormal==before.spareNormal);
  const auto checkpoint=a.SaveState();T::KeyedRandomStream resumed(seed,41,72,0);check.Require(resumed.RestoreState(checkpoint).ok());
  J::FluxSample direct,restored;check.Require(J::SampleMaterialFlux(input,&a,&direct).ok()&&J::SampleMaterialFlux(input,&resumed,&restored).ok()&&SameSample(direct,restored)&&SameRandom(a.SaveState(),resumed.SaveState()));
  // Other source/particle keys cannot perturb this event's future. Full
  // four-momenta and a pointwise positive moving-plane speed are compared,
  // rather than only checking that a seed label changed in a report.
  T::KeyedRandomStream other(seed,42,72,0);J::FluxSample ignored;check.Require(J::SampleMaterialFlux(input,&other,&ignored).ok());
  check.Require(!SameSample(first,ignored));
  T::KeyedRandomStream changed(seed^UINT64_C(1),41,72,0);changed.Normal01();J::FluxSample changedSample;
  check.Require(J::SampleMaterialFlux(input,&changed,&changedSample).ok()&&!SameSample(first,changedSample));
  const Wide mc=static_cast<Wide>(Mass)*J::SpeedOfLightMPerS;
  Wide p2=0.0L;for(double p:first.simulationMomentumKgMPerS)p2+=(p/mc)*(p/mc);
  check.Require(std::abs((first.simulationGamma*first.simulationGamma-1.0L)/p2-1.0L)<=2e-14L);
  const Wide gamma=first.simulationGamma;
  Wide directNormal=0.0L,norm=std::sqrt(14.0L);
  for(unsigned i=0;i<3;++i)directNormal+=(J::SpeedOfLightMPerS*(first.simulationMomentumKgMPerS[i]/mc)/gamma-input.bulkVelocityMPerS[i])*(input.inwardNormal[i]/norm);
  check.Require(directNormal>0.0L&&std::abs(directNormal/first.relativeInwardSpeedMPerS-1.0L)<2e-14L);
  return check.Finish("six-draw keyed event resumes complete cached RNG state and full boosted momentum without rank/order-dependent source randomness",true,first.restKineticGamma+first.simulationGamma);
}
} // namespace

std::vector<Descriptor> MaxwellJuttnerFluxDescriptors() {
  const char* names[]={"Proper Maxwell-Juttner flux normalization","Invariant moving material-plane flux measure",
      "Exact radial mixture and cosine angular distribution","Maxwellian low-temperature flux limit",
      "Atomic source-domain and unrepresentable-tail rejection","Keyed full-momentum source event lifecycle"};
  const TestCallback callbacks[]={Normalization,InvariantMeasure,RadialAndAngular,ColdMoments,AtomicRejection,KeyedLifecycle};
  std::vector<Descriptor> result;
  for(unsigned i=0;i<6;++i) {
    Descriptor descriptor;descriptor.id="MJ0"+std::to_string(i+1);descriptor.name=names[i];descriptor.description=names[i];
    descriptor.group="maxwell-juttner-flux";descriptor.initialization=InitializationLevel::None;descriptor.runtime=RuntimeClass::Routine;
    const bool stochastic=i==2||i==3||i==5;
    descriptor.seedCapability=stochastic?SeedCapability::Context:SeedCapability::Deterministic;descriptor.repeatSupported=true;
    descriptor.seedPolicy=stochastic?"TestSeed(default20261117), stable particle/operator/event keys; no legacy PIC RNG":"deterministic quadrature/reference/domain checks; no draws";
    descriptor.supportedBuildModes="dependency-light/native identical full-particle material-plane callbacks; native guiding-center source bridge remains unqualified";
    descriptor.stateIsolation="private input/output and staged keyed streams; no source scheduler, native allocation, geometry, density-provider or wave-owner mutation";
    descriptor.callback=callbacks[i];result.push_back(descriptor);
  }
  return result;
}
} }
