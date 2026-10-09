#include "perpendicular_diffusion.h"

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

namespace PD = SEP::PerpendicularDiffusion;

namespace {

constexpr double Pi = 3.141592653589793238462643383279502884;

struct TestState {
  int passed = 0;
  int failed = 0;
};

bool Near(double actual, double expected, double relative,
          double absolute = 1.0e-14) {
  return std::isfinite(actual) &&
      std::fabs(actual - expected) <=
          absolute + relative * std::fabs(expected);
}

void Check(TestState* state, const std::string& name, bool condition) {
  if (condition) {
    ++state->passed;
    std::cout << "PASS " << name << '\n';
  } else {
    ++state->failed;
    std::cerr << "FAIL " << name << '\n';
  }
}

PD::ParticleState Particle() {
  // Synthetic charged particle. Tests compare dimensionless closures, so no
  // species calibration is implied by this deliberately simple state.
  PD::ParticleState p;
  p.massKg = 1.0;
  p.chargeC = 1.0;
  p.momentumKgMPerS = 0.1 * PD::SpeedOfLightMPerS;
  p.mu = 0.25;
  return p;
}

PD::LocalState OrderedState() {
  PD::LocalState local;
  local.meanFieldT = std::array<double,3>{{0.0,0.0,1.0}};
  local.positionM = {{2.0,0.0,0.0}};
  local.turbulence.sampleFingerprint = "synthetic-sample";
  return local;
}

void SupplyParallel(PD::LocalState* local, double kappa) {
  PD::ParallelInput input;
  input.kappaM2PerS = kappa;
  input.modelId = "independent_fixture";
  input.equationVersion = "test-equation";
  input.sampleFingerprint = "synthetic-sample";
  local->parallelDependency = input;
}

PD::ModelConfiguration Configure(const std::string& model,
                                 std::initializer_list<PD::InputParameter> p,
                                 TestState* state) {
  PD::ModelConfiguration configuration;
  const PD::Status status = PD::BuildConfiguration(model, p, &configuration);
  Check(state, "configure/" + model, status.ok());
  return configuration;
}

void TestRegistryAndParser(TestState* state) {
  Check(state,"registry/stable-count",PD::ModelRegistry().size()==46);
  PD::ModelConfiguration configuration;
  const PD::Status gated=PD::BuildConfiguration("qlt_perp",{},&configuration);
  Check(state,"registry/source-gate",gated.code==PD::StatusCode::SourceGate);
  bool everyGate=true;
  for(const PD::ModelDescriptor& descriptor:PD::ModelRegistry())
    if(!descriptor.executable)
      everyGate=everyGate&&
          PD::BuildConfiguration(descriptor.stableId,{},&configuration).code==
              PD::StatusCode::SourceGate;
  Check(state,"registry/all-source-gates",everyGate);
  const PD::Status unknown=PD::BuildConfiguration("not_a_model",{},&configuration);
  Check(state,"registry/unknown",unknown.code==PD::StatusCode::UnsupportedModel);
  const PD::Status duplicate=PD::BuildConfiguration("constant_kappa_perp",
      {{"kappa_perp_m2_per_s","1"},{"kappa_perp_m2_per_s","2"}},
      &configuration);
  Check(state,"parser/duplicate",duplicate.code==PD::StatusCode::InvalidConfiguration);
  const PD::Status unknownKey=PD::BuildConfiguration("constant_kappa_perp",
      {{"kappa_perp_m2_per_s","1"},{"invented","2"}},&configuration);
  Check(state,"parser/unknown-key",unknownKey.code==PD::StatusCode::InvalidConfiguration);
  const PD::ModelResult malformed=PD::Evaluate(Particle(),OrderedState(),
                                               PD::ModelConfiguration{});
  Check(state,"parser/direct-malformed",malformed.status.code==PD::StatusCode::MissingInput&&
      !malformed.coefficients);
  const PD::Status missingCalibration=PD::BuildConfiguration("iso_fit_kuhlen_2025",
      {{"A","1"},{"rho_star","1"},{"s_kappa","1"},{"C_K","1"},
       {"z1_m","1"},{"z2_m","2"}},&configuration);
  Check(state,"parser/missing-calibration",
        missingCalibration.code==PD::StatusCode::MissingCalibration);

  PD::Status active=PD::ConfigureActiveModel("constant_kappa_perp",
      {{"kappa_perp_m2_per_s","2"}});
  const PD::ModelFunction first=PD::ActiveModelFunction;
  active=PD::ConfigureActiveModel("constant_lambda_perp",{{"lambda_perp_m","3"}});
  const PD::ModelFunction second=PD::ActiveModelFunction;
  const PD::Status failedUpdate=PD::ConfigureActiveModel("constant_lambda_perp",
      {{"lambda_perp_m","-1"}});
  Check(state,"manager/model-specific-pointer",active.ok()&&first!=second);
  Check(state,"manager/transactional-failure",!failedUpdate.ok()&&
      PD::ActiveModelFunction==second&&
      PD::GetActiveConfiguration().model==PD::ModelId::ConstantLambdaPerp);
}

void TestSpectra(TestState* state) {
  const double s=5.0/3.0;
  Check(state,"spectrum/C",Near(PD::SpectrumC(s),0.118862354635443,2.0e-14));
  double moment=0.0;
  PD::Status status=PD::SmoothTwoDSpectralMoment(s,3.0,-1,0.8,2.0,&moment);
  // Independent beta-integral simplification for q=3 is kept explicit in
  // the test; it does not call the production length helper.
  const double expected=PD::SpectrumD(s,3.0)/Pi*0.8*2.0*
      std::tgamma(1.5)*std::tgamma(s/2.0)/std::tgamma((s+3.0)/2.0);
  Check(state,"spectrum/moment-Iminus1",status.ok()&&Near(moment,expected,3.0e-14));
  status=PD::SmoothTwoDSpectralMoment(s,0.0,-2,1.0,1.0,&moment);
  Check(state,"spectrum/divergent-moment",status.code==PD::StatusCode::DivergentMoment);
  double ultra=0.0,integral=0.0;
  status=PD::SmoothSpectrumLengths(s,3.0,1.0,&ultra,&integral);
  Check(state,"spectrum/lengths",status.ok()&&
      Near(ultra,0.577350269189626,2.0e-14)&&
      Near(integral,0.746834200222187,2.0e-14));
}

void TestFieldLines(TestState* state) {
  PD::LocalState local=OrderedState();
  local.turbulence.geometry=PD::GeometryKind::CompositeSlab2D;
  local.turbulence.slabVarianceT2=0.2;
  local.turbulence.twoDVarianceT2=0.8;
  local.turbulence.slabBendoverLengthM=1.0;
  local.turbulence.twoDBendoverLengthM=1.0;
  local.turbulence.inertialIndex=5.0/3.0;
  local.turbulence.energyRangeIndex=3.0;
  PD::ModelConfiguration configuration=Configure("fl_composite",{},state);
  const PD::ModelResult result=PD::Evaluate(Particle(),local,configuration);
  Check(state,"field-line/F4",result.status.ok()&&result.fieldLineM&&
      Near(*result.fieldLineM,0.404394481,3.0e-9));

  local.turbulence.energyRangeIndex=1.0;
  const PD::ModelResult divergent=PD::Evaluate(Particle(),local,configuration);
  Check(state,"field-line/q-gate",divergent.status.code==PD::StatusCode::DivergentMoment&&
      !divergent.fieldLineM);

  PD::LocalState zero=OrderedState();
  zero.turbulence.geometry=PD::GeometryKind::Pure2D;
  zero.turbulence.twoDVarianceT2=0.0;
  PD::ModelConfiguration twoD=Configure("fl_2d",{},state);
  const PD::ModelResult zeroResult=PD::Evaluate(Particle(),zero,twoD);
  Check(state,"field-line/zero-without-undefined-shape",zeroResult.status.ok()&&
      zeroResult.fieldLineM&&*zeroResult.fieldLineM==0.0);
}

void TestPrescribedAndTensor(TestState* state) {
  PD::LocalState local=OrderedState();
  SupplyParallel(&local,8.0);
  PD::ModelConfiguration constant=Configure("constant_kappa_perp",
      {{"kappa_perp_m2_per_s","2.5"}},state);
  const PD::ModelResult result=PD::Evaluate(Particle(),local,constant);
  Check(state,"prescribed/constant",result.status.ok()&&result.coefficients&&
      Near(*result.coefficients->perpendicular1M2PerS,2.5,0.0)&&
      Near(*result.coefficients->parallelM2PerS,8.0,0.0));
  std::array<std::array<double,3>,3> tensor;
  const PD::Status assembled=PD::AssembleSymmetricTensor(*result.coefficients,&tensor);
  Check(state,"tensor/axisymmetric",assembled.ok()&&Near(tensor[0][0],2.5,1.0e-14)&&
      Near(tensor[1][1],2.5,1.0e-14)&&Near(tensor[2][2],8.0,1.0e-14));

  PD::ModelConfiguration pitch=Configure("pitch_angle_perp",
      {{"D0_m2_per_s","3"},{"shape","abs_mu"}},state);
  const PD::ModelResult pitchResult=PD::Evaluate(Particle(),local,pitch);
  Check(state,"prescribed/pitch-angle",pitchResult.status.ok()&&
      pitchResult.pitchAngleM2PerS&&Near(*pitchResult.pitchAngleM2PerS,1.5,0.0));

  PD::ModelConfiguration ratio=Configure("ratio_kappa",{{"eta_kappa","0.125"}},state);
  const PD::ModelResult ratioResult=PD::Evaluate(Particle(),local,ratio);
  Check(state,"prescribed/ratio",ratioResult.status.ok()&&
      Near(*ratioResult.coefficients->perpendicular1M2PerS,1.0,1.0e-14));
}

void TestKernels(TestState* state) {
  double exact=0.0,rational=0.0;
  PD::Status status=PD::EvaluateImplicitSlabKernel(1.0,false,&exact);
  PD::Status status2=PD::EvaluateImplicitSlabKernel(1.0,true,&rational);
  Check(state,"kernel/U6",status.ok()&&Near(exact,0.242127843858688,2.0e-13));
  Check(state,"kernel/U8-distinct",status2.ok()&&Near(rational,1.0/3.0,1.0e-15)&&
      !Near(exact,rational,1.0e-2));
  status=PD::EvaluateImplicitSlabKernel(1.0e4,false,&exact);
  const double asymptotic=1.0/(2.0e8)-3.0/(4.0e16)+15.0/(8.0e24);
  Check(state,"kernel/U7-large",status.ok()&&Near(exact,asymptotic,2.0e-15,1.0e-25));
}

void TestDirectClosures(TestState* state) {
  PD::ParticleState particle=Particle();
  PD::ParticleKinematics kin;
  PD::ComputeParticleKinematics(particle,&kin);

  // The chosen B0=ell=v-normalized states are exactly the dimensionless
  // Section 17 fixtures after kappa_parallel=v*lambda_parallel/3.
  PD::LocalState flpd=OrderedState();
  flpd.turbulence.geometry=PD::GeometryKind::Pure2D;
  flpd.turbulence.twoDVarianceT2=1.0;
  flpd.turbulence.twoDBendoverLengthM=1.0;
  flpd.turbulence.inertialIndex=5.0/3.0;
  flpd.turbulence.energyRangeIndex=3.0;
  SupplyParallel(&flpd,kin.speedMPerS*0.1/3.0);
  PD::ModelConfiguration flpdConfig=Configure("flpd_complete",{},state);
  const PD::ModelResult flpdResult=PD::Evaluate(particle,flpd,flpdConfig);
  const double flpdRatio=flpdResult.status.ok()?
      *flpdResult.coefficients->perpendicular1M2PerS/
           *flpdResult.coefficients->parallelM2PerS:0.0;
  const bool flpdPass=flpdResult.status.ok()&&
      Near(flpdRatio,
           0.0898707818,3.0e-7);
  if(!flpdPass) std::cerr<<"  FLPD actual="<<flpdRatio<<" status="
                         <<flpdResult.status.detail<<'\n';
  Check(state,"closure/FLPD-D2",flpdPass);

  PD::LocalState enlgc=OrderedState();
  enlgc.turbulence.geometry=PD::GeometryKind::Pure2D;
  enlgc.turbulence.twoDVarianceT2=0.8;
  enlgc.turbulence.twoDBendoverLengthM=1.0;
  enlgc.turbulence.inertialIndex=5.0/3.0;
  enlgc.turbulence.energyRangeIndex=0.0;
  SupplyParallel(&enlgc,kin.speedMPerS/3.0);
  PD::ModelConfiguration enlgcConfig=Configure("enlgc_2d",{},state);
  const PD::ModelResult enlgcResult=PD::Evaluate(particle,enlgc,enlgcConfig);
  const double enlgcLambda=3.0**enlgcResult.coefficients->perpendicular1M2PerS/
                           kin.speedMPerS;
  const bool enlgcPass=enlgcResult.status.ok()&&Near(enlgcLambda,0.263129509,3.0e-7);
  if(!enlgcPass) std::cerr<<"  ENLGC actual="<<enlgcLambda<<" status="
                          <<enlgcResult.status.detail<<'\n';
  Check(state,"closure/ENLGC-N5",enlgcPass);

  PD::LocalState rbd=OrderedState();
  rbd.turbulence.geometry=PD::GeometryKind::Pure2D;
  rbd.turbulence.twoDVarianceT2=0.4;
  rbd.turbulence.twoDBendoverLengthM=1.0;
  SupplyParallel(&rbd,kin.speedMPerS*100.0/3.0);
  PD::ModelConfiguration rbdConfig=Configure("rbd_bc",
      {{"a_squared","0.3333333333333333333"},
       {"area_energy_index","2"},{"area_inertial_index","1.6666666666666667"}},state);
  const PD::ModelResult rbdResult=PD::Evaluate(particle,rbd,rbdConfig);
  const double rbdValue=rbdResult.status.ok()?
      *rbdResult.coefficients->perpendicular1M2PerS/kin.speedMPerS:0.0;
  const bool rbdPass=rbdResult.status.ok()&&Near(rbdValue,0.0946606068,4.0e-7);
  if(!rbdPass) std::cerr<<"  RBD actual="<<rbdValue<<" status="
                        <<rbdResult.status.detail<<'\n';
  Check(state,"closure/RBD-B3",rbdPass);

  // Every implicit closure must take its analytic zero-parallel branch before
  // constructing Q or a logarithm. Geometry is still validated first.
  SupplyParallel(&flpd,0.0);
  const PD::ModelResult zero=PD::Evaluate(particle,flpd,flpdConfig);
  Check(state,"closure/zero-parallel",zero.status.ok()&&
      *zero.coefficients->perpendicular1M2PerS==0.0);

  PD::LocalState composite=OrderedState();
  composite.turbulence.geometry=PD::GeometryKind::CompositeSlab2D;
  composite.turbulence.slabVarianceT2=0.2;
  composite.turbulence.twoDVarianceT2=0.8;
  composite.turbulence.slabBendoverLengthM=1.0;
  composite.turbulence.twoDBendoverLengthM=1.0;
  composite.turbulence.inertialIndex=5.0/3.0;
  composite.turbulence.energyRangeIndex=3.0;
  SupplyParallel(&composite,kin.speedMPerS/3.0);
  PD::ModelConfiguration nlgc=Configure("nlgc",{{"a_squared","0.3333333333333333"}},state);
  const PD::ModelResult nlgcResult=PD::Evaluate(particle,composite,nlgc);
  const double nlgcBound=composite.parallelDependency->kappaM2PerS*
      (0.2+0.8)/2.0/3.0;
  Check(state,"closure/NLGC-bound",nlgcResult.status.ok()&&
      *nlgcResult.coefficients->perpendicular1M2PerS>0.0&&
      *nlgcResult.coefficients->perpendicular1M2PerS<=nlgcBound*(1.0+1.0e-8));

  PD::ModelConfiguration unlt=Configure("unlt",{{"a_squared","1"}},state);
  const PD::ModelResult unltResult=PD::Evaluate(particle,composite,unlt);
  Check(state,"closure/UNLT-bound",unltResult.status.ok()&&
      *unltResult.coefficients->perpendicular1M2PerS>0.0&&
      *unltResult.coefficients->perpendicular1M2PerS<=
          composite.parallelDependency->kappaM2PerS*0.8/2.0*(1.0+1.0e-8));

  PD::ModelConfiguration exact=Configure("implicit_slab_exact_2016",{},state);
  PD::ModelConfiguration rational=Configure("implicit_slab_rational_2016",{},state);
  const PD::ModelResult exactResult=PD::Evaluate(particle,composite,exact);
  const PD::ModelResult rationalResult=PD::Evaluate(particle,composite,rational);
  Check(state,"closure/implicit-distinct",exactResult.status.ok()&&rationalResult.status.ok()&&
      *exactResult.coefficients->perpendicular1M2PerS>0.0&&
      *rationalResult.coefficients->perpendicular1M2PerS>0.0&&
      !Near(*exactResult.coefficients->perpendicular1M2PerS,
            *rationalResult.coefficients->perpendicular1M2PerS,1.0e-6));

  PD::LocalState slab=OrderedState();
  slab.turbulence.geometry=PD::GeometryKind::PureSlab;
  slab.turbulence.slabVarianceT2=0.2;
  slab.turbulence.slabBendoverLengthM=1.0;
  slab.turbulence.inertialIndex=5.0/3.0;
  SupplyParallel(&slab,kin.speedMPerS/3.0);
  PD::ModelConfiguration slabDiagnostic=Configure("nlgc_slab_kernel_diagnostic",
      {{"a_squared","0.3333333333333333"}},state);
  const PD::ModelResult diagnostic=PD::Evaluate(particle,slab,slabDiagnostic);
  Check(state,"closure/NLGC-slab-diagnostic",diagnostic.status.ok()&&
      diagnostic.quality==PD::Quality::Diagnostic&&
      *diagnostic.coefficients->perpendicular1M2PerS>0.0);
  const PD::ModelResult prohibited=PD::Evaluate(particle,slab,nlgc);
  Check(state,"closure/NLGC-pure-slab-rejected",
        prohibited.status.code==PD::StatusCode::IncompatibleGeometry&&
        !prohibited.coefficients);
}

void TestClosedAndDiagnostics(TestState* state) {
  PD::ParticleState particle=Particle(); PD::ParticleKinematics kin;
  PD::ComputeParticleKinematics(particle,&kin);
  PD::LocalState local=OrderedState();
  SupplyParallel(&local,kin.speedMPerS/3.0);
  local.fieldLineCoefficientM=0.1;
  local.fieldLineModelId="fixture";
  local.fieldLineSampleFingerprint="synthetic-sample";
  PD::ModelConfiguration closed=Configure("composite_closed_2019",
      {{"perpendicular_length_m","1"},
       {"length_profile","parameterized_tagged_length"}},state);
  const PD::ModelResult result=PD::Evaluate(particle,local,closed);
  const double lambda=3.0**result.coefficients->perpendicular1M2PerS/kin.speedMPerS;
  Check(state,"closed/B6",result.status.ok()&&Near(lambda,0.00885427379,3.0e-8));

  PD::ModelConfiguration compound=Configure("compound_diffusive_lines",
      {{"age_s","4"}},state);
  const PD::ModelResult anomalous=PD::Evaluate(particle,local,compound);
  const double expected=4.0*0.1*std::sqrt((*local.parallelDependency).kappaM2PerS*4.0/Pi);
  Check(state,"diagnostic/F8",anomalous.status.ok()&&anomalous.particleMoments&&
      anomalous.diffusionRegime==PD::DiffusionRegime::NoNormalDiffusion&&
      Near(anomalous.particleMoments->rawMsdM2[0],expected,2.0e-14));
}

void TestPairAndTable(TestState* state) {
  PD::LocalState local=OrderedState();
  local.turbulence.geometry=PD::GeometryKind::CompositeSlab2D;
  local.turbulence.slabVarianceT2=0.2;
  local.turbulence.twoDVarianceT2=0.8;
  local.turbulence.slabBendoverLengthM=1.0;
  local.turbulence.twoDBendoverLengthM=0.1;
  PD::ParticleState particle=Particle();
  particle.momentumKgMPerS=0.035848041610665;
  PD::ModelConfiguration pair=Configure("nlgce_f_2014",{},state);
  const PD::ModelResult pairResult=PD::Evaluate(particle,local,pair);
  Check(state,"paired/NLGCE-owner",pairResult.status.ok()&&pairResult.coefficients&&
      pairResult.coefficients->parallelOwner==PD::DependencyOwner::PairedBackend&&
      pairResult.coefficients->parallelM2PerS&&
      *pairResult.coefficients->parallelM2PerS>0.0&&
      *pairResult.coefficients->perpendicular1M2PerS>0.0);

  PD::ModelConfiguration table=Configure("tabulated_perp",
      {{"axis","time"},{"axis_values_SI","0,10,20"},
       {"values_m2_per_s","2,4,8"},{"generation_identity","fixture-v1"},
       {"table_checksum","fixture-sha256"},{"boundary_policy","reject"},
       {"interpolation","linear"},{"zero_policy","linear_explicit_zero"}},state);
  PD::LocalState tabState=OrderedState(); tabState.timeS=5.0;
  const PD::ModelResult tableResult=PD::Evaluate(Particle(),tabState,table);
  Check(state,"table/linear",tableResult.status.ok()&&
      Near(*tableResult.coefficients->perpendicular1M2PerS,3.0,1.0e-14));
  tabState.timeS=21.0;
  const PD::ModelResult outside=PD::Evaluate(Particle(),tabState,table);
  Check(state,"table/reject-boundary",outside.status.code==PD::StatusCode::OutsideModelDomain&&
      !outside.coefficients);
}

void TestFramesAndIsotropy(TestState* state) {
  PD::LocalState unequal=OrderedState();
  SupplyParallel(&unequal,10.0);
  unequal.heliocentricColatitudeRad=Pi/2.0;
  unequal.perpendicularAxis1=std::array<double,3>{{1.0,0.0,0.0}};
  unequal.perpendicularAxis2=std::array<double,3>{{0.0,1.0,0.0}};
  PD::ModelConfiguration nwu=Configure("nwu_ratio_polar",
      {{"eta_r","0.02"},{"eta_theta","0.01"},
       {"polar_enhancement","3"},{"theta_F_rad","0.6"},
       {"width_per_rad","2"},{"angular_convention","folded_radian"}},state);
  const PD::ModelResult oriented=PD::Evaluate(Particle(),unequal,nwu);
  Check(state,"frame/oriented-unequal",oriented.status.ok()&&oriented.coefficients&&
      oriented.coefficients->frame.kind==PD::FrameKind::OrientedUnequal&&
      *oriented.coefficients->perpendicular1M2PerS!=
      *oriented.coefficients->perpendicular2M2PerS);

  PD::LocalState isotropic=OrderedState();
  isotropic.meanFieldT=std::array<double,3>{{0.0,0.0,0.0}};
  isotropic.turbulence.geometry=PD::GeometryKind::Isotropic3D;
  isotropic.turbulence.totalVarianceT2=1.0;
  PD::ParticleState particle=Particle();
  particle.momentumKgMPerS=0.01;
  PD::ModelConfiguration kappa0=Configure("iso_kappa0_snodin_2016",
      {{"outer_scale_m","1"},{"spectrum","sharp_kolmogorov"},
       {"domain_policy","reject"}},state);
  const PD::ModelResult isoResult=PD::Evaluate(particle,isotropic,kappa0);
  Check(state,"frame/zero-mean-isotropic",isoResult.status.ok()&&
      isoResult.coefficients&&
      isoResult.coefficients->frame.kind==PD::FrameKind::Isotropic&&
      !isoResult.coefficients->frame.b&&isoResult.coefficients->isotropicM2PerS);
}

void TestFitsAndDrifts(TestState* state) {
  PD::LocalState isotropic=OrderedState();
  isotropic.turbulence.geometry=PD::GeometryKind::Isotropic3D;
  isotropic.turbulence.totalVarianceT2=1.0;
  PD::ParticleState candiaParticle=Particle();
  candiaParticle.momentumKgMPerS=0.1;
  PD::ModelConfiguration candia=Configure("iso_fit_candia_roulet_2004",
      {{"outer_scale_m","1"},{"spectrum","kolmogorov"},
       {"domain_policy","reject"}},state);
  const PD::ModelResult candiaResult=PD::Evaluate(candiaParticle,isotropic,candia);
  Check(state,"fit/Candia-pair",candiaResult.status.ok()&&candiaResult.coefficients&&
      candiaResult.coefficients->parallelOwner==PD::DependencyOwner::PairedBackend&&
      *candiaResult.coefficients->parallelM2PerS>0.0);

  PD::ParticleState snodinParticle=Particle();
  snodinParticle.momentumKgMPerS=std::sqrt(2.0)*0.01;
  PD::ModelConfiguration snodin=Configure("iso_fit_snodin_2016",
      {{"outer_scale_m","1"},{"domain_policy","reject"}},state);
  const PD::ModelResult snodinResult=PD::Evaluate(snodinParticle,isotropic,snodin);
  Check(state,"fit/Snodin-pair",snodinResult.status.ok()&&snodinResult.coefficients&&
      *snodinResult.coefficients->parallelM2PerS>0.0&&
      *snodinResult.coefficients->perpendicular1M2PerS>0.0);

  PD::LocalState kuhlenState=OrderedState();
  kuhlenState.turbulence.geometry=PD::GeometryKind::Isotropic3D;
  kuhlenState.turbulence.totalVarianceT2=0.25;
  kuhlenState.turbulence.correlationLengthM=1.0;
  PD::ParticleState kuhlenParticle=Particle();
  kuhlenParticle.momentumKgMPerS=std::sqrt(1.25)*0.01;
  PD::ModelConfiguration kuhlen=Configure("iso_fit_kuhlen_2025",
      {{"A","1"},{"rho_star","0.5"},{"s_kappa","1"},{"C_K","0.5"},
       {"z1_m","1.5"},{"z2_m","5"},{"gamma_K","-0.5"},
       {"transverse_correlation_length_m","0.1"},
       {"calibration_id","synthetic-equation-test"},
       {"root_selection","first_upward_crossing"}},state);
  const PD::ModelResult kuhlenResult=PD::Evaluate(kuhlenParticle,kuhlenState,kuhlen);
  Check(state,"fit/Kuhlen-calibrated-algorithm",kuhlenResult.status.ok()&&
      kuhlenResult.coefficients&&*kuhlenResult.coefficients->perpendicular1M2PerS>0.0&&
      kuhlenResult.numerical.controls.at("decorrelation_time_s")>0.0);

  PD::LocalState classicalState=OrderedState();
  SupplyParallel(&classicalState,10.0);
  PD::ModelConfiguration classical=Configure("classical_scattering",{},state);
  const PD::ModelResult classicalResult=PD::Evaluate(Particle(),classicalState,classical);
  Check(state,"fit/classical-pair-and-Hall",classicalResult.status.ok()&&
      classicalResult.coefficients&&classicalResult.signedHallM2PerS&&
      *classicalResult.signedHallM2PerS>0.0);

  PD::ModelConfiguration drift=Configure("drift_rigidity_reduction",
      {{"K_A0","0.8"},{"rigidity_A_V","1e15"}},state);
  const PD::ModelResult positive=PD::Evaluate(Particle(),classicalState,drift);
  PD::ParticleState negativeParticle=Particle(); negativeParticle.chargeC=-1.0;
  const PD::ModelResult negative=PD::Evaluate(negativeParticle,classicalState,drift);
  Check(state,"drift/charge-sign",positive.status.ok()&&negative.status.ok()&&
      positive.signedHallM2PerS&&negative.signedHallM2PerS&&
      Near(*positive.signedHallM2PerS,-*negative.signedHallM2PerS,1.0e-14));
}

}  // namespace

int main() {
  TestState state;
  TestRegistryAndParser(&state);
  TestSpectra(&state);
  TestFieldLines(&state);
  TestPrescribedAndTensor(&state);
  TestKernels(&state);
  TestDirectClosures(&state);
  TestClosedAndDiagnostics(&state);
  TestPairAndTable(&state);
  TestFramesAndIsotropy(&state);
  TestFitsAndDrifts(&state);
  std::cout << "perpendicular_diffusion: " << state.passed << " passed, "
            << state.failed << " failed\n";
  return state.failed==0 ? 0 : 1;
}
