#include "parallel_diffusion.h"

#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <string>
#include <thread>
#include <vector>

namespace PD = SEP::ParallelDiffusion;

namespace {

constexpr double ElementaryChargeC = 1.602176634e-19;
constexpr double ProtonMassKg = 1.67262192369e-27;
constexpr double FieldT = 1.0e-9;
constexpr double SlabLengthM = 1.0e9;

double Number(const char* text) {
  char* end = nullptr;
  const double value = std::strtod(text, &end);
  if (!end || *end != '\0' || !std::isfinite(value)) {
    std::cerr << "invalid numeric argument: " << text << '\n';
    std::exit(2);
  }
  return value;
}

PD::ParticleState Particle(double rLOverSlabLength) {
  PD::ParticleState particle;
  particle.massKg = ProtonMassKg;
  particle.chargeC = ElementaryChargeC;
  particle.momentumKgMPerS = rLOverSlabLength * SlabLengthM *
      ElementaryChargeC * FieldT;
  particle.nucleonCount = 1.0;
  return particle;
}

PD::LocalState State(double fs, double eps2, double rho) {
  PD::LocalState state;
  state.meanFieldT = std::array<double, 3>{{0.0, 0.0, FieldT}};
  state.turbulence.emplace();
  state.turbulence->slabVarianceT2 = fs * eps2 * FieldT * FieldT;
  state.turbulence->twoDVarianceT2 = (1.0 - fs) * eps2 * FieldT * FieldT;
  state.turbulence->slabBendoverLengthM = SlabLengthM;
  state.turbulence->twoDBendoverLengthM = SlabLengthM / rho;
  state.turbulence->inertialIndex = 5.0 / 3.0;
  return state;
}

int Print(const PD::ParallelResult& result) {
  if (!result.status.ok()) std::cerr << result.status.detail << '\n';
  std::cout << std::setprecision(17) << static_cast<int>(result.status.code);
  if (result.lambdaParallelM.has_value())
    std::cout << ' ' << *result.lambdaParallelM / SlabLengthM;
  else
    std::cout << " nan";
  if (result.perpendicular.has_value())
    std::cout << ' ' << result.perpendicular->lambdaPerpendicularM /
        SlabLengthM;
  else
    std::cout << " nan";
  if (result.nonlinearMaxLogResidual.has_value())
    std::cout << ' ' << *result.nonlinearMaxLogResidual;
  else
    std::cout << " nan";
  std::cout << '\n';
  return result.status.ok() ? 0 : 1;
}

}  // namespace

int main(int argc, char** argv) {
  if (argc < 2) return 2;
  const std::string mode = argv[1];
  if (mode == "selfcheck" && argc == 2) {
    bool pass = true;
    const PD::ParticleState particle = Particle(0.01);

    // Equation (23) must reproduce the target lambda even though D0 is
    // recomputed from the independently integrated shape normalization.
    PD::ModelConfiguration shape;
    shape.model = PD::ModelId::PrescribedLambdaMuShape;
    shape.prescribedLambdaMu.amplitudeMode =
        PD::PitchAngleAmplitudeMode::TargetLambda;
    shape.prescribedLambdaMu.qMu = 5.0 / 3.0;
    shape.prescribedLambdaMu.hMu = 0.01;
    shape.prescribedLambdaMu.targetLambdaM = 2.5e9;
    const PD::ParallelResult shaped = PD::Evaluate(particle, PD::LocalState(), shape);
    double dMu = 0.0;
    const PD::Status pitch = PD::EvaluatePitchAngleDiffusion(
        0.4, particle, PD::LocalState(), shape, &dMu);
    pass = pass && shaped.status.ok() && pitch.ok() && dMu > 0.0 &&
        std::fabs(*shaped.lambdaParallelM / 2.5e9 - 1.0) < 1.0e-12;

    // A constant canonical spectrum over a finite interval with explicitly
    // zero tails integrates analytically to P*(k_max-k_min). The evaluator
    // accepts the matching declaration and distinguishes a mismatched one as
    // inconsistent spectrum rather than missing resonance coverage.
    PD::ModelConfiguration supplied;
    supplied.model = PD::ModelId::QltSlabSpectrum;
    supplied.qltSlab.spectrum.form = PD::SpectrumForm::SuppliedLogLog;
    supplied.qltSlab.spectrum.wavenumberRadPerM = {1.0e-8, 1.0e-6};
    supplied.qltSlab.spectrum.powerT2M = {1.0e-12, 1.0e-12};
    supplied.qltSlab.spectrum.declaredVarianceT2 = 9.9e-19;
    supplied.qltSlab.spectrum.sourceIdentity = "constant_fixture";
    supplied.qltSlab.spectrum.lowKPolicy = PD::SpectrumTailPolicy::Zero;
    supplied.qltSlab.spectrum.highKPolicy = PD::SpectrumTailPolicy::Zero;
    PD::LocalState suppliedState;
    suppliedState.meanFieldT = std::array<double, 3>{{0.0, 0.0, FieldT}};
    pass = pass && PD::EvaluatePitchAngleDiffusion(
        0.5, particle, suppliedState, supplied, &dMu).ok() && dMu > 0.0;
    supplied.qltSlab.spectrum.declaredVarianceT2 *= 2.0;
    pass = pass && PD::EvaluatePitchAngleDiffusion(
        0.5, particle, suppliedState, supplied, &dMu).code ==
        PD::StatusCode::InconsistentSpectrum;

    // Equation (35), including its explicit s=1 logarithmic normalization,
    // is checked against arithmetic assembled independently of the production
    // spectrum evaluator.  This point resonates in the dissipation branch.
    PD::ModelConfiguration multirange;
    multirange.model = PD::ModelId::QltSlabSpectrum;
    multirange.qltSlab.spectrum.form = PD::SpectrumForm::Multirange;
    multirange.qltSlab.spectrum.energyRangeIndex = 0.0;
    multirange.qltSlab.spectrum.dissipationIndex = 1.5;
    multirange.qltSlab.spectrum.dissipationWavenumberRadPerM = 1.0e-8;
    multirange.qltSlab.numerical.relativeTolerance = 1.0e-7;
    multirange.qltSlab.numerical.maximumRefinements = 24;
    PD::LocalState multirangeState = State(1.0, 0.04, 1.0);
    multirangeState.turbulence->inertialIndex = 1.0;
    PD::ParticleKinematics multirangeKinematics;
    PD::ComputeParticleKinematics(particle, &multirangeKinematics);
    const double mu = 0.5;
    const double xD = 10.0;
    const double xResonant = 1.0 / (0.01 * mu);
    const double normalization = 1.0 + std::log(xD) + 2.0;
    const double shapeAtResonance =
        std::pow(xD, 0.5) * std::pow(xResonant, -1.5);
    const double canonicalPower = 0.04 * FieldT * FieldT * SlabLengthM *
        shapeAtResonance / normalization;
    const double gyrofrequency = ElementaryChargeC * FieldT /
        (multirangeKinematics.gamma * ProtonMassKg);
    const double expectedDmu = std::acos(-1.0) * gyrofrequency *
        gyrofrequency * (1.0 - mu * mu) * canonicalPower /
        (4.0 * FieldT * FieldT * multirangeKinematics.speedMPerS * mu);
    pass = pass && PD::EvaluatePitchAngleDiffusion(
        mu, particle, multirangeState, multirange, &dMu).ok() &&
        std::fabs(dMu / expectedDmu - 1.0) < 2.0e-14;
    pass = pass && PD::Evaluate(particle, multirangeState, multirange).status.ok();
    multirange.qltSlab.spectrum.dissipationWavenumberRadPerM = 5.0e-10;
    pass = pass && PD::EvaluatePitchAngleDiffusion(
        mu, particle, multirangeState, multirange, &dMu).code ==
        PD::StatusCode::InvalidBackground;
    multirange.qltSlab.spectrum.dissipationWavenumberRadPerM = 1.0e-8;
    multirange.qltSlab.spectrum.dissipationIndex = 2.0;
    pass = pass && PD::Evaluate(particle, multirangeState, multirange).status.code ==
        PD::StatusCode::InfiniteMeanFreePath;

    // Section 8.7 conversion factors must all reproduce the same canonical
    // one-sided total-transverse sample.  The deliberately unequal component
    // and signed-k values also verify that the general converters do not
    // assume axisymmetry or evenness.  Coordinate mappings independently
    // check the cycles/radian and frozen-flow Jacobians.
    double canonical = 0.0, mappedK = 0.0;
    const double canonicalFixture = 12.0;
    pass = pass && PD::ConvertOneSidedComponentsToCanonical(
        5.0, 7.0, &canonical).ok() && canonical == canonicalFixture;
    pass = pass && PD::ConvertTwoSidedComponentsToCanonical(
        2.0, 3.0, 3.0, 4.0, &canonical).ok() &&
        canonical == canonicalFixture;
    pass = pass && PD::ConvertEvenTwoSidedTotalToCanonical(
        6.0, &canonical).ok() && canonical == canonicalFixture;
    const double twoPi = 2.0 * std::acos(-1.0);
    pass = pass && PD::ConvertOneSidedCyclesPerMToCanonical(
        2.0, twoPi * canonicalFixture, &mappedK, &canonical).ok() &&
        std::fabs(mappedK / (2.0 * twoPi) - 1.0) < 1.0e-15 &&
        std::fabs(canonical / canonicalFixture - 1.0) < 1.0e-15;
    PD::FrozenFlowMapping frozenFlow;
    frozenFlow.samplingVelocityProjectionMPerS = 4.0;
    frozenFlow.assumptionIdentity = "mathematical_frozen_flow_fixture";
    pass = pass && PD::ConvertFrozenFlowFrequencyToCanonical(
        2.0, twoPi * canonicalFixture / 4.0, frozenFlow,
        &mappedK, &canonical).ok() &&
        std::fabs(mappedK / std::acos(-1.0) - 1.0) < 1.0e-15 &&
        std::fabs(canonical / canonicalFixture - 1.0) < 1.0e-15;
    pass = pass && PD::ConvertQinZhangSlabComponentToCanonical(
        3.0, &canonical).ok() && canonical == canonicalFixture;
    pass = pass && PD::ConvertQinZhangTwoDComponentToCanonical(
        3.0, &canonical).ok() && canonical == canonicalFixture;
    frozenFlow.assumptionIdentity.clear();
    pass = pass && PD::ConvertFrozenFlowFrequencyToCanonical(
        2.0, 1.0, frozenFlow, &mappedK, &canonical).code ==
        PD::StatusCode::InvalidConfiguration;

    // Equation (57b): an unhalved Elsasser sum is twice canonical Z^2;
    // converting it must equal the direct kinetic-plus-magnetic convention.
    double direct = 0.0, halfElsasser = 0.0, elsasser = 0.0;
    double specificEnergy = 0.0;
    pass = pass && PD::ConvertTurbulenceMoment(
        4.0, 5.0e-20, -0.2, 1.25663706212e-6,
        PD::MomentEnergyConvention::KineticPlusMagneticVariance,
        PD::ResidualEnergyConvention::KineticMinusMagnetic, &direct).ok();
    pass = pass && PD::ConvertTurbulenceMoment(
        8.0, 5.0e-20, 0.2, 1.25663706212e-6,
        PD::MomentEnergyConvention::ElsasserSum,
        PD::ResidualEnergyConvention::MagneticMinusKinetic, &elsasser).ok();
    pass = pass && PD::ConvertTurbulenceMoment(
        4.0, 5.0e-20, -0.2, 1.25663706212e-6,
        PD::MomentEnergyConvention::HalfElsasserSum,
        PD::ResidualEnergyConvention::KineticMinusMagnetic,
        &halfElsasser).ok();
    pass = pass && PD::ConvertTurbulenceMoment(
        2.0, 5.0e-20, -0.2, 1.25663706212e-6,
        PD::MomentEnergyConvention::SpecificTotalFluctuationEnergy,
        PD::ResidualEnergyConvention::KineticMinusMagnetic,
        &specificEnergy).ok();
    pass = pass && std::fabs(direct / elsasser - 1.0) < 1.0e-15 &&
        std::fabs(direct / halfElsasser - 1.0) < 1.0e-15 &&
        std::fabs(direct / specificEnergy - 1.0) < 1.0e-15;

    // A local-state moment provider is read for every call.  Doubling its
    // canonical energy at fixed density/residual/slab fraction doubles the
    // slab variance and therefore halves exact QLT lambda; a stale adapter
    // snapshot or hidden cache would fail this identity.
    PD::ModelConfiguration momentAdapter;
    momentAdapter.model = PD::ModelId::TurbulenceAdapter;
    momentAdapter.turbulenceAdapter.momentSource = PD::MomentSource::LocalState;
    momentAdapter.turbulenceAdapter.energyConvention =
        PD::MomentEnergyConvention::HalfElsasserSum;
    momentAdapter.turbulenceAdapter.residualConvention =
        PD::ResidualEnergyConvention::KineticMinusMagnetic;
    momentAdapter.turbulenceAdapter.vacuumPermeabilityHPerM =
        1.25663706212e-6;
    momentAdapter.turbulenceAdapter.closure =
        PD::AdapterClosure::QltSlabSpectrum;
    PD::LocalState momentState = State(1.0, 0.04, 1.0);
    momentState.turbulence->densityKgPerM3 = 5.0e-20;
    momentState.turbulence->providerMomentM2PerS2 = 4.0;
    momentState.turbulence->residualEnergy = -0.2;
    momentState.turbulence->slabFraction = 1.0;
    const PD::ParallelResult firstMoment =
        PD::Evaluate(particle, momentState, momentAdapter);
    momentState.turbulence->providerMomentM2PerS2 = 8.0;
    const PD::ParallelResult secondMoment =
        PD::Evaluate(particle, momentState, momentAdapter);
    pass = pass && firstMoment.status.ok() && secondMoment.status.ok() &&
        std::fabs(*firstMoment.lambdaParallelM /
                  *secondMoment.lambdaParallelM - 2.0) < 1.0e-13 &&
        firstMoment.provenance.requestedModelId == "turbulence_adapter" &&
        firstMoment.provenance.evaluatedModelId == "qlt_slab_spectrum";

    // A two-dimensional table checks row-major multilinear log interpolation.
    PD::ModelConfiguration table;
    table.model = PD::ModelId::TabulatedParallel;
    table.table.storedCoefficient = PD::StoredCoefficient::LambdaParallel;
    table.table.axes = {PD::TableAxis::Rigidity, PD::TableAxis::Time};
    PD::ParticleKinematics kin;
    PD::ComputeParticleKinematics(particle, &kin);
    table.table.axisSI = {{kin.rigidityV / 2.0, kin.rigidityV * 2.0},
                          {0.0, 10.0}};
    table.table.coefficientSI = {1.0e9, 4.0e9, 4.0e9, 16.0e9};
    table.table.generationIdentity = "mathematical_fixture";
    table.table.timeInterpolation = PD::TimeInterpolation::Linear;
    PD::LocalState tableState;
    tableState.timeS = 5.0;
    const PD::ParallelResult tabulated = PD::Evaluate(particle, tableState, table);
    pass = pass && tabulated.status.ok() &&
        std::fabs(*tabulated.lambdaParallelM / 5.0e9 - 1.0) < 1.0e-14;
    table.table.timeInterpolation = PD::TimeInterpolation::StepPrevious;
    const PD::ParallelResult stepped = PD::Evaluate(particle, tableState, table);
    pass = pass && stepped.status.ok() &&
        std::fabs(*stepped.lambdaParallelM / 2.0e9 - 1.0) < 1.0e-14;
    tableState.timeS = 11.0;
    pass = pass && PD::Evaluate(particle, tableState, table).status.code ==
        PD::StatusCode::OutsideModelDomain;

    // Batch shape errors are transactional; a valid batch preserves each
    // point's individual failure rather than failing unrelated points.
    std::vector<PD::ParallelResult> batch(1);
    pass = pass && PD::EvaluateBatch({particle}, {}, shape, &batch).code ==
        PD::StatusCode::InvalidConfiguration && batch.size() == 1;
    PD::ParticleState invalid = particle;
    invalid.chargeC = 0.0;
    pass = pass && PD::EvaluateBatch({particle, invalid},
        {PD::LocalState(), PD::LocalState()}, shape, &batch).ok() &&
        batch.size() == 2 && batch[0].status.ok() &&
        batch[1].status.code == PD::StatusCode::InvalidParticle;

    // Direct evaluation is reentrant: every worker receives immutable copies
    // of one mathematical state and must reproduce the scalar bit-for-bit.
    // The startup-only active manager is intentionally not reconfigured while
    // these calls run, matching the documented lifetime contract.
    const double expectedThreadLambda =
        shaped.lambdaParallelM.value_or(-1.0);
    std::vector<double> threadedLambda(8, 0.0);
    std::vector<std::thread> workers;
    for (std::size_t worker = 0; worker < threadedLambda.size(); ++worker) {
      workers.emplace_back([&, worker]() {
        for (int repetition = 0; repetition < 32; ++repetition) {
          const PD::ParallelResult value =
              PD::Evaluate(particle, PD::LocalState(), shape);
          if (!value.status.ok() || !value.lambdaParallelM.has_value() ||
              *value.lambdaParallelM != expectedThreadLambda) {
            threadedLambda[worker] = -1.0;
            return;
          }
          threadedLambda[worker] = *value.lambdaParallelM;
        }
      });
    }
    for (std::thread& worker : workers) worker.join();
    for (double value : threadedLambda)
      pass = pass && value == expectedThreadLambda;

    // With no cache, a changed local turbulence snapshot must immediately
    // change the exact QLT result.  Batch and scalar evaluation must agree for
    // each point rather than accidentally reusing the first state.
    PD::ModelConfiguration stateSensitive;
    stateSensitive.model = PD::ModelId::QltSlabSpectrum;
    stateSensitive.qltSlab.spectrum.form = PD::SpectrumForm::SmoothBendover;
    const PD::LocalState weaker = State(1.0, 0.04, 1.0);
    const PD::LocalState stronger = State(1.0, 0.08, 1.0);
    const PD::ParallelResult weakerScalar =
        PD::Evaluate(particle, weaker, stateSensitive);
    const PD::ParallelResult strongerScalar =
        PD::Evaluate(particle, stronger, stateSensitive);
    pass = pass && weakerScalar.status.ok() && strongerScalar.status.ok() &&
        std::fabs(*weakerScalar.lambdaParallelM /
                  *strongerScalar.lambdaParallelM - 2.0) < 1.0e-13;
    pass = pass && PD::EvaluateBatch({particle, particle}, {weaker, stronger},
        stateSensitive, &batch).ok() && batch.size() == 2 &&
        *batch[0].lambdaParallelM == *weakerScalar.lambdaParallelM &&
        *batch[1].lambdaParallelM == *strongerScalar.lambdaParallelM;

    // Every parameterized first-release backend is reachable through its own
    // fail-closed semantic reader.  These are parser-neutral assignments; the
    // deferred D14 host syntax is deliberately not exercised here.
    PD::ModelConfiguration parsed;
    pass = pass && PD::BuildConfiguration("qlt_slab_spectrum",
        {{"spectrum_form", "smooth_bendover"}}, &parsed).ok();
    pass = pass && PD::BuildConfiguration("qlt_slab_spectrum",
        {{"spectrum_form", "multirange"}, {"energy_range_index", "0"},
         {"dissipation_index", "1.5"},
         {"dissipation_wavenumber_rad_per_m", "1e-8"}}, &parsed).ok();
    pass = pass && PD::BuildConfiguration("qlt_slab_inertial", {}, &parsed).ok();
    pass = pass && PD::BuildConfiguration("prescribed_lambda_mu_shape",
        {{"amplitude_mode", "target_lambda"}, {"q_mu", "1.6666666666666667"},
         {"h_mu", "0.01"}, {"target_lambda_m", "2.5e9"}}, &parsed).ok();
    pass = pass && PD::BuildConfiguration("broadened_slab",
        {{"spectrum_form", "smooth_bendover"},
         {"kernel", "gaussian"}, {"width0_per_s", "0.1"}}, &parsed).ok();
    for (const char* model : {"nlpa_given_perp", "nlgc_e", "nlgce_n",
                              "nlgce_f_2014"})
      pass = pass && PD::BuildConfiguration(model, {}, &parsed).ok();
    pass = pass && PD::BuildConfiguration("turbulence_adapter",
        {{"energy_convention", "elsasser_sum"},
         {"residual_energy_convention", "magnetic_minus_kinetic"},
         {"moment_source", "input_parameters"},
         {"underlying_closure", "qlt_slab_spectrum"},
         {"vacuum_permeability_H_per_m", "1.25663706212e-6"},
         {"provider_moment_m2_per_s2", "8"}, {"residual_energy", "0.2"},
         {"slab_fraction", "0.2"}, {"spectrum_form", "smooth_bendover"}},
        &parsed).ok();
    pass = pass && PD::BuildConfiguration("wave_spectrum_adapter",
        {{"propagation", "balanced"},
         {"polarization", "transverse_axisymmetric"}, {"frame", "plasma"},
         {"underlying_closure", "qlt_slab_spectrum"},
         {"spectrum_form", "smooth_bendover"}}, &parsed).ok();
    pass = pass && PD::BuildConfiguration("wave_spectrum_adapter",
        {{"propagation", "imbalanced"},
         {"polarization", "transverse_axisymmetric"}, {"frame", "plasma"},
         {"underlying_closure", "qlt_slab_spectrum"},
         {"spectrum_form", "smooth_bendover"}}, &parsed).code ==
        PD::StatusCode::InvalidConfiguration;
    pass = pass && PD::BuildConfiguration("tabulated_parallel",
        {{"stored_quantity", "lambda_parallel"}, {"axes", "rigidity,time"},
         {"axis_0_values_SI", "1,2"}, {"axis_1_values_SI", "0,1"},
         {"coefficient_values_SI", "1,2,2,4"},
         {"time_rule", "linear"},
         {"generation_identity", "fixture"}}, &parsed).ok();
    pass = pass && PD::BuildConfiguration("nlgce_f_2014",
        {{"unknown", "1"}}, &parsed).code ==
        PD::StatusCode::InvalidConfiguration;

    std::cout << "0 " << (pass ? 1 : 0) << " nan nan\n";
    return pass ? 0 : 1;
  }
  if (mode == "qlt" && argc == 4) {
    const double r = Number(argv[2]);
    const double eps2 = Number(argv[3]);
    PD::ModelConfiguration configuration;
    configuration.model = PD::ModelId::QltSlabSpectrum;
    configuration.qltSlab.spectrum.form = PD::SpectrumForm::SmoothBendover;
    configuration.qltSlab.numerical.relativeTolerance = 1.0e-11;
    return Print(PD::Evaluate(Particle(r), State(1.0, eps2, 1.0),
                              configuration));
  }
  if ((mode == "nlgce_f" || mode == "nlgce_n" || mode == "nlgc_e") &&
      argc == 6) {
    const double r = Number(argv[2]), fs = Number(argv[3]);
    const double eps2 = Number(argv[4]), rho = Number(argv[5]);
    PD::ModelConfiguration configuration;
    configuration.model = mode == "nlgce_f" ? PD::ModelId::NlgceF2014 :
        (mode == "nlgce_n" ? PD::ModelId::NlgceN : PD::ModelId::NlgcE);
    configuration.nonlinear.numerical.relativeTolerance = 1.0e-11;
    configuration.nonlinear.numerical.maximumRefinements = 22;
    return Print(PD::Evaluate(Particle(r), State(fs, eps2, rho),
                              configuration));
  }
  if (mode == "nlpa" && argc == 7) {
    const double r = Number(argv[2]), fs = Number(argv[3]);
    const double eps2 = Number(argv[4]), rho = Number(argv[5]);
    const double kxOverVEll = Number(argv[6]);
    PD::ModelConfiguration configuration;
    configuration.model = PD::ModelId::NlpaGivenPerp;
    configuration.nonlinear.numerical.relativeTolerance = 1.0e-11;
    configuration.nonlinear.numerical.maximumRefinements = 22;
    PD::LocalState state = State(fs, eps2, rho);
    PD::ParticleKinematics kin;
    PD::ParticleState particle = Particle(r);
    PD::ComputeParticleKinematics(particle, &kin);
    state.suppliedKappaPerpendicularM2PerS =
        kxOverVEll * kin.speedMPerS * SlabLengthM;
    state.suppliedPerpendicularModelId = "mathematical_fixture";
    return Print(PD::Evaluate(particle, state, configuration));
  }
  if ((mode == "broadened" || mode == "broadened_dmu") && argc == 6) {
    const double r = Number(argv[2]), eps2 = Number(argv[3]);
    const std::string kernel = argv[4];
    const double widthOverOmega = Number(argv[5]);
    PD::ParticleState particle = Particle(r);
    PD::ParticleKinematics kin;
    PD::ComputeParticleKinematics(particle, &kin);
    const double omega = std::fabs(particle.chargeC) * FieldT /
        (kin.gamma * particle.massKg);
    PD::ModelConfiguration configuration;
    configuration.model = PD::ModelId::BroadenedSlab;
    configuration.broadenedSlab.spectrum.form = PD::SpectrumForm::SmoothBendover;
    configuration.broadenedSlab.kernel = kernel == "gaussian"
        ? PD::BroadeningKernel::Gaussian
        : PD::BroadeningKernel::LorentzianConstant;
    configuration.broadenedSlab.width0PerS = widthOverOmega * omega;
    configuration.broadenedSlab.numerical.relativeTolerance = 1.0e-7;
    configuration.broadenedSlab.numerical.maximumRefinements = 22;
    PD::LocalState state = State(1.0, eps2, 1.0);
    if (mode == "broadened_dmu") {
      double d = 0.0;
      const PD::Status status = PD::EvaluatePitchAngleDiffusion(
          0.5, particle, state, configuration, &d);
      if (!status.ok()) std::cerr << status.detail << '\n';
      std::cout << static_cast<int>(status.code) << ' ' <<
          std::setprecision(17) << d << " nan nan\n";
      return status.ok() ? 0 : 1;
    }
    return Print(PD::Evaluate(particle, state, configuration));
  }
  std::cerr << "unknown advanced-test operation or argument count\n";
  return 2;
}
