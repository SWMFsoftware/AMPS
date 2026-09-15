#include "sep_coefficient_physics.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace SEP {
namespace Transport {
namespace CoefficientPhysics {
namespace {

const double Pi = 3.141592653589793238462643383279502884;

PitchAngleResult Error(StatusCode code, const std::string& message) {
  PitchAngleResult result;
  result.status = Status::Error(code, message);
  result.valueState = ValueState::Invalid;
  return result;
}

bool NearlyEqual(double left, double right) {
  const double scale = std::max(1.0, std::max(std::fabs(left), std::fabs(right)));
  return std::fabs(left - right) <=
      32.0 * std::numeric_limits<double>::epsilon() * scale;
}

double ScaledWavenumber(double reference, double referenceRadius,
                        double radius, double exponent) {
  return reference * std::pow(referenceRadius / radius, exponent);
}

double FlorinskiyPower(double kAbsPerM, double varianceT2,
                       const FlorinskiyParameters& p) {
  if (!(kAbsPerM >= 0.0) || !(varianceT2 >= 0.0) ||
      varianceT2 == 0.0) return 0.0;
  const double halfIndex = 0.5 * p.spectralIndex;
  const double k0 = std::sqrt(Pi) * std::tgamma(halfIndex) /
      (std::tgamma(halfIndex - 0.5) * p.parallelCorrelationLengthM);
  const double normalization = 2.0 * std::tgamma(halfIndex) /
      (std::sqrt(Pi) * std::tgamma(halfIndex - 0.5) * k0);
  const double ratio = kAbsPerM / k0;
  return normalization * varianceT2 *
      std::pow(1.0 + ratio * ratio, -halfIndex);
}

double FlorinskiyValue(const LocalInputView& input,
                       const FlorinskiyParameters& parameters,
                       const SpeciesProperties& species,
                       double speedMPerS, double mu) {
  const ScalarResult omega = GyrofrequencyRadPerS(
      input.magneticFieldT, species);
  if (!omega.status.ok()) return std::numeric_limits<double>::quiet_NaN();
  const double parallelSpeed = speedMPerS * mu;
  const double plusDenominator =
      std::fabs(parallelSpeed - input.alfvenSpeedMPerS);
  const double minusDenominator =
      std::fabs(parallelSpeed + input.alfvenSpeedMPerS);

  // For spectralIndex>1, P(Omega/t)/t tends to zero as t tends to zero.
  // Encoding that analytic limit avoids Omega/0 and preserves the physical
  // no-resonance limit at exact wave-frame co-motion.
  const double plus = plusDenominator == 0.0 ? 0.0 :
      FlorinskiyPower(omega.value / plusDenominator,
                      input.deltaBPlus2T2, parameters) / plusDenominator;
  const double minus = minusDenominator == 0.0 ? 0.0 :
      FlorinskiyPower(omega.value / minusDenominator,
                      input.deltaBMinus2T2, parameters) / minusDenominator;
  return Pi * omega.value * omega.value * (1.0 - mu * mu) *
      (plus + minus) /
      (4.0 * input.magneticFieldT * input.magneticFieldT);
}

struct IntegrandState {
  bool failed = false;
  bool gap = false;
  Status failure;
  std::size_t evaluations = 0;
};

double Integrand(double mu, const PitchAngleFunction& dmumu,
                 IntegrandState* state) {
  ++state->evaluations;
  const PitchAngleResult sample = dmumu(mu);
  if (!sample.status.ok() || sample.valueState != ValueState::Finite ||
      !std::isfinite(sample.dMuMuPerS) || sample.dMuMuPerS < 0.0) {
    state->failed = true;
    state->failure = sample.status.ok()
        ? Status::Error(StatusCode::InvalidCoefficient,
                        "Dmumu quadrature callback returned an invalid state")
        : sample.status;
    return 0.0;
  }
  const double shape = std::max(0.0, 1.0 - mu * mu);
  if (sample.dMuMuPerS == 0.0) {
    // Endpoint zeros have zero measure and are sampled just inside the interval
    // by the caller.  Any remaining zero therefore represents a finite-width
    // or interior resonance gap and cannot be regularized numerically.
    if (shape > 64.0 * std::numeric_limits<double>::epsilon()) state->gap = true;
    return 0.0;
  }
  return shape * shape / sample.dMuMuPerS;
}

double AdaptiveSimpson(const PitchAngleFunction& function,
                       double left, double right,
                       double fLeft, double fMid, double fRight,
                       double whole, double tolerance, int depth,
                       IntegrandState* state, double* error) {
  const double mid = 0.5 * (left + right);
  const double leftMid = 0.5 * (left + mid);
  const double rightMid = 0.5 * (mid + right);
  const double fLeftMid = Integrand(leftMid, function, state);
  const double fRightMid = Integrand(rightMid, function, state);
  if (state->failed || state->gap) return 0.0;
  const double leftIntegral = (mid - left) *
      (fLeft + 4.0 * fLeftMid + fMid) / 6.0;
  const double rightIntegral = (right - mid) *
      (fMid + 4.0 * fRightMid + fRight) / 6.0;
  const double refined = leftIntegral + rightIntegral;
  const double localError = std::fabs(refined - whole) / 15.0;
  if (depth <= 0 || localError <= tolerance) {
    *error += localError;
    return refined + (refined - whole) / 15.0;
  }
  return AdaptiveSimpson(function, left, mid, fLeft, fLeftMid, fMid,
                         leftIntegral, 0.5 * tolerance, depth - 1,
                         state, error) +
         AdaptiveSimpson(function, mid, right, fMid, fRightMid, fRight,
                         rightIntegral, 0.5 * tolerance, depth - 1,
                         state, error);
}

}  // namespace

Status ValidateSpecies(const SpeciesProperties& species) {
  if (species.modelSpecies < 0 || species.name.empty() ||
      !std::isfinite(species.signedChargeC) || species.signedChargeC == 0.0 ||
      !std::isfinite(species.restMassKg) || species.restMassKg <= 0.0 ||
      !std::isfinite(species.nucleonCount) || species.nucleonCount < 0.0) {
    return Status::Error(StatusCode::InvalidArgument,
        "species requires an ID, name, nonzero signed charge, positive rest "
        "mass, and non-negative nucleon count");
  }
  return Status::Ok();
}

ScalarResult GyrofrequencyRadPerS(double magneticFieldT,
                                  const SpeciesProperties& species) {
  ScalarResult result;
  const Status valid = ValidateSpecies(species);
  if (!valid.ok() || !std::isfinite(magneticFieldT) || magneticFieldT <= 0.0) {
    result.status = valid.ok()
        ? Status::Error(StatusCode::InvalidArgument,
                        "gyrofrequency requires B>0 tesla")
        : valid;
    return result;
  }
  result.value = std::fabs(species.signedChargeC) * magneticFieldT /
      species.restMassKg;
  result.status = std::isfinite(result.value) ? Status::Ok() :
      Status::Error(StatusCode::InvalidCoefficient,
                    "gyrofrequency overflowed");
  return result;
}

ScalarResult LarmorRadiusM(double perpendicularMomentumKgMPerS,
                           double magneticFieldT,
                           const SpeciesProperties& species) {
  ScalarResult result;
  const Status valid = ValidateSpecies(species);
  if (!valid.ok() || !std::isfinite(perpendicularMomentumKgMPerS) ||
      perpendicularMomentumKgMPerS < 0.0 ||
      !std::isfinite(magneticFieldT) || magneticFieldT <= 0.0) {
    result.status = valid.ok()
        ? Status::Error(StatusCode::InvalidArgument,
                        "Larmor radius requires p_perp>=0 and B>0")
        : valid;
    return result;
  }
  result.value = perpendicularMomentumKgMPerS /
      (std::fabs(species.signedChargeC) * magneticFieldT);
  result.status = std::isfinite(result.value) ? Status::Ok() :
      Status::Error(StatusCode::InvalidCoefficient,
                    "Larmor radius overflowed");
  return result;
}

ScalarResult RigidityVolt(double momentumKgMPerS,
                          const SpeciesProperties& species,
                          double speedOfLightMPerS) {
  ScalarResult result;
  const Status valid = ValidateSpecies(species);
  if (!valid.ok() || !std::isfinite(momentumKgMPerS) ||
      momentumKgMPerS < 0.0 || !std::isfinite(speedOfLightMPerS) ||
      speedOfLightMPerS <= 0.0) {
    result.status = valid.ok()
        ? Status::Error(StatusCode::InvalidArgument,
                        "rigidity requires p>=0 and c>0")
        : valid;
    return result;
  }
  result.value = momentumKgMPerS * speedOfLightMPerS /
      std::fabs(species.signedChargeC);
  result.status = std::isfinite(result.value) ? Status::Ok() :
      Status::Error(StatusCode::InvalidCoefficient, "rigidity overflowed");
  return result;
}

Status ValidateSpectrum(const SpectrumParameters& s) {
  if (!std::isfinite(s.referenceRadiusM) || s.referenceRadiusM <= 0.0 ||
      !std::isfinite(s.kMinAtReferencePerM) ||
      !std::isfinite(s.kMaxAtReferencePerM) ||
      s.kMinAtReferencePerM <= 0.0 ||
      !(s.kMaxAtReferencePerM > s.kMinAtReferencePerM) ||
      !std::isfinite(s.kMinRadialExponent) ||
      !std::isfinite(s.kMaxRadialExponent) ||
      !std::isfinite(s.spectralIndex) || s.spectralIndex <= 1.0) {
    return Status::Error(StatusCode::InvalidArgument,
        "spectrum requires positive ordered bounds, positive reference "
        "radius, finite radial exponents, and spectral index greater than one");
  }
  return Status::Ok();
}

Status ValidateLocalInput(const LocalInputView& input) {
  if (input.source.empty() || input.representation.empty() ||
      input.generation == 0 || input.checksum == 0 ||
      !std::isfinite(input.heliocentricRadiusM) ||
      input.heliocentricRadiusM <= 0.0 ||
      !std::isfinite(input.magneticFieldT) || input.magneticFieldT <= 0.0 ||
      !std::isfinite(input.deltaB2T2) || input.deltaB2T2 < 0.0 ||
      !std::isfinite(input.deltaBPlus2T2) || input.deltaBPlus2T2 < 0.0 ||
      !std::isfinite(input.deltaBMinus2T2) || input.deltaBMinus2T2 < 0.0 ||
      !std::isfinite(input.alfvenSpeedMPerS) ||
      input.alfvenSpeedMPerS < 0.0) {
    return Status::Error(StatusCode::InvalidCoefficient,
        "source-bound coefficient input is incomplete or non-physical");
  }
  return Status::Ok();
}

PitchAngleResult EvaluateConstantDmumu(double configuredDmumuPerS,
                                       double mu) {
  PitchAngleResult result;
  if (!std::isfinite(configuredDmumuPerS) || configuredDmumuPerS < 0.0 ||
      !std::isfinite(mu) || std::fabs(mu) > 1.0) {
    return Error(StatusCode::InvalidCoefficient,
                 "constant Dmumu requires a finite non-negative rate and |mu|<=1");
  }
  result.dMuMuPerS = configuredDmumuPerS;
  result.dDmuMuDmuPerS = 0.0;
  result.valueState = ValueState::Finite;
  result.status = Status::Ok();
  return result;
}

PitchAngleResult EvaluateJokipiiSlab(
    const LocalInputView& input, const SpectrumParameters& spectrum,
    const SpeciesProperties& species, double speedMPerS, double mu) {
  const Status inputStatus = ValidateLocalInput(input);
  const Status spectrumStatus = ValidateSpectrum(spectrum);
  const Status speciesStatus = ValidateSpecies(species);
  if (!inputStatus.ok()) return Error(inputStatus.code, inputStatus.message);
  if (!spectrumStatus.ok())
    return Error(spectrumStatus.code, spectrumStatus.message);
  if (!speciesStatus.ok()) return Error(speciesStatus.code, speciesStatus.message);
  if (!std::isfinite(speedMPerS) || speedMPerS <= 0.0 ||
      !std::isfinite(mu) || std::fabs(mu) > 1.0) {
    return Error(StatusCode::InvalidArgument,
                 "Jokipii kernel requires v>0 and |mu|<=1");
  }

  PitchAngleResult result;
  result.valueState = ValueState::Finite;
  const ScalarResult omega = GyrofrequencyRadPerS(input.magneticFieldT, species);
  if (!omega.status.ok()) return Error(omega.status.code, omega.status.message);
  const double kMin = ScaledWavenumber(
      spectrum.kMinAtReferencePerM, spectrum.referenceRadiusM,
      input.heliocentricRadiusM, spectrum.kMinRadialExponent);
  const double kMax = ScaledWavenumber(
      spectrum.kMaxAtReferencePerM, spectrum.referenceRadiusM,
      input.heliocentricRadiusM, spectrum.kMaxRadialExponent);
  if (!std::isfinite(kMin) || !std::isfinite(kMax) || !(kMax > kMin)) {
    return Error(StatusCode::InvalidCoefficient,
                 "local resonance spectrum has invalid bounds");
  }

  const double absMu = std::fabs(mu);
  if (absMu == 0.0 || input.deltaB2T2 == 0.0) {
    result.status = Status::Ok();
    return result;
  }
  const double resonantSpeed = speedMPerS * absMu;
  const double kRes = omega.value / resonantSpeed;
  if (kRes < kMin || kRes > kMax) {
    result.status = Status::Ok();
    return result;
  }
  if (NearlyEqual(kRes, kMin) || NearlyEqual(kRes, kMax)) {
    result.derivativeState = DerivativeState::Nondifferentiable;
    result.status = Status::Error(StatusCode::NondifferentiableCoefficient,
        "Jokipii resonance lies exactly on a finite spectral-band boundary");
    return result;
  }

  const double index = spectrum.spectralIndex;
  const double normalizationDenominator =
      std::pow(kMin, 1.0 - index) - std::pow(kMax, 1.0 - index);
  if (!(normalizationDenominator > 0.0) ||
      !std::isfinite(normalizationDenominator)) {
    return Error(StatusCode::InvalidCoefficient,
                 "power-law spectrum normalization is invalid");
  }
  const double normalization = (index - 1.0) * input.deltaB2T2 /
      normalizationDenominator;
  const double powerT2M = normalization * std::pow(kRes, -index);
  const double amplitude = (Pi / 2.0) * omega.value * omega.value *
      powerT2M /
      (input.magneticFieldT * input.magneticFieldT * resonantSpeed);
  result.dMuMuPerS = amplitude * (1.0 - mu * mu);

  // Because k_res is proportional to 1/|mu| and P is proportional to k^-q,
  // the smooth-branch logarithmic derivative is
  // -2mu/(1-mu^2)+(q-1)/mu.  At |mu|=1 the first form is singular although D
  // has a finite one-sided derivative; use that limit explicitly.
  if (absMu == 1.0) {
    result.dDmuMuDmuPerS = mu > 0.0 ? -2.0 * amplitude : 2.0 * amplitude;
    result.derivativeState = DerivativeState::OneSidedLimit;
  }
  else {
    result.dDmuMuDmuPerS = result.dMuMuPerS *
        (-2.0 * mu / (1.0 - mu * mu) + (index - 1.0) / mu);
  }
  if (!std::isfinite(result.dMuMuPerS) || result.dMuMuPerS < 0.0 ||
      !std::isfinite(result.dDmuMuDmuPerS)) {
    return Error(StatusCode::InvalidCoefficient,
                 "Jokipii coefficient or derivative is not finite");
  }
  result.status = Status::Ok();
  return result;
}

PitchAngleResult EvaluateFlorinskiySlab(
    const LocalInputView& input, const FlorinskiyParameters& parameters,
    const SpeciesProperties& species, double speedMPerS, double mu) {
  const Status inputStatus = ValidateLocalInput(input);
  const Status speciesStatus = ValidateSpecies(species);
  if (!inputStatus.ok()) return Error(inputStatus.code, inputStatus.message);
  if (!speciesStatus.ok()) return Error(speciesStatus.code, speciesStatus.message);
  if (!std::isfinite(parameters.spectralIndex) ||
      parameters.spectralIndex <= 1.0 ||
      !std::isfinite(parameters.parallelCorrelationLengthM) ||
      parameters.parallelCorrelationLengthM <= 0.0 ||
      !std::isfinite(speedMPerS) || speedMPerS <= 0.0 ||
      !std::isfinite(mu) || std::fabs(mu) > 1.0) {
    return Error(StatusCode::InvalidArgument,
                 "Florinskiy kernel received invalid parameters, speed, or mu");
  }

  const double plusFrame = speedMPerS * mu - input.alfvenSpeedMPerS;
  const double minusFrame = speedMPerS * mu + input.alfvenSpeedMPerS;
  if (plusFrame == 0.0 || minusFrame == 0.0) {
    PitchAngleResult result = Error(StatusCode::NondifferentiableCoefficient,
        "Florinskiy dynamic resonance is nondifferentiable at wave-frame co-motion");
    result.derivativeState = DerivativeState::Nondifferentiable;
    return result;
  }

  PitchAngleResult result;
  result.valueState = ValueState::Finite;
  result.dMuMuPerS = FlorinskiyValue(input, parameters, species, speedMPerS, mu);
  const double hNominal = 1.0e-5;
  if (mu <= -1.0 + hNominal) {
    const double h = std::min(hNominal, 0.25);
    const double f0 = FlorinskiyValue(input, parameters, species, speedMPerS, mu);
    const double f1 = FlorinskiyValue(input, parameters, species, speedMPerS, mu + h);
    const double f2 = FlorinskiyValue(input, parameters, species, speedMPerS, mu + 2.0 * h);
    result.dDmuMuDmuPerS = (-3.0 * f0 + 4.0 * f1 - f2) / (2.0 * h);
    result.derivativeState = DerivativeState::OneSidedLimit;
  }
  else if (mu >= 1.0 - hNominal) {
    const double h = std::min(hNominal, 0.25);
    const double f0 = FlorinskiyValue(input, parameters, species, speedMPerS, mu);
    const double f1 = FlorinskiyValue(input, parameters, species, speedMPerS, mu - h);
    const double f2 = FlorinskiyValue(input, parameters, species, speedMPerS, mu - 2.0 * h);
    result.dDmuMuDmuPerS = (3.0 * f0 - 4.0 * f1 + f2) / (2.0 * h);
    result.derivativeState = DerivativeState::OneSidedLimit;
  }
  else {
    const double distance = std::min(1.0 - mu, mu + 1.0);
    const double h = std::min(hNominal, 0.25 * distance);
    const double plus = FlorinskiyValue(
        input, parameters, species, speedMPerS, mu + h);
    const double minus = FlorinskiyValue(
        input, parameters, species, speedMPerS, mu - h);
    result.dDmuMuDmuPerS = (plus - minus) / (2.0 * h);
  }
  if (!std::isfinite(result.dMuMuPerS) || result.dMuMuPerS < 0.0 ||
      !std::isfinite(result.dDmuMuDmuPerS)) {
    return Error(StatusCode::InvalidCoefficient,
                 "Florinskiy coefficient or bounded derivative is invalid");
  }
  result.status = Status::Ok();
  return result;
}

MeanFreePathResult EvaluateCorrelationMeanFreePath(
    const LocalInputView& input, const SpeciesProperties& species,
    double momentumKgMPerS, double correlationLengthAtReferenceM,
    double referenceRadiusM) {
  MeanFreePathResult result;
  const Status inputStatus = ValidateLocalInput(input);
  const Status speciesStatus = ValidateSpecies(species);
  if (!inputStatus.ok()) { result.status = inputStatus; return result; }
  if (!speciesStatus.ok()) { result.status = speciesStatus; return result; }
  if (!std::isfinite(momentumKgMPerS) || momentumKgMPerS < 0.0 ||
      !std::isfinite(correlationLengthAtReferenceM) ||
      correlationLengthAtReferenceM <= 0.0 ||
      !std::isfinite(referenceRadiusM) || referenceRadiusM <= 0.0) {
    result.status = Status::Error(StatusCode::InvalidArgument,
        "correlation mean free path requires p>=0 and positive scales");
    return result;
  }
  if (input.deltaB2T2 == 0.0) {
    result.valueState = ValueState::Ballistic;
    result.lambdaParallelM = std::numeric_limits<double>::infinity();
    result.status = Status::Ok();
    return result;
  }
  const double correlationLengthM = correlationLengthAtReferenceM *
      input.heliocentricRadiusM / referenceRadiusM;
  const ScalarResult larmor = LarmorRadiusM(
      std::sqrt(2.0 / 3.0) * momentumKgMPerS,
      input.magneticFieldT, species);
  if (!larmor.status.ok()) { result.status = larmor.status; return result; }
  const double ratio2 = input.deltaB2T2 /
      (input.magneticFieldT * input.magneticFieldT);
  result.lambdaParallelM = correlationLengthM / ratio2 *
      std::pow(larmor.value / correlationLengthM, 1.0 / 3.0);
  result.valueState = ValueState::Finite;
  result.status = std::isfinite(result.lambdaParallelM) &&
      result.lambdaParallelM > 0.0 ? Status::Ok() :
      Status::Error(StatusCode::InvalidCoefficient,
                    "species-aware correlation mean free path is invalid");
  return result;
}

SpatialDiffusionResult IntegrateSpatialDiffusion(
    double speedMPerS, const PitchAngleFunction& dmumu,
    const SpatialQuadratureConfiguration& configuration) {
  SpatialDiffusionResult result;
  if (!std::isfinite(speedMPerS) || speedMPerS <= 0.0 || !dmumu ||
      !std::isfinite(configuration.absoluteToleranceM2PerS) ||
      configuration.absoluteToleranceM2PerS < 0.0 ||
      !std::isfinite(configuration.relativeTolerance) ||
      configuration.relativeTolerance <= 0.0 ||
      configuration.maximumRecursion < 1 ||
      configuration.maximumRecursion > 60) {
    result.status = Status::Error(StatusCode::InvalidArgument,
                                  "invalid spatial quadrature configuration");
    return result;
  }

  IntegrandState state;
  double integral = 0.0;
  double integralError = 0.0;
  // Move only the two measure-zero physical endpoints inward.  mu=0 remains an
  // exact sample so a ninety-degree resonance gap is detected, not skipped.
  const double endpoint = 1.0 - 1.0e-12;
  const double intervals[3] = {-endpoint, 0.0, endpoint};
  const double prefactor = speedMPerS * speedMPerS / 8.0;
  const double integralAbsTol = configuration.absoluteToleranceM2PerS /
      prefactor;
  for (int part = 0; part < 2; ++part) {
    const double left = intervals[part];
    const double right = intervals[part + 1];
    const double mid = 0.5 * (left + right);
    const double fLeft = Integrand(left, dmumu, &state);
    const double fMid = Integrand(mid, dmumu, &state);
    const double fRight = Integrand(right, dmumu, &state);
    if (state.failed || state.gap) break;
    const double whole = (right - left) *
        (fLeft + 4.0 * fMid + fRight) / 6.0;
    const double tolerance = 0.5 * std::max(
        integralAbsTol,
        configuration.relativeTolerance * std::fabs(whole));
    integral += AdaptiveSimpson(dmumu, left, right, fLeft, fMid, fRight,
        whole, tolerance, configuration.maximumRecursion, &state,
        &integralError);
    if (state.failed || state.gap) break;
  }
  result.evaluations = state.evaluations;
  if (state.failed) {
    result.status = state.failure;
    return result;
  }
  if (state.gap) {
    if (configuration.gapPolicy == ResonanceGapPolicy::Ballistic) {
      result.valueState = ValueState::Ballistic;
      result.kappaParallelM2PerS = std::numeric_limits<double>::infinity();
      result.status = Status::Ok();
    }
    else {
      result.valueState = ValueState::Unavailable;
      result.status = Status::Error(StatusCode::UnresolvedCoefficient,
          "Dmumu has an interior resonance gap; no numerical floor is configured");
    }
    return result;
  }
  result.kappaParallelM2PerS = prefactor * integral;
  result.estimatedAbsoluteErrorM2PerS = prefactor * integralError;
  result.valueState = ValueState::Finite;
  result.status = std::isfinite(result.kappaParallelM2PerS) &&
      result.kappaParallelM2PerS >= 0.0 ? Status::Ok() :
      Status::Error(StatusCode::InvalidCoefficient,
                    "adaptive spatial-diffusion quadrature failed");
  return result;
}

}  // namespace CoefficientPhysics
}  // namespace Transport
}  // namespace SEP
