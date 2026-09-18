#include "sep_physics_extensions.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>

namespace SEP {
namespace PhysicsExtensions {
namespace {

Transport::Status Error(Transport::StatusCode code,
                        const std::string& message) {
  return Transport::Status::Error(code, message);
}

bool Finite(double value) { return std::isfinite(value); }

double Sum(const std::vector<double>& values) {
  long double result = 0.0L;
  for (std::size_t i = 0; i < values.size(); ++i) result += values[i];
  return static_cast<double>(result);
}

bool StrictPositiveGrid(const std::vector<double>& values) {
  if (values.empty()) return false;
  for (std::size_t i = 0; i < values.size(); ++i)
    if (!Finite(values[i]) || values[i] <= 0.0 ||
        (i && values[i] <= values[i - 1])) return false;
  return true;
}

}  // namespace

Transport::ScalarResult ParkerItoDriftMPerS(
    const ParkerGeometryInput& input) {
  Transport::ScalarResult result;
  if (!Finite(input.plasmaAdvectionMPerS) ||
      !Finite(input.kappaParallelM2PerS) ||
      !Finite(input.dKappaDsMPerS) ||
      !Finite(input.dLnAreaDsPerM) || input.kappaParallelM2PerS < 0.0) {
    result.status = Error(Transport::StatusCode::InvalidArgument,
                          "Parker geometry coefficients are invalid");
    return result;
  }
  // For volume density n, the one-dimensional conservative flux contains
  // A*kappa*dn/ds.  Rewriting the Fokker-Planck equation for the conserved
  // walker density q=A*n adds kappa*d(ln A)/ds to the Ito drift.  A walker
  // already defined per arc length has no additional area drift.
  result.value = input.plasmaAdvectionMPerS + input.dKappaDsMPerS;
  if (input.measure == ParkerMeasure::PerVolume)
    result.value += input.kappaParallelM2PerS * input.dLnAreaDsPerM;
  result.status = Transport::Status::Ok();
  return result;
}

Transport::ScalarResult WalkerToPhysicalDensity(
    ParkerMeasure measure, double walkerDensity, double tubeAreaM2) {
  Transport::ScalarResult result;
  if (!Finite(walkerDensity) || walkerDensity < 0.0 ||
      !Finite(tubeAreaM2) || tubeAreaM2 <= 0.0) {
    result.status = Error(Transport::StatusCode::InvalidArgument,
                          "walker density conversion requires positive area");
    return result;
  }
  result.value = measure == ParkerMeasure::PerArcLength
      ? walkerDensity / tubeAreaM2 : walkerDensity;
  result.status = Transport::Status::Ok();
  return result;
}

ResonanceMetadata SolveDynamicResonance(
    const DynamicResonanceInput& input) {
  ResonanceMetadata result;
  result.branch = input.branch;
  result.harmonic = input.harmonic;
  if (!Finite(input.particleSpeedMPerS) || input.particleSpeedMPerS < 0.0 ||
      !Finite(input.mu) || std::fabs(input.mu) > 1.0 ||
      !Finite(input.chargeC) || input.chargeC == 0.0 ||
      !Finite(input.massKg) || input.massKg <= 0.0 ||
      !Finite(input.signedMagneticFieldT) ||
      input.signedMagneticFieldT == 0.0 ||
      !Finite(input.alfvenSpeedMPerS) || input.alfvenSpeedMPerS < 0.0 ||
      !Finite(input.speedOfLightMPerS) || input.speedOfLightMPerS <= 0.0 ||
      input.particleSpeedMPerS >= input.speedOfLightMPerS ||
      input.harmonic == 0 || (input.branch != -1 && input.branch != 1) ||
      !StrictPositiveGrid(input.binCentersPerM)) {
    result.status = Error(Transport::StatusCode::InvalidArgument,
                          "dynamic resonance input is invalid");
    result.reason = result.status.message;
    return result;
  }
  const double beta = input.particleSpeedMPerS / input.speedOfLightMPerS;
  const double gamma = 1.0 / std::sqrt(1.0 - beta * beta);
  const double gyrofrequency = input.chargeC * input.signedMagneticFieldT /
                               (gamma * input.massKg);
  const double parallelSpeed = input.particleSpeedMPerS * input.mu;
  const double denominator = input.branch * input.alfvenSpeedMPerS -
                             parallelSpeed;
  const double scale = std::max(input.particleSpeedMPerS,
                                input.alfvenSpeedMPerS);
  if (std::fabs(denominator) <=
      32.0 * std::numeric_limits<double>::epsilon() *
      std::max(1.0, scale)) {
    result.resonanceStatus = ResonanceStatus::SingularComoving;
    result.status = Error(Transport::StatusCode::UnresolvedCoefficient,
                          "particle is comoving with the selected wave branch");
    result.reason = result.status.message;
    return result;
  }
  result.signedWaveNumberPerM =
      input.harmonic * gyrofrequency / denominator;
  result.polarization = result.signedWaveNumberPerM >= 0.0 ? 1 : -1;
  const double k = std::fabs(result.signedWaveNumberPerM);
  if (!Finite(k) || k == 0.0) {
    result.resonanceStatus = ResonanceStatus::NoRoot;
    result.status = Error(Transport::StatusCode::UnresolvedCoefficient,
                          "resonance relation has no finite nonzero root");
    result.reason = result.status.message;
    return result;
  }
  if (k < input.binCentersPerM.front() || k > input.binCentersPerM.back()) {
    result.resonanceStatus = ResonanceStatus::OutsideBand;
    result.status = Error(Transport::StatusCode::OutOfDomain,
                          "resonant wavenumber lies outside the spectral band");
    result.reason = result.status.message;
    return result;
  }
  std::vector<double>::const_iterator upper = std::lower_bound(
      input.binCentersPerM.begin(), input.binCentersPerM.end(), k);
  if (upper == input.binCentersPerM.begin()) {
    result.lowerBin = result.upperBin = 0;
    result.lowerWeight = 1.0;
  } else if (upper == input.binCentersPerM.end()) {
    result.lowerBin = result.upperBin = input.binCentersPerM.size() - 1;
    result.lowerWeight = 1.0;
  } else {
    result.upperBin = static_cast<std::size_t>(
        upper - input.binCentersPerM.begin());
    result.lowerBin = result.upperBin - 1;
    const double logLower = std::log(input.binCentersPerM[result.lowerBin]);
    const double logUpper = std::log(input.binCentersPerM[result.upperBin]);
    result.upperWeight = (std::log(k) - logLower) / (logUpper - logLower);
    result.lowerWeight = 1.0 - result.upperWeight;
  }
  result.resonanceStatus = ResonanceStatus::Resolved;
  result.status = Transport::Status::Ok();
  result.reason = "signed low-frequency cyclotron resonance resolved";
  return result;
}

PitchAngleCoefficient ApplyNinetyDegreeClosure(
    double mu, double slabD, double slabDerivative,
    const NinetyDegreeClosure& closure) {
  PitchAngleCoefficient result;
  if (!Finite(mu) || std::fabs(mu) > 1.0 || !Finite(slabD) || slabD < 0.0 ||
      !Finite(slabDerivative) || !Finite(closure.amplitudePerS) ||
      closure.amplitudePerS < 0.0 || !Finite(closure.halfWidthMu) ||
      closure.halfWidthMu <= 0.0 || closure.halfWidthMu > 1.0 ||
      (closure.amplitudePerS > 0.0 && closure.provenance.empty())) {
    result.status = Error(Transport::StatusCode::InvalidCoefficient,
                          "ninety-degree scattering closure is invalid");
    return result;
  }
  if (closure.amplitudePerS == 0.0) {
    result.dMuMuPerS = slabD;
    result.derivativePerS = slabDerivative;
    result.status = Transport::Status::Ok();
    return result;
  }
  const double u = mu / closure.halfWidthMu;
  const double gaussian = std::exp(-u * u);
  const double envelope = 1.0 - mu * mu;
  const double addition = closure.amplitudePerS * envelope * gaussian;
  const double derivative = closure.amplitudePerS * gaussian *
      (-2.0 * mu - 2.0 * mu * envelope /
                         (closure.halfWidthMu * closure.halfWidthMu));
  result.dMuMuPerS = slabD + addition;
  result.derivativePerS = slabDerivative + derivative;
  result.insideClosureRegion = std::fabs(mu) <= closure.halfWidthMu;
  result.status = Transport::Status::Ok();
  return result;
}

Transport::ScalarResult WaveEnergyJ(const WaveActionCell& cell,
                                    WaveInvariant invariant) {
  Transport::ScalarResult result;
  if (!Finite(cell.authoritativeValue) || cell.authoritativeValue < 0.0 ||
      !Finite(cell.volumeM3) || cell.volumeM3 <= 0.0 ||
      !Finite(cell.intrinsicFrequencyRadPerS) ||
      cell.intrinsicFrequencyRadPerS <= 0.0) {
    result.status = Error(Transport::StatusCode::InvalidArgument,
                          "wave invariant cell is invalid");
    return result;
  }
  result.value = invariant == WaveInvariant::IntegratedEnergy
      ? cell.authoritativeValue
      : cell.authoritativeValue * cell.intrinsicFrequencyRadPerS;
  result.status = Transport::Status::Ok();
  return result;
}

WaveActionUpdate ApplyGeometricConservation(
    const WaveActionCell& oldCell, double newVolumeM3,
    double newFrequency, WaveInvariant invariant) {
  WaveActionUpdate result;
  const Transport::ScalarResult oldEnergy = WaveEnergyJ(oldCell, invariant);
  if (!oldEnergy.status.ok() || !Finite(newVolumeM3) || newVolumeM3 <= 0.0 ||
      !Finite(newFrequency) || newFrequency <= 0.0) {
    result.status = oldEnergy.status.ok()
        ? Error(Transport::StatusCode::InvalidArgument,
                "new wave geometry/frequency is invalid")
        : oldEnergy.status;
    return result;
  }
  result.cell = oldCell;
  result.cell.volumeM3 = newVolumeM3;
  result.cell.intrinsicFrequencyRadPerS = newFrequency;
  if (invariant == WaveInvariant::IntegratedEnergy) {
    // Preserving uniform energy density under a pure mesh-volume change is the
    // discrete geometric-conservation-law target for this representation.
    result.cell.authoritativeValue = oldCell.authoritativeValue *
                                     newVolumeM3 / oldCell.volumeM3;
  } else {
    // Wave action is authoritative and remains unchanged for a reversible
    // adiabatic geometric update; the induced energy change is background work.
    result.cell.authoritativeValue = oldCell.authoritativeValue;
  }
  const Transport::ScalarResult newEnergy = WaveEnergyJ(result.cell, invariant);
  result.backgroundWorkJ = newEnergy.value - oldEnergy.value;
  result.status = newEnergy.status;
  return result;
}

SpectralCascadeResult AdvanceConservativeCascade(
    const std::vector<double>& plus, const std::vector<double>& minus,
    const std::vector<double>& plusFlux,
    const std::vector<double>& minusFlux,
    double dtS, double electronFraction) {
  SpectralCascadeResult result;
  const std::size_t n = plus.size();
  if (n == 0 || minus.size() != n || plusFlux.size() != n + 1 ||
      minusFlux.size() != n + 1 || !Finite(dtS) || dtS < 0.0 ||
      !Finite(electronFraction) || electronFraction < 0.0 ||
      electronFraction > 1.0) {
    result.status = Error(Transport::StatusCode::InvalidArgument,
                          "spectral cascade grid or heat partition is invalid");
    return result;
  }
  for (std::size_t i = 0; i < n; ++i)
    if (!Finite(plus[i]) || plus[i] < 0.0 || !Finite(minus[i]) || minus[i] < 0.0) {
      result.status = Error(Transport::StatusCode::InvalidParticleState,
                            "spectral cascade energy is non-finite or negative");
      return result;
    }
  for (std::size_t i = 0; i <= n; ++i)
    if (!Finite(plusFlux[i]) || !Finite(minusFlux[i]) ||
        plusFlux[i] < 0.0 || minusFlux[i] < 0.0) {
      result.status = Error(Transport::StatusCode::InvalidCoefficient,
                            "cascade interface flux must be finite and forward");
      return result;
    }
  result.plusEnergyJ = plus;
  result.minusEnergyJ = minus;
  for (std::size_t i = 0; i < n; ++i) {
    result.plusEnergyJ[i] += dtS * (plusFlux[i] - plusFlux[i + 1]);
    result.minusEnergyJ[i] += dtS * (minusFlux[i] - minusFlux[i + 1]);
    if (result.plusEnergyJ[i] < -1.0e-13 || result.minusEnergyJ[i] < -1.0e-13) {
      result.status = Error(Transport::StatusCode::StepUnderflow,
                            "cascade step violates positivity; subcycle required");
      return result;
    }
    result.plusEnergyJ[i] = std::max(0.0, result.plusEnergyJ[i]);
    result.minusEnergyJ[i] = std::max(0.0, result.minusEnergyJ[i]);
  }
  const double dissipated = dtS * (plusFlux[n] + minusFlux[n]);
  result.electronHeatJ = electronFraction * dissipated;
  result.ionHeatJ = (1.0 - electronFraction) * dissipated;
  const double initial = Sum(plus) + Sum(minus);
  const double final = Sum(result.plusEnergyJ) + Sum(result.minusEnergyJ);
  const double injection = dtS * (plusFlux[0] + minusFlux[0]);
  result.closureResidualJ = final + dissipated - initial - injection;
  result.status = Transport::Status::Ok();
  return result;
}

Transport::Status ValidateProfile(const AnalyticProfile& profile) {
  if (profile.id.empty() || profile.version.empty() || profile.units.empty() ||
      profile.provenance.empty() || profile.knots.size() < 2)
    return Error(Transport::StatusCode::InvalidArgument,
                 "analytic profile metadata or knot table is incomplete");
  for (std::size_t i = 0; i < profile.knots.size(); ++i) {
    const ProfileKnot& k = profile.knots[i];
    if (!Finite(k.radiusM) || k.radiusM <= 0.0 || !Finite(k.value) ||
        !Finite(k.derivativePerM) ||
        (i && k.radiusM <= profile.knots[i - 1].radiusM))
      return Error(Transport::StatusCode::InvalidArgument,
                   "analytic profile knots must be finite and ordered");
  }
  return Transport::Status::Ok();
}

ProfileValue EvaluateProfile(const AnalyticProfile& profile, double radiusM) {
  ProfileValue result;
  result.status = ValidateProfile(profile);
  if (!result.status.ok()) return result;
  if (!Finite(radiusM) || radiusM <= 0.0) {
    result.status = Error(Transport::StatusCode::InvalidArgument,
                          "profile radius must be finite and positive");
    return result;
  }
  const ProfileKnot& first = profile.knots.front();
  const ProfileKnot& last = profile.knots.back();
  if (radiusM < first.radiusM || radiusM > last.radiusM) {
    if (profile.extrapolation == ProfileExtrapolation::Reject) {
      result.status = Error(Transport::StatusCode::OutOfDomain,
                            "profile radius is outside its calibrated domain");
      return result;
    }
    const ProfileKnot& endpoint = radiusM < first.radiusM ? first : last;
    result.extrapolated = true;
    if (profile.extrapolation == ProfileExtrapolation::ConstantEndpoint) {
      result.value = endpoint.value;
      result.derivativePerM = 0.0;
    } else {
      if (endpoint.value == 0.0) {
        result.status = Error(Transport::StatusCode::InvalidCoefficient,
                              "power-law extrapolation requires nonzero endpoint");
        return result;
      }
      const double exponent = endpoint.radiusM * endpoint.derivativePerM /
                              endpoint.value;
      result.value = endpoint.value * std::pow(radiusM / endpoint.radiusM,
                                               exponent);
      result.derivativePerM = exponent * result.value / radiusM;
    }
    result.status = Transport::Status::Ok();
    return result;
  }
  std::size_t upper = 1;
  while (upper < profile.knots.size() &&
         profile.knots[upper].radiusM < radiusM) ++upper;
  if (upper == profile.knots.size()) upper = profile.knots.size() - 1;
  const ProfileKnot& a = profile.knots[upper - 1];
  const ProfileKnot& b = profile.knots[upper];
  const double h = b.radiusM - a.radiusM;
  const double t = (radiusM - a.radiusM) / h;
  const double h00 = 2.0*t*t*t - 3.0*t*t + 1.0;
  const double h10 = t*t*t - 2.0*t*t + t;
  const double h01 = -2.0*t*t*t + 3.0*t*t;
  const double h11 = t*t*t - t*t;
  result.value = h00*a.value + h10*h*a.derivativePerM +
                 h01*b.value + h11*h*b.derivativePerM;
  const double dh00 = (6.0*t*t - 6.0*t) / h;
  const double dh10 = 3.0*t*t - 4.0*t + 1.0;
  const double dh01 = (-6.0*t*t + 6.0*t) / h;
  const double dh11 = 3.0*t*t - 2.0*t;
  result.derivativePerM = dh00*a.value + dh10*a.derivativePerM +
                          dh01*b.value + dh11*b.derivativePerM;
  result.status = Transport::Status::Ok();
  return result;
}

std::string ProfileManifest(const AnalyticProfile& profile) {
  if (!ValidateProfile(profile).ok()) return std::string();
  std::ostringstream out;
  out << "SEP_ANALYTIC_PROFILE 1\n"
      << "id " << profile.id << "\nversion " << profile.version
      << "\nunits " << profile.units << "\nprovenance "
      << profile.provenance << "\nextrapolation "
      << static_cast<int>(profile.extrapolation) << "\nknots "
      << profile.knots.size() << '\n' << std::setprecision(17);
  for (std::size_t i = 0; i < profile.knots.size(); ++i)
    out << profile.knots[i].radiusM << ' ' << profile.knots[i].value << ' '
        << profile.knots[i].derivativePerM << '\n';
  return out.str();
}

}  // namespace PhysicsExtensions
}  // namespace SEP
