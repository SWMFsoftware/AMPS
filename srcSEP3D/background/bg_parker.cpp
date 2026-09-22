#include "bg_parker.h"

#include "sep_background_snapshot.h"

#include <array>
#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>

namespace SEP3D {
namespace Background {
namespace {

constexpr double kMu0 = 4.0e-7 * Core::Const::kPi;

Core::Status Invalid(const std::string& message) {
  return Core::Status(Core::StatusCode::InvalidInput, message);
}

std::uint64_t Digest64(const std::string& text) {
  // FNV-1a is used only as a compact in-cell tag; the complete versioned
  // fingerprint string remains in SnapshotMetadata and the run manifest.
  std::uint64_t result = 1469598103934665603ULL;
  for (unsigned char value : text) {
    result ^= value;
    result *= 1099511628211ULL;
  }
  return result;
}

Core::Vec3 AxisCrossBasisColumn(const Core::Vec3& axis, int column) {
  if (column == 0) return {0.0, axis.z, -axis.y};
  if (column == 1) return {-axis.z, 0.0, axis.x};
  return {axis.y, -axis.x, 0.0};
}

double Component(const Core::Vec3& value, int index) {
  return index == 0 ? value.x : (index == 1 ? value.y : value.z);
}

Core::ParkerSpiralGeometry Geometry(const ParkerConfiguration& configuration) {
  Core::ParkerSpiralGeometry geometry;
  geometry.sourceRadiusM = configuration.sourceRadiusM;
  geometry.sourceLongitudeRad = configuration.sourceLongitudeRad;
  geometry.sourceColatitudeRad = configuration.sourceColatitudeRad;
  geometry.solarWindSpeedMPerS = configuration.solarWindSpeedMPerS;
  geometry.solarRotationRateRadPerS =
      configuration.solarRotationRateRadPerS;
  geometry.rotationAxis = configuration.rotationAxis;
  return geometry;
}

// Translate the provider-neutral SI record into the one canonical SWCME
// ambient state.  SWCME accepts total |B| at one AU, whereas srcSEP3D's Parker
// cross-check records radial Br at an arbitrary reference radius.  The
// conversion below is the exact inverse of SWCME's documented one-AU Parker
// normalization and therefore preserves both conventions without maintaining
// a second magnetic-field model.
swcme::solarwind::ConfigSI SwcmeConfiguration(
    const ParkerConfiguration& configuration) {
  swcme::solarwind::ConfigSI result;
  result.V_sw_m_s = configuration.solarWindSpeedMPerS;
  result.T_K = configuration.temperatureK;
  result.gamma_ad = configuration.adiabaticIndex;
  result.thermodynamic_closure = configuration.thermodynamicClosure;
  result.alpha_to_proton_ratio = configuration.alphaToProtonRatio;
  result.electron_T_K = configuration.electronTemperatureK;
  result.alpha_T_K = configuration.alphaTemperatureK;
  result.parker_radial_polarity = configuration.magneticPolarity;
  result.reference_sin_theta = configuration.referenceSinColatitude;
  result.solar_rotation_rate_rad_s =
      configuration.solarRotationRateRadPerS;
  result.parker_source_radius_m = configuration.sourceRadiusM;

  const double brAtOneAuT = configuration.radialFieldAtReferenceT *
      std::pow(configuration.referenceRadiusM / Core::Const::AU, 2);
  const double referencePitch = configuration.solarRotationRateRadPerS *
      (Core::Const::AU - configuration.sourceRadiusM) /
      configuration.solarWindSpeedMPerS *
      configuration.referenceSinColatitude;
  result.B1AU_T = brAtOneAuT *
      std::sqrt(1.0 + referencePitch * referencePitch);

  // SWCME's Leblanc helper is parameterized by n_e(1 AU).  A provider record
  // may state the same normalization at another explicit radius, so first
  // prepare a unit-one-AU profile and use its linear scaling to recover the
  // equivalent n_e(1 AU).  This is exact for the complete C2/r^2+C4/r^4+C6/r^6
  // law; replacing it with an r^-2 conversion would discard the near-Sun
  // Leblanc terms that this integration is meant to preserve.
  result.n1AU_m3 = 1.0;
  const swcme::solarwind::PreparedState unitDensity =
      swcme::solarwind::prepare(result);
  const double unitAtReference = swcme::solarwind::density_m3(
      unitDensity, configuration.densityReferenceRadiusM);
  result.n1AU_m3 =
      configuration.numberDensityAtReferenceM3 / unitAtReference;
  return result;
}

}  // namespace

AnalyticParkerProvider::AnalyticParkerProvider(
    const ParkerConfiguration& configuration)
    : configuration_(configuration) {
  metadata_.provider = ProviderKind::AnalyticParker;
  metadata_.ownership = StorageOwnership::ModelOwned;
  metadata_.coordinateFrame = configuration_.coordinateFrame;
  metadata_.providerIdentity = CanonicalName();
  metadata_.configurationFingerprint =
      SEP::Background::FingerprintConfiguration(ResolvedManifest());
  configurationDigest_ = Digest64(metadata_.configurationFingerprint);
}

Core::Status AnalyticParkerProvider::Validate() const {
  const double values[] = {
      configuration_.sourceRadiusM, configuration_.referenceRadiusM,
      configuration_.sourceLongitudeRad,
      configuration_.sourceColatitudeRad,
      configuration_.radialFieldAtReferenceT,
      configuration_.numberDensityAtReferenceM3,
      configuration_.densityReferenceRadiusM, configuration_.temperatureK,
      configuration_.adiabaticIndex, configuration_.alphaToProtonRatio,
      configuration_.electronTemperatureK,
      configuration_.alphaTemperatureK,
      configuration_.referenceSinColatitude,
      configuration_.solarWindSpeedMPerS,
      configuration_.solarRotationRateRadPerS,
      configuration_.validityCadenceS};
  for (double value : values) {
    if (!std::isfinite(value)) return Invalid("Parker configuration contains a non-finite value");
  }
  if (configuration_.sourceRadiusM <= 0.0 ||
      configuration_.referenceRadiusM <= configuration_.sourceRadiusM ||
      configuration_.radialFieldAtReferenceT <= 0.0 ||
      configuration_.numberDensityAtReferenceM3 <= 0.0 ||
      configuration_.densityReferenceRadiusM <= configuration_.sourceRadiusM ||
      configuration_.temperatureK <= 0.0 ||
      configuration_.adiabaticIndex <= 1.0 ||
      configuration_.alphaToProtonRatio < 0.0 ||
      configuration_.electronTemperatureK <= 0.0 ||
      configuration_.alphaTemperatureK <= 0.0 ||
      configuration_.referenceSinColatitude < 0.0 ||
      configuration_.referenceSinColatitude > 1.0 ||
      configuration_.solarWindSpeedMPerS <= 0.0 ||
      configuration_.solarRotationRateRadPerS < 0.0 ||
      configuration_.validityCadenceS <= 0.0 ||
      configuration_.sourceColatitudeRad < 0.0 ||
      configuration_.sourceColatitudeRad > Core::Const::kPi ||
      (configuration_.magneticPolarity != 1 &&
       configuration_.magneticPolarity != -1) ||
      configuration_.coordinateFrame.empty()) {
    return Invalid("Parker configuration is outside its physical range");
  }
  if (configuration_.thermodynamicClosure !=
          swcme::solarwind::ThermodynamicClosure::ProtonOnly &&
      configuration_.thermodynamicClosure !=
          swcme::solarwind::ThermodynamicClosure::MultiSpecies) {
    return Invalid("Parker thermodynamic closure is not a supported SWCME closure");
  }
  const double axisNorm = configuration_.rotationAxis.Norm();
  if (!std::isfinite(axisNorm) || axisNorm <= 0.0) {
    return Invalid("Parker rotation axis must be finite and non-zero");
  }
  const Core::Status geometry = Core::ValidateParkerGeometry(
      Geometry(configuration_));
  if (!geometry.ok()) return geometry;
  return Core::Status::OK();
}

Core::Status AnalyticParkerProvider::Prepare(double timeS) {
  const Core::Status valid = Validate();
  if (!valid.ok()) return valid;
  if (!std::isfinite(timeS)) return Invalid("Parker preparation time is not finite");
  // Build the complete candidate cache only after all divisors and physical
  // ranges have been validated.  Assigning it before metadata publication
  // keeps Prepare transactional if future SWCME preparation adds diagnostics.
  const swcme::solarwind::PreparedState candidateSolarWind =
      swcme::solarwind::prepare(SwcmeConfiguration(configuration_));
  if (!std::isfinite(candidateSolarWind.Br1AU_T) ||
      !std::isfinite(candidateSolarWind.C2) ||
      !std::isfinite(candidateSolarWind.C4) ||
      !std::isfinite(candidateSolarWind.C6)) {
    return Invalid("SWCME Parker/Leblanc preparation produced a non-finite state");
  }
  SnapshotMetadata candidate = metadata_;
  candidate.epochS = timeS;
  candidate.validFromS = timeS;
  candidate.validUntilS = timeS + configuration_.validityCadenceS;
  candidate.generation = metadata_.generation + 1;
  solarWind_ = candidateSolarWind;
  metadata_ = candidate;
  prepared_ = true;
  return Core::Status::OK();
}

const SnapshotMetadata* AnalyticParkerProvider::PreparedMetadata() const {
  return prepared_ ? &metadata_ : nullptr;
}

BackgroundSample AnalyticParkerProvider::Evaluate(
    const Core::Vec3& positionM) const {
  BackgroundSample sample;
  if (!prepared_) {
    sample.status = Core::Status(Core::StatusCode::SnapshotUnavailable,
                                 "Parker provider is not prepared");
    return sample;
  }
  const double radius = positionM.Norm();
  if (!std::isfinite(radius) || radius < configuration_.sourceRadiusM) {
    sample.status = Core::Status(Core::StatusCode::BackgroundInvalid,
                                 "position is inside the Parker source surface");
    return sample;
  }

  const Core::Vec3 rHat = positionM / radius;
  const Core::Vec3 axis = configuration_.rotationAxis.Normalized();
  const Core::Vec3 axisCrossX = axis.Cross(positionM);
  const std::array<double, 3> axisArray = {{axis.x, axis.y, axis.z}};
  const std::array<double, 3> radialArray = {{rHat.x, rHat.y, rHat.z}};
  const std::array<double, 3> canonicalField =
      swcme::solarwind::parker_field_cartesian(
          solarWind_, axisArray, radialArray, radius);
  sample.B = {canonicalField[0], canonicalField[1], canonicalField[2]};

  // The closed Cartesian derivative below differentiates the exact same
  // vector returned by swcme::solarwind::parker_field_cartesian.  Br1AU is
  // signed, so polarity enters once here and never alters mesh geometry.
  const double coefficient = solarWind_.Br1AU_T *
      Core::Const::AU * Core::Const::AU;
  const double winding = solarWind_.solar_rotation_rate_rad_s /
                         solarWind_.V_sw_m_s;
  const double inverseR3 = 1.0 / (radius * radius * radius);
  sample.absB = sample.B.Norm();
  if (!(sample.absB > 0.0) || !std::isfinite(sample.absB)) {
    sample.status = Core::Status(Core::StatusCode::BackgroundInvalid,
                                 "Parker field magnitude is invalid");
    return sample;
  }
  sample.bHat = sample.B / sample.absB;

  // Closed Cartesian derivative.  Writing the spiral as
  // B=C[x/r^3-k(r-r0)(a x x)/r^3] avoids singular spherical basis vectors at
  // the rotation axis and makes div(B)=0 to round-off.
  const double inverseR5 = inverseR3 / (radius * radius);
  const double f = (radius - configuration_.sourceRadiusM) * inverseR3;
  const double dfScale = -2.0 / std::pow(radius, 4) +
                         3.0 * configuration_.sourceRadiusM * inverseR5;
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      const double xi = Component(positionM, i);
      const double xj = Component(positionM, j);
      const double radialDerivative = (i == j ? inverseR3 : 0.0) -
                                      3.0 * xi * xj * inverseR5;
      const Core::Vec3 crossColumn = AxisCrossBasisColumn(axis, j);
      const double spiralDerivative = dfScale * xj *
          Component(axisCrossX, i) + f * Component(crossColumn, i);
      sample.gradB(i, j) = coefficient *
          (radialDerivative - winding * spiralDerivative);
    }
  }

  Core::Vec3 gradAbsB;
  for (int j = 0; j < 3; ++j) {
    double value = 0.0;
    for (int i = 0; i < 3; ++i)
      value += Component(sample.B, i) * sample.gradB(i, j);
    if (j == 0) gradAbsB.x = value / sample.absB;
    else if (j == 1) gradAbsB.y = value / sample.absB;
    else gradAbsB.z = value / sample.absB;
  }
  sample.divBhat = sample.gradB.Trace() / sample.absB -
      sample.B.Dot(gradAbsB) / (sample.absB * sample.absB);
  const double dLnBds = sample.bHat.Dot(gradAbsB) / sample.absB;
  sample.focusingLenM = std::fabs(dLnBds) > 0.0
      ? -1.0 / dLnBds : std::numeric_limits<double>::infinity();
  Core::Tensor3 gradBhat;
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      gradBhat(i, j) = sample.gradB(i, j) / sample.absB -
          Component(sample.B, i) * Component(gradAbsB, j) /
              (sample.absB * sample.absB);
  sample.curvature = gradBhat.Apply(sample.bHat);

  sample.U = solarWind_.V_sw_m_s * rHat;
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      sample.gradU(i, j) = solarWind_.V_sw_m_s / radius *
          ((i == j ? 1.0 : 0.0) -
           Component(rHat, i) * Component(rHat, j));
  sample.divU = sample.gradU.Trace();
  sample.fieldAlignedStrain =
      sample.gradU.DoubleContract(sample.bHat, sample.bHat);
  // SWCME defines n as electron density.  thermodynamic_state() then applies
  // the selected proton-only or charge-neutral electron/proton/alpha closure;
  // in particular, Alfvén speed uses the resulting mass density rather than
  // assuming rho=m_p*n_e in a multi-species plasma.
  sample.numberDensityM3 =
      swcme::solarwind::density_m3(solarWind_, radius);
  const swcme::solarwind::ThermodynamicState thermodynamics =
      swcme::solarwind::thermodynamic_state(
          solarWind_, sample.numberDensityM3);
  sample.temperatureK = solarWind_.T_K;
  sample.pressurePa = thermodynamics.pressure_Pa;
  sample.alfvenSpeedMpS = sample.absB /
      std::sqrt(kMu0 * thermodynamics.mass_density_kg_m3);
  sample.generation = metadata_.generation;
  sample.configurationDigest = configurationDigest_;
  sample.valid = true;
  sample.status = Core::Status::OK();
  return sample;
}

std::string AnalyticParkerProvider::ResolvedManifest() const {
  std::ostringstream out;
  out << std::setprecision(17) << std::scientific
      << "parker-provider-swcme-v2"
      << ";source_m=" << configuration_.sourceRadiusM
      << ";source_lon_rad=" << configuration_.sourceLongitudeRad
      << ";source_colat_rad=" << configuration_.sourceColatitudeRad
      << ";reference_m=" << configuration_.referenceRadiusM
      << ";Br_ref_T=" << configuration_.radialFieldAtReferenceT
      << ";n_ref_m-3=" << configuration_.numberDensityAtReferenceM3
      << ";n_reference_m=" << configuration_.densityReferenceRadiusM
      << ";temperature_K=" << configuration_.temperatureK
      << ";gamma=" << configuration_.adiabaticIndex
      << ";thermodynamic_closure="
      << swcme::solarwind::thermodynamic_closure_name(
             configuration_.thermodynamicClosure)
      << ";alpha_to_proton=" << configuration_.alphaToProtonRatio
      << ";electron_temperature_K="
      << configuration_.electronTemperatureK
      << ";alpha_temperature_K=" << configuration_.alphaTemperatureK
      << ";reference_sin_colatitude="
      << configuration_.referenceSinColatitude
      << ";wind_m_s=" << configuration_.solarWindSpeedMPerS
      << ";omega_rad_s=" << configuration_.solarRotationRateRadPerS
      << ";axis=" << configuration_.rotationAxis.x << ','
      << configuration_.rotationAxis.y << ',' << configuration_.rotationAxis.z
      << ";polarity=" << configuration_.magneticPolarity
      << ";cadence_s=" << configuration_.validityCadenceS
      << ";frame=" << configuration_.coordinateFrame;
  return out.str();
}

ProviderCapabilities AnalyticParkerProvider::Capabilities() const {
  ProviderCapabilities capabilities;
  capabilities.hasAnalyticGradB = true;
  capabilities.hasAnalyticDivBhat = true;
  capabilities.hasAnalyticCurvature = true;
  capabilities.hasAnalyticDivU = true;
  capabilities.hasFieldAlignedStrain = true;
  capabilities.hasPlasmaState = true;
  capabilities.supportsBatchEval = true;
  return capabilities;
}

}  // namespace Background
}  // namespace SEP3D
