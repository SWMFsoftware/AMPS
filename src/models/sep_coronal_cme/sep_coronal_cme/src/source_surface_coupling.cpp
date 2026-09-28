#include "sep_coronal_cme/source_surface_coupling.h"

#include "sep_coronal_cme/constants.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <utility>

namespace SEP { namespace CoronalCME { namespace {

bool Finite(double value) { return std::isfinite(value); }

double Clamp(double value, double lower, double upper) {
  return std::max(lower, std::min(upper, value));
}

double FactorialRatio(int lower, int upper) {
  // Form lower!/upper! as a product so high-degree normalization does not
  // overflow merely because the two factorials are individually enormous.
  double ratio = 1.0;
  for (int value = lower + 1; value <= upper; ++value) ratio /= value;
  return ratio;
}

double AssociatedLegendre(int degree, int order, double x) {
  double pmm = 1.0;
  if (order > 0) {
    const double root = std::sqrt(std::max(0.0, 1.0 - x * x));
    double factor = 1.0;
    for (int index = 1; index <= order; ++index) {
      pmm *= -factor * root;
      factor += 2.0;
    }
  }
  if (degree == order) return pmm;
  double current = x * (2 * order + 1) * pmm;
  if (degree == order + 1) return current;
  double previous = pmm;
  for (int value = order + 2; value <= degree; ++value) {
    const double next = ((2 * value - 1) * x * current -
        (value + order - 1) * previous) / (value - order);
    previous = current;
    current = next;
  }
  return current;
}

struct Basis {
  double value = 0.0;
  double theta = 0.0;
  double theta2 = 0.0;
  double phi = 0.0;
  double thetaPhi = 0.0;
  double phi2 = 0.0;
};

Basis RealBasis(int degree, int order, bool sine, double theta, double phi) {
  const double x = std::cos(theta);
  const double sinTheta = std::sin(theta);
  const double safeSin = std::abs(sinTheta) > 1.0e-12 ? sinTheta : 1.0e-12;
  const double p = AssociatedLegendre(degree, order, x);
  const double previous = degree > order
      ? AssociatedLegendre(degree - 1, order, x) : 0.0;
  const double pTheta = (degree * x * p - (degree + order) * previous) /
      safeSin;
  const double pTheta2 = -x / safeSin * pTheta -
      (degree * (degree + 1.0) - order * order / (safeSin * safeSin)) * p;
  const double normalization = std::sqrt(
      (2.0 * degree + 1.0) / (4.0 * Constants::kPi) *
      (order == 0 ? 1.0 : 2.0) *
      FactorialRatio(degree - order, degree + order));
  const double angle = order * phi;
  const double trig = sine ? std::sin(angle) : std::cos(angle);
  const double trigPhi = order *
      (sine ? std::cos(angle) : -std::sin(angle));
  return {normalization * p * trig,
          normalization * pTheta * trig,
          normalization * pTheta2 * trig,
          normalization * p * trigPhi,
          normalization * pTheta * trigPhi,
          -order * order * normalization * p * trig};
}

Vec3 Scale(Vec3 value, double scale) { return scale * value; }

}  // namespace

Core::Result<FiniteShellScs> FiniteShellScs::Create(
    double interfaceRadiusM, double radialOuterRadiusM,
    std::vector<HarmonicCoefficient> coefficients) {
  if (!(Finite(interfaceRadiusM) && Finite(radialOuterRadiusM) &&
        interfaceRadiusM > 0.0 && radialOuterRadiusM > interfaceRadiusM)) {
    return Core::Result<FiniteShellScs>::Failure(
        Core::StatusCode::InvalidConfiguration,
        "finite-shell SCS requires 0<R_i<R_scs");
  }
  std::map<std::pair<int, int>, bool> seen;
  bool hasMonopole = false;
  int maximumDegree = 0;
  double monopoleAmplitude = 0.0;
  for (const auto& coefficient : coefficients) {
    if (coefficient.degree < 0 || coefficient.order < 0 ||
        coefficient.order > coefficient.degree ||
        !Finite(coefficient.cosineT) || !Finite(coefficient.sineT) ||
        (coefficient.order == 0 && coefficient.sineT != 0.0) ||
        !seen.emplace(std::make_pair(coefficient.degree, coefficient.order),
                      true).second) {
      return Core::Result<FiniteShellScs>::Failure(
          Core::StatusCode::InvalidConfiguration,
          "invalid or duplicate SCS real-harmonic coefficient");
    }
    if (coefficient.degree == 0) {
      hasMonopole = coefficient.cosineT > 0.0;
      monopoleAmplitude = coefficient.cosineT;
    }
    maximumDegree = std::max(maximumDegree, coefficient.degree);
  }
  if (!hasMonopole) {
    return Core::Result<FiniteShellScs>::Failure(
        Core::StatusCode::InvalidConfiguration,
        "unsigned SCS boundary requires a positive degree-zero mode");
  }

  // The constrained fit itself is normally produced by the preprocessing
  // layer.  The shared provider nevertheless verifies its defining invariant
  // on a deterministic grid oversampling every retained angular wavelength.
  // A negative lobe is rejected as one invalid asset; it is never clipped at
  // individual cells, which would change flux and introduce grid-dependent
  // divergence.  The tolerance covers only roundoff relative to the monopole.
  const int thetaSamples = std::max(16, 8 * (maximumDegree + 1));
  const int phiSamples = std::max(32, 16 * (maximumDegree + 1));
  const double negativityTolerance = 1.0e-12 * monopoleAmplitude /
      std::sqrt(4.0 * Constants::kPi);
  for (int i = 0; i <= thetaSamples; ++i) {
    const double theta = Constants::kPi * i / thetaSamples;
    for (int j = 0; j < phiSamples; ++j) {
      const double phi = 2.0 * Constants::kPi * j / phiSamples;
      double boundary = 0.0;
      for (const auto& coefficient : coefficients) {
        boundary += coefficient.cosineT *
            RealBasis(coefficient.degree, coefficient.order, false,
                      theta, phi).value;
        if (coefficient.order > 0) {
          boundary += coefficient.sineT *
              RealBasis(coefficient.degree, coefficient.order, true,
                        theta, phi).value;
        }
      }
      if (boundary < -negativityTolerance) {
        return Core::Result<FiniteShellScs>::Failure(
            Core::StatusCode::DataIntegrityFailure,
            "SCS constrained inner-boundary fit has a negative lobe");
      }
    }
  }
  FiniteShellScs result;
  result.interfaceRadiusM_ = interfaceRadiusM;
  result.outerRadiusM_ = radialOuterRadiusM;
  result.coefficients_ = std::move(coefficients);
  return Core::Result<FiniteShellScs>::Success(std::move(result));
}

Core::Result<SphericalField> FiniteShellScs::Evaluate(
    double radiusM, double thetaRad, double phiRad) const {
  if (!(Finite(radiusM) && Finite(thetaRad) && Finite(phiRad) &&
        radiusM >= interfaceRadiusM_ && radiusM <= outerRadiusM_ &&
        thetaRad >= 0.0 && thetaRad <= Constants::kPi)) {
    return Core::Result<SphericalField>::Failure(
        Core::StatusCode::OutOfDomain,
        "SCS query is outside the finite spherical shell");
  }
  const double sinTheta = std::sin(thetaRad);
  if (std::abs(sinTheta) < 1.0e-10) {
    return Core::Result<SphericalField>::Failure(
        Core::StatusCode::OutOfDomain,
        "spherical SCS components are singular at a coordinate pole");
  }
  SphericalField field;
  for (const auto& coefficient : coefficients_) {
    const int degree = coefficient.degree;
    double radialPotential = 0.0;
    double radialDerivative = 0.0;
    double radialSecondDerivative = 0.0;
    if (degree == 0) {
      // Psi_00=d(1/r-1/R_scs), d=h_00 R_i^2.  The outer value is
      // the fixed monopole gauge; the magnetic field is gauge independent.
      const double d = interfaceRadiusM_ * interfaceRadiusM_;
      radialPotential = d * (1.0 / radiusM - 1.0 / outerRadiusM_);
      radialDerivative = -d / (radiusM * radiusM);
      radialSecondDerivative = 2.0 * d /
          (radiusM * radiusM * radiusM);
    } else {
      const double outerPower = std::pow(outerRadiusM_, 2 * degree + 1);
      const double denominator = degree *
          std::pow(interfaceRadiusM_, degree - 1) +
          (degree + 1.0) * outerPower *
          std::pow(interfaceRadiusM_, -degree - 2);
      const double cPerBoundaryAmplitude = -1.0 / denominator;
      radialPotential = cPerBoundaryAmplitude *
          (std::pow(radiusM, degree) -
           outerPower * std::pow(radiusM, -degree - 1));
      radialDerivative = cPerBoundaryAmplitude *
          (degree * std::pow(radiusM, degree - 1) +
           (degree + 1.0) * outerPower *
               std::pow(radiusM, -degree - 2));
      radialSecondDerivative = cPerBoundaryAmplitude *
          (degree * (degree - 1.0) * std::pow(radiusM, degree - 2) -
           (degree + 1.0) * (degree + 2.0) * outerPower *
               std::pow(radiusM, -degree - 3));
    }

    for (int part = 0; part < (coefficient.order == 0 ? 1 : 2); ++part) {
      const double amplitude = part == 0 ? coefficient.cosineT
                                         : coefficient.sineT;
      const Basis basis = RealBasis(degree, coefficient.order, part == 1,
                                    thetaRad, phiRad);
      field.brT += -amplitude * radialDerivative * basis.value;
      field.bThetaT += -amplitude * radialPotential / radiusM * basis.theta;
      field.bPhiT += -amplitude * radialPotential /
          (radiusM * sinTheta) * basis.phi;

      field.derivative[0][0] +=
          -amplitude * radialSecondDerivative * basis.value;
      field.derivative[0][1] +=
          -amplitude * radialDerivative * basis.theta;
      field.derivative[0][2] +=
          -amplitude * radialDerivative * basis.phi;
      const double radialTangential = radialDerivative / radiusM -
          radialPotential / (radiusM * radiusM);
      field.derivative[1][0] +=
          -amplitude * radialTangential * basis.theta;
      field.derivative[1][1] +=
          -amplitude * radialPotential / radiusM * basis.theta2;
      field.derivative[1][2] +=
          -amplitude * radialPotential / radiusM * basis.thetaPhi;
      field.derivative[2][0] += -amplitude * radialTangential /
          sinTheta * basis.phi;
      field.derivative[2][1] += -amplitude * radialPotential / radiusM *
          (basis.thetaPhi / sinTheta -
           basis.phi * std::cos(thetaRad) / (sinTheta * sinTheta));
      field.derivative[2][2] += -amplitude * radialPotential /
          (radiusM * sinTheta) * basis.phi2;
    }
  }
  return Core::Result<SphericalField>::Success(field);
}

Core::Result<double> ScsHarmonicAttenuation(
    int degree, double outerToInnerRadius) {
  if (degree < 1 || !Finite(outerToInnerRadius) ||
      outerToInnerRadius <= 1.0) {
    return Core::Result<double>::Failure(
        Core::StatusCode::InvalidConfiguration,
        "SCS attenuation requires degree>=1 and R_scs/R_i>1");
  }
  const double numerator = (2.0 * degree + 1.0) *
      std::pow(outerToInnerRadius, degree + 1);
  const double denominator = degree + (degree + 1.0) *
      std::pow(outerToInnerRadius, 2 * degree + 1);
  return Core::Result<double>::Success(numerator / denominator);
}

Core::Result<ScsSpectralDiagnostics> EvaluateScsSpectrum(
    const std::vector<HarmonicCoefficient>& coefficients,
    double outerToInnerRadius) {
  double monopolePower = 0.0;
  int maximumDegree = 0;
  ScsSpectralDiagnostics result;
  for (const auto& coefficient : coefficients) {
    if (coefficient.degree < 0 || coefficient.order < 0 ||
        coefficient.order > coefficient.degree ||
        !Finite(coefficient.cosineT) || !Finite(coefficient.sineT)) {
      return Core::Result<ScsSpectralDiagnostics>::Failure(
          Core::StatusCode::InvalidConfiguration,
          "invalid coefficient in SCS spectral diagnostic");
    }
    maximumDegree = std::max(maximumDegree, coefficient.degree);
  }
  result.attenuationByDegree.assign(maximumDegree + 1, 1.0);
  for (int degree = 1; degree <= maximumDegree; ++degree) {
    const auto attenuation = ScsHarmonicAttenuation(
        degree, outerToInnerRadius);
    if (!attenuation.ok()) {
      return Core::Result<ScsSpectralDiagnostics>::Failure(
          attenuation.status.code, attenuation.status.message);
    }
    result.attenuationByDegree[degree] = attenuation.value;
  }
  double zonalOuterPower = 0.0;
  for (const auto& coefficient : coefficients) {
    const double power = coefficient.cosineT * coefficient.cosineT +
        coefficient.sineT * coefficient.sineT;
    if (coefficient.degree == 0) {
      monopolePower += power;
      continue;
    }
    result.innerNonMonopolePowerT2 += power;
    const double outerPower = power * std::pow(
        result.attenuationByDegree[coefficient.degree], 2);
    result.outerNonMonopolePowerT2 += outerPower;
    if (coefficient.order == 0) zonalOuterPower += outerPower;
  }
  const double denominator = monopolePower + zonalOuterPower;
  if (!(denominator > 0.0)) {
    return Core::Result<ScsSpectralDiagnostics>::Failure(
        Core::StatusCode::InvalidState,
        "SCS outer zonal fraction has no positive denominator");
  }
  result.outerZonalFraction = zonalOuterPower / denominator;
  return Core::Result<ScsSpectralDiagnostics>::Success(std::move(result));
}

Core::Result<OneSidedMagneticField> RestoreSector(
    Vec3 unsignedFieldT, MagneticSector sector, bool onIdealSheet,
    bool sideWasSupplied, InterfaceSide side) {
  if (!(Finite(unsignedFieldT.x) && Finite(unsignedFieldT.y) &&
        Finite(unsignedFieldT.z)) || Norm(unsignedFieldT) == 0.0) {
    return Core::Result<OneSidedMagneticField>::Failure(
        Core::StatusCode::InvalidState,
        "sector restoration requires a finite nonzero unsigned field");
  }
  if (onIdealSheet && !sideWasSupplied) {
    return Core::Result<OneSidedMagneticField>::Failure(
        Core::StatusCode::OutOfDomain,
        "ideal HCS query requires an explicit one-sided trace");
  }
  const double sign = sector == MagneticSector::Positive ? 1.0 : -1.0;
  return Core::Result<OneSidedMagneticField>::Success(
      {Scale(unsignedFieldT, sign), sector, side});
}

Core::Status ValidateIdealHcsTransport(
    IdealHcsTransport requested, bool qualifiedFiniteSheetProvider) {
  if (requested == IdealHcsTransport::CrossSectorOrDrift &&
      !qualifiedFiniteSheetProvider) {
    return Core::Status::Failure(Core::StatusCode::UnsupportedCapability,
        "cross-sector/HCS-drift transport requires a qualified finite sheet");
  }
  return Core::Status::Success();
}

Core::Status ValidateCouplingMode(
    CouplingMode mode, bool productionIntent, double interfaceRadiusM,
    double scsOuterRadiusM, double transitionWidthM) {
  if (!(interfaceRadiusM > 0.0 && scsOuterRadiusM >= interfaceRadiusM &&
        transitionWidthM >= 0.0)) {
    return Core::Status::Failure(Core::StatusCode::InvalidConfiguration,
        "coupling radii and transition width are invalid");
  }
  if (productionIntent && mode != CouplingMode::ProductionResolvedTransition) {
    return Core::Status::Failure(Core::StatusCode::UnsupportedCapability,
        "production requires the resolved PFSS/SCS transition");
  }
  if (mode == CouplingMode::ProductionResolvedTransition &&
      (!(scsOuterRadiusM > interfaceRadiusM) || transitionWidthM <= 0.0 ||
       transitionWidthM > scsOuterRadiusM - interfaceRadiusM)) {
    return Core::Status::Failure(Core::StatusCode::InvalidConfiguration,
        "resolved transition must have 0<w<=R_scs-R_i");
  }
  if ((mode == CouplingMode::SharpInterfaceVerification ||
       mode == CouplingMode::NoScsVerification) && transitionWidthM != 0.0) {
    return Core::Status::Failure(Core::StatusCode::InvalidConfiguration,
        "verification-only sharp/no-SCS coupling requires zero width");
  }
  if (mode == CouplingMode::NoScsVerification &&
      interfaceRadiusM != scsOuterRadiusM) {
    return Core::Status::Failure(Core::StatusCode::InvalidConfiguration,
        "no-SCS verification derives R_i=R_scs=R_b");
  }
  return Core::Status::Success();
}

Core::Result<InterfaceMagneticBalance> EvaluateMagneticInterface(
    Vec3 pfssFieldT, Vec3 scsFieldT, Vec3 outwardNormal,
    double normalAbsoluteToleranceT) {
  const Vec3 normal = Unit(outwardNormal);
  if (Norm(normal) == 0.0 || normalAbsoluteToleranceT < 0.0) {
    return Core::Result<InterfaceMagneticBalance>::Failure(
        Core::StatusCode::InvalidConfiguration,
        "magnetic interface requires a normal and nonnegative tolerance");
  }
  InterfaceMagneticBalance result;
  result.normalJumpT = Dot(scsFieldT - pfssFieldT, normal);
  result.surfaceCurrentAPerM =
      Cross(normal, scsFieldT - pfssFieldT) / Constants::kVacuumPermeabilityHPerM;
  const double product = Norm(pfssFieldT) * Norm(scsFieldT);
  if (product > 0.0) {
    result.kinkAngleRad = std::acos(Clamp(
        Dot(pfssFieldT, scsFieldT) / product, -1.0, 1.0));
  }
  if (std::abs(result.normalJumpT) > normalAbsoluteToleranceT) {
    return Core::Result<InterfaceMagneticBalance>::Failure(
        Core::StatusCode::InvalidState,
        "PFSS/SCS conservative normal-flux trace does not close");
  }
  return Core::Result<InterfaceMagneticBalance>::Success(result);
}

Core::Result<TransitionBlend> BlendVectorPotentialFields(
    double radiusM, double innerRadiusM, double outerRadiusM,
    Vec3 radialUnit, Vec3 pfssFieldT, Vec3 scsFieldT,
    Vec3 pfssPotentialTm, Vec3 scsPotentialTm) {
  if (!(Finite(radiusM) && innerRadiusM > 0.0 && outerRadiusM > innerRadiusM &&
        radiusM >= innerRadiusM && radiusM <= outerRadiusM &&
        Norm(radialUnit) > 0.0)) {
    return Core::Result<TransitionBlend>::Failure(
        Core::StatusCode::InvalidConfiguration,
        "transition query requires R_i<=r<=R_i+w and a radial unit vector");
  }
  const double width = outerRadiusM - innerRadiusM;
  const double z = (radiusM - innerRadiusM) / width;
  const double chi = 10.0 * z * z * z - 15.0 * std::pow(z, 4) +
      6.0 * std::pow(z, 5);
  const double dChi = (30.0 * z * z - 60.0 * z * z * z +
      30.0 * std::pow(z, 4)) / width;
  const Vec3 gradientChi = dChi * Unit(radialUnit);
  const Vec3 field = (1.0 - chi) * pfssFieldT + chi * scsFieldT +
      Cross(gradientChi, scsPotentialTm - pfssPotentialTm);
  return Core::Result<TransitionBlend>::Success({chi, dChi, field});
}

Core::Status ValidateSignedPotentialQualification(
    const SignedPotentialQualification& qualification) {
  if (!qualification.fluxBalanced) {
    return Core::Status::Failure(Core::StatusCode::InvalidState,
        "signed Mie potential requires zero net magnetic flux");
  }
  if (!qualification.fixedMieGauge) {
    return Core::Status::Failure(Core::StatusCode::InvalidState,
        "PFSS and SCS potentials require the common fixed Mie gauge");
  }
  if (!qualification.sheetTangentialTraceContinuous) {
    return Core::Status::Failure(Core::StatusCode::InvalidState,
        "piecewise potential has an incompatible sheet tangential trace");
  }
  return Core::Status::Success();
}

Core::Result<LongitudeMap> IntegrateLongitudeMap(
    double sourceLongitudeRad, double sourceRadiusM, double radiusM,
    int radialSteps, const WindingRate& windingRatePerM,
    const WindingRate& longitudeDerivativePerM,
    double minimumForwardJacobian) {
  if (!(Finite(sourceLongitudeRad) && Finite(sourceRadiusM) && Finite(radiusM) &&
        sourceRadiusM > 0.0 && radiusM > 0.0 && radialSteps > 0 &&
        minimumForwardJacobian > 0.0 && windingRatePerM &&
        longitudeDerivativePerM)) {
    return Core::Result<LongitudeMap>::Failure(
        Core::StatusCode::InvalidConfiguration,
        "longitude map requires finite radii, callbacks, steps, and A_phi_min>0");
  }
  double phi = sourceLongitudeRad;
  double aPhi = 1.0;
  const double step = (radiusM - sourceRadiusM) / radialSteps;
  auto derivative = [&](double r, double mappedPhi, double a) {
    return std::pair<double, double>{-windingRatePerM(r, mappedPhi),
        -longitudeDerivativePerM(r, mappedPhi) * a};
  };
  for (int index = 0; index < radialSteps; ++index) {
    const double r = sourceRadiusM + index * step;
    const auto k1 = derivative(r, phi, aPhi);
    const auto k2 = derivative(r + 0.5 * step,
        phi + 0.5 * step * k1.first,
        aPhi + 0.5 * step * k1.second);
    const auto k3 = derivative(r + 0.5 * step,
        phi + 0.5 * step * k2.first,
        aPhi + 0.5 * step * k2.second);
    const auto k4 = derivative(r + step, phi + step * k3.first,
                               aPhi + step * k3.second);
    phi += step * (k1.first + 2.0 * k2.first + 2.0 * k3.first +
                   k4.first) / 6.0;
    aPhi += step * (k1.second + 2.0 * k2.second + 2.0 * k3.second +
                    k4.second) / 6.0;
    if (!Finite(phi) || !Finite(aPhi) || aPhi <= minimumForwardJacobian) {
      return Core::Result<LongitudeMap>::Failure(
          Core::StatusCode::InvalidState,
          "rotating-footpoint map folded or crossed its A_phi floor");
    }
  }
  return Core::Result<LongitudeMap>::Success(
      {sourceLongitudeRad, phi, aPhi, 1.0 / aPhi,
       aPhi - minimumForwardJacobian});
}

Core::Result<ParkerMappedState> MapParkerState(
    const LongitudeMap& map, double sourceRadiusM, double radiusM,
    double thetaRad, double sourceRadialFieldT, double sourceDensityKgM3,
    double sourceRadialSpeedMPerS, double radialSpeedMPerS,
    double windingRatePerM) {
  if (!(sourceRadiusM > 0.0 && radiusM >= sourceRadiusM &&
        thetaRad >= 0.0 && thetaRad <= Constants::kPi &&
        map.forwardJacobian > 0.0 && map.inverseJacobian > 0.0 &&
        sourceDensityKgM3 > 0.0 && sourceRadialSpeedMPerS > 0.0 &&
        radialSpeedMPerS > 0.0 && Finite(sourceRadialFieldT) &&
        Finite(windingRatePerM))) {
    return Core::Result<ParkerMappedState>::Failure(
        Core::StatusCode::InvalidConfiguration,
        "Parker mapping requires a valid unfolded map and positive plasma state");
  }
  ParkerMappedState result;
  const double radialScale = std::pow(sourceRadiusM / radiusM, 2);
  const double br = radialScale * sourceRadialFieldT * map.inverseJacobian;
  const double bPhi = -radiusM * std::sin(thetaRad) * windingRatePerM * br;
  result.magneticSphericalT = {br, 0.0, bPhi};
  result.densityKgM3 = sourceDensityKgM3 *
      sourceRadialSpeedMPerS / radialSpeedMPerS * radialScale *
      map.inverseJacobian;
  result.fieldAlignedSpeedMPerS = radialSpeedMPerS *
      Norm(result.magneticSphericalT) / std::abs(br);
  result.mappedMassFluxKgPerSPerSr = radiusM * radiusM *
      result.densityKgM3 * radialSpeedMPerS;
  return Core::Result<ParkerMappedState>::Success(result);
}

Core::Result<PlasmaSheetNormalization> BuildPlasmaSheetNormalization(
    double neutralLineDistanceRad, double centralContrast,
    double angularWidthRad, PlasmaSheetAuthority authority,
    bool fixedTemperature) {
  if (!(Finite(neutralLineDistanceRad) && centralContrast >= 1.0 &&
        angularWidthRad > 0.0 && fixedTemperature)) {
    return Core::Result<PlasmaSheetNormalization>::Failure(
        Core::StatusCode::InvalidConfiguration,
        "plasma sheet requires C>=1, width>0, and fixed-temperature EOS");
  }
  const double ratio = neutralLineDistanceRad / angularWidthRad;
  const double contrast = 1.0 + (centralContrast - 1.0) *
      std::exp(-ratio * ratio);
  PlasmaSheetNormalization result;
  result.contrast = contrast;
  if (authority == PlasmaSheetAuthority::BaseDensity) {
    // Scaling both rho and p preserves p/rho and therefore temperature.
    result.baseDensityMultiplier = contrast;
    result.basePressureMultiplier = contrast;
  } else {
    // An outer/mass-loading authority is modified exactly once; changing the
    // base density as well would overdetermine the transonic tube solution.
    result.massLoadingMultiplier = contrast;
  }
  return Core::Result<PlasmaSheetNormalization>::Success(result);
}

Core::Result<TransitionDiagnostics> EvaluateTransitionDiagnostics(
    Vec3 plusFieldT, Vec3 minusFieldT, Vec3 sheetNormal,
    double patchAreaM2, double openFluxWb, double minimumFieldT,
    double jumpSupportThreshold, double normalTraceTolerance) {
  const Vec3 normal = Unit(sheetNormal);
  const double plusMagnitude = Norm(plusFieldT);
  const double minusMagnitude = Norm(minusFieldT);
  if (!(Norm(normal) > 0.0 && patchAreaM2 > 0.0 && openFluxWb > 0.0 &&
        minimumFieldT > 0.0 && jumpSupportThreshold >= 0.0 &&
        jumpSupportThreshold <= 1.0 && normalTraceTolerance >= 0.0)) {
    return Core::Result<TransitionDiagnostics>::Failure(
        Core::StatusCode::InvalidConfiguration,
        "D9 transition diagnostic received invalid physical bounds");
  }
  if (plusMagnitude < minimumFieldT || minusMagnitude < minimumFieldT) {
    return Core::Result<TransitionDiagnostics>::Failure(
        Core::StatusCode::InvalidState,
        "transition direction is inapplicable in a weak-field trace");
  }
  const double normalPlus = Dot(plusFieldT, normal);
  const double normalMinus = Dot(minusFieldT, normal);
  TransitionDiagnostics result;
  result.normalTraceRelativeError = std::abs(normalPlus - normalMinus) *
      patchAreaM2 / openFluxWb;
  if (result.normalTraceRelativeError > normalTraceTolerance) {
    return Core::Result<TransitionDiagnostics>::Failure(
        Core::StatusCode::InvalidState,
        "one-sided transition normal traces exceed the mortar tolerance");
  }
  const double mortarNormal = 0.5 * (normalPlus + normalMinus);
  result.absoluteCrossingFluxWb = std::abs(mortarNormal) * patchAreaM2;
  result.netCrossingFluxWb = mortarNormal * patchAreaM2;
  result.jumpSupport = Norm(plusFieldT - minusFieldT) /
      (plusMagnitude + minusMagnitude);
  result.incidencePlus = normalPlus / plusMagnitude;
  result.incidenceMinus = normalMinus / minusMagnitude;
  if (result.jumpSupport >= jumpSupportThreshold) {
    result.antipodalityApplicable = true;
    result.antipodalityDefectRad = std::acos(Clamp(
        -Dot(plusFieldT, minusFieldT) /
            (plusMagnitude * minusMagnitude), -1.0, 1.0));
  }
  return Core::Result<TransitionDiagnostics>::Success(result);
}

Core::Result<BudgetRatio> EvaluateBudgetRatio(
    double excluded, double counterfactual, double maximumFraction) {
  if (!(Finite(excluded) && Finite(counterfactual) &&
        Finite(maximumFraction) && excluded >= 0.0 &&
        counterfactual >= 0.0 && excluded <= counterfactual &&
        maximumFraction >= 0.0 && maximumFraction <= 1.0)) {
    return Core::Result<BudgetRatio>::Failure(
        Core::StatusCode::InvalidConfiguration,
        "budget measures require 0<=excluded<=counterfactual and a unit bound");
  }
  BudgetRatio result;
  result.numerator = excluded;
  result.denominator = counterfactual;
  if (counterfactual == 0.0) {
    result.validity = MeasureValidity::InapplicableZeroDenominator;
    return Core::Result<BudgetRatio>::Success(result);
  }
  result.ratio = excluded / counterfactual;
  result.withinBound = result.ratio <= maximumFraction;
  return Core::Result<BudgetRatio>::Success(result);
}

BudgetRatio EvaluatePointClearance(bool intersectsClearance) {
  BudgetRatio result;
  result.validity = intersectsClearance
      ? MeasureValidity::RejectedClearance
      : MeasureValidity::InapplicablePointMeasure;
  result.withinBound = !intersectsClearance;
  return result;
}

Core::Result<OpenFluxLifecycle> BuildOpenFluxLifecycle(
    double unscaledFluxWb, double targetFluxWb, double rebuiltPassBFluxWb,
    double relativeTolerance, const std::string& constructionChecksum,
    const std::string& qualificationChecksum) {
  if (!(unscaledFluxWb > 0.0 && targetFluxWb > 0.0 &&
        rebuiltPassBFluxWb > 0.0 && relativeTolerance >= 0.0) ||
      constructionChecksum.empty() || qualificationChecksum.empty() ||
      constructionChecksum == qualificationChecksum) {
    return Core::Result<OpenFluxLifecycle>::Failure(
        Core::StatusCode::DataIntegrityFailure,
        "open-flux lifecycle requires positive fluxes and disjoint checksums");
  }
  OpenFluxLifecycle result;
  result.passAScale = targetFluxWb / unscaledFluxWb;
  result.passBFluxWb = rebuiltPassBFluxWb;
  result.constructionChecksum = constructionChecksum;
  result.qualificationChecksum = qualificationChecksum;
  result.qualificationPassed = std::abs(rebuiltPassBFluxWb - targetFluxWb) <=
      relativeTolerance * targetFluxWb;
  if (!result.qualificationPassed) {
    return Core::Result<OpenFluxLifecycle>::Failure(
        Core::StatusCode::InvalidState,
        "rebuilt Pass-B composite failed independent open-flux qualification");
  }
  return Core::Result<OpenFluxLifecycle>::Success(std::move(result));
}

} }  // namespace SEP::CoronalCME
