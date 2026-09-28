#ifndef SEP_CORONAL_CME_SOURCE_SURFACE_COUPLING_H
#define SEP_CORONAL_CME_SOURCE_SURFACE_COUPLING_H

#include "sep_coronal_cme/pfss_harmonics.h"
#include "sep_coronal_cme/vector_math.h"
#include "sep_status.h"

#include <cstdint>
#include <functional>
#include <string>
#include <vector>

namespace SEP { namespace CoronalCME {

// Stage 3's finite-shell Schatten current-sheet (SCS) representation uses the
// same real, orthonormal spherical-harmonic convention as PFSS.  Unlike PFSS,
// degree zero is mandatory: taking the absolute value of the inner radial
// field produces a nonzero unsigned magnetic flux through the shell.
class FiniteShellScs {
 public:
  static Core::Result<FiniteShellScs> Create(
      double interfaceRadiusM, double radialOuterRadiusM,
      std::vector<HarmonicCoefficient> unsignedInnerBoundary);

  Core::Result<SphericalField> Evaluate(double radiusM, double thetaRad,
                                        double phiRad) const;
  double InterfaceRadiusM() const noexcept { return interfaceRadiusM_; }
  double OuterRadiusM() const noexcept { return outerRadiusM_; }
  const std::vector<HarmonicCoefficient>& Coefficients() const noexcept {
    return coefficients_;
  }

 private:
  double interfaceRadiusM_ = 0.0;
  double outerRadiusM_ = 0.0;
  std::vector<HarmonicCoefficient> coefficients_;
};

// Analytic attenuation of one non-monopole SCS harmonic relative to the
// monopole between the inner and outer shell surfaces (Section 6.2).
Core::Result<double> ScsHarmonicAttenuation(int degree,
                                            double outerToInnerRadius);

struct ScsSpectralDiagnostics {
  double innerNonMonopolePowerT2 = 0.0;
  double outerNonMonopolePowerT2 = 0.0;
  double outerZonalFraction = 0.0;
  std::vector<double> attenuationByDegree;
};
Core::Result<ScsSpectralDiagnostics> EvaluateScsSpectrum(
    const std::vector<HarmonicCoefficient>& coefficients,
    double outerToInnerRadius);

enum class MagneticSector : int { Negative = -1, Positive = 1 };
enum class InterfaceSide { Minus, Plus };

struct OneSidedMagneticField {
  Vec3 valueT;
  MagneticSector sector = MagneticSector::Positive;
  InterfaceSide side = InterfaceSide::Plus;
};

// The sector is categorical.  A query on the mathematical HCS must provide a
// side; this API never interpolates +1 and -1 into an artificial weak field.
Core::Result<OneSidedMagneticField> RestoreSector(
    Vec3 unsignedFieldT, MagneticSector sector, bool onIdealSheet,
    bool sideWasSupplied, InterfaceSide side = InterfaceSide::Plus);

enum class IdealHcsTransport { FieldAlignedNoCrossing, CrossSectorOrDrift };
Core::Status ValidateIdealHcsTransport(IdealHcsTransport requested,
                                       bool qualifiedFiniteSheetProvider);

enum class CouplingMode { ProductionResolvedTransition,
                          SharpInterfaceVerification,
                          NoScsVerification };
Core::Status ValidateCouplingMode(CouplingMode mode, bool productionIntent,
                                  double interfaceRadiusM,
                                  double scsOuterRadiusM,
                                  double transitionWidthM);

struct InterfaceMagneticBalance {
  double normalJumpT = 0.0;
  Vec3 surfaceCurrentAPerM;
  double kinkAngleRad = 0.0;
};
Core::Result<InterfaceMagneticBalance> EvaluateMagneticInterface(
    Vec3 pfssFieldT, Vec3 scsFieldT, Vec3 outwardNormal,
    double normalAbsoluteToleranceT);

// Quintic transition data.  The field is formed by taking the curl of the
// blended vector potential, not by directly blending B.  The cross term is
// essential for both endpoint matching and analytical div(B)=0.
struct TransitionBlend {
  double chi = 0.0;
  double dChiDrPerM = 0.0;
  Vec3 magneticFieldT;
};
Core::Result<TransitionBlend> BlendVectorPotentialFields(
    double radiusM, double innerRadiusM, double outerRadiusM,
    Vec3 radialUnit, Vec3 pfssFieldT, Vec3 scsFieldT,
    Vec3 pfssPotentialTm, Vec3 scsPotentialTm);

struct SignedPotentialQualification {
  bool fluxBalanced = false;
  bool fixedMieGauge = false;
  bool sheetTangentialTraceContinuous = false;
};
Core::Status ValidateSignedPotentialQualification(
    const SignedPotentialQualification& qualification);

// A longitude characteristic carries its derivative A_phi=dPhi/dalpha.
// J_phi=1/A_phi is computed only after the strictly positive fold guard has
// passed.  K and dK/dphi are evaluated in the same rotating-frame convention.
struct LongitudeMap {
  double sourceLongitudeRad = 0.0;
  double mappedLongitudeRad = 0.0;
  double forwardJacobian = 1.0;
  double inverseJacobian = 1.0;
  double foldMargin = 0.0;
};
using WindingRate = std::function<double(double, double)>;
Core::Result<LongitudeMap> IntegrateLongitudeMap(
    double sourceLongitudeRad, double sourceRadiusM, double radiusM,
    int radialSteps, const WindingRate& windingRatePerM,
    const WindingRate& longitudeDerivativePerM, double minimumForwardJacobian);

struct ParkerMappedState {
  Vec3 magneticSphericalT;  // components (B_r,B_theta,B_phi)
  double densityKgM3 = 0.0;
  double fieldAlignedSpeedMPerS = 0.0;
  double mappedMassFluxKgPerSPerSr = 0.0;
};
Core::Result<ParkerMappedState> MapParkerState(
    const LongitudeMap& map, double sourceRadiusM, double radiusM,
    double thetaRad, double sourceRadialFieldT, double sourceDensityKgM3,
    double sourceRadialSpeedMPerS, double radialSpeedMPerS,
    double windingRatePerM);

enum class PlasmaSheetAuthority { BaseDensity, MassPerMagneticFlux };
struct PlasmaSheetNormalization {
  double contrast = 1.0;
  double baseDensityMultiplier = 1.0;
  double basePressureMultiplier = 1.0;
  double massLoadingMultiplier = 1.0;
};
Core::Result<PlasmaSheetNormalization> BuildPlasmaSheetNormalization(
    double neutralLineDistanceRad, double centralContrast,
    double angularWidthRad, PlasmaSheetAuthority authority,
    bool fixedTemperature);

struct TransitionDiagnostics {
  double absoluteCrossingFluxWb = 0.0;
  double netCrossingFluxWb = 0.0;
  double normalTraceRelativeError = 0.0;
  double jumpSupport = 0.0;
  double incidencePlus = 0.0;
  double incidenceMinus = 0.0;
  double antipodalityDefectRad = 0.0;
  bool antipodalityApplicable = false;
};
Core::Result<TransitionDiagnostics> EvaluateTransitionDiagnostics(
    Vec3 plusFieldT, Vec3 minusFieldT, Vec3 sheetNormal,
    double patchAreaM2, double openFluxWb, double minimumFieldT,
    double jumpSupportThreshold, double normalTraceTolerance);

enum class MeasureValidity { Valid, InapplicableZeroDenominator,
                             InapplicablePointMeasure, RejectedClearance };
struct BudgetRatio {
  double numerator = 0.0;
  double denominator = 0.0;
  double ratio = 0.0;
  MeasureValidity validity = MeasureValidity::Valid;
  bool withinBound = false;
};
Core::Result<BudgetRatio> EvaluateBudgetRatio(double excluded,
                                              double counterfactual,
                                              double maximumFraction);
BudgetRatio EvaluatePointClearance(bool intersectsClearance);

enum class CalibrationRole { Construction, Qualification };
struct OpenFluxLifecycle {
  double passAScale = 1.0;
  double passBFluxWb = 0.0;
  bool qualificationPassed = false;
  std::string constructionChecksum;
  std::string qualificationChecksum;
};
Core::Result<OpenFluxLifecycle> BuildOpenFluxLifecycle(
    double unscaledFluxWb, double targetFluxWb, double rebuiltPassBFluxWb,
    double relativeTolerance, const std::string& constructionChecksum,
    const std::string& qualificationChecksum);

} }  // namespace SEP::CoronalCME

#endif  // SEP_CORONAL_CME_SOURCE_SURFACE_COUPLING_H
