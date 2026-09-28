#ifndef SEP_CORONAL_CME_MODEL_CONFIGURATION_H
#define SEP_CORONAL_CME_MODEL_CONFIGURATION_H

#include "sep_status.h"

#include <cstdint>
#include <map>
#include <string>
#include <vector>

namespace SEP {
namespace CoronalCME {

// Schema-5 enums are explicit and closed.  Values are never inferred from a
// neighboring selector, because that would make a typo silently change the
// governing equation or its reference frame.
enum class RunIntent { ProductionShockInjection, AnalyticVerification };
enum class TransportModel {
  BallisticVerification,
  Parker,
  FocusedPitchAngleDiffusion,
  FocusedDiscreteScattering
};
enum class TransportFrame { Inertial, RigidCorotating };
enum class SolarRotationModel { Rigid, LatitudeDependentVerification };
enum class WindModel { FluxTubePolytropic, EmpiricalKinematic };
enum class WindEnergyClosure { Isothermal, Polytropic, EmpiricalProfile };
enum class ClosedFieldModel { IsothermalHydrostatic, PolytropicHydrostatic };
enum class InterfaceRepresentation { SharpOneSided, FiniteWidthVolume };
enum class InterfacePolicy {
  DiagnosticKinematic,
  BoundedApproximation,
  StationaryTangentialDiscontinuity
};
enum class DataUseRole { None, Construction, Qualification, WithheldValidation };

// Only Stage-0--2 physical authorities are materialized here.  The complete
// raw assignment map is retained separately, allowing the schema to freeze
// later-stage selectors without falsely constructing unavailable providers.
struct ModelConfiguration {
  int schemaVersion = 0;
  RunIntent intent = RunIntent::AnalyticVerification;
  TransportModel transport = TransportModel::BallisticVerification;
  TransportFrame transportFrame = TransportFrame::Inertial;
  double startTimeS = 0.0;
  double endTimeS = 0.0;
  std::uint64_t campaignSeed = 0;

  double solarRadiusM = 0.0;
  double outerRadiusM = 0.0;
  double pfssOuterBoundaryRadiusM = 0.0;
  double currentSheetInterfaceRadiusM = 0.0;
  double currentSheetOuterRadiusM = 0.0;

  SolarRotationModel rotationModel = SolarRotationModel::Rigid;
  double siderealRotationRateRadPerS = 0.0;
  WindModel windModel = WindModel::FluxTubePolytropic;
  WindEnergyClosure windEnergyClosure = WindEnergyClosure::Polytropic;
  ClosedFieldModel closedFieldModel = ClosedFieldModel::IsothermalHydrostatic;
  InterfaceRepresentation interfaceRepresentation =
      InterfaceRepresentation::SharpOneSided;
  InterfacePolicy interfacePolicy = InterfacePolicy::DiagnosticKinematic;

  // These three indices have independent physical roles and are deliberately
  // stored separately even when a verification deck gives them equal values.
  double gammaWind = 0.0;
  double gammaAdiabatic = 0.0;
  double gammaClosed = 0.0;
  bool includeElectronMassInDensity = false;

  // Canonical dotted-key assignments are normalized before hashing.  Paths
  // are provenance only; checksums and physical selectors define identity.
  std::map<std::string, std::string> assignments;
  std::string physicsFingerprint;
};

enum class ParseDisposition { ResolvedSchema5, LegacyPassThrough };

struct VersionedConfiguration {
  ParseDisposition disposition = ParseDisposition::ResolvedSchema5;
  ModelConfiguration schema5;
  int legacySchemaVersion = 0;
  std::string legacyBytes;
  std::string fingerprint;
};

// Returns a stable lower-case spelling used in diagnostics and manifests.
const char* Name(RunIntent value) noexcept;
const char* Name(TransportModel value) noexcept;
const char* Name(TransportFrame value) noexcept;
const char* Name(SolarRotationModel value) noexcept;
const char* Name(WindModel value) noexcept;
const char* Name(WindEnergyClosure value) noexcept;
const char* Name(ClosedFieldModel value) noexcept;
const char* Name(InterfaceRepresentation value) noexcept;
const char* Name(InterfacePolicy value) noexcept;

}  // namespace CoronalCME
}  // namespace SEP

#endif  // SEP_CORONAL_CME_MODEL_CONFIGURATION_H
