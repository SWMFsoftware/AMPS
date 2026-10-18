// Strict file format and evaluation for the optional H2O activity residual.
//
// The spherical-harmonic mathematics lives in ActivityEnhancementBasis.* and
// is specified by HARMONIC_BASIS.md.  This component adds the versioned file
// contract needed by the forward solver.  Keeping parsing separate makes it
// possible to reject an incompatible fitted product before any propagation.
#ifndef ACTIVITY_ENHANCEMENT_H
#define ACTIVITY_ENHANCEMENT_H

#include <string>
#include <vector>

#include "Vec3.h"

namespace ActivityEnhancement {

constexpr const char* ModelId="real_spherical_harmonics_v1";
constexpr const char* FileFormatId="liouville_cops.activity_enhancement";
constexpr int FileFormatVersion=1;
constexpr const char* AngleConvention="theta=acos(z/r);phi=atan2(y,x)";
constexpr const char* Ordering="increasing_l_then_m_minus_l_to_plus_l";
constexpr const char* ClippingRule="max(0,legacy_linear+enhancement)";
constexpr const char* ScaleInterpretation="dimensionless_multiplier_of_normalized_basis";

// All values below are either validated file contents or explicit runtime
// state.  Coefficients follow the normative increasing-l/increasing-m order.
// MaximumDegree=4 is the zero-length no-enhancement control; otherwise the
// represented range is exactly degrees 5..MaximumDegree.
struct Configuration {
  std::string Model="none";
  std::string FilePath;
  std::string Format=FileFormatId;
  int FormatVersion=FileFormatVersion;
  std::string BasisId;
  std::string Species="H2O";
  std::string BodyFrame;
  std::string AngleDefinition=AngleConvention;
  std::string CoefficientOrdering=Ordering;
  std::string Clipping=ClippingRule;
  int MinimumDegree=5;
  int MaximumDegree=4;
  double Scale=1.0;
  std::string ScaleMeaning=ScaleInterpretation;
  std::string LegacyActivityIdentity;
  std::string Regularization;
  double Lambda=0.0;
  std::vector<double> Coefficients;

  bool Enabled() const {return Model==ModelId;}
  bool HasResidualModes() const {return Enabled()&&MaximumDegree>=5;}
};

// Parse and validate every field and coefficient in a version-1 file.
// Throws std::runtime_error on missing, duplicate, unknown, nonfinite, or
// incompatible content.  The caller separately verifies body-frame and
// legacy-map identity against the active forward model.
Configuration ReadFile(const std::string& filename);

//Reject a coefficient product fitted for a different mesh/SPICE frame or a
//different legacy map.  Exact string equality is intentional: both identities
//are provenance tokens, not approximate physical values.
void ValidateCompatibility(const Configuration& configuration,
                           const std::string& bodyFrame,
                           const std::string& legacyActivityIdentity);

// Evaluate DeltaA = Scale * sum b_lm Y_lm at a mesh/body-frame direction.
// A degree-4 control returns exactly zero.  Invalid Cartesian vectors are
// rejected by the shared basis evaluator.
double EvaluateResidual(const Configuration& configuration,const Vec3& position);

// Apply the one-and-only clipping operation required by the amendment.
// Optional outputs expose the two linear components for diagnostics/tests.
double Apply(const Configuration& configuration,const Vec3& position,
             double legacyLinear,double* residual=nullptr,
             double* totalLinear=nullptr);

} // namespace ActivityEnhancement

#endif
