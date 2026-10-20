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
#include <cstddef>
#include <limits>

#include "Vec3.h"

class cMesh;

namespace ActivityEnhancement {

constexpr const char* ModelId="real_spherical_harmonics_v1";
constexpr const char* FileFormatId="liouville_cops.activity_enhancement";
constexpr const char* FacetModelId="per_facet_residual_v1";
constexpr const char* FacetFileFormatId="liouville_cops.activity_enhancement_facets";
constexpr const char* CompositeBasisId="real_spherical_harmonics_plus_declared_facet_templates_v1";
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
  //Explicit residuals are used only to replay a Python fit containing one or
  //more localized mesh-bound templates.  The reader validates every saved
  //centroid and area against the active mesh before this vector is accepted.
  std::vector<double> FacetResiduals;

  bool Enabled() const {return Model==ModelId||Model==FacetModelId;}
  bool IsHarmonic() const {return Model==ModelId;}
  bool IsFacetResidual() const {return Model==FacetModelId;}
  bool HasResidualModes() const {
    return (IsHarmonic()&&MaximumDegree>=5)||
           (IsFacetResidual()&&!FacetResiduals.empty());
  }
};

// Parse and validate every field and coefficient in a version-1 file.
// Throws std::runtime_error on missing, duplicate, unknown, nonfinite, or
// incompatible content.  The caller separately verifies body-frame and
// legacy-map identity against the active forward model.
Configuration ReadFile(const std::string& filename);

//Read an explicit per-facet residual and prove that face count, order,
//centroids and areas match the active mesh exactly.  This is intentionally a
//replay format, not a portable spatial interpolation format.
Configuration ReadFacetFile(const std::string& filename,const ::cMesh& mesh);

//Reject a coefficient product fitted for a different mesh/SPICE frame or a
//different legacy map.  Exact string equality is intentional: both identities
//are provenance tokens, not approximate physical values.
void ValidateCompatibility(const Configuration& configuration,
                           const std::string& bodyFrame,
                           const std::string& legacyActivityIdentity);

// Evaluate DeltaA = Scale * sum b_lm Y_lm at a mesh/body-frame direction.
// A degree-4 control returns exactly zero.  Invalid Cartesian vectors are
// rejected by the shared basis evaluator.
double EvaluateResidual(const Configuration& configuration,const Vec3& position,
                        std::size_t faceIndex=std::numeric_limits<std::size_t>::max());

// Apply the one-and-only clipping operation required by the amendment.
// Optional outputs expose the two linear components for diagnostics/tests.
double Apply(const Configuration& configuration,const Vec3& position,
             double legacyLinear,double* residual=nullptr,
             double* totalLinear=nullptr,
             std::size_t faceIndex=std::numeric_limits<std::size_t>::max());

} // namespace ActivityEnhancement

#endif
