#ifndef SRCSEP_UTIL_SEP_MATERIAL_TUBE_KINEMATICS_H
#define SRCSEP_UTIL_SEP_MATERIAL_TUBE_KINEMATICS_H
#include "sep_common_header_path.h"
#include SRCSEP_SEP_COMMON_HEADER(sep_transport_common.h)
#include <array>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

namespace SEP { namespace FieldLine { namespace MaterialTube {
using Vector=std::array<double,3>;
using Matrix=std::array<double,9>; // Row-major F_ij=partial X_i/partial a_j.

// F and Fdot refer to one fixed material reference, not successive unrelated
// local charts. J is explicitly supplied as an independent declaration and
// checked against det(F). No map, rate, origin or density is inferred.
struct Deformation {
  Matrix gradient{{0,0,0,0,0,0,0,0,0}};
  Matrix ratePerS{{0,0,0,0,0,0,0,0,0}};
  double declaredJacobian=0;
};
struct PointInput {
  Deformation deformation;
  Vector referenceMagneticFieldT{{0,0,0}},referenceAreaVectorM2{{0,0,0}};
  double referenceMassDensityKgM3=0;
};
struct PointState {
  double jacobian=0,jacobianRatePerS=0,massDensityKgM3=0,massDensityRateKgM3PerS=0;
  Matrix cofactor{{0,0,0,0,0,0,0,0,0}},velocityGradientPerS{{0,0,0,0,0,0,0,0,0}};
  Vector magneticFieldT{{0,0,0}},magneticFieldRateTPerS{{0,0,0}},areaVectorM2{{0,0,0}};
  double referenceFluxWb=0,currentFluxWb=0;
};
struct PointResult {Transport::Status status;PointState state;};
// Cauchy induction, source-free mass conservation and Nanson area mapping:
// B=F*Bref/J, rho=rhoref/J, dA=cof(F)*dAref. Pointwise nonaffine maps are
// allowed if the caller supplies their actual compatible F/Fdot; this does
// not establish a nonaffine cell interpolation, quadrature or native adapter.
// Both wide construction and stored binary64 B/area must retain signed flux
// within512*eps64 relative to resolved nonzero referenceflux. For exactzero,
// the absolute criterion is512*eps64*||Bref||*||dAref||, not a flux floor.
// Error-free signed-dot expansions explicitly require nearest rounding and
// normal product/residual exponent headroom; unresolved precision rejects.
// Reported currentflux derives from encoded B/area. Failure publishes no state.
PointResult PushForwardPoint(const PointInput& input);

enum class FaceRole { LowerCap,UpperCap,Lateral };
struct ReferenceFace {
  std::uint64_t id=0;
  FaceRole role=FaceRole::Lateral;
  // These are INTEGRALS over the actual reference surface, which may be
  // curved: A_j=integral n_j dA [m2], M_ij=integral a_i*n_j dA [m3].
  // A may vanish for a complete closed curved side; no fictitious planar
  // area/normal is assigned to it. The moment still contributes to GCL.
  Vector orientedAreaM2{{0,0,0}};
  Matrix firstMomentM3{{0,0,0,0,0,0,0,0,0}};
  // Signed material magnetic flux is supplied from the reference B field.
  // Cauchy/Nanson imply its invariance pointwise, including varying Bref.
  // This integral is not reconstructed as B_mid dot integratedArea.
  double magneticFluxWb=0;
};
struct ReferenceCell {
  std::uint64_t fieldLine=0,segment=0;
  double volumeM3=0,materialMassKg=0;
  std::vector<ReferenceFace> faces;
  std::string identity,provenance;
};
struct AffineSample {
  double timeS=0;
  Deformation deformation;
  Vector translationM{{0,0,0}},translationVelocityMPerS{{0,0,0}};
};
struct TimeQuadraturePoint {double weightS=0;AffineSample sample;};
enum class MapCapability { Unspecified,GlobalAffine,NonAffineUnavailable };
struct QuadratureBudget {
  bool declared=false;
  // Explicit NUMERICAL absolute error allowances, not physical cutoffs.
  // Zero asks for roundoff-only closure. The module measures residuals and
  // never changes a swept face or endpoint volume to make a check pass.
  double volumeM3=0,deformation=0,translationM=0;
  std::string provenance;
};
struct IntervalInput {
  MapCapability capability=MapCapability::Unspecified;
  AffineSample initial,final;
  std::vector<TimeQuadraturePoint> quadrature;
  QuadratureBudget budget;
  std::string operationProvenance;
};
struct FaceGeometry {
  std::uint64_t id=0;FaceRole role=FaceRole::Lateral;
  Vector initialAreaM2{{0,0,0}},finalAreaM2{{0,0,0}};
  Matrix initialFirstMomentM3{{0,0,0,0,0,0,0,0,0}},finalFirstMomentM3{{0,0,0,0,0,0,0,0,0}};
  double outwardGridSweptVolumeM3=0,materialMagneticFluxWb=0;
};
struct GeometryCandidate {
  std::uint64_t fieldLine=0,segment=0;
  double initialTimeS=0,finalTimeS=0,initialVolumeM3=0,finalVolumeM3=0;
  double materialMassKg=0,initialMeanDensityKgM3=0,finalMeanDensityKgM3=0;
  std::vector<FaceGeometry> faces;
  double capSweptVolumeM3=0,lateralSweptVolumeM3=0;
  double gclResidualM3=0,gclAllowedResidualM3=0;
  double maximumDeformationResidual=0,maximumTranslationResidualM=0;
  std::string referenceIdentity,referenceProvenance,operationProvenance,quadratureProvenance;
  bool withinDeclaredQuadratureBudget=false,qualifiedForProduction=false;
};
struct GeometryResult {
  Transport::Status status;
  std::shared_ptr<const GeometryCandidate> candidate;
};
// For X(a,t)=F(t)*a+c(t), face moments integrate U dot dA exactly in SPACE:
// integral U_i dA_i = Fdot_ij*cof(F)_ik*M_jk + cdot_i*cof(F)_ik*A_k.
// TIME quadrature is explicit and checked against independent endpoint J/F/c.
// Closed-reference area/moment identities and all cap+side contributions are
// checked first. True material mass is retained, never rho_midpoint*V guessed.
// Failure yields no candidate; inputs, native arrays and physical owners are
// never mutated. General nonaffine cell maps remain explicitly unsupported.
GeometryResult PrepareAffineInterval(const ReferenceCell& reference,const IntervalInput& interval);
} } }
#endif
