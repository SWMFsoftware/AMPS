#ifndef SEP_CORONA_SWCME_SHEATH_MODEL_H
#define SEP_CORONA_SWCME_SHEATH_MODEL_H

#include "sep_corona_swcme/surface_shock.h"

#include <array>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

namespace SEP { namespace CoronaSwcme {

// Material coordinates are immutable: two surface coordinates identify the
// shock element and crossingTimeS identifies its unique admission event.
// They are background-fluid labels, not particle/source identities.
struct SheathMaterialLabel {
  double polarRad = 0.0;
  double azimuthRad = 0.0;
  double crossingTimeS = 0.0;
  std::uint64_t shockPatchLineage = 0;
};

struct SheathMappedState {
  SheathMaterialLabel label;
  CoronalCME::Vec3 positionM;
  CoronalCME::MhdPrimitiveState primitive;
  // Columns of F_rel=A(t) A0^-1 in the inertial Cartesian basis.  A0 is the
  // centre basis of the finite shock-born cohort prism documented below; it
  // is not an identity assumed in (theta,phi,tau) labels.
  std::array<CoronalCME::Vec3,3> deformationColumns;
  // d x/d(theta,phi,tau), with units m/rad, m/rad and m/s. These columns
  // make the spatial inverse an actual material-coordinate solve rather than
  // a nearest-cell assignment.
  std::array<CoronalCME::Vec3,3> labelDerivativeColumns;
  std::array<CoronalCME::Vec3,3> referenceDerivativeColumns;
  // Physical volume metrics |det A| in m^3 rad^-2 s^-1.  Multiplication by
  // dtheta*dphi*dtau gives the curved cohort volume; no unweighted radial
  // column or current-area surrogate is used.
  double referenceVolumeDensityM3PerRad2S = 0.0;
  double currentVolumeDensityM3PerRad2S = 0.0;
  double jacobian = 0.0;
  double upstreamRelativeNormalSpeedMPerS = 0.0;
  double downstreamRelativeNormalSpeedMPerS = 0.0;
  double birthAreaDensityM2PerRad2 = 0.0;
  double shockResidual = 0.0;
  std::string eventIdentity;
};

enum class SheathCellDisposition {
  Retained,
  SolarExit,
  OuterExit
};

// Sheath-side value on the material rear boundary.  The selected zero-volume
// startup makes this the globally oldest cohort (tau=start), not a fixed
// fraction of the shock and not a surface re-created at local fast activation.
struct SheathContactState {
  SheathMaterialLabel label;
  CoronalCME::Vec3 positionM;
  CoronalCME::Vec3 outwardNormal;
  CoronalCME::Vec3 velocityMPerS;
  CoronalCME::MhdPrimitiveState primitive;
  CoronalCME::Vec3 tractionPa;
  double normalSpeedMPerS = 0.0;
  double relativeMassFluxKgM2S = 0.0;
  double normalMagneticFieldT = 0.0;
  double areaDensityM2PerRad2 = 0.0;
  std::string eventIdentity;
};

struct SheathAdmissionCell {
  SheathMaterialLabel label;
  double intervalS = 0.0;
  // Retained cells use the inventory epoch. Exited cells use their first
  // event-located boundary time and are never extrapolated until a later fold.
  double evaluationTimeS = 0.0;
  double birthAreaM2 = 0.0;
  double admittedMassKg = 0.0;
  double currentVolumeM3 = 0.0; // volume at evaluationTimeS
  SheathCellDisposition disposition = SheathCellDisposition::Retained;
  SheathMappedState state;
};

struct SheathInventory {
  double epochS = 0.0;
  std::uint64_t backgroundGeneration = 0;
  std::string eventIdentity;
  std::vector<SheathAdmissionCell> cells;
  double admittedMassKg = 0.0;
  double initialMassKg = 0.0;
  double retainedMassKg = 0.0;
  double solarExitMassKg = 0.0;
  double outerExitMassKg = 0.0;
  double fastAreaTimeM2S = 0.0;
  double subfastAreaTimeM2S = 0.0;
  double maximumAdmissionMassResidual = 0.0;
  double maximumInventoryMassResidual = 0.0;
  double maximumInventoryMassAbsoluteResidualKg = 0.0;
  double minimumJacobian = 0.0;
  bool zeroVolumeStartup = false;
};

// Level A+ diagnostic only: this class is the superseded experimental BG3D-4
// shock-fed map.  Production BG3D-4 is specified by
// docs/BG3D4_PISTON_CLOSURE.md as an ejecta/contact-driven per-ray Lagrangian
// finite-volume model.  Do not extend this class into that solver or publish
// its oldest cohort as the production contact: the two closures have different
// unknowns and authorities.  In particular, this class treats the prescribed
// front and its current RH state as inputs, whereas the replacement prescribes
// one ejecta contact and obtains compression/shock motion from the plasma
// evolution.
//
// For a parcel born at (q,tau), this diagnostic defines
// D(q,t)=U2(q,t)-d_t X_front(q,t), age=t-tau and
// L(age)=kappa*age+(1-kappa)*T*(1-exp(-age/T)).  The current map is
//   x(q,tau,t)=X_front(q,t)+L(age) D(q,t).
// L(0)=0 and L'(0)=1 give x=X_front and U=U2 at admission.  Kappa and T are
// empirical closure inputs; conservation identities do not make this velocity
// law a momentum/energy solution.  In particular it remains unqualified until
// a compatible finite initial state and one cross-stage contact authority are
// supplied.  A finite cell centred on tau owns the
// regular reference prism X0(q,s)=X_front(q,tau)-s deficit(q,tau), where s is
// the within-cohort time coordinate.  At s=0 its basis is
// [d_theta X_front,d_phi X_front,-deficit], whose determinant is the curved
// shock area density times w2.  A(t), A0, F_rel=A A0^-1 and J_rel are retained
// explicitly; rho, p and B use rho2/J_rel, p2*J_rel^-gamma and
// F_rel*B2/J_rel.  These Cauchy/adiabatic identities verify mass, ideal
// induction and entropy along the chosen deformation; they do not determine
// the deformation from momentum or energy.  The measured non-refining force
// residual is therefore model discrepancy, not evidence that this diagnostic
// is a fluid solution.
//
// The use of D(q,t) at query time also means an old parcel depends on a valid
// *current* fast-shock RH state.  If a patch becomes sub-fast, a future
// candidate is rejected and committed cells remain readable; this class does
// not provide the replacement model's supported linear compression region.
class ShockFedSheathModel final {
 public:
  static Core::Result<std::shared_ptr<ShockFedSheathModel>> Create(
      std::shared_ptr<const EventConfiguration> event,
      std::shared_ptr<const AmbientModel> ambient);

  Core::Result<SheathMappedState> Evaluate(
      const SheathMaterialLabel& label,double epochS) const;

  Core::Result<SheathContactState> EvaluateContact(
      double polarRad,double azimuthRad,double epochS) const;

  // Requires a committed inventory at the same epoch. Retained cell centers
  // seed a damped Newton solve; toleranceM is a physical Cartesian distance.
  Core::Result<SheathMappedState> EvaluateAtPosition(
      CoronalCME::Vec3 positionM,double epochS,double toleranceM,
      int maximumIterations=20) const;

  // Read a stored material cell without re-running shock admission at a later
  // candidate time.  This is deliberately limited to the committed epoch:
  // unsupported future evolution must reject transactionally, not fabricate a
  // sub-fast compression state or make already committed inventory disappear.
  Core::Result<SheathMappedState> QueryCommitted(
      const SheathMaterialLabel& label) const;

  Core::Result<std::shared_ptr<const SheathInventory>> PrepareInventory(
      double epochS,std::uint64_t backgroundGeneration,int polarCells,
      int azimuthCells,int admissionTimeCells);

  std::shared_ptr<const SheathInventory> Current() const { return current_; }

 private:
  std::shared_ptr<const EventConfiguration> event_;
  std::shared_ptr<const AmbientModel> ambient_;
  std::shared_ptr<const SheathInventory> current_;
};

const char* Name(SheathCellDisposition) noexcept;

} } // namespace SEP::CoronaSwcme

#endif
