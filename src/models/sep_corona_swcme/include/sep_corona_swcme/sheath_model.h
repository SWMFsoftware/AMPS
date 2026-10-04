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
  // Columns of F=d x(current)/d x(birth), in the inertial Cartesian basis.
  std::array<CoronalCME::Vec3,3> deformationColumns;
  // d x/d(theta,phi,tau), with units m/rad, m/rad and m/s. These columns
  // make the spatial inverse an actual material-coordinate solve rather than
  // a nearest-cell assignment.
  std::array<CoronalCME::Vec3,3> labelDerivativeColumns;
  double jacobian = 0.0;
  double upstreamRelativeNormalSpeedMPerS = 0.0;
  double downstreamRelativeNormalSpeedMPerS = 0.0;
  double birthAreaDensityM2PerRad2 = 0.0;
  double shockResidual = 0.0;
  std::string eventIdentity;
};

enum class SheathCellDisposition {
  Retained,
  ContactExit,
  SolarExit,
  OuterExit
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
  double retainedMassKg = 0.0;
  double contactExitMassKg = 0.0;
  double solarExitMassKg = 0.0;
  double outerExitMassKg = 0.0;
  double fastAreaTimeM2S = 0.0;
  double subfastAreaTimeM2S = 0.0;
  double maximumAdmissionMassResidual = 0.0;
  double maximumInventoryMassResidual = 0.0;
  double minimumJacobian = 0.0;
};

// BG3D-4 analytical shock-fed sheath.  For a parcel born at (q,tau),
//   x(q,tau,t)=X_front(q,t)+(t-tau)[U2(q,tau)-d_t X_front(q,tau)].
// Thus x=X_front and U=U2 at admission.  The label-space derivative supplies
// F and J relative to the finite crossing-volume basis; rho, p and B then use
// rho2/J, p2*J^-gamma and F*B2/J.  This is a prescribed map, not an MHD solve.
class ShockFedSheathModel final {
 public:
  static Core::Result<std::shared_ptr<ShockFedSheathModel>> Create(
      std::shared_ptr<const EventConfiguration> event,
      std::shared_ptr<const AmbientModel> ambient);

  Core::Result<SheathMappedState> Evaluate(
      const SheathMaterialLabel& label,double epochS) const;

  // Requires a committed inventory at the same epoch. Retained cell centers
  // seed a damped Newton solve; toleranceM is a physical Cartesian distance.
  Core::Result<SheathMappedState> EvaluateAtPosition(
      CoronalCME::Vec3 positionM,double epochS,double toleranceM,
      int maximumIterations=20) const;

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
