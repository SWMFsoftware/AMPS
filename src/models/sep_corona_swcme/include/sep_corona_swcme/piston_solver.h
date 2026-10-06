#ifndef SEP_CORONA_SWCME_PISTON_SOLVER_H
#define SEP_CORONA_SWCME_PISTON_SOLVER_H

#include "sep_corona_swcme/cme_event.h"

#include <array>
#include <functional>
#include <memory>
#include <vector>

namespace SEP { namespace CoronaSwcme {

enum class PistonTubeGeometry { Planar, RadialSpherical };

// Numerical controls for the first production-solver increment.  The planar
// geometry is deliberately explicit: it qualifies the moving material
// boundary, conservative pressure/viscous work and shock capture before the
// same discretization is extended to A=DeltaOmega*r^2 and magnetic stresses.
// All quantities are SI in production; tests may choose a self-consistent
// nondimensional system because the Euler equations are scale invariant.
struct PlanarPistonInput {
  PistonTubeGeometry geometry = PistonTubeGeometry::Planar;
  double startS = 0.0;
  double endS = 0.0;
  double leftPositionM = 0.0;
  double columnLengthM = 0.0;
  double areaM2 = 1.0;
  // Physical solid angle for RadialSpherical.  The exact shell volume is
  // DeltaOmega*(r_R^3-r_L^3)/3 and face area is DeltaOmega*r^2.
  double solidAngleSr = 0.0;
  double initialDensityKgM3 = 0.0;
  double initialPressurePa = 0.0;
  double initialVelocityMPerS = 0.0;
  // Signed transverse field.  Planar ideal induction preserves B_t/rho;
  // density compression therefore supplies magnetic pressure without an
  // independently evolved, potentially inconsistent field variable.
  double initialTransverseMagneticFieldT = 0.0;
  // Second component in a fixed right-handed tangent basis attached to the
  // ray.  Both components obey the same radial frozen-flux equation but must
  // remain distinct: replacing them by a magnitude would lose the signed IMF
  // direction needed by later 3-D assembly and shock-obliquity diagnostics.
  double initialTransverseMagneticField2T = 0.0;
  // Ray-parallel IMF.  In spherical geometry each material cell preserves
  // B_r r^2.  It contributes magnetic energy and moving-boundary Maxwell
  // work but no bulk radial Lorentz force in the selected 1-D reduction.
  double initialRadialMagneticFieldT = 0.0;
  double magneticPermeabilityNPerA2 = 1.25663706212e-6;
  double gammaAdiabatic = 0.0;
  int cells = 0;
  double cfl = 0.0;
  double quadraticViscosity = 0.0;
  double linearViscosity = 0.0;
  // The linear VNR term is disabled for infinitesimal compression when
  // -DeltaU/c is below this dimensionless sensor.  This prevents a nominal
  // shock stabilizer from becoming leading-order damping of the linear
  // sub-fast acceptance wave.  The quadratic term remains compression-only.
  double linearViscosityActivation = 0.0;
  double shockThreshold = 0.0;
};

using PistonHistory = std::function<Core::Result<CoronalCME::KinematicValue>(
    double timeS)>;

// Eulerian volume sources evaluated at the current material position.  The
// radial acceleration is the *net* permitted body acceleration
// -GM_sun/r^2+f_amb.  Keeping that sum in the integrator avoids cancellation
// of two large numbers, while the producer retains and reports its separate
// gravity and ambient-maintaining terms.  Heating is per physical volume.
// The two magnetic entries are d[B_t/(rho r)]/dt in spherical geometry and
// d(B_t/rho)/dt in planar geometry, in the same tangent basis as the state.
struct PistonVolumeSource {
  double radialAccelerationMPerS2 = 0.0;
  double heatingWPerM3 = 0.0;
  std::array<double,2> transverseInvariantRate{};
};

using PistonSourceHistory = std::function<Core::Result<PistonVolumeSource>(
    double positionM,double timeS)>;

struct PlanarCellState {
  double centerM = 0.0;
  double widthM = 0.0;
  double massKg = 0.0;
  double densityKgM3 = 0.0;
  double pressurePa = 0.0;
  double specificInternalEnergyJPerKg = 0.0;
  double artificialPressurePa = 0.0;
  double transverseMagneticFieldT = 0.0;
  double transverseMagneticField2T = 0.0;
  double radialMagneticFieldT = 0.0;
  double totalPressurePa = 0.0;
  double velocityMPerS = 0.0;
};

struct PlanarShockState {
  bool present = false;
  bool statesAvailable = false;
  // Detector evidence is retained because an outermost Q/P excursion can be
  // a weak leading precursor rather than the settled piston shock.  A change
  // in zone count or selected-zone strength under refinement is therefore a
  // shock-identity failure, not merely a radius truncation error.
  int shockZoneCount = 0;
  double maximumArtificialPressureRatio = 0.0;
  double selectedZoneMaximumArtificialPressureRatio = 0.0;
  double radiusM = 0.0;
  double speedMPerS = 0.0;
  double compressionRatio = 0.0;
  double upstreamDensityKgM3 = 0.0;
  double upstreamPressurePa = 0.0;
  double upstreamVelocityMPerS = 0.0;
  double upstreamTransverseMagneticFieldT = 0.0;
  double upstreamTransverseMagneticField2T = 0.0;
  double upstreamRadialMagneticFieldT = 0.0;
  double downstreamDensityKgM3 = 0.0;
  double downstreamPressurePa = 0.0;
  double downstreamVelocityMPerS = 0.0;
  double downstreamTransverseMagneticFieldT = 0.0;
  double downstreamTransverseMagneticField2T = 0.0;
  double downstreamRadialMagneticFieldT = 0.0;
  int firstShockCell = -1;
  int lastShockCell = -1;
};

struct PlanarEnergyLedger {
  double initialEnergyJ = 0.0;
  double currentEnergyJ = 0.0;
  double pistonWorkJ = 0.0;
  double outerBoundaryWorkJ = 0.0;
  double bodyForceWorkJ = 0.0;
  double volumeHeatingJ = 0.0;
  double magneticSourceWorkJ = 0.0;
  double radialMagneticBoundaryWorkJ = 0.0;
  double appendedEnergyJ = 0.0;
  double residualJ = 0.0;
};

struct PistonInitialState {
  std::vector<double> nodeVelocityMPerS;
  std::vector<double> cellDensityKgM3;
  std::vector<double> cellPressurePa;
  std::vector<double> cellTransverseMagneticFieldT;
  std::vector<double> cellTransverseMagneticField2T;
  std::vector<double> cellRadialMagneticFieldT;
};

struct PistonAppendState {
  // Physical HCI positions include the old outer node as element zero.  The
  // remaining nodes and all cell primitives describe newly admitted ambient
  // material ordered outward.  No interpolation/remap of committed cells is
  // performed.
  std::vector<double> nodePositionM;
  PistonInitialState state;
  PistonHistory outerBoundary;
};

class PlanarPistonSolver final {
 public:
  static Core::Result<std::unique_ptr<PlanarPistonSolver>> Create(
      PlanarPistonInput input,PistonHistory piston,
      PistonSourceHistory source=PistonSourceHistory{},
      PistonHistory outerBoundary=PistonHistory{});
  // Manufactured references and restarts may supply a complete compatible
  // material state.  Production ambient-only startup continues to call
  // Create(); this overload must not be used to hide a missing CME prehistory.
  static Core::Result<std::unique_ptr<PlanarPistonSolver>> CreateInitialized(
      PlanarPistonInput input,PistonHistory piston,PistonInitialState state,
      PistonSourceHistory source=PistonSourceHistory{},
      PistonHistory outerBoundary=PistonHistory{});

  // Advance is transactional: integration occurs on a private copy and is
  // committed only if every cell remains finite with positive volume,
  // density, pressure and internal energy.  A rejected step cannot corrupt a
  // state that future material queries depend on.
  Core::Status AdvanceTo(double targetTimeS);

  // Extend the material domain transactionally at the current epoch.  The
  // shared first node must match the committed outer material surface in
  // position and velocity.  The exact energy introduced by new cell masses
  // and endpoint half masses is recorded as append/boundary transport; failed
  // input leaves the committed inventory and outer history unchanged.
  Core::Status AppendAmbient(PistonAppendState appended);

  // Deep private candidate for multi-ray epoch transactions.  std::function
  // authorities and immutable material arrays are copied; callers commit by
  // replacing the owning pointer only after every ray candidate succeeds.
  std::unique_ptr<PlanarPistonSolver> Clone() const;

  Core::Result<std::vector<PlanarCellState>> Cells() const;
  Core::Result<PlanarShockState> DetectShock() const;
  Core::Result<PlanarEnergyLedger> EnergyLedger() const;

  double TimeS() const noexcept { return timeS_; }
  const PlanarPistonInput& Input() const noexcept { return input_; }
  std::vector<double> NodePositionsM() const;
  const std::vector<double>& NodeVelocitiesMPerS() const noexcept { return u_; }
  std::size_t CellCount() const noexcept { return cellMass_.size(); }

 private:
  PlanarPistonInput input_;
  PistonHistory piston_,outerBoundary_;
  PistonSourceHistory source_;
  double timeS_ = 0.0;
  // x_ is stored in the frame translating with the initial ambient velocity.
  // Widths are therefore differences of stationary coordinates in the exact
  // well-balanced solution, avoiding loss of small cell widths when a large
  // common heliocentric translation is repeatedly accumulated.
  std::vector<double> cellMass_,x_,u_,e_;
  std::array<std::vector<double>,2> transverseInvariant_;
  std::vector<double> radialFluxInvariant_;
  double initialEnergyJ_ = 0.0;
  double pistonWorkJ_ = 0.0;
  double outerBoundaryWorkJ_ = 0.0;
  double bodyForceWorkJ_ = 0.0;
  double volumeHeatingJ_ = 0.0;
  double magneticSourceWorkJ_ = 0.0;
  double radialMagneticBoundaryWorkJ_ = 0.0;
  double appendedEnergyJ_ = 0.0;
};

} } // namespace SEP::CoronaSwcme

#endif
