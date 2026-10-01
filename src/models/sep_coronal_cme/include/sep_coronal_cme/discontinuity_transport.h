#ifndef SEP_CORONAL_CME_DISCONTINUITY_TRANSPORT_H
#define SEP_CORONAL_CME_DISCONTINUITY_TRANSPORT_H

#include "sep_coronal_cme/mhd_jump_solver.h"
#include "sep_coronal_cme/particle_source.h"
#include "sep_coronal_cme/turbulence_transport.h"
#include "sep_status.h"
#include <cstdint>
#include <memory>
#include <set>
#include <string>
#include <tuple>
#include <vector>

namespace SEP { namespace CoronalCME {

// These identities are intentionally disjoint. A finite exterior HCS never
// qualifies a separatrix or the composite PFSS/SCS transition sheet.
enum class TransportSurfaceKind { IdealHcs, FiniteHcs, Separatrix,
                                  PfssScsTransition, Shock };
struct TransportSurfaceIdentity {
  TransportSurfaceKind kind = TransportSurfaceKind::FiniteHcs;
  std::string stableId;
  std::uint64_t generation = 0;
};
struct MovingPlane {
  TransportSurfaceIdentity identity;
  Vec3 originM;
  Vec3 normal;                    // must already be a unit vector
  double epochS = 0.0;
  double normalSpeedMPerS = 0.0;
};
struct CrossingEvent {
  TransportSurfaceIdentity surface;
  double timeS = 0.0;
  double fraction = 0.0;
  Vec3 positionM;
  InterfaceSide incoming = InterfaceSide::Minus;
  InterfaceSide outgoing = InterfaceSide::Plus;
};
struct SegmentEvent {
  bool crossed = false;
  CrossingEvent event;
};
// Exact signed-distance root on a straight space-time segment. A segment
// starting on the sheet does not create a second event; an ending root belongs
// to the incoming segment. Same-side and tangent segments are not crossings.
Core::Result<SegmentEvent> LocatePlaneCrossing(
    const MovingPlane& plane, Vec3 beginM, Vec3 endM,
    double beginTimeS, double endTimeS);
std::vector<CrossingEvent> OrderCrossingEvents(
    std::vector<CrossingEvent> events, double coincidenceToleranceS);

struct FiniteHcsParameters {
  MovingPlane sheet;
  Vec3 tangent;                  // unit, perpendicular to sheet.normal
  double halfThicknessM = 0.0;  // physical HCS scale, not an interface width
  double fieldMagnitudeT = 0.0;
  double massDensityKgM3 = 0.0;
  double pressurePa = 0.0;
  double outwardWaveEnergyJPerM3 = 0.0;
  double inwardWaveEnergyJPerM3 = 0.0;
};
struct FiniteHcsState {
  MhdPrimitiveState plasma;
  DirectionalWaveState waves;
  double signedDistanceM = 0.0;
  bool signedWaveLabelsValid = false;
};
struct OrbitState {
  Vec3 positionM;
  Vec3 momentumKgMPerS;          // inertial full vector, never a gyro-average
  double timeS = 0.0;
};
struct OrbitControls {
  double maximumGyroAngleRad = 0.05;
  double maximumThicknessFraction = 0.05;
  double eventTimeToleranceS = 1.0e-10;
  // Zero uses the analytic field. Positive spacing samples/normalizes the
  // same force-free orientation on a 1-D mesh for convergence studies.
  double fieldGridSpacingM = 0.0;
  std::uint64_t maximumSubsteps = 1000000;
};
struct HcsOrbitResult {
  OrbitState state;
  std::vector<CrossingEvent> events;
  std::uint64_t substeps = 0;
};

// Exact stationary, force-free planar rotational sheet:
// B=B0[tanh(d/a)t + sech(d/a)(n x t)], rho,p constant, u=E=0.
// B.n=0, div B=0, |B|=B0 and J x B=0. This supplied, qualified
// full-orbit family is exterior only. No unqualified curved-sheet or
// PFSS/SCS-overlap construction is inferred from these parameters.
class FiniteHcsSheet {
 public:
  static Core::Result<std::shared_ptr<const FiniteHcsSheet>> Create(
      const FiniteHcsParameters& parameters);
  Core::Result<FiniteHcsState> Evaluate(Vec3 positionM) const;
  Core::Result<HcsOrbitResult> Advance(
      const OrbitState& initial, double massKg, double chargeC,
      double durationS, const OrbitControls& controls = {}) const;
  const FiniteHcsParameters& Parameters() const noexcept { return parameters_; }
 private:
  FiniteHcsParameters parameters_;
};

struct SheathParameters {
  MovingPlane shock;
  MhdPrimitiveState upstream;
  double gammaAdiabatic = 5.0 / 3.0;
  double thicknessAtEpochM = 0.0;
  double thicknessRateMPerS = 0.0;
  double beginTimeS = 0.0, endTimeS = 0.0;
  double auditAreaM2 = 0.0;
  double conservationTolerance = 1.0e-8;
  // Sector belongs to the background, not the sign of B dot shock normal.
  MagneticSector magneticSector = MagneticSector::Positive;
  double maximumPassiveWavePressureFraction = 0.01;
  double outwardWaveEnergyJPerM3 = 0.0;
  double inwardWaveEnergyJPerM3 = 0.0;
};
struct SheathAudit {
  double divergenceRelative = 0.0;
  double massRelative = 0.0;
  double normalFluxRelative = 0.0;
  double energyRelative = 0.0;
  double waveFluxRelative = 0.0;
  double massInventoryRateKgPerS = 0.0;
  double energyInventoryRateW = 0.0;
  bool passed = false;
};
struct ShockCrossingState {
  FourMomentum inertial;
  FourMomentum localPlasma;
  double pitchAngleCosine = 0.0;
  double shockFrameEnergyJ = 0.0;
  bool applied = false;
};

struct SheathOrbitResult {
  OrbitState state;
  std::vector<CrossingEvent> events;
  std::vector<ShockCrossingState> crossingStates;
  std::uint64_t substeps = 0;
};

// Particle lineage is part of event identity. The same geometric event cannot
// be committed again after AMR subdivision or a repeated callback. A restart
// owner may serialize Keys() and reconstruct the ledger through Restore().
class CrossingLedger {
 public:
  using Key = std::tuple<std::string, std::string, std::uint64_t>;
  bool Contains(const Key& key) const { return keys_.count(key) != 0; }
  const std::set<Key>& Keys() const noexcept { return keys_; }
  Core::Status Restore(const std::set<Key>& keys);
 private:
  friend class PlanarShockSheath;
  std::set<Key> keys_;
};

// A finite moving planar control volume, not a collection of independently
// glued ellipsoid patches. Every downstream point is the SAME solved RH
// state; B.n is continuous at the front, tangential side fluxes cancel, and
// a moving rear boundary has explicitly accounted mass/energy outflow.
// Passive waves transmit their shock-frame normal energy flux with zero
// reflection; both characteristic speeds must be nonsingular. Wave feedback
// on the MHD jump and electrostatic shock potentials are outside this family.
class PlanarShockSheath {
 public:
  static Core::Result<std::shared_ptr<const PlanarShockSheath>> Create(
      const SheathParameters& parameters);
  Core::Result<MhdPrimitiveState> Evaluate(Vec3 positionM, double timeS,
      bool onFront = false, InterfaceSide side = InterfaceSide::Plus) const;
  Core::Result<SheathAudit> Audit(int normalCells, double timeS) const;
  Core::Result<ShockCrossingState> Cross(
      const std::string& particleId, const CrossingEvent& event,
      const FourMomentum& inertial, CrossingLedger* ledger) const;
  // Uniform one-sided ideal-MHD fields are integrated with a full-vector
  // relativistic Boris map; a front event splits the map BEFORE the downstream
  // field is used. Ledger changes publish only when the entire step succeeds.
  Core::Result<SheathOrbitResult> Advance(
      const std::string& particleId, const OrbitState& initial,
      double massKg, double chargeC, double durationS, CrossingLedger* ledger,
      const OrbitControls& controls = {}) const;
  const MhdShockSolution& Jump() const noexcept { return jump_; }
  const DirectionalWaveState& DownstreamWaves() const noexcept { return waves_; }
  const SheathParameters& Parameters() const noexcept { return parameters_; }
 private:
  SheathParameters parameters_;
  MhdShockSolution jump_;
  DirectionalWaveState waves_;
};

Core::Status ValidateDiscontinuityCapabilities(
    bool requestHcs, bool requestSheath,
    const std::shared_ptr<const FiniteHcsSheet>& hcs,
    const std::shared_ptr<const PlanarShockSheath>& sheath,
    bool requestCompositeTransition = false);

} } // namespace SEP::CoronalCME
#endif
