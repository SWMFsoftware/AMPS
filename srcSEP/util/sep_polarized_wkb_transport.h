#ifndef SRCSEP_UTIL_SEP_POLARIZED_WKB_TRANSPORT_H
#define SRCSEP_UTIL_SEP_POLARIZED_WKB_TRANSPORT_H
#include "sep_polarized_wave_state.h"
#include "sep_exact_nonnegative_sum.h"

namespace SEP { namespace Turbulence { namespace Polarized { namespace WKB {

using Vector=std::array<double,3>;
using Tensor=std::array<double,9>; // G_ij=partial U_i/partial x_j, row major.
using FourBins=std::array<std::vector<double>,4>;
using FourTotals=std::array<double,4>;
enum class Side { Lower,Upper };
enum class GeometryMode { Unspecified,MaterialTube,ControlledExplicitAle };

// A declared external source, not a second evolving wave state. Values are
// INITIAL magnetic J/m^3 per finite bin, held fixed thereafter on incoming
// characteristics only. Capture requires revision1; restart must retain this
// object rather than reconstructing it from an evolved interior spectrum.
struct InitialReservoir {
  Grid grid;
  FourBins magneticDensityJPerM3;
  std::uint64_t fieldLine=0;
  Side side=Side::Lower;
  double initialEpochS=0;
  std::string initialSnapshotIdentity,configurationFingerprint,provenance,identity;
  // identity binds all declared values against accidental mutation/codec
  // corruption. It is not proof of historical origin: capture/restart must
  // come from the enclosing trusted immutable initial-reservoir owner.
};
struct ReservoirResult {Transport::Status status;InitialReservoir reservoir;};
ReservoirResult CaptureInitialReservoir(const Snapshot& initial,std::uint64_t fieldLine,Side side);

// Explicit time-integrated face geometry/kinematics. The caller supplies the
// actual quadrature; this module never invents a material/Eulerian tube map.
// Both endpoint states and the known no-reversal flag constrain flux direction.
// Splitting at a characteristic sign change is the caller's responsibility.
struct FaceInterval {
  // End-cap normals point toward increasing field-line label. Material caps
  // can tilt relative to B; area*normal, not A_perp*b, defines geometric and
  // pressure-work fluxes. The Alfvén flux uses the actual b dot normal.
  Vector initialNormal{{0,0,0}},finalNormal{{0,0,0}};
  Vector initialMagneticFieldT{{0,0,0}},finalMagneticFieldT{{0,0,0}};
  Vector initialPlasmaVelocityMPerS{{0,0,0}},finalPlasmaVelocityMPerS{{0,0,0}};
  Vector initialGridVelocityMPerS{{0,0,0}},finalGridVelocityMPerS{{0,0,0}};
  double initialAreaM2=0,finalAreaM2=0,initialMassDensityKgM3=0,finalMassDensityKgM3=0;
  double plasmaSweptVolumeM3=0; // integral A*(U dot normal) dt.
  double gridSweptVolumeM3=0;   // integral A*(v_grid dot normal) dt.
  double alfvenSweptVolumeM3=0;// integral A*v_A*(b dot normal) dt.
  double relativeMassTransportKg=0; // integral rho*A*(U-v_grid) dot normal dt.
  double relativeFlowOverAlfvenAreaTimeM2S=0; // integral A*(U-v_grid)dot(normal)/v_A dt.
  double signedAreaTimeM2S=0; // integral A*(b dot normal) dt.
  bool integralsKnown=false,noCharacteristicReversalKnown=false;
};

// Every flux is integrated over the interval and oriented OUT of the domain.
// Action is TOTAL-wave action [J*s], epsilon_total=2*epsilon_magnetic.
// Pressure-work flux [J] is U*P_wave, P_wave=W_B; it is separate from intrinsic
// transported energy and cannot be counted again as local compression work.
struct BoundaryFlux {
  FourBins outwardMagneticJ,outwardActionJs,outwardPressureWorkJ;
  std::array<bool,2> incoming{{false,false}};
};
struct BoundaryFluxResult {Transport::Status status;BoundaryFlux flux;};
BoundaryFluxResult EvaluateBoundaryFlux(const InitialReservoir& reservoir,
    const FourBins& interiorMagneticDensityJPerM3,const FaceInterval& face,
    double vacuumPermeabilityHPerM);

// Endpoints are prescribed physical inputs. The integrated full velocity
// tensor gives theta*dt and b.G.b*dt, not d(U dot b)/ds. The two scalar
// advection integrals convert material ideal-MHD derivatives into derivatives
// following the explicitly declared grid. Fixed B direction is an explicit
// controlled-interval capability; rotating directions need resolved quadrature.
struct CellInterval {
  Vector initialMagneticFieldT{{0,0,0}},finalMagneticFieldT{{0,0,0}};
  double finalVolumeM3=0,finalMassDensityKgM3=0;
  // Integrated fluid mass is supplied independently. The B/rho samples used
  // for v_A may be midpoint values, so rho_midpoint*V is not an exact mass.
  double initialCellMassKg=0,finalCellMassKg=0;
  Tensor integratedVelocityGradient{{0,0,0,0,0,0,0,0,0}};
  Vector integratedAlfvenSpeedGradient{{0,0,0}}; // integral grad(v_A) dt, dimensionless.
  double relativeFlowDotGradLogB=0; // integral (U-v_grid).grad ln|B| dt.
  double relativeFlowDotGradLogRho=0;
  double lateralGridSweptVolumeM3=0; // complete side-wall geometric contribution.
  Transport::ResonantCollision::DomainContract finalDomain;
  bool inputsKnown=false,fixedFieldDirectionKnown=false,lateralTubeMaterialKnown=false;
};

// Only the spectral range is fixed, not the physical waves. Power leaving it
// is retained in this external account at the crossing energy/action. It is
// neither deleted nor relabeled thermal heating. Zero exterior incoming power
// is an explicit first-case policy. The account is bound to each state revision.
// The enclosing future transaction must retain this account immutably and
// checkpoint it with the wave snapshot; a supplied value is not an independent
// proof of historical fluxes. This prepare-only module creates no second owner.
struct ExternalSpectralAccount {
  FourTotals lowerMagneticJ{{0,0,0,0}},upperMagneticJ{{0,0,0,0}};
  FourTotals lowerActionJs{{0,0,0,0}},upperActionJs{{0,0,0,0}};
  // Exact integer accumulators below are authoritative. The floating displays
  // retain a signed remainder and explicit diagnostic rounding bound; their
  // two-float sum is not an exact representation of an arbitrary stream.
  FourTotals lowerMagneticRoundoffJ{{0,0,0,0}},upperMagneticRoundoffJ{{0,0,0,0}};
  FourTotals lowerActionRoundoffJs{{0,0,0,0}},upperActionRoundoffJs{{0,0,0,0}};
  FourTotals lowerMagneticDiagnosticBoundJ{{0,0,0,0}},upperMagneticDiagnosticBoundJ{{0,0,0,0}};
  FourTotals lowerActionDiagnosticBoundJs{{0,0,0,0}},upperActionDiagnosticBoundJs{{0,0,0,0}};
  std::array<Transport::ExactNonnegativeSum,4> lowerMagneticExact,upperMagneticExact,lowerActionExact,upperActionExact;
  std::string originSnapshotIdentity,lastSnapshotIdentity;
};
struct ExternalAccountResult {Transport::Status status;ExternalSpectralAccount account;};
ExternalAccountResult BeginExternalAccount(const Snapshot& initial);

struct TubeInterval {
  std::uint64_t fieldLine=0;
  std::vector<CellInterval> cells; // all contiguous native segment IDs, ordered from0.
  std::vector<FaceInterval> faces; // cells+1 shared faces, increasing line coordinate.
  InitialReservoir lower,upper;
};
struct Interval {
  double durationS=0;
  std::vector<TubeInterval> tubes; // exactly cover the snapshot, no silent skipped lines.
  EvolutionMetadata target;
  ExternalSpectralAccount externalBefore;
  GeometryMode geometryMode=GeometryMode::Unspecified;
  bool zeroExteriorSpectralIncoming=false;
  bool reflectionEnabled=false,cascadeEnabled=false;
};

struct ComponentLedger {
  // Ordinary floating diagnostic sums over cells/bins. They are numerical
  // balance/convergence evidence, not exact historical accounts. Each accepted
  // cell's escaped increment reaches ExternalSpectralAccount integer authority
  // before these sums, so loss of a low diagnostic bit cannot delete a wave.
  double magneticBeforeJ=0,magneticAfterJ=0;
  double spatialIncomingMagneticJ=0,spatialOutgoingMagneticJ=0;
  double lowerSpectralOutgoingMagneticJ=0,upperSpectralOutgoingMagneticJ=0;
  double magneticCompressionWorkJ=0; // signed work ON magnetic waves.
  double totalBackgroundCompressionReceiptJ=0; // -2*magneticCompressionWorkJ.
  double pressureWorkBoundaryOutJ=0;
  double lateralPressureWorkOutJ=0; // Included in pressureWorkBoundaryOutJ.
  // Pressure is frozen at the donor stage on caps and material side walls.
  // This first-order pressure/splitting quadrature needs its own convergence
  // before the receipt is interpreted as -integral U.grad(P) dt.
  double totalBackgroundMechanicalReceiptJ=0; // compression minus ALL-face U*P.
  double actionBeforeJs=0,actionAfterJs=0,spatialIncomingActionJs=0,spatialOutgoingActionJs=0;
  double lowerSpectralOutgoingActionJs=0,upperSpectralOutgoingActionJs=0;
  double energyRoundoffResidualJ=0;
  // Finite-volume Ebin/Delta k reconstruction does not preserve its 1/k action
  // moment exactly under remapping. This is exposed for mesh/time convergence,
  // never assigned to physical work, cascade, heating or an exterior source.
  double actionDiscretizationResidualJs=0;
};
struct StageResult {
  Transport::Status status;
  EvolutionCandidateResult candidate;
  std::array<ComponentLedger,4> ledger;
  ExternalSpectralAccount externalAfter;
  double maximumSpatialDonorFraction=0;
  bool qualifiedForProduction=false;
};
struct SpectralCellResult {
  Transport::Status status;
  FourBins magneticEnergyAfterJ;
  std::array<ComponentLedger,4> ledger;
};
// A separately testable source/ray stage. All three integrated coefficients
// are dimensionless; gradAlongB is integral d_s(v_A) dt. Frozen/common-profile
// histories are exact. Other histories require resolved operator substeps.
SpectralCellResult PrepareSpectralWork(const Grid& grid,const FourBins& energyBeforeJ,
    double alfvenSpeedBeforeMPerS,double alfvenSpeedAfterMPerS,
    double integratedDivergence,double integratedParallelStrain,double integratedGradAlfvenAlongB);

// Spatial transport uses conservative donor-cell fluxes and an explicit CFL
// rejection; geometry obeys the supplied discrete GCL and mass continuity.
// The source/ray substep pushes the actual Ebin/Delta k reconstruction along
// k=k0*exp(-integral(S+sigma*d_s v_A)dt), with magnetic work exp(-integral theta dt/2).
// Frozen/common-profile coefficients are exact in that substep; spatial/source
// splitting and projection require independent convergence. All edits are
// private candidates, and no native owner, particle or input is mutated.
StageResult PrepareInterval(const Snapshot& before,const Interval& interval);

} } } }
#endif
