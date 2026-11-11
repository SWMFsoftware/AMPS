#ifndef SEP_CORONA_SWCME_CROSS_RAY_TRANSPORT_H
#define SEP_CORONA_SWCME_CROSS_RAY_TRANSPORT_H

#include "sep_coronal_cme/vector_math.h"
#include "sep_status.h"

namespace SEP { namespace CoronaSwcme { namespace CrossRay {

// Bounded conservative transport REFERENCE kernels, independent of PIC/MPI.
// Production uses the selected analytical/semi-empirical material closure,
// without an MHD runtime; these fluxes supply independent conservation checks.
// They are not yet an evolved cross-ray mesh, constrained transport, piston
// coupling, or production G1 provider. See CONSERVATIVE_CROSS_RAY_TRANSPORT.md.
// Every vector is Cartesian in the same Sun-centred inertial HCI frame.
struct Thermodynamics {
  double gammaAdiabatic=0;
  double magneticPermeabilityNPerA2=0;
};

struct Primitive {
  double densityKgM3=0;
  CoronalCME::Vec3 velocityMPerS;
  double pressurePa=0;
  CoronalCME::Vec3 magneticFieldT;
};

struct ConservedDensity {
  double massDensityKgM3=0;
  CoronalCME::Vec3 momentumDensityKgM2S;
  double totalEnergyDensityJPerM3=0;
  CoronalCME::Vec3 magneticFieldT;
};

struct FaceGeometry {
  // Unit normal points from the left cell into the right cell. normalSpeed
  // is the physical face velocity dotted with that normal, in m/s, not the
  // fluid's relative speed. A unit-normal roundoff allowance is numerical
  // only; invalid normals are rejected instead of silently normalised.
  CoronalCME::Vec3 normalLeftToRight;
  double normalSpeedMPerS=0;
};

struct Flux {
  double massKgM2S=0;
  CoronalCME::Vec3 momentumPa;
  double energyWPerM2=0;
  CoronalCME::Vec3 inductionTMPerS;
};

struct FaceFlux {
  Flux value;
  double leftWaveMPerS=0;
  double rightWaveMPerS=0;
  double maximumRelativeSignalMPerS=0;
  FaceGeometry geometry;
};

struct Cell {
  double volumeM3=0;
  ConservedDensity density;
};

struct ExtensiveTransfer {
  // Positive is from left to right. The B entry is a volume integral of
  // cell-centred B, not a face magnetic flux or a proof of div(B)=0.
  double massKg=0;
  CoronalCME::Vec3 momentumKgMPerS;
  double totalEnergyJ=0;
  CoronalCME::Vec3 magneticVolumeIntegralTM3;
};

struct ExchangeCandidate {
  Cell left,right;
  ExtensiveTransfer transfer;
};

// E=p/(gamma-1)+rho|v|^2/2+|B|^2/(2mu0), SI J/m^3. Invalid/nonfinite
// thermodynamics or nonpositive density/pressure reject; no floor is applied.
Core::Result<ConservedDensity> ToConserved(
    const Primitive& state,const Thermodynamics& gas);
Core::Result<Primitive> ToPrimitive(
    const ConservedDensity& state,const Thermodynamics& gas);

// Canonical ideal-MHD fast speed along a unit HCI normal. Unlike the legacy
// radial reduced speed, this includes both normal and transverse magnetic
// components. Near equal acoustic/Alfven speeds the discriminant is evaluated
// as a sum of nonnegative terms instead of subtracting nearly equal squares.
Core::Result<double> FastSpeedMPerS(const Primitive& state,
    CoronalCME::Vec3 normal,const Thermodynamics& gas);

// HLL flux with waves measured relative to the face and the complete ALE
// correction F_n - w_n U. Left/right B.n must agree within floating-point
// dot-product error; no field projection, monopole source or cleaning occurs.
// This is a bounded two-wave solver, not an unconditional positivity theorem
// or a qualification of rotational/compound MHD discontinuities.
Core::Result<FaceFlux> HllFlux(const Primitive& left,const Primitive& right,
    const FaceGeometry& geometry,const Thermodynamics& gas);

// Apply ONE already evaluated shared flux to both stationary fixed-volume
// cells with opposite signs. Both private candidates must have positive
// density/thermal pressure before the immutable inputs can be superseded.
// Moving-face use returns UnsupportedCapability until an actual moving mesh
// satisfies geometric conservation; an ALE flux alone is insufficient.
Core::Result<ExchangeCandidate> TransferStationaryFace(
    const Cell& left,const Cell& right,const FaceFlux& flux,
    double areaM2,double stepS,const Thermodynamics& gas);

} } } // namespace SEP::CoronaSwcme::CrossRay

#endif
