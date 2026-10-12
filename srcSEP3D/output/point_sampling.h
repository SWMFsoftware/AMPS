// ============================================================================
// srcSEP3D point sampling: channel grid, accumulation, normalization and the
// Tecplot product of the 'sampling: points' input block.
//
// This header is AMPS-, PIC- and MPI-independent.  The AMPS adapter
// (amps/amps_point_sampler.cpp) owns the particle walk, the MPI reduction and
// the registration with the AMPS sampling cycle; everything numerical lives
// here so that it is unit tested against hand-computed fixtures (stage-1
// tests PS3D01-PS3D05) without a mesh.
//
// Sampling volume.  A sampling point samples exactly the macroparticles that
// AMPS links to the mesh cell containing the point, and divides by that
// cell's volume V [m^3] (user decision 2026-10-11).  Particles are therefore
// counted where AMPS stores them, at the resolution of the mesh.
//
// Two reference frames (user decision 2026-10-11, "both").  srcSEP3D stores
// the scalar momentum p and the pitch cosine mu of a particle in the LOCAL
// PLASMA FRAME (transport/focused_transport.h: the guiding centre moves with
// U + mu v b).  Every product is written twice:
//   plasma frame    E, v from p; v_par = mu v (mu relative to the particle's
//                   local magnetic direction, as stored by the mover);
//   inertial frame  the plasma-frame velocity vector v' (the AMPS particle
//                   velocity, rebuilt by the mover from p, mu and the
//                   gyrophase) is Lorentz-boosted with the bulk velocity U of
//                   the sampled cell:
//                     u = [v'/g + U + g/(1+g) (U.v'/c^2) U] / (1 + U.v'/c^2),
//                     g = 1/sqrt(1 - U^2/c^2),
//                   E = m c^2 (gamma(u) - 1), v = |u|, v_par = u.b (b = the
//                   cell's magnetic direction).  This is what an instrument
//                   at rest in the heliocentric frame measures.
// A particle is assigned to an energy channel separately in each frame, so a
// row of the product refers to particles with that frame's energy in the
// channel.
//
// Quantities.  For a point P, species s, frame f and channel k, with
// N = the number of sampling calls in the AMPS output interval
// (the AMPS LastSampleLength) and the sums running over those calls and over the
// particles in the cell whose frame-f kinetic energy lies in channel k:
//
//   density                 n_k   = (1/N) sum w / V                  [m^-3]
//   omnidirectional flux    F_k   = (1/N) sum w v / V                [m^-2 s^-1]
//   field-aligned flux      F|_k  = (1/N) sum w v_par / V            [m^-2 s^-1]
//                           (positive along +B; net streaming along b)
//   differential intensity  j_k   = F_k / (4 pi dE_k)                [m^-2 s^-1 sr^-1 MeV^-1]
//                           (the isotropic-equivalent intensity; exact
//                           for an isotropic distribution, otherwise the
//                           solid-angle average of the intensity)
//   interval fluence        Phi_k = sum w v dt / V                   [m^-2]
//                           (sum over calls of the instantaneous flux times
//                           the species time step dt of that call)
//   cumulative fluence      the running sum of Phi_k since the first sample
//                           of the run, and its intensity Phi_k/(4 pi dE_k)
//
// w is the particle's statistical weight (global species weight times the
// individual correction).  The product is written in the units the user
// selected: cm^-3, cm^-2 s^-1, cm^-2 s^-1 sr^-1 MeV^-1, cm^-2 and
// cm^-2 sr^-1 MeV^-1 (factors 1e-6 and 1e-4 from the SI sums).  Energies are
// kinetic energies per particle in MeV.
//
// The number of macroparticle samples in each channel (summed over the N
// calls) is written with every row so that a reader can judge the Monte
// Carlo noise; no uncertainty model is implied.  Particles outside the
// channel range are not dropped silently: their density, flux and sample
// counts below the lowest and above the highest edge are written as zone
// auxiliary data.
// ============================================================================

#ifndef SEP3D_OUTPUT_POINT_SAMPLING_H
#define SEP3D_OUTPUT_POINT_SAMPLING_H

#include "../core/sep3d_types.h"

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace SEP3D {
namespace Output {
namespace PointSampling {

// One sampling location.  'origin' records where it was defined
// ("input.in:41" or "points.txt:7") for diagnostics.
struct PointDefinition {
  std::string name;
  Core::Vec3 positionM;
  std::string origin;
  // Centre and edge lengths [m] of the AMPS cell whose particles the point
  // samples, filled when the point is bound to the mesh; written into the
  // product so a reader knows the actual sampling volume.
  bool cellKnown = false;
  Core::Vec3 cellCenterM;
  Core::Vec3 cellSizeM;
};

// Channel spacing of an energy grid.
enum class EnergyScale { Linear, Logarithmic };

// The energy channels of one sampled species.  'edgesMeV' has
// channels + 1 strictly increasing kinetic-energy edges [MeV].
struct SpeciesChannels {
  int ampsIndex = -1;       // AMPS species index (compiled SpeciesList)
  std::string symbol;       // AMPS chemical symbol, e.g. "H_PLUS"
  double massKg = 0.0;      // rest mass [kg]
  EnergyScale scale = EnergyScale::Logarithmic;
  std::vector<double> edgesMeV;
  std::size_t channels() const {
    return edgesMeV.empty() ? 0 : edgesMeV.size() - 1;
  }
};

// Build the channel edges [MeV].  Linear: E_i = E_min + i (E_max-E_min)/n.
// Logarithmic: E_i = E_min (E_max/E_min)^(i/n).  The end points are set
// exactly to the configured values (no round-off drift at E_max).  Requires
// 0 < E_min < E_max (E_min may be 0 for a linear grid) and n >= 1.
Core::Status BuildEnergyEdges(double minimumMeV, double maximumMeV,
                              unsigned channels, EnergyScale scale,
                              std::vector<double>* edgesMeV);

// Channel of a kinetic energy [MeV]: 0..n-1, kBelowRange for E < E_min,
// kAboveRange for E > E_max.  A value equal to an interior edge belongs to
// the upper channel; E == E_max belongs to the last channel, so a closed
// range [E_min, E_max] is covered.  Binary search: O(log n).
constexpr int kBelowRange = -1;
constexpr int kAboveRange = -2;
int ChannelIndex(const std::vector<double>& edgesMeV, double energyMeV);

// Representative channel energy written in the product: the arithmetic
// midpoint for a linear grid, the geometric mean for a logarithmic grid.
double ChannelCenterMeV(const SpeciesChannels& species, std::size_t channel);

// ---- Relativistic kinematics (SI) -----------------------------------------
//
// Kinetic energy from the momentum magnitude, written without cancellation
// for p << m c:  E = p^2 c^2 / (sqrt(p^2 c^2 + m^2 c^4) + m c^2)   [J].
double KineticEnergyFromMomentumJ(double momentumKgMPerS, double massKg);
// Speed from the momentum magnitude: v = p c^2 / sqrt(p^2 c^2 + m^2 c^4).
double SpeedFromMomentumMPerS(double momentumKgMPerS, double massKg);
// Kinetic energy from a speed below c, without cancellation:
// gamma - 1 = beta^2 / (sqrt(1-beta^2) (1 + sqrt(1-beta^2))).
double KineticEnergyFromSpeedJ(double speedMPerS, double massKg);
// Relativistic velocity addition: velocity in the inertial frame of a
// particle moving with 'velocityInMovingFrame' in a frame that moves with
// 'frameVelocity' (both |.| < c).  See the header comment for the formula.
Core::Vec3 BoostVelocity(const Core::Vec3& velocityInMovingFrame,
                         const Core::Vec3& frameVelocity);

// Joule <-> MeV for kinetic energies (1 MeV = 1e6 e J, CODATA 2018 e).
constexpr double kJoulePerMeV = 1.0e6 * Core::Const::e;

// One sampled macroparticle, reduced to what the accumulation needs.
struct ParticleInput {
  double weight = 0.0;                // statistical weight
  double massKg = 0.0;                // rest mass [kg]
  double momentumKgMPerS = 0.0;       // plasma-frame |p| [kg m/s]
  double mu = 0.0;                    // plasma-frame pitch cosine
  Core::Vec3 plasmaVelocityDirection; // direction of the plasma-frame velocity
                                      // (AMPS particle velocity; any length)
};

// The frame-resolved kinematics of one particle (frame 0 = plasma,
// frame 1 = inertial).
constexpr int kPlasmaFrame = 0;
constexpr int kInertialFrame = 1;
constexpr int kFrames = 2;
struct ParticleKinematics {
  double energyMeV[kFrames] = {0.0, 0.0};
  double speedMPerS[kFrames] = {0.0, 0.0};
  double parallelSpeedMPerS[kFrames] = {0.0, 0.0};  // v.b, signed
};

// Evaluate both frames.  'bulkVelocityMPerS' and 'magneticDirection' are the
// background of the sampled cell (magneticDirection must be a unit vector).
// Fails for non-finite input, a non-positive mass or weight, |U| >= c, or a
// zero velocity direction with a non-zero momentum.
Core::Status EvaluateKinematics(const ParticleInput& particle,
                                const Core::Vec3& bulkVelocityMPerS,
                                const Core::Vec3& magneticDirection,
                                ParticleKinematics* kinematics);

// ---- Accumulator layout --------------------------------------------------
//
// All sums of one process live in one flat std::vector<double> so the AMPS
// adapter reduces them with a single MPI_Reduce(MPI_SUM).  For every
// (point, species, frame) block the layout is
//   channels x {density, flux, field-aligned flux, fluence, samples}
//   followed by {below density, below flux, below samples,
//                above density, above flux, above samples}.
// 'samples' are macroparticle counts (summed over calls) stored as doubles;
// they are exact integers below 2^53.
enum ChannelQuantity : std::size_t {
  kDensity = 0,          // sum w / V                [m^-3]
  kFlux = 1,             // sum w v / V              [m^-2 s^-1]
  kFieldAlignedFlux = 2, // sum w v_par / V          [m^-2 s^-1]
  kFluence = 3,          // sum w v dt / V           [m^-2]
  kSamples = 4,          // macroparticle samples
  kChannelQuantities = 5
};
enum TailQuantity : std::size_t {
  kBelowDensity = 0, kBelowFlux = 1, kBelowSamples = 2,
  kAboveDensity = 3, kAboveFlux = 4, kAboveSamples = 5,
  kTailQuantities = 6
};

class Layout {
 public:
  Layout() = default;
  Layout(std::size_t points, const std::vector<SpeciesChannels>& species);

  std::size_t points() const { return points_; }
  std::size_t species() const { return speciesStride_.size(); }
  std::size_t size() const { return points_ * pointStride_; }

  // Offset of channel 'channel', quantity 'q' of (point, slot, frame).
  std::size_t Channel(std::size_t point, std::size_t slot, int frame,
                      std::size_t channel, std::size_t q) const;
  // Offset of tail quantity 'q' of (point, slot, frame).
  std::size_t Tail(std::size_t point, std::size_t slot, int frame,
                   std::size_t q) const;

 private:
  std::size_t FrameBase(std::size_t point, std::size_t slot, int frame) const;
  std::size_t points_ = 0;
  std::size_t pointStride_ = 0;
  std::vector<std::size_t> channels_;       // per slot
  std::vector<std::size_t> speciesOffset_;  // per slot, within a point
  std::vector<std::size_t> speciesStride_;  // per slot, both frames
};

// Add one particle's contribution of one sampling call.  'cellVolumeM3' is
// the AMPS cell volume (> 0) and 'timeStepS' the species time step of the
// call (> 0, used only by the fluence).
void Accumulate(const Layout& layout, const SpeciesChannels& species,
                std::size_t point, std::size_t slot,
                double weight, double cellVolumeM3, double timeStepS,
                const ParticleKinematics& kinematics,
                std::vector<double>* sums);

// Add the interval fluence of 'sums' to the running cumulative fluence
// ('cumulative' has layout.size() entries; only the kFluence slots are used).
void FoldFluence(const Layout& layout,
                 const std::vector<SpeciesChannels>& species,
                 const std::vector<double>& sums,
                 std::vector<double>* cumulative);

// Products written by the 'products' key of the input block.
struct ProductSelection {
  bool density = false;
  bool flux = false;      // F and F|
  bool spectrum = false;  // j
  bool fluence = false;   // interval and cumulative fluence (+ intensity)
};
// Parse "flux, density, spectrum, fluence" (any order, comma separated,
// each at most once, at least one).
Core::Status ParseProducts(const std::string& text, ProductSelection* products);

// Interval metadata written into the product.
struct IntervalInfo {
  long outputFileNumber = 0;    // AMPS DataOutputFileNumber of this file
  long sampleCalls = 0;         // N = the AMPS LastSampleLength
  double firstSampleTimeS = 0;  // simulation time of the first call
  double lastSampleTimeS = 0;   // simulation time of the last call
};

// Format the complete Tecplot ASCII product of one interval: one zone per
// (point, species), one row per channel.  'sums' are the reduced interval
// sums (all processes), 'cumulative' the cumulative fluence after
// FoldFluence() of this interval.  Requires info.sampleCalls >= 1.
std::string FormatTecplot(const Layout& layout,
                          const std::vector<PointDefinition>& points,
                          const std::vector<SpeciesChannels>& species,
                          const ProductSelection& products,
                          const IntervalInfo& info,
                          const std::vector<double>& sums,
                          const std::vector<double>& cumulative);

// Parse a point file: one point per line, "x y z [name]" with x, y, z in
// metres in the AMPS heliocentric Cartesian frame.  '!' and '#' start a
// comment; blank lines are ignored.  A point without a name is named by its
// zero-based index among the points of the file ("0", "1", ...), matching the
// zero-based numbering of AMPS sample points.  Errors name
// 'displayName:line'.
Core::Status ParsePointFile(const std::string& text,
                            const std::string& displayName,
                            std::vector<PointDefinition>* points);

}  // namespace PointSampling
}  // namespace Output
}  // namespace SEP3D

#endif  // SEP3D_OUTPUT_POINT_SAMPLING_H
