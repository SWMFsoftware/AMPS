#ifndef SEP_COMMON_PITCH_ANGLE_DIFFUSION_PITCH_ANGLE_KERNELS_H
#define SEP_COMMON_PITCH_ANGLE_DIFFUSION_PITCH_ANGLE_KERNELS_H

// ============================================================================
// Pure pitch-angle scattering kernels and the independent finite-volume
// reference operator (PITCH_ANGLE_SCATTERING_MODEL.md v1.2, Sections 16.3,
// 16.4, 16.7, 16.9 and 16.12; roadmap stages PS04 and PS06).
//
// Ownership.  These functions realize only the scattering operator (21) (and,
// for the Ito kicks, add a host-supplied deterministic pitch-angle rate once,
// as Eq. (25) prescribes).  They consume no random-number generator: the
// caller draws the standard normal or uniform variates from its own keyed
// stream and passes them in, so RNG ownership, reproducibility across ranks
// and restart state stay with the host (Section 17, "RNG contract").  They
// never stream particles in space; Section 16.12 requires the host to
// integrate ds = v mu dt along the subcycled trajectory.
//
// Boundary at mu = +-1 (Section 16.4, U-60).  The exact process does not
// reach the endpoints when D vanishes linearly there (Feller test, derived),
// but a discrete step can leave [-1, 1].  The Ito kicks therefore apply the
// mirror rule Eq. (161) repeatedly until the value is inside ("folding",
// implemented in closed form as a period-4 triangular wave) and report the
// number of reflections.  The random rotation stays in [-1, 1] by
// construction (Cauchy-Schwarz) and needs no rule.
// ============================================================================

#include "pitch_angle_diffusion.h"

#include <cstddef>
#include <vector>

namespace SEP {
namespace PitchAngleDiffusion {
namespace Kernels {

// Repeated mirror reflection about mu = +-1 (Eq. 161 applied until the value
// is in [-1, 1]).  O(1) for any finite input; NaN/inf are returned unchanged
// so the caller's state check rejects them.
double FoldPitchCosine(double mu, unsigned* reflections = nullptr);

struct KickResult {
  Status status;
  double mu = 0.0;               // folded result
  double driftIncrement = 0.0;   // (dmu/dt + D') dt
  double noiseIncrement = 0.0;   // sqrt(2 D dt) zeta (+ Milstein term)
  unsigned reflections = 0;
};

// Euler-Maruyama, Eq. (157):
//   mu_{n+1} = mu_n + (dmu/dt + D') dt + sqrt(2 D dt) zeta,
// all coefficients at mu_n (Ito).  One-step moments (F-PUSH-04a/b):
//   E[dmu] = A dt,  E[dmu^2] = A^2 dt^2 + 2 D dt.
KickResult EulerMaruyamaKick(double mu, double dMuMuPerS,
                             double dDmuMuDmuPerS,
                             double deterministicMuRatePerS, double dtS,
                             double standardNormal);

// Milstein, Eq. (159): Euler-Maruyama plus (1/2) D' dt (zeta^2 - 1), since
// b b' = D' for b = sqrt(2D) wherever D > 0 (Section 16.3, F-PUSH-03).
//   E[dmu^2] = A^2 dt^2 + 2 D dt + (1/2) D'^2 dt^2.
KickResult MilsteinKick(double mu, double dMuMuPerS, double dDmuMuDmuPerS,
                        double deterministicMuRatePerS, double dtS,
                        double standardNormal);

// Small-angle random rotation for isotropic scattering of frequency nu,
// Eqs. (167)-(168) (Agueda 2008; Wijsen et al. 2019, MC-ROT): theta^2 is
// exponential with mean 2 nu dt, the azimuth phi uniform, and
//   mu' = mu cos(theta) + sqrt(1 - mu^2) sin(theta) cos(phi).
// uniformTheta in (0, 1] gives theta^2 = -2 nu dt ln(uniformTheta);
// uniformPhi in [0, 1) gives phi = 2 pi uniformPhi.  Restricted to the
// isotropic law D = (nu/2)(1 - mu^2) (Section 16.11); its exact one-step
// moments are Eq. (169).
KickResult RandomRotationKick(double mu, double nuPerS, double dtS,
                              double uniformTheta, double uniformPhi);

// ---------------------------------------------------------------------------
// Finite-volume reference operator (Eq. 178, Section 16.12).
//
// Cells [e_i, e_{i+1}] with e_0 = -1 and e_n = +1.  Neighbouring cell centres
// are connected by the reciprocal-integral transmissibility
//
//     T_{i,i+1} = [ Int_{c_i}^{c_{i+1}} dmu / D(mu) ]^(-1),
//
// never by sampling D at the interface, so an integrable internal resistance
// (PA-QFORM, H = 0) keeps a finite cross-zero conductance and a nonintegrable
// one (PA-EPS, eps = 0) disconnects the hemispheres.  With stationary masses
// m_i (cell widths for pure scattering) the continuous-time Markov chain has
// rates Q_{i->j} = T_ij / m_i (rows are origin states), Q_ii = -sum_j Q_ij,
// conserves number, and satisfies detailed balance m_i Q_ij = m_j Q_ji.
// No endpoint face exists: zero flux at mu = +-1.
struct FvGenerator {
  std::vector<double> edges;           // n+1 strictly increasing, -1 .. 1
  std::vector<double> centres;         // n cell centres
  std::vector<double> stationaryMass;  // n stationary masses (unnormalized)
  std::vector<double> conductance;     // n-1 values T_{i,i+1}
  std::vector<double> rateLeft;        // Q_{i->i-1}
  std::vector<double> rateRight;       // Q_{i->i+1}
  std::vector<double> rateDiagonal;    // Q_ii
  double quadratureRelativeError = 0.0;  // largest estimate over faces
  std::string internalExtension;  // "connected" or "disconnected at mu = 0"
};

// Pure scattering of `model` with amplitude D0 [s^-1] (EvaluateShape
// conventions: H = 0 and eps = 0 are admitted).  Stationary masses are the
// cell widths (uniform stationary density).
Status BuildFvGenerator(ModelId model, const ShapeParameters& shape,
                        double amplitudeD0PerS,
                        const std::vector<double>& edges, FvGenerator* out);

// Generic assembly from conductances and stationary masses (used, for
// example, by the constant-focusing reference of Section 16.12 with
// T^rho = [Int (D rho_+)^(-1)]^(-1) and m_i = Int_cell rho_+).
Status BuildFvGeneratorFromConductances(const std::vector<double>& edges,
                                        const std::vector<double>& masses,
                                        const std::vector<double>& conductance,
                                        FvGenerator* out);

// Exact time propagation of a mass vector, p(t) = p(0) exp(Q t), through the
// eigen-decomposition of the symmetric tridiagonal matrix
// S = M^(-1/2) K M^(-1/2), K = M Q (implicit QL with Wilkinson shifts).
struct FvPropagator {
  std::size_t n = 0;
  std::vector<double> sqrtMass;
  std::vector<double> eigenvalues;   // n, all <= 0 up to roundoff
  std::vector<double> eigenvectors;  // n x n, column k is vector k (row-major)
};
Status PrepareFvPropagator(const FvGenerator& generator, FvPropagator* out);
Status PropagateFv(const FvPropagator& propagator,
                   const std::vector<double>& initialMass, double timeS,
                   std::vector<double>* mass);

// Asymptotic (Green-Kubo) variance rate of the additive functional
// Int f(mu(t)) dt under the reversible chain: <f, (-Q)^(-1) f>_pi with pi the
// normalized stationary mass and f the cell values (f is centred by the
// function).  Solved exactly through the flux recursion
// T_{i,i+1} (phi_i - phi_{i+1}) = sum_{k<=i} m_k f_k (derived), so no matrix
// is inverted.  For f = mu and constant v, Var[z]/(2t) -> v^2 times it.
Status FvAsymptoticVariance(const FvGenerator& generator,
                            const std::vector<double>& cellValues,
                            double* variance);

// Legendre polynomial P_l(x) (three-term recurrence) and the exact cell
// average (1/w) Int_cell P_l = (P_{l+1} - P_{l-1})|/((2l+1) w), l >= 1.
double Legendre(int l, double x);
double LegendreCellAverage(int l, double left, double right);

}  // namespace Kernels
}  // namespace PitchAngleDiffusion
}  // namespace SEP

#endif  // SEP_COMMON_PITCH_ANGLE_DIFFUSION_PITCH_ANGLE_KERNELS_H
