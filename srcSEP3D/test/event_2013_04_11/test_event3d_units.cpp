// ============================================================================
// srcSEP3D/test/event_2013_04_11/test_event3d_units.cpp
//
// Unit tests E3D01-E3D06 of the srcSEP3D components added for the 3-D
// version of the 2013-04-11 real-event test (README.md, sections 2 and 5).
//
// WHY THESE TESTS EXIST:
//   The nightly Table tests (SEP3D--event-2013-04-11--*) compare whole runs
//   with AMPS references.  They detect a change of the results, but not an
//   error that was already present when the references were made.  These
//   tests compare each new component with an analytic result on a small
//   input, without AMPS, PIC or MPI:
//
//     E3D01  computation frame                core/computation_frame.h
//     E3D02  inertial terms of the focused    transport/focused_transport.*
//            movers (Coriolis/centrifugal
//            option)
//     E3D03  ephemeris table and its          output/observer_trajectory.*
//            interpolation
//     E3D04  field-line tracing and the       mesh/observer_flux_tube.*
//            observer flux tube
//     E3D05  absorbing-shock test             adapters/shock_return.*
//     E3D06  1-D-equivalent tube-line source  adapters/tube_line_source.*
//
//   E3D07 (the SPICE body-observer path) needs a CSPICE toolkit and lives in
//   test_spice_ephemeris.cpp.
//
// HOW IT RUNS:
//   run_unit_tests.sh compiles the canonical sources together with this file
//   under AddressSanitizer and UndefinedBehaviorSanitizer and runs it.  Each
//   check prints one "PASS <id> ..." or "FAIL <id> ..." line; the program
//   prints a summary and exits 1 when any check failed, 0 otherwise.
//
// TOLERANCES:
//   Deterministic checks use round-off tolerances.  The three Monte Carlo
//   checks (E3D05d, E3D06c, E3D06f) draw from keyed random streams
//   (transport/keyed_random.h), so every run draws the same numbers and the
//   result is reproducible; their tolerances are stated in standard errors of
//   the estimate at the comment of each check.
// ============================================================================
#include "adapters/shock_return.h"
#include "adapters/tube_line_source.h"
#include "core/computation_frame.h"
#include "mesh/observer_flux_tube.h"
#include "output/observer_trajectory.h"
#include "transport/focused_transport.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

using namespace SEP3D;
using Core::Vec3;

namespace {
// Number of failed checks; main() returns non-zero when it is positive.
int failures = 0;
// Print one PASS/FAIL line for a check and count a failure.
void Check(bool condition, const std::string& what) {
  std::printf("%s %s\n", condition ? "PASS" : "FAIL", what.c_str());
  if (!condition) ++failures;
}
// |a - b| <= absolute + relative * max(|a|, |b|): a mixed relative/absolute
// comparison, so that a value expected to be zero can use the absolute part.
bool Near(double a, double b, double relative, double absolute = 0.0) {
  return std::fabs(a - b) <= absolute + relative * std::max(std::fabs(a), std::fabs(b));
}
// OMEGA is the event's solar rotation rate [rad/s], the rotation rate of the
// corotating frame of the 2013-04-11 decks (README.md section 2.1).
const double AU = Core::Const::AU, OMEGA = 2.86533e-6;

// ---- E3D01 computation frame ------------------------------------------------
// The corotating frame at the event's rate OMEGA with reference epoch 100 s.
// At t = 5e4 s the frame has turned by OMEGA (t - 100 s) = 0.143 rad, so every
// rotation below is a real one.
//  a  frame -> HCI -> frame position round trip.  Tolerance 1e-4 m at
//     |x| = 1.25e11 m: a few units of round-off (~1e-15 relative).
//  b  a rigidly corotating flow U_hci = Omega x x is at rest in the frame
//     (U' = U - Omega x x).  The frame velocity subtracts the same cross
//     product, so the result is zero to round-off; the limit 1e-9 m/s is
//     tiny against the corotation speed there (3.5e5 m/s).
//  c  the velocity gradient of rigid rotation, grad(Omega x x), whose only
//     nonzero entries are d(Omega x x)_x/dy = -Omega and
//     d(Omega x x)_y/dx = +Omega, is removed by FromInertialVelocityGradient.
//     The removal cancels the two entries exactly, so the limit 1e-20 1/s
//     (against Omega ~ 3e-6 1/s) is a test for exact zero.
//  d  for an arbitrary gradient g and unit vector b, div U = tr(grad U) and
//     the field-aligned strain bb:grad U are the same in both frames: the
//     rotation is orthogonal, and the removed rigid-rotation part is
//     antisymmetric (zero trace, zero bb contraction).  These are the two
//     invariants behind the statement that the transport coefficients do not
//     depend on the frame (README.md section 2.1).  Round-off tolerance:
//     relative 1e-12, absolute 1e-22 1/s.
void TestFrame() {
  const auto frame = Core::ComputationFrame::Corotating(OMEGA, 100.0);
  const Vec3 x(1.2e11, -3.0e10, 2.0e10);
  const double t = 5.0e4;
  const Vec3 back = frame.ToInertialPosition(frame.FromInertialPosition(x, t), t);
  Check((back - x).Norm() < 1e-4, "E3D01a position round trip");
  // A rigidly corotating flow is at rest in the frame.
  const Vec3 corotation = frame.AngularVelocity().Cross(x);
  Check(frame.FromInertialVelocity(corotation, x, t).Norm() < 1e-9,
        "E3D01b corotation flow is at rest in the frame");
  // grad(Omega x x) vanishes in the frame; div U and bb:grad U are invariant.
  Core::ComputationFrame::Tensor rigid{0, -OMEGA, 0, OMEGA, 0, 0, 0, 0, 0};
  const auto removed = frame.FromInertialVelocityGradient(rigid, t);
  double maximum = 0;
  for (double v : removed) maximum = std::max(maximum, std::fabs(v));
  Check(maximum < 1e-20, "E3D01c grad(Omega x x) is removed");
  Core::ComputationFrame::Tensor g{1e-5, 2e-6, -3e-6, 4e-7, -2e-5, 5e-6, 1e-6, 7e-7, 3e-6};
  const auto gf = frame.FromInertialVelocityGradient(g, t);
  const double div = g[0] + g[4] + g[8], divf = gf[0] + gf[4] + gf[8];
  const Vec3 b = Vec3(0.3, -0.5, 0.81).Normalized();
  const Vec3 bf = frame.FromInertialVector(b, t);
  double strain = 0, strainf = 0;
  const double bv[3] = {b.x, b.y, b.z}, bfv[3] = {bf.x, bf.y, bf.z};
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j) {
      strain += bv[i] * bv[j] * g[3 * i + j];
      strainf += bfv[i] * bfv[j] * gf[3 * i + j];
    }
  Check(Near(div, divf, 1e-12, 1e-22) && Near(strain, strainf, 1e-12, 1e-22),
        "E3D01d div U and bb:grad U are frame invariant");
}

// ---- E3D02 inertial terms: a uniform force along b ---------------------------
// The optional rotating-frame terms of the focused movers
// (application.sep.transport.frame_inertial_forces, README.md section 2.6),
// with every other focused-transport coefficient zero:
//
//   dmu/dt = -(1 - mu^2) A,   d ln p/dt = -mu A,   A = (b.a)/v [1/s],
//
// a the frame acceleration of the plasma.  Together they are the parallel
// component of a uniform force -p A = -gamma m (b.a) along b:
//
//   d(p mu)/dt     = p dmu/dt + mu dp/dt = -p A,
//   d ln p_perp/dt = d ln p/dt - mu/(1 - mu^2) dmu/dt = 0.
//
// A is held constant (the speed v that defines it is fixed), and the two
// equations are integrated over 1 s with 1000 midpoint steps of 1 ms.
//  a  p_par changes by -p0 A x 1 s.  Over the second, ln p changes by only
//     ~mu A ~ 1e-6, so replacing p(t) by p0 in the integral of -p A changes
//     p_par by only ~5e-12 of p_par0.  The tolerance, 1e-6 of p_par, is about 13 % of the
//     change being checked (p0 A = 7.5e-6 p_par0).
//  b  p_perp is conserved to 1e-9 relative.  An error in either term would
//     change p_perp by order mu A x 1 s ~ 1e-6 relative.
//  c  a default FocusedLocalState (all coefficients zero, the inertial term
//     defaulting to zero) gives exactly zero rates: the option adds nothing
//     unless it is set.
void TestInertialForce() {
  Transport::FocusedLocalState local;
  const double v = 5.0e7, A = 3.0e-6;   // A = b.a / v [1/s]
  local.frameAccelerationOverSpeedPerS = A;
  const double mu0 = 0.4, p0 = 1.0e-20, dt = 1.0e-3;
  double mu = mu0, lnp = std::log(p0);
  // Midpoint (second-order Runge-Kutta) integration over 1 s: both updates
  // use the rates at the half-step pitch m1.  l1, the log-momentum rate at
  // the start of the step, is computed but enters multiplied by zero.
  for (int i = 0; i < 1000; ++i) {   // midpoint integration over 1 s
    const double m1 = mu + 0.5 * dt * Transport::FocusedPitchDriftPerS(mu, v, local);
    const double l1 = Transport::FocusedLogMomentumRatePerS(mu, local);
    lnp += dt * Transport::FocusedLogMomentumRatePerS(m1, local) + 0.0 * l1;
    mu += dt * Transport::FocusedPitchDriftPerS(m1, v, local);
  }
  const double p = std::exp(lnp);
  // dp_par/dt = -p A, p_perp conserved (nonrelativistic speed held fixed).
  const double ppar0 = p0 * mu0, pperp0 = p0 * std::sqrt(1 - mu0 * mu0);
  Check(Near(p * mu, ppar0 - p0 * A * 1.0, 1e-6),
        "E3D02a inertial term changes p_par by -gamma m (b.a) dt");
  Check(Near(p * std::sqrt(1 - mu * mu), pperp0, 1e-9),
        "E3D02b inertial term conserves p_perp");
  Transport::FocusedLocalState none;
  Check(Transport::FocusedPitchDriftPerS(0.3, v, none) == 0.0 &&
        Transport::FocusedLogMomentumRatePerS(0.3, none) == 0.0,
        "E3D02c zero inertial term reproduces the released coefficients");
}

// ---- E3D03 ephemeris interpolation -------------------------------------------
// The body-observer table path (ephemeris = table, README.md section 3).
// x(t) is a cubic polynomial and v(t) = dx/dt; four hourly samples of it are
// written in the checked-in table format.  Cubic Hermite interpolation from
// the positions and velocities at both ends of an interval is exact for a
// cubic, so:
//  a  the table parses into its four samples;
//  b  at t = 1234.5 s (inside the second interval, [0, 3600] s) the
//     interpolated position and velocity equal the polynomial to round-off:
//     1e-3 m at |x| ~ 1e11 m and 1e-9 m/s at |v| ~ 3e4 m/s (both ~1e-14
//     relative);
//  c  t = 9000 s, after the last sample (7200 s), is refused: there is no
//     extrapolation;
//  d  a header other than time_s,x_m,y_m,z_m,vx_m_per_s,vy_m_per_s,vz_m_per_s
//     is refused.
void TestEphemeris() {
  // A cubic trajectory is reproduced exactly by cubic Hermite interpolation.
  auto x = [](double t) { return Vec3(1e11 + 3e4 * t + 2e-1 * t * t - 1e-6 * t * t * t,
                                      -2e10 + 1e3 * t, 5e9 - 7e-2 * t * t); };
  auto v = [](double t) { return Vec3(3e4 + 4e-1 * t - 3e-6 * t * t, 1e3, -1.4e-1 * t); };
  std::string text = "# test\ntime_s,x_m,y_m,z_m,vx_m_per_s,vy_m_per_s,vz_m_per_s\n";
  for (double t : {-3600.0, 0.0, 3600.0, 7200.0}) {
    char row[256];
    const Vec3 p = x(t), w = v(t);
    std::snprintf(row, sizeof(row), "%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n",
                  t, p.x, p.y, p.z, w.x, w.y, w.z);
    text += row;
  }
  std::vector<Output::EphemerisSample> samples;
  Check(Output::ParseEphemerisTable(text, "test", &samples).ok() && samples.size() == 4,
        "E3D03a ephemeris table parses");
  Vec3 p, w;
  const bool ok = Output::InterpolateEphemeris(samples, 1234.5, &p, &w).ok();
  Check(ok && (p - x(1234.5)).Norm() < 1e-3 && (w - v(1234.5)).Norm() < 1e-9,
        "E3D03b cubic Hermite reproduces a cubic trajectory");
  Check(!Output::InterpolateEphemeris(samples, 9000.0, &p, &w).ok(),
        "E3D03c no extrapolation outside the table");
  std::vector<Output::EphemerisSample> bad;
  Check(!Output::ParseEphemerisTable("time_s,x_m\n0,1\n", "bad", &bad).ok(),
        "E3D03d a malformed header is rejected");
}

// ---- E3D04 field-line tracing and the observer flux tube ----------------------
// ParkerField: an analytic Parker spiral in the corotating frame, the
// FrameFieldFunction that the tube construction consumes.  B_r = B0/r^2 with
// B0 = 5 nT (1 au)^2, and B_phi/B_r = -Omega rho/U_sw with U_sw = 400 km/s
// and rho the cylindrical radius.  The field is static in the frame, and in
// the equatorial plane its lines are the spirals dphi/dr = -Omega/U_sw.  The
// whole field is one open sector (+1: outward is +B).
const double USW = 4.0e5, B0 = 5.0e-9 * AU * AU;   // B_r = B0 / r^2
Core::Status ParkerField(const Vec3& x, Vec3* b, int* sector) {
  const double r = x.Norm(), rho = std::hypot(x.x, x.y);
  const Vec3 er = x / r, ephi = Vec3(-x.y, x.x, 0.0) / rho;
  *b = (er - ephi * (OMEGA * rho / USW)) * (B0 / (r * r));
  *sector = 1;
  return Core::Status::OK();
}

//  a  one line from (1 au, 0, 0) sunward (sign -1) to r = 0.1 au, RK4 with
//     arc step <= 1e9 m.  The end lies on r = 0.1 au to round-off (1e-12
//     relative; the last step is cut onto the sphere by bisection), and its
//     azimuth equals the analytic phi = Omega (0.9 au)/U_sw = 0.96 rad to
//     1e-6 relative.  A wrong winding sense would give phi < 0.
//
// The tube (b-e): one observer "probe" with collection radius rho = 0.02 au,
// sampled five times 900 s apart while it moves 0.03 au in azimuth along
// r = 1 au in the equatorial plane (angles 0 to -0.03 rad).  The reference
// sphere is r = 1 au (the radius of the first sample), with lattice spacing
// s = 0.005 au; lines run from 0.05 au to 1.5 au.
//  b  the tube is built and holds lines;
//  c  every kept line passes, at the radius of some sample centre, within
//     rho + s_k of that centre (the keep rule of BuildObserverFluxTube, s_k
//     the local line spacing there), and the radius of its points increases
//     strictly from the inner sphere outward (the TracedTubeLine contract);
//  d  the total flux of the kept lines equals B_r(1 au) times the area of the
//     region the collection sphere sweeps on r = 1 au: a stadium of
//     half-width rho and length 0.03 au, area pi rho^2 + 2 rho (0.03 au).  The
//     bounds shrink the half-width to rho - s and widen it to rho + 2 s: the
//     lattice resolves the edge only to one cell, a line is kept out to
//     rho + s_k from a centre, and a kept line carries the flux of its whole
//     cell.  (The five sample discs, 0.0075 au apart, fill the stadium up to
//     scallops of 3.6e-4 au, well inside one spacing.)
//  e  the stored local spacing at the first point of a line is
//     sqrt(Phi/|B|), the flux-conservation width, to round-off.
void TestTube() {
  std::vector<Vec3> line;
  const Core::Status traced = Mesh::TraceFieldLineToRadius(
      ParkerField, Vec3(AU, 0, 0), -1, 0.1 * AU, 1e9, 100000, &line);
  const Vec3 end = line.empty() ? Vec3() : line.back();
  const double phi = std::atan2(end.y, end.x);
  Check(traced.ok() && Near(end.Norm(), 0.1 * AU, 1e-12) &&
        Near(phi, OMEGA * 0.9 * AU / USW, 1e-6),
        "E3D04a RK4 trace follows the analytic Parker spiral to 1e-6");
  Mesh::TubeObserverPath path;
  path.id = "probe";
  path.collectionRadiusM = 0.02 * AU;
  for (int k = 0; k < 5; ++k) {
    // The observer moves 0.03 au in azimuth over its stay.
    const double angle = -0.0075 * k;
    path.timesS.push_back(900.0 * k);
    path.centresM.push_back(Vec3(std::cos(angle), std::sin(angle), 0.0) * AU);
  }
  Mesh::ObserverTubeControls controls;
  controls.innerRadiusM = 0.05 * AU;
  controls.outerRadiusM = 1.5 * AU;
  controls.traceStepM = 2e9;
  controls.seedSpacingM = 0.005 * AU;
  Mesh::ObserverFluxTube tube;
  const Core::Status built = Mesh::BuildObserverFluxTube({path}, controls, ParkerField, &tube);
  Check(built.ok() && !tube.lines.empty(), "E3D04b observer flux tube is traced");
  bool within = true, monotone = true;
  double flux = 0;
  for (const auto& l : tube.lines) {
    flux += l.fluxWb;
    double best = 1e300;
    for (const Vec3& c : path.centresM) {
      Vec3 at;
      double spacing = 0;
      if (Mesh::PointAtRadius(l, c.Norm(), &at, &spacing))
        best = std::min(best, (at - c).Norm() - path.collectionRadiusM - spacing);
    }
    within = within && best <= 0.0;
    for (std::size_t i = 1; i < l.pointsM.size(); ++i)
      monotone = monotone && l.pointsM[i].Norm() > l.pointsM[i - 1].Norm();
  }
  Check(within && monotone, "E3D04c every kept line meets a collection sphere; radius monotone");
  // Flux through the swept region at 1 au: a stadium of half-width rho and
  // length 0.03 au, within the lattice-cell rounding (one spacing).
  const double rho = path.collectionRadiusM, s = controls.seedSpacingM;
  const double inner = (3.14159265358979 * (rho - s) * (rho - s) + 2 * (rho - s) * 0.03 * AU) * B0 / (AU * AU);
  const double outer = (3.14159265358979 * (rho + 2 * s) * (rho + 2 * s) + 2 * (rho + 2 * s) * 0.03 * AU) * B0 / (AU * AU);
  Check(flux > inner && flux < outer, "E3D04d tube flux matches the swept collection region");
  // Local spacing follows flux conservation: s(r) = sqrt(Phi/|B|).
  const auto& first = tube.lines.front();
  Vec3 b;
  int sector = 0;
  ParkerField(first.pointsM.front(), &b, &sector);
  Check(Near(first.localSpacingM.front(), std::sqrt(first.fluxWb / b.Norm()), 1e-12),
        "E3D04e local line spacing is sqrt(Phi/|B|)");
}

// ---- E3D05 absorbing shock -------------------------------------------------------
// TestShockReturn (README.md section 2.5) on a heliocentric sphere of radius
// 1e10 m.  g(tau) = |x(tau)| - R(tau) is the particle's height above the
// front along the substep, tau in [0, 1]; the test reports a hit and the
// fraction tau of the substep at which the particle reaches the front.
//  a  front speed 1e6 m/s, particle at rest 3e7 m upstream, dt = 60 s: the
//     front advances 6e7 m, g = 3e7 (1 - 2 tau), and the hit fraction
//     g0/(g0 - g1) is 0.5 (round-off tolerance).
//  b  static front; a straight chord between two upstream ends
//     (|x| = 1e10 + 1.15e7 m) whose midpoint lies 1e6 m below the front.  With
//     no diffusion (the focused movers' test) g is convex, its minimum is
//     negative, and the first root lies before the midpoint: 0 < fraction < 0.5.
//  c  radial motion away from a static front is never absorbed.
//  d  the Parker (diffusive) test with both ends upstream, g0 = 2e6 m and
//     g1 = 3e6 m, and kappa_r dt = 1e13 m^2: the Brownian bridge hits with
//     probability exp(-g0 g1/(kappa_r dt)) = exp(-0.6) = 0.549.  Each of the
//     200000 trials has its own keyed stream (campaign seed 7, particle id
//     i + 1, purpose ShockReturnBridge); the hit frequency must lie within 4
//     binomial standard errors (sigma = 1.1e-3) of that probability.  The
//     check covers the sampling of the formula that shock_return.cpp
//     implements (the 1-D boundary's bridge, shock_return.h); it does not
//     simulate Brownian paths to re-derive the formula.
void TestShockReturn() {
  Adapters::ExpandingShock shock;
  shock.radiusAtStepStartM = 1.0e10;
  shock.radialSpeedMPerS = 1.0e6;
  shock.active = true;
  shock.generation = 1;
  // A particle at rest 3e7 m outside is overtaken during a 60 s substep.
  const Vec3 x0(1.0e10 + 3.0e7, 0, 0);
  auto hit = Adapters::TestShockReturn(x0, x0, 60.0, shock, 0.0, nullptr);
  Check(hit.hit && Near(hit.fraction, 0.5, 1e-12), "E3D05a front overtakes a particle");
  // A chord that dips below the front between two upstream ends.
  shock.radialSpeedMPerS = 0.0;
  const Vec3 a(1.0e10 - 1.0e6, -5e8, 0), b2(1.0e10 - 1.0e6, 5e8, 0);
  hit = Adapters::TestShockReturn(a, b2, 1.0, shock, 0.0, nullptr);
  Check(hit.hit && hit.fraction > 0 && hit.fraction < 0.5, "E3D05b grazing chord is caught");
  hit = Adapters::TestShockReturn(Vec3(2e10, 0, 0), Vec3(2.1e10, 0, 0), 1.0, shock, 0.0, nullptr);
  Check(!hit.hit, "E3D05c upstream motion is not absorbed");
  // Brownian bridge: hit frequency exp(-g0 g1 / (kappa dt)).
  const double g0 = 2e6, g1 = 3e6, kdt = 1.0e13;
  int hits = 0;
  const int trials = 200000;
  for (int i = 0; i < trials; ++i) {
    Transport::RandomKey key;
    key.campaignSeed = 7;
    key.particleId = static_cast<std::uint64_t>(i) + 1;
    key.purpose = Transport::RandomPurpose::ShockReturnBridge;
    Transport::KeyedRandomStream stream(key);
    if (Adapters::TestShockReturn(Vec3(1e10 + g0, 0, 0), Vec3(1e10 + g1, 0, 0), 1.0,
                                  shock, kdt, &stream).hit) ++hits;
  }
  const double expected = std::exp(-g0 * g1 / kdt);
  const double sigma = std::sqrt(expected * (1 - expected) / trials);
  Check(std::fabs(hits / static_cast<double>(trials) - expected) < 4 * sigma,
        "E3D05d Brownian-bridge crossing probability");
}

// ---- E3D06 tube-line source on an analytic shock ---------------------------------
// UniformSphere: an analytic ShockSurface (adapters/shock_surface.h) for the
// source test.  A heliocentric sphere R(t) = R0 + V t (R0 = 3e9 m,
// V = 1.5e6 m/s) that is an accepted shock everywhere, with a radial upstream
// field of constant magnitude B = 1e-5 T, upstream proton density
// n = 1e9 m^-3, radial upstream wind U = 1e5 m/s (inflow w1 = V - U) and
// compression ratio r = 3.  Intersect reports no crossing: the tube-line
// source locates the front and the release surface on a line through Level
// and ReleaseLevel (LineCrossing in tube_line_source.cpp), so Intersect is not
// exercised here.
class UniformSphere final : public Adapters::ShockSurface {
 public:
  double R0 = 3.0e9, V = 1.5e6, n = 1.0e9, U = 1.0e5, B = 1.0e-5, r = 3.0;
  double Radius(double t) const { return R0 + V * t; }
  std::string GeometryName() const override { return "test-sphere"; }
  std::string Identity() const override { return "test"; }
  double SupportBeginS() const override { return 0; }
  double SupportEndS() const override { return 1e6; }
  Core::Status Level(const Vec3& x, double t, double* level) const override {
    *level = x.Norm() - Radius(t);
    return Core::Status::OK();
  }
  Core::Status Intersect(const Vec3&, const Vec3&, double,
                         std::vector<Adapters::ShockSegmentCrossing>* c) const override {
    c->clear();
    return Core::Status::OK();
  }
  Core::Status FrontAt(const Vec3& x, double t, Adapters::ShockFrontState* s) const override {
    const Vec3 e = x.Normalized();
    s->supported = s->accepted = true;
    s->positionM = e * Radius(t);
    s->outwardNormal = e;
    s->normalSpeedMPerS = V;
    s->inflowMPerS = V - U;
    s->compressionRatio = r;
    s->upstreamProtonDensityM3 = n;
    s->upstreamMagneticFieldT = e * B;
    s->upstreamVelocityMPerS = e * U;
    return Core::Status::OK();
  }
  Core::Status ReleaseLevel(const Vec3& x, double t, double f, double* level) const override {
    *level = x.Norm() - (1 + f) * Radius(t);
    return Core::Status::OK();
  }
  double MaximumNormalSpeedMPerS(double, double) const override { return V; }
};

// The tube: three straight radial lines (directions 0.01 rad apart) from
// 1e9 m to 3e11 m with fluxes 1e9, 2e9 and 3e9 Wb, and with the local spacing
// sqrt(Phi/|B|) of the field |B| = B (3e9 m / r)^2 that the plasma callback
// also returns.  One source step of 1 s with N = 40000 samples:
//  a  the net rate of line 1 is eta n1 w1 Phi/|B1.n| with w1 = V - U, the 1-D
//     NetProtonRate (README.md section 2.4); round-off tolerance;
//  b  every stratified sample yields a candidate (the front meets every line
//     throughout the step, so no sample is empty);
//  c  the sum of the statistical weights equals sum_l R_l dt.  R_l is
//     constant here, so the sum is sum_l R_l dt times the mean importance
//     weight, whose expectation is 1.  For this spectrum (q = 4.5,
//     pmax/pmin = 105) the weight has a standard deviation of about 1.6, so
//     the mean of 40000 draws has a standard error of about 0.8 %; the 3 %
//     tolerance is about 4 standard errors;
//  d  each particle starts on the release sphere (1 + f) R(t), f = 0.3, at
//     its own birth epoch t = (1 - remainingFirstStepFraction) dt, also after
//     the transverse jitter (which keeps the radius); tolerance 1e-9 relative;
//  e  with the focused pitch law every particle has v mu + U - V > 0 (here
//     b.n = 1): a positive normal velocity relative to the front speed V,
//     which is the advection term U.n - V_n the 1-D law uses (the variable
//     `advection` in GenerateTubeLineSource).  The bound -1e-9 v allows
//     round-off at the edge of the admitted mu range;
//  f  the weighted mean of ln(p/pmin) equals that of dN/dp ~ p^k, k = 2 - q,
//     q = 3r/(r - 1) = 4.5, truncated to [pmin, pmax]:
//       <ln(p/pmin)> = [e^(aL) (L - 1/a) + 1/a] / (e^(aL) - 1),
//       a = k + 1,  L = ln(pmax/pmin),
//     about 0.66.  The self-normalised weighted mean of 40000 draws has a
//     standard error of about 0.66 % of that value, so the 2 % tolerance is
//     about 3 standard errors.
void TestSource() {
  UniformSphere shock;
  Mesh::ObserverFluxTube tube;
  for (int k = 0; k < 3; ++k) {
    Mesh::TracedTubeLine line;
    line.fluxWb = 1.0e9 * (k + 1);
    line.sector = 1;
    const Vec3 e = Vec3(1.0, 0.01 * k, 0.0).Normalized();
    for (int i = 0; i <= 400; ++i) {
      const double radius = 1.0e9 * std::pow(300.0, i / 400.0);
      line.pointsM.push_back(e * radius);
      line.localSpacingM.push_back(std::sqrt(line.fluxWb / (shock.B * std::pow(3.0e9 / radius, 2))));
    }
    tube.lines.push_back(line);
  }
  double rate = 0;
  Adapters::TubeLineNetRate(tube.lines[1], shock, 0.0, 1e-3, &rate);
  const double expected = 1e-3 * shock.n * (shock.V - shock.U) * 2.0e9 / shock.B;
  Check(Near(rate, expected, 1e-12), "E3D06a line rate = eta n1 w1 Phi/|B.n| (1-D NetProtonRate)");
  Adapters::TubeLineSourceRequest q;
  q.tube = &tube;
  q.shock = &shock;
  q.plasma = [&](const Vec3& x, double, Vec3* b, Vec3* u) {
    const Vec3 e = x.Normalized();
    *b = e * (shock.B * std::pow(3.0e9 / x.Norm(), 2));
    *u = e * shock.U;
    return Core::Status::OK();
  };
  q.species = 0;
  q.massKg = Core::Const::m_p;
  q.injectionEfficiency = 1e-3;
  q.releaseFraction = 0.3;
  q.minimumEnergyJ = 0.02e6 * Core::Const::e;
  q.maximumEnergyJ = 200e6 * Core::Const::e;
  q.logUniformMomentum = true;
  q.focusedPitchLaw = true;
  q.samplesPerStep = 40000;
  q.campaignSeed = 11;
  q.step = 3;
  q.stepBeginS = 0.0;
  q.dtS = 1.0;   // short step: the rate is constant in this test
  const auto result = Adapters::GenerateTubeLineSource(q);
  double total = 0;
  for (const auto& l : tube.lines) {
    double r = 0;
    Adapters::TubeLineNetRate(l, shock, 0.0, 1e-3, &r);
    total += r;
  }
  Check(result.status.ok() && result.particles.size() == q.samplesPerStep,
        "E3D06b one candidate per stratified sample");
  Check(Near(result.representedParticles, total * q.dtS, 0.03),
        "E3D06c log-uniform importance weights are unbiased");
  bool onRelease = true, outward = true;
  for (const auto& injected : result.particles) {
    const double t = (1.0 - injected.remainingFirstStepFraction) * q.dtS;
    onRelease = onRelease && Near(injected.particle.positionM.Norm(), 1.3 * shock.Radius(t), 1e-9);
    const double speed = Core::Const::c * injected.particle.momentumKgMPerS /
        std::hypot(injected.particle.momentumKgMPerS, Core::Const::m_p * Core::Const::c);
    outward = outward && speed * injected.particle.mu + shock.U - shock.V > -1e-9 * speed;
  }
  Check(onRelease, "E3D06d particles are released on the (1+f) R_sh sphere");
  Check(outward, "E3D06e focused pitch law gives positive normal flux through the release surface");
  // The weighted DSA spectrum: <ln p> of dN/dp ~ p^(2-q), q = 3r/(r-1).
  const double qd = 3 * shock.r / (shock.r - 1), k = 2 - qd;
  const double pmin = std::sqrt(q.minimumEnergyJ * (q.minimumEnergyJ + 2 * Core::Const::m_p * Core::Const::c * Core::Const::c)) / Core::Const::c;
  const double pmax = std::sqrt(q.maximumEnergyJ * (q.maximumEnergyJ + 2 * Core::Const::m_p * Core::Const::c * Core::Const::c)) / Core::Const::c;
  double sw = 0, swl = 0;
  for (const auto& injected : result.particles) {
    sw += injected.particle.statisticalWeight;
    swl += injected.particle.statisticalWeight * std::log(injected.particle.momentumKgMPerS / pmin);
  }
  // Analytic mean of ln(p/pmin) for the truncated power law of index k.
  const double L = std::log(pmax / pmin), a = k + 1;
  const double mean = (std::exp(a * L) * (L - 1 / a) + 1 / a) / (std::exp(a * L) - 1);
  Check(Near(swl / sw, mean, 0.02), "E3D06f weighted spectrum follows the local DSA law");
}

}  // namespace

int main() {
  TestFrame();
  TestInertialForce();
  TestEphemeris();
  TestTube();
  TestShockReturn();
  TestSource();
  std::printf("%s: %d failure(s)\n", failures ? "FAILED" : "PASSED", failures);
  return failures ? 1 : 0;
}
