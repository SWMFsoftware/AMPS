// ============================================================================
// srcSEP3D/mesh/observer_flux_tube.cpp
//
// Field-line tracing and construction of the observer flux tube.  The
// geometry, the formulas and the frame conventions are described in
// observer_flux_tube.h; the comments here explain the individual steps.
// ============================================================================

#include "observer_flux_tube.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>
#include <set>
#include <utility>

namespace SEP3D {
namespace Mesh {
namespace {

Core::Status Invalid(const std::string& message) {
  return Core::Status(Core::StatusCode::InvalidInput, message);
}

// The unit tangent sign * sector * B/|B| at x: the right-hand side of the
// field-line ODE dx/ds.
//   sector  in/out: 0 on the first call of a trace, which sets it to the
//           sector at x.  Afterwards every evaluation must return the same
//           sector, because an open line never changes polarity; a change
//           means that the trace has crossed a sector boundary (the
//           heliospheric current sheet) or a closed region.
// Errors: closed field at x, a change of sector, or a zero or nonfinite |B|.
// The field's own failure status is returned unchanged.
Core::Status Direction(const FrameFieldFunction& field, const Core::Vec3& x,
                       int sign, int* sector, Core::Vec3* direction) {
  Core::Vec3 b;
  int s = 0;
  const Core::Status status = field(x, &b, &s);
  if (!status.ok()) return status;
  const double magnitude = b.Norm();
  if (s == 0)
    return Invalid("field line reaches a closed-field region");
  if (*sector != 0 && s != *sector)
    return Invalid("field line changes open-field sector");
  if (!(magnitude > 0.0) || !std::isfinite(magnitude))
    return Invalid("field line reaches a zero or nonfinite magnetic field");
  *sector = s;
  *direction = b * (static_cast<double>(sign * s) / magnitude);
  return Core::Status::OK();
}

// One classical fourth-order Runge-Kutta step of arc length ds [m] for
// dx/ds = Direction(x):
//   k1 = d(x),  k2 = d(x + k1 ds/2),  k3 = d(x + k2 ds/2),  k4 = d(x + k3 ds),
//   x_next = x + (k1 + 2 k2 + 2 k3 + k4) ds / 6.
// All four stages must lie on open field of the same sector.  The ODE is in
// arc length, so its right-hand side is a unit vector and the step length is
// independent of |B|.
Core::Status Rk4Step(const FrameFieldFunction& field, const Core::Vec3& x,
                     double ds, int sign, int* sector, Core::Vec3* next) {
  Core::Vec3 k1, k2, k3, k4;
  Core::Status status = Direction(field, x, sign, sector, &k1);
  if (!status.ok()) return status;
  status = Direction(field, x + k1 * (0.5 * ds), sign, sector, &k2);
  if (!status.ok()) return status;
  status = Direction(field, x + k2 * (0.5 * ds), sign, sector, &k3);
  if (!status.ok()) return status;
  status = Direction(field, x + k3 * ds, sign, sector, &k4);
  if (!status.ok()) return status;
  *next = x + (k1 + k2 * 2.0 + k3 * 2.0 + k4) * (ds / 6.0);
  return Core::Status::OK();
}

// One FNV-1a pass over `size` bytes: for each byte, hash ^= byte and then
// hash *= 1099511628211 (the 64-bit FNV prime), modulo 2^64.  The digest is
// only compared for equality (across MPI ranks and in the run log), so the
// raw bytes of the doubles are hashed as they are.
std::uint64_t Fnv(std::uint64_t hash, const void* data, std::size_t size) {
  const unsigned char* bytes = static_cast<const unsigned char*>(data);
  for (std::size_t i = 0; i < size; ++i) {
    hash ^= bytes[i];
    hash *= 1099511628211ULL;
  }
  return hash;
}

// Angle [rad] between two nonzero vectors, acos of the normalized dot
// product.  The cosine is clamped to [-1, 1] against roundoff.  This is used
// for the candidate window of step 2, whose margin is a whole lattice
// spacing or more, so acos's loss of accuracy at nearly parallel vectors is
// immaterial there.
double AngleBetween(const Core::Vec3& a, const Core::Vec3& b) {
  const double c = a.Dot(b) / (a.Norm() * b.Norm());
  return std::acos(std::max(-1.0, std::min(1.0, c)));
}

}  // namespace

Core::Status TraceFieldLineToRadius(const FrameFieldFunction& field,
                                    const Core::Vec3& start, int sign,
                                    double targetRadiusM, double traceStepM,
                                    std::size_t maximumSteps,
                                    std::vector<Core::Vec3>* points) {
  if (!points || (sign != 1 && sign != -1) || !std::isfinite(targetRadiusM) ||
      targetRadiusM <= 0.0 || !std::isfinite(traceStepM) || traceStepM <= 0.0)
    return Invalid("invalid field-line trace request");
  points->assign(1, start);
  Core::Vec3 x = start;
  double r = x.Norm();
  const bool outward = sign > 0;
  // Already at or beyond the target in the requested direction: the line
  // is the single point `start`.
  if (outward ? r >= targetRadiusM : r <= targetRadiusM)
    return Core::Status::OK();
  int sector = 0;   // set by the first Direction() call, then held fixed
  for (std::size_t step = 0; step < maximumSteps; ++step) {
    // The arc step is at most 2 % of the local heliocentric radius, so the
    // step shrinks towards the Sun, where the field geometry varies on
    // smaller length scales.  trace_step_m bounds it at large r.
    const double ds = std::min(traceStepM, 0.02 * r);
    Core::Vec3 next;
    Core::Status status = Rk4Step(field, x, ds, sign, &sector, &next);
    if (!status.ok()) return status;
    const double rn = next.Norm();
    // Every step must move strictly in the requested radial direction.  This
    // admits only lines that are monotone in r (the Parker regime), so each
    // radius occurs once on the line (PointAtRadius, the keep test of
    // BuildObserverFluxTube).  A line that turns back, e.g. a near-closed
    // loop, is an error rather than a truncated line.
    if (outward ? !(rn > r) : !(rn < r))
      return Invalid("field line is not monotone in heliocentric radius at r = " +
                     std::to_string(r) + " m");
    const bool crossed = outward ? rn >= targetRadiusM : rn <= targetRadiusM;
    if (!crossed) {
      x = next;
      r = rn;
      points->push_back(x);
      continue;
    }
    // This step crosses the target sphere.  Bisect the length h in [0, ds]
    // of the last step: each trial is a full RK4 step of length h from x,
    // so the landing point stays on the RK4 solution, and `landing` holds the
    // shortest trial found so far that reaches the sphere.  60 halvings take
    // the bracket width to ds * 2^-60, below double-precision resolution.
    double low = 0.0, high = ds;
    Core::Vec3 landing = next;
    for (int iteration = 0; iteration < 60; ++iteration) {
      const double middle = 0.5 * (low + high);
      Core::Vec3 trial;
      status = Rk4Step(field, x, middle, sign, &sector, &trial);
      if (!status.ok()) return status;
      const double rt = trial.Norm();
      const bool beyond = outward ? rt >= targetRadiusM : rt <= targetRadiusM;
      if (beyond) {
        high = middle;
        landing = trial;
      } else {
        low = middle;
      }
    }
    // Remove the residual radial error exactly (it is below the bisection
    // width); the landing point is then on the sphere to roundoff.  This
    // rescales an already converged point by a factor 1 + O(2^-60); it never
    // replaces a point that the trace failed to reach.
    landing = landing * (targetRadiusM / landing.Norm());
    points->push_back(landing);
    return Core::Status::OK();
  }
  // The budget ran out before the sphere.  This is an error, so a field
  // that winds without making radial progress cannot hang the build.
  return Invalid("field line did not reach r = " + std::to_string(targetRadiusM) +
                 " m within the step budget");
}

bool PointAtRadius(const TracedTubeLine& line, double radiusM,
                   Core::Vec3* pointM, double* localSpacingM) {
  const auto& p = line.pointsM;
  if (p.size() < 2) return false;
  if (radiusM < p.front().Norm() || radiusM > p.back().Norm()) return false;
  // Bisection on the strictly increasing point radii: at exit
  // |p[low]| <= r <= |p[high]| with high = low + 1.
  std::size_t low = 0, high = p.size() - 1;
  while (high - low > 1) {
    const std::size_t middle = (low + high) / 2;
    if (p[middle].Norm() <= radiusM) low = middle; else high = middle;
  }
  // The fraction is linear in radius, and the returned point lies on the
  // chord between the bracketing points.  |x| is convex along a chord, so
  // the point's radius is at most r.  The deficit is second order in the
  // chord's transverse extent, which the RK4 step bound keeps small.
  const double r0 = p[low].Norm(), r1 = p[high].Norm();
  const double f = r1 > r0 ? (radiusM - r0) / (r1 - r0) : 0.0;
  if (pointM) *pointM = p[low] + (p[high] - p[low]) * f;
  if (localSpacingM)
    *localSpacingM = line.localSpacingM[low] +
        f * (line.localSpacingM[high] - line.localSpacingM[low]);
  return true;
}

Core::Status BuildObserverFluxTube(const std::vector<TubeObserverPath>& paths,
                                   const ObserverTubeControls& controls,
                                   const FrameFieldFunction& field,
                                   ObserverFluxTube* tube) {
  // ---- Input validation --------------------------------------------------
  if (!tube || paths.empty() || !field)
    return Invalid("observer flux tube needs observers, a field and an output");
  if (!(controls.innerRadiusM > 0.0) ||
      !(controls.outerRadiusM > controls.innerRadiusM) ||
      !(controls.traceStepM > 0.0) || !(controls.seedSpacingM > 0.0))
    return Invalid("observer flux tube controls are invalid");
  for (const TubeObserverPath& path : paths) {
    if (path.centresM.empty() || path.centresM.size() != path.timesS.size() ||
        !(path.collectionRadiusM > 0.0))
      return Invalid("observer '" + path.id + "' has no tube samples");
    for (const Core::Vec3& c : path.centresM) {
      // The whole collection sphere [r - rho, r + rho] must lie strictly
      // inside the transport shell, so every traced line covers the
      // sample's radius and the sphere never touches a boundary.
      const double r = c.Norm();
      if (!(r - path.collectionRadiusM > controls.innerRadiusM) ||
          !(r + path.collectionRadiusM < controls.outerRadiusM))
        return Invalid("observer '" + path.id + "' collection sphere leaves "
                       "the radial transport shell");
    }
  }
  // Build into a local result.  *tube is written only after every step has
  // succeeded, so a failure leaves the caller's tube untouched.
  ObserverFluxTube result;
  // ---- 1. Reference sphere and seed lattice --------------------------------
  // R0 is the radius of the first tube observer at the start of its stay.
  // The caller's default s (half the requested cell size) is evaluated at
  // the same observer position (main_lib.cpp).
  const double R0 = paths.front().centresM.front().Norm();
  const double s = controls.seedSpacingM;
  result.referenceRadiusM = R0;
  result.seedSpacingM = s;

  // Colatitude rings: the nominal width s/R0 is rounded so that an integer
  // number of rings tiles [0, pi] exactly (ringWidth <= s/R0).  Ring i covers
  // theta in [i, i + 1] * ringWidth.  ringCells(i) = round(2 pi R0
  // sin(theta_mid) / s) longitude cells of equal width 2 pi / n_i (at least
  // 1), so every cell is about s x s and the polar rings are not
  // oversubdivided.  A cell is addressed by (ring i, longitude index j), and
  // its centre is (theta, phi) = ((i + 0.5) ringWidth, (j + 0.5) 2 pi / n_i).
  const double dTheta = s / R0;
  const int rings = std::max(1, static_cast<int>(std::ceil(Core::Const::kPi / dTheta)));
  const double ringWidth = Core::Const::kPi / rings;
  auto ringCells = [&](int i) {
    const double theta = (i + 0.5) * ringWidth;
    return std::max(1, static_cast<int>(std::lround(
        2.0 * Core::Const::kPi * R0 * std::sin(theta) / s)));
  };

  // ---- 2. Candidate cells around the R0 image of every observer sample ----
  // A std::set removes cells found from several samples, and it iterates in
  // sorted (i, j) order.  That makes the kept lines, their order and the
  // digest independent of the order of the samples' discovery and identical
  // on every rank.
  std::set<std::pair<int, int>> candidates;
  for (const TubeObserverPath& path : paths) {
    for (const Core::Vec3& c : path.centresM) {
      ++result.observerSamples;
      // q: where the field line through the sample centre c crosses R0.  A
      // sample within 1 m of R0 is used as it is; any other sample is traced
      // outward (when R0 > |c|) or sunward (when R0 < |c|) to R0.
      Core::Vec3 q = c;
      const double rc = c.Norm();
      if (std::fabs(rc - R0) > 1.0) {
        std::vector<Core::Vec3> leg;
        const Core::Status traced = TraceFieldLineToRadius(
            field, c, R0 > rc ? +1 : -1, R0, controls.traceStepM,
            controls.maximumStepsPerLine, &leg);
        if (!traced.ok())
          return Invalid("observer '" + path.id + "' centre line: " +
                         traced.message);
        q = leg.back();
      }
      // Image of the collection sphere on R0.  By flux conservation an area
      // A at c maps to A |B(c)|/|B(q)| at q, so a transverse length rho
      // maps to rho * m with m = sqrt(|B(c)|/|B(q)|).  (The sector output
      // is not used here; the traces below check it.)
      Core::Vec3 bc, bq;
      int sector = 0;
      Core::Status status = field(c, &bc, &sector);
      if (!status.ok()) return status;
      status = field(q, &bq, &sector);
      if (!status.ok()) return status;
      const double mapping = std::sqrt(bc.Norm() / bq.Norm());
      // Half-angle of the candidate window around q.  The 2 s margin is
      // larger than the image of the keep test's s_k term (at most about s
      // on R0, since Phi ~ |B_r| s^2 <= |B| s^2 there) plus the half-diagonal
      // of a lattice cell (about 0.7 s).  m is one isotropic factor; see
      // ASSUMPTIONS AND LIMITS in the header.
      const double alpha = (path.collectionRadiusM * mapping + 2.0 * s) / R0;
      const double thetaQ = std::acos(std::max(-1.0, std::min(1.0, q.z / q.Norm())));
      // Rings that can hold a cell centre within alpha of q, plus one ring
      // on each side.  Every cell of those rings is then tested by its exact
      // great-circle angle to q (near a pole the window spans all
      // longitudes, which the full ring scan handles).
      const int iLow = std::max(0, static_cast<int>(std::floor((thetaQ - alpha) / ringWidth)) - 1);
      const int iHigh = std::min(rings - 1, static_cast<int>(std::ceil((thetaQ + alpha) / ringWidth)) + 1);
      for (int i = iLow; i <= iHigh; ++i) {
        const int n = ringCells(i);
        const double theta = (i + 0.5) * ringWidth;
        for (int j = 0; j < n; ++j) {
          const double phi = (j + 0.5) * 2.0 * Core::Const::kPi / n;
          const Core::Vec3 p(std::sin(theta) * std::cos(phi),
                             std::sin(theta) * std::sin(phi), std::cos(theta));
          if (AngleBetween(p, q) <= alpha) candidates.insert({i, j});
        }
      }
    }
  }
  result.candidateLines = candidates.size();

  // ---- 3. Trace each candidate and keep it when it meets a sphere --------
  for (const auto& cell : candidates) {
    const int i = cell.first, j = cell.second, n = ringCells(i);
    const double theta = (i + 0.5) * ringWidth;
    const double phi = (j + 0.5) * 2.0 * Core::Const::kPi / n;
    const Core::Vec3 unit(std::sin(theta) * std::cos(phi),
                          std::sin(theta) * std::sin(phi), std::cos(theta));
    const Core::Vec3 p = unit * R0;   // seed point: the cell centre on R0
    Core::Vec3 b;
    int sector = 0;
    Core::Status status = field(p, &b, &sector);
    if (!status.ok()) return status;
    // A seed on closed field is skipped: a closed loop cannot connect the
    // inner sphere to an observer.  A candidate that starts on open field
    // but enters closed field along its trace stops the build (below).
    if (sector == 0) continue;
    // The exact area of the cell (spherical zone, divided into n equal
    // longitude cells) times |B_r| at its centre: the flux element Phi.
    const double area = R0 * R0 *
        (std::cos(i * ringWidth) - std::cos((i + 1) * ringWidth)) *
        (2.0 * Core::Const::kPi / n);
    const double flux = std::fabs(b.Dot(unit)) * area;
    // A cell with no radial flux (B tangent to R0) carries no element.
    if (!(flux > 0.0)) continue;
    // Both legs start at the seed.  Each leg fixes its own sector at p, and
    // both legs see the same sector because they start at the same point.
    std::vector<Core::Vec3> sunward, outward;
    status = TraceFieldLineToRadius(field, p, -1, controls.innerRadiusM,
                                    controls.traceStepM,
                                    controls.maximumStepsPerLine, &sunward);
    if (!status.ok())
      return Invalid("tube candidate line (sunward): " + status.message);
    status = TraceFieldLineToRadius(field, p, +1, controls.outerRadiusM,
                                    controls.traceStepM,
                                    controls.maximumStepsPerLine, &outward);
    if (!status.ok())
      return Invalid("tube candidate line (outward): " + status.message);
    // Join the legs into one inner-to-outer line: the reversed sunward leg
    // (inner sphere ... p) and then the outward leg without its first point
    // (p again).  Both legs are strictly monotone in r, so the joined line is
    // monotone too.
    TracedTubeLine line;
    line.fluxWb = flux;
    line.sector = sector;
    line.referencePointM = p;
    line.pointsM.assign(sunward.rbegin(), sunward.rend());
    line.pointsM.insert(line.pointsM.end(), outward.begin() + 1, outward.end());
    // Local transverse size of the element, s(x) = sqrt(Phi/|B(x)|): with
    // the cross-section A(x) = Phi/|B(x)|, s is the side of a square of
    // that area.  Direction() along the trace guaranteed |B| > 0 at the RK4
    // stages; the field is evaluated again at the stored points here.
    line.localSpacingM.reserve(line.pointsM.size());
    for (const Core::Vec3& x : line.pointsM) {
      Core::Vec3 bx;
      int sx = 0;
      status = field(x, &bx, &sx);
      if (!status.ok()) return status;
      line.localSpacingM.push_back(std::sqrt(flux / bx.Norm()));
    }
    // ---- Keep test ---------------------------------------------------------
    // The line is kept when, at the radius of some observer sample c, its
    // point x(|c|) satisfies |x - c| <= rho + s(x): the element of
    // transverse size s around the line then reaches the collection sphere.
    // The separation is measured at equal radius.  It is never smaller than
    // the minimum distance from c to the line.  When c is offset from the
    // line within the plane that holds the line and the radial direction,
    // it is larger by up to 1/cos(psi), psi being the line's angle to the
    // radial direction (about 45 degrees for the Parker spiral at 1 au).
    bool keep = false;
    for (const TubeObserverPath& path : paths) {
      for (const Core::Vec3& c : path.centresM) {
        Core::Vec3 x;
        double spacing = 0.0;
        if (!PointAtRadius(line, c.Norm(), &x, &spacing)) continue;
        if ((x - c).Norm() <= path.collectionRadiusM + spacing) {
          keep = true;
          break;
        }
      }
      if (keep) break;
    }
    if (keep) result.lines.push_back(std::move(line));
  }
  if (result.lines.empty())
    return Invalid("observer flux tube contains no field line; check the "
                   "observer trajectories and seed_spacing_m");
  // ---- Identity digest -----------------------------------------------------
  // FNV-1a over the raw bytes of every point of every kept line, in line
  // order.  The seed 1469598103934665603 is the value used by the other
  // digests of this code base (e.g. background/bg_parker.cpp).  It is not
  // the published FNV-1a 64-bit offset basis 14695981039346656037.  Any
  // fixed seed works here, because the digest is only compared for
  // equality (across MPI ranks in main_lib.cpp, and in the run log).
  std::uint64_t hash = 1469598103934665603ULL;
  for (const TracedTubeLine& line : result.lines)
    for (const Core::Vec3& x : line.pointsM) {
      const double v[3] = {x.x, x.y, x.z};
      hash = Fnv(hash, v, sizeof(v));
    }
  result.digest = hash;
  *tube = std::move(result);
  return Core::Status::OK();
}

}  // namespace Mesh
}  // namespace SEP3D
