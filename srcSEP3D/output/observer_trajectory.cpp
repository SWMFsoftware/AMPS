// ============================================================================
// srcSEP3D/output/observer_trajectory.cpp
//
// The four trajectory implementations (frame-fixed, linear HCI, ephemeris
// table, SPICE), the ephemeris-table reader and its cubic Hermite
// interpolation.  Kinds, frames and units: observer_trajectory.h.
// ============================================================================

#include "observer_trajectory.h"

#include "spice_ephemeris.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>

// SHA-256 of the shared coronal-CME library (linked into every srcSEP3D
// build; declared in its private src/sha256.h).  It returns the 64-digit
// lower-case hex digest.  The configuration lower-cases the deck's
// ephemeris_sha256 (runtime/configuration_io.cpp), so the two compare
// directly.  The declaration is repeated here because that header is
// private to the library.
namespace SEP { namespace CoronalCME { namespace Internal {
std::string Sha256Hex(const std::string& bytes);
} } }

namespace SEP3D {
namespace Output {
namespace {

Core::Status Invalid(const std::string& message) {
  return Core::Status(Core::StatusCode::InvalidInput, message);
}

// An observer fixed in the computation frame (fixed-cartesian,
// fixed-heliographic, field-connected, spherical-shell).  In a corotating
// run it stays on one field line, like a 1-D srcSEP observer.  Its support
// is unbounded.
class FrameFixedTrajectory final : public ObserverTrajectory {
 public:
  explicit FrameFixedTrajectory(const Core::Vec3& position) : x_(position) {}
  Core::Status StateHci(double, Core::Vec3* x, Core::Vec3* v) const override {
    // Only the frame position is meaningful; ObserverFramePosition uses it.
    if (x) *x = x_;
    if (v) *v = Core::Vec3();
    return Core::Status::OK();
  }
  double BeginS() const override { return -1e300; }
  double EndS() const override { return 1e300; }
  bool FrameFixed() const override { return true; }
  std::string Description() const override { return "frame-fixed"; }

 private:
  Core::Vec3 x_;
};

// moving-cartesian: x_hci(t) = x0 + v t in the inertial HCI frame (t model
// time [s], x0 the position at t = 0 [m], v [m/s]); a point plus the
// velocity with which it moves.  Unbounded support.
class LinearTrajectory final : public ObserverTrajectory {
 public:
  LinearTrajectory(const Core::Vec3& x0, const Core::Vec3& v) : x0_(x0), v_(v) {}
  Core::Status StateHci(double t, Core::Vec3* x, Core::Vec3* v) const override {
    if (x) *x = x0_ + v_ * t;
    if (v) *v = v_;
    return Core::Status::OK();
  }
  double BeginS() const override { return -1e300; }
  double EndS() const override { return 1e300; }
  bool FrameFixed() const override { return false; }
  std::string Description() const override { return "moving-cartesian (HCI)"; }

 private:
  Core::Vec3 x0_, v_;
};

// body with ephemeris = table: the parsed, checksum-verified samples (at
// least two, strictly increasing times) and Hermite interpolation between
// them.  The support is the table's time range, and a query outside it is
// NotFound (no extrapolation).
class TableTrajectory final : public ObserverTrajectory {
 public:
  TableTrajectory(std::vector<EphemerisSample> samples, std::string origin)
      : samples_(std::move(samples)), origin_(std::move(origin)) {}
  Core::Status StateHci(double t, Core::Vec3* x, Core::Vec3* v) const override {
    return InterpolateEphemeris(samples_, t, x, v);
  }
  double BeginS() const override { return samples_.front().timeS; }
  double EndS() const override { return samples_.back().timeS; }
  bool FrameFixed() const override { return false; }
  std::string Description() const override { return "ephemeris table " + origin_; }

 private:
  std::vector<EphemerisSample> samples_;
  std::string origin_;
};

// body with ephemeris = spice: every query asks CSPICE for the state at
// et = referenceEt + t (spice_ephemeris.h) and rotates it into the model
// HCI frame.  referenceEt [s] is str2et(spice_reference_utc), computed once
// by PrepareSpiceEphemeris when the trajectory is created.
class SpiceTrajectory final : public ObserverTrajectory {
 public:
  SpiceTrajectory(SpiceEphemerisRequest request, double referenceEt)
      : request_(std::move(request)), referenceEt_(referenceEt) {}
  Core::Status StateHci(double t, Core::Vec3* x, Core::Vec3* v) const override {
    return SpiceEphemerisState(request_, referenceEt_, t, x, v);
  }
  // SPICE coverage is checked by every query; the support is the kernels'.
  double BeginS() const override { return -1e300; }
  double EndS() const override { return 1e300; }
  bool FrameFixed() const override { return false; }
  std::string Description() const override {
    return "SPICE " + request_.target + " wrt " + request_.center + " in " +
           request_.frame;
  }

 private:
  SpiceEphemerisRequest request_;
  double referenceEt_ = 0.0;
};

}  // namespace

Core::Status ParseEphemerisTable(const std::string& text,
                                 const std::string& origin,
                                 std::vector<EphemerisSample>* samples) {
  if (!samples) return Invalid("ephemeris output is null");
  samples->clear();
  std::istringstream lines(text);
  std::string line;
  bool header = false;
  std::size_t number = 0;
  while (std::getline(lines, line)) {
    ++number;
    if (!line.empty() && line.back() == '\r') line.pop_back();
    if (line.empty() || line[0] == '#') continue;
    if (!header) {
      if (line != "time_s,x_m,y_m,z_m,vx_m_per_s,vy_m_per_s,vz_m_per_s")
        return Invalid(origin + ":" + std::to_string(number) +
                       ": ephemeris header must be time_s,x_m,y_m,z_m,"
                       "vx_m_per_s,vy_m_per_s,vz_m_per_s");
      header = true;
      continue;
    }
    // Split on ','.  Up to 8 cells are read, so a row with too many columns
    // is detected (count == 8) without reading the rest of the line.  Each
    // cell must be one complete finite number (strtod must consume it all).
    std::stringstream row(line);
    std::string cell;
    double value[7];
    int count = 0;
    while (std::getline(row, cell, ',') && count < 8) {
      char* end = nullptr;
      const double parsed = std::strtod(cell.c_str(), &end);
      if (end == cell.c_str() || *end != '\0' || !std::isfinite(parsed))
        return Invalid(origin + ":" + std::to_string(number) +
                       ": ephemeris value '" + cell + "' is not a finite number");
      if (count < 7) value[count] = parsed;
      ++count;
    }
    if (count != 7)
      return Invalid(origin + ":" + std::to_string(number) +
                     ": ephemeris row must have 7 columns");
    EphemerisSample sample;
    sample.timeS = value[0];
    sample.positionM = {value[1], value[2], value[3]};
    sample.velocityMPerS = {value[4], value[5], value[6]};
    if (!samples->empty() && !(sample.timeS > samples->back().timeS))
      return Invalid(origin + ":" + std::to_string(number) +
                     ": ephemeris times must strictly increase");
    samples->push_back(sample);
  }
  if (samples->size() < 2)
    return Invalid(origin + ": ephemeris table needs at least two samples");
  return Core::Status::OK();
}

Core::Status InterpolateEphemeris(const std::vector<EphemerisSample>& samples,
                                  double t, Core::Vec3* positionM,
                                  Core::Vec3* velocityMPerS) {
  if (samples.size() < 2 || !std::isfinite(t))
    return Invalid("invalid ephemeris interpolation request");
  if (t < samples.front().timeS || t > samples.back().timeS)
    return Core::Status(Core::StatusCode::NotFound,
                        "model time " + std::to_string(t) +
                            " s is outside the ephemeris table [" +
                            std::to_string(samples.front().timeS) + ", " +
                            std::to_string(samples.back().timeS) + "] s");
  // Bracketing interval [a, b]: b is the first sample with time > t.  At
  // t = last time there is no such sample, and the last interval is used
  // (u = 1).  At t = first time, b is the second sample (u = 0).  So a is
  // never before the first sample.
  auto upper = std::upper_bound(
      samples.begin(), samples.end(), t,
      [](double value, const EphemerisSample& s) { return value < s.timeS; });
  if (upper == samples.end()) --upper;
  const EphemerisSample& b = *upper;
  const EphemerisSample& a = *(upper - 1);
  const double h = b.timeS - a.timeS, u = (t - a.timeS) / h;
  // Cubic Hermite basis on u in [0, 1] (h00, h10, h01, h11) and the basis
  // derivatives d/du (d00, d10, d01, d11):
  //   x(u)   = h00 x_a + h10 h v_a + h01 x_b + h11 h v_b,
  //   dx/dt  = (1/h) dx/du = (d00 x_a + d01 x_b)/h + d10 v_a + d11 v_b.
  // The tabulated velocities enter as slopes dx/du = h v, so x(t) matches
  // both the positions and the velocities at the two samples.
  const double h00 = (1 + 2 * u) * (1 - u) * (1 - u), h10 = u * (1 - u) * (1 - u);
  const double h01 = u * u * (3 - 2 * u), h11 = u * u * (u - 1);
  const double d00 = 6 * u * u - 6 * u, d10 = 3 * u * u - 4 * u + 1;
  const double d01 = -6 * u * u + 6 * u, d11 = 3 * u * u - 2 * u;
  if (positionM)
    *positionM = a.positionM * h00 + a.velocityMPerS * (h10 * h) +
                 b.positionM * h01 + b.velocityMPerS * (h11 * h);
  if (velocityMPerS)
    *velocityMPerS = (a.positionM * d00 + b.positionM * d01) / h +
                     a.velocityMPerS * d10 + b.velocityMPerS * d11;
  return Core::Status::OK();
}

Core::Status CreateObserverTrajectory(
    const RuntimeModel::ObserverOptions& observer,
    ObserverTrajectoryPointer* output) {
  using RuntimeModel::ObserverKind;
  using RuntimeModel::ObserverEphemerisSource;
  if (!output) return Invalid("observer trajectory output is null");
  // A point and a velocity: the straight HCI line x0 + v t.
  if (observer.kind == ObserverKind::MovingCartesian) {
    *output = std::make_shared<LinearTrajectory>(observer.positionM,
                                                 observer.velocityMPerS);
    return Core::Status::OK();
  }
  // Every other non-body kind is fixed in the computation frame at
  // observer.positionM as configured.  (For spherical-shell,
  // output/observer_runtime.cpp places the sphere at shellRadiusM along the
  // direction of positionM; this trajectory does not apply that.)
  if (observer.kind != ObserverKind::Body) {
    *output = std::make_shared<FrameFixedTrajectory>(observer.positionM);
    return Core::Status::OK();
  }
  // Ephemeris table.  The checksum is computed over the exact bytes that are
  // then parsed (binary read, no newline translation), so the table that is
  // used is the one the deck pins with ephemeris_sha256.
  if (observer.ephemerisSource == ObserverEphemerisSource::Table) {
    std::ifstream input(observer.ephemerisResolvedPath, std::ios::binary);
    if (!input)
      return Invalid("cannot read ephemeris table '" +
                     observer.ephemerisResolvedPath + "' of observer '" +
                     observer.id + "'");
    std::ostringstream bytes;
    bytes << input.rdbuf();
    const std::string text = bytes.str();
    const std::string digest = SEP::CoronalCME::Internal::Sha256Hex(text);
    if (digest != observer.ephemerisSha256)
      return Invalid("ephemeris table '" + observer.ephemerisResolvedPath +
                     "' has SHA-256 " + digest + ", the deck declares " +
                     observer.ephemerisSha256);
    std::vector<EphemerisSample> samples;
    const Core::Status parsed =
        ParseEphemerisTable(text, observer.ephemerisResolvedPath, &samples);
    if (!parsed.ok()) return parsed;
    *output = std::make_shared<TableTrajectory>(std::move(samples),
                                                observer.ephemerisResolvedPath);
    return Core::Status::OK();
  }
  // SPICE.  The deck's rotation is already validated as a proper orthonormal
  // matrix (runtime/run_configuration.cpp).  The kernels are loaded and the
  // reference UTC converted here, so a missing kernel or a build without
  // SPICE fails when the trajectory is created.
  if (observer.ephemerisSource == ObserverEphemerisSource::Spice) {
    SpiceEphemerisRequest request;
    request.target = observer.spiceTarget;
    request.center = observer.spiceCenter;
    request.frame = observer.spiceFrame;
    request.kernels = observer.spiceKernels;
    request.referenceUtc = observer.spiceReferenceUtc;
    request.toModelRotation = observer.spiceToModelRotation;
    double referenceEt = 0.0;
    const Core::Status prepared = PrepareSpiceEphemeris(request, &referenceEt);
    if (!prepared.ok())
      return Core::Status(prepared.code, "observer '" + observer.id +
                                             "': " + prepared.message);
    *output = std::make_shared<SpiceTrajectory>(std::move(request), referenceEt);
    return Core::Status::OK();
  }
  return Invalid("body observer '" + observer.id + "' has no ephemeris");
}

Core::Status ObserverFramePosition(const ObserverTrajectory& trajectory,
                                   const Core::ComputationFrame& frame,
                                   double epochS, Core::Vec3* positionM) {
  if (!positionM) return Invalid("observer position output is null");
  // A frame-fixed trajectory already returns its frame position.  An
  // inertial one returns x_hci(t), which is rotated into the frame:
  // x_frame = Rz(-Omega (t - t_ref)) x_hci (the identity for an HCI run).
  Core::Vec3 x;
  const Core::Status state = trajectory.StateHci(epochS, &x, nullptr);
  if (!state.ok()) return state;
  *positionM = trajectory.FrameFixed() ? x : frame.FromInertialPosition(x, epochS);
  return Core::Status::OK();
}

}  // namespace Output
}  // namespace SEP3D
