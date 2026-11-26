// ============================================================================
// srcSEP3D/core/computation_frame.h
//
// Computation frame of a srcSEP3D run.
//
// srcSEP3D solves particle transport in one of two Cartesian frames, both
// centred on the Sun:
//
//   hci         the heliocentric inertial frame of the background provider;
//   corotating  the frame that rotates rigidly with the Sun about the +z axis
//               of HCI at the angular rate Omega.  The two frames coincide at
//               the reference epoch t_ref.
//
// SELECTION
//   mesh.domain.coordinate_frame              HCI | HCI-corotating
//   mesh.domain.frame_rotation_rate_rad_per_s Omega [rad/s] (corotating only)
//   mesh.domain.frame_reference_epoch_s       t_ref [model elapsed s]
// (srcSEP3D/CONFIGURATION.md, "Corotating frame, observer flux tube, ...").
// RunOptions::ComputationFrameValue() (runtime/run_configuration.*) builds a
// ComputationFrame from these keys; Canonical() enters the run's physics
// fingerprint.  With the reduced shock-front background, Omega must equal the
// event's ambient.rotation_rate_rad_per_s (checked in
// adapters/shock_front_background_adapter.cpp), so the ambient is static in
// the frame.
//
// TRANSFORMATIONS
// With theta(t) = Omega (t - t_ref) and Rz(theta) the right-handed rotation
// about +z,
//
//               | cos(theta)  -sin(theta)  0 |
//   Rz(theta) = | sin(theta)   cos(theta)  0 |
//               |     0            0       1 |
//
//   x_hci = Rz(theta) x_rot,                 x_rot = Rz(-theta) x_hci,
//   U_rot = Rz(-theta) (U_hci - Omega x x_hci),
//   B_rot = Rz(-theta) B_hci                 (and any other true vector),
//   grad-tensors rotate as T_rot = Rz(-theta) T_hci Rz(theta), except that the
//   velocity gradient also loses the rigid-rotation part:
//   grad U_rot = Rz(-theta) (grad U_hci - grad(Omega x x)) Rz(theta).
//
// Positions are heliocentric: a caller whose mesh has an application origin
// offset subtracts that origin before calling, and adds it back afterwards
// (see adapters/shock_front_background_adapter.cpp and
// adapters/reduced_front_shock_surface.cpp).  Epochs are the model's elapsed
// SI seconds.
//
// The velocity transformation is Galilean (no relativistic velocity
// addition).  The frame speed |Omega x x| = Omega r is 4.3e5 m/s at 1 au for
// Omega = 2.86533e-6 rad/s, and would reach c only at r = c/Omega ~ 700 au.
//
// WHY THE COROTATING FRAME
// (srcSEP3D/test/event_2013_04_11/README.md, section 2.1.)  SEPs move along
// field lines, and a field line rooted in the corona turns with the Sun.  In
// HCI the guiding centre is carried across the spiral by the perpendicular
// flow at exactly the pattern rate, so the population rotates.  AMPS does not
// move its mesh, so a domain fixed along the field line is a domain fixed in
// the corotating frame.  For a background that is steady in that frame (the
// v1.3 PFSS / open-tube / Parker ambient) the field, the field lines and the
// active tube are all static; only observers move.
//
// TRANSPORT IN THE FRAME
// The focused and Parker transport coefficients that srcSEP3D uses are
// unchanged by the change of frame: div(Omega x x) = 0 and
// bb:grad(Omega x x) = 0 because grad(Omega x x) is antisymmetric, so div U and
// bb:grad U are frame invariant, and the plasma-frame momentum is the same in
// both frames.  The optional inertial-force terms of the focused movers
// (Coriolis and centrifugal, application.sep.transport.frame_inertial_forces)
// are the only place where the frame enters the transport equations.  Those
// terms are assembled in adapters/transport_adapter.cpp from
// AngularVelocity(); see README.md of the event test, section 2.6.
// ============================================================================

#ifndef SEP3D_CORE_COMPUTATION_FRAME_H
#define SEP3D_CORE_COMPUTATION_FRAME_H

#include "sep3d_types.h"

#include <array>
#include <cmath>
#include <cstdio>
#include <string>

namespace SEP3D {
namespace Core {

// Hci: the provider's heliocentric inertial frame (identity transformation).
// Corotating: rigid rotation about HCI +z at Omega, equal to HCI at t_ref.
enum class ComputationFrameKind { Hci, Corotating };

// Lower-case frame-kind name for diagnostics and receipts.
inline const char* ComputationFrameName(ComputationFrameKind kind) {
  return kind == ComputationFrameKind::Corotating ? "corotating" : "hci";
}

// Value type describing the run's frame and converting points, vectors,
// velocities and gradient tensors between HCI and the frame.  Cheap to copy;
// every conversion is a pure function of its arguments and the epoch.  In the
// inertial frame every conversion returns its input unchanged.
class ComputationFrame {
 public:
  // The default-constructed frame is HCI.
  ComputationFrame() = default;
  static ComputationFrame Inertial() { return ComputationFrame(); }
  // Corotating frame with rate rotationRateRadPerS [rad/s] about HCI +z (a
  // negative rate rotates clockwise seen from +z) that equals HCI at
  // referenceEpochS [model elapsed s].  Not validated here: call Validate().
  static ComputationFrame Corotating(double rotationRateRadPerS,
                                     double referenceEpochS) {
    ComputationFrame frame;
    frame.kind_ = ComputationFrameKind::Corotating;
    frame.omega_ = rotationRateRadPerS;
    frame.referenceEpochS_ = referenceEpochS;
    return frame;
  }

  ComputationFrameKind Kind() const { return kind_; }
  bool IsInertial() const { return kind_ == ComputationFrameKind::Hci; }
  // Omega [rad/s]; exactly 0 for the inertial frame.
  double RotationRateRadPerS() const { return IsInertial() ? 0.0 : omega_; }
  // t_ref [model elapsed s] at which the frame equals HCI.
  double ReferenceEpochS() const { return referenceEpochS_; }
  // Angular velocity vector of the frame relative to HCI [rad/s].  It is
  // (0, 0, Omega) in both HCI and frame components, because the rotation is
  // about the z axis itself.
  Vec3 AngularVelocity() const { return {0.0, 0.0, RotationRateRadPerS()}; }

  // OK for the inertial frame.  A corotating frame needs a finite, nonzero
  // rate (either sign) and a finite reference epoch; otherwise InvalidInput.
  Status Validate() const {
    if (IsInertial()) return Status::OK();
    if (!std::isfinite(omega_) || omega_ == 0.0 ||
        !std::isfinite(referenceEpochS_))
      return Status(StatusCode::InvalidInput,
                    "corotating frame requires a finite nonzero rotation rate "
                    "and a finite reference epoch");
    return Status::OK();
  }

  // Rotation angle theta(t) = Omega (t - t_ref) [rad] of the frame at epoch t
  // [model elapsed s]; zero at the reference epoch and for the inertial frame.
  // The angle is not reduced modulo 2 pi (std::cos/std::sin reduce it).
  double AngleRad(double epochS) const {
    return IsInertial() ? 0.0 : omega_ * (epochS - referenceEpochS_);
  }

  // Frame point -> HCI point, x_hci = Rz(theta) x_rot.  Heliocentric [m].
  Vec3 ToInertialPosition(const Vec3& x, double epochS) const {
    return Rotate(x, AngleRad(epochS));
  }
  // HCI point -> frame point, x_rot = Rz(-theta) x_hci.  Heliocentric [m].
  Vec3 FromInertialPosition(const Vec3& x, double epochS) const {
    return Rotate(x, -AngleRad(epochS));
  }
  // A true vector (B, a normal, a displacement) from HCI into the frame:
  // rotation only, v_rot = Rz(-theta) v_hci.  Not for velocities of matter,
  // which also lose the frame velocity (FromInertialVelocity).
  Vec3 FromInertialVector(const Vec3& v, double epochS) const {
    return Rotate(v, -AngleRad(epochS));
  }
  // Inverse of FromInertialVector: v_hci = Rz(theta) v_rot.
  Vec3 ToInertialVector(const Vec3& v, double epochS) const {
    return Rotate(v, AngleRad(epochS));
  }
  // Velocity of matter at HCI point xHci [m]: subtract the frame velocity
  // there and rotate, U_rot = Rz(-theta) (U_hci - Omega x x_hci) [m/s].
  // For the corotating solar wind of the reduced ambient this is the
  // field-aligned flow U' of the event test README (section 2.1).
  Vec3 FromInertialVelocity(const Vec3& velocityHci, const Vec3& xHci,
                            double epochS) const {
    if (IsInertial()) return velocityHci;
    return Rotate(velocityHci - AngularVelocity().Cross(xHci),
                  -AngleRad(epochS));
  }
  // Inverse of FromInertialVelocity for matter at frame point xFrame [m]:
  // U_hci = Rz(theta) U_rot + Omega x x_hci, with x_hci = Rz(theta) xFrame.
  Vec3 ToInertialVelocity(const Vec3& velocityFrame, const Vec3& xFrame,
                          double epochS) const {
    if (IsInertial()) return velocityFrame;
    const Vec3 xHci = ToInertialPosition(xFrame, epochS);
    return Rotate(velocityFrame, AngleRad(epochS)) +
           AngularVelocity().Cross(xHci);
  }

  // Row-major 3x3 tensor T_ij = d v_i / d x_j of a true vector field,
  // element (i, j) at index 3 i + j.
  using Tensor = std::array<double, 9>;
  // Gradient of a true vector field (e.g. grad B [T/m]) from HCI into the
  // frame.  With R = Rz(-theta), v_rot = R v_hci and x_hci = R^T x_rot, the
  // chain rule gives d v_rot / d x_rot = R (d v_hci / d x_hci) R^T, which is
  // the product computed below (first R T, then (R T) R^T).
  Tensor FromInertialTensor(const Tensor& t, double epochS) const {
    if (IsInertial()) return t;
    const double a = -AngleRad(epochS);
    const double c = std::cos(a), s = std::sin(a);
    const double r[9] = {c, -s, 0, s, c, 0, 0, 0, 1};   // Rz(a)
    Tensor tmp{}, out{};
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j) {
        double sum = 0.0;
        for (int k = 0; k < 3; ++k) sum += r[3 * i + k] * t[3 * k + j];
        tmp[3 * i + j] = sum;
      }
    for (int i = 0; i < 3; ++i)                 // tmp * Rz(a)^T
      for (int j = 0; j < 3; ++j) {
        double sum = 0.0;
        for (int k = 0; k < 3; ++k) sum += tmp[3 * i + k] * r[3 * j + k];
        out[3 * i + j] = sum;
      }
    return out;
  }
  // Velocity gradient [1/s]: also remove the rigid-rotation part
  // grad(Omega x x) before rotating the tensor into the frame,
  //
  //                     |   0    -Omega   0 |
  //   grad(Omega x x) = | Omega     0     0 |   (T_ij = d(Omega x x)_i / d x_j).
  //                     |   0       0     0 |
  //
  // This tensor is antisymmetric (so div U and bb:grad U do not change) and
  // commutes with Rz, so subtracting it before or after the rotation gives
  // the same result.
  Tensor FromInertialVelocityGradient(const Tensor& gradU, double epochS) const {
    if (IsInertial()) return gradU;
    Tensor g = gradU;
    // Omega x x = (-Omega y, Omega x, 0):  d/dy of x-component = -Omega,
    // d/dx of y-component = +Omega.
    g[1] += omega_;   // (i=x, j=y)
    g[3] -= omega_;   // (i=y, j=x)
    return FromInertialTensor(g, epochS);
  }

  // Inertial accelerations per unit mass [m/s^2] felt in the frame by a body
  // at frame position x [m] moving with frame velocity v [m/s] (Coriolis
  // -2 Omega x v and centrifugal -Omega x (Omega x x)); zero in the inertial
  // frame.  The focused movers' optional inertial terms do not call these
  // helpers: adapters/transport_adapter.cpp assembles the plasma-frame
  // acceleration 2 Omega x U' (+ Omega x (Omega x x)) from AngularVelocity().
  Vec3 CoriolisAcceleration(const Vec3& velocityFrame) const {
    return AngularVelocity().Cross(velocityFrame) * -2.0;
  }
  Vec3 CentrifugalAcceleration(const Vec3& xFrame) const {
    const Vec3 w = AngularVelocity();
    return w.Cross(w.Cross(xFrame)) * -1.0;
  }

  // Canonical text for fingerprints and receipts: "frame=hci", or the rate
  // and reference epoch at %.17g, which round-trips every double exactly, so
  // two runs share a fingerprint only when their frames are bit-identical.
  std::string Canonical() const;

 private:
  // v rotated by angleRad about +z, i.e. Rz(angleRad) v.  A zero angle returns
  // v unchanged, so the conversions are the bit-exact identity at t_ref.
  static Vec3 Rotate(const Vec3& v, double angleRad) {
    if (angleRad == 0.0) return v;
    const double c = std::cos(angleRad), s = std::sin(angleRad);
    return {c * v.x - s * v.y, s * v.x + c * v.y, v.z};
  }
  ComputationFrameKind kind_ = ComputationFrameKind::Hci;
  double omega_ = 0.0;            // Omega [rad/s]; unused when kind_ is Hci
  double referenceEpochS_ = 0.0;  // t_ref [model elapsed s]
};

inline std::string ComputationFrame::Canonical() const {
  if (IsInertial()) return "frame=hci";
  char buffer[160];
  std::snprintf(buffer, sizeof(buffer),
                "frame=corotating;omega_rad_per_s=%.17g;reference_epoch_s=%.17g",
                omega_, referenceEpochS_);
  return buffer;
}

}  // namespace Core
}  // namespace SEP3D

#endif  // SEP3D_CORE_COMPUTATION_FRAME_H
