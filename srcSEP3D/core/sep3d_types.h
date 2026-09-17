// ============================================================================
// srcSEP3D/core/sep3d_types.h
//
// Fundamental types, physical constants, and status codes used throughout the
// srcSEP3D application.
//
// LAYER: L0 (numerical core).
//   This header must NEVER include pic.h, mpi.h, or any AMPS/MPI symbol.
//   The make check-layering target enforces this at build time.
//   Rationale: L0 files must compile into the Stage-1 standalone test binary
//   that links neither AMPS nor MPI.  Keeping this header free of those
//   dependencies is what makes unit-testing the numerical kernels possible
//   without building the entire AMPS framework.
//
// PHYSICAL CONSTANTS:
//   All values are copied from the SWCME library headers by value (not by
//   reference to SWCME symbols) so that this header has no dependency on
//   SWCME.  The values are verified in test GEO03 to be bitwise equal to
//   SWCME's copies.  Any discrepancy there is a bug here, not in SWCME.
//
// ============================================================================

#ifndef SEP3D_TYPES_H
#define SEP3D_TYPES_H

#include <cmath>
#include <cstdint>
#include <string>

namespace SEP3D {
namespace Core {

// ----------------------------------------------------------------------------
// Physical constants — SI units throughout.
//
// Values are consistent with IAU 2012 / NIST CODATA 2018.
// GEO03 asserts bitwise equality with SWCME's copies of these constants.
// ----------------------------------------------------------------------------
namespace Const {

  // Speed of light [m/s]
  constexpr double c      = 2.99792458e+08;
  // Elementary charge [C]
  constexpr double e      = 1.602176634e-19;
  // Proton mass [kg]
  constexpr double m_p    = 1.67262192369e-27;
  // Electron mass [kg]
  constexpr double m_e    = 9.1093837015e-31;
  // Boltzmann constant [J/K]
  constexpr double k_B    = 1.380649e-23;
  // Gravitational constant [m^3 kg^-1 s^-2]
  constexpr double G      = 6.674e-11;
  // Solar mass [kg]
  constexpr double M_sun  = 1.989e+30;
  // Solar radius [m]
  constexpr double R_sun  = 6.957e+08;
  // 1 AU [m]
  constexpr double AU     = 1.495978707e+11;
  // Solar rotation rate [rad/s]  (sidereal, equatorial)
  constexpr double Omega_sun = 2.865e-06;
  // Standard solar wind speed [m/s]  (used as a default; overridable via CLI)
  constexpr double V_sw_default = 4.0e+05;
  // Pi.  The k-prefix is required at the AMPS boundary: general/constants.h
  // defines a legacy function-like global token named Pi as a preprocessor
  // macro.  A namespace cannot protect an identifier from macro expansion, so
  // `constexpr double Pi` breaks every AMPS translation unit that includes
  // pic.h before this header (including cell-centered interpolation).
  constexpr double kPi    = 3.14159265358979323846;

} // namespace Const


// ----------------------------------------------------------------------------
// Status and StatusCode
//
// Every function that can fail returns a Status rather than a bare bool or
// int.  The StatusCode enum provides fine-grained diagnostics that can be
// logged, counted, and branched on without string parsing.
//
// StatusCode taxonomy (mirrors SEP::Transport::StatusCode in srcSEP so that
// the two applications share the same diagnostic vocabulary):
//
//   OK             - normal completion
//   Ballistic      - coefficient provider returned the λ=+∞ zero-rate state;
//                    the particle propagates freely this step
//   StepUnderflow  - the computed substep fell below the declared floor;
//                    the mover records the event and takes the floor step
//   ReservedFeature- an option that is specified but not yet implemented was
//                    encountered; treated as a hard error at parse time
//   InvalidInput   - a numerical argument is NaN, Inf, or out of range
//   NotFound       - a lookup (e.g. AMR node search) returned no result
//   BackgroundInvalid - the background validity flag is clear at this point
//   DomainExit     - the particle left the computational domain
//   InnerBoundary  - the particle crossed the inner sphere r < r_in
//   Error          - generic unrecoverable error; message carries details
// ----------------------------------------------------------------------------
enum class StatusCode : int {
  OK              =  0,
  Ballistic       =  1,
  StepUnderflow   =  2,
  ReservedFeature =  3,
  InvalidInput    =  4,
  NotFound        =  5,
  BackgroundInvalid = 6,
  DomainExit      =  7,
  InnerBoundary   =  8,
  Error           = -1
};

struct Status {
  StatusCode code    = StatusCode::OK;
  std::string message;                  // human-readable; empty when code==OK

  // Convenience constructor for the common success case
  Status() = default;
  explicit Status(StatusCode c, std::string msg = {})
    : code(c), message(std::move(msg)) {}

  // Returns true iff the status represents successful completion or the
  // special-case ballistic state (which is a valid physics outcome, not
  // an error).
  bool ok()       const { return code == StatusCode::OK; }
  bool ballistic() const { return code == StatusCode::Ballistic; }
  bool usable()   const { return ok() || ballistic(); }

  static Status OK()     { return Status{StatusCode::OK}; }
  static Status Ballistic(std::string msg = {}) {
    return Status{StatusCode::Ballistic, std::move(msg)};
  }
  static Status Error(std::string msg) {
    return Status{StatusCode::Error, std::move(msg)};
  }
  static Status Reserved(std::string feature) {
    return Status{StatusCode::ReservedFeature,
                  "reserved, not implemented in this release: " + feature};
  }
};


// ----------------------------------------------------------------------------
// Vec3 — a plain three-component double vector.
//
// Deliberately not a general linear-algebra type: only the operations
// actually used in the transport kernels are provided.  This keeps the
// file free of template machinery and makes the dependency chain trivial.
// ----------------------------------------------------------------------------
struct Vec3 {
  double x = 0.0, y = 0.0, z = 0.0;

  Vec3() = default;
  Vec3(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}

  // Construct from a C array (the interface used by AMPS internals)
  explicit Vec3(const double* p) : x(p[0]), y(p[1]), z(p[2]) {}

  // Write back into a C array
  void CopyTo(double* p) const { p[0]=x; p[1]=y; p[2]=z; }

  double  Dot(const Vec3& b)   const { return x*b.x + y*b.y + z*b.z; }
  double  NormSq()             const { return Dot(*this); }
  double  Norm()               const { return std::sqrt(NormSq()); }

  Vec3 Cross(const Vec3& b) const {
    return {y*b.z - z*b.y,
            z*b.x - x*b.z,
            x*b.y - y*b.x};
  }

  // Returns the unit vector, or the zero vector if |this| < eps.
  // eps defaults to 1e-300 to catch true zeros without false positives
  // from small physical fields.
  Vec3 Normalized(double eps = 1e-300) const {
    double n = Norm();
    if (n < eps) return {};
    return {x/n, y/n, z/n};
  }

  Vec3 operator+(const Vec3& b) const { return {x+b.x, y+b.y, z+b.z}; }
  Vec3 operator-(const Vec3& b) const { return {x-b.x, y-b.y, z-b.z}; }
  Vec3 operator*(double s)      const { return {x*s,   y*s,   z*s};   }
  Vec3 operator/(double s)      const { return {x/s,   y/s,   z/s};   }
  Vec3& operator+=(const Vec3& b){ x+=b.x; y+=b.y; z+=b.z; return *this; }
  Vec3& operator*=(double s)     { x*=s;   y*=s;   z*=s;   return *this; }

  bool operator==(const Vec3& b) const {
    return x == b.x && y == b.y && z == b.z;
  }
};

inline Vec3 operator*(double s, const Vec3& v) { return v * s; }


// ----------------------------------------------------------------------------
// Tensor3 — a 3x3 matrix stored row-major: T[i][j] = dF_i/dx_j.
//
// Used for the magnetic-field gradient ∂B_i/∂x_j and the velocity gradient
// ∂U_i/∂x_j.  Row-major storage is chosen to match the C multi-dimensional
// array layout, which makes the AMPS central-difference stencil loops
// natural.
// ----------------------------------------------------------------------------
struct Tensor3 {
  double m[3][3] = {};          // m[row][col], zero-initialised

  // Element access
  double  operator()(int i, int j) const { return m[i][j]; }
  double& operator()(int i, int j)       { return m[i][j]; }

  // Trace: Σ m[i][i]
  double Trace() const { return m[0][0] + m[1][1] + m[2][2]; }

  // Double contraction with a pair of unit vectors: b·T·b = Σ b_i T_ij b_j
  // Used to compute b̂b̂:∇U = the field-aligned strain.
  double DoubleContract(const Vec3& a, const Vec3& b) const {
    double result = 0.0;
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
        result += a.x*(i==0) * m[i][j] * (b.x*(j==0));
    // Compact form:
    const double* av = &a.x;
    const double* bv = &b.x;
    result = 0.0;
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
        result += av[i] * m[i][j] * bv[j];
    return result;
  }

  // Matrix-vector product: T·v
  Vec3 Apply(const Vec3& v) const {
    return { m[0][0]*v.x + m[0][1]*v.y + m[0][2]*v.z,
             m[1][0]*v.x + m[1][1]*v.y + m[1][2]*v.z,
             m[2][0]*v.x + m[2][1]*v.y + m[2][2]*v.z };
  }
};


// ----------------------------------------------------------------------------
// SpeciesProperties — per-species physical parameters.
//
// Stored once at initialisation.  The mover reads it on every substep so
// the layout is kept compact.  Charge is signed (negative for electrons).
// ----------------------------------------------------------------------------
struct SpeciesProperties {
  int    index       = 0;       // species index in the AMPS species table
  double massKg      = 0.0;     // rest mass [kg]
  double chargeCoul  = 0.0;     // signed electric charge [C]
  int    nucleons    = 0;       // nucleon count (for rigidity computation)
  std::string name;             // human-readable, e.g. "h+", "e-", "he2+"
};


// ----------------------------------------------------------------------------
// ParticleMotionOutcome — semantic result of one production mover call.
//
// IMPORTANT: this enum is deliberately *not* assigned AMPS integer values.
// The prototype copied values into the numerical core and they did not match
// pic.h: MotionFinished was recorded as zero although AMPS defines it as
// three, and a srcSEP3D extension collided with that real AMPS value.  Such a
// mismatch can delete or retain the wrong particle while still compiling.
//
// The only conversion to the AMPS ABI is now in
// amps/amps_mover_status.h.  That adapter includes pic.h, uses the actual
// _PARTICLE_* macros, and contains compile-time assertions.  Keeping the
// semantic enum here lets L0 remain independent of AMPS and prevents future
// physics/status additions from becoming accidental mover return codes.
// ----------------------------------------------------------------------------
enum class ParticleMotionOutcome {
  Advanced,
  LeftDomain
};

} // namespace Core
} // namespace SEP3D

#endif // SEP3D_TYPES_H
