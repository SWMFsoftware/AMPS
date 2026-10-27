#ifndef _SRC_GEOSPACE_GPU_ENGINE_CONFIG_H_
#define _SRC_GEOSPACE_GPU_ENGINE_CONFIG_H_

//======================================================================================
// EngineConfig.h -- portability layer of the batched backward-trajectory engine
//======================================================================================
//
// The engine source is a single code base that is compiled in three ways:
//
//   (1) g++ -fopenmp            CPU build.  Kernels are OpenMP host loops.  This is the
//                               build that is linked into the default AMPS executable
//                               and that reproduces the legacy tracer bit for bit in double/SI.
//   (2) nvc++ -acc=gpu          GPU build.  Kernels are OpenACC compute regions; every
//                               function reached from a kernel is a header-only inline
//                               template, so nvc++ generates "acc routine seq" device
//                               versions implicitly (verified with nvc++ 24.5).
//   (3) nvc++ -acc=multicore    OpenACC on CPU cores (optional).
//
// No device code uses virtual functions, exceptions, std::vector, std::string or
// dynamic allocation.  Host-only magnetic-field models (Geopack/Tsyganenko) are reached
// only through the host executor, never from a device kernel.
//
// PRECISION AND UNITS
// -------------------
// The trajectory state, the field sample and the mover arithmetic use the compile-time
// type EG_REAL together with a unit system:
//
//   GEOSPACE_PRECISION_FLOAT undefined (default)
//       Real=double, SI units.  Every expression is a literal port of the legacy code
//       code; all unit conversions are compile-time identities, so the engine result is
//       bitwise identical to the legacy tracer when compiled with the same compiler and flags.
//
//   GEOSPACE_PRECISION_FLOAT defined
//       Real=float in normalized geospace units: length Re, speed c, mass amu,
//       charge e, time Re/c, momentum amu*c, magnetic field amu*c/(e*Re) (~488 nT).
//       In SI, |p|^2 and (m c)^2 of protons/electrons underflow single precision
//       (|p|^2 ~ 1e-40..1e-43); in these units all state magnitudes are O(1e-4..1e4).
//       Time, path length, trap-detector diagnostics and boundary products are kept in
//       double precision.  See srcGeospace/doc/PRECISION.md.
//
// The AMPS PIC::Units normalization used by the field solver (mSI=amu, lSI=1e-7 e^2/amu,
// uSI=1e6 m/s) was evaluated for this purpose and rejected: it produces |x|^2~6e52 at
// 60 Re and |B|^2~2e-50 at 1 nT, i.e. it overflows/underflows float itself.
//======================================================================================

#include <cmath>
#include <cstdint>
#include <limits>

#if defined(_OPENACC)
  #define EG_HAVE_OPENACC 1
#else
  #define EG_HAVE_OPENACC 0
#endif

#define EG_PRAGMA_(x) _Pragma(#x)
#define EG_PRAGMA(x) EG_PRAGMA_(x)

// Inline marker for functions callable from kernels.  nvc++ generates device versions
// of inline/template functions implicitly; the macro documents the device contract.
#define EG_HD inline

namespace Geospace {
namespace GPU {

#if defined(GEOSPACE_PRECISION_FLOAT)
typedef float Real;
#define EG_REAL_NAME "float"
#else
typedef double Real;
#define EG_REAL_NAME "double"
#endif

//--------------------------------------------------------------------------------------
// Unit systems.  Each provides SI->internal conversion factors; the SI system is the
// identity and uses `if constexpr` at every conversion so no rounding is introduced.
//--------------------------------------------------------------------------------------
struct UnitsSI {
  static constexpr bool kSI=true;
  static constexpr double kLength_m=1.0;          // internal length unit [m]
  static constexpr double kSpeed_m_s=1.0;         // internal speed unit [m/s]
  static constexpr double kMass_kg=1.0;
  static constexpr double kCharge_C=1.0;
  static constexpr double kTime_s=1.0;
  static constexpr double kMomentum_SI=1.0;
  static constexpr double kField_T=1.0;
  static constexpr double kSpeedOfLight=299792458.0; // internal value of c
};

struct UnitsGeo {
  static constexpr bool kSI=false;
  static constexpr double kLength_m=6371.0E3;                 // _EARTH__RADIUS_
  static constexpr double kSpeed_m_s=299792458.0;             // c
  static constexpr double kMass_kg=1.66053906660E-27;         // _AMU_
  static constexpr double kCharge_C=1.602176565E-19;          // ElectronCharge (AMPS)
  static constexpr double kTime_s=kLength_m/kSpeed_m_s;
  static constexpr double kMomentum_SI=kMass_kg*kSpeed_m_s;
  static constexpr double kField_T=kMass_kg*kSpeed_m_s/(kCharge_C*kLength_m);
  static constexpr double kSpeedOfLight=1.0;
};

#if defined(GEOSPACE_PRECISION_FLOAT)
typedef UnitsGeo Units;
#else
typedef UnitsSI Units;
#endif

// SI -> internal conversions.  For UnitsSI these return the argument unchanged.
template <class U> EG_HD double ToLength(double v_m) {
  if constexpr (U::kSI) return v_m; else return v_m/U::kLength_m;
}
template <class U> EG_HD double ToSpeed(double v) {
  if constexpr (U::kSI) return v; else return v/U::kSpeed_m_s;
}
template <class U> EG_HD double ToTime(double v_s) {
  if constexpr (U::kSI) return v_s; else return v_s/U::kTime_s;
}
template <class U> EG_HD double ToMass(double v_kg) {
  if constexpr (U::kSI) return v_kg; else return v_kg/U::kMass_kg;
}
template <class U> EG_HD double ToCharge(double v_C) {
  if constexpr (U::kSI) return v_C; else return v_C/U::kCharge_C;
}
template <class U> EG_HD double ToMomentum(double v) {
  if constexpr (U::kSI) return v; else return v/U::kMomentum_SI;
}
template <class U> EG_HD double ToField(double v_T) {
  if constexpr (U::kSI) return v_T; else return v_T/U::kField_T;
}
// internal -> SI
template <class U> EG_HD double FromLength(double v) {
  if constexpr (U::kSI) return v; else return v*U::kLength_m;
}
template <class U> EG_HD double FromTime(double v) {
  if constexpr (U::kSI) return v; else return v*U::kTime_s;
}
template <class U> EG_HD double FromMomentum(double v) {
  if constexpr (U::kSI) return v; else return v*U::kMomentum_SI;
}
template <class U> EG_HD double FromField(double v) {
  if constexpr (U::kSI) return v; else return v*U::kField_T;
}

// Sentinel used by the legacy code for "no limit" distances (1e300 / 1e299).  float cannot
// represent it, so the normalized build uses a large finite float value instead.
template <class R> struct Sentinel;
template <> struct Sentinel<double> {
  static constexpr double kBig=1.0e300;
  static constexpr double kBigLimit=1.0e299;
};
template <> struct Sentinel<float> {
  static constexpr float kBig=1.0e30f;
  static constexpr float kBigLimit=1.0e29f;
};

//--------------------------------------------------------------------------------------
// Exact equivalents of std::min/std::max/std::isfinite usable in device code.  std::min
// returns (b<a)?b:a and std::max returns (a<b)?b:a; keeping these definitions preserves
// NaN propagation semantics of the legacy expressions.
//--------------------------------------------------------------------------------------
template <class T> EG_HD T emin(T a,T b) { return (b<a) ? b : a; }
template <class T> EG_HD T emax(T a,T b) { return (a<b) ? b : a; }
EG_HD bool efinite(double v) { return std::isfinite(v); }
EG_HD bool efinite(float v) { return std::isfinite(v); }

// std::fmod for device code.  nvc++ 24.5 has no device implementation of fmod/fmodf
// (-gpu=nordc: "ptxas fatal: Unresolved extern function 'fmod'"; default rdc: "nvlink
// error: Undefined reference to 'fmod'" when the device code is linked), while
// remainder() is available.  Both fmod and remainder are exact operations:
// remainder(x,y)=x-n*y with n=nearest integer to x/y, |r|<=|y|/2.  When r has the
// opposite sign of x, fmod(x,y)=r+sign(x)*|y| is exactly representable, so the IEEE
// addition returns it exactly; the result is bitwise identical to std::fmod for all
// inputs, including zeros, infinities and NaN.  Host-compiled (g++) code keeps
// std::fmod.
template <class T> EG_HD T efmod(T x,T y) {
#if EG_HAVE_OPENACC
  const T r=std::remainder(x,y);
  if (x>=T(0)) return (r<T(0)) ? r+std::fabs(y) : r;
  return (r>T(0)) ? r-std::fabs(y) : r;
#else
  return std::fmod(x,y);
#endif
}

} // namespace GPU
} // namespace Geospace

#endif
