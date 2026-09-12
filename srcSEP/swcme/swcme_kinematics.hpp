#pragma once

// ============================================================================
// swcme_kinematics.hpp
// ----------------------------------------------------------------------------
// Shared CME/shock-apex kinematics used by BOTH the 1-D and 3-D SWCME models.
//
// Why this file exists
// --------------------
// The original SWCME code carried separate DBM formulas in swcme1d.hpp and
// swcme3d.cpp.  Those copies had diverged in two important limiting cases:
//   * the 1-D slow-CME branch clipped V0-Vsw to zero and instantaneously placed
//     a slow CME at the ambient wind speed;
//   * the 3-D implementation used the fast-CME algebra for negative V0-Vsw and
//     divided by Gamma when Gamma=0, after which a finite-value fallback could
//     hide the invalid radius.
//
// This header provides one mathematically sign-aware, numerically stable
// implementation so every dimensional wrapper obtains exactly the same apex
// radius and speed for the same physical configuration.
//
// Supported modes
// ---------------
//   BALLISTIC   R(t)=R0+V0*t, V(t)=V0.
//   DBM         Constant-background drag-based model
//                 d(Delta V)/dt = -Gamma DeltaV |DeltaV|,
//               valid for both fast and slow CMEs.
//   DATA_DRIVEN Monotonic PCHIP interpolation of observed/model height-time
//               knots.  The derivative of the same interpolant is returned as
//               the apex speed, so radius and speed remain kinematically
//               consistent.
//
// Units
// -----
// All quantities in this common layer are SI: meters, seconds, m/s and 1/m.
// The 1-D and 3-D public parameter structures remain free to expose convenient
// heliophysics units and convert once when building KinematicsConfig.
// ============================================================================

#include "swcme_solarwind.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace swcme {
namespace kinematics {

enum class Mode { Ballistic = 0, DBM = 1, DataDriven = 2 };

// The default data-driven policy deliberately rejects queries outside the
// supplied measurement interval.  Extrapolation changes the physical model
// and therefore must be selected explicitly by the caller.
enum class ExtrapolationPolicy { OutsideTime = 0, Ballistic = 1 };

enum class Status {
  Ok = 0,
  InvalidInput,
  OutsideTime,
  OutsideDomain
};

inline const char* status_name(Status s) {
  switch (s) {
    case Status::Ok: return "OK";
    case Status::InvalidInput: return "INVALID_INPUT";
    case Status::OutsideTime: return "OUTSIDE_TIME";
    case Status::OutsideDomain: return "OUTSIDE_DOMAIN";
  }
  return "UNKNOWN";
}

struct State {
  Status status = Status::InvalidInput;
  double radius_m = std::numeric_limits<double>::quiet_NaN();
  double speed_m_s = std::numeric_limits<double>::quiet_NaN();
};

struct Config {
  Mode mode = Mode::DBM;

  // BALLISTIC/DBM parameters.  The reference radius must lie inside the
  // analytical solar-wind domain, V0 may be zero for a stationary front, and
  // the ambient wind is strictly outward.  KIN09 applies the same physical
  // postconditions to every extrapolated state.
  double r0_m = 0.0;
  double V0_m_s = 0.0;
  double Vsw_m_s = 0.0;
  double Gamma_m_inv = 0.0;

  // DATA_DRIVEN height-time knots.  Radius must be nondecreasing and time must
  // be strictly increasing.  A nondecreasing table is sufficient for PCHIP;
  // flat intervals are allowed and produce zero local propagation speed.  All
  // radii must lie at or above solarwind::MIN_RADIUS_M.
  std::vector<double> data_time_s;
  std::vector<double> data_radius_m;
  ExtrapolationPolicy extrapolation = ExtrapolationPolicy::OutsideTime;
};

// ----------------------------------------------------------------------------
// Input validation is retained locally as a defense-in-depth boundary for
// direct users of this common component.  The public 1-D/3-D validators provide
// field-specific diagnostics first; this lower layer independently prevents a
// malformed SI config from being clipped or interpreted as an extrapolation.
// ----------------------------------------------------------------------------
inline bool valid_scalar(double x) { return std::isfinite(x); }

inline Status validate(const Config& c) {
  if (!valid_scalar(c.r0_m) || !valid_scalar(c.V0_m_s) ||
      !valid_scalar(c.Vsw_m_s) || !valid_scalar(c.Gamma_m_inv)) {
    return Status::InvalidInput;
  }
  if (c.r0_m < swcme::solarwind::MIN_RADIUS_M || c.V0_m_s < 0.0 ||
      !(c.Vsw_m_s > 0.0) || c.Gamma_m_inv < 0.0) {
    return Status::InvalidInput;
  }

  if (c.mode != Mode::DataDriven) {
    return Status::Ok;
  }

  if (c.data_time_s.size() < 2 ||
      c.data_time_s.size() != c.data_radius_m.size()) {
    return Status::InvalidInput;
  }
  for (std::size_t i=0; i<c.data_time_s.size(); ++i) {
    if (!valid_scalar(c.data_time_s[i]) || !valid_scalar(c.data_radius_m[i]) ||
        c.data_radius_m[i] < swcme::solarwind::MIN_RADIUS_M) {
      return Status::InvalidInput;
    }
    if (i>0) {
      if (!(c.data_time_s[i] > c.data_time_s[i-1])) {
        return Status::InvalidInput;
      }
      if (c.data_radius_m[i] < c.data_radius_m[i-1]) {
        return Status::InvalidInput;
      }
    }
  }
  return Status::Ok;
}

// Convert a candidate formula result into the public kinematic-state contract.
// Input validation and physical-domain validation are intentionally distinct:
// malformed parameters return INVALID_INPUT, whereas a finite-time extension
// of a valid model that leaves r>=1.05 R_sun, overflows, or reverses sunward
// returns OUTSIDE_DOMAIN.  Failed states retain NaN payloads so callers cannot
// accidentally consume a stale/partial radius while ignoring the status.
inline State domain_checked_state(double radius_m,double speed_m_s) {
  State out;
  if (!std::isfinite(radius_m) || !std::isfinite(speed_m_s) ||
      radius_m<swcme::solarwind::MIN_RADIUS_M || speed_m_s<0.0) {
    out.status=Status::OutsideDomain;
    return out;
  }
  out.status=Status::Ok;
  out.radius_m=radius_m;
  out.speed_m_s=speed_m_s;
  return out;
}

// ----------------------------------------------------------------------------
// Ballistic reference/operating mode.
// ----------------------------------------------------------------------------
inline State ballistic_state(const Config& c, double t_s) {
  State out;
  if (validate(c) != Status::Ok || !std::isfinite(t_s)) {
    out.status = Status::InvalidInput;
    return out;
  }
  // Negative time is a legitimate backward extrapolation until the analytical
  // trajectory crosses the configured radial domain.  fma performs the one
  // multiplication/addition with a single rounding and makes overflow visible
  // to the common postcondition instead of returning an unchecked OK state.
  return domain_checked_state(std::fma(c.V0_m_s,t_s,c.r0_m),c.V0_m_s);
}

// ----------------------------------------------------------------------------
// Sign-aware drag-based model (DBM).
//
// Let DeltaV0=V0-Vsw and a=|DeltaV0|.  The exact solution of
//   d(DeltaV)/dt = -Gamma DeltaV |DeltaV|
// is
//   DeltaV(t)=DeltaV0/(1+Gamma a t),
//   R(t)=R0+Vsw t + sign(DeltaV0)/Gamma * log(1+Gamma a t).
// The absolute value is essential: it makes drag decelerate fast CMEs and
// accelerate slow CMEs toward Vsw.
//
// Numerics near Gamma=0
// ---------------------
// Gamma=0 has an exact ballistic solution and is handled explicitly.  For very
// small x=Gamma*a*t, evaluating log1p(x)/Gamma is already stable, but the ratio
// still contains an avoidable 0/0 sensitivity as Gamma approaches machine
// underflow.  We therefore use a short Taylor series for |x|<1e-8:
//   log(1+x)/Gamma = a t [1-x/2+x^2/3-x^3/4+x^4/5+...].
// This makes the DBM-to-ballistic transition continuous without a tolerance-
// sized jump in radius.
// ----------------------------------------------------------------------------
inline State dbm_state(const Config& c, double t_s) {
  State out;
  if (validate(c) != Status::Ok || !std::isfinite(t_s)) {
    out.status = Status::InvalidInput;
    return out;
  }

  const double dv0 = c.V0_m_s - c.Vsw_m_s;
  const double a = std::abs(dv0);

  // No initial speed difference means drag is irrelevant even when Gamma>0.
  if (c.Gamma_m_inv == 0.0 || a == 0.0) {
    return domain_checked_state(
        std::fma(c.V0_m_s,t_s,c.r0_m),c.V0_m_s);
  }

  const double x = c.Gamma_m_inv*a*t_s;
  // Backward DBM continuation exists only while 1+Gamma*|dV0|*t is positive.
  // At/beyond that pole the closed form is not a physical branch.  Overflow
  // of x likewise means the requested finite time cannot be represented by a
  // finite model state, so both cases are explicit domain failures.
  if (!std::isfinite(x) || !(1.0+x>0.0)) {
    out.status = Status::OutsideDomain;
    return out;
  }

  const double denom = 1.0 + x;
  const double dv = dv0/denom;

  const double sign = (dv0 > 0.0) ? 1.0 : -1.0;
  double drag_distance = 0.0;
  if (std::abs(x) < 1.0e-8) {
    // log1p(x)/Gamma = a*t * log1p(x)/x.  Expanding the final ratio avoids
    // loss of continuity when Gamma is tiny while preserving the exact
    // ballistic limit as x->0.
    const double x2=x*x, x3=x2*x, x4=x3*x;
    const double log1p_over_x = 1.0 - 0.5*x + x2/3.0 - 0.25*x3 + 0.2*x4;
    drag_distance = sign*a*t_s*log1p_over_x;
  } else {
    drag_distance = sign*std::log1p(x)/c.Gamma_m_inv;
  }

  const double radius_m=std::fma(c.Vsw_m_s,t_s,c.r0_m)+drag_distance;
  const double speed_m_s=c.Vsw_m_s+dv;
  return domain_checked_state(radius_m,speed_m_s);
}

// ----------------------------------------------------------------------------
// Monotone PCHIP helpers.
//
// These implement the Fritsch-Carlson/Fritsch-Butland monotone cubic Hermite
// construction.  The interpolation passes exactly through every supplied
// height-time knot and, for a nondecreasing radius table, does not introduce
// cubic overshoot or a physically spurious sunward segment.
// ----------------------------------------------------------------------------
inline double pchip_endpoint_slope(double h0, double h1,
                                   double d0, double d1) {
  double m = ((2.0*h0 + h1)*d0 - h0*d1)/(h0+h1);
  if (m*d0 <= 0.0) return 0.0;
  if ((d0*d1 < 0.0) && (std::abs(m) > std::abs(3.0*d0))) {
    return 3.0*d0;
  }
  return m;
}

inline std::vector<double> pchip_slopes(const std::vector<double>& x,
                                        const std::vector<double>& y) {
  const std::size_t n=x.size();
  std::vector<double> m(n,0.0);
  if (n<2 || y.size()!=n) return m;

  if (n==2) {
    const double d=(y[1]-y[0])/(x[1]-x[0]);
    m[0]=d; m[1]=d;
    return m;
  }

  std::vector<double> h(n-1), d(n-1);
  for (std::size_t i=0;i+1<n;++i) {
    h[i]=x[i+1]-x[i];
    d[i]=(y[i+1]-y[i])/h[i];
  }

  m[0]=pchip_endpoint_slope(h[0],h[1],d[0],d[1]);
  for (std::size_t k=1;k+1<n;++k) {
    if (d[k-1]==0.0 || d[k]==0.0 || d[k-1]*d[k]<=0.0) {
      m[k]=0.0;
    } else {
      const double w1=2.0*h[k]+h[k-1];
      const double w2=h[k]+2.0*h[k-1];
      m[k]=(w1+w2)/(w1/d[k-1]+w2/d[k]);
    }
  }
  m[n-1]=pchip_endpoint_slope(h[n-2],h[n-3],d[n-2],d[n-3]);

  // For a nondecreasing radius table, tiny negative endpoint slopes can occur
  // only from roundoff in nearly flat data.  Clamping those numerical remnants
  // to zero preserves the monotone physical contract without masking a truly
  // decreasing input table (which validate() rejects above).
  for (double& v:m) {
    const double scale=std::max(1.0,std::abs(v));
    if (v<0.0 && std::abs(v)<=64.0*std::numeric_limits<double>::epsilon()*scale) {
      v=0.0;
    }
  }
  return m;
}

inline State data_driven_state(const Config& c, double t_s) {
  State out;
  if (validate(c) != Status::Ok || !std::isfinite(t_s)) {
    out.status = Status::InvalidInput;
    return out;
  }

  const std::vector<double>& x=c.data_time_s;
  const std::vector<double>& y=c.data_radius_m;
  const std::vector<double> m=pchip_slopes(x,y);
  // Extrapolation is never implicit.  OUTSIDE_TIME is the science-safe default;
  // BALLISTIC extends the nearest endpoint with the PCHIP endpoint derivative.
  if (t_s < x.front()) {
    if (c.extrapolation==ExtrapolationPolicy::OutsideTime) {
      out.status=Status::OutsideTime;
      return out;
    }
    return domain_checked_state(
        std::fma(m.front(),t_s-x.front(),y.front()),m.front());
  }
  if (t_s > x.back()) {
    if (c.extrapolation==ExtrapolationPolicy::OutsideTime) {
      out.status=Status::OutsideTime;
      return out;
    }
    return domain_checked_state(
        std::fma(m.back(),t_s-x.back(),y.back()),m.back());
  }

  // Exact knot queries return the supplied radius exactly, while the speed is
  // the derivative of the same monotone interpolant at that knot.
  auto upper=std::lower_bound(x.begin(),x.end(),t_s);
  if (upper!=x.end() && *upper==t_s) {
    const std::size_t k=static_cast<std::size_t>(upper-x.begin());
    return domain_checked_state(y[k],m[k]);
  }

  const std::size_t i=static_cast<std::size_t>(upper-x.begin()-1);
  const double h=x[i+1]-x[i];
  const double q=(t_s-x[i])/h;
  const double q2=q*q, q3=q2*q;

  const double h00= 2.0*q3-3.0*q2+1.0;
  const double h10= q3-2.0*q2+q;
  const double h01=-2.0*q3+3.0*q2;
  const double h11= q3-q2;
  out.radius_m=h00*y[i]+h10*h*m[i]+h01*y[i+1]+h11*h*m[i+1];

  const double dh00=(6.0*q2-6.0*q)/h;
  const double dh10=(3.0*q2-4.0*q+1.0);
  const double dh01=(-6.0*q2+6.0*q)/h;
  const double dh11=(3.0*q2-2.0*q);
  out.speed_m_s=dh00*y[i]+dh10*m[i]+dh01*y[i+1]+dh11*m[i+1];

  // PCHIP should be nonnegative for nondecreasing data.  Remove only tiny
  // roundoff-level negatives; anything larger is a numerical failure and is
  // surfaced as INVALID_INPUT rather than silently repaired.
  const double speed_scale=std::max({1.0,std::abs(m[i]),std::abs(m[i+1])});
  if (out.speed_m_s<0.0) {
    if (std::abs(out.speed_m_s)<=128.0*std::numeric_limits<double>::epsilon()*speed_scale) {
      out.speed_m_s=0.0;
    } else {
      out.status=Status::InvalidInput;
      return out;
    }
  }

  return domain_checked_state(out.radius_m,out.speed_m_s);
}

inline State evaluate(const Config& c, double t_s) {
  if (c.mode==Mode::Ballistic) return ballistic_state(c,t_s);
  if (c.mode==Mode::DBM) return dbm_state(c,t_s);
  return data_driven_state(c,t_s);
}

} // namespace kinematics
} // namespace swcme
