#ifndef SWCME_SHOCK_HPP
#define SWCME_SHOCK_HPP

#include "swcme_constants.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>

// ============================================================================
// swcme_shock.hpp
// ----------------------------------------------------------------------------
// Shared ideal-MHD forward-shock jump solver used by both the 1-D and 3-D
// SWCME models.
//
// The solver deliberately separates three questions that were conflated in the
// legacy implementation:
//   (1) does a geometric CME/shock surface exist in this direction?;
//   (2) is the surface moving super-fast relative to the upstream normal flow,
//       so that a physical forward fast shock exists?;
//   (3) if a shock exists, what downstream primitive state satisfies the
//       ideal-MHD Rankine-Hugoniot conditions?
//
// A compression floor is NOT part of this calculation.  In particular, a
// sheath-shaping parameter may never manufacture a shock when M_fast <= 1.
//
// Frame convention
// ----------------
// The supplied upstream velocity is in the heliocentric/Sun frame.  The shock
// is transformed to a frame moving at V_sh,n along the outward unit normal n.
// The tangential velocity of the mathematical surface is a reparameterization
// freedom and is set to zero for the jump calculation.  Thus
//
//   u_1 = V_1 - V_sh,n n .
//
// For an outward forward shock u_1n < 0, while the positive inflow speed used
// for Mach-number classification is U_1n = -u_1n = V_sh,n - V_1.n.
//
// Solver reduction
// ----------------
// For a trial compression r=rho_2/rho_1, mass conservation fixes
// u_2n=u_1n/r.  Tangential momentum and tangential electric-field continuity
// form a 2x2 linear system for B_2t and u_2t.  Normal momentum then gives p_2.
// The remaining total-energy flux condition is one scalar nonlinear equation
// in r.  The trivial r=1 root is divided out numerically by solving
// energy_residual/(r-1)=0 on the physical interval 1<r<r_max, where
// r_max=(gamma+1)/(gamma-1).  A bracketed bisection search is used so the
// production result does not depend on an unconstrained Newton initial guess.
// ============================================================================

namespace swcme {
namespace shock {

using Vec3 = std::array<double, 3>;

struct PrimitiveState {
  double rho_kg_m3 = 0.0;
  double pressure_Pa = 0.0;
  Vec3 velocity_m_s{{0.0, 0.0, 0.0}};
  Vec3 magnetic_T{{0.0, 0.0, 0.0}};
};

struct JumpResult {
  // has_shock answers the physical shock-existence question.  A valid CME
  // surface may exist while has_shock=false if the normal relative speed is
  // sub-fast.
  bool has_shock = false;

  // solver_converged is true for a no-shock state (no nonlinear solve was
  // needed) and for a successfully solved physical fast shock.  It is false
  // only when a super-fast state was identified but no admissible RH root was
  // found.
  bool solver_converged = true;

  double compression = 1.0;
  double theta_Bn_rad = 0.0;
  double fast_speed_m_s = 0.0;
  double fast_mach = 0.0;
  double shock_normal_speed_m_s = 0.0;
  double upstream_inflow_normal_m_s = 0.0;
  int root_iterations = 0;

  PrimitiveState upstream;
  PrimitiveState downstream;

  // Dimensionless conservation diagnostics, evaluated independently after
  // reconstruction of the accepted downstream state.  They are useful both
  // for validation and for rejecting numerically ill-conditioned roots.
  double mass_residual = 0.0;
  double normal_B_residual = 0.0;
  double electric_residual = 0.0;
  double momentum_residual = 0.0;
  double energy_residual = 0.0;
  double entropy_ratio = 1.0;  // (p/rho^gamma)_2 / (p/rho^gamma)_1
};

namespace detail {

inline double dot(const Vec3& a, const Vec3& b) {
  return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

inline Vec3 add(const Vec3& a, const Vec3& b) {
  return {{a[0]+b[0], a[1]+b[1], a[2]+b[2]}};
}

inline Vec3 subtract(const Vec3& a, const Vec3& b) {
  return {{a[0]-b[0], a[1]-b[1], a[2]-b[2]}};
}

inline Vec3 scale(const Vec3& a, double s) {
  return {{s*a[0], s*a[1], s*a[2]}};
}

inline double norm(const Vec3& a) {
  return std::sqrt(std::max(0.0, dot(a,a)));
}

inline Vec3 normalized(const Vec3& a) {
  const double n = norm(a);
  return (n>0.0 && std::isfinite(n)) ? scale(a,1.0/n)
                                     : Vec3{{0.0,0.0,0.0}};
}

inline Vec3 tangential(const Vec3& a, const Vec3& n) {
  return subtract(a, scale(n, dot(a,n)));
}

inline Vec3 cross(const Vec3& a, const Vec3& b) {
  return {{a[1]*b[2]-a[2]*b[1],
           a[2]*b[0]-a[0]*b[2],
           a[0]*b[1]-a[1]*b[0]}};
}

inline double max_component_abs(const Vec3& a) {
  return std::max({std::abs(a[0]),std::abs(a[1]),std::abs(a[2])});
}

inline double total_energy_flux_normal(const PrimitiveState& state,
                                       const Vec3& shock_frame_velocity,
                                       const Vec3& n,
                                       double gamma) {
  const double mu0 = constants::VACUUM_PERMEABILITY_N_A2;
  const double un = dot(shock_frame_velocity,n);
  const double Bn = dot(state.magnetic_T,n);
  const double kinetic = 0.5*state.rho_kg_m3*dot(shock_frame_velocity,
                                                 shock_frame_velocity);
  const double enthalpy = gamma/(gamma-1.0)*state.pressure_Pa;
  const double magnetic = dot(state.magnetic_T,state.magnetic_T)/mu0;
  return un*(kinetic+enthalpy+magnetic)
       - (Bn/mu0)*dot(shock_frame_velocity,state.magnetic_T);
}

struct Candidate {
  bool valid = false;
  double energy_residual = 0.0;
  PrimitiveState downstream;
  Vec3 u2{{0.0,0.0,0.0}};
};

inline Candidate candidate_for_compression(const PrimitiveState& upstream,
                                           const Vec3& n,
                                           double Vsh_n,
                                           double gamma,
                                           double compression) {
  Candidate out;
  if (!(compression>1.0) || !(upstream.rho_kg_m3>0.0) ||
      !(upstream.pressure_Pa>0.0) || !(gamma>1.0)) return out;

  const double mu0 = constants::VACUUM_PERMEABILITY_N_A2;
  const Vec3 shock_velocity = scale(n,Vsh_n);
  const Vec3 u1 = subtract(upstream.velocity_m_s,shock_velocity);
  const double u1n = dot(u1,n);
  const Vec3 u1t = subtract(u1,scale(n,u1n));
  const double Bn = dot(upstream.magnetic_T,n);
  const Vec3 B1t = subtract(upstream.magnetic_T,scale(n,Bn));
  const double mass_flux = upstream.rho_kg_m3*u1n;

  const double rho2 = compression*upstream.rho_kg_m3;
  const double u2n = u1n/compression;

  // Tangential ideal-MHD jump conditions, solved component-by-component:
  //
  //   u2n B2t - Bn u2t = u1n B1t - Bn u1t             (E_t continuity)
  //   m u2t - (Bn/mu0) B2t = m u1t - (Bn/mu0) B1t     (tangential momentum)
  //
  // The vector coefficients are scalars, so the same 2x2 inverse applies to
  // both tangential basis directions.  A strictly parallel fixture has zero
  // tangential vectors; treating that exact limit analytically avoids the
  // removable Alfvénic determinant singularity without changing its physics.
  const Vec3 C = subtract(scale(B1t,u1n),scale(u1t,Bn));
  const Vec3 D = subtract(scale(u1t,mass_flux),scale(B1t,Bn/mu0));

  Vec3 B2t{{0.0,0.0,0.0}};
  Vec3 u2t{{0.0,0.0,0.0}};
  const double Bscale = std::max(norm(upstream.magnetic_T),1.0e-30);
  const double uscale = std::max(norm(u1),1.0);
  const bool strict_parallel = norm(B1t)<=64.0*std::numeric_limits<double>::epsilon()*Bscale &&
                               norm(u1t)<=64.0*std::numeric_limits<double>::epsilon()*uscale;

  if (!strict_parallel) {
    const double det = u2n*mass_flux - Bn*Bn/mu0;
    const double det_scale = std::max({std::abs(u2n*mass_flux),
                                       std::abs(Bn*Bn/mu0),1.0e-300});
    // Near this determinant the ideal-MHD tangential system is singular.  Do
    // not hide the conditioning by replacing it with a large finite number;
    // the bracket scan simply skips this mathematically singular candidate.
    if (!std::isfinite(det) || std::abs(det)<=1.0e-12*det_scale) return out;
    B2t = scale(add(scale(C,mass_flux),scale(D,Bn)),1.0/det);
    u2t = scale(add(scale(D,u2n),scale(C,Bn/mu0)),1.0/det);
  }

  const double B1t2 = dot(B1t,B1t);
  const double B2t2 = dot(B2t,B2t);
  const double p2 = upstream.rho_kg_m3*u1n*u1n + upstream.pressure_Pa
                  + B1t2/(2.0*mu0)
                  - rho2*u2n*u2n - B2t2/(2.0*mu0);
  if (!(p2>0.0) || !std::isfinite(p2)) return out;

  const Vec3 u2 = add(scale(n,u2n),u2t);
  PrimitiveState downstream;
  downstream.rho_kg_m3 = rho2;
  downstream.pressure_Pa = p2;
  downstream.velocity_m_s = add(u2,shock_velocity);
  downstream.magnetic_T = add(scale(n,Bn),B2t);

  const double F1 = total_energy_flux_normal(upstream,u1,n,gamma);
  const double F2 = total_energy_flux_normal(downstream,u2,n,gamma);
  if (!std::isfinite(F1) || !std::isfinite(F2)) return out;

  out.valid = true;
  out.energy_residual = F2-F1;
  out.downstream = downstream;
  out.u2 = u2;
  return out;
}

inline double fast_mode_speed(const PrimitiveState& upstream,
                              const Vec3& n,
                              double gamma) {
  if (!(upstream.rho_kg_m3>0.0) || !(upstream.pressure_Pa>=0.0) ||
      !(gamma>1.0)) return 0.0;
  const double B2 = dot(upstream.magnetic_T,upstream.magnetic_T);
  const double vA2 = B2/(constants::VACUUM_PERMEABILITY_N_A2*upstream.rho_kg_m3);
  const double cs2 = gamma*upstream.pressure_Pa/upstream.rho_kg_m3;
  double cosBn = 0.0;
  if (B2>0.0) cosBn = dot(upstream.magnetic_T,n)/std::sqrt(B2);
  cosBn = std::max(-1.0,std::min(1.0,cosBn));
  const double a = vA2+cs2;
  const double disc = std::max(0.0,a*a-4.0*vA2*cs2*cosBn*cosBn);
  return std::sqrt(std::max(0.0,0.5*(a+std::sqrt(disc))));
}

inline double relative_residual(double a, double b, double scale_floor) {
  return std::abs(a-b)/std::max({std::abs(a),std::abs(b),scale_floor});
}

}  // namespace detail

inline JumpResult solve_ideal_mhd_fast_shock(const PrimitiveState& upstream,
                                             const Vec3& outward_normal,
                                             double shock_normal_speed_m_s,
                                             double gamma) {
  JumpResult result;
  result.upstream = upstream;
  result.downstream = upstream;
  result.shock_normal_speed_m_s = shock_normal_speed_m_s;

  const Vec3 n = detail::normalized(outward_normal);
  if (detail::norm(n)==0.0 || !(upstream.rho_kg_m3>0.0) ||
      !(upstream.pressure_Pa>0.0) || !(gamma>1.0) ||
      !std::isfinite(shock_normal_speed_m_s)) {
    result.solver_converged = false;
    return result;
  }

  const double Bmag = detail::norm(upstream.magnetic_T);
  if (Bmag>0.0) {
    double c = std::abs(detail::dot(upstream.magnetic_T,n))/Bmag;
    c = std::max(0.0,std::min(1.0,c));
    result.theta_Bn_rad = std::acos(c);
  } else {
    // theta_Bn is undefined for B=0.  Zero is used only as a harmless scalar
    // diagnostic; shock existence still follows the hydrodynamic fast/sound
    // speed limit returned by fast_mode_speed().
    result.theta_Bn_rad = 0.0;
  }

  result.fast_speed_m_s = detail::fast_mode_speed(upstream,n,gamma);
  const double V1n = detail::dot(upstream.velocity_m_s,n);
  const double U1n = shock_normal_speed_m_s - V1n;
  result.upstream_inflow_normal_m_s = U1n;
  result.fast_mach = (result.fast_speed_m_s>0.0) ? U1n/result.fast_speed_m_s : 0.0;

  // The CME surface and the fast shock are deliberately separate concepts.
  // Equality at M_fast=1 is classified as no shock; no empirical compression
  // floor is permitted to override this criterion.
  const double mach_tol = 64.0*std::numeric_limits<double>::epsilon();
  if (!(U1n>0.0) || !(result.fast_mach>1.0+mach_tol)) {
    result.has_shock = false;
    result.compression = 1.0;
    return result;
  }

  result.has_shock = true;
  const double rmax = (gamma+1.0)/(gamma-1.0);
  const double eps_r = 1.0e-9;
  const int scan_points = 800;

  // Solve E(r)/(r-1)=0 instead of E(r)=0 so the always-present trivial
  // no-jump solution at r=1 cannot be selected as the physical fast-shock
  // branch.  Quadratic spacing clusters trial points near r=1, which is
  // important when M_fast is only slightly above unity.
  bool have_previous = false;
  double r_prev=0.0, q_prev=0.0;
  bool have_bracket = false;
  double r_lo=0.0, r_hi=0.0, q_lo=0.0;

  for (int i=0;i<=scan_points;++i) {
    const double x=static_cast<double>(i)/static_cast<double>(scan_points);
    const double delta=eps_r + (rmax-1.0-eps_r)*x*x;
    const double r=1.0+delta;
    const detail::Candidate c=detail::candidate_for_compression(upstream,n,
                                                                 shock_normal_speed_m_s,
                                                                 gamma,r);
    if (!c.valid) continue;
    const double q=c.energy_residual/(r-1.0);
    if (!std::isfinite(q)) continue;

    if (have_previous && (q==0.0 || q_prev==0.0 || (q_prev<0.0)!=(q<0.0))) {
      // Keep the outermost physical bracket found by the scan.  The trivial
      // branch was divided out; choosing the outermost remaining root selects
      // the compressive fast-shock branch for the supported solar-wind regime
      // and avoids lower-compression intermediate/switch branches when they
      // occur in pathological parameter sets.
      have_bracket=true;
      r_lo=r_prev; q_lo=q_prev;
      r_hi=r;
    }
    have_previous=true;
    r_prev=r; q_prev=q;
  }

  if (!have_bracket) {
    result.solver_converged=false;
    return result;
  }

  detail::Candidate accepted;
  double r_mid=0.0;
  for (int iter=0;iter<120;++iter) {
    r_mid=0.5*(r_lo+r_hi);
    const detail::Candidate c=detail::candidate_for_compression(upstream,n,
                                                                 shock_normal_speed_m_s,
                                                                 gamma,r_mid);
    if (!c.valid) {
      // A singular candidate inside the bracket is uncommon for the selected
      // fast branch.  Contract toward the side with the smaller interval; if
      // this persists the post-solve admissibility checks will reject it.
      r_hi=r_mid;
      continue;
    }
    const double q_mid=c.energy_residual/(r_mid-1.0);
    accepted=c;
    result.root_iterations=iter+1;

    if (std::abs(r_hi-r_lo)<=2.0e-13*std::max(1.0,r_mid)) break;
    if (q_mid==0.0) break;
    if ((q_lo<0.0)!=(q_mid<0.0)) {
      r_hi=r_mid;
    } else {
      r_lo=r_mid; q_lo=q_mid;
    }
  }

  const double compression=0.5*(r_lo+r_hi);
  accepted=detail::candidate_for_compression(upstream,n,shock_normal_speed_m_s,
                                               gamma,compression);
  if (!accepted.valid) {
    result.solver_converged=false;
    return result;
  }

  result.compression=compression;
  result.downstream=accepted.downstream;

  // ---------------- Independent conservation diagnostics -----------------
  // Recompute every ideal-MHD invariant from the final primitive states rather
  // than reusing candidate_for_compression() intermediates.  This makes these
  // residuals meaningful validation data rather than a restatement of the
  // scalar root function.
  const double mu0=constants::VACUUM_PERMEABILITY_N_A2;
  const Vec3 shock_velocity=detail::scale(n,shock_normal_speed_m_s);
  const Vec3 u1=detail::subtract(upstream.velocity_m_s,shock_velocity);
  const Vec3 u2=detail::subtract(result.downstream.velocity_m_s,shock_velocity);
  const double u1n=detail::dot(u1,n), u2n=detail::dot(u2,n);
  const double B1n=detail::dot(upstream.magnetic_T,n);
  const double B2n=detail::dot(result.downstream.magnetic_T,n);

  const double m1=upstream.rho_kg_m3*u1n;
  const double m2=result.downstream.rho_kg_m3*u2n;
  result.mass_residual=detail::relative_residual(m1,m2,1.0e-300);
  result.normal_B_residual=detail::relative_residual(B1n,B2n,
                                                     std::max(Bmag,1.0e-30));

  const Vec3 Et1=detail::tangential(detail::cross(u1,upstream.magnetic_T),n);
  const Vec3 Et2=detail::tangential(detail::cross(u2,result.downstream.magnetic_T),n);
  result.electric_residual=detail::norm(detail::subtract(Et1,Et2)) /
      std::max({detail::norm(Et1),detail::norm(Et2),1.0e-30});

  auto momentum_flux=[&](const PrimitiveState& s,const Vec3& u)->Vec3 {
    const double un=detail::dot(u,n);
    const double Bn=detail::dot(s.magnetic_T,n);
    const double magnetic_pressure=detail::dot(s.magnetic_T,s.magnetic_T)/(2.0*mu0);
    Vec3 f=detail::scale(u,s.rho_kg_m3*un);
    f=detail::add(f,detail::scale(n,s.pressure_Pa+magnetic_pressure));
    f=detail::subtract(f,detail::scale(s.magnetic_T,Bn/mu0));
    return f;
  };
  const Vec3 mom1=momentum_flux(upstream,u1);
  const Vec3 mom2=momentum_flux(result.downstream,u2);
  result.momentum_residual=detail::norm(detail::subtract(mom1,mom2)) /
      std::max({detail::norm(mom1),detail::norm(mom2),1.0e-30});

  const double F1=detail::total_energy_flux_normal(upstream,u1,n,gamma);
  const double F2=detail::total_energy_flux_normal(result.downstream,u2,n,gamma);
  result.energy_residual=detail::relative_residual(F1,F2,1.0e-300);

  const double K1=upstream.pressure_Pa/std::pow(upstream.rho_kg_m3,gamma);
  const double K2=result.downstream.pressure_Pa/
                  std::pow(result.downstream.rho_kg_m3,gamma);
  result.entropy_ratio=(K1>0.0)? K2/K1 : 0.0;

  // A mathematical root is not accepted unless it is a compressive,
  // entropy-increasing forward shock and satisfies the conserved quantities at
  // a level tighter than the public validation thresholds.  This prevents a
  // formally converged but wrong branch from reaching SEP source physics.
  const bool admissible = result.compression>1.0 && result.compression<=rmax*(1.0+1e-12) &&
                          result.downstream.rho_kg_m3>upstream.rho_kg_m3 &&
                          result.downstream.pressure_Pa>0.0 &&
                          result.entropy_ratio>=1.0-1.0e-10 &&
                          result.mass_residual<=1.0e-9 &&
                          result.normal_B_residual<=1.0e-10 &&
                          result.electric_residual<=1.0e-8 &&
                          result.momentum_residual<=1.0e-8 &&
                          result.energy_residual<=1.0e-8;
  result.solver_converged=admissible;
  return result;
}

}  // namespace shock
}  // namespace swcme

#endif  // SWCME_SHOCK_HPP
