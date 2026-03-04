//======================================================================================
// GridlessParticleMovers.cpp
//======================================================================================
// See GridlessParticleMovers.h for detailed documentation and design notes.
//======================================================================================

#include "GridlessParticleMovers.h"

//------------------------------------------------------------------------------
// Classic relativistic Boris pusher (magnetic field only)
//------------------------------------------------------------------------------
// This is the same algorithm historically used in the gridless cutoff and density
// solvers, moved here for reuse.
//
// Key details:
//   - We integrate momentum p (not velocity).
//   - Relativistic gamma is computed from |p|.
//   - In the absence of E, the Boris push reduces to a pure rotation of p around B.
//   - Position is advanced using the updated momentum.
void BorisStep(V3& x, V3& p, double q_C, double m0_kg, double dt,
               const IGridlessFieldEvaluator& field) {
  const double p2 = dot(p,p);
  const double mc = m0_kg*SpeedOfLight;
  const double gamma = std::sqrt(1.0 + p2/(mc*mc));

  V3 B; field.GetB_T(x,B);

  // Boris rotation parameters.
  // t = (q dt / (2 gamma m)) * B
  V3 t = mul((q_C*dt)/(2.0*gamma*m0_kg), B);
  const double t2 = dot(t,t);
  V3 s = mul(2.0/(1.0+t2), t);

  // p' = p + p x t
  V3 p_prime = add(p, cross(p, t));
  // p+ = p + p' x s
  V3 p_plus  = add(p, cross(p_prime, s));
  p = p_plus;

  // Update position using v = p / (gamma m)
  const double p2n = dot(p,p);
  const double gamman = std::sqrt(1.0 + p2n/(mc*mc));
  V3 vnew = mul(1.0/(gamman*m0_kg), p);
  x = add(x, mul(dt, vnew));
}

//------------------------------------------------------------------------------
// Higher-accuracy Boris variant: midpoint B sampling
//------------------------------------------------------------------------------
// Motivation:
//   In strongly varying magnetic fields, sampling B only at x_n introduces a
//   first-order error in how the rotation axis is chosen over the step.
//   A simple and effective improvement is to sample B at an approximate midpoint
//   position x_{n+1/2}.
//
// Algorithm sketch:
//   1) compute v_n from p_n
//   2) x_half = x_n + 0.5*dt*v_n
//   3) evaluate B_half = B(x_half)
//   4) perform Boris rotation using B_half
//   5) advance position with the updated momentum
//
// Notes:
//   - Still no E-field.
//   - Cost: one extra B evaluation per step.
//   - Benefit: reduced phase/trajectory error for a given dt.
void BorisMidpointStep(V3& x, V3& p, double q_C, double m0_kg, double dt,
                       const IGridlessFieldEvaluator& field) {
  const double p2 = dot(p,p);
  const double mc = m0_kg*SpeedOfLight;
  const double gamma = std::sqrt(1.0 + p2/(mc*mc));

  // v_n = p_n/(gamma m)
  const V3 v_n = mul(1.0/(gamma*m0_kg), p);
  const V3 x_half = add(x, mul(0.5*dt, v_n));

  // B sampled at midpoint position
  V3 B_half; field.GetB_T(x_half, B_half);

  V3 t = mul((q_C*dt)/(2.0*gamma*m0_kg), B_half);
  const double t2 = dot(t,t);
  V3 s = mul(2.0/(1.0+t2), t);

  V3 p_prime = add(p, cross(p, t));
  V3 p_plus  = add(p, cross(p_prime, s));
  p = p_plus;

  const double p2n = dot(p,p);
  const double gamman = std::sqrt(1.0 + p2n/(mc*mc));
  V3 vnew = mul(1.0/(gamman*m0_kg), p);
  x = add(x, mul(dt, vnew));
}

//------------------------------------------------------------------------------
// Dispatcher
//------------------------------------------------------------------------------
void StepParticle(MoverType mover,
                  V3& x_m, V3& p_SI,
                  double q_C, double m0_kg,
                  double dt,
                  const IGridlessFieldEvaluator& field) {
  switch (mover) {
    case MoverType::BORIS:
      BorisStep(x_m, p_SI, q_C, m0_kg, dt, field);
      break;
    case MoverType::BORIS_MIDPOINT:
      BorisMidpointStep(x_m, p_SI, q_C, m0_kg, dt, field);
      break;
    default:
      // Defensive default: use classic Boris.
      BorisStep(x_m, p_SI, q_C, m0_kg, dt, field);
      break;
  }
}
