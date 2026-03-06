//======================================================================================
// GridlessParticleMovers.cpp
//======================================================================================
// See GridlessParticleMovers.h for detailed documentation and design notes.
//======================================================================================

#include "GridlessParticleMovers.h"

#include <algorithm>
#include <cctype>

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
void BorisStep(V3& x_m, V3& p_SI, double q_C, double m0_kg, double dt,
               const IGridlessFieldEvaluator& field) {

  V3 B_T; field.GetB_T(x_m,B_T);

  auto TcppGammaFromMomentum = [](const V3& p_SI, double m0_kg) -> double {
  const double p2 = dot(p_SI, p_SI);
  const double mc = m0_kg * SpeedOfLight;
  return std::sqrt(1.0 + p2 / (mc * mc));
  };
  

  const double gam = TcppGammaFromMomentum(p_SI, m0_kg);
  const V3 v = mul(1.0/(gam*m0_kg), p_SI);

  // half drift
  x_m = add(x_m, mul(0.5*dt, v));

  // Boris rotation in momentum
  const V3 t = mul((q_C*dt)/(2.0*gam*m0_kg), B_T);
  const double t2 = dot(t,t);
  const V3 s = mul(2.0/(1.0+t2), t);

  const V3 p_minus = p_SI;
  const V3 p_prime = add(p_minus, cross(p_minus, t));
  const V3 p_plus  = add(p_minus, cross(p_prime, s));
  p_SI = p_plus;

  // second half drift with updated momentum
  const double gam2 = TcppGammaFromMomentum(p_SI, m0_kg);
  const V3 v2 = mul(1.0/(gam2*m0_kg), p_SI);
  x_m = add(x_m, mul(0.5*dt, v2));
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
// Process-wide default mover
//------------------------------------------------------------------------------
// See header for intended precedence and usage.
MoverType gDefaultMover = MoverType::BORIS;

void SetDefaultMoverType(MoverType m) {
  gDefaultMover = m;
}

MoverType GetDefaultMoverType() {
  return gDefaultMover;
}

//------------------------------------------------------------------------------
// Parsing / printing helpers
//------------------------------------------------------------------------------
static std::string ToUpperCopy(std::string s) {
  std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c){return std::toupper(c);});
  return s;
}

bool ParseMoverType(const std::string& s, MoverType& out) {
  const std::string u = ToUpperCopy(s);
  if (u=="BORIS") { out = MoverType::BORIS; return true; }
  if (u=="BORIS_MIDPOINT" || u=="MIDPOINT" || u=="BORIS-MIDPOINT") {
    out = MoverType::BORIS_MIDPOINT; return true;
  }
  return false;
}

const char* MoverTypeToString(MoverType m) {
  switch (m) {
    case MoverType::BORIS: return "BORIS";
    case MoverType::BORIS_MIDPOINT: return "BORIS_MIDPOINT";
    default: return "BORIS";
  }
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
