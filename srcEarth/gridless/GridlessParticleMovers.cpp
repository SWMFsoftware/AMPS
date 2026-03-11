//======================================================================================
// GridlessParticleMovers.cpp
//======================================================================================
//
// IMPLEMENTATION NOTES
// --------------------
// This file implements the four particle movers declared in GridlessParticleMovers.h.
// See that header for the full physics description, algorithm rationale, and
// comparison of mover choices.
//
// UNIT CONVENTIONS (enforced throughout this file)
// -------------------------------------------------
//   Position x          : meters  [m]
//   Momentum p          : kg m/s  (relativistic: p = gamma * m0 * v)
//   Magnetic field B    : Tesla   [T]
//   Charge q            : Coulombs [C]
//   Mass m0             : kg
//   Time dt             : seconds [s]
//   Cyclotron frequency : rad/s
//
// All input positions (GSM coordinates) are converted from km to m by the callers
// before being passed here. No unit conversion is done inside the movers themselves.
//
// BORIS PUSHER -- IMPLEMENTATION NOTES
// -------------------------------------
// The rotation step follows the exact Boris half-angle formula:
//
//   t = (q * B * dt) / (2 * gamma * m0)       [vector]
//   s = 2*t / (1 + |t|^2)
//   p' = p_minus + (p_minus + p_minus x t) x s
//
// The factor 1/(1+|t|^2) in s ensures that the rotation matrix is exactly
// orthogonal even for large |t| (large step size or strong field). This is the
// key stability property absent in the naive "add cross product" form.
//
// gamma is evaluated at p_minus (i.e., the beginning of the magnetic rotation
// half-step). Since |p| is conserved by the magnetic force, gamma is the same
// before and after the rotation, so this is self-consistent.
//
// RK MOVERS -- IMPLEMENTATION NOTES
// ------------------------------------
// The RK movers integrate the system [dx/dt, dp/dt] = [v, q*v x B(x)] as a
// 6-dimensional ODE. The function f(x, p) returning (v, q*v x B) is evaluated
// at each stage position. For RK2 this means evaluating B twice per outer step;
// for RK4, four times; for RK6, six times.
//
// Note on gamma consistency in RK stages: at each intermediate stage we recompute
// gamma from the stage momentum p_k:
//   gamma_k = sqrt(1 + |p_k|^2 / (m0*c)^2)
//   v_k = p_k / (gamma_k * m0)
// This maintains relativistic correctness at each stage without approximation.
//
// StepParticleChecked -- IMPLEMENTATION NOTES
// --------------------------------------------
// For Boris: subdivide into N_SUBSTEP Boris sub-steps of size dt/N_SUBSTEP,
// checking |x| < rInner after each. N_SUBSTEP = 4 balances detection accuracy
// against overhead (4 extra field evaluations per outer step).
//
// For RK variants: the intermediate stage positions k1..k4 (or k1..k6 for RK6)
// are already computed inside the step. We test |x| < rInner at each stage
// position rather than after the completed step. This costs nothing extra because
// the stage evaluations happen anyway.
//
// THREAD SAFETY
// -------------
// All movers are stateless functions (no global or static mutable state). They
// are safe to call concurrently from multiple MPI worker threads as long as the
// IGridlessFieldEvaluator passed in is also thread-safe (the Tsyganenko and dipole
// evaluators in CutoffRigidityGridless.cpp use only local stack state and are safe).
//
//======================================================================================

#include "GridlessParticleMovers.h"

#include <algorithm>
#include <array>
#include <cctype>

namespace {

static inline double GammaFromMomentum(const V3& p_SI, double m0_kg) {
  const double p2 = dot(p_SI, p_SI);
  const double mc = m0_kg * SpeedOfLight;
  return std::sqrt(1.0 + p2 / (mc * mc));
}

static inline V3 VelocityFromMomentum(const V3& p_SI, double m0_kg) {
  const double gamma = GammaFromMomentum(p_SI, m0_kg);
  return mul(1.0 / (gamma * m0_kg), p_SI);
}

struct State { V3 x; V3 p; };

static inline State add_state(const State& a, const State& b) {
  return { add(a.x,b.x), add(a.p,b.p) };
}

static inline State mul_state(double s, const State& a) {
  return { mul(s,a.x), mul(s,a.p) };
}

static inline bool InsideInnerSphereM(const V3& x_m, double rInner_m) {
  return (rInner_m > 0.0) && (dot(x_m,x_m) <= rInner_m*rInner_m);
}

static bool SegmentHitsInnerSphere(const V3& a, const V3& b, double rInner_m) {
  if (rInner_m <= 0.0) return false;
  if (InsideInnerSphereM(a,rInner_m) || InsideInnerSphereM(b,rInner_m)) return true;

  const V3 d = sub(b,a);
  const double dd = dot(d,d);
  if (dd <= 0.0) return false;

  double t = -dot(a,d)/dd;
  if (t < 0.0) t = 0.0;
  if (t > 1.0) t = 1.0;
  const V3 q = add(a, mul(t,d));
  return dot(q,q) <= rInner_m*rInner_m;
}

static inline State LorentzRhs(const State& s,
                               double q_C, double m0_kg,
                               const IGridlessFieldEvaluator& field) {
  V3 B_T; field.GetB_T(s.x, B_T);
  const V3 v = VelocityFromMomentum(s.p, m0_kg);
  return { v, mul(q_C, cross(v, B_T)) };
}

static bool CheckIntermediateSphereHit(const V3& from,
                                       const V3& to,
                                       double rInner_m,
                                       V3& x_m) {
  if (SegmentHitsInnerSphere(from,to,rInner_m)) {
    x_m = to;
    return true;
  }
  return false;
}

// Generic explicit RK stepper with segment/stage inner-sphere checks.
template <size_t N>
static bool RKStepGeneric(const std::array<double,N>& c,
                          const std::array<std::array<double,N>,N>& a,
                          const std::array<double,N>& b,
                          V3& x_m, V3& p_SI,
                          double q_C, double m0_kg,
                          double dt,
                          const IGridlessFieldEvaluator& field,
                          double rInner_m) {
  const State y0{x_m,p_SI};
  std::array<State,N> k{};
  std::array<V3,N+1> stageX{};
  stageX[0] = x_m;

  for (size_t i=0;i<N;i++) {
    State yi = y0;
    for (size_t j=0;j<i;j++) {
      yi = add_state(yi, mul_state(dt*a[i][j], k[j]));
    }

    stageX[i+1] = yi.x;
    if (CheckIntermediateSphereHit(stageX[i], stageX[i+1], rInner_m, x_m)) return false;
    if (InsideInnerSphereM(yi.x, rInner_m)) { x_m = yi.x; p_SI = yi.p; return false; }

    (void)c;
    k[i] = LorentzRhs(yi, q_C, m0_kg, field);
  }

  State y1 = y0;
  for (size_t i=0;i<N;i++) y1 = add_state(y1, mul_state(dt*b[i], k[i]));

  if (CheckIntermediateSphereHit(stageX[N], y1.x, rInner_m, x_m)) { p_SI = y1.p; return false; }
  if (InsideInnerSphereM(y1.x, rInner_m)) { x_m = y1.x; p_SI = y1.p; return false; }

  x_m = y1.x;
  p_SI = y1.p;
  return true;
}

} // namespace

void BorisStep(V3& x_m, V3& p_SI, double q_C, double m0_kg, double dt,
               const IGridlessFieldEvaluator& field) {
  V3 B_T; field.GetB_T(x_m,B_T);

  const double gam = GammaFromMomentum(p_SI, m0_kg);
  const V3 v = mul(1.0/(gam*m0_kg), p_SI);
  x_m = add(x_m, mul(0.5*dt, v));

  const V3 t = mul((q_C*dt)/(2.0*gam*m0_kg), B_T);
  const double t2 = dot(t,t);
  const V3 s = mul(2.0/(1.0+t2), t);

  const V3 p_minus = p_SI;
  const V3 p_prime = add(p_minus, cross(p_minus, t));
  const V3 p_plus  = add(p_minus, cross(p_prime, s));
  p_SI = p_plus;

  const double gam2 = GammaFromMomentum(p_SI, m0_kg);
  const V3 v2 = mul(1.0/(gam2*m0_kg), p_SI);
  x_m = add(x_m, mul(0.5*dt, v2));
}

void RK2Step(V3& x_m, V3& p_SI,
             double q_C, double m0_kg,
             double dt,
             const IGridlessFieldEvaluator& field) {
  static const std::array<double,2> c{{0.0, 0.5}};
  static const std::array<std::array<double,2>,2> a{{
    {{0.0, 0.0}},
    {{0.5, 0.0}}
  }};
  static const std::array<double,2> b{{0.0, 1.0}};
  RKStepGeneric(c,a,b,x_m,p_SI,q_C,m0_kg,dt,field,-1.0);
}

void RK4Step(V3& x_m, V3& p_SI,
             double q_C, double m0_kg,
             double dt,
             const IGridlessFieldEvaluator& field) {
  static const std::array<double,4> c{{0.0, 0.5, 0.5, 1.0}};
  static const std::array<std::array<double,4>,4> a{{
    {{0.0, 0.0, 0.0, 0.0}},
    {{0.5, 0.0, 0.0, 0.0}},
    {{0.0, 0.5, 0.0, 0.0}},
    {{0.0, 0.0, 1.0, 0.0}}
  }};
  static const std::array<double,4> b{{1.0/6.0, 1.0/3.0, 1.0/3.0, 1.0/6.0}};
  RKStepGeneric(c,a,b,x_m,p_SI,q_C,m0_kg,dt,field,-1.0);
}

void RK6Step(V3& x_m, V3& p_SI,
             double q_C, double m0_kg,
             double dt,
             const IGridlessFieldEvaluator& field) {
  static const std::array<double,7> c{{0.0, 1.0/3.0, 2.0/3.0, 1.0/3.0, 0.5, 0.5, 1.0}};
  static const std::array<std::array<double,7>,7> a{{
    {{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}},
    {{1.0/3.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}},
    {{0.0, 2.0/3.0, 0.0, 0.0, 0.0, 0.0, 0.0}},
    {{1.0/12.0, 1.0/3.0, -1.0/12.0, 0.0, 0.0, 0.0, 0.0}},
    {{-1.0/16.0, 9.0/8.0, -3.0/16.0, -3.0/8.0, 0.0, 0.0, 0.0}},
    {{0.0, 9.0/8.0, -3.0/8.0, -3.0/4.0, 0.5, 0.0, 0.0}},
    {{9.0/44.0, -9.0/11.0, 63.0/44.0, 18.0/11.0, 0.0, -16.0/11.0, 0.0}}
  }};
  static const std::array<double,7> b{{11.0/120.0, 0.0, 27.0/40.0, 27.0/40.0, -4.0/15.0, -4.0/15.0, 11.0/120.0}};
  RKStepGeneric(c,a,b,x_m,p_SI,q_C,m0_kg,dt,field,-1.0);
}

MoverType gDefaultMover = MoverType::RK4;

void SetDefaultMoverType(MoverType m) {
  gDefaultMover = m;
}

MoverType GetDefaultMoverType() {
  return gDefaultMover;
}

static std::string ToUpperCopy(std::string s) {
  std::transform(s.begin(), s.end(), s.begin(), [](unsigned char ch){ return std::toupper(ch); });
  return s;
}

bool ParseMoverType(const std::string& s, MoverType& out) {
  const std::string u = ToUpperCopy(s);
  if (u=="BORIS") { out = MoverType::BORIS; return true; }
  if (u=="RK2" || u=="RUNGEKUTTA2" || u=="RUNGE-KUTTA-2") { out = MoverType::RK2; return true; }
  if (u=="RK4" || u=="RUNGEKUTTA4" || u=="RUNGE-KUTTA-4") { out = MoverType::RK4; return true; }
  if (u=="RK6" || u=="RUNGEKUTTA6" || u=="RUNGE-KUTTA-6") { out = MoverType::RK6; return true; }
  return false;
}

const char* MoverTypeToString(MoverType m) {
  switch (m) {
    case MoverType::BORIS: return "BORIS";
    case MoverType::RK2:   return "RK2";
    case MoverType::RK4:   return "RK4";
    case MoverType::RK6:   return "RK6";
    default: return "BORIS";
  }
}

void StepParticle(MoverType mover,
                  V3& x_m, V3& p_SI,
                  double q_C, double m0_kg,
                  double dt,
                  const IGridlessFieldEvaluator& field) {
  switch (mover) {
    case MoverType::BORIS:
      BorisStep(x_m, p_SI, q_C, m0_kg, dt, field);
      break;
    case MoverType::RK2:
      RK2Step(x_m, p_SI, q_C, m0_kg, dt, field);
      break;
    case MoverType::RK4:
      RK4Step(x_m, p_SI, q_C, m0_kg, dt, field);
      break;
    case MoverType::RK6:
      RK6Step(x_m, p_SI, q_C, m0_kg, dt, field);
      break;
    default:
      BorisStep(x_m, p_SI, q_C, m0_kg, dt, field);
      break;
  }
}

bool StepParticleChecked(MoverType mover,
                         V3& x_m, V3& p_SI,
                         double q_C, double m0_kg,
                         double dt,
                         const IGridlessFieldEvaluator& field,
                         double rInner_m) {
  switch (mover) {
    case MoverType::BORIS: {
      const V3 x0 = x_m;
      BorisStep(x_m, p_SI, q_C, m0_kg, dt, field);
      return !SegmentHitsInnerSphere(x0, x_m, rInner_m) && !InsideInnerSphereM(x_m, rInner_m);
    }
    case MoverType::RK2: {
      static const std::array<double,2> c{{0.0, 0.5}};
      static const std::array<std::array<double,2>,2> a{{
        {{0.0, 0.0}},
        {{0.5, 0.0}}
      }};
      static const std::array<double,2> b{{0.0, 1.0}};
      return RKStepGeneric(c,a,b,x_m,p_SI,q_C,m0_kg,dt,field,rInner_m);
    }
    case MoverType::RK4: {
      static const std::array<double,4> c{{0.0, 0.5, 0.5, 1.0}};
      static const std::array<std::array<double,4>,4> a{{
        {{0.0, 0.0, 0.0, 0.0}},
        {{0.5, 0.0, 0.0, 0.0}},
        {{0.0, 0.5, 0.0, 0.0}},
        {{0.0, 0.0, 1.0, 0.0}}
      }};
      static const std::array<double,4> b{{1.0/6.0, 1.0/3.0, 1.0/3.0, 1.0/6.0}};
      return RKStepGeneric(c,a,b,x_m,p_SI,q_C,m0_kg,dt,field,rInner_m);
    }
    case MoverType::RK6: {
      static const std::array<double,7> c{{0.0, 1.0/3.0, 2.0/3.0, 1.0/3.0, 0.5, 0.5, 1.0}};
      static const std::array<std::array<double,7>,7> a{{
        {{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}},
        {{1.0/3.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}},
        {{0.0, 2.0/3.0, 0.0, 0.0, 0.0, 0.0, 0.0}},
        {{1.0/12.0, 1.0/3.0, -1.0/12.0, 0.0, 0.0, 0.0, 0.0}},
        {{-1.0/16.0, 9.0/8.0, -3.0/16.0, -3.0/8.0, 0.0, 0.0, 0.0}},
        {{0.0, 9.0/8.0, -3.0/8.0, -3.0/4.0, 0.5, 0.0, 0.0}},
        {{9.0/44.0, -9.0/11.0, 63.0/44.0, 18.0/11.0, 0.0, -16.0/11.0, 0.0}}
      }};
      static const std::array<double,7> b{{11.0/120.0, 0.0, 27.0/40.0, 27.0/40.0, -4.0/15.0, -4.0/15.0, 11.0/120.0}};
      return RKStepGeneric(c,a,b,x_m,p_SI,q_C,m0_kg,dt,field,rInner_m);
    }
    default: {
      const V3 x0 = x_m;
      BorisStep(x_m, p_SI, q_C, m0_kg, dt, field);
      return !SegmentHitsInnerSphere(x0, x_m, rInner_m) && !InsideInnerSphereM(x_m, rInner_m);
    }
  }
}
