//======================================================================================
// CutoffRigidityGridless.cpp
//======================================================================================
// IMPLEMENTATION SUMMARY
//
// 1) Direct access to Tsyganenko models
//    The standard AMPS interfaces (src/interface/T96Interface.* and
//    src/interface/T05Interface.*) guard model calls by compile-time
//    _PIC_COUPLER_MODE_. For the gridless tool we need to choose T96/T05 at
//    runtime. Therefore we call the underlying Fortran entry points directly:
//      - t96_01_ (T96)
//      - t04_s_  (T05)
//
//    Both models are called with X,Y,Z in GSM [Re] and return BX,BY,BZ in nT.
//    We convert to Tesla and add the internal IGRF field obtained from
//    Geopack::IGRF::GetMagneticField().
//
// 2) Particle tracing
//    - Relativistic Boris pusher (magnetic field only).
//    - State variables are kept in SI units:
//        x [m], p=gamma m v [kg m/s]
//    - Rigidity R [GV] is converted to momentum magnitude:
//        p = (R * 1e9 * |q|) / c
//
// 3) Classification
//    Allowed if the particle escapes the rectangular domain boundary without
//    hitting the inner loss sphere. If it times out, we classify as forbidden
//    (conservative).
//
//======================================================================================

#include "pic.h"
#include "CutoffRigidityGridless.h"
#include "DipoleInterface.h" 

//--------------------------------------------------------------------------------------
// Shared particle movers (extracted)
//--------------------------------------------------------------------------------------
// Per user request, particle movers used by gridless cutoff and gridless density/spectrum
// solvers must be identical. We therefore include the shared mover module.
//
// This header provides:
//   - struct V3 and common vector ops (add/mul/dot/cross/norm/unit)
//   - interface IGridlessFieldEvaluator
//   - enum class MoverType { BORIS, BORIS_MIDPOINT }
//   - ::StepParticle(...) dispatcher and individual mover implementations
//
// IMPORTANT:
//   This file historically contained local copies of V3 and BorisStep. Those historical
//   implementations are archived below inside a large comment block (comments preserved) so future
//   developers can compare behavior if needed.
#include "GridlessParticleMovers.h"

#include <cstdio>
#include <cmath>
#include <stdexcept>
#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <numeric>
#include <mpi.h> 

//--------------------------------------------------------------------------------------
// OPTIONAL SPICE SUPPORT (FRAME TRANSFORMS)
//--------------------------------------------------------------------------------------
// This solver keeps the actual particle tracing in GSM (because Tsyganenko models are
// defined in GSM), but for *visualization / directional sampling maps* it is often
// desirable to label directions in another global frame (e.g., SM).
//
// IMPORTANT DESIGN CHOICE:
//   - We intentionally DO NOT use Geopack for coordinate transformations here.
//   - If a direction-frame transform is needed, we do it via SPICE pxform.
//   - Geopack remains IGRF-only (see cFieldEvaluator).
//
// To enable SPICE transforms, compile with -DAMPS_USE_SPICE and link CSPICE.
// The runtime must furnish kernels that define the requested frames.
// Per user request:
//   - GSM frame name in SPICE is "GSM".
//   - Solar Magnetic frame name is "SM".
//   - If you later want Earth-fixed labeling, you may use "GCE" (if present in your FK).
#ifndef _NO_SPICE_CALLS_
#include "SpiceUsr.h"
#endif

#include "constants.h"
#include "constants.PlanetaryData.h"
#include "GeopackInterface.h"
#include "DipoleInterface.h"

extern "C" {
  void t96_01_(int*,double*,double*,double*,double*,double*,double*,double*,double*);
  void t04_s_(int*,double*,double*,double*,double*,double*,double*,double*,double*);
}

// Compute Stormer vertical-cutoff coefficient R0(M) in GV, using the *same* dipole moment
// normalization as the dipole B-field implementation.
//
// Dipole B uses:  B = (mu0/4pi) * ( 3 r (m·r)/r^5 - m/r^3 )
// with m = M * m_hat, and M = momentScale_Me * M_E_Am2.
//
// Stormer vertical cutoff: Rc = R0(M) * cos^4(lambda) / r^2, with
// R0(M) [GV] = 0.299792458 * (1/4) * (mu0/4pi) * M / Re^2
static inline double StormerVerticalCoeff_GV(double momentScale_Me,
                                             double Re_m /* _EARTH__RADIUS_ */) {
  constexpr double mu0_over_4pi = 1.0e-7;          // SI
  constexpr double c_to_GV_per_Tm = 0.299792458;   // (c/1e9) converts T·m -> GV
  constexpr double M_E_Am2 = Earth::GridlessMode::Dipole::M_E_Am2;              // your code constant (modern-era representative)

  const double M = momentScale_Me * M_E_Am2;       // A·m^2
  const double Brho_Tm = (mu0_over_4pi * M) / (Re_m * Re_m); // T·m scale at equator, r=Re
  const double R0_GV = c_to_GV_per_Tm * 0.25 * Brho_Tm;      // Stormer vertical factor 1/4
  return R0_GV;
}



namespace {

/*
//--------------------------------------------------------------------------------------
// Legacy local V3 + vector ops (kept for reference; DO NOT DELETE)
//--------------------------------------------------------------------------------------
// NOTE:
//   These definitions were historically embedded in CutoffRigidityGridless.cpp.
//   They have been superseded by the shared module GridlessParticleMovers.h so that
//   cutoff rigidity and density/spectrum tools cannot diverge.
//   We keep the original code here (disabled) because it is valuable for regression
//   archaeology and for quickly comparing behavior if a future change raises questions.
struct V3 { double x,y,z; };
static inline V3 add(const V3&a,const V3&b){return {a.x+b.x,a.y+b.y,a.z+b.z};}
static inline V3 mul(double s,const V3&a){return {s*a.x,s*a.y,s*a.z};}
static inline double dot(const V3&a,const V3&b){return a.x*b.x+a.y*b.y+a.z*b.z;}
static inline V3 cross(const V3&a,const V3&b){return {a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x};}
static inline double norm(const V3&a){return std::sqrt(dot(a,a));}
static inline V3 unit(const V3&a){double n=norm(a); return (n>0)?mul(1.0/n,a):V3{0,0,0};}

*/

//--------------------------------------------------------------------------------------
// Small 3x3 rotation helper
//--------------------------------------------------------------------------------------
struct Mat3 { double a[3][3]; };
static inline Mat3 Identity3() {
  Mat3 R{};
  R.a[0][0]=1.0; R.a[0][1]=0.0; R.a[0][2]=0.0;
  R.a[1][0]=0.0; R.a[1][1]=1.0; R.a[1][2]=0.0;
  R.a[2][0]=0.0; R.a[2][1]=0.0; R.a[2][2]=1.0;
  return R;
}

static inline V3 Apply(const Mat3& R, const V3& v) {
  return {
    R.a[0][0]*v.x + R.a[0][1]*v.y + R.a[0][2]*v.z,
    R.a[1][0]*v.x + R.a[1][1]*v.y + R.a[1][2]*v.z,
    R.a[2][0]*v.x + R.a[2][1]*v.y + R.a[2][2]*v.z
  };
}

//--------------------------------------------------------------------------------------
// SPICE frame transform helper
//--------------------------------------------------------------------------------------
// Returns a rotation matrix that maps vectors from `fromFrame` to `toFrame`.
// The matrix is time-dependent; we use prm.field.epoch (a SPICE-parsable time string).
//
// Behavior when SPICE is not available or the transform fails:
//   - ok=false
//   - returns identity (so the program can still run, but the labeled frame is not real)
//
// NOTE:
//   This helper is used ONLY for direction labeling / sampling maps. The physics tracing
//   remains in GSM.
static inline Mat3 GetSpiceRotationOrIdentity(const char* fromFrame,
                                              const char* toFrame,
                                              const std::string& epoch,
                                              bool& ok) {
  ok = false;
  Mat3 R = Identity3();

#ifndef _NO_SPICE_CALLS_
  try {
    SpiceDouble et = 0.0;
    str2et_c(epoch.c_str(), &et);

    SpiceDouble m[3][3];
    pxform_c(fromFrame, toFrame, et, m);

    for (int i=0;i<3;i++) for (int j=0;j<3;j++) R.a[i][j] = m[i][j];
    ok = true;
  }
  catch (...) {
    // CSPICE is C, but some builds may still throw through wrappers; keep this as
    // a defensive guard. If it fails we fall back to identity.
    ok = false;
  }
#else
  (void)fromFrame; (void)toFrame; (void)epoch;
#endif

  return R;
}

static inline double MomentumFromKineticEnergy_MeV(double E_MeV,double m0_kg) {
  const double E_J = E_MeV * 1.0e6 * ElectronCharge;
  return Relativistic::Energy2Momentum(E_J,m0_kg);
}

static inline double KineticEnergyFromMomentum_MeV(double p,double m0_kg) {
  const double E_J = Relativistic::Momentum2Energy(p,m0_kg);
  return E_J / (1.0e6 * ElectronCharge);
}

static inline double MomentumFromRigidity_GV(double R_GV,double q_C_abs) {
  return (R_GV*1.0e9*q_C_abs)/SpeedOfLight;
}

static inline double RigidityFromMomentum_GV(double p,double q_C_abs) {
  return (q_C_abs>0.0) ? (p*SpeedOfLight/q_C_abs/1.0e9) : 0.0;
}

class cFieldEvaluator : public IGridlessFieldEvaluator {
public:
  explicit cFieldEvaluator(const EarthUtil::AmpsParam& p) : prm(p) {
    // For Tsyganenko models we need Geopack initialization.
    //
    // IMPORTANT NOTE (coordinate transforms vs magnetic field evaluation):
    //   - In this gridless cutoff tool Geopack is used **ONLY** to evaluate the
    //     internal IGRF field via Geopack::IGRF::GetMagneticField().
    //   - We DO NOT use Geopack for any coordinate transformations in this file.
    //     All coordinate / direction transforms (e.g., SM<->GSM, GEO/GCE<->GSM)
    //     are intended to be handled by SPICE (pxform) when enabled.
    //
    // Historical comment (kept): older revisions mentioned "coordinate setup" here
    // because Geopack commonly provides GEO<->GSM helpers in other projects. That
    // is NOT the design in this solver; keep Geopack confined to IGRF.
    // For the analytic DIPOLE model we do NOT call Geopack at all; this avoids
    // link-time dependencies and keeps the dipole test self-contained.
    if (Model()!="DIPOLE") {
      Geopack::Init(prm.field.epoch.c_str(),"GSM");
    }

    // Configure analytic dipole parameters (used only when FIELD_MODEL=DIPOLE).
    Earth::GridlessMode::Dipole::SetMomentScale(prm.field.dipoleMoment_Me);
    Earth::GridlessMode::Dipole::SetTiltDeg(prm.field.dipoleTilt_deg);

    for (int i=0;i<11;i++) PARMOD[i]=0.0;
    PS = 0.170481; // same default as interfaces

    PARMOD[0]=prm.field.pdyn_nPa;
    PARMOD[1]=prm.field.dst_nT;
    PARMOD[2]=prm.field.imfBy_nT;
    PARMOD[3]=prm.field.imfBz_nT;

    if (Model()=="T05") {
      for (int i=0;i<6;i++) PARMOD[4+i]=prm.field.w[i];
    }
  }

  std::string Model() const { return prm.field.model; }

  void GetB_T(const V3& x_m, V3& B_T) const override {
    // Dipole branch: internal field only, analytic, no IGRF/external field.
    if (Model()=="DIPOLE") {
      double x_arr[3]={x_m.x,x_m.y,x_m.z};
      double b_arr[3];
      Earth::GridlessMode::Dipole::GetB_Tesla(x_arr,b_arr);
      B_T.x=b_arr[0]; B_T.y=b_arr[1]; B_T.z=b_arr[2];
      return;
    }

    double x_arr[3]={x_m.x,x_m.y,x_m.z};
    double b_int[3];
    Geopack::IGRF::GetMagneticField(b_int,x_arr);

    double xRe[3]={x_m.x/_EARTH__RADIUS_, x_m.y/_EARTH__RADIUS_, x_m.z/_EARTH__RADIUS_};
    double b_ext_nT[3]={0,0,0};
    int IOPT=0;

    if (Model()=="T96") {
      t96_01_(&IOPT,const_cast<double*>(PARMOD),const_cast<double*>(&PS),xRe+0,xRe+1,xRe+2,b_ext_nT+0,b_ext_nT+1,b_ext_nT+2);
    }
    else if (Model()=="T05") {
      t04_s_(&IOPT,const_cast<double*>(PARMOD),const_cast<double*>(&PS),xRe+0,xRe+1,xRe+2,b_ext_nT+0,b_ext_nT+1,b_ext_nT+2);
    }
    else {
      throw std::runtime_error("Unsupported FIELD_MODEL in gridless solver: "+Model()+" (supported: T96,T05,DIPOLE)");
    }

    B_T.x = b_int[0] + b_ext_nT[0]*_NANO_;
    B_T.y = b_int[1] + b_ext_nT[1]*_NANO_;
    B_T.z = b_int[2] + b_ext_nT[2]*_NANO_;
  }

private:
  const EarthUtil::AmpsParam& prm;
  double PARMOD[11];
  double PS;
};


//--------------------------------------------------------------------------------------
// Particle movers (integrators)
//--------------------------------------------------------------------------------------
// Historically this solver used a single mover: a relativistic Boris pusher evaluated
// with B(x_n) (field sampled at the start of the step). That is robust and widely used,
// but for strongly spatially varying B (dipole / Tsyganenko) its accuracy is limited by
// the fact that the field is assumed constant over the whole step.
//
// RECENT UPDATE (requested):
//   - Keep classic BORIS as an option.
//   - Add a more accurate option without a large cost increase.
//   - Mover selection will be wired to the input parser later; for now it is controlled
//     by a single file-scope variable (gMoverType) and a small dispatcher.
//
// The new option implemented here is BORIS_MIDPOINT:
//   - Still uses the Boris rotation structure (good long-term stability).
//   - Samples the magnetic field at a *midpoint position* x_{n+1/2}:
//         x_{n+1/2} = x_n + (dt/2) * v(x_n, p_n)
//     and uses B(x_{n+1/2}) for the rotation.
//   - This reduces error when B varies over the step (common in cutoff problems).
//
// IMPORTANT:
//   - We do NOT remove the original Boris implementation. It remains available.
//   - Until the parser is updated, gMoverType defaults to BORIS (no behavior change).
//--------------------------------------------------------------------------------------

/*
//--------------------------------------------------------------------------------------
// Legacy in-file mover selection + implementations (kept for reference; DO NOT DELETE)
//--------------------------------------------------------------------------------------
// NOTE:
//   These movers (MoverType enum, gMoverType, BorisStep/BorisMidpointStep/StepParticle)
//   were originally implemented directly in this file.
//   They have been superseded by the shared module GridlessParticleMovers.{h,cpp}.
//
//   We keep the original code here (disabled) for two reasons:
//     (1) Preserve historical comments and implementation details (per user request).
//     (2) Allow future developers to compare behavior if a regression is suspected.
//
//   ACTIVE CODE PATH:
//     - See the small "Mover selection (active)" section below, which selects a mover
//       and calls ::::StepParticle(...) from GridlessParticleMovers.
enum class MoverType {
  BORIS = 0,          // Classic Boris with B(x_n)
  BORIS_MIDPOINT = 1  // Boris rotation using B(x_{n+1/2})
};

// TODO(parser): expose this via a new input keyword, e.g. CUTOFF_MOVER BORIS|BORIS_MIDPOINT.
// For now we keep it as a compile-/edit-time selection.
static MoverType gMoverType = MoverType::BORIS;

static inline void BorisStep(V3& x, V3& p, double q_C, double m0_kg, double dt,
                             const cFieldEvaluator& field) {
  const double p2 = dot(p,p);
  const double mc = m0_kg*SpeedOfLight;
  const double gamma = std::sqrt(1.0 + p2/(mc*mc));

  V3 B; field.GetB_T(x,B);
  V3 t = mul((q_C*dt)/(2.0*gamma*m0_kg), B);
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


//--------------------------------------------------------------------------------------
// More accurate mover: Boris with midpoint-B sampling
//--------------------------------------------------------------------------------------
// This is a "minimal disruption" accuracy upgrade:
//   1) compute v_n from p_n and gamma_n
//   2) predict midpoint position x_{n+1/2} = x_n + (dt/2) v_n
//   3) evaluate B at x_{n+1/2}
//   4) perform the standard Boris rotation using that B
//   5) advance x using the updated momentum p_{n+1} (as in classic Boris)
//
// Notes:
//   - For magnetic-field-only motion, |p| (and hence gamma) should remain constant;
//     numerical drift is reduced because B is sampled at a position more representative
//     of the step.
//   - Cost: one extra field evaluation per step.
//--------------------------------------------------------------------------------------
static inline void BorisMidpointStep(V3& x, V3& p, double q_C, double m0_kg, double dt,
                                     const cFieldEvaluator& field) {
  const double p2 = dot(p,p);
  const double mc = m0_kg*SpeedOfLight;
  const double gamma = std::sqrt(1.0 + p2/(mc*mc));

  // v_n from p_n
  const V3 v_n = mul(1.0/(gamma*m0_kg), p);

  // Midpoint position predictor
  const V3 x_half = add(x, mul(0.5*dt, v_n));

  // Sample B at the midpoint position
  V3 B; field.GetB_T(x_half,B);

  // Standard Boris rotation using B(x_{n+1/2})
  V3 t = mul((q_C*dt)/(2.0*gamma*m0_kg), B);
  const double t2 = dot(t,t);
  V3 s = mul(2.0/(1.0+t2), t);

  V3 p_prime = add(p, cross(p, t));
  V3 p_plus  = add(p, cross(p_prime, s));
  p = p_plus;

  // Position update with updated momentum (same as classic Boris)
  const double p2n = dot(p,p);
  const double gamman = std::sqrt(1.0 + p2n/(mc*mc));
  V3 vnew = mul(1.0/(gamman*m0_kg), p);
  x = add(x, mul(dt, vnew));
}

//--------------------------------------------------------------------------------------
// Unified mover dispatcher (single hook point for future parser integration)
//--------------------------------------------------------------------------------------
// Keep all stepping logic behind this function so that:
//   - TraceAllowed() does not need to know which mover is active.
//   - Adding future movers (e.g., Vay, RK4) is localized.
//--------------------------------------------------------------------------------------
static inline void ::StepParticle(V3& x, V3& p, double q_C, double m0_kg, double dt,
                                const cFieldEvaluator& field) {
  switch (gMoverType) {
    case MoverType::BORIS:
      BorisStep(x,p,q_C,m0_kg,dt,field);
      break;
    case MoverType::BORIS_MIDPOINT:
      BorisMidpointStep(x,p,q_C,m0_kg,dt,field);
      break;
    default:
      // Defensive fallback: if gMoverType ever holds an unknown value,
      // use the classic Boris mover.
      BorisStep(x,p,q_C,m0_kg,dt,field);
      break;
  }
}
*/



//--------------------------------------------------------------------------------------
// Domain geometry checks
//--------------------------------------------------------------------------------------
// INPUT UNITS POLICY
//   For the RoR-style AMPS_PARAM inputs used by this gridless mode we interpret all
//   geometry distances as **kilometers** (km): DOMAIN_* bounds, R_INNER, and POINT
//   coordinates. The Tsyganenko Fortran models require coordinates in Re, so we
//   pre-convert the domain bounds to Re for fast in-loop checks.

struct DomainBoxRe {
  double xMin, xMax, yMin, yMax, zMin, zMax, rInner;
};

static inline DomainBoxRe ToDomainBoxRe(const EarthUtil::DomainBox& bKm) {
  const double km2Re = 1000.0/_EARTH__RADIUS_;
  DomainBoxRe r;
  r.xMin   = bKm.xMin*km2Re;  r.xMax   = bKm.xMax*km2Re;
  r.yMin   = bKm.yMin*km2Re;  r.yMax   = bKm.yMax*km2Re;
  r.zMin   = bKm.zMin*km2Re;  r.zMax   = bKm.zMax*km2Re;
  r.rInner = bKm.rInner*km2Re;
  return r;
}

static inline bool InsideBoxRe(const V3& xRe,const DomainBoxRe& b) {
  return (xRe.x>=b.xMin && xRe.x<=b.xMax &&
          xRe.y>=b.yMin && xRe.y<=b.yMax &&
          xRe.z>=b.zMin && xRe.z<=b.zMax);
}

static inline bool LostInnerSphere(const V3& xRe,double rInnerRe) {
  return (std::sqrt(xRe.x*xRe.x + xRe.y*xRe.y + xRe.z*xRe.z) <= rInnerRe);
}

//------------------------------------------------------------------------------
// Automatic variable time-step selection for gridless cutoff tracing
//------------------------------------------------------------------------------
// DEVELOPMENT NOTE
//   The original gridless cutoff implementation used a fixed trace step
//   (dt = prm.numerics.dtTrace_s) for every trajectory, everywhere in the
//   domain. This becomes problematic when particles enter strong-field regions:
//   the local gyroperiod shrinks, but the fixed dt does not, so the orbit can be
//   under-resolved and the cutoff classification becomes sensitive to timestep.
//
//   To improve robustness without making all runs uniformly expensive, we select
//   dt automatically at each step using the *current* particle state and local B.
//   The user-provided DT_TRACE value is preserved as a hard upper bound.
//
// CONTROLLER DESIGN (cheap, deterministic, conservative)
//   We choose dt = min( dt_user_cap,
//                       dt_gyro,
//                       dt_boundary,
//                       time_remaining )
//   where:
//     dt_gyro     limits the Boris rotation angle per step (resolves gyration)
//     dt_boundary limits travel distance per step relative to nearest stopping
//                 surface (box face or inner loss sphere) to reduce overshoot
//                 near classification boundaries.
//
//   This is intentionally not a full adaptive error-estimate / reject-retry
//   integrator. The cutoff solver calls TraceAllowed() many times; a lightweight
//   controller gives a better cost/benefit tradeoff for this application.
//
// IMPLEMENTATION DETAILS
//   - B is sampled once here to estimate local gyrofrequency; BorisStep() will
//     sample B again for the actual push. This extra field evaluation is the
//     price of adaptivity and is usually worth it when fixed dt is too large.
//   - If |B| is small, the gyro constraint becomes inactive and dt is set by
//     the user cap and/or boundary-distance cap.
//   - A small floor avoids zero or denormal dt and guarantees forward progress.
static inline double SelectAdaptiveDt(const EarthUtil::AmpsParam& prm,
                                      const cFieldEvaluator& field,
                                      const V3& x,
                                      const V3& p,
                                      double q_C,
                                      double m0_kg,
                                      const DomainBoxRe& boxRe,
                                      double timeRemaining_s) {
  // DT_TRACE from input remains the *maximum* allowed step in the adaptive mode.
  double dt = prm.numerics.dtTrace_s;
  if (timeRemaining_s < dt) dt = timeRemaining_s;

  // Compute relativistic gamma and speed from momentum p = gamma m v.
  const double p2 = dot(p,p);
  const double mc = m0_kg*SpeedOfLight;
  const double gamma = std::sqrt(1.0 + p2/(mc*mc));
  const double pMag = std::sqrt(std::max(0.0,p2));
  const double vMag = (gamma>0.0 && m0_kg>0.0) ? pMag/(gamma*m0_kg) : 0.0;

  // Local field magnitude for gyrofrequency estimate.
  V3 B; field.GetB_T(x,B);
  const double Bmag = norm(B);

  // (1) Gyro-angle criterion: keep omega_c * dt below a conservative limit.
  //     This directly addresses the "dt too big" issue in strong magnetic field.
  //     0.15 rad (~8.6 deg) is a practical compromise between accuracy and cost.
  const double dphiMax = 0.15;
  if (Bmag>0.0) {
    const double omega = std::fabs(q_C)*Bmag/(gamma*m0_kg);
    if (omega>0.0) dt = std::min(dt, dphiMax/omega);
  }

  // (2) Geometry-aware travel criterion: avoid very large jumps near stopping
  //     surfaces (box faces and inner loss sphere) where one oversized step can
  //     change the loss/escape classification.
  V3 xRe{ x.x/_EARTH__RADIUS_, x.y/_EARTH__RADIUS_, x.z/_EARTH__RADIUS_ };
  const double rRe = std::sqrt(xRe.x*xRe.x + xRe.y*xRe.y + xRe.z*xRe.z);

  double dInner_m = 1.0e300;
  if (rRe > boxRe.rInner) dInner_m = (rRe - boxRe.rInner)*_EARTH__RADIUS_;

  double dBox_m = 1.0e300;
  if (InsideBoxRe(xRe,boxRe)) {
    const double dxp = (boxRe.xMax - xRe.x)*_EARTH__RADIUS_;
    const double dxm = (xRe.x - boxRe.xMin)*_EARTH__RADIUS_;
    const double dyp = (boxRe.yMax - xRe.y)*_EARTH__RADIUS_;
    const double dym = (xRe.y - boxRe.yMin)*_EARTH__RADIUS_;
    const double dzp = (boxRe.zMax - xRe.z)*_EARTH__RADIUS_;
    const double dzm = (xRe.z - boxRe.zMin)*_EARTH__RADIUS_;
    dBox_m = std::min(std::min(std::min(dxp,dxm),std::min(dyp,dym)),std::min(dzp,dzm));
  }

  // Limit per-step travel distance to a fraction of the nearest termination
  // surface distance. This is a CFL-like geometric criterion, not a formal LTE.
  const double dNear_m = std::min(dInner_m,dBox_m);
  const double fDist = 0.20; // allow <=20% of nearest-boundary distance per step
  if (vMag>0.0 && dNear_m<1.0e299) dt = std::min(dt, fDist*dNear_m/vMag);

  // Floor for numerical robustness and guaranteed forward progress.
  // We still clamp by timeRemaining_s so we never exceed the time cap.
  const double dtFloor = std::max(1.0e-12, 1.0e-9*std::max(prm.numerics.dtTrace_s,1.0));
  dt = std::max(dtFloor, dt);
  if (timeRemaining_s > 0.0) dt = std::min(dt, timeRemaining_s);

  return dt;
}

}

static bool TraceAllowed(const EarthUtil::AmpsParam& prm,
                         const cFieldEvaluator& field,
                         const V3& x0_m,
                         const V3& v0_unit,
                         double R_GV) {
  const double q = prm.species.charge_e * ElectronCharge;
  const double qabs = std::fabs(q);
  const double m0 = prm.species.mass_amu * _AMU_;

  const double pMag = MomentumFromRigidity_GV(R_GV,qabs);
  V3 p = mul(pMag, v0_unit);
  V3 x = x0_m;

  // Convert domain bounds from km (input) to Re for checks.
  const DomainBoxRe boxRe = ToDomainBoxRe(prm.domain);

  // Adaptive integration bookkeeping.
  // We keep both a physical-time cap and a hard step-count cap. The former
  // preserves the intended tracing horizon from the input file; the latter is a
  // safety valve for pathological orbits in weak-field / trapped configurations.
  double tTrace_s = 0.0;
  int nSteps = 0;

  //----------------------------------------------------------------------------------
  // Per-trajectory time limit
  //----------------------------------------------------------------------------------
  // By default, we use the global time cap from #NUMERICAL MAX_TRACE_TIME.
  //
  // New feature (requested):
  //   If the input file provides CUTOFF_MAX_TRAJ_TIME > 0 in #CUTOFF_RIGIDITY,
  //   we use that tighter limit *only for cutoff-rigidity tracing*.
  //
  // Rationale:
  //   Cutoff scans can backtrace many trajectories; a small fraction may be
  //   quasi-trapped and would otherwise run up to the global cap, dominating
  //   runtime. This knob allows you to control that cost specifically for
  //   cutoff calculations.
  //
  // IMPORTANT: This does not change any other AMPS workflows.
  const double maxTraceTime_s_effective =
    (prm.cutoff.maxTrajTime_s > 0.0) ? prm.cutoff.maxTrajTime_s : prm.numerics.maxTraceTime_s;

  // Main trace loop with automatic dt selection. The geometric classification
  // checks are intentionally evaluated *before* the push so that starting exactly
  // outside the box (allowed) or inside the loss sphere (forbidden) is handled
  // consistently and independently of the step size.
  while (nSteps < prm.numerics.maxSteps && tTrace_s < maxTraceTime_s_effective) {
    V3 xRe{ x.x/_EARTH__RADIUS_, x.y/_EARTH__RADIUS_, x.z/_EARTH__RADIUS_ };
    if (LostInnerSphere(xRe,boxRe.rInner)) return false;
    if (!InsideBoxRe(xRe,boxRe)) return true;

    const double timeRemaining_s = maxTraceTime_s_effective - tTrace_s;

    // Automatic dt selection based on local gyrofrequency and proximity to the
    // termination surfaces. DT_TRACE remains a hard upper bound for backward
    // compatibility with existing input files.
    const double dt = SelectAdaptiveDt(prm,field,x,p,q,m0,boxRe,timeRemaining_s);

    // One relativistic Boris push with the selected local step.
    // NOTE: This used to call BorisStep() unconditionally.
    //       We now dispatch through ::StepParticle() so that alternative movers
    //       (e.g., BorisMidpointStep) can be selected without touching the trace loop.
    ::StepParticle(gDefaultMover, x, p,q,m0,dt,field);

    tTrace_s += dt;
    ++nSteps;
  }

  // Conservative fallback: if no escape is observed before time/step caps are
  // reached, classify trajectory as not allowed in the current rigidity bracket.
  return false;
}

static std::vector<V3> BuildDirGrid(int nZenith,int nAz) {
  std::vector<V3> dirs;
  dirs.reserve(nZenith*nAz);

  for (int i=0;i<nZenith;i++) {
    double mu = -1.0 + (2.0*(i+0.5))/nZenith;
    double theta = std::acos(std::max(-1.0,std::min(1.0,mu)));
    double st = std::sin(theta);

    for (int j=0;j<nAz;j++) {
      double phi = (2.0*M_PI)*(j+0.5)/nAz;
      V3 v{ st*std::cos(phi), st*std::sin(phi), std::cos(theta) };
      dirs.push_back(unit(v));
    }
  }

  return dirs;
}

//======================================================================================
// LOCAL "GSM-LIKE" LON/LAT FRAME FOR DIRECTIONAL SKY MAPS (POINT-CENTERED)
//======================================================================================
// Requested behavior:
//   - POINTS are specified in GSM Cartesian coordinates.
//   - When building a *directional* cutoff rigidity map (Rc vs arrival direction),
//     use a "GSM-type" frame that is centered on the point.
//
// What this means in practice:
//   We define a local orthonormal basis (E,N,U) at the point, analogous to the
//   usual East/North/Up frame on a sphere, but with the "north reference" taken
//   from the *global GSM Z axis*.
//
// Definitions (all vectors in GSM):
//   U = r_hat = unit(position)
//   E = unit( Z_GSM x U )
//   N = U x E
// where Z_GSM=(0,0,1). This yields:
//   - U: local "up" (radially outward from Earth center)
//   - N: local "north" (closest to +Z_GSM while tangent to the sphere)
//   - E: local "east"  (completes a right-handed triad)
//
// IMPORTANT:
//   - This is NOT geographic (GEO). It is a *local* frame that remains tied to
//     the GSM axes.
//   - Near the GSM poles (U nearly parallel to Z_GSM), the cross product becomes
//     ill-conditioned. We include a deterministic fallback to X_GSM=(1,0,0).
//
// FUTURE REPLACEMENT HOOK:
//   If you later decide to define directional-map lon/lat differently (e.g.,
//   using a magnetic-aligned frame, GEO local frame, or some event-specific
//   reference), this is the ONLY code you need to swap.
//======================================================================================

struct LocalENUFrame {
  V3 E; // local east
  V3 N; // local north
  V3 U; // local up (radial)
};

static inline LocalENUFrame BuildLocalENU_GSM(const V3& x0_m) {
  const V3 U = unit(x0_m);

  // Primary reference axis: global GSM +Z.
  const V3 Zgsm{0.0,0.0,1.0};
  V3 E = cross(Zgsm, U);

  // If near singular (point near GSM pole), fall back to GSM +X.
  if (norm(E) < 1.0e-12) {
    const V3 Xgsm{1.0,0.0,0.0};
    E = cross(Xgsm, U);
  }

  E = unit(E);
  const V3 N = unit(cross(U, E));
  return {E,N,U};
}

// Convert local (lon,lat) on the sky to a GSM unit direction vector.
//
// Convention (documented for Tecplot maps):
//   - lon_deg = 0   points toward local +N
//   - lon_deg = +90 points toward local +E
//   - lat_deg = 0   is in the local tangent plane
//   - lat_deg = +90 points toward +U (radially outward)
//
// This is a standard spherical parameterization where lon is azimuth measured
// from +N toward +E, and lat is elevation toward +U.
static inline V3 LocalLonLatToDir_GSM(const LocalENUFrame& fr, double lon_deg, double lat_deg) {
  const double lon = lon_deg * M_PI/180.0;
  const double lat = lat_deg * M_PI/180.0;

  const double cl = std::cos(lat);
  const double sl = std::sin(lat);

  // Local components in (E,N,U) basis (see convention above).
  const double dE = cl * std::sin(lon);
  const double dN = cl * std::cos(lon);
  const double dU = sl;

  // Convert to GSM Cartesian.
  const V3 d = add(add(mul(dE, fr.E), mul(dN, fr.N)), mul(dU, fr.U));
  return unit(d);
}

static double CutoffAtPoint_GV(const EarthUtil::AmpsParam& prm,
                               const cFieldEvaluator& field,
                               const V3& x0_m,
                               const std::vector<V3>& dirs,
                               double Rmin_GV,
                               double Rmax_GV,
                               int maxIter=30) {
  double Rc=-1.0;

  for (const auto& d : dirs) {
    // Backtrace convention
    V3 v0 = mul(-1.0, d);

    bool alo = TraceAllowed(prm,field,x0_m,v0,Rmin_GV);
    bool ahi = TraceAllowed(prm,field,x0_m,v0,Rmax_GV);

    if (alo && ahi) {
      Rc = (Rc<0.0) ? Rmin_GV : std::min(Rc,Rmin_GV);
      continue;
    }
    if (!ahi) continue;

    double lo=Rmin_GV, hi=Rmax_GV;
    for (int it=0;it<maxIter;it++) {
      double mid=0.5*(lo+hi);
      bool a = TraceAllowed(prm,field,x0_m,v0,mid);
      if (a) hi=mid; else lo=mid;
    }

    Rc = (Rc<0.0) ? hi : std::min(Rc,hi);
  }

  return Rc;
}

static void WriteTecplotPoints(const std::vector<EarthUtil::Vec3>& points,
                               const std::vector<double>& Rc,
                               const std::vector<double>& Emin) {
  FILE* f=std::fopen("cutoff_gridless_points.dat","w");
  if (!f) throw std::runtime_error("Cannot write Tecplot file: cutoff_gridless_points.dat");

  std::fprintf(f,"TITLE=\"Cutoff Rigidity (Gridless)\"\n");

  //----------------------------------------------------------------------------------
  // POINTS MODE OUTPUT: adding Lon/Lat labels
  //----------------------------------------------------------------------------------
  // Requested change:
  //   When OUTPUT_MODE=POINTS, also provide a Lon/Lat map of the cutoff rigidity.
  //
  // IMPORTANT CLARIFICATION (because POINTS are in GSM):
  //   - points[i].(x,y,z) are interpreted as Cartesian coordinates in the *GSM* frame.
  //   - Therefore the Lon/Lat we compute here are *GSM spherical* Lon/Lat labels:
  //       lon_GSM = atan2(y,x)
  //       lat_GSM = atan2(z, sqrt(x^2+y^2))
  //   - These are NOT geographic lon/lat (GEO). Converting to geographic requires
  //     a time-dependent GEO/GCE<->GSM transform (recommended via SPICE pxform) and is intentionally
  //     NOT done here to keep the output consistent with your GSM input.
  //
  // FUTURE REPLACEMENT HOOK:
  //   If later you decide you want geographic lon/lat, replace ONLY the two lines
  //   that compute lon_deg/lat_deg below with your preferred GEO mapping.
  //----------------------------------------------------------------------------------
  std::fprintf(f,"VARIABLES=\"id\",\"x\",\"y\",\"z\",\"lon_deg\",\"lat_deg\",\"Rc_GV\",\"Emin_MeV\"\n");
  std::fprintf(f,"ZONE T=\"points\" I=%zu F=POINT\n", points.size());

  for (size_t i=0;i<points.size();i++) {
    const double x = points[i].x;
    const double y = points[i].y;
    const double z = points[i].z;

    // GSM spherical labels (degrees). See note above.
    const double lon_deg = std::atan2(y,x) * 180.0/M_PI;
    const double lat_deg = std::atan2(z, std::sqrt(x*x + y*y)) * 180.0/M_PI;

    std::fprintf(f,"%zu %e %e %e %e %e %e %e\n",
                 i, x,y,z, lon_deg, lat_deg, Rc[i],Emin[i]);
  }

  std::fclose(f);
}



//======================================================================================
// DIPOLE ANALYTIC REFERENCE (POINTS)
//======================================================================================
// For a centered dipole, Størmer theory gives an analytic expression for the
// *vertical* cutoff rigidity. A widely used approximation (Earth-normalized) is:
//
//   Rv(λ,r) ≈ 14.9 * cos^4(λ) / (r/Re)^2   [GV]
//
// where λ is magnetic latitude and r is geocentric radius.
//
// In our verification mode we compute λ with respect to the (optionally tilted)
// dipole axis m_hat used by DipoleInterface:
//   sin(λ) = m_hat · r_hat
//
// We also scale the constant linearly with DIPOLE_MOMENT (multiples of M_E).
//
// This provides a simple analytic benchmark to compare against the numerically
// evaluated cutoff (which may include directional averaging/search).
//======================================================================================
static void WriteTecplotPoints_DipoleAnalyticCompare(const EarthUtil::AmpsParam& prm,
                                                     const std::vector<EarthUtil::Vec3>& points,
                                                     const std::vector<double>& Rc_num_GV) {
  FILE* f=std::fopen("cutoff_gridless_points_dipole_compare.dat","w");
  if (!f) throw std::runtime_error("Cannot write Tecplot file: cutoff_gridless_points_dipole_compare.dat");

  std::fprintf(f,"TITLE=\"Dipole Cutoff Rigidity: Numeric vs Analytic Vertical\"\n");
  std::fprintf(f,"VARIABLES=\"id\",\"x\",\"y\",\"z\",\"Rc_num_GV\",\"Rc_vert_GV\",\"rel_err\"\n");
  std::fprintf(f,"ZONE T=\"points\" I=%zu F=POINT\n", points.size());
// Unit dipole axis used by DipoleInterface.
  const double mx = Earth::GridlessMode::Dipole::gParams.m_hat[0];
  const double my = Earth::GridlessMode::Dipole::gParams.m_hat[1];
  const double mz = Earth::GridlessMode::Dipole::gParams.m_hat[2];

  const double R0 = StormerVerticalCoeff_GV(prm.field.dipoleMoment_Me, _EARTH__RADIUS_);

  for (size_t i=0;i<points.size();i++) {
    const double x_m = points[i].x*1000.0;
    const double y_m = points[i].y*1000.0;
    const double z_m = points[i].z*1000.0;
    const double r_m = std::sqrt(x_m*x_m + y_m*y_m + z_m*z_m);
    const double rhatx = x_m/r_m;
    const double rhaty = y_m/r_m;
    const double rhatz = z_m/r_m;

    const double sinLam = mx*rhatx + my*rhaty + mz*rhatz;
    const double cosLam = std::sqrt(std::max(0.0, 1.0 - sinLam*sinLam));
    const double rRe = r_m/_EARTH__RADIUS_;

    const double Rc_vert = R0 * prm.field.dipoleMoment_Me * std::pow(cosLam,4) / (rRe*rRe);
    const double Rc_num = Rc_num_GV[i];

    const double rel = (Rc_vert>0.0 && Rc_num>0.0) ? (Rc_num-Rc_vert)/Rc_vert : 0.0;

        std::fprintf(f,"%zu %e %e %e %e %e %e \n", i, points[i].x,points[i].y,points[i].z,Rc_num, Rc_vert, rel); 
  }

  std::fclose(f);
}

static void WriteTecplotShells(const std::vector<double>& shellAlt_km,
                               double res_deg,
                               const std::vector< std::vector<double> >& RcShell,
                               const std::vector< std::vector<double> >& EminShell) {
  // One file with multiple Tecplot zones, one per altitude.
  FILE* f=std::fopen("cutoff_gridless_shells.dat","w");
  if (!f) throw std::runtime_error("Cannot write Tecplot file: cutoff_gridless_shells.dat");

  std::fprintf(f,"TITLE=\"Cutoff Rigidity (Gridless Shells)\"\n");
  std::fprintf(f,"VARIABLES=\"lon_deg\",\"lat_deg\",\"Rc_GV\",\"Emin_MeV\"\n");

  const int nLon = static_cast<int>(std::floor(360.0/res_deg + 0.5));
  const int nLat = static_cast<int>(std::floor(180.0/res_deg + 0.5)) + 1; // include poles

  for (size_t s=0;s<shellAlt_km.size();s++) {
    const double alt=shellAlt_km[s];
    std::fprintf(f,"ZONE T=\"alt_km=%g\" I=%d J=%d F=POINT\n", alt, nLon, nLat);

    // k = i + nLon*j ordering
    for (int j=0;j<nLat;j++) {
      double lat=-90.0 + res_deg*j;
      if (lat>90.0) lat=90.0;

      for (int i=0;i<nLon;i++) {
        double lon = res_deg*i;
        int k=i+nLon*j;

        std::fprintf(f,"%e %e %e %e\n", lon, lat, 
          RcShell[s][k], EminShell[s][k]);
      }
    }
  }

  std::fclose(f);
}


//======================================================================================
// Dipole analytic vs numeric cutoff comparison (SHELLS mode)
//
// Requested update:
//   Provide the same style of benchmark output as cutoff_gridless_points_dipole_compare.dat,
//   but for SHELLS mode, formatted as a Tecplot multi-zone file analogous to
//   cutoff_gridless_shells.dat.
//
// File produced (only in nightly test dipole case):
//   cutoff_gridless_shells_dipole_compare.dat
//
// VARIABLES (per node):
//   lon_deg, lat_deg : shell grid coordinates (degrees) matching cutoff_gridless_shells.dat
//   x_km,y_km,z_km    : Cartesian location of the shell node in the solver coordinate system
//                      (currently interpreted as GSM; for DIPOLE nightly tests this matches
//                      the analytic dipole frame used by DipoleInterface).
//   Rc_num_GV         : numerically computed cutoff rigidity (what the solver produced)
//   Rc_vert_GV        : analytic vertical Størmer cutoff approximation for a dipole
//   rel_err           : (Rc_num - Rc_vert)/Rc_vert
//
// IMPORTANT NOTE ABOUT "VERTICAL":
//   The analytic formula is for the *vertical* cutoff in an ideal dipole.
//   If the numerical solver is configured for ISOTROPIC sampling, then Rc_num_GV will
//   generally be <= the vertical cutoff, and discrepancies are expected. For the
//   cleanest regression comparison, use CUTOFF_SAMPLING VERTICAL in the nightly dipole case.
//======================================================================================
static void WriteTecplotShells_DipoleAnalyticCompare(const EarthUtil::AmpsParam& prm,
                                                     const std::vector<double>& shellAlt_km,
                                                     double res_deg,
                                                     const std::vector< std::vector<double> >& RcShell) {
  FILE* f=std::fopen("cutoff_gridless_shells_dipole_compare.dat","w");
  if (!f) throw std::runtime_error("Cannot write Tecplot file: cutoff_gridless_shells_dipole_compare.dat");

  std::fprintf(f,"TITLE=\"Dipole Cutoff Rigidity (Shells): Numeric vs Analytic Vertical\"\n");
  std::fprintf(f,"VARIABLES=\"lon_deg\",\"lat_deg\",\"x_km\",\"y_km\",\"z_km\",\"Rc_num_GV\",\"Rc_vert_GV\",\"rel_err\"\n");

  const int nLon = static_cast<int>(std::floor(360.0/res_deg + 0.5));
  const int nLat = static_cast<int>(std::floor(180.0/res_deg + 0.5)) + 1; // include poles

  // Unit dipole axis used by DipoleInterface.
  const double mx = Earth::GridlessMode::Dipole::gParams.m_hat[0];
  const double my = Earth::GridlessMode::Dipole::gParams.m_hat[1];
  const double mz = Earth::GridlessMode::Dipole::gParams.m_hat[2];

  for (size_t s=0;s<shellAlt_km.size();s++) {
    const double alt_km = shellAlt_km[s];
    std::fprintf(f,"ZONE T=\"alt_km=%g\" I=%d J=%d F=POINT\n", alt_km, nLon, nLat);

    // Geocentric radius of this shell in meters.
    const double r_m = _EARTH__RADIUS_ + alt_km*1000.0;

    const double R0 = StormerVerticalCoeff_GV(prm.field.dipoleMoment_Me, _EARTH__RADIUS_);

    for (int j=0;j<nLat;j++) {
      double lat_deg = -90.0 + res_deg*j;
      if (lat_deg>90.0) lat_deg=90.0;

      const double lat = lat_deg*Pi/180.0;
      const double clat = std::cos(lat);
      const double slat = std::sin(lat);

      for (int i=0;i<nLon;i++) {
        const double lon_deg = res_deg*i;
        const double lon = lon_deg*Pi/180.0;

        const double clon = std::cos(lon);
        const double slon = std::sin(lon);

        // Cartesian position (km) consistent with shell writer:
        const double x_km = (r_m*clat*clon)/1000.0;
        const double y_km = (r_m*clat*slon)/1000.0;
        const double z_km = (r_m*slat)/1000.0;

        // Compute magnetic latitude with respect to the dipole axis m_hat used by DipoleInterface.
        const double rhatx = clat*clon;
        const double rhaty = clat*slon;
        const double rhatz = slat;

        const double sinLam = mx*rhatx + my*rhaty + mz*rhatz;
        const double cosLam = std::sqrt(std::max(0.0, 1.0 - sinLam*sinLam));
        const double rRe = r_m/_EARTH__RADIUS_;

        const double Rc_vert = R0 * prm.field.dipoleMoment_Me * std::pow(cosLam,4) / (rRe*rRe);

        const int k=i+nLon*j;
        const double Rc_num = RcShell[s][(size_t)k];
        const double rel = (Rc_vert>0.0 && Rc_num>0.0) ? (Rc_num-Rc_vert)/Rc_vert : 0.0;

        std::fprintf(f,"%e %e %e %e %e %e %e %e\n",
                     lon_deg, lat_deg, x_km, y_km, z_km, Rc_num, Rc_vert, rel);
      }
    }
  }

  std::fclose(f);
}


//--------------------------------------------------------------------------------------
// Directional cutoff rigidity sky-map writer (POINTS mode)
//--------------------------------------------------------------------------------------
// This writes a Tecplot POINT zone on a (lon,lat) grid.
//
// UPDATED (per request):
//   The directional map is now parameterized in the Solar Magnetic (SM) frame
//   using **global spherical** lon/lat. The actual tracing remains in GSM:
//     dir_SM(lon,lat) -> (SPICE pxform) -> dir_GSM -> trace
//
// Historical note (kept):
//   Earlier versions used a *local GSM-like ENU* sky coordinate system centered
//   on the injection point (BuildLocalENU_GSM / LocalLonLatToDir_GSM). That legacy
//   mapping is still present in the source as a reference/fallback.
//
// Lon/Lat meaning here is explicitly *directional* (arrival direction), in SM:
//   - lon_deg: global spherical longitude in SM (atan2(dy,dx)), [0,360)
//   - lat_deg: global spherical latitude  in SM (asin(dz)),   [-90,90]
//
// This is NOT GEO lon/lat. It is a sky-direction coordinate system tied to the
// chosen *direction-labeling* frame (SM in the current implementation).
//--------------------------------------------------------------------------------------
static void WriteTecplotDirectionalMap_Point(const std::string& fileName,
                                             int pointId,
                                             const EarthUtil::Vec3& point_km,
                                             double lonRes_deg,
                                             double latRes_deg,
                                             int nLon,
                                             int nLat,
                                             const std::vector<double>& RcCell,
                                             double qabs,
                                             double m0_kg) {
  FILE* f = std::fopen(fileName.c_str(),"w");
  if (!f) throw std::runtime_error("Cannot write Tecplot file: "+fileName);

  std::fprintf(f,"TITLE=\"Directional cutoff rigidity sky-map (POINT %d)\"\n", pointId);
  std::fprintf(f,"VARIABLES=\"lon_deg\",\"lat_deg\",\"Rc_GV\",\"Emin_MeV\"\n");
  std::fprintf(f,"ZONE T=\"point=%d x_km=%g y_km=%g z_km=%g\" I=%d J=%d F=POINT\n",
               pointId, point_km.x, point_km.y, point_km.z, nLon, nLat);

  // Cell ordering: k = iLon + nLon*jLat (same as used throughout this solver).
  for (int j=0; j<nLat; j++) {
    double lat = -90.0 + latRes_deg * j;
    if (lat > 90.0) lat = 90.0;

    for (int i=0; i<nLon; i++) {
      const double lon = lonRes_deg * i;
      const int k = i + nLon*j;
      const double rc = RcCell[(size_t)k];

      double Emin = -1.0;
      if (rc > 0.0) {
        const double pCut = MomentumFromRigidity_GV(rc, qabs);
        Emin = KineticEnergyFromMomentum_MeV(pCut, m0_kg);
      }

      std::fprintf(f,"%e %e %e %e\n", lon, lat, rc, Emin);
    }
  }

  std::fclose(f);
}


namespace Earth {
namespace GridlessMode {

int RunCutoffRigidity(const EarthUtil::AmpsParam& prm) {
  //====================================================================================
  // MPI PARALLEL EXECUTION MODEL (TRAJECTORY-BASED DYNAMIC SCHEDULING)
  //
  // Why we parallelize by "trajectory" rather than by "location":
  //   - A "location" (a point or a shell grid node) requires evaluating many candidate
  //     backtraced trajectories (different launch directions) and, for each direction,
  //     a cutoff search (multiple trace runs at different rigidities).
  //   - The wall-time of an individual trajectory evaluation can vary drastically:
  //       * some trajectories quickly escape the domain,
  //       * some hit the loss sphere quickly,
  //       * some remain trapped and run close to the time limit (slow).
  //   - If we parallelize only by locations and the number of locations is small
  //     (e.g., a few points), we may have fewer tasks than cores and waste compute.
  //
  // Strategy implemented here:
  //   - Define a "task" as:  (location_id, direction_id)
  //       -> compute the cutoff rigidity for that single direction at that location
  //          (internally: endpoint checks + bisection).
  //   - Total number of tasks = N_locations * N_directions, typically >> #cores.
  //   - Use a master/worker dynamic scheduler:
  //       * Rank 0 is the master: it feeds tasks to workers as they become idle and
  //         accumulates results for each location.
  //       * Ranks 1..(size-1) are workers: each runs expensive trajectory tracing for
  //         the task assigned, and returns the result.
  //
  // Benefits:
  //   - Excellent load balance even with highly variable trajectory wall times.
  //   - Scales well even when N_locations is small, because we have many tasks.
  //
  // Important note about MPI initialization:
  //   - You requested "assume MPI is always on", so we include mpi.h unconditionally
  //     and always compile MPI code.
  //   - However, depending on how AMPS is launched, MPI may or may not have already
  //     been initialized by the time we reach here.
  //   - To be robust for local debugging, we call MPI_Initialized() and, if needed,
  //     initialize MPI here. This allows running without mpirun in many MPI stacks
  //     (single-process MPI_COMM_WORLD), while still behaving normally under mpirun.
  //====================================================================================

  //----------------------------
  // MPI runtime initialization
  //----------------------------
  int mpiInitialized = 0;
  MPI_Initialized(&mpiInitialized);

  bool mpiInitByThisModule = false;
  if (!mpiInitialized) {
    int argc_dummy = 0;
    char** argv_dummy = nullptr;
    MPI_Init(&argc_dummy, &argv_dummy);
    mpiInitByThisModule = true;
  }

  int mpiRank = 0, mpiSize = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);

  //----------------------------------------------------------------------------
  // Precompute species + rigidity bracket derived from input energy bracket
  //----------------------------------------------------------------------------
  const double qabs = std::fabs(prm.species.charge_e * ElectronCharge);
  const double m0   = prm.species.mass_amu * _AMU_;

  const double pMin = MomentumFromKineticEnergy_MeV(prm.cutoff.eMin_MeV,m0);
  const double pMax = MomentumFromKineticEnergy_MeV(prm.cutoff.eMax_MeV,m0);
  const double Rmin = RigidityFromMomentum_GV(pMin,qabs);
  const double Rmax = RigidityFromMomentum_GV(pMax,qabs);

  if (!(Rmax>Rmin) || !(Rmax>0.0)) {
    throw std::runtime_error("Invalid cutoff energy bracket in input; cannot compute rigidity range");
  }

  //----------------------------------------------------------------------------
  // Background field evaluator (T96/T05 + IGRF) is local per rank.
  //
  // NOTE:
  //   This object typically holds model parameters and performs calls into the
  //   Tsyganenko/Geopack routines. It is safer to have one instance per rank to
  //   avoid any accidental shared state between ranks (MPI ranks are separate
  //   processes, but some libraries may still have hidden global state).
  //----------------------------------------------------------------------------
  const cFieldEvaluator field(prm);

  //----------------------------------------------------------------------------
  // Direction grid (kept identical to your current serial "meaningful results")
  //----------------------------------------------------------------------------
  const int nZenith=24;
  const int nAz=48;
  const std::vector<V3> dirs = BuildDirGrid(nZenith,nAz);

  //----------------------------------------------------------------------------------
  // Cutoff sampling mode (VERTICAL vs ISOTROPIC)
  //----------------------------------------------------------------------------------
  // Requested feature:
  //   Allow computing either:
  //     - VERTICAL cutoff: single arrival direction (toward Earth) per point.
  //     - ISOTROPIC cutoff: min over a pre-defined sky sampling grid.
  //
  // Implementation detail:
  //   We keep the existing dirs grid for isotropic sampling. For vertical, we
  //   do not use dirs[]; instead, we compute the vertical direction directly
  //   from the location position vector at runtime.
  const std::string samplingMode = EarthUtil::ToUpper(prm.cutoff.sampling);
  const bool samplingVertical = (samplingMode=="VERTICAL");
  const bool samplingIsotropic = (samplingMode=="ISOTROPIC" || samplingMode.empty());

  if (!samplingVertical && !samplingIsotropic) {
    throw std::runtime_error("Unsupported CUTOFF_SAMPLING: '"+prm.cutoff.sampling+"' (use VERTICAL or ISOTROPIC)");
  }

  //----------------------------------------------------------------------------
  // Print a run summary (rank 0 only), and flush to avoid buffered stdout issues.
  //----------------------------------------------------------------------------
  if (mpiRank==0) {
    std::cout << "================ Gridless cutoff rigidity ================\n";
    std::cout << "Run ID          : " << prm.runId << "\n";
    std::cout << "Mode            : GRIDLESS\n";
    std::cout << "Field model     : " << prm.field.model << "\n";
    std::cout << "Epoch           : " << prm.field.epoch << "\n";
    std::cout << "Species         : " << prm.species.name << " (q=" << prm.species.charge_e
              << " e, m=" << prm.species.mass_amu << " amu)\n";
    std::cout << "Rigidity bracket: [" << Rmin << ", " << Rmax << "] GV\n";
    std::cout << "CUTOFF_SAMPLING : " << (samplingVertical ? "VERTICAL" : "ISOTROPIC") << "\n";
    if (!samplingVertical) {
      std::cout << "Directions grid : " << dirs.size() << " (nZenith=" << nZenith << ", nAz=" << nAz << ")\n";
    } else {
      std::cout << "Directions grid : (not used for VERTICAL)\n";
    }
    std::cout << "Directional map : " << (prm.cutoff.directionalMap ? "ON" : "OFF") << "\n";
    if (prm.cutoff.directionalMap) {
      std::cout << "  DIRMAP_LON_RES: " << prm.cutoff.dirMapLonRes_deg << " deg\n";
      std::cout << "  DIRMAP_LAT_RES: " << prm.cutoff.dirMapLatRes_deg << " deg\n";
    }
    const DomainBoxRe boxRe = ToDomainBoxRe(prm.domain);
    std::cout << "Domain box (km) : x[" << prm.domain.xMin << "," << prm.domain.xMax << "] "
              << "y[" << prm.domain.yMin << "," << prm.domain.yMax << "] "
              << "z[" << prm.domain.zMin << "," << prm.domain.zMax << "] "
              << "rInner=" << prm.domain.rInner << "\n";
    std::cout << "Domain box (Re) : x[" << boxRe.xMin << "," << boxRe.xMax << "] "
              << "y[" << boxRe.yMin << "," << boxRe.yMax << "] "
              << "z[" << boxRe.zMin << "," << boxRe.zMax << "] "
              << "rInner=" << boxRe.rInner << "\n";
    std::cout << "dtTrace max [s] : " << prm.numerics.dtTrace_s << "  (adaptive upper bound)\n";
    std::cout << "MPI ranks       : " << mpiSize
              << " (trajectory-based dynamic scheduling)\n";
    std::cout << "==========================================================\n";
    std::cout.flush();
  }

  //====================================================================================
  // Helper: compute cutoff for a SINGLE (location, direction) task.
  //
  // The logic is the same as inside CutoffAtPoint_GV(), but extracted so we can
  // distribute per-direction work across MPI ranks.
  //
  // Return value:
  //   - If direction never becomes allowed up to Rmax -> return -1 (no cutoff / forbidden).
  //   - Else return the cutoff rigidity for this direction (>= Rmin).
  //====================================================================================
  auto CutoffForDirection_GV = [&](const V3& x0_m, const V3& dir_unit, double Rmin_GV, double Rmax_GV) -> double {
    // Backtrace convention: initial velocity points opposite to the desired arrival direction.
    const V3 v0 = mul(-1.0, dir_unit);

    // Quick endpoint classification:
    const bool alo = TraceAllowed(prm,field,x0_m,v0,Rmin_GV);
    const bool ahi = TraceAllowed(prm,field,x0_m,v0,Rmax_GV);

    // If already allowed at Rmin, the cutoff for this direction is at/below Rmin.
    if (alo && ahi) return Rmin_GV;

    // If still forbidden at Rmax, there is no allowed trajectory in the bracket.
    if (!ahi) return -1.0;

    // Otherwise: bracket exists (forbidden at Rmin, allowed at Rmax). Refine by bisection.
    double lo=Rmin_GV, hi=Rmax_GV;
    const int maxIter=24; // keep identical to serial behavior for reproducibility

    for (int it=0; it<maxIter; it++) {
      const double mid=0.5*(lo+hi);
      const bool a = TraceAllowed(prm,field,x0_m,v0,mid);
      if (a) hi=mid; else lo=mid;
    }

    // "hi" is our conservative estimate of the first allowed rigidity for this direction.
    return hi;
  };

  //====================================================================================
  // Location indexing
  //
  // POINTS mode:
  //   locationId = iPoint in [0, nPoints)
  //
  // SHELLS mode:
  //   For each shell altitude s, there is a structured lon/lat grid of size nPts=nLon*nLat.
  //   We flatten the global list of locations as:
  //     locationId = s*nPts + k,  where k in [0, nPts)
  //
  // This scheme lets workers reconstruct the 3D coordinate from (locationId) using only
  // prm.output.* data (which is available on every rank).
  //====================================================================================
  const bool isPoints = (prm.output.mode=="POINTS");
  const bool isShells = (prm.output.mode=="SHELLS");

  if (!isPoints && !isShells) {
    throw std::runtime_error("Unsupported OUTPUT_MODE for gridless cutoff solver: "+prm.output.mode);
  }

  // Shell grid geometry (only used in SHELLS mode).
  const double d_deg = isShells ? prm.output.shellRes_deg : 0.0;
  // Number of shells (altitude surfaces) requested in SHELLS mode.
  // NOTE: We define this here (next to other shell geometry quantities) so that
  //       progress reporting and per-shell completion tracking can use it without
  //       relying on any later declarations.
  const int nShells = isShells ? static_cast<int>(prm.output.shellAlt_km.size()) : 0;
  const int nLon = isShells ? static_cast<int>(std::floor(360.0/d_deg + 0.5)) : 0;
  const int nLat = isShells ? static_cast<int>(std::floor(180.0/d_deg + 0.5)) + 1 : 0;
  const int nPtsShell = isShells ? (nLon*nLat) : 0;

  // Total number of locations in this run.
  const int nLoc =
    isPoints ? static_cast<int>(prm.output.points.size())
             : static_cast<int>(prm.output.shellAlt_km.size()) * nPtsShell;

  //----------------------------------------------------------------------------------
  // Number of sampling directions for the *primary* cutoff result
  //----------------------------------------------------------------------------------
  // ISOTROPIC: use the full dirs[] grid.
  // VERTICAL : use a single direction computed from the location position.
  const int nDirSampling = samplingVertical ? 1 : static_cast<int>(dirs.size());

  //----------------------------------------------------------------------------------
  // Optional directional sky-map configuration (POINTS mode only)
  //----------------------------------------------------------------------------------
  // Requested feature: MPI-parallel sky-map computation.
  //
  // Notes:
  //   - We currently enable this only for OUTPUT_MODE=POINTS because it is
  //     primarily a diagnostic product per injection point.
  //   - Extending to SHELLS is possible but would generate extremely large
  //     outputs; if you need it, we can add a dedicated output mode.
  const bool doDirMap = (prm.cutoff.directionalMap && isPoints);

  // Directional map grid dimensions. We interpret resolutions in degrees.
  // - lon: [0,360) in steps of lonRes
  // - lat: [-90,90] inclusive in steps of latRes
  const double lonRes_deg = prm.cutoff.dirMapLonRes_deg;
  const double latRes_deg = prm.cutoff.dirMapLatRes_deg;

  int nLonMap = 0;
  int nLatMap = 0;
  int nDirMapCells = 0;
  if (doDirMap) {
    if (!(lonRes_deg>0.0) || !(latRes_deg>0.0)) {
      throw std::runtime_error("DIRMAP_LON_RES and DIRMAP_LAT_RES must be > 0 when DIRECTIONAL_MAP=T");
    }

    nLonMap = static_cast<int>(std::floor(360.0/lonRes_deg + 0.5));
    nLatMap = static_cast<int>(std::floor(180.0/latRes_deg + 0.5)) + 1; // include poles
    nDirMapCells = nLonMap * nLatMap;
  }

  //----------------------------------------------------------------------------------
  // Directional map frame transform: SM -> GSM (SPICE)
  //----------------------------------------------------------------------------------
  // The directional cutoff sky-map is parameterized in the Solar Magnetic frame (SM)
  // using **global spherical** lon/lat (not the older local ENU definition).
  //
  // For each sky-map cell (lon_SM, lat_SM):
  //   1) build a unit direction vector in SM:
  //        d_SM = (cos(lat)*cos(lon), cos(lat)*sin(lon), sin(lat))
  //   2) rotate into GSM using SPICE pxform:
  //        d_GSM = R_SM->GSM(epoch) * d_SM
  //   3) run the cutoff search/tracing in GSM (Tsyganenko expects GSM).
  //
  // WHY THIS FRAME:
  //   This matches common literature plots where the Earth-shadow for a point on
  //   +X (sunward) appears around lon~180, lat~0.
  //
  // FALLBACK BEHAVIOR:
  //   If SPICE is not enabled or the SM->GSM transform is not available, we fall
  //   back to identity. In that case the map is effectively labeled in GSM but
  //   still written as lon/lat.
  //----------------------------------------------------------------------------------
  bool spice_ok_sm2gsm = false;
  const Mat3 R_sm2gsm = doDirMap ? GetSpiceRotationOrIdentity("SM","GSM",prm.field.epoch,spice_ok_sm2gsm)
                                 : Identity3();

  // Inform the user early if SPICE transforms are not active.
  // We keep this as a WARNING (not a hard error) because some development
  // environments may compile without CSPICE. In that case the map still
  // computes, but lon/lat are effectively interpreted in GSM.
  if (mpiRank==0 && doDirMap && !spice_ok_sm2gsm) {
    std::cout << "[gridless][warning] SPICE SM->GSM transform not available. "
              << "Directional maps will fall back to identity (SM treated as GSM). "
              << "Enable AMPS_USE_SPICE and furnish kernels defining frames SM and GSM.\n";
    std::cout.flush();
  }

  //====================================================================================
  // MASTER/WORKER message protocol
  //
  // We keep messages minimal to reduce MPI overhead:
  //   Task:   {int locationId; int dirId;}
  //   Result: {int locationId; double RcDir;}
  //
  // The worker recomputes x0_m from locationId and reads dir vector from dirs[dirId].
  //====================================================================================
  // Task kinds:
  //   TASK_SAMPLING: compute Rc for one sampling direction (isotropic grid or vertical)
  //   TASK_DIRMAP  : compute Rc for one directional-map cell (lon/lat sky-map)
  enum : int {
    TASK_SAMPLING = 0,
    TASK_DIRMAP   = 1
  };

  // Task message:
  //   type : TASK_SAMPLING or TASK_DIRMAP
  //   loc  : locationId (for both task types)
  //   idx  :
  //          - TASK_SAMPLING: dirId (ignored for VERTICAL)
  //          - TASK_DIRMAP  : cellId on the lon/lat grid
  struct TaskMsg { int type; int loc; int idx; };

  // Result message mirrors task identification so the master can reduce/store.
  struct ResultMsg { int type; int loc; int idx; double rc; };

  const int TAG_TASK   = 1001;
  const int TAG_RESULT = 1002;
  const int TAG_STOP   = 1003;

  //====================================================================================
  // Rank-0 progress reporting
  //
  // We print infrequently (time-throttled) to avoid turning stdout into a scalability
  // bottleneck. The "done" count here refers to completed TRAJECTORY tasks, not locations.
  //====================================================================================
  auto nowSeconds = []() -> double {
    return MPI_Wtime();
  };

  // Wall-clock reference time used for ETA estimation in the progress bar.
  // IMPORTANT: This must be defined *before* the progress lambda below; otherwise
  //            the lambda will refer to an out-of-scope identifier and fail to compile.
  const double tStart = nowSeconds();

//====================================================================================
// Progress reporting (MASTER ONLY)
//
// IMPORTANT CONTEXT (why this is more complicated than a simple "location done" bar):
//   In this trajectory-parallel MPI design we distribute work by (locationId,dirId)
//   tasks, because:
//     * Different directions/rigidities take very different walltime (escape quickly vs
//       long trapping vs near-boundary grazing trajectories).
//     * The number of search locations (points/shell nodes) can be smaller than the
//       number of MPI ranks, but the number of trajectories is usually much larger:
//           Ntasks = Nlocations * Ndirections
//
//   As a result, the natural "unit of work" that advances smoothly is TASKS, not
//   LOCATIONS. A location is only "complete" once all its direction tasks return.
//
//   To make runtime behavior transparent (and to avoid the confusion you observed),
//   the progress bar prints BOTH:
//     - completed locations:   locDone / nLoc
//     - completed tasks:       doneTasks / totalTasks
//
//   This makes it obvious that the scheduler is working even when locDone stays at 0
//   for some time (because many direction tasks must finish before the first location
//   can be reduced to a final cutoff).
//
// OUTPUT DESIGN:
//   - Only rank 0 prints, to keep stdout readable.
//   - Print at most once per second (stdout can become a scalability bottleneck).
//   - For SHELLS mode, we additionally show "SHELL i/N alt=..." for the first not-yet-
//     completed shell, so the output resembles the original serial progress format.
//====================================================================================
auto maybePrintProgress = [&](long long doneTasks, long long totalTasks,
                              int locDone,
                              const std::vector<int>& locDonePerShell) {
  if (mpiRank!=0) return; // MASTER ONLY

  static double tLast = -1.0;
  const double t = nowSeconds();
  if (tLast < 0.0) tLast = t;
  if (t - tLast < 1.0) return; // throttle
  tLast = t;

  // --- Compute fraction based on TASKS, because tasks advance smoothly ---
  const double frac = (totalTasks>0) ? (double(doneTasks)/double(totalTasks)) : 1.0;

  // --- ETA based on tasks ---
  const double dt = t - tStart;
  const double rate = (dt>0.0) ? (double(doneTasks)/dt) : 0.0;
  double eta_s = -1.0;
  if (rate>0.0 && totalTasks>doneTasks) eta_s = double(totalTasks-doneTasks)/rate;

  auto fmt_hms = [&](double s)->std::string{
    if (s < 0.0) return std::string("--:--:--");
    long long is = (long long)std::llround(s);
    long long hh = is/3600; is-=hh*3600;
    long long mm = is/60;   is-=mm*60;
    long long ss = is;
    char buf[64];
    std::snprintf(buf,sizeof(buf),"%02lld:%02lld:%02lld",hh,mm,ss);
    return std::string(buf);
  };

  const int barW = 36;
  const int filled = (int)std::floor(frac*barW + 0.5);

  // Header label: POINTS or SHELL i/N with altitude
  if (isPoints) {
    std::cout << "[POINTS] ";
  } else {
    // Identify the first incomplete shell (for user-friendly output)
    int curShell = 0;
    for (int s=0;s<nShells;s++) {
      if (locDonePerShell[(size_t)s] < nPtsShell) { curShell = s; break; }
      curShell = s; // if all done, prints last shell
    }
    const double alt_km = prm.output.shellAlt_km[(size_t)curShell];
    std::cout << "[SHELL " << (curShell+1) << "/" << nShells << " alt=" << alt_km << "km] ";
  }

  // Prefix with rank to make it unambiguous that only rank 0 is printing.
  std::cout << "[rank " << mpiRank << "] ";

  std::cout << "[";
  for (int i=0;i<barW;i++) std::cout << (i<filled ? "#" : "-");
  std::cout << "] ";

  std::cout << std::fixed;
  std::cout.precision(1);
  std::cout << (frac*100.0) << "%  ";

  // Show BOTH locations and tasks
  std::cout << "(Loc " << locDone << "/" << nLoc << ", "
            << "Task " << doneTasks << "/" << totalTasks << ")  "
            << "ETA " << fmt_hms(eta_s) << "\n";
  std::cout.flush();
};

  //====================================================================================
  // Helper: reconstruct the starting position x0_m [m] from a flattened locationId.
  //
  // IMPORTANT: This matches your current "meaningful results" conventions:
  //   - POINTS are interpreted as GSM kilometers in the input file.
  //   - SHELLS positions are constructed in GSM Cartesian from lon/lat/alt.
  //====================================================================================
  auto LocationToX0m = [&](int locationId) -> V3 {
    if (isPoints) {
      const auto& P = prm.output.points[(size_t)locationId];
      return { P.x*1000.0, P.y*1000.0, P.z*1000.0 };
    }

    // SHELLS
    const int s = locationId / nPtsShell;
    const int k = locationId - s*nPtsShell;

    const int iLon = k % nLon;
    const int jLat = k / nLon;

    double lon = d_deg * iLon;
    double lat = -90.0 + d_deg * jLat;
    if (lat > 90.0) lat = 90.0;

    const double alt_km = prm.output.shellAlt_km[(size_t)s];
    const double r_m = (_RADIUS_(_EARTH_) + alt_km*1000.0);

    const double lonRad = lon*M_PI/180.0;
    const double latRad = lat*M_PI/180.0;
    const double cl = std::cos(latRad);

    return { r_m*cl*std::cos(lonRad), r_m*cl*std::sin(lonRad), r_m*std::sin(latRad) };
  };

  //====================================================================================
  // Storage for final results (computed on master, then written to Tecplot).
  //
  // NOTE:
  //   We only need per-location minimum cutoff (min over directions). We do not store
  //   per-direction cutoffs to keep memory bounded for large shell grids.
  //====================================================================================
  std::vector<double> RcMin;
  std::vector<double> EminMin;

  // Directional sky-map storage (POINTS only). Flattened as:
  //   RcDirMap[ pointId*nDirMapCells + cellId ]
  // where cellId = iLon + nLonMap*jLat.
  //
  // We store this only on the master (rank 0). It can be large, but for POINTS
  // mode the number of points is typically modest.
  std::vector<double> RcDirMap;

  if (mpiRank==0) {
    RcMin.assign((size_t)nLoc, -1.0);
    EminMin.assign((size_t)nLoc, -1.0);

    if (doDirMap) {
      RcDirMap.assign((size_t)prm.output.points.size() * (size_t)nDirMapCells, -1.0);
    }
  }

  //----------------------------------------------------------------------------------
  // Total number of independent tasks
  //----------------------------------------------------------------------------------
  // We now have up to two task families:
  //   (A) Primary cutoff sampling tasks (always):
  //       - ISOTROPIC: nDirSampling = size(dirs)
  //       - VERTICAL : nDirSampling = 1
  //       Total = nLoc * nDirSampling
  //
  //   (B) Optional directional sky-map tasks (POINTS only, if enabled):
  //       Total = nPoints * nDirMapCells
  //
  // We schedule both families through the same dynamic master/worker scheduler.
  const long long totalSamplingTasks = (long long)nLoc * (long long)nDirSampling;
  const long long totalDirMapTasks   = (doDirMap ? (long long)prm.output.points.size() * (long long)nDirMapCells : 0LL);
  const long long totalTasks = totalSamplingTasks + totalDirMapTasks;

  // Each worker counts how many trajectory-tasks it actually computed.
  // This is used at the end to *prove* tasks were distributed (no duplicated work).
  long long myTasksProcessed = 0;

  //====================================================================================
  // WORKER LOOP
  //
  // Workers wait for tasks, compute the cutoff for one direction, and send back results.
  //====================================================================================
  auto WorkerLoop = [&]() {
    while (true) {
      MPI_Status st;
      TaskMsg task;

      MPI_Recv(&task, (int)sizeof(TaskMsg), MPI_BYTE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &st);

      if (st.MPI_TAG == TAG_STOP) {
        // Master tells us there is no more work.
        break;
      }

      // Defensive: only TAG_TASK is expected here.
      if (st.MPI_TAG != TAG_TASK) {
        continue;
      }


      // Count this task. Each task corresponds to one (locationId,dirId) trajectory.
      // If scheduling is correct, the sum of these counters over ranks 1..N-1 should
      // equal totalTasks.
      myTasksProcessed++; 

      // Reconstruct start position (needed for both task types).
      const V3 x0_m = LocationToX0m(task.loc);

      // Compute direction cutoff depending on task type.
      double rc = -1.0;

      if (task.type == TASK_SAMPLING) {
        // Primary cutoff sampling direction.
        //
        // ISOTROPIC: use the precomputed direction grid (dirs[idx]).
        // VERTICAL : ignore idx and compute the local vertical arrival direction
        //            toward Earth: d = -unit(r0).
        V3 dir;
        if (samplingVertical) {
          dir = mul(-1.0, unit(x0_m));
        } else {
          dir = dirs[(size_t)task.idx];
        }

        rc = CutoffForDirection_GV(x0_m, dir, Rmin, Rmax);
      }
      else if (task.type == TASK_DIRMAP) {
        // Directional sky-map cell.
        //
        // UPDATED (per request): parameterize the sky-map in the global SM frame,
        // then transform the direction into GSM using SPICE.
        //
        // IMPORTANT: do NOT confuse "labeling / sampling" frame with the physics
        // tracing frame. The particle tracing is always performed in GSM because
        // the external field models (T96/T05) are defined in GSM.
        //
        // Lon/Lat meaning here is explicitly *directional* (arrival direction):
        //   - lon_SM_deg = atan2(dy,dx)  mapped to [0,360)
        //   - lat_SM_deg = asin(dz)      mapped to [-90,90]
        //
        // We sample these angles on a rectangular lon/lat grid and convert
        // (lon_SM,lat_SM) -> d_SM using the standard spherical parameterization.
        // Then we rotate d_SM into GSM:
        //    d_GSM = R_SM->GSM(epoch) * d_SM
        // and run the cutoff search using d_GSM.
        //
        // FUTURE REPLACEMENT HOOK:
        //   If later you prefer a different plotting frame (e.g., GCE, GEO),
        //   change ONLY the SPICE frame names in the rotation definition above
        //   (R_sm2gsm), and keep the math below identical.
        //
        // Cell indexing convention:
        //   cellId = iLon + nLonMap*jLat
        // where:
        //   iLon in [0,nLonMap)
        //   jLat in [0,nLatMap)
        //   lon_deg = lonRes_deg * iLon
        //   lat_deg = -90 + latRes_deg * jLat  (clamped to +90 for the last row)
        const int cellId = task.idx;
        const int iLon = cellId % nLonMap;
        const int jLat = cellId / nLonMap;

        double lon_deg = lonRes_deg * iLon;
        double lat_deg = -90.0 + latRes_deg * jLat;
        if (lat_deg > 90.0) lat_deg = 90.0;

        // Build unit direction in SM from global spherical lon/lat.
        // Note: lon_deg is sampled on [0,360) and lat_deg on [-90,90].
        const double lon = lon_deg * M_PI/180.0;
        const double lat = lat_deg * M_PI/180.0;
        const double cl  = std::cos(lat);
        const V3 dir_sm { cl*std::cos(lon), cl*std::sin(lon), std::sin(lat) };

        // Transform the direction to GSM (for tracing).
        const V3 dir_gsm = unit(Apply(R_sm2gsm, dir_sm));

        // Compute cutoff for this arrival direction (in GSM).
        rc = CutoffForDirection_GV(x0_m, dir_gsm, Rmin, Rmax);
      }

      ResultMsg res{ task.type, task.loc, task.idx, rc };
      MPI_Send(&res, (int)sizeof(ResultMsg), MPI_BYTE, 0, TAG_RESULT, MPI_COMM_WORLD);
    }
  };

  //====================================================================================
  // MASTER SCHEDULER
  //
  // Rank 0 feeds tasks to workers as they become available. This is the key piece
  // that keeps the cluster busy when trajectory costs vary widely.
  //====================================================================================
  auto MasterScheduler = [&]() {
    // Edge case: running with a single rank -> just do serial execution.
    if (mpiSize==1) {
      if (mpiRank==0) std::cout << "[gridless][MPI] size==1 -> serial fallback.\n";

      for (int loc=0; loc<nLoc; loc++) {
        const V3 x0_m = LocationToX0m(loc);

        double rcMin = -1.0;
        for (int dId=0; dId<nDirSampling; dId++) {
          V3 dir;
          if (samplingVertical) {
            // Local vertical arrival direction toward Earth.
            dir = mul(-1.0, unit(x0_m));
          } else {
            dir = dirs[(size_t)dId];
          }

          const double rc = CutoffForDirection_GV(x0_m, dir, Rmin, Rmax);
          if (rc>0.0) rcMin = (rcMin<0.0) ? rc : std::min(rcMin, rc);
        }

        RcMin[(size_t)loc] = rcMin;
        if (rcMin>0.0) {
          const double pCut = MomentumFromRigidity_GV(rcMin, qabs);
          EminMin[(size_t)loc] = KineticEnergyFromMomentum_MeV(pCut, m0);
        }

        // Optional directional sky-map (POINTS only).
        if (doDirMap) {
          const int pointId = loc;
          const LocalENUFrame fr = BuildLocalENU_GSM(x0_m);
          (void)fr; // silence unused warning in case doDirMap is compiled out later

          for (int cell=0; cell<nDirMapCells; cell++) {
            const int iLon = cell % nLonMap;
            const int jLat = cell / nLonMap;
            double lon_deg = lonRes_deg * iLon;
            double lat_deg = -90.0 + latRes_deg * jLat;
            if (lat_deg > 90.0) lat_deg = 90.0;

            const V3 dirMap = LocalLonLatToDir_GSM(BuildLocalENU_GSM(x0_m), lon_deg, lat_deg);
            const double rc = CutoffForDirection_GV(x0_m, dirMap, Rmin, Rmax);
            RcDirMap[(size_t)pointId*(size_t)nDirMapCells + (size_t)cell] = rc;
          }
        }
      }

      return;
    }


// Next task to issue (linearized over both task families).
long long nextTask = 0;
long long doneTasks = 0;

// Decode a linear task index into a TaskMsg.
// This isolates the scheduling logic so we can extend task families without
// rewriting the master/worker protocol.
auto DecodeTask = [&](long long taskId) -> TaskMsg {
  if (taskId < totalSamplingTasks) {
    const int loc = (int)(taskId / nDirSampling);
    const int idx = (int)(taskId - (long long)loc*nDirSampling);
    return TaskMsg{ TASK_SAMPLING, loc, idx };
  }

  // Directional sky-map tasks start after sampling tasks.
  const long long rem = taskId - totalSamplingTasks;
  const int pointId = (int)(rem / nDirMapCells);
  const int cellId  = (int)(rem - (long long)pointId*nDirMapCells);
  return TaskMsg{ TASK_DIRMAP, pointId, cellId };
};

//----------------------------------------------------------------------------------
// Location-level completion tracking (for user-visible progress reporting)
//
// We schedule work as (locationId,dirId) tasks, but users naturally think in terms
// of "how many points/shell nodes are finished?".
//
// A location is COMPLETE only when *all* its direction tasks have returned.
// We track that with a simple countdown initialized to nDir for each location.
//
// For SHELLS, we also track completion per shell index so the progress line can
// report "SHELL i/N alt=..." similar to the legacy serial output.
//----------------------------------------------------------------------------------
// Location completion tracking:
// If directional maps are enabled, we treat a location as "complete" only when
// BOTH the primary cutoff sampling tasks AND all its map cells are finished.
const int tasksPerLocation = nDirSampling + ((doDirMap && isPoints) ? nDirMapCells : 0);
std::vector<int> locRemain((size_t)nLoc, tasksPerLocation);
int locDone = 0;
std::vector<int> locDonePerShell;
if (!isPoints) locDonePerShell.assign((size_t)nShells, 0);

    // 1) Prime the workers with one task each (or until we run out of tasks).
    const int nWorkers = mpiSize - 1;
    for (int w=1; w<=nWorkers; w++) {
      if (nextTask >= totalTasks) break;

      TaskMsg t = DecodeTask(nextTask);
      MPI_Send(&t, (int)sizeof(TaskMsg), MPI_BYTE, w, TAG_TASK, MPI_COMM_WORLD);
      nextTask++;
    }

    // 2) Receive results and keep issuing new tasks until completion.
    while (doneTasks < totalTasks) {
      MPI_Status st;
      ResultMsg res;

      MPI_Recv(&res, (int)sizeof(ResultMsg), MPI_BYTE, MPI_ANY_SOURCE, TAG_RESULT, MPI_COMM_WORLD, &st);
      doneTasks++;

      // Update master-side storage.
      if (res.type == TASK_SAMPLING) {
        // Update per-location minimum cutoff (isotropic = min over directions;
        // vertical = min over a single direction).
        if (res.rc > 0.0) {
          double& cur = RcMin[(size_t)res.loc];
          cur = (cur<0.0) ? res.rc : std::min(cur, res.rc);
        }
      }
      else if (res.type == TASK_DIRMAP) {
        // Store Rc for this directional map cell.
        if (doDirMap) {
          RcDirMap[(size_t)res.loc*(size_t)nDirMapCells + (size_t)res.idx] = res.rc;
        }
      }

      //--------------------------------------------------------------------------------
      // Mark this (location,dir) task as completed.
      //
      // A location is 'done' only after ALL directions have been processed.
      // Tracking this lets the progress bar show both Task and Location completion.
      //--------------------------------------------------------------------------------
      {
        const int loc = res.loc;
        if (loc >= 0 && loc < nLoc) {
          locRemain[(size_t)loc]--;
          if (locRemain[(size_t)loc] == 0) {
            locDone++;
            if (!isPoints) {
              const int s = loc / nPtsShell;
              if (s >= 0 && s < nShells) locDonePerShell[(size_t)s]++;
            }
          }
        }
      }

      // Print progress occasionally (throttled).
      maybePrintProgress(doneTasks, totalTasks, locDone, locDonePerShell);

      // Send a new task to the worker that just returned, or stop it if we're done.
      if (nextTask < totalTasks) {
        TaskMsg t = DecodeTask(nextTask);
        MPI_Send(&t, (int)sizeof(TaskMsg), MPI_BYTE, st.MPI_SOURCE, TAG_TASK, MPI_COMM_WORLD);
        nextTask++;
      } else {
        // No more tasks left -> stop this worker.
        TaskMsg stopMsg{ -1, -1, -1 };
        MPI_Send(&stopMsg, (int)sizeof(TaskMsg), MPI_BYTE, st.MPI_SOURCE, TAG_STOP, MPI_COMM_WORLD);
      }
    }

    // 3) Ensure any workers that never received a STOP also get one.
    //    (This can happen when totalTasks < #workers).
    for (int w=1; w<mpiSize; w++) {
      // We attempt a nonblocking "probe" for safety; if a worker is still waiting,
      // it will receive the stop message. If it already exited, the message is harmless.
      TaskMsg stopMsg{ -1, -1, -1 };
      MPI_Send(&stopMsg, (int)sizeof(TaskMsg), MPI_BYTE, w, TAG_STOP, MPI_COMM_WORLD);
    }

    // 4) Convert RcMin -> EminMin (rank 0 only) after all *sampling* tasks
    //    have been reduced.
    for (int loc=0; loc<nLoc; loc++) {
      const double rc = RcMin[(size_t)loc];
      if (rc>0.0) {
        const double pCut = MomentumFromRigidity_GV(rc,qabs);
        EminMin[(size_t)loc] = KineticEnergyFromMomentum_MeV(pCut,m0);
      }
    }
  };

  //====================================================================================
  // Run: master or worker
  //====================================================================================
  if (mpiRank==0) {
    MasterScheduler();
  } else {
    WorkerLoop();
  }

  

      //====================================================================================
      // POST-RUN DIAGNOSTIC: Verify that tasks were truly distributed (no duplicated work)
      //
      // Symptom you observed:
      //   "All processes output the progress bar" and the intervals between progress prints
      //   increase over time. The second effect is normal when some tasks become long:
      //   the master blocks in MPI_Recv waiting for a worker to finish, so it cannot print
      //   at a fixed cadence. The first effect is usually stdout interleaving/tagging by
      //   the MPI launcher, or accidental worker printing.
      //
      // The definitive check is to count how many tasks each rank actually executed.
      // Each worker increments myTasksProcessed once per received TAG_TASK message.
      //
      // We gather these counts on rank 0 and print:
      //   - per-rank counts
      //   - sum over workers
      // and we compare sum to totalTasks.
      //
      // If scheduling is correct:
      //   sum_{r=1..size-1} myTasksProcessed[r] == totalTasks
      // (rank 0 does not compute tasks in this design).
      //====================================================================================
      std::vector<long long> taskCounts;
      if (mpiRank == 0) taskCounts.assign((size_t)mpiSize, 0);

      MPI_Gather(&myTasksProcessed, 1, MPI_LONG_LONG,
                 (mpiRank==0 ? taskCounts.data() : nullptr), 1, MPI_LONG_LONG,
                 0, MPI_COMM_WORLD);

      if (mpiRank == 0) {
        long long sumWorkers = 0;
        long long minW = (mpiSize>1 ? taskCounts[1] : 0);
        long long maxW = (mpiSize>1 ? taskCounts[1] : 0);

        for (int r=1; r<mpiSize; ++r) {
          sumWorkers += taskCounts[(size_t)r];
          minW = std::min(minW, taskCounts[(size_t)r]);
          maxW = std::max(maxW, taskCounts[(size_t)r]);
        }

        std::cout << "[gridless][MPI] Task distribution check:\n";
        std::cout << "  totalTasks (expected) = " << totalTasks << "\n";
        std::cout << "  sum(worker tasks)     = " << sumWorkers << "\n";
        if (mpiSize > 1) {
          std::cout << "  per-worker min/avg/max = " << minW
                    << " / " << (double(sumWorkers)/double(mpiSize-1))
                    << " / " << maxW << "\n";
        }
        for (int r=0; r<mpiSize; ++r) {
          std::cout << "    rank " << r << ": " << taskCounts[(size_t)r] << " tasks\n";
        }

        if (sumWorkers != totalTasks) {
          std::cout << "[gridless][MPI][WARNING] sum(worker tasks) != totalTasks.\n"
                    << "  This can happen only if:\n"
                    << "   - rank 0 also computed tasks (not in this design), or\n"
                    << "   - tasks were dropped/duplicated due to a scheduler bug.\n"
                    << "  Investigate immediately if this persists.\n";
        }
        std::cout.flush();
      }

  //====================================================================================
  // Output (rank 0 only)
  //====================================================================================
  if (mpiRank==0) {
    if (isPoints) {
      // Repackage into per-point arrays for Tecplot writer.
      std::vector<double> Rc((size_t)nLoc), Emin((size_t)nLoc);
      for (int i=0;i<nLoc;i++) { Rc[(size_t)i]=RcMin[(size_t)i]; Emin[(size_t)i]=EminMin[(size_t)i]; }

      // Preserve the original per-point console summary (optional; can be large).
      for (size_t i=0;i<prm.output.points.size();i++) {
        const auto& P = prm.output.points[i];
        std::cout << "Point " << i << " (" << P.x << "," << P.y << "," << P.z << ")"
                  << " -> Rc=" << Rc[i] << " GV, Emin=" << Emin[i] << " MeV\n";
      }

      WriteTecplotPoints(prm.output.points,Rc,Emin);

      //--------------------------------------------------------------------------------
      // Optional directional sky-maps (MPI-parallelized)
      //--------------------------------------------------------------------------------
      // If DIRECTIONAL_MAP=T, we write one Tecplot file per point:
      //   cutoff_gridless_dir_map_point_XXXX.dat
      //
      // Each file is a (lon,lat) grid in the Solar Magnetic (SM) frame using
      // global spherical lon/lat. Directions are built in SM and then rotated
      // into GSM via SPICE for the actual tracing.
      //
      // Historical note (kept): older revisions used a local GSM-like ENU frame
      // centered on each point (see BuildLocalENU_GSM / LocalLonLatToDir_GSM).
      // That mapping is still present in the file as a reference/fallback.
      //
      // IMPORTANT:
      //   The maps are computed in parallel by MPI tasks and stored in RcDirMap.
      //--------------------------------------------------------------------------------
      if (doDirMap) {
        for (size_t ip=0; ip<prm.output.points.size(); ip++) {
          char fname[256];
          std::snprintf(fname, sizeof(fname), "cutoff_gridless_dir_map_point_%04zu.dat", ip);

          // Slice out this point's cell array.
          std::vector<double> RcCell((size_t)nDirMapCells, -1.0);
          const size_t base = ip*(size_t)nDirMapCells;
          for (int k=0; k<nDirMapCells; k++) RcCell[(size_t)k] = RcDirMap[base + (size_t)k];

          WriteTecplotDirectionalMap_Point(fname,
                                           (int)ip,
                                           prm.output.points[ip],
                                           lonRes_deg,
                                           latRes_deg,
                                           nLonMap,
                                           nLatMap,
                                           RcCell,
                                           qabs,
                                           m0);
        }
        std::cout << "Wrote Tecplot: cutoff_gridless_dir_map_point_####.dat (" << prm.output.points.size() << " files)\n";
      }

// If the background model is an analytic dipole, also write an analytic
// Størmer vertical-cutoff reference for quick verification.
//
// NOTE: We only emit this comparison file in nightly-test mode to keep
// production outputs unchanged and to avoid extra I/O for regular runs.
// In nightly test mode, produce an analytic-vs-numeric comparison for the DIPOLE case.
// This avoids extra I/O for regular runs.
#if _PIC_NIGHTLY_TEST_MODE_ == _PIC_MODE_ON_ 
if (EarthUtil::ToUpper(prm.field.model)=="DIPOLE") {
  WriteTecplotPoints_DipoleAnalyticCompare(prm,prm.output.points,Rc);
}
#endif
      std::cout << "Wrote Tecplot: cutoff_gridless_points.dat\n";
    } else {
      // SHELLS: RcMin/EminMin are flattened [s*nPtsShell + k].
      if (prm.output.shellAlt_km.empty()) {
        throw std::runtime_error("OUTPUT_MODE=SHELLS but SHELL_ALTS_KM list is empty");
      }

      std::vector< std::vector<double> > RcShell(prm.output.shellAlt_km.size());
      std::vector< std::vector<double> > EminShell(prm.output.shellAlt_km.size());

      for (size_t s=0; s<prm.output.shellAlt_km.size(); s++) {
        RcShell[s].assign((size_t)nPtsShell, -1.0);
        EminShell[s].assign((size_t)nPtsShell, -1.0);

        for (int k=0;k<nPtsShell;k++) {
          const int locId = (int)(s*nPtsShell + k);
          RcShell[s][(size_t)k]   = RcMin[(size_t)locId];
          EminShell[s][(size_t)k] = EminMin[(size_t)locId];
        }

        std::cout << "Shell alt=" << prm.output.shellAlt_km[s] << " km done.\n";
      }

      WriteTecplotShells(prm.output.shellAlt_km,prm.output.shellRes_deg,RcShell,EminShell);
      std::cout << "Wrote Tecplot: cutoff_gridless_shells.dat\n";

      // In nightly test mode, produce an analytic-vs-numeric comparison for the DIPOLE case
      // in SHELLS mode. This is analogous to cutoff_gridless_points_dipole_compare.dat, but
      // formatted as a multi-zone shells file (I/J grid per altitude).
#if _PIC_NIGHTLY_TEST_MODE_ == _PIC_MODE_ON_
      if (EarthUtil::ToUpper(prm.field.model)=="DIPOLE") {
        WriteTecplotShells_DipoleAnalyticCompare(prm,prm.output.shellAlt_km,prm.output.shellRes_deg,RcShell);
        std::cout << "Wrote Tecplot: cutoff_gridless_shells_dipole_compare.dat\n";
      }
#endif
    }

    if (prm.output.coords!="GSM") {
      std::cout << "[gridless] NOTE: OUTPUT_COORDS=" << prm.output.coords
                << ". This prototype interprets positions as GSM.\n";
    }

    std::cout.flush();
  }

  //----------------------------
  // MPI finalize (only if we initialized it here)
  //----------------------------
  if (mpiInitByThisModule) {
    int mpiFinalized = 0;
    MPI_Finalized(&mpiFinalized);
    if (!mpiFinalized) MPI_Finalize();
  }

  return 0;
}


}
}
