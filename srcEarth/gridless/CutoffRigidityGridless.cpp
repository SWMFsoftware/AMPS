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

#include "specfunc.h"
#include "CutoffRigidityGridless.h"

#include <cstdio>
#include <cmath>
#include <stdexcept>
#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <numeric>
#include <mpi.h> 

#include "constants.h"
#include "constants.PlanetaryData.h"
#include "GeopackInterface.h"

extern "C" {
  void t96_01_(int*,double*,double*,double*,double*,double*,double*,double*,double*);
  void t04_s_(int*,double*,double*,double*,double*,double*,double*,double*,double*);
}

namespace {

struct V3 { double x,y,z; };
static inline V3 add(const V3&a,const V3&b){return {a.x+b.x,a.y+b.y,a.z+b.z};}
static inline V3 mul(double s,const V3&a){return {s*a.x,s*a.y,s*a.z};}
static inline double dot(const V3&a,const V3&b){return a.x*b.x+a.y*b.y+a.z*b.z;}
static inline V3 cross(const V3&a,const V3&b){return {a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x};}
static inline double norm(const V3&a){return std::sqrt(dot(a,a));}
static inline V3 unit(const V3&a){double n=norm(a); return (n>0)?mul(1.0/n,a):V3{0,0,0};}

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

class cFieldEvaluator {
public:
  explicit cFieldEvaluator(const EarthUtil::AmpsParam& p) : prm(p) {
    Geopack::Init(prm.field.epoch.c_str(),"GSM");

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

  void GetB_T(const V3& x_m, V3& B_T) const {
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
      throw std::runtime_error("Unsupported FIELD_MODEL in gridless solver: "+Model()+" (supported: T96,T05)");
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

  // Main trace loop with automatic dt selection. The geometric classification
  // checks are intentionally evaluated *before* the push so that starting exactly
  // outside the box (allowed) or inside the loss sphere (forbidden) is handled
  // consistently and independently of the step size.
  while (nSteps < prm.numerics.maxSteps && tTrace_s < prm.numerics.maxTraceTime_s) {
    V3 xRe{ x.x/_EARTH__RADIUS_, x.y/_EARTH__RADIUS_, x.z/_EARTH__RADIUS_ };
    if (LostInnerSphere(xRe,boxRe.rInner)) return false;
    if (!InsideBoxRe(xRe,boxRe)) return true;

    const double timeRemaining_s = prm.numerics.maxTraceTime_s - tTrace_s;

    // Automatic dt selection based on local gyrofrequency and proximity to the
    // termination surfaces. DT_TRACE remains a hard upper bound for backward
    // compatibility with existing input files.
    const double dt = 100*SelectAdaptiveDt(prm,field,x,p,q,m0,boxRe,timeRemaining_s);

    // One relativistic Boris push with the selected local step.
    BorisStep(x,p,q,m0,dt,field);

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

static double CutoffAtPoint_GV(const EarthUtil::AmpsParam& prm,
                               const cFieldEvaluator& field,
                               const V3& x0_m,
                               const std::vector<V3>& dirs,
                               double Rmin_GV,
                               double Rmax_GV,
                               int maxIter=24) {
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
  std::fprintf(f,"VARIABLES=\"id\",\"x\",\"y\",\"z\",\"Rc_GV\",\"Emin_MeV\"\n");
  std::fprintf(f,"ZONE T=\"points\" I=%zu F=POINT\n", points.size());

  for (size_t i=0;i<points.size();i++) {
    std::fprintf(f,"%zu %e %e %e %e %e\n", i, points[i].x,points[i].y,points[i].z, Rc[i],Emin[i]);
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
  std::fprintf(f,"VARIABLES=\"lon_deg\",\"lat_deg\",\"alt_km\",\"Rc_GV\",\"Emin_MeV\"\n");

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

        std::fprintf(f,"%e %e %e %e %e\n", lon, lat, alt,
          RcShell[s][k], EminShell[s][k]);
      }
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
    std::cout << "Directions grid : " << dirs.size() << " (nZenith=" << nZenith << ", nAz=" << nAz << ")\n";
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

  const int nDir = static_cast<int>(dirs.size());

  //====================================================================================
  // MASTER/WORKER message protocol
  //
  // We keep messages minimal to reduce MPI overhead:
  //   Task:   {int locationId; int dirId;}
  //   Result: {int locationId; double RcDir;}
  //
  // The worker recomputes x0_m from locationId and reads dir vector from dirs[dirId].
  //====================================================================================
  struct TaskMsg { int loc; int dir; };
  struct ResultMsg { int loc; double rc; };

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

  if (mpiRank==0) {
    RcMin.assign((size_t)nLoc, -1.0);
    EminMin.assign((size_t)nLoc, -1.0);
  }

  // Total number of independent tasks.
  const long long totalTasks = (long long)nLoc * (long long)nDir;

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

      // Reconstruct start position and direction.
      const V3 x0_m = LocationToX0m(task.loc);
      const V3 dir  = dirs[(size_t)task.dir];

      // Compute direction cutoff.
      const double rc = CutoffForDirection_GV(x0_m, dir, Rmin, Rmax);

      ResultMsg res{ task.loc, rc };
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
        for (int dId=0; dId<nDir; dId++) {
          const double rc = CutoffForDirection_GV(x0_m, dirs[(size_t)dId], Rmin, Rmax);
          if (rc>0.0) rcMin = (rcMin<0.0) ? rc : std::min(rcMin, rc);
        }

        RcMin[(size_t)loc] = rcMin;
        if (rcMin>0.0) {
          const double pCut = MomentumFromRigidity_GV(rcMin, qabs);
          EminMin[(size_t)loc] = KineticEnergyFromMomentum_MeV(pCut, m0);
        }
      }

      return;
    }


// Next task to issue (linearized over loc then dir).
long long nextTask = 0;
long long doneTasks = 0;

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
std::vector<int> locRemain((size_t)nLoc, nDir);
int locDone = 0;
std::vector<int> locDonePerShell;
if (!isPoints) locDonePerShell.assign((size_t)nShells, 0);

    // 1) Prime the workers with one task each (or until we run out of tasks).
    const int nWorkers = mpiSize - 1;
    for (int w=1; w<=nWorkers; w++) {
      if (nextTask >= totalTasks) break;

      const int loc = (int)(nextTask / nDir);
      const int dir = (int)(nextTask - (long long)loc*nDir);
      TaskMsg t{ loc, dir };

      MPI_Send(&t, (int)sizeof(TaskMsg), MPI_BYTE, w, TAG_TASK, MPI_COMM_WORLD);
      nextTask++;
    }

    // 2) Receive results and keep issuing new tasks until completion.
    while (doneTasks < totalTasks) {
      MPI_Status st;
      ResultMsg res;

      MPI_Recv(&res, (int)sizeof(ResultMsg), MPI_BYTE, MPI_ANY_SOURCE, TAG_RESULT, MPI_COMM_WORLD, &st);
      doneTasks++;

      // Update per-location minimum cutoff.
      if (res.rc > 0.0) {
        double& cur = RcMin[(size_t)res.loc];
        cur = (cur<0.0) ? res.rc : std::min(cur, res.rc);
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
        const int loc = (int)(nextTask / nDir);
        const int dir = (int)(nextTask - (long long)loc*nDir);
        TaskMsg t{ loc, dir };

        MPI_Send(&t, (int)sizeof(TaskMsg), MPI_BYTE, st.MPI_SOURCE, TAG_TASK, MPI_COMM_WORLD);
        nextTask++;
      } else {
        // No more tasks left -> stop this worker.
        TaskMsg stopMsg{ -1, -1 };
        MPI_Send(&stopMsg, (int)sizeof(TaskMsg), MPI_BYTE, st.MPI_SOURCE, TAG_STOP, MPI_COMM_WORLD);
      }
    }

    // 3) Ensure any workers that never received a STOP also get one.
    //    (This can happen when totalTasks < #workers).
    for (int w=1; w<mpiSize; w++) {
      // We attempt a nonblocking "probe" for safety; if a worker is still waiting,
      // it will receive the stop message. If it already exited, the message is harmless.
      TaskMsg stopMsg{ -1, -1 };
      MPI_Send(&stopMsg, (int)sizeof(TaskMsg), MPI_BYTE, w, TAG_STOP, MPI_COMM_WORLD);
    }

    // 4) Convert RcMin -> EminMin (rank 0 only) after all tasks have been reduced.
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
