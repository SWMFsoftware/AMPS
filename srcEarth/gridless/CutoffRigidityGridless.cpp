//======================================================================================
// CutoffRigidityGridless.cpp
//======================================================================================
// GRIDLESS CUTOFF RIGIDITY SOLVER (MPI DYNAMIC SCHEDULING VERSION)
//
// OVERVIEW
//   This implementation computes geomagnetic cutoff rigidity in a "gridless"
//   mode: instead of interpolating fields from a precomputed 3D grid, it
//   evaluates the background magnetic field directly at the particle position
//   during orbit integration. The field is composed of:
//     (1) internal IGRF field via Geopack::IGRF::GetMagneticField()
//     (2) external Tsyganenko model selected at runtime from the input file
//         (currently T96 or T05)
//
// WHY THIS MODULE EXISTS
//   The legacy AMPS Earth workflow for cutoff rigidity is tightly integrated
//   with the PIC / AMR data structures and historically selected Tsyganenko
//   models through compile-time switches (_PIC_COUPLER_MODE_). For the new
//   "gridless" workflow we need:
//     - runtime selection of background model (T96/T05 from input),
//     - a lightweight executable path driven by AMPS_PARAM files,
//     - straightforward extension toward flux calculations later.
//
// NUMERICAL ALGORITHM (PER LOCATION)
//   1) Build a fixed angular sampling grid of launch / backtrace directions.
//   2) For each direction, classify trajectories as allowed/forbidden for a
//      given rigidity R by integrating the relativistic Lorentz equation
//      (magnetic field only) with a Boris pusher.
//   3) For each direction, bracket and refine the transition using bisection
//      in rigidity between [Rmin, Rmax] derived from the input energy range.
//   4) The location cutoff rigidity Rc is the minimum directional cutoff over
//      all sampled directions.
//   5) Convert Rc to minimum kinetic energy (Emin) for reporting.
//
// TRAJECTORY CLASSIFICATION (ALLOWED / FORBIDDEN)
//   - "Allowed"    : particle exits the rectangular domain box before entering
//                    the inner loss sphere (R_INNER).
//   - "Forbidden"  : particle reaches the inner sphere, or exceeds the
//                    integration step/time cap (conservative classification).
//
// UNITS AND CONVERSIONS
//   Input parser stores geometry in kilometers (km):
//     DOMAIN_* bounds, R_INNER, POINT coordinates, shell altitudes.
//   This solver converts as needed:
//     - km -> m   : particle state and integration (SI units)
//     - km -> Re  : boundary checks and Tsyganenko model coordinates
//   Tsyganenko models are called with GSM position in Re and return B_ext in nT.
//   IGRF is returned in Tesla by Geopack interface (as used in AMPS Earth code);
//   the external field is converted nT -> T before summation.
//
// MPI PARALLELIZATION STRATEGY
//   The computational cost per location is highly variable because orbit
//   lifetime depends strongly on launch direction, rigidity, and location.
//   Some trajectories escape quickly while others approach the integration cap.
//   A static partition (equal block of locations per rank) leads to severe
//   load imbalance. To address this, this file implements a dynamic master/
//   worker scheduler:
//
//     Rank 0 (master):
//       - creates tasks (one location per task),
//       - dispatches work to workers on demand,
//       - collects results and immediately reassigns the next task,
//       - writes Tecplot output after all tasks complete.
//
//     Ranks 1..N-1 (workers):
//       - receive a task,
//       - compute Rc/Emin for that location,
//       - send result back,
//       - repeat until STOP tag is received.
//
//   Task granularity is intentionally set to ONE location per task because the
//   variance in trajectory count / runtime within a location is already large;
//   finer granularity (per direction) would increase MPI traffic and bookkeeping.
//   This one-location dynamic scheme is a good balance for shell/point scans.
//
// SERIAL FALLBACK
//   If AMPS is built/run without MPI or with a single rank, the same code path
//   executes in serial with identical physics and output format.
//
// OUTPUT
//   Tecplot ASCII files are written by rank 0:
//     - cutoff_gridless_points.dat
//     - cutoff_gridless_shells.dat
//
// EXTENSIBILITY NOTES
//   - The cutoff kernel is isolated in EvaluateLocationCutoff() and can be
//     reused later for flux / transmission-function calculations.
//   - The direction grid and bisection settings are currently fixed constants;
//     these can be moved to input-file controls in a future revision.
//======================================================================================

#include "specfunc.h"
#include "CutoffRigidityGridless.h"

#include <cstdio>
#include <cmath>
#include <stdexcept>
#include <iostream>
#include <vector>

#include "constants.h"
#include "constants.PlanetaryData.h"
#include "GeopackInterface.h"

#include <mpi.h>

extern "C" {
  void t96_01_(int*,double*,double*,double*,double*,double*,double*,double*,double*);
  void t04_s_(int*,double*,double*,double*,double*,double*,double*,double*,double*);
}

namespace {

//------------------------------------------------------------------------------
// Rank-0 text progress bar (MPI dynamic scheduling)
//
// WHY THIS EXISTS
// ---------------
// Orbit integration cost varies strongly from one search point to another:
// some trajectories escape quickly, others hit the inner sphere quickly, and
// some consume almost the full trace budget before classification. Because the
// MPI scheduler is dynamic (master/worker), the most meaningful and cheapest
// progress metric is the number of COMPLETED TASKS received by rank 0.
//
// DESIGN
// ------
// - Rank 0 only (master): no extra MPI messages are needed.
// - Progress is task-based, not step-based (robust under variable trajectory cost).
// - Output is throttled to avoid flooding stdout / scheduler logs.
// - Uses MPI_Wtime() (wall clock) for ETA estimation.
//
// NOTE
// ----
// In SHELLS mode the progress bar is instantiated PER SHELL ALTITUDE, but the
// progress unit is still a single SEARCH POINT (shell grid location / task).
//------------------------------------------------------------------------------
struct MPIRank0ProgressBar {
  int total_tasks = 1;
  int done_tasks  = 0;
  int update_count = 0;  // Track number of Print() calls for periodic newline
  double t_start = 0.0;
  double t_last_print = -1.0;
  int last_percent_printed = -1;
  const char* prefix = "[Rank 0] Progress";

  static const int kBarWidth = 32;
  static const int kNewlineEveryNUpdates = 1;  // Force newline periodically

  void Start(int total, const char* pfx = "[Rank 0] Progress") {
    total_tasks = (total > 0) ? total : 1;
    done_tasks = 0;
    update_count = 0;
    prefix = pfx ? pfx : "[Rank 0] Progress";
    t_start = MPI_Wtime();
    t_last_print = -1.0;
    last_percent_printed = -1;
    Print(true);
  }

  void Advance(int delta = 1) {
    done_tasks += delta;
    if (done_tasks > total_tasks) done_tasks = total_tasks;
    if (done_tasks < 0) done_tasks = 0;
    Print(false);
  }

  void Finish() {
    done_tasks = total_tasks;
    Print(true);
    std::fprintf(stderr, "\n");
    std::fflush(stderr);
  }

  void Print(bool force) {
    const double now = MPI_Wtime();
    const double elapsed = std::max(1.0e-12, now - t_start);
    const double frac = std::min(1.0, std::max(0.0,
      static_cast<double>(done_tasks) / static_cast<double>(total_tasks)));
    const int pct = static_cast<int>(frac * 100.0);

    // Throttle updates: print immediately if forced; otherwise update when
    // integer percent changes or after ~1 second has elapsed.
    if (!force) {
      const bool too_soon = (t_last_print >= 0.0) && ((now - t_last_print) < 1.0);
      const bool same_pct = (pct == last_percent_printed);
      if (too_soon && same_pct) return;
    }

    const double rate = static_cast<double>(done_tasks) / elapsed; // tasks/sec
    const int remain = std::max(0, total_tasks - done_tasks);
    const double eta_sec = (rate > 1.0e-12) ? (static_cast<double>(remain) / rate) : 0.0;

    const int eta_h = static_cast<int>(eta_sec / 3600.0);
    const int eta_m = (static_cast<int>(eta_sec) / 60) % 60;
    const int eta_s = static_cast<int>(eta_sec) % 60;

    int fill = static_cast<int>(frac * static_cast<double>(kBarWidth) + 0.5);
    if (fill < 0) fill = 0;
    if (fill > kBarWidth) fill = kBarWidth;

    char bar[kBarWidth + 1];
    for (int i = 0; i < kBarWidth; ++i) bar[i] = (i < fill ? '#' : '-');
    bar[kBarWidth] = '\0';

    // Pad with ~20 spaces at end to clear any trailing characters from longer previous output
    std::fprintf(stderr,
      "\r%s [%s] %6.2f%%  (%d/%d tasks)  ETA %02d:%02d:%02d                    ",
      prefix, bar, 100.0 * frac, done_tasks, total_tasks, eta_h, eta_m, eta_s);
    
    // In severely buffered environments (batch schedulers, heavy redirection),
    // even explicit fflush(stderr) with _IONBF mode may not show up until
    // process termination. Force a newline every N updates to guarantee SOMETHING
    // appears, even if it's not as pretty as a single-line updating bar.
    ++update_count;
    const bool force_newline = (update_count % kNewlineEveryNUpdates == 0) || force;
    if (force_newline) {
      std::fprintf(stderr, "\n");
    }
    
    std::fflush(stderr);

    t_last_print = now;
    last_percent_printed = pct;
  }
};


// Small local 3-vector utility used to keep the implementation self-contained.
// We intentionally avoid introducing heavier AMPS vector dependencies here
// because this module is meant to be lightweight and easy to port/refactor.
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

// Field evaluator wrapper:
//   - initializes Geopack once using the requested epoch/frame,
//   - stores Tsyganenko parameters (PARMOD, PS),
//   - evaluates B = B_IGRF + B_Tsyganenko at arbitrary x(t) during tracing.
// This object is created per rank. Workers do not communicate field values;
// each rank evaluates fields locally to minimize MPI traffic.
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

// Relativistic Boris pusher (magnetic field only).
// The Boris scheme is used because it is robust and volume-preserving for
// Lorentz motion, and handles gyration much better than naïve explicit Euler.
// State variables:
//   x [m]                  : position
//   p = gamma m v [kg m/s] : relativistic momentum
// Inputs:
//   q_C, m0_kg, dt, field evaluator
// Notes:
//   - No electric field is included in this prototype (E=0).
//   - We update momentum first, then drift position with v_{n+1}.
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

// Integrate one trajectory and classify it as allowed / forbidden.
// This is the core expensive kernel called repeatedly by the cutoff search.
// Runtime variability originates here: depending on the initial direction and
// rigidity, the orbit may terminate quickly (loss/escape) or run near the cap.
// This is exactly why the MPI layer uses dynamic scheduling.
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

  const double dt = prm.numerics.dtTrace_s;
  const int nStepsTimeCap = static_cast<int>(prm.numerics.maxTraceTime_s/dt);
  const int nSteps = std::min(prm.numerics.maxSteps, std::max(1,nStepsTimeCap));

  // Convert domain bounds from km (input) to Re for checks.
  const DomainBoxRe boxRe = ToDomainBoxRe(prm.domain);

    // Main trajectory integration loop. Each iteration performs:
  //   (1) boundary/loss checks at current position,
  //   (2) one Boris push step.
  // Boundary checks are done before stepping so immediate out-of-domain /
  // inner-sphere starting conditions are classified consistently.
  for (int i=0;i<nSteps;i++) {
    V3 xRe{ x.x/_EARTH__RADIUS_, x.y/_EARTH__RADIUS_, x.z/_EARTH__RADIUS_ };
    if (LostInnerSphere(xRe,boxRe.rInner)) return false;
    if (!InsideBoxRe(xRe,boxRe)) return true;

    BorisStep(x,p,q,m0,dt,field);
  }

  return false;
}

// Build an approximately uniform angular sampling over the sphere using
// uniform bins in mu=cos(theta) and uniform azimuth bins. Sampling in mu
// (instead of theta) avoids over-weighting the poles.
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

// Compute cutoff rigidity at one location by scanning all directions and
// taking the minimum directional cutoff. For each direction:
//   - evaluate at Rmin and Rmax,
//   - skip if no allowed orbit exists in the bracket,
//   - if both are allowed then directional cutoff <= Rmin,
//   - otherwise refine the allowed/forbidden transition by bisection.
// This produces a practical directional cutoff estimate for the chosen grid.
static double CutoffAtPoint_GV(const EarthUtil::AmpsParam& prm,
                               const cFieldEvaluator& field,
                               const V3& x0_m,
                               const std::vector<V3>& dirs,
                               double Rmin_GV,
                               double Rmax_GV,
                               int maxIter=24) {
  double Rc=-1.0;

    // Loop over direction samples. This loop is intentionally local and serial
  // within a task to keep MPI task granularity coarse (one location per task).
  // Coarse tasks reduce communication overhead and still balance well because
  // cost variance across locations is high.
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
        // Bisection search for the directional cutoff within [lo,hi].
    // Invariant after each iteration:
    //   lo -> forbidden (or less allowed),
    //   hi -> allowed.
    // We stop after maxIter fixed iterations for deterministic cost.
    for (int it=0;it<maxIter;it++) {
      double mid=0.5*(lo+hi);
      bool a = TraceAllowed(prm,field,x0_m,v0,mid);
      if (a) hi=mid; else lo=mid;
    }

    Rc = (Rc<0.0) ? hi : std::min(Rc,hi);
  }

  return Rc;
}

// Tecplot writer for explicit point list mode.
// Coordinates are written exactly as provided by the parser (km in current
// gridless convention), together with computed Rc and Emin.
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

// Tecplot writer for shell mode.
// We emit one zone per altitude and keep the logical structured indexing
// (I=nLon, J=nLat) so maps are easy to plot as contours on lon/lat grids.
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

//--------------------------------------------------------------------------------------
// MPI helpers (dynamic scheduling for variable-cost trajectory tracing)
//--------------------------------------------------------------------------------------
// WHY DYNAMIC SCHEDULING?
//   Trajectory integration wall-time is highly non-uniform. Different injection
//   points, pitch angles, and rigidities can terminate quickly (escape/loss) or
//   consume the full time cap. A static partition causes severe load imbalance.
//   We therefore use a master/worker scheduler with one location per task.
//
//   This approach is well-suited for cutoff-rigidity orbit tracing because task
//   cost variance is dominated by orbit lifetime and bisection convergence, both
//   of which vary strongly with position and direction. Dynamic scheduling keeps
//   workers busy and minimizes idle walltime.
//--------------------------------------------------------------------------------------

//--------------------------------------------------------------------------------------
// MPI runtime guard
//--------------------------------------------------------------------------------------
// The gridless cutoff solver is often tested locally as `./amps -mode gridless ...`
// without `mpirun`. Some MPI implementations support this as a valid single-rank run,
// but only if MPI is initialized before any communicator calls. The earlier version of
// this file called MPI_Comm_rank/MPI_Comm_size directly, which can abort when MPI_Init
// has not yet been called by the surrounding AMPS code.
//
// To make the solver robust in all launch modes, we:
//   1) query whether MPI is already initialized,
//   2) initialize it locally only if needed,
//   3) finalize only if *this module* initialized it.
//
// This preserves compatibility with the broader AMPS MPI lifecycle while allowing a
// safe serial fallback when the executable is run without mpirun.
struct cGridlessMpiRuntime {
  int rank = 0;
  int size = 1;
  bool initializedByThisModule = false;
};

static inline cGridlessMpiRuntime InitGridlessMpiRuntime() {
  cGridlessMpiRuntime rt;
  int mpiInitialized = 0;
  MPI_Initialized(&mpiInitialized);

  if (!mpiInitialized) {
    int argc_dummy = 0;
    char** argv_dummy = nullptr;
    MPI_Init(&argc_dummy, &argv_dummy);
    rt.initializedByThisModule = true;
  }

  MPI_Comm_rank(MPI_COMM_WORLD, &rt.rank);
  MPI_Comm_size(MPI_COMM_WORLD, &rt.size);
  return rt;
}

static inline void FinalizeGridlessMpiRuntime(cGridlessMpiRuntime& rt) {
  if (!rt.initializedByThisModule) return;
  int mpiFinalized = 0;
  MPI_Finalized(&mpiFinalized);
  if (!mpiFinalized) MPI_Finalize();
  rt.initializedByThisModule = false;
}

// Thin wrapper used by both serial and MPI code paths.
// It keeps the expensive per-location computation in one place so master/worker
// and serial loops share exactly the same physics kernel and output conversion.
static inline void EvaluateLocationCutoff(const EarthUtil::AmpsParam& prm,
                                          const cFieldEvaluator& field,
                                          const std::vector<V3>& dirs,
                                          double Rmin,double Rmax,
                                          const V3& x0_m,
                                          double qabs,double m0,
                                          double& rc_out,double& Emin_out) {
  rc_out = CutoffAtPoint_GV(prm,field,x0_m,dirs,Rmin,Rmax);
  Emin_out = -1.0;
  if (rc_out>0.0) {
    const double pCut = MomentumFromRigidity_GV(rc_out,qabs);
    Emin_out = KineticEnergyFromMomentum_MeV(pCut,m0);
  }
}


enum { GRIDLESS_MPI_TAG_TASK=5001, GRIDLESS_MPI_TAG_RESULT=5002, GRIDLESS_MPI_TAG_STOP=5003 };
struct PointTaskPayload { int idx; double x_m,y_m,z_m; };
struct PointResultPayload { int idx; double Rc_GV,Emin_MeV; };
struct ShellTaskPayload { int idx; double lon_deg,lat_deg,r_m; };
struct ShellResultPayload { int idx; double Rc_GV,Emin_MeV; };

static inline V3 Sph2CartDeg(double lonDeg,double latDeg,double r_m) {
  const double lon=lonDeg*M_PI/180.0, lat=latDeg*M_PI/180.0;
  const double cl=std::cos(lat);
  return { r_m*cl*std::cos(lon), r_m*cl*std::sin(lon), r_m*std::sin(lat)};
}

// Worker loop for explicit point tasks.
// The worker blocks waiting for either:
//   - GRIDLESS_MPI_TAG_TASK : compute one point and send result
//   - GRIDLESS_MPI_TAG_STOP : terminate worker loop
// Any unexpected tag is treated as a hard error to avoid silent corruption.
static void RunPointWorkerLoopMPI(const EarthUtil::AmpsParam& prm,const cFieldEvaluator& field,
                                  const std::vector<V3>& dirs,double Rmin,double Rmax,double qabs,double m0) {
  while (true) {
    PointTaskPayload task{}; MPI_Status st;
    MPI_Recv(&task,sizeof(task),MPI_BYTE,0,MPI_ANY_TAG,MPI_COMM_WORLD,&st);
    if (st.MPI_TAG==GRIDLESS_MPI_TAG_STOP) break;
    if (st.MPI_TAG!=GRIDLESS_MPI_TAG_TASK) throw std::runtime_error("Unexpected MPI tag (POINT worker)");
    PointResultPayload out{}; out.idx=task.idx;
    EvaluateLocationCutoff(prm,field,dirs,Rmin,Rmax,V3{task.x_m,task.y_m,task.z_m},qabs,m0,out.Rc_GV,out.Emin_MeV);
    MPI_Send(&out,sizeof(out),MPI_BYTE,0,GRIDLESS_MPI_TAG_RESULT,MPI_COMM_WORLD);
  }
}

// Worker loop for shell-grid tasks (lon/lat at a given radius).
// Same protocol as point worker, but payload carries spherical coordinates.
static void RunShellWorkerLoopMPI(const EarthUtil::AmpsParam& prm,const cFieldEvaluator& field,
                                  const std::vector<V3>& dirs,double Rmin,double Rmax,double qabs,double m0) {
  while (true) {
    ShellTaskPayload task{}; MPI_Status st;
    MPI_Recv(&task,sizeof(task),MPI_BYTE,0,MPI_ANY_TAG,MPI_COMM_WORLD,&st);
    if (st.MPI_TAG==GRIDLESS_MPI_TAG_STOP) break;
    if (st.MPI_TAG!=GRIDLESS_MPI_TAG_TASK) throw std::runtime_error("Unexpected MPI tag (SHELL worker)");
    ShellResultPayload out{}; out.idx=task.idx;
    EvaluateLocationCutoff(prm,field,dirs,Rmin,Rmax,Sph2CartDeg(task.lon_deg,task.lat_deg,task.r_m),qabs,m0,out.Rc_GV,out.Emin_MeV);
    MPI_Send(&out,sizeof(out),MPI_BYTE,0,GRIDLESS_MPI_TAG_RESULT,MPI_COMM_WORLD);
  }
}

}

namespace Earth {
namespace GridlessMode {

// Top-level driver called from srcEarth/main.cpp when '-mode gridless' is
// selected. This function:
//   1) Initializes MPI rank/size info (if enabled),
//   2) Builds the rigidity bracket from input energy limits,
//   3) Initializes field evaluator and direction grid,
//   4) Executes POINTS or SHELLS workflow (serial or MPI dynamic scheduling),
//   5) Writes Tecplot output on rank 0.
int RunCutoffRigidity(const EarthUtil::AmpsParam& prm) {
  cGridlessMpiRuntime mpiRt = InitGridlessMpiRuntime();
  
  // CRITICAL: Unbuffer stderr IMMEDIATELY, before any other output.
  // Some environments (batch schedulers, redirected streams) will buffer stderr
  // aggressively. We need unbuffered mode + explicit fflush to get real-time
  // progress visibility. We also unbuffer stdout in case there's stream interaction.
  if (mpiRt.rank == 0) {
    setvbuf(stderr, nullptr, _IONBF, 0);  // Unbuffered mode for stderr
    setvbuf(stdout, nullptr, _IOLBF, 0);  // Line-buffered mode for stdout
    // Immediate feedback that unbuffering is active:
    std::fprintf(stderr, "[Rank 0] Progress tracking enabled (unbuffered stderr)\n");
    std::fflush(stderr);
  }
  
  try {
  const int mpiRank = mpiRt.rank;
  const int mpiSize = mpiRt.size;
  const double qabs = std::fabs(prm.species.charge_e * ElectronCharge);
  const double m0   = prm.species.mass_amu * _AMU_;

  // Convert user-specified kinetic-energy bracket to rigidity bracket.
  // The cutoff search itself is done in rigidity because rigidity is the
  // natural control parameter for geomagnetic transmission/cutoff work.
  const double pMin = MomentumFromKineticEnergy_MeV(prm.cutoff.eMin_MeV,m0);
  const double pMax = MomentumFromKineticEnergy_MeV(prm.cutoff.eMax_MeV,m0);
  const double Rmin = RigidityFromMomentum_GV(pMin,qabs);
  const double Rmax = RigidityFromMomentum_GV(pMax,qabs);

  if (!(Rmax>Rmin) || !(Rmax>0.0)) {
    throw std::runtime_error("Invalid cutoff energy bracket in input; cannot compute rigidity range");
  }

  // Direction grid controls are fixed in this prototype. We declare them before
  // the summary so the summary can report the intended angular resolution even
  // before field initialization and direction-vector construction occurs.
  const int nZenith=24;
  const int nAz=48;

  // Print a run summary early (before field initialization) so the user gets
  // immediate feedback even if later initialization fails (e.g., unsupported
  // field model, Geopack/Tsyganenko setup issue, bad runtime parameters).
  if (mpiRank==0) {
  std::cout << "================ Gridless cutoff rigidity ================\n";
  std::cout << "Run ID          : " << prm.runId << "\n";
  std::cout << "Mode            : GRIDLESS\n";
  std::cout << "Field model     : " << prm.field.model << "\n";
  std::cout << "Epoch           : " << prm.field.epoch << "\n";
  std::cout << "Species         : " << prm.species.name << " (q=" << prm.species.charge_e
            << " e, m=" << prm.species.mass_amu << " amu)\n";
  std::cout << "Rigidity bracket: [" << Rmin << ", " << Rmax << "] GV\n";
  std::cout << "Directions grid : " << (nZenith*nAz) << " (nZenith=" << nZenith << ", nAz=" << nAz << ")\n";
  const DomainBoxRe boxRe = ToDomainBoxRe(prm.domain);
  std::cout << "Domain box (km) : x[" << prm.domain.xMin << "," << prm.domain.xMax << "] "
            << "y[" << prm.domain.yMin << "," << prm.domain.yMax << "] "
            << "z[" << prm.domain.zMin << "," << prm.domain.zMax << "] "
            << "rInner=" << prm.domain.rInner << "\n";
  std::cout << "Domain box (Re) : x[" << boxRe.xMin << "," << boxRe.xMax << "] "
            << "y[" << boxRe.yMin << "," << boxRe.yMax << "] "
            << "z[" << boxRe.zMin << "," << boxRe.zMax << "] "
            << "rInner=" << boxRe.rInner << "\n";
  std::cout << "dtTrace [s]     : " << prm.numerics.dtTrace_s << "\n";
  std::cout << "MPI ranks        : " << mpiSize << "\n";
  std::cout << "Scheduling       : " << ((mpiSize>1) ? "dynamic master/worker (per-location tasks)" : "serial") << "\n";
  std::cout << "==========================================================\n";
  std::cout.flush();
  }

  cFieldEvaluator field(prm);

  // Direction grid (prototype constants). These values control angular
  // resolution of the directional cutoff search and therefore the accuracy/
  // cost tradeoff. They can be exposed in the input file in a later revision.
  std::vector<V3> dirs = BuildDirGrid(nZenith,nAz);

  if (prm.output.mode=="POINTS") {
    std::vector<double> Rc(prm.output.points.size(),-1.0);
    std::vector<double> Emin(prm.output.points.size(),-1.0);

    if (mpiSize>1) {
      if (mpiRank==0) {
        size_t nextTask = 0, nDone = 0;
                // Send one point task to a worker. Returns false when the global task
        // queue is exhausted. This lambda is used both for initial filling of
        // workers and for refill-after-result dynamic scheduling.
        auto sendPointTask = [&](int dest)->bool {
          if (nextTask >= prm.output.points.size()) return false;
          const auto& P = prm.output.points[nextTask];
          PointTaskPayload t{}; t.idx = static_cast<int>(nextTask);
          t.x_m = P.x*1000.0; t.y_m = P.y*1000.0; t.z_m = P.z*1000.0;
          MPI_Send(&t,sizeof(t),MPI_BYTE,dest,GRIDLESS_MPI_TAG_TASK,MPI_COMM_WORLD);
          ++nextTask; return true;
        };
        for (int r=1; r<mpiSize; ++r) {
          if (!sendPointTask(r)) MPI_Send(nullptr,0,MPI_BYTE,r,GRIDLESS_MPI_TAG_STOP,MPI_COMM_WORLD);
        }
        // Rank-0 task-level progress for POINTS mode. One task = one search point.
        MPIRank0ProgressBar pbarPoints;
        pbarPoints.Start(static_cast<int>(prm.output.points.size()), "[Rank 0] POINTS");
                // Dynamic scheduling loop:
        //   receive a completed result from any worker,
        //   store it,
        //   immediately refill that same worker with the next task (or STOP).
        while (nDone < prm.output.points.size()) {
          PointResultPayload out{}; MPI_Status st;
          MPI_Recv(&out,sizeof(out),MPI_BYTE,MPI_ANY_SOURCE,GRIDLESS_MPI_TAG_RESULT,MPI_COMM_WORLD,&st);
          if (out.idx>=0 && static_cast<size_t>(out.idx)<Rc.size()) { Rc[out.idx]=out.Rc_GV; Emin[out.idx]=out.Emin_MeV; }
          ++nDone;
          // Progress unit is one completed search point (location task).
          pbarPoints.Advance(1);
          if (!sendPointTask(st.MPI_SOURCE)) MPI_Send(nullptr,0,MPI_BYTE,st.MPI_SOURCE,GRIDLESS_MPI_TAG_STOP,MPI_COMM_WORLD);
        }
        pbarPoints.Finish();
      }
      else {
        RunPointWorkerLoopMPI(prm,field,dirs,Rmin,Rmax,qabs,m0);
      }
            // Broadcast final arrays so all ranks leave this function with the same
      // state (useful for debugging and for future extensions that may perform
      // collective post-processing after cutoff computation).
      MPI_Bcast(Rc.data(), static_cast<int>(Rc.size()), MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Bcast(Emin.data(), static_cast<int>(Emin.size()), MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }
    else
    {
      // Serial fallback: same per-location kernel, no MPI scheduling.
      // IMPORTANT: We still provide a rank-0 progress bar here (with one MPI rank,
      // rank 0 is the only rank). This gives user feedback for long serial runs and
      // avoids the "no progress visible until Ctrl-C" issue when combined with
      // explicit stdout flushing in MPIRank0ProgressBar.
      MPIRank0ProgressBar pbarPointsSerial;
      pbarPointsSerial.Start(static_cast<int>(prm.output.points.size()), "[Rank 0] POINTS");
      for (size_t i=0;i<prm.output.points.size();i++) {
        const auto& P = prm.output.points[i];
        V3 x0_m{ P.x*1000.0, P.y*1000.0, P.z*1000.0 };
        EvaluateLocationCutoff(prm,field,dirs,Rmin,Rmax,x0_m,qabs,m0,Rc[i],Emin[i]);
        // Progress unit is one completed search point (location task).
        pbarPointsSerial.Advance(1);
      }
      pbarPointsSerial.Finish();
    }

    if (mpiRank==0) {
      WriteTecplotPoints(prm.output.points,Rc,Emin);
      std::cout << "Wrote Tecplot: cutoff_gridless_points.dat\n";
      std::cout.flush();
    }
  }
  else if (prm.output.mode=="SHELLS") {
    if (prm.output.shellAlt_km.empty()) {
      throw std::runtime_error("OUTPUT_MODE=SHELLS but SHELL_ALTS_KM list is empty");
    }

    const double d = prm.output.shellRes_deg;
    const int nLon = static_cast<int>(std::floor(360.0/d + 0.5));
    const int nLat = static_cast<int>(std::floor(180.0/d + 0.5)) + 1;
    const int nPts = nLon*nLat;

    std::vector< std::vector<double> > RcShell(prm.output.shellAlt_km.size());
    std::vector< std::vector<double> > EminShell(prm.output.shellAlt_km.size());

    // Process shell altitudes one-by-one. Keeping altitude as the outer loop
    // limits memory footprint and lets us reuse the same worker code with a
    // simple task payload (lon, lat, r).
    for (size_t s=0;s<prm.output.shellAlt_km.size();s++) {
      const double alt_km=prm.output.shellAlt_km[s];
      // Radius of the shell in meters (Earth radius + altitude).
      const double r_m = (_RADIUS_(_EARTH_) + alt_km*1000.0);

      RcShell[s].assign(nPts,-1.0);
      EminShell[s].assign(nPts,-1.0);

      if (mpiSize>1) {
        if (mpiRank==0) {
          int nextK=0, nDone=0;
                    // Send one shell-grid location (lon/lat at current shell radius).
          // Task index k maps to (i,j) using the Tecplot-compatible ordering.
          auto sendShellTask = [&](int dest)->bool {
            if (nextK>=nPts) return false;
            const int k=nextK++;
            const int i=k % nLon;
            const int j=k / nLon;
            double lat=-90.0 + d*j; if (lat>90.0) lat=90.0;
            double lon=d*i;
            ShellTaskPayload t{}; t.idx=k; t.lon_deg=lon; t.lat_deg=lat; t.r_m=r_m;
            MPI_Send(&t,sizeof(t),MPI_BYTE,dest,GRIDLESS_MPI_TAG_TASK,MPI_COMM_WORLD);
            return true;
          };
          for (int r=1; r<mpiSize; ++r) {
            if (!sendShellTask(r)) MPI_Send(nullptr,0,MPI_BYTE,r,GRIDLESS_MPI_TAG_STOP,MPI_COMM_WORLD);
          }
          // Rank-0 task-level progress for this shell altitude.
          // IMPORTANT: progress granularity is ONE SEARCH POINT (shell grid cell),
          // not one altitude. This provides useful feedback even for large shells.
          char shellPfx[160];
          std::snprintf(shellPfx, sizeof(shellPfx), "[Rank 0] SHELL %zu/%zu alt=%.3f km", s+1, prm.output.shellAlt_km.size(), alt_km);
          MPIRank0ProgressBar pbarShell;
          pbarShell.Start(nPts, shellPfx);
                    // Same dynamic scheduling pattern as POINTS mode, but for shell
          // grid cells of the current altitude.
          while (nDone<nPts) {
            ShellResultPayload out{}; MPI_Status st;
            MPI_Recv(&out,sizeof(out),MPI_BYTE,MPI_ANY_SOURCE,GRIDLESS_MPI_TAG_RESULT,MPI_COMM_WORLD,&st);
            if (out.idx>=0 && out.idx<nPts) { RcShell[s][out.idx]=out.Rc_GV; EminShell[s][out.idx]=out.Emin_MeV; }
            ++nDone;
            // Progress unit is one completed shell-grid search point.
            pbarShell.Advance(1);
            if (!sendShellTask(st.MPI_SOURCE)) MPI_Send(nullptr,0,MPI_BYTE,st.MPI_SOURCE,GRIDLESS_MPI_TAG_STOP,MPI_COMM_WORLD);
          }
          pbarShell.Finish();
        }
        else {
          RunShellWorkerLoopMPI(prm,field,dirs,Rmin,Rmax,qabs,m0);
        }
        MPI_Bcast(RcShell[s].data(), nPts, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(EminShell[s].data(), nPts, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      }
      else
      {
        // Serial fallback for shell grid traversal.
        // We show progress per SEARCH POINT (shell grid cell), matching the MPI
        // dynamic-scheduling progress granularity. This gives fine-grained feedback
        // during large shell runs even when only one rank is used.
        char shellPfx[160];
        std::snprintf(shellPfx, sizeof(shellPfx), "[Rank 0] SHELL %zu/%zu alt=%.3f km", s+1, prm.output.shellAlt_km.size(), alt_km);
        MPIRank0ProgressBar pbarShellSerial;
        pbarShellSerial.Start(nPts, shellPfx);
        for (int j=0;j<nLat;j++) {
          double lat=-90.0 + d*j; if (lat>90.0) lat=90.0;
          for (int i=0;i<nLon;i++) {
            double lon=d*i; int k=i+nLon*j;
            EvaluateLocationCutoff(prm,field,dirs,Rmin,Rmax,Sph2CartDeg(lon,lat,r_m),qabs,m0,RcShell[s][k],EminShell[s][k]);
            // Progress unit is one completed shell-grid search point.
            pbarShellSerial.Advance(1);
          }
        }
        pbarShellSerial.Finish();
      }

      if (mpiRank==0) { std::cout << "Shell alt=" << alt_km << " km done.\n"; std::cout.flush(); }
    }

    if (mpiRank==0) {
      WriteTecplotShells(prm.output.shellAlt_km,prm.output.shellRes_deg,RcShell,EminShell);
      std::cout << "Wrote Tecplot: cutoff_gridless_shells.dat\n";
      std::cout.flush();
    }
  }
  else {
    throw std::runtime_error("Unsupported OUTPUT_MODE for gridless cutoff solver: "+prm.output.mode);
  }

  if (mpiRank==0 && prm.output.coords!="GSM") {
    std::cout << "[gridless] NOTE: OUTPUT_COORDS=" << prm.output.coords
              << ". This prototype interprets positions as GSM.\n";
    std::cout.flush();
  }

  if (mpiSize>1) MPI_Barrier(MPI_COMM_WORLD);
  FinalizeGridlessMpiRuntime(mpiRt);
  return 0;
  }
  catch (...) {
    FinalizeGridlessMpiRuntime(mpiRt);
    throw;
  }
}

}
}
