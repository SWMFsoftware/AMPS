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
  const double qabs = std::fabs(prm.species.charge_e * ElectronCharge);
  const double m0   = prm.species.mass_amu * _AMU_;

  const double pMin = MomentumFromKineticEnergy_MeV(prm.cutoff.eMin_MeV,m0);
  const double pMax = MomentumFromKineticEnergy_MeV(prm.cutoff.eMax_MeV,m0);
  const double Rmin = RigidityFromMomentum_GV(pMin,qabs);
  const double Rmax = RigidityFromMomentum_GV(pMax,qabs);

  if (!(Rmax>Rmin) || !(Rmax>0.0)) {
    throw std::runtime_error("Invalid cutoff energy bracket in input; cannot compute rigidity range");
  }

  cFieldEvaluator field(prm);

  // Direction grid (prototype constants)
  const int nZenith=24;
  const int nAz=48;
  std::vector<V3> dirs = BuildDirGrid(nZenith,nAz);

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
  std::cout << "dt integrator   : variable dt (gyro-angle + boundary-distance caps)\n";
  std::cout << "==========================================================\n";

  if (prm.output.mode=="POINTS") {
    std::vector<double> Rc(prm.output.points.size(),-1.0);
    std::vector<double> Emin(prm.output.points.size(),-1.0);

    for (size_t i=0;i<prm.output.points.size();i++) {
      const auto& P = prm.output.points[i];

      // IMPORTANT (UNITS): POINT coordinates are interpreted as GSM kilometers.
      V3 x0_m{ P.x*1000.0, P.y*1000.0, P.z*1000.0 };

      double rc = CutoffAtPoint_GV(prm,field,x0_m,dirs,Rmin,Rmax);
      Rc[i]=rc;
      if (rc>0.0) {
        double pCut = MomentumFromRigidity_GV(rc,qabs);
        Emin[i] = KineticEnergyFromMomentum_MeV(pCut,m0);
      }

      std::cout << "Point " << i << " (" << P.x << "," << P.y << "," << P.z << ")"
                << " -> Rc=" << Rc[i] << " GV, Emin=" << Emin[i] << " MeV\n";
    }

    WriteTecplotPoints(prm.output.points,Rc,Emin);
    std::cout << "Wrote Tecplot: cutoff_gridless_points.dat\n";
  }
  else if (prm.output.mode=="SHELLS") {
    if (prm.output.shellAlt_km.empty()) {
      throw std::runtime_error("OUTPUT_MODE=SHELLS but SHELL_ALTS_KM list is empty");
    }

    const double d = prm.output.shellRes_deg;
    const int nLon = static_cast<int>(std::floor(360.0/d + 0.5));
    const int nLat = static_cast<int>(std::floor(180.0/d + 0.5)) + 1;
    const int nPts = nLon*nLat;

    auto sph2cart = [](double lonDeg,double latDeg,double r_m)->V3{
      const double lon=lonDeg*M_PI/180.0;
      const double lat=latDeg*M_PI/180.0;
      const double cl=std::cos(lat);
      return { r_m*cl*std::cos(lon), r_m*cl*std::sin(lon), r_m*std::sin(lat)};
    };

    std::vector< std::vector<double> > RcShell(prm.output.shellAlt_km.size());
    std::vector< std::vector<double> > EminShell(prm.output.shellAlt_km.size());

    for (size_t s=0;s<prm.output.shellAlt_km.size();s++) {
      const double alt_km=prm.output.shellAlt_km[s];
      const double r_m = (_RADIUS_(_EARTH_) + alt_km*1000.0);

      RcShell[s].assign(nPts,-1.0);
      EminShell[s].assign(nPts,-1.0);

      for (int j=0;j<nLat;j++) {
        double lat=-90.0 + d*j;
        if (lat>90.0) lat=90.0;

        for (int i=0;i<nLon;i++) {
          double lon=d*i;
          int k=i+nLon*j;

          V3 x0 = sph2cart(lon,lat,r_m);
          double rc = CutoffAtPoint_GV(prm,field,x0,dirs,Rmin,Rmax);
          RcShell[s][k]=rc;
          if (rc>0.0) {
            double pCut = MomentumFromRigidity_GV(rc,qabs);
            EminShell[s][k] = KineticEnergyFromMomentum_MeV(pCut,m0);
          }
        }
      }

      std::cout << "Shell alt=" << alt_km << " km done.\n";
    }

    WriteTecplotShells(prm.output.shellAlt_km,prm.output.shellRes_deg,RcShell,EminShell);
    std::cout << "Wrote Tecplot: cutoff_gridless_shells.dat\n";
  }
  else {
    throw std::runtime_error("Unsupported OUTPUT_MODE for gridless cutoff solver: "+prm.output.mode);
  }

  if (prm.output.coords!="GSM") {
    std::cout << "[gridless] NOTE: OUTPUT_COORDS=" << prm.output.coords
              << ". This prototype interprets positions as GSM.\n";
  }

  return 0;
}

}
}
