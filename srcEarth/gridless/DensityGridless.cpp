//======================================================================================
// DensityGridless.cpp
//======================================================================================
// THEORY (what is computed)
// -------------------------------------------------------------------------------------
// We compute the *local energetic-particle spectrum* and the corresponding *total number
// density* at user-specified spatial locations (OUTPUT_MODE = POINTS).
//
// 1) Directional transmissivity T(E)
//    For a given kinetic energy E and a given location x0, we sample a fixed grid of
//    directions on the unit sphere. For each direction we backtrace a test particle in
//    the prescribed magnetic field (IGRF + Tsyganenko external field, evaluated on the
//    fly). A trajectory is classified as:
//      - ALLOWED   : it escapes the outer domain box before hitting the inner loss sphere
//      - FORBIDDEN : it hits the inner sphere or times out before escaping
//
//    The transmissivity is the *fraction of allowed directions*:
//        T(E; x0) = N_allowed(E; x0) / N_dirs
//
//    This is the same "allowed vs forbidden" concept used in the cutoff rigidity tool,
//    but here we evaluate it over an energy grid rather than searching for a single
//    cutoff.
//
// 2) Local spectrum folding
//    The boundary spectrum J_b(E) is defined in the #SPECTRUM section of AMPS_PARAM.in
//    and is evaluated via the global ::gSpectrum object:
//        J_b(E) = gSpectrum.GetSpectrum(E[J])
//
//    The local spectrum is modeled as
//        J_loc(E; x0) = T(E; x0) * J_b(E)
//
//    i.e., the magnetosphere acts as an energy-dependent filter with transmissivity T.
//    This is a deliberately simple (but reproducible) model that matches how cutoff
//    rigidity calculations are commonly interpreted.
//
// 3) Total number density
//    Under an isotropic assumption, differential flux J(E) [m^-2 s^-1 sr^-1 J^-1] is
//    related to differential number density n_E [m^-3 J^-1] by:
//        J(E) = (v(E) / 4π) * n_E(E)
//    therefore
//        n_E(E) = 4π * J(E) / v(E)
//
//    The total number density in the energy range [Emin, Emax] is:
//        n = ∫ n_E(E) dE = 4π ∫_{Emin}^{Emax} J_loc(E) / v(E) dE
//
//    We evaluate the integral numerically using trapezoidal integration on an energy
//    grid defined by #DENSITY_SPECTRUM:
//        DS_EMIN, DS_EMAX, DS_NINTERVALS, DS_ENERGY_SPACING (LOG or LINEAR).
//
// IMPLEMENTATION OVERVIEW (how it is computed)
// -------------------------------------------------------------------------------------
// - Parse and validate inputs using EarthUtil::ParseAmpsParamFile().
// - Build the energy grid:
//     Npoints = DS_NINTERVALS + 1
//     LOG spacing: E[i] = Emin * (Emax/Emin)^(i/(Npoints-1))
//     LINEAR     : E[i] = Emin + (Emax-Emin)*i/(Npoints-1)
// - Build a direction grid identical to the cutoff tool (nZenith=24, nAz=48).
// - For each point x0:
//     For each energy E:
//       - Convert energy -> rigidity R(GV)
//       - For each direction d:
//           backtrace using TraceAllowed(prm, field, x0, v0=-d, R)
//       - T(E) = allowed/Ndirs
//     Compute J_loc(E) and integrate for density.
// - MPI: dynamic master/worker scheduling over POINT indices.
//     Each task is one location; the worker returns T(E) array + integrated density.
// - Rank 0 writes two Tecplot ASCII files:
//     (1) gridless_points_density.dat
//     (2) gridless_points_spectrum.dat
//======================================================================================

#include "DensityGridless.h"

#include "CutoffRigidityGridless.h" // for namespace consistency only (no symbols used)

#include "../boundary/spectrum.h"  // ::gSpectrum and spectrum metadata writers

#include <mpi.h>
#include <cmath>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <algorithm>
#include <iostream>

// --- Field model dependencies (match CutoffRigidityGridless.cpp) ----------------------
// We intentionally use the *same* access pattern as the gridless cutoff solver to avoid
// link-time mismatches across platforms/compilers:
//   - GeopackInterface provides C++ wrappers for initializing Geopack state for a given
//     epoch and evaluating IGRF in GSM.
//   - Tsyganenko external fields are accessed via the Fortran entry points (t96_01_,
//     t04_s_) with the same signatures as in the cutoff code.
//
// IMPORTANT:
//   Do NOT call raw Geopack Fortran symbols such as recalc_ or igrf_gsm_ here. Those
//   may not be linked in your AMPS build depending on how Geopack is wrapped.
//   Using GeopackInterface is the portable approach already proven by CutoffRigidity.

#include "constants.h"
#include "constants.PlanetaryData.h"
#include "GeopackInterface.h"

namespace {

//--------------------------------------------------------------------------------------
// Small 3D vector helper (private to this file).
// We duplicate the minimal utilities from CutoffRigidityGridless.cpp intentionally so
// this module is self-contained.
//--------------------------------------------------------------------------------------

struct V3 { double x,y,z; };

static inline V3 add(const V3&a,const V3&b){return {a.x+b.x,a.y+b.y,a.z+b.z};}
static inline V3 sub(const V3&a,const V3&b){return {a.x-b.x,a.y-b.y,a.z-b.z};}
static inline V3 mul(double s,const V3&a){return {s*a.x,s*a.y,s*a.z};}
static inline double dot(const V3&a,const V3&b){return a.x*b.x+a.y*b.y+a.z*b.z;}
static inline double norm(const V3&a){return std::sqrt(dot(a,a));}
static inline V3 unit(const V3&a){double n=norm(a); return (n>0)?mul(1.0/n,a):V3{0,0,0};}

//--------------------------------------------------------------------------------------
// Constants
//--------------------------------------------------------------------------------------

static constexpr double C_LIGHT = 299792458.0;
static constexpr double QE = 1.602176634e-19;
static constexpr double AMU = 1.66053906660e-27;
static constexpr double MEV_TO_J = 1.0e6 * QE;

static constexpr double RE_KM = 6371.2;

struct DomainBoxRe {
  double xMin,xMax,yMin,yMax,zMin,zMax,rInner;
};

static DomainBoxRe ToDomainBoxRe(const EarthUtil::DomainBox& d) {
  DomainBoxRe b;
  b.xMin = d.xMin/RE_KM; b.xMax = d.xMax/RE_KM;
  b.yMin = d.yMin/RE_KM; b.yMax = d.yMax/RE_KM;
  b.zMin = d.zMin/RE_KM; b.zMax = d.zMax/RE_KM;
  b.rInner= d.rInner/RE_KM;
  return b;
}

//--------------------------------------------------------------------------------------
// Field evaluator and Boris pusher
//--------------------------------------------------------------------------------------
// We reuse the same magnetic field stack as gridless cutoff:
//   - internal: IGRF via Geopack
//   - external: Tsyganenko T96 or T05
//
// NOTE: The actual calls into Geopack/Tsyganenko are already implemented in the
// cutoff module; to keep this module independent, we replicate the evaluator code
// pattern here.
//--------------------------------------------------------------------------------------

extern "C" {
  // Tsyganenko Fortran entry points (same declarations as in CutoffRigidityGridless.cpp)
  void t96_01_(int*,double*,double*,double*,double*,double*,double*,double*,double*);
  void t04_s_ (int*,double*,double*,double*,double*,double*,double*,double*,double*);
}

// Convert epoch "YYYY-MM-DDTHH:MM" or similar to geopack time inputs.
static void ParseEpochToGeopack(const std::string& epoch,
                                int& iyear,int& iday,int& ihour,int& imin,int& isec) {
  // Very lightweight: expect at least YYYY-MM-DD.
  if (epoch.size()<10) throw std::runtime_error("Invalid epoch string: '"+epoch+"'");
  iyear = std::stoi(epoch.substr(0,4));
  int imon = std::stoi(epoch.substr(5,2));
  int idom = std::stoi(epoch.substr(8,2));

  ihour = 0; imin = 0; isec = 0;
  if (epoch.size()>=16) {
    ihour = std::stoi(epoch.substr(11,2));
    imin  = std::stoi(epoch.substr(14,2));
  }

  // Convert month/day to day-of-year.
  static const int mdays_norm[12] = {31,28,31,30,31,30,31,31,30,31,30,31};
  bool leap = ((iyear%4==0 && iyear%100!=0) || (iyear%400==0));
  int id = 0;
  for (int m=1;m<imon;m++) {
    id += mdays_norm[m-1];
    if (leap && m==2) id += 1;
  }
  id += idom;
  iday = id;
}

class cFieldEvaluator {
public:
  explicit cFieldEvaluator(const EarthUtil::AmpsParam& prm) : prm_(prm) {
    // Exactly as in CutoffRigidityGridless.cpp: initialize Geopack in GSM.
    Geopack::Init(prm_.field.epoch.c_str(),"GSM");

    for (int i=0;i<11;i++) PARMOD_[i]=0.0;
    PS_ = 0.170481; // default dipole tilt used by existing interfaces

    PARMOD_[0]=prm_.field.pdyn_nPa;
    PARMOD_[1]=prm_.field.dst_nT;
    PARMOD_[2]=prm_.field.imfBy_nT;
    PARMOD_[3]=prm_.field.imfBz_nT;

    if (EarthUtil::ToUpper(prm_.field.model)=="T05") {
      for (int i=0;i<6;i++) PARMOD_[4+i]=prm_.field.w[i];
    }
  }

  std::string Model() const { return prm_.field.model; }

  // Evaluate B at a position expressed in SI meters (GSM). Output is Tesla.
  void GetB_T(const V3& x_m, V3& B_T) const {
    // Internal field from IGRF via Geopack wrapper.
    double x_arr[3]={x_m.x,x_m.y,x_m.z};
    double b_int[3];
    Geopack::IGRF::GetMagneticField(b_int,x_arr);

    // External field from Tsyganenko (expects coordinates in Re and returns nT).
    double xRe[3]={x_m.x/_EARTH__RADIUS_, x_m.y/_EARTH__RADIUS_, x_m.z/_EARTH__RADIUS_};
    double b_ext_nT[3]={0,0,0};
    int IOPT=0;

    if (EarthUtil::ToUpper(Model())=="T96") {
      t96_01_(&IOPT,const_cast<double*>(PARMOD_),const_cast<double*>(&PS_),
              xRe+0,xRe+1,xRe+2,b_ext_nT+0,b_ext_nT+1,b_ext_nT+2);
    }
    else if (EarthUtil::ToUpper(Model())=="T05") {
      t04_s_(&IOPT,const_cast<double*>(PARMOD_),const_cast<double*>(&PS_),
             xRe+0,xRe+1,xRe+2,b_ext_nT+0,b_ext_nT+1,b_ext_nT+2);
    }
    else {
      throw std::runtime_error("Unsupported FIELD_MODEL for gridless density: '"+Model()+"' (supported: T96,T05)");
    }

    B_T.x = b_int[0] + b_ext_nT[0]*_NANO_;
    B_T.y = b_int[1] + b_ext_nT[1]*_NANO_;
    B_T.z = b_int[2] + b_ext_nT[2]*_NANO_;
  }

private:
  const EarthUtil::AmpsParam& prm_;
  double PARMOD_[11];
  double PS_;
};

// Relativistic Boris push (magnetic-only).
static void BorisStep(V3& x_m, V3& p_SI, double q_C, double m0_kg,
                      double dt, const cFieldEvaluator& field) {
  // Field evaluation is in SI meters, consistent with CutoffRigidityGridless.
  V3 B; field.GetB_T(x_m,B);

  // gamma = sqrt(1 + (p/(m c))^2)
  const double p2 = dot(p_SI,p_SI);
  const double mc = m0_kg*C_LIGHT;
  const double gamma = std::sqrt(1.0 + p2/(mc*mc));

  // t = (q dt / (2 gamma m)) B
  const double coef = (q_C*dt)/(2.0*gamma*m0_kg);
  const V3 t = mul(coef,B);
  const double t2 = dot(t,t);
  const V3 s = mul(2.0/(1.0+t2), t);

  // p- = p
  V3 pminus = p_SI;

  // p' = p- + p- x t
  V3 pprime { pminus.x + (pminus.y*t.z - pminus.z*t.y),
              pminus.y + (pminus.z*t.x - pminus.x*t.z),
              pminus.z + (pminus.x*t.y - pminus.y*t.x) };

  // p+ = p- + p' x s
  V3 pplus { pminus.x + (pprime.y*s.z - pprime.z*s.y),
             pminus.y + (pprime.z*s.x - pprime.x*s.z),
             pminus.z + (pprime.x*s.y - pprime.y*s.x) };

  p_SI = pplus;

  // Update position: v = p/(gamma m)
  const double p2n = dot(p_SI,p_SI);
  const double gamman = std::sqrt(1.0 + p2n/(mc*mc));
  const V3 v = mul(1.0/(gamman*m0_kg), p_SI);
  x_m = add(x_m, mul(dt, v));
}

// Adaptive dt selection: copy the conservative approach from cutoff.
static double ChooseDt(const EarthUtil::AmpsParam& prm,
                       double qabs_C,double m0_kg,
                       const V3& p_SI,const V3& B_T,
                       double dtMax) {
  // Limit the gyro phase advance per step.
  const double Bmag = norm(B_T);
  if (Bmag>0.0) {
    const double mc = m0_kg*C_LIGHT;
    const double p2 = dot(p_SI,p_SI);
    const double gamma = std::sqrt(1.0 + p2/(mc*mc));
    const double omega = qabs_C*Bmag/(gamma*m0_kg); // rad/s
    if (omega>0.0) {
      const double dphiMax = 0.15;
      dtMax = std::min(dtMax, dphiMax/omega);
    }
  }
  // Also clamp by user dtTrace_s.
  return std::min(dtMax, prm.numerics.dtTrace_s);
}

static bool InBoxRe(const DomainBoxRe& b, const V3& xRe) {
  return (xRe.x>=b.xMin && xRe.x<=b.xMax &&
          xRe.y>=b.yMin && xRe.y<=b.yMax &&
          xRe.z>=b.zMin && xRe.z<=b.zMax);
}

static bool HitInnerSphereRe(const DomainBoxRe& b, const V3& xRe) {
  return (norm(xRe) <= b.rInner);
}

// Convert kinetic energy [J] to rigidity [GV].
static double RigidityFromEnergy_GV(double E_J, double qabs_C, double m0_kg) {
  // Total energy: E_tot = E_kin + m c^2
  const double mc2 = m0_kg*C_LIGHT*C_LIGHT;
  const double Etot = E_J + mc2;
  const double pc = std::sqrt(std::max(0.0, Etot*Etot - mc2*mc2)); // p*c
  const double p = pc / C_LIGHT;
  const double R_V = (p*C_LIGHT)/qabs_C; // = p c / q
  return R_V * 1.0e-9;
}

// Relativistic speed from kinetic energy.
static double SpeedFromEnergy(double E_J, double m0_kg) {
  const double mc2 = m0_kg*C_LIGHT*C_LIGHT;
  const double gamma = 1.0 + E_J/mc2;
  if (gamma<=1.0) return 0.0;
  const double beta2 = 1.0 - 1.0/(gamma*gamma);
  return C_LIGHT*std::sqrt(std::max(0.0,beta2));
}

// Trace a single backtraced trajectory and classify allowed vs forbidden.
static bool TraceAllowed(const EarthUtil::AmpsParam& prm,
                         const cFieldEvaluator& field,
                         const V3& x0_m,
                         const V3& v0_unit,
                         double R_GV) {
  const DomainBoxRe box = ToDomainBoxRe(prm.domain);
  const double qabs = std::abs(prm.species.charge_e)*QE;
  const double m0 = prm.species.mass_amu*AMU;

  // Convert rigidity to momentum magnitude p = (R*1e9 * q)/c
  const double p_mag = (R_GV*1.0e9*qabs)/C_LIGHT;
  V3 p = mul(p_mag, v0_unit);
  V3 x = x0_m;

  double t=0.0;
  int nSteps=0;

  while (t < prm.numerics.maxTraceTime_s && nSteps < prm.numerics.maxSteps) {
    // classification checks in Re
    V3 xRe { x.x/1000.0/RE_KM, x.y/1000.0/RE_KM, x.z/1000.0/RE_KM };
    if (!InBoxRe(box,xRe)) return true;       // escaped -> allowed
    if (HitInnerSphereRe(box,xRe)) return false; // hit inner sphere -> forbidden

	    // dt selection
	    // IMPORTANT: evaluate B at the *current Cartesian position in meters* using the
	    // same background-field access path as CutoffRigidityGridless.
	    // (xRe is only used for geometry checks; the field evaluator expects SI meters.)
	    V3 B; field.GetB_T(x, B);
    double dt = ChooseDt(prm,qabs,m0,p,B,prm.numerics.dtTrace_s);
    if (!(dt>0.0)) dt = prm.numerics.dtTrace_s;

    BorisStep(x,p,prm.species.charge_e*QE,m0,dt,field);
    t += dt;
    ++nSteps;
  }

  // conservative: timeout -> forbidden
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
      V3 d{ st*std::cos(phi), st*std::sin(phi), std::cos(theta) };
      dirs.push_back(unit(d));
    }
  }
  return dirs;
}

// Select a deterministic subset of directions from the full direction grid.
//
// Why do we need this?
//   The density/spectrum workflow can be requested with a large energy grid.
//   Tracing the full direction grid at every energy can become expensive.
//   DS_MAX_PARTICLES provides a cap on the *total* number of trajectories per
//   observation point across all energies.
//
// Policy:
//   Let N_E be the number of energy points. We choose
//     N_dir_use = min(N_dir_default, floor(DS_MAX_PARTICLES / N_E))
//   and trace only those directions for each energy.
//
// Determinism:
//   We avoid random sampling to keep results reproducible.
static std::vector<V3> SelectDirectionsDeterministic(const std::vector<V3>& dirs, int nUse) {
  if (nUse <= 0 || nUse >= (int)dirs.size()) return dirs;

  std::vector<V3> out;
  out.reserve(nUse);

  const double step = double(dirs.size()) / double(nUse);
  for (int k=0; k<nUse; ++k) {
    const int idx = std::min((int)dirs.size()-1, (int)std::floor(k*step));
    out.push_back(dirs[idx]);
  }

  return out;
}

static std::vector<double> BuildEnergyGrid_MeV(const EarthUtil::AmpsParam& prm) {
  const int n = prm.densitySpectrum.nPoints();
  std::vector<double> E(n);
  const double Emin = prm.densitySpectrum.Emin_MeV;
  const double Emax = prm.densitySpectrum.Emax_MeV;
  if (n<2) throw std::runtime_error("Energy grid requires at least 2 points");

  if (prm.densitySpectrum.spacing==EarthUtil::DensitySpectrumParam::Spacing::LOG) {
    const double r = Emax/Emin;
    for (int i=0;i<n;i++) {
      double a = double(i)/double(n-1);
      E[i] = Emin*std::pow(r,a);
    }
  }
  else {
    for (int i=0;i<n;i++) {
      double a = double(i)/double(n-1);
      E[i] = Emin + (Emax-Emin)*a;
    }
  }
  return E;
}

// Simple trapezoid integral on non-uniform grid.
static double Trapz(const std::vector<double>& x, const std::vector<double>& f) {
  if (x.size()!=f.size() || x.size()<2) return 0.0;
  double s=0.0;
  for (size_t i=0;i+1<x.size();i++) {
    s += 0.5*(f[i]+f[i+1])*(x[i+1]-x[i]);
  }
  return s;
}

}

namespace Earth {
namespace GridlessMode {

int RunDensityAndSpectrumPoints(const EarthUtil::AmpsParam& prm) {
  int mpiRank=0, mpiSize=1;
  MPI_Comm_rank(MPI_COMM_WORLD,&mpiRank);
  MPI_Comm_size(MPI_COMM_WORLD,&mpiSize);

  if (EarthUtil::ToUpper(prm.output.mode)!="POINTS") {
    throw std::runtime_error("DENSITY_SPECTRUM currently supports OUTPUT_MODE=POINTS only");
  }

  // Build shared grids.
  const int nZenith=24;
  const int nAz=48;
  const std::vector<V3> dirs = BuildDirGrid(nZenith,nAz);
  const std::vector<double> E_MeV = BuildEnergyGrid_MeV(prm);
  const int nE = (int)E_MeV.size();

  // Enforce DS_MAX_PARTICLES (if provided) by limiting the number of directions
  // traced at each energy. See SelectDirectionsDeterministic() comment block.
  int nDirsUse = (int)dirs.size();
  if (prm.densitySpectrum.maxParticlesPerPoint > 0 && nE > 0) {
    nDirsUse = std::max(1, prm.densitySpectrum.maxParticlesPerPoint / nE);
    nDirsUse = std::min(nDirsUse, (int)dirs.size());
  }
  const std::vector<V3> dirsUse = SelectDirectionsDeterministic(dirs, nDirsUse);

  const cFieldEvaluator field(prm);

  // Rank 0 summary.
  if (mpiRank==0) {
    std::cout << "================ Gridless density & spectrum ================\n";
    std::cout << "Run ID          : " << prm.runId << "\n";
    std::cout << "Field model     : " << prm.field.model << "\n";
    std::cout << "Epoch           : " << prm.field.epoch << "\n";
    std::cout << "Species         : " << prm.species.name << " (q=" << prm.species.charge_e
              << " e, m=" << prm.species.mass_amu << " amu)\n";
    std::cout << "Energy grid     : [" << prm.densitySpectrum.Emin_MeV << ", " << prm.densitySpectrum.Emax_MeV
              << "] MeV, Npoints=" << nE << " (" << (prm.densitySpectrum.spacing==EarthUtil::DensitySpectrumParam::Spacing::LOG?"LOG":"LINEAR")
              << ")\n";
    std::cout << "Directions grid : " << dirsUse.size() << " (nZenith=" << nZenith << ", nAz=" << nAz << ")";
    if (prm.densitySpectrum.maxParticlesPerPoint > 0) {
      std::cout << " [capped by DS_MAX_PARTICLES=" << prm.densitySpectrum.maxParticlesPerPoint << "]";
    }
    std::cout << "\n";
    std::cout << "MPI ranks       : " << mpiSize << "\n";
    std::cout << "============================================================\n";
    std::cout.flush();
  }

  const int nPoints = (int)prm.output.points.size();

  // Output buffers on master.
  std::vector<double> density_m3(nPoints, 0.0);
  std::vector< std::vector<double> > T_byPoint;
  if (mpiRank==0) T_byPoint.assign(nPoints, std::vector<double>(nE, 0.0));

  // Master/worker scheduling over point indices.
  if (mpiSize==1) {
    // Serial path.
    for (int ip=0; ip<nPoints; ++ip) {
      const EarthUtil::Vec3 pk = prm.output.points[ip];
      const V3 x0_m { pk.x*1000.0, pk.y*1000.0, pk.z*1000.0 };

      std::vector<double> T(nE,0.0);
      for (int ie=0; ie<nE; ++ie) {
        const double Ej = E_MeV[ie]*MEV_TO_J;
        const double Rgv = RigidityFromEnergy_GV(Ej, std::abs(prm.species.charge_e)*QE, prm.species.mass_amu*AMU);
        int allowed=0;
        for (const auto& d : dirsUse) {
          if (TraceAllowed(prm,field,x0_m, mul(-1.0,d), Rgv)) ++allowed;
        }
        T[ie] = double(allowed)/double(dirsUse.size());
      }

      // Density integral
      std::vector<double> integrand(nE,0.0);
      std::vector<double> EjGrid(nE,0.0);
      for (int ie=0; ie<nE; ++ie) {
        const double Ej = E_MeV[ie]*MEV_TO_J;
        EjGrid[ie]=Ej;
        const double v = SpeedFromEnergy(Ej, prm.species.mass_amu*AMU);
        const double Jb = ::gSpectrum.GetSpectrum(Ej); // per Joule
        const double Jloc = T[ie]*Jb;
        integrand[ie] = (v>0.0) ? (4.0*M_PI*Jloc/v) : 0.0;
      }
      density_m3[ip] = Trapz(EjGrid, integrand);
      T_byPoint[ip] = std::move(T);
    }
  }
  else {
    // Parallel dynamic scheduling.
    const int TAG_TASK=100, TAG_RES=200;
    if (mpiRank==0) {
      int next=0;
      // seed workers
      for (int r=1;r<mpiSize;r++) {
        int idx = (next<nPoints)? next++ : -1;
        MPI_Send(&idx,1,MPI_INT,r,TAG_TASK,MPI_COMM_WORLD);
      }

      int done=0;
      while (done<nPoints) {
        // Receive header: point index + density.
        struct { int idx; double dens; } hdr;
        MPI_Status st;
        MPI_Recv(&hdr,sizeof(hdr),MPI_BYTE,MPI_ANY_SOURCE,TAG_RES,MPI_COMM_WORLD,&st);

        // Receive T array
        std::vector<double> T(nE,0.0);
        MPI_Recv(T.data(),nE,MPI_DOUBLE,st.MPI_SOURCE,TAG_RES+1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

        density_m3[hdr.idx]=hdr.dens;
        T_byPoint[hdr.idx]=std::move(T);
        done++;

        // send next task
        int idx = (next<nPoints)? next++ : -1;
        MPI_Send(&idx,1,MPI_INT,st.MPI_SOURCE,TAG_TASK,MPI_COMM_WORLD);
      }
    }
    else {
      while (true) {
        int idx=-1;
        MPI_Recv(&idx,1,MPI_INT,0,TAG_TASK,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        if (idx<0) break;

        const EarthUtil::Vec3 pk = prm.output.points[idx];
        const V3 x0_m { pk.x*1000.0, pk.y*1000.0, pk.z*1000.0 };

        std::vector<double> T(nE,0.0);
        for (int ie=0; ie<nE; ++ie) {
          const double Ej = E_MeV[ie]*MEV_TO_J;
          const double Rgv = RigidityFromEnergy_GV(Ej, std::abs(prm.species.charge_e)*QE, prm.species.mass_amu*AMU);
          int allowed=0;
          for (const auto& d : dirsUse) {
            if (TraceAllowed(prm,field,x0_m, mul(-1.0,d), Rgv)) ++allowed;
          }
          T[ie] = double(allowed)/double(dirsUse.size());
        }

        // Density integral
        std::vector<double> integrand(nE,0.0);
        std::vector<double> EjGrid(nE,0.0);
        for (int ie=0; ie<nE; ++ie) {
          const double Ej = E_MeV[ie]*MEV_TO_J;
          EjGrid[ie]=Ej;
          const double v = SpeedFromEnergy(Ej, prm.species.mass_amu*AMU);
          const double Jb = ::gSpectrum.GetSpectrum(Ej);
          const double Jloc = T[ie]*Jb;
          integrand[ie] = (v>0.0) ? (4.0*M_PI*Jloc/v) : 0.0;
        }
        const double dens = Trapz(EjGrid, integrand);

        struct { int idx; double dens; } hdr{idx,dens};
        MPI_Send(&hdr,sizeof(hdr),MPI_BYTE,0,TAG_RES,MPI_COMM_WORLD);
        MPI_Send(T.data(),nE,MPI_DOUBLE,0,TAG_RES+1,MPI_COMM_WORLD);
      }
    }
  }

  // Rank 0 writes Tecplot outputs.
  if (mpiRank==0) {
    // File 1: densities
    {
      std::ofstream out("gridless_points_density.dat");
      out << "TITLE=\"Gridless energetic particle density (POINTS)\"\n";
      out << "VARIABLES=\"X_km\" \"Y_km\" \"Z_km\" \"N_m^-3\" \"N_cm^-3\"\n";
      out << "ZONE T=\"density\" I=" << nPoints << " F=POINT\n";
      for (int ip=0; ip<nPoints; ++ip) {
        const auto& p0 = prm.output.points[ip];
        const double n_m3 = density_m3[ip];
        const double n_cm3 = n_m3*1.0e-6;
        out << p0.x << " " << p0.y << " " << p0.z << " " << n_m3 << " " << n_cm3 << "\n";
      }
    }

    // File 2: spectra
    {
      std::ofstream out("gridless_points_spectrum.dat");
      out << "TITLE=\"Gridless energetic particle spectrum (POINTS)\"\n";
      out << "VARIABLES=\"E_MeV\" \"T\" \"J_boundary_perMeV\" \"J_local_perMeV\"\n";

      // Record spectrum definition for provenance.
      for (int ip=0; ip<nPoints; ++ip) {
        std::ostringstream zt;
        zt << "P" << ip;
        out << "ZONE T=\"" << zt.str() << "\" I=" << nE << " F=POINT\n";
        for (int ie=0; ie<nE; ++ie) {
          const double Ej = E_MeV[ie]*MEV_TO_J;
          const double T = T_byPoint[ip][ie];
          const double Jb_perMeV = ::gSpectrum.GetSpectrumPerMeV(E_MeV[ie]);
          const double Jloc_perMeV = T*Jb_perMeV;
          out << E_MeV[ie] << " " << T << " " << Jb_perMeV << " " << Jloc_perMeV << "\n";
        }
      }
    }
  }

  return 0;
}

}
}
