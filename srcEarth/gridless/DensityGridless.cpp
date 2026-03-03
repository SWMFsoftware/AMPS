//======================================================================================
// DensityGridless.cpp
//======================================================================================
/**
 * \file DensityGridless.cpp
 *
 * GRIDLESS ENERGETIC-PARTICLE DENSITY + SPECTRUM (POINTS MODE)
 * -----------------------------------------------------------
 *
 * This module computes, at a set of user-provided observation points, (i) the *local*
 * differential energy spectrum and (ii) the corresponding *total energetic-particle number
 * density* implied by that local spectrum.
 *
 * The implementation is intentionally "gridless": the magnetic field is evaluated on the fly
 * using the same background-field access pattern as the gridless cutoff-rigidity solver
 * (IGRF + Tsyganenko external field), and particle trajectories are integrated directly in SI.
 *
 * -----------------------------------------------------------------------------
 * 1. INPUTS (from AMPS_PARAM.in)
 * -----------------------------------------------------------------------------
 * The calculation is activated by:
 *   #CALCULATION_MODE
 *     CALC_TARGET       DENSITY_SPECTRUM
 *     FIELD_EVAL_METHOD GRIDLESS
 *
 * Energy grid controls (this section is REQUIRED for CALC_TARGET=DENSITY_SPECTRUM):
 *   #DENSITY_SPECTRUM
 *     DS_EMIN           [MeV/n]
 *     DS_EMAX           [MeV/n]
 *     DS_NINTERVALS     number of energy intervals (Npoints = NINTERVALS + 1)
 *     DS_ENERGY_SPACING LOG | LINEAR
 *     DS_MAX_PARTICLES  (optional) cap on total trajectory evaluations per point
 *     DS_MAX_TRAJ_TIME  (optional) cap on integration time per trajectory [s]
 *
 * Boundary spectrum:
 *   #SPECTRUM  ... (parsed into amps::gSpectrum, see boundary/spectrum.*)
 * The function amps::gSpectrum.GetSpectrum(E[J]) must return the *boundary differential
 * intensity* J_b(E) in units of [m^-2 s^-1 sr^-1 J^-1] (or an equivalent consistent system).
 *
 * Geometry / stopping conditions:
 *   - outer boundary: rectangular domain box (escape = ALLOWED)
 *   - inner boundary: loss sphere (atmosphere proxy; hit = FORBIDDEN)
 *   - time/step caps: from #NUMERICAL plus DS_MAX_TRAJ_TIME override for this mode
 *
 * -----------------------------------------------------------------------------
 * 2. THEORY (what is computed and why)
 * -----------------------------------------------------------------------------
 * We assume that far outside the magnetosphere the energetic particle population is an
 * (approximately) *isotropic* distribution characterized by a boundary differential intensity
 * J_b(E). The magnetosphere acts as a direction-dependent filter: for a particle at an
 * observation point x0 and kinetic energy E, only some asymptotic directions connect to the
 * boundary without intersecting the loss surface.
 *
 * 2.1 Directional transmissivity
 * For each observation point x0 and energy E we sample N_dirs directions {\hat{u}_k} and
 * classify each direction via backtracing. Define:
 *
 *   A_k(E; x0) = 1  if direction k is ALLOWED (escapes outer box)
 *             = 0  otherwise (hits loss sphere or times out).
 *
 * The (isotropic) transmissivity is approximated by a simple directional average:
 *
 *   T(E; x0) \approx (1/N_dirs) \sum_{k=1}^{N_dirs} A_k(E; x0) .
 *
 * Notes:
 *   - This is the same physical classifier used in cutoff-rigidity: "allowed vs forbidden".
 *   - The direction set is deterministic for reproducibility.
 *   - If DS_MAX_PARTICLES limits work, we deterministically subsample directions.
 *
 * 2.2 Local spectrum
 * Under the isotropy assumption, the local differential intensity is modeled as:
 *
 *   J_loc(E; x0) = T(E; x0) * J_b(E).
 *
 * This corresponds to a "grey" (direction-averaged) transmissivity applied to the boundary
 * spectrum.
 *
 * 2.3 Number density from differential intensity
 * Let n(E) = dn/dE be the differential number density [m^-3 J^-1]. For an isotropic
 * distribution, the differential intensity and density are related by:
 *
 *   J(E) = (v(E) / 4\pi) * n(E),
 *
 * where v(E) is the particle speed at kinetic energy E. Therefore:
 *
 *   n(E) = (4\pi / v(E)) * J(E),
 *   n_tot = \int_{Emin}^{Emax} n(E) dE = 4\pi \int_{Emin}^{Emax} J(E)/v(E) dE .
 *
 * We apply this to the local spectrum:
 *
 *   n_tot(x0) = 4\pi \int_{Emin}^{Emax} J_loc(E; x0)/v(E) dE
 *             = 4\pi \int_{Emin}^{Emax} T(E; x0) * J_b(E) / v(E) dE .
 *
 * Numerical quadrature uses the trapezoidal rule on the user-specified energy grid.
 *
 * 2.4 Relativistic kinematics (kinetic energy -> speed)
 * For a particle of rest mass m0:
 *   gamma = 1 + E/(m0 c^2)
 *   beta  = sqrt(1 - 1/gamma^2)
 *   v     = beta * c
 *
 * This is used consistently in both trajectory integration (Boris pusher in momentum form)
 * and in the density integral.
 *
 * -----------------------------------------------------------------------------
 * 3. OUTPUTS
 * -----------------------------------------------------------------------------
 * Two Tecplot ASCII files are written (rank 0):
 *
 *  (1) gridless_points_density.dat
 *      Variables: X_km Y_km Z_km N_m3 N_cm3
 *
 *  (2) gridless_points_spectrum.dat
 *      One ZONE per observation point.
 *      Variables: E_MeV  T  J_boundary_perMeV  J_local_perMeV
 *
 * -----------------------------------------------------------------------------
 * 4. IMPLEMENTATION OVERVIEW (algorithmic steps)
 * -----------------------------------------------------------------------------
 * For each observation point x0:
 *   (a) Build energy grid {E_i} with LOG or LINEAR spacing in [Emin,Emax].
 *   (b) For each E_i:
 *       - choose direction subset (possibly limited by DS_MAX_PARTICLES)
 *       - for each direction, backtrace and classify allowed/forbidden
 *       - compute T(E_i) as allowed fraction
 *       - compute boundary spectrum J_b(E_i) = gSpectrum.GetSpectrum(E_i[J])
 *       - compute local spectrum  J_loc(E_i) = T(E_i)*J_b(E_i)
 *   (c) Integrate n_tot = 4*pi * sum_i trapezoid( J_loc(E_i)/v(E_i) ).
 *   (d) Output: density per point and per-point spectrum ZONE.
 *
 * The heavy lifting per trajectory is identical in spirit to the cutoff-rigidity gridless
 * solver: same field evaluation, same geometry checks, same time-step constraints.
 */
#include "pic.h"
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
#include "DipoleInterface.h"

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
    // Exactly as in CutoffRigidityGridless.cpp: Geopack is required only for
    // Tsyganenko models (IGRF + coordinate setup). For the analytic DIPOLE model
    // we keep the calculation self-contained and do not call Geopack at all.
    if (EarthUtil::ToUpper(Model())!="DIPOLE") {
      Geopack::Init(prm_.field.epoch.c_str(),"GSM");
    }

    // Configure analytic dipole parameters (used only when FIELD_MODEL=DIPOLE).
    Earth::GridlessMode::Dipole::SetMomentScale(prm_.field.dipoleMoment_Me);
    Earth::GridlessMode::Dipole::SetTiltDeg(prm_.field.dipoleTilt_deg);

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
    // Dipole branch: internal field only, analytic.
    if (EarthUtil::ToUpper(Model())=="DIPOLE") {
      double x_arr[3]={x_m.x,x_m.y,x_m.z};
      double b_arr[3];
      Earth::GridlessMode::Dipole::GetB_Tesla(x_arr,b_arr);
      B_T.x=b_arr[0]; B_T.y=b_arr[1]; B_T.z=b_arr[2];
      return;
    }

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
      throw std::runtime_error("Unsupported FIELD_MODEL for gridless density: '"+Model()+"' (supported: T96,T05,DIPOLE)");
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

//-------------------------------------------------------------------------------------------------
// Convenience overloads used by the dipole analytic-comparison utilities.
//-------------------------------------------------------------------------------------------------
// In the core tight-loop trajectory integrator we keep the low-level helpers:
//
//   RigidityFromEnergy_GV(E_J, qabs_C, m0_kg)
//   SpeedFromEnergy(E_J, m0_kg)
//
// because they avoid any dependency on a higher-level parameter object.
//
// Some of the analytic reference writers (used only in nightly test mode) naturally have access to
// the full parsed run configuration (EarthUtil::AmpsParam).  To keep those call sites readable, we
// provide small wrappers that extract the charge and rest mass from prm.species and forward to the
// low-level implementations.
//
// Units:
//   - E_J      : kinetic energy [J]
//   - charge_e : particle charge in units of elementary charge (e)
//   - mass_amu : particle rest mass in atomic mass units (amu)
//   - Returns rigidity [GV] and speed [m/s].
//
// NOTE: These overloads intentionally live in the same anonymous namespace as the low-level helpers
// so the compiler can inline them and so we do not export any new symbols.
static double RigidityFromEnergy_GV(double E_J, const EarthUtil::AmpsParam& prm) {
  const double qabs_C = std::abs(prm.species.charge_e) * QE;
  const double m0_kg  = prm.species.mass_amu * AMU;
  return RigidityFromEnergy_GV(E_J, qabs_C, m0_kg);
}

static double SpeedFromEnergy(double E_J, const EarthUtil::AmpsParam& prm) {
  const double m0_kg = prm.species.mass_amu * AMU;
  return SpeedFromEnergy(E_J, m0_kg);
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

  // Per-trajectory integration time limit.
  //
  // Rationale:
  //   Density/spectrum calculations can trace a large number of trajectories
  //   (N_energy * N_dir per point). A dedicated DS_MAX_TRAJ_TIME provides a
  //   convenient knob to cap work per trajectory without changing the cutoff
  //   tool's global #NUMERICAL MAX_TRACE_TIME setting.
  //
  // Policy:
  //   - If DS_MAX_TRAJ_TIME > 0: use it.
  //   - Else: fall back to #NUMERICAL MAX_TRACE_TIME.
  const double maxTrajTime_s = (prm.densitySpectrum.maxTrajTime_s > 0.0)
                                 ? prm.densitySpectrum.maxTrajTime_s
                                 : prm.numerics.maxTraceTime_s;

  while (t < maxTrajTime_s && nSteps < prm.numerics.maxSteps) {
    // classification checks in Re
    V3 xRe { x.x/1000.0/RE_KM, x.y/1000.0/RE_KM, x.z/1000.0/RE_KM };
    if (!InBoxRe(box,xRe)) return true;       // escaped -> allowed
    if (HitInnerSphereRe(box,xRe)) return false; // hit inner sphere -> forbidden

	    // dt selection
	    // IMPORTANT: evaluate B at the *current Cartesian position in meters* using the
	    // same background-field access path as CutoffRigidityGridless.
	    // (xRe is only used for geometry checks; the field evaluator expects SI meters.)
	    V3 B; field.GetB_T(x, B);
    double dt = 100*ChooseDt(prm,qabs,m0,p,B,prm.numerics.dtTrace_s);
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

//====================================================================================
// Progress bar helper for density/spectrum calculations
//
// DESIGN RATIONALE
//   The cutoff-rigidity module already has a sophisticated progress bar; this module
//   (density/spectrum) previously had none, making long runs opaque to the user.
//
//   We follow the same visual conventions as the cutoff progress bar:
//     - ASCII bar [####----] with percentage
//     - Completed/total task counts
//     - ETA based on elapsed walltime and current rate
//     - Time-throttled printing (at most once per second) to avoid stdout overhead
//     - Master (rank 0) only prints; workers never touch stdout
//
//   The progress tracker is a lightweight struct so it can be instantiated locally
//   in each driver function without polluting the file-level namespace.
//====================================================================================
struct ProgressBar {
  double tStart;          // wall-clock reference
  double tLastPrint;      // last time we printed
  int    mpiRank;         // only rank 0 prints
  std::string label;      // e.g. "[POINTS]" or "[SHELL 2/5 alt=450km]"

  ProgressBar(int rank, const std::string& lbl)
    : tStart(MPI_Wtime()), tLastPrint(-1.0), mpiRank(rank), label(lbl) {}

  void setLabel(const std::string& lbl) { label = lbl; }

  // Print progress if rank==0 and at least 1 second has elapsed since last print.
  // done/total are the number of completed vs total work items.
  void update(long long done, long long total) {
    if (mpiRank != 0) return;

    const double tNow = MPI_Wtime();
    if (tLastPrint >= 0.0 && (tNow - tLastPrint) < 1.0) return;
    tLastPrint = tNow;

    const double frac = (total > 0) ? double(done) / double(total) : 1.0;

    // ETA
    const double elapsed = tNow - tStart;
    const double rate = (elapsed > 0.0) ? double(done) / elapsed : 0.0;
    double eta_s = -1.0;
    if (rate > 0.0 && done < total) eta_s = double(total - done) / rate;

    // Format HH:MM:SS
    auto fmt_hms = [](double s) -> std::string {
      if (s < 0.0) return std::string("--:--:--");
      long long is = (long long)std::llround(s);
      long long hh = is / 3600; is -= hh * 3600;
      long long mm = is / 60;   is -= mm * 60;
      long long ss = is;
      char buf[64];
      std::snprintf(buf, sizeof(buf), "%02lld:%02lld:%02lld", hh, mm, ss);
      return std::string(buf);
    };

    const int barW = 36;
    const int filled = (int)std::floor(frac * barW + 0.5);

    std::cout << label << " [rank 0] [";
    for (int i = 0; i < barW; i++) std::cout << (i < filled ? "#" : "-");
    std::cout << "] ";

    std::cout << std::fixed;
    std::cout.precision(1);
    std::cout << (frac * 100.0) << "%  ";
    std::cout << "(" << done << "/" << total << ")  "
              << "ETA " << fmt_hms(eta_s) << "\n";
    std::cout.flush();
  }

  // Force a final 100% line (always prints regardless of throttle).
  void finish(long long total) {
    if (mpiRank != 0) return;
    tLastPrint = -1.0; // force print
    update(total, total);
  }
};

}



//======================================================================================
// DIPOLE ANALYTIC REFERENCE (DENSITY, POINTS)
//======================================================================================
// This helper writes a Tecplot file comparing the numerically computed energetic
// particle number density (from the full backtracing transmissivity calculation)
// against a simple *analytic* dipole benchmark.
//
// Benchmark model (hard-cutoff approximation):
//   1) Use the vertical Størmer cutoff rigidity Rv(λ,r) as a sharp threshold.
//   2) Assume transmissivity T(E)=1 for R(E)≥Rv and 0 otherwise.
//   3) Fold the boundary spectrum: J_loc(E)=T(E) J_b(E).
//   4) Convert intensity to number density via n=4π\int J_loc(E)/v(E) dE.
//
// This intentionally ignores the penumbra and direction-dependent cutoffs. It is
// nonetheless an excellent regression reference for unit/geometry conversions and
// the numerical quadrature in the density tool.
//======================================================================================
static void WriteTecplotPoints_DipoleAnalyticDensityCompare(const EarthUtil::AmpsParam& prm,
                                                            const std::vector<EarthUtil::Vec3>& points,
                                                            const std::vector<double>& n_num_m3);

namespace Earth {
namespace GridlessMode {

//======================================================================================
// INTERNAL DRIVER: POINTS MODE
//======================================================================================
// This is the original implementation: a list of explicit observation points is
// processed (with MPI dynamic scheduling if available) and two Tecplot files are
// written:
//   - gridless_points_density.dat
//   - gridless_points_spectrum.dat
//
// The logic is kept as a standalone helper so we can reuse the same per-point kernel
// for SHELLS mode without mixing the output logic.
static int RunDensityAndSpectrum_POINTS(const EarthUtil::AmpsParam& prm) {
  int mpiRank=0, mpiSize=1;
  MPI_Comm_rank(MPI_COMM_WORLD,&mpiRank);
  MPI_Comm_size(MPI_COMM_WORLD,&mpiSize);

  if (EarthUtil::ToUpper(prm.output.mode)!="POINTS") {
    throw std::runtime_error("Internal error: RunDensityAndSpectrum_POINTS called for non-POINTS mode");
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
    ProgressBar prog(mpiRank, "[POINTS]");

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

      prog.update((long long)(ip + 1), (long long)nPoints);
    }

    prog.finish((long long)nPoints);
  }
  else {
    // Parallel dynamic scheduling.
    const int TAG_TASK=100, TAG_RES=200;
    if (mpiRank==0) {
      ProgressBar prog(mpiRank, "[POINTS]");

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

        prog.update((long long)done, (long long)nPoints);

        // send next task
        int idx = (next<nPoints)? next++ : -1;
        MPI_Send(&idx,1,MPI_INT,st.MPI_SOURCE,TAG_TASK,MPI_COMM_WORLD);
      }

      prog.finish((long long)nPoints);
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

    //------------------------------------------------------------------------------------------
    // Nightly test mode: write an analytic-vs-numeric density comparison for the DIPOLE field.
    //
    // IMPORTANT:
    //  - This comparison is intentionally generated only in nightly mode to avoid extra I/O in
    //    production runs.
    //  - The analytic reference used here is a "hard-cutoff" benchmark derived from the vertical
    //    Størmer cutoff approximation in a centered dipole (with optional dipole tilt/moment).
    //  - The numerical density is the one already computed above (density_m3).
    //------------------------------------------------------------------------------------------
#if _PIC_NIGHTLY_TEST_MODE_ == _PIC_MODE_ON_ 
    if (EarthUtil::ToUpper(prm.field.model)=="DIPOLE") {
      // Writes: density_gridless_points_dipole_compare.dat
      WriteTecplotPoints_DipoleAnalyticDensityCompare(prm,prm.output.points,density_m3);
    }
#endif

  }

  return 0;
}

//======================================================================================
// INTERNAL DRIVER: SHELLS MODE
//======================================================================================
// For each requested shell altitude, we discretize the spherical surface into a
// deterministic lon/lat grid and evaluate the density model at each grid point.
//
// OUTPUT (per shell altitude): a single Tecplot file with multiple ZONEs:
//   - ZONE 1: total density integrated over the full energy band
//   - ZONE i: density contribution from energy channel i-2
//
// Each ZONE is written in POINT format with the same spatial coordinates, so the
// zones can be loaded together and visualized as separate scalar fields.
//
// IMPORTANT ABOUT "ENERGY CHANNELS"
//   The energy grid contains Npoints = DS_NINTERVALS+1 values {E_i}.
//   Channel k corresponds to interval [E_k, E_{k+1}] and its density contribution is
//     n_k = 4*pi * \int_{E_k}^{E_{k+1}} J_loc(E)/v(E) dE
//   which we approximate by a single trapezoid over the two endpoints.
static std::vector<EarthUtil::Vec3> BuildShellPoints_km(double alt_km, double res_deg) {
  // Radius of the shell in km.
  const double r_km = RE_KM + alt_km;

  // We intentionally match the Tecplot shell-grid convention used by
  // CutoffRigidityGridless.cpp so that shell outputs from different tools
  // are directly comparable.
  //
  //  - Longitude: periodic; generate nLon samples on [0, 360) degrees.
  //    Do NOT include 360 deg to avoid duplicating the seam.
  //  - Latitude: include both poles; generate nLat samples on [-90, 90] degrees.
  //
  // Using floor( ... + 0.5 ) matches the cutoff-rigidity writer.
  const int nLon = std::max(1, (int)std::floor(360.0 / res_deg + 0.5));
  const int nLat = std::max(2, (int)std::floor(180.0 / res_deg + 0.5) + 1);

  std::vector<EarthUtil::Vec3> pts;
  pts.reserve((size_t)nLat * (size_t)nLon);

  for (int ilat = 0; ilat < nLat; ++ilat) {
    // Match cutoff rigidity: lat = -90 + res_deg*j (with the last sample clamped at +90)
    double lat_deg = -90.0 + res_deg * ilat;
    if (lat_deg > 90.0) lat_deg = 90.0;
    const double lat = lat_deg * M_PI / 180.0;
    const double clat = std::cos(lat);
    const double slat = std::sin(lat);
    for (int ilon = 0; ilon < nLon; ++ilon) {
      // Match cutoff rigidity: lon = res_deg*i, i in [0, nLon)
      const double lon_deg = res_deg * ilon;
      const double lon = lon_deg * M_PI / 180.0;
      const double clon = std::cos(lon);
      const double slon = std::sin(lon);
      // GSM Cartesian convention used by the rest of the tool: X sunward, Z north.
      // This spherical parameterization is purely geometric; it assumes the shell is
      // centered at Earth's center in GSM coordinates.
      EarthUtil::Vec3 p;
      p.x = r_km * clat * clon;
      p.y = r_km * clat * slon;
      p.z = r_km * slat;
      pts.push_back(p);
    }
  }
  return pts;
}

static int RunDensityAndSpectrum_SHELLS(const EarthUtil::AmpsParam& prm) {
  int mpiRank=0, mpiSize=1;
  MPI_Comm_rank(MPI_COMM_WORLD,&mpiRank);
  MPI_Comm_size(MPI_COMM_WORLD,&mpiSize);

  if (EarthUtil::ToUpper(prm.output.mode)!="SHELLS") {
    throw std::runtime_error("Internal error: RunDensityAndSpectrum_SHELLS called for non-SHELLS mode");
  }
  if (prm.output.shellAlt_km.empty()) {
    throw std::runtime_error("OUTPUT_MODE=SHELLS requires at least one altitude in SHELL_ALTITUDES");
  }

  // Shared direction and energy grids.
  const int nZenith=24;
  const int nAz=48;
  const std::vector<V3> dirs = BuildDirGrid(nZenith,nAz);
  const std::vector<double> E_MeV = BuildEnergyGrid_MeV(prm);
  const int nE = (int)E_MeV.size();
  const int nIntervals = std::max(0, prm.densitySpectrum.nIntervals);
  if (nE != nIntervals + 1) {
    throw std::runtime_error("Internal inconsistency: energy grid size != DS_NINTERVALS+1");
  }

  int nDirsUse = (int)dirs.size();
  if (prm.densitySpectrum.maxParticlesPerPoint > 0 && nE > 0) {
    nDirsUse = std::max(1, prm.densitySpectrum.maxParticlesPerPoint / nE);
    nDirsUse = std::min(nDirsUse, (int)dirs.size());
  }
  const std::vector<V3> dirsUse = SelectDirectionsDeterministic(dirs, nDirsUse);

  const cFieldEvaluator field(prm);

  if (mpiRank==0) {
    std::cout << "================ Gridless density & spectrum (SHELLS) ================\n";
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
    std::cout << "Shells          : " << prm.output.shellAlt_km.size() << " altitude(s), res=" << prm.output.shellRes_deg << " deg\n";
    std::cout << "MPI ranks       : " << mpiSize << "\n";
    std::cout << "======================================================================\n";
    std::cout.flush();
  }

  // Process each shell independently: this keeps memory bounded and also matches
  // the requested output semantics (one file per shell).
  for (double alt_km : prm.output.shellAlt_km) {
    // Build the shell point list in km and convert to meters when tracing.
    const std::vector<EarthUtil::Vec3> shellPts_km = BuildShellPoints_km(alt_km, prm.output.shellRes_deg);
    const int nPts = (int)shellPts_km.size();

    // Sanity check: the shell point list must be a complete structured lon/lat grid
    // with the exact same (nLon,nLat) definition used by the cutoff-rigidity tool.
    const int nLon_expected = std::max(1, (int)std::floor(360.0 / prm.output.shellRes_deg + 0.5));
    const int nLat_expected = std::max(2, (int)std::floor(180.0 / prm.output.shellRes_deg + 0.5) + 1);
    if (nPts != nLon_expected * nLat_expected) {
      std::ostringstream oss;
      oss << "Internal error: shell point grid size mismatch: got nPts=" << nPts
          << ", expected nLon*nLat=" << (nLon_expected * nLat_expected)
          << " (nLon=" << nLon_expected << ", nLat=" << nLat_expected
          << ", res_deg=" << prm.output.shellRes_deg << ")";
      throw std::runtime_error(oss.str());
    }

    // Storage on rank 0.
    std::vector<double> nTot_m3;
    std::vector< std::vector<double> > nChan_m3; // [interval][point]
    if (mpiRank==0) {
      nTot_m3.assign(nPts, 0.0);
      nChan_m3.assign(nIntervals, std::vector<double>(nPts, 0.0));
    }

    // MPI scheduling over shell point indices (same pattern as POINTS).
    if (mpiSize==1) {
      // Build shell-specific label: "[SHELL 2/5 alt=450km]"
      {
        int shellIdx = 0;
        for (size_t si = 0; si < prm.output.shellAlt_km.size(); ++si) {
          if (std::fabs(prm.output.shellAlt_km[si] - alt_km) < 0.01) { shellIdx = (int)si; break; }
        }
        std::ostringstream lbl;
        lbl << "[SHELL " << (shellIdx + 1) << "/" << prm.output.shellAlt_km.size()
            << " alt=" << alt_km << "km]";
        ProgressBar prog(mpiRank, lbl.str());

        for (int ip=0; ip<nPts; ++ip) {
          const auto& pk = shellPts_km[ip];
          const V3 x0_m{pk.x*1000.0, pk.y*1000.0, pk.z*1000.0};

          // Compute transmissivity T(E) on the energy grid.
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

          // Build per-energy integrand g(E)=4*pi*J_loc(E)/v(E) and integrate.
          std::vector<double> g(nE,0.0);
          std::vector<double> EjGrid(nE,0.0);
          for (int ie=0; ie<nE; ++ie) {
            const double Ej = E_MeV[ie]*MEV_TO_J;
            EjGrid[ie]=Ej;
            const double v = SpeedFromEnergy(Ej, prm.species.mass_amu*AMU);
            const double Jb = ::gSpectrum.GetSpectrum(Ej);
            const double Jloc = T[ie]*Jb;
            g[ie] = (v>0.0) ? (4.0*M_PI*Jloc/v) : 0.0;
          }
          const double nTot = Trapz(EjGrid, g);
          nTot_m3[ip] = nTot;

          // Energy-channel contributions: one trapezoid per interval.
          for (int ic=0; ic<nIntervals; ++ic) {
            const double dE = EjGrid[ic+1] - EjGrid[ic];
            const double nC = 0.5*(g[ic] + g[ic+1]) * dE;
            nChan_m3[ic][ip] = nC;
          }

          prog.update((long long)(ip + 1), (long long)nPts);
        }

        prog.finish((long long)nPts);
      }
    }
    else {
      // Parallel: rank 0 distributes point indices.
      const int TAG_TASK=3100, TAG_RES=3200;
      if (mpiRank==0) {
        // Build shell-specific label.
        int shellIdx = 0;
        for (size_t si = 0; si < prm.output.shellAlt_km.size(); ++si) {
          if (std::fabs(prm.output.shellAlt_km[si] - alt_km) < 0.01) { shellIdx = (int)si; break; }
        }
        std::ostringstream lbl;
        lbl << "[SHELL " << (shellIdx + 1) << "/" << prm.output.shellAlt_km.size()
            << " alt=" << alt_km << "km]";
        ProgressBar prog(mpiRank, lbl.str());

        int next=0;
        for (int r=1;r<mpiSize;r++) {
          int idx = (next<nPts)? next++ : -1;
          MPI_Send(&idx,1,MPI_INT,r,TAG_TASK,MPI_COMM_WORLD);
        }
        int done=0;
        while (done<nPts) {
          // Header: point index + total density
          struct { int idx; double nTot; } hdr;
          MPI_Status st;
          MPI_Recv(&hdr,sizeof(hdr),MPI_BYTE,MPI_ANY_SOURCE,TAG_RES,MPI_COMM_WORLD,&st);

          // Channel densities
          std::vector<double> buf(nIntervals,0.0);
          MPI_Recv(buf.data(),nIntervals,MPI_DOUBLE,st.MPI_SOURCE,TAG_RES+1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

          nTot_m3[hdr.idx]=hdr.nTot;
          for (int ic=0; ic<nIntervals; ++ic) nChan_m3[ic][hdr.idx] = buf[ic];
          done++;

          prog.update((long long)done, (long long)nPts);

          int idx = (next<nPts)? next++ : -1;
          MPI_Send(&idx,1,MPI_INT,st.MPI_SOURCE,TAG_TASK,MPI_COMM_WORLD);
        }

        prog.finish((long long)nPts);
      }
      else {
        while (true) {
          int idx=-1;
          MPI_Recv(&idx,1,MPI_INT,0,TAG_TASK,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
          if (idx<0) break;

          const auto& pk = shellPts_km[idx];
          const V3 x0_m{pk.x*1000.0, pk.y*1000.0, pk.z*1000.0};

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

          std::vector<double> g(nE,0.0);
          std::vector<double> EjGrid(nE,0.0);
          for (int ie=0; ie<nE; ++ie) {
            const double Ej = E_MeV[ie]*MEV_TO_J;
            EjGrid[ie]=Ej;
            const double v = SpeedFromEnergy(Ej, prm.species.mass_amu*AMU);
            const double Jb = ::gSpectrum.GetSpectrum(Ej);
            const double Jloc = T[ie]*Jb;
            g[ie] = (v>0.0) ? (4.0*M_PI*Jloc/v) : 0.0;
          }
          const double nTot = Trapz(EjGrid, g);
          std::vector<double> nCh(nIntervals,0.0);
          for (int ic=0; ic<nIntervals; ++ic) {
            const double dE = EjGrid[ic+1] - EjGrid[ic];
            nCh[ic] = 0.5*(g[ic] + g[ic+1]) * dE;
          }

          struct { int idx; double nTot; } hdr{idx,nTot};
          MPI_Send(&hdr,sizeof(hdr),MPI_BYTE,0,TAG_RES,MPI_COMM_WORLD);
          MPI_Send(nCh.data(),nIntervals,MPI_DOUBLE,0,TAG_RES+1,MPI_COMM_WORLD);
        }
      }
    }

    // Rank 0 output: one Tecplot file per shell.
    if (mpiRank==0) {
      std::ostringstream fn;
      fn << "gridless_shell_" << (int)std::round(alt_km) << "km_density_channels.dat";
      std::ofstream out(fn.str());
      out << "TITLE=\"Gridless energetic particle density (SHELL alt=" << alt_km << " km)\"\n";
      // Match the structured (I,J) Tecplot layout used by the cutoff-rigidity tool.
      // We write lon/lat grids (not X/Y/Z) because Tecplot structured grids are
      // naturally indexed in 2D. The underlying physical shell is still a sphere
      // centered at Earth; lon/lat are just a convenient parameterization.
      out << "VARIABLES=\"lon_deg\" \"lat_deg\" \"N_m^-3\" \"N_cm^-3\"\n";

      // Grid dimensions (must match BuildShellPoints_km).
      const int nLon = std::max(1, (int)std::floor(360.0 / prm.output.shellRes_deg + 0.5));
      const int nLat = std::max(2, (int)std::floor(180.0 / prm.output.shellRes_deg + 0.5) + 1);

      // Zone 1: total density on the structured lon/lat grid.
      // Ordering must match cutoff rigidity: k = i + nLon*j.
      out << "ZONE T=\"TotalDensity\" I=" << nLon << " J=" << nLat << " F=POINT\n";
      for (int j=0; j<nLat; ++j) {
        double lat = -90.0 + prm.output.shellRes_deg * j;
        if (lat > 90.0) lat = 90.0;
        for (int i=0; i<nLon; ++i) {
          const int k = i + nLon * j;
          const double lon = prm.output.shellRes_deg * i;
          const double n_m3 = nTot_m3[k];
          const double n_cm3 = n_m3 * 1.0e-6;
          out << lon << " " << lat << " " << n_m3 << " " << n_cm3 << "\n";
        }
      }

      // Zones 2..: per-channel densities.
      // Each channel corresponds to the energy interval [E_i, E_{i+1}] on the
      // DS energy grid. We store each channel as its own Tecplot zone so that
      // Tecplot users can toggle visibility / compute derived quantities per band.
      for (int ic=0; ic<nIntervals; ++ic) {
        const double Elo = E_MeV[ic];
        const double Ehi = E_MeV[ic+1];
        out << "ZONE T=\"Chan" << ic << "_" << Elo << "_" << Ehi << "MeV\" I=" << nLon << " J=" << nLat << " F=POINT\n";
        for (int j=0; j<nLat; ++j) {
          double lat = -90.0 + prm.output.shellRes_deg * j;
          if (lat > 90.0) lat = 90.0;
          for (int i=0; i<nLon; ++i) {
            const int k = i + nLon * j;
            const double lon = prm.output.shellRes_deg * i;
            const double n_m3 = nChan_m3[ic][k];
            const double n_cm3 = n_m3 * 1.0e-6;
            out << lon << " " << lat << " " << n_m3 << " " << n_cm3 << "\n";
          }
        }
      }
    }
  }

  return 0;
}

//======================================================================================
// PUBLIC ENTRY POINT
//======================================================================================
int RunDensityAndSpectrum(const EarthUtil::AmpsParam& prm) {
  const std::string mode = EarthUtil::ToUpper(prm.output.mode);
  if (mode=="POINTS") return RunDensityAndSpectrum_POINTS(prm);
  if (mode=="SHELLS") return RunDensityAndSpectrum_SHELLS(prm);
  throw std::runtime_error("Unsupported OUTPUT_MODE for DENSITY_SPECTRUM: '"+prm.output.mode+"' (supported: POINTS,SHELLS)");
}

}
}


  // NOTE:
  // This file historically accumulated two independent implementations of the
  // dipole analytic density comparison Tecplot writer with the same function name.
  // Unqualified calls then became ambiguous at compile time.
  //
  // We keep BOTH implementations (for traceability/regression comparisons), but
  // the older/alternate implementation below is renamed with a _v2 suffix so the
  // primary implementation above remains the one used by the code path.

static void WriteTecplotPoints_DipoleAnalyticDensityCompare(const EarthUtil::AmpsParam& prm,
                                                            const std::vector<EarthUtil::Vec3>& points,
                                                            const std::vector<double>& n_num_m3) {
  FILE* f=std::fopen("density_gridless_points_dipole_compare.dat","w");
  if (!f) throw std::runtime_error("Cannot write Tecplot file: density_gridless_points_dipole_compare.dat");

  std::fprintf(f,"TITLE=\"Dipole Density: Numeric vs Analytic Hard-Cutoff\"\n"); 
  std::fprintf(f,"VARIABLES=\"id\",\"x\",\"y\",\"z\",\"n_num_m^-3\",\"n_ana_m^-3\",\"rel_err\"\n"); 
  std::fprintf(f,"ZONE T=\"points\" I=%zu F=POINT\n",points.size()); 

  const double mx = Earth::GridlessMode::Dipole::gParams.m_hat[0];
  const double my = Earth::GridlessMode::Dipole::gParams.m_hat[1];
  const double mz = Earth::GridlessMode::Dipole::gParams.m_hat[2];

  const int nPts = prm.densitySpectrum.nPoints();
  const double Emin = prm.densitySpectrum.Emin_MeV;
  const double Emax = prm.densitySpectrum.Emax_MeV;

  for (size_t ip=0; ip<points.size(); ++ip) {
    const double x_m = points[ip].x*1000.0;
    const double y_m = points[ip].y*1000.0;
    const double z_m = points[ip].z*1000.0;
    const double r_m = std::sqrt(x_m*x_m + y_m*y_m + z_m*z_m);
    const double rhatx = x_m/r_m;
    const double rhaty = y_m/r_m;
    const double rhatz = z_m/r_m;

    const double sinLam = mx*rhatx + my*rhaty + mz*rhatz;
    const double cosLam = std::sqrt(std::max(0.0, 1.0 - sinLam*sinLam));
    const double rRe = r_m/_EARTH__RADIUS_;

    const double Rv_GV = 14.9 * prm.field.dipoleMoment_Me * std::pow(cosLam,4) / (rRe*rRe);

    double n_ana = 0.0;
    double Eprev_J = 0.0;
    double gprev = 0.0;

    for (int i=0;i<nPts;i++) {
      double Ei_MeV;
      const double a = (nPts==1)?0.0:double(i)/(nPts-1);
      if (prm.densitySpectrum.spacing==EarthUtil::DensitySpectrumParam::Spacing::LOG) {
        Ei_MeV = Emin*std::pow(Emax/Emin, a);
      }
      else {
        Ei_MeV = Emin + (Emax-Emin)*a;
      }

      const double E_J = Ei_MeV * 1.0e6 * 1.602176634e-19;

      const double R_GV = RigidityFromEnergy_GV(E_J, prm);
      const double T = (R_GV >= Rv_GV) ? 1.0 : 0.0;

      const double Jb = gSpectrum.GetSpectrum(E_J);
      const double Jloc = T * Jb;

      const double v = SpeedFromEnergy(E_J, prm);
      const double g = 4.0*M_PI * Jloc / v;

      if (i>0) {
        const double dE = E_J - Eprev_J;
        n_ana += 0.5*(g+gprev)*dE;
      }

      Eprev_J = E_J;
      gprev = g;
    }

    const double rel = (n_ana>0.0) ? (n_num_m3[ip]-n_ana)/n_ana : 0.0;

    std::fprintf(f,"%zu %e %e %e %e %e %e\n", ip, points[ip].x,points[ip].y,points[ip].z,n_num_m3[ip], n_ana, rel); 
  }

  std::fclose(f);
}
