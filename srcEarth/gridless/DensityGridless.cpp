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

 * 2.5 Anisotropic boundary spectrum (ANISOTROPIC mode)
 * When DS_BOUNDARY_MODE = ANISOTROPIC the above is replaced by:
 *
 *   T_aniso(E; x0) = (1/N_dirs) * sum_{k=1}^{N_dirs} A_k * f_aniso_k
 *
 * where f_aniso_k = f_PAD(cos_alpha_k) * f_spatial(x_exit_k) is the anisotropy
 * weight for trajectory k, evaluated by EvalAnisotropyFactor (AnisotropicSpectrum.h).
 * The exit state (cos_alpha_k, x_exit_k) is retrieved from TraceAllowedSharedEx.
 *
 * The final flux and density formulas (steps 2.2 and 2.3) are unchanged:
 *   J_loc(E; x0) = T_aniso(E; x0) * J_b_iso(E)
 *   n_tot(x0)    = 4*pi * integral J_loc(E)/v(E) dE
 *
 * -----------------------------------------------------------------------------
 * 5. ComputeT_atEnergy HELPER
 * -----------------------------------------------------------------------------
 * A file-local static function that encapsulates the inner direction-tracing loop
 * for a single energy point at a single observation location. Both the ISOTROPIC
 * and ANISOTROPIC branches, and both the serial and parallel code paths, call the
 * same function.
 *
 * Signature:
 *   static double ComputeT_atEnergy(
 *       const EarthUtil::AmpsParam& prm,
 *       const V3& x0_m,
 *       double Rgv,
 *       const std::vector<V3>& dirs,
 *       double maxTrajTime_s,
 *       bool doAnisotropic,
 *       const EarthUtil::AnisotropyParam& anisoPar);
 *
 * Operation:
 *   double weightSum = 0.0;
 *   for each direction d in dirs:
 *     if doAnisotropic:
 *       TrajectoryExitState exit = {};
 *       allowed = TraceAllowedSharedEx(prm, x0, d, Rgv, &exit, maxTrajTime_s)
 *       if allowed:
 *         weightSum += EvalAnisotropyFactor(anisoPar, exit.cosAlpha, exit.x_exit_m)
 *     else:
 *       allowed = TraceAllowedShared(prm, x0, d, Rgv, maxTrajTime_s)
 *       if allowed:
 *         weightSum += 1.0
 *   return weightSum / dirs.size()
 *
 * The rationale for this extraction is described in DensityGridless.h: prior to
 * this refactor the inner loop appeared four times (serial/parallel x POINTS/SHELLS),
 * making the code fragile. ComputeT_atEnergy is the single authoritative copy.
 *
 * -----------------------------------------------------------------------------
 * 6. MPI SCHEDULING FOR DENSITY/SPECTRUM
 * -----------------------------------------------------------------------------
 * Density mode uses *point-level* dynamic scheduling (one observation point per
 * MPI task), analogous to the cutoff solver but with heavier payloads:
 *
 *   ResultMsg for density:
 *     { int pointIdx, double density, double T[nEnergyPoints] }
 *
 * Each worker:
 *   (a) Receives a TaskMsg containing a point index.
 *   (b) Evaluates ComputeT_atEnergy at all nEnergyPoints energies for that point.
 *   (c) Integrates n_tot using the trapezoidal rule.
 *   (d) Sends { pointIdx, n_tot, T[0..N-1] } back to master rank 0.
 *   (e) Requests the next task.
 *
 * Master rank 0:
 *   (a) Dispatches tasks in point-index order (not spatial order) to minimize
 *       wait time for the first result.
 *   (b) Collects results and stores them in a per-point array indexed by pointIdx.
 *   (c) After all results are in, writes both output files in the original point
 *       order regardless of which worker processed which point.
 *
 * Because the T[] array can be large (nEnergyPoints doubles per point), the result
 * message size is variable. MPI_Send with MPI_DOUBLE handles this transparently.
 *
 * SHELLS mode uses the same protocol with the observation point being a (shell, lon, lat)
 * tuple encoded as a single integer index.
 *
 * -----------------------------------------------------------------------------
 * 7. PROGRESS REPORTING
 * -----------------------------------------------------------------------------
 * Rank 0 prints a progress line to stdout each time it collects a result:
 *   [DensityGridless] Point k/N done.  n_tot = X.XXe+YY m^-3
 *
 * This is the only console output from the density solver; all scientific output
 * goes to Tecplot files.
 *
 * -----------------------------------------------------------------------------
 * 8. KNOWN LIMITATIONS
 * -----------------------------------------------------------------------------
 * (a) Coordinate transform for SHELLS mode: the observation point grid is built in
 *     GSM spherical coordinates (altitude + GSM lon/lat). For T96/T05 field models
 *     this is correct. For DIPOLE field model the tilt of the dipole relative to
 *     GSM Z means that a Cartesian-distance-based inner sphere check in GSM is not
 *     exactly aligned with the dipole loss cone. This is acceptable for the current
 *     validation use cases but should be documented.
 *
 * (b) The trapezoidal rule is first-order accurate in energy spacing. For strongly
 *     curved transmissivity T(E) (near the cutoff edge), the spectrum output can
 *     benefit from using a finer energy grid (increase DS_NINTERVALS). The density
 *     integral is less sensitive because T(E)*J_b(E)/v(E) is smoother than T(E) alone.
 *
 * (c) The anisotropy weighting in T_aniso is not normalised to preserve n_tot from
 *     the isotropic case. Users who need a normalised comparison should post-process
 *     (divide by the isotropic n_tot from a separate run with ISOTROPIC mode).
 */
 
#include "pic.h"
#include "DensityGridless.h"

#include "CutoffRigidityGridless.h" // shared trajectory classifier used directly by density
#include "GridlessParticleMovers.h"   // shared V3 helpers / mover vocabulary used across gridless modules
#include "AnisotropicSpectrum.h"      // EvalAnisotropyFactor for ANISOTROPIC boundary mode

#include "../boundary/spectrum.h"  // ::gSpectrum and spectrum metadata writers

#include <mpi.h>
#include <cmath>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>

#include <algorithm>
#include <iostream>

#ifdef _OPENMP
#include <omp.h>
#endif

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
// OpenMP helper notes
//--------------------------------------------------------------------------------------
// The density solver now uses OpenMP at two levels of embarrassingly parallel work:
//   (1) over directions for a *single* energy point inside ComputeT_atEnergy(), and
//   (2) over the energy grid for a *single* observation point in the POINTS/SHELLS
//       drivers.
//
// However, enabling both levels simultaneously would create nested parallel regions.
// In practice that usually hurts performance because the total number of worker threads
// becomes roughly (threads over energy) x (threads over directions), which leads to
// oversubscription, cache thrash, and much higher OpenMP runtime overhead.
//
// To avoid that, the inner direction loop checks whether the current thread is already
// inside an outer OpenMP parallel region. If yes, the direction loop falls back to its
// serial implementation for that energy point. If not, the direction loop is free to
// use OpenMP itself. This gives the following behavior:
//   - serial/MPI-only build + OpenMP enabled: directions can be threaded
//   - energy-loop threading active: directions stay serial inside each energy task
//   - OpenMP disabled at compile time: everything reduces to ordinary serial loops
//
// The helper below hides the omp_in_parallel() call so the rest of the code can stay
// readable and can still compile cleanly when _OPENMP is not defined.
static inline bool DensityGridless_IsInsideOpenMPParallelRegion() {
#ifdef _OPENMP
  return omp_in_parallel() != 0;
#else
  return false;
#endif
}

//--------------------------------------------------------------------------------------
// Small 3D vector helper note
//--------------------------------------------------------------------------------------
// The density module now intentionally uses the *same* V3 type and vector helper
// functions that are declared in GridlessParticleMovers.h.  This avoids having one
// V3/add/mul/unit implementation in cutoff and another in density, which in turn makes
// it easier to pass trajectory state between modules without any translation layer.

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
  if (epoch.size()<10) { std::ostringstream _exit_msg; _exit_msg << "Invalid epoch string: '"+epoch+"'"; exit(__LINE__,__FILE__,_exit_msg.str().c_str()); };
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

//--------------------------------------------------------------------------------------
// Shared tracing policy
//--------------------------------------------------------------------------------------
// IMPORTANT ARCHITECTURAL CHANGE
// --------------------------------
// The density/spectrum solver no longer owns a private field evaluator, private Boris
// pusher, private adaptive-dt selector, or private TraceAllowed(...) function.  All of
// those pieces now live in the cutoff-rigidity module and are accessed through the
// exported helper Earth::GridlessMode::TraceAllowedShared(...).
//
// Why this is better:
//   1. One executable now contains exactly one authoritative definition of how a
//      backtraced trajectory is advanced and classified.
//   2. Any future fix in the cutoff tracer (mover bug, dt refinement, inner-sphere
//      crossing logic, geometry handling, etc.) automatically benefits density runs.
//   3. The old density-only bug "dt = 100 * ChooseDt(...)" becomes impossible because
//      density no longer has its own step-size logic at all.
//
// The functions below (energy grid generation, transmissivity averaging, density
// integration, Tecplot output) remain density-specific.  Only the actual single-particle
// trajectory classification has been centralized.

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
  if (n<2) exit(__LINE__,__FILE__,"Energy grid requires at least 2 points");

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
// FluxIntegrateTotal  —  omnidirectional integral flux over the full energy grid
//====================================================================================
// Computes:
//   F_tot(x0) = 4π ∫_{Emin}^{Emax} J_loc(E; x0) dE
//             = 4π ∫_{Emin}^{Emax} T(E; x0) · J_b(E) dE          [m^-2 s^-1]
//
// PHYSICS vs DENSITY
//   Number density:      n   = 4π ∫ J_loc(E)/v(E) dE   [m^-3]
//   Omnidirectional flux: F  = 4π ∫ J_loc(E)    dE    [m^-2 s^-1]
//
//   The only difference is the presence/absence of the 1/v(E) factor.
//   For relativistic protons in the GeV range, v ≈ c, so n ≈ F/c.
//   At non-relativistic energies (E << mp c^2 = 938 MeV) the 1/v factor
//   matters significantly: equal fluxes give more density at low energy.
//
// UNITS
//   Input:
//     E_MeV    : energy grid nodes [MeV]
//     T        : transmissivity at each node (dimensionless, 0–1)
//     GetSpectrum(E_J): boundary intensity J_b [m^-2 s^-1 sr^-1 J^-1]
//
//   The integrand is:  4π sr  ×  T(E) × J_b(E) [m^-2 s^-1 sr^-1 J^-1]  ×  dE [J]
//   Result:  [m^-2 s^-1]
//
// ANALYTIC CHECK (power-law J_b = J0*(E/E0)^{-γ}, γ≠1, constant T)
//   F_tot = 4π · T · J0 · E0^γ / (γ−1) · ( Emin^{1-γ} − Emax^{1-γ} )
//   For γ=2: F_tot = 4π · T · J0 · E0² · ( 1/Emin − 1/Emax )
//   This is evaluated analytically in test_density_analytic.py for comparison.
//====================================================================================
static double FluxIntegrateTotal(
    const std::vector<double>& E_MeV,
    const std::vector<double>& T) {

  const int n = (int)E_MeV.size();
  if (n < 2) return 0.0;
  double s = 0.0;
  for (int i = 0; i+1 < n; ++i) {
    // Convert MeV -> J for the boundary-spectrum call.
    const double Ei_J   = E_MeV[i]   * MEV_TO_J;
    const double Ei1_J  = E_MeV[i+1] * MEV_TO_J;
    const double Jloc_i  = T[i]   * ::gSpectrum.GetSpectrum(Ei_J);
    const double Jloc_i1 = T[i+1] * ::gSpectrum.GetSpectrum(Ei1_J);
    const double dE_J    = Ei1_J - Ei_J;
    s += 0.5 * (Jloc_i + Jloc_i1) * dE_J;
  }
  return 4.0 * M_PI * s;  // [m^-2 s^-1]
}

//====================================================================================
// FluxIntegrateChannel  —  omnidirectional integral flux in a user-defined channel
//====================================================================================
// Computes:
//   F_ch(x0) = 4π ∫_{E1}^{E2} T(E; x0) · J_b(E) dE               [m^-2 s^-1]
//
// ALGORITHM
//   The energy grid {E_i} is the solver's internal grid from BuildEnergyGrid_MeV.
//   The channel [E1, E2] may be narrower than, identical to, or wider than the grid.
//
//   Step 1: Find interior grid points:  i_lo ≤ i ≤ i_hi  where E_i ∈ (E1, E2).
//   Step 2: Build an augmented sub-grid that includes E1 and E2 as endpoints.
//           T is linearly interpolated at E1 and E2 if they do not coincide with
//           grid nodes. Linear interpolation is appropriate because T(E) is
//           smooth on scales larger than one grid interval.
//   Step 3: Integrate the augmented sub-grid with the trapezoidal rule.
//
// BOUNDARY CASES
//   - If E1 >= Emax or E2 <= Emin (channel completely outside the grid): F_ch = 0.
//   - If E1 == E2 (zero-width channel): F_ch = 0 (no trapezoid formed).
//   - If only one grid node falls inside [E1, E2]: a single trapezoid from E1 to E2
//     is formed using interpolated T at both endpoints.
//
// LINEAR INTERPOLATION OF T AT CHANNEL BOUNDARY
//   Given two adjacent grid nodes (E_a, T_a) and (E_b, T_b) with E_a < x < E_b:
//     T(x) ≈ T_a + (T_b - T_a) * (x - E_a) / (E_b - E_a)
//   J_b(x) is evaluated directly (not interpolated) at x via gSpectrum.GetSpectrum.
//   This is exact for analytic spectra; for tabulated spectra the same interpolation
//   scheme as gSpectrum uses applies.
//====================================================================================
static double FluxIntegrateChannel(
    const std::vector<double>& E_MeV,  // solver energy grid [MeV]
    const std::vector<double>& T,       // transmissivity at each grid node
    double E1_MeV,                      // channel lower bound [MeV]
    double E2_MeV) {                    // channel upper bound [MeV]

  const int n = (int)E_MeV.size();
  if (n < 2 || E1_MeV >= E2_MeV) return 0.0;

  // Clip channel to the grid range.
  const double Eg_lo = E_MeV.front();
  const double Eg_hi = E_MeV.back();
  const double lo = std::max(E1_MeV, Eg_lo);
  const double hi = std::min(E2_MeV, Eg_hi);
  if (lo >= hi) return 0.0;

  // Linear interpolation of T at an arbitrary energy x, given that x falls
  // inside the grid (lo >= Eg_lo, hi <= Eg_hi is already guaranteed above).
  auto interpT = [&](double x) -> double {
    // Find the first grid index with E_MeV[k] >= x.
    auto it = std::lower_bound(E_MeV.begin(), E_MeV.end(), x);
    if (it == E_MeV.end()) return T.back();
    const int k = (int)(it - E_MeV.begin());
    if (k == 0) return T[0];
    const double Ea = E_MeV[k-1], Eb = E_MeV[k];
    const double Ta = T[k-1],      Tb = T[k];
    if (Eb <= Ea) return Ta;
    return Ta + (Tb - Ta) * (x - Ea) / (Eb - Ea);
  };

  // Collect all grid nodes strictly inside (lo, hi).
  struct Pt { double E_MeV; double T_val; };
  std::vector<Pt> pts;
  pts.push_back({lo, interpT(lo)});
  for (int i = 0; i < n; ++i) {
    if (E_MeV[i] > lo && E_MeV[i] < hi)
      pts.push_back({E_MeV[i], T[i]});
  }
  pts.push_back({hi, interpT(hi)});

  // Trapezoidal integration of J_loc(E) = T(E) * J_b(E) over the sub-grid.
  double s = 0.0;
  for (size_t k = 0; k+1 < pts.size(); ++k) {
    const double Ei_J  = pts[k].E_MeV   * MEV_TO_J;
    const double Ei1_J = pts[k+1].E_MeV * MEV_TO_J;
    const double Jloc_i  = pts[k].T_val   * ::gSpectrum.GetSpectrum(Ei_J);
    const double Jloc_i1 = pts[k+1].T_val * ::gSpectrum.GetSpectrum(Ei1_J);
    s += 0.5 * (Jloc_i + Jloc_i1) * (Ei1_J - Ei_J);
  }
  return 4.0 * M_PI * s;  // [m^-2 s^-1]
}

//====================================================================================
// ComputeT_atEnergy  —  unified transmissivity helper for ISOTROPIC and ANISOTROPIC
//====================================================================================
// Returns the effective transmissivity T(E; x0) at a single energy point for the
// given observation position x0_m and direction set dirs.
//
// ISOTROPIC mode (doAnisotropic=false):
//   T = N_allowed / N_dirs
//   Each allowed trajectory contributes weight 1.
//
// ANISOTROPIC mode (doAnisotropic=true):
//   T = (1/N_dirs) * sum_k [ A_k * f_PAD(cos_alpha_k) * f_spatial(x_exit_k) ]
//   Each allowed trajectory k contributes f_aniso_k from EvalAnisotropyFactor().
//   This requires TraceAllowedSharedEx() to return the exit state per trajectory.
//
// In both modes the rest of the density computation is identical:
//   J_loc(E) = J_b_iso(E) * T(E)
//   n_tot = 4*pi * integral J_loc(E)/v(E) dE
//
// Arguments:
//   prm           : full parsed parameter block
//   x0_m          : observation point [m] in GSM
//   Rgv           : particle rigidity [GV] at this energy
//   dirs          : direction sample set (full grid or subsampled)
//   maxTrajTime_s : per-trajectory time cap (<= 0 means use global cap)
//   doAnisotropic : false -> isotropic branch, true -> anisotropic branch
//   anisoPar      : anisotropy parameters (only used when doAnisotropic=true)
//====================================================================================
static double ComputeT_atEnergy(const EarthUtil::AmpsParam& prm,
                                 const V3& x0_m,
                                 double Rgv,
                                 const std::vector<V3>& dirs,
                                 double maxTrajTime_s,
                                 bool doAnisotropic,
                                 const EarthUtil::AnisotropyParam& anisoPar) {
  const double x0_arr[3] = {x0_m.x, x0_m.y, x0_m.z};
  double weightSum = 0.0;

  // The direction loop is embarrassingly parallel because each backtraced trajectory
  // is completely independent from the others: one launch direction produces exactly
  // one allowed/forbidden classification (and, in ANISOTROPIC mode, one independent
  // anisotropy weight). The only shared quantity is the accumulated scalar weightSum,
  // which is handled with an OpenMP reduction.
  //
  // IMPORTANT: the loop is executed in parallel only when we are *not* already inside
  // another OpenMP region. This prevents nested parallelism when the caller has already
  // parallelized the outer energy loop for the same observation point.
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(dirs, prm, x0_arr, Rgv, maxTrajTime_s, doAnisotropic, anisoPar) reduction(+:weightSum) if(!DensityGridless_IsInsideOpenMPParallelRegion() && (int)dirs.size() > 1) schedule(dynamic)
#endif
  for (int idir = 0; idir < (int)dirs.size(); ++idir) {
    const auto& d = dirs[idir];

    // Backtrace convention: launch opposite to the desired arrival direction.
    const V3 vTry = mul(-1.0, d);
    const double v0_arr[3] = {vTry.x, vTry.y, vTry.z};

    if (doAnisotropic) {
      // ANISOTROPIC branch: call the extended tracer to get exit state.
      Earth::GridlessMode::TrajectoryExitState exitSt;
      if (Earth::GridlessMode::TraceAllowedSharedEx(
              prm, x0_arr, v0_arr, Rgv, &exitSt, maxTrajTime_s)) {
        weightSum += EvalAnisotropyFactor(anisoPar, exitSt.cosAlpha, exitSt.x_exit_m);
      }
    } else {
      // ISOTROPIC branch: boolean classification, weight = 1 for allowed.
      if (Earth::GridlessMode::TraceAllowedShared(
              prm, x0_arr, v0_arr, Rgv, maxTrajTime_s)) {
        weightSum += 1.0;
      }
    }
  }

  return weightSum / static_cast<double>(dirs.size());
}

//====================================================================================
// ComputeWeightBlockAtEnergy  —  partial direction-block helper for fine MPI balance
//====================================================================================
// This helper evaluates only a contiguous sub-range of the sampled arrival directions
// for one fixed observation point and one fixed energy / rigidity.
//
// WHY THIS EXISTS
//   The original MPI work decomposition assigned one entire point, then later one whole
//   (point,energy) pair, to a worker. That is still sometimes too coarse because the
//   actual expensive physics work is the set of independent traced trajectories over the
//   sampled directions. In difficult regions some direction traces can be much slower
//   than others, so a worker that owns a "bad" full energy sample can still become a
//   straggler.
//
//   To make the scheduling even finer without exploding the MPI message count, we group
//   a *small batch* of directions (typically 10–50 trajectories) into one MPI task.
//   Rank 0 later accumulates the returned partial weights from all batches belonging to
//   the same (point,energy) pair and divides by the total number of sampled directions.
//
// MATHEMATICAL EQUIVALENCE
//   This is an exact repartitioning of the same sum already used by ComputeT_atEnergy():
//
//      weightSum(E,x0) = sum_over_dirs w_k
//      T(E,x0)         = weightSum / N_dirs
//
//   We simply compute the sum in chunks:
//
//      weightSum = sum_over_blocks [ sum_over_dirs_in_block w_k ]
//
//   Therefore the final transmissivity is bitwise-identical up to normal floating-point
//   reordering effects that are already unavoidable in parallel reductions.
//====================================================================================
static double ComputeWeightBlockAtEnergy(const EarthUtil::AmpsParam& prm,
                                         const V3& x0_m,
                                         double Rgv,
                                         const std::vector<V3>& dirs,
                                         int idirBegin,
                                         int idirEnd,
                                         double maxTrajTime_s,
                                         bool doAnisotropic,
                                         const EarthUtil::AnisotropyParam& anisoPar) {
  const double x0_arr[3] = {x0_m.x, x0_m.y, x0_m.z};
  double weightSum = 0.0;

  // Each trajectory in the requested sub-range is independent. We intentionally do not
  // open another OpenMP region here because the MPI work unit is already deliberately
  // small (10–50 trajectories). Creating nested thread teams for such small batches
  // would usually add more overhead than benefit and would also complicate reasoning
  // about the intended MPI-vs-OpenMP work split.
  for (int idir = idirBegin; idir < idirEnd; ++idir) {
    const auto& d = dirs[(size_t)idir];
    const V3 vTry = mul(-1.0, d);
    const double v0_arr[3] = {vTry.x, vTry.y, vTry.z};

    if (doAnisotropic) {
      Earth::GridlessMode::TrajectoryExitState exitSt;
      if (Earth::GridlessMode::TraceAllowedSharedEx(
              prm, x0_arr, v0_arr, Rgv, &exitSt, maxTrajTime_s)) {
        weightSum += EvalAnisotropyFactor(anisoPar, exitSt.cosAlpha, exitSt.x_exit_m);
      }
    }
    else {
      if (Earth::GridlessMode::TraceAllowedShared(
              prm, x0_arr, v0_arr, Rgv, maxTrajTime_s)) {
        weightSum += 1.0;
      }
    }
  }

  return weightSum;
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
namespace Earth {
namespace GridlessMode {

void WriteTecplotPoints_DipoleAnalyticCompare(const EarthUtil::AmpsParam& prm,
                                                     const std::vector<EarthUtil::Vec3>& points,
                                                     const std::vector<double>& n_num_m3,
                                                     const std::vector< std::vector<double> >& T_byPoint,
                                                     const std::vector<double>& flux_tot_m2s1,
                                                     const std::vector< std::vector<double> >& flux_ch);

void WriteTecplotShells_DipoleAnalyticCompare(const EarthUtil::AmpsParam& prm,
                                                     double alt_km,
                                                     double res_deg,
                                                     const std::vector<EarthUtil::Vec3>& shellPts_km,
                                                     const std::vector<double>& n_num_m3,
                                                     const std::vector< std::vector<double> >& T_byPoint,
                                                     const std::vector<double>& flux_tot_m2s1,
                                                     const std::vector< std::vector<double> >& flux_ch);

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
    exit(__LINE__,__FILE__,"Internal error: RunDensityAndSpectrum_POINTS called for non-POINTS mode");
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

  // Branch selection: resolved once per run (CLI override already applied by caller).
  const bool doAnisotropic = (EarthUtil::ToUpper(prm.densitySpectrum.boundaryMode) == "ANISOTROPIC");
  if (mpiRank==0) {
    std::cout << "Boundary mode    : " << (doAnisotropic ? "ANISOTROPIC" : "ISOTROPIC") << "\n";
    if (doAnisotropic) {
      std::cout << "  PAD model      : " << prm.anisotropy.padModel
                << " (n=" << prm.anisotropy.padExponent << ")\n";
      std::cout << "  Spatial model  : " << prm.anisotropy.spatialModel;
      if (EarthUtil::ToUpper(prm.anisotropy.spatialModel)=="DAYSIDE_NIGHTSIDE")
        std::cout << " (day=" << prm.anisotropy.daysideFactor
                  << ", night=" << prm.anisotropy.nightsideFactor << ")";
      std::cout << "\n";
    }
    std::cout.flush();
  }

  // Output buffers on master.
  std::vector<double> density_m3(nPoints, 0.0);
  std::vector<double> flux_tot_m2s1(nPoints, 0.0);
  const int nCh = (int)prm.fluxChannels.size();
  // flux_ch[ic][ip] = integral flux in channel ic at point ip  [m^-2 s^-1]
  std::vector< std::vector<double> > flux_ch(nCh, std::vector<double>(nPoints, 0.0));
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

      // Each energy point can be processed independently: rigidity depends only on the
      // current grid energy, while transmissivity T(E) is computed from a self-contained
      // set of backtraced directions for that energy. Therefore the energy loop is a
      // natural OpenMP target.
      //
      // We use schedule(dynamic) because the tracing cost is not uniform across energy:
      // some energies escape quickly, while others may spend much longer near the cutoff
      // penumbra before being classified. Dynamic scheduling improves load balance among
      // threads for such irregular work.
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(T, E_MeV, prm, x0_m, dirsUse, doAnisotropic, nE) if(nE > 1) schedule(dynamic)
#endif
      for (int ie=0; ie<nE; ++ie) {
        const double Ej = E_MeV[ie]*MEV_TO_J;
        const double Rgv = RigidityFromEnergy_GV(Ej, std::abs(prm.species.charge_e)*QE, prm.species.mass_amu*AMU);
        const double maxTraceTime_s = (prm.densitySpectrum.maxTrajTime_s > 0.0)
                                        ? prm.densitySpectrum.maxTrajTime_s
                                        : -1.0;
        // ComputeT_atEnergy handles both ISOTROPIC and ANISOTROPIC branches.
        T[ie] = ComputeT_atEnergy(prm, x0_m, Rgv, dirsUse, maxTraceTime_s,
                                  doAnisotropic, prm.anisotropy);
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
      // ---- Integral flux (total and per channel) ----
      // F_tot = 4pi * int T(E)*Jb(E) dE  [m^-2 s^-1]  (no 1/v weight)
      flux_tot_m2s1[ip] = FluxIntegrateTotal(E_MeV, T);
      for (int ic = 0; ic < nCh; ++ic) {
        flux_ch[ic][ip] = FluxIntegrateChannel(
            E_MeV, T,
            prm.fluxChannels[ic].E1_MeV,
            prm.fluxChannels[ic].E2_MeV);
      }
      T_byPoint[ip] = std::move(T);

      prog.update((long long)(ip + 1), (long long)nPoints);
    }

    prog.finish((long long)nPoints);
  }
  else {
    // Parallel dynamic scheduling.
    //
    // LOAD-BALANCING IMPROVEMENT
    // --------------------------
    // The original MPI implementation used one *point* as the scheduling unit. That is
    // simple, but it can leave noticeable load imbalance because the actual expensive
    // work is not "one point"; it is "one backtraced transmissivity evaluation at one
    // point and one energy". Different energies can have very different tracing costs,
    // especially near transmissivity transitions / cutoff-like behavior, so assigning a
    // whole point to one worker can make that worker busy for much longer than others.
    //
    // Here we switch to an even finer scheduler whose MPI work unit is a *small
    // block of traced trajectories* for one fixed (point,energy) pair.
    //
    // TOTAL KNOWN WORK
    //   For POINTS mode the total number of traced trajectories is known exactly at the
    //   start of the run:
    //
    //      totalTraj = nPoints * nE * nDirs
    //
    //   where nDirs is the number of sampled arrival directions actually used in the
    //   calculation (full grid or user-requested subsample).
    //
    // WHY BLOCK THE DIRECTIONS
    //   Sending one MPI message per single trajectory would maximize balance but would
    //   also flood the code with very small MPI messages. Instead we batch a modest
    //   number of trajectories (here 32, which lies in the requested 10–50 range) into
    //   one work item. This keeps the granularity fine enough to smooth out stragglers
    //   while still amortizing MPI latency over a meaningful amount of physics work.
    //
    // LINEAR TASK SPACE
    //   A task id now corresponds to:
    //
    //      taskId -> (point index, energy index, direction-block index)
    //
    //   The worker computes the partial weight sum over only that direction block and
    //   returns it to rank 0. Rank 0 accumulates block sums until all directions for an
    //   energy are complete, forms T(E)=weightSum/N_dirs, and once all energies of one
    //   point are complete performs the cheap density / flux reductions.
    const int TAG_TASK=100, TAG_RES_META=200, TAG_RES_VAL=201;
    const int kTrajBatch = 32; // Deliberately in the requested 10–50 range.
    const int nDirs = (int)dirsUse.size();
    const int nDirBlocks = std::max(1, (nDirs + kTrajBatch - 1) / kTrajBatch);
    if (mpiRank==0) {
      ProgressBar prog(mpiRank, "[POINTS]");

      const long long totalTasks = (long long)nPoints * (long long)nE * (long long)nDirBlocks;
      const long long totalTraj  = (long long)nPoints * (long long)nE * (long long)nDirs;
      long long nextTask = 0;
      long long doneTraj = 0;

      std::vector<std::vector<double>> partialWeight((size_t)nPoints, std::vector<double>((size_t)nE, 0.0));
      std::vector<std::vector<int>>    dirsDone     ((size_t)nPoints, std::vector<int>((size_t)nE, 0));
      std::vector<int> energiesDone((size_t)nPoints, 0);

      for (int r=1; r<mpiSize; r++) {
        long long taskId = (nextTask < totalTasks) ? nextTask++ : -1;
        MPI_Send(&taskId, 1, MPI_LONG_LONG, r, TAG_TASK, MPI_COMM_WORLD);
      }

      while (doneTraj < totalTraj) {
        // Receive the result in two explicitly typed messages instead of one raw struct.
        //
        // WHY THIS IS SAFER THAN MPI_BYTE + C++ STRUCT
        //   The earlier implementation sent a locally defined C++ struct as a blob of bytes.
        //   That usually works on homogeneous systems, but it relies on compiler-specific
        //   padding / layout rules and therefore is more fragile than necessary.
        //
        //   Here the worker sends:
        //     1) a fixed-length integer metadata packet: [point index, energy index,
        //        first direction, one-past-last direction]
        //     2) one double containing the partial accumulated weight for that block
        //
        //   This keeps the protocol explicit, typed, and independent of struct packing.
        int meta[4];
        double weightSum = 0.0;
        MPI_Status st;
        MPI_Recv(meta, 4, MPI_INT, MPI_ANY_SOURCE, TAG_RES_META, MPI_COMM_WORLD, &st);
        MPI_Recv(&weightSum, 1, MPI_DOUBLE, st.MPI_SOURCE, TAG_RES_VAL, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        const int idx       = meta[0];
        const int ie        = meta[1];
        const int idirBegin = meta[2];
        const int idirEnd   = meta[3];

        const int batchCount = idirEnd - idirBegin;
        doneTraj += (long long)batchCount;

        partialWeight[(size_t)idx][(size_t)ie] += weightSum;
        dirsDone[(size_t)idx][(size_t)ie] += batchCount;

        if (dirsDone[(size_t)idx][(size_t)ie] == nDirs) {
          T_byPoint[(size_t)idx][(size_t)ie] = partialWeight[(size_t)idx][(size_t)ie] / (double)nDirs;
          energiesDone[(size_t)idx]++;

          if (energiesDone[(size_t)idx] == nE) {
            const std::vector<double>& T = T_byPoint[(size_t)idx];

            std::vector<double> integrand(nE,0.0);
            std::vector<double> EjGrid(nE,0.0);
            for (int ie=0; ie<nE; ++ie) {
              const double Ej = E_MeV[ie]*MEV_TO_J;
              EjGrid[ie] = Ej;
              const double v = SpeedFromEnergy(Ej, prm.species.mass_amu*AMU);
              const double Jb = ::gSpectrum.GetSpectrum(Ej); // boundary spectrum per Joule
              const double Jloc = T[ie] * Jb;
              integrand[ie] = (v>0.0) ? (4.0*M_PI*Jloc/v) : 0.0;
            }

            density_m3[(size_t)idx] = Trapz(EjGrid, integrand);
            flux_tot_m2s1[(size_t)idx] = FluxIntegrateTotal(E_MeV, T);
            for (int ic = 0; ic < nCh; ++ic) {
              flux_ch[ic][(size_t)idx] = FluxIntegrateChannel(
                  E_MeV, T,
                  prm.fluxChannels[ic].E1_MeV,
                  prm.fluxChannels[ic].E2_MeV);
            }
          }
        }

        prog.update(doneTraj, totalTraj);

        long long taskId = (nextTask < totalTasks) ? nextTask++ : -1;
        MPI_Send(&taskId, 1, MPI_LONG_LONG, st.MPI_SOURCE, TAG_TASK, MPI_COMM_WORLD);
      }

      prog.finish(totalTraj);
    }
    else {
      while (true) {
        long long taskId = -1;
        MPI_Recv(&taskId, 1, MPI_LONG_LONG, 0, TAG_TASK, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        if (taskId < 0) break;

        const long long perPoint = (long long)nE * (long long)nDirBlocks;
        const int idx = (int)(taskId / perPoint);
        const long long rem1 = taskId - (long long)idx * perPoint;
        const int ie = (int)(rem1 / (long long)nDirBlocks);
        const int iblock = (int)(rem1 - (long long)ie * (long long)nDirBlocks);
        const int idirBegin = iblock * kTrajBatch;
        const int idirEnd = std::min(idirBegin + kTrajBatch, nDirs);

        const EarthUtil::Vec3 pk = prm.output.points[idx];
        const V3 x0_m { pk.x*1000.0, pk.y*1000.0, pk.z*1000.0 };

        const double Ej = E_MeV[ie]*MEV_TO_J;
        const double Rgv = RigidityFromEnergy_GV(Ej, std::abs(prm.species.charge_e)*QE,
                                                 prm.species.mass_amu*AMU);
        const double maxTraceTime_s = (prm.densitySpectrum.maxTrajTime_s > 0.0)
                                        ? prm.densitySpectrum.maxTrajTime_s
                                        : -1.0;

        const double weightSum = ComputeWeightBlockAtEnergy(prm, x0_m, Rgv, dirsUse,
                                                            idirBegin, idirEnd,
                                                            maxTraceTime_s,
                                                            doAnisotropic,
                                                            prm.anisotropy);

        // Send the result back in two typed messages rather than as a raw struct.
        // The metadata packet fully identifies the finished block, and the separate
        // floating-point message carries the accumulated partial weight.
        int meta[4] = {idx, ie, idirBegin, idirEnd};
        MPI_Send(meta, 4, MPI_INT, 0, TAG_RES_META, MPI_COMM_WORLD);
        MPI_Send(&weightSum, 1, MPI_DOUBLE, 0, TAG_RES_VAL, MPI_COMM_WORLD);
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

    //--------------------------------------------------------------------------
    // File 3: integral flux (total and per user-defined channel)
    //--------------------------------------------------------------------------
    // PHYSICS
    //   Total omnidirectional flux [m^-2 s^-1]:
    //     F_tot(x0) = 4π ∫_{Emin}^{Emax} T(E;x0) · J_b(E) dE
    //
    //   Per-channel flux for channel [E1,E2] [m^-2 s^-1]:
    //     F_ch(x0)  = 4π ∫_{E1}^{E2}  T(E;x0) · J_b(E) dE
    //
    //   Analytic reference (power-law J_b = J0*(E/E0)^{-γ}, γ≠1, T constant):
    //     F = 4π · T · J0 · E0^γ / (γ−1) · ( E1^{1−γ} − E2^{1−γ} )
    //     For γ=2: F = 4π · T · J0 · E0² · ( 1/E1 − 1/E2 )
    //     This is the benchmark evaluated in test_density_analytic.py §9.
    //
    // OUTPUT FORMAT
    //   Variables:  X_km  Y_km  Z_km  F_tot_m2s1  [F_NAME_m2s1 ...]
    //   One ZONE per energy channel plus the total.
    //   Separate zones allow Tecplot to toggle channel visibility.
    //--------------------------------------------------------------------------
    {
      std::ofstream out("gridless_points_flux.dat");
      out << "TITLE=\"Gridless omnidirectional integral flux (POINTS)\"\n";

      // Build variable list header dynamically.
      out << "VARIABLES=\"X_km\" \"Y_km\" \"Z_km\" \"F_tot_m2s1\"";
      for (int ic = 0; ic < nCh; ++ic) {
        out << " \"F_" << prm.fluxChannels[ic].name << "_m2s1\"";
      }
      out << "\n";

      // Single zone: all points.
      out << "ZONE T=\"flux\" I=" << nPoints << " F=POINT\n";
      // Auxiliary data: channel definitions for post-processing scripts.
      out << "AUXDATA Emin_MeV=\"" << prm.densitySpectrum.Emin_MeV << "\"\n";
      out << "AUXDATA Emax_MeV=\"" << prm.densitySpectrum.Emax_MeV << "\"\n";
      for (int ic = 0; ic < nCh; ++ic) {
        out << "AUXDATA CH_" << prm.fluxChannels[ic].name
            << "=\"" << prm.fluxChannels[ic].E1_MeV
            << "_" << prm.fluxChannels[ic].E2_MeV << "_MeV\"\n";
      }

      for (int ip = 0; ip < nPoints; ++ip) {
        const auto& p0 = prm.output.points[ip];
        out << p0.x << " " << p0.y << " " << p0.z
            << " " << flux_tot_m2s1[ip];
        for (int ic = 0; ic < nCh; ++ic) {
          out << " " << flux_ch[ic][ip];
        }
        out << "\n";
      }
    }

    //------------------------------------------------------------------------------------------
    // Nightly test mode: write an analytic-vs-numeric density comparison for the DIPOLE field.
    //------------------------------------------------------------------------------------------
#if _PIC_NIGHTLY_TEST_MODE_ == _PIC_MODE_ON_
    if (EarthUtil::ToUpper(prm.field.model)=="DIPOLE") {
      WriteTecplotPoints_DipoleAnalyticCompare(prm, prm.output.points, density_m3, T_byPoint, flux_tot_m2s1, flux_ch);
    }
#endif

  }

  MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
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
    exit(__LINE__,__FILE__,"Internal error: RunDensityAndSpectrum_SHELLS called for non-SHELLS mode");
  }
  if (prm.output.shellAlt_km.empty()) {
    exit(__LINE__,__FILE__,"OUTPUT_MODE=SHELLS requires at least one altitude in SHELL_ALTITUDES");
  }

  // Shared direction and energy grids.
  const int nZenith=24;
  const int nAz=48;
  const std::vector<V3> dirs = BuildDirGrid(nZenith,nAz);
  const std::vector<double> E_MeV = BuildEnergyGrid_MeV(prm);
  const int nE = (int)E_MeV.size();
  const int nIntervals = std::max(0, prm.densitySpectrum.nIntervals);
  if (nE != nIntervals + 1) {
    exit(__LINE__,__FILE__,"Internal inconsistency: energy grid size != DS_NINTERVALS+1");
  }

  int nDirsUse = (int)dirs.size();
  if (prm.densitySpectrum.maxParticlesPerPoint > 0 && nE > 0) {
    nDirsUse = std::max(1, prm.densitySpectrum.maxParticlesPerPoint / nE);
    nDirsUse = std::min(nDirsUse, (int)dirs.size());
  }
  const std::vector<V3> dirsUse = SelectDirectionsDeterministic(dirs, nDirsUse);

  // Branch selection (mirrors POINTS driver).
  const bool doAnisotropic = (EarthUtil::ToUpper(prm.densitySpectrum.boundaryMode) == "ANISOTROPIC");

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
    std::cout << "Boundary mode   : " << (doAnisotropic ? "ANISOTROPIC" : "ISOTROPIC") << "\n";
    if (doAnisotropic) {
      std::cout << "  PAD model     : " << prm.anisotropy.padModel
                << " (n=" << prm.anisotropy.padExponent << ")\n";
      std::cout << "  Spatial model : " << prm.anisotropy.spatialModel;
      if (EarthUtil::ToUpper(prm.anisotropy.spatialModel)=="DAYSIDE_NIGHTSIDE")
        std::cout << " (day=" << prm.anisotropy.daysideFactor
                  << ", night=" << prm.anisotropy.nightsideFactor << ")";
      std::cout << "\n";
    }
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
      exit(__LINE__,__FILE__,oss.str().c_str());
    }

    // Storage on rank 0.
    std::vector<double> nTot_m3;
    std::vector< std::vector<double> > nChan_m3; // [interval][point]
    std::vector<double> flux_tot_m2s1;
    const int nFluxCh = (int)prm.fluxChannels.size();
    std::vector< std::vector<double> > flux_ch;
    std::vector< std::vector<double> > T_byPoint;
    if (mpiRank==0) {
      nTot_m3.assign(nPts, 0.0);
      nChan_m3.assign(nIntervals, std::vector<double>(nPts, 0.0));
      flux_tot_m2s1.assign(nPts, 0.0);
      flux_ch.assign(nFluxCh, std::vector<double>(nPts, 0.0));
      T_byPoint.assign(nPts, std::vector<double>(nE, 0.0));
    }

    // MPI scheduling over shell point indices (same pattern as POINTS).
    if (mpiSize==1) {
      // Serial path for the current shell.  Even though the MPI scheduler below now
      // works at a much finer granularity (direction blocks for a fixed point/energy
      // pair), the serial case should remain simple and avoid any MPI-specific state.
      // We therefore keep the original serial semantics here: loop over shell points on
      // rank 0, and for each point compute the full T(E) array before doing the cheap
      // density / flux reductions.
      int shellIdx = 0;
      for (size_t si = 0; si < prm.output.shellAlt_km.size(); ++si) {
        if (std::fabs(prm.output.shellAlt_km[si] - alt_km) < 0.01) { shellIdx = (int)si; break; }
      }
      std::ostringstream lbl;
      lbl << "[SHELL " << (shellIdx + 1) << "/" << prm.output.shellAlt_km.size()
          << " alt=" << alt_km << "km]";
      ProgressBar prog(mpiRank, lbl.str());

      for (int idx=0; idx<nPts; ++idx) {
        const auto& pk = shellPts_km[idx];
        const V3 x0_m{pk.x*1000.0, pk.y*1000.0, pk.z*1000.0};

        std::vector<double> T(nE,0.0);

        // As in POINTS mode, the shell serial path can still exploit OpenMP across
        // energies because each transmissivity evaluation T(E) is independent for a
        // fixed spatial location.
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(T, E_MeV, prm, x0_m, dirsUse, doAnisotropic, nE) if(nE > 1) schedule(dynamic)
#endif
        for (int ie=0; ie<nE; ++ie) {
          const double Ej = E_MeV[ie]*MEV_TO_J;
          const double Rgv = RigidityFromEnergy_GV(Ej, std::abs(prm.species.charge_e)*QE,
                                                   prm.species.mass_amu*AMU);
          const double maxTraceTime_s = (prm.densitySpectrum.maxTrajTime_s > 0.0)
                                          ? prm.densitySpectrum.maxTrajTime_s
                                          : -1.0;
          T[ie] = ComputeT_atEnergy(prm, x0_m, Rgv, dirsUse, maxTraceTime_s,
                                    doAnisotropic, prm.anisotropy);
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
        nTot_m3[(size_t)idx] = Trapz(EjGrid, g);
        for (int ic=0; ic<nIntervals; ++ic) {
          const double dE = EjGrid[ic+1] - EjGrid[ic];
          nChan_m3[ic][(size_t)idx] = 0.5*(g[ic] + g[ic+1]) * dE;
        }
        flux_tot_m2s1[(size_t)idx] = FluxIntegrateTotal(E_MeV, T);
        for (int icf=0; icf<nFluxCh; ++icf) {
          flux_ch[icf][(size_t)idx] = FluxIntegrateChannel(E_MeV, T,
                                                           prm.fluxChannels[icf].E1_MeV,
                                                           prm.fluxChannels[icf].E2_MeV);
        }
        T_byPoint[(size_t)idx] = std::move(T);

        prog.update((long long)(idx + 1), (long long)nPts);
      }

      prog.finish((long long)nPts);
    }
    else {
      // Parallel dynamic scheduling for SHELLS mode.
      //
      // The work unit mirrors the new POINTS-mode scheduler: one MPI task covers a
      // small batch of traced directions for a fixed (shell point, energy) pair.
      // This provides finer balancing than assigning a whole shell point to one rank,
      // while still avoiding the excessive MPI overhead that one-message-per-trajectory
      // would create.
      const int TAG_TASK=100, TAG_RES_META=200, TAG_RES_VAL=201;
      const int kTrajBatch = 32; // Deliberately in the requested 10–50 range.
      const int nDirs = (int)dirsUse.size();
      const int nDirBlocks = std::max(1, (nDirs + kTrajBatch - 1) / kTrajBatch);

      if (mpiRank==0) {
        int shellIdx = 0;
        for (size_t si = 0; si < prm.output.shellAlt_km.size(); ++si) {
          if (std::fabs(prm.output.shellAlt_km[si] - alt_km) < 0.01) { shellIdx = (int)si; break; }
        }
        std::ostringstream lbl;
        lbl << "[SHELL " << (shellIdx + 1) << "/" << prm.output.shellAlt_km.size()
            << " alt=" << alt_km << "km]";
        ProgressBar prog(mpiRank, lbl.str());

        const long long totalTasks = (long long)nPts * (long long)nE * (long long)nDirBlocks;
        const long long totalTraj  = (long long)nPts * (long long)nE * (long long)nDirs;
        long long nextTask = 0;
        long long doneTraj = 0;

        std::vector<std::vector<double>> partialWeight((size_t)nPts, std::vector<double>((size_t)nE, 0.0));
        std::vector<std::vector<int>>    dirsDone     ((size_t)nPts, std::vector<int>((size_t)nE, 0));
        std::vector<int> energiesDone((size_t)nPts, 0);

        for (int r=1; r<mpiSize; r++) {
          long long taskId = (nextTask < totalTasks) ? nextTask++ : -1;
          MPI_Send(&taskId, 1, MPI_LONG_LONG, r, TAG_TASK, MPI_COMM_WORLD);
        }

        while (doneTraj < totalTraj) {
          // Receive shell-mode worker results in two typed messages.  This avoids
          // sending a raw packed C++ struct over MPI and therefore avoids any hidden
          // dependence on compiler-inserted padding bytes.
          int meta[4];
          double weightSum = 0.0;
          MPI_Status st;
          MPI_Recv(meta, 4, MPI_INT, MPI_ANY_SOURCE, TAG_RES_META, MPI_COMM_WORLD, &st);
          MPI_Recv(&weightSum, 1, MPI_DOUBLE, st.MPI_SOURCE, TAG_RES_VAL, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

          const int idx       = meta[0];
          const int ie        = meta[1];
          const int idirBegin = meta[2];
          const int idirEnd   = meta[3];

          const int batchCount = idirEnd - idirBegin;
          doneTraj += (long long)batchCount;
          partialWeight[(size_t)idx][(size_t)ie] += weightSum;
          dirsDone[(size_t)idx][(size_t)ie] += batchCount;

          if (dirsDone[(size_t)idx][(size_t)ie] == nDirs) {
            T_byPoint[(size_t)idx][(size_t)ie] = partialWeight[(size_t)idx][(size_t)ie] / (double)nDirs;
            energiesDone[(size_t)idx]++;

            if (energiesDone[(size_t)idx] == nE) {
              const std::vector<double>& T = T_byPoint[(size_t)idx];
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
              nTot_m3[(size_t)idx] = Trapz(EjGrid, g);
              for (int ic=0; ic<nIntervals; ++ic) {
                const double dE = EjGrid[ic+1] - EjGrid[ic];
                nChan_m3[ic][(size_t)idx] = 0.5*(g[ic] + g[ic+1]) * dE;
              }
              flux_tot_m2s1[(size_t)idx] = FluxIntegrateTotal(E_MeV, T);
              for (int icf=0; icf<nFluxCh; ++icf) {
                flux_ch[icf][(size_t)idx] = FluxIntegrateChannel(E_MeV, T,
                                                                     prm.fluxChannels[icf].E1_MeV,
                                                                     prm.fluxChannels[icf].E2_MeV);
              }
            }
          }

          prog.update(doneTraj, totalTraj);

          long long taskId = (nextTask < totalTasks) ? nextTask++ : -1;
          MPI_Send(&taskId, 1, MPI_LONG_LONG, st.MPI_SOURCE, TAG_TASK, MPI_COMM_WORLD);
        }

        prog.finish(totalTraj);
      }
      else {
        while (true) {
          long long taskId = -1;
          MPI_Recv(&taskId, 1, MPI_LONG_LONG, 0, TAG_TASK, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
          if (taskId < 0) break;

          const long long perPoint = (long long)nE * (long long)nDirBlocks;
          const int idx = (int)(taskId / perPoint);
          const long long rem1 = taskId - (long long)idx * perPoint;
          const int ie  = (int)(rem1 / (long long)nDirBlocks);
          const int iblock = (int)(rem1 - (long long)ie * (long long)nDirBlocks);
          const int idirBegin = iblock * kTrajBatch;
          const int idirEnd = std::min(idirBegin + kTrajBatch, nDirs);

          const auto& pk = shellPts_km[idx];
          const V3 x0_m{pk.x*1000.0, pk.y*1000.0, pk.z*1000.0};

          const double Ej = E_MeV[ie]*MEV_TO_J;
          const double Rgv = RigidityFromEnergy_GV(Ej, std::abs(prm.species.charge_e)*QE,
                                                   prm.species.mass_amu*AMU);
          const double maxTraceTime_s = (prm.densitySpectrum.maxTrajTime_s > 0.0)
                                          ? prm.densitySpectrum.maxTrajTime_s
                                          : -1.0;
          const double weightSum = ComputeWeightBlockAtEnergy(prm, x0_m, Rgv, dirsUse,
                                                              idirBegin, idirEnd,
                                                              maxTraceTime_s,
                                                              doAnisotropic,
                                                              prm.anisotropy);

          // Send shell-mode results as two typed messages rather than a raw struct.
          int meta[4] = {idx, ie, idirBegin, idirEnd};
          MPI_Send(meta, 4, MPI_INT, 0, TAG_RES_META, MPI_COMM_WORLD);
          MPI_Send(&weightSum, 1, MPI_DOUBLE, 0, TAG_RES_VAL, MPI_COMM_WORLD);
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

#if _PIC_NIGHTLY_TEST_MODE_ == _PIC_MODE_ON_
      if (EarthUtil::ToUpper(prm.field.model)=="DIPOLE") {
        WriteTecplotShells_DipoleAnalyticCompare(prm, alt_km, prm.output.shellRes_deg,
                                                 shellPts_km, nTot_m3, T_byPoint,
                                                 flux_tot_m2s1, flux_ch);
      }
#endif
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
  { std::ostringstream _exit_msg; _exit_msg << "Unsupported OUTPUT_MODE for DENSITY_SPECTRUM: '"+prm.output.mode+"'  (supported: POINTS,SHELLS)"; exit(__LINE__,__FILE__,_exit_msg.str().c_str()); }

  return -1; //the code should not get to this point -- it is here just to make the compiler happy
}



  // NOTE:
  // This file historically accumulated two independent implementations of the
  // dipole analytic density comparison Tecplot writer with the same function name.
  // Unqualified calls then became ambiguous at compile time.
  //
  // We keep BOTH implementations (for traceability/regression comparisons), but
  // the older/alternate implementation below is renamed with a _v2 suffix so the
  // primary implementation above remains the one used by the code path.


//======================================================================================
// Dipole analytic / semi-analytic comparison helpers for density, spectrum, and flux
//======================================================================================
// NOTE:
// The helper block below was previously wrapped in an anonymous namespace.
// That made the later namespace-qualified writer definitions invalid because
// a definition of Earth::GridlessMode::... must appear at global scope or
// inside a namespace that encloses Earth::GridlessMode.  We keep all comments
// and helper code intact, but remove the anonymous namespace wrapper so the
// qualified writer definitions below are legal and link correctly.
//

struct DipoleAnalyticReference {
  double T_geo{0.0};
  double T_weighted{0.0};
  double Rc_vert_GV{0.0};
  std::vector<double> T_ref;
  double density_m3{0.0};
  double flux_tot_m2s1{0.0};
  std::vector<double> flux_ch_m2s1;
};

static double StormerVerticalCutoff_GV(const EarthUtil::AmpsParam& prm, const EarthUtil::Vec3& p_km) {
  const double x_m = p_km.x*1000.0;
  const double y_m = p_km.y*1000.0;
  const double z_m = p_km.z*1000.0;
  const double r_m = std::sqrt(x_m*x_m + y_m*y_m + z_m*z_m);
  if (r_m <= 0.0) return 0.0;

  const double mx = Earth::GridlessMode::Dipole::gParams.m_hat[0];
  const double my = Earth::GridlessMode::Dipole::gParams.m_hat[1];
  const double mz = Earth::GridlessMode::Dipole::gParams.m_hat[2];
  const double sinLam = (mx*x_m + my*y_m + mz*z_m) / r_m;
  const double cosLam = std::sqrt(std::max(0.0, 1.0 - sinLam*sinLam));
  const double rRe = r_m / _EARTH__RADIUS_;
  return 14.9 * prm.field.dipoleMoment_Me * std::pow(cosLam, 4) / (rRe*rRe);
}

static bool RayBoxExit_m(const EarthUtil::DomainBox& box_km,
                         const double x0_m[3], const double u[3],
                         double& tExit_m, double xExit_m[3]) {
  const double bmin[3] = {box_km.xMin*1000.0, box_km.yMin*1000.0, box_km.zMin*1000.0};
  const double bmax[3] = {box_km.xMax*1000.0, box_km.yMax*1000.0, box_km.zMax*1000.0};

  double tBest = 1.0e300;
  for (int d=0; d<3; ++d) {
    if (std::abs(u[d]) < 1.0e-14) continue;
    for (int side=0; side<2; ++side) {
      const double plane = (side==0) ? bmin[d] : bmax[d];
      const double t = (plane - x0_m[d]) / u[d];
      if (t <= 0.0 || t >= tBest) continue;

      const double x[3] = {x0_m[0] + t*u[0], x0_m[1] + t*u[1], x0_m[2] + t*u[2]};
      bool inside = true;
      for (int q=0; q<3; ++q) {
        if (q==d) continue;
        if (x[q] < bmin[q]-1.0e-9 || x[q] > bmax[q]+1.0e-9) { inside = false; break; }
      }
      if (inside) {
        tBest = t;
        xExit_m[0] = x[0]; xExit_m[1] = x[1]; xExit_m[2] = x[2];
      }
    }
  }

  if (tBest >= 1.0e299) return false;
  tExit_m = tBest;
  return true;
}

static bool RayHitsInnerSphereBeforeExit_m(const EarthUtil::DomainBox& box_km,
                                           const double x0_m[3], const double u[3],
                                           double& tHit_m) {
  double xExit_m[3];
  double tExit_m = 0.0;
  if (!RayBoxExit_m(box_km, x0_m, u, tExit_m, xExit_m)) return false;

  const double rInner_m = box_km.rInner * 1000.0;
  const double b = 2.0*(x0_m[0]*u[0] + x0_m[1]*u[1] + x0_m[2]*u[2]);
  const double c = x0_m[0]*x0_m[0] + x0_m[1]*x0_m[1] + x0_m[2]*x0_m[2] - rInner_m*rInner_m;
  const double disc = b*b - 4.0*c;
  if (disc < 0.0) return false;

  const double sq = std::sqrt(std::max(0.0, disc));
  const double t1 = 0.5*(-b - sq);
  const double t2 = 0.5*(-b + sq);

  double tMin = 1.0e300;
  if (t1 > 0.0) tMin = std::min(tMin, t1);
  if (t2 > 0.0) tMin = std::min(tMin, t2);
  if (tMin >= 1.0e299) return false;
  if (tMin < tExit_m) { tHit_m = tMin; return true; }
  return false;
}

static void DipoleBhat(const double x_m[3], double bhat[3]) {
  double B[3];
  Earth::GridlessMode::Dipole::GetB_Tesla(x_m, B);
  const double bmag = std::sqrt(B[0]*B[0] + B[1]*B[1] + B[2]*B[2]);
  if (bmag <= 0.0) {
    bhat[0] = bhat[1] = 0.0; bhat[2] = 1.0;
  }
  else {
    bhat[0] = B[0]/bmag; bhat[1] = B[1]/bmag; bhat[2] = B[2]/bmag;
  }
}

//======================================================================================
// ComputeStraightLineOpenFractions
//======================================================================================
// PURPOSE
//   Construct the *analytic / semi-analytic open-sky factor* used by the built-in
//   DIPOLE nightly benchmark for density, spectrum, and flux.
//
// WHY THIS EXISTS
//   The full numerical solver obtains transmissivity T(E;x0) by integrating the exact
//   Lorentz equation in the chosen magnetic field model and classifying each trial
//   direction as ALLOWED or FORBIDDEN.  That is the physically richer solution, but
//   for regression testing we also need a simpler reference that:
//
//     1) is deterministic,
//     2) is fast,
//     3) can be evaluated inside the C++ code itself,
//     4) is sensitive to geometry mistakes, unit mistakes, and spectrum-folding bugs.
//
//   The benchmark therefore separates the problem into two conceptually distinct pieces:
//
//     (A) a purely geometric / directional "open-sky" factor,
//     (B) an energy-dependent dipole cutoff factor.
//
//   The product of (A) and (B) becomes the analytic reference transmissivity T_ref(E).
//
// PHYSICAL MEANING OF THE TWO OUTPUTS
//   T_geo
//     Fraction of sampled asymptotic directions that can leave the computational box
//     without intersecting the inner loss sphere, *ignoring magnetic bending*.
//     In the isotropic + uniform case, this is simply the fraction of the sky not
//     occulted by the planet/loss sphere.
//
//   T_weighted
//     Same open fraction, but with each open direction weighted by the boundary
//     anisotropy model:
//
//       T_weighted = (1/Ndir) * Σ_open w_k ,
//
//     where
//
//       w_k = f_PAD(cos α_k) · f_spatial(x_exit,k) .
//
//     In isotropic/uniform mode, w_k = 1 and therefore T_weighted = T_geo.
//
// STRAIGHT-LINE GEOMETRY
//   For the isotropic/uniform case we can derive the open fraction analytically.
//   Let r = |x0| be the observer radius and r_in be the inner loss-sphere radius.
//   The blocked directions form a cone subtending the loss sphere.  If ψ is the
//   half-angle of that cone, then
//
//       sin ψ = r_in / r ,
//
//   and the blocked solid angle is
//
//       Ω_blocked = 2π (1 - cos ψ) .
//
//   The open fraction is then
//
//       T_geo = 1 - Ω_blocked / (4π)
//             = 1/2 (1 + cos ψ)
//             = 1/2 (1 + sqrt(1 - (r_in/r)^2)) .
//
//   That is exactly what is used below when the boundary mode is isotropic/uniform.
//
// ANISOTROPIC CASE
//   Once the boundary population is anisotropic, a closed-form solid-angle integral is
//   generally no longer available because the weight depends on:
//
//     • pitch angle α between the asymptotic velocity direction and the local magnetic
//       field direction at the exit point,
//     • possible spatial asymmetry (e.g. dayside/nightside weighting) through x_exit.
//
//   The code therefore evaluates a *semi-analytic* straight-line reference:
//     • it still ignores magnetic bending,
//     • it still uses line-of-sight ray geometry to decide whether a direction is open,
//     • but it computes the anisotropy weight explicitly for each sampled open direction.
//
//   This is still far cheaper and simpler than the full backtracing solution while being
//   sensitive to mistakes in anisotropy implementation and spectrum folding.
//
// IMPORTANT LIMITATION
//   This helper is *not* meant to be the best physical model of access in a dipole.
//   It is only a benchmark reference.  The full numerical solution remains the source
//   of truth for actual science calculations.
//======================================================================================
static void ComputeStraightLineOpenFractions(const EarthUtil::AmpsParam& prm,
                                             const EarthUtil::Vec3& p_km,
                                             double& T_geo,
                                             double& T_weighted) {
  const std::string bmode = EarthUtil::ToUpper(prm.densitySpectrum.boundaryMode);
  const bool doAniso = (bmode == "ANISOTROPIC");
  const std::string padModel = EarthUtil::ToUpper(prm.anisotropy.padModel);
  const std::string spatialModel = EarthUtil::ToUpper(prm.anisotropy.spatialModel);

  if (!doAniso || (padModel == "ISOTROPIC" && spatialModel == "UNIFORM")) {
    const double r_m = std::sqrt(p_km.x*p_km.x + p_km.y*p_km.y + p_km.z*p_km.z) * 1000.0;
    const double rInner_m = prm.domain.rInner * 1000.0;
    if (r_m <= rInner_m) {
      T_geo = 0.0;
      T_weighted = 0.0;
    }
    else {
      const double s = rInner_m / r_m;
      T_geo = 0.5 * (1.0 + std::sqrt(std::max(0.0, 1.0 - s*s)));
      T_weighted = T_geo;
    }
    return;
  }

  const std::vector<V3> dirs = BuildDirGrid(24, 48);
  const double x0_m[3] = {p_km.x*1000.0, p_km.y*1000.0, p_km.z*1000.0};
  double sumW = 0.0;
  int nAllowed = 0;

  for (const auto& d : dirs) {
    const double u[3] = {d.x, d.y, d.z};
    // Each sampled direction u is interpreted as an *asymptotic arrival direction*.
    // In the straight-line benchmark there is no curvature: the trajectory is just the
    // ray x(s)=x0+s*u.  If that ray intersects the inner loss sphere before it exits the
    // outer box, the direction is considered blocked.  Otherwise it is open.
    double tHit = 0.0;
    if (RayHitsInnerSphereBeforeExit_m(prm.domain, x0_m, u, tHit)) continue;

    double tExit = 0.0;
    double xExit[3];
    if (!RayBoxExit_m(prm.domain, x0_m, u, tExit, xExit)) continue;

    ++nAllowed;
    double bhat[3];
    DipoleBhat(xExit, bhat);
    // Pitch-angle cosine used by the anisotropy model:
    //
    //   cos(alpha) = û · b̂_exit ,
    //
    // where û is the asymptotic particle direction at the outer boundary and b̂_exit is
    // the unit magnetic-field direction of the dipole at the exit point.  This is the
    // physically natural quantity for many heliospheric and SEP boundary distributions,
    // because they are organized relative to the magnetic field rather than to the GSM
    // axes themselves.
    const double cosAlpha = u[0]*bhat[0] + u[1]*bhat[1] + u[2]*bhat[2];
    // EvalAnisotropyFactor() combines:
    //   • pitch-angle dependence f_PAD(cos alpha)
    //   • optional spatial asymmetry f_spatial(x_exit)
    // into one multiplicative weight.  By averaging that weight over *open* directions
    // and dividing by the total number of sampled directions, we obtain the effective
    // anisotropy-weighted transmissivity used by the benchmark.
    sumW += EvalAnisotropyFactor(prm.anisotropy, cosAlpha, xExit);
  }

  T_geo = dirs.empty() ? 0.0 : double(nAllowed) / double(dirs.size());
  T_weighted = dirs.empty() ? 0.0 : sumW / double(dirs.size());
}

//======================================================================================
// BuildDipoleAnalyticReference
//======================================================================================
// PURPOSE
//   Build the full nightly-test reference for one observation point in DIPOLE mode.
//
// REFERENCE MODEL
//   The benchmark transmissivity is approximated as
//
//     T_ref(E; x0) = T_open(x0) · H( R(E) - Rc_vert(x0) ) ,
//
//   where
//
//     T_open(x0)   = T_weighted from ComputeStraightLineOpenFractions(),
//     H(...)       = Heaviside step function,
//     R(E)         = particle rigidity corresponding to kinetic energy E,
//     Rc_vert(x0)  = vertical Størmer cutoff rigidity for a dipole.
//
//   In words:
//
//     • the point sees only a fraction of the asymptotic sky;
//     • among those directions, access is additionally suppressed below the vertical
//       Størmer cutoff;
//     • above the cutoff, all open directions are treated as equally accessible in the
//       isotropic case, or anisotropy-weighted in the anisotropic case.
//
//   This is intentionally simpler than the full numerical solution.  It discards
//   penumbra structure and direction-dependent cutoffs.  That simplification is exactly
//   what makes the reference robust for regression testing.
//
// FROM TRANSMISSIVITY TO OBSERVABLES
//   Once T_ref(E) is known, the benchmark uses the *same* physics definitions as the
//   production solver:
//
//     Local differential intensity:
//       J_loc(E; x0) = T_ref(E; x0) · J_b(E)
//
//     Number density:
//       n(x0) = 4π ∫ J_loc(E; x0) / v(E) dE
//
//     Omnidirectional integral flux:
//       F_tot(x0) = 4π ∫ J_loc(E; x0) dE
//
//     Channel flux for user channel [E1,E2]:
//       F_ch(x0) = 4π ∫_{E1}^{E2} J_loc(E; x0) dE
//
//   Therefore the comparison isolates the *transmissivity* model while keeping all
//   spectrum, relativistic, and quadrature machinery consistent between analytic and
//   numerical paths.
//
// IMPLEMENTATION NOTES
//   • E_MeV is the exact same energy grid as used by the solver.
//   • RigidityFromEnergy_GV() converts energy -> rigidity using species charge and mass.
//   • SpeedFromEnergy() provides v(E) for the density conversion.
//   • FluxIntegrateTotal() and FluxIntegrateChannel() are reused so that the analytic
//     and numerical outputs differ only in T(E), not in integration details.
//======================================================================================
static DipoleAnalyticReference BuildDipoleAnalyticReference(const EarthUtil::AmpsParam& prm,
                                                            const EarthUtil::Vec3& p_km,
                                                            const std::vector<double>& E_MeV) {
  DipoleAnalyticReference ref;
  ref.flux_ch_m2s1.assign(prm.fluxChannels.size(), 0.0);
  ComputeStraightLineOpenFractions(prm, p_km, ref.T_geo, ref.T_weighted);
  ref.Rc_vert_GV = StormerVerticalCutoff_GV(prm, p_km);

  ref.T_ref.resize(E_MeV.size(), 0.0);
  std::vector<double> Ej(E_MeV.size(), 0.0), g(E_MeV.size(), 0.0);
  for (size_t i=0; i<E_MeV.size(); ++i) {
    // The reference is evaluated on the *same* discrete energy nodes as the production
    // solver.  This is deliberate: if analytic and numerical results differ, we want the
    // difference to reflect transmissivity physics, not different quadrature grids.
    const double E_J = E_MeV[i] * MEV_TO_J;
    Ej[i] = E_J;
    const double R_GV = RigidityFromEnergy_GV(E_J, std::abs(prm.species.charge_e)*QE, prm.species.mass_amu*AMU);
    // Hard-cutoff approximation:
    //   below Rc_vert -> inaccessible  -> T_ref = 0
    //   above Rc_vert -> accessible    -> T_ref = open-sky factor
    //
    // This is the spectrum-space analogue of the traditional Størmer cutoff benchmark:
    // the detailed penumbra is collapsed into a step function in rigidity.
    const double Tref = (R_GV >= ref.Rc_vert_GV) ? ref.T_weighted : 0.0;
    ref.T_ref[i] = Tref;
    const double v = SpeedFromEnergy(E_J, prm.species.mass_amu*AMU);
    const double Jb = ::gSpectrum.GetSpectrum(E_J);
    const double Jloc = Tref * Jb;
    g[i] = (v > 0.0) ? (4.0*M_PI*Jloc/v) : 0.0;
  }
  ref.density_m3 = Trapz(Ej, g);
  ref.flux_tot_m2s1 = FluxIntegrateTotal(E_MeV, ref.T_ref);
  for (size_t ic=0; ic<prm.fluxChannels.size(); ++ic) {
    ref.flux_ch_m2s1[ic] = FluxIntegrateChannel(E_MeV, ref.T_ref,
                                                prm.fluxChannels[ic].E1_MeV,
                                                prm.fluxChannels[ic].E2_MeV);
  }
  return ref;
}

//======================================================================================
// WriteTecplotPoints_DipoleAnalyticCompare
//======================================================================================
// WHAT IS WRITTEN
//   Three benchmark files are produced for explicit observation points:
//
//     1) density_gridless_points_dipole_compare.dat
//     2) spectrum_gridless_points_dipole_compare.dat
//     3) flux_gridless_points_dipole_compare.dat
//
//   Together they answer three distinct verification questions:
//
//     • Density file:
//         "Does the integrated local phase-space population agree with the dipole
//          benchmark after applying transmissivity and 1/v weighting?"
//
//     • Spectrum file:
//         "Does the numerical T(E) curve produce the expected filtered spectrum
//          J_loc(E)=T(E)J_b(E), and where in energy do deviations appear?"
//
//     • Flux file:
//         "Does the same T(E) fold into the correct total and channel-integrated
//          omnidirectional fluxes?"
//
// WHY SEPARATE FILES
//   Density, spectrum, and flux emphasize different parts of the pipeline.  A bug can
//   affect one and not the others:
//
//     • a rigidity/cutoff bug mainly appears in T(E) and therefore in the spectrum file;
//     • a speed conversion bug mainly appears in density, because density contains 1/v;
//     • a channel-boundary interpolation bug mainly appears in the flux file.
//
// RELATIVE ERROR DEFINITIONS
//   The comparison files use the convention
//
//       rel_err = (numeric - analytic) / analytic ,
//
//   when the analytic quantity is non-zero.  If the analytic reference is exactly zero,
//   the code falls back to storing the raw numerical value.  This avoids division by zero
//   while still flagging spurious numerical leakage into forbidden regions.
//
// UNITS
//   Density file:
//     n_num_m^-3, n_ana_m^-3
//
//   Spectrum file:
//     J_boundary_perMeV, J_local_num_perMeV, J_local_ana_perMeV
//
//   Flux file:
//     Ftot_num_m2s1, Ftot_ana_m2s1, plus one triple per user-defined channel.
//======================================================================================

void WriteTecplotPoints_DipoleAnalyticCompare(const EarthUtil::AmpsParam& prm,
                                                     const std::vector<EarthUtil::Vec3>& points,
                                                     const std::vector<double>& n_num_m3,
                                                     const std::vector< std::vector<double> >& T_byPoint,
                                                     const std::vector<double>& flux_tot_m2s1,
                                                     const std::vector< std::vector<double> >& flux_ch) {
  const std::vector<double> E_MeV = BuildEnergyGrid_MeV(prm);

  {
    FILE* f = std::fopen("density_gridless_points_dipole_compare.dat", "w");
    if (!f) exit(__LINE__,__FILE__,"Cannot write Tecplot file: density_gridless_points_dipole_compare.dat");
    std::fprintf(f,"TITLE=\"Dipole Density (POINTS): Numeric vs Analytic/Semi-analytic\"\n");
    std::fprintf(f,"VARIABLES=\"id\",\"x_km\",\"y_km\",\"z_km\",\"T_geo_ref\",\"T_open_ref\",\"Rc_vert_GV\",\"n_num_m^-3\",\"n_ana_m^-3\",\"rel_err\"\n");
    std::fprintf(f,"ZONE T=\"points\" I=%zu F=POINT\n", points.size());
    for (size_t ip=0; ip<points.size(); ++ip) {
      DipoleAnalyticReference ref = BuildDipoleAnalyticReference(prm, points[ip], E_MeV);
      const double rel = (std::abs(ref.density_m3) > 0.0) ? (n_num_m3[ip] - ref.density_m3)/ref.density_m3 : 0.0;
      std::fprintf(f, "%zu %e %e %e %e %e %e %e %e %e\n", ip,
                   points[ip].x, points[ip].y, points[ip].z,
                   ref.T_geo, ref.T_weighted, ref.Rc_vert_GV,
                   n_num_m3[ip], ref.density_m3, rel);
    }
    std::fclose(f);
  }

  {
    FILE* f = std::fopen("spectrum_gridless_points_dipole_compare.dat", "w");
    if (!f) exit(__LINE__,__FILE__,"Cannot write Tecplot file: spectrum_gridless_points_dipole_compare.dat");
    std::fprintf(f,"TITLE=\"Dipole Spectrum (POINTS): Numeric vs Analytic/Semi-analytic\"\n");
    std::fprintf(f,"VARIABLES=\"E_MeV\",\"T_num\",\"T_ana\",\"J_boundary_perMeV\",\"J_local_num_perMeV\",\"J_local_ana_perMeV\",\"rel_err_T\",\"rel_err_Jloc\"\n");
    for (size_t ip=0; ip<points.size(); ++ip) {
      DipoleAnalyticReference ref = BuildDipoleAnalyticReference(prm, points[ip], E_MeV);
      std::fprintf(f, "ZONE T=\"P%zu\" I=%zu F=POINT\n", ip, E_MeV.size());
      std::fprintf(f, "AUXDATA X_km=\"%g\"\nAUXDATA Y_km=\"%g\"\nAUXDATA Z_km=\"%g\"\n", points[ip].x, points[ip].y, points[ip].z);
      for (size_t ie=0; ie<E_MeV.size(); ++ie) {
        const double Tnum = T_byPoint[ip][ie];
        const double Tana = ref.T_ref[ie];
        const double JbMeV = ::gSpectrum.GetSpectrumPerMeV(E_MeV[ie]);
        const double Jnum = Tnum * JbMeV;
        const double Jana = Tana * JbMeV;
        const double relT = (std::abs(Tana) > 0.0) ? (Tnum - Tana)/Tana : Tnum;
        const double relJ = (std::abs(Jana) > 0.0) ? (Jnum - Jana)/Jana : Jnum;
        std::fprintf(f, "%e %e %e %e %e %e %e %e\n", E_MeV[ie], Tnum, Tana, JbMeV, Jnum, Jana, relT, relJ);
      }
    }
    std::fclose(f);
  }

  {
    FILE* f = std::fopen("flux_gridless_points_dipole_compare.dat", "w");
    if (!f) exit(__LINE__,__FILE__,"Cannot write Tecplot file: flux_gridless_points_dipole_compare.dat");
    std::fprintf(f,"TITLE=\"Dipole Integral Flux (POINTS): Numeric vs Analytic/Semi-analytic\"\n");
    std::fprintf(f,"VARIABLES=\"id\",\"x_km\",\"y_km\",\"z_km\",\"T_geo_ref\",\"T_open_ref\",\"Rc_vert_GV\",\"Ftot_num_m2s1\",\"Ftot_ana_m2s1\",\"Ftot_rel_err\"");
    for (const auto& ch : prm.fluxChannels) {
      std::fprintf(f, " \"F_%s_num_m2s1\" \"F_%s_ana_m2s1\" \"F_%s_rel_err\"", ch.name.c_str(), ch.name.c_str(), ch.name.c_str());
    }
    std::fprintf(f, "\n");
    std::fprintf(f,"ZONE T=\"points\" I=%zu F=POINT\n", points.size());
    for (size_t ip=0; ip<points.size(); ++ip) {
      DipoleAnalyticReference ref = BuildDipoleAnalyticReference(prm, points[ip], E_MeV);
      const double relTot = (std::abs(ref.flux_tot_m2s1) > 0.0) ? (flux_tot_m2s1[ip] - ref.flux_tot_m2s1)/ref.flux_tot_m2s1 : flux_tot_m2s1[ip];
      std::fprintf(f, "%zu %e %e %e %e %e %e %e %e %e", ip,
                   points[ip].x, points[ip].y, points[ip].z,
                   ref.T_geo, ref.T_weighted, ref.Rc_vert_GV,
                   flux_tot_m2s1[ip], ref.flux_tot_m2s1, relTot);
      for (size_t ic=0; ic<prm.fluxChannels.size(); ++ic) {
        const double Fnum = flux_ch[ic][ip];
        const double Fana = ref.flux_ch_m2s1[ic];
        const double rel = (std::abs(Fana) > 0.0) ? (Fnum - Fana)/Fana : Fnum;
        std::fprintf(f, " %e %e %e", Fnum, Fana, rel);
      }
      std::fprintf(f, "\n");
    }
    std::fclose(f);
  }
}

//======================================================================================
// WriteTecplotShells_DipoleAnalyticCompare
//======================================================================================
// SHELL BENCHMARK PHILOSOPHY
//   In SHELLS mode the same physics as POINTS mode is evaluated, but now on a structured
//   longitude/latitude grid at one or more fixed altitudes.  This turns the benchmark from
//   a pointwise regression check into a *map-based* regression check.
//
//   That matters because many implementation mistakes are geometric and only become obvious
//   when viewed spatially:
//
//     • wrong longitude convention,
//     • latitude indexing mistakes,
//     • GSM sign errors,
//     • altitude/radius conversion errors,
//     • north/south asymmetries caused by dipole geometry,
//     • dayside/nightside asymmetry in anisotropic runs.
//
// FILE ORGANIZATION
//   Instead of writing one file per altitude, the benchmark appends one Tecplot ZONE per
//   altitude into a common comparison file.  This matches the style already used by the
//   cutoff nightly benchmark and makes side-by-side loading easier.
//
//   Three files are written:
//
//     • density_gridless_shells_dipole_compare.dat
//     • flux_gridless_shells_dipole_compare.dat
//     • spectrum_gridless_shells_dipole_compare.dat
//
//   The spectrum file is necessarily larger because each shell point becomes its own
//   energy-dependent ZONE.
//
// SPECTRUM PHYSICS IN SHELLS MODE
//   For every shell grid point x0 and every solver energy E_i:
//
//       T_num(E_i;x0)    = full backtracing transmissivity
//       T_ana(E_i;x0)    = dipole benchmark transmissivity
//       J_num(E_i;x0)    = T_num(E_i;x0) · J_b(E_i)
//       J_ana(E_i;x0)    = T_ana(E_i;x0) · J_b(E_i)
//
//   so the comparison remains fully consistent with the POINTS-mode definition.
//
// NOTE ON ANISOTROPY
//   If DS_BOUNDARY_MODE = ANISOTROPIC, the benchmark still remains meaningful:
//   T_open_ref becomes the anisotropy-weighted open-sky fraction from the straight-line
//   reference model, while the numerical solution uses the exact backtraced exit state.
//   The shell maps are therefore especially useful for diagnosing whether the anisotropic
//   weighting introduces the expected large-scale spatial patterns.
//======================================================================================
void WriteTecplotShells_DipoleAnalyticCompare(const EarthUtil::AmpsParam& prm,
                                                     double alt_km,
                                                     double res_deg,
                                                     const std::vector<EarthUtil::Vec3>& shellPts_km,
                                                     const std::vector<double>& n_num_m3,
                                                     const std::vector< std::vector<double> >& T_byPoint,
                                                     const std::vector<double>& flux_tot_m2s1,
                                                     const std::vector< std::vector<double> >& flux_ch) {
  const std::vector<double> E_MeV = BuildEnergyGrid_MeV(prm);
  const int nLon = std::max(1, (int)std::floor(360.0 / res_deg + 0.5));
  const int nLat = std::max(2, (int)std::floor(180.0 / res_deg + 0.5) + 1);

  {
    static bool firstDensityShellWrite = true;
    const bool newFile = firstDensityShellWrite;
    firstDensityShellWrite = false;
    FILE* f = std::fopen("density_gridless_shells_dipole_compare.dat", newFile ? "w" : "a");
    if (!f) exit(__LINE__,__FILE__,"Cannot write Tecplot file: density_gridless_shells_dipole_compare.dat");
    if (newFile) {
      std::fprintf(f,"TITLE=\"Dipole Density (SHELLS): Numeric vs Analytic/Semi-analytic\"\n");
      std::fprintf(f,"VARIABLES=\"lon_deg\",\"lat_deg\",\"x_km\",\"y_km\",\"z_km\",\"T_geo_ref\",\"T_open_ref\",\"Rc_vert_GV\",\"n_num_m^-3\",\"n_ana_m^-3\",\"rel_err\"\n");
    }
    std::fprintf(f,"ZONE T=\"alt_km=%g\" I=%d J=%d F=POINT\n", alt_km, nLon, nLat);
    for (int j=0; j<nLat; ++j) {
      double lat_deg = -90.0 + res_deg*j; if (lat_deg > 90.0) lat_deg = 90.0;
      for (int i=0; i<nLon; ++i) {
        const int k = i + nLon*j;
        const double lon_deg = res_deg*i;
        DipoleAnalyticReference ref = BuildDipoleAnalyticReference(prm, shellPts_km[k], E_MeV);
        const double rel = (std::abs(ref.density_m3) > 0.0) ? (n_num_m3[k] - ref.density_m3)/ref.density_m3 : 0.0;
        std::fprintf(f, "%e %e %e %e %e %e %e %e %e %e %e\n", lon_deg, lat_deg,
                     shellPts_km[k].x, shellPts_km[k].y, shellPts_km[k].z,
                     ref.T_geo, ref.T_weighted, ref.Rc_vert_GV,
                     n_num_m3[k], ref.density_m3, rel);
      }
    }
    std::fclose(f);
  }

  {
    static bool firstFluxShellWrite = true;
    const bool newFile = firstFluxShellWrite;
    firstFluxShellWrite = false;
    FILE* f = std::fopen("flux_gridless_shells_dipole_compare.dat", newFile ? "w" : "a");
    if (!f) exit(__LINE__,__FILE__,"Cannot write Tecplot file: flux_gridless_shells_dipole_compare.dat");
    if (newFile) {
      std::fprintf(f,"TITLE=\"Dipole Integral Flux (SHELLS): Numeric vs Analytic/Semi-analytic\"\n");
      std::fprintf(f,"VARIABLES=\"lon_deg\",\"lat_deg\",\"x_km\",\"y_km\",\"z_km\",\"T_geo_ref\",\"T_open_ref\",\"Rc_vert_GV\",\"Ftot_num_m2s1\",\"Ftot_ana_m2s1\",\"Ftot_rel_err\"");
      for (const auto& ch : prm.fluxChannels) {
        std::fprintf(f, " \"F_%s_num_m2s1\" \"F_%s_ana_m2s1\" \"F_%s_rel_err\"", ch.name.c_str(), ch.name.c_str(), ch.name.c_str());
      }
      std::fprintf(f, "\n");
    }
    std::fprintf(f,"ZONE T=\"alt_km=%g\" I=%d J=%d F=POINT\n", alt_km, nLon, nLat);
    for (int j=0; j<nLat; ++j) {
      double lat_deg = -90.0 + res_deg*j; if (lat_deg > 90.0) lat_deg = 90.0;
      for (int i=0; i<nLon; ++i) {
        const int k = i + nLon*j;
        const double lon_deg = res_deg*i;
        DipoleAnalyticReference ref = BuildDipoleAnalyticReference(prm, shellPts_km[k], E_MeV);
        const double relTot = (std::abs(ref.flux_tot_m2s1) > 0.0) ? (flux_tot_m2s1[k] - ref.flux_tot_m2s1)/ref.flux_tot_m2s1 : flux_tot_m2s1[k];
        std::fprintf(f, "%e %e %e %e %e %e %e %e %e %e %e", lon_deg, lat_deg,
                     shellPts_km[k].x, shellPts_km[k].y, shellPts_km[k].z,
                     ref.T_geo, ref.T_weighted, ref.Rc_vert_GV,
                     flux_tot_m2s1[k], ref.flux_tot_m2s1, relTot);
        for (size_t ic=0; ic<prm.fluxChannels.size(); ++ic) {
          const double Fnum = flux_ch[ic][k];
          const double Fana = ref.flux_ch_m2s1[ic];
          const double rel = (std::abs(Fana) > 0.0) ? (Fnum - Fana)/Fana : Fnum;
          std::fprintf(f, " %e %e %e", Fnum, Fana, rel);
        }
        std::fprintf(f, "\n");
      }
    }
    std::fclose(f);
  }

  {
    static bool firstSpectrumShellWrite = true;
    const bool newFile = firstSpectrumShellWrite;
    firstSpectrumShellWrite = false;
    FILE* f = std::fopen("spectrum_gridless_shells_dipole_compare.dat", newFile ? "w" : "a");
    if (!f) exit(__LINE__,__FILE__,"Cannot write Tecplot file: spectrum_gridless_shells_dipole_compare.dat");
    if (newFile) {
      std::fprintf(f,"TITLE=\"Dipole Spectrum (SHELLS): Numeric vs Analytic/Semi-analytic\"\n");
      std::fprintf(f,"VARIABLES=\"E_MeV\",\"T_num\",\"T_ana\",\"J_boundary_perMeV\",\"J_local_num_perMeV\",\"J_local_ana_perMeV\",\"rel_err_T\",\"rel_err_Jloc\"\n");
    }
    for (int j=0; j<nLat; ++j) {
      double lat_deg = -90.0 + res_deg*j; if (lat_deg > 90.0) lat_deg = 90.0;
      for (int i=0; i<nLon; ++i) {
        const int k = i + nLon*j;
        const double lon_deg = res_deg*i;
        DipoleAnalyticReference ref = BuildDipoleAnalyticReference(prm, shellPts_km[k], E_MeV);
        std::fprintf(f, "ZONE T=\"alt_%g_lon_%g_lat_%g\" I=%zu F=POINT\n", alt_km, lon_deg, lat_deg, E_MeV.size());
        std::fprintf(f, "AUXDATA ALT_km=\"%g\"\nAUXDATA LON_deg=\"%g\"\nAUXDATA LAT_deg=\"%g\"\nAUXDATA X_km=\"%g\"\nAUXDATA Y_km=\"%g\"\nAUXDATA Z_km=\"%g\"\n",
                     alt_km, lon_deg, lat_deg, shellPts_km[k].x, shellPts_km[k].y, shellPts_km[k].z);
        for (size_t ie=0; ie<E_MeV.size(); ++ie) {
          const double Tnum = T_byPoint[k][ie];
          const double Tana = ref.T_ref[ie];
          const double JbMeV = ::gSpectrum.GetSpectrumPerMeV(E_MeV[ie]);
          const double Jnum = Tnum * JbMeV;
          const double Jana = Tana * JbMeV;
          const double relT = (std::abs(Tana) > 0.0) ? (Tnum - Tana)/Tana : Tnum;
          const double relJ = (std::abs(Jana) > 0.0) ? (Jnum - Jana)/Jana : Jnum;
          std::fprintf(f, "%e %e %e %e %e %e %e %e\n", E_MeV[ie], Tnum, Tana, JbMeV, Jnum, Jana, relT, relJ);
        }
      }
    }
    std::fclose(f);
  }
}


} // namespace GridlessMode
} // namespace Earth
