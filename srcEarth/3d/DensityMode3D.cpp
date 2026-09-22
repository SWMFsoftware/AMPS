//======================================================================================
// DensityMode3D.cpp
//======================================================================================
//
// MODE3D MESH-FIELD DENSITY + FLUX FROM BACKWARD TRANSMISSIVITY
//
// This module adds to standalone/SWMF-coupled `-mode 3d` the same *type* of
// energetic-particle density and omnidirectional flux calculation implemented in
// gridless/DensityGridless.cpp:
//
//   1. At an observation location x0 and energy E, launch a deterministic set of
//      backtraced arrival directions.
//   2. Classify each trajectory as ALLOWED if it escapes the outer model box before
//      hitting the inner loss sphere, or FORBIDDEN otherwise.
//   3. Form resolved-only T plus conservative lower/upper values for unresolved
//      numerical outcomes; an all-unresolved nominal value is NaN rather than zero.
//   4. Fold T(E;x0) with the boundary spectrum J_b(E) to obtain
//        local spectrum: J_loc(E;x0) = T(E;x0) J_b(E),
//        flux:          F = 4*pi * int J_loc(E) dE,
//        density:       n = 4*pi * int J_loc(E)/v(E) dE.
//
// The intentional difference from gridless mode is the field evaluator.  Gridless mode
// calls Tsyganenko/dipole routines directly for every trajectory step.  Mode3D calls
// Earth::Mode3D::TraceTrajectoryMesh(), which uses the AMPS AMR mesh field already
// prepared by Mode3D.cpp (standalone snapshots) or Mode3DForwardSWMF.cpp (coupled SWMF
// snapshots).  This makes the density/flux products consistent with the mesh-backed
// cutoff-rigidity and directional-map products written from the same snapshot.
//
// This is NOT the same as 3d_forward/Density3D.cpp.  The 3d_forward module samples a
// forward Monte-Carlo particle population in AMPS cells.  This module computes the
// gridless-style backward-access transmissivity and then integrates a prescribed
// boundary spectrum.
//======================================================================================

#include "DensityMode3D.h"
#include "CutoffRigidityMode3D.h"  // TraceAllowedMesh/Ex and TrajectoryExitState
#include "Mode3D.h"
#include "Mode3DParallel.h"

#include "pic.h"
#include "Earth.h"

#include "../gridless/GridlessParticleMovers.h"  // V3 helpers
#include "../gridless/AnisotropicSpectrum.h"     // EvalAnisotropyFactor
#include "../boundary/spectrum.h"                // ::gSpectrum
#include "../util/amps_param_parser.h"
#include "../util/FluxNumerics.h"                // common units, grids, quadrature, access accounting

#include "constants.h"
#include "constants.PlanetaryData.h"

#ifndef _NO_SPICE_CALLS_
#include "SpiceUsr.h"
#endif

#include <mpi.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#include <algorithm>
#include <atomic>
#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <mutex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

namespace {

//--------------------------------------------------------------------------------------
// Physical constants and small helpers
//--------------------------------------------------------------------------------------
static constexpr double QE       = Earth::FluxNumerics::kElementaryCharge_C;
static constexpr double AMU      = Earth::FluxNumerics::kAtomicMassUnit_kg;
static constexpr double MEV_TO_J = Earth::FluxNumerics::kMeVToJ;

static inline bool InsideOpenMPParallelRegion_() {
#ifdef _OPENMP
  return omp_in_parallel() != 0;
#else
  return false;
#endif
}

// Direct std::thread density workers run outside OpenMP.  Guard nested OpenMP loops
// with this thread-local flag so THREADS mode does not oversubscribe by launching
// OpenMP teams inside every std::thread worker.
static thread_local bool gInsideDirectDensityWorker_ = false;

struct DirectDensityWorkerScope_ {
  DirectDensityWorkerScope_()  { gInsideDirectDensityWorker_ = true; }
  ~DirectDensityWorkerScope_() { gInsideDirectDensityWorker_ = false; }
};

static inline bool InsideDirectDensityWorker_() {
  return gInsideDirectDensityWorker_;
}

using DensityParallelBackend_ = Earth::Mode3D::ParallelBackend;

static const char* DensityParallelBackendName_(DensityParallelBackend_ backend) {
  return Earth::Mode3D::ParallelBackendName(backend);
}

static DensityParallelBackend_ ResolveDensityParallelBackend_(const EarthUtil::AmpsParam& prm) {
  return Earth::Mode3D::ResolveParallelBackend(prm,"Mode3D density/flux");
}

static int ResolveDensityThreadCount_(const EarthUtil::AmpsParam& prm,
                                      DensityParallelBackend_ backend) {
  return Earth::Mode3D::ResolveParallelThreadCount(prm,backend);
}

static void ApplyWideAffinityForDirectDensityThreadsOnce_(DensityParallelBackend_ backend,
                                                          int densityThreadCount) {
  Earth::Mode3D::ApplyWideAffinityForDirectThreadsOnce(backend,densityThreadCount,
                                                       "Mode3D density/flux");
}

// Convert kinetic energy [J] to rigidity [GV].  This is copied intentionally from the
// shared numerical module so identical energies produce identical rigidity bins
// independent of whether the field backend is direct or mesh based.
static double RigidityFromEnergy_GV_(double E_J, double qabs_C, double m0_kg) {
  return Earth::FluxNumerics::RigidityFromEnergyGV(E_J,qabs_C,m0_kg);
}
//--------------------------------------------------------------------------------------
// Output naming
//--------------------------------------------------------------------------------------
static std::string gDensityOutputFileSuffix;

static std::string DensityOutputFileName_(const char* stem) {
  return std::string(stem) + gDensityOutputFileSuffix + ".dat";
}

static std::string FormatEnergyBoundForName_(double E_MeV) {
  const double rounded = std::round(E_MeV);
  if (std::fabs(E_MeV-rounded) < 1.0e-9) return std::to_string((long long)rounded);
  std::ostringstream os;
  os.setf(std::ios::fixed);
  os << std::setprecision(6) << E_MeV;
  std::string s = os.str();
  while (!s.empty() && s.back()=='0') s.pop_back();
  if (!s.empty() && s.back()=='.') s.pop_back();
  std::replace(s.begin(), s.end(), '.', 'p');
  return s;
}

//--------------------------------------------------------------------------------------
// Energy and direction grids
//--------------------------------------------------------------------------------------
static std::vector<double> BuildEnergyGrid_MeV_(const EarthUtil::AmpsParam& prm) {
  const int nLegacy = prm.densitySpectrum.nPoints();
  if (nLegacy < 2) exit(__LINE__,__FILE__,"Mode3D density requires DS_NINTERVALS >= 1");

  const double Emin = prm.densitySpectrum.Emin_MeV;
  const double Emax = prm.densitySpectrum.Emax_MeV;
  const std::string mode=EarthUtil::ToUpper(prm.densitySpectrum.transmissionMode);
  const bool scan=(mode=="SCAN" || mode=="ADAPTIVE");
  try {
    return Earth::FluxNumerics::BuildEnergyGridMeV(
        Emin,Emax,nLegacy,
        prm.densitySpectrum.spacing==EarthUtil::DensitySpectrumParam::Spacing::LOG
          ? Earth::FluxNumerics::EnergySpacing::Log
          : Earth::FluxNumerics::EnergySpacing::Linear,
        scan,prm.densitySpectrum.transmissionScanN,
        prm.densitySpectrum.transmissionMaxN,
        std::fabs(prm.species.charge_e)*QE,prm.species.mass_amu*AMU);
  }
  catch (const std::exception& e) { exit(__LINE__,__FILE__,e.what()); }
  return std::vector<double>();
}

// Keep the gridless density convention (24 zenith x 48 azimuth) for backend-to-backend
// comparability.  The point set is deterministic, so repeated runs and different MPI
// decompositions produce reproducible transmissivity values.
static std::vector<V3> BuildDirGrid_(int nZenith,int nAz) {
  const std::vector<Earth::FluxNumerics::DirectionSample> shared=
      Earth::FluxNumerics::BuildEqualSolidAngleDirections(nZenith,nAz);
  std::vector<V3> dirs;
  dirs.reserve(shared.size());
  for (std::size_t i=0;i<shared.size();++i)
    dirs.push_back(V3{shared[i].x,shared[i].y,shared[i].z});
  return dirs;
}

// Deterministic subsampling used when DS_MAX_PARTICLES limits total work per point.
// Rather than taking the first N directions (which would bias the sky coverage), choose
// approximately uniformly spaced indices across the full direction list.
static std::vector<V3> SelectDirectionsDeterministic_(const std::vector<V3>& all, int nUse) {
  return Earth::FluxNumerics::SelectDeterministic(all,nUse);
}

//--------------------------------------------------------------------------------------
// Location mapping
//--------------------------------------------------------------------------------------
static int ShellNLon_(double res_deg) {
  return std::max(1, (int)std::floor(360.0/res_deg + 0.5));
}
static int ShellNLat_(double res_deg) {
  return std::max(2, (int)std::floor(180.0/res_deg + 0.5) + 1);
}

static V3 ShellLocationGSM_m_(const EarthUtil::AmpsParam& prm,
                              int shellIndex,int iLon,int jLat,
                              int nLon,int nLat,double res_deg) {
  (void)nLon; (void)nLat;
  double lon = res_deg * (double)iLon;
  double lat = -90.0 + res_deg * (double)jLat;
  if (lat > 90.0) lat = 90.0;

  const double alt_km = prm.output.shellAlt_km[(std::size_t)shellIndex];
  const double r_m    = _RADIUS_(_EARTH_) + alt_km*1000.0;
  const double lonRad = lon*M_PI/180.0;
  const double latRad = lat*M_PI/180.0;
  const double cl     = std::cos(latRad);

  const V3 xFixed{r_m*cl*std::cos(lonRad),
                  r_m*cl*std::sin(lonRad),
                  r_m*std::sin(latRad)};

  // DIPOLE validation runs historically use the spherical grid directly in GSM.
  // For external-field snapshots, match the cutoff solver and rotate Earth-fixed
  // shell labels into GSM with SPICE when available.
  if (EarthUtil::ToUpper(prm.field.model) == "DIPOLE") return xFixed;

#ifndef _NO_SPICE_CALLS_
  // SPICE itself and this small rotation cache are shared state.  Protect this
  // location transform so direct std::thread density workers cannot race while
  // updating cachedEpoch/rot or calling SPICE routines.
  static std::mutex spiceRotationMutex;
  std::lock_guard<std::mutex> lock(spiceRotationMutex);

  static std::string cachedEpoch;
  static SpiceDouble rot[3][3];
  if (cachedEpoch != prm.field.epoch) {
    cachedEpoch = prm.field.epoch;
    SpiceDouble et;
    str2et_c(prm.field.epoch.c_str(), &et);
    pxform_c("ITRF93", "GSM", et, rot);
  }
  SpiceDouble xIn[3]  = {xFixed.x, xFixed.y, xFixed.z};
  SpiceDouble xOut[3] = {0.0,0.0,0.0};
  mxv_c(rot, xIn, xOut);
  return V3{xOut[0],xOut[1],xOut[2]};
#else
  return xFixed;
#endif
}

static V3 LocationByIndex_m_(const EarthUtil::AmpsParam& prm,
                             int loc,int nLon,int nLat,double res_deg,int nPtsShell) {
  const std::string mode = EarthUtil::ToUpper(prm.output.mode);
  if (mode == "POINTS" || mode == "TRAJECTORY") {
    const auto& p0 = prm.output.points[(std::size_t)loc];
    return V3{p0.x*1000.0,p0.y*1000.0,p0.z*1000.0};
  }

  const int shellIndex = loc / nPtsShell;
  const int k          = loc - shellIndex*nPtsShell;
  const int iLon       = k % nLon;
  const int jLat       = k / nLon;
  return ShellLocationGSM_m_(prm,shellIndex,iLon,jLat,nLon,nLat,res_deg);
}

//--------------------------------------------------------------------------------------
// Flux integration
//--------------------------------------------------------------------------------------
static double FluxIntegrateInterval_(const std::vector<double>& E_MeV,
                                     const std::vector<double>& T,
                                     double E1_MeV,double E2_MeV) {
  return Earth::FluxNumerics::IntegrateFlux(E_MeV,T,E1_MeV,E2_MeV,
      [](double E_J) { return ::gSpectrum.GetSpectrum(E_J); });
}


//--------------------------------------------------------------------------------------
// Transmission diagnostics derived from T(E)
//--------------------------------------------------------------------------------------
using TransmissionDiagnostics_=Earth::FluxNumerics::TransmissionDiagnostics;

static TransmissionDiagnostics_ ComputeTransmissionDiagnostics_(
    const EarthUtil::AmpsParam& prm,
    const std::vector<double>& E_MeV,
    const std::vector<double>& T) {
  return Earth::FluxNumerics::ComputeTransmissionDiagnostics(
      E_MeV,T,std::fabs(prm.species.charge_e)*QE,prm.species.mass_amu*AMU);
}

//--------------------------------------------------------------------------------------
// Transmissivity at one energy/location using the Mode3D mesh tracer
//--------------------------------------------------------------------------------------
static Earth::FluxNumerics::AccessAccumulator ComputeT_atEnergy_Mode3D_(const EarthUtil::AmpsParam& prm,
                                        const V3& x0_m,
                                        double Rgv,
                                        const std::vector<V3>& dirs,
                                        double maxTrajTime_s,
                                        bool doAnisotropic,
                                        const EarthUtil::AnisotropyParam& anisoPar) {
  Earth::FluxNumerics::AccessAccumulator total;
  if (dirs.empty()) return total;

  const double x0_arr[3]={x0_m.x,x0_m.y,x0_m.z};
  std::vector<double> weights(dirs.size(),0.0);
  std::vector<int> terminations(dirs.size(),
      static_cast<int>(Earth::GridlessMode::TrajectoryTermination::NumericalFailure));
  std::vector<int> retried(dirs.size(),0);

#ifdef _OPENMP
#pragma omp parallel for default(none) shared(dirs,prm,x0_arr,Rgv,maxTrajTime_s,doAnisotropic,anisoPar,weights,terminations,retried) if(!InsideOpenMPParallelRegion_() && !InsideDirectDensityWorker_() && (int)dirs.size() > 1) schedule(dynamic)
#endif
  for (int idir=0; idir<(int)dirs.size(); ++idir) {
    const V3& arrivalDir=dirs[(std::size_t)idir];
    const V3 vTry=mul(-1.0,arrivalDir);
    const double v0_arr[3]={vTry.x,vTry.y,vTry.z};

    auto result=Earth::Mode3D::TraceTrajectoryMesh(
        prm,x0_arr,v0_arr,Rgv,doAnisotropic,maxTrajTime_s);
    if (!result.resolved() && prm.densitySpectrum.retryUnresolved) {
      EarthUtil::AmpsParam retryPrm=prm;
      retryPrm.numerics.dtTrace_s=std::max(1.0e-12,0.5*prm.numerics.dtTrace_s);
      retryPrm.numerics.maxSteps=(prm.numerics.maxSteps<=std::numeric_limits<int>::max()/2)
        ? 2*prm.numerics.maxSteps : std::numeric_limits<int>::max();
      const double retryBase=(maxTrajTime_s>0.0)
          ? maxTrajTime_s
          : ((prm.cutoff.maxTrajTime_s>0.0)
             ? prm.cutoff.maxTrajTime_s : prm.numerics.maxTraceTime_s);
      const double retryTime=2.0*retryBase;
      result=Earth::Mode3D::TraceTrajectoryMesh(
          retryPrm,x0_arr,v0_arr,Rgv,doAnisotropic,retryTime);
      retried[(std::size_t)idir]=1;
    }

    terminations[(std::size_t)idir]=static_cast<int>(result.termination);
    if (!result.resolved() || !result.allowed()) continue;
    weights[(std::size_t)idir]=doAnisotropic
        ? EvalAnisotropyFactor(anisoPar,result.exitState.cosAlpha,result.exitState.x_exit_m)
        : 1.0;
  }

  // Aggregate serially after the parallel trace loop.  This makes every count and
  // termination category deterministic and avoids a custom OpenMP reduction for the
  // structured result.
  for (std::size_t i=0;i<dirs.size();++i) {
    const Earth::GridlessMode::TrajectoryTermination termination=
        static_cast<Earth::GridlessMode::TrajectoryTermination>(terminations[i]);
    total.Record(termination,weights[i],retried[i]!=0);
  }
  return total;
}

struct DensityResultBuffers {
  std::vector<double> density_m3,density_lower_m3,density_upper_m3;
  std::vector<double> flux_total_m2s1,flux_total_lower_m2s1,flux_total_upper_m2s1;
  std::vector<double> T_flat,T_lower_flat,T_upper_flat,unresolved_fraction_flat;
  std::vector<double> flux_ch_flat,flux_ch_lower_flat,flux_ch_upper_flat;
  std::vector<int> sampled_flat,resolved_flat,allowed_flat,retried_flat;
  std::vector<int> termination_flat; // [loc][energy][TrajectoryTermination]
};

static double MaximumAllowedWeight_(const EarthUtil::AmpsParam& prm,bool anisotropic) {
  if (!anisotropic) return 1.0;
  if (EarthUtil::ToUpper(prm.anisotropy.spatialModel)=="DAYSIDE_NIGHTSIDE")
    return std::max(0.0,std::max(prm.anisotropy.daysideFactor,prm.anisotropy.nightsideFactor));
  return 1.0; // all implemented PAD models have maximum one
}

static DensityResultBuffers ComputeAllLocations_(const EarthUtil::AmpsParam& prm,
                                                 int nLoc,
                                                 int nLon,int nLat,double res_deg,int nPtsShell,
                                                 const std::vector<double>& E_MeV,
                                                 const std::vector<V3>& dirsUse) {
  int mpiRank=0, mpiSize=1;
  MPI_Comm_rank(MPI_GLOBAL_COMMUNICATOR,&mpiRank);
  MPI_Comm_size(MPI_GLOBAL_COMMUNICATOR,&mpiSize);

  const int nE  = (int)E_MeV.size();
  const int nCh = (int)prm.fluxChannels.size();

  DensityResultBuffers local;
  local.density_m3.assign((std::size_t)nLoc,0.0);
  local.density_lower_m3.assign((std::size_t)nLoc,0.0);
  local.density_upper_m3.assign((std::size_t)nLoc,0.0);
  local.flux_total_m2s1.assign((std::size_t)nLoc,0.0);
  local.flux_total_lower_m2s1.assign((std::size_t)nLoc,0.0);
  local.flux_total_upper_m2s1.assign((std::size_t)nLoc,0.0);
  local.T_flat.assign((std::size_t)nLoc*(std::size_t)nE,0.0);
  local.T_lower_flat.assign((std::size_t)nLoc*(std::size_t)nE,0.0);
  local.T_upper_flat.assign((std::size_t)nLoc*(std::size_t)nE,0.0);
  local.unresolved_fraction_flat.assign((std::size_t)nLoc*(std::size_t)nE,0.0);
  local.flux_ch_flat.assign((std::size_t)nCh*(std::size_t)nLoc,0.0);
  local.flux_ch_lower_flat.assign((std::size_t)nCh*(std::size_t)nLoc,0.0);
  local.flux_ch_upper_flat.assign((std::size_t)nCh*(std::size_t)nLoc,0.0);
  local.sampled_flat.assign((std::size_t)nLoc*(std::size_t)nE,0);
  local.resolved_flat.assign((std::size_t)nLoc*(std::size_t)nE,0);
  local.allowed_flat.assign((std::size_t)nLoc*(std::size_t)nE,0);
  local.retried_flat.assign((std::size_t)nLoc*(std::size_t)nE,0);
  local.termination_flat.assign((std::size_t)nLoc*(std::size_t)nE*
                                (std::size_t)Earth::FluxNumerics::kTerminationCount,0);

  //====================================================================================
  // MPI location scheduling
  //====================================================================================
  // As in the Mode3D cutoff solver, every MPI rank has access to the compact global
  // B/E arrays.  The inter-rank scheduler therefore operates on global
  // observation-location indices rather than on rank-owned mesh subdomains.
  //
  // DYNAMIC uses an MPI one-sided atomic work queue: the rank/main thread fetches a
  // chunk of global locations and then the intra-rank backend computes that chunk.  MPI
  // is never called from std::thread workers, so MPI_THREAD_MULTIPLE is not required.
  // BLOCK_CYCLIC and STATIC keep deterministic fallback schedules for debugging.
  //====================================================================================
  const Earth::Mode3D::MpiScheduler mpiScheduler =
      Earth::Mode3D::ResolveMpiScheduler(prm,"Mode3D density/flux");

  std::vector<int> rankWorkList;
  if (mpiScheduler == Earth::Mode3D::MpiScheduler::BLOCK_CYCLIC) {
    rankWorkList.reserve((std::size_t)((nLoc + mpiSize - 1) / mpiSize));
    for (int loc=mpiRank; loc<nLoc; loc+=mpiSize) rankWorkList.push_back(loc);
  }
  else if (mpiScheduler == Earth::Mode3D::MpiScheduler::STATIC) {
    const int begin = (int)((static_cast<long long>(nLoc) * mpiRank) / mpiSize);
    const int end   = (int)((static_cast<long long>(nLoc) * (mpiRank+1)) / mpiSize);
    rankWorkList.reserve((std::size_t)std::max(0,end-begin));
    for (int loc=begin; loc<end; ++loc) rankWorkList.push_back(loc);
  }
  const int nLocalStatic = (int)rankWorkList.size();

  const bool doAnisotropic = (EarthUtil::ToUpper(prm.densitySpectrum.boundaryMode) == "ANISOTROPIC");
  const double maximumAllowedWeight=MaximumAllowedWeight_(prm,doAnisotropic);
  const double qabs_C = std::fabs(prm.species.charge_e)*QE;
  const double m0_kg  = prm.species.mass_amu*AMU;

  const DensityParallelBackend_ densityBackend = ResolveDensityParallelBackend_(prm);
  const int densityThreadCount = ResolveDensityThreadCount_(prm,densityBackend);
  const long long mpiDynamicChunk = Earth::Mode3D::ResolveMpiDynamicChunk(
      prm,densityThreadCount,static_cast<long long>(nLoc));

  // If the direct std::thread backend is selected, repair/widen the MPI-rank CPU
  // affinity before any density worker threads are created.  This reproduces the
  // manual `taskset -apc` operation that is needed on systems where the MPI runtime
  // pins each rank to one CPU.  The function is a no-op for OPENMP/SERIAL backends.
  ApplyWideAffinityForDirectDensityThreadsOnce_(densityBackend,densityThreadCount);

#ifdef _OPENMP
  if (densityBackend == DensityParallelBackend_::OPENMP && densityThreadCount > 0) {
    omp_set_num_threads(densityThreadCount);
  }
#endif

  //====================================================================================
  // Global progress-bar bookkeeping
  //====================================================================================
  // The older implementation printed only rank-0 local progress.  That was misleading in
  // MPI runs because rank 0 can finish its slab while other ranks are still tracing.  The
  // new implementation follows the Mode3D cutoff solver:
  //
  //   1. split each rank's local slab into the same number of synchronized batches;
  //   2. after each batch, all ranks participate in MPI_Allreduce;
  //   3. rank 0 renders one global progress line based on the sum over all ranks.
  //
  // The progress counter is intentionally coarse-grained at the location level.  A full
  // location contains many trajectory traces: N_energy * N_direction.  We report both
  // completed locations and an approximate task count where one task is one
  // energy-direction trace.  The task count gives a realistic ETA scale without requiring
  // MPI calls inside the inner trajectory loops or inside OpenMP regions.
  //====================================================================================
  const std::string outputMode = EarthUtil::ToUpper(prm.output.mode);
  const bool isPoints = (outputMode == "POINTS" || outputMode == "TRAJECTORY");
  const bool isShells = (outputMode == "SHELLS");
  const int nShells = isShells ? (int)prm.output.shellAlt_km.size() : 0;

  const long long tasksPerLocation =
      std::max(1LL, (long long)nE * (long long)std::max<std::size_t>(1,dirsUse.size()));
  const long long totalLocationsGlobal = (long long)nLoc;
  const long long totalTasksGlobal     = totalLocationsGlobal * tasksPerLocation;

  long long doneLocationsLocal  = 0;
  long long doneLocationsGlobal = 0;
  long long doneTasksLocal      = 0;
  long long doneTasksGlobal     = 0;

  std::vector<int> locDonePerShellLocal((std::size_t)std::max(nShells,0),0);
  std::vector<int> locDonePerShellGlobal((std::size_t)std::max(nShells,0),0);
  std::vector<int> locTotalPerShellLocal((std::size_t)std::max(nShells,0),0);
  std::vector<int> locTotalPerShellGlobal((std::size_t)std::max(nShells,0),0);

  if (isShells && nShells > 0) {
    if (mpiScheduler == Earth::Mode3D::MpiScheduler::DYNAMIC) {
      // In dynamic mode no rank has a predetermined shell subset.  The denominator is
      // a property of the global shell grid itself: each shell contains nPtsShell
      // locations.  Fill it directly and avoid an unnecessary collective.
      for (int s=0; s<nShells; ++s) locTotalPerShellGlobal[(std::size_t)s] = nPtsShell;
    }
    else {
      for (int localIdx=0; localIdx<nLocalStatic; ++localIdx) {
        const int loc = rankWorkList[(std::size_t)localIdx];
        const int shellIdx = loc / std::max(1,nPtsShell);
        if (shellIdx>=0 && shellIdx<nShells) locTotalPerShellLocal[(std::size_t)shellIdx]++;
      }

      // Rank 0 needs the denominator for every shell in order to print meaningful
      // per-shell progress.  The allreduce is outside the compute loop and is used only
      // for diagnostics; it does not participate in the physical calculation.
      MPI_Allreduce(locTotalPerShellLocal.data(), locTotalPerShellGlobal.data(),
                    nShells, MPI_INT, MPI_SUM, MPI_GLOBAL_COMMUNICATOR);
    }
  }

  auto mode3d_density_now_seconds = []() -> double { return MPI_Wtime(); };
  const double progressStartTime = mode3d_density_now_seconds();
  double progressLastPrintTime = -1.0;

  auto mode3d_density_fmt_hms = [](double s) -> std::string {
    if (s < 0.0) return std::string("--:--:--");
    long long is = (long long)std::llround(s);
    long long hh = is/3600; is-=hh*3600;
    long long mm = is/60;   is-=mm*60;
    long long ss = is;
    char buf[64];
    std::snprintf(buf,sizeof(buf),"%02lld:%02lld:%02lld",hh,mm,ss);
    return std::string(buf);
  };

  auto maybePrintProgress = [&](long long doneLocations,
                                long long doneTasks,
                                const std::vector<int>& shellDoneGlobal,
                                bool forcePrint) {
    if (mpiRank != 0) return;

    const double t = mode3d_density_now_seconds();
    if (!forcePrint) {
      if (progressLastPrintTime < 0.0) progressLastPrintTime = t;
      if (t - progressLastPrintTime < 1.0) return;
    }
    progressLastPrintTime = t;

    const double frac = (totalTasksGlobal > 0)
        ? (double(doneTasks)/double(totalTasksGlobal)) : 1.0;
    const double dt = t - progressStartTime;
    const double rate = (dt > 0.0) ? (double(doneTasks)/dt) : 0.0;
    double eta_s = -1.0;
    if (rate > 0.0 && totalTasksGlobal > doneTasks)
      eta_s = double(totalTasksGlobal-doneTasks)/rate;

    const int barW = 36;
    int filled = (int)std::floor(frac*barW + 0.5);
    if (filled < 0) filled = 0;
    if (filled > barW) filled = barW;

    std::ostringstream line;

    if (isPoints) {
      if (outputMode == "TRAJECTORY") line << "[Mode3D density TRAJECTORY] ";
      else line << "[Mode3D density POINTS] ";
    }
    else {
      line << "[Mode3D density SHELLS " << nShells << " zones";
      if (nShells == 1) {
        line << " alt=" << prm.output.shellAlt_km[0] << "km";
      }
      else if (nShells > 1 && nShells <= 4) {
        line << " alt=";
        for (int s=0; s<nShells; ++s) {
          if (s) line << ",";
          line << prm.output.shellAlt_km[(std::size_t)s];
        }
        line << "km";
      }
      else if (nShells > 4) {
        line << " alt=" << prm.output.shellAlt_km.front()
             << ".." << prm.output.shellAlt_km.back() << "km";
      }
      line << "] ";
    }

    line << "[rank 0/global over " << mpiSize << " MPI ranks] ";
    line << "[";
    for (int i=0; i<barW; ++i) line << (i<filled ? "#" : "-");
    line << "] ";

    line.setf(std::ios::fixed);
    line.precision(1);
    line << (frac*100.0) << "%  ";

    line << "(Loc " << doneLocations << "/" << totalLocationsGlobal
         << ", Task " << doneTasks << "/" << totalTasksGlobal;

    if (isShells && nShells > 0) {
      line << "; Shells ";
      if (nShells <= 4) {
        for (int s=0; s<nShells; ++s) {
          if (s) line << ", ";
          const int shellDone  = shellDoneGlobal[(std::size_t)s];
          const int shellTotal = locTotalPerShellGlobal[(std::size_t)s];
          const double pct = (shellTotal > 0)
              ? 100.0*double(shellDone)/double(shellTotal) : 100.0;
          line << (s+1) << ":" << shellDone << "/" << shellTotal << " " << pct << "%";
        }
      }
      else {
        int nCompleteShells = 0;
        int slowestShell = -1;
        double slowestFrac = 2.0;
        for (int s=0; s<nShells; ++s) {
          const int shellDone  = shellDoneGlobal[(std::size_t)s];
          const int shellTotal = locTotalPerShellGlobal[(std::size_t)s];
          const double shellFrac = (shellTotal > 0) ? double(shellDone)/double(shellTotal) : 1.0;
          if (shellFrac >= 1.0) nCompleteShells++;
          if (shellFrac < slowestFrac) {
            slowestFrac = shellFrac;
            slowestShell = s;
          }
        }
        line << nCompleteShells << "/" << nShells << " complete";
        if (slowestShell >= 0) {
          const int shellDone  = shellDoneGlobal[(std::size_t)slowestShell];
          const int shellTotal = locTotalPerShellGlobal[(std::size_t)slowestShell];
          line << ", slowest " << (slowestShell+1) << ":"
               << shellDone << "/" << shellTotal << " " << (100.0*slowestFrac) << "%";
        }
      }
    }

    line << ")  ETA " << mode3d_density_fmt_hms(eta_s) << "\n";
    std::cout << line.str();
    std::cout.flush();
  };

  auto computeOneLocation = [&](int loc) {
    const V3 x0_m = LocationByIndex_m_(prm,loc,nLon,nLat,res_deg,nPtsShell);
    std::vector<double> T((std::size_t)nE,0.0),TLower((std::size_t)nE,0.0),
                        TUpper((std::size_t)nE,0.0);
    std::vector<Earth::FluxNumerics::AccessAccumulator> blocks((std::size_t)nE);

    // The OpenMP loop is inside a C++ lambda.  With default(none), GCC treats
    // variables captured by the enclosing lambda as closure members and may not
    // accept them in the OpenMP data-sharing clauses by their original names.
    // Use local pointer/value aliases declared inside this lambda body, and
    // reference only those aliases inside the parallel region.
    const int nE_local = nE;
    const std::vector<double>* E_MeV_ptr = &E_MeV;
    const EarthUtil::AmpsParam* prm_ptr = &prm;
    const std::vector<V3>* dirsUse_ptr = &dirsUse;
    const EarthUtil::AnisotropyParam* aniso_ptr = &prm.anisotropy;
    const bool doAnisotropic_local = doAnisotropic;
    const double qabs_C_local = qabs_C;
    const double m0_kg_local = m0_kg;
    const double maxTraceTime_s = (prm.densitySpectrum.maxTrajTime_s > 0.0)
                                ? prm.densitySpectrum.maxTrajTime_s
                                : -1.0;

#ifdef _OPENMP
#pragma omp parallel for default(none) \
  shared(blocks,x0_m,E_MeV_ptr,prm_ptr,dirsUse_ptr,aniso_ptr) \
  firstprivate(nE_local,doAnisotropic_local,qabs_C_local,m0_kg_local,maxTraceTime_s) \
  if(!InsideDirectDensityWorker_() && nE_local > 1) schedule(dynamic)
#endif
    for (int ie=0; ie<nE_local; ++ie) {
      const double Ej  = (*E_MeV_ptr)[(std::size_t)ie]*MEV_TO_J;
      const double Rgv = RigidityFromEnergy_GV_(Ej,qabs_C_local,m0_kg_local);
      blocks[(std::size_t)ie] = ComputeT_atEnergy_Mode3D_(*prm_ptr,x0_m,Rgv,*dirsUse_ptr,
                                                          maxTraceTime_s,doAnisotropic_local,
                                                          *aniso_ptr);
    }

    for (int ie=0; ie<nE; ++ie) {
      const Earth::FluxNumerics::AccessEstimate estimate=
          Earth::FluxNumerics::ResolveAccess(blocks[(std::size_t)ie],maximumAllowedWeight);
      T[(std::size_t)ie]=estimate.nominal;
      TLower[(std::size_t)ie]=estimate.lower;
      TUpper[(std::size_t)ie]=estimate.upper;
      const std::size_t flat=(std::size_t)loc*(std::size_t)nE+(std::size_t)ie;
      local.T_flat[flat]=estimate.nominal;
      local.T_lower_flat[flat]=estimate.lower;
      local.T_upper_flat[flat]=estimate.upper;
      local.unresolved_fraction_flat[flat]=estimate.unresolvedFraction;
      local.sampled_flat[flat]=blocks[(std::size_t)ie].sampled;
      local.resolved_flat[flat]=blocks[(std::size_t)ie].resolved;
      local.allowed_flat[flat]=blocks[(std::size_t)ie].allowed;
      local.retried_flat[flat]=blocks[(std::size_t)ie].retried;
      for (int it=0;it<Earth::FluxNumerics::kTerminationCount;++it)
        local.termination_flat[flat*(std::size_t)Earth::FluxNumerics::kTerminationCount+
                               (std::size_t)it]=blocks[(std::size_t)ie].terminationCounts[(std::size_t)it];
    }

    const auto spectrum=[](double E_J) { return ::gSpectrum.GetSpectrum(E_J); };
    local.density_m3[(std::size_t)loc]=
        Earth::FluxNumerics::IntegrateDensity(E_MeV,T,m0_kg,spectrum);
    local.density_lower_m3[(std::size_t)loc]=
        Earth::FluxNumerics::IntegrateDensity(E_MeV,TLower,m0_kg,spectrum);
    local.density_upper_m3[(std::size_t)loc]=
        Earth::FluxNumerics::IntegrateDensity(E_MeV,TUpper,m0_kg,spectrum);
    local.flux_total_m2s1[(std::size_t)loc] = FluxIntegrateInterval_(E_MeV,T,E_MeV.front(),E_MeV.back());
    local.flux_total_lower_m2s1[(std::size_t)loc] = FluxIntegrateInterval_(E_MeV,TLower,E_MeV.front(),E_MeV.back());
    local.flux_total_upper_m2s1[(std::size_t)loc] = FluxIntegrateInterval_(E_MeV,TUpper,E_MeV.front(),E_MeV.back());
    for (int ic=0; ic<nCh; ++ic) {
      local.flux_ch_flat[(std::size_t)ic*(std::size_t)nLoc + (std::size_t)loc] =
          FluxIntegrateInterval_(E_MeV,T,prm.fluxChannels[(std::size_t)ic].E1_MeV,
                                      prm.fluxChannels[(std::size_t)ic].E2_MeV);
      local.flux_ch_lower_flat[(std::size_t)ic*(std::size_t)nLoc + (std::size_t)loc] =
          FluxIntegrateInterval_(E_MeV,TLower,prm.fluxChannels[(std::size_t)ic].E1_MeV,
                                           prm.fluxChannels[(std::size_t)ic].E2_MeV);
      local.flux_ch_upper_flat[(std::size_t)ic*(std::size_t)nLoc + (std::size_t)loc] =
          FluxIntegrateInterval_(E_MeV,TUpper,prm.fluxChannels[(std::size_t)ic].E1_MeV,
                                           prm.fluxChannels[(std::size_t)ic].E2_MeV);
    }
  };

  auto computeGlobalRange = [&](int begin, int end) {
    if (end <= begin) return;

    if (densityBackend == DensityParallelBackend_::THREADS && densityThreadCount > 1) {
      const int nWork = end - begin;
      const int nWorkers = std::max(1,std::min(densityThreadCount,nWork));
      // Thread-safe dynamic work queue over the contiguous global-location chunk fetched
      // by this MPI rank.  The atomic is local to this rank and assigns each global
      // location to exactly one std::thread worker.  MPI is not called inside workers.
      std::atomic<int> nextLoc(begin);
      std::vector<std::thread> workers;
      workers.reserve((std::size_t)nWorkers);

      for (int iw=0; iw<nWorkers; ++iw) {
        workers.emplace_back([&]() {
          DirectDensityWorkerScope_ workerScope;
          for (;;) {
            const int loc = nextLoc.fetch_add(1,std::memory_order_relaxed);
            if (loc >= end) break;
            computeOneLocation(loc);
          }
        });
      }

      for (std::thread& worker : workers) worker.join();
      return;
    }

    const bool suppressNestedOpenMP =
        (densityBackend == DensityParallelBackend_::SERIAL ||
         densityBackend == DensityParallelBackend_::THREADS);
    if (suppressNestedOpenMP) {
      DirectDensityWorkerScope_ serialScope;
      for (int loc=begin; loc<end; ++loc) computeOneLocation(loc);
      return;
    }

    for (int loc=begin; loc<end; ++loc) computeOneLocation(loc);
  };

  auto computeWorkListRange = [&](int begin, int end) {
    if (end <= begin) return;

    if (densityBackend == DensityParallelBackend_::THREADS && densityThreadCount > 1) {
      const int nWork = end - begin;
      const int nWorkers = std::max(1,std::min(densityThreadCount,nWork));
      // Thread-safe dynamic work queue over this rank's deterministic work list.
      // The MPI work distribution is fixed for STATIC/BLOCK_CYCLIC, but the thread
      // scheduler remains dynamic so one long trajectory does not pin one worker while
      // other workers sit idle.
      std::atomic<int> nextLocalIdx(begin);
      std::vector<std::thread> workers;
      workers.reserve((std::size_t)nWorkers);

      for (int iw=0; iw<nWorkers; ++iw) {
        workers.emplace_back([&]() {
          DirectDensityWorkerScope_ workerScope;
          for (;;) {
            const int localIdx = nextLocalIdx.fetch_add(1,std::memory_order_relaxed);
            if (localIdx >= end) break;
            const int loc = rankWorkList[(std::size_t)localIdx];
            computeOneLocation(loc);
          }
        });
      }

      for (std::thread& worker : workers) worker.join();
      return;
    }

    const bool suppressNestedOpenMP =
        (densityBackend == DensityParallelBackend_::SERIAL ||
         densityBackend == DensityParallelBackend_::THREADS);
    if (suppressNestedOpenMP) {
      DirectDensityWorkerScope_ serialScope;
      for (int localIdx=begin; localIdx<end; ++localIdx) {
        computeOneLocation(rankWorkList[(std::size_t)localIdx]);
      }
      return;
    }

    for (int localIdx=begin; localIdx<end; ++localIdx) {
      computeOneLocation(rankWorkList[(std::size_t)localIdx]);
    }
  };

  auto accountCompletedGlobalRange = [&](int begin, int end) {
    const long long nBatchLocations = (long long)std::max(0,end-begin);
    doneLocationsLocal += nBatchLocations;
    doneTasksLocal     += nBatchLocations * tasksPerLocation;
    if (isShells && nShells > 0) {
      for (int loc=begin; loc<end; ++loc) {
        const int shellIdx = loc / std::max(1,nPtsShell);
        if (shellIdx>=0 && shellIdx<nShells) locDonePerShellLocal[(std::size_t)shellIdx]++;
      }
    }
  };

  auto accountCompletedWorkListRange = [&](int begin, int end) {
    const long long nBatchLocations = (long long)std::max(0,end-begin);
    doneLocationsLocal += nBatchLocations;
    doneTasksLocal     += nBatchLocations * tasksPerLocation;
    if (isShells && nShells > 0) {
      for (int localIdx=begin; localIdx<end; ++localIdx) {
        const int loc = rankWorkList[(std::size_t)localIdx];
        const int shellIdx = loc / std::max(1,nPtsShell);
        if (shellIdx>=0 && shellIdx<nShells) locDonePerShellLocal[(std::size_t)shellIdx]++;
      }
    }
  };

  maybePrintProgress(0,0,locDonePerShellGlobal,true);

  if (mpiScheduler == Earth::Mode3D::MpiScheduler::DYNAMIC) {
    // Two-level dynamic scheduling: each rank dynamically fetches chunks from the MPI
    // RMA counter, and the selected intra-rank backend dynamically computes locations
    // within that chunk.  No progress collectives are placed inside the work loop, so
    // ranks are never forced to wait at a synchronization point between chunks.
    Earth::Mode3D::DynamicMpiLocationScheduler scheduler(
        MPI_GLOBAL_COMMUNICATOR,
        static_cast<long long>(nLoc),
        mpiDynamicChunk,
        "Mode3D density/flux");

    for (;;) {
      const long long chunkStartLL = scheduler.FetchNextChunkStart();
      if (chunkStartLL >= static_cast<long long>(nLoc)) break;

      const long long chunkEndLL = std::min(
          chunkStartLL + scheduler.ChunkSize(),
          static_cast<long long>(nLoc));

      const int chunkStart = static_cast<int>(chunkStartLL);
      const int chunkEnd   = static_cast<int>(chunkEndLL);

      computeGlobalRange(chunkStart,chunkEnd);
      accountCompletedGlobalRange(chunkStart,chunkEnd);
    }

    MPI_Allreduce(&doneLocationsLocal,&doneLocationsGlobal,
                  1,MPI_LONG_LONG,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
    MPI_Allreduce(&doneTasksLocal,&doneTasksGlobal,
                  1,MPI_LONG_LONG,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
    if (isShells && nShells > 0) {
      MPI_Allreduce(locDonePerShellLocal.data(),locDonePerShellGlobal.data(),
                    nShells,MPI_INT,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
    }
    maybePrintProgress(doneLocationsGlobal,doneTasksGlobal,locDonePerShellGlobal,true);
  }
  else {
    int nProgressBatches = std::max(1,std::min(nLoc,200));
    if (densityBackend == DensityParallelBackend_::THREADS && densityThreadCount > 1) {
      nProgressBatches = 1;
    }

    for (int ibatch=0; ibatch<nProgressBatches; ++ibatch) {
      const int localBatchBegin = (int)((static_cast<long long>(nLocalStatic) * ibatch) / nProgressBatches);
      const int localBatchEnd   = (int)((static_cast<long long>(nLocalStatic) * (ibatch+1)) / nProgressBatches);

      if (localBatchEnd > localBatchBegin) {
        computeWorkListRange(localBatchBegin,localBatchEnd);
        accountCompletedWorkListRange(localBatchBegin,localBatchEnd);
      }

      MPI_Allreduce(&doneLocationsLocal,&doneLocationsGlobal,
                    1,MPI_LONG_LONG,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
      MPI_Allreduce(&doneTasksLocal,&doneTasksGlobal,
                    1,MPI_LONG_LONG,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
      if (isShells && nShells > 0) {
        MPI_Allreduce(locDonePerShellLocal.data(),locDonePerShellGlobal.data(),
                      nShells,MPI_INT,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
      }

      maybePrintProgress(doneLocationsGlobal,doneTasksGlobal,locDonePerShellGlobal,
                         ibatch == nProgressBatches-1);
    }
  }

  DensityResultBuffers global;
  global.density_m3.assign((std::size_t)nLoc,0.0);
  global.density_lower_m3.assign((std::size_t)nLoc,0.0);
  global.density_upper_m3.assign((std::size_t)nLoc,0.0);
  global.flux_total_m2s1.assign((std::size_t)nLoc,0.0);
  global.flux_total_lower_m2s1.assign((std::size_t)nLoc,0.0);
  global.flux_total_upper_m2s1.assign((std::size_t)nLoc,0.0);
  global.T_flat.assign((std::size_t)nLoc*(std::size_t)nE,0.0);
  global.T_lower_flat.assign((std::size_t)nLoc*(std::size_t)nE,0.0);
  global.T_upper_flat.assign((std::size_t)nLoc*(std::size_t)nE,0.0);
  global.unresolved_fraction_flat.assign((std::size_t)nLoc*(std::size_t)nE,0.0);
  global.flux_ch_flat.assign((std::size_t)nCh*(std::size_t)nLoc,0.0);
  global.flux_ch_lower_flat.assign((std::size_t)nCh*(std::size_t)nLoc,0.0);
  global.flux_ch_upper_flat.assign((std::size_t)nCh*(std::size_t)nLoc,0.0);
  global.sampled_flat.assign((std::size_t)nLoc*(std::size_t)nE,0);
  global.resolved_flat.assign((std::size_t)nLoc*(std::size_t)nE,0);
  global.allowed_flat.assign((std::size_t)nLoc*(std::size_t)nE,0);
  global.retried_flat.assign((std::size_t)nLoc*(std::size_t)nE,0);
  global.termination_flat.assign((std::size_t)nLoc*(std::size_t)nE*
                                 (std::size_t)Earth::FluxNumerics::kTerminationCount,0);

  MPI_Allreduce(local.density_m3.data(), global.density_m3.data(), nLoc, MPI_DOUBLE, MPI_SUM, MPI_GLOBAL_COMMUNICATOR);
  MPI_Allreduce(local.density_lower_m3.data(), global.density_lower_m3.data(), nLoc, MPI_DOUBLE, MPI_SUM, MPI_GLOBAL_COMMUNICATOR);
  MPI_Allreduce(local.density_upper_m3.data(), global.density_upper_m3.data(), nLoc, MPI_DOUBLE, MPI_SUM, MPI_GLOBAL_COMMUNICATOR);
  MPI_Allreduce(local.flux_total_m2s1.data(), global.flux_total_m2s1.data(), nLoc, MPI_DOUBLE, MPI_SUM, MPI_GLOBAL_COMMUNICATOR);
  MPI_Allreduce(local.flux_total_lower_m2s1.data(), global.flux_total_lower_m2s1.data(), nLoc, MPI_DOUBLE, MPI_SUM, MPI_GLOBAL_COMMUNICATOR);
  MPI_Allreduce(local.flux_total_upper_m2s1.data(), global.flux_total_upper_m2s1.data(), nLoc, MPI_DOUBLE, MPI_SUM, MPI_GLOBAL_COMMUNICATOR);
  MPI_Allreduce(local.T_flat.data(), global.T_flat.data(), nLoc*nE, MPI_DOUBLE, MPI_SUM, MPI_GLOBAL_COMMUNICATOR);
  MPI_Allreduce(local.T_lower_flat.data(), global.T_lower_flat.data(), nLoc*nE, MPI_DOUBLE, MPI_SUM, MPI_GLOBAL_COMMUNICATOR);
  MPI_Allreduce(local.T_upper_flat.data(), global.T_upper_flat.data(), nLoc*nE, MPI_DOUBLE, MPI_SUM, MPI_GLOBAL_COMMUNICATOR);
  MPI_Allreduce(local.unresolved_fraction_flat.data(), global.unresolved_fraction_flat.data(), nLoc*nE, MPI_DOUBLE, MPI_SUM, MPI_GLOBAL_COMMUNICATOR);
  MPI_Allreduce(local.sampled_flat.data(),global.sampled_flat.data(),nLoc*nE,MPI_INT,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
  MPI_Allreduce(local.resolved_flat.data(),global.resolved_flat.data(),nLoc*nE,MPI_INT,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
  MPI_Allreduce(local.allowed_flat.data(),global.allowed_flat.data(),nLoc*nE,MPI_INT,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
  MPI_Allreduce(local.retried_flat.data(),global.retried_flat.data(),nLoc*nE,MPI_INT,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
  MPI_Allreduce(local.termination_flat.data(),global.termination_flat.data(),
                nLoc*nE*Earth::FluxNumerics::kTerminationCount,MPI_INT,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
  if (nCh*nLoc > 0) {
    MPI_Allreduce(local.flux_ch_flat.data(), global.flux_ch_flat.data(), nCh*nLoc, MPI_DOUBLE, MPI_SUM, MPI_GLOBAL_COMMUNICATOR);
    MPI_Allreduce(local.flux_ch_lower_flat.data(), global.flux_ch_lower_flat.data(), nCh*nLoc, MPI_DOUBLE, MPI_SUM, MPI_GLOBAL_COMMUNICATOR);
    MPI_Allreduce(local.flux_ch_upper_flat.data(), global.flux_ch_upper_flat.data(), nCh*nLoc, MPI_DOUBLE, MPI_SUM, MPI_GLOBAL_COMMUNICATOR);
  }

  return global;
}

//--------------------------------------------------------------------------------------
// Tecplot writers
//--------------------------------------------------------------------------------------
static void WritePointOutputs_(const EarthUtil::AmpsParam& prm,
                               const std::vector<double>& E_MeV,
                               const DensityResultBuffers& res) {
  const int nLoc = (int)prm.output.points.size();
  const int nE   = (int)E_MeV.size();
  const int nCh  = (int)prm.fluxChannels.size();

  {
    std::ofstream out(DensityOutputFileName_("mode3d_points_density"));
    out << "TITLE=\"Mode3D mesh-field energetic particle density\"\n";
    out << "VARIABLES=\"X_km\" \"Y_km\" \"Z_km\" \"N_m^-3\" \"N_lower_m^-3\" \"N_upper_m^-3\" "
        << "\"N_cm^-3\" \"N_lower_cm^-3\" \"N_upper_cm^-3\" "
        << "\"Rc_lower_GV\" \"Rc_effective_GV\" \"Rc_upper_GV\" \"PenumbraWidth_GV\" \"T_high\"\n";
    out << "ZONE T=\"density\" I=" << nLoc << " F=POINT\n";
    for (int i=0;i<nLoc;i++) {
      const auto& p0 = prm.output.points[(std::size_t)i];
      const double n_m3  = res.density_m3[(std::size_t)i];
      const double n_cm3 = n_m3*1.0e-6;
      std::vector<double> Tloc((std::size_t)nE,0.0);
      for (int ie=0; ie<nE; ++ie)
        Tloc[(std::size_t)ie] = res.T_flat[(std::size_t)i*(std::size_t)nE + (std::size_t)ie];
      const TransmissionDiagnostics_ td = ComputeTransmissionDiagnostics_(prm,E_MeV,Tloc);
      out << p0.x << " " << p0.y << " " << p0.z << " " << n_m3
          << " " << res.density_lower_m3[(std::size_t)i]
          << " " << res.density_upper_m3[(std::size_t)i] << " " << n_cm3
          << " " << res.density_lower_m3[(std::size_t)i]*1.0e-6
          << " " << res.density_upper_m3[(std::size_t)i]*1.0e-6
          << " " << td.RcLower_GV << " " << td.RcEffective_GV << " " << td.RcUpper_GV
          << " " << td.PenumbraWidth_GV << " " << td.THigh << "\n";
    }
  }

  {
    std::ofstream out(DensityOutputFileName_("mode3d_points_spectrum"));
    out << "TITLE=\"Mode3D mesh-field local energetic particle spectrum\"\n";
    out << "VARIABLES=\"E_MeV\" \"T\" \"T_lower\" \"T_upper\" "
        << "\"unresolved_fraction\" \"N_sampled\" \"N_resolved\" \"N_allowed\" "
        << "\"J_boundary_perMeV\" \"J_local_perMeV\" \"J_local_lower_perMeV\" \"J_local_upper_perMeV\"\n";
    for (int loc=0; loc<nLoc; ++loc) {
      out << "ZONE T=\"loc_" << std::setw(6) << std::setfill('0') << loc << std::setfill(' ')
          << "\" I=" << nE << " F=POINT\n";
      for (int ie=0; ie<nE; ++ie) {
        const double T = res.T_flat[(std::size_t)loc*(std::size_t)nE + (std::size_t)ie];
        const std::size_t flat=(std::size_t)loc*(std::size_t)nE+(std::size_t)ie;
        const double lower=res.T_lower_flat[flat],upper=res.T_upper_flat[flat];
        const double Jb_perMeV = ::gSpectrum.GetSpectrumPerMeV(E_MeV[(std::size_t)ie]);
        out << E_MeV[(std::size_t)ie] << " " << T << " " << lower << " " << upper
            << " " << res.unresolved_fraction_flat[flat] << " " << res.sampled_flat[flat]
            << " " << res.resolved_flat[flat] << " " << res.allowed_flat[flat]
            << " " << Jb_perMeV << " " << T*Jb_perMeV
            << " " << lower*Jb_perMeV << " " << upper*Jb_perMeV << "\n";
      }
    }
  }

  {
    std::ofstream out(DensityOutputFileName_("mode3d_points_flux"));
    out << "TITLE=\"Mode3D mesh-field omnidirectional integral flux\"\n";
    out << "VARIABLES=\"X_km\" \"Y_km\" \"Z_km\" \"F_tot_m2s1\" \"F_tot_lower_m2s1\" \"F_tot_upper_m2s1\"";
    for (int ic=0; ic<nCh; ++ic) out << " \"F_" << prm.fluxChannels[(std::size_t)ic].name << "_m2s1\""
        << " \"F_" << prm.fluxChannels[(std::size_t)ic].name << "_lower_m2s1\""
        << " \"F_" << prm.fluxChannels[(std::size_t)ic].name << "_upper_m2s1\"";
    out << "\n";
    out << "ZONE T=\"flux\" I=" << nLoc << " F=POINT\n";
    for (int i=0;i<nLoc;i++) {
      const auto& p0 = prm.output.points[(std::size_t)i];
      out << p0.x << " " << p0.y << " " << p0.z << " " << res.flux_total_m2s1[(std::size_t)i]
          << " " << res.flux_total_lower_m2s1[(std::size_t)i]
          << " " << res.flux_total_upper_m2s1[(std::size_t)i];
      for (int ic=0; ic<nCh; ++ic) {
        const std::size_t flat=(std::size_t)ic*(std::size_t)nLoc+(std::size_t)i;
        out << " " << res.flux_ch_flat[flat] << " " << res.flux_ch_lower_flat[flat]
            << " " << res.flux_ch_upper_flat[flat];
      }
      out << "\n";
    }
  }
}

static void WriteShellOutputs_(const EarthUtil::AmpsParam& prm,
                               int nLon,int nLat,double res_deg,int nPtsShell,
                               const std::vector<double>& E_MeV,
                               const DensityResultBuffers& res) {
  const int nE = (int)E_MeV.size();
  const int nShells = (int)prm.output.shellAlt_km.size();
  const int nCh     = (int)prm.fluxChannels.size();
  const int nLoc    = nShells*nPtsShell;

  for (int s=0; s<nShells; ++s) {
    const std::string altLabel = FormatEnergyBoundForName_(prm.output.shellAlt_km[(std::size_t)s]);
    const std::string stem = "mode3d_shell_" + altLabel + "km_density_flux";
    std::ofstream out(DensityOutputFileName_(stem.c_str()));

    out << "TITLE=\"Mode3D mesh-field density and flux shell alt="
        << prm.output.shellAlt_km[(std::size_t)s] << " km\"\n";
    out << "VARIABLES=\"Lon_deg\" \"Lat_deg\" \"N_m^-3\" \"N_lower_m^-3\" \"N_upper_m^-3\" "
        << "\"N_cm^-3\" \"F_tot_m2s1\" \"F_tot_lower_m2s1\" \"F_tot_upper_m2s1\" "
        << "\"Rc_lower_GV\" \"Rc_effective_GV\" \"Rc_upper_GV\" \"PenumbraWidth_GV\" \"T_high\"";
    for (int ic=0; ic<nCh; ++ic) out << " \"F_" << prm.fluxChannels[(std::size_t)ic].name << "_m2s1\""
        << " \"F_" << prm.fluxChannels[(std::size_t)ic].name << "_lower_m2s1\""
        << " \"F_" << prm.fluxChannels[(std::size_t)ic].name << "_upper_m2s1\"";
    out << "\n";
    out << "ZONE T=\"shell\" I=" << nLon << " J=" << nLat << " F=POINT\n";

    for (int j=0;j<nLat;j++) {
      double lat = -90.0 + res_deg*(double)j;
      if (lat > 90.0) lat = 90.0;
      for (int i=0;i<nLon;i++) {
        const double lon = res_deg*(double)i;
        const int loc = s*nPtsShell + j*nLon + i;
        const double n_m3 = res.density_m3[(std::size_t)loc];
        std::vector<double> Tloc((std::size_t)nE,0.0);
        for (int ie=0; ie<nE; ++ie)
          Tloc[(std::size_t)ie] = res.T_flat[(std::size_t)loc*(std::size_t)nE + (std::size_t)ie];
        const TransmissionDiagnostics_ td = ComputeTransmissionDiagnostics_(prm,E_MeV,Tloc);
        out << lon << " " << lat << " " << n_m3
            << " " << res.density_lower_m3[(std::size_t)loc]
            << " " << res.density_upper_m3[(std::size_t)loc] << " " << n_m3*1.0e-6
            << " " << res.flux_total_m2s1[(std::size_t)loc]
            << " " << res.flux_total_lower_m2s1[(std::size_t)loc]
            << " " << res.flux_total_upper_m2s1[(std::size_t)loc]
            << " " << td.RcLower_GV << " " << td.RcEffective_GV << " " << td.RcUpper_GV
            << " " << td.PenumbraWidth_GV << " " << td.THigh;
        for (int ic=0; ic<nCh; ++ic) {
          const std::size_t flat=(std::size_t)ic*(std::size_t)nLoc+(std::size_t)loc;
          out << " " << res.flux_ch_flat[flat] << " " << res.flux_ch_lower_flat[flat]
              << " " << res.flux_ch_upper_flat[flat];
        }
        out << "\n";
      }
    }
  }
}

static double WriteTerminationSummary_(const EarthUtil::AmpsParam& prm,
                                       const std::vector<double>& E_MeV,
                                       int nLoc,const DensityResultBuffers& res) {
  const int nE=static_cast<int>(E_MeV.size());
  double maximumUnresolved=0.0;
  std::ofstream out;
  if (prm.densitySpectrum.saveTerminationSummary) {
    out.open(DensityOutputFileName_("mode3d_termination_summary"));
    out << std::setprecision(std::numeric_limits<double>::max_digits10);
    out << "TITLE=\"Mode3D trajectory termination summary\"\n";
    out << "VARIABLES=\"location_index\" \"E_MeV\" \"N_sampled\" \"N_retried\" "
        << "\"N_resolved\" \"N_allowed\" \"T\" \"T_lower\" \"T_upper\" "
        << "\"unresolved_fraction\"";
    for (int it=0;it<Earth::FluxNumerics::kTerminationCount;++it) {
      const Earth::GridlessMode::TrajectoryTermination termination=
          static_cast<Earth::GridlessMode::TrajectoryTermination>(it);
      out << " \"N_" << Earth::GridlessMode::TrajectoryTerminationName(termination) << "\"";
    }
    out << "\nZONE T=\"termination\" I=" << nLoc*nE << " F=POINT\n";
  }
  for (int loc=0;loc<nLoc;++loc) {
    for (int ie=0;ie<nE;++ie) {
      const std::size_t flat=(std::size_t)loc*(std::size_t)nE+(std::size_t)ie;
      maximumUnresolved=std::max(maximumUnresolved,res.unresolved_fraction_flat[flat]);
      if (!out.is_open()) continue;
      out << loc << " " << E_MeV[(std::size_t)ie] << " " << res.sampled_flat[flat]
          << " " << res.retried_flat[flat] << " " << res.resolved_flat[flat]
          << " " << res.allowed_flat[flat] << " " << res.T_flat[flat]
          << " " << res.T_lower_flat[flat] << " " << res.T_upper_flat[flat]
          << " " << res.unresolved_fraction_flat[flat];
      for (int it=0;it<Earth::FluxNumerics::kTerminationCount;++it)
        out << " " << res.termination_flat[flat*(std::size_t)Earth::FluxNumerics::kTerminationCount+
                                            (std::size_t)it];
      out << "\n";
    }
  }
  return maximumUnresolved;
}

} // anonymous namespace

namespace Earth {
namespace Mode3D {

void SetDensityOutputFileSuffix(const std::string& suffix) {
  gDensityOutputFileSuffix = suffix;
}

int RunDensityAndFlux(const EarthUtil::AmpsParam& prm) {
  int mpiRank=0, mpiSize=1;
  MPI_Comm_rank(MPI_GLOBAL_COMMUNICATOR,&mpiRank);
  MPI_Comm_size(MPI_GLOBAL_COMMUNICATOR,&mpiSize);

  // Start a new sample-weighted magnetic-field accuracy diagnostic for this density/flux
  // calculation.  ResetDipoleMagneticFieldErrorStatistics() enables collection only for
  // mesh-backed DIPOLE runs, so no extra analytic evaluations are made for other models.
  ResetDipoleMagneticFieldErrorStatistics(prm);

  const std::string outputMode = EarthUtil::ToUpper(prm.output.mode);
  if (outputMode!="POINTS" && outputMode!="TRAJECTORY" && outputMode!="SHELLS") {
    exit(__LINE__,__FILE__,"Mode3D density/flux supports OUTPUT_MODE POINTS, TRAJECTORY, or SHELLS");
  }

  // Time-dependent boundary spectra are selected once per magnetic-field snapshot.  This
  // matches the snapshot contract: every call to RunDensityAndFlux() sees one already
  // assembled compact field snapshot and should fold it with the boundary spectrum at
  // the same epoch.  TRAJECTORY per-sample epochs are intentionally not used to mutate
  // gSpectrum because the global arrays are not rebuilt per spacecraft sample.
  ::gSpectrum.SetEvaluationEpochUTC(prm.field.epoch);

  const int nZenith = 24;
  const int nAz     = 48;
  const std::vector<V3> dirsFull = BuildDirGrid_(nZenith,nAz);
  const std::vector<double> E_MeV = BuildEnergyGrid_MeV_(prm);
  const int nE = (int)E_MeV.size();

  int nDirsUse = (int)dirsFull.size();
  if (prm.densitySpectrum.maxParticlesPerPoint > 0 && nE > 0) {
    nDirsUse = std::max(1, prm.densitySpectrum.maxParticlesPerPoint / nE);
    nDirsUse = std::min(nDirsUse, (int)dirsFull.size());
  }
  const std::vector<V3> dirsUse = SelectDirectionsDeterministic_(dirsFull,nDirsUse);

  int nLon=1,nLat=1,nPtsShell=1,nLoc=0;
  double res_deg = prm.output.shellRes_deg;
  if (outputMode=="SHELLS") {
    if (prm.output.shellAlt_km.empty()) exit(__LINE__,__FILE__,"OUTPUT_MODE=SHELLS requires SHELL_ALT_KM");
    if (!(res_deg > 0.0)) exit(__LINE__,__FILE__,"SHELL_RES_DEG must be > 0");
    nLon = ShellNLon_(res_deg);
    nLat = ShellNLat_(res_deg);
    nPtsShell = nLon*nLat;
    nLoc = nPtsShell*(int)prm.output.shellAlt_km.size();
  }
  else {
    nLoc = (int)prm.output.points.size();
  }

  if (nLoc <= 0) {
    exit(__LINE__,__FILE__,"Mode3D density/flux: no observation points or shell cells are defined");
  }

  const DensityParallelBackend_ densityBackendForBanner = ResolveDensityParallelBackend_(prm);
  const int densityThreadsForBanner = ResolveDensityThreadCount_(prm,densityBackendForBanner);
  const Earth::Mode3D::MpiScheduler mpiSchedulerForBanner =
      Earth::Mode3D::ResolveMpiScheduler(prm,"Mode3D density/flux");
  const long long mpiDynamicChunkForBanner = Earth::Mode3D::ResolveMpiDynamicChunk(
      prm,densityThreadsForBanner,static_cast<long long>(nLoc));

  if (mpiRank==0) {
    std::cout << "================ Mode3D mesh density & flux ================\n";
    std::cout << "Field model     : " << prm.field.model << "\n";
    std::cout << "Epoch           : " << prm.field.epoch << "\n";
    std::cout << "Species         : " << prm.species.name << " (q=" << prm.species.charge_e
              << " e, m=" << prm.species.mass_amu << " amu)\n";
    std::cout << "Energy grid     : [" << prm.densitySpectrum.Emin_MeV << ", "
              << prm.densitySpectrum.Emax_MeV << "] MeV, Npoints=" << nE << "\n";
    std::cout << "Transmission    : " << EarthUtil::ToUpper(prm.densitySpectrum.transmissionMode);
    if (EarthUtil::ToUpper(prm.densitySpectrum.transmissionMode)!="DIRECT") {
      std::cout << " (log-rigidity scan";
      if (prm.densitySpectrum.transmissionScanN > 0)
        std::cout << ", requested N=" << prm.densitySpectrum.transmissionScanN;
      if (prm.densitySpectrum.transmissionRefineN > 0)
        std::cout << ", refine N=" << prm.densitySpectrum.transmissionRefineN << " reserved";
      std::cout << ")";
    }
    std::cout << "\n";
    std::cout << "Directions      : " << dirsUse.size() << " / " << dirsFull.size()
              << " (" << nZenith << "x" << nAz << ")";
    if (prm.densitySpectrum.maxParticlesPerPoint > 0)
      std::cout << " [DS_MAX_PARTICLES=" << prm.densitySpectrum.maxParticlesPerPoint << "]";
    std::cout << "\n";
    std::cout << "Boundary mode   : " << prm.densitySpectrum.boundaryMode << "\n";
    std::cout << "Output mode     : " << outputMode << ", N_locations=" << nLoc << "\n";
    std::cout << "MPI ranks       : " << mpiSize << " (compact global field arrays)\n";
    std::cout << "Density backend : " << DensityParallelBackendName_(densityBackendForBanner)
              << ", threads/MPI rank=" << densityThreadsForBanner << "\n";
    std::cout << "MPI scheduler   : " << Earth::Mode3D::MpiSchedulerName(mpiSchedulerForBanner) << "\n";
    if (mpiSchedulerForBanner == Earth::Mode3D::MpiScheduler::DYNAMIC) {
      std::cout << "MPI dyn chunk   : " << mpiDynamicChunkForBanner
                << " global location(s) per atomic fetch\n";
    }
    std::cout << "Progress bar    : ON\n";
    if (!gDensityOutputFileSuffix.empty()) std::cout << "Output suffix   : " << gDensityOutputFileSuffix << "\n";
    std::cout << "============================================================\n";
    std::cout.flush();
  }

  DensityResultBuffers res = ComputeAllLocations_(prm,nLoc,nLon,nLat,res_deg,nPtsShell,E_MeV,dirsUse);

  if (mpiRank==0) {
    const double maximumUnresolved=WriteTerminationSummary_(prm,E_MeV,nLoc,res);
    std::cout << "[mode3d-density][termination] max unresolved fraction = "
              << std::setprecision(17) << maximumUnresolved << " (tolerance="
              << prm.densitySpectrum.unresolvedTolerance << ")\n";
    if (maximumUnresolved>prm.densitySpectrum.unresolvedTolerance) {
      std::cout << "[mode3d-density][termination][WARNING] unresolved fraction exceeds tolerance; "
                << "nominal values are not replaced by physical zero and bounds remain in output.\n";
      if (prm.densitySpectrum.failOnUnresolved)
        exit(__LINE__,__FILE__,"Mode3D density unresolved fraction exceeds DS_UNRESOLVED_TOL");
    }
    if (outputMode=="POINTS" || outputMode=="TRAJECTORY") {
      WritePointOutputs_(prm,E_MeV,res);
      std::cout << "Wrote: " << DensityOutputFileName_("mode3d_points_density") << "\n";
      std::cout << "Wrote: " << DensityOutputFileName_("mode3d_points_spectrum") << "\n";
      std::cout << "Wrote: " << DensityOutputFileName_("mode3d_points_flux") << "\n";
    }
    else {
      WriteShellOutputs_(prm,nLon,nLat,res_deg,nPtsShell,E_MeV,res);
      std::cout << "Wrote: " << prm.output.shellAlt_km.size()
                << " Mode3D shell density/flux file(s).\n";
    }
    std::cout.flush();
  }

  // Every TraceAllowedMesh() call has returned, so all per-trajectory field evaluators
  // have merged their local DIPOLE samples.  Perform one global reduction and print the
  // mean/max interpolation error and the location of the maximum before returning.
  ReportDipoleMagneticFieldErrorStatistics("Mode3D density and flux");

  MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
  return 0;
}

} // namespace Mode3D
} // namespace Earth
