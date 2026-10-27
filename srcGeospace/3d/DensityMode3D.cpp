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
// Geospace::Mode3D::TraceTrajectoryMesh(), which uses the AMPS AMR mesh field already
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
#include "GlobalMagneticField.h"

#include "pic.h"
#include "Earth.h"

#include "../gridless/GridlessParticleMovers.h"  // V3 helpers
#include "../gridless/AnisotropicSpectrum.h"     // EvalAnisotropyFactor
#include "../boundary/spectrum.h"                // ::gSpectrum
#include "../util/amps_param_parser.h"
#include "../util/FluxNumerics.h"                // common units, grids, quadrature, access accounting
#include "../util/BoundaryProducts.h"            // shared Step-6 product integrator
#include "../util/SWMFCoupledProductsContract.h" // Step-11 identity and manifest accounting
#include "../util/CommonPhysicsRelease.h"        // Step-12 cross-build physics identity

// Batched engine: batched trajectory engine (fiber driver support).  See the
// "BATCHED TRAJECTORY ENGINE" section below for how this module uses it.
#include "../gpu/BatchTracer.h"
#include "../gpu/DriverSupport.h"
#include "../gpu/Fiber.h"

// Mesh residency for engine-routed Mode3D traces.  gpu/Mode3DMeshExport.h belongs to the
// Mode3D tracer port (TraceTrajectoryMesh -> engine, kFieldMesh requests); its contract
// is that every Mode3D driver calls Mode3DMesh::EnsureUploaded() on the main thread
// before its first trajectory fiber starts, because BatchTracer::SetMesh() cannot run
// while trajectories are in flight.  The include is optional so this file builds with
// and without that port; see RunDensityAndFlux().
#if defined(__has_include)
#if __has_include("../gpu/Mode3DMeshExport.h")
#include "../gpu/Mode3DMeshExport.h"
#define SRC_EARTH_DENSITY_MODE3D_HAVE_MESH_EXPORT_ 1
#endif
#endif

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
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
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
static constexpr double QE       = Geospace::FluxNumerics::kElementaryCharge_C;
static constexpr double AMU      = Geospace::FluxNumerics::kAtomicMassUnit_kg;

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

using DensityParallelBackend_ = Geospace::Mode3D::ParallelBackend;

static const char* DensityParallelBackendName_(DensityParallelBackend_ backend) {
  return Geospace::Mode3D::ParallelBackendName(backend);
}

static DensityParallelBackend_ ResolveDensityParallelBackend_(const GeospaceUtil::AmpsParam& prm) {
  return Geospace::Mode3D::ResolveParallelBackend(prm,"Mode3D density/flux");
}

static int ResolveDensityThreadCount_(const GeospaceUtil::AmpsParam& prm,
                                      DensityParallelBackend_ backend) {
  return Geospace::Mode3D::ResolveParallelThreadCount(prm,backend);
}

static void ApplyWideAffinityForDirectDensityThreadsOnce_(DensityParallelBackend_ backend,
                                                          int densityThreadCount) {
  Geospace::Mode3D::ApplyWideAffinityForDirectThreadsOnce(backend,densityThreadCount,
                                                       "Mode3D density/flux");
}

//======================================================================================
// BATCHED TRAJECTORY ENGINE
//======================================================================================
// With AMPS_EARTH_ENGINE unset or BATCH (Geospace::GPU::DriverBatchEnabled()), the
// rank-local worker pools of this module (OpenMP teams, direct std::thread workers, or
// the serial loop) are replaced by cooperative fibers running on the MPI rank's main
// thread (gpu/Fiber.h, gpu/DriverSupport.h).  Every trajectory requested by a fiber
// through Geospace::Mode3D::TraceTrajectoryMesh() is queued, the fiber is suspended, and
// the batched engine integrates all queued trajectories of the rank together (an OpenMP
// host kernel or an OpenACC GPU kernel).  (Where the Mode3D tracer is not routed to the
// engine, a request made in a fiber is simply traced synchronously by the legacy
// tracer; the driver structure and the results are the same, only without overlap.)
// With AMPS_EARTH_ENGINE=LEGACY no fiber is ever created, Fiber::InFiber() is false
// everywhere, and every original OpenMP/std::thread/serial path below executes
// unchanged.
//
// Converted loops (all keep their original code as the LEGACY branch):
//
//   ComputeT_atEnergy_Mode3D_      direction loop     -> FiberLoop over directions
//   computeOneLocation (lambda)    energy loop        -> FiberLoop over energies
//   computeGlobalRange /
//   computeWorkListRange           location workers   -> RunDriverFibers over the
//                                                        rank's location stream
//   computeSharedEpochTaskRange /
//   computeSharedEpochTaskList     (location,energy)  -> RunDriverFibers over the
//                                  task workers          rank's task stream
//
// The nesting is therefore  location/task fibers -> energy fibers -> direction fibers,
// so a single rank keeps (locations x energies x directions) trajectories in flight
// instead of one per OS thread.  The fiber budget (AMPS_EARTH_FIBERS) bounds the
// number of live fibers; when it is exhausted Fiber::ParallelFor runs a nested loop
// inline in the calling fiber, which is still correct (just less concurrent).
//
// WHY THE RESULTS ARE BITWISE IDENTICAL
//   * Every trajectory is traced by the same TraceDensityDirection_Mode3D_() with the
//     same seed, parameter block, budgets, numerical retry and unresolved-extension
//     policy; only the moment at which the call returns changes.
//   * Per-direction results are written into per-direction slots and folded into the
//     AccessAccumulator by the original serial loop in increasing direction index, so
//     the order-dependent floating-point weight sum and the integer counters
//     (sampled/resolved/allowed/retried/termination histogram) are unchanged.
//   * Per-energy accumulators are written into per-energy slots and the per-location
//     product integration (EvaluateProductSet_) runs after the energy loop exactly as
//     before.  Each (location,energy) or location result is ASSIGNED to its own slot of
//     the rank-local buffers; no slot is written by two tasks and nothing is summed
//     across tasks inside a rank, so the order in which tasks finish is irrelevant.
//   * The MPI assignment of work is unchanged for STATIC and BLOCK_CYCLIC.  DYNAMIC
//     remains dynamic (fibers fetch chunks from the same RMA counter on demand); since
//     a slot is written by exactly one rank and the MPI_SUM reductions add only exact
//     zeros from the other ranks, the reduced buffers do not depend on which rank
//     computed which task -- the same property the original DYNAMIC scheduler relies on.
//
// FIBER / THREAD-SAFETY ASSUMPTIONS
//   * All fiber code (including the MPI RMA chunk fetches of DynamicMpiLocationScheduler
//     and the progress printing in the service hook) runs on the rank's main thread, so
//     MPI is never called from another OS thread (MPI_THREAD_MULTIPLE is not needed) and
//     the plain (non-atomic) counters/slots written by fibers are race free.
//   * A fiber may be suspended only inside a trajectory request.  No mutex of this file
//     is held across a trajectory call (the diagnostic mutex in
//     TraceDensityDirection_Mode3D_ is taken only around the stderr write, and the SPICE
//     shell-coordinate mutex is used only before the workers start), so fibers cannot
//     deadlock on one OS thread.
//   * The nested-OpenMP guards (InsideOpenMPParallelRegion_/InsideDirectDensityWorker_)
//     are not consulted on the batch path: the batch branches are selected by
//     Fiber::InFiber(), which is false inside any OpenMP region, and the batch drivers
//     never install DirectDensityWorkerScope_ on the main thread.
//   * This module has no thread_local caches or omp_get_thread_num()-indexed scratch
//     arrays on the trajectory path; the Mode3D field evaluator is constructed per
//     trajectory call inside Geospace::Mode3D::TraceTrajectoryMesh().
//======================================================================================

// Number of host OpenMP threads the engine may use for this rank: the number of OS
// threads the original backend would have used for the same trajectory loops.  SERIAL
// (and an OPENMP request in a build without OpenMP) used one thread.
static int DensityBatchHostThreads_(DensityParallelBackend_ backend,int densityThreadCount) {
  if (backend==DensityParallelBackend_::SERIAL) return 1;
#ifndef _OPENMP
  if (backend==DensityParallelBackend_::OPENMP) return 1;
#endif
  return std::max(1,densityThreadCount);
}

// One rank's stream of global work ids for the batch drivers (locations in the
// whole-location mode, flattened (location,energy) ids in the shared-epoch mode).
//
// It yields exactly the ids the original scheduler gives this rank: the precomputed
// STATIC/BLOCK_CYCLIC list in list order, or DYNAMIC chunks fetched from the RMA counter
// in fetch order.  Unlike gpu/TaskStream.h it uses the communicator of this module
// (MPI_GLOBAL_COMMUNICATOR, a singleton in EPOCH_PARALLEL TRAJECTORY runs) instead of
// MPI_COMM_WORLD.  Next() is called by fibers, i.e. on the main thread, so the RMA
// fetch is a main-thread MPI call.  Stop() makes every later Next() return -1; the batch
// drivers call it when a fiber fails so the remaining fibers do not start new work.
class DensityBatchWorkStream_ {
 public:
  // DYNAMIC: chunks of [0,nGlobal) fetched from `dynamicScheduler` (not owned).
  DensityBatchWorkStream_(Geospace::Mode3D::DynamicMpiLocationScheduler* dynamicScheduler,
                          long long nGlobal)
      : dynamic_(dynamicScheduler), total_(nGlobal) {}
  // STATIC/BLOCK_CYCLIC: the explicit list of this rank's ids, in processing order.
  template <class Id>
  explicit DensityBatchWorkStream_(const std::vector<Id>& ids) : total_(0) {
    list_.assign(ids.begin(),ids.end());
  }

  long long Next() {
    if (stopped_) return -1;
    if (dynamic_!=nullptr) {
      if (next_>=end_) {
        if (exhausted_) return -1;
        const long long start=dynamic_->FetchNextChunkStart();
        if (start>=total_) { exhausted_=true; return -1; }
        next_=start;
        end_=std::min(start+dynamic_->ChunkSize(),total_);
      }
      return next_++;
    }
    if (listPos_>=list_.size()) return -1;
    return list_[listPos_++];
  }

  void Stop() { stopped_=true; }

 private:
  Geospace::Mode3D::DynamicMpiLocationScheduler* dynamic_=nullptr;
  long long total_=0;
  long long next_=0,end_=0;
  bool exhausted_=false;
  std::vector<long long> list_;
  std::size_t listPos_=0;
  bool stopped_=false;
};

// Convert kinetic energy [J] to rigidity [GV].  This is copied intentionally from the
// shared numerical module so identical energies produce identical rigidity bins
// independent of whether the field backend is direct or mesh based.
static double RigidityFromEnergy_GV_(double E_J, double qabs_C, double m0_kg) {
  return Geospace::FluxNumerics::RigidityFromEnergyGV(E_J,qabs_C,m0_kg);
}
//--------------------------------------------------------------------------------------
// Output naming
//--------------------------------------------------------------------------------------
static std::string gDensityOutputFileSuffix;
static std::vector<std::string> gLastDensityFluxArtifactFiles;
static Geospace::SWMFCoupledProducts::ProductRunSummary gLastDensityFluxRunSummary;

static std::string DensityOutputFileName_(const char* stem) {
  return std::string(stem) + gDensityOutputFileSuffix + ".dat";
}

static std::ofstream OpenDensityArtifact_(const std::string& fileName) {
  std::ofstream out(fileName.c_str(),std::ios::out|std::ios::trunc);
  if (!out)
    throw std::runtime_error("Cannot open Mode3D density/flux artifact: "+fileName);
  return out;
}

static void CloseAndRecordDensityArtifact_(std::ofstream& out,
                                           const std::string& fileName) {
  // ostream errors can be delayed until the userspace buffer is flushed.  A coupled
  // snapshot must never advertise PASS merely because construction of an ofstream
  // succeeded, so explicitly flush and close before enrolling the file in the Step-11
  // transaction.  Duplicate enrollment is rejected rather than silently hidden.
  out.flush();
  if (!out)
    throw std::runtime_error("Failed while flushing Mode3D product: "+fileName);
  out.close();
  if (out.fail())
    throw std::runtime_error("Failed while closing Mode3D product: "+fileName);
  if (std::find(gLastDensityFluxArtifactFiles.begin(),
                gLastDensityFluxArtifactFiles.end(),fileName)!=
      gLastDensityFluxArtifactFiles.end())
    throw std::runtime_error("Mode3D product was recorded twice: "+fileName);
  gLastDensityFluxArtifactFiles.push_back(fileName);
}

static Geospace::SWMFCoupledProducts::ProductControl BuildProductControl_(
    const GeospaceUtil::AmpsParam& prm) {
  namespace CP=Geospace::SWMFCoupledProducts;
  CP::ProductControl control;
  control.outputMode=GeospaceUtil::ToUpper(prm.output.mode);
  control.speciesName=prm.species.name;
  control.charge_e=prm.species.charge_e;
  control.mass_amu=prm.species.mass_amu;
  control.boundaryMode=GeospaceUtil::ToUpper(prm.densitySpectrum.boundaryMode);
  control.transmissionMode=GeospaceUtil::ToUpper(prm.densitySpectrum.transmissionMode);
  control.minimumEnergy_MeV=prm.densitySpectrum.Emin_MeV;
  control.maximumEnergy_MeV=prm.densitySpectrum.Emax_MeV;
  control.energyIntervals=prm.densitySpectrum.nIntervals;
  control.transmissionScanPoints=prm.densitySpectrum.transmissionScanN;
  control.maximumParticlesPerPoint=prm.densitySpectrum.maxParticlesPerPoint;
  control.directionNMu=prm.densitySpectrum.directionNMu;
  control.directionNPhi=prm.densitySpectrum.directionNPhi;
  control.energySpacing=(prm.densitySpectrum.spacing==
      GeospaceUtil::DensitySpectrumParam::Spacing::LOG) ? "LOG" : "LINEAR";
  control.spectrumType=prm.particleSpectrum.typeName.empty()
      ? std::string("UNKNOWN") : prm.particleSpectrum.typeName;
  control.energyBasis=(::gSpectrum.EnergyCoordinateBasis()==
      Geospace::BoundaryProducts::EnergyBasis::PerNucleon)
      ? "PER_NUCLEON" : "PER_PARTICLE";
  control.spectrumMassNumber=::gSpectrum.MassNumber();
  control.intensityUnit=::gSpectrum.IntensityUnitLabel();
  control.spectrumRelativeUncertainty=::gSpectrum.RelativeUncertainty();
  for (std::map<std::string,std::string>::const_iterator it=prm.spectrum.begin();
       it!=prm.spectrum.end();++it)
    control.spectrumKeyValues.push_back(*it);
  control.spectrumTableEnergy_MeV=::gSpectrum.TableEnergy_MeV();
  control.spectrumTableIntensityPerMeV=::gSpectrum.TableFlux_perMeV();
  for (std::size_t i=0;i<prm.fluxChannels.size();++i) {
    CP::EnergyChannelDefinition channel;
    channel.name=prm.fluxChannels[i].name;
    channel.lower_MeV=prm.fluxChannels[i].E1_MeV;
    channel.upper_MeV=prm.fluxChannels[i].E2_MeV;
    control.channels.push_back(channel);
  }
  for (std::size_t i=0;i<prm.detectorResponses.size();++i) {
    CP::DetectorResponseDefinition response;
    response.name=prm.detectorResponses[i].name;
    response.lower_MeV=prm.detectorResponses[i].E1_MeV;
    response.upper_MeV=prm.detectorResponses[i].E2_MeV;
    response.geometricFactor_m2_sr=prm.detectorResponses[i].geometricFactor_m2_sr;
    control.detectorResponses.push_back(response);
  }
  control.coordinateFrame=(control.outputMode=="TRAJECTORY")
      ? prm.output.trajFrame : prm.output.coords;
  if (control.outputMode=="TRAJECTORY") {
    // Preserve ephemeris time and position together.  The flattened points alone are
    // insufficient: two spacecraft samples can occupy the same cell at different UTCs
    // while selecting different field/spectrum snapshots in a later campaign.
    for (std::size_t it=0;it<prm.output.trajectories.size();++it) {
      const GeospaceUtil::SpacecraftTrajectory& trajectory=prm.output.trajectories[it];
      for (std::size_t is=0;is<trajectory.samples.size();++is) {
        CP::ObservationDefinition observation;
        observation.epochUTC=trajectory.samples[is].timeUTC;
        observation.x_km=trajectory.samples[is].xGSM_m.x/1000.0;
        observation.y_km=trajectory.samples[is].xGSM_m.y/1000.0;
        observation.z_km=trajectory.samples[is].xGSM_m.z/1000.0;
        control.observations.push_back(observation);
      }
    }
  }
  else if (control.outputMode=="POINTS") {
    for (std::size_t i=0;i<prm.output.points.size();++i) {
      CP::ObservationDefinition observation;
      observation.x_km=prm.output.points[i].x;
      observation.y_km=prm.output.points[i].y;
      observation.z_km=prm.output.points[i].z;
      control.observations.push_back(observation);
    }
  }
  control.shellAltitude_km=prm.output.shellAlt_km;
  control.shellResolution_deg=prm.output.shellRes_deg;
  control.shellGeometry=prm.output.shellGeometry;
  CP::ValidateProductControl(control);
  return control;
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
static std::vector<double> BuildEnergyGrid_MeV_(const GeospaceUtil::AmpsParam& prm) {
  const int nLegacy = prm.densitySpectrum.nPoints();
  if (nLegacy < 2) exit(__LINE__,__FILE__,"Mode3D density requires DS_NINTERVALS >= 1");

  const double Emin = prm.densitySpectrum.Emin_MeV;
  const double Emax = prm.densitySpectrum.Emax_MeV;
  const std::string mode=GeospaceUtil::ToUpper(prm.densitySpectrum.transmissionMode);
  const bool scan=(mode=="SCAN" || mode=="ADAPTIVE");
  try {
    return Geospace::BoundaryProducts::BuildEnergyCoordinateGridMeV(
        Emin,Emax,nLegacy,
        prm.densitySpectrum.spacing==GeospaceUtil::DensitySpectrumParam::Spacing::LOG
          ? Geospace::FluxNumerics::EnergySpacing::Log
          : Geospace::FluxNumerics::EnergySpacing::Linear,
        scan,prm.densitySpectrum.transmissionScanN,
        prm.densitySpectrum.transmissionMaxN,
        std::fabs(prm.species.charge_e)*QE,prm.species.mass_amu*AMU,
        ::gSpectrum.Units());
  }
  catch (const std::exception& e) { exit(__LINE__,__FILE__,e.what()); }
  return std::vector<double>();
}

// Keep the gridless density convention (24 zenith x 48 azimuth) for backend-to-backend
// comparability.  The point set is deterministic, so repeated runs and different MPI
// decompositions produce reproducible transmissivity values.
static std::vector<V3> BuildDirGrid_(int nZenith,int nAz) {
  const std::vector<Geospace::FluxNumerics::DirectionSample> shared=
      Geospace::FluxNumerics::BuildEqualSolidAngleDirections(nZenith,nAz);
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
  return Geospace::FluxNumerics::SelectDeterministic(all,nUse);
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

static V3 ShellLocationGSM_m_(const GeospaceUtil::AmpsParam& prm,
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
  if (GeospaceUtil::ToUpper(prm.field.model) == "DIPOLE") return xFixed;

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

static V3 LocationByIndex_m_(const GeospaceUtil::AmpsParam& prm,
                             int loc,int nLon,int nLat,double res_deg,int nPtsShell) {
  const std::string mode = GeospaceUtil::ToUpper(prm.output.mode);
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

// Convert the backend-neutral access curve into every Step-6 scalar and spectral
// product. Mode3D (including the SWMF coupled caller) supplies only A(E) and its
// unresolved bounds; all unit conversions, boundary uncertainty, energy-channel
// clipping, and detector folding are performed by the same production kernel used by
// the gridless backend.
static Geospace::BoundaryProducts::ProductSet EvaluateProductSet_(
    const GeospaceUtil::AmpsParam& prm,const std::vector<double>& E_MeV,
    const std::vector<double>& T,const std::vector<double>& TLower,
    const std::vector<double>& TUpper) {
  namespace BP=Geospace::BoundaryProducts;
  std::vector<BP::EnergyChannel> channels;
  channels.reserve(prm.fluxChannels.size());
  for (const GeospaceUtil::EnergyChannel& channel:prm.fluxChannels) {
    channels.push_back(BP::EnergyChannel(
        channel.name,channel.E1_MeV,channel.E2_MeV));
  }

  std::vector<BP::DetectorResponse> responses;
  responses.reserve(prm.detectorResponses.size());
  for (const GeospaceUtil::DetectorResponseChannel& input:prm.detectorResponses) {
    BP::DetectorResponse response;
    response.name=input.name;
    response.energy_MeV={input.E1_MeV,input.E2_MeV};
    response.relativeResponse={1.0,1.0};
    response.geometricFactor_m2_sr=input.geometricFactor_m2_sr;
    responses.push_back(response);
  }

  return BP::EvaluateIsotropicProducts(
      E_MeV,T,TLower,TUpper,prm.species.mass_amu*AMU,
      [](double E_J) { return ::gSpectrum.GetSpectrum(E_J); },channels,responses,
      ::gSpectrum.RelativeUncertainty(),::gSpectrum.Units());
}


//--------------------------------------------------------------------------------------
// Transmission diagnostics derived from T(E)
//--------------------------------------------------------------------------------------
using TransmissionDiagnostics_=Geospace::FluxNumerics::TransmissionDiagnostics;

static TransmissionDiagnostics_ ComputeTransmissionDiagnostics_(
    const GeospaceUtil::AmpsParam& prm,
    const std::vector<double>& E_MeV,
    const std::vector<double>& T) {
  // The diagnostic helper converts energy to rigidity and therefore requires total
  // particle kinetic energy.  Output/integration remain in the declared coordinate.
  std::vector<double> particleEnergy_MeV(E_MeV.size(),0.0);
  const Geospace::BoundaryProducts::SpectrumUnits units=::gSpectrum.Units();
  for (std::size_t i=0;i<E_MeV.size();++i)
    particleEnergy_MeV[i]=units.ParticleEnergyMeV(E_MeV[i]);
  return Geospace::FluxNumerics::ComputeTransmissionDiagnostics(
      particleEnergy_MeV,T,std::fabs(prm.species.charge_e)*QE,
      prm.species.mass_amu*AMU);
}

//--------------------------------------------------------------------------------------
// Transmissivity at one energy/location using the Mode3D mesh tracer
//--------------------------------------------------------------------------------------
static Geospace::GridlessMode::TrajectoryResult TraceDensityDirection_Mode3D_(
    const GeospaceUtil::AmpsParam& prm,
    const double x0_arr[3],
    const double v0_arr[3],
    double Rgv,
    double maxTrajTime_s,
    bool captureExit,
    int directionIndex) {
  using Geospace::GridlessMode::TrajectoryResult;

  // Opt-in forensics retain the exact deterministic seed and every attempted
  // budget before the strict product gate can abort. No MPI or spectrum mutation
  // is permitted here: this function also runs inside OpenMP workers. One locked
  // write per record prevents rank-local worker output from interleaving.
  //
  // Batched engine: in batch mode this function runs inside fibers on the rank's main
  // thread, and a fiber may be suspended inside Geospace::Mode3D::TraceTrajectoryMesh()
  // while the engine integrates the trajectory.  The diagnostic mutex below is taken
  // only around the stderr write (never across a trajectory call), so suspended fibers
  // can never hold it.  The retry/extension sequence is evaluated per call on the
  // fiber's own stack, so interleaving with other directions cannot change it.
  static const bool diagnostics=[]() {
    const char* value=std::getenv("AMPS_DENSITY_TRACE_DIAGNOSTICS");
    return value && std::string(value)=="1";
  }();
  auto record=[&](const TrajectoryResult& out,const char* stage,int pass,
                  double budget,int stepBudget,bool force) {
    // Full attempt histories are opt-in because a normal SA02 epoch contains thousands
    // of trajectories.  The final unresolved state is always emitted, however, so a
    // strict production failure retains the exact rejected trap sub-gates in amps.log.
    if (!diagnostics && !force) return;
    std::ostringstream line;
    line << std::setprecision(17)
         << "[density-trace] epoch=" << prm.field.epoch
         << " direction=" << directionIndex << " rigidity_GV=" << Rgv
         << " particle_energy_MeV=" << Geospace::FluxNumerics::EnergyFromRigidityMeV(
                Rgv,std::fabs(prm.species.charge_e)*QE,prm.species.mass_amu*AMU)
         << " x_m=" << x0_arr[0] << ',' << x0_arr[1] << ',' << x0_arr[2]
         << " backtrace_unit=" << v0_arr[0] << ',' << v0_arr[1] << ',' << v0_arr[2]
         << " stage=" << stage << " pass=" << pass
         << " budget_s=" << budget << " step_budget=" << stepBudget
         << " termination=" << Geospace::GridlessMode::TrajectoryTerminationName(out.termination)
         << " trace_s=" << out.traceTime_s << " steps=" << out.steps
         << " retries=" << out.retryCount << " mirrors=" << out.mirrorPoints
         << " bounces=" << out.bounceCycles << " drift_turns=" << out.driftRevolutions
         << " momentum_spread=" << out.momentumRelativeSpread
         << " drift_radius_change_m=" << out.driftMeanRadiusChange_m
         << " drift_radius_span_m=" << out.driftMeanRadiusSpan_m
         << " drift_stable_comparisons=" << out.driftStableComparisons
         << " drift_recurrence_lag_turns=" << out.driftRecurrenceLagRevolutions
         << " minimum_outer_margin_m=" << out.minimumOuterMargin_m
         << " bounce_envelope_change_m=" << out.lastBounceEnvelopeChange_m
         << " drift_comparison_available=" << (out.driftComparisonAvailable ? 1 : 0)
         << " drift_previous_coverage=" << out.driftPreviousCoverageFraction
         << " drift_current_coverage=" << out.driftCurrentCoverageFraction
         << " drift_common_bins=" << out.driftCommonBins
         << " drift_matched_bins=" << out.driftMatchedBins
         << " drift_matched_fraction=" << out.driftMatchedBinFraction
         << " drift_radial_mismatches=" << out.driftRadialMismatchBins
         << " drift_latitude_mismatches=" << out.driftLatitudeMismatchBins
         << " drift_pitch_mismatches=" << out.driftPitchMismatchBins
         << " drift_radial_residual_norm=" << out.driftMaximumNormalizedRadialResidual
         << " drift_latitude_residual_norm=" << out.driftMaximumNormalizedLatitudeResidual
         << " drift_pitch_residual_norm=" << out.driftMaximumNormalizedPitchResidual
         << " drift_secular_ok=" << (out.driftSecularDriftPassed ? 1 : 0)
         << " drift_profiles_recurred=" << (out.driftProfilesRecurred ? 1 : 0) << '\n';
    static std::mutex diagnosticMutex;
    std::lock_guard<std::mutex> lock(diagnosticMutex);
    std::cerr << line.str();
    // Emit the compact Poincare series only for the final unresolved attempt.  The
    // states are recorded only when AMPS_TRAP_TURN_DIAGNOSTICS=1, so ordinary SA02
    // output and memory use are unchanged.  Keeping the selector here (where energy
    // and direction are known) avoids concurrent per-step files and guarantees each
    // line carries the exact deterministic trajectory seed needed for reproduction.
    if (force) {
      for (const auto& state:out.driftTurnStates) {
        std::cerr << std::setprecision(17)
             << "[trap-turn] epoch=" << prm.field.epoch
             << " direction=" << directionIndex
             << " rigidity_GV=" << Rgv
             << " particle_energy_MeV=" << Geospace::FluxNumerics::EnergyFromRigidityMeV(
                    Rgv,std::fabs(prm.species.charge_e)*QE,
                    prm.species.mass_amu*AMU)
             << " pass=" << pass
             << " revolution=" << state.revolution
             << " trace_s=" << state.traceTime_s
             << " accumulated_angle_rad=" << state.accumulatedAngle_rad
             << " x_m=" << state.x_m[0] << ',' << state.x_m[1] << ',' << state.x_m[2]
             << " p_unit=" << state.p_unit[0] << ',' << state.p_unit[1] << ','
                            << state.p_unit[2]
             << " b_unit=" << state.b_unit[0] << ',' << state.b_unit[1] << ','
                            << state.b_unit[2]
             << " radius_m=" << state.radius_m
             << " latitude=" << state.latitude
             << " pitch_cos2=" << state.pitchCos2
             << " outer_margin_m=" << state.outerMargin_m
             << " momentum_spread=" << state.momentumRelativeSpread << '\n';
      }
    }
    std::cerr.flush();
  };

  Geospace::Trajectory::RetryPolicy retryPolicy;
  retryPolicy.maximumNumericalRetries=prm.densitySpectrum.retryUnresolved ? 1 : 0;
  retryPolicy.unresolvedExtensionPasses=prm.densitySpectrum.retryUnresolved
      ? prm.densitySpectrum.unresolvedExtensionPasses : 0;
  retryPolicy.unresolvedExtensionFactor=
      prm.densitySpectrum.unresolvedExtensionFactor;
  const double primaryTime_s=(maxTrajTime_s>0.0)
      ? maxTrajTime_s
      : ((prm.densitySpectrum.maxTrajTime_s>0.0)
          ? prm.densitySpectrum.maxTrajTime_s : prm.numerics.maxTraceTime_s);

  // Keep numerical recovery and trace-time convergence as two auditable mechanisms.
  // INVALID_TIME_STEP/INVALID_FIELD/NUMERICAL_FAILURE may receive the established
  // smaller-step retry.  A valid TIME_LIMIT or STEP_LIMIT retains the original
  // integrator settings and is eligible only for the explicit larger-budget passes
  // below.  This prevents a timeout from being silently converted to a physical zero.
  auto RunOneBudget=[&](const GeospaceUtil::AmpsParam& attemptPrm,
                        double requestedTime_s,
                        double& actualTimeBudget_s) -> TrajectoryResult {
    actualTimeBudget_s=requestedTime_s;
    auto out=Geospace::Mode3D::TraceTrajectoryMesh(
        attemptPrm,x0_arr,v0_arr,Rgv,captureExit,requestedTime_s);
    if (!Geospace::Trajectory::ShouldRetryNumerical(out,0,retryPolicy)) return out;

    GeospaceUtil::AmpsParam retryPrm=attemptPrm;
    retryPrm.numerics.dtTrace_s=std::max(
        1.0e-12,retryPolicy.numericalDtScale*attemptPrm.numerics.dtTrace_s);
    retryPrm.numerics.maxSteps=Geospace::Trajectory::ScaledStepBudget(
        attemptPrm.numerics.maxSteps,retryPolicy.numericalStepScale,
        out,retryPolicy.numericalTimeScale*requestedTime_s);
    actualTimeBudget_s=retryPolicy.numericalTimeScale*requestedTime_s;
    out=Geospace::Mode3D::TraceTrajectoryMesh(
        retryPrm,x0_arr,v0_arr,Rgv,captureExit,actualTimeBudget_s);
    out.retryCount=1;
    return out;
  };

  double usedBudget_s=primaryTime_s;
  int finalStepBudget=prm.numerics.maxSteps;
  TrajectoryResult result=RunOneBudget(prm,primaryTime_s,usedBudget_s);
  record(result,"primary",0,usedBudget_s,prm.numerics.maxSteps,false);
  const int primaryTerminationCode=static_cast<int>(result.termination);
  const double primaryTraceTime_s=result.traceTime_s;
  result.primaryTerminationCode=primaryTerminationCode;
  result.primaryTraceTime_s=primaryTraceTime_s;
  result.traceExtensionCount=0;
  result.initialTraceLimit_s=primaryTime_s;
  result.finalTraceLimit_s=usedBudget_s;

  // Restart only unresolved limit outcomes from the same phase-space seed.  A true
  // continuation would be faster but would have to serialize private mover, field-cache,
  // and trap-detector state.  Restarting makes every extended result exactly reproducible
  // by an independent run configured with the same final limit and frozen field snapshot.
  for (int pass=1;pass<=retryPolicy.unresolvedExtensionPasses;++pass) {
    if (!Geospace::Trajectory::ShouldExtendUnresolved(result,pass-1,retryPolicy)) break;
    const double factor=std::pow(retryPolicy.unresolvedExtensionFactor,
                                 static_cast<double>(pass));
    const double targetTime_s=Geospace::Trajectory::ExtensionTimeBudget(
        primaryTime_s,pass,retryPolicy);
    if (!(targetTime_s>primaryTime_s) || !std::isfinite(targetTime_s)) break;

    GeospaceUtil::AmpsParam extensionPrm=prm;
    extensionPrm.numerics.maxSteps=Geospace::Trajectory::ScaledStepBudget(
        prm.numerics.maxSteps,factor,result,targetTime_s);
    double extensionBudget_s=targetTime_s;
    TrajectoryResult extended=RunOneBudget(
        extensionPrm,targetTime_s,extensionBudget_s);
    record(extended,"extension",pass,extensionBudget_s,extensionPrm.numerics.maxSteps,false);
    extended.primaryTerminationCode=primaryTerminationCode;
    extended.primaryTraceTime_s=primaryTraceTime_s;
    extended.traceExtensionCount=pass;
    extended.initialTraceLimit_s=primaryTime_s;
    extended.finalTraceLimit_s=extensionBudget_s;
    result=extended;
    finalStepBudget=extensionPrm.numerics.maxSteps;
  }
  if (!result.resolved())
    record(result,"final",result.traceExtensionCount,result.finalTraceLimit_s,
           finalStepBudget,true);
  return result;
}

static Geospace::FluxNumerics::AccessAccumulator ComputeT_atEnergy_Mode3D_(const GeospaceUtil::AmpsParam& prm,
                                        const V3& x0_m,
                                        double Rgv,
                                        const std::vector<V3>& dirs,
                                        double maxTrajTime_s,
                                        bool doAnisotropic,
                                        const GeospaceUtil::AnisotropyParam& anisoPar) {
  Geospace::FluxNumerics::AccessAccumulator total;
  if (dirs.empty()) return total;

  const double x0_arr[3]={x0_m.x,x0_m.y,x0_m.z};
  std::vector<double> weights(dirs.size(),0.0);
  std::vector<int> terminations(dirs.size(),
      static_cast<int>(Geospace::GridlessMode::TrajectoryTermination::NumericalFailure));
  std::vector<int> retried(dirs.size(),0);

  if (Geospace::GPU::Fiber::InFiber() && dirs.size()>1) {
    //----------------------------------------------------------------------------------
    // Batch mode: mirror of the OpenMP direction loop below.
    //
    // Each direction runs in its own child fiber of the calling (energy/task) fiber, so
    // all directions of this energy are in the engine's trajectory pool together.  The
    // body is the original loop body verbatim: it traces the same seed with the same
    // retry/extension policy and writes only this direction's slots of weights[],
    // terminations[] and retried[].  The serial aggregation after the branch is shared
    // with the LEGACY path and folds the slots in increasing direction index, so the
    // AccessAccumulator (including its floating-point weight sum) is bitwise identical
    // to a serial or OpenMP evaluation.  The direction index passed to the tracer is
    // the same deterministic index used for diagnostics in the original loop.
    //
    // This branch replaces the OpenMP team only where the original would have used one
    // (InFiber() is false inside OpenMP regions and in LEGACY runs, and a one-direction
    // list is serial in both paths).  When the original loop would have been serial
    // because an enclosing OpenMP/std::thread worker owned the energy or location, the
    // batch drivers have replaced that worker by a fiber, and fanning the directions
    // out here is what keeps the engine pool full.
    //----------------------------------------------------------------------------------
    Geospace::GPU::FiberLoop(0,static_cast<long long>(dirs.size()),[&](long long index,int) {
      const int idir=static_cast<int>(index);
      const V3& arrivalDir=dirs[(std::size_t)idir];
      const V3 vTry=mul(-1.0,arrivalDir);
      const double v0_arr[3]={vTry.x,vTry.y,vTry.z};

      auto result=TraceDensityDirection_Mode3D_(
          prm,x0_arr,v0_arr,Rgv,maxTrajTime_s,doAnisotropic,idir);
      retried[(std::size_t)idir]=
          (result.retryCount>0 || result.traceExtensionCount>0) ? 1 : 0;

      terminations[(std::size_t)idir]=static_cast<int>(result.termination);
      if (!result.resolved() || !result.allowed()) return;
      weights[(std::size_t)idir]=doAnisotropic
          ? EvalAnisotropyFactor(anisoPar,result.exitState.cosAlpha,result.exitState.x_exit_m)
          : 1.0;
    });
  }
  else {
    // Original (legacy) direction loop (LEGACY engine, or a single direction).
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(dirs,prm,x0_arr,Rgv,maxTrajTime_s,doAnisotropic,anisoPar,weights,terminations,retried) if(!InsideOpenMPParallelRegion_() && !InsideDirectDensityWorker_() && (int)dirs.size() > 1) schedule(dynamic)
#endif
  for (int idir=0; idir<(int)dirs.size(); ++idir) {
    const V3& arrivalDir=dirs[(std::size_t)idir];
    const V3 vTry=mul(-1.0,arrivalDir);
    const double v0_arr[3]={vTry.x,vTry.y,vTry.z};

    auto result=TraceDensityDirection_Mode3D_(
        prm,x0_arr,v0_arr,Rgv,maxTrajTime_s,doAnisotropic,idir);
    retried[(std::size_t)idir]=
        (result.retryCount>0 || result.traceExtensionCount>0) ? 1 : 0;

    terminations[(std::size_t)idir]=static_cast<int>(result.termination);
    if (!result.resolved() || !result.allowed()) continue;
    weights[(std::size_t)idir]=doAnisotropic
        ? EvalAnisotropyFactor(anisoPar,result.exitState.cosAlpha,result.exitState.x_exit_m)
        : 1.0;
  }
  }

  // Aggregate serially after the parallel trace loop.  This makes every count and
  // termination category deterministic and avoids a custom OpenMP reduction for the
  // structured result.
  for (std::size_t i=0;i<dirs.size();++i) {
    const Geospace::GridlessMode::TrajectoryTermination termination=
        static_cast<Geospace::GridlessMode::TrajectoryTermination>(terminations[i]);
    total.Record(termination,weights[i],retried[i]!=0);
  }
  return total;
}

struct DensityResultBuffers {
  std::vector<double> density_m3,density_lower_m3,density_upper_m3;
  std::vector<double> flux_total_m2s1,flux_total_lower_m2s1,flux_total_upper_m2s1;
  std::vector<double> flux_planar_m2s1,flux_planar_lower_m2s1,flux_planar_upper_m2s1;
  std::vector<double> T_flat,T_lower_flat,T_upper_flat,unresolved_fraction_flat;
  std::vector<double> flux_ch_flat,flux_ch_lower_flat,flux_ch_upper_flat;
  std::vector<double> detector_rate_flat,detector_rate_lower_flat,detector_rate_upper_flat;
  std::vector<int> sampled_flat,resolved_flat,allowed_flat,retried_flat;
  std::vector<int> termination_flat; // [loc][energy][TrajectoryTermination]
};

static Geospace::SWMFCoupledProducts::ProductRunSummary BuildRunSummary_(
    const GeospaceUtil::AmpsParam& prm,int nLoc,int nE,int nDirections,
    const DensityResultBuffers& result) {
  namespace CP=Geospace::SWMFCoupledProducts;
  CP::ProductRunSummary summary;
  summary.control=BuildProductControl_(prm);
  summary.spectrumEvaluationEpochUTC=prm.field.epoch;
  summary.activeSpectrumTableEpochUTC=::gSpectrum.ActiveTableEpochUTC();
  // Analytic spectra do not own a table row.  Recording the requested evaluation UTC
  // in both fields avoids an empty provenance token while still allowing tabulated
  // spectra to report the exact/interpolated source rows selected by cSpectrum.
  if (summary.activeSpectrumTableEpochUTC.empty())
    summary.activeSpectrumTableEpochUTC=prm.field.epoch;
  summary.spectrumTemporalStatus=Geospace::BoundaryProducts::TemporalStatusName(
      ::gSpectrum.LastTemporalStatus());
  summary.spectrumTemporalGap=::gSpectrum.LastTemporalSelectionCrossedGap();
  summary.spectrumTemporalFraction=::gSpectrum.LastTemporalInterpolationFraction();
  summary.locationCount=nLoc;
  summary.energyCount=nE;
  summary.directionCount=nDirections;
  summary.terminationCounts.assign(
      static_cast<std::size_t>(Geospace::FluxNumerics::kTerminationCount),0);
  summary.maximumUnresolvedFraction=0.0;
  summary.unresolvedTolerance=prm.densitySpectrum.unresolvedTolerance;

  const std::size_t nBins=static_cast<std::size_t>(nLoc)*
                          static_cast<std::size_t>(nE);
  for (std::size_t flat=0;flat<nBins;++flat) {
    summary.sampled+=result.sampled_flat[flat];
    summary.retried+=result.retried_flat[flat];
    summary.resolved+=result.resolved_flat[flat];
    summary.allowed+=result.allowed_flat[flat];
    summary.maximumUnresolvedFraction=std::max(
        summary.maximumUnresolvedFraction,result.unresolved_fraction_flat[flat]);
    for (int it=0;it<Geospace::FluxNumerics::kTerminationCount;++it)
      summary.terminationCounts[static_cast<std::size_t>(it)]+=
          result.termination_flat[
              flat*static_cast<std::size_t>(Geospace::FluxNumerics::kTerminationCount)+
              static_cast<std::size_t>(it)];
  }
  summary.valid=true;
  // Artifact completeness is validated only after every writer has closed.  At this
  // stage validate the physics/provenance/count portion independently so a closure
  // defect is distinguishable from a bad trajectory accounting defect.
  CP::ValidateRunSummary(summary,false);
  return summary;
}

static double MaximumAllowedWeight_(const GeospaceUtil::AmpsParam& prm,bool anisotropic) {
  if (!anisotropic) return 1.0;
  namespace BP=Geospace::BoundaryProducts;
  const std::string pad=GeospaceUtil::ToUpper(prm.anisotropy.padModel);
  const BP::PadModel padModel=pad=="SINALPHA_N" ? BP::PadModel::SinAlphaN :
      (pad=="COSALPHA_N" ? BP::PadModel::CosAlphaN :
       (pad=="BIDIRECTIONAL" ? BP::PadModel::Bidirectional : BP::PadModel::Isotropic));
  const BP::NormalizationMode padNorm=
      GeospaceUtil::ToUpper(prm.anisotropy.padNormalization)=="RAW"
      ? BP::NormalizationMode::Raw : BP::NormalizationMode::UnitMean;
  const BP::SpatialModel spatial=
      GeospaceUtil::ToUpper(prm.anisotropy.spatialModel)=="DAYSIDE_NIGHTSIDE"
      ? BP::SpatialModel::DaysideNightside : BP::SpatialModel::Uniform;
  const BP::NormalizationMode spatialNorm=
      GeospaceUtil::ToUpper(prm.anisotropy.spatialNormalization)=="RAW"
      ? BP::NormalizationMode::Raw : BP::NormalizationMode::UnitMean;
  return BP::MaximumPadWeight(padModel,prm.anisotropy.padExponent,padNorm)*
         BP::MaximumSpatialWeight(spatial,prm.anisotropy.daysideFactor,
                                  prm.anisotropy.nightsideFactor,spatialNorm);
}

static DensityResultBuffers ComputeAllLocations_(const GeospaceUtil::AmpsParam& prm,
                                                 int nLoc,
                                                 int nLon,int nLat,double res_deg,int nPtsShell,
                                                 const std::vector<double>& E_MeV,
                                                 const std::vector<V3>& dirsUse) {
  int mpiRank=0, mpiSize=1;
  MPI_Comm_rank(MPI_GLOBAL_COMMUNICATOR,&mpiRank);
  MPI_Comm_size(MPI_GLOBAL_COMMUNICATOR,&mpiSize);

  const int nE  = (int)E_MeV.size();
  const int nCh = (int)prm.fluxChannels.size();
  const int nDetector = (int)prm.detectorResponses.size();

  DensityResultBuffers local;
  local.density_m3.assign((std::size_t)nLoc,0.0);
  local.density_lower_m3.assign((std::size_t)nLoc,0.0);
  local.density_upper_m3.assign((std::size_t)nLoc,0.0);
  local.flux_total_m2s1.assign((std::size_t)nLoc,0.0);
  local.flux_total_lower_m2s1.assign((std::size_t)nLoc,0.0);
  local.flux_total_upper_m2s1.assign((std::size_t)nLoc,0.0);
  local.flux_planar_m2s1.assign((std::size_t)nLoc,0.0);
  local.flux_planar_lower_m2s1.assign((std::size_t)nLoc,0.0);
  local.flux_planar_upper_m2s1.assign((std::size_t)nLoc,0.0);
  local.T_flat.assign((std::size_t)nLoc*(std::size_t)nE,0.0);
  local.T_lower_flat.assign((std::size_t)nLoc*(std::size_t)nE,0.0);
  local.T_upper_flat.assign((std::size_t)nLoc*(std::size_t)nE,0.0);
  local.unresolved_fraction_flat.assign((std::size_t)nLoc*(std::size_t)nE,0.0);
  local.flux_ch_flat.assign((std::size_t)nCh*(std::size_t)nLoc,0.0);
  local.flux_ch_lower_flat.assign((std::size_t)nCh*(std::size_t)nLoc,0.0);
  local.flux_ch_upper_flat.assign((std::size_t)nCh*(std::size_t)nLoc,0.0);
  local.detector_rate_flat.assign((std::size_t)nDetector*(std::size_t)nLoc,0.0);
  local.detector_rate_lower_flat.assign((std::size_t)nDetector*(std::size_t)nLoc,0.0);
  local.detector_rate_upper_flat.assign((std::size_t)nDetector*(std::size_t)nLoc,0.0);
  local.sampled_flat.assign((std::size_t)nLoc*(std::size_t)nE,0);
  local.resolved_flat.assign((std::size_t)nLoc*(std::size_t)nE,0);
  local.allowed_flat.assign((std::size_t)nLoc*(std::size_t)nE,0);
  local.retried_flat.assign((std::size_t)nLoc*(std::size_t)nE,0);
  local.termination_flat.assign((std::size_t)nLoc*(std::size_t)nE*
                                (std::size_t)Geospace::FluxNumerics::kTerminationCount,0);

  //====================================================================================
  // MPI work scheduling
  //====================================================================================
  // As in the Mode3D cutoff solver, every MPI rank has access to the compact global
  // B/E arrays.  The inter-rank scheduler therefore operates on global
  // observation-location indices rather than on rank-owned mesh subdomains.
  //
  // The work-unit granularity is selected below from OUTPUT_MODE:
  //
  //   TRAJECTORY
  //     Mode3D::Run has already assigned one epoch/location transaction to this rank
  //     and installed a singleton MPI_GLOBAL_COMMUNICATOR.  Threads share that one
  //     immutable rank-local field and divide its energy/direction traces.
  //
  //   POINTS / SHELLS
  //     Every location in the call shares one epoch and compact read-only field.  The
  //     MPI scheduler therefore flattens (location,energy) pairs into global work IDs.
  //     This is finer than assigning a whole location to a rank and remains useful
  //     when there are fewer locations than MPI processes.
  //
  // DYNAMIC uses an MPI one-sided atomic queue. MPI is called only by each rank's main
  // thread, never by std::thread/OpenMP workers, so MPI_THREAD_MULTIPLE is unnecessary.
  // BLOCK_CYCLIC and STATIC remain deterministic fallback schedulers for validation.
  //
  // Batched engine: in batch mode the same schedulers feed fibers (see
  // runDensityBatchStream below).  Fibers run on the main thread, so the DYNAMIC chunk
  // fetches they perform are still main-thread MPI calls.
  //====================================================================================
  const Geospace::Mode3D::MpiScheduler mpiScheduler =
      Geospace::Mode3D::ResolveMpiScheduler(prm,"Mode3D density/flux");

  std::vector<int> rankWorkList;
  if (mpiScheduler == Geospace::Mode3D::MpiScheduler::BLOCK_CYCLIC) {
    rankWorkList.reserve((std::size_t)((nLoc + mpiSize - 1) / mpiSize));
    for (int loc=mpiRank; loc<nLoc; loc+=mpiSize) rankWorkList.push_back(loc);
  }
  else if (mpiScheduler == Geospace::Mode3D::MpiScheduler::STATIC) {
    const int begin = (int)((static_cast<long long>(nLoc) * mpiRank) / mpiSize);
    const int end   = (int)((static_cast<long long>(nLoc) * (mpiRank+1)) / mpiSize);
    rankWorkList.reserve((std::size_t)std::max(0,end-begin));
    for (int loc=begin; loc<end; ++loc) rankWorkList.push_back(loc);
  }
  const int nLocalStatic = (int)rankWorkList.size();

  const bool doAnisotropic = (GeospaceUtil::ToUpper(prm.densitySpectrum.boundaryMode) == "ANISOTROPIC");
  const double maximumAllowedWeight=MaximumAllowedWeight_(prm,doAnisotropic);
  const double qabs_C = std::fabs(prm.species.charge_e)*QE;
  const double m0_kg  = prm.species.mass_amu*AMU;

  const DensityParallelBackend_ densityBackend = ResolveDensityParallelBackend_(prm);
  const int densityThreadCount = ResolveDensityThreadCount_(prm,densityBackend);
  const long long mpiDynamicChunk = Geospace::Mode3D::ResolveMpiDynamicChunk(
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

  // Batched engine: batch-engine selection for this call.  The decision depends only on
  // AMPS_EARTH_ENGINE, which is identical on every rank of a run, so all ranks take the
  // same branch below and enter the same sequence of collectives.  The engine's host
  // executor uses the thread count of the selected legacy backend (see
  // DensityBatchHostThreads_); the backend itself is otherwise not used in batch mode.
  const bool densityBatchEngine = Geospace::GPU::DriverBatchEnabled();
  const int densityBatchHostThreads =
      DensityBatchHostThreads_(densityBackend,densityThreadCount);

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
  const std::string outputMode = GeospaceUtil::ToUpper(prm.output.mode);
  const bool isPoints = (outputMode == "POINTS" || outputMode == "TRAJECTORY");
  const bool isShells = (outputMode == "SHELLS");
  const bool sharedEpochFineTasks =
      (outputMode == "POINTS" || outputMode == "SHELLS") && mpiSize > 1;
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
    if (mpiScheduler == Geospace::Mode3D::MpiScheduler::DYNAMIC) {
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

  // Resolve coordinates once before worker threads start.  SHELLS conversion can call
  // SPICE, whose legacy coordinate-transform cache is protected but need not be paid
  // once per energy.  The resulting Cartesian coordinates are immutable and safe for
  // all MPI-scheduled energy tasks and rank-local workers.
  std::vector<V3> locationCoordinates((std::size_t)nLoc);
  for (int loc=0;loc<nLoc;++loc)
    locationCoordinates[(std::size_t)loc]=
        LocationByIndex_m_(prm,loc,nLon,nLat,res_deg,nPtsShell);

  auto storeEnergyAccumulator = [&](int loc,int ie,
                                    const Geospace::FluxNumerics::AccessAccumulator& block) {
    const Geospace::FluxNumerics::AccessEstimate estimate=
        Geospace::FluxNumerics::ResolveAccess(block,maximumAllowedWeight);
    const std::size_t flat=(std::size_t)loc*(std::size_t)nE+(std::size_t)ie;
    local.T_flat[flat]=estimate.nominal;
    local.T_lower_flat[flat]=estimate.lower;
    local.T_upper_flat[flat]=estimate.upper;
    local.unresolved_fraction_flat[flat]=estimate.unresolvedFraction;
    local.sampled_flat[flat]=block.sampled;
    local.resolved_flat[flat]=block.resolved;
    local.allowed_flat[flat]=block.allowed;
    local.retried_flat[flat]=block.retried;
    for (int it=0;it<Geospace::FluxNumerics::kTerminationCount;++it)
      local.termination_flat[
          flat*(std::size_t)Geospace::FluxNumerics::kTerminationCount+(std::size_t)it]=
          block.terminationCounts[(std::size_t)it];
  };

  auto computeOneLocation = [&](int loc) {
    const V3 x0_m = locationCoordinates[(std::size_t)loc];
    std::vector<double> T((std::size_t)nE,0.0),TLower((std::size_t)nE,0.0),
                        TUpper((std::size_t)nE,0.0);
    std::vector<Geospace::FluxNumerics::AccessAccumulator> blocks((std::size_t)nE);

    // The OpenMP loop is inside a C++ lambda.  With default(none), GCC treats
    // variables captured by the enclosing lambda as closure members and may not
    // accept them in the OpenMP data-sharing clauses by their original names.
    // Use local pointer/value aliases declared inside this lambda body, and
    // reference only those aliases inside the parallel region.
    const int nE_local = nE;
    const std::vector<double>* E_MeV_ptr = &E_MeV;
    const GeospaceUtil::AmpsParam* prm_ptr = &prm;
    const std::vector<V3>* dirsUse_ptr = &dirsUse;
    const GeospaceUtil::AnisotropyParam* aniso_ptr = &prm.anisotropy;
    const bool doAnisotropic_local = doAnisotropic;
    const double qabs_C_local = qabs_C;
    const double m0_kg_local = m0_kg;
    // Capture only immutable unit-conversion metadata before the parallel region.
    // gSpectrum's temporal selection is advanced outside worker regions; sharing this
    // value makes that separation explicit and satisfies OpenMP default(none).
    const Geospace::BoundaryProducts::SpectrumUnits spectrumUnits=::gSpectrum.Units();
    const double maxTraceTime_s = (prm.densitySpectrum.maxTrajTime_s > 0.0)
                                ? prm.densitySpectrum.maxTrajTime_s
                                : -1.0;

    if (Geospace::GPU::Fiber::InFiber() && nE_local > 1) {
      //--------------------------------------------------------------------------------
      // Batch mode: mirror of the OpenMP energy loop below.
      //
      // The calling location fiber spawns one child fiber per energy; each traces its
      // direction set through ComputeT_atEnergy_Mode3D_ (which fans the directions out
      // again) and assigns only blocks[ie].  The body is the original loop body with
      // the same inputs (rigidity from the same SpectrumUnits copy and constants), and
      // the per-energy accumulators are consumed by the unchanged serial loop below in
      // increasing energy index, so T/TLower/TUpper, the stored counters and the
      // integrated products are bitwise identical to the LEGACY result.
      //--------------------------------------------------------------------------------
      Geospace::GPU::FiberLoop(0,nE_local,[&](long long index,int) {
        const int ie=static_cast<int>(index);
        const double Ej=spectrumUnits.ParticleEnergyJ(
            (*E_MeV_ptr)[(std::size_t)ie]);
        const double Rgv = RigidityFromEnergy_GV_(Ej,qabs_C_local,m0_kg_local);
        blocks[(std::size_t)ie] = ComputeT_atEnergy_Mode3D_(*prm_ptr,x0_m,Rgv,*dirsUse_ptr,
                                                            maxTraceTime_s,doAnisotropic_local,
                                                            *aniso_ptr);
      });
    }
    else {
      // Original (legacy) energy loop (LEGACY engine, or a single energy).
#ifdef _OPENMP
#pragma omp parallel for default(none) \
  shared(blocks,x0_m,E_MeV_ptr,prm_ptr,dirsUse_ptr,aniso_ptr,spectrumUnits) \
  firstprivate(nE_local,doAnisotropic_local,qabs_C_local,m0_kg_local,maxTraceTime_s) \
  if(!InsideDirectDensityWorker_() && nE_local > 1) schedule(dynamic)
#endif
    for (int ie=0; ie<nE_local; ++ie) {
      // E_MeV is a spectrum coordinate.  Per-nucleon inputs must be converted to
      // total particle kinetic energy before the trajectory rigidity is calculated.
      const double Ej=spectrumUnits.ParticleEnergyJ(
          (*E_MeV_ptr)[(std::size_t)ie]);
      const double Rgv = RigidityFromEnergy_GV_(Ej,qabs_C_local,m0_kg_local);
      blocks[(std::size_t)ie] = ComputeT_atEnergy_Mode3D_(*prm_ptr,x0_m,Rgv,*dirsUse_ptr,
                                                          maxTraceTime_s,doAnisotropic_local,
                                                          *aniso_ptr);
    }
    }

    for (int ie=0; ie<nE; ++ie) {
      const Geospace::FluxNumerics::AccessEstimate estimate=
          Geospace::FluxNumerics::ResolveAccess(blocks[(std::size_t)ie],maximumAllowedWeight);
      T[(std::size_t)ie]=estimate.nominal;
      TLower[(std::size_t)ie]=estimate.lower;
      TUpper[(std::size_t)ie]=estimate.upper;
      storeEnergyAccumulator(loc,ie,blocks[(std::size_t)ie]);
    }

    // All derived quantities must come from this one product set.  Computing density,
    // integral flux, and channels independently allowed unit or clipping behavior to
    // drift between the gridless, Mode3D, and SWMF-coupled paths.
    const Geospace::BoundaryProducts::ProductSet products=
        EvaluateProductSet_(prm,E_MeV,T,TLower,TUpper);
    local.density_m3[(std::size_t)loc]=products.numberDensity_m3.nominal;
    local.density_lower_m3[(std::size_t)loc]=products.numberDensity_m3.lower;
    local.density_upper_m3[(std::size_t)loc]=products.numberDensity_m3.upper;
    local.flux_total_m2s1[(std::size_t)loc]=products.omnidirectionalFlux_m2_s.nominal;
    local.flux_total_lower_m2s1[(std::size_t)loc]=products.omnidirectionalFlux_m2_s.lower;
    local.flux_total_upper_m2s1[(std::size_t)loc]=products.omnidirectionalFlux_m2_s.upper;
    local.flux_planar_m2s1[(std::size_t)loc]=products.oneWayPlanarFlux_m2_s.nominal;
    local.flux_planar_lower_m2s1[(std::size_t)loc]=products.oneWayPlanarFlux_m2_s.lower;
    local.flux_planar_upper_m2s1[(std::size_t)loc]=products.oneWayPlanarFlux_m2_s.upper;
    for (int ic=0; ic<nCh; ++ic) {
      const std::size_t flat=(std::size_t)ic*(std::size_t)nLoc+(std::size_t)loc;
      local.flux_ch_flat[flat]=products.channelFlux_m2_s[(std::size_t)ic].nominal;
      local.flux_ch_lower_flat[flat]=products.channelFlux_m2_s[(std::size_t)ic].lower;
      local.flux_ch_upper_flat[flat]=products.channelFlux_m2_s[(std::size_t)ic].upper;
    }
    for (int id=0;id<nDetector;++id) {
      const std::size_t flat=(std::size_t)id*(std::size_t)nLoc+(std::size_t)loc;
      local.detector_rate_flat[flat]=products.detectorRate_s[(std::size_t)id].nominal;
      local.detector_rate_lower_flat[flat]=products.detectorRate_s[(std::size_t)id].lower;
      local.detector_rate_upper_flat[flat]=products.detectorRate_s[(std::size_t)id].upper;
    }
  };

  // Compute exactly one flattened shared-epoch work item.  A work item owns one
  // (location,energy) result slot, including its termination histogram, so separate
  // MPI ranks and threads never write the same slot.  Directions remain an internal
  // unit of ComputeT_atEnergy_Mode3D_: when this function is invoked inside an OpenMP
  // or direct-thread worker, the nested-parallel guard evaluates them serially; when
  // called outside a worker, that kernel may use its own OpenMP direction team.
  const Geospace::BoundaryProducts::SpectrumUnits sharedEpochSpectrumUnits=
      ::gSpectrum.Units();
  auto computeOneSharedEpochEnergy = [&](long long taskId) {
    const int loc=static_cast<int>(taskId/static_cast<long long>(nE));
    const int ie=static_cast<int>(taskId%static_cast<long long>(nE));
    if (loc<0 || loc>=nLoc || ie<0 || ie>=nE)
      throw std::runtime_error("Mode3D density decoded an invalid shared-epoch task ID.");

    const double Ej=sharedEpochSpectrumUnits.ParticleEnergyJ(E_MeV[(std::size_t)ie]);
    const double Rgv=RigidityFromEnergy_GV_(Ej,qabs_C,m0_kg);
    const double maxTraceTime_s=(prm.densitySpectrum.maxTrajTime_s>0.0)
        ? prm.densitySpectrum.maxTrajTime_s : -1.0;
    const Geospace::FluxNumerics::AccessAccumulator block=
        ComputeT_atEnergy_Mode3D_(prm,locationCoordinates[(std::size_t)loc],Rgv,
                                  dirsUse,maxTraceTime_s,doAnisotropic,prm.anisotropy);
    storeEnergyAccumulator(loc,ie,block);
  };

  // Batched engine: the four rank-local worker-pool lambdas below (computeGlobalRange,
  // computeWorkListRange, computeSharedEpochTaskRange, computeSharedEpochTaskList) are
  // the LEGACY drivers.  In batch mode they are not called; runDensityBatchStream runs
  // the same per-item functions (computeOneLocation / computeOneSharedEpochEnergy) on
  // fibers instead.
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

  // Execute a contiguous range of flattened (location,energy) tasks.  The OpenMP and
  // direct-thread variants deliberately parallelize this outer shared-epoch task
  // space.  ComputeT_atEnergy_Mode3D_ detects the enclosing worker region and avoids
  // creating a nested team for directions, preventing oversubscription.
  auto computeSharedEpochTaskRange = [&](long long begin,long long end) {
    if (end<=begin) return;

    if (densityBackend==DensityParallelBackend_::THREADS && densityThreadCount>1) {
      std::atomic<long long> nextTask(begin);
      const long long nWork=end-begin;
      const int nWorkers=std::max(
          1,std::min(densityThreadCount,static_cast<int>(
              std::min<long long>(nWork,std::numeric_limits<int>::max()))));
      std::vector<std::thread> workers;
      workers.reserve((std::size_t)nWorkers);
      for (int iw=0;iw<nWorkers;++iw) {
        workers.emplace_back([&]() {
          DirectDensityWorkerScope_ workerScope;
          for (;;) {
            const long long taskId=nextTask.fetch_add(1,std::memory_order_relaxed);
            if (taskId>=end) break;
            computeOneSharedEpochEnergy(taskId);
          }
        });
      }
      for (std::thread& worker:workers) worker.join();
      return;
    }

    if (densityBackend==DensityParallelBackend_::OPENMP && densityThreadCount>1) {
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,1)
      for (long long taskId=begin;taskId<end;++taskId)
        computeOneSharedEpochEnergy(taskId);
      return;
#endif
    }

    DirectDensityWorkerScope_ serialScope;
    for (long long taskId=begin;taskId<end;++taskId)
      computeOneSharedEpochEnergy(taskId);
  };

  auto computeSharedEpochTaskList = [&](const std::vector<long long>& work) {
    if (work.empty()) return;

    if (densityBackend==DensityParallelBackend_::THREADS && densityThreadCount>1) {
      std::atomic<std::size_t> nextIndex(0);
      const int nWorkers=std::max(
          1,std::min(densityThreadCount,static_cast<int>(work.size())));
      std::vector<std::thread> workers;
      workers.reserve((std::size_t)nWorkers);
      for (int iw=0;iw<nWorkers;++iw) {
        workers.emplace_back([&]() {
          DirectDensityWorkerScope_ workerScope;
          for (;;) {
            const std::size_t index=nextIndex.fetch_add(1,std::memory_order_relaxed);
            if (index>=work.size()) break;
            computeOneSharedEpochEnergy(work[index]);
          }
        });
      }
      for (std::thread& worker:workers) worker.join();
      return;
    }

    if (densityBackend==DensityParallelBackend_::OPENMP && densityThreadCount>1) {
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,1)
      for (long long index=0;index<static_cast<long long>(work.size());++index)
        computeOneSharedEpochEnergy(work[(std::size_t)index]);
      return;
#endif
    }

    DirectDensityWorkerScope_ serialScope;
    for (const long long taskId:work) computeOneSharedEpochEnergy(taskId);
  };

  //------------------------------------------------------------------------------------
  // Batch mode: rank-local fiber driver shared by the location and the
  // (location,energy) schedulers below.  It replaces the std::thread / OpenMP / serial
  // worker pools of computeGlobalRange, computeWorkListRange,
  // computeSharedEpochTaskRange and computeSharedEpochTaskList.
  //
  //   * `nFibers` worker fibers pull ids from `stream` (this rank's ids, in the original
  //     scheduler order) and run body(id), which is the unchanged per-location or
  //     per-(location,energy) work function.  The trajectories of all fibers and of
  //     their nested energy/direction fibers are integrated together by the engine.
  //   * The engine's host executor gets densityBatchHostThreads threads, i.e. the
  //     thread count the legacy backend would have used for these workers.
  //   * body() writes only its own result slots and the plain counters it owns; all
  //     fibers run on this (main) thread, so no atomics or locks are needed.
  //   * service() runs on the main thread between engine pool advances (never inside
  //     a fiber), which is where live progress is printed.
  //   * If body() throws, the stream is stopped so the other fibers finish only their
  //     current item, and RunDriverFibers rethrows the first exception after all fibers
  //     have stopped (the LEGACY std::thread/OpenMP workers would have terminated the
  //     process instead; the serial LEGACY path propagated the exception the same way).
  //------------------------------------------------------------------------------------
  auto runDensityBatchStream=[&](DensityBatchWorkStream_& stream,
                                 long long nItemsUpperBound,
                                 const std::function<void(long long)>& body,
                                 const std::function<void()>& service,
                                 const char* label) {
    if (nItemsUpperBound<=0) return;
    const int nFibers=Geospace::GPU::DriverFiberCount(nItemsUpperBound,densityBatchHostThreads);
    Geospace::GPU::RunDriverFibers(nFibers,densityBatchHostThreads,[&](int) {
      for (;;) {
        const long long id=stream.Next();
        if (id<0) break;
        try {
          body(id);
        }
        catch (...) {
          stream.Stop();
          throw;
        }
      }
    },service,label);
  };

  maybePrintProgress(0,0,locDonePerShellGlobal,true);

  if (sharedEpochFineTasks) {
    // POINTS/SHELLS share one immutable field, so ranks may safely own different
    // energies of the same location.  The global result reduction below reassembles
    // complete access curves before density and flux integration.
    const long long nEnergyTasks=static_cast<long long>(nLoc)*
                                 static_cast<long long>(nE);
    long long completedEnergyTasksLocal=0;

    if (densityBatchEngine) {
      //--------------------------------------------------------------------------------
      // Batch mode: mirror of the DYNAMIC / BLOCK_CYCLIC / STATIC branches
      // below.
      //
      // The rank computes exactly the task ids the original scheduler assigns to it:
      // the same BLOCK_CYCLIC/STATIC lists, or DYNAMIC chunks of the same size from an
      // RMA counter constructed collectively on the same communicator with the same
      // label.  The only difference for DYNAMIC is that fibers fetch the next chunk as
      // soon as they need work instead of after the whole previous chunk has finished,
      // so the trajectory pool is not drained at chunk boundaries.  Each task runs the
      // unchanged computeOneSharedEpochEnergy(), which assigns the (location,energy)
      // slot it owns; completedEnergyTasksLocal counts finished tasks exactly like the
      // original branches, so the completeness collective below is unchanged.
      //
      // Progress: the original prints only the initial and the final (collective)
      // state on this path because rank-local counts are not global; batch mode keeps
      // that and installs no service hook.
      //--------------------------------------------------------------------------------
      std::unique_ptr<Geospace::Mode3D::DynamicMpiLocationScheduler> dynamicScheduler;
      std::vector<long long> work;
      if (mpiScheduler==Geospace::Mode3D::MpiScheduler::DYNAMIC) {
        const long long fineChunk=Geospace::Mode3D::ResolveMpiDynamicChunk(
            prm,densityThreadCount,nEnergyTasks);
        dynamicScheduler.reset(new Geospace::Mode3D::DynamicMpiLocationScheduler(
            MPI_GLOBAL_COMMUNICATOR,nEnergyTasks,fineChunk,
            "Mode3D density shared-epoch (location,energy) scheduler"));
      }
      else if (mpiScheduler==Geospace::Mode3D::MpiScheduler::BLOCK_CYCLIC) {
        for (long long taskId=mpiRank;taskId<nEnergyTasks;taskId+=mpiSize)
          work.push_back(taskId);
      }
      else {
        const long long begin=(nEnergyTasks*mpiRank)/mpiSize;
        const long long end=(nEnergyTasks*(mpiRank+1))/mpiSize;
        work.reserve((std::size_t)std::max(0LL,end-begin));
        for (long long taskId=begin;taskId<end;++taskId) work.push_back(taskId);
      }
      DensityBatchWorkStream_ stream=dynamicScheduler
          ? DensityBatchWorkStream_(dynamicScheduler.get(),nEnergyTasks)
          : DensityBatchWorkStream_(work);
      runDensityBatchStream(
          stream,
          dynamicScheduler ? nEnergyTasks : static_cast<long long>(work.size()),
          [&](long long taskId) {
            computeOneSharedEpochEnergy(taskId);
            ++completedEnergyTasksLocal;
          },
          nullptr,"Mode3D density shared-epoch (location,energy)");
    }
    else if (mpiScheduler==Geospace::Mode3D::MpiScheduler::DYNAMIC) {
      const long long fineChunk=Geospace::Mode3D::ResolveMpiDynamicChunk(
          prm,densityThreadCount,nEnergyTasks);
      Geospace::Mode3D::DynamicMpiLocationScheduler scheduler(
          MPI_GLOBAL_COMMUNICATOR,nEnergyTasks,fineChunk,
          "Mode3D density shared-epoch (location,energy) scheduler");
      for (;;) {
        const long long begin=scheduler.FetchNextChunkStart();
        if (begin>=nEnergyTasks) break;
        const long long end=std::min(begin+scheduler.ChunkSize(),nEnergyTasks);
        computeSharedEpochTaskRange(begin,end);
        completedEnergyTasksLocal+=end-begin;
      }
    }
    else {
      std::vector<long long> work;
      if (mpiScheduler==Geospace::Mode3D::MpiScheduler::BLOCK_CYCLIC) {
        for (long long taskId=mpiRank;taskId<nEnergyTasks;taskId+=mpiSize)
          work.push_back(taskId);
      }
      else {
        const long long begin=(nEnergyTasks*mpiRank)/mpiSize;
        const long long end=(nEnergyTasks*(mpiRank+1))/mpiSize;
        work.reserve((std::size_t)std::max(0LL,end-begin));
        for (long long taskId=begin;taskId<end;++taskId) work.push_back(taskId);
      }
      computeSharedEpochTaskList(work);
      completedEnergyTasksLocal=static_cast<long long>(work.size());
    }

    long long completedEnergyTasksGlobal=0;
    MPI_Allreduce(&completedEnergyTasksLocal,&completedEnergyTasksGlobal,1,
                  MPI_LONG_LONG,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
    if (completedEnergyTasksGlobal!=nEnergyTasks)
      throw std::runtime_error(
          "Mode3D density shared-epoch scheduler did not complete every (location,energy) task.");

    // A location is reported complete only after the collective proves that all of its
    // energy tasks are present. This keeps the existing progress schema honest rather
    // than claiming completion when merely one energy was assigned.
    doneLocationsGlobal=totalLocationsGlobal;
    doneTasksGlobal=totalTasksGlobal;
    if (isShells && nShells>0)
      for (int s=0;s<nShells;++s)
        locDonePerShellGlobal[(std::size_t)s]=nPtsShell;
    maybePrintProgress(doneLocationsGlobal,doneTasksGlobal,locDonePerShellGlobal,true);
  }
  else if (densityBatchEngine) {
    //----------------------------------------------------------------------------------
    // Batch mode, whole-location work units (single-rank POINTS/SHELLS and
    // every TRAJECTORY call): mirror of the DYNAMIC branch and of the batched
    // STATIC/BLOCK_CYCLIC branch below.
    //
    // Work assignment: DYNAMIC fetches chunks of mpiDynamicChunk locations from an RMA
    // counter built exactly like the original one (same communicator, size, chunk and
    // label); STATIC/BLOCK_CYCLIC process rankWorkList in list order.  Each location
    // runs the unchanged computeOneLocation() (energy and direction loops fan out as
    // nested fibers) and is then accounted with accountCompletedGlobalRange(loc,loc+1),
    // which adds exactly what the original chunk/batch accounting adds for it.
    //
    // Progress: the whole rank-local work list is one fiber batch, as the original
    // THREADS backend already does (nProgressBatches=1), because splitting it into
    // ~200 synchronized batches would drain the trajectory pool at every batch
    // boundary.  Consequently one set of progress collectives follows the batch, on
    // every rank (the batch decision is the same on all ranks), so the collective
    // sequence stays matched.  When this call runs on a single rank the rank-local
    // counters are the global ones and the service hook prints live, throttled
    // progress from the main thread between pool advances; with several ranks only the
    // final collective progress is printed, as in the original DYNAMIC branch.
    //----------------------------------------------------------------------------------
    std::unique_ptr<Geospace::Mode3D::DynamicMpiLocationScheduler> dynamicScheduler;
    if (mpiScheduler == Geospace::Mode3D::MpiScheduler::DYNAMIC) {
      dynamicScheduler.reset(new Geospace::Mode3D::DynamicMpiLocationScheduler(
          MPI_GLOBAL_COMMUNICATOR,
          static_cast<long long>(nLoc),
          mpiDynamicChunk,
          "Mode3D density/flux"));
    }
    DensityBatchWorkStream_ stream=dynamicScheduler
        ? DensityBatchWorkStream_(dynamicScheduler.get(),static_cast<long long>(nLoc))
        : DensityBatchWorkStream_(rankWorkList);
    auto progressService=[&]() {
      if (mpiSize!=1) return;
      maybePrintProgress(doneLocationsLocal,doneTasksLocal,locDonePerShellLocal,false);
    };
    runDensityBatchStream(
        stream,
        dynamicScheduler ? static_cast<long long>(nLoc) : static_cast<long long>(nLocalStatic),
        [&](long long id) {
          const int loc=static_cast<int>(id);
          computeOneLocation(loc);
          accountCompletedGlobalRange(loc,loc+1);
        },
        progressService,"Mode3D density/flux");

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
  else if (mpiScheduler == Geospace::Mode3D::MpiScheduler::DYNAMIC) {
    // Two-level dynamic scheduling: each rank dynamically fetches chunks from the MPI
    // RMA counter, and the selected intra-rank backend dynamically computes locations
    // within that chunk.  No progress collectives are placed inside the work loop, so
    // ranks are never forced to wait at a synchronization point between chunks.
    Geospace::Mode3D::DynamicMpiLocationScheduler scheduler(
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
  global.flux_planar_m2s1.assign((std::size_t)nLoc,0.0);
  global.flux_planar_lower_m2s1.assign((std::size_t)nLoc,0.0);
  global.flux_planar_upper_m2s1.assign((std::size_t)nLoc,0.0);
  global.T_flat.assign((std::size_t)nLoc*(std::size_t)nE,0.0);
  global.T_lower_flat.assign((std::size_t)nLoc*(std::size_t)nE,0.0);
  global.T_upper_flat.assign((std::size_t)nLoc*(std::size_t)nE,0.0);
  global.unresolved_fraction_flat.assign((std::size_t)nLoc*(std::size_t)nE,0.0);
  global.flux_ch_flat.assign((std::size_t)nCh*(std::size_t)nLoc,0.0);
  global.flux_ch_lower_flat.assign((std::size_t)nCh*(std::size_t)nLoc,0.0);
  global.flux_ch_upper_flat.assign((std::size_t)nCh*(std::size_t)nLoc,0.0);
  global.detector_rate_flat.assign((std::size_t)nDetector*(std::size_t)nLoc,0.0);
  global.detector_rate_lower_flat.assign((std::size_t)nDetector*(std::size_t)nLoc,0.0);
  global.detector_rate_upper_flat.assign((std::size_t)nDetector*(std::size_t)nLoc,0.0);
  global.sampled_flat.assign((std::size_t)nLoc*(std::size_t)nE,0);
  global.resolved_flat.assign((std::size_t)nLoc*(std::size_t)nE,0);
  global.allowed_flat.assign((std::size_t)nLoc*(std::size_t)nE,0);
  global.retried_flat.assign((std::size_t)nLoc*(std::size_t)nE,0);
  global.termination_flat.assign((std::size_t)nLoc*(std::size_t)nE*
                                 (std::size_t)Geospace::FluxNumerics::kTerminationCount,0);

  // Whole-location scheduling computes integrated products on the owning rank and
  // therefore retains the historical reductions. Shared-epoch fine scheduling first
  // reduces the disjoint energy rows below, then integrates each complete access curve
  // identically on every rank. Integrating before the energy reduction would combine
  // incomplete curves and is scientifically invalid.
  if (!sharedEpochFineTasks) {
    MPI_Allreduce(local.density_m3.data(), global.density_m3.data(), nLoc, MPI_DOUBLE, MPI_SUM, MPI_GLOBAL_COMMUNICATOR);
    MPI_Allreduce(local.density_lower_m3.data(), global.density_lower_m3.data(), nLoc, MPI_DOUBLE, MPI_SUM, MPI_GLOBAL_COMMUNICATOR);
    MPI_Allreduce(local.density_upper_m3.data(), global.density_upper_m3.data(), nLoc, MPI_DOUBLE, MPI_SUM, MPI_GLOBAL_COMMUNICATOR);
    MPI_Allreduce(local.flux_total_m2s1.data(), global.flux_total_m2s1.data(), nLoc, MPI_DOUBLE, MPI_SUM, MPI_GLOBAL_COMMUNICATOR);
    MPI_Allreduce(local.flux_total_lower_m2s1.data(), global.flux_total_lower_m2s1.data(), nLoc, MPI_DOUBLE, MPI_SUM, MPI_GLOBAL_COMMUNICATOR);
    MPI_Allreduce(local.flux_total_upper_m2s1.data(), global.flux_total_upper_m2s1.data(), nLoc, MPI_DOUBLE, MPI_SUM, MPI_GLOBAL_COMMUNICATOR);
    MPI_Allreduce(local.flux_planar_m2s1.data(), global.flux_planar_m2s1.data(), nLoc, MPI_DOUBLE, MPI_SUM, MPI_GLOBAL_COMMUNICATOR);
    MPI_Allreduce(local.flux_planar_lower_m2s1.data(), global.flux_planar_lower_m2s1.data(), nLoc, MPI_DOUBLE, MPI_SUM, MPI_GLOBAL_COMMUNICATOR);
    MPI_Allreduce(local.flux_planar_upper_m2s1.data(), global.flux_planar_upper_m2s1.data(), nLoc, MPI_DOUBLE, MPI_SUM, MPI_GLOBAL_COMMUNICATOR);
  }
  MPI_Allreduce(local.T_flat.data(), global.T_flat.data(), nLoc*nE, MPI_DOUBLE, MPI_SUM, MPI_GLOBAL_COMMUNICATOR);
  MPI_Allreduce(local.T_lower_flat.data(), global.T_lower_flat.data(), nLoc*nE, MPI_DOUBLE, MPI_SUM, MPI_GLOBAL_COMMUNICATOR);
  MPI_Allreduce(local.T_upper_flat.data(), global.T_upper_flat.data(), nLoc*nE, MPI_DOUBLE, MPI_SUM, MPI_GLOBAL_COMMUNICATOR);
  MPI_Allreduce(local.unresolved_fraction_flat.data(), global.unresolved_fraction_flat.data(), nLoc*nE, MPI_DOUBLE, MPI_SUM, MPI_GLOBAL_COMMUNICATOR);
  MPI_Allreduce(local.sampled_flat.data(),global.sampled_flat.data(),nLoc*nE,MPI_INT,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
  MPI_Allreduce(local.resolved_flat.data(),global.resolved_flat.data(),nLoc*nE,MPI_INT,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
  MPI_Allreduce(local.allowed_flat.data(),global.allowed_flat.data(),nLoc*nE,MPI_INT,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
  MPI_Allreduce(local.retried_flat.data(),global.retried_flat.data(),nLoc*nE,MPI_INT,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
  MPI_Allreduce(local.termination_flat.data(),global.termination_flat.data(),
                nLoc*nE*Geospace::FluxNumerics::kTerminationCount,MPI_INT,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
  if (!sharedEpochFineTasks && nCh*nLoc > 0) {
    MPI_Allreduce(local.flux_ch_flat.data(), global.flux_ch_flat.data(), nCh*nLoc, MPI_DOUBLE, MPI_SUM, MPI_GLOBAL_COMMUNICATOR);
    MPI_Allreduce(local.flux_ch_lower_flat.data(), global.flux_ch_lower_flat.data(), nCh*nLoc, MPI_DOUBLE, MPI_SUM, MPI_GLOBAL_COMMUNICATOR);
    MPI_Allreduce(local.flux_ch_upper_flat.data(), global.flux_ch_upper_flat.data(), nCh*nLoc, MPI_DOUBLE, MPI_SUM, MPI_GLOBAL_COMMUNICATOR);
  }
  if (!sharedEpochFineTasks && nDetector*nLoc > 0) {
    MPI_Allreduce(local.detector_rate_flat.data(),global.detector_rate_flat.data(),
                  nDetector*nLoc,MPI_DOUBLE,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
    MPI_Allreduce(local.detector_rate_lower_flat.data(),global.detector_rate_lower_flat.data(),
                  nDetector*nLoc,MPI_DOUBLE,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
    MPI_Allreduce(local.detector_rate_upper_flat.data(),global.detector_rate_upper_flat.data(),
                  nDetector*nLoc,MPI_DOUBLE,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
  }

  if (sharedEpochFineTasks) {
    // Every global access row is now complete. Derive density, total/planar flux,
    // channels, and detector responses from that one assembled product set, exactly as
    // the whole-location path does. All ranks execute this deterministic inexpensive
    // integration so later root-only output and collective validation see identical
    // buffers without another large reduction.
    for (int loc=0;loc<nLoc;++loc) {
      std::vector<double> T((std::size_t)nE,0.0);
      std::vector<double> TLower((std::size_t)nE,0.0);
      std::vector<double> TUpper((std::size_t)nE,0.0);
      for (int ie=0;ie<nE;++ie) {
        const std::size_t flat=(std::size_t)loc*(std::size_t)nE+(std::size_t)ie;
        T[(std::size_t)ie]=global.T_flat[flat];
        TLower[(std::size_t)ie]=global.T_lower_flat[flat];
        TUpper[(std::size_t)ie]=global.T_upper_flat[flat];
      }

      const Geospace::BoundaryProducts::ProductSet products=
          EvaluateProductSet_(prm,E_MeV,T,TLower,TUpper);
      global.density_m3[(std::size_t)loc]=products.numberDensity_m3.nominal;
      global.density_lower_m3[(std::size_t)loc]=products.numberDensity_m3.lower;
      global.density_upper_m3[(std::size_t)loc]=products.numberDensity_m3.upper;
      global.flux_total_m2s1[(std::size_t)loc]=products.omnidirectionalFlux_m2_s.nominal;
      global.flux_total_lower_m2s1[(std::size_t)loc]=products.omnidirectionalFlux_m2_s.lower;
      global.flux_total_upper_m2s1[(std::size_t)loc]=products.omnidirectionalFlux_m2_s.upper;
      global.flux_planar_m2s1[(std::size_t)loc]=products.oneWayPlanarFlux_m2_s.nominal;
      global.flux_planar_lower_m2s1[(std::size_t)loc]=products.oneWayPlanarFlux_m2_s.lower;
      global.flux_planar_upper_m2s1[(std::size_t)loc]=products.oneWayPlanarFlux_m2_s.upper;
      for (int ic=0;ic<nCh;++ic) {
        const std::size_t flat=(std::size_t)ic*(std::size_t)nLoc+(std::size_t)loc;
        global.flux_ch_flat[flat]=products.channelFlux_m2_s[(std::size_t)ic].nominal;
        global.flux_ch_lower_flat[flat]=products.channelFlux_m2_s[(std::size_t)ic].lower;
        global.flux_ch_upper_flat[flat]=products.channelFlux_m2_s[(std::size_t)ic].upper;
      }
      for (int id=0;id<nDetector;++id) {
        const std::size_t flat=(std::size_t)id*(std::size_t)nLoc+(std::size_t)loc;
        global.detector_rate_flat[flat]=products.detectorRate_s[(std::size_t)id].nominal;
        global.detector_rate_lower_flat[flat]=products.detectorRate_s[(std::size_t)id].lower;
        global.detector_rate_upper_flat[flat]=products.detectorRate_s[(std::size_t)id].upper;
      }
    }
  }

  return global;
}

//--------------------------------------------------------------------------------------
// Tecplot writers
//--------------------------------------------------------------------------------------
template<class Stream>
static void WriteStep6ProductMetadata_(Stream& out) {
  const char* basis=(::gSpectrum.EnergyCoordinateBasis()==
      Geospace::BoundaryProducts::EnergyBasis::PerNucleon)
      ? "PER_NUCLEON" : "PER_PARTICLE";
  out << "AUXDATA STEP6_CHARACTERISTIC_MAPPING=\"STATIC_MAGNETIC\"\n"
      << "AUXDATA SPECTRUM_ENERGY_BASIS=\"" << basis << "\"\n"
      << "AUXDATA SPECTRUM_MASS_NUMBER=\"" << ::gSpectrum.MassNumber() << "\"\n"
      << "AUXDATA SPECTRUM_INTENSITY_UNIT=\"" << ::gSpectrum.IntensityUnitLabel() << "\"\n"
      << "AUXDATA SPECTRUM_RELATIVE_UNCERTAINTY=\""
      << ::gSpectrum.RelativeUncertainty() << "\"\n"
      << "AUXDATA SPECTRUM_TEMPORAL_STATUS=\""
      << Geospace::BoundaryProducts::TemporalStatusName(::gSpectrum.LastTemporalStatus())
      << "\"\n"
      << "AUXDATA SPECTRUM_TEMPORAL_GAP=\""
      << (::gSpectrum.LastTemporalSelectionCrossedGap() ? 1 : 0) << "\"\n"
      << "AUXDATA SPECTRUM_TEMPORAL_FRACTION=\""
      << ::gSpectrum.LastTemporalInterpolationFraction() << "\"\n"
      << "AUXDATA PLANAR_FLUX_CONVENTION=\"ISOTROPIC_EQUIVALENT_PI_J\"\n";
}

static std::string DensityTecplotAuxValue_(std::string value) {
  std::replace(value.begin(),value.end(),'"','\'');
  std::replace(value.begin(),value.end(),'\n',' ');
  std::replace(value.begin(),value.end(),'\r',' ');
  return value;
}

template<class Stream>
static void WriteStep11ProductMetadata_(Stream& out,
                                        const GeospaceUtil::AmpsParam& prm) {
  namespace CP=Geospace::SWMFCoupledProducts;
  const Geospace::Field::SnapshotMetadata& metadata=
      Geospace::Mode3D::GlobalMagneticField::CurrentSnapshotMetadata();
  const CP::ProductControl& control=gLastDensityFluxRunSummary.control;

  // These fields deliberately duplicate the coupled manifest's identity in every
  // numeric artifact.  A file copied away from its manifest remains self-describing,
  // and the strict replay comparator can reject a stale field, spectrum table,
  // channel definition, or instrument response before looking at numeric rows.
  out << "AUXDATA COMMON_PHYSICS_TAG=\""
      << Geospace::CommonPhysicsRelease::kTag << "\"\n"
      << "AUXDATA COMMON_PHYSICS_SCOPE=\""
      << Geospace::CommonPhysicsRelease::kScope << "\"\n"
      << "AUXDATA PHASE_1_INTERPRETATION=\"INSTANTANEOUS_QUASI_STATIC\"\n"
      << "AUXDATA SNAPSHOT_ID=\""
      << DensityTecplotAuxValue_(metadata.snapshotId) << "\"\n"
      << "AUXDATA SNAPSHOT_EPOCH_UTC=\""
      << DensityTecplotAuxValue_(metadata.epochUTC) << "\"\n"
      << "AUXDATA SNAPSHOT_MESH_REVISION=\""
      << DensityTecplotAuxValue_(
             Geospace::Mode3D::GlobalMagneticField::CurrentMeshRevision()) << "\"\n"
      << "AUXDATA SNAPSHOT_CONTENT_FINGERPRINT=\""
      << DensityTecplotAuxValue_(
             Geospace::Mode3D::GlobalMagneticField::CurrentContentFingerprint()) << "\"\n"
      << "AUXDATA BOUNDARY_SPECTRUM_EVALUATION_EPOCH_UTC=\""
      << DensityTecplotAuxValue_(
             gLastDensityFluxRunSummary.spectrumEvaluationEpochUTC) << "\"\n"
      << "AUXDATA ACTIVE_SPECTRUM_TABLE_EPOCH_UTC=\""
      << DensityTecplotAuxValue_(
             gLastDensityFluxRunSummary.activeSpectrumTableEpochUTC) << "\"\n"
      << "AUXDATA PRODUCT_CONTROL_FINGERPRINT=\""
      << CP::ProductControlFingerprint(control) << "\"\n"
      << "AUXDATA BOUNDARY_SPECTRUM_FINGERPRINT=\""
      << CP::BoundarySpectrumFingerprint(control) << "\"\n"
      << "AUXDATA CHANNEL_SCHEMA_FINGERPRINT=\""
      << CP::ChannelSchemaFingerprint(control) << "\"\n"
      << "AUXDATA DETECTOR_RESPONSE_FINGERPRINT=\""
      << CP::DetectorResponseFingerprint(control) << "\"\n"
      << "AUXDATA OBSERVATION_STATE_FINGERPRINT=\""
      << CP::ObservationStateFingerprint(control) << "\"\n"
      << "AUXDATA MAXIMUM_UNRESOLVED_FRACTION=\""
      << gLastDensityFluxRunSummary.maximumUnresolvedFraction << "\"\n"
      << "AUXDATA RESPONSE_WEIGHTED_UNRESOLVED_UPPER_BOUND=\""
      << gLastDensityFluxRunSummary.maximumUnresolvedFraction << "\"\n"
      << "AUXDATA UNRESOLVED_TOLERANCE=\""
      << gLastDensityFluxRunSummary.unresolvedTolerance << "\"\n";

  const std::string boundary=GeospaceUtil::ToUpper(prm.domain.boundaryType.empty()
      ? std::string("BOX") : prm.domain.boundaryType);
  out << "AUXDATA OUTER_BOUNDARY_POLICY=\"" << boundary << "\"\n";
  if (boundary=="SHUE") {
    const Geospace::SWMFCoupledAccess::ShueParameters shue=
        Geospace::SWMFCoupledAccess::ResolveShueParameters(
            prm.domain.shueR0Token,prm.domain.shueAlphaToken,
            prm.field.pdyn_nPa,prm.field.imfBz_nT,_EARTH__RADIUS_,
            1000.0*prm.domain.xMin);
    out << "AUXDATA SHUE_R0_RE=\"" << shue.r0_Re << "\"\n"
        << "AUXDATA SHUE_ALPHA=\"" << shue.alpha << "\"\n"
        << "AUXDATA SHUE_TAIL_CAP_X_M=\"" << shue.tailCapX_m << "\"\n";
  }
}

template<class Stream>
static void WriteProductMetadata_(Stream& out,const GeospaceUtil::AmpsParam& prm) {
  WriteStep6ProductMetadata_(out);
  WriteStep11ProductMetadata_(out,prm);
}

static void WritePointOutputs_(const GeospaceUtil::AmpsParam& prm,
                               const std::vector<double>& E_MeV,
                               const DensityResultBuffers& res) {
  const int nLoc = (int)prm.output.points.size();
  const int nE   = (int)E_MeV.size();
  const int nCh  = (int)prm.fluxChannels.size();
  const int nDetector = (int)prm.detectorResponses.size();

  {
    const std::string fileName=DensityOutputFileName_("mode3d_points_density");
    std::ofstream out=OpenDensityArtifact_(fileName);
    out << std::setprecision(std::numeric_limits<double>::max_digits10);
    out << "TITLE=\"Mode3D mesh-field energetic particle density\"\n";
    WriteProductMetadata_(out,prm);
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
    CloseAndRecordDensityArtifact_(out,fileName);
  }

  {
    const std::string fileName=DensityOutputFileName_("mode3d_points_spectrum");
    std::ofstream out=OpenDensityArtifact_(fileName);
    out << std::setprecision(std::numeric_limits<double>::max_digits10);
    out << "TITLE=\"Mode3D mesh-field local energetic particle spectrum\"\n";
    WriteProductMetadata_(out,prm);
    out << "VARIABLES=\"E_MeV\" \"T\" \"T_lower\" \"T_upper\" "
        << "\"unresolved_fraction\" \"N_sampled\" \"N_resolved\" \"N_allowed\" "
        << "\"J_boundary_perMeV\" \"J_local_perMeV\" \"J_local_lower_perMeV\" \"J_local_upper_perMeV\" "
        << "\"J_boundary_lower_perMeV\" \"J_boundary_upper_perMeV\" "
        << "\"J_omni_perMeV\" \"J_omni_lower_perMeV\" \"J_omni_upper_perMeV\" "
        << "\"J_planar_perMeV\" \"J_planar_lower_perMeV\" \"J_planar_upper_perMeV\"\n";
    for (int loc=0; loc<nLoc; ++loc) {
      out << "ZONE T=\"loc_" << std::setw(6) << std::setfill('0') << loc << std::setfill(' ')
          << "\" I=" << nE << " F=POINT\n";
      for (int ie=0; ie<nE; ++ie) {
        const double T = res.T_flat[(std::size_t)loc*(std::size_t)nE + (std::size_t)ie];
        const std::size_t flat=(std::size_t)loc*(std::size_t)nE+(std::size_t)ie;
        const double lower=res.T_lower_flat[flat],upper=res.T_upper_flat[flat];
        const Geospace::BoundaryProducts::Bounds Jb=
            ::gSpectrum.GetSpectrumPerMeVBounds(E_MeV[(std::size_t)ie]);
        const Geospace::BoundaryProducts::Bounds local=
            Geospace::BoundaryProducts::MapBoundaryIntensity(
                Geospace::BoundaryProducts::Bounds(T,lower,upper),Jb,
                Geospace::BoundaryProducts::CharacteristicMapping::StaticMagnetic);
        out << E_MeV[(std::size_t)ie] << " " << T << " " << lower << " " << upper
            << " " << res.unresolved_fraction_flat[flat] << " " << res.sampled_flat[flat]
            << " " << res.resolved_flat[flat] << " " << res.allowed_flat[flat]
            << " " << Jb.nominal << " " << local.nominal
            << " " << local.lower << " " << local.upper
            << " " << Jb.lower << " " << Jb.upper
            << " " << 4.0*Geospace::FluxNumerics::kPi*local.nominal
            << " " << 4.0*Geospace::FluxNumerics::kPi*local.lower
            << " " << 4.0*Geospace::FluxNumerics::kPi*local.upper
            << " " << Geospace::FluxNumerics::kPi*local.nominal
            << " " << Geospace::FluxNumerics::kPi*local.lower
            << " " << Geospace::FluxNumerics::kPi*local.upper << "\n";
      }
    }
    CloseAndRecordDensityArtifact_(out,fileName);
  }

  {
    const std::string fileName=DensityOutputFileName_("mode3d_points_flux");
    std::ofstream out=OpenDensityArtifact_(fileName);
    out << std::setprecision(std::numeric_limits<double>::max_digits10);
    out << "TITLE=\"Mode3D mesh-field omnidirectional integral flux\"\n";
    WriteProductMetadata_(out,prm);
    out << "VARIABLES=\"X_km\" \"Y_km\" \"Z_km\" \"F_tot_m2s1\" \"F_tot_lower_m2s1\" \"F_tot_upper_m2s1\"";
    for (int ic=0; ic<nCh; ++ic) out << " \"F_" << prm.fluxChannels[(std::size_t)ic].name << "_m2s1\""
        << " \"F_" << prm.fluxChannels[(std::size_t)ic].name << "_lower_m2s1\""
        << " \"F_" << prm.fluxChannels[(std::size_t)ic].name << "_upper_m2s1\"";
    out << " \"F_planar_m2s1\" \"F_planar_lower_m2s1\" \"F_planar_upper_m2s1\"";
    for (int id=0;id<nDetector;++id) {
      out << " \"R_" << prm.detectorResponses[(std::size_t)id].name << "_s1\""
          << " \"R_" << prm.detectorResponses[(std::size_t)id].name << "_lower_s1\""
          << " \"R_" << prm.detectorResponses[(std::size_t)id].name << "_upper_s1\"";
    }
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
      out << " " << res.flux_planar_m2s1[(std::size_t)i]
          << " " << res.flux_planar_lower_m2s1[(std::size_t)i]
          << " " << res.flux_planar_upper_m2s1[(std::size_t)i];
      for (int id=0;id<nDetector;++id) {
        const std::size_t flat=(std::size_t)id*(std::size_t)nLoc+(std::size_t)i;
        out << " " << res.detector_rate_flat[flat]
            << " " << res.detector_rate_lower_flat[flat]
            << " " << res.detector_rate_upper_flat[flat];
      }
      out << "\n";
    }
    CloseAndRecordDensityArtifact_(out,fileName);
  }
}

static void WriteShellOutputs_(const GeospaceUtil::AmpsParam& prm,
                               int nLon,int nLat,double res_deg,int nPtsShell,
                               const std::vector<double>& E_MeV,
                               const DensityResultBuffers& res) {
  const int nE = (int)E_MeV.size();
  const int nShells = (int)prm.output.shellAlt_km.size();
  const int nCh     = (int)prm.fluxChannels.size();
  const int nDetector = (int)prm.detectorResponses.size();
  const int nLoc    = nShells*nPtsShell;

  for (int s=0; s<nShells; ++s) {
    const std::string altLabel = FormatEnergyBoundForName_(prm.output.shellAlt_km[(std::size_t)s]);
    const std::string stem = "mode3d_shell_" + altLabel + "km_density_flux";
    const std::string fileName=DensityOutputFileName_(stem.c_str());
    std::ofstream out=OpenDensityArtifact_(fileName);
    out << std::setprecision(std::numeric_limits<double>::max_digits10);

    out << "TITLE=\"Mode3D mesh-field density and flux shell alt="
        << prm.output.shellAlt_km[(std::size_t)s] << " km\"\n";
    WriteProductMetadata_(out,prm);
    out << "VARIABLES=\"Lon_deg\" \"Lat_deg\" \"N_m^-3\" \"N_lower_m^-3\" \"N_upper_m^-3\" "
        << "\"N_cm^-3\" \"F_tot_m2s1\" \"F_tot_lower_m2s1\" \"F_tot_upper_m2s1\" "
        << "\"Rc_lower_GV\" \"Rc_effective_GV\" \"Rc_upper_GV\" \"PenumbraWidth_GV\" \"T_high\"";
    for (int ic=0; ic<nCh; ++ic) out << " \"F_" << prm.fluxChannels[(std::size_t)ic].name << "_m2s1\""
        << " \"F_" << prm.fluxChannels[(std::size_t)ic].name << "_lower_m2s1\""
        << " \"F_" << prm.fluxChannels[(std::size_t)ic].name << "_upper_m2s1\"";
    out << " \"F_planar_m2s1\" \"F_planar_lower_m2s1\" \"F_planar_upper_m2s1\"";
    for (int id=0;id<nDetector;++id) {
      out << " \"R_" << prm.detectorResponses[(std::size_t)id].name << "_s1\""
          << " \"R_" << prm.detectorResponses[(std::size_t)id].name << "_lower_s1\""
          << " \"R_" << prm.detectorResponses[(std::size_t)id].name << "_upper_s1\"";
    }
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
        out << " " << res.flux_planar_m2s1[(std::size_t)loc]
            << " " << res.flux_planar_lower_m2s1[(std::size_t)loc]
            << " " << res.flux_planar_upper_m2s1[(std::size_t)loc];
        for (int id=0;id<nDetector;++id) {
          const std::size_t flat=(std::size_t)id*(std::size_t)nLoc+(std::size_t)loc;
          out << " " << res.detector_rate_flat[flat]
              << " " << res.detector_rate_lower_flat[flat]
              << " " << res.detector_rate_upper_flat[flat];
        }
        out << "\n";
      }
    }
    CloseAndRecordDensityArtifact_(out,fileName);

    // A shell spectrum is a new Step-6 artifact, so it need not alter the historical
    // combined density/flux schema.  One structured zone per energy makes every
    // reported integral independently reconstructable from emitted rows.
    const std::string spectrumStem="mode3d_shell_"+altLabel+"km_spectrum";
    const std::string spectrumFileName=DensityOutputFileName_(spectrumStem.c_str());
    std::ofstream spectrumOut=OpenDensityArtifact_(spectrumFileName);
    spectrumOut << std::setprecision(std::numeric_limits<double>::max_digits10);
    spectrumOut << "TITLE=\"Mode3D shell differential spectra\"\n";
    WriteProductMetadata_(spectrumOut,prm);
    spectrumOut << "VARIABLES=\"Lon_deg\" \"Lat_deg\" \"E_MeV\" \"T\" "
                << "\"T_lower\" \"T_upper\" \"unresolved_fraction\" "
                << "\"J_boundary_perMeV\" \"J_local_perMeV\" "
                << "\"J_local_lower_perMeV\" \"J_local_upper_perMeV\" "
                << "\"J_omni_perMeV\" \"J_omni_lower_perMeV\" \"J_omni_upper_perMeV\" "
                << "\"J_planar_perMeV\" \"J_planar_lower_perMeV\" \"J_planar_upper_perMeV\"\n";
    for (int ie=0;ie<nE;++ie) {
      const Geospace::BoundaryProducts::Bounds boundary=
          ::gSpectrum.GetSpectrumPerMeVBounds(E_MeV[(std::size_t)ie]);
      spectrumOut << "ZONE T=\"E_" << E_MeV[(std::size_t)ie]
                  << "MeV\" I=" << nLon << " J=" << nLat << " F=POINT\n";
      for (int j=0;j<nLat;++j) {
        double lat=-90.0+res_deg*(double)j;
        if (lat>90.0) lat=90.0;
        for (int i=0;i<nLon;++i) {
          const double lon=res_deg*(double)i;
          const int loc=s*nPtsShell+j*nLon+i;
          const std::size_t flat=(std::size_t)loc*(std::size_t)nE+(std::size_t)ie;
          const Geospace::BoundaryProducts::Bounds local=
              Geospace::BoundaryProducts::MapBoundaryIntensity(
                  Geospace::BoundaryProducts::Bounds(res.T_flat[flat],
                      res.T_lower_flat[flat],res.T_upper_flat[flat]),boundary,
                  Geospace::BoundaryProducts::CharacteristicMapping::StaticMagnetic);
          spectrumOut << lon << " " << lat << " " << E_MeV[(std::size_t)ie]
                      << " " << res.T_flat[flat] << " " << res.T_lower_flat[flat]
                      << " " << res.T_upper_flat[flat]
                      << " " << res.unresolved_fraction_flat[flat]
                      << " " << boundary.nominal
                      << " " << local.nominal << " " << local.lower << " " << local.upper
                      << " " << 4.0*Geospace::FluxNumerics::kPi*local.nominal
                      << " " << 4.0*Geospace::FluxNumerics::kPi*local.lower
                      << " " << 4.0*Geospace::FluxNumerics::kPi*local.upper
                      << " " << Geospace::FluxNumerics::kPi*local.nominal
                      << " " << Geospace::FluxNumerics::kPi*local.lower
                      << " " << Geospace::FluxNumerics::kPi*local.upper << "\n";
        }
      }
    }
    CloseAndRecordDensityArtifact_(spectrumOut,spectrumFileName);
  }
}

static void WriteTerminationSummary_(const GeospaceUtil::AmpsParam& prm,
                                     const std::vector<double>& E_MeV,
                                     int nLoc,const DensityResultBuffers& res) {
  const int nE=static_cast<int>(E_MeV.size());
  std::ofstream out;
  std::string fileName;
  if (prm.densitySpectrum.saveTerminationSummary) {
    fileName=DensityOutputFileName_("mode3d_termination_summary");
    out=OpenDensityArtifact_(fileName);
    out << std::setprecision(std::numeric_limits<double>::max_digits10);
    out << "TITLE=\"Mode3D trajectory termination summary\"\n";
    WriteProductMetadata_(out,prm);
    out << "VARIABLES=\"location_index\" \"E_MeV\" \"N_sampled\" \"N_retried\" "
        << "\"N_resolved\" \"N_allowed\" \"T\" \"T_lower\" \"T_upper\" "
        << "\"unresolved_fraction\"";
    for (int it=0;it<Geospace::FluxNumerics::kTerminationCount;++it) {
      const Geospace::GridlessMode::TrajectoryTermination termination=
          static_cast<Geospace::GridlessMode::TrajectoryTermination>(it);
      out << " \"N_" << Geospace::GridlessMode::TrajectoryTerminationName(termination) << "\"";
    }
    out << "\nZONE T=\"termination\" I=" << nLoc*nE << " F=POINT\n";
  }
  for (int loc=0;loc<nLoc;++loc) {
    for (int ie=0;ie<nE;++ie) {
      const std::size_t flat=(std::size_t)loc*(std::size_t)nE+(std::size_t)ie;
      if (!out.is_open()) continue;
      out << loc << " " << E_MeV[(std::size_t)ie] << " " << res.sampled_flat[flat]
          << " " << res.retried_flat[flat] << " " << res.resolved_flat[flat]
          << " " << res.allowed_flat[flat] << " " << res.T_flat[flat]
          << " " << res.T_lower_flat[flat] << " " << res.T_upper_flat[flat]
          << " " << res.unresolved_fraction_flat[flat];
      for (int it=0;it<Geospace::FluxNumerics::kTerminationCount;++it)
        out << " " << res.termination_flat[flat*(std::size_t)Geospace::FluxNumerics::kTerminationCount+
                                            (std::size_t)it];
      out << "\n";
    }
  }
  if (out.is_open()) CloseAndRecordDensityArtifact_(out,fileName);
}

} // anonymous namespace

namespace Geospace {
namespace Mode3D {

void SetDensityOutputFileSuffix(const std::string& suffix) {
  gDensityOutputFileSuffix = suffix;
}

Geospace::SWMFCoupledProducts::ProductControl DescribeDensityFluxProductControl(
    const GeospaceUtil::AmpsParam& prm) {
  return BuildProductControl_(prm);
}

int RunDensityAndFlux(const GeospaceUtil::AmpsParam& prm) {
  int mpiRank=0, mpiSize=1;
  MPI_Comm_rank(MPI_GLOBAL_COMMUNICATOR,&mpiRank);
  MPI_Comm_size(MPI_GLOBAL_COMMUNICATOR,&mpiSize);

  // Begin a new output transaction before validation or trajectory work.  If this call
  // fails, callers cannot accidentally retrieve the closed files or PASS-ready summary
  // from an earlier standalone/coupled field epoch.
  gLastDensityFluxArtifactFiles.clear();
  gLastDensityFluxRunSummary=
      Geospace::SWMFCoupledProducts::ProductRunSummary();

  // Start a new sample-weighted magnetic-field accuracy diagnostic for this density/flux
  // calculation.  ResetDipoleMagneticFieldErrorStatistics() enables collection only for
  // mesh-backed DIPOLE runs, so no extra analytic evaluations are made for other models.
  ResetDipoleMagneticFieldErrorStatistics(prm);

#ifdef SRC_EARTH_DENSITY_MODE3D_HAVE_MESH_EXPORT_
  // Batched engine: make the published compact field generation and AMR tree resident in
  // the batched trajectory engine before ComputeAllLocations_() starts any trajectory
  // fiber (the engine cannot replace its mesh while trajectories are in flight, and a
  // stale mesh inside a fiber is fatal in the engine-routed Mode3D tracer).
  // Standalone Mode3D::Run() already uploads after freezing the snapshot; this call is
  // idempotent (keyed by the snapshot generation) and covers the SWMF-coupled caller.
  // Rank-local (no MPI), main thread, outside fibers.  No-op for
  // AMPS_EARTH_ENGINE=LEGACY and for FORCE_ANALYTIC field evaluation, which never
  // samples the mesh.
  if (Geospace::GPU::BatchTracer::Enabled() && !prm.mode3d.forceAnalyticMagneticField)
    Geospace::GPU::Mode3DMesh::EnsureUploaded("Mode3D density/flux");
#endif

  const std::string outputMode = GeospaceUtil::ToUpper(prm.output.mode);
  if (outputMode!="POINTS" && outputMode!="TRAJECTORY" && outputMode!="SHELLS") {
    exit(__LINE__,__FILE__,"Mode3D density/flux supports OUTPUT_MODE POINTS, TRAJECTORY, or SHELLS");
  }

  // Time-dependent boundary spectra are selected once per magnetic-field snapshot.  This
  // matches the snapshot contract: every call to RunDensityAndFlux() sees one already
  // assembled compact field snapshot and should fold it with the boundary spectrum at
  // the same epoch.  TRAJECTORY per-sample epochs are intentionally not used to mutate
  // gSpectrum because the global arrays are not rebuilt per spacecraft sample.
  if (prm.densitySpectrum.spectrumEpochOffsetActive)
    ::gSpectrum.SetEvaluationEpochUTCOffset(
        prm.field.epoch,prm.densitySpectrum.spectrumEpochOffset_s);
  else
    ::gSpectrum.SetEvaluationEpochUTC(prm.field.epoch);

  const int nZenith = prm.densitySpectrum.directionNMu;
  const int nAz     = prm.densitySpectrum.directionNPhi;
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
  const Geospace::Mode3D::MpiScheduler mpiSchedulerForBanner =
      Geospace::Mode3D::ResolveMpiScheduler(prm,"Mode3D density/flux");
  const bool sharedEpochFineTasksForBanner=
      (outputMode=="POINTS" || outputMode=="SHELLS") && mpiSize>1;
  const long long mpiWorkItemsForBanner=sharedEpochFineTasksForBanner
      ? static_cast<long long>(nLoc)*static_cast<long long>(nE)
      : static_cast<long long>(nLoc);
  const long long mpiDynamicChunkForBanner = Geospace::Mode3D::ResolveMpiDynamicChunk(
      prm,densityThreadsForBanner,mpiWorkItemsForBanner);

  if (mpiRank==0) {
    std::cout << "================ Mode3D mesh density & flux ================\n";
    std::cout << "Field model     : " << prm.field.model << "\n";
    std::cout << "Epoch           : " << prm.field.epoch << "\n";
    std::cout << "Species         : " << prm.species.name << " (q=" << prm.species.charge_e
              << " e, m=" << prm.species.mass_amu << " amu)\n";
    std::cout << "Energy grid     : [" << prm.densitySpectrum.Emin_MeV << ", "
              << prm.densitySpectrum.Emax_MeV << "] MeV, Npoints=" << nE << "\n";
    std::cout << "Transmission    : " << GeospaceUtil::ToUpper(prm.densitySpectrum.transmissionMode);
    if (GeospaceUtil::ToUpper(prm.densitySpectrum.transmissionMode)!="DIRECT") {
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
    // Batched engine: in batch mode report the trajectory driver (see ComputeAllLocations_).
    // The line is additional and is not printed by LEGACY runs, whose banner is unchanged.
    if (Geospace::GPU::DriverBatchEnabled()) {
      std::cout << "Traj. engine    : BATCH ("
                << Geospace::GPU::BatchTracer::Get().BackendName() << "), "
                << DensityBatchHostThreads_(densityBackendForBanner,densityThreadsForBanner)
                << " host thread(s)/MPI rank; fibers replace the "
                << DensityParallelBackendName_(densityBackendForBanner) << " workers\n";
    }
    std::cout << "MPI work unit   : "
              << (sharedEpochFineTasksForBanner
                    ? "flattened (location,energy) shared-epoch task"
                    : (outputMode=="TRAJECTORY"
                         ? "rank-owned epoch/location; threads divide energy/directions"
                         : "location (single MPI rank)"))
              << "\n";
    std::cout << "MPI scheduler   : " << Geospace::Mode3D::MpiSchedulerName(mpiSchedulerForBanner) << "\n";
    if (mpiSchedulerForBanner == Geospace::Mode3D::MpiScheduler::DYNAMIC) {
      std::cout << "MPI dyn chunk   : " << mpiDynamicChunkForBanner
                << (sharedEpochFineTasksForBanner
                      ? " global (location,energy) task(s) per atomic fetch\n"
                      : " global location(s) per atomic fetch\n");
    }
    std::cout << "Progress bar    : ON\n";
    if (!gDensityOutputFileSuffix.empty()) std::cout << "Output suffix   : " << gDensityOutputFileSuffix << "\n";
    std::cout << "============================================================\n";
    std::cout.flush();
  }

  DensityResultBuffers res = ComputeAllLocations_(prm,nLoc,nLon,nLat,res_deg,nPtsShell,E_MeV,dirsUse);

  int outputSucceeded=1;
  std::string outputError;
  if (mpiRank==0) {
    try {
    // Construct and validate the trajectory accounting before the first product file
    // is opened.  All writers below then embed exactly the same summary/fingerprints;
    // no writer is allowed to reconstruct provenance independently.
    gLastDensityFluxRunSummary=BuildRunSummary_(
        prm,nLoc,nE,static_cast<int>(dirsUse.size()),res);
    WriteTerminationSummary_(prm,E_MeV,nLoc,res);
    const double maximumUnresolved=
        gLastDensityFluxRunSummary.maximumUnresolvedFraction;
    std::cout << "[mode3d-density][termination] max unresolved fraction = "
              << std::setprecision(17) << maximumUnresolved << " (tolerance="
              << prm.densitySpectrum.unresolvedTolerance << ")\n";
    if (maximumUnresolved>prm.densitySpectrum.unresolvedTolerance) {
      std::cout << "[mode3d-density][termination][WARNING] unresolved fraction exceeds tolerance; "
                << "nominal values are not replaced by physical zero and bounds remain in output.\n";
      if (prm.densitySpectrum.failOnUnresolved) {
        // The generic AMPS exit helper forwards the packed function code to
        // MPI_Abort.  Its process-global default is zero, which made a fatal scientific
        // gate look successful to mpiexec.  Set a nonzero code immediately before this
        // rank-0 fatal path; the SA01 artifact gate remains an independent safeguard.
        SetExitErrorCode(__LINE__,1);
        exit(__LINE__,__FILE__,"Mode3D density unresolved fraction exceeds DS_UNRESOLVED_TOL");
      }
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
    // Publish the artifact inventory only after every stream has been explicitly
    // closed.  The coupled bridge will apply the stronger complete-artifact and
    // unresolved-tolerance checks before it writes its PASS manifest.
    gLastDensityFluxRunSummary.artifacts=gLastDensityFluxArtifactFiles;
    Geospace::SWMFCoupledProducts::ValidateRunSummary(
        gLastDensityFluxRunSummary,false);
    std::cout.flush();
    }
    catch (const std::exception& error) {
      outputSucceeded=0;
      outputError=error.what();
    }
    catch (...) {
      outputSucceeded=0;
      outputError="unknown Mode3D density/flux output exception";
    }
  }

  // Output is root-owned but the caller and all remaining diagnostics are collective.
  // Propagate a close/manifest-accounting exception before any rank enters the next
  // reduction or barrier; otherwise rank zero could unwind while its peers wait
  // forever.  This synchronization changes no numerical gate or product value.
  MPI_Bcast(&outputSucceeded,1,MPI_INT,0,MPI_GLOBAL_COMMUNICATOR);
  int outputErrorLength=(mpiRank==0) ? static_cast<int>(outputError.size()) : 0;
  MPI_Bcast(&outputErrorLength,1,MPI_INT,0,MPI_GLOBAL_COMMUNICATOR);
  if (outputErrorLength<0)
    throw std::runtime_error("invalid collective Mode3D output error length");
  if (mpiRank!=0)
    outputError.assign(static_cast<std::size_t>(outputErrorLength),'\0');
  if (outputErrorLength>0)
    MPI_Bcast(&outputError[0],outputErrorLength,MPI_CHAR,0,
              MPI_GLOBAL_COMMUNICATOR);
  if (!outputSucceeded)
    throw std::runtime_error("Mode3D density/flux output transaction failed: "+
                             outputError);

  // Every TraceAllowedMesh() call has returned, so all per-trajectory field evaluators
  // have merged their local DIPOLE samples.  Perform one global reduction and print the
  // mean/max interpolation error and the location of the maximum before returning.
  ReportDipoleMagneticFieldErrorStatistics("Mode3D density and flux");

  MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
  return 0;
}

std::vector<std::string> GetLastDensityFluxArtifactFiles() {
  return gLastDensityFluxArtifactFiles;
}

Geospace::SWMFCoupledProducts::ProductRunSummary GetLastDensityFluxRunSummary() {
  return gLastDensityFluxRunSummary;
}

} // namespace Mode3D
} // namespace Geospace
