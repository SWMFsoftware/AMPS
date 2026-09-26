//======================================================================================
// Mode3DForwardSWMF.cpp
//======================================================================================
//
// SWMF-coupled initialization bridge for the Earth 3-D forward energetic-particle
// branch.
//
// DESIGN — TWO-PHASE INITIALIZATION
// -----------------------------------
// The original single amps_init() was called at the END of main_lib.cpp::amps_init(),
// after both amps_init_mesh() and amps_init() had already executed.  That ordering
// caused six initialization gaps relative to the standalone 3d_forward path:
//
//   Gap 1  Earth::ModelMode not set before amps_init_mesh()
//          → main_lib.cpp took the cutoff-rigidity branch instead of the
//            BoundaryInjection/MAIN_LIB_GEO branch.
//   Gap 2  cDensity3D::ConfigureEnergyGrid() called too late
//          → PIC::Mesh::initCellSamplingDataBuffer() allocated per-cell buffers
//            using the wrong (compile-time default) nEnergyBins.
//   Gap 3  Earth::Init_AfterParser() never called
//          → only reachable via MAIN_LIB_GEO::amps_init_mesh(), which was
//            never entered because of Gap 1.
//   Gap 4  BoundingBoxInjection::InitDirectionIMF() not called
//          → main_lib.cpp::amps_init() calls it in the BoundaryInjectionMode
//            branch, but that branch was skipped because of Gap 1.
//   Gap 5  PIC::Mover::ProcessOutsideDomainParticles left pointing at the
//          cutoff-rigidity handler set by the cutoff-rigidity amps_init_mesh()
//          path; should be null for the forward mode.
//   Gap 6  BoundingBoxInjection::SetPrm() not called before InitDirectionIMF().
//
// The fix splits initialization into two hooks:
//
//   amps_pre_init()  — called BEFORE amps_init_mesh() in main_lib.cpp
//      Fixes Gaps 1, 2, 3, 4, 6 by setting the mode and energy grid before
//      the mesh is built, and by registering prm before InitDirectionIMF().
//
//   amps_init()      — called at the END of main_lib.cpp::amps_init()
//      Fixes Gap 5.  Continues to do all post-mesh forward-mode setup,
//      reusing the prm parsed by amps_pre_init() (no second file read).
//
// The only intended physical difference from standalone 3d_forward remains the
// source of the background state: in _PIC_COUPLER_MODE__SWMF_ builds B/u come from
// PIC::CPLR rather than Tsyganenko/DATAFILE. Step-9 backward products are magnetic-
// only by default; derived E=-u x B is explicitly labelled experimental.
//
//======================================================================================

#include "Mode3DForwardSWMF.h"

#include "../3d_forward/Mode3DForward.h"
#include "../3d_forward/ForwardParticleMovers.h"
#include "../3d_forward/Density3D.h"
#include "../3d_forward/SphereFlux3D.h"
#include "../3d/Mode3D.h"
#include "../3d/CutoffRigidityMode3D.h"
#include "../3d/DensityMode3D.h"
#include "../3d/GlobalMagneticField.h"
#include "../gridless/DipoleInterface.h"
#include "../boundary/spectrum.h"
#include "../util/amps_param_parser.h"
#include "../util/SWMFSnapshotContract.h"
#include "../Earth.h"

#include "pic.h"

#include <iostream>
#include <stdexcept>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <vector>
#include <limits>
#include <cstring>
#include <memory>

#ifndef _NO_SPICE_CALLS_
#include "SpiceUsr.h"
#endif

namespace Earth {
namespace Mode3DForwardSWMF {

//======================================================================================
// Shared parsed parameters
//======================================================================================
// Populated once by amps_pre_init() and consumed read-only by amps_init().
// Using a file-scope static avoids re-parsing AMPS_PARAM.in and guarantees both
// hooks operate on an identical parameter set.
#if _PIC_COUPLER_MODE_ == _PIC_COUPLER_MODE__SWMF_
static EarthUtil::AmpsParam s_prm;
static bool                 s_pre_initialized = false;
static long int             s_cutoff_call_counter = 0;

// Suffix used by the most recent SWMF-coupled cutoff-rigidity output.
//
// The cutoff solver itself writes several products through
// Earth::Mode3D::SetCutoffOutputFileSuffix().  main_lib.cpp also writes the
// coupled AMPS mesh snapshot after the cutoff calculation finishes.  Keeping
// the last suffix here lets that mesh snapshot reuse exactly the same stamp as
// the cutoff products, so a single coupled snapshot produces a clearly matched
// file set.
static std::string          s_last_cutoff_output_suffix;
static bool                 s_last_cutoff_output_suffix_valid = false;
#endif


namespace {

void AnalyticDipoleMagneticField_(double *xIn,double *bOut) {
  Earth::GridlessMode::Dipole::GetB_Tesla(xIn,bOut);
}

bool TargetRequestsCutoff_(const EarthUtil::AmpsParam& prm) {
  const std::string t = EarthUtil::ToUpper(prm.calc.target);
  return t.find("CUTOFF") != std::string::npos || t=="ALL" || t=="BOTH";
}

bool TargetRequestsDensityFlux_(const EarthUtil::AmpsParam& prm) {
  const std::string t = EarthUtil::ToUpper(prm.calc.target);
  return t.find("DENSITY") != std::string::npos || t.find("FLUX") != std::string::npos ||
         t=="ALL" || t=="BOTH";
}

bool IsBackwardProductTarget_(const EarthUtil::AmpsParam& prm) {
  // The parser default is CUTOFF_RIGIDITY for standalone backward compatibility.
  // In coupled SWMF runs, keep requiring an explicit CALC_TARGET so old coupled
  // 3d_forward inputs that omit #CALCULATION_MODE keep their historical forward
  // particle-injection behavior.  Once CALC_TARGET is explicit, allow cutoff,
  // density/flux, or both.
  return prm.calc.targetExplicit &&
         (TargetRequestsCutoff_(prm) || TargetRequestsDensityFlux_(prm));
}

bool IsForwardTarget_(const EarthUtil::AmpsParam& prm) {
  return !IsBackwardProductTarget_(prm);
}

// Return the physical simulation-time interval, in seconds, between two
// expensive SWMF-coupled backward-product calculations.  Backward products are
// the particle-backtracking diagnostics handled by Mode3D in a coupled PT run:
// cutoff rigidity, density/flux, or both.  The same helper is used by two
// pieces of code that must stay synchronized:
//
//   1. main_lib.cpp::amps_time_step(), through the public
//      GetCoupledCalculationCadenceSeconds(), uses this value to decide whether
//      the current SWMF callback should run a new backtracking calculation or
//      be skipped until the next requested physical time.
//
//   2. ConfigureBackwardProductGlobalTimeStep_(), below, writes the same value
//      into PIC::ParticleWeightTimeStep::GlobalTimeStep[0].  In SWMF-coupled
//      runs that array is also visible to the generic AMPS/PT time-step
//      infrastructure.  For cutoff/density backtracking there is no forward
//      particle population whose CFL/gyro time step should control the PT
//      component.  Therefore species 0 is assigned the diagnostic cadence
//      itself: one PT time step corresponds to one requested cutoff/density
//      calculation interval.
//
// Units: #TEMPORAL/FIELD_UPDATE_DT is specified in minutes in AMPS_PARAM.in;
// this helper returns seconds.  If the input does not provide a positive value,
// preserve the historical coupled-diagnostic default of one hour.
double BackwardProductCadenceSecondsFromPrm_(const EarthUtil::AmpsParam& prm) {
  if (prm.temporal.fieldUpdateDt_min > 0.0) {
    return 60.0*prm.temporal.fieldUpdateDt_min;
  }

  return 3600.0;
}

// Configure the AMPS/PIC global time step used by the SWMF-coupled Mode3D
// backtracking products.  This function is intentionally private to the SWMF
// bridge because standalone Mode3D does not use PIC::SimulationTime to schedule
// a live coupled component; standalone snapshot looping is driven directly by
// the parsed input file.
//
// Why GlobalTimeStep[0] must be set here:
//   The cutoff and density/flux products are not advanced by a normal forward
//   particle push.  They are instantaneous diagnostics computed by backward
//   tracing many test trajectories through the current meshed magnetic field.
//   In a live SWMF-coupled run the physically meaningful interval is therefore
//   the interval between two diagnostic evaluations, not a particle CFL time
//   step.  Assigning GlobalTimeStep[0] to that interval makes the PT component's
//   species-0 time step consistent with the cadence gate that actually launches
//   the calculations.
//
// Only species zero is modified because the backtracking products currently use
// the proton/species-0 conventions throughout the Mode3D cutoff and density
// modules.  Forward 3d_forward_swmf runs are not affected; they overwrite
// GlobalTimeStep[species] later with the true forward-particle integration step.
//
// The function returns the value written so callers can include it in diagnostic
// messages without recomputing or risking a mismatch.
double ConfigureBackwardProductGlobalTimeStep_(const EarthUtil::AmpsParam& prm,
                                               bool verbose) {
  const double cadence_s = BackwardProductCadenceSecondsFromPrm_(prm);

  if (!std::isfinite(cadence_s) || cadence_s <= 0.0) {
    exit(__LINE__,__FILE__,
         "[Mode3DForwardSWMF] Invalid coupled Mode3D backward-product cadence. "
         "#TEMPORAL/FIELD_UPDATE_DT must define a positive interval in minutes, "
         "or be omitted to use the one-hour default.");
  }

  PIC::ParticleWeightTimeStep::GlobalTimeStep[0] = cadence_s;

  if (verbose && PIC::ThisThread == 0) {
    std::cout << "[Mode3DForwardSWMF] Set "
              << "PIC::ParticleWeightTimeStep::GlobalTimeStep[0]="
              << cadence_s
              << " s for SWMF-coupled Mode3D backward products "
              << "(cutoff and/or density/flux cadence).\n";
    std::cout.flush();
  }

  return cadence_s;
}

void ApplyParsedDomainForCutoff_(const EarthUtil::AmpsParam& prm) {
  Earth::Mode3D::ParsedDomainActive = true;
  Earth::Mode3D::ParsedDomainMin[0] = prm.domain.xMin * 1000.0;
  Earth::Mode3D::ParsedDomainMin[1] = prm.domain.yMin * 1000.0;
  Earth::Mode3D::ParsedDomainMin[2] = prm.domain.zMin * 1000.0;
  Earth::Mode3D::ParsedDomainMax[0] = prm.domain.xMax * 1000.0;
  Earth::Mode3D::ParsedDomainMax[1] = prm.domain.yMax * 1000.0;
  Earth::Mode3D::ParsedDomainMax[2] = prm.domain.zMax * 1000.0;
}

double CurrentCutoffOutputTimeSeconds_() {
  // Use the simulation clock maintained by AMPS/SWMF coupling, not a locally
  // reconstructed time.  The previous implementation advanced an internal
  // elapsed-time estimate by the configured output/coupling cadence after every
  // cutoff call.  That was vulnerable to drift whenever the actual SWMF coupling
  // schedule differed from the nominal cadence, or when the first/last snapshot
  // was not exactly aligned with that cadence.
  //
  // PIC::SimulationTime::TimeCounter is the authoritative AMPS simulation time
  // visible to the PT component.  Using it here makes the time part of the file
  // stamp correspond to the exact coupled snapshot being processed.
  const double tSim_s = PIC::SimulationTime::TimeCounter;

  // File names should remain well formed even if a bad clock value is exposed
  // during early initialization or error recovery.  Clamp only for the purpose
  // of naming the file; the physical solver state itself is not modified here.
  if (!std::isfinite(tSim_s)) return 0.0;

  return std::max(0.0,tSim_s);
}

std::string FormatCutoffOutputSuffix_(long int callIndex) {
  const double tSim_s = CurrentCutoffOutputTimeSeconds_();

  std::ostringstream ss;

  // Keep the call counter in the suffix even though the physical time is now
  // taken from PIC::SimulationTime::TimeCounter.  The counter guarantees unique
  // names if two coupling callbacks occur at the same simulation time, while
  // the TimeCounter part preserves the physically meaningful snapshot time.
  ss << ".swmf_n" << std::setw(6) << std::setfill('0') << callIndex
     << "_t" << std::setw(12) << std::setfill('0')
     << std::fixed << std::setprecision(3) << tSim_s << "s";

  return ss.str();
}

#if _PIC_COUPLER_MODE_ == _PIC_COUPLER_MODE__SWMF_
void RequireCollectiveStatusWriteSuccess_(int rootSucceeded,
                                          std::string rootError) {
  // Product callbacks are collective. If rank zero alone failed while creating the
  // status artifact, other ranks could continue into field collectives and deadlock.
  // Broadcast both outcome and diagnostic so all ranks take the same failure path.
  MPI_Bcast(&rootSucceeded,1,MPI_INT,0,MPI_GLOBAL_COMMUNICATOR);
  int messageLength=(PIC::ThisThread==0) ?
      static_cast<int>(rootError.size()) : 0;
  MPI_Bcast(&messageLength,1,MPI_INT,0,MPI_GLOBAL_COMMUNICATOR);
  if (messageLength<0)
    throw std::runtime_error("Invalid collective SWMF status error length");
  if (PIC::ThisThread!=0)
    rootError.assign(static_cast<std::size_t>(messageLength),'\0');
  if (messageLength>0)
    MPI_Bcast(&rootError[0],messageLength,MPI_CHAR,0,MPI_GLOBAL_COMMUNICATOR);
  if (!rootSucceeded) {
    std::string message="SWMF product-status write failed";
    if (!rootError.empty()) message+=": "+rootError;
    throw std::runtime_error(message);
  }
}

std::string CurrentCoupledEpochUTC_() {
  const std::string reference=s_prm.field.epoch.empty() ?
      std::string("UNSPECIFIED") : s_prm.field.epoch;
  const double offset_s=CurrentCutoffOutputTimeSeconds_();
#ifdef _NO_SPICE_CALLS_
  if (offset_s!=0.0)
    throw std::runtime_error(
        "Step-9 SWMF snapshot epoch requires SPICE to add nonzero PT simulation time");
  return reference;
#else
  SpiceDouble et=0.0;
  str2et_c(reference.c_str(),&et);
  et+=static_cast<SpiceDouble>(offset_s);
  char utc[64];
  // Retain nanosecond text precision so a sub-millisecond PT clock is not rounded to
  // another boundary-spectrum/ephemeris epoch. The exact double simulation time is
  // also stored separately in the content fingerprint.
  et2utc_c(et,"ISOC",9,static_cast<SpiceInt>(sizeof(utc)),utc);
  return std::string(utc)+"Z";
#endif
}

unsigned long long StableControlStringHash_(const std::string& value,
                                             unsigned long long seed) {
  // Local FNV-1a is sufficient here: the value is compared only among ranks in this
  // executable, not persisted as scientific provenance. Include an explicit terminator
  // so concatenated selectors cannot alias ("ab"+"c" versus "a"+"bc").
  unsigned long long hash=seed;
  for (std::string::const_iterator it=value.begin();it!=value.end();++it) {
    hash^=static_cast<unsigned char>(*it);
    hash*=1099511628211ULL;
  }
  hash^=0xffU;
  hash*=1099511628211ULL;
  return hash;
}

void ValidateLiveSWMFStateCoherence_() {
  // A coupled callback is expected to expose one logical MHD state on every rank.
  // Boolean readiness alone is insufficient: rank-local clocks or buffer layouts could
  // differ and still all be "ready".  Compare the exact state selectors collectively
  // before any owner cell is copied so the compact arrays cannot combine epochs.
  const double localTime=PIC::SimulationTime::TimeCounter;
  double minTime=0.0,maxTime=0.0;
  MPI_Allreduce((void*)&localTime,&minTime,1,MPI_DOUBLE,MPI_MIN,
                MPI_GLOBAL_COMMUNICATOR);
  MPI_Allreduce((void*)&localTime,&maxTime,1,MPI_DOUBLE,MPI_MAX,
                MPI_GLOBAL_COMMUNICATOR);
  if (!std::isfinite(localTime) || !std::isfinite(minTime) ||
      !std::isfinite(maxTime) || localTime<0.0 || minTime!=maxTime)
    throw std::runtime_error(
        "SWMF PT simulation time is negative, non-finite, or inconsistent across MPI ranks");

  // A receive may be repeated at the same epoch (restart/retry), but an older receive
  // must never replace a newer one in the same process. Equality is intentionally
  // allowed so restart identity remains reproducible.
  static double lastAcceptedTime_s=-1.0;
  if (lastAcceptedTime_s>=0.0 && localTime<lastAcceptedTime_s)
    throw std::runtime_error(
        "stale SWMF receive time is older than the last accepted field snapshot");

  const long int selectors[3]={
    PIC::CPLR::SWMF::MagneticFieldOffset,
    PIC::CPLR::SWMF::BulkVelocityOffset,
    s_cutoff_call_counter};
  for (int n=0;n<3;++n) {
    long int minimum=0,maximum=0;
    MPI_Allreduce((void*)&selectors[n],&minimum,1,MPI_LONG,MPI_MIN,
                  MPI_GLOBAL_COMMUNICATOR);
    MPI_Allreduce((void*)&selectors[n],&maximum,1,MPI_LONG,MPI_MAX,
                  MPI_GLOBAL_COMMUNICATOR);
    if (minimum!=maximum)
      throw std::runtime_error(
          "SWMF B/u offsets or backward-product generation differ across MPI ranks");
  }

  // Domain faces and the electric-field release mode are physical selectors, not just
  // input cosmetics. Exact collective agreement prevents ranks from hashing/tracing
  // different boxes or mixing magnetic-only and experimental-E metadata.
  const double domainValues[6]={
    s_prm.domain.xMin,s_prm.domain.xMax,s_prm.domain.yMin,
    s_prm.domain.yMax,s_prm.domain.zMin,s_prm.domain.zMax};
  for (int n=0;n<6;++n) {
    double minimum=0.0,maximum=0.0;
    MPI_Allreduce((void*)&domainValues[n],&minimum,1,MPI_DOUBLE,MPI_MIN,
                  MPI_GLOBAL_COMMUNICATOR);
    MPI_Allreduce((void*)&domainValues[n],&maximum,1,MPI_DOUBLE,MPI_MAX,
                  MPI_GLOBAL_COMMUNICATOR);
    if (!std::isfinite(domainValues[n]) || minimum!=maximum)
      throw std::runtime_error(
          "SWMF validity domain differs across MPI ranks");
  }
  const int localModes[2]={
    s_prm.field.swmfDerivedElectricFieldExperimental ? 1 : 0,
    s_prm.field.swmfSnapshotExport ? 1 : 0};
  for (int n=0;n<2;++n) {
    int minimum=0,maximum=0;
    MPI_Allreduce((void*)&localModes[n],&minimum,1,MPI_INT,MPI_MIN,
                  MPI_GLOBAL_COMMUNICATOR);
    MPI_Allreduce((void*)&localModes[n],&maximum,1,MPI_INT,MPI_MAX,
                  MPI_GLOBAL_COMMUNICATOR);
    if (minimum!=maximum)
      throw std::runtime_error(
          "SWMF electric-field or snapshot-export mode differs across MPI ranks");
  }

  // Differing product targets would send ranks into different solver collectives, and
  // differing epoch/export strings would create contradictory provenance. Compare a
  // stable hash of all three before gathering the field or branching on the target.
  unsigned long long localControlHash=14695981039346656037ULL;
  localControlHash=StableControlStringHash_(s_prm.field.epoch,localControlHash);
  localControlHash=StableControlStringHash_(
      s_prm.field.swmfSnapshotExportPrefix,localControlHash);
  localControlHash=StableControlStringHash_(s_prm.calc.target,localControlHash);
  unsigned long long minimumControlHash=0,maximumControlHash=0;
  MPI_Allreduce((void*)&localControlHash,&minimumControlHash,1,
                MPI_UNSIGNED_LONG_LONG,MPI_MIN,MPI_GLOBAL_COMMUNICATOR);
  MPI_Allreduce((void*)&localControlHash,&maximumControlHash,1,
                MPI_UNSIGNED_LONG_LONG,MPI_MAX,MPI_GLOBAL_COMMUNICATOR);
  if (minimumControlHash!=maximumControlHash)
    throw std::runtime_error(
        "SWMF epoch, export prefix, or product target differs across MPI ranks");
  lastAcceptedTime_s=localTime;
}

Earth::Field::SnapshotMetadata BuildSWMFFieldSnapshotMetadata_() {
  Earth::Field::SnapshotMetadata metadata;
  metadata.sourceId="PIC::CPLR:SWMF";
  metadata.modelName="SWMF";

  // Convert the configured absolute reference plus the authoritative PT clock to a
  // real ISO UTC.  Standalone replay and observation/ephemeris selection must consume
  // a parseable physical epoch, not an implementation string such as "base+PT...".
  metadata.epochUTC=CurrentCoupledEpochUTC_();
  metadata.frame=Earth::Field::CoordinateFrame::GSM;
  metadata.interpolation=s_prm.field.swmfDerivedElectricFieldExperimental ?
      Earth::Field::InterpolationMode::CellCenteredLinearDerivedElectric :
      Earth::Field::InterpolationMode::CellCenteredLinear;
  metadata.magneticFieldAvailable=true;
  metadata.electricFieldAvailable=
      s_prm.field.swmfDerivedElectricFieldExperimental;
  metadata.immutableDuringBatch=true;
  metadata.valid=true;

  metadata.domain.enabled=true;
  metadata.domain.minimum_m[0]=1000.0*s_prm.domain.xMin;
  metadata.domain.minimum_m[1]=1000.0*s_prm.domain.yMin;
  metadata.domain.minimum_m[2]=1000.0*s_prm.domain.zMin;
  metadata.domain.maximum_m[0]=1000.0*s_prm.domain.xMax;
  metadata.domain.maximum_m[1]=1000.0*s_prm.domain.yMax;
  metadata.domain.maximum_m[2]=1000.0*s_prm.domain.zMax;

  // This is a provisional identity used only while validating the gather request.  The
  // compact-array publisher replaces it with the Step-9 content-derived identity after
  // B/u have been reduced and validated.  Thus restarts and different MPI layouts get
  // the same final ID for the same epoch/data, while the call number remains only an
  // output-name discriminator.
  std::ostringstream state;
  state << std::setprecision(17)
        << "call=" << s_cutoff_call_counter
        << "|time=" << CurrentCutoffOutputTimeSeconds_()
        << "|Boffset=" << PIC::CPLR::SWMF::MagneticFieldOffset
        << "|Voffset=" << PIC::CPLR::SWMF::BulkVelocityOffset
        << "|xmin=" << s_prm.domain.xMin << "|xmax=" << s_prm.domain.xMax
        << "|ymin=" << s_prm.domain.yMin << "|ymax=" << s_prm.domain.yMax
        << "|zmin=" << s_prm.domain.zMin << "|zmax=" << s_prm.domain.zMax;
  metadata.snapshotId=Earth::Field::MakeSnapshotId(
      metadata.sourceId,metadata.epochUTC,state.str());
  return metadata;
}

std::string CoupledSnapshotExportFile_(const std::string& suffix) {
  return s_prm.field.swmfSnapshotExport ?
      s_prm.field.swmfSnapshotExportPrefix+suffix+".csv" : std::string();
}

std::string CoupledSnapshotStatusFile_(const std::string& suffix) {
  return "swmf_product_status"+suffix+".json";
}

void WriteCoupledSnapshotStatus_(const std::string& fileName,
                                 const std::string& status,
                                 const Earth::Field::SnapshotMetadata& metadata,
                                 const std::string& suffix,
                                 const std::string& exportedSnapshot,
                                 const std::string& message) {
  int writeSucceeded=1;
  std::string writeError;
  if (PIC::ThisThread==0) {
    try {
      Earth::SWMFSnapshot::WriteProductStatus(
          fileName,
          Earth::SWMFSnapshot::BuildProductStatusJson(
              status,metadata.snapshotId,metadata.epochUTC,suffix,exportedSnapshot,
              TargetRequestsCutoff_(s_prm),TargetRequestsDensityFlux_(s_prm),message));
    }
    catch (const std::exception& error) {
      writeSucceeded=0;
      writeError=error.what();
    }
    catch (...) {
      writeSucceeded=0;
      writeError="unknown exception";
    }
  }
  RequireCollectiveStatusWriteSuccess_(writeSucceeded,writeError);
}
#endif

//======================================================================================
// SWMF-coupled cutoff global-field assembly is implemented by the shared
// Earth::Mode3D::GlobalMagneticField helper.  Both standalone and coupled paths reset
// and assign node->Temp_ID, gather owner interior cells, and create compact global
// arrays. No nonlocal AMR blocks or ghost-cell state vectors are allocated. For SWMF,
// B and u are read from their authoritative offsets. The released OFF mode keeps E
// unavailable; only explicit EXPERIMENTAL mode derives E=-u x B. In both modes u is
// retained so export/replay and pointwise diagnostics compare the same physical state.
//======================================================================================

} // anonymous namespace

//======================================================================================
// amps_pre_init  —  MUST be called before amps_init_mesh()
//======================================================================================
//
// Performs the subset of 3d_forward initialization that must precede the AMPS
// mesh-build phase.  main_lib.cpp::amps_init_mesh() calls this at its very start
// under a #if _PIC_COUPLER_MODE__SWMF_ guard.
//
// After this function returns:
//   • Earth::ModelMode == BoundaryInjectionMode
//     → amps_init_mesh() routes to MAIN_LIB_GEO::amps_init_mesh(), which calls
//       Earth::Init_AfterParser() (fixes Gap 3).
//   • cDensity3D::nEnergyBins/Emin_J/Emax_J reflect the input file
//     → initCellSamplingDataBuffer() allocates correctly-sized per-cell buffers
//       (fixes Gap 2).
//   • BoundingBoxInjection::prm is set
//     → main_lib.cpp::amps_init() can call InitDirectionIMF() correctly
//       (fixes Gaps 4 and 6).
//
void amps_pre_init() {
#if _PIC_COUPLER_MODE_ != _PIC_COUPLER_MODE__SWMF_
  return;
#else
  const char* inputFile = "AMPS_PARAM.in";

  if (PIC::ThisThread == 0)
    std::cout << "[Mode3DForwardSWMF] amps_pre_init: reading '" << inputFile << "'.\n";

  // Parse once.  amps_init() will use s_prm directly — no second read.
  s_prm = EarthUtil::ParseAmpsParamFile(inputFile);

  // Mark the parameter object as SWMF-coupled: forces field access through
  // PIC::CPLR instead of Tsyganenko/DATAFILE in the shared 3-D helpers.
  s_prm.field.model                       = "SWMF";
  s_prm.calc.fieldEvalMethod              = "SWMF";
  s_prm.mode3d.forceAnalyticMagneticField = false;

  if (IsBackwardProductTarget_(s_prm)) {
    // Coupled Mode3D backward-product mode.  Keep the regular AMPS/SWMF mesh path,
    // but route per-coupling-step work to the Mode3D backward-product driver instead
    // of the 3d_forward injector.  The mesh tracer will sample the current SWMF
    // fields through the globalized PIC::CPLR magnetic-field data.
    Earth::ModelMode = Earth::CutoffRigidityMode;
    ApplyParsedDomainForCutoff_(s_prm);

    s_pre_initialized = true;

    if (PIC::ThisThread == 0)
      std::cout << "[Mode3DForwardSWMF] amps_pre_init complete:"
                << " ModelMode=CutoffRigidityMode"
                << ", field source=PIC::CPLR/SWMF"
                << ", output will be timestamped per amps_time_step().\n";
    return;
  }

  // Historical coupled 3d_forward mode.
  // ── Fix Gap 1 ────────────────────────────────────────────────────────────
  // Set Earth::ModelMode before amps_init_mesh() is called.
  // main_lib.cpp::amps_init_mesh() checks this at its entry:
  //   if (Earth::ModelMode == BoundaryInjectionMode)  →  MAIN_LIB_GEO path
  //   else                                            →  cutoff-rigidity path
  // The cutoff-rigidity path misses Earth::Init_AfterParser() and allocates
  // GeospaceFlag cell data that 3d_forward does not use.
  Earth::ModelMode = Earth::BoundaryInjectionMode;

  // ── Fix Gap 2 ────────────────────────────────────────────────────────────
  // Configure the cDensity3D energy grid before amps_init_mesh() calls
  // PIC::Mesh::initCellSamplingDataBuffer().  That function queries
  // cDensity3D::RequestSamplingData (already pushed into
  // PIC::IndividualModelSampling::RequestSamplingData by the first line of
  // amps_init_mesh()) to calculate each cell's sampling-buffer byte count:
  //     bytes = nEnergyBins * sizeof(double)
  // If nEnergyBins is still at its compile-time default here, every cell's
  // buffer will be under-allocated.
  Earth::Mode3DForward::cDensity3D::ConfigureEnergyGrid(s_prm);

  // ── Fix Gaps 4 & 6 ───────────────────────────────────────────────────────
  // Register prm with the bounding-box injector before amps_init_mesh() and
  // amps_init() run.  main_lib.cpp::amps_init() calls
  // Earth::BoundingBoxInjection::InitDirectionIMF() in the BoundaryInjectionMode
  // branch (now reachable because of the Gap-1 fix above), and InitDirectionIMF()
  // needs a valid prm to evaluate the background field at the domain boundary.
  // In the SWMF build ConfigureBackgroundFieldModel() is a no-op, so
  // InitDirectionIMF() falls through to the PIC::CPLR path or the safe default.
  Earth::BoundingBoxInjection::SetPrm(s_prm);

  s_pre_initialized = true;

  if (PIC::ThisThread == 0)
    std::cout << "[Mode3DForwardSWMF] amps_pre_init complete:"
              << " ModelMode=BoundaryInjectionMode"
              << ", nEnergyBins=" << Earth::Mode3DForward::cDensity3D::nEnergyBins
              << ".\n";
#endif
}


//======================================================================================
// amps_init  —  called at the END of main_lib.cpp::amps_init()
//======================================================================================
//
// Completes the SWMF-coupled 3d_forward runtime after the full AMPS mesh and PIC
// core have been initialized.  Uses s_prm cached by amps_pre_init() — AMPS_PARAM.in
// is not re-read.
//
void amps_init() {
#if _PIC_COUPLER_MODE_ != _PIC_COUPLER_MODE__SWMF_
  return;
#else
  static bool initialized = false;
  if (initialized) {
    if (PIC::ThisThread == 0)
      std::cout << "[Mode3DForwardSWMF] Coupled 3d_forward runtime was already "
                   "initialized; skipping duplicate initialization.\n";
    return;
  }
  initialized = true;

  // Guard: amps_pre_init() must have run first.
  if (!s_pre_initialized)
    exit(__LINE__, __FILE__,
         "[Mode3DForwardSWMF] amps_init() called before amps_pre_init(). "
         "Ensure amps_pre_init() is called at the start of amps_init_mesh().");

  // Use the prm already parsed in amps_pre_init() — do not re-read the file.
  const EarthUtil::AmpsParam& prm = s_prm;

  if (IsBackwardProductTarget_(prm)) {
    // Coupled Mode3D backward-product mode needs only the generic AMPS/PIC
    // initialization already performed in main_lib.cpp before this hook is
    // called.  Do not register forward-injection callbacks.  The requested
    // cutoff and/or density/flux calculation is launched from amps_cutoff_time_step()
    // each time the cadence gate accepts an SWMF/PT callback and the current MHD
    // fields have been imported.
    Earth::ModelMode = Earth::CutoffRigidityMode;
    ApplyParsedDomainForCutoff_(prm);

    for (int s=0; s<PIC::nTotalSpecies; ++s) {
      PIC::ParticleWeightTimeStep::SetGlobalParticleWeight(s,1.0);
    }

    // Synchronize the generic AMPS/PIC species-0 time step with the physical
    // cadence of the SWMF-coupled backtracking diagnostics.  The cadence gate in
    // main_lib.cpp uses exactly the same value through
    // GetCoupledCalculationCadenceSeconds(), so the stored time step, skip/run
    // decision, and output timestamps are all based on one common interval.
    const double backwardProductCadence_s =
        ConfigureBackwardProductGlobalTimeStep_(prm,/*verbose=*/true);

    PIC::Mover::BackwardTimeIntegrationMode = _PIC_MODE_OFF_;

    if (PIC::ThisThread == 0)
      std::cout << "[Mode3DForwardSWMF] Initialized SWMF-coupled Mode3D "
                << "backward-product runtime. Calculation cadence="
                << backwardProductCadence_s
                << " s; each accepted amps_time_step() call writes "
                << "timestamped cutoff/density outputs as requested by CALC_TARGET.\n";
    return;
  }

  if (PIC::ThisThread == 0)
    std::cout << "[Mode3DForwardSWMF] amps_init: completing coupled 3d_forward setup.\n";

  // Route AMPS particle motion through Earth::ParticleMover() →
  // Earth3DForward::MoverManager().  Earth::ModelMode was already set to
  // BoundaryInjectionMode in amps_pre_init(); confirm it is still correct.
  Earth::ModelMode = Earth::BoundaryInjectionMode;

  if (!Earth::Earth3DForward::SetMoverByName(prm.mode3dForward.particleMover)) {
    throw std::runtime_error(
        "Unknown 3d_forward_swmf particle mover '" +
        prm.mode3dForward.particleMover +
        "'. Valid movers are BORIS, RK4, GC/GC4, and HYBRID.");
  }

  Earth::Earth3DForward::ConfigureFieldEvaluation(prm);

  Earth::Mode3DForward::sInjectionEnergyDistribution =
      Earth::Mode3DForward::ParseInjectionEnergyDistribution(
          prm.mode3dForward.injectionEnergyDistribution);

#if _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ != _INDIVIDUAL_PARTICLE_WEIGHT_ON_
  if (Earth::Mode3DForward::sInjectionEnergyDistribution ==
      Earth::Mode3DForward::InjectionEnergyDistribution::LOG_UNIFORM) {
    throw std::runtime_error(
        "3d_forward_swmf LOG_UNIFORM injection requires AMPS to be compiled with "
        "_INDIVIDUAL_PARTICLE_WEIGHT_MODE_ == _INDIVIDUAL_PARTICLE_WEIGHT_ON_.");
  }
#endif

  Earth::Mode3DForward::sInitializeParticleTrajectories =
      prm.mode3dForward.initializeParticleTrajectories;
  Earth::Mode3DForward::sMaxParticleTrajectories = prm.mode3dForward.nParticleTrajectories;

  // Register the 3d_forward source-rate and particle-injection callbacks.
  // These override any BoundaryInjectionMode callbacks set by main_lib.cpp::amps_init()
  // (LocalBlockInjectionRate / InjectionProcessor) with the correct 3d_forward ones.
  PIC::ParticleWeightTimeStep::UserDefinedExtraSourceRate =
      Earth::Mode3DForward::BoundaryInjectionSourceRate;
  PIC::BC::UserDefinedParticleInjectionFunction =
      Earth::Mode3DForward::InjectParticles;

  // In an SWMF-coupled build this call is a deliberate no-op: it leaves the field
  // source neutral so all field access goes through PIC::CPLR/SWMF.
  Earth::Mode3DForward::ConfigureBackgroundFieldModel(prm);

  // The 3d_forward injector does not use the IMF direction vector b to construct
  // particle velocities (it samples the full cosine-weighted inward hemisphere per
  // face).  Set a safe default so diagnostic messages that print b are well-defined.
  // Note: BoundingBoxInjection::SetPrm() was already called in amps_pre_init() and
  // InitDirectionIMF() was called by main_lib.cpp::amps_init() for BoundaryInjectionMode.
  Earth::BoundingBoxInjection::b[0] = 1.0;
  Earth::BoundingBoxInjection::b[1] = 0.0;
  Earth::BoundingBoxInjection::b[2] = 0.0;

#if _PIC_PARTICLE_TRACKER_MODE_ == _PIC_MODE_ON_
  Earth::ParticleTracker::ConfigureRuntimeTrajectoryTracking(
      prm.mode3dForward.initializeParticleTrajectories,
      prm.mode3dForward.nParticleTrajectories,
      true);
#else
  if (Earth::Mode3DForward::sInitializeParticleTrajectories && PIC::ThisThread == 0) {
    std::cerr << "[Mode3DForwardSWMF] WARNING: trajectory initialization was "
              << "requested, but AMPS was compiled with "
              << "_PIC_PARTICLE_TRACKER_MODE_ != _PIC_MODE_ON_. No trajectory "
              << "records will be initialized.\n";
  }
#endif

  // cDensity3D::ConfigureEnergyGrid() was already called in amps_pre_init().
  // Calling it again here would be harmless but wasteful; the values are unchanged.

  // Enable the AMPS per-particle sampling dispatcher and register the 3d_forward
  // density sampler.  The duplicate-check handles the case where the generic
  // startup path already inserted the callback (cutoff-rigidity amps_init_mesh path).
  Earth::Sampling::ParticleData::SamplingMode = true;
  bool densitySamplerAlreadyRegistered = false;
  for (auto fn : Earth::Sampling::ParticleData::SampleParticleDataCallbacks) {
    if (fn == Mode3DForward::cDensity3D::SampleParticleData) {
      densitySamplerAlreadyRegistered = true;
      break;
    }
  }
  if (!densitySamplerAlreadyRegistered) {
    Earth::Sampling::ParticleData::SampleParticleDataCallbacks.push_back(
        Mode3DForward::cDensity3D::SampleParticleData);
  }

  // Initialize the spectrum before evaluating the source rate and particle weight.
  // 3d_forward convention: #SPECTRUM supplies the shape/normalization and
  // #DENSITY_3D supplies the actual simulated energy range.
  const auto forwardSpectrum =
      Earth::Mode3DForward::BuildForwardInjectionSpectrumMap(prm);
  InitGlobalSpectrumFromKeyValueMap(forwardSpectrum);

  Earth::Mode3DForward::sSpecies = 0;
  Earth::Mode3DForward::sDt      = Earth::Mode3DForward::EvaluateTimeStep(prm);
  PIC::ParticleWeightTimeStep::GlobalTimeStep[Earth::Mode3DForward::sSpecies] =
      Earth::Mode3DForward::sDt;

  Earth::Mode3DForward::sNParticlesPerIter = prm.mode3dForward.nParticlesPerIter;
  PIC::ParticleWeightTimeStep::maxReferenceInjectedParticleNumber =
      Earth::Mode3DForward::sNParticlesPerIter;
  PIC::ParticleWeightTimeStep::initParticleWeight_ConstantWeight(
      Earth::Mode3DForward::sSpecies);
  Earth::Mode3DForward::sParticleWeight =
      PIC::ParticleWeightTimeStep::GlobalParticleWeight[Earth::Mode3DForward::sSpecies];

  Earth::Mode3DForward::InitAbsorptionSphere(prm);

  // ── Fix Gap 5 ────────────────────────────────────────────────────────────
  // The cutoff-rigidity path in main_lib.cpp::amps_init_mesh() sets
  //   PIC::Mover::ProcessOutsideDomainParticles =
  //       Earth::CutoffRigidity::ProcessOutsideDomainParticles
  // With the amps_pre_init() fix, amps_init_mesh() now takes the
  // MAIN_LIB_GEO path, which does NOT set this pointer (correct for
  // 3d_forward).  This reset is kept as an explicit guard in case another
  // code path sets it between amps_pre_init() and here.
  PIC::Mover::ProcessOutsideDomainParticles = nullptr;

  Mode3DForward::cDensity3D::Init(prm);
  Mode3DForward::cSphereFlux3D::Init(
      prm, Earth::Mode3DForward::sAbsorptionSphere, Earth::Mode3DForward::sDt);
  Earth::Mode3DForward::InitBoundaryInjectionTable();

  PIC::Mover::BackwardTimeIntegrationMode = _PIC_MODE_OFF_;
  PIC::SamplingMode = _RESTART_SAMPLING_MODE_;

  if (PIC::ThisThread == 0) {
    std::cout << "[Mode3DForwardSWMF] Initialized SWMF-coupled 3d_forward runtime. "
              << "Field source=PIC::CPLR/SWMF"
              << ", mover="       << Earth::Earth3DForward::GetMoverName()
              << ", nParticlesPerIter=" << Earth::Mode3DForward::sNParticlesPerIter
              << ", particleWeight="    << Earth::Mode3DForward::sParticleWeight
              << ", dt="               << Earth::Mode3DForward::sDt
              << ", energySampling="
              << Earth::Mode3DForward::InjectionEnergyDistributionName(
                     Earth::Mode3DForward::sInjectionEnergyDistribution)
              << ", nEnergyBins=" << Mode3DForward::cDensity3D::nEnergyBins
              << "\n";
  }
#endif
}


bool IsForwardMode() {
#if _PIC_COUPLER_MODE_ != _PIC_COUPLER_MODE__SWMF_
  return false;
#else
  return s_pre_initialized && IsForwardTarget_(s_prm);
#endif
}

bool IsCutoffRigidityMode() {
#if _PIC_COUPLER_MODE_ != _PIC_COUPLER_MODE__SWMF_
  return false;
#else
  return s_pre_initialized && IsBackwardProductTarget_(s_prm);
#endif
}


double GetCoupledCalculationCadenceSeconds() {
#if _PIC_COUPLER_MODE_ != _PIC_COUPLER_MODE__SWMF_
  return 0.0;
#else
  // Keep the public cadence query bit-for-bit consistent with the value written
  // to PIC::ParticleWeightTimeStep::GlobalTimeStep[0] in amps_init().  Do not
  // duplicate the FIELD_UPDATE_DT parsing/default logic here; otherwise the PT
  // scheduler time step and the cutoff/density cadence gate could drift apart
  // after a future input-format edit.
  return BackwardProductCadenceSecondsFromPrm_(s_prm);
#endif
}

bool ReadyForBackwardProductCalculation(bool verbose) {
#if _PIC_COUPLER_MODE_ != _PIC_COUPLER_MODE__SWMF_
  return false;
#else
  int preInitLocal = s_pre_initialized ? 1 : 0;
  int meshLocal    = (PIC::Mesh::mesh != NULL) ? 1 : 0;
  int bOffsetLocal = (PIC::CPLR::SWMF::MagneticFieldOffset >= 0) ? 1 : 0;
  int vOffsetLocal = (PIC::CPLR::SWMF::BulkVelocityOffset >= 0) ? 1 : 0;
  int coupledLocal = PIC::CPLR::SWMF::FirstCouplingOccured ? 1 : 0;

  int preInitGlobal = 0;
  int meshGlobal    = 0;
  int bOffsetGlobal = 0;
  int vOffsetGlobal = 0;
  int coupledGlobal = 0;

  MPI_Allreduce((void*)&preInitLocal,&preInitGlobal,1,MPI_INT,MPI_MIN,
                MPI_GLOBAL_COMMUNICATOR);
  MPI_Allreduce((void*)&meshLocal,&meshGlobal,1,MPI_INT,MPI_MIN,
                MPI_GLOBAL_COMMUNICATOR);
  MPI_Allreduce((void*)&bOffsetLocal,&bOffsetGlobal,1,MPI_INT,MPI_MIN,
                MPI_GLOBAL_COMMUNICATOR);
  MPI_Allreduce((void*)&vOffsetLocal,&vOffsetGlobal,1,MPI_INT,MPI_MIN,
                MPI_GLOBAL_COMMUNICATOR);
  MPI_Allreduce((void*)&coupledLocal,&coupledGlobal,1,MPI_INT,MPI_MIN,
                MPI_GLOBAL_COMMUNICATOR);

  const bool ready =
      (preInitGlobal == 1) &&
      (meshGlobal    == 1) &&
      (bOffsetGlobal == 1) &&
      (vOffsetGlobal == 1) &&
      (coupledGlobal == 1);

  if (!ready) {
    static bool reportedWaiting = false;

    if (verbose && PIC::ThisThread == 0 && !reportedWaiting) {
      std::cout << "[Mode3DForwardSWMF] Waiting to calculate SWMF-coupled "
                << "Mode3D backward products until the first complete coupling "
                << "snapshot is available: "
                << "preInit=" << preInitGlobal
                << ", mesh=" << meshGlobal
                << ", magneticFieldOffset=" << bOffsetGlobal
                << ", bulkVelocityOffset=" << vOffsetGlobal
                << ", firstCoupling=" << coupledGlobal
                << ".\n";
      std::cout.flush();
    }

    reportedWaiting = true;
    return false;
  }

  return true;
#endif
}

// Public preparation hook called immediately before the SWMF-coupled cutoff solver.
//
// This routine intentionally performs all expensive globalisation work once per cutoff
// snapshot, not inside the particle mover.  After it returns, field access during the
// cutoff calculation is a local memory lookup/interpolation on every MPI rank.
void PrepareGlobalSWMFCoupledMagneticFieldForCutoff(bool verbose) {
#if _PIC_COUPLER_MODE_ != _PIC_COUPLER_MODE__SWMF_
  return;
#else
  if (!s_pre_initialized) {
    exit(__LINE__,__FILE__,
         "[Mode3DForwardSWMF] PrepareGlobalSWMFCoupledMagneticFieldForCutoff() called before amps_pre_init().");
  }

  if (PIC::Mesh::mesh==NULL) {
    exit(__LINE__,__FILE__,
         "[Mode3DForwardSWMF] PrepareGlobalSWMFCoupledMagneticFieldForCutoff() called before PIC::Mesh::mesh is initialized.");
  }

  if (PIC::CPLR::SWMF::MagneticFieldOffset<0) {
    exit(__LINE__,__FILE__,
         "[Mode3DForwardSWMF] PrepareGlobalSWMFCoupledMagneticFieldForCutoff() called before the SWMF magnetic-field buffer is allocated.");
  }

  if (PIC::CPLR::SWMF::FirstCouplingOccured==false) {
    exit(__LINE__,__FILE__,
         "[Mode3DForwardSWMF] PrepareGlobalSWMFCoupledMagneticFieldForCutoff() called before the first SWMF coupling data receive.");
  }

  if (PIC::CPLR::SWMF::BulkVelocityOffset<0) {
    exit(__LINE__,__FILE__,
         "[Mode3DForwardSWMF] SWMF bulk-velocity buffer is unavailable; "
         "cannot validate/export the required frozen B/u snapshot.");
  }

  ValidateLiveSWMFStateCoherence_();

  // Assemble compact fields from authoritative owner cells.  The globally replicated
  // AMR tree supplies geometry and node->Temp_ID indexing; row-stencil interpolation
  // later accesses these arrays without requiring remote block allocation.
  Earth::Mode3D::GlobalMagneticField::AssembleCellCenteredFieldsForCutoff(
      "Mode3DForwardSWMF",
      PIC::CPLR::SWMF::MagneticFieldOffset,
      -1,
      PIC::CPLR::SWMF::BulkVelocityOffset,
      BuildSWMFFieldSnapshotMetadata_(),
      verbose,
      CurrentCutoffOutputTimeSeconds_());
#endif
}

// Optional debug override.  It is deliberately not called by default in the coupled
// cutoff path: production SWMF-coupled cutoff now uses
// PrepareGlobalSWMFCoupledMagneticFieldForCutoff().  Call this manually when a known
// analytic dipole is desired for debugging the cutoff machinery independently of SWMF.
void RedefineSWMFCoupledMagneticFieldToAnalyticDipole() {
#if _PIC_COUPLER_MODE_ != _PIC_COUPLER_MODE__SWMF_
  return;
#else
  if (!s_pre_initialized) {
    exit(__LINE__,__FILE__,
         "[Mode3DForwardSWMF] RedefineSWMFCoupledMagneticFieldToAnalyticDipole() called before amps_pre_init().");
  }

  if (PIC::Mesh::mesh==NULL) {
    exit(__LINE__,__FILE__,
         "[Mode3DForwardSWMF] RedefineSWMFCoupledMagneticFieldToAnalyticDipole() called before PIC::Mesh::mesh is initialized.");
  }

  if (PIC::CPLR::SWMF::MagneticFieldOffset<0) {
    exit(__LINE__,__FILE__,
         "[Mode3DForwardSWMF] RedefineSWMFCoupledMagneticFieldToAnalyticDipole() called before the SWMF magnetic-field buffer is allocated.");
  }

  Earth::GridlessMode::Dipole::SetMomentScale(s_prm.field.dipoleMoment_Me);
  Earth::GridlessMode::Dipole::SetTiltDeg(s_prm.field.dipoleTilt_deg);

  // Replace only the compact global B array.  The helper evaluates cell-center
  // coordinates directly from tree geometry and does not allocate remote blocks.
  const long int nCells =
      Earth::Mode3D::GlobalMagneticField::RedefineGlobalMagneticField(
          "Mode3DForwardSWMF",
          AnalyticDipoleMagneticField_,
          false);

  if (PIC::ThisThread == 0) {
    std::cout << "[Mode3DForwardSWMF] Replaced SWMF-coupled cell-centered B field "
              << "with analytic dipole:"
              << " DIPOLE_MOMENT=" << s_prm.field.dipoleMoment_Me
              << ", DIPOLE_TILT=" << s_prm.field.dipoleTilt_deg << " deg"
              << ", updatedGlobalCells=" << nCells << ".\n";
    std::cout.flush();
  }
#endif
}

void amps_cutoff_time_step() {
#if _PIC_COUPLER_MODE_ != _PIC_COUPLER_MODE__SWMF_
  return;
#else
  if (!s_pre_initialized) {
    exit(__LINE__,__FILE__,
         "[Mode3DForwardSWMF] amps_cutoff_time_step() called before amps_pre_init().");
  }

  if (!IsBackwardProductTarget_(s_prm)) return;

  const long int callIndex = s_cutoff_call_counter;
  const double   tSim_s    = CurrentCutoffOutputTimeSeconds_();
  const std::string suffix = FormatCutoffOutputSuffix_(callIndex);

  // Save the suffix before the cutoff calculation starts.  The cutoff products
  // and the diagnostic coupled-mesh dump written by main_lib.cpp must carry the
  // same stamp, because they describe the same SWMF/PT coupled snapshot.
  s_last_cutoff_output_suffix = suffix;
  s_last_cutoff_output_suffix_valid = true;

  Earth::Mode3D::SetCutoffOutputFileSuffix(suffix);
  Earth::Mode3D::SetDensityOutputFileSuffix(suffix);

  if (PIC::ThisThread == 0) {
    std::cout << "[Mode3DForwardSWMF] cutoff snapshot " << callIndex
              << ": t_sim=" << tSim_s << " s, suffix='" << suffix << "'.\n";
    std::cout.flush();
  }

  const std::string exportedSnapshot=CoupledSnapshotExportFile_(suffix);
  const std::string statusFile=CoupledSnapshotStatusFile_(suffix);
  Earth::Field::SnapshotMetadata batchFieldMetadata;
  batchFieldMetadata.sourceId="PIC::CPLR:SWMF";
  batchFieldMetadata.modelName="SWMF";
  batchFieldMetadata.epochUTC=s_prm.field.epoch;

  try {
    // Fail closed before gathering or tracing.  Several lower-level AMPS consistency
    // guards terminate through exit() rather than C++ exceptions; prewriting FAILED
    // ensures such a termination still leaves an explicit non-passing artifact for
    // this requested coupling snapshot.  PASS is written only at the end below.
    WriteCoupledSnapshotStatus_(
        statusFile,"FAILED",batchFieldMetadata,suffix,exportedSnapshot,
        "snapshot calculation started but has not completed all requested products");

    // Replace the fail-safe reference metadata with the precise coupled UTC/domain/
    // mode record. The final content-derived ID is installed only after B/u assembly.
    batchFieldMetadata=BuildSWMFFieldSnapshotMetadata_();

    // The cutoff solver can run trajectories across the entire AMR domain on each
    // MPI rank.  Assemble the distributed SWMF fields into compact global arrays
    // before starting any backward tracing.
    PrepareGlobalSWMFCoupledMagneticFieldForCutoff(true);
    std::shared_ptr<Earth::Field::IFieldProvider> batchFieldProvider=
        Earth::Mode3D::GlobalMagneticField::CurrentFieldProvider();
    Earth::Field::SnapshotRequest batchFieldRequest;
    batchFieldRequest.epochUTC=
        Earth::Mode3D::GlobalMagneticField::CurrentSnapshotMetadata().epochUTC;
    batchFieldRequest.requestId=suffix;
    const std::shared_ptr<const Earth::Field::IFieldSnapshot> batchFieldSnapshot=
        batchFieldProvider->CreateSnapshot(batchFieldRequest);
    batchFieldMetadata=batchFieldSnapshot->Metadata();

    // Freeze before export or trajectories. SWMF serializes PT callbacks, naturally
    // queueing the next receive until this callback returns; the compact-array lease
    // also rejects any accidental reassembly initiated from inside this batch.
    Earth::Mode3D::GlobalMagneticField::BeginFrozenFieldBatch(
        batchFieldMetadata.snapshotId);

    // Export only after compact-array assembly and identity publication.  The writer
    // consumes the frozen arrays, not the mutable coupler buffers, and verifies their
    // content fingerprint against batchFieldMetadata before rank 0 writes the file.
    if (!exportedSnapshot.empty()) {
      const std::string exportedFingerprint=
          Earth::Mode3D::GlobalMagneticField::ExportCurrentSWMFSnapshot(
              exportedSnapshot);
      if (exportedFingerprint!=
          Earth::Mode3D::GlobalMagneticField::CurrentContentFingerprint())
        throw std::runtime_error(
            "exported SWMF snapshot fingerprint differs from the active field state");
    }
    // Use a value-owned parameter block at the absolute epoch published with the
    // frozen field.  Earlier code kept the reference epoch and added PT time only for
    // the boundary spectrum; that could synchronize the spectrum while leaving frame
    // rotations/ephemerides at the wrong time.  Step 9 binds every consumer to the
    // same actual UTC, so no secondary offset is applied.
    EarthUtil::AmpsParam batchProductPrm=s_prm;
    batchProductPrm.field.epoch=batchFieldMetadata.epochUTC;
    batchProductPrm.densitySpectrum.spectrumEpochOffsetActive=false;
    batchProductPrm.densitySpectrum.spectrumEpochOffset_s=0.0;

    // Run the products requested by CALC_TARGET. Both products share the compact SWMF
    // B snapshot (and explicitly experimental E, if enabled), so cutoff, directional
    // access, spectra, and flux are synchronized under one identity and suffix. The
    // released default remains magnetic-only; no test or observation gate is relaxed.
    if (TargetRequestsCutoff_(s_prm)) {
      Earth::Mode3D::RunCutoffRigidity(batchProductPrm,true);
      Earth::Field::RequireSameSnapshot(
          batchFieldMetadata,
          Earth::Mode3D::GlobalMagneticField::CurrentSnapshotMetadata(),
          "SWMF-coupled cutoff product");
    }
    if (TargetRequestsDensityFlux_(s_prm)) {
      Earth::Mode3D::RunDensityAndFlux(batchProductPrm);
      Earth::Field::RequireSameSnapshot(
          batchFieldMetadata,
          Earth::Mode3D::GlobalMagneticField::CurrentSnapshotMetadata(),
          "SWMF-coupled density/flux/spectrum product");
    }

    Earth::Mode3D::GlobalMagneticField::EndFrozenFieldBatch(
        batchFieldMetadata.snapshotId);

    WriteCoupledSnapshotStatus_(
        statusFile,"PASS",batchFieldMetadata,suffix,exportedSnapshot,
        "all requested products completed from one immutable SWMF B/u snapshot");
  }
  catch (const std::exception& e) {
    // FAILED was written collectively before field assembly. Do not start another MPI
    // collective here: a rank-local product exception may reach this catch while peers
    // are still inside the product, and a status broadcast would then deadlock. Rank
    // zero may enrich the already-failing artifact on a best-effort basis; failure to
    // do so never masks the original validation/physics error.
    if (PIC::ThisThread==0) {
      try {
        Earth::SWMFSnapshot::WriteProductStatus(
            statusFile,
            Earth::SWMFSnapshot::BuildProductStatusJson(
                "FAILED",batchFieldMetadata.snapshotId,batchFieldMetadata.epochUTC,
                suffix,exportedSnapshot,TargetRequestsCutoff_(s_prm),
                TargetRequestsDensityFlux_(s_prm),e.what()));
      }
      catch (const std::exception&) {}
    }
    std::ostringstream msg;
    msg << "[Mode3DForwardSWMF] SWMF-coupled backward-product calculation failed: "
        << e.what();
    const std::string text = msg.str();
    exit(__LINE__,__FILE__,text.c_str());
  }

  ++s_cutoff_call_counter;
#endif
}

std::string GetLastCutoffOutputFileName(const char* stem,const char* extension) {
#if _PIC_COUPLER_MODE_ != _PIC_COUPLER_MODE__SWMF_
  // This helper is meaningful only for SWMF-coupled PT runs.  In non-SWMF
  // builds, return the conventional unstamped name so accidental calls remain
  // harmless and backward compatible.
  return std::string(stem ? stem : "") + std::string(extension ? extension : "");
#else
  // main_lib.cpp calls this immediately after amps_cutoff_time_step().  At that
  // point s_last_cutoff_output_suffix should contain the exact suffix already
  // passed to Earth::Mode3D::SetCutoffOutputFileSuffix() for the cutoff products.
  // Reusing the cached suffix is intentional: it prevents the AMPS mesh dump from
  // observing a later TimeCounter value if the code path is ever extended to do
  // additional work between cutoff output and mesh output.
  const std::string safeStem      = (stem      != nullptr) ? stem      : "";
  const std::string safeExtension = (extension != nullptr) ? extension : "";

  if (!s_last_cutoff_output_suffix_valid) {
    // Fallback for defensive robustness.  The normal cutoff branch always sets
    // the suffix before this helper is used, but returning a valid file name is
    // safer than failing during diagnostic output.
    return safeStem + safeExtension;
  }

  return safeStem + s_last_cutoff_output_suffix + safeExtension;
#endif
}

} // namespace Mode3DForwardSWMF
} // namespace Earth
