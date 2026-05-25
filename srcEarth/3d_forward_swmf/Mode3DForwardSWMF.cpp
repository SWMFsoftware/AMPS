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
// source of background B/E fields: in _PIC_COUPLER_MODE__SWMF_ builds they come
// from PIC::CPLR rather than Tsyganenko/DATAFILE.
//
//======================================================================================

#include "Mode3DForwardSWMF.h"

#include "../3d_forward/Mode3DForward.h"
#include "../3d_forward/ForwardParticleMovers.h"
#include "../3d_forward/Density3D.h"
#include "../3d_forward/SphereFlux3D.h"
#include "../boundary/spectrum.h"
#include "../util/amps_param_parser.h"
#include "../Earth.h"

#include "pic.h"

#include <iostream>
#include <stdexcept>

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
#endif


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
  // PIC::CPLR instead of Tsyganenko/DATAFILE in all shared 3d_forward helpers.
  s_prm.field.model                       = "SWMF";
  s_prm.calc.fieldEvalMethod              = "SWMF";
  s_prm.mode3d.forceAnalyticMagneticField = false;

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

} // namespace Mode3DForwardSWMF
} // namespace Earth
