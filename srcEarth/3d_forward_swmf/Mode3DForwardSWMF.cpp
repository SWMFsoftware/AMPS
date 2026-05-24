//======================================================================================
// Mode3DForwardSWMF.cpp
//======================================================================================
//
// SWMF-coupled initialization bridge for the Earth 3-D forward energetic-particle
// branch.
//
// DESIGN
// ------
// This file deliberately exposes only Earth::Mode3DForwardSWMF::amps_init().  The
// function is a middleman between the SWMF-coupled AMPS lifecycle and the shared
// srcEarth/3d_forward implementation.  It does not call a separate
// coupled-runtime wrapper.  Instead, it performs the coupled initialization
// directly and reuses the existing 3d_forward helpers/callbacks for the actual work:
//
//   * Earth3DForward::SetMoverByName() / ConfigureFieldEvaluation()
//   * Mode3DForward::BoundaryInjectionSourceRate()
//   * Mode3DForward::InjectParticles()
//   * Mode3DForward::BuildForwardInjectionSpectrumMap()
//   * Mode3DForward::EvaluateTimeStep()
//   * Mode3DForward::InitAbsorptionSphere()
//   * Mode3DForward::InitBoundaryInjectionTable()
//   * cDensity3D and cSphereFlux3D samplers
//
// In a coupled SWMF run the standalone command line is not available to this AMPS
// component.  Therefore the input file name is fixed here to AMPS_PARAM.in.
//
// The only intended physical difference from standalone 3d_forward is the source of
// background magnetic/electric fields.  In _PIC_COUPLER_MODE__SWMF_ builds, the
// forward movers obtain B and E from PIC::CPLR, which dispatches to the live SWMF
// coupler state, rather than from standalone analytic or DATAFILE-populated fields.
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

void amps_init() {
#if _PIC_COUPLER_MODE_ != _PIC_COUPLER_MODE__SWMF_
  // This function is meaningful only in the SWMF-coupled build.  Keep a compiled
  // no-op in other builds so the source can still be included unconditionally.
  return;
#else
  static bool initialized = false;
  if (initialized) {
    if (PIC::ThisThread == 0) {
      std::cout << "[Mode3DForwardSWMF] Coupled 3d_forward runtime was already "
                   "initialized; skipping duplicate initialization.\n";
    }
    return;
  }
  initialized = true;

  // In standalone AMPS runs the input file is normally selected through the command
  // line.  In a coupled SWMF run that argument line is not available to this AMPS
  // component, so the coupled middleman must use a fixed conventional input name.
  const char* inputFile = "AMPS_PARAM.in";

  if (PIC::ThisThread == 0) {
    std::cout << "[Mode3DForwardSWMF] Reading coupled-run parameters from '"
              << inputFile << "'.\n";
  }

  EarthUtil::AmpsParam prm = EarthUtil::ParseAmpsParamFile(inputFile);

  // Keep the parameter object self-describing for diagnostics and for the shared
  // 3d_forward mover configuration.  The actual field source is enforced by the
  // SWMF build branch in ForwardParticleMovers.cpp/ElectricField.cpp: B and E are
  // obtained from PIC::CPLR, not from the standalone DATAFILE/Tsyganenko buffers.
  prm.field.model = "SWMF";
  prm.calc.fieldEvalMethod = "SWMF";
  prm.mode3d.forceAnalyticMagneticField = false;

  // Route AMPS particle motion through Earth::ParticleMover() ->
  // Earth3DForward::MoverManager().  This follows standalone 3d_forward and lets
  // the normal AMPS time-step call advance particles with the selected forward mover.
  Earth::ModelMode = Earth::BoundaryInjectionMode;

  if (!Earth::Earth3DForward::SetMoverByName(prm.mode3dForward.particleMover)) {
    throw std::runtime_error(
        "Unknown 3d_forward_swmf particle mover '" +
        prm.mode3dForward.particleMover +
        "'. Valid movers are BORIS, RK4, GC/GC4, and HYBRID.");
  }

  Earth::Earth3DForward::ConfigureFieldEvaluation(prm);

  Earth::Mode3DForward::sInjectionEnergyDistribution =
      Earth::Mode3DForward::ParseInjectionEnergyDistribution(prm.mode3dForward.injectionEnergyDistribution);

#if _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ != _INDIVIDUAL_PARTICLE_WEIGHT_ON_
  if (Earth::Mode3DForward::sInjectionEnergyDistribution == Earth::Mode3DForward::InjectionEnergyDistribution::LOG_UNIFORM) {
    throw std::runtime_error(
        "3d_forward_swmf LOG_UNIFORM injection requires AMPS to be compiled with "
        "_INDIVIDUAL_PARTICLE_WEIGHT_MODE_ == _INDIVIDUAL_PARTICLE_WEIGHT_ON_.");
  }
#endif

  Earth::Mode3DForward::sInitializeParticleTrajectories =
      prm.mode3dForward.initializeParticleTrajectories;
  Earth::Mode3DForward::sMaxParticleTrajectories = prm.mode3dForward.nParticleTrajectories;

  // Register the 3d_forward source-rate and particle-injection callbacks directly.
  // These replace any older BoundaryInjectionMode callbacks installed earlier by the
  // generic Earth initialization path.
  PIC::ParticleWeightTimeStep::UserDefinedExtraSourceRate = Earth::Mode3DForward::BoundaryInjectionSourceRate;
  PIC::BC::UserDefinedParticleInjectionFunction = Earth::Mode3DForward::InjectParticles;

  // In an SWMF-coupled build this helper intentionally does not initialize T96/T05/TA16
  // or fill DATAFILE field buffers.  It leaves the field source neutral so the actual
  // field access goes through PIC::CPLR/SWMF.
  Earth::Mode3DForward::ConfigureBackgroundFieldModel(prm);

  // Keep legacy boundary-injection diagnostics well-defined.  The current forward
  // injector samples the inward cosine-weighted hemisphere relative to each domain
  // face and does not use this vector to construct particle velocities, but some
  // generic messages still print it.
  Earth::BoundingBoxInjection::b[0] = 1.0;
  Earth::BoundingBoxInjection::b[1] = 0.0;
  Earth::BoundingBoxInjection::b[2] = 0.0;
  Earth::BoundingBoxInjection::SetPrm(prm);

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

  // Apply #DENSITY_3D settings before local density/sphere-flux initialization.
  // Note: the per-cell sampling buffer allocation is still controlled by the AMPS
  // mesh lifecycle.  This call keeps all runtime/output routines synchronized with
  // the coupled-run input file and is the latest safe point available from
  // main_lib.cpp::amps_init().
  Mode3DForward::cDensity3D::ConfigureEnergyGrid(prm);

  // Enable the AMPS per-particle sampling dispatcher and register the 3d_forward
  // density sampler, avoiding duplicate registrations if the generic startup path
  // already inserted the callback.
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
  const auto forwardSpectrum = Earth::Mode3DForward::BuildForwardInjectionSpectrumMap(prm);
  InitGlobalSpectrumFromKeyValueMap(forwardSpectrum);

  Earth::Mode3DForward::sSpecies = 0;
  Earth::Mode3DForward::sDt = Earth::Mode3DForward::EvaluateTimeStep(prm);
  PIC::ParticleWeightTimeStep::GlobalTimeStep[Earth::Mode3DForward::sSpecies] = Earth::Mode3DForward::sDt;

  Earth::Mode3DForward::sNParticlesPerIter = prm.mode3dForward.nParticlesPerIter;
  PIC::ParticleWeightTimeStep::maxReferenceInjectedParticleNumber =
      Earth::Mode3DForward::sNParticlesPerIter;
  PIC::ParticleWeightTimeStep::initParticleWeight_ConstantWeight(Earth::Mode3DForward::sSpecies);
  Earth::Mode3DForward::sParticleWeight = PIC::ParticleWeightTimeStep::GlobalParticleWeight[Earth::Mode3DForward::sSpecies];

  Earth::Mode3DForward::InitAbsorptionSphere(prm);
  Mode3DForward::cDensity3D::Init(prm);
  Mode3DForward::cSphereFlux3D::Init(prm, Earth::Mode3DForward::sAbsorptionSphere, Earth::Mode3DForward::sDt);
  Earth::Mode3DForward::InitBoundaryInjectionTable();

  PIC::Mover::BackwardTimeIntegrationMode = _PIC_MODE_OFF_;
  PIC::SamplingMode = _RESTART_SAMPLING_MODE_;

  if (PIC::ThisThread == 0) {
    std::cout << "[Mode3DForwardSWMF] Initialized SWMF-coupled 3d_forward runtime. "
              << "Field source=PIC::CPLR/SWMF, mover="
              << Earth::Earth3DForward::GetMoverName()
              << ", nParticlesPerIter=" << Earth::Mode3DForward::sNParticlesPerIter
              << ", particleWeight=" << Earth::Mode3DForward::sParticleWeight
              << ", dt=" << Earth::Mode3DForward::sDt
              << ", energySampling="
              << Earth::Mode3DForward::InjectionEnergyDistributionName(Earth::Mode3DForward::sInjectionEnergyDistribution)
              << "\n";
  }
#endif
}

} // namespace Mode3DForwardSWMF
} // namespace Earth
