#include "sep.h"

#include <cmath>
#include <ostream>

// srcSEP is the field-line transport application.  Reject an incompatible
// enclosing AMPS configuration at compile time, before a particle can be
// attached to an AMR mesh cell and bypass the production adapter contract.
#if _PIC_FIELD_LINE_MODE_ != _PIC_MODE_ON_
#error "srcSEP Step 5 requires PIC field-line mode"
#endif

#if _PIC_PARTICLE_LIST_ATTACHING_ != _PIC_PARTICLE_LIST_ATTACHING_FL_SEGMENT_
#error "srcSEP Step 5 requires field-line-segment particle attachment"
#endif

namespace {

SEP::Mover::ProductionMover g_selected_mover =
    SEP::Mover::ProductionMover::FocusedTransportDiffusion;
using MoverImplementation = int (*)(
    long int, double, cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>*);
MoverImplementation g_selected_implementation =
    SEP::ParticleMover_FocusedTransport_Dmumu;

MoverImplementation ImplementationFor(SEP::Mover::ProductionMover mover) {
  switch (mover) {
    case SEP::Mover::ProductionMover::Parker:
      return SEP::ParticleMover_Parker;
    case SEP::Mover::ProductionMover::FocusedTransportDiffusion:
      return SEP::ParticleMover_FocusedTransport_Dmumu;
    case SEP::Mover::ProductionMover::FocusedTransportMeanFreePath:
      return SEP::ParticleMover_FocusedTransport_EventDriven;
  }
  exit(__LINE__,__FILE__,"unknown production mover");
  return nullptr;
}

const char* PitchAngleDiffusionProviderName() {
  const SEP::Diffusion::fGetPitchAngleDiffusionCoefficient provider =
      SEP::Diffusion::GetPitchAngleDiffusionCoefficient;
  if (provider == static_cast<SEP::Diffusion::fGetPitchAngleDiffusionCoefficient>(
                      &SEP::Diffusion::Constant::GetPitchAngleDiffusionCoefficient))
    return "Constant::GetPitchAngleDiffusionCoefficient";
  if (provider == static_cast<SEP::Diffusion::fGetPitchAngleDiffusionCoefficient>(
                      &SEP::Diffusion::Roux2004AJ::GetPitchAngleDiffusionCoefficient))
    return "Roux2004AJ::GetPitchAngleDiffusionCoefficient";
  if (provider == static_cast<SEP::Diffusion::fGetPitchAngleDiffusionCoefficient>(
                      &SEP::Diffusion::Qin2013AJ::GetPitchAngleDiffusionCoefficient))
    return "Qin2013AJ::GetPitchAngleDiffusionCoefficient";
  if (provider == static_cast<SEP::Diffusion::fGetPitchAngleDiffusionCoefficient>(
                      &SEP::Diffusion::Borovokov_2019_ARXIV::GetPitchAngleDiffusionCoefficient))
    return "Borovokov_2019_ARXIV::GetPitchAngleDiffusionCoefficient";
  if (provider == static_cast<SEP::Diffusion::fGetPitchAngleDiffusionCoefficient>(
                      &SEP::Diffusion::Jokopii1966AJ::GetPitchAngleDiffusionCoefficient))
    return "Jokopii1966AJ::GetPitchAngleDiffusionCoefficient";
  if (provider == static_cast<SEP::Diffusion::fGetPitchAngleDiffusionCoefficient>(
                      &SEP::Diffusion::Florinskiy::GetPitchAngleDiffusionCoefficient))
    return "Florinskiy::GetPitchAngleDiffusionCoefficient";
  return "custom function pointer";
}

const char* MeanFreePathProviderName() {
  switch (SEP::Scattering::MeanFreePathMode) {
    case SEP::Scattering::MeanFreePathMode_QLT: return "QLT";
    case SEP::Scattering::MeanFreePathMode_QLT1: return "QLT1";
    case SEP::Scattering::MeanFreePathMode_Tenishev2005AIAA:
      return "Tenishev2005AIAA";
    case SEP::Scattering::MeanFreePathMode_Chen2024AA: return "Chen2024AA";
  }
  return "unknown";
}

}  // namespace

void SEP::Mover::SelectProductionMover(ProductionMover mover) {
  // Resolve through one audited, private mapping.  Step 14 removes the public
  // mutable function pointer so callers cannot bypass state validation or
  // manufacture a mover that is absent from the three-entry registry.
  g_selected_mover = mover;
  g_selected_implementation = ImplementationFor(mover);
}

SEP::Mover::ProductionMover SEP::Mover::CurrentProductionMover() {
  return g_selected_mover;
}

const SEP::Mover::MoverCapabilities& SEP::Mover::CurrentCapabilities() {
  return Describe(g_selected_mover).capabilities;
}

int SEP::Mover::DispatchProductionMover(
    long int ptr,
    double dtTotal,
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* startNode) {
  if (ptr < 0 || !std::isfinite(dtTotal) || dtTotal < 0.0) {
    exit(__LINE__,__FILE__,"invalid particle handle or mover time interval");
  }

  PIC::ParticleBuffer::byte* particle =
      PIC::ParticleBuffer::GetParticleDataPointer(ptr);
  if (particle == NULL) {
    exit(__LINE__,__FILE__,"production mover received a null particle record");
  }

  const int field_line_id = PIC::ParticleBuffer::GetFieldLineId(particle);
  const int species = PIC::ParticleBuffer::GetI(particle);
  const double field_line_coordinate =
      PIC::ParticleBuffer::GetFieldLineCoord(particle);
  const double v_parallel_m_s = PIC::ParticleBuffer::GetVParallel(particle);
  const double v_normal_m_s = PIC::ParticleBuffer::GetVNormal(particle);

  if (PIC::FieldLine::FieldLinesAll == NULL ||
      field_line_id < 0 || field_line_id >= PIC::FieldLine::nFieldLine ||
      species < 0 || species >= PIC::nTotalSpecies ||
      !std::isfinite(field_line_coordinate) ||
      !std::isfinite(v_parallel_m_s) || !std::isfinite(v_normal_m_s) ||
      PIC::FieldLine::FieldLinesAll[field_line_id].GetSegment(
          field_line_coordinate) == NULL) {
    exit(__LINE__,__FILE__,
         "particle violates production field-line representation contract");
  }

  return g_selected_implementation(ptr,dtTotal,startNode);
}

void SEP::Mover::PrintRuntimeConfiguration(std::ostream& out) {
  const Descriptor& descriptor = Describe(g_selected_mover);
  out << "SEP production mover configuration:\n"
      << "  canonical name:               " << descriptor.canonicalName << "\n"
      << "  field-line attachment:        required\n"
      << "  particle pitch-angle state:   "
      << (descriptor.capabilities.usesPitchAngleState ? "required" : "averaged")
      << "\n"
      << "  wave streaming source:        "
      << (descriptor.capabilities.accumulatesWaveStreaming ? "manager input" : "none")
      << "\n"
      << "  direct wave-state evolution:  "
      << (descriptor.capabilities.evolvesWaveStateDirectly ? "yes" : "no")
      << "\n"
      << "  coefficient contract:         "
      << CoefficientContractName(descriptor.capabilities.coefficientContract)
      << "\n";

  // Report the actual provider behind the selected contract.  Parker Dxx is
  // calculated from the active Dmumu provider in diffusion_dxx.cpp, whereas
  // fte-mfp chooses one of the explicit mean-free-path closures.
  switch (descriptor.capabilities.coefficientContract) {
    case CoefficientContract::SpatialDiffusion:
      out << "  coefficient provider:         Dxx integral of "
          << PitchAngleDiffusionProviderName() << "\n";
      break;
    case CoefficientContract::PitchAngleDiffusion:
      out << "  coefficient provider:         "
          << PitchAngleDiffusionProviderName() << "\n";
      break;
    case CoefficientContract::MeanFreePath:
      out << "  coefficient provider:         "
          << MeanFreePathProviderName() << "\n";
      break;
  }

  if (descriptor.capabilities.coefficientContract ==
          CoefficientContract::SpatialDiffusion ||
      descriptor.capabilities.coefficientContract ==
          CoefficientContract::PitchAngleDiffusion) {
    out << "  Dmumu derivative:              "
        << (SEP::Diffusion::PitchAngleDifferentialMode ==
                    SEP::Diffusion::PitchAngleDifferentialModeAnalytical
                ? "analytical" : "numerical")
        << "\n"
        << "  constant Dmumu [1/s]:          "
        << SEP::Diffusion::ConstPitchAngleDiffusionValue << "\n"
        << "  special-mu limiter:            "
        << (SEP::Diffusion::LimitSpecialMuPointsMode ==
                    SEP::Diffusion::LimitSpecialMuPointsModeOn
                ? "on" : "off")
        << "\n"
        << "  special-mu distance:           "
        << SEP::Diffusion::LimitSpecialMuPointsDistance << "\n";

    if (SEP::Diffusion::GetPitchAngleDiffusionCoefficient ==
        static_cast<SEP::Diffusion::fGetPitchAngleDiffusionCoefficient>(
            &SEP::Diffusion::Jokopii1966AJ::GetPitchAngleDiffusionCoefficient)) {
      out << "  Jokipii spectrum mode:         "
          << (SEP::Diffusion::Jokopii1966AJ::Mode ==
                      SEP::Diffusion::Jokopii1966AJ::_awsom
                  ? "AWSoM" : "fractional")
          << "\n"
          << "  Jokipii fraction/power:        "
          << SEP::Diffusion::Jokopii1966AJ::FractionValue << " / "
          << SEP::Diffusion::Jokopii1966AJ::FractionPowerIndex << "\n"
          << "  Jokipii k reference [1/m]:     "
          << SEP::Diffusion::Jokopii1966AJ::k_ref_min << " .. "
          << SEP::Diffusion::Jokopii1966AJ::k_ref_max << " at r="
          << SEP::Diffusion::Jokopii1966AJ::k_ref_R << " m\n";
    }
  } else {
    out << "  MFP lower-Larmor limiter:      "
        << (SEP::LimitMeanFreePath ? "on" : "off") << "\n";

    switch (SEP::Scattering::MeanFreePathMode) {
      case SEP::Scattering::MeanFreePathMode_QLT:
        out << "  QLT k range at 1 AU [1/m]:     "
            << QLT::k_min_1AU << " .. " << QLT::k_max_1AU << "\n";
        break;
      case SEP::Scattering::MeanFreePathMode_QLT1:
        out << "  QLT1 deltaB/B:                 "
            << SEP::Transport::Coefficient::ActiveConfiguration().
                   prescribedDeltaBOverB << "\n"
            << "  QLT1 correlation length [m]:   "
            << SEP::Transport::Coefficient::ActiveConfiguration().
                   correlationLengthAt1AuM << "\n";
        break;
      case SEP::Scattering::MeanFreePathMode_Tenishev2005AIAA:
        out << "  Tenishev lambda0 [m]:          "
            << SEP::Scattering::Tenishev2005AIAA::lambda0 << "\n"
            << "  Tenishev energy exponent:      "
            << SEP::Scattering::Tenishev2005AIAA::alpha << "\n"
            << "  Tenishev radial exponent:      "
            << SEP::Scattering::Tenishev2005AIAA::beta << "\n";
        break;
      case SEP::Scattering::MeanFreePathMode_Chen2024AA:
        out << "  Chen Dxx closure:              5.16e14 r_AU^1.17 E_keV^0.71 m2/s\n";
        break;
    }
  }

  const SEP::Transport::Coefficient::Configuration& coefficients =
      SEP::Transport::Coefficient::ActiveConfiguration();
  const SEP::Transport::NumericalTolerances& tolerances =
      SEP::Transport::ActiveNumericalTolerances();
  out << "  source/representation policy:  "
      << SEP::Transport::Coefficient::SourceName(coefficients.source)
      << " / "
      << SEP::Transport::Coefficient::PitchAngleName(coefficients.pitchAngle)
      << "\n"
      << "  amplitude/gap policy:          "
      << SEP::Transport::Coefficient::TurbulenceAmplitudePolicyName(
             coefficients.amplitudePolicy) << " / "
      << SEP::Transport::Coefficient::ResonanceGapPolicyName(
             coefficients.resonanceGapPolicy) << "\n"
      << "  adaptive quadrature abs/rel:   "
      << coefficients.spatialQuadrature.absoluteToleranceM2PerS << " / "
      << coefficients.spatialQuadrature.relativeTolerance << "\n"
      << "  coefficient fingerprint:      "
      << SEP::Transport::Coefficient::ConfigurationFingerprint(coefficients)
      << "\n"
      << "  mover geometry/stochastic:     "
      << tolerances.geometryFraction << " / "
      << tolerances.stochasticPitchRms << "\n"
      << "  mover deterministic tolerance: "
      << tolerances.deterministicRelativeTolerance << "\n"
      << "  mover cooling/focusing:        "
      << tolerances.coolingLogChange << " / "
      << tolerances.focusingPitchChange << "\n"
      << "  mover shock/min-step [s]:      "
      << tolerances.shockFraction << " / " << tolerances.minimumStepS
      << "\n";
}
