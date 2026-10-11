#ifndef SRCSEP_UTIL_SEP_INITIALIZATION_H
#define SRCSEP_UTIL_SEP_INITIALIZATION_H

#include "sep_common_header_path.h"
#include SRCSEP_SEP_COMMON_HEADER(sep_transport_common.h)
#include "sep_production_mover.h"
#include "sep_turbulence_core.h"

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace SEP {
namespace Initialization {

// Small AMPS-independent vector used while parsing and validating startup
// input.  Keeping PIC and MPI out of this layer lets the complete input and
// mesh-resolution laws run before AMPS allocates global state.
struct Vec3 {
  double x = 0.0;
  double y = 0.0;
  double z = 0.0;
  Vec3() {}
  Vec3(double xValue, double yValue, double zValue)
      : x(xValue), y(yValue), z(zValue) {}
};

// Raw provider assignment transported to the canonical SWCME1D resolver.
// This startup layer deliberately does not interpret model keys or units: the
// model-owned resolver remains the single scientific authority.  Keeping line
// provenance here lets a malformed physical value identify its input record.
struct SwcmeAssignment {
  std::string key;
  std::string value;
  std::size_t line = 0;
};

// One record of the physics: parallel_diffusion block body, stored
// as "key = value" (comments already removed by the generic reader) because
// the shared parallel-diffusion library, not this reader, owns the block's
// model selection and key schema.  The line number is the 1-based physical
// line in the input file that holds the record.
struct ParallelDiffusionLine {
  std::size_t line = 0;
  std::string text;
};

enum class RefinementProfile { Linear, PowerLaw, Smoothstep };
enum class TubeRadiusMode { PhysicalConstant, ConstantAngularWidth };
enum class EnergyChannelSpacing { Logarithmic, Linear };

// A one-dimensional observer is identified by its radius along the embedded
// field line (shared observer kind field-line-radius, decision D5 of
// src/models/sep_common/sep_input_layout.h).  Each repeatable
// '#begin observer: ID' block owns a complete spectrum grid; observers
// therefore need not share energy bounds or bin spacing.
struct ObserverConfiguration {
  std::string id;
  double heliocentricRadiusM = 0.0;  // observer.ID.heliocentric_radius_m [m]
  double minimumEnergyJ = 0.0;       // kinetic energy [J]
  double maximumEnergyJ = 0.0;       // kinetic energy [J]
  unsigned energyChannels = 0;       // observer.ID.energy_bins
  EnergyChannelSpacing energySpacing = EnergyChannelSpacing::Logarithmic;
  unsigned pitchAngleBins = 0;
};

// SI-only startup record of srcSEP, read from the shared SEP content schema 6
// (src/INPUT_FILE.md section 5; vocabulary table
// src/models/sep_common/sep_input_layout.h).  The spatial mesh supports the
// embedded field line; particles themselves remain on the one-dimensional PIC
// field-line representation.  Every value is SI.
struct Configuration {
  // run.time_step_s and source.shock.samples_per_step: explicit numerical
  // choices, never derived from mesh resolution or an unrelated boundary
  // source.  samples_per_step is the macroparticle count of every active
  // field-line source event.
  double timeStepS = 0.0;
  std::uint64_t macroparticlesPerStep = 0;
  // run.maximum_time_steps (optional): the main-loop step budget.  Absent
  // means the historical driver default; --total-iterations overrides it.
  bool maximumTimeStepsSpecified = false;
  std::uint64_t maximumTimeSteps = 0;
  // species.macroparticle_weight: common base statistical weight installed
  // for every species in AMPS' compiled SpeciesList.  Species identity/mass/
  // charge never come from this post-compile file; injection applies any
  // species-specific correction via the already validated source table.
  double particleWeight = 0.0;

  // One-dimensional observers (kind field-line-radius), installed into the
  // retained field-line sampler by heliocentric radius.
  std::vector<ObserverConfiguration> observers;

  // output.initialization_{mesh,field_line,data}_tecplot_file.  Products are
  // written only after the corresponding AMPS mesh/field line exists.
  // Explicit paths avoid an undocumented working-directory convention and are
  // part of the startup fingerprint.
  std::string meshTecplotFile;
  std::string fieldLineTecplotFile;
  std::string dataTecplotFile;

  Vec3 parkerOriginM;
  Vec3 parkerInitialPointM;
  double parkerLengthM = 0.0;
  std::uint64_t parkerPointCount = 0;
  // Ambient wind speed [m/s] and solar rotation rate [rad/s] of the Parker
  // winding.  Decision D4: physics: swcme owns both, so they are not written
  // in the file; the parser copies them from the canonical SWCME1D resolver
  // (ambient.wind_speed and the fixed SWCME1D Carrington rate).
  double solarWindSpeedMPerS = 0.0;
  double solarRotationRateRadPerS = 0.0;

  double innerRadiusM = 0.0;
  double outerRadiusM = 0.0;
  double globalCellSizeM = 0.0;
  double minimumCellSizeM = 0.0;
  unsigned maximumMeshLevel = 0;

  bool solarRefinementEnabled = true;
  double solarSurfaceCellSizeM = 0.0;
  double solarTransitionOuterRadiusM = 0.0;
  RefinementProfile solarProfile = RefinementProfile::Smoothstep;
  double solarExponent = 1.0;

  bool tubeRefinementEnabled = true;
  double tubeReferenceRadiusM = 0.0;
  double tubeRadiusAtReferenceM = 0.0;
  TubeRadiusMode tubeRadiusMode = TubeRadiusMode::ConstantAngularWidth;
  double tubeCenterCellSizeM = 0.0;
  RefinementProfile tubeProfile = RefinementProfile::Smoothstep;
  double tubeExponent = 1.0;

  // Every canonical SWCME setting: the physics: swcme block plus the
  // source.* assignments synthesized from source: shock (decision D4).  The
  // application forwards this layer unchanged to SW1DAdapter::Configure().
  std::vector<SwcmeAssignment> swcmeAssignments;

  // application.sep.transport = parker | focused-diffusion |
  // focused-scattering selects the production mover parker | fte-dmumu |
  // fte-mfp unless --particle-mover is given (the command line keeps its
  // historical higher precedence; main.cpp applies the file value only when
  // the CLI is silent).  particleMoverLine is the 1-based input line.
  bool particleMoverSpecified = false;
  Mover::ProductionMover particleMover = Mover::ProductionMover::Parker;
  std::size_t particleMoverLine = 0;

  // application.sep.turbulence: srcSEP's self-consistent Alfven-wave
  // turbulence (decision D2).  Every value has a command-line option that
  // overrides it: main.cpp layers defaults < this record < command line by
  // re-parsing the command line on top of SEP::Util::CLI::ApplyInputFile-
  // Turbulence().  shockSourceEfficiency and shockSourcePlusFraction are the
  // fraction of shock energy given to waves and its E+ share;
  // diagnosticCadenceSteps is the spectrum-diagnostic cadence (0 = off).
  bool turbulenceSpecified = false;
  Turbulence::Configuration turbulence;
  double shockSourceEfficiency = 0.0;
  double shockSourcePlusFraction = 0.0;
  std::uint64_t diagnosticCadenceSteps = 0;

  // Optional physics: parallel_diffusion block.  This layer only records
  // that the section exists and its raw body; it does not interpret model
  // keys or units.  The srcSEP parallel-diffusion adapter
  // (adapters/parallel_diffusion_adapter.h) hands the body to the library
  // parser, cross-checks it against --spatial-diffusion-provider, installs
  // the selected model, and copies the resolved model ID and library
  // configuration fingerprint into the two strings below for the startup
  // fingerprint.
  bool hasParallelDiffusionSection = false;
  std::vector<ParallelDiffusionLine> parallelDiffusionLines;
  std::string parallelDiffusionModelId;
  std::string parallelDiffusionConfigurationFingerprint;
};

// Parse a complete srcSEP input written in the generic AMPS runtime-input
// format (src/INPUT_FILE.md) with the content schema 6 shared with srcSEP3D
// (application.sep.schema_version = 6; vocabulary table
// src/models/sep_common/sep_input_layout.h).  Unknown blocks and keys are
// errors, srcSEP3D-only settings are typed "not supported by srcSEP" errors,
// and no required value is silently inherited from a production default.
Transport::Status ParseText(const std::string& text, Configuration* result);
Transport::Status LoadFile(const std::string& path, Configuration* result);
Transport::Status Validate(const Configuration& configuration);

// Retarget all declared initialization products to one command-line-selected
// directory while preserving the leaf names declared in [output].  This is a
// pure configuration transform; the MPI-aware runtime creates parent
// directories immediately before collective output begins.
Transport::Status ApplyOutputDirectoryOverride(
    const std::string& directory, Configuration* configuration);

// Install exactly once before amps_init_mesh().  Absence means legacy mode and
// intentionally preserves every historical hard-coded srcSEP mesh choice.
Transport::Status Install(const Configuration& configuration);
bool HasActive();
const Configuration& Active();

// Stable configuration identity printed before AMPS initialization.
std::string Fingerprint(const Configuration& configuration);

// Generate an outward Parker curve at uniform arc-length spacing. pointCount
// includes both endpoints. The midpoint tangent method is second-order in the
// geometry while preserving each requested segment length to roundoff.
Transport::Status BuildParkerLine(const Configuration& configuration,
                                  std::vector<Vec3>* points);

// Write the finite one-dimensional transport mesh as a Tecplot ordered POINT
// zone.  The writer validates/builds the complete curve before opening the
// destination, so invalid physics never leaves a partial visualization file.
Transport::Status WriteParkerLineTecplot(
    const Configuration& configuration, const std::string& path);

double TubeRadiusM(double heliocentricRadiusM,
                   const Configuration& configuration);
double TubeDistanceM(const Vec3& positionM,
                     const Configuration& configuration);
double RequestedCellSizeM(const Vec3& positionM,
                          const Configuration& configuration);

const char* Name(RefinementProfile profile);
const char* Name(TubeRadiusMode mode);
const char* Name(EnergyChannelSpacing spacing);

}  // namespace Initialization
}  // namespace SEP

#endif  // SRCSEP_UTIL_SEP_INITIALIZATION_H
