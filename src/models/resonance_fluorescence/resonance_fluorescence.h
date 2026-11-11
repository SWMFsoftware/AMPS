#ifndef _AMPS_RESONANCE_FLUORESCENCE_H_
#define _AMPS_RESONANCE_FLUORESCENCE_H_

#include <cstddef>

// Geometry-, PIC-, MPI-, and application-independent resonance-fluorescence
// models.  Public inputs are SI.  Model names identify a transition and a
// published physical treatment; a species name alone is never sufficient to
// select a g factor.
namespace AMPS {
namespace ResonanceFluorescence {

enum class Model {
  HydrogenLymanAlphaIKL2020=1,
  HydroxylAX00SchleicherAhearn1988Unquenched=2
};

enum class Status {
  Ok=0,
  UnknownModel,
  InvalidHeliocentricDistance,
  InvalidRadialVelocity,
  RadialVelocityOutOfRange,
  MissingSolarLymanAlphaFlux,
  InvalidModelResult
};

struct Conditions {
  // Positive velocity is radially away from the Sun.
  double radialVelocityMPS;
  double heliocentricDistanceM;

  // Required by the IKL H model and ignored by the Lowell OH model.
  // The public unit is photons m^-2 s^-1 at 1 AU.
  double solarLymanAlphaFlux1AUPhotonM2S;
};

struct Result {
  Status status;

  // NaN means that the selected source does not define the quantity without
  // an additional convention.  In particular, Lowell supplies an integrated
  // OH-band energy efficiency; converting it to photons requires an explicit
  // effective wavelength selected by the application.
  double photonRateS1;
  double emittedPowerW;
  double radiationPressureRatio;
};

struct ModelMetadata {
  Model model;
  const char *id;
  const char *species;
  const char *transition;
  const char *citation;
  const char *sourceUrl;
  const char *referenceDataSha256;
  double referenceDistanceM;
  double minimumRadialVelocityMPS;
  double maximumRadialVelocityMPS;
  bool requiresSolarLymanAlphaFlux;
  bool photonRateAvailable;
  bool emittedPowerAvailable;
  bool radiationPressureRatioAvailable;
};

// Strict name conversion used by application input adapters.  Unknown names
// fail rather than silently selecting a default model.
bool ModelFromName(const char *name,Model *model);
const char *ModelName(Model model);
const char *StatusMessage(Status status);
const ModelMetadata *GetMetadata(Model model);

// Evaluate the selected model and apply heliocentric inverse-square scaling
// exactly once.  No velocity clamping or extrapolation is performed.
Result Evaluate(Model model,const Conditions &conditions);

// Validate immutable compiled reference tables and metadata.  Applications
// should call this once during startup; focused tests exercise it directly.
bool ValidateReferenceData(const char **error=NULL);

} // namespace ResonanceFluorescence
} // namespace AMPS

#endif
