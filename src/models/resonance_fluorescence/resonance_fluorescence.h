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
  HydroxylAX00SchleicherAhearn1988Unquenched=2,
  // Atomic oxygen 130.4-nm triplet (130.2 + 130.5 + 130.6 nm).  All three
  // models assume absorbers in the 2p4 3P2 ground level and return the rate
  // of triplet photons per atom; see README.md for the physical scope.
  //   direct pumping of 3s 3S1 at 130.217 nm, Killen et al. (2022) table;
  OxygenI1304DirectKVB2022GroundJ2=3,
  //   direct pumping at 130.217 nm, SOHO/SUMER disk-centre atlas profile;
  OxygenI1304DirectSUMER2001GroundJ2=4,
  //   solar Lyman-beta pumping of 3d 3D (102.576 nm) and the 1128.7-nm ->
  //   844.6-nm cascade into the triplet (Bowen fluorescence).
  OxygenI1304LymanBetaBowenSUMER2001GroundJ2=5
};

enum class Status {
  Ok=0,
  UnknownModel,
  InvalidHeliocentricDistance,
  InvalidRadialVelocity,
  RadialVelocityOutOfRange,
  MissingSolarLymanAlphaFlux,
  MissingSolarOxygen1304Flux,
  MissingSolarLymanBetaFlux,
  InvalidModelResult
};

struct Conditions {
  // Positive velocity is radially away from the Sun.
  double radialVelocityMPS;
  double heliocentricDistanceM;

  // Required by the IKL H model and ignored by the Lowell OH model.
  // The public unit is photons m^-2 s^-1 at 1 AU.
  double solarLymanAlphaFlux1AUPhotonM2S;

  // Required by the two direct O I models: the line-integrated solar photon
  // flux of the complete O I 130.4-nm triplet at 1 AU, photons m^-2 s^-1, on
  // the SORCE SOLSTICE absolute scale of Killen et al. (2022).
  double solarOxygen1304TripletFlux1AUPhotonM2S;

  // Required by the O I Bowen model: solar flux in 102.4-102.8 nm (H
  // Lyman-beta and its wings) at 1 AU, photons m^-2 s^-1.
  double solarLymanBetaFlux1AUPhotonM2S;
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
  bool requiresSolarOxygen1304Flux;
  bool requiresSolarLymanBetaFlux;
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
