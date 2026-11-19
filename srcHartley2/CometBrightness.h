#ifndef _GENERIC_COMET_BRIGHTNESS_H_
#define _GENERIC_COMET_BRIGHTNESS_H_

#include <string>

namespace Comet {
namespace Brightness {

enum {MaxEmissionLines=16,MetadataStringLength=160};

namespace GridMode {
enum Type {Linear=0,Logarithmic=1};
}
namespace OpticalDepthDiagnostic {
enum Type {None=0,LineCenter=1};
}
namespace GFactor {
// Fixed preserves the historical input path.  The other values are resolved
// through the reusable core resonance-fluorescence module during startup.
enum Type {Fixed=0,HydrogenLymanAlphaIKL2020=1,
    HydroxylAX00SchleicherAhearn1988Unquenched=2,
    OxygenI1304DirectKVB2022GroundJ2=3,
    OxygenI1304DirectSUMER2001GroundJ2=4,
    OxygenI1304LymanBetaBowenSUMER2001GroundJ2=5};
}
namespace VelocityDomainPolicy {
enum Type {Abort=0,ExcludeAndReport=1};
}

// All values below are generated from explicit comet input records.  The
// source defaults disable output, preserving every existing Hartley/Halley
// input that does not opt in to brightness diagnostics.
extern bool Enabled;
extern double MinimumImpactParameter;
extern double MaximumImpactParameter;
extern long int ImpactParameterCount;
extern int ImpactParameterGridMode;
extern int EmissionLineCount;
extern int EmissionSpecies[];
extern double PhotonGFactor[];
extern int GFactorModel[];
extern double EffectivePhotonWavelength[];
// Full comet heliocentric velocity in the brightness Cartesian frame.  The
// source default is deliberately (0,0,0); an input adapter or another comet
// application may install a different value through the public API below.
extern double CometHeliocentricVelocity[3];
extern double HeliocentricRadialDirection[3];
extern double LineOfSightDirection[3];
extern double ImpactParameterDirection[3];
extern bool ObservationGeometryConfigured;
extern double BulkHeliocentricRadialVelocity;
extern double SolarLymanAlphaFlux1AU;
// 1-AU solar photon fluxes [photons m^-2 s^-1] required only by the O I
// 130.4-nm models: the complete O I triplet (direct models) and the
// 102.4-102.8 nm band containing H Lyman-beta (Bowen model).
extern double SolarOxygen1304TripletFlux1AU;
extern double SolarLymanBetaFlux1AU;
extern int OutOfRangeVelocityPolicy;
extern int TauDiagnosticMode[];
extern double LineCenterCrossSection[];
extern double TauFlagThreshold[];

bool Validate(std::string *error=NULL);

// Runtime geometry API. Directions are normalized by the setter and the LOS
// and impact directions must be orthogonal. The radial direction points away
// from the Sun; the LOS direction points from the comet toward the observer.
bool SetCometHeliocentricVelocity(const double velocity[3],
    std::string *error=NULL);
bool SetObservationGeometry(const double heliocentricRadialDirection[3],
    const double lineOfSightDirection[3],
    const double impactParameterDirection[3],std::string *error=NULL);

// Runtime comet receipts configure the same arrays that legacy compile-time
// inputs generate. ResetRuntimeConfiguration is legal only before sampler
// registration. SetEmissionMetadata retains transition/provenance strings in
// fixed compiled storage, so runtime selection cannot change memory layout.
bool ResetRuntimeConfiguration(std::string *error=NULL);
bool SetEmissionMetadata(int line,const char *transition,
    const char *gFactorProvenance,const char *tauProvenance,
    std::string *error=NULL);

// Register the particle-resolved sampler after radial blocks and cell measures
// exist. PIC calls SampleParticleResolvedBrightness() every sampling iteration
// and OutputSampledBrightness() at the same completed-buffer boundary used by
// ordinary pic.* output. The output callback performs the MPI reduction,
// writes Tecplot files, and clears the external accumulation buffers.
void InitializeSampling();
void SampleParticleResolvedBrightness();
void OutputSampledBrightness(int outputIndex);

} // namespace Brightness
} // namespace Comet

#endif
