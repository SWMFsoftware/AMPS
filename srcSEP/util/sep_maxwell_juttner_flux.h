#ifndef SRCSEP_UTIL_SEP_MAXWELL_JUTTNER_FLUX_H
#define SRCSEP_UTIL_SEP_MAXWELL_JUTTNER_FLUX_H

#include "sep_common_header_path.h"
#include SRCSEP_SEP_COMMON_HEADER(sep_transport_common.h)
#include <array>

namespace SEP { namespace Sources { namespace MaxwellJuttner {

// Exact SI defining constants. Temperature and density below are PROPER
// rest-frame quantities; an imported electron density or lab density is not an
// implicit substitute. This kernel does not decide a provider's density basis.
constexpr double SpeedOfLightMPerS=299792458.0;
constexpr double BoltzmannJPerK=1.380649e-23;
// Full-output events must retain their intended positive material crossing
// after SI momenta/gamma are encoded in binary64. This is a relative precision
// requirement on the thermal crossing, never a floor based on bulk speed.
constexpr double CrossingRelativeResolution=2.0e-12;

struct MaterialBoundaryInput {
  double properNumberDensityPerM3=0.0,properTemperatureK=0.0,restMassKg=0.0;
  std::array<double,3> bulkVelocityMPerS{{0.0,0.0,0.0}};
  std::array<double,3> boundaryVelocityMPerS{{0.0,0.0,0.0}};
  // Direction INTO the simulation, expressed in the simulation frame. Its
  // magnitude is immaterial: this API computes flux per physical unit area,
  // not a source count or a material-face area evolution law.
  std::array<double,3> inwardNormal{{1.0,0.0,0.0}};
};

enum class NormalizationMethod {ColdZero,BoundedAsymptotic,LibraryBessel};
struct FluxNormalization {
  double dimensionlessTemperature=0.0; // k_B*T/(m*c^2)
  double meanRestSpeedMPerS=0.0;
  double transformedNormalMagnitude=0.0;
  double numberFluxPerM2S=0.0; // per simulation-frame instantaneous area/time
  NormalizationMethod method=NormalizationMethod::ColdZero;
  // A mathematical relative truncation bound for the large-1/theta series.
  // It is zero for the library branch, which uses the C++ special function;
  // zero does not assert exact library/roundoff error. Independent quadrature
  // validation sets the supported numerical accuracy of that branch.
  double asymptoticRelativeRemainderBound=0.0;
};

struct FluxSample {
  // Keep gamma'-1 separately: adding one destroys cold kinetic energy below
  // machine epsilon. Momentum is full three-component momentum, not velocity.
  double restKineticGamma=0.0,restGamma=1.0,simulationGamma=1.0;
  std::array<double,3> restMomentumKgMPerS{{0.0,0.0,0.0}};
  std::array<double,3> simulationMomentumKgMPerS{{0.0,0.0,0.0}};
  double relativeInwardSpeedMPerS=0.0;
};

// The scalar distribution f and d^3p/E invariance give a material-plane flux
// F=n0*|a|*<v'>/4. The selected domain is exactly v_boundary=U. A nonmaterial
// plane has a different constant term/angular law and is explicitly rejected.
// Output replacement is atomic on every error, including numeric overflow.
Transport::Status EvaluateMaterialFlux(const MaterialBoundaryInput& input,
                                       FluxNormalization* output);

// Exact conditional FULL-PARTICLE crossing distribution. It is not a native
// guiding-center source: projecting a tilted-cap gyrophase-dependent flux onto
// (p_parallel,p_perp) requires a separate derived boundary matching law.
// Uniform slots are mixture choice, three exponentials, cosine and azimuth.
// All six slots must be strictly interior (0,1); shape two ignores the third
// exponential. This fixed draw budget is part of keyed-event reproducibility.
Transport::Status SampleMaterialFluxFromUniforms(const MaterialBoundaryInput& input,
    const std::array<double,6>& openUniforms,FluxSample* output);

// Stage the complete keyed generator (including its unused normal cache) and
// sample privately. Only a representable, valid result commits RNG and output;
// numerical failure must never become a tail rejection/redraw or consumed event.
Transport::Status SampleMaterialFlux(const MaterialBoundaryInput& input,
    Transport::KeyedRandomStream* random,FluxSample* output);

} } }
#endif
