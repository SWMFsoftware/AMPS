#ifndef SRCSEP_ADAPTERS_REDUCED_SHOCK_BACKGROUND_ADAPTER_H
#define SRCSEP_ADAPTERS_REDUCED_SHOCK_BACKGROUND_ADAPTER_H

#include <array>
#include <cstdint>
#include <memory>
#include <string>

namespace SEP { namespace CoronaSwcme { namespace ShockFront {
class Provider;
struct Epoch;
} } }

namespace SEP { namespace ReducedShock {

// AMPS-independent SI sample returned by the shared reduced provider.  This
// type is intentionally small: srcSEP owns the conversion into field-line
// vertex storage, while all coronal/Parker plasma and IMF equations remain in
// src/models/sep_corona_swcme.  ``numberDensityM3`` is electron number density
// and ``pressurePa`` is the provider's total thermal pressure.
struct AmbientSample {
  std::array<double,3> magneticFieldT{{0.0,0.0,0.0}};
  std::array<double,3> velocityMPerS{{0.0,0.0,0.0}};
  double numberDensityM3 = 0.0;
  double pressurePa = 0.0;
  double protonTemperatureK = 0.0;
};

// Metadata that binds native field-line storage to one immutable provider
// epoch.  Generation follows the event's declared background cadence rather
// than the number of calls made by the application.
struct EpochMetadata {
  double epochS = 0.0;
  double validUntilS = 0.0;
  std::uint64_t generation = 0;
  std::string eventIdentity;
};

struct FrontSummary {
  double apexRadiusM = 0.0;
  double acceptedShockAreaM2 = 0.0;
  bool apexShockAccepted = false;
};

// Configure the process-wide srcSEP adapter from one event file and its
// repository-relative assets.  Configuration is transactional: a failed
// candidate does not replace an already configured provider.  The selected
// event must be HCI and particle_mode=disabled because this coupling stage is
// background-only.
bool Configure(const std::string& eventPath, std::string* error);
bool Enabled();

// Prepare a complete front/ambient epoch and commit it only after all shared
// geometry, ambient, and RH checks succeed.  A repeated call at the same epoch
// is idempotent; time reversal and two different epochs mapping to one cadence
// generation are rejected because either would make snapshot provenance
// ambiguous.
bool Prepare(double epochS, std::string* error);

// Query the undisturbed ambient authority at an HCI Cartesian position.  The
// reduced model deliberately has no spatial downstream CME volume, so this
// routine never paints a local RH jump behind the front into the field line.
bool EvaluateAmbient(const std::array<double,3>& positionM,
                     AmbientSample* sample, std::string* error);

const EpochMetadata& Metadata();
double BackgroundCadenceS();
const std::string& EventPath();
const std::string& EventIdentity();
std::shared_ptr<const CoronaSwcme::ShockFront::Epoch> FrontEpoch();
FrontSummary CurrentFrontSummary();

} } // namespace SEP::ReducedShock

#endif
