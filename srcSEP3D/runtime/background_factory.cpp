#include "background_factory.h"
#include "../background/bg_parker.h"
#include "../background/bg_swcme.h"
#include "swcme3d_input.hpp"
#include <map>
#include <mutex>
namespace SEP3D { namespace RuntimeModel {
namespace {
// Function-local statics avoid cross-translation-unit initialization order.
// The registry is process-local: MPI ranks must register the same models
// before acquisition, while the AMPS publication gate checks model identity.
std::map<std::string,BackgroundFactory>& Registry() {
  static std::map<std::string,BackgroundFactory> registry;return registry;
}
std::mutex& RegistryMutex() { static std::mutex mutex;return mutex; }
Core::Status Invalid(const std::string& text) {
  return Core::Status(Core::StatusCode::ConfigurationConflict,text);
}
// Translate every frozen ambient option, including composition and frame.
// SWCME and the inner continuation must share these values; copying only B/U
// would hide density, pressure or rotation differences at the handoff.
Background::ParkerConfiguration Parker(const RunConfiguration3D& configuration) {
  Background::ParkerConfiguration parker;
  const RuntimeModel::ParkerPhysicsOptions& configured =
    configuration.options().parker;
  parker.sourceRadiusM = configured.sourceRadiusM;
  parker.sourceLongitudeRad = configured.sourceLongitudeRad;
  parker.sourceColatitudeRad = configured.sourceColatitudeRad;
  parker.referenceRadiusM = configured.referenceRadiusM;
  parker.radialFieldAtReferenceT = configured.radialFieldAtReferenceT;
  parker.numberDensityAtReferenceM3 = configured.numberDensityAtReferenceM3;
  parker.densityReferenceRadiusM = configured.densityReferenceRadiusM;
  parker.temperatureK = configured.temperatureK;
  parker.adiabaticIndex = configured.adiabaticIndex;
  parker.thermodynamicClosure =
    configured.thermodynamicClosure ==
        RuntimeModel::SolarWindThermodynamicClosure::MultiSpecies
      ? swcme::solarwind::ThermodynamicClosure::MultiSpecies
      : swcme::solarwind::ThermodynamicClosure::ProtonOnly;
  parker.alphaToProtonRatio = configured.alphaToProtonRatio;
  parker.electronTemperatureK = configured.electronTemperatureK;
  parker.alphaTemperatureK = configured.alphaTemperatureK;
  parker.referenceSinColatitude = configured.referenceSinColatitude;
  parker.solarWindSpeedMPerS = configured.solarWindSpeedMPerS;
  parker.solarRotationRateRadPerS = configured.solarRotationRateRadPerS;
  parker.rotationAxis = configured.rotationAxis;
  parker.magneticPolarity = configured.magneticPolarity;
  parker.validityCadenceS = configured.validityCadenceS;
  parker.coordinateFrame = configured.coordinateFrame;
  return parker;
}
}
bool BackgroundAuthorityMatches(BackgroundAuthority a,Background::ProviderKind p) {
  return (a==BackgroundAuthority::AnalyticParker && p==Background::ProviderKind::AnalyticParker) ||
    (a==BackgroundAuthority::Swcme && p==Background::ProviderKind::Swcme) ||
    (a==BackgroundAuthority::Swmf && p==Background::ProviderKind::SwmfAwsom) ||
    (a==BackgroundAuthority::RuntimeModel && p==Background::ProviderKind::RuntimeModel);
}
Core::Status RegisterBackgroundModel(const std::string& id,BackgroundFactory factory) {
  if(id.empty()||!factory||id=="analytic-parker"||id=="swcme"||id=="swmf")
    return Invalid("background factory name is empty, reserved, or callback is null");
  std::lock_guard<std::mutex> lock(RegistryMutex());
  if(!Registry().emplace(id,std::move(factory)).second)
    return Invalid("background model already registered: "+id);
  return Core::Status::OK();
}
Core::Status CreateBackgroundProvider(const RunConfiguration3D& configuration,
    std::shared_ptr<Background::BackgroundProvider>* output) {
  if(!output)return Invalid("background provider output is null");
  const auto& options=configuration.options();
  std::shared_ptr<Background::BackgroundProvider> candidate;
  Core::Status status;
  try {
    if(options.background==BackgroundAuthority::AnalyticParker)
      candidate.reset(new Background::AnalyticParkerProvider(Parker(configuration)));
    else if(options.background==BackgroundAuthority::Swcme) {
      // Re-run the canonical resolver on retained assignments rather than
      // reconstructing a second set of CME parameters from application defaults.
      // Both the normalized manifest and fingerprint must match the shock's
      // frozen identity before this provider can sample a single mesh point.
      std::vector<swcme::input3d::Assignment> assignments;
      for(const auto& raw:options.swcmeAssignments) {
        swcme::input3d::Assignment a;a.key=raw.key;a.value=raw.value;
        a.line=raw.line;a.origin="frozen AMPS background";assignments.push_back(a);
      }
      const auto resolved=swcme::input3d::Resolve(assignments);
      if(!resolved.ok())return Invalid("cannot resolve SWCME background: "+resolved.status.message);
      if(resolved.configuration.fingerprint!=options.swcmeConfigurationFingerprint ||
         resolved.configuration.normalized_manifest!=options.swcmeResolvedManifest)
        return Invalid("SWCME background identity differs from frozen shock configuration");
      candidate.reset(new Background::SwcmeBackgroundProvider(resolved.configuration,
          Parker(configuration),options.requestedTimeStepS,
          options.requestedTimeStepS*options.backgroundCadenceSteps,options.coordinateOriginM));
    } else if(options.background==BackgroundAuthority::RuntimeModel) {
      BackgroundFactory factory;
      // Copy the callable while locked so registration cannot invalidate the
      // lookup. Construct outside the lock to permit dependency construction
      // and avoid holding a global registry lock through expensive model work.
      { std::lock_guard<std::mutex> lock(RegistryMutex());
        const auto i=Registry().find(options.backgroundModelId);
        if(i==Registry().end())return Invalid("unregistered background model: "+options.backgroundModelId);
        factory=i->second;
      }
      // Factories may return a typed failure; successful null results are also
      // rejected below. Neither path replaces the caller's existing provider.
      status=factory(configuration,&candidate);if(!status.ok())return status;
    } else return Invalid("this authority requires a host-installed snapshot");
    if(!candidate)return Invalid("background factory returned a null provider");
    status=candidate->Validate();if(!status.ok())return status;
  } catch(const std::exception& e) { return Invalid(std::string("background factory: ")+e.what()); }
  // This assignment is the only publication point in the factory. Preparing
  // epochs and publishing mesh fields remain separate lifecycle operations.
  *output=std::move(candidate);return Core::Status::OK();
}
} }
