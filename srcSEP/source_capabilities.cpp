#include "source_capabilities.h"
#include "sep.h"
#include <cmath>

namespace SEP { namespace Sources {
Transport::Status CartesianBirthPreflight(CartesianSource source) {
  const char* name=nullptr;
  switch(source) {
    case CartesianSource::GalacticCosmicRayBox:name="legacy Cartesian GCR box source";break;
    case CartesianSource::ThermalOuterBox:name="legacy Cartesian outer Maxwellian source";break;
    case CartesianSource::SphericalInnerBoundary:
      return Transport::Status::Error(Transport::StatusCode::UnsupportedConfiguration,
          "legacy rank-local spherical source has no declared surface-to-field-line normalization and persisted collective birth scheduler; explicit field-line source semantics are required before enabling it");
    case CartesianSource::SphericalShock:name="legacy Cartesian spherical shock source";break;
  }
  if(!name)return Transport::Status::Error(Transport::StatusCode::InvalidArgument,
      "unknown legacy source capability");
  return Transport::Status::Error(Transport::StatusCode::UnsupportedConfiguration,
      std::string(name)+" has no declared physical field-line attachment and persisted collective birth scheduler; three-dimensional source capability belongs in srcSEP3D (see srcSEP3D/MIGRATION_MANIFEST.md); no sampler port is implied");
}
Transport::Status MaxwellianFieldLinePreflight() {
  return Transport::Status::Error(Transport::StatusCode::UnsupportedConfiguration,
      "SEP beginning source has not implemented the selected drifting relativistic Maxwell-Juttner flux with explicit proper-density/temperature and bulk/boundary velocity/normal inputs, exact local-plasma conversion and purpose-keyed sampling; available birth-window bookkeeping does not qualify the retained generic classical pre-drift helper");
}
Transport::Status ThermalPopulationPreflight() {
  // New source-time metadata cannot repair a distribution/weighting law. In
  // particular do not introduce an unapproved speed cutoff to make a classical
  // unbounded Maxwellian look like a relativistic source.
  return Transport::Status::Error(Transport::StatusCode::UnsupportedConfiguration,
      "optional thermal solar-wind population has not implemented the selected relativistic Maxwell-Juttner distribution, explicit proper-density/temperature basis, normalization and purpose-keyed sampler; birth metadata cannot qualify the retained classical helper");
}
Transport::Status ThermalBeginningPreflight() {
  return Transport::Status::Error(Transport::StatusCode::UnsupportedConfiguration,
      "optional thermal solar-wind beginning source has not implemented the selected drifting relativistic Maxwell-Juttner incoming-flux law with explicit proper density/temperature, bulk/boundary velocity, normal and purpose-keyed sampling; the retained classical helper is unsupported");
}
Transport::Status MovingShockDeparturePreflight(double departureEpochS,
                                               double trajectoryQueryEpochS) {
  if(!std::isfinite(departureEpochS)||!std::isfinite(trajectoryQueryEpochS))
    return Transport::Status::Error(Transport::StatusCode::InvalidArgument,
        "moving-shock source departure/query epochs must be finite [s]");
  if(departureEpochS!=trajectoryQueryEpochS)
    return Transport::Status::Error(Transport::StatusCode::UnsupportedConfiguration,
        "analytic moving-shock source view epoch differs from its frozen step-start birth; an explicit matching departure view is required before event reservation or particle creation");
  return Transport::Status::Ok();
}
Transport::Status ValidateConfiguredSources() {
#if _DOMAIN_GEOMETRY_ == _DOMAIN_GEOMETRY_BOX_ && _PIC_EARTH_GCR__MODE_ == _PIC_MODE_ON_
  if(SEP::BoundingBoxInjection::BoundaryInjectionMode)
    return CartesianBirthPreflight(CartesianSource::GalacticCosmicRayBox);
#endif
  if(::SEP::FieldLine::InjectionParameters::InjectionMomentumModel==
      ::SEP::FieldLine::InjectionParameters::_background_sw_temperature)
    return MaxwellianFieldLinePreflight();
  return Transport::Status::Ok();
}
} }
