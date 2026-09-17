#include "sep_species_source.h"

#include <cmath>
#include <set>

namespace SEP {
namespace Transport {
namespace SpeciesSource {
namespace {

std::vector<Configuration> gActiveConfiguration;

Status ValidateEntry(const Configuration& entry) {
  const Status speciesStatus =
      CoefficientPhysics::ValidateSpecies(entry.species);
  if (!speciesStatus.ok()) return speciesStatus;
  if (!std::isfinite(entry.abundanceFraction) ||
      entry.abundanceFraction < 0.0 ||
      !std::isfinite(entry.injectionEfficiency) ||
      entry.injectionEfficiency < 0.0 ||
      !std::isfinite(entry.spectralIndex) || entry.spectralIndex <= 1.0) {
    return Status::Error(StatusCode::InvalidArgument,
        "species source requires non-negative abundance and efficiency and "
        "a finite spectral index greater than one");
  }
  if (entry.energyConvention == EnergyConvention::PerNucleon &&
      !(entry.species.nucleonCount > 0.0)) {
    return Status::Error(StatusCode::UnsupportedConfiguration,
        "energy per nucleon requires a positive species nucleon count");
  }
  return Status::Ok();
}

}  // namespace

const char* EnergyConventionName(EnergyConvention convention) {
  return convention == EnergyConvention::TotalKinetic
      ? "total-kinetic-energy" : "kinetic-energy-per-nucleon";
}

Status ValidateAndNormalize(const std::vector<Configuration>& requested,
                            std::vector<Configuration>* normalized) {
  if (normalized == NULL) {
    return Status::Error(StatusCode::InvalidArgument,
                         "normalized species-source output is null");
  }
  if (requested.empty()) {
    return Status::Error(StatusCode::InvalidArgument,
                         "an explicit species source table cannot be empty");
  }

  std::set<int> identifiers;
  long double abundanceSum = 0.0L;
  for (std::size_t i = 0; i < requested.size(); ++i) {
    const Status status = ValidateEntry(requested[i]);
    if (!status.ok()) return status;
    if (!identifiers.insert(requested[i].species.modelSpecies).second) {
      return Status::Error(StatusCode::InvalidArgument,
                           "species source table contains a duplicate ID");
    }
    abundanceSum +=
        static_cast<long double>(requested[i].abundanceFraction);
  }
  if (!(abundanceSum > 0.0L) || !std::isfinite(abundanceSum)) {
    return Status::Error(StatusCode::InvalidArgument,
                         "species abundance sum must be finite and positive");
  }

  *normalized = requested;
  double accumulated = 0.0;
  for (std::size_t i = 0; i < normalized->size(); ++i) {
    Configuration& entry = (*normalized)[i];
    if (i + 1 == normalized->size()) {
      // Assign the residual to the last species so downstream number-source
      // accounting closes in the same double-precision summation order used
      // here.  Validation above guarantees a positive mathematical total.
      entry.abundanceFraction = 1.0 - accumulated;
    }
    else {
      entry.abundanceFraction = static_cast<double>(
          static_cast<long double>(entry.abundanceFraction) / abundanceSum);
      accumulated += entry.abundanceFraction;
    }
    if (entry.abundanceFraction < 0.0 ||
        !std::isfinite(entry.abundanceFraction)) {
      return Status::Error(StatusCode::InvalidArgument,
                           "normalized species abundance is invalid");
    }
  }
  return Status::Ok();
}

Status SetActiveConfiguration(
    const std::vector<Configuration>& requested) {
  std::vector<Configuration> normalized;
  const Status status = ValidateAndNormalize(requested, &normalized);
  if (!status.ok()) return status;
  gActiveConfiguration.swap(normalized);
  return Status::Ok();
}

bool HasActiveConfiguration() { return !gActiveConfiguration.empty(); }

Status FindActiveConfiguration(int modelSpecies, Configuration* result) {
  if (result == NULL) {
    return Status::Error(StatusCode::InvalidArgument,
                         "species-source result is null");
  }
  for (std::size_t i = 0; i < gActiveConfiguration.size(); ++i) {
    if (gActiveConfiguration[i].species.modelSpecies == modelSpecies) {
      *result = gActiveConfiguration[i];
      return Status::Ok();
    }
  }
  return Status::Error(StatusCode::UnsupportedConfiguration,
                       "injected species has no configured source entry");
}

ScalarResult TotalKineticEnergyJ(const Configuration& configuration,
                                 double configuredEnergyJ) {
  ScalarResult result;
  const Status status = ValidateEntry(configuration);
  if (!status.ok() || !std::isfinite(configuredEnergyJ) ||
      configuredEnergyJ < 0.0) {
    result.status = status.ok()
        ? Status::Error(StatusCode::InvalidArgument,
                        "configured source energy must be finite and non-negative")
        : status;
    return result;
  }
  const double multiplier =
      configuration.energyConvention == EnergyConvention::PerNucleon
          ? configuration.species.nucleonCount : 1.0;
  result.value = configuredEnergyJ * multiplier;
  result.status = std::isfinite(result.value)
      ? Status::Ok()
      : Status::Error(StatusCode::InvalidCoefficient,
                      "species source total energy overflowed");
  return result;
}

}  // namespace SpeciesSource
}  // namespace Transport
}  // namespace SEP
