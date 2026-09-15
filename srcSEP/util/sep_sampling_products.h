#ifndef SEP_UTIL_SEP_SAMPLING_PRODUCTS_H
#define SEP_UTIL_SEP_SAMPLING_PRODUCTS_H

#include "sep_transport_common.h"

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace SEP {
namespace SamplingCore {

enum class InvalidParticlePolicy { Exclude, Fail };
enum class BinClassification { InRange, Underflow, Overflow, Invalid };
enum class Product {
  NumberDensityPerEnergy,
  OmnidirectionalDifferentialIntensity,
  DirectionalCrossingFlux,
  NormalizedShape
};

struct InvalidCounters {
  std::uint64_t accepted = 0;
  std::uint64_t invalidKinematics = 0;
  std::uint64_t invalidSpecies = 0;
  std::uint64_t missingBackground = 0;
  std::uint64_t underflow = 0;
  std::uint64_t overflow = 0;
};

struct Kinematics {
  Transport::Status status;
  double speedMPerS = 0.0;
  double momentumKgMPerS = 0.0;
  double mu = 0.0;
  double kineticEnergyJ = 0.0;
  double larmorRadiusM = 0.0;
};

// Invalid states are never repaired by diagnostics.  The caller can exclude
// them with counters or promote the same status to a fatal run error.
Kinematics EvaluateKinematics(double vParallelMPerS, double vPerpendicularMPerS,
                              double massKg, double chargeC,
                              double magneticFieldT, double speedOfLightMPerS);

struct BinGrid {
  std::vector<double> edges;
};

struct BinResult {
  BinClassification classification = BinClassification::Invalid;
  std::size_t index = 0;
};

BinResult Classify(const BinGrid& grid, double value);

struct WeightedBins {
  BinGrid grid;
  std::vector<double> sumWeight;
  std::vector<double> sumWeightSquared;
  InvalidCounters counters;
};

Transport::Status Initialize(const BinGrid& grid, WeightedBins* bins);
Transport::Status Accumulate(WeightedBins* bins, double value, double weight);
double EffectiveSampleSize(const WeightedBins& bins, std::size_t index);
double StandardErrorOfWeightedCount(const WeightedBins& bins,
                                    std::size_t index);

struct NormalizationContext {
  double volumeM3 = 0.0;
  double areaM2 = 0.0;
  double durationS = 0.0;
  double representativeSpeedMPerS = 0.0;
};

struct ProductValues {
  Transport::Status status;
  std::vector<double> value;
  std::vector<double> standardError;
  std::string units;
};

// Convert raw weighted counts without changing the accumulator.  Absolute
// products require the corresponding volume/area/time normalization; the
// normalized shape is explicitly dimensionless probability per coordinate.
ProductValues BuildProduct(const WeightedBins& bins, Product product,
                           const NormalizationContext& context);

struct ProductMetadata {
  std::string schemaVersion = "srcsep-sampling-v2";
  Product product = Product::NumberDensityPerEnergy;
  std::string units;
  std::string geometryFingerprint;
  std::string backgroundFingerprint;
  std::string coefficientFingerprint;
  std::string runFingerprint;
  std::string species;
  double snapshotEpochS = 0.0;
};

Transport::Status ValidateMetadata(const ProductMetadata& metadata);
std::string ProductName(Product product);

}  // namespace SamplingCore
}  // namespace SEP

#endif  // SEP_UTIL_SEP_SAMPLING_PRODUCTS_H
