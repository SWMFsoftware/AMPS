#ifndef SRCMOON_MOON_SURFACE_TEMPERATURE_H
#define SRCMOON_MOON_SURFACE_TEMPERATURE_H

#include <map>
#include <string>
#include <vector>

class cInternalSphericalData;

namespace Exosphere {
namespace SurfaceInteraction {
enum class AccommodationSpeedLaw;
}
}

namespace Moon { namespace Runtime { struct Configuration; } }

namespace Moon {
namespace SurfaceTemperature {

enum class Mode { SimpleCosine, DivinerMap, ThermalModel };

// Diagnostics returned by the observation-free periodic initialization used
// by I17.  Values describe the final pair of cycle boundaries.  They are
// written to the linked receipt so convergence cannot be inferred merely from
// successful process completion.
struct SpinUpMetrics {
  long long completedCycles = 0;
  double maximumNodeTemperatureChangeK = 0.0;
  double maximumColumnEnthalpyChangeJm2 = 0.0;
  double cycleStartEphemerisTimeS = 0.0;
  double cycleEndEphemerisTimeS = 0.0;
};

struct Sample {
  double latitudeDeg = 0.0;
  double longitudeEastDeg = 0.0;
  double localTimeHour = 0.0;
  double temperatureK = 0.0;
  double uncertaintyK = 0.0;
  std::string qualityFlag;
  std::string productId;
  std::string sourceRecord;
};

// Strict loader/interpolator for the processed D03 product.  The processed
// CSV is evidence, not configuration: it must retain source record/product
// IDs and quality flags.  Interpolation is permitted only inside a complete
// eight-corner cell of GOOD samples; gaps and rejected records stay gaps.
class DivinerTable {
 public:
  bool LoadCsv(const std::string& path, std::string* error);
  bool Interpolate(double latitudeDeg, double longitudeEastDeg,
                   double localTimeHour, double* temperatureK,
                   double* uncertaintyK, std::string* error) const;
  bool empty() const { return samples_.empty(); }
  std::size_t size() const { return samples_.size(); }

 private:
  struct Key {
    std::size_t latitude = 0;
    std::size_t longitude = 0;
    std::size_t localTime = 0;
    bool operator<(const Key& other) const;
  };

  std::vector<double> latitudesDeg_;
  std::vector<double> longitudesDeg_;
  std::vector<double> localTimesHour_;
  std::map<Key, Sample> samples_;
};

const char* ModeName(Mode mode);

// Production adapter shared by impacts, thermal desorption, and the epoch
// update.  Initialization occurs only after the logical AMPS surface has been
// allocated.  Thereafter particle threads only read the state; Advance() is
// called by the synchronized application time-step hook before particles move.
bool InitializeProductionState(const Runtime::Configuration& configuration,
                               cInternalSphericalData* logicalSurface,
                               std::string* error);
bool ProductionStateInitialized();
double ProductionTemperatureK(const double* positionLsoM);
bool AdvanceProductionState(double ephemerisTimeBeginS,
                            double ephemerisTimeEndS,
                            std::string* error);
// Largest absolute surface-temperature change [K] in any one internal
// thermal substep of the most recent AdvanceProductionState call.  This is an
// update-continuity diagnostic, not an observational residual.
double LastAdvanceMaximumSurfaceTemperatureChangeK();

// Exact planar-circle overlap for the approved apparent-disc definition.  The
// three inputs are apparent angular radii/separation [rad].  The Sun is
// uniformly bright and Earth is an opaque solid disc.  Atmospheric light and
// solar limb darkening are therefore absent by construction, not hidden in a
// fitted factor.  Invalid geometry returns NaN.
double UniformSolarDiscVisibleFraction(double solarAngularRadiusRad,
                                       double earthAngularRadiusRad,
                                       double centreSeparationRad);

// Replay the configured observation-free synodic forcing interval until the
// entire production thermal state is periodic at both declared tolerances.
// The routine has no arbitrary cycle limit.  Allocation/resource exhaustion
// is handled by the application's fatal diagnostics and the linked runner's
// process-resource monitor.
bool SpinUpProductionState(SpinUpMetrics* metrics,std::string* error);

// Runtime hooks installed into the normal Exosphere inventory/source path.
// The rate is atoms s^-1 for one logical surface element.  Unsupported ions
// fail at impact rather than inheriting a neutral residence law silently.
double ProductionThermalDesorptionRate(int species, int surfaceElement,
                                       void* sphereDataPointer);
bool ProductionStickingProbability(int species, double temperatureK,
                                   double* probability, std::string* error);
bool ProductionAccommodationCoefficient(int species, double temperatureK,
                                        double* coefficient,
                                        Exosphere::SurfaceInteraction::
                                            AccommodationSpeedLaw* speedLaw,
                                        std::string* error);
bool UseProductionFluxBalancedInventory();
bool WriteProductionInventoryCheckpoint(double ephemerisTimeS,
                                        std::string* error);

}  // namespace SurfaceTemperature
}  // namespace Moon

#endif  // SRCMOON_MOON_SURFACE_TEMPERATURE_H
