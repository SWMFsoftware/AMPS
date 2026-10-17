#ifndef SRCMOON_MOON_SURFACE_TEMPERATURE_H
#define SRCMOON_MOON_SURFACE_TEMPERATURE_H

#include <array>
#include <cstddef>
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

// Complete read-only state of one production logical surface element at one
// SPICE epoch.  Cartesian vectors are in MOON_ME_DE421; distances are metres.
// Terrain visibility and finite-disc Earth occultation remain separate so a
// downstream plot cannot mistake an eclipse, the geometric nightside, and a
// topographic horizon for the same physical mechanism.
struct SurfaceEnvironmentSample {
  std::size_t surfaceElement = 0;
  int zenithIndex = 0;
  int azimuthIndex = 0;
  double latitudeDeg = 0.0;
  double longitudeEastDeg = 0.0;
  double localSolarTimeHour = 0.0;
  std::array<double,3> surfacePositionBodyFixedM{{0.0,0.0,0.0}};
  std::array<double,3> terrainNormalBodyFixed{{0.0,0.0,0.0}};
  double solarIncidenceCosine = 0.0;
  bool sunAboveLocalHorizon = false;
  bool terrainSunVisible = false;
  bool terrainOccluded = false;
  double solarDiscVisibleFraction = 0.0;
  double directBeamFraction = 0.0;
  double heliocentricDistanceAu = 0.0;
  double absorbedSolarFluxWm2 = 0.0;
  double surfaceTemperatureK = 0.0;
  double lastEnergyResidualJm2 = 0.0;
  int lastNonlinearIterations = 0;
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

// Sample every logical production surface element without changing its
// thermal state.  The returned illumination is evaluated at ephemerisTimeS;
// temperatures are the state reached by the caller's preceding initialization
// and AdvanceProductionState()/SpinUpProductionState() sequence.  This is a
// model-output service, not an observational comparison or I15/I17 result.
bool SampleProductionSurfaceEnvironment(
    double ephemerisTimeS, std::vector<SurfaceEnvironmentSample>* samples,
    std::string* error);

// Return all finite-volume layer-centre depths [m] and temperatures [K] for
// one production logical surface element.  The surface boundary temperature
// itself is carried in SurfaceEnvironmentSample and is not duplicated as a
// zero-depth finite-volume cell.
bool ProductionThermalProfile(std::size_t surfaceElement,
                              std::vector<double>* layerDepthM,
                              std::vector<double>* layerTemperatureK,
                              std::string* error);

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
