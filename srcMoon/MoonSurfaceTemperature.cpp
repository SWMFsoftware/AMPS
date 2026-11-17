#include "pic.h"

#include "MoonSurfaceTemperature.h"

#include "Moon.h"
#include "MoonInput.h"
#include "MoonSurfacePhysics.h"
#include "MoonSurfaceEnvironment.h"
#include "MoonThermalModel.h"
#include "MoonDivinerObservation.h"

#include <algorithm>
#include <array>
#include <cerrno>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <fcntl.h>
#include <limits>
#include <set>
#include <sstream>
#include <tuple>
#include <unistd.h>

#ifndef _NO_SPICE_CALLS_
#include "SpiceUsr.h"
#endif

namespace Moon {
namespace SurfaceTemperature {
namespace {

struct ProductionState {
  bool initialized = false;
  const Runtime::Configuration* configuration = nullptr;
  cInternalSphericalData* logicalSurface = nullptr;
  Mode mode = Mode::SimpleCosine;
  DivinerTable diviner;
  std::vector<Thermal::Column> columns;
  std::vector<double> temperaturesK;
  std::vector<unsigned char> geometricShadow;
  double lastAdvanceMaximumSurfaceTemperatureChangeK = 0.0;
  SurfacePhysics::SpeciesParameters sodium;
  SurfacePhysics::SpeciesParameters argon;
  double lastCheckpointEtS = -std::numeric_limits<double>::infinity();
  double photonFluxEpochS=std::numeric_limits<double>::quiet_NaN();
  std::vector<double> photonFluxFactors;
  double prescribedTemperatureEpochS=std::numeric_limits<double>::quiet_NaN();
  ProductionThermalInitializationInfo thermalInitialization;
  bool thermalInitializationAvailable=false;
  bool initialThermalStateStillApplied=false;
  double currentThermalEpochS=std::numeric_limits<double>::quiet_NaN();
  // A completed spin-up can be reused only while these exact columns have
  // not evolved. An epoch comparison alone is insufficient: replaying an
  // interval can end at the same timestamp with different physical columns.
  bool completedSpinupAtCurrentState=false;
  SpinUpMetrics completedSpinup;
};

ProductionState gProduction;

constexpr double kElectronVoltJ = 1.602176634e-19;
constexpr double kPi = 3.141592653589793238462643383279502884;

bool PositiveFinite(double value) {
  return std::isfinite(value) && value > 0.0;
}

// Small fixed-size vector helpers are kept local so the production thermal
// driver uses one explicit Cartesian convention (metres for surface points,
// kilometres for raw SPICE states) without depending on an AMPS vector type
// whose units are implicit at the call site.
double Dot(const std::array<double,3>& first,
           const std::array<double,3>& second) {
  return first[0]*second[0]+first[1]*second[1]+first[2]*second[2];
}

double Norm(const std::array<double,3>& value) {
  return std::sqrt(Dot(value,value));
}

std::array<double,3> TransformPosition(const double transform[6][6],
                                       const double* position) {
  std::array<double,3> result{{0.0,0.0,0.0}};
  for (int row=0;row<3;++row) for (int column=0;column<3;++column)
    result[row]+=transform[row][column]*position[column];
  return result;
}

double NormalizeDegrees(double value) {
  value=std::fmod(value,360.0);
  return value<0.0 ? value+360.0 : value;
}

double LocalSolarTimeHour(const std::array<double,3>& positionBody,
                          const std::array<double,3>& sunBody) {
  // Local apparent solar time is 12 h on the body-fixed meridian containing
  // the Sun.  Positive-east longitude therefore advances local time by
  // 24 h/360 deg relative to the subsolar longitude.
  const double longitude=std::atan2(positionBody[1],positionBody[0])*180.0/kPi;
  const double subsolarLongitude=
      std::atan2(sunBody[1],sunBody[0])*180.0/kPi;
  return NormalizeDegrees(longitude-subsolarLongitude+180.0)/15.0;
}

#ifndef _NO_SPICE_CALLS_
// Geometry shared by every surface element at one TDB epoch.  SPICE states
// are retained in their native kilometres until translated to an actual
// terrain point below; the unit Sun direction is body-fixed and dimensionless.
struct EpochIlluminationGeometry {
  std::array<double,3> sunStateKm{{0.0,0.0,0.0}};
  std::array<double,3> earthStateKm{{0.0,0.0,0.0}};
  std::array<double,3> unitSunBodyFixed{{0.0,0.0,0.0}};
  double heliocentricDistanceAu=0.0;
  double ephemerisTimeS=0.0;
};

// Element-local quantities used both by the thermal time integrator and by
// the read-only surface-environment exporter.  Keeping one implementation is
// essential: a visually plausible demonstration must not use a shadow law
// that differs from the production thermal boundary condition.
struct ElementIlluminationGeometry {
  std::array<double,3> surfacePositionM{{0.0,0.0,0.0}};
  std::array<double,3> terrainNormal{{0.0,0.0,0.0}};
  double solarIncidenceCosine=0.0;
  bool sunAboveLocalHorizon=false;
  bool terrainSunVisible=false;
  bool terrainOccluded=false;
  double solarDiscVisibleFraction=0.0;
  double directBeamFraction=0.0;
  double absorbedSolarFluxWm2=0.0;
  double heliocentricDistanceAu=0.0;
  double pointSunDistanceM=0.0;
  double projectedSolarIncidence=0.0;
};

bool EvaluateEpochIlluminationGeometry(double ephemerisTimeS,
    EpochIlluminationGeometry* geometry,std::string* error) {
  if (geometry==nullptr || error==nullptr || !std::isfinite(ephemerisTimeS))
    return false;
  SpiceDouble sunState[6],earthState[6],lightTime;
  spkezr_c("SUN",ephemerisTimeS,Frames::BodyFixed,
      gProduction.configuration->illuminationAberrationCorrection.c_str(),ObjectName,sunState,&lightTime);
  spkezr_c("EARTH",ephemerisTimeS,Frames::BodyFixed,
      gProduction.configuration->illuminationAberrationCorrection.c_str(),ObjectName,earthState,&lightTime);
  const double sunLength=std::sqrt(sunState[0]*sunState[0]+
      sunState[1]*sunState[1]+sunState[2]*sunState[2]);
  if (!(sunLength>0.0) || !std::isfinite(sunLength)) {
    *error="SPICE returned a zero or non-finite Moon-to-Sun vector for "
        "thermal forcing";
    return false;
  }
  for (int dimension=0;dimension<3;++dimension) {
    geometry->sunStateKm[dimension]=sunState[dimension];
    geometry->earthStateKm[dimension]=earthState[dimension];
    geometry->unitSunBodyFixed[dimension]=sunState[dimension]/sunLength;
  }
  geometry->heliocentricDistanceAu=sunLength/1.495978707e8;
  geometry->ephemerisTimeS=ephemerisTimeS;
  return true;
}

bool EvaluateElementIlluminationGeometry(std::size_t element,
    const EpochIlluminationGeometry& epoch,
    ElementIlluminationGeometry* geometry,std::string* error,
    bool includeAbsorbedFlux=true) {
  if (geometry==nullptr || error==nullptr || gProduction.logicalSurface==nullptr ||
      element>=gProduction.temperaturesK.size()) return false;
  const Runtime::Configuration& configuration=*gProduction.configuration;
  geometry->heliocentricDistanceAu=epoch.heliocentricDistanceAu;
  double radialNormal[3];
  gProduction.logicalSurface->GetSurfaceElementNormal(
      radialNormal,static_cast<int>(element));
  geometry->surfacePositionM={{_RADIUS_(_MOON_)*radialNormal[0],
      _RADIUS_(_MOON_)*radialNormal[1],
      _RADIUS_(_MOON_)*radialNormal[2]}};
  geometry->terrainNormal={{radialNormal[0],radialNormal[1],radialNormal[2]}};
  std::size_t terrainFace=std::numeric_limits<std::size_t>::max();
  const auto* terrain=SurfaceEnvironment::ProductionTerrain();
  if (terrain!=nullptr) {
    SurfaceEnvironment::RayHit radialHit;
    const std::array<double,3> centre{{0.0,0.0,0.0}};
    if (!terrain->IntersectRay(centre,geometry->terrainNormal,0.0,
            std::numeric_limits<double>::infinity(),terrainFace,
            &radialHit,error)) return false;
    if (!radialHit.hit) {
      *error="thermal column radial direction misses production terrain";
      return false;
    }
    geometry->surfacePositionM=radialHit.positionM;
    geometry->terrainNormal=radialHit.outwardNormal;
    terrainFace=radialHit.face;
  }

  // Earth eclipse is independent of terrain blocking.  Translate Moon-centred
  // SPICE states [km] to the actual terrain point [m] before calculating
  // apparent angular radii, preserving the production parallax correction.
  std::array<double,3> pointToSun,pointToEarth;
  for (int dimension=0;dimension<3;++dimension) {
    pointToSun[dimension]=epoch.sunStateKm[dimension]-
        geometry->surfacePositionM[dimension]/1000.0;
    pointToEarth[dimension]=epoch.earthStateKm[dimension]-
        geometry->surfacePositionM[dimension]/1000.0;
  }
  if (configuration.illuminationAberrationCorrection!="NONE") {
    // Apparent reception corrections must use the actual rotating surface
    // observer. Correcting a Moon-centred state and then subtracting a site
    // is not a topocentric light-time/stellar-aberration calculation.
    std::array<double,3> sunM,earthM;
    if (!SurfaceEnvironment::ApparentSolarSystemGeometryAtPoint(epoch.ephemerisTimeS,
        geometry->surfacePositionM,configuration.illuminationAberrationCorrection,
        &sunM,&earthM,error)) return false;
    for (int k=0;k<3;++k) {pointToSun[k]=sunM[k]/1000;pointToEarth[k]=earthM[k]/1000;}
  }
  const double sunDistanceKm=Norm(pointToSun);
  const double earthDistanceKm=Norm(pointToEarth);
  if (!(sunDistanceKm>configuration.thermalSolarDiscRadiusM/1000.0) ||
      !(earthDistanceKm>configuration.thermalOpaqueEarthRadiusM/1000.0)) {
    *error="surface point lies inside a configured eclipse body radius";
    return false;
  }
  geometry->pointSunDistanceM=1000*sunDistanceKm;
  geometry->heliocentricDistanceAu=sunDistanceKm/1.495978707e8;
  std::array<double,3> unitSun=pointToSun;
  for (double& value:unitSun) value/=sunDistanceKm;
  geometry->solarIncidenceCosine=std::max(-1.0,std::min(1.0,Dot(geometry->terrainNormal,unitSun)));
  geometry->sunAboveLocalHorizon=geometry->solarIncidenceCosine>0;
  geometry->terrainSunVisible=geometry->sunAboveLocalHorizon;
  if (terrain && geometry->terrainSunVisible && !terrain->IsSunVisible(
      geometry->surfacePositionM,geometry->terrainNormal,unitSun,
      configuration.illuminationRayOriginOffsetM,terrainFace,&geometry->terrainSunVisible,error)) return false;
  geometry->terrainOccluded=geometry->sunAboveLocalHorizon && !geometry->terrainSunVisible;
  geometry->solarDiscVisibleFraction=1.0;
  if (earthDistanceKm<sunDistanceKm) {
    const double directionCosine=std::max(-1.0,std::min(1.0,
        Dot(pointToSun,pointToEarth)/(sunDistanceKm*earthDistanceKm)));
    const double separation=std::acos(directionCosine);
    const double solarRadius=std::asin(
        configuration.thermalSolarDiscRadiusM/(1000.0*sunDistanceKm));
    const double earthRadius=std::asin(
        configuration.thermalOpaqueEarthRadiusM/(1000.0*earthDistanceKm));
    geometry->solarDiscVisibleFraction=UniformSolarDiscVisibleFraction(
        solarRadius,earthRadius,separation);
    if (!std::isfinite(geometry->solarDiscVisibleFraction)) {
      *error="finite-disc Earth-eclipse geometry is invalid";
      return false;
    }
  }
  geometry->directBeamFraction=geometry->terrainSunVisible ?
      geometry->solarDiscVisibleFraction:0.0;
  geometry->projectedSolarIncidence=std::max(0.0,geometry->solarIncidenceCosine)*geometry->directBeamFraction;
  if (configuration.illuminationTerrainSolarModel=="uniform-disc") {
    SurfaceEnvironment::SolarDiscRequest request;
    request.surfacePositionM=geometry->surfacePositionM;
    request.outwardNormal=geometry->terrainNormal;
    for (int k=0;k<3;++k) {
      request.pointToSunM[k]=1000*pointToSun[k];
      request.pointToEarthM[k]=1000*pointToEarth[k];
    }
    request.solarRadiusM=configuration.thermalSolarDiscRadiusM;
    request.opaqueEarthRadiusM=configuration.thermalOpaqueEarthRadiusM;
    request.rayOriginOffsetM=configuration.illuminationRayOriginOffsetM;
    request.sourceFace=terrainFace;
    request.radialPoints=configuration.illuminationDiscRadialPoints;
    request.azimuthPoints=configuration.illuminationDiscAzimuthPoints;
    SurfaceEnvironment::SolarDiscSample disc;
    std::function<double(double)> absorptivity;
    if (includeAbsorbedFlux) {
      if (element>=gProduction.columns.size()) {
        *error="absorbed forcing requires the initialized production thermal column";return false;
      }
      const auto* column=&gProduction.columns[element];
      absorptivity=[column](double cosine) {return 1-column->DirectionalBondAlbedo(std::acos(cosine));};
    }
    if (!SurfaceEnvironment::EvaluateSolarDisc(terrain,request,&disc,error,absorptivity)) return false;
    // Joint masks are integrated ray by ray, including directional albedo;
    // incidence/area/occurrence remain distinct diagnostics. Use the actual
    // point-to-Sun distance for the inverse-square irradiance in this mode.
    geometry->terrainSunVisible=disc.sampledAnyTerrainVisible;
    geometry->terrainOccluded=geometry->sunAboveLocalHorizon &&
        !disc.sampledAnyTerrainVisible;
    geometry->solarDiscVisibleFraction=disc.earthVisibleAreaFraction;
    geometry->directBeamFraction=disc.jointlyVisibleAreaFraction;
    geometry->projectedSolarIncidence=disc.projectedIncidence;
    const double distanceAu=sunDistanceKm/1.495978707e8;
    geometry->heliocentricDistanceAu=distanceAu;
    geometry->solarIncidenceCosine=std::max(-1.0,std::min(1.0,
        Dot(geometry->terrainNormal,pointToSun)/sunDistanceKm));
    geometry->sunAboveLocalHorizon=geometry->solarIncidenceCosine>0;
    geometry->terrainOccluded=geometry->sunAboveLocalHorizon &&
        disc.terrainVisibleAreaFraction<1;
    geometry->absorbedSolarFluxWm2=includeAbsorbedFlux ? configuration.thermal.solarConstantWm2/
        (distanceAu*distanceAu)*disc.absorbedProjectedIncidence:0;
  } else if (includeAbsorbedFlux) {
    geometry->absorbedSolarFluxWm2=gProduction.columns[element].AbsorbedSolarFluxWm2(
        geometry->solarIncidenceCosine,geometry->directBeamFraction,
        geometry->heliocentricDistanceAu);
  }
  if (!std::isfinite(geometry->absorbedSolarFluxWm2)) {
    *error="non-finite production absorbed solar flux";
    return false;
  }
  return true;
}

// One prescribed boundary function serves actual impacts and representative
// logical-cell source temperatures. Coordinates are body-fixed metres and
// the supplied Sun/Earth states are at one explicit TDB epoch; no global
// end/midpoint orbit state is borrowed during a split-step update.
bool PrescribedTemperatureAtBodyPoint(const std::array<double,3>& body,
    const EpochIlluminationGeometry& epoch,double* kineticK,std::string* error,
    bool* mapCoverageUnavailable=nullptr) {
  if (kineticK==nullptr || error==nullptr) return false;
  if (mapCoverageUnavailable) *mapCoverageUnavailable=false;
  const double radius=Norm(body);
  if (!(radius>0) || !std::isfinite(radius)) {
    *error="prescribed temperature requires a finite body-fixed surface point";return false;
  }
  const auto& configuration=*gProduction.configuration;
  double temperature=0;
  if (gProduction.mode==Mode::SimpleCosine) {
    // This explicitly retained debug boundary keeps its original radial
    // cosine and cylindrical Earth shadow. It is not a predictive eclipse
    // or terrain heat equation. Source cells and impacts now use it alike.
    const double cosine=Dot(body,epoch.unitSunBodyFixed)/radius;
    std::array<double,3> sunM,earthM;
    for (int k=0;k<3;++k) {
      sunM[k]=1000*epoch.sunStateKm[k];earthM[k]=1000*epoch.earthStateKm[k];
    }
    const bool night=cosine<0 || Moon::EarthShadowCheck(body.data(),sunM.data(),earthM.data());
    temperature=night ? configuration.simpleCosineNightTemperatureK:
        configuration.simpleCosineNightTemperatureK+
        configuration.simpleCosineDayAmplitudeK*std::pow(std::min(1.0,cosine),configuration.simpleCosineExponent);
  } else if (gProduction.mode==Mode::DivinerMap) {
    if (configuration.divinerMapKineticInterpretation!="isothermal-grey") {
      *error="Diviner brightness cannot supply kinetic temperature without explicit isothermal-grey interpretation";
      return false;
    }
    const double latitude=std::asin(std::max(-1.0,std::min(1.0,body[2]/radius)))*180/kPi;
    const double longitude=NormalizeDegrees(std::atan2(body[1],body[0])*180/kPi);
    double brightnessK=0,uncertaintyK=0;
    if (!gProduction.diviner.Interpolate(latitude,longitude,
            LocalSolarTimeHour(body,epoch.unitSunBodyFixed),&brightnessK,&uncertaintyK,error)) {
      // The table was structurally validated on load. A missing interpolation
      // support is real observational coverage, not an invented climate value.
      if (mapCoverageUnavailable) *mapCoverageUnavailable=true;
      return false;
    }
    // Approved isothermal grey-facet MODEL interpretation. The archived
    // brightness and its uncertainty are never changed or relabelled as
    // measured kinetic temperature; anisothermality remains a scope gate.
    temperature=DivinerObservation::GreyKineticTemperatureK(
        brightnessK,configuration.thermal.infraredEmissivity);
  } else {
    *error="prescribed temperature helper cannot replace a thermal-model state";return false;
  }
  if (!std::isfinite(temperature) || temperature<=0) {
    *error="prescribed boundary produced nonpositive/nonfinite kinetic temperature";return false;
  }
  *kineticK=temperature;
  return true;
}

bool RefreshPrescribedTemperatureCache(double epochS,std::string* error) {
  if (gProduction.mode==Mode::ThermalModel) return true;
  if (gProduction.prescribedTemperatureEpochS==epochS) return true;
  EpochIlluminationGeometry epoch;
  if (!EvaluateEpochIlluminationGeometry(epochS,&epoch,error)) return false;
  std::vector<double> temperatures(gProduction.temperaturesK.size(),
      std::numeric_limits<double>::quiet_NaN());
  for (std::size_t element=0;element<temperatures.size();++element) {
    ElementIlluminationGeometry geometry;
    if (!EvaluateElementIlluminationGeometry(element,epoch,&geometry,error,false)) return false;
    bool missingMapSupport=false;
    if (!PrescribedTemperatureAtBodyPoint(geometry.surfacePositionM,epoch,
            &temperatures[element],error,&missingMapSupport)) {
      if (!missingMapSupport) return false;
      // Empty inventory cells need not have a measured map temperature.
      // Leave the NaN coverage marker: any later populated source or kinetic
      // query fails explicitly. In particular, do not copy the initialization
      // temperature into real archive gaps just to make I16 initialize.
      error->clear();
    }
  }
  // Publish a complete epoch atomically with respect to subsequent serial
  // source dispatch. Particle callbacks only read this immutable vector.
  gProduction.temperaturesK=std::move(temperatures);
  gProduction.prescribedTemperatureEpochS=epochS;
  return true;
}
#endif

const SurfacePhysics::SpeciesParameters* ParametersForSpecies(int species) {
#if _NA_SPEC_ >= 0
  if (species==_NA_SPEC_) return &gProduction.sodium;
#endif
#if _AR_SPEC_ >= 0
  if (species==_AR_SPEC_) return &gProduction.argon;
#endif
  return nullptr;
}

std::string CheckpointPath(const std::string& base,const char* species) {
  return base+"."+species;
}

bool ParseFinite(const std::string& text, double* value) {
  errno = 0;
  char* end = nullptr;
  const double parsed = std::strtod(text.c_str(), &end);
  if (errno == ERANGE || end == text.c_str() || *end != '\0' ||
      !std::isfinite(parsed)) return false;
  *value = parsed;
  return true;
}

// Processed products deliberately use a restricted CSV dialect: fields may
// not contain commas or embedded quoting.  Rejecting such records prevents a
// malformed provenance string from shifting scientific columns silently.
std::vector<std::string> SplitCsv(const std::string& line) {
  std::vector<std::string> fields;
  std::stringstream stream(line);
  std::string field;
  while (std::getline(stream, field, ',')) fields.push_back(field);
  if (!line.empty() && line.back() == ',') fields.emplace_back();
  return fields;
}

double NormalizeLongitude(double longitudeDeg) {
  longitudeDeg = std::fmod(longitudeDeg, 360.0);
  if (longitudeDeg < 0.0) longitudeDeg += 360.0;
  return longitudeDeg;
}

double NormalizeLocalTime(double hour) {
  hour = std::fmod(hour, 24.0);
  if (hour < 0.0) hour += 24.0;
  return hour;
}

bool BoundingLinear(const std::vector<double>& axis, double value,
                    std::size_t* lower, std::size_t* upper, double* fraction) {
  if (axis.empty() || value < axis.front() || value > axis.back()) return false;
  const auto high = std::lower_bound(axis.begin(), axis.end(), value);
  if (high != axis.end() && *high == value) {
    *lower = *upper = static_cast<std::size_t>(high - axis.begin());
    *fraction = 0.0;
    return true;
  }
  if (high == axis.begin() || high == axis.end()) return false;
  *upper = static_cast<std::size_t>(high - axis.begin());
  *lower = *upper - 1;
  *fraction = (value - axis[*lower]) / (axis[*upper] - axis[*lower]);
  return true;
}

// Periodic axes use a synthetic wrap interval from the final node to the first
// node plus one period.  This preserves continuity at 0/360 longitude and
// 0/24 local time without duplicating rows in the processed evidence file.
bool BoundingPeriodic(const std::vector<double>& axis, double value,
                      double period, std::size_t* lower, std::size_t* upper,
                      double* fraction) {
  if (axis.empty()) return false;
  const auto high = std::lower_bound(axis.begin(), axis.end(), value);
  if (high != axis.end() && *high == value) {
    *lower = *upper = static_cast<std::size_t>(high - axis.begin());
    *fraction = 0.0;
    return true;
  }
  if (high == axis.end()) {
    *lower = axis.size() - 1;
    *upper = 0;
    *fraction = (value - axis.back()) /
        (axis.front() + period - axis.back());
  } else if (high == axis.begin()) {
    *lower = axis.size() - 1;
    *upper = 0;
    *fraction = (value + period - axis.back()) /
        (axis.front() + period - axis.back());
  } else {
    *upper = static_cast<std::size_t>(high - axis.begin());
    *lower = *upper - 1;
    *fraction = (value - axis[*lower]) / (axis[*upper] - axis[*lower]);
  }
  return std::isfinite(*fraction) && *fraction >= 0.0 && *fraction <= 1.0;
}

}  // namespace

bool DivinerTable::Key::operator<(const Key& other) const {
  return std::tie(latitude, longitude, localTime) <
      std::tie(other.latitude, other.longitude, other.localTime);
}

bool DivinerTable::LoadCsv(const std::string& path, std::string* error) {
  if (error == nullptr) return false;
  error->clear();
  std::ifstream input(path);
  if (!input) {
    *error = "cannot open processed Diviner table '" + path + "'";
    return false;
  }
  std::string line;
  if (!std::getline(input, line)) {
    *error = "processed Diviner table is empty";
    return false;
  }
  const std::string requiredHeader =
      "latitude_deg,longitude_east_deg,local_time_hour,temperature_k,"
      "uncertainty_k,quality_flag,product_id,source_record";
  if (line != requiredHeader) {
    *error = "processed Diviner table header does not match the frozen D03 "
             "schema";
    return false;
  }

  std::vector<Sample> rows;
  std::set<double> latitudes, longitudes, localTimes;
  std::size_t lineNumber = 1;
  while (std::getline(input, line)) {
    ++lineNumber;
    if (line.empty()) continue;
    const std::vector<std::string> fields = SplitCsv(line);
    Sample sample;
    // The Diviner GCP archive does not publish a per-bin uncertainty.  An
    // empty field therefore means "not supplied", not zero uncertainty.
    // Preserve that distinction as NaN while still rejecting malformed
    // non-empty values.  Callers may use the temperature, but validation code
    // must not construct an uncertainty-weighted residual from the NaN.
    sample.uncertaintyK = std::numeric_limits<double>::quiet_NaN();
    const bool uncertaintyValid = fields.size() == 8 &&
        (fields[4].empty() || ParseFinite(fields[4], &sample.uncertaintyK));
    if (fields.size() != 8 ||
        !ParseFinite(fields[0], &sample.latitudeDeg) ||
        !ParseFinite(fields[1], &sample.longitudeEastDeg) ||
        !ParseFinite(fields[2], &sample.localTimeHour) ||
        !ParseFinite(fields[3], &sample.temperatureK) ||
        !uncertaintyValid ||
        sample.latitudeDeg < -90.0 || sample.latitudeDeg > 90.0 ||
        sample.longitudeEastDeg < 0.0 || sample.longitudeEastDeg >= 360.0 ||
        sample.localTimeHour < 0.0 || sample.localTimeHour >= 24.0 ||
        sample.temperatureK <= 0.0 ||
        (std::isfinite(sample.uncertaintyK) && sample.uncertaintyK < 0.0) ||
        fields[5].empty() || fields[6].empty() || fields[7].empty()) {
      std::ostringstream message;
      message << "invalid processed Diviner record at line " << lineNumber;
      *error = message.str();
      return false;
    }
    sample.qualityFlag = fields[5];
    sample.productId = fields[6];
    sample.sourceRecord = fields[7];
    rows.push_back(sample);
    latitudes.insert(sample.latitudeDeg);
    longitudes.insert(sample.longitudeEastDeg);
    localTimes.insert(sample.localTimeHour);
  }
  if (rows.empty()) {
    *error = "processed Diviner table contains no records";
    return false;
  }

  std::vector<double> latitudeAxis(latitudes.begin(), latitudes.end());
  std::vector<double> longitudeAxis(longitudes.begin(), longitudes.end());
  std::vector<double> localTimeAxis(localTimes.begin(), localTimes.end());
  std::map<Key, Sample> indexed;
  for (const Sample& sample : rows) {
    const Key key{
        static_cast<std::size_t>(std::lower_bound(latitudeAxis.begin(),
            latitudeAxis.end(), sample.latitudeDeg) - latitudeAxis.begin()),
        static_cast<std::size_t>(std::lower_bound(longitudeAxis.begin(),
            longitudeAxis.end(), sample.longitudeEastDeg) -
            longitudeAxis.begin()),
        static_cast<std::size_t>(std::lower_bound(localTimeAxis.begin(),
            localTimeAxis.end(), sample.localTimeHour) - localTimeAxis.begin())};
    if (!indexed.emplace(key, sample).second) {
      *error = "processed Diviner table contains a duplicate grid coordinate";
      return false;
    }
  }
  latitudesDeg_ = std::move(latitudeAxis);
  longitudesDeg_ = std::move(longitudeAxis);
  localTimesHour_ = std::move(localTimeAxis);
  samples_ = std::move(indexed);
  return true;
}

bool DivinerTable::Interpolate(double latitudeDeg, double longitudeEastDeg,
                               double localTimeHour, double* temperatureK,
                               double* uncertaintyK,
                               std::string* error) const {
  if (error == nullptr || temperatureK == nullptr || uncertaintyK == nullptr)
    return false;
  error->clear();
  if (samples_.empty() || !std::isfinite(latitudeDeg) ||
      !std::isfinite(longitudeEastDeg) || !std::isfinite(localTimeHour)) {
    *error = "invalid or unloaded Diviner interpolation request";
    return false;
  }
  longitudeEastDeg = NormalizeLongitude(longitudeEastDeg);
  localTimeHour = NormalizeLocalTime(localTimeHour);
  std::size_t lat[2], lon[2], time[2];
  double latitudeFraction = 0.0, longitudeFraction = 0.0,
      timeFraction = 0.0;
  if (!BoundingLinear(latitudesDeg_, latitudeDeg, &lat[0], &lat[1],
                      &latitudeFraction) ||
      !BoundingPeriodic(longitudesDeg_, longitudeEastDeg, 360.0,
                        &lon[0], &lon[1], &longitudeFraction) ||
      !BoundingPeriodic(localTimesHour_, localTimeHour, 24.0,
                        &time[0], &time[1], &timeFraction)) {
    *error = "Diviner request lies outside the processed latitude support";
    return false;
  }

  double weightedTemperature = 0.0;
  double weightedVariance = 0.0;
  double weightSum = 0.0;
  bool allUncertaintiesAvailable = true;
  for (int i = 0; i < 2; ++i) for (int j = 0; j < 2; ++j)
    for (int k = 0; k < 2; ++k) {
      // When the request lies exactly on an axis node, the duplicate branch
      // receives zero weight and is not required.  Thus a real gap adjacent to
      // an exact sample does not incorrectly invalidate that sample.
      const double weight = (i ? latitudeFraction : 1.0-latitudeFraction) *
          (j ? longitudeFraction : 1.0-longitudeFraction) *
          (k ? timeFraction : 1.0-timeFraction);
      if (weight == 0.0) continue;
      const auto found = samples_.find(Key{lat[i], lon[j], time[k]});
      if (found == samples_.end() || found->second.qualityFlag != "GOOD") {
        *error = "Diviner interpolation cell contains missing or non-GOOD "
                 "evidence; no gap filling is permitted";
        return false;
      }
      weightedTemperature += weight * found->second.temperatureK;
      if (std::isfinite(found->second.uncertaintyK)) {
        weightedVariance += weight * weight *
            found->second.uncertaintyK * found->second.uncertaintyK;
      }
      else {
        // Never replace an archive metadata gap with a numerical zero.  One
        // missing corner is enough to make the interpolated uncertainty
        // unavailable even though the temperature interpolation remains
        // physically well defined.
        allUncertaintiesAvailable = false;
      }
      weightSum += weight;
    }
  if (!(weightSum > 0.0)) {
    *error = "Diviner interpolation produced zero total weight";
    return false;
  }
  *temperatureK = weightedTemperature / weightSum;
  *uncertaintyK = allUncertaintiesAvailable
      ? std::sqrt(weightedVariance) / weightSum
      : std::numeric_limits<double>::quiet_NaN();
  return true;
}

const char* ModeName(Mode mode) {
  switch (mode) {
    case Mode::SimpleCosine: return "simple-cosine";
    case Mode::DivinerMap: return "diviner-map";
    case Mode::ThermalModel: return "thermal-model";
  }
  return "unknown";
}

double UniformSolarDiscVisibleFraction(double solarRadius,
                                       double earthRadius,
                                       double separation) {
  if (!PositiveFinite(solarRadius) || !PositiveFinite(earthRadius) ||
      !std::isfinite(separation) || separation<0.0)
    return std::numeric_limits<double>::quiet_NaN();

  // The apparent discs are defined as circles in angular image coordinates.
  // The following standard two-circle intersection is analytical: there is no
  // ray-grid resolution, Monte-Carlo noise, or tunable transition width.  It
  // is exact for the approved uniform projected-disc model; omitted solar
  // limb darkening and Earth-atmosphere light are separate declared physics.
  if (separation>=solarRadius+earthRadius) return 1.0;
  if (earthRadius>=separation+solarRadius) return 0.0;
  if (solarRadius>=separation+earthRadius) {
    const double covered=(earthRadius*earthRadius)/
        (solarRadius*solarRadius);
    return std::max(0.0,std::min(1.0,1.0-covered));
  }

  const auto clamped=[](double value) {
    return std::max(-1.0,std::min(1.0,value));
  };
  const double d2=separation*separation;
  const double s2=solarRadius*solarRadius;
  const double e2=earthRadius*earthRadius;
  const double solarSector=s2*std::acos(clamped(
      (d2+s2-e2)/(2.0*separation*solarRadius)));
  const double earthSector=e2*std::acos(clamped(
      (d2+e2-s2)/(2.0*separation*earthRadius)));
  const double radical=std::max(0.0,
      (-separation+solarRadius+earthRadius)*
      ( separation+solarRadius-earthRadius)*
      ( separation-solarRadius+earthRadius)*
      ( separation+solarRadius+earthRadius));
  const double overlap=solarSector+earthSector-0.5*std::sqrt(radical);
  return std::max(0.0,std::min(1.0,
      1.0-overlap/(kPi*s2)));
}

bool InitializeProductionState(const Runtime::Configuration& configuration,
                               cInternalSphericalData* logicalSurface,
                               std::string* error) {
  if (error==nullptr) return false;
  error->clear();
  if (!configuration.m2EnvironmentEnabled) {
    // A caller may deliberately return to the baseline in the same process.
    // Release only OUR provider before its immutable geometry cache is
    // destroyed; a foreign application's optional core hook is not ours to
    // reset. A fresh baseline process already has the historical NULL hook.
    auto& photonProvider=Exosphere::SourceProcesses::PhotonStimulatedDesorption::RuntimePhotonFluxFactor;
    if (photonProvider==ProductionPhotonFluxFactor) photonProvider=nullptr;
    gProduction=ProductionState{};
    return true;
  }
#ifdef _NO_SPICE_CALLS_
  // All approved M2 modes require an epoch-dependent LSO-to-body transform;
  // the thermal mode additionally requires Sun/Earth ephemerides.  Refuse the
  // configuration instead of inventing a frozen illumination geometry.
  *error="M2 surface/environment state requires a SPICE-enabled build";
  return false;
#endif
  if (logicalSurface==nullptr) {
    *error="M2 surface state requires an allocated logical AMPS surface";
    return false;
  }

  ProductionState candidate;
  candidate.configuration=&configuration;
  candidate.logicalSurface=logicalSurface;
  if (configuration.surfaceTemperatureMode=="simple-cosine")
    candidate.mode=Mode::SimpleCosine;
  else if (configuration.surfaceTemperatureMode=="diviner-map")
    candidate.mode=Mode::DivinerMap;
  else if (configuration.surfaceTemperatureMode=="thermal-model")
    candidate.mode=Mode::ThermalModel;
  else {
    *error="unsupported M2 surface-temperature mode";
    return false;
  }

  if (candidate.mode==Mode::DivinerMap &&
      !candidate.diviner.LoadCsv(configuration.divinerTemperatureTableFile,
                                 error)) return false;

  const std::size_t count=static_cast<std::size_t>(
      logicalSurface->GetTotalSurfaceElementsNumber());
  if (count==0) {
    *error="M2 logical surface contains no inventory elements";
    return false;
  }
  candidate.temperaturesK.assign(count,configuration.thermalInitialTemperatureK);
  candidate.geometricShadow.assign(count,0);
  if (candidate.mode==Mode::ThermalModel) {
    candidate.columns.resize(count);
    for (std::size_t element=0;element<count;++element)
      if (!candidate.columns[element].Initialize(configuration.thermal,
              configuration.thermalInitialTemperatureK,error)) return false;
  }

  candidate.sodium.species=SurfacePhysics::Species::Sodium;
  candidate.sodium.residenceLaw=
      SurfacePhysics::ResidenceLaw::FirstOrderArrhenius;
  candidate.sodium.attemptFrequencyHz=configuration.sodiumAttemptFrequencyHz;
  candidate.sodium.activationEnergyJ=
      configuration.sodiumActivationEnergyEv*kElectronVoltJ;
  candidate.sodium.stickingTemperatureK=
      configuration.sodiumStickingTemperatureK;
  candidate.sodium.stickingProbability=
      configuration.sodiumStickingProbability;
  candidate.sodium.citation=configuration.sodiumSurfaceInteractionReference;
  candidate.argon.species=SurfacePhysics::Species::Argon;
  candidate.argon.residenceLaw=
      SurfacePhysics::ResidenceLaw::ModifiedArrhenius;
  candidate.argon.activationEnergyJmol=configuration.argonDesorptionEnergyJmol;
  candidate.argon.stickingTemperatureK=
      configuration.argonStickingTemperatureK;
  candidate.argon.stickingProbability=
      configuration.argonStickingProbability;
  candidate.argon.citation=configuration.argonSurfaceInteractionReference;

  // M2 starts with no adsorbed material unless the receipt selects a
  // versioned checkpoint.  This prevents the legacy compiled 3e16 m^-2
  // abundance from becoming an undocumented initial condition.
  for (int species : std::array<int,2>{{_NA_SPEC_,_AR_SPEC_}}) {
    if (species<0 || species>=PIC::nTotalSpecies) continue;
    for (std::size_t element=0;element<count;++element)
      logicalSurface->SurfaceElementPopulation[species][element]=0.0;
  }
  if (configuration.surfaceInventoryInitialState=="checkpoint") {
    for (const auto& entry : std::array<std::pair<int,const char*>,2>{{
             {_NA_SPEC_,"Na"},{_AR_SPEC_,"Ar"}}}) {
      if (entry.first<0 || entry.first>=PIC::nTotalSpecies) continue;
      SurfacePhysics::Inventory inventory;
      double epoch=0.0;
      if (!inventory.ReadCheckpoint(CheckpointPath(
              configuration.surfaceInventoryCheckpointFile,entry.second),
              entry.second,&epoch,error)) return false;
      if (inventory.cells().size()!=count) {
        *error="surface checkpoint logical-grid size does not match this run";
        return false;
      }
      for (std::size_t element=0;element<count;++element)
        logicalSurface->SurfaceElementPopulation[entry.first][element]=
            inventory.cells()[element].adsorbedAtoms;
      candidate.lastCheckpointEtS=epoch;
    }
  }
  else if (configuration.surfaceInventoryInitialState!="zero" &&
          configuration.surfaceInventoryInitialState!="prescribed-uniform") {
    *error="unsupported M2 surface-inventory initial state";
    return false;
  }

  candidate.initialized=true;
  gProduction=std::move(candidate);
#ifndef _NO_SPICE_CALLS_
  const auto failPublishedInitialization=[]() {
    // Spin-up needs the candidate published while it exercises the ordinary
    // production driver. A failed trial must nevertheless leave no usable
    // "uniform" proof/export behind for a caller that mishandles false.
    auto& provider=Exosphere::SourceProcesses::PhotonStimulatedDesorption::RuntimePhotonFluxFactor;
    if (provider==ProductionPhotonFluxFactor) provider=nullptr;
    gProduction=ProductionState{};
    return false;
  };
#endif
  // Explicit initialization precedes the first weight/source-rate query;
  // J2000 zero is not an implicit source epoch. Later particle steps update
  // this through the same BeginProductionParticleStep entry point.
#ifndef _NO_SPICE_CALLS_
  SpiceDouble startEt=0;
  str2et_c(configuration.simulationStartUtc.c_str(),&startEt);
  if (!std::isfinite(startEt)) { *error="invalid M2 source start epoch";return failPublishedInitialization(); }
  if (gProduction.mode==Mode::ThermalModel) {
    gProduction.currentThermalEpochS=startEt;
    gProduction.thermalInitialization.mode="declared-uniform";
    gProduction.thermalInitialization.ephemerisTimeS=startEt;
    gProduction.thermalInitializationAvailable=true;
    gProduction.initialThermalStateStillApplied=true;
    if (configuration.thermalInitializationMode=="periodic-spinup") {
      SpiceDouble cycleEndEt=0;
      str2et_c(configuration.thermalSpinupCycleEndUtc.c_str(),&cycleEndEt);
      // UTC spelling can differ while representing the same instant. Enforce
      // identity after SPICE conversion, before a potentially long spin-up.
      if (!std::isfinite(cycleEndEt) || cycleEndEt!=startEt) {
        *error="periodic-spinup requires thermal_spinup_cycle_end_utc and simulation_start_utc to represent exactly the same SPICE TDB epoch; declare the intended physical interval explicitly";
        return failPublishedInitialization();
      }
      SpinUpMetrics spinup;
      if (!SpinUpProductionState(&spinup,error)) return failPublishedInitialization();
    } else if (configuration.thermalInitializationMode!="declared-uniform") {
      *error="unsupported predictive thermal initialization mode";return failPublishedInitialization();
    }
  } else if (configuration.thermalInitializationMode!="declared-uniform") {
    *error="periodic-spinup initialization requires the predictive thermal-model";
    return failPublishedInitialization();
  }
  if (!RefreshPrescribedTemperatureCache(startEt,error)) return failPublishedInitialization();
  SurfaceEnvironment::BeginProductionParticleStep(startEt);
#endif
  return true;
}

bool ProductionStateInitialized() { return gProduction.initialized; }

bool GetProductionThermalInitialization(ProductionThermalInitializationInfo* info,
                                       std::string* error) {
  if (info==nullptr || error==nullptr) return false;
  error->clear();
  if (!gProduction.initialized || gProduction.mode!=Mode::ThermalModel ||
      !gProduction.thermalInitializationAvailable) {
    *error="applied thermal initialization is unavailable outside an initialized predictive thermal model";
    return false;
  }
  *info=gProduction.thermalInitialization;
  return true;
}

bool WriteProductionInitialThermalStateCsv(const std::string& path,
                                          std::string* error) {
  if (error==nullptr) return false;
  error->clear();
  if (!gProduction.initialized || gProduction.mode!=Mode::ThermalModel ||
      !gProduction.thermalInitializationAvailable ||
      !gProduction.initialThermalStateStillApplied ||
      gProduction.currentThermalEpochS!=gProduction.thermalInitialization.ephemerisTimeS ||
      gProduction.columns.empty()) {
    *error="initial thermal-state export requires the actual applied t0 predictive columns before evolution; reinitialize in a new run rather than exporting a later state";
    return false;
  }
  // Validate the full state before creating an evidence file. Layer centres
  // are positive depth [m]; the distinct zero-depth surface has no heat mass.
  for (const auto& column:gProduction.columns) {
    const auto& depths=column.LayerDepthM();
    const auto& temperatures=column.LayerTemperatureK();
    if (!PositiveFinite(column.SurfaceTemperatureK()) || depths.empty() ||
        depths.size()!=temperatures.size()) {
      *error="initial thermal-state export encountered incomplete production column arrays";return false;
    }
    double previousDepth=0;
    for (std::size_t layer=0;layer<depths.size();++layer) {
      if (!PositiveFinite(depths[layer]) || depths[layer]<=previousDepth ||
          !PositiveFinite(temperatures[layer])) {
        *error="initial thermal-state export encountered invalid depth/kinetic temperature";return false;
      }
      previousDepth=depths[layer];
    }
  }
  if (path.empty()) { *error="initial thermal-state export path is empty";return false; }
  // O_EXCL enforces the no-overwrite contract atomically, even if two faulty
  // launchers accidentally name the same path. Only the native rank-zero
  // adapter calls this writer; no MPI-rank or later-run replacement is allowed.
  const int descriptor=::open(path.c_str(),O_WRONLY|O_CREAT|O_EXCL,0666);
  if (descriptor<0) {
    *error="cannot exclusively create initial thermal-state CSV '"+path+"' (it may already exist); use a new campaign path and inspect permissions/free disk";
    return false;
  }
  FILE* output=::fdopen(descriptor,"w");
  if (output==nullptr) {
    ::close(descriptor);
    *error="cannot open the exclusively created initial thermal-state CSV stream; preserve the partial artifact and use a new path";
    return false;
  }
  bool good=std::fprintf(output,"logical_surface_element,depth_index,depth_m,kinetic_temperature_K\n")>0;
  for (std::size_t element=0;good && element<gProduction.columns.size();++element) {
    const auto& column=gProduction.columns[element];
    good=std::fprintf(output,"%zu,-1,0,%.17g\n",element,column.SurfaceTemperatureK())>0;
    const auto& depths=column.LayerDepthM();
    const auto& temperatures=column.LayerTemperatureK();
    for (std::size_t layer=0;good && layer<depths.size();++layer)
      good=std::fprintf(output,"%zu,%zu,%.17g,%.17g\n",element,layer,depths[layer],temperatures[layer])>0;
  }
  if (std::fclose(output)!=0) good=false;
  if (!good) {
    *error="writing initial thermal-state CSV failed; preserve the partial file, inspect free disk/quota and use a new campaign path";
    return false;
  }
  return true;
}

bool RefreshProductionPhotonFluxCache(double sourceBeginEtS,std::string* error) {
  if (error==nullptr) return false;
  error->clear();
  if (!gProduction.initialized || !std::isfinite(sourceBeginEtS)) {
    *error="M2 photon-flux cache requires initialized state and an explicit finite source epoch";return false;
  }
  if (gProduction.photonFluxEpochS==sourceBeginEtS &&
      gProduction.photonFluxFactors.size()==gProduction.temperaturesK.size()) return true;
#ifdef _NO_SPICE_CALLS_
  *error="terrain-aware PSD forcing requires SPICE";return false;
#else
  EpochIlluminationGeometry epoch;
  if (!EvaluateEpochIlluminationGeometry(sourceBeginEtS,&epoch,error)) return false;
  std::vector<double> factors(gProduction.temperaturesK.size());
  for (std::size_t element=0;element<factors.size();++element) {
    ElementIlluminationGeometry geometry;
    if (!EvaluateElementIlluminationGeometry(element,epoch,&geometry,error,false)) return false;
    // The core's published photon normalization uses the AMPS AU constant.
    // Preserve that normalization; the actual point distance is SI metres.
    const double distanceM=geometry.pointSunDistanceM;
    factors[element]=geometry.projectedSolarIncidence*std::pow(_AU_/distanceM,2);
    if (!std::isfinite(factors[element]) || factors[element]<0) {
      *error="M2 direct-photon factor is negative/nonfinite";return false;
    }
  }
  gProduction.photonFluxFactors=std::move(factors);
  gProduction.photonFluxEpochS=sourceBeginEtS;
  return true;
#endif
}

double ProductionPhotonFluxFactor(int surfaceElement,void* logicalSurface) {
  if (!gProduction.initialized || logicalSurface!=gProduction.logicalSurface ||
      surfaceElement<0 || static_cast<std::size_t>(surfaceElement)>=gProduction.photonFluxFactors.size() ||
      gProduction.photonFluxEpochS!=SurfaceEnvironment::ProductionSourceEpochS())
    exit(__LINE__,__FILE__,"PSD source geometry cache is unavailable or stale; initialize at the declared source-BEGIN epoch before source rate/CDF evaluation");
  return gProduction.photonFluxFactors[surfaceElement];
}

double LastAdvanceMaximumSurfaceTemperatureChangeK() {
  return gProduction.initialized ?
      gProduction.lastAdvanceMaximumSurfaceTemperatureChangeK:
      std::numeric_limits<double>::quiet_NaN();
}

double ProductionTemperatureK(const double* positionLsoM) {
  if (!gProduction.initialized || positionLsoM==nullptr)
    return std::numeric_limits<double>::quiet_NaN();
#ifdef _NO_SPICE_CALLS_
  // InitializeProductionState() prevents this branch in a supported run.
  // Keep a defensive NaN rather than interpreting LSO coordinates as lunar
  // body-fixed coordinates.
  return std::numeric_limits<double>::quiet_NaN();
#else
  SpiceDouble transform[6][6];
  sxform_c(Frames::SolarOrbital,Frames::BodyFixed,OrbitalMotion::et,transform);
  const auto body=TransformPosition(transform,positionLsoM);
  long int zenith=0,azimuth=0;
  gProduction.logicalSurface->GetSurfaceElementProjectionIndex(
      const_cast<double*>(body.data()),zenith,azimuth);
  const std::size_t element=static_cast<std::size_t>(
      gProduction.logicalSurface->GetLocalSurfaceElementNumber(zenith,azimuth));
  if (element>=gProduction.temperaturesK.size())
    return std::numeric_limits<double>::quiet_NaN();

  if (gProduction.mode==Mode::ThermalModel)
    return gProduction.temperaturesK[element];
  EpochIlluminationGeometry epoch;
  double temperature=0;
  std::string error;
  if (!EvaluateEpochIlluminationGeometry(OrbitalMotion::et,&epoch,&error) ||
      !PrescribedTemperatureAtBodyPoint(body,epoch,&temperature,&error))
    exit(__LINE__,__FILE__,("M2 prescribed kinetic-temperature lookup failed: "+error).c_str());
  return temperature;
#endif
}

bool AdvanceProductionState(double ephemerisTimeBeginS,
                            double ephemerisTimeEndS,
                            std::string* error) {
  if (error==nullptr) return false;
  error->clear();
  if (!gProduction.initialized) return true;
#ifdef _NO_SPICE_CALLS_
  *error="M2 thermal-state advancement requires a SPICE-enabled build";
  return false;
#else
  const double duration=ephemerisTimeEndS-ephemerisTimeBeginS;
  if (!(duration>=0.0) || !std::isfinite(duration)) {
    *error="invalid M2 thermal-state epoch interval";
    return false;
  }
  gProduction.lastAdvanceMaximumSurfaceTemperatureChangeK=0.0;
  // Prescribed boundaries have no evolving subsurface state, but their
  // kinetic source temperatures still follow the explicit end epoch. Keep
  // desorption/cohort diagnostics synchronized even when duration is zero.
  if (gProduction.mode!=Mode::ThermalModel)
    return RefreshPrescribedTemperatureCache(ephemerisTimeEndS,error);
  if (duration==0.0) return true;
  // Any real evolution invalidates reuse/export of the previous exact state,
  // including an interval whose end numerically equals an older cached epoch.
  gProduction.completedSpinupAtCurrentState=false;
  gProduction.initialThermalStateStillApplied=false;
  const Runtime::Configuration& configuration=*gProduction.configuration;
  const int substeps=static_cast<int>(std::ceil(duration/
      configuration.thermalTimeStepS));
  const double timeStep=duration/substeps;
  for (int step=0;step<substeps;++step) {
    const double epoch=ephemerisTimeBeginS+(step+0.5)*timeStep;
    EpochIlluminationGeometry epochGeometry;
    if (!EvaluateEpochIlluminationGeometry(epoch,&epochGeometry,error))
      return false;
    for (std::size_t element=0;element<gProduction.columns.size();++element) {
      ElementIlluminationGeometry elementGeometry;
      if (!EvaluateElementIlluminationGeometry(element,epochGeometry,
              &elementGeometry,error)) return false;
      const double previousSurfaceTemperatureK=
          gProduction.columns[element].SurfaceTemperatureK();
      if (!gProduction.columns[element].AdvanceWithAbsorbedFlux(timeStep,
              elementGeometry.absorbedSolarFluxWm2,error)) return false;
      gProduction.temperaturesK[element]=
          gProduction.columns[element].SurfaceTemperatureK();
      gProduction.lastAdvanceMaximumSurfaceTemperatureChangeK=std::max(
          gProduction.lastAdvanceMaximumSurfaceTemperatureChangeK,
          std::abs(gProduction.temperaturesK[element]-
                   previousSurfaceTemperatureK));
      gProduction.geometricShadow[element]=
          elementGeometry.terrainSunVisible ? 0 : 1;
    }
  }
  gProduction.currentThermalEpochS=ephemerisTimeEndS;
  return true;
#endif
}

bool SampleProductionSurfaceEnvironment(double ephemerisTimeS,
    std::vector<SurfaceEnvironmentSample>* samples,std::string* error) {
  if (samples==nullptr || error==nullptr) return false;
  samples->clear();
  error->clear();
  if (!gProduction.initialized || gProduction.mode!=Mode::ThermalModel ||
      gProduction.logicalSurface==nullptr) {
    *error="surface-environment sampling requires an initialized production "
        "thermal model";
    return false;
  }
#ifdef _NO_SPICE_CALLS_
  (void)ephemerisTimeS;
  *error="surface-environment sampling requires a SPICE-enabled build";
  return false;
#else
  EpochIlluminationGeometry epochGeometry;
  if (!EvaluateEpochIlluminationGeometry(ephemerisTimeS,&epochGeometry,error))
    return false;
  samples->reserve(gProduction.columns.size());
  for (std::size_t element=0;element<gProduction.columns.size();++element) {
    ElementIlluminationGeometry forcing;
    if (!EvaluateElementIlluminationGeometry(element,epochGeometry,
            &forcing,error)) return false;
    SurfaceEnvironmentSample sample;
    sample.surfaceElement=element;
    gProduction.logicalSurface->GetSurfaceElementIndex(
        sample.zenithIndex,sample.azimuthIndex,static_cast<int>(element));
    sample.surfacePositionBodyFixedM=forcing.surfacePositionM;
    sample.terrainNormalBodyFixed=forcing.terrainNormal;
    const double radius=Norm(forcing.surfacePositionM);
    if (!(radius>0.0)) {
      *error="sampled production terrain point has invalid radius";
      return false;
    }
    sample.latitudeDeg=std::asin(forcing.surfacePositionM[2]/radius)*180.0/kPi;
    sample.longitudeEastDeg=NormalizeDegrees(std::atan2(
        forcing.surfacePositionM[1],forcing.surfacePositionM[0])*180.0/kPi);
    sample.localSolarTimeHour=LocalSolarTimeHour(
        forcing.surfacePositionM,epochGeometry.unitSunBodyFixed);
    sample.solarIncidenceCosine=forcing.solarIncidenceCosine;
    sample.sunAboveLocalHorizon=forcing.sunAboveLocalHorizon;
    sample.terrainSunVisible=forcing.terrainSunVisible;
    sample.terrainOccluded=forcing.terrainOccluded;
    sample.solarDiscVisibleFraction=forcing.solarDiscVisibleFraction;
    sample.directBeamFraction=forcing.directBeamFraction;
    sample.heliocentricDistanceAu=forcing.heliocentricDistanceAu;
    sample.absorbedSolarFluxWm2=forcing.absorbedSolarFluxWm2;
    sample.surfaceTemperatureK=
        gProduction.columns[element].SurfaceTemperatureK();
    sample.lastEnergyResidualJm2=
        gProduction.columns[element].LastEnergyResidualJm2();
    sample.lastNonlinearIterations=
        gProduction.columns[element].LastNonlinearIterations();
    if (!std::isfinite(sample.absorbedSolarFluxWm2) ||
        !PositiveFinite(sample.surfaceTemperatureK)) {
      *error="sampled production thermal state is non-finite";
      return false;
    }
    samples->push_back(sample);
  }
  return true;
#endif
}

bool ProductionThermalProfile(std::size_t surfaceElement,
    std::vector<double>* layerDepthM,std::vector<double>* layerTemperatureK,
    std::string* error) {
  if (layerDepthM==nullptr || layerTemperatureK==nullptr || error==nullptr)
    return false;
  layerDepthM->clear();
  layerTemperatureK->clear();
  error->clear();
  if (!gProduction.initialized || gProduction.mode!=Mode::ThermalModel ||
      surfaceElement>=gProduction.columns.size()) {
    *error="thermal-profile request is outside the initialized production "
        "surface";
    return false;
  }
  *layerDepthM=gProduction.columns[surfaceElement].LayerDepthM();
  *layerTemperatureK=
      gProduction.columns[surfaceElement].LayerTemperatureK();
  if (layerDepthM->empty() || layerDepthM->size()!=layerTemperatureK->size()) {
    *error="production thermal column has inconsistent layer arrays";
    return false;
  }
  return true;
}

bool SpinUpProductionState(SpinUpMetrics* metrics,std::string* error) {
  if (metrics==nullptr || error==nullptr) return false;
  *metrics=SpinUpMetrics{};
  error->clear();
  if (!gProduction.initialized || gProduction.mode!=Mode::ThermalModel) {
    *error="periodic spin-up requires an initialized production thermal model";
    return false;
  }
#ifdef _NO_SPICE_CALLS_
  *error="periodic thermal spin-up requires a SPICE-enabled build";
  return false;
#else
  const Runtime::Configuration& configuration=*gProduction.configuration;
  SpiceDouble cycleStart=0.0,cycleEnd=0.0;
  str2et_c(configuration.thermalSpinupCycleStartUtc.c_str(),&cycleStart);
  str2et_c(configuration.thermalSpinupCycleEndUtc.c_str(),&cycleEnd);
  if (!(cycleEnd>cycleStart) || !std::isfinite(cycleStart) ||
      !std::isfinite(cycleEnd)) {
    *error="thermal spin-up UTC interval is not a positive CSPICE interval";
    return false;
  }
  metrics->cycleStartEphemerisTimeS=cycleStart;
  metrics->cycleEndEphemerisTimeS=cycleEnd;
  if (gProduction.completedSpinupAtCurrentState &&
      gProduction.currentThermalEpochS==cycleEnd &&
      gProduction.completedSpinup.cycleStartEphemerisTimeS==cycleStart &&
      gProduction.completedSpinup.cycleEndEphemerisTimeS==cycleEnd) {
    // Auto-initialization already solved this exact state. I17/I20 may ask
    // explicitly for its diagnostics without paying for another spin-up or
    // changing the initial columns used by sources and the proof exporter.
    *metrics=gProduction.completedSpinup;
    return true;
  }

  // Compare every finite-volume node at consecutive cycle boundaries.  A
  // surface-only check can converge while the deeper heat reservoir still
  // remembers the arbitrary initial temperature and biases a later eclipse.
  // The state copy is one additional production-column array, independent of
  // cycle count; std::new_handler supplies the approved diagnostic if this
  // allocation exhausts the job memory.
  for (;;) {
    const std::vector<Thermal::Column> previous=gProduction.columns;
    if (!AdvanceProductionState(cycleStart,cycleEnd,error)) return false;
    ++metrics->completedCycles;
    double maximumTemperatureChange=0.0;
    double maximumEnthalpyChange=0.0;
    for (std::size_t element=0;element<gProduction.columns.size();++element) {
      maximumTemperatureChange=std::max(maximumTemperatureChange,std::abs(
          gProduction.columns[element].SurfaceTemperatureK()-
          previous[element].SurfaceTemperatureK()));
      const auto& currentLayers=gProduction.columns[element].LayerTemperatureK();
      const auto& previousLayers=previous[element].LayerTemperatureK();
      if (currentLayers.size()!=previousLayers.size()) {
        *error="thermal spin-up changed the configured layer count";
        return false;
      }
      for (std::size_t layer=0;layer<currentLayers.size();++layer)
        maximumTemperatureChange=std::max(maximumTemperatureChange,
            std::abs(currentLayers[layer]-previousLayers[layer]));
      maximumEnthalpyChange=std::max(maximumEnthalpyChange,std::abs(
          gProduction.columns[element].TotalSensibleEnthalpyJm2()-
          previous[element].TotalSensibleEnthalpyJm2()));
    }
    metrics->maximumNodeTemperatureChangeK=maximumTemperatureChange;
    metrics->maximumColumnEnthalpyChangeJm2=maximumEnthalpyChange;
    if (maximumTemperatureChange<=
            configuration.thermalSpinupTemperatureToleranceK &&
        maximumEnthalpyChange<=
            configuration.thermalSpinupEnergyToleranceJm2) {
      gProduction.completedSpinup=*metrics;
      gProduction.completedSpinupAtCurrentState=true;
      gProduction.thermalInitialization.mode="periodic-spinup";
      gProduction.thermalInitialization.ephemerisTimeS=cycleEnd;
      gProduction.thermalInitialization.spinup=*metrics;
      gProduction.thermalInitializationAvailable=true;
      gProduction.initialThermalStateStillApplied=true;
      return true;
    }
  }
#endif
}

double ProductionElementTemperatureK(int element) {
  if (!gProduction.initialized || element<0 ||
      static_cast<std::size_t>(element)>=gProduction.temperaturesK.size())
    exit(__LINE__,__FILE__,"invalid logical element in passive temperature diagnostic");
  if (!std::isfinite(gProduction.temperaturesK[element]) || gProduction.temperaturesK[element]<=0)
    exit(__LINE__,__FILE__,"logical element has no supported kinetic surface temperature; preserve the Diviner coverage gap and use a qualified predictive thermal boundary or restrict the model support explicitly");
  return gProduction.temperaturesK[element];
}

double ProductionRetentionTimescaleS() {
  if (!gProduction.initialized)
    exit(__LINE__,__FILE__,"retention timescale requires initialized production surface state");
  return gProduction.configuration->coldTrapRetentionTimescaleS;
}

bool SampleProductionElementKinetics(int element,int species,
    double* residenceTimeS,bool* thermallyRetained,std::string* error) {
  if (error==nullptr || residenceTimeS==nullptr || thermallyRetained==nullptr)
    return false;
  error->clear();
  if (!gProduction.initialized || element<0 ||
      static_cast<std::size_t>(element)>=gProduction.temperaturesK.size()) {
    *error="surface kinetics diagnostic needs an initialized, available logical element";
    return false;
  }
  const auto* parameters=ParametersForSpecies(species);
  if (parameters==nullptr) {
    *error="surface kinetics diagnostic requested an unconfigured species";
    return false;
  }
  if (!std::isfinite(gProduction.temperaturesK[element]) || gProduction.temperaturesK[element]<=0) {
    *error="surface kinetics unavailable at a real Diviner-map coverage gap; no initialization temperature is substituted";
    return false;
  }
  const double residence=SurfacePhysics::ResidenceTimeS(
      *parameters,gProduction.temperaturesK[element]);
  const double timescale=gProduction.configuration->coldTrapRetentionTimescaleS;
  if (std::isnan(residence) || residence<0 || !std::isfinite(timescale) || timescale<=0) {
    *error="invalid production residence time or configured retention-reporting timescale";
    return false;
  }
  *residenceTimeS=residence;
  *thermallyRetained=residence>=timescale;
  return true;
}

double ProductionThermalDesorptionRate(int species,int surfaceElement,
                                       void* sphereDataPointer) {
  if (!gProduction.initialized) return 0.0;
  const SurfacePhysics::SpeciesParameters* parameters=
      ParametersForSpecies(species);
  if (parameters==nullptr) return 0.0;
  if (surfaceElement<0 || static_cast<std::size_t>(surfaceElement)>=
          gProduction.temperaturesK.size() || sphereDataPointer==nullptr)
    exit(__LINE__,__FILE__,"invalid M2 thermal-desorption surface element");
  auto* sphere=static_cast<cInternalSphericalData*>(sphereDataPointer);
  const double population=
      sphere->SurfaceElementPopulation[species][surfaceElement];
  if (!std::isfinite(population) || population<0)
    exit(__LINE__,__FILE__,"M2 thermal desorption requires finite nonnegative physical surface stock [atoms]");
  // An unoccupied cell has exactly zero release even where the measured map
  // has no coverage. Positive stock must never silently use an invented or
  // stale initialization temperature to evaluate its residence law.
  if (population==0) return 0;
  if (!std::isfinite(gProduction.temperaturesK[surfaceElement]) ||
      gProduction.temperaturesK[surfaceElement]<=0)
    exit(__LINE__,__FILE__,"M2 thermal desorption has positive stock in a Diviner temperature coverage gap; use a qualified predictive thermal mode or explicitly restrict sources/model support; no gap filling is performed");
  const double residence=SurfacePhysics::ResidenceTimeS(
      *parameters,gProduction.temperaturesK[surfaceElement]);
  if (population<0.0 || std::isnan(residence) || residence<=0.0)
    exit(__LINE__,__FILE__,"invalid M2 surface population/residence time");
  if (std::isinf(residence)) return 0.0;
  return population/residence;
}

bool ProductionStickingProbability(int species,double temperatureK,
                                   double* probability,std::string* error) {
  if (probability==nullptr || error==nullptr) return false;
  error->clear();
  if (!gProduction.initialized) {
    *error="M2 production surface state is not initialized";
    return false;
  }
#if _NA_PLUS_SPEC_ >= 0
  if (species==_NA_PLUS_SPEC_) {
    *error="Na+ surface interaction is unconfigured; no neutral sticking law "
           "is applied to ions";
    return false;
  }
#endif
#if _AR_PLUS_SPEC_ >= 0
  if (species==_AR_PLUS_SPEC_) {
    *error="Ar+ surface interaction is unconfigured; no neutral sticking law "
           "is applied to ions";
    return false;
  }
#endif
  const SurfacePhysics::SpeciesParameters* parameters=
      ParametersForSpecies(species);
  if (parameters==nullptr) {
    *probability=0.0;
    return true;
  }
  *probability=SurfacePhysics::StickingProbability(*parameters,temperatureK);
  if (!std::isfinite(*probability) || *probability<0.0 || *probability>1.0) {
    *error="configured M2 sticking law returned an invalid probability";
    return false;
  }
  return true;
}

bool ProductionAccommodationCoefficient(int species,double temperatureK,
                                        double* coefficient,
                                        Exosphere::SurfaceInteraction::
                                            AccommodationSpeedLaw* speedLaw,
                                        std::string* error) {
  if (coefficient==nullptr || speedLaw==nullptr || error==nullptr) return false;
  error->clear();
  if (!gProduction.initialized || !std::isfinite(temperatureK) ||
      temperatureK<=0.0) {
    *error="M2 accommodation requested without a valid production thermal "
           "state";
    return false;
  }
#if _AR_SPEC_ >= 0
  if (species==_AR_SPEC_) {
    *coefficient=gProduction.configuration->argonAccommodationCoefficient;
    const std::string& model=
        gProduction.configuration->argonAccommodationSpeedModel;
    *speedLaw=model=="kinetic-energy" ?
        Exosphere::SurfaceInteraction::AccommodationSpeedLaw::KineticEnergy :
        Exosphere::SurfaceInteraction::AccommodationSpeedLaw::LinearSpeed;
    return model=="kinetic-energy" || model=="linear-speed";
  }
#endif
#if _NA_SPEC_ >= 0
  if (species==_NA_SPEC_) {
    // Tenishev et al. (2013), Eq. 14, defines omega as the fraction used to
    // interpolate incident and thermal kinetic energies.  The input still
    // names both the coefficient and mathematical model so an investigator
    // can select a different published law without recompiling AMPS.
    *coefficient=gProduction.configuration->sodiumAccommodationCoefficient;
    const std::string& model=
        gProduction.configuration->sodiumAccommodationSpeedModel;
    *speedLaw=model=="kinetic-energy" ?
        Exosphere::SurfaceInteraction::AccommodationSpeedLaw::KineticEnergy :
        Exosphere::SurfaceInteraction::AccommodationSpeedLaw::LinearSpeed;
    return model=="kinetic-energy" || model=="linear-speed";
  }
#endif
  *error="surface accommodation is configured only for M2 neutral Na and Ar";
  return false;
}

bool UseProductionFluxBalancedInventory() {
  // Tenishev et al. (2013) Case 1 prescribes a time-independent uniform Na
  // surface abundance and permanently removes returned particles.  It is not
  // the finite flux-balanced M2 inventory used by migration/cold-trap runs.
  return gProduction.initialized &&
      gProduction.configuration->surfaceInventoryInitialState!=
          "prescribed-uniform";
}

bool WriteProductionInventoryCheckpoint(double ephemerisTimeS,
                                        std::string* error) {
  if (error==nullptr) return false;
  error->clear();
  if (!gProduction.initialized) return true;
  const Runtime::Configuration& configuration=*gProduction.configuration;
  // MoonInput normalizes the literal `none` to an empty string.  Test the
  // normalized representation here; otherwise an intentionally disabled
  // checkpoint would create the relative files `.Na` and `.Ar` in the run
  // directory, silently defeating the run-receipt setting.
  if (configuration.surfaceInventoryCheckpointFile.empty()) return true;
  if (std::isfinite(gProduction.lastCheckpointEtS) &&
      ephemerisTimeS-gProduction.lastCheckpointEtS<
          configuration.surfaceInventoryCheckpointIntervalS) return true;
  const std::size_t count=gProduction.temperaturesK.size();
  for (const auto& entry : std::array<std::pair<int,const char*>,2>{{
           {_NA_SPEC_,"Na"},{_AR_SPEC_,"Ar"}}}) {
    if (entry.first<0 || entry.first>=PIC::nTotalSpecies) continue;
    SurfacePhysics::Inventory inventory;
    if (!inventory.Initialize(count,0.0,error)) return false;
    for (std::size_t element=0;element<count;++element) {
      double residence=0;
      bool retained=false;
      if (!SampleProductionElementKinetics(static_cast<int>(element),entry.first,
              &residence,&retained,error)) return false;
      if (!inventory.SetCellState(element,
              gProduction.logicalSurface->SurfaceElementPopulation
                  [entry.first][element],
              gProduction.geometricShadow[element]!=0,
              retained,
              error)) return false;
    }
    if (!inventory.WriteCheckpoint(CheckpointPath(
            configuration.surfaceInventoryCheckpointFile,entry.second),
            entry.second,ephemerisTimeS,error)) return false;
  }
  gProduction.lastCheckpointEtS=ephemerisTimeS;
  return true;
}

}  // namespace SurfaceTemperature
}  // namespace Moon
