#include "ForwardModelAMPS.h"

#include "Enceladus.h"
#include "constants.h"

#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <utility>
#include <vector>

namespace Enceladus {
namespace ForwardModel {
namespace AMPS {
namespace {

namespace fs = std::filesystem;

bool gInstalled = false;
bool gNumericsFinalized = false;
Configuration gConfiguration;
SourceDistribution gDistribution;
NumericalPlan gPlan;
std::map<std::string, int> gSpeciesIndex;
std::uint64_t gInjectionStep = 0;
std::vector<unsigned long long> gLocalInjected;
std::vector<unsigned long long> gLocalPhotolyticLoss;
std::vector<unsigned long long> gLocalSurfaceLoss;
std::vector<unsigned long long> gLocalDomainEscape;
thread_local long int gSurfaceAbsorptionParticle = -1;
constexpr double kBoltzmann = 1.380649e-23;

// A flyby point is an observation location, not a moving simulation
// particle. Positions and ellipsoidal altitude are supplied by the frozen
// CSPICE trajectory in the declared body-fixed frame; AMPS evaluates its
// completed, time-averaged cell sampling buffer at each position.
struct FlybyPoint {
  double timeFromReferenceS = 0.0;
  double positionM[3] = {0.0, 0.0, 0.0};
  double altitudeM = 0.0;
};

std::vector<FlybyPoint> gFlybyTrajectory;
std::vector<FlybyPoint> gFlybyRefinementTrajectory;

struct ObservationPoint {
  double timeS = 0.0;
  double positionM[3] = {0.0, 0.0, 0.0};
  double observedValue = 0.0;
  double observedSigma = 0.0;
  double rayDirection[3] = {0.0, 0.0, 0.0};
  double observerCoordinateM = 0.0;
};

struct LoadedObservation {
  ObservationRequest request;
  std::vector<ObservationPoint> points;
};

std::vector<LoadedObservation> gObservations;

void Increment(std::vector<unsigned long long>* counters, int species) {
#if _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
#pragma omp atomic update
#endif
  ++(*counters)[species];
}

[[noreturn]] void Stop(const std::string& reason) {
  std::cerr << "[srcEnceladus] FATAL: " << reason << '\n';
  MPI_Abort(MPI_GLOBAL_COMMUNICATOR, 71);
  std::abort();
}

std::uint64_t Mix(std::uint64_t value) {
  value += 0x9e3779b97f4a7c15ULL;
  value = (value ^ (value >> 30U)) * 0xbf58476d1ce4e5b9ULL;
  value = (value ^ (value >> 27U)) * 0x94d049bb133111ebULL;
  return value ^ (value >> 31U);
}

class DeterministicRandom {
 public:
  explicit DeterministicRandom(std::uint64_t seed) : state_(seed) {}

  double Uniform() {
    state_ = Mix(state_);
    // Open interval (0,1): the half-unit offset makes log(U) finite and the
    // 53-bit projection is exactly reproducible on every MPI rank.
    return (static_cast<double>(state_ >> 11U) + 0.5) /
           9007199254740992.0;
  }

 private:
  std::uint64_t state_;
};

double Norm(const Vec3& value) {
  return std::sqrt(value.x*value.x + value.y*value.y + value.z*value.z);
}

Vec3 Normalize(const Vec3& value) {
  const double norm = Norm(value);
  return Vec3{value.x/norm, value.y/norm, value.z/norm};
}

double Dot(const Vec3& a, const Vec3& b) {
  return a.x*b.x + a.y*b.y + a.z*b.z;
}

std::vector<std::string> SplitCsv(const std::string& line) {
  std::vector<std::string> fields;
  std::stringstream stream(line);
  std::string field;
  while (std::getline(stream, field, ',')) fields.push_back(field);
  return fields;
}

bool ParseCsvFinite(const std::string& text, double* value) {
  if (text.empty()) return false;
  char* end = nullptr;
  errno = 0;
  const double parsed = std::strtod(text.c_str(), &end);
  if (errno == ERANGE || end == text.c_str() || *end != '\0' ||
      !std::isfinite(parsed)) return false;
  *value = parsed;
  return true;
}

std::vector<std::string> SplitCsvQuoted(const std::string& line) {
  std::vector<std::string> fields;
  std::string field;
  bool quoted = false;
  for (std::size_t i = 0; i < line.size(); ++i) {
    const char c = line[i];
    if (c == '"') {
      if (quoted && i+1 < line.size() && line[i+1] == '"') {
        field.push_back('"');
        ++i;
      } else {
        quoted = !quoted;
      }
    } else if (c == ',' && !quoted) {
      fields.push_back(field);
      field.clear();
    } else {
      field.push_back(c);
    }
  }
  if (quoted) return {};
  fields.push_back(field);
  return fields;
}

Status CsvColumns(const std::string& header,
                  std::map<std::string, std::size_t>* columns) {
  columns->clear();
  const std::vector<std::string> names = SplitCsvQuoted(header);
  if (names.empty()) return Status::Error("CSV header is empty or malformed");
  for (std::size_t i = 0; i < names.size(); ++i)
    if (names[i].empty() || !columns->emplace(names[i], i).second)
      return Status::Error("CSV header contains an empty or duplicate column");
  return Status::Ok();
}

Status RequireCsvColumn(const std::map<std::string, std::size_t>& columns,
                        const std::string& name, std::size_t* index) {
  const auto found = columns.find(name);
  if (found == columns.end())
    return Status::Error("CSV is missing required column '" + name + "'");
  *index = found->second;
  return Status::Ok();
}

Status LoadObservation(const ObservationRequest& request,
                       LoadedObservation* output) {
  output->request = request;
  output->points.clear();
  std::ifstream input(request.dataFile);
  if (!input)
    return Status::Error("cannot open observation data_file: " +
                         request.dataFile);
  std::string line;
  if (!std::getline(input, line))
    return Status::Error("observation CSV is empty: " + request.dataFile);
  std::map<std::string, std::size_t> columns;
  Status status = CsvColumns(line, &columns);
  if (!status.success) return status;

  std::vector<std::string> required;
  if (request.instrument == ObservationInstrument::Inms) {
    required = {request.timeColumn, request.xColumn, request.yColumn,
                request.zColumn};
    if (!request.comparisonSpecies.empty()) {
      required.push_back(request.valueColumn);
      required.push_back(request.uncertaintyColumn);
    } else if (!request.contextObservable.empty()) {
      required.push_back(request.contextValueColumn);
      required.push_back(request.contextUncertaintyColumn);
    }
  } else {
    required = {"event_id", "time_source_s", "time_from_min_altitude_s",
                "tangent_x_m", "tangent_y_m", "tangent_z_m",
                "optical_depth_raw", "poisson_sigma_optical_depth"};
  }
  std::map<std::string, std::size_t> index;
  for (const std::string& name : required) {
    std::size_t column = 0;
    status = RequireCsvColumn(columns, name, &column);
    if (!status.success)
      return Status::Error(request.id + ": " + status.message);
    index[name] = column;
  }

  std::vector<double> rayTimes;
  double previousTime = -std::numeric_limits<double>::infinity();
  std::size_t rowNumber = 1;
  while (std::getline(input, line)) {
    ++rowNumber;
    if (line.empty()) continue;
    const std::vector<std::string> fields = SplitCsvQuoted(line);
    if (fields.size() != columns.size())
      return Status::Error(request.id + " observation row " +
          std::to_string(rowNumber) + " has " +
          std::to_string(fields.size()) + " fields; expected " +
          std::to_string(columns.size()));
    ObservationPoint point;
    if (request.instrument == ObservationInstrument::Inms) {
      if (!ParseCsvFinite(fields[index[request.timeColumn]], &point.timeS) ||
          !ParseCsvFinite(fields[index[request.xColumn]],
                          &point.positionM[0]) ||
          !ParseCsvFinite(fields[index[request.yColumn]],
                          &point.positionM[1]) ||
          !ParseCsvFinite(fields[index[request.zColumn]],
                          &point.positionM[2]))
        return Status::Error(request.id + " observation row " +
            std::to_string(rowNumber) +
            " has invalid time or IAU_ENCELADUS position");
      if (!request.comparisonSpecies.empty()) {
        if (!ParseCsvFinite(fields[index[request.valueColumn]],
                            &point.observedValue) ||
            !ParseCsvFinite(fields[index[request.uncertaintyColumn]],
                            &point.observedSigma) ||
            point.observedValue < 0.0 || !(point.observedSigma > 0.0))
          return Status::Error(request.id + " observation row " +
              std::to_string(rowNumber) +
              " requires nonnegative density and positive uncertainty");
      } else if (!request.contextObservable.empty()) {
        if (!ParseCsvFinite(fields[index[request.contextValueColumn]],
                            &point.observedValue) ||
            !ParseCsvFinite(fields[index[request.contextUncertaintyColumn]],
                            &point.observedSigma) ||
            point.observedValue < 0.0 || !(point.observedSigma > 0.0))
          return Status::Error(request.id + " observation row " +
              std::to_string(rowNumber) +
              " requires nonnegative native observable and positive uncertainty");
      }
    } else {
      if (fields[index["event_id"]] != request.eventId)
        return Status::Error(request.id + " observation row " +
            std::to_string(rowNumber) + " has event_id '" +
            fields[index["event_id"]] + "', expected '" + request.eventId +
            "'");
      double rayTime = 0.0;
      if (!ParseCsvFinite(fields[index["time_source_s"]], &rayTime) ||
          !ParseCsvFinite(fields[index["time_from_min_altitude_s"]],
                          &point.timeS) ||
          !ParseCsvFinite(fields[index["tangent_x_m"]],
                          &point.positionM[0]) ||
          !ParseCsvFinite(fields[index["tangent_y_m"]],
                          &point.positionM[1]) ||
          !ParseCsvFinite(fields[index["tangent_z_m"]],
                          &point.positionM[2]) ||
          !ParseCsvFinite(fields[index["optical_depth_raw"]],
                          &point.observedValue) ||
          !ParseCsvFinite(fields[index["poisson_sigma_optical_depth"]],
                          &point.observedSigma) ||
          !(point.observedSigma > 0.0))
        return Status::Error(request.id + " observation row " +
            std::to_string(rowNumber) + " has invalid UVIS data");
      rayTimes.push_back(rayTime);
    }
    if (!(point.timeS > previousTime))
      return Status::Error(request.id +
          " observation times must be strictly increasing");
    previousTime = point.timeS;
    const double radius = std::sqrt(
        point.positionM[0]*point.positionM[0] +
        point.positionM[1]*point.positionM[1] +
        point.positionM[2]*point.positionM[2]);
    if (!(radius > _RADIUS_(_TARGET_)))
      return Status::Error(request.id + " observation row " +
          std::to_string(rowNumber) +
          " lies on or inside the reference sphere");
    output->points.push_back(point);
  }
  if (output->points.size() < 2)
    return Status::Error(request.id + " requires at least two observation rows");
  if (request.instrument == ObservationInstrument::Inms)
    return Status::Ok();

  if (request.uvisRayGeometry ==
      UvisRayGeometry::InfiniteChordTangentTrackNormal) {
    if (output->points.size() < 3)
      return Status::Error(request.id +
          " tangent-track ray construction requires at least three rows");
    // The SwRI reconstruction defines the stellar line of sight as the unit
    // normal to both the body-centred tangent radius r and its measured track
    // derivative dr/dt.  A three-point Lagrange derivative retains nonuniform
    // sample timing; the ray sign is made continuous but has no effect on the
    // two-sided chord integral.
    const std::size_t count = output->points.size();
    for (std::size_t i = 0; i < count; ++i) {
      std::size_t a = i-1, b = i, c = i+1;
      if (i == 0) { a = 0; b = 1; c = 2; }
      if (i+1 == count) { a = count-3; b = count-2; c = count-1; }
      const double t = rayTimes[i];
      const double ta = rayTimes[a], tb = rayTimes[b], tc = rayTimes[c];
      const double wa = (2.0*t-tb-tc)/((ta-tb)*(ta-tc));
      const double wb = (2.0*t-ta-tc)/((tb-ta)*(tb-tc));
      const double wc = (2.0*t-ta-tb)/((tc-ta)*(tc-tb));
      double velocity[3];
      for (int d = 0; d < 3; ++d)
        velocity[d] = wa*output->points[a].positionM[d] +
                      wb*output->points[b].positionM[d] +
                      wc*output->points[c].positionM[d];
      const double* r = output->points[i].positionM;
      double* direction = output->points[i].rayDirection;
      direction[0] = r[1]*velocity[2]-r[2]*velocity[1];
      direction[1] = r[2]*velocity[0]-r[0]*velocity[2];
      direction[2] = r[0]*velocity[1]-r[1]*velocity[0];
      double norm = std::sqrt(direction[0]*direction[0] +
                              direction[1]*direction[1] +
                              direction[2]*direction[2]);
      if (!(norm > 0.0) || !std::isfinite(norm))
        return Status::Error(request.id + " has a degenerate tangent track");
      for (double& component : output->points[i].rayDirection)
        component /= norm;
      if (i > 0) {
        double continuity = 0.0;
        for (int d = 0; d < 3; ++d)
          continuity += output->points[i-1].rayDirection[d]*direction[d];
        if (continuity < 0.0)
          for (double& component : output->points[i].rayDirection)
            component = -component;
      }
    }
    return Status::Ok();
  }

  std::ifstream observer(request.observerGeometryFile);
  if (!observer)
    return Status::Error("cannot open observer_geometry_file: " +
                         request.observerGeometryFile);
  if (!std::getline(observer, line))
    return Status::Error("observer geometry CSV is empty");
  std::map<std::string, std::size_t> observerColumns;
  status = CsvColumns(line, &observerColumns);
  if (!status.success) return status;
  const std::vector<std::string> observerRequired = {
      "time_source_s", "observer_coordinate_from_tangent_m",
      "ray_ux_table", "ray_uy_table", "ray_uz_table"};
  std::map<std::string, std::size_t> oi;
  for (const std::string& name : observerRequired) {
    std::size_t column = 0;
    status = RequireCsvColumn(observerColumns, name, &column);
    if (!status.success) return status;
    oi[name] = column;
  }
  std::size_t row = 0;
  while (std::getline(observer, line)) {
    if (line.empty()) continue;
    if (row >= output->points.size())
      return Status::Error(request.id + " observer geometry has extra rows");
    const std::vector<std::string> fields = SplitCsvQuoted(line);
    if (fields.size() != observerColumns.size())
      return Status::Error(request.id + " observer geometry row is malformed");
    double time = 0.0;
    ObservationPoint& point = output->points[row];
    if (!ParseCsvFinite(fields[oi["time_source_s"]], &time) ||
        std::abs(time-rayTimes[row]) > 1.0e-5 ||
        !ParseCsvFinite(fields[oi["observer_coordinate_from_tangent_m"]],
                        &point.observerCoordinateM) ||
        !ParseCsvFinite(fields[oi["ray_ux_table"]], &point.rayDirection[0]) ||
        !ParseCsvFinite(fields[oi["ray_uy_table"]], &point.rayDirection[1]) ||
        !ParseCsvFinite(fields[oi["ray_uz_table"]], &point.rayDirection[2]))
      return Status::Error(request.id +
          " observer geometry time/direction is invalid or misaligned");
    double norm = 0.0;
    for (double component : point.rayDirection) norm += component*component;
    norm = std::sqrt(norm);
    if (!(norm > 0.0))
      return Status::Error(request.id + " observer ray direction is zero");
    for (double& component : point.rayDirection) component /= norm;
    if (point.observerCoordinateM < 0.0) {
      double radius2 = 0.0, projection = 0.0;
      for (int d = 0; d < 3; ++d) {
        radius2 += point.positionM[d]*point.positionM[d];
        projection += point.rayDirection[d]*point.positionM[d];
      }
      for (int d = 0; d < 3; ++d)
        point.rayDirection[d] -=
            projection/radius2*point.positionM[d];
      norm = 0.0;
      for (double component : point.rayDirection) norm += component*component;
      norm = std::sqrt(norm);
      if (!(norm > 0.0))
        return Status::Error(request.id + " projected observer ray is zero");
      for (double& component : point.rayDirection) component /= norm;
    } else {
      double radial = 0.0;
      for (int d = 0; d < 3; ++d)
        radial += point.positionM[d]*point.rayDirection[d];
      if (radial < 0.0)
        for (double& component : point.rayDirection) component = -component;
    }
    ++row;
  }
  if (row != output->points.size())
    return Status::Error(request.id +
        " observer geometry row count does not match UVIS data");
  return Status::Ok();
}

Status LoadFlybyTrajectory(const Configuration& configuration,
                           std::vector<FlybyPoint>* output) {
  output->clear();
  std::ifstream input(configuration.flybyTrajectoryFile);
  if (!input)
    return Status::Error("cannot open flyby_trajectory_file: " +
                         configuration.flybyTrajectoryFile);
  std::string line;
  if (!std::getline(input, line))
    return Status::Error("flyby trajectory is empty: " +
                         configuration.flybyTrajectoryFile);
  const std::string expectedHeader =
      "flyby,time_from_reference_s,utc,spacecraft_x_m,spacecraft_y_m,"
      "spacecraft_z_m,subpoint_x_m,subpoint_y_m,subpoint_z_m,altitude_m,"
      "west_lon_planetocentric_deg,lat_planetocentric_deg,"
      "west_lon_planetographic_deg,lat_planetographic_deg";
  if (line != expectedHeader)
    return Status::Error(
        "flyby trajectory header does not match the frozen SPICE schema");
  std::size_t row = 1;
  double previousTime = -std::numeric_limits<double>::infinity();
  while (std::getline(input, line)) {
    ++row;
    if (line.empty()) continue;
    const std::vector<std::string> fields = SplitCsv(line);
    if (fields.size() != 14)
      return Status::Error("flyby trajectory row " + std::to_string(row) +
                           " does not have 14 CSV fields");
    if (fields[0] != configuration.flybyId)
      return Status::Error("flyby trajectory row " + std::to_string(row) +
                           " has flyby id '" + fields[0] + "', expected '" +
                           configuration.flybyId + "'");
    FlybyPoint point;
    if (!ParseCsvFinite(fields[1], &point.timeFromReferenceS) ||
        !ParseCsvFinite(fields[3], &point.positionM[0]) ||
        !ParseCsvFinite(fields[4], &point.positionM[1]) ||
        !ParseCsvFinite(fields[5], &point.positionM[2]) ||
        !ParseCsvFinite(fields[9], &point.altitudeM) ||
        point.altitudeM < 0.0)
      return Status::Error("flyby trajectory row " + std::to_string(row) +
                           " has invalid SI time, position, or altitude");
    if (!(point.timeFromReferenceS > previousTime))
      return Status::Error("flyby trajectory times must be strictly increasing");
    previousTime = point.timeFromReferenceS;
    output->push_back(point);
  }
  if (output->size() < 2)
    return Status::Error("flyby trajectory requires at least two points");
  return Status::Ok();
}

Vec3 Cross(const Vec3& a, const Vec3& b) {
  return Vec3{a.y*b.z-a.z*b.y, a.z*b.x-a.x*b.z, a.x*b.y-a.y*b.x};
}

int Species(const std::string& symbol) {
  const auto found = gSpeciesIndex.find(symbol);
  if (found == gSpeciesIndex.end()) Stop("compiled species table has no " + symbol);
  return found->second;
}

bool IsSimulated(const std::string& symbol) {
  return std::find(gConfiguration.simulatedSpecies.begin(),
                   gConfiguration.simulatedSpecies.end(), symbol) !=
         gConfiguration.simulatedSpecies.end();
}

bool IsChemistrySpecies(const std::string& symbol) {
  return std::find(gConfiguration.chemistrySpecies.begin(),
                   gConfiguration.chemistrySpecies.end(), symbol) !=
         gConfiguration.chemistrySpecies.end();
}

std::string ParticleSymbol(int species) {
  if (species < 0 || species >= PIC::nTotalSpecies)
    Stop("particle has an invalid compiled species index");
  return PIC::MolecularData::GetChemSymbol(species);
}

double MinimumAllocatedCellSizeLocal() {
  double minimum = std::numeric_limits<double>::infinity();
  PIC::DomainBlockDecomposition::UpdateBlockTable();
  for (unsigned int index = 0;
       index < PIC::DomainBlockDecomposition::nLocalBlocks; ++index) {
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node =
        PIC::DomainBlockDecomposition::BlockTable[index];
    if (node == nullptr || node->block == nullptr) continue;
    const double widths[3] = {
        (node->xmax[0]-node->xmin[0])/_BLOCK_CELLS_X_,
        (node->xmax[1]-node->xmin[1])/_BLOCK_CELLS_Y_,
        (node->xmax[2]-node->xmin[2])/_BLOCK_CELLS_Z_};
    for (double width : widths)
      if (width > 0.0 && std::isfinite(width)) minimum = std::min(minimum, width);
  }
  return minimum;
}

double SourcePhysicalRate(const Source& source, const std::string& species) {
  const auto found = source.strength.find(species);
  if (found == source.strength.end()) return 0.0;
  if (gDistribution.mode == SourceMode::Absolute) return found->second;
  return found->second * gPlan.physicalSourceRatePerS.at(species);
}

const Source& SelectSource(const std::string& species, double target) {
  double cumulative = 0.0;
  for (const Source& source : gDistribution.sources) {
    cumulative += SourcePhysicalRate(source, species);
    if (target < cumulative) return source;
  }
  // target is strictly below the total, so only accumulated roundoff can
  // reach this fallback. Return the last positive source, never a zero-rate
  // component.
  for (auto source = gDistribution.sources.rbegin();
       source != gDistribution.sources.rend(); ++source)
    if (SourcePhysicalRate(*source, species) > 0.0) return *source;
  Stop("source selection has no positive component for " + species);
}

void SampleSourceGeometry(const Source& source, DeterministicRandom* random,
                          Vec3* position, Vec3* normal, Vec3* axis) {
  if (source.geometry == SourceGeometry::PointJet) {
    *position = PointJetPositionM(source, _RADIUS_(_TARGET_));
    *normal = PointJetNormal(source);
    *axis = PointJetAxis(source);
    return;
  }
  std::vector<double> weights;
  const Status status = StripeNodeWeights(source, &weights);
  if (!status.success) Stop(status.message);
  const double selection = random->Uniform();
  double cumulative = 0.0;
  std::size_t node = weights.size() - 1;
  for (std::size_t index = 0; index < weights.size(); ++index) {
    cumulative += weights[index];
    if (selection < cumulative) {
      node = index;
      break;
    }
  }
  *position = source.verticesM[node];
  *normal = Normalize(*position);
  *axis = *normal;
}

Vec3 SampleTenishevVelocity(const Source& source, const Vec3& normal,
                            const Vec3& axis, double massKg,
                            DeterministicRandom* random) {
  // A model particle represents an emitted molecule, so it samples the
  // production flux, not the local phase-space density.  Equations (3)-(7)
  // of the standalone law give the joint event density
  //   t^3 exp[-t^2+2*s*c*t-s^2] d(t) d(Omega),
  // where t=v/sqrt(2 k_B T/m), s=u/v_th, and c=v_hat.axis.  The Tenishev
  // angular factor is unity over v_hat.normal>0.  Directions are therefore
  // proposed uniformly on that hemisphere and accepted by K3(c,s)/K3(1,s);
  // the conditional speed is inverted from the analytic K3 CDF.
  const double thermalSpeed =
      std::sqrt(2.0*kBoltzmann*source.temperatureK/massKg);
  const double speedRatio = source.bulkSpeedMPerS/thermalSpeed;
  double maximumMoment = 0.0;
  Status status = Tenishev2010FluxMoment3(1.0, speedRatio, &maximumMoment);
  if (!status.success) Stop(source.id + ": " + status.message);

  const Vec3 reference = std::abs(normal.z) < 0.9
      ? Vec3{0.0,0.0,1.0} : Vec3{1.0,0.0,0.0};
  const Vec3 tangent0 = Normalize(Cross(reference, normal));
  const Vec3 tangent1 = Cross(normal, tangent0);
  for (int attempt = 0; attempt < 100000; ++attempt) {
    const double mu = random->Uniform();
    const double phi = 2.0*Pi*random->Uniform();
    const double transverse = std::sqrt(std::max(0.0, 1.0-mu*mu));
    const Vec3 direction{
        mu*normal.x + transverse*(std::cos(phi)*tangent0.x +
                                  std::sin(phi)*tangent1.x),
        mu*normal.y + transverse*(std::cos(phi)*tangent0.y +
                                  std::sin(phi)*tangent1.y),
        mu*normal.z + transverse*(std::cos(phi)*tangent0.z +
                                  std::sin(phi)*tangent1.z)};
    const double cosine = Dot(direction, axis);
    double moment = 0.0;
    status = Tenishev2010FluxMoment3(cosine, speedRatio, &moment);
    if (!status.success) Stop(source.id + ": " + status.message);
    if (random->Uniform()*maximumMoment > moment) continue;
    double speedOverThermal = 0.0;
    status = Tenishev2010SpeedQuantile(
        cosine, speedRatio, random->Uniform(), &speedOverThermal);
    if (!status.success) Stop(source.id + ": " + status.message);
    const double speed = thermalSpeed*speedOverThermal;
    return Vec3{speed*direction.x, speed*direction.y, speed*direction.z};
  }
  Stop("Tenishev flux sampler did not converge for " + source.id);
}

void LinkParticle(long int particle, long int& firstParticleCell) {
  PIC::ParticleBuffer::SetNext(firstParticleCell, particle);
  PIC::ParticleBuffer::SetPrev(-1, particle);
  if (firstParticleCell != -1)
    PIC::ParticleBuffer::SetPrev(particle, firstParticleCell);
  firstParticleCell = particle;
}

std::vector<const Reaction*> Channels(const std::string& parent) {
  std::vector<const Reaction*> result;
  for (const std::string& id : gConfiguration.selectedReactions) {
    const Reaction* reaction = FindReaction(id);
    if (reaction != nullptr && reaction->parent == parent)
      result.push_back(reaction);
  }
  return result;
}

double LocalReactionRate(const Reaction& reaction, const double* x,
                         const double* neutralVelocity,
                         cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node,
                         double* plasmaVelocity) {
  if (reaction.kind == ReactionKind::Photo) {
    double rate = 0.0;
    const Status status = EvaluateReactionRate(
        reaction, gConfiguration.solarActivity, 0.0, 0.0, 0.0, &rate);
    if (!status.success) Stop(reaction.id + ": " + status.message);
    return rate;
  }

  PIC::InterpolationRoutines::CellCentered::Linear::InitStencil(
      const_cast<double*>(x), node);
  const double density =
      PIC::CPLR::GetBackgroundPlasmaNumberDensity();
  PIC::CPLR::GetBackgroundPlasmaVelocity(plasmaVelocity);
  double relative2 = 0.0;
  for (int dimension = 0; dimension < 3; ++dimension) {
    const double delta = plasmaVelocity[dimension] - neutralVelocity[dimension];
    relative2 += delta*delta;
  }
  const auto fraction =
      gConfiguration.plasmaIonFraction.find(reaction.plasmaReactant);
  if (fraction == gConfiguration.plasmaIonFraction.end())
    Stop("no ion fraction for " + reaction.plasmaReactant);
  double rate = 0.0;
  const Status status = EvaluateReactionRate(
      reaction, gConfiguration.solarActivity, density, fraction->second,
      std::sqrt(relative2), &rate);
  if (!status.success) Stop(reaction.id + ": " + status.message);
  return rate;
}

void FlybySamplingNoOp() {}

void OutputFlybyDensity(int dataOutputFileNumber) {
  if (!gConfiguration.flybyDensityOutput) return;
  if (PIC::LastSampleLength != gConfiguration.sampleLengthSteps)
    Stop("flyby density output received AMPS sample length " +
         std::to_string(PIC::LastSampleLength) + ", expected " +
         std::to_string(gConfiguration.sampleLengthSteps));

  const std::size_t pointCount = gFlybyTrajectory.size();
  const std::size_t speciesCount = gConfiguration.simulatedSpecies.size();
  std::vector<double> localDensity(pointCount*speciesCount, 0.0);
  std::vector<double> localDomain(pointCount, 0.0);
  std::vector<double> localCellMaxEdge(pointCount, 0.0);

  // Only the rank owning a leaf evaluates that point. The completed AMPS
  // number-density datum is the physical-particle weight divided by cell
  // volume and by LastSampleLength [m^-3]. This build selects the native
  // cell-centered constant stencil, so the value is the containing cell's
  // sampling-window mean. Recording the cell edge makes spatial resolution an
  // auditable part of every observation row rather than an inferred property.
  for (std::size_t pointIndex = 0; pointIndex < pointCount; ++pointIndex) {
    double position[3] = {
        gFlybyTrajectory[pointIndex].positionM[0],
        gFlybyTrajectory[pointIndex].positionM[1],
        gFlybyTrajectory[pointIndex].positionM[2]};
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node =
        PIC::Mesh::mesh->findTreeNode(position);
    if (node == nullptr || !node->IsUsedInCalculationFlag ||
        node->Thread != PIC::ThisThread || node->block == nullptr) continue;
    localDomain[pointIndex] = 1.0;
    localCellMaxEdge[pointIndex] = std::max(
        (node->xmax[0]-node->xmin[0])/_BLOCK_CELLS_X_,
        std::max((node->xmax[1]-node->xmin[1])/_BLOCK_CELLS_Y_,
                 (node->xmax[2]-node->xmin[2])/_BLOCK_CELLS_Z_));
    double coefficients[PIC::InterpolationRoutines::nMaxStencilLength];
    PIC::Mesh::cDataCenterNode* stencil[
        PIC::InterpolationRoutines::nMaxStencilLength];
    const int stencilLength = PIC::Mesh::GetCenterNodesInterpolationCoefficients(
        position, coefficients, stencil, node,
        PIC::InterpolationRoutines::nMaxStencilLength);
    if (stencilLength <= 0)
      Stop("flyby density interpolation produced an empty stencil");
    for (std::size_t speciesSlot = 0; speciesSlot < speciesCount;
         ++speciesSlot) {
      const int species = Species(gConfiguration.simulatedSpecies[speciesSlot]);
      double density = 0.0;
      for (int cell = 0; cell < stencilLength; ++cell)
        density += coefficients[cell]*stencil[cell]->GetNumberDensity(species);
      if (!std::isfinite(density) || density < 0.0)
        Stop("flyby density interpolation produced a non-finite or negative "
             "number density");
      localDensity[pointIndex*speciesCount+speciesSlot] = density;
    }
  }

  std::vector<double> globalDensity(pointCount*speciesCount, 0.0);
  std::vector<double> globalDomain(pointCount, 0.0);
  std::vector<double> globalCellMaxEdge(pointCount, 0.0);
  MPI_Reduce(localDensity.data(), globalDensity.data(),
             static_cast<int>(localDensity.size()), MPI_DOUBLE, MPI_SUM, 0,
             MPI_GLOBAL_COMMUNICATOR);
  MPI_Reduce(localDomain.data(), globalDomain.data(),
             static_cast<int>(localDomain.size()), MPI_DOUBLE, MPI_SUM, 0,
             MPI_GLOBAL_COMMUNICATOR);
  MPI_Reduce(localCellMaxEdge.data(), globalCellMaxEdge.data(),
             static_cast<int>(localCellMaxEdge.size()), MPI_DOUBLE, MPI_SUM, 0,
             MPI_GLOBAL_COMMUNICATOR);
  if (PIC::ThisThread != 0) return;

  std::size_t inDomainCount = 0;
  for (double owners : globalDomain) {
    if (owners == 1.0) ++inDomainCount;
    else if (owners != 0.0)
      Stop("flyby trajectory point has more than one AMR owner");
  }
  const std::string path = std::string(PIC::OutputDataFileDirectory) +
      "/flyby_density." + gConfiguration.flybyId + ".dat";
  if (dataOutputFileNumber > 0 && !fs::is_regular_file(path))
    Stop("flyby Tecplot file is missing before append: " + path);
  std::ofstream output(path, dataOutputFileNumber == 0
      ? (std::ios::out | std::ios::trunc)
      : (std::ios::out | std::ios::app));
  if (!output) Stop("cannot open flyby Tecplot output: " + path);
  output << std::scientific << std::setprecision(17);
  if (dataOutputFileNumber == 0) {
    output << "TITLE=\"AMPS local neutral density along "
           << gConfiguration.flybyId << "\"\n"
           << "# IAU_ENCELADUS spacecraft positions from frozen SPICE CSV: "
           << gConfiguration.flybyTrajectoryFile << "\n"
           << "# Density is the completed AMPS sampling-window mean. "
              "Rows outside the AMR domain have in_domain=0 and zero "
              "placeholder densities; exclude them from physical analysis.\n"
           << "VARIABLES=\"time_from_reference_s\",\"spacecraft_x_m\","
              "\"spacecraft_y_m\",\"spacecraft_z_m\",\"altitude_m\","
              "\"in_domain\",\"cell_max_edge_m\"";
    for (const std::string& symbol : gConfiguration.simulatedSpecies)
      output << ",\"n_" << symbol << "_m-3\"";
    output << ",\"n_total_primary_m-3\"\n";
  }
  const double solutionTimeS =
      (dataOutputFileNumber+1.0)*gConfiguration.sampleLengthSteps*
      gPlan.globalTimeStepS;
  output << "ZONE T=\"output=" << dataOutputFileNumber
         << ", iterations="
         << (dataOutputFileNumber+1)*gConfiguration.sampleLengthSteps
         << ", sample_steps=" << PIC::LastSampleLength
         << "\" I=" << pointCount
         << " DATAPACKING=POINT STRANDID=1 SOLUTIONTIME="
         << solutionTimeS << "\n";
  for (std::size_t pointIndex = 0; pointIndex < pointCount; ++pointIndex) {
    const FlybyPoint& point = gFlybyTrajectory[pointIndex];
    output << point.timeFromReferenceS << ' ' << point.positionM[0] << ' '
           << point.positionM[1] << ' ' << point.positionM[2] << ' '
           << point.altitudeM << ' ' << globalDomain[pointIndex] << ' '
           << globalCellMaxEdge[pointIndex];
    double total = 0.0;
    for (std::size_t speciesSlot = 0; speciesSlot < speciesCount;
         ++speciesSlot) {
      const double density =
          globalDensity[pointIndex*speciesCount+speciesSlot];
      output << ' ' << density;
      total += density;
    }
    output << ' ' << total << '\n';
  }
  output.close();
  if (!output) Stop("failed while writing flyby Tecplot output: " + path);
  std::cout << "[srcEnceladus] flyby_density_output=" << dataOutputFileNumber
            << " flyby=" << gConfiguration.flybyId
            << " trajectory_points=" << pointCount
            << " in_domain_points=" << inDomainCount
            << " sample_length_steps=" << PIC::LastSampleLength
            << " file=" << path << '\n';
}

void ObservationColumnIntegrand(
    double* result, int resultLength, double* position,
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node) {
  const std::size_t speciesCount = gConfiguration.simulatedSpecies.size();
  if (resultLength != static_cast<int>(speciesCount))
    Stop("observation column integrand received an inconsistent vector size");
  double coefficients[PIC::InterpolationRoutines::nMaxStencilLength];
  PIC::Mesh::cDataCenterNode* stencil[
      PIC::InterpolationRoutines::nMaxStencilLength];
  const int stencilLength = PIC::Mesh::GetCenterNodesInterpolationCoefficients(
      position, coefficients, stencil, node,
      PIC::InterpolationRoutines::nMaxStencilLength);
  if (stencilLength <= 0)
    Stop("UVIS column integration produced an empty density stencil");
  for (std::size_t slot = 0; slot < speciesCount; ++slot) {
    const int species = Species(gConfiguration.simulatedSpecies[slot]);
    double density = 0.0;
    for (int cell = 0; cell < stencilLength; ++cell)
      density += coefficients[cell]*stencil[cell]->GetNumberDensity(species);
    if (!std::isfinite(density) || density < 0.0)
      Stop("UVIS column integration produced a non-finite or negative density");
    result[slot] = density;
  }
}

double RayPathLengthInDomain(double* origin, double* direction) {
  double length = 0.0, start[3], finish[3];
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* startNode = nullptr;
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* finishNode = nullptr;
  if (!PIC::ColumnIntegration::FindIntegrationLimits(
          origin, direction, length, start, startNode, finish, finishNode))
    return 0.0;
  return length;
}

void WriteInmsObservation(const LoadedObservation& observation,
                          int dataOutputFileNumber) {
  const std::size_t pointCount = observation.points.size();
  const std::size_t speciesCount = gConfiguration.simulatedSpecies.size();
  std::vector<double> localDensity(pointCount*speciesCount, 0.0);
  std::vector<double> localDomain(pointCount, 0.0);
  std::vector<double> localCellMaxEdge(pointCount, 0.0);
  for (std::size_t pointIndex = 0; pointIndex < pointCount; ++pointIndex) {
    double position[3];
    for (int d = 0; d < 3; ++d)
      position[d] = observation.points[pointIndex].positionM[d];
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node =
        PIC::Mesh::mesh->findTreeNode(position);
    if (node == nullptr || !node->IsUsedInCalculationFlag ||
        node->Thread != PIC::ThisThread || node->block == nullptr) continue;
    localDomain[pointIndex] = 1.0;
    localCellMaxEdge[pointIndex] = std::max(
        (node->xmax[0]-node->xmin[0])/_BLOCK_CELLS_X_,
        std::max((node->xmax[1]-node->xmin[1])/_BLOCK_CELLS_Y_,
                 (node->xmax[2]-node->xmin[2])/_BLOCK_CELLS_Z_));
    double coefficients[PIC::InterpolationRoutines::nMaxStencilLength];
    PIC::Mesh::cDataCenterNode* stencil[
        PIC::InterpolationRoutines::nMaxStencilLength];
    const int stencilLength = PIC::Mesh::GetCenterNodesInterpolationCoefficients(
        position, coefficients, stencil, node,
        PIC::InterpolationRoutines::nMaxStencilLength);
    if (stencilLength <= 0)
      Stop(observation.request.id +
           " INMS interpolation produced an empty stencil");
    for (std::size_t slot = 0; slot < speciesCount; ++slot) {
      const int species = Species(gConfiguration.simulatedSpecies[slot]);
      double density = 0.0;
      for (int cell = 0; cell < stencilLength; ++cell)
        density += coefficients[cell]*stencil[cell]->GetNumberDensity(species);
      if (!std::isfinite(density) || density < 0.0)
        Stop(observation.request.id +
             " INMS interpolation produced an invalid density");
      localDensity[pointIndex*speciesCount+slot] = density;
    }
  }
  std::vector<double> density(pointCount*speciesCount, 0.0);
  std::vector<double> domain(pointCount, 0.0);
  std::vector<double> cellMaxEdge(pointCount, 0.0);
  MPI_Reduce(localDensity.data(), density.data(),
             static_cast<int>(density.size()), MPI_DOUBLE, MPI_SUM, 0,
             MPI_GLOBAL_COMMUNICATOR);
  MPI_Reduce(localDomain.data(), domain.data(), static_cast<int>(domain.size()),
             MPI_DOUBLE, MPI_SUM, 0, MPI_GLOBAL_COMMUNICATOR);
  MPI_Reduce(localCellMaxEdge.data(), cellMaxEdge.data(),
             static_cast<int>(cellMaxEdge.size()), MPI_DOUBLE, MPI_SUM, 0,
             MPI_GLOBAL_COMMUNICATOR);
  if (PIC::ThisThread != 0) return;

  const std::string path = std::string(PIC::OutputDataFileDirectory) +
      "/observation." + observation.request.id + ".dat";
  if (dataOutputFileNumber > 0 && !fs::is_regular_file(path))
    Stop("INMS observation file is missing before append: " + path);
  std::ofstream output(path, dataOutputFileNumber == 0
      ? (std::ios::out | std::ios::trunc)
      : (std::ios::out | std::ios::app));
  if (!output) Stop("cannot open INMS observation output: " + path);
  output << std::scientific << std::setprecision(17);
  if (dataOutputFileNumber == 0) {
    output << "TITLE=\"AMPS INMS comparison " << observation.request.eventId
           << "\"\n# request_id=" << observation.request.id
           << " event_id=" << observation.request.eventId
           << " data_file=" << observation.request.dataFile << "\n"
           << "# model-only rows have no observed density columns; this is "
              "intentional when no response-qualified H2O/CO2 product exists.\n"
           << (observation.request.contextObservable.empty()
                   ? std::string()
                   : std::string("# context_observable=") +
                         observation.request.contextObservable +
                         " is native instrument data; it is not converted or "
                         "scored against model density.\n")
           << "VARIABLES=\"time_s\",\"spacecraft_x_m\","
              "\"spacecraft_y_m\",\"spacecraft_z_m\",\"altitude_m\","
              "\"in_domain\",\"cell_max_edge_m\"";
    for (const std::string& symbol : gConfiguration.simulatedSpecies)
      output << ",\"model_n_" << symbol << "_m-3\"";
    if (!observation.request.comparisonSpecies.empty())
      output << ",\"observed_n_" << observation.request.comparisonSpecies
             << "_m-3\",\"observed_sigma_"
             << observation.request.comparisonSpecies << "_m-3\"";
    else if (!observation.request.contextObservable.empty())
      output << ",\"observed_context_value\",\"observed_context_sigma\"";
    output << '\n';
  }
  const double solutionTimeS =
      (dataOutputFileNumber+1.0)*gConfiguration.sampleLengthSteps*
      gPlan.globalTimeStepS;
  output << "ZONE T=\"output=" << dataOutputFileNumber
         << ", iterations="
         << (dataOutputFileNumber+1)*gConfiguration.sampleLengthSteps
         << ", sample_steps=" << PIC::LastSampleLength
         << "\" I=" << pointCount
         << " DATAPACKING=POINT STRANDID=1 SOLUTIONTIME="
         << solutionTimeS << '\n';
  std::size_t covered = 0;
  for (std::size_t pointIndex = 0; pointIndex < pointCount; ++pointIndex) {
    const ObservationPoint& point = observation.points[pointIndex];
    const double radius = std::sqrt(
        point.positionM[0]*point.positionM[0] +
        point.positionM[1]*point.positionM[1] +
        point.positionM[2]*point.positionM[2]);
    if (domain[pointIndex] == 1.0) ++covered;
    else if (domain[pointIndex] != 0.0)
      Stop(observation.request.id + " point has more than one AMR owner");
    output << point.timeS << ' ' << point.positionM[0] << ' '
           << point.positionM[1] << ' ' << point.positionM[2] << ' '
           << radius-_RADIUS_(_TARGET_) << ' ' << domain[pointIndex] << ' '
           << cellMaxEdge[pointIndex];
    for (std::size_t slot = 0; slot < speciesCount; ++slot)
      output << ' ' << density[pointIndex*speciesCount+slot];
    if (!observation.request.comparisonSpecies.empty())
      output << ' ' << point.observedValue << ' ' << point.observedSigma;
    else if (!observation.request.contextObservable.empty())
      output << ' ' << point.observedValue << ' ' << point.observedSigma;
    output << '\n';
  }
  output.close();
  if (!output) Stop("failed while writing INMS observation output: " + path);
  std::cout << "[srcEnceladus] observation_output=" << dataOutputFileNumber
            << " instrument=INMS request=" << observation.request.id
            << " event=" << observation.request.eventId
            << " rows=" << pointCount << " in_domain=" << covered
            << " file=" << path << '\n';
}

void WriteUvisObservation(const LoadedObservation& observation,
                          int dataOutputFileNumber) {
  const std::size_t pointCount = observation.points.size();
  const std::size_t speciesCount = gConfiguration.simulatedSpecies.size();
  std::vector<double> columns(pointCount*speciesCount, 0.0);
  std::vector<double> pathLengths(pointCount, 0.0);
  for (std::size_t pointIndex = 0; pointIndex < pointCount; ++pointIndex) {
    const ObservationPoint& point = observation.points[pointIndex];
    double direction[3], origin[3];
    for (int d = 0; d < 3; ++d) {
      direction[d] = point.rayDirection[d];
      origin[d] = point.positionM[d];
    }
    std::vector<double> result(speciesCount, 0.0);
    if (observation.request.uvisRayGeometry ==
        UvisRayGeometry::InfiniteChordTangentTrackNormal) {
      std::vector<double> reverse(speciesCount, 0.0);
      pathLengths[pointIndex] = RayPathLengthInDomain(origin, direction);
      PIC::ColumnIntegration::GetCoulumnIntegral(
          result.data(), static_cast<int>(speciesCount), origin, direction,
          ObservationColumnIntegrand);
      for (double& component : direction) component = -component;
      pathLengths[pointIndex] += RayPathLengthInDomain(origin, direction);
      PIC::ColumnIntegration::GetCoulumnIntegral(
          reverse.data(), static_cast<int>(speciesCount), origin, direction,
          ObservationColumnIntegrand);
      for (std::size_t slot = 0; slot < speciesCount; ++slot)
        result[slot] += reverse[slot];
    } else {
      if (point.observerCoordinateM < 0.0)
        for (int d = 0; d < 3; ++d)
          origin[d] += point.observerCoordinateM*direction[d];
      pathLengths[pointIndex] = RayPathLengthInDomain(origin, direction);
      PIC::ColumnIntegration::GetCoulumnIntegral(
          result.data(), static_cast<int>(speciesCount), origin, direction,
          ObservationColumnIntegrand);
    }
    for (std::size_t slot = 0; slot < speciesCount; ++slot) {
      if (!std::isfinite(result[slot]) || result[slot] < 0.0)
        Stop(observation.request.id + " produced an invalid column density");
      columns[pointIndex*speciesCount+slot] = result[slot];
    }
  }
  if (PIC::ThisThread != 0) return;
  const std::string path = std::string(PIC::OutputDataFileDirectory) +
      "/observation." + observation.request.id + ".dat";
  if (dataOutputFileNumber > 0 && !fs::is_regular_file(path))
    Stop("UVIS observation file is missing before append: " + path);
  std::ofstream output(path, dataOutputFileNumber == 0
      ? (std::ios::out | std::ios::trunc)
      : (std::ios::out | std::ios::app));
  if (!output) Stop("cannot open UVIS observation output: " + path);
  output << std::scientific << std::setprecision(17);
  if (dataOutputFileNumber == 0) {
    output << "TITLE=\"AMPS UVIS comparison " << observation.request.eventId
           << "\"\n# request_id=" << observation.request.id
           << " event_id=" << observation.request.eventId
           << " ray_geometry=" << Name(observation.request.uvisRayGeometry)
           << " data_file=" << observation.request.dataFile << "\n"
           << "# published_h2o_column_m-2="
           << observation.request.publishedH2OColumnM2
           << " published_h2o_column_sigma_m-2="
           << observation.request.publishedH2OColumnSigmaM2 << "\n"
           << "# Optical depth is a shape proxy. The comparison plot fits "
              "one common column-to-optical-depth coefficient and one signed "
              "baseline; it does not convert each row into an absolute column.\n"
           << "VARIABLES=\"time_from_min_altitude_s\",\"tangent_x_m\","
              "\"tangent_y_m\",\"tangent_z_m\","
              "\"integrated_path_in_domain_m\"";
    for (const std::string& symbol : gConfiguration.simulatedSpecies)
      output << ",\"model_column_" << symbol << "_m-2\"";
    output << ",\"observed_optical_depth\","
              "\"observed_sigma_optical_depth\"\n";
  }
  const double solutionTimeS =
      (dataOutputFileNumber+1.0)*gConfiguration.sampleLengthSteps*
      gPlan.globalTimeStepS;
  output << "ZONE T=\"output=" << dataOutputFileNumber
         << ", iterations="
         << (dataOutputFileNumber+1)*gConfiguration.sampleLengthSteps
         << ", sample_steps=" << PIC::LastSampleLength
         << "\" I=" << pointCount
         << " DATAPACKING=POINT STRANDID=1 SOLUTIONTIME="
         << solutionTimeS << '\n';
  for (std::size_t pointIndex = 0; pointIndex < pointCount; ++pointIndex) {
    const ObservationPoint& point = observation.points[pointIndex];
    output << point.timeS << ' ' << point.positionM[0] << ' '
           << point.positionM[1] << ' ' << point.positionM[2] << ' '
           << pathLengths[pointIndex];
    for (std::size_t slot = 0; slot < speciesCount; ++slot)
      output << ' ' << columns[pointIndex*speciesCount+slot];
    output << ' ' << point.observedValue << ' ' << point.observedSigma << '\n';
  }
  output.close();
  if (!output) Stop("failed while writing UVIS observation output: " + path);
  std::cout << "[srcEnceladus] observation_output=" << dataOutputFileNumber
            << " instrument=UVIS request=" << observation.request.id
            << " event=" << observation.request.eventId
            << " rays=" << pointCount << " file=" << path << '\n';
}

void OutputRequestedObservations(int dataOutputFileNumber) {
  if (gObservations.empty()) return;
  if (PIC::LastSampleLength != gConfiguration.sampleLengthSteps)
    Stop("observation output received an unexpected AMPS sample length");
  for (const LoadedObservation& observation : gObservations) {
    if (observation.request.instrument == ObservationInstrument::Inms)
      WriteInmsObservation(observation, dataOutputFileNumber);
    else
      WriteUvisObservation(observation, dataOutputFileNumber);
  }
}

void OutputAllObservationProducts(int dataOutputFileNumber) {
  OutputFlybyDensity(dataOutputFileNumber);
  OutputRequestedObservations(dataOutputFileNumber);
  if (!gConfiguration.volumeFieldOutput) {
    // PIC switches completed/collecting sampling buffers before invoking any
    // output routine. Generic volume output normally clears the newly active
    // collecting buffer inside its per-species loop. With every volume file
    // disabled, clear it here only after all observation callbacks have read
    // the immutable completed buffer, preserving independent sampling windows.
    for (cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node =
             PIC::Mesh::mesh->ParallelNodesDistributionList[PIC::ThisThread];
         node != nullptr; node = node->nextNodeThisThread) {
      if (node->block == nullptr) continue;
      for (int k = 0; k < _BLOCK_CELLS_Z_; ++k) {
        for (int j = 0; j < _BLOCK_CELLS_Y_; ++j) {
          for (int i = 0; i < _BLOCK_CELLS_X_; ++i) {
            const int local = _getCenterNodeLocalNumber(i, j, k);
            PIC::Mesh::flushCollectingSamplingBuffer(
                node->block->GetCenterNode(local));
          }
        }
      }
    }
  }
}

}  // namespace

bool InstallFromInput(const std::string& inputFile, std::string* error) {
  if (gInstalled) {
    if (error != nullptr) *error = "Enceladus runtime configuration is already installed";
    return false;
  }
  Configuration configuration;
  Status status = ParseConfiguration(inputFile, &configuration);
  if (!status.success) {
    if (error != nullptr) *error = status.message;
    return false;
  }
  SourceDistribution distribution;
  status = ParseSourceDistribution(configuration.sourceDistributionFile,
                                   &distribution);
  if (!status.success) {
    if (error != nullptr) *error = status.message;
    return false;
  }
  status = ValidateConfiguration(configuration, distribution);
  if (!status.success) {
    if (error != nullptr) *error = status.message;
    return false;
  }
  std::vector<FlybyPoint> flybyTrajectory;
  if (configuration.flybyDensityOutput) {
    status = LoadFlybyTrajectory(configuration, &flybyTrajectory);
    if (!status.success) {
      if (error != nullptr) *error = status.message;
      return false;
    }
  }
  std::vector<LoadedObservation> observations;
  observations.reserve(configuration.observations.size());
  for (const ObservationRequest& request : configuration.observations) {
    LoadedObservation observation;
    status = LoadObservation(request, &observation);
    if (!status.success) {
      if (error != nullptr) *error = request.id + ": " + status.message;
      return false;
    }
    observations.push_back(std::move(observation));
  }
  if (std::abs(distribution.stripeReferenceRadiusM-_RADIUS_(_TARGET_)) >
      1.0e-6) {
    if (error != nullptr) {
      std::ostringstream message;
      message << "source receipt stripe_reference_radius_m="
              << std::setprecision(17) << distribution.stripeReferenceRadiusM
              << " does not match the configured spherical boundary "
              << _RADIUS_(_TARGET_) << " m";
      *error = message.str();
    }
    return false;
  }
  // Stripe coordinates are declared as physical surface source locations, not
  // generic points outside the body. Require agreement with the configured
  // spherical internal boundary to a 1 micrometre serialization tolerance.
  // This tolerance only absorbs decimal/binary roundoff; it is far too small
  // to hide a topographic or body-shape discrepancy.
  constexpr double surfaceSerializationToleranceM = 1.0e-6;
  for (const Source& source : distribution.sources) {
    if (source.geometry != SourceGeometry::PolylineStripe) continue;
    for (std::size_t vertex = 0; vertex < source.verticesM.size(); ++vertex) {
      const double radius = Norm(source.verticesM[vertex]);
      const double mismatch = std::abs(radius-_RADIUS_(_TARGET_));
      if (mismatch > surfaceSerializationToleranceM) {
        if (error != nullptr) {
          std::ostringstream message;
          message << source.id << " vertex " << vertex << " has radius "
                  << std::setprecision(17) << radius
                  << " m and does not lie on the configured spherical boundary "
                  << _RADIUS_(_TARGET_) << " m (absolute mismatch " << mismatch
                  << " m exceeds the " << surfaceSerializationToleranceM
                  << " m serialization tolerance) despite the declared "
                  << Name(distribution.stripeSurfaceMapping) << " mapping";
          *error = message.str();
        }
        return false;
      }
    }
  }
  std::error_code filesystemError;
  if (configuration.plasmaBackgroundMode != PlasmaBackgroundMode::None &&
      (!fs::is_directory(configuration.plasmaCacheDirectory, filesystemError) ||
       filesystemError)) {
    if (error != nullptr)
      *error = "plasma_cache_directory is not an accessible directory: " +
               configuration.plasmaCacheDirectory;
    return false;
  }
  if (configuration.plasmaBackgroundMode ==
      PlasmaBackgroundMode::CacheOrIces) {
    filesystemError.clear();
    if (!fs::is_directory(
            configuration.icesInstallationDirectory, filesystemError) ||
        filesystemError) {
      if (error != nullptr)
        *error = "ices_installation_directory is not an accessible directory: " +
                 configuration.icesInstallationDirectory;
      return false;
    }
  }
  if ((configuration.plasmaBackgroundMode != PlasmaBackgroundMode::None &&
       configuration.plasmaCacheDirectory.size() >= _MAX_STRING_LENGTH_PIC_) ||
      configuration.icesInstallationDirectory.size() >=
          _MAX_STRING_LENGTH_PIC_ ||
      configuration.icesCase.size() >= _MAX_STRING_LENGTH_PIC_) {
    if (error != nullptr) *error = "configured plasma/ICES path or case is too long";
    return false;
  }
  if (configuration.plasmaBackgroundMode != PlasmaBackgroundMode::None)
    std::snprintf(PIC::CPLR::DATAFILE::path, _MAX_STRING_LENGTH_PIC_, "%s",
                  configuration.plasmaCacheDirectory.c_str());
  if (configuration.plasmaBackgroundMode ==
      PlasmaBackgroundMode::CacheOrIces)
    PIC::CPLR::DATAFILE::ICES::SetLocationICES(
        configuration.icesInstallationDirectory.c_str());
  gConfiguration = std::move(configuration);
  gDistribution = std::move(distribution);
  gFlybyTrajectory = std::move(flybyTrajectory);
  gObservations = std::move(observations);
  gFlybyRefinementTrajectory.clear();
  if (gConfiguration.flybyMeshRefinement) {
    for (const FlybyPoint& point : gFlybyTrajectory)
      if (point.timeFromReferenceS >=
              gConfiguration.flybyRefinementTimeMinS &&
          point.timeFromReferenceS <=
              gConfiguration.flybyRefinementTimeMaxS)
        gFlybyRefinementTrajectory.push_back(point);
    if (gFlybyRefinementTrajectory.size() < 2) {
      if (error != nullptr)
        *error = "flyby refinement interval contains fewer than two trajectory points";
      gFlybyTrajectory.clear();
      gObservations.clear();
      gConfiguration = Configuration{};
      gDistribution = SourceDistribution{};
      return false;
    }
  }
  gInstalled = true;
  return true;
}

bool IsInstalled() { return gInstalled; }

double RequestedMeshResolution(const double* position,
                               double legacyResolutionM) {
  if (!gInstalled || !gConfiguration.flybyMeshRefinement ||
      gFlybyRefinementTrajectory.size() < 2)
    return legacyResolutionM;

  // The refinement region is the union of radius-R capsules around successive
  // straight segments of the frozen body-fixed trajectory. For a segment
  // a--b, the closest point is a+clamp(((x-a).(b-a))/|b-a|^2,0,1)(b-a).
  // Coordinates and R are metres in IAU_ENCELADUS. This geometrical AMR
  // selector changes no source, force, chemistry, or spacecraft ephemeris.
  const double radius2 = gConfiguration.flybyRefinementRadiusM *
                         gConfiguration.flybyRefinementRadiusM;
  bool inside = false;
  for (std::size_t index = 1;
       index < gFlybyRefinementTrajectory.size(); ++index) {
    const FlybyPoint& left = gFlybyRefinementTrajectory[index-1];
    const FlybyPoint& right = gFlybyRefinementTrajectory[index];
    double delta[3], relative[3];
    double length2 = 0.0, projection = 0.0;
    for (int dimension = 0; dimension < 3; ++dimension) {
      delta[dimension] = right.positionM[dimension]-left.positionM[dimension];
      relative[dimension] = position[dimension]-left.positionM[dimension];
      length2 += delta[dimension]*delta[dimension];
      projection += relative[dimension]*delta[dimension];
    }
    if (!(length2 > 0.0)) continue;
    const double fraction = std::max(0.0, std::min(1.0, projection/length2));
    double distance2 = 0.0;
    for (int dimension = 0; dimension < 3; ++dimension) {
      const double separation = relative[dimension]-fraction*delta[dimension];
      distance2 += separation*separation;
    }
    if (distance2 <= radius2) {
      inside = true;
      break;
    }
  }
  if (!inside) return legacyResolutionM;

  // meshAMR compares its requested resolution with the Euclidean diagonal of
  // a Cartesian cell, whereas the input parameter is the more user-facing
  // maximum edge width. Multiplication by sqrt(3) is therefore the exact
  // isotropic conversion. Dyadic AMR may produce a smaller actual edge.
  const double requestedDiagonal = std::sqrt(3.0) *
      gConfiguration.flybyRefinementTargetCellSizeM;
  return std::min(legacyResolutionM, requestedDiagonal);
}

const Configuration& InstalledConfiguration() {
  if (!gInstalled) Stop("runtime configuration was not installed");
  return gConfiguration;
}

const SourceDistribution& InstalledDistribution() {
  if (!gInstalled) Stop("source distribution was not installed");
  return gDistribution;
}

const NumericalPlan& InstalledNumericalPlan() {
  if (!gNumericsFinalized) Stop("numerical plan was not finalized");
  return gPlan;
}

void FinalizeNumericsAfterMeshAllocation() {
  if (!gInstalled) return;
  if (gNumericsFinalized) Stop("numerics were finalized more than once");
  gSpeciesIndex.clear();
  for (int species = 0; species < PIC::nTotalSpecies; ++species) {
    const std::string symbol = PIC::MolecularData::GetChemSymbol(species);
    if (!gSpeciesIndex.emplace(symbol, species).second)
      Stop("compiled species symbol occurs more than once: " + symbol);
  }
  gLocalInjected.assign(PIC::nTotalSpecies, 0);
  gLocalPhotolyticLoss.assign(PIC::nTotalSpecies, 0);
  gLocalSurfaceLoss.assign(PIC::nTotalSpecies, 0);
  gLocalDomainEscape.assign(PIC::nTotalSpecies, 0);
  for (const std::string& symbol : gConfiguration.simulatedSpecies)
    (void)Species(symbol);

  // The build contains a capability union of neutral, charged, electron, and
  // dust species, but this runtime stage activates only selected neutral
  // primaries. Generic AMPS defaults to writing a full mesh file for every
  // compiled species, including inactive zero fields. On the refined E3 mesh
  // each such file is hundreds of MB. When volume_field_output is true, retain
  // generic fields only for explicitly simulated species. When false, the
  // registered observation callback reads the same completed buffer and then
  // resets the next collecting buffer without serializing the AMR volume.
  // Storage, mover, chemistry, sampling, and observation values are unchanged.
  if (PIC::Sampling::SaveOutputDataFile == nullptr)
    Stop("generic AMPS output-selection table is not initialized");
  for (int species = 0; species < PIC::nTotalSpecies; ++species)
    PIC::Sampling::SaveOutputDataFile[species] =
        gConfiguration.volumeFieldOutput && IsSimulated(ParticleSymbol(species));

  // The analytic recurrence used by the exact hybrid-8 source sampler is
  // independently qualified for a=s*c>=-2. Check the entire allowed outward
  // hemisphere, not only whichever random direction is drawn first.
  for (const std::string& symbol : gConfiguration.injectedSpecies) {
    const int species = Species(symbol);
    const double mass = PIC::MolecularData::GetMass(species);
    for (const Source& source : gDistribution.sources) {
      const auto strength = source.strength.find(symbol);
      if (strength == source.strength.end() || strength->second == 0.0) continue;
      const Vec3 normal = source.geometry == SourceGeometry::PointJet
          ? PointJetNormal(source) : Normalize(source.verticesM.front());
      const Vec3 axis = source.geometry == SourceGeometry::PointJet
          ? PointJetAxis(source) : normal;
      const double gamma = std::max(0.0, std::min(1.0, Dot(normal, axis)));
      const double minimumCosine = -std::sqrt(std::max(0.0, 1.0-gamma*gamma));
      const double thermal =
          std::sqrt(2.0*kBoltzmann*source.temperatureK/mass);
      if (source.bulkSpeedMPerS/thermal*minimumCosine < -2.0)
        Stop(source.id + " exceeds the qualified Tenishev sampler range for " +
             symbol + " (minimum speed_ratio*direction_cosine < -2)");
    }
  }

  const double localMinimum = MinimumAllocatedCellSizeLocal();
  double globalMinimum = 0.0;
  MPI_Allreduce(&localMinimum, &globalMinimum, 1, MPI_DOUBLE, MPI_MIN,
                MPI_GLOBAL_COMMUNICATOR);
  const Status status = BuildNumericalPlan(
      gConfiguration, gDistribution, globalMinimum, &gPlan);
  if (!status.success) Stop(status.message);

  double fallbackWeight = std::numeric_limits<double>::infinity();
  for (const auto& item : gPlan.statisticalWeight)
    if (item.second > 0.0) fallbackWeight = std::min(fallbackWeight, item.second);
  if (!std::isfinite(fallbackWeight)) Stop("no positive particle weight was planned");
  for (int species = 0; species < PIC::nTotalSpecies; ++species) {
    const std::string symbol = ParticleSymbol(species);
    const auto weight = gPlan.statisticalWeight.find(symbol);
    PIC::ParticleWeightTimeStep::GlobalParticleWeight[species] =
        weight == gPlan.statisticalWeight.end() ? fallbackWeight : weight->second;
    PIC::ParticleWeightTimeStep::GlobalTimeStep[species] = gPlan.globalTimeStepS;
  }
  PIC::ParticleWeightTimeStep::GlobalTimeStepInitialized = true;
  PIC::ParticleWeightTimeStep::LocalTimeStep = GlobalTimeStep;
  PIC::ParticleWeightTimeStep::UserDefinedExtraSourceRate = PlannedSourceRate;
  PIC::Sampling::SampleTimeInterval =
      gConfiguration.sampleLengthSteps * gPlan.globalTimeStepS;
  PIC::RequiredSampleLength = gConfiguration.sampleLengthSteps;

  int splittingMode = PIC::ParticleSplitting::_disactivated;
  if (gConfiguration.populationControl == PopulationControl::VelocityShift)
    splittingMode = PIC::ParticleSplitting::_VelocityShift;
  else if (gConfiguration.populationControl == PopulationControl::Scatter)
    splittingMode = PIC::ParticleSplitting::_Scatter;
  PIC::ParticleSplitting::SetMode(splittingMode);
  if (splittingMode != PIC::ParticleSplitting::_disactivated) {
    PIC::ParticleSplitting::SetParam(
        gConfiguration.splitPositionShiftCellDiagonalFraction,
        gConfiguration.splitVelocityShiftThermalSpeedFraction,
        gConfiguration.minParticlesPerCell,
        gConfiguration.maxParticlesPerCell,
        gConfiguration.nonuniformSplitPositionShift);
  }
  // Pointer-mode dispatch is selected in input/enceladus.input. The wrapper
  // retains Boris gravity and internal-sphere handling and records domain
  // escape. In plasma-free mode LoadCoupledPlasmaBackground explicitly zeros
  // every allocated background field before any particle can be moved.
  PIC::Mover::SetUserDefinedParticleMover(&MoveParticle);
  PIC::BC::UserDefinedParticleInjectionFunction = InjectSurfaceParticles;
  if (gConfiguration.flybyDensityOutput || !gObservations.empty() ||
      !gConfiguration.volumeFieldOutput)
    PIC::Sampling::ExternalSamplingLocalVariables::RegisterSamplingRoutine(
        FlybySamplingNoOp, OutputAllObservationProducts);
  gNumericsFinalized = true;
  PrintReceipt();
}

void LoadCoupledPlasmaBackground() {
  if (!gInstalled) Stop("runtime configuration was not installed");
  if (gConfiguration.plasmaBackgroundMode == PlasmaBackgroundMode::None) {
    // Coupler capability remains compiled for other runs, so Boris will still
    // interpolate its allocated center-node arrays even for q=0 neutrals.
    // Initialize the complete data-file background region to exact zeros:
    // E=B=u=n=T=p=0 in their native SI fields. This makes "none" a defined
    // physical state and prevents uninitialized memory from influencing the
    // mover, while gravity remains active through the application callback.
    if (PIC::CPLR::DATAFILE::CenterNodeAssociatedDataOffsetBegin < 0 ||
        PIC::CPLR::DATAFILE::nTotalBackgroundVariables <= 0)
      Stop("plasma-free mode cannot locate allocated background fields");
    PIC::Mesh::SetCenterNodeAssociatedDataValue(
        0.0, PIC::CPLR::DATAFILE::nTotalBackgroundVariables,
        PIC::CPLR::DATAFILE::CenterNodeAssociatedDataOffsetBegin);
    MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
    return;
  }
  constexpr const char* cacheIdentity = "ENCELADUS-BATSRUS";
  if (PIC::CPLR::DATAFILE::BinaryFileExists(cacheIdentity)) {
    PIC::CPLR::DATAFILE::LoadBinaryFile(cacheIdentity);
    return;
  }
  if (gConfiguration.plasmaBackgroundMode ==
      PlasmaBackgroundMode::CacheOnly)
    Stop("plasma_background_mode=cache_only but no mesh-matched "
         "ENCELADUS-BATSRUS cache exists in " +
         gConfiguration.plasmaCacheDirectory);

  // The legacy ICES adapter extracts values at the allocated AMR cell
  // centers, then the generic reader populates coupled density, flow, E, B,
  // pressure, and temperature fields. Extraction is an explicit fallback;
  // cache_only can therefore guarantee that a production run never launches
  // an external executable unexpectedly.
#ifdef _ICES_CREATE_COORDINATE_LIST_
#if _ICES_CREATE_COORDINATE_LIST_ == _PIC_MODE_ON_
  PIC::CPLR::DATAFILE::ICES::createCellCenterCoordinateList();
  PIC::CPLR::DATAFILE::ICES::retriveSWMFdata(
      gConfiguration.icesCase.c_str());
#else
  Stop("cache_or_ices requires _ICES_CREATE_COORDINATE_LIST_=_PIC_MODE_ON_");
#endif
#else
  Stop("cache_or_ices requires _ICES_CREATE_COORDINATE_LIST_ at compile time");
#endif
  PIC::CPLR::DATAFILE::ICES::readSWMFdata(1.0);
  PIC::CPLR::DATAFILE::SaveBinaryFile(cacheIdentity);
}

double GlobalTimeStep(int, cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>*) {
  return InstalledNumericalPlan().globalTimeStepS;
}

double PlannedSourceRate(int species) {
  if (!gNumericsFinalized) return 0.0;
  const auto found = gPlan.physicalSourceRatePerS.find(ParticleSymbol(species));
  return found == gPlan.physicalSourceRatePerS.end() ? 0.0 : found->second;
}

long int InjectSurfaceParticles() {
  if (!gNumericsFinalized) return 0;
  const double dt = gPlan.globalTimeStepS;
  unsigned long long localInjected = 0;
  for (const std::string& symbol : gConfiguration.injectedSpecies) {
    const int species = Species(symbol);
    const double physicalRate = gPlan.physicalSourceRatePerS.at(symbol);
    const double weight = gPlan.statisticalWeight.at(symbol);
    const std::uint64_t requested =
        gConfiguration.primaryMacroparticlesPerStep.at(symbol);
    DeterministicRandom random(Mix(gConfiguration.randomSeed) ^
        Mix(gInjectionStep) ^ Mix(static_cast<std::uint64_t>(species)));
    // All ranks replay the identical event sequence and only the AMR owner
    // creates each particle. Consequently the global count is exactly the
    // requested integer N for each primary on every step. Uniform independent
    // birth times avoid an artificial synchronized sheet within the step.
    for (std::uint64_t event = 0; event < requested; ++event) {
      const double eventTime = random.Uniform()*dt;
      const Source& source = SelectSource(symbol, random.Uniform()*physicalRate);
      Vec3 position, normal, axis;
      SampleSourceGeometry(source, &random, &position, &normal, &axis);
      const Vec3 velocity = SampleTenishevVelocity(
          source, normal, axis, PIC::MolecularData::GetMass(species), &random);
      // Move the birth point 2^-40 radii outward. This is far below any AMR
      // scale but avoids classifying an exactly-on-sphere point as an incoming
      // collision through roundoff.
      const double radius = Norm(position);
      const double offset = std::ldexp(radius, -40);
      position.x += offset*normal.x;
      position.y += offset*normal.y;
      position.z += offset*normal.z;
      double x[3] = {position.x, position.y, position.z};
      cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node =
          PIC::Mesh::mesh->findTreeNode(x);
      if (node == nullptr) Stop("source point is outside the AMR tree: " + source.id);
      if (node->Thread != PIC::ThisThread) continue;
      if (node->block == nullptr) Stop("owner of source point has no allocated block");
      const long int particle = PIC::ParticleBuffer::GetNewParticle();
      PIC::ParticleBuffer::byte* data =
          PIC::ParticleBuffer::GetParticleDataPointer(particle);
      double v[3] = {velocity.x, velocity.y, velocity.z};
      PIC::ParticleBuffer::SetX(x, data);
      PIC::ParticleBuffer::SetV(v, data);
      PIC::ParticleBuffer::SetI(species, data);
      PIC::ParticleBuffer::SetIndividualStatWeightCorrection(1.0, data);
      ++localInjected;
      Increment(&gLocalInjected, species);
      ++PIC::BC::nInjectedParticles[species];
      PIC::BC::ParticleProductionRate[species] += weight/dt;
      PIC::BC::ParticleMassProductionRate[species] +=
          weight*PIC::MolecularData::GetMass(species)/dt;
      PIC::Mover::DispatchParticleMover(particle, dt-eventTime, node);
    }
  }
  ++gInjectionStep;
  if (localInjected > static_cast<unsigned long long>(
          std::numeric_limits<long int>::max()))
    Stop("rank-local injected particle count overflow");
  return static_cast<long int>(localInjected);
}

int MoveParticle(long int particle, double timeStep,
                 cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node) {
  PIC::ParticleBuffer::byte* data =
      PIC::ParticleBuffer::GetParticleDataPointer(particle);
  const int species = PIC::ParticleBuffer::GetI(data);
  gSurfaceAbsorptionParticle = -1;
  const int result = PIC::Mover::Boris(particle, timeStep, node);
  // Boris uses _PARTICLE_LEFT_THE_DOMAIN_ for both deletion at the outer AMR
  // boundary and deletion by an internal-sphere callback.  The callback is
  // synchronous on this mover thread, so the thread-local particle marker
  // separates the mutually exclusive sinks without changing mover physics.
  if (result == _PARTICLE_LEFT_THE_DOMAIN_ &&
      gSurfaceAbsorptionParticle != particle)
    Increment(&gLocalDomainEscape, species);
  gSurfaceAbsorptionParticle = -1;
  return result;
}

void ProcessParticle(long int particle, long int& firstParticleCell,
                     cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node) {
  PIC::ParticleBuffer::byte* parentData =
      PIC::ParticleBuffer::GetParticleDataPointer(particle);
  const int parentIndex = PIC::ParticleBuffer::GetI(parentData);
  const std::string parent = ParticleSymbol(parentIndex);
  if (!gNumericsFinalized || !IsSimulated(parent) ||
      !IsChemistrySpecies(parent)) {
    LinkParticle(particle, firstParticleCell);
    return;
  }
  const std::vector<const Reaction*> channels = Channels(parent);
  if (channels.empty()) {
    LinkParticle(particle, firstParticleCell);
    return;
  }

  double x[3], parentVelocity[3];
  PIC::ParticleBuffer::GetX(x, parentData);
  PIC::ParticleBuffer::GetV(parentVelocity, parentData);
  double totalRate = 0.0;
  for (const Reaction* reaction : channels) {
    double plasmaVelocity[3] = {};
    const double rate = LocalReactionRate(
        *reaction, x, parentVelocity, node, plasmaVelocity);
    totalRate += rate;
  }
  const double dt = gPlan.globalTimeStepS;
  if (!(totalRate > 0.0) || rnd() >= -std::expm1(-totalRate*dt)) {
    LinkParticle(particle, firstParticleCell);
    return;
  }
  // Application schema 2 represents photolysis solely as loss of the neutral
  // parent. The event probability is exact for constant rates over this
  // global step, 1-exp[-dt sum(nu_r)]. Physical channel products remain table
  // metadata but no daughter particle, charge, momentum, or weight is created.
  Increment(&gLocalPhotolyticLoss, parentIndex);
  PIC::ParticleBuffer::DeleteParticle(particle);
}

int AbsorbParticleAtSurface(int species, long int particle, double*, double*,
                            double&, void*, void*) {
  // Generic PIC deletes the particle after this disposition is returned.
  // Deleting it here as well would release the particle-buffer slot twice.
  gSurfaceAbsorptionParticle = particle;
  Increment(&gLocalSurfaceLoss, species);
  return _PARTICLE_DELETED_ON_THE_FACE_;
}

void TotalParticleAcceleration(
    double* acceleration, int species, long int particle, double* position,
    double* velocity, cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node) {
  if (!gInstalled) {
    ElectricallyChargedDust::TotalGrainAcceleration(
        acceleration, species, particle, position, velocity, node);
    return;
  }
  // The coupled electric/magnetic force for charged species is owned by the
  // selected PIC mover. This callback adds only Enceladus point-mass gravity
  // in the body-centred frame, a=-GM r/r^3 [m s^-2].
  const double radius2 = position[0]*position[0] + position[1]*position[1] +
                         position[2]*position[2];
  const double radius = std::sqrt(radius2);
  if (!(radius > 0.0)) Stop("particle acceleration requested at r=0");
  const double factor = -GravityConstant*_MASS_(_TARGET_)/(radius2*radius);
  for (int dimension = 0; dimension < 3; ++dimension)
    acceleration[dimension] = factor*position[dimension];
}

int TransformationIndicator(
    double* position, double* velocity, int species, long int particle,
    PIC::ParticleBuffer::byte* particleData, double& dt, bool& limitTimeStep,
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node) {
  if (gInstalled) {
    limitTimeStep = false;
    return _GENERIC_PARTICLE_TRANSFORMATION_CODE__NO_TRANSFORMATION_;
  }
  return ElectricallyChargedDust::DustChargingProcessorIndicator(
      position, velocity, species, particle, particleData, dt, limitTimeStep,
      node);
}

int TransformationProcessor(
    double* initialPosition, double* finalPosition, double* velocity,
    int& species, long int particle, PIC::ParticleBuffer::byte* particleData,
    double dt, cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node) {
  if (gInstalled)
    return _GENERIC_PARTICLE_TRANSFORMATION_CODE__NO_TRANSFORMATION_;
  return ElectricallyChargedDust::DustChargingProcessor_Implicit_SecondOrder(
      initialPosition, finalPosition, velocity, species, particle,
      particleData, dt, node);
}

void PrintReceipt() {
  if (PIC::ThisThread != 0 || !gNumericsFinalized) return;
  std::ostringstream out;
  out << std::scientific << std::setprecision(17)
      << "[srcEnceladus] forward_model source_case=" << gDistribution.sourceCase
      << " source_mode=" << Name(gDistribution.mode)
      << " global_dt_s=" << gPlan.globalTimeStepS
      << " limiting_species=" << gPlan.limitingSpecies
      << " sample_length_steps=" << gConfiguration.sampleLengthSteps
      << " population_control=" << Name(gConfiguration.populationControl)
      << " reaction_products=" << Name(gConfiguration.reactionProducts)
      << " stripe_surface_mapping="
      << Name(gDistribution.stripeSurfaceMapping)
      << " plasma_background_mode="
      << Name(gConfiguration.plasmaBackgroundMode)
      << " total_iterations=" << gConfiguration.totalIterations
      << " volume_field_output="
      << (gConfiguration.volumeFieldOutput ? "on" : "off")
      << " flyby_density_output="
      << (gConfiguration.flybyDensityOutput ? "on" : "off")
      << " observation_requests=" << gObservations.size()
      << " plasma_back_coupling=off\n";
  if (gConfiguration.flybyDensityOutput)
    out << "[srcEnceladus] flyby_id=" << gConfiguration.flybyId
        << " flyby_trajectory_points=" << gFlybyTrajectory.size()
        << " flyby_trajectory_file=" << gConfiguration.flybyTrajectoryFile
        << '\n';
  if (gConfiguration.flybyMeshRefinement)
    out << "[srcEnceladus] flyby_mesh_refinement=on"
        << " target_cell_size_m="
        << gConfiguration.flybyRefinementTargetCellSizeM
        << " tube_radius_m=" << gConfiguration.flybyRefinementRadiusM
        << " time_min_s=" << gConfiguration.flybyRefinementTimeMinS
        << " time_max_s=" << gConfiguration.flybyRefinementTimeMaxS
        << " refinement_trajectory_points="
        << gFlybyRefinementTrajectory.size() << '\n';
  for (const LoadedObservation& observation : gObservations)
    out << "[srcEnceladus] observation_id=" << observation.request.id
        << " instrument=" << Name(observation.request.instrument)
        << " event_id=" << observation.request.eventId
        << " rows=" << observation.points.size()
        << " comparison_species="
        << (observation.request.comparisonSpecies.empty()
                ? "none" : observation.request.comparisonSpecies)
        << (observation.request.contextObservable.empty()
                ? std::string()
                : std::string(" context_observable=") +
                      observation.request.contextObservable)
        << (observation.request.instrument == ObservationInstrument::Uvis
                ? std::string(" ray_geometry=") +
                      Name(observation.request.uvisRayGeometry)
                : std::string())
        << '\n';
  if (gConfiguration.plasmaBackgroundMode != PlasmaBackgroundMode::None)
    out << "[srcEnceladus] plasma_cache_directory="
        << gConfiguration.plasmaCacheDirectory << '\n';
  for (const auto& item : gPlan.statisticalWeight)
    out << "[srcEnceladus] species=" << item.first
        << " rate_s-1=" << gPlan.physicalSourceRatePerS.at(item.first)
        << " expected_macroparticles_per_step="
        << gPlan.expectedMacroparticlesPerStep.at(item.first)
        << " statistical_weight=" << item.second << '\n';
  for (const std::string& reaction : gConfiguration.selectedReactions)
    out << "[srcEnceladus] selected_reaction=" << reaction << '\n';
  std::cout << out.str();
}

void PrintFinalReceipt(long completedIterations) {
  if (!gNumericsFinalized) return;
  unsigned long long minimumInjectionSteps = gInjectionStep;
  unsigned long long maximumInjectionSteps = gInjectionStep;
  MPI_Allreduce(MPI_IN_PLACE, &minimumInjectionSteps, 1,
                MPI_UNSIGNED_LONG_LONG, MPI_MIN, MPI_GLOBAL_COMMUNICATOR);
  MPI_Allreduce(MPI_IN_PLACE, &maximumInjectionSteps, 1,
                MPI_UNSIGNED_LONG_LONG, MPI_MAX, MPI_GLOBAL_COMMUNICATOR);
  const long int resident = PIC::ParticleBuffer::GetTotalParticleNumber();
  std::vector<unsigned long long> residentBySpecies(PIC::nTotalSpecies, 0);
  PIC::DomainBlockDecomposition::UpdateBlockTable();
  for (unsigned int blockIndex = 0;
       blockIndex < PIC::DomainBlockDecomposition::nLocalBlocks; ++blockIndex) {
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node =
        PIC::DomainBlockDecomposition::BlockTable[blockIndex];
    if (node == nullptr || node->block == nullptr) continue;
    const int cellCount = _BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_;
    for (int cell = 0; cell < cellCount; ++cell) {
      for (long int particle = node->block->FirstCellParticleTable[cell];
           particle != -1; particle = PIC::ParticleBuffer::GetNext(particle)) {
        const int species = PIC::ParticleBuffer::GetI(particle);
        if (species < 0 || species >= PIC::nTotalSpecies)
          Stop("resident particle has an invalid compiled species index");
        ++residentBySpecies[species];
      }
    }
  }
  MPI_Allreduce(MPI_IN_PLACE, residentBySpecies.data(), PIC::nTotalSpecies,
                MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_GLOBAL_COMMUNICATOR);
  std::ostringstream out;
  if (PIC::ThisThread == 0)
    out << "[srcEnceladus] final_receipt completed_iterations="
        << completedIterations << " injection_step_min="
        << minimumInjectionSteps << " injection_step_max="
        << maximumInjectionSteps << " resident_macroparticles="
        << resident << '\n';
  for (const std::string& symbol : gConfiguration.injectedSpecies) {
    const int species = Species(symbol);
    unsigned long long counters[4] = {
        gLocalInjected[species], gLocalPhotolyticLoss[species],
        gLocalSurfaceLoss[species], gLocalDomainEscape[species]};
    MPI_Allreduce(MPI_IN_PLACE, counters, 4, MPI_UNSIGNED_LONG_LONG, MPI_SUM,
                  MPI_GLOBAL_COMMUNICATOR);
    if (PIC::ThisThread == 0) {
      const unsigned long long expected =
          static_cast<unsigned long long>(completedIterations) *
          gConfiguration.primaryMacroparticlesPerStep.at(symbol);
      const long long ledgerResidual =
          static_cast<long long>(counters[0]) -
          static_cast<long long>(counters[1]) -
          static_cast<long long>(counters[2]) -
          static_cast<long long>(counters[3]) -
          static_cast<long long>(residentBySpecies[species]);
      out << "[srcEnceladus] final_species=" << symbol
          << " injected_macroparticles=" << counters[0]
          << " expected_injected_macroparticles=" << expected
          << " photolytic_parent_losses=" << counters[1]
          << " surface_losses=" << counters[2]
          << " domain_escapes=" << counters[3]
          << " resident_macroparticles=" << residentBySpecies[species]
          << " ledger_residual=" << ledgerResidual << '\n';
    }
  }
  if (PIC::ThisThread == 0) std::cout << out.str();
}

}  // namespace AMPS
}  // namespace ForwardModel
}  // namespace Enceladus
