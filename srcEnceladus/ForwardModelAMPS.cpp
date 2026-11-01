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
  // each such file is hundreds of MB, so retain generic fields only for the
  // explicitly simulated species. This changes output selection only: storage,
  // mover, chemistry, and sampling behavior for active particles are untouched.
  if (PIC::Sampling::SaveOutputDataFile == nullptr)
    Stop("generic AMPS output-selection table is not initialized");
  for (int species = 0; species < PIC::nTotalSpecies; ++species)
    PIC::Sampling::SaveOutputDataFile[species] =
        IsSimulated(ParticleSymbol(species));

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
  if (gConfiguration.flybyDensityOutput)
    PIC::Sampling::ExternalSamplingLocalVariables::RegisterSamplingRoutine(
        FlybySamplingNoOp, OutputFlybyDensity);
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
      << " flyby_density_output="
      << (gConfiguration.flybyDensityOutput ? "on" : "off")
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
