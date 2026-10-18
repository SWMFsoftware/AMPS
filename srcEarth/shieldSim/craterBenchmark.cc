/*
 * shieldSimCraterBenchmark
 * ========================
 *
 * A deliberately bounded Geant4 benchmark of the CRaTER foreground response.
 * It transports only zenith-side rays whose unperturbed trajectories cross the
 * published D1 and D6 apertures.  It is not a CRaTER engineering mass model:
 * housing, boards, dead layers, thermal blanket, spacecraft, lunar albedo, and
 * side-entry/background particles are absent.  Every output repeats the scope
 * token "foreground_only_axial" so a plot cannot reasonably be presented as
 * a full-instrument calculation.
 *
 * The axial dimensions, material compositions, detector aperture, geometry
 * factor, and electronics thresholds come from Spence et al. (2010).  The
 * D5/D6 selection logic and 20 keV/um splice come from Looper et al. (2020)
 * and its supporting information.  BON2020 spectra and the observed D5/D6
 * spectrum are checksum-pinned public fixtures selected by the input case.
 */

#include "BuildProvenance.hh"
#include "CraterBenchmark.hh"
#include "StrictJson.hh"

#include <G4Box.hh>
#include <G4Element.hh>
#include <G4EmStandardPhysics_option4.hh>
#include <G4Event.hh>
#include <G4IonTable.hh>
#include <G4LogicalVolume.hh>
#include <G4Material.hh>
#include <G4NistManager.hh>
#include <G4ParticleDefinition.hh>
#include <G4ParticleGun.hh>
#include <G4ParticleTable.hh>
#include <G4PhysListFactory.hh>
#include <G4PVPlacement.hh>
#include <G4RunManager.hh>
#include <G4Step.hh>
#include <G4SystemOfUnits.hh>
#include <G4ThreeVector.hh>
#include <G4Tubs.hh>
#include <G4VModularPhysicsList.hh>
#include <G4VPhysicalVolume.hh>
#include <G4VUserDetectorConstruction.hh>
#include <G4VUserPrimaryGeneratorAction.hh>
#include <G4UserEventAction.hh>
#include <G4UserSteppingAction.hh>
#include <G4Version.hh>
#include <G4ios.hh>
#include <Randomize.hh>

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using shield_sim::CraterBenchmarkCase;
using shield_sim::CraterBonSpectrum;
using shield_sim::CraterDeposits;
using shield_sim::CraterObservation;

constexpr double kSecondsPerDay = 86400.0;
constexpr double kSquareCentimetersPerSquareMeter = 10000.0;

struct WeightedBins {
  std::vector<double> sum;
  std::vector<double> sumSquares;

  explicit WeightedBins(std::size_t size) : sum(size, 0.0), sumSquares(size, 0.0) {}
  void Add(std::size_t bin, double value) {
    sum.at(bin) += value;
    sumSquares.at(bin) += value * value;
  }
};

struct DetectorPlanes {
  double stackFrontMm = 0.0;
  double stackBackMm = 0.0;
  double d1FrontMm = 0.0;
  double d6FrontMm = 0.0;
};

struct ParallelConfig {
  std::int64_t shardIndex = 0;
  std::int64_t shardCount = 1;
  std::int64_t localHistories = 0;
  std::int64_t effectiveSeed = 0;
};

std::int64_t ParseNonnegativeInteger(const std::string& text,
                                     const std::string& name) {
  std::size_t consumed = 0;
  long long value = 0;
  try { value = std::stoll(text, &consumed); }
  catch (const std::exception&) { throw std::runtime_error(name+" is not an integer"); }
  if (consumed != text.size() || value < 0)
    throw std::runtime_error(name+" must be a nonnegative integer");
  return static_cast<std::int64_t>(value);
}

ParallelConfig ResolveParallelConfig(const CraterBenchmarkCase& definition,
                                     int argc, char** argv) {
  bool commandLineParallel = false;
  ParallelConfig result;
  for (int i = 2; i < argc; ++i) {
    const std::string argument = argv[i];
    const std::string indexPrefix = "--shard-index=";
    const std::string countPrefix = "--shard-count=";
    if (argument.rfind(indexPrefix, 0) == 0) {
      result.shardIndex = ParseNonnegativeInteger(argument.substr(indexPrefix.size()),
                                                  "shard index");
      commandLineParallel = true;
    } else if (argument.rfind(countPrefix, 0) == 0) {
      result.shardCount = ParseNonnegativeInteger(argument.substr(countPrefix.size()),
                                                  "shard count");
      commandLineParallel = true;
    } else {
      throw std::runtime_error("unknown CRaTER benchmark argument '"+argument+"'");
    }
  }
  if (!commandLineParallel) {
    // No MPI library is required: standard launchers publish rank and size in
    // the process environment.  Supporting Open MPI and PMI covers mpiexec,
    // MPICH, and common scheduler launchers without coupling Geant4 to MPI.
    const char* rank = std::getenv("OMPI_COMM_WORLD_RANK");
    const char* size = std::getenv("OMPI_COMM_WORLD_SIZE");
    if (!rank || !size) {
      rank = std::getenv("PMI_RANK");
      size = std::getenv("PMI_SIZE");
    }
    if ((rank == nullptr) != (size == nullptr))
      throw std::runtime_error("MPI environment provides only one of rank/size");
    if (rank && size) {
      result.shardIndex = ParseNonnegativeInteger(rank, "MPI rank");
      result.shardCount = ParseNonnegativeInteger(size, "MPI size");
    }
  }
  if (result.shardCount < 1 || result.shardIndex >= result.shardCount)
    throw std::runtime_error("require 0 <= shard_index < shard_count");
  if (result.shardCount > definition.histories)
    throw std::runtime_error("shard_count cannot exceed the requested histories");
  const std::int64_t quotient = definition.histories / result.shardCount;
  const std::int64_t remainder = definition.histories % result.shardCount;
  result.localHistories = quotient + (result.shardIndex < remainder ? 1 : 0);
  if (definition.randomSeed > std::numeric_limits<std::int64_t>::max()-result.shardIndex)
    throw std::runtime_error("random_seed + shard_index overflows");
  // Unique deterministic seed identifiers give each rank a distinct Geant4
  // pseudo-random stream.  Statistical independence is the documented RNG
  // property; this does not claim a formal proof that finite stream segments
  // can never share a value.
  result.effectiveSeed = definition.randomSeed + result.shardIndex;
  return result;
}

DetectorPlanes AxialPlanes() {
  DetectorPlanes result;
  double total = 0.0;
  for (const auto& element : shield_sim::PublishedCraterAxialStack())
    total += element.thicknessMm;
  result.stackFrontMm = -0.5 * total;
  result.stackBackMm = 0.5 * total;
  double front = result.stackFrontMm;
  for (const auto& element : shield_sim::PublishedCraterAxialStack()) {
    if (element.id == "D1") result.d1FrontMm = front;
    if (element.id == "D6") result.d6FrontMm = front;
    front += element.thicknessMm;
  }
  return result;
}

std::size_t Locate(const std::vector<double>& edges, double value) {
  if (value < edges.front() || value > edges.back()) return edges.size();
  if (value == edges.back()) return edges.size() - 2;
  return static_cast<std::size_t>(
      std::upper_bound(edges.begin(), edges.end(), value) - edges.begin() - 1);
}

double StandardErrorOfMean(double sum, double sumSquares, std::int64_t histories) {
  if (histories < 2) return 0.0;
  const double n = static_cast<double>(histories);
  const double centered = std::max(0.0, sumSquares - sum * sum / n);
  return std::sqrt(centered / (n * (n - 1.0)));
}

shield_sim::Json NumberOrNull(double value) {
  return std::isfinite(value) ? shield_sim::Json(value) : shield_sim::Json(nullptr);
}

std::string SafeFileToken(std::string value) {
  for (char& c : value)
    if (!std::isalnum(static_cast<unsigned char>(c)) && c != '_' && c != '-') c = '_';
  return value;
}

class BenchmarkState {
public:
  BenchmarkState(CraterBenchmarkCase definition, CraterBonSpectrum source,
                 CraterObservation observation, ParallelConfig parallel)
      : definition(std::move(definition)), source(std::move(source)),
        observation(std::move(observation)),
        parallel(std::move(parallel)),
        edges(shield_sim::CraterObservationBinEdges(this->observation.centersKevPerUm)),
        ideal(edges.size() - 1), afterCuts(edges.size() - 1) {}

  void BeginEvent() {
    depositedEnergy.fill(0.0);
    sourceWeight = 0.0;
  }

  void ScoreDetector(int detectorIndex, double energy) {
    depositedEnergy.at(static_cast<std::size_t>(detectorIndex)) += energy;
  }

  void EndEvent() {
    ++completedHistories;
    CraterDeposits normalized;
    static constexpr std::array<double, 6> thicknessUm =
        {148.0, 1000.0, 149.0, 1000.0, 149.0, 1000.0};
    for (std::size_t i = 0; i < 6; ++i)
      normalized.kevPerUm[i] = depositedEnergy[i] / keV / thicknessUm[i];

    const double observable = shield_sim::CraterCombinedObservable(normalized);
    if (observable < 0.0) return;
    ++coincidenceEvents;
    const std::size_t bin = Locate(edges, observable);
    if (bin >= ideal.sum.size()) {
      if (observable < edges.front()) ++underflowEvents;
      else ++overflowEvents;
      return;
    }
    ideal.Add(bin, sourceWeight);
    if (shield_sim::CraterCombinedPublishedCut(normalized)) {
      ++cutAcceptedEvents;
      afterCuts.Add(bin, sourceWeight);
    }
  }

  void WriteResults() const {
    namespace fs = std::filesystem;
    fs::create_directories(definition.outputDirectory);
    std::string stem = SafeFileToken(definition.caseId) + "__foreground_only_axial";
    if (parallel.shardCount > 1) {
      std::ostringstream suffix;
      suffix << "__shard_" << std::setw(4) << std::setfill('0') << parallel.shardIndex
             << "_of_" << std::setw(4) << parallel.shardCount << "_partial";
      stem += suffix.str();
    }
    const fs::path jsonPath = fs::path(definition.outputDirectory) / (stem + "_result.json");
    const fs::path csvPath = fs::path(definition.outputDirectory) / (stem + "_spectrum.csv");

    shield_sim::Json::Array bins;
    std::ofstream csv(csvPath);
    if (!csv) throw std::runtime_error("cannot create CRaTER spectrum CSV");
    csv << "lower_keV_per_um,upper_keV_per_um,center_keV_per_um,"
           "ideal_foreground_per_cm2_sr_day_keV_per_um,ideal_mc_standard_error,"
           "foreground_after_cuts_per_cm2_sr_day_keV_per_um,cuts_mc_standard_error,"
           "observed_per_cm2_sr_day_keV_per_um\n";

    const double n = static_cast<double>(completedHistories);
    std::size_t assessed = 0;
    std::size_t meetingTarget = 0;
    double idealIntegrated = 0.0;
    double cutIntegrated = 0.0;
    double observedIntegrated = 0.0;
    double logResidualSum = 0.0;
    std::size_t logResidualCount = 0;
    for (std::size_t i = 0; i < ideal.sum.size(); ++i) {
      const double width = edges[i + 1] - edges[i];
      // BON2020 is per m2 s sr.  Multiplication by seconds/day and division by
      // cm2/m2 gives the units printed on Looper et al. (2020), Figure 4.
      const double conversion = kSecondsPerDay /
          (kSquareCentimetersPerSquareMeter * width);
      const double idealValue = ideal.sum[i] / n * conversion;
      const double idealError = StandardErrorOfMean(
          ideal.sum[i], ideal.sumSquares[i], completedHistories) * conversion;
      const double cutValue = afterCuts.sum[i] / n * conversion;
      const double cutError = StandardErrorOfMean(
          afterCuts.sum[i], afterCuts.sumSquares[i], completedHistories) * conversion;
      const double observedValue = observation.differentialIntensity[i];
      idealIntegrated += idealValue * width;
      cutIntegrated += cutValue * width;
      observedIntegrated += observedValue * width;
      if (observation.centersKevPerUm[i] >= 0.5 &&
          observation.centersKevPerUm[i] <= 10.0 && cutValue > 0.0) {
        ++assessed;
        if (cutError / cutValue <= definition.targetRelativeStatisticalUncertainty)
          ++meetingTarget;
      }
      if (cutValue > 0.0 && observedValue > 0.0) {
        const double residual = std::log(cutValue / observedValue);
        logResidualSum += residual * residual;
        ++logResidualCount;
      }
      bins.emplace_back(shield_sim::Json::Object{
          {"lower_keV_per_um", edges[i]},
          {"upper_keV_per_um", edges[i + 1]},
          {"center_keV_per_um", observation.centersKevPerUm[i]},
          {"ideal_foreground_per_cm2_sr_day_keV_per_um", idealValue},
          {"ideal_mc_standard_error", idealError},
          {"foreground_after_cuts_per_cm2_sr_day_keV_per_um", cutValue},
          {"cuts_mc_standard_error", cutError},
          {"observed_per_cm2_sr_day_keV_per_um", observedValue},
          {"ideal_raw_weight_sum", ideal.sum[i]},
          {"ideal_raw_weight_sum_squares", ideal.sumSquares[i]},
          {"cuts_raw_weight_sum", afterCuts.sum[i]},
          {"cuts_raw_weight_sum_squares", afterCuts.sumSquares[i]},
      });
      csv << std::scientific << std::setprecision(12)
          << edges[i] << ',' << edges[i + 1] << ','
          << observation.centersKevPerUm[i] << ','
          << idealValue << ',' << idealError << ','
          << cutValue << ',' << cutError << ',' << observedValue << '\n';
    }

    shield_sim::Json::Array geometry;
    double cumulative = 0.0;
    for (const auto& element : shield_sim::PublishedCraterAxialStack()) {
      cumulative += element.thicknessMm;
      geometry.emplace_back(shield_sim::Json::Object{
          {"id", element.id}, {"material", element.material},
          {"thickness_mm", element.thicknessMm},
          {"cumulative_back_face_mm", cumulative},
          {"sensitive", element.sensitive},
          {"detector_number", element.detectorNumber},
      });
    }
    shield_sim::Json::Array ions;
    for (const auto& ion : shield_sim::PublishedCraterPrimaryIons())
      ions.emplace_back(shield_sim::Json::Object{
          {"Z", ion.z}, {"A", ion.a}, {"symbol", ion.symbol}});

    const double geometryFactorM2Sr =
        shield_sim::kCraterD1D6GeometryFactorCm2Sr / kSquareCentimetersPerSquareMeter;
    const double idealRatePerDay = idealIntegrated *
        shield_sim::kCraterD1D6GeometryFactorCm2Sr;
    const double cutRatePerDay = cutIntegrated *
        shield_sim::kCraterD1D6GeometryFactorCm2Sr;
    const double logRmse = logResidualCount
        ? std::sqrt(logResidualSum / static_cast<double>(logResidualCount))
        : std::numeric_limits<double>::quiet_NaN();

    shield_sim::Json document(shield_sim::Json::Object{
      {"schema_version", "shieldSim.crater_foreground_axial_result.v1"},
      {"case_id", definition.caseId},
      {"benchmark_scope", shield_sim::kCraterBenchmarkScope},
      {"validation_status",
       parallel.shardCount > 1
         ? "partial_shard_requires_merge_not_stage8_acceptance_or_full_instrument_validation"
         : "diagnostic_benchmark_not_stage8_acceptance_or_full_instrument_validation"},
      {"result_type", parallel.shardCount > 1 ? "partial_shard" : "complete_serial"},
      {"parallel", shield_sim::Json::Object{
          {"shard_index", static_cast<double>(parallel.shardIndex)},
          {"shard_count", static_cast<double>(parallel.shardCount)},
          {"local_histories", static_cast<double>(completedHistories)},
          {"global_requested_histories", static_cast<double>(definition.histories)},
          {"base_random_seed", static_cast<double>(definition.randomSeed)},
          {"effective_random_seed", static_cast<double>(parallel.effectiveSeed)},
          {"merge_status", parallel.shardCount > 1 ? "required" : "not_applicable"},
      }},
      {"observable",
       "measured_energy_deposit_divided_by_detector_thickness_not_unrestricted_LET"},
      {"comparison_units", "particles / (cm2 sr day (keV/um))"},
      {"source", shield_sim::Json::Object{
          {"path", definition.sourcePath}, {"sha256", definition.sourceSha256},
          {"model_version", source.modelVersion}, {"run_dates", source.runDates},
          {"modulation_potential_MV", source.modulationPotentialMv},
          {"interpolation", definition.sourceInterpolation},
          {"transport_energy_domain",
           "10-50000 MeV/n for H; 10-10000 MeV/n for Z=2..28"},
          {"importance_sampling",
           "uniform_Z_1_to_28_and_uniform_log_energy_with_exact_event_weights"},
      }},
      {"observation", shield_sim::Json::Object{
          {"path", definition.observationPath},
          {"sha256", definition.observationSha256},
          {"column", definition.observationColumn},
          {"uncertainty_status", "not_provided_in_public_v1_3_table"},
      }},
      {"transport", shield_sim::Json::Object{
          {"geant4_version", G4Version}, {"physics_list", definition.physicsList},
          {"em_physics", definition.emPhysics},
          {"production_cut_um", definition.productionCutUm},
          {"histories", static_cast<double>(completedHistories)},
          {"global_requested_histories", static_cast<double>(definition.histories)},
          {"random_seed", static_cast<double>(parallel.effectiveSeed)},
          {"code_revision", shield_sim_build::GitRevision},
          {"source_tree_state", shield_sim_build::GitState},
      }},
      {"geometry", shield_sim::Json::Object{
          {"reference", "Spence_et_al_2010_Tables_2_to_5"},
          {"detector_diameter_mm", 2.0 * shield_sim::kCraterDetectorRadiusMm},
          {"D1_D6_geometry_factor_cm2_sr", shield_sim::kCraterD1D6GeometryFactorCm2Sr},
          {"full_angle_field_of_view_deg", shield_sim::kCraterFullAngleFieldOfViewDeg},
          {"axial_elements", shield_sim::Json(std::move(geometry))},
          {"non_detector_lateral_boundary",
           "computationally_wide_axial_slabs_not_claimed_as_published_housing_geometry"},
      }},
      {"primary_ions", shield_sim::Json(std::move(ions))},
      {"event_selection", shield_sim::Json::Object{
          {"coincidence", "D2_and_D4_and_(D5_or_D6)_using_published_approximate_LLDs"},
          {"splice", "D6_when_D6_normalized_deposit_le_20_keV_per_um_else_D5"},
          {"published_background_cuts_output", true},
          {"incident_phase_space",
           "zenith_rays_conditioned_on_unperturbed_D1_and_D6_aperture_crossing"},
      }},
      {"event_counts", shield_sim::Json::Object{
          {"coincidence", static_cast<double>(coincidenceEvents)},
          {"after_published_cuts", static_cast<double>(cutAcceptedEvents)},
          {"histogram_underflow", static_cast<double>(underflowEvents)},
          {"histogram_overflow", static_cast<double>(overflowEvents)},
      }},
      {"accuracy", shield_sim::Json::Object{
          {"requested_relative_statistical_uncertainty",
           definition.targetRelativeStatisticalUncertainty},
          {"assessment_range_keV_per_um", "0.5_to_10"},
          {"nonzero_bins_assessed", static_cast<double>(assessed)},
          {"bins_meeting_target", static_cast<double>(meetingTarget)},
          {"status", assessed > 0 && meetingTarget == assessed
              ? "requested_MC_precision_met_in_all_nonzero_assessed_bins"
              : "requested_MC_precision_not_met_in_all_nonzero_assessed_bins"},
      }},
      {"diagnostic_comparison", shield_sim::Json::Object{
          {"ideal_integral_per_cm2_sr_day", idealIntegrated},
          {"after_cuts_integral_per_cm2_sr_day", cutIntegrated},
          {"observed_integral_per_cm2_sr_day", observedIntegrated},
          {"ideal_telescope_rate_per_day", idealRatePerDay},
          {"after_cuts_telescope_rate_per_day", cutRatePerDay},
          {"published_geometry_factor_m2_sr", geometryFactorM2Sr},
          {"log_space_RMSE_after_cuts_vs_observation", NumberOrNull(logRmse)},
          {"acceptance_decision", "not_defined_no_observational_covariance_or_preregistered_threshold"},
      }},
      {"limitations", shield_sim::Json::Array{
          "foreground_only_no_side_entry_or_out_of_aperture_background",
          "no_sensor_housing_electronics_board_spacecraft_or_lunar_albedo",
          "no_dead_layers_or_thermal_blanket_because_public_dimensions_are_unavailable",
          "published_axial_non_detector_diameters_are_only_described_as_greater_than_35_mm",
          "approximate_published_LLDs_not_time_resolved_operational_thresholds",
          "no_ADC_gain_offset_resolution_saturation_or_D5_8_percent_empirical_gain_adjustment",
          "BON2020_411.143_MV_is_not_the_published_BON2010_417_MV_input",
          "current_Geant4_version_is_not_the_published_10.5_version",
          "observational_table_has_no_bin_covariance_or_uncertainty",
          "observational_bin_edges_are_geometric_midpoint_reconstructions_from_released_centers",
      }},
      {"bins", shield_sim::Json(std::move(bins))},
    });
    std::ofstream json(jsonPath);
    if (!json) throw std::runtime_error("cannot create CRaTER result JSON");
    json << document.DumpCanonical() << '\n';
    G4cout << "CRaTER foreground-only result: " << jsonPath << G4endl;
    G4cout << "CRaTER foreground-only spectrum: " << csvPath << G4endl;
  }

  CraterBenchmarkCase definition;
  CraterBonSpectrum source;
  CraterObservation observation;
  ParallelConfig parallel;
  std::vector<double> edges;
  WeightedBins ideal;
  WeightedBins afterCuts;
  std::array<double, 6> depositedEnergy{};
  std::array<G4LogicalVolume*, 6> detectorVolumes{};
  double sourceWeight = 0.0;
  std::int64_t completedHistories = 0;
  std::int64_t coincidenceEvents = 0;
  std::int64_t cutAcceptedEvents = 0;
  std::int64_t underflowEvents = 0;
  std::int64_t overflowEvents = 0;
};

class CraterDetectorConstruction final : public G4VUserDetectorConstruction {
public:
  explicit CraterDetectorConstruction(BenchmarkState& state) : state_(state) {}

  G4VPhysicalVolume* Construct() override {
    auto* nist = G4NistManager::Instance();
    auto* vacuum = nist->FindOrBuildMaterial("G4_Galactic");
    auto* silicon = nist->FindOrBuildMaterial("G4_Si");

    // A-150 mass fractions and density are from Spence et al. (2010), Table 2.
    // The percentages sum to exactly 100%; no residual element is invented.
    auto* a150 = new G4Material("CRATER_A150_TEP", 1.127*g/cm3, 6);
    a150->AddElement(nist->FindOrBuildElement("H"),  0.1033);
    a150->AddElement(nist->FindOrBuildElement("C"),  0.7693);
    a150->AddElement(nist->FindOrBuildElement("N"),  0.0330);
    a150->AddElement(nist->FindOrBuildElement("O"),  0.0693);
    a150->AddElement(nist->FindOrBuildElement("F"),  0.0114);
    a150->AddElement(nist->FindOrBuildElement("Ca"), 0.0137);

    // Flight endcap Al 6061-T6 composition and density are from Table 3.
    auto* aluminum = new G4Material("CRATER_Al6061_T6", 2.70*g/cm3, 9);
    aluminum->AddElement(nist->FindOrBuildElement("Al"), 0.97200);
    aluminum->AddElement(nist->FindOrBuildElement("Cr"), 0.00195);
    aluminum->AddElement(nist->FindOrBuildElement("Cu"), 0.00275);
    aluminum->AddElement(nist->FindOrBuildElement("Fe"), 0.00455);
    aluminum->AddElement(nist->FindOrBuildElement("Mg"), 0.01000);
    aluminum->AddElement(nist->FindOrBuildElement("Mn"), 0.00075);
    aluminum->AddElement(nist->FindOrBuildElement("Si"), 0.00600);
    aluminum->AddElement(nist->FindOrBuildElement("Ti"), 0.00075);
    aluminum->AddElement(nist->FindOrBuildElement("Zn"), 0.00125);

    const DetectorPlanes planes = AxialPlanes();
    auto* worldSolid = new G4Box("CRATER_world_solid", 60.0*mm, 60.0*mm,
                                 (0.5*(planes.stackBackMm-planes.stackFrontMm)+10.0)*mm);
    auto* worldLogical = new G4LogicalVolume(worldSolid, vacuum, "CRATER_world_logical");
    auto* world = new G4PVPlacement(nullptr, G4ThreeVector(), worldLogical,
                                    "CRATER_world_physical", nullptr, false, 0, true);

    double front = planes.stackFrontMm;
    for (const auto& element : shield_sim::PublishedCraterAxialStack()) {
      const double center = front + 0.5 * element.thicknessMm;
      front += element.thicknessMm;
      if (element.material == "vacuum") continue;
      G4Material* material = element.material == "silicon" ? silicon
                             : element.material == "A150_TEP" ? a150 : aluminum;
      G4VSolid* solid = nullptr;
      if (element.sensitive) {
        solid = new G4Tubs((element.id+"_solid").c_str(), 0.0,
            shield_sim::kCraterDetectorRadiusMm*mm,
            0.5*element.thicknessMm*mm, 0.0, 360.0*deg);
      } else {
        // Published sources give only the axial thickness and say these parts
        // are slightly wider than 35 mm.  A wide slab preserves the published
        // foreground areal density without claiming an unpublished diameter.
        solid = new G4Box((element.id+"_axial_slab_solid").c_str(),
                          50.0*mm, 50.0*mm, 0.5*element.thicknessMm*mm);
      }
      auto* logical = new G4LogicalVolume(
          solid, material, (element.id+"_logical").c_str());
      new G4PVPlacement(nullptr, G4ThreeVector(0.0,0.0,center*mm), logical,
                        (element.id+"_physical").c_str(), worldLogical,
                        false, element.detectorNumber, true);
      if (element.sensitive)
        state_.detectorVolumes.at(static_cast<std::size_t>(element.detectorNumber-1)) = logical;
    }
    return world;
  }

private:
  BenchmarkState& state_;
};

class CraterPrimaryGenerator final : public G4VUserPrimaryGeneratorAction {
public:
  explicit CraterPrimaryGenerator(BenchmarkState& state)
      : state_(state), gun_(std::make_unique<G4ParticleGun>(1)) {}

  void GeneratePrimaries(G4Event* event) override {
    state_.BeginEvent();
    const auto& ions = shield_sim::PublishedCraterPrimaryIons();
    std::size_t species = static_cast<std::size_t>(G4UniformRand() * ions.size());
    if (species == ions.size()) species = ions.size() - 1;
    const auto& ion = ions[species];

    const double lower = 10.0;
    const double requestedUpper = ion.z == 1 ? 50000.0 : 10000.0;
    const double upper = std::min(requestedUpper,
                                  state_.source.energyMeVPerNucleon.back());
    if (lower < state_.source.energyMeVPerNucleon.front() || !(upper > lower))
      throw std::runtime_error("BON2020 table does not cover the registered transport domain");
    const double logRange = std::log(upper / lower);
    const double energyPerNucleon = lower * std::exp(G4UniformRand() * logRange);
    const double flux = shield_sim::InterpolateCraterFlux(
        state_.source, ion.z, energyPerNucleon);
    // p(Z,E) = (1/28)/(E ln(Emax/Emin)).  Dividing J(E) by this proposal
    // density gives an unbiased contribution to integral J(E)dE summed over Z.
    state_.sourceWeight = flux * static_cast<double>(ions.size()) *
                          energyPerNucleon * logRange;

    G4ParticleDefinition* particle = nullptr;
    if (ion.z == 1 && ion.a == 1)
      particle = G4ParticleTable::GetParticleTable()->FindParticle("proton");
    else if (ion.z == 2 && ion.a == 4)
      particle = G4ParticleTable::GetParticleTable()->FindParticle("alpha");
    else
      particle = G4IonTable::GetIonTable()->GetIon(ion.z, ion.a, 0.0);
    if (!particle) throw std::runtime_error("Geant4 could not create a registered CRaTER ion");

    const DetectorPlanes planes = AxialPlanes();
    const double separation = planes.d6FrontMm - planes.d1FrontMm;
    const double thetaMax = 0.5 * shield_sim::kCraterFullAngleFieldOfViewDeg * deg;
    G4ThreeVector d1;
    G4ThreeVector direction;
    // Generate the isotropic plane-crossing phase-space measure (p(mu)=2mu),
    // then condition it on straight-line intersection of both D1 and D6.  The
    // published D1-D6 geometry factor supplies the rate normalization.
    for (;;) {
      const double radius = shield_sim::kCraterDetectorRadiusMm * std::sqrt(G4UniformRand());
      const double positionPhi = CLHEP::twopi * G4UniformRand();
      d1 = G4ThreeVector(radius*std::cos(positionPhi)*mm,
                         radius*std::sin(positionPhi)*mm,
                         planes.d1FrontMm*mm);
      const double cos2Min = std::pow(std::cos(thetaMax), 2);
      const double mu = std::sqrt(cos2Min + (1.0-cos2Min)*G4UniformRand());
      const double phi = CLHEP::twopi * G4UniformRand();
      const double transverse = std::sqrt(std::max(0.0, 1.0-mu*mu));
      direction = G4ThreeVector(transverse*std::cos(phi),
                                transverse*std::sin(phi), mu);
      const double dx = separation * direction.x() / direction.z();
      const double dy = separation * direction.y() / direction.z();
      const double x6 = d1.x()/mm + dx;
      const double y6 = d1.y()/mm + dy;
      if (x6*x6 + y6*y6 <= shield_sim::kCraterDetectorRadiusMm *
                             shield_sim::kCraterDetectorRadiusMm) break;
    }
    const double backDistance = (planes.d1FrontMm - planes.stackFrontMm) /
                                direction.z();
    const G4ThreeVector sourcePosition = d1 - backDistance*mm*direction -
                                         0.1*um*direction;
    gun_->SetParticleDefinition(particle);
    gun_->SetParticleEnergy(energyPerNucleon * ion.a * MeV);
    gun_->SetParticlePosition(sourcePosition);
    gun_->SetParticleMomentumDirection(direction);
    gun_->GeneratePrimaryVertex(event);
  }

private:
  BenchmarkState& state_;
  std::unique_ptr<G4ParticleGun> gun_;
};

class CraterSteppingAction final : public G4UserSteppingAction {
public:
  explicit CraterSteppingAction(BenchmarkState& state) : state_(state) {}
  void UserSteppingAction(const G4Step* step) override {
    const auto* volume = step->GetPreStepPoint()->GetPhysicalVolume();
    if (!volume) return;
    auto* logical = volume->GetLogicalVolume();
    for (std::size_t i = 0; i < state_.detectorVolumes.size(); ++i) {
      if (logical == state_.detectorVolumes[i]) {
        state_.ScoreDetector(static_cast<int>(i), step->GetTotalEnergyDeposit());
        return;
      }
    }
  }
private:
  BenchmarkState& state_;
};

class CraterEventAction final : public G4UserEventAction {
public:
  explicit CraterEventAction(BenchmarkState& state) : state_(state) {}
  void EndOfEventAction(const G4Event*) override { state_.EndEvent(); }
private:
  BenchmarkState& state_;
};

void PrintScope(const BenchmarkState& state) {
  G4cout << "============ CRaTER benchmark ============" << G4endl
         << "Scope     : foreground_only_axial" << G4endl
         << "Claim     : diagnostic benchmark; not a full CRaTER/LRO model" << G4endl
         << "Source    : " << state.source.modelVersion << G4endl
         << "Interval  : " << state.source.runDates << G4endl
         << "Observation: " << state.observation.column << G4endl
         << "Physics   : Shielding + G4EmStandardPhysics_option4" << G4endl
         << "Cut       : " << state.definition.productionCutUm << " um" << G4endl
         << "Histories : " << state.definition.histories << G4endl
         << "Shard     : " << state.parallel.shardIndex << "/"
         << state.parallel.shardCount << " (" << state.parallel.localHistories
         << " local histories)" << G4endl
         << "Output    : " << state.definition.outputDirectory << G4endl
         << "==========================================" << G4endl;
}

}  // namespace

int main(int argc, char** argv) {
  try {
    if (argc < 2) {
      std::cerr << "usage: shieldSimCraterBenchmark CASE.json "
                   "[--shard-index=N --shard-count=N]\n";
      return 2;
    }
    CraterBenchmarkCase definition = shield_sim::LoadCraterBenchmarkCase(argv[1]);
    CraterBonSpectrum source = shield_sim::LoadCraterBonSpectrum(
        definition.sourcePath, definition.sourceSha256);
    CraterObservation observation = shield_sim::LoadCraterObservation(
        definition.observationPath, definition.observationSha256,
        definition.observationColumn);
    const ParallelConfig parallel = ResolveParallelConfig(definition, argc, argv);
    BenchmarkState state(std::move(definition), std::move(source),
                         std::move(observation), parallel);
    PrintScope(state);

    G4Random::setTheSeed(state.parallel.effectiveSeed);
    auto runManager = std::make_unique<G4RunManager>();
    runManager->SetUserInitialization(new CraterDetectorConstruction(state));
    G4PhysListFactory factory;
    auto* physics = factory.GetReferencePhysList(state.definition.physicsList);
    if (!physics) throw std::runtime_error("Geant4 cannot construct Shielding physics");
    // The 2020 paper names Shielding_EMZ.  In current Geant4 this is assembled
    // explicitly by replacing Shielding's EM constructor with option4.
    physics->ReplacePhysics(new G4EmStandardPhysics_option4());
    physics->SetDefaultCutValue(state.definition.productionCutUm * um);
    physics->SetVerboseLevel(0);
    runManager->SetUserInitialization(physics);
    runManager->SetUserAction(new CraterPrimaryGenerator(state));
    runManager->SetUserAction(new CraterEventAction(state));
    runManager->SetUserAction(new CraterSteppingAction(state));
    runManager->Initialize();

    std::int64_t remaining = state.parallel.localHistories;
    while (remaining > 0) {
      const std::int64_t batch = std::min<std::int64_t>(remaining,
          std::numeric_limits<G4int>::max());
      runManager->BeamOn(static_cast<G4int>(batch));
      remaining -= batch;
    }
    state.WriteResults();
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "CRaTER benchmark error: " << error.what() << '\n';
    return 1;
  }
}
