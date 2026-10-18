#include "Stage7Qualification.hh"

#include "ContentHash.hh"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <tuple>

namespace shield_sim {
namespace {

constexpr double kSymmetryTolerance = 1.0e-12;

double RequiredNumber(const Json::Object& object, const std::string& key,
                      const std::string& context) {
  const auto found = object.find(key);
  if (found == object.end()) throw JsonError(context + " missing '" + key + "'");
  const double value = found->second.AsNumber(context + "." + key);
  if (!std::isfinite(value)) throw JsonError(context + "." + key + " must be finite");
  return value;
}

std::string RequiredString(const Json::Object& object, const std::string& key,
                           const std::string& context) {
  const auto found = object.find(key);
  if (found == object.end()) throw JsonError(context + " missing '" + key + "'");
  const std::string value = found->second.AsString(context + "." + key);
  if (value.empty()) throw JsonError(context + "." + key + " must not be empty");
  return value;
}

std::size_t RequiredPositiveInteger(const Json::Object& object,
                                    const std::string& key,
                                    const std::string& context) {
  const double value = RequiredNumber(object, key, context);
  if (value <= 0.0 || std::floor(value) != value ||
      value > static_cast<double>(std::numeric_limits<std::size_t>::max()))
    throw JsonError(context + "." + key + " must be a positive integer");
  return static_cast<std::size_t>(value);
}

std::vector<double> NumberArray(const Json& value, const std::string& context) {
  std::vector<double> result;
  const auto& array = value.AsArray(context);
  result.reserve(array.size());
  for (std::size_t i = 0; i < array.size(); ++i) {
    const double number = array[i].AsNumber(context + "[" + std::to_string(i) + "]");
    if (!std::isfinite(number)) throw JsonError(context + " entries must be finite");
    result.push_back(number);
  }
  return result;
}

Json DoubleArray(const std::vector<double>& values) {
  Json::Array result;
  result.reserve(values.size());
  for (const double value : values) result.emplace_back(value);
  return Json(std::move(result));
}

bool LeapYear(int year) {
  return (year % 4 == 0 && year % 100 != 0) || year % 400 == 0;
}

int DaysInMonth(int year, int month) {
  static const int days[] = {31,28,31,30,31,30,31,31,30,31,30,31};
  return days[month - 1] + (month == 2 && LeapYear(year) ? 1 : 0);
}

int ParseDigits(const std::string& text, std::size_t begin, std::size_t count,
                const std::string& context) {
  int value = 0;
  for (std::size_t i = begin; i < begin + count; ++i) {
    if (text[i] < '0' || text[i] > '9')
      throw JsonError(context + " must use YYYY-MM-DDTHH:MM:SSZ");
    value = 10 * value + text[i] - '0';
  }
  return value;
}

// Howard Hinnant's civil-date transform, expressed here to avoid locale and
// process-time-zone dependencies in mktime/timegm.
std::int64_t DaysFromCivil(int year, unsigned month, unsigned day) {
  year -= month <= 2;
  const int era = (year >= 0 ? year : year - 399) / 400;
  const unsigned yearOfEra = static_cast<unsigned>(year - era * 400);
  const unsigned shiftedMonth = static_cast<unsigned>(
      static_cast<int>(month) + (month > 2 ? -3 : 9));
  const unsigned dayOfYear = (153 * shiftedMonth + 2) / 5 + day - 1;
  const unsigned dayOfEra = yearOfEra * 365 + yearOfEra / 4 - yearOfEra / 100 + dayOfYear;
  return static_cast<std::int64_t>(era) * 146097 + static_cast<int>(dayOfEra) - 719468;
}

void ValidateStrictlyIncreasingPositive(const std::vector<double>& values,
                                        const std::string& context) {
  if (values.empty()) throw JsonError(context + " must not be empty");
  for (std::size_t i = 0; i < values.size(); ++i) {
    if (!(values[i] > 0.0) || !std::isfinite(values[i]))
      throw JsonError(context + " entries must be finite and positive");
    if (i && !(values[i] > values[i - 1]))
      throw JsonError(context + " must be strictly increasing");
  }
}

double RelativeChange(double current, double previous) {
  const double scale = std::max(std::abs(current), std::abs(previous));
  if (scale == 0.0) return 0.0;
  return std::abs(current - previous) / scale;
}

std::vector<std::vector<double>> ParseSquareMatrix(const Json& value,
                                                    std::size_t expected,
                                                    const std::string& context) {
  const auto& rows = value.AsArray(context);
  if (rows.size() != expected)
    throw JsonError(context + " must have one row per exposure interval");
  std::vector<std::vector<double>> matrix(expected, std::vector<double>(expected));
  for (std::size_t i = 0; i < expected; ++i) {
    const auto parsed = NumberArray(rows[i], context + "[" + std::to_string(i) + "]");
    if (parsed.size() != expected) throw JsonError(context + " must be square");
    matrix[i] = parsed;
  }
  return matrix;
}

void ValidateCovariance(const std::vector<std::vector<double>>& matrix,
                        const std::string& context) {
  const std::size_t n = matrix.size();
  double scale = 1.0;
  for (std::size_t i = 0; i < n; ++i) {
    if (matrix[i][i] < 0.0) throw JsonError(context + " has a negative diagonal");
    scale = std::max(scale, std::abs(matrix[i][i]));
    for (std::size_t j = 0; j < n; ++j) {
      const double tolerance = kSymmetryTolerance *
          std::max({1.0, std::abs(matrix[i][j]), std::abs(matrix[j][i])});
      if (std::abs(matrix[i][j] - matrix[j][i]) > tolerance)
        throw JsonError(context + " must be symmetric");
    }
  }

  // A tolerant LDL^T factorization rejects a materially negative pivot while
  // allowing exact/near-exact rank deficiency common to fully correlated
  // systematic components.
  std::vector<std::vector<double>> lower(n, std::vector<double>(n));
  std::vector<double> diagonal(n, 0.0);
  const double negativeTolerance = 1.0e-12 * scale;
  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t j = 0; j < i; ++j) {
      double residual = matrix[i][j];
      for (std::size_t k = 0; k < j; ++k)
        residual -= lower[i][k] * diagonal[k] * lower[j][k];
      if (std::abs(diagonal[j]) <= negativeTolerance) {
        if (std::abs(residual) > 10.0 * negativeTolerance)
          throw JsonError(context + " is not positive semidefinite");
        lower[i][j] = 0.0;
      } else {
        lower[i][j] = residual / diagonal[j];
      }
    }
    double pivot = matrix[i][i];
    for (std::size_t k = 0; k < i; ++k)
      pivot -= lower[i][k] * lower[i][k] * diagonal[k];
    if (pivot < -negativeTolerance)
      throw JsonError(context + " is not positive semidefinite");
    diagonal[i] = std::max(0.0, pivot);
    lower[i][i] = 1.0;
  }
}

double QuadraticForm(const std::vector<double>& weights,
                     const std::vector<std::vector<double>>& covariance) {
  double value = 0.0;
  for (std::size_t i = 0; i < weights.size(); ++i)
    for (std::size_t j = 0; j < weights.size(); ++j)
      value += weights[i] * covariance[i][j] * weights[j];
  return std::max(0.0, value);
}

struct IntervalRate {
  std::string seriesId;
  std::string quantity;
  std::string rateUnit;
  std::string integratedUnit;
  double value = 0.0;
};

struct ExposureInterval {
  std::string id;
  std::string startUtc;
  std::string endUtc;
  std::int64_t start = 0;
  std::int64_t end = 0;
  std::string originalTimeCoordinate;
  std::string sourceId;
  std::string sourceHash;
  std::string responseId;
  std::string responseHash;
  std::string shieldStackId;
  std::string absorberId;
  std::string attitudeId;
  std::string qualityState;
  double coverageFraction = 0.0;
  std::vector<IntervalRate> rates;
};

struct SeriesDefinition {
  std::string id;
  std::string quantity;
  std::string rateUnit;
  std::string integratedUnit;
};

Json NullOrNumber(bool available, double value) {
  return available ? Json(value) : Json(nullptr);
}

} // namespace

double InterpolateResponseInLogEnergy(const std::vector<double>& energyMeV,
                                      const std::vector<double>& response,
                                      double queryEnergyMeV) {
  if (energyMeV.size() != response.size() || energyMeV.empty())
    throw std::invalid_argument("response energy/value arrays must be nonempty and equal-sized");
  ValidateStrictlyIncreasingPositive(energyMeV, "response_energy_nodes_MeV");
  if (!(queryEnergyMeV > 0.0) || !std::isfinite(queryEnergyMeV))
    throw std::invalid_argument("query energy must be finite and positive");
  for (const double value : response)
    if (!std::isfinite(value)) throw std::invalid_argument("response values must be finite");
  const double tolerance = 1.0e-12;
  if (queryEnergyMeV < energyMeV.front() * (1.0 - tolerance) ||
      queryEnergyMeV > energyMeV.back() * (1.0 + tolerance))
    throw std::out_of_range("response query requires prohibited energy extrapolation");
  if (queryEnergyMeV <= energyMeV.front()) return response.front();
  if (queryEnergyMeV >= energyMeV.back()) return response.back();
  const auto upper = std::upper_bound(energyMeV.begin(), energyMeV.end(), queryEnergyMeV);
  const std::size_t high = static_cast<std::size_t>(upper - energyMeV.begin());
  const std::size_t low = high - 1;
  const double fraction = (std::log(queryEnergyMeV) - std::log(energyMeV[low])) /
                          (std::log(energyMeV[high]) - std::log(energyMeV[low]));
  return response[low] + fraction * (response[high] - response[low]);
}

AdaptiveResponseFoldResult FoldPiecewiseConstantSpectrum(
    const std::vector<double>& sourceEdgesMeV,
    const std::vector<double>& differentialSource,
    const std::vector<double>& responseEnergyMeV,
    const std::vector<double>& response,
    double relativeTolerance,
    std::size_t maximumRefinements) {
  if (sourceEdgesMeV.size() != differentialSource.size() + 1 || differentialSource.empty())
    throw std::invalid_argument("source edges must contain exactly one more value than source bins");
  ValidateStrictlyIncreasingPositive(sourceEdgesMeV, "source_bin_edges_MeV");
  ValidateStrictlyIncreasingPositive(responseEnergyMeV, "response_energy_nodes_MeV");
  if (responseEnergyMeV.size() != response.size())
    throw std::invalid_argument("response arrays must have equal size");
  if (sourceEdgesMeV.front() < responseEnergyMeV.front() ||
      sourceEdgesMeV.back() > responseEnergyMeV.back())
    throw std::out_of_range("source support requires prohibited response extrapolation");
  for (const double density : differentialSource)
    if (!std::isfinite(density) || density < 0.0)
      throw std::invalid_argument("differential source values must be finite and nonnegative");
  if (!(relativeTolerance > 0.0) || !std::isfinite(relativeTolerance) ||
      relativeTolerance > 0.02)
    throw std::invalid_argument("Stage-7 relative quadrature tolerance must be in (0,0.02]");
  if (maximumRefinements == 0 || maximumRefinements > 30)
    throw std::invalid_argument("maximum refinements must be in [1,30]");

  AdaptiveResponseFoldResult result;
  std::size_t subdivisions = 1;
  double previous = 0.0;
  bool havePrevious = false;
  for (std::size_t refinement = 0; refinement <= maximumRefinements; ++refinement) {
    double folded = 0.0;
    std::size_t evaluations = 0;
    for (std::size_t bin = 0; bin < differentialSource.size(); ++bin) {
      const double width = (sourceEdgesMeV[bin + 1] - sourceEdgesMeV[bin]) /
                           static_cast<double>(subdivisions);
      for (std::size_t part = 0; part < subdivisions; ++part) {
        const double energy = sourceEdgesMeV[bin] +
            (static_cast<double>(part) + 0.5) * width;
        folded += differentialSource[bin] * width *
                  InterpolateResponseInLogEnergy(responseEnergyMeV, response, energy);
        ++evaluations;
      }
    }
    result.value = folded;
    result.previousValue = havePrevious ? previous : folded;
    result.finalRelativeChange = havePrevious ? RelativeChange(folded, previous) : 0.0;
    result.subdivisionsPerSourceBin = subdivisions;
    result.refinements = refinement;
    result.responseEvaluations = evaluations;
    if (havePrevious && result.finalRelativeChange <= relativeTolerance) {
      result.converged = true;
      return result;
    }
    havePrevious = true;
    previous = folded;
    subdivisions *= 2;
  }
  return result;
}

std::vector<NeymanAllocation> AllocateNeymanHistories(
    const std::vector<NeymanStratum>& strata,
    std::size_t totalHistories,
    std::size_t minimumPerStratum) {
  if (strata.empty()) throw std::invalid_argument("Neyman allocation requires at least one stratum");
  if (minimumPerStratum == 0)
    throw std::invalid_argument("minimum histories per stratum must be positive");
  if (minimumPerStratum > totalHistories / strata.size())
    throw std::invalid_argument("total histories cannot satisfy the per-stratum minimum");
  std::set<std::string> ids;
  std::vector<double> importance;
  importance.reserve(strata.size());
  double totalImportance = 0.0;
  for (const auto& stratum : strata) {
    if (stratum.id.empty() || !ids.insert(stratum.id).second)
      throw std::invalid_argument("stratum IDs must be nonempty and unique");
    if (!std::isfinite(stratum.sourceWeight) || stratum.sourceWeight < 0.0 ||
        !std::isfinite(stratum.pilotStandardDeviation) || stratum.pilotStandardDeviation < 0.0)
      throw std::invalid_argument("Neyman weights and pilot deviations must be finite and nonnegative");
    importance.push_back(std::abs(stratum.sourceWeight) * stratum.pilotStandardDeviation);
    totalImportance += importance.back();
  }

  const std::size_t remaining = totalHistories - minimumPerStratum * strata.size();
  std::vector<NeymanAllocation> result;
  result.reserve(strata.size());
  std::vector<std::pair<double,std::size_t>> remainders;
  std::size_t assigned = minimumPerStratum * strata.size();
  for (std::size_t i = 0; i < strata.size(); ++i) {
    const double fraction = totalImportance > 0.0
        ? importance[i] / totalImportance : 1.0 / static_cast<double>(strata.size());
    const double exactExtra = fraction * static_cast<double>(remaining);
    const std::size_t extra = static_cast<std::size_t>(std::floor(exactExtra));
    result.push_back({strata[i].id, minimumPerStratum + extra});
    assigned += extra;
    remainders.emplace_back(exactExtra - static_cast<double>(extra), i);
  }
  std::sort(remainders.begin(), remainders.end(), [&](const auto& first, const auto& second) {
    if (first.first != second.first) return first.first > second.first;
    return strata[first.second].id < strata[second.second].id;
  });
  for (std::size_t i = 0; assigned < totalHistories; ++i, ++assigned)
    ++result[remainders[i].second].histories;
  return result;
}

std::int64_t ParseUtcTimestamp(const std::string& utc) {
  const std::string context = "UTC timestamp '" + utc + "'";
  if (utc.size() != 20 || utc[4] != '-' || utc[7] != '-' || utc[10] != 'T' ||
      utc[13] != ':' || utc[16] != ':' || utc[19] != 'Z')
    throw JsonError(context + " must use YYYY-MM-DDTHH:MM:SSZ");
  const int year = ParseDigits(utc, 0, 4, context);
  const int month = ParseDigits(utc, 5, 2, context);
  const int day = ParseDigits(utc, 8, 2, context);
  const int hour = ParseDigits(utc, 11, 2, context);
  const int minute = ParseDigits(utc, 14, 2, context);
  const int second = ParseDigits(utc, 17, 2, context);
  if (year < 1970 || month < 1 || month > 12 || day < 1 ||
      day > DaysInMonth(year, month) || hour > 23 || minute > 59 || second > 59)
    throw JsonError(context + " is not a supported Gregorian UTC instant");
  return 86400 * DaysFromCivil(year, static_cast<unsigned>(month), static_cast<unsigned>(day)) +
         3600 * hour + 60 * minute + second;
}

Json FoldSpectrumDocument(const Json& input) {
  const auto& root = input.AsObject("response quadrature");
  RejectUnknownFields(root, {"schema_version","source_bin_edges_MeV",
      "differential_source","response_energy_nodes_MeV","response_values",
      "source_density_unit","response_unit","output_unit","relative_tolerance",
      "maximum_refinements","extrapolation_policy"}, "response quadrature");
  if (RequiredString(root, "schema_version", "response quadrature") !=
      "shieldSim.response_quadrature.v1")
    throw JsonError("response quadrature schema_version must be shieldSim.response_quadrature.v1");
  if (RequiredString(root, "extrapolation_policy", "response quadrature") != "refuse")
    throw JsonError("Stage-7 response extrapolation_policy must be 'refuse'");
  const auto sourceEdges = NumberArray(root.at("source_bin_edges_MeV"),
                                       "response quadrature.source_bin_edges_MeV");
  const auto source = NumberArray(root.at("differential_source"),
                                  "response quadrature.differential_source");
  const auto responseEnergy = NumberArray(root.at("response_energy_nodes_MeV"),
                                          "response quadrature.response_energy_nodes_MeV");
  const auto response = NumberArray(root.at("response_values"),
                                    "response quadrature.response_values");
  const double tolerance = RequiredNumber(root, "relative_tolerance", "response quadrature");
  const std::size_t maximum = RequiredPositiveInteger(root, "maximum_refinements",
                                                       "response quadrature");
  const auto folded = FoldPiecewiseConstantSpectrum(sourceEdges, source, responseEnergy,
                                                     response, tolerance, maximum);
  Json::Object result{
    {"converged", folded.converged},
    {"extrapolation_used", false},
    {"final_relative_change", folded.finalRelativeChange},
    {"input_hash", Sha256(input.DumpCanonical())},
    {"maximum_relative_tolerance", tolerance},
    {"output_unit", RequiredString(root, "output_unit", "response quadrature")},
    {"previous_value", folded.previousValue},
    {"refinements", static_cast<double>(folded.refinements)},
    {"response_evaluations", static_cast<double>(folded.responseEvaluations)},
    {"schema_version", "shieldSim.response_quadrature_result.v1"},
    {"subdivisions_per_source_bin", static_cast<double>(folded.subdivisionsPerSourceBin)},
    {"value", folded.value},
    {"validation_status", folded.converged
        ? "stage7_quadrature_converged_no_extrapolation"
        : "not_assessable_quadrature_did_not_converge"}
  };
  return Json(std::move(result));
}

Json AllocateHistoriesDocument(const Json& input) {
  const auto& root = input.AsObject("event allocation");
  RejectUnknownFields(root, {"schema_version","policy","total_histories",
      "minimum_per_stratum","strata"}, "event allocation");
  if (RequiredString(root, "schema_version", "event allocation") !=
      "shieldSim.neyman_allocation.v1")
    throw JsonError("event allocation schema_version must be shieldSim.neyman_allocation.v1");
  if (RequiredString(root, "policy", "event allocation") !=
      "pilot_neyman_no_transport_biasing")
    throw JsonError("event allocation policy must be pilot_neyman_no_transport_biasing");
  const std::size_t total = RequiredPositiveInteger(root, "total_histories", "event allocation");
  const std::size_t minimum = RequiredPositiveInteger(root, "minimum_per_stratum",
                                                       "event allocation");
  std::vector<NeymanStratum> strata;
  for (const auto& item : root.at("strata").AsArray("event allocation.strata")) {
    const auto& object = item.AsObject("event allocation.strata[]");
    RejectUnknownFields(object, {"stratum_id","source_weight","pilot_standard_deviation"},
                        "event allocation.strata[]");
    strata.push_back({RequiredString(object, "stratum_id", "event allocation.strata[]"),
                      RequiredNumber(object, "source_weight", "event allocation.strata[]"),
                      RequiredNumber(object, "pilot_standard_deviation",
                                     "event allocation.strata[]")});
  }
  const auto allocation = AllocateNeymanHistories(strata, total, minimum);
  Json::Array rows;
  for (std::size_t i = 0; i < allocation.size(); ++i)
    rows.emplace_back(Json::Object{
      {"histories", static_cast<double>(allocation[i].histories)},
      {"pilot_standard_deviation", strata[i].pilotStandardDeviation},
      {"source_weight", strata[i].sourceWeight},
      {"stratum_id", allocation[i].id}});
  return Json(Json::Object{
    {"allocation", Json(std::move(rows))},
    {"input_hash", Sha256(input.DumpCanonical())},
    {"minimum_per_stratum", static_cast<double>(minimum)},
    {"physics_biasing", false},
    {"policy", "pilot_neyman_no_transport_biasing"},
    {"schema_version", "shieldSim.neyman_allocation_result.v1"},
    {"total_histories", static_cast<double>(total)}});
}

Json IntegrateTimeResolvedExposure(const Json& input) {
  const auto& root = input.AsObject("time-resolved exposure");
  RejectUnknownFields(root, {"schema_version","case_id","time_system",
      "source_qualification","integration_policy","anomaly_timestamps","intervals",
      "uncertainty"}, "time-resolved exposure");
  if (RequiredString(root, "schema_version", "time-resolved exposure") !=
      "shieldSim.time_resolved_exposure.v1")
    throw JsonError("time-resolved exposure schema_version must be shieldSim.time_resolved_exposure.v1");
  const std::string caseId = RequiredString(root, "case_id", "time-resolved exposure");
  if (RequiredString(root, "time_system", "time-resolved exposure") != "UTC")
    throw JsonError("time-resolved exposure time_system must be UTC");

  const auto& qualification = root.at("source_qualification").AsObject(
      "time-resolved exposure.source_qualification");
  RejectUnknownFields(qualification, {"environment","environment_completeness",
      "goes_transport_status","transport_admission","limitations"},
      "time-resolved exposure.source_qualification");
  const std::string environment = RequiredString(qualification, "environment",
      "time-resolved exposure.source_qualification");
  const std::string completeness = RequiredString(qualification, "environment_completeness",
      "time-resolved exposure.source_qualification");
  const std::string goesStatus = RequiredString(qualification, "goes_transport_status",
      "time-resolved exposure.source_qualification");
  const std::string admission = RequiredString(qualification, "transport_admission",
      "time-resolved exposure.source_qualification");
  const std::set<std::string> allowedEnvironments{
      "ams02_a1_proton_component", "goes16_sgps_directional_ingestion",
      "synthetic_verification"};
  if (!allowedEnvironments.count(environment))
    throw JsonError("Stage-7 source environment is not admitted by the approved qualification boundary");
  if (environment == "ams02_a1_proton_component" &&
      completeness != "incomplete_not_full_gcr")
    throw JsonError("AMS-02 A=1 qualification input must be labeled incomplete_not_full_gcr");
  if (environment == "goes16_sgps_directional_ingestion" &&
      (goesStatus != "disabled_missing_boundary_mapping" || admission != "disabled"))
    throw JsonError("GOES transport must remain disabled until boundary mapping is supplied");
  const auto& limitations = qualification.at("limitations").AsArray(
      "time-resolved exposure.source_qualification.limitations");
  if (limitations.empty()) throw JsonError("source qualification must state at least one limitation");
  for (std::size_t i = 0; i < limitations.size(); ++i)
    if (limitations[i].AsString("source qualification limitation").empty())
      throw JsonError("source qualification limitations must not be empty");

  const auto& policy = root.at("integration_policy").AsObject(
      "time-resolved exposure.integration_policy");
  RejectUnknownFields(policy, {"gap_policy","cadence_policy","configuration_policy",
      "extrapolation_policy"}, "time-resolved exposure.integration_policy");
  if (RequiredString(policy, "gap_policy", "integration policy") != "preserve_no_fill" ||
      RequiredString(policy, "cadence_policy", "integration policy") !=
          "native_intervals_no_implicit_rebin" ||
      RequiredString(policy, "configuration_policy", "integration policy") !=
          "split_before_fold" ||
      RequiredString(policy, "extrapolation_policy", "integration policy") != "refuse")
    throw JsonError("time-resolved exposure uses a non-approved integration policy");

  Json::Array anomalyOutput;
  std::set<std::string> anomalyIds;
  for (const auto& item : root.at("anomaly_timestamps").AsArray(
           "time-resolved exposure.anomaly_timestamps")) {
    const auto& object = item.AsObject("time-resolved exposure.anomaly_timestamps[]");
    RejectUnknownFields(object, {"timestamp_id","utc"},
                        "time-resolved exposure.anomaly_timestamps[]");
    const std::string id = RequiredString(object, "timestamp_id", "anomaly timestamp");
    const std::string utc = RequiredString(object, "utc", "anomaly timestamp");
    if (!anomalyIds.insert(id).second) throw JsonError("duplicate anomaly timestamp_id '" + id + "'");
    const std::int64_t seconds = ParseUtcTimestamp(utc);
    anomalyOutput.emplace_back(Json::Object{{"timestamp_id", id}, {"unix_seconds", static_cast<double>(seconds)},
                                            {"utc", utc}});
  }

  std::vector<ExposureInterval> intervals;
  std::map<std::string,SeriesDefinition> series;
  std::set<std::string> intervalIds;
  const auto& intervalArray = root.at("intervals").AsArray("time-resolved exposure.intervals");
  if (intervalArray.empty()) throw JsonError("time-resolved exposure intervals must not be empty");
  for (std::size_t index = 0; index < intervalArray.size(); ++index) {
    const std::string context = "time-resolved exposure.intervals[" + std::to_string(index) + "]";
    const auto& object = intervalArray[index].AsObject(context);
    RejectUnknownFields(object, {"interval_id","start_utc","end_utc",
        "original_time_coordinate","source_id","source_hash","response_id","response_hash",
        "shield_stack_id","absorber_id","attitude_id","quality_state","coverage_fraction",
        "rates"}, context);
    ExposureInterval parsed;
    parsed.id = RequiredString(object, "interval_id", context);
    if (!intervalIds.insert(parsed.id).second) throw JsonError("duplicate interval_id '" + parsed.id + "'");
    parsed.startUtc = RequiredString(object, "start_utc", context);
    parsed.endUtc = RequiredString(object, "end_utc", context);
    parsed.start = ParseUtcTimestamp(parsed.startUtc);
    parsed.end = ParseUtcTimestamp(parsed.endUtc);
    if (parsed.end <= parsed.start) throw JsonError(context + " end must be after start");
    if (!intervals.empty() && parsed.start != intervals.back().end)
      throw JsonError(context + " must begin at the preceding end; gaps require explicit gap intervals");
    parsed.originalTimeCoordinate = RequiredString(object, "original_time_coordinate", context);
    parsed.sourceId = RequiredString(object, "source_id", context);
    parsed.sourceHash = RequiredString(object, "source_hash", context);
    parsed.responseId = RequiredString(object, "response_id", context);
    parsed.responseHash = RequiredString(object, "response_hash", context);
    parsed.shieldStackId = RequiredString(object, "shield_stack_id", context);
    parsed.absorberId = RequiredString(object, "absorber_id", context);
    parsed.attitudeId = RequiredString(object, "attitude_id", context);
    parsed.qualityState = RequiredString(object, "quality_state", context);
    parsed.coverageFraction = RequiredNumber(object, "coverage_fraction", context);
    if (parsed.coverageFraction < 0.0 || parsed.coverageFraction > 1.0)
      throw JsonError(context + ".coverage_fraction must be in [0,1]");
    if (parsed.qualityState == "valid") {
      if (parsed.coverageFraction != 1.0)
        throw JsonError(context + " valid interval requires coverage_fraction=1");
    } else if (parsed.qualityState == "partial") {
      if (!(parsed.coverageFraction > 0.0 && parsed.coverageFraction < 1.0))
        throw JsonError(context + " partial interval requires coverage strictly between zero and one");
    } else if (parsed.qualityState == "gap") {
      if (parsed.coverageFraction != 0.0)
        throw JsonError(context + " gap interval requires coverage_fraction=0");
    } else {
      throw JsonError(context + ".quality_state must be valid, partial, or gap");
    }
    const auto& rates = object.at("rates").AsArray(context + ".rates");
    if (parsed.qualityState == "gap" && !rates.empty())
      throw JsonError(context + " gap interval must not contain calculated rates");
    if (parsed.qualityState != "gap" && rates.empty())
      throw JsonError(context + " non-gap interval must contain calculated rates");
    std::set<std::string> intervalSeries;
    for (const auto& rateValue : rates) {
      const auto& rate = rateValue.AsObject(context + ".rates[]");
      RejectUnknownFields(rate, {"series_id","quantity","rate_unit","integrated_unit","value"},
                          context + ".rates[]");
      IntervalRate parsedRate;
      parsedRate.seriesId = RequiredString(rate, "series_id", context + ".rates[]");
      parsedRate.quantity = RequiredString(rate, "quantity", context + ".rates[]");
      parsedRate.rateUnit = RequiredString(rate, "rate_unit", context + ".rates[]");
      parsedRate.integratedUnit = RequiredString(rate, "integrated_unit", context + ".rates[]");
      parsedRate.value = RequiredNumber(rate, "value", context + ".rates[]");
      if (!intervalSeries.insert(parsedRate.seriesId).second)
        throw JsonError(context + " contains duplicate series_id '" + parsedRate.seriesId + "'");
      const SeriesDefinition definition{parsedRate.seriesId, parsedRate.quantity,
                                        parsedRate.rateUnit, parsedRate.integratedUnit};
      const auto existing = series.find(parsedRate.seriesId);
      if (existing == series.end()) series.emplace(parsedRate.seriesId, definition);
      else if (std::tie(existing->second.quantity, existing->second.rateUnit,
                        existing->second.integratedUnit) !=
               std::tie(definition.quantity, definition.rateUnit, definition.integratedUnit))
        throw JsonError("series_id '" + parsedRate.seriesId + "' changes quantity or units");
      parsed.rates.push_back(std::move(parsedRate));
    }
    intervals.push_back(std::move(parsed));
  }

  const std::size_t count = intervals.size();
  std::map<std::string,std::vector<double>> ratesBySeries;
  std::vector<double> timeWeights(count, 0.0);
  for (const auto& item : series) ratesBySeries[item.first] = std::vector<double>(count, 0.0);
  Json::Array intervalOutput;
  std::set<std::string> configurationSegments;
  for (std::size_t i = 0; i < count; ++i) {
    const auto& interval = intervals[i];
    const double duration = static_cast<double>(interval.end - interval.start);
    const double validDuration = duration * interval.coverageFraction;
    timeWeights[i] = validDuration;
    Json::Array contributions;
    for (const auto& rate : interval.rates) {
      ratesBySeries[rate.seriesId][i] = rate.value;
      contributions.emplace_back(Json::Object{
        {"integrated_unit", rate.integratedUnit},
        {"interval_exposure", rate.value * validDuration},
        {"quantity", rate.quantity}, {"rate", rate.value},
        {"rate_unit", rate.rateUnit}, {"series_id", rate.seriesId}});
    }
    if (interval.qualityState != "gap") {
      const std::string segment = interval.shieldStackId + "|" + interval.absorberId + "|" +
                                  interval.attitudeId + "|" + interval.responseHash;
      configurationSegments.insert(segment);
    }
    intervalOutput.emplace_back(Json::Object{
      {"absorber_id", interval.absorberId}, {"attitude_id", interval.attitudeId},
      {"coverage_fraction", interval.coverageFraction},
      {"duration_s", duration}, {"end_utc", interval.endUtc},
      {"interval_exposures", Json(std::move(contributions))}, {"interval_id", interval.id},
      {"original_time_coordinate", interval.originalTimeCoordinate},
      {"quality_state", interval.qualityState}, {"response_hash", interval.responseHash},
      {"response_id", interval.responseId}, {"shield_stack_id", interval.shieldStackId},
      {"source_hash", interval.sourceHash}, {"source_id", interval.sourceId},
      {"start_utc", interval.startUtc}, {"valid_duration_s", validDuration}});
  }

  const auto& uncertainty = root.at("uncertainty").AsObject("time-resolved exposure.uncertainty");
  RejectUnknownFields(uncertainty, {"status","limitation","independence_basis",
      "covariance_components","sampled_ensembles"}, "time-resolved exposure.uncertainty");
  const std::string uncertaintyStatus = RequiredString(uncertainty, "status", "uncertainty");
  if (uncertaintyStatus != "complete" && uncertaintyStatus != "incomplete")
    throw JsonError("uncertainty.status must be complete or incomplete");
  const std::string limitation = RequiredString(uncertainty, "limitation", "uncertainty");
  const std::string independence = RequiredString(uncertainty, "independence_basis", "uncertainty");

  struct ComponentResult { std::string id, category, provenance; double variance = 0.0; };
  std::map<std::string,std::vector<ComponentResult>> componentResults;
  std::set<std::string> componentIds;
  for (const auto& item : uncertainty.at("covariance_components").AsArray(
           "uncertainty.covariance_components")) {
    const auto& object = item.AsObject("uncertainty.covariance_components[]");
    RejectUnknownFields(object, {"component_id","category","series_id","provenance",
        "rate_unit_squared","interval_rate_covariance"}, "uncertainty.covariance_components[]");
    const std::string id = RequiredString(object, "component_id", "covariance component");
    if (!componentIds.insert(id).second) throw JsonError("duplicate covariance component_id '" + id + "'");
    const std::string category = RequiredString(object, "category", "covariance component");
    if (category != "source" && category != "response_matrix_statistical" &&
        category != "numerical" && category != "physics_reference")
      throw JsonError("unsupported Stage-7 covariance component category '" + category + "'");
    const std::string seriesId = RequiredString(object, "series_id", "covariance component");
    if (!series.count(seriesId)) throw JsonError("covariance component references unknown series_id");
    // Parentheses are part of the contract: (Gy/s)^2 is a rate covariance,
    // whereas the unparenthesized Gy/s^2 denotes a physically different unit.
    const std::string expectedUnit = "(" + series.at(seriesId).rateUnit + ")^2";
    if (RequiredString(object, "rate_unit_squared", "covariance component") != expectedUnit)
      throw JsonError("covariance component rate_unit_squared must be '" + expectedUnit + "'");
    const auto matrix = ParseSquareMatrix(object.at("interval_rate_covariance"), count,
                                          "covariance component.interval_rate_covariance");
    ValidateCovariance(matrix, "covariance component.interval_rate_covariance");
    componentResults[seriesId].push_back({id, category,
      RequiredString(object, "provenance", "covariance component"),
      QuadraticForm(timeWeights, matrix)});
  }

  struct EnsembleResult { std::string id, method, provenance; double mean = 0.0, sd = 0.0; std::size_t n = 0; };
  std::map<std::string,std::vector<EnsembleResult>> ensembleResults;
  std::set<std::string> ensembleIds;
  for (const auto& item : uncertainty.at("sampled_ensembles").AsArray(
           "uncertainty.sampled_ensembles")) {
    const auto& object = item.AsObject("uncertainty.sampled_ensembles[]");
    RejectUnknownFields(object, {"ensemble_id","category","series_id","method","provenance",
        "rate_samples"}, "uncertainty.sampled_ensembles[]");
    const std::string id = RequiredString(object, "ensemble_id", "sampled ensemble");
    if (!ensembleIds.insert(id).second) throw JsonError("duplicate ensemble_id '" + id + "'");
    if (RequiredString(object, "category", "sampled ensemble") != "source")
      throw JsonError("Stage-7 sampled ensembles currently represent source uncertainty only");
    const std::string seriesId = RequiredString(object, "series_id", "sampled ensemble");
    if (!series.count(seriesId)) throw JsonError("sampled ensemble references unknown series_id");
    const auto& samples = object.at("rate_samples").AsArray("sampled ensemble.rate_samples");
    if (samples.size() < 2) throw JsonError("sampled ensemble requires at least two samples");
    std::vector<double> exposures;
    exposures.reserve(samples.size());
    for (std::size_t sample = 0; sample < samples.size(); ++sample) {
      const auto row = NumberArray(samples[sample], "sampled ensemble.rate_samples[]");
      if (row.size() != count) throw JsonError("each sampled ensemble row must match interval count");
      double exposure = 0.0;
      for (std::size_t i = 0; i < count; ++i) exposure += timeWeights[i] * row[i];
      exposures.push_back(exposure);
    }
    double mean = 0.0;
    for (const double value : exposures) mean += value;
    mean /= static_cast<double>(exposures.size());
    double sumSquares = 0.0;
    for (const double value : exposures) sumSquares += (value - mean) * (value - mean);
    ensembleResults[seriesId].push_back({id,
      RequiredString(object, "method", "sampled ensemble"),
      RequiredString(object, "provenance", "sampled ensemble"), mean,
      std::sqrt(sumSquares / static_cast<double>(exposures.size() - 1)), exposures.size()});
  }
  if (uncertaintyStatus == "complete") {
    if (componentResults.empty() || ensembleResults.empty())
      throw JsonError("complete Stage-7 uncertainty requires covariance and sampled-ensemble paths");
    for (const auto& item : series)
      if (!componentResults.count(item.first) || !ensembleResults.count(item.first))
        throw JsonError("complete Stage-7 uncertainty must cover every series");
  }

  Json::Array integrated;
  for (const auto& item : series) {
    const auto& definition = item.second;
    double value = 0.0;
    for (std::size_t i = 0; i < count; ++i)
      value += timeWeights[i] * ratesBySeries[definition.id][i];
    double variance = 0.0;
    Json::Array components;
    for (const auto& component : componentResults[definition.id]) {
      variance += component.variance;
      components.emplace_back(Json::Object{
        {"category", component.category}, {"component_id", component.id},
        {"integrated_standard_uncertainty", std::sqrt(component.variance)},
        {"integrated_variance", component.variance}, {"provenance", component.provenance}});
    }
    Json::Array ensembles;
    for (const auto& ensemble : ensembleResults[definition.id])
      ensembles.emplace_back(Json::Object{
        {"ensemble_id", ensemble.id}, {"integrated_mean", ensemble.mean},
        {"integrated_standard_deviation", ensemble.sd}, {"method", ensemble.method},
        {"provenance", ensemble.provenance}, {"samples", static_cast<double>(ensemble.n)}});
    const bool quantified = !componentResults[definition.id].empty();
    integrated.emplace_back(Json::Object{
      {"accuracy_status", uncertaintyStatus == "complete" ? "assessable" : "not_assessable"},
      {"covariance_components", Json(std::move(components))},
      {"integrated_unit", definition.integratedUnit},
      {"quantity", definition.quantity}, {"sampled_ensembles", Json(std::move(ensembles))},
      {"series_id", definition.id},
      {"standard_uncertainty", NullOrNumber(quantified, std::sqrt(variance))},
      {"uncertainty_status", uncertaintyStatus}, {"value", value}});
  }

  Json::Array segments;
  for (const auto& segment : configurationSegments) segments.emplace_back(segment);
  return Json(Json::Object{
    {"anomaly_timestamps", Json(std::move(anomalyOutput))}, {"case_id", caseId},
    {"configuration_segments", Json(std::move(segments))},
    {"input_hash", Sha256(input.DumpCanonical())},
    {"integrated_exposures", Json(std::move(integrated))},
    {"integration_policy", root.at("integration_policy")},
    {"interval_results", Json(std::move(intervalOutput))},
    {"schema_version", "shieldSim.time_resolved_exposure_result.v1"},
    {"source_qualification", root.at("source_qualification")}, {"time_system", "UTC"},
    {"uncertainty", Json(Json::Object{{"independence_basis", independence},
        {"limitation", limitation}, {"status", uncertaintyStatus}})},
    {"validation_status", "stage7_software_framework_incomplete_physical_gcr_goes_qualification"}});
}

} // namespace shield_sim
