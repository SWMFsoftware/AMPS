#include "Stage7Qualification.hh"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace fs = std::filesystem;

namespace {
int failures = 0;

void Check(bool condition, const std::string& message) {
  if (!condition) {
    ++failures;
    std::cerr << "FAIL: " << message << '\n';
  }
}

void CheckNear(double actual, double expected, double relativeTolerance,
               const std::string& message) {
  const double scale = std::max({1.0, std::abs(actual), std::abs(expected)});
  Check(std::abs(actual - expected) <= relativeTolerance * scale,
        message + " actual=" + std::to_string(actual) +
        " expected=" + std::to_string(expected));
}

void CheckThrows(const std::function<void()>& action, const std::string& message) {
  try {
    action();
    Check(false, message + " did not throw");
  } catch (const std::exception&) {
  }
}

shield_sim::Json ValidExposure() {
  return shield_sim::Json::Parse(R"json(
{
  "schema_version":"shieldSim.time_resolved_exposure.v1",
  "case_id":"stage7_piecewise_verification",
  "time_system":"UTC",
  "source_qualification":{
    "environment":"synthetic_verification",
    "environment_completeness":"synthetic_only",
    "goes_transport_status":"disabled_missing_boundary_mapping",
    "transport_admission":"verification_test_only",
    "limitations":["synthetic rates test integration and uncertainty algebra only"]
  },
  "integration_policy":{
    "gap_policy":"preserve_no_fill",
    "cadence_policy":"native_intervals_no_implicit_rebin",
    "configuration_policy":"split_before_fold",
    "extrapolation_policy":"refuse"
  },
  "anomaly_timestamps":[{"timestamp_id":"telemetry_marker","utc":"2020-02-29T00:00:17Z"}],
  "intervals":[
    {
      "interval_id":"i0","start_utc":"2020-02-29T00:00:00Z","end_utc":"2020-02-29T00:00:10Z",
      "original_time_coordinate":"synthetic_start_0","source_id":"source-a","source_hash":"sha256:source-a",
      "response_id":"response-a","response_hash":"sha256:response-a","shield_stack_id":"al-hdpe",
      "absorber_id":"silicon","attitude_id":"fixed-a","quality_state":"valid","coverage_fraction":1.0,
      "rates":[{"series_id":"silicon/TID","quantity":"TID","rate_unit":"proxy/s","integrated_unit":"proxy","value":1.0}]
    },
    {
      "interval_id":"i1","start_utc":"2020-02-29T00:00:10Z","end_utc":"2020-02-29T00:00:25Z",
      "original_time_coordinate":"synthetic_start_10","source_id":"source-a","source_hash":"sha256:source-a",
      "response_id":"response-a","response_hash":"sha256:response-a","shield_stack_id":"al-hdpe",
      "absorber_id":"silicon","attitude_id":"fixed-a","quality_state":"partial","coverage_fraction":0.5,
      "rates":[{"series_id":"silicon/TID","quantity":"TID","rate_unit":"proxy/s","integrated_unit":"proxy","value":2.0}]
    },
    {
      "interval_id":"gap","start_utc":"2020-02-29T00:00:25Z","end_utc":"2020-02-29T00:00:35Z",
      "original_time_coordinate":"explicit_missing_support","source_id":"source-gap","source_hash":"sha256:gap",
      "response_id":"none","response_hash":"sha256:none","shield_stack_id":"unknown-during-gap",
      "absorber_id":"silicon","attitude_id":"unknown-during-gap","quality_state":"gap","coverage_fraction":0.0,
      "rates":[]
    },
    {
      "interval_id":"i3","start_utc":"2020-02-29T00:00:35Z","end_utc":"2020-02-29T00:00:55Z",
      "original_time_coordinate":"synthetic_start_35","source_id":"source-b","source_hash":"sha256:source-b",
      "response_id":"response-b","response_hash":"sha256:response-b","shield_stack_id":"hdpe-al",
      "absorber_id":"silicon","attitude_id":"fixed-b","quality_state":"valid","coverage_fraction":1.0,
      "rates":[{"series_id":"silicon/TID","quantity":"TID","rate_unit":"proxy/s","integrated_unit":"proxy","value":3.0}]
    }
  ],
  "uncertainty":{
    "status":"complete",
    "limitation":"synthetic covariance and ensemble are exact test fixtures",
    "independence_basis":"source and response test estimators are generated independently",
    "covariance_components":[
      {
        "component_id":"source_covariance","category":"source","series_id":"silicon/TID",
        "provenance":"analytic synthetic diagonal covariance","rate_unit_squared":"(proxy/s)^2",
        "interval_rate_covariance":[[0.01,0,0,0],[0,0.04,0,0],[0,0,0,0],[0,0,0,0.09]]
      },
      {
        "component_id":"response_mc_covariance","category":"response_matrix_statistical","series_id":"silicon/TID",
        "provenance":"independent synthetic response estimates","rate_unit_squared":"(proxy/s)^2",
        "interval_rate_covariance":[[0.0025,0,0,0],[0,0.01,0,0],[0,0,0,0],[0,0,0,0.0225]]
      }
    ],
    "sampled_ensembles":[
      {
        "ensemble_id":"source_symmetric_three_sample","category":"source","series_id":"silicon/TID",
        "method":"preregistered symmetric verification ensemble","provenance":"analytic synthetic samples",
        "rate_samples":[[0.9,1.8,0,2.7],[1.0,2.0,0,3.0],[1.1,2.2,0,3.3]]
      }
    ]
  }
}
)json", "stage7-valid-exposure");
}

} // namespace

int main(int argc, char** argv) {
  if (argc != 2) {
    std::cerr << "usage: stage7_core_tests EVIDENCE_DIR\n";
    return 2;
  }
  const fs::path evidence = argv[1];
  fs::create_directories(evidence);

  // The interpolation oracle is exactly affine in log(E), so the midpoint at
  // E=2 must be reproduced without a tolerance associated with transport.
  const std::vector<double> responseEnergy{1.0, 2.0, 4.0};
  const std::vector<double> response{
      3.0, 3.0 + 0.5 * std::log(2.0), 3.0 + 0.5 * std::log(4.0)};
  CheckNear(shield_sim::InterpolateResponseInLogEnergy(responseEnergy, response, std::sqrt(2.0)),
            3.0 + 0.25 * std::log(2.0), 1.0e-14, "log-energy interpolation");
  CheckThrows([&] { (void)shield_sim::InterpolateResponseInLogEnergy(
      responseEnergy, response, 0.99); }, "response extrapolation rejection");

  const auto fold = shield_sim::FoldPiecewiseConstantSpectrum(
      {1.0, 4.0}, {2.0}, responseEnergy, response, 1.0e-6, 16);
  const double analyticFold = 2.0 * (3.0 * (4.0 - 1.0) +
      0.5 * ((4.0 * std::log(4.0) - 4.0) - (1.0 * std::log(1.0) - 1.0)));
  Check(fold.converged, "adaptive spectrum fold converged");
  Check(fold.finalRelativeChange <= 1.0e-6, "adaptive fold satisfies requested refinement");
  CheckNear(fold.value, analyticFold, 3.0e-6, "adaptive fold analytic integral");
  CheckThrows([&] { (void)shield_sim::FoldPiecewiseConstantSpectrum(
      {0.5, 4.0}, {1.0}, responseEnergy, response, 0.02, 8); },
      "spectrum fold extrapolation rejection");

  const auto allocation = shield_sim::AllocateNeymanHistories(
      {{"low",0.5,1.0},{"rare",0.3,2.0},{"zero_pilot",0.2,0.0}}, 100, 10);
  Check(allocation.size() == 3 && allocation[0].histories == 42 &&
        allocation[1].histories == 48 && allocation[2].histories == 10,
        "deterministic Neyman allocation with rare-stratum minimum");
  CheckThrows([&] { (void)shield_sim::AllocateNeymanHistories(
      {{"a",1.0,1.0},{"b",1.0,1.0}}, 10, 6); }, "infeasible Neyman minimum rejection");

  Check(shield_sim::ParseUtcTimestamp("1970-01-01T00:00:00Z") == 0,
        "UTC Unix epoch");
  Check(shield_sim::ParseUtcTimestamp("2020-02-29T00:00:55Z") -
        shield_sim::ParseUtcTimestamp("2020-02-29T00:00:00Z") == 55,
        "UTC leap-day interval");
  CheckThrows([&] { (void)shield_sim::ParseUtcTimestamp("2020-02-29T00:00:00-06:00"); },
              "numeric time-zone rejection");
  CheckThrows([&] { (void)shield_sim::ParseUtcTimestamp("2019-02-29T00:00:00Z"); },
              "invalid Gregorian date rejection");

  const auto valid = ValidExposure();
  const auto integrated = shield_sim::IntegrateTimeResolvedExposure(valid);
  const auto& integratedRows = integrated.At("integrated_exposures").AsArray();
  Check(integratedRows.size() == 1, "one integrated exposure series");
  const auto& row = integratedRows[0].AsObject();
  CheckNear(row.at("value").AsNumber(), 85.0, 1.0e-14,
            "constant/piecewise/partial/gap exposure integral");
  CheckNear(row.at("standard_uncertainty").AsNumber(), std::sqrt(49.0625), 1.0e-14,
            "covariance propagation through interval durations");
  const auto& ensembles = row.at("sampled_ensembles").AsArray();
  CheckNear(ensembles[0].At("integrated_mean").AsNumber(), 85.0, 1.0e-14,
            "sampled ensemble integrated mean");
  CheckNear(ensembles[0].At("integrated_standard_deviation").AsNumber(), 8.5, 1.0e-14,
            "sampled ensemble integrated deviation");
  Check(integrated.At("configuration_segments").AsArray().size() == 2,
        "configuration changes remain distinct");
  Check(integrated.At("anomaly_timestamps").AsArray().size() == 1,
        "anomaly timestamp retained separately");

  auto implicitGap = valid;
  implicitGap.At("intervals").AsArray()[1].At("end_utc") = "2020-02-29T00:00:20Z";
  CheckThrows([&] { (void)shield_sim::IntegrateTimeResolvedExposure(implicitGap); },
              "implicit time gap rejection");
  auto badOrder = valid;
  badOrder.At("intervals").AsArray()[1].At("start_utc") = "2020-02-29T00:00:09Z";
  CheckThrows([&] { (void)shield_sim::IntegrateTimeResolvedExposure(badOrder); },
              "overlap/out-of-order interval rejection");
  auto badCovariance = valid;
  auto& matrix = badCovariance.At("uncertainty").At("covariance_components").AsArray()[0]
      .At("interval_rate_covariance").AsArray();
  matrix[0].AsArray()[1] = 1.0;
  matrix[1].AsArray()[0] = 1.0;
  CheckThrows([&] { (void)shield_sim::IntegrateTimeResolvedExposure(badCovariance); },
              "indefinite source covariance rejection");
  auto falseGcr = valid;
  falseGcr.At("source_qualification").At("environment") = "ams02_a1_proton_component";
  falseGcr.At("source_qualification").At("environment_completeness") = "complete_gcr";
  CheckThrows([&] { (void)shield_sim::IntegrateTimeResolvedExposure(falseGcr); },
              "incomplete AMS component cannot be labeled full GCR");
  auto falseGoes = valid;
  falseGoes.At("source_qualification").At("environment") =
      "goes16_sgps_directional_ingestion";
  falseGoes.At("source_qualification").At("environment_completeness") = "directional_sep";
  falseGoes.At("source_qualification").At("transport_admission") = "verification_test_only";
  CheckThrows([&] { (void)shield_sim::IntegrateTimeResolvedExposure(falseGoes); },
              "GOES transport admission without boundary mapping rejection");

  shield_sim::Json evidenceDocument(shield_sim::Json::Object{
    {"comparisons", shield_sim::Json(shield_sim::Json::Array{
      shield_sim::Json(shield_sim::Json::Object{{"actual",fold.value},{"expected",analyticFold},
        {"id","adaptive_log_energy_response_fold"},{"unit","arbitrary response-source product"}}),
      shield_sim::Json(shield_sim::Json::Object{{"actual",row.at("value")},{"expected",85.0},
        {"id","piecewise_time_integral_with_explicit_gap"},{"unit","proxy"}}),
      shield_sim::Json(shield_sim::Json::Object{{"actual",row.at("standard_uncertainty")},
        {"expected",std::sqrt(49.0625)},{"id","interval_covariance_propagation"},{"unit","proxy"}}),
      shield_sim::Json(shield_sim::Json::Object{{"actual",ensembles[0].At("integrated_standard_deviation")},
        {"expected",8.5},{"id","sampled_source_ensemble_propagation"},{"unit","proxy"}})
    })},
    {"quadrature_final_relative_change",fold.finalRelativeChange},
    {"schema_version","shieldSim.validation_evidence.v1"},
    {"test_case","stage7_core"},
    {"validation_status","stage7_software_framework_incomplete_physical_gcr_goes_qualification"}});
  std::ofstream output(evidence / "stage7_core.json", std::ios::binary | std::ios::trunc);
  output << evidenceDocument.DumpCanonical() << '\n';

  if (failures) {
    std::cerr << failures << " Stage-7 core test(s) failed\n";
    return 1;
  }
  std::cout << "All Stage-7 core tests passed\n";
  return 0;
}
