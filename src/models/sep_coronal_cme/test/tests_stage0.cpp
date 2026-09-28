#include "test_framework.h"

#include "sep_coronal_cme/model_configuration.h"

#include <map>
#include <string>
#include <vector>

namespace SCCMTest {
namespace {
using SEP::Core::StatusCode;

void CFG3D12() {
  const std::string good = Fixture();
  const auto parsed = ParseGood(good);
  Require(parsed.schema5.schemaVersion == 5, "schema-5 dispatch failed");
  RequireRejected(good + "\n[run]\nunknown_key = 1\n", "unknown key");
  RequireRejected(ReplaceOnce(good, "schema_version = 5",
                              "schema_version = 5\nschema_version = 5"),
                  "duplicate version");
  RequireRejected(ReplaceOnce(good, "end_time_s = 86400",
                              "end_time_s = 86400 s"),
                  "inline SI unit");
}

void CFG3D13() {
  for (int version = 1; version <= 4; ++version) {
    const std::string bytes = "# legacy bytes stay opaque\n[run]\n"
        "schema_version = " + std::to_string(version) +
        "\nlegacy application syntax = untouched\n";
    const auto result = ParseGood(bytes);
    Require(result.disposition ==
                SEP::CoronalCME::ParseDisposition::LegacyPassThrough,
            "legacy deck was reinterpreted");
    Require(result.legacyBytes == bytes, "legacy bytes changed");
    Require(!result.fingerprint.empty(), "legacy fingerprint missing");
  }
}

void CFG3D14() {
  const std::string good = Fixture();
  RequireRejected(ReplaceOnce(good, "outer_radius_m = 2.0871e10",
                              "outer_radius_m = 1.0e9"),
                  "radius ordering");
  RequireRejected(ReplaceOnce(good, "adiabatic_index = 1.6666666666666667",
                              "adiabatic_index = 1"),
                  "adiabatic EOS");
  RequireRejected(ReplaceOnce(good, "interface_radius_m = 1.53054e9",
                              "interface_radius_m = 3e9"),
                  "SCS geometry");
}

void CFG3D15() {
  const std::string good = Fixture();
  std::string noSheet = ReplaceOnce(good, "model = finite-shell-schatten",
                                    "model = none");
  RequireRejected(noSheet, "inactive SCS values are nonzero");
  RequireRejected(ReplaceOnce(good, "closed_polytropic_index = 0",
                              "closed_polytropic_index = 1.1"),
                  "inactive closed gamma");
}

void CFG3D16() {
  const auto parsed = ParseGood(Fixture());
  const auto& values = parsed.schema5.assignments;
  const std::vector<std::string> authorities = {
      "run.transport", "pfss.coefficients_file", "solar_wind.model",
      "open_closed_interface.policy", "plasma_eos.electron_mass_in_density"};
  for (const std::string& authority : authorities)
    Require(values.count(authority) == 1u, "missing selector: " + authority);
  RequireRejected(ReplaceOnce(Fixture(), "transport = ballistic-verification",
                              "transport = invented-fallback"),
                  "unknown transport must not fall back");
}

void CFG3D17() {
  const auto configuration = ParseGood(Fixture()).schema5;
  const auto stage0 = SEP::CoronalCME::CheckCapabilityAvailability(configuration, 0);
  const auto stage1 = SEP::CoronalCME::CheckCapabilityAvailability(configuration, 1);
  const auto stage2 = SEP::CoronalCME::CheckCapabilityAvailability(configuration, 2);
  Require(stage0.code == StatusCode::NotImplemented, "Stage-1 PFSS activated early");
  Require(stage1.code == StatusCode::NotImplemented, "Stage-2 wind activated early");
  Require(stage2.code == StatusCode::NotImplemented, "Stage-3 SCS fell back");
}

void CFG3D18() {
  const auto base = ParseGood(Fixture());
  const auto changed = ParseGood(ReplaceOnce(Fixture(),
      "electron_mass_in_density = neglect-recorded", "electron_mass_in_density = include"));
  Require(base.fingerprint != changed.fingerprint,
          "electron mass discriminant omitted from identity");
}

void CFG3D19() {
  const std::string good = Fixture();
  RequireRejected(ReplaceOnce(good, "input_rate_convention = sidereal",
                              "input_rate_convention = synodic"),
                  "synodic rate without ephemeris");
  RequireRejected(ReplaceOnce(good, "differential_rotation_coefficients_file = none",
                              "differential_rotation_coefficients_file = active.dat"),
                  "duplicate rigid/differential authority");
}

void CFG3D20() {
  const std::string good = Fixture();
  RequireRejected(ReplaceOnce(good,
      "maximum_outer_zonal_nonmonopole_power_fraction = 0",
      "maximum_outer_zonal_nonmonopole_power_fraction = 0.1"),
      "diagnostic branch with active bounds");
  std::string production = ReplaceOnce(good, "intent = analytic-verification",
                                       "intent = production-shock-injection");
  production = ReplaceOnce(production, "radialization_gate = diagnostic-only",
                            "radialization_gate = outer-zonal-power-and-latitude-flatness");
  production = ReplaceOnce(production,
      "maximum_outer_zonal_nonmonopole_power_fraction = 0",
      "maximum_outer_zonal_nonmonopole_power_fraction = 0.02");
  production = ReplaceOnce(production,
      "latitude_minimum_unmasked_longitude_fraction = 0",
      "latitude_minimum_unmasked_longitude_fraction = 0.8");
  production = ReplaceOnce(production,
      "maximum_unsigned_radial_flux_rms_fraction = 0",
      "maximum_unsigned_radial_flux_rms_fraction = 0.1");
  production = ReplaceOnce(production,
      "maximum_unsigned_radial_flux_p95_to_p05_ratio = 0",
      "maximum_unsigned_radial_flux_p95_to_p05_ratio = 2");
  ParseGood(production);
}

void RST3D04() {
  std::map<std::string, std::string> base = {
      {"pfss.coefficients_file", "a.dat"},
      {"pfss.coefficients_file.content_checksum", "sha-a"},
      {"open_closed_interface.policy", "diagnostic-kinematic"},
      {"comment", "first"}};
  const std::string fingerprint = SEP::CoronalCME::ComputePhysicsFingerprint(base);
  base["pfss.coefficients_file"] = "/different/path/a.dat";
  base["comment"] = "second";
  Require(fingerprint == SEP::CoronalCME::ComputePhysicsFingerprint(base),
          "path/comment changed physics identity");
  base["pfss.coefficients_file.content_checksum"] = "sha-b";
  Require(fingerprint != SEP::CoronalCME::ComputePhysicsFingerprint(base),
          "asset checksum did not change identity");
}

void RST3D05() {
  const auto base = ParseGood(Fixture());
  const std::vector<std::pair<std::string, std::string>> mutations = {
      {"campaign_seed_u64 = 731993", "campaign_seed_u64 = 731994"},
      {"outer_radius_m = 2.0871e10", "outer_radius_m = 2.1e10"},
      {"adiabatic_index = 1.6666666666666667", "adiabatic_index = 1.4"},
      {"maximum_degree = 8", "maximum_degree = 9"}};
  for (const auto& mutation : mutations) {
    const auto changed = ParseGood(ReplaceOnce(Fixture(), mutation.first,
                                               mutation.second));
    Require(base.fingerprint != changed.fingerprint,
            "restart category absent from fingerprint: " + mutation.first);
    Require(!SEP::CoronalCME::CheckRestartIdentity(
                 base.fingerprint, changed.fingerprint, mutation.first).ok(),
            "restart mismatch accepted");
  }
}

}  // namespace

void RegisterStage0(Registry* tests) {
  (*tests)["CFG3D12"] = CFG3D12;
  (*tests)["CFG3D13"] = CFG3D13;
  (*tests)["CFG3D14"] = CFG3D14;
  (*tests)["CFG3D15"] = CFG3D15;
  (*tests)["CFG3D16"] = CFG3D16;
  (*tests)["CFG3D17"] = CFG3D17;
  (*tests)["CFG3D18"] = CFG3D18;
  (*tests)["CFG3D19"] = CFG3D19;
  (*tests)["CFG3D20"] = CFG3D20;
  (*tests)["RST3D04"] = RST3D04;
  (*tests)["RST3D05"] = RST3D05;
}

}  // namespace SCCMTest
