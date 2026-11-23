#include "sep_stage13_validation.h"

#include "amps_input.h"
#include "swcme_constants.hpp"

#include <cmath>
#include <filesystem>
#include <iomanip>
#include <locale>
#include <ostream>
#include <set>
#include <sstream>

namespace SEP::Stage13 {
namespace {

namespace Input = AMPS::Input;
namespace fs = std::filesystem;

constexpr const char* Profile = "swcme3d-sphere-shock-only";
constexpr double InnerRadiusSolar = 1.05;
constexpr double OuterRadiusAu = 5.0;
constexpr const char* OuterBoundary = "escape-no-incoming-sep";
constexpr const char* ShockReturn = "delete";

const std::set<std::string> roles = {
    "coronal-density", "early-shock-formation", "near-earth-spectra",
    "sunward-return", "outer-transport", "larger-domain-control"};

bool ReadNumber(const Input::Assignment& item, double& value) {
  // Numbers in this deck are dimensionless multiples of the named IAU unit.
  // A simulator or the report must never interpret them as metres or arc length.
  std::istringstream input(item.value);
  input.imbue(std::locale::classic());
  return bool(input >> value) && input.eof() && std::isfinite(value);
}

Input::Status CheckSelected(const Input::Block& block) {
  auto status = Input::RequireOnlyBlocks(block, {}, "the selected SEP validation case");
  if (!status.ok()) return status;

  const std::set<std::string> required = {
      "schema", "profile", "domain.inner_radius_solar", "domain.outer_radius_au",
      "boundary.outer", "shock.return", "case.role", "case.manifest"};
  const std::set<std::string> optional = {
      "case.continuity_reference", "case.control_outer_radius_au"};
  for (const auto& item : block.assignments) {
    if (required.count(item.key) == 0 && optional.count(item.key) == 0)
      return Input::ErrorAt(item.where, "unknown Stage 13 key '" + item.key + "'");
  }
  for (const auto& key : required) {
    const auto* item = block.Find(key);
    if (item == nullptr || item->value.empty())
      return Input::ErrorAt(block.begin, "missing Stage 13 key '" + key + "'");
  }

  const auto require = [&](const char* key, const char* expected) {
    const auto& item = *block.Find(key);
    return item.value == expected ? Input::Status::OK() :
        Input::ErrorAt(item.where, "Stage 13 key '" + item.key + "' requires '" + expected + "'");
  };
  for (const auto& check : {require("schema", "1"), require("profile", Profile),
                            require("boundary.outer", OuterBoundary),
                            require("shock.return", ShockReturn)})
    if (!check.ok()) return check;

  double inner = 0.0, outer = 0.0;
  const auto& innerItem = *block.Find("domain.inner_radius_solar");
  const auto& outerItem = *block.Find("domain.outer_radius_au");
  if (!ReadNumber(innerItem, inner) || inner != InnerRadiusSolar)
    return Input::ErrorAt(innerItem.where, "Stage 13 inner particle boundary must be 1.05 solar radii");
  if (!ReadNumber(outerItem, outer) || outer != OuterRadiusAu)
    return Input::ErrorAt(outerItem.where, "Stage 13 outer particle boundary must be 5 AU");

  const auto& role = *block.Find("case.role");
  if (roles.count(role.value) == 0)
    return Input::ErrorAt(role.where, "unknown Stage 13 comparison role '" + role.value + "'");
  const auto* control = block.Find("case.control_outer_radius_au");
  if (role.value == "larger-domain-control" && control == nullptr)
    return Input::ErrorAt(role.where, "a larger-domain control requires an explicit control outer radius");
  if (control != nullptr) {
    double controlRadius = 0.0;
    if (role.value != "larger-domain-control" ||
        !ReadNumber(*control, controlRadius) || controlRadius <= OuterRadiusAu)
      return Input::ErrorAt(control->where, "control radius requires the larger-domain-control role and a value above 5 AU");
  }
  const auto* continuity = block.Find("case.continuity_reference");
  if (continuity != nullptr && (role.value != "near-earth-spectra" || continuity->value != "OV01"))
    return Input::ErrorAt(continuity->where, "only the near-Earth spectra case can select the OV01 continuity reference");
  return Input::Status::OK();
}

}  // namespace

int RunValidationDeck(const std::string& path, const std::string& caseName,
                      bool listOnly, bool validateOnly, std::ostream& out,
                      std::ostream& error) {
  Input::Document document;
  auto status = Input::ParseFile(path, &document);
  if (status.ok()) status = Input::RequireNoAssignments(document.root);
  const auto* application = document.root.FindBlock("application", "sep-validation");
  if (status.ok() && (application == nullptr || document.root.blocks.size() != 1))
    status = Input::Status::Error("requires exactly one '#begin application: sep-validation' block");
  if (status.ok()) status = Input::RequireOnlyBlocks(*application, {"case"}, "SEP validation");
  if (status.ok() && application->blocks.empty())
    status = Input::ErrorAt(application->begin, "the SEP validation deck requires named case blocks");
  Input::Block selected;
  if (status.ok() && listOnly) {
    // Resolve each name using the same generic overlay as selection. Listing
    // does not touch the observation files or require a configured model.
    for (const auto& entry : application->blocks) {
      status = Input::ResolveCase(*application, entry.name, &selected);
      if (status.ok()) status = CheckSelected(selected);
      if (!status.ok()) break;
    }
    if (status.ok()) {
      for (const auto& entry : application->blocks) out << entry.name << '\n';
      return 0;
    }
  } else if (status.ok()) {
    status = Input::ResolveCase(*application, caseName, &selected);
    if (status.ok()) status = CheckSelected(selected);
  }
  if (!status.ok()) {
    error << "ERROR: SEP validation runtime input: " << status.message << '\n';
    return 2;
  }

  // Includes preserve each record's file location. This matters when the
  // campaign is submitted from a separate working directory or batch stage.
  const auto& item = *selected.Find("case.manifest");
  std::error_code ec;
  const auto manifest = fs::absolute(fs::path(item.where.file).parent_path() / item.value, ec)
      .lexically_normal();
  if (ec || !fs::is_regular_file(manifest, ec) || ec) {
    error << "ERROR: SEP validation runtime input: " << item.where.Prefix()
          << "cannot read selected observation manifest " << manifest.string() << '\n';
    return 2;
  }

  // One line per field gives the campaign runner unambiguous provenance even
  // for absolute paths with spaces. Use the same IAU constants as SWCME.
  out << "stage13.case=" << caseName << '\n'
      << "stage13.profile=" << Profile << '\n'
      << "stage13.role=" << selected.Find("case.role")->value << '\n'
      << "stage13.manifest=" << manifest.string() << '\n'
      << std::setprecision(17)
      << "stage13.physical_inner_radius_m=" << InnerRadiusSolar * swcme::constants::SOLAR_RADIUS_M << '\n'
      << "stage13.physical_outer_radius_m=" << OuterRadiusAu * swcme::constants::AU_M << '\n'
      << "stage13.outer_boundary=" << OuterBoundary << '\n'
      << "stage13.shock_return=" << ShockReturn << '\n';
  if (const auto* control = selected.Find("case.control_outer_radius_au")) {
    double radius = 0.0;
    ReadNumber(*control, radius); // Already validated, in astronomical units.
    out << "stage13.control_outer_radius_m=" << radius * swcme::constants::AU_M << '\n';
  }
  if (const auto* continuity = selected.Find("case.continuity_reference"))
    out << "stage13.continuity_reference=" << continuity->value << '\n';

  if (validateOnly) {
    out << "stage13.runtime_input=PASS\n";
    return 0;
  }
  // Native validation must never fall through to ARC/Parker defaults, start
  // time stepping on an unqualified background, or reuse OV01 as evidence
  // for the new spherical model. Preserve a machine-readable terminal state.
  out << "stage13.comparison=SKIP\n"
      << "stage13.reason=Native spherical observation operators and reviewed event settings are required.\n";
  return 77;
}

}  // namespace SEP::Stage13
